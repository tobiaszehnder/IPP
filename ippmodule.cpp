/**
 * Defines the "ipp" Python module and translates the calls to the Ipp class.
 */
#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <csignal>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

#include "ipp.h"

namespace {

void signalHandler(int signal);
class AbortSignalHandler;

AbortSignalHandler* currentAbortSignalHandler(nullptr);

class AbortSignalHandler {
    // Registers signal handlers for SIGINT and SIGTERM and calls ipp->cancel()
    // upon receiving them.
    // The previous signal handlers are re-installed upon destruction.
public:
    explicit AbortSignalHandler(Ipp* ipp)
        : ipp_(ipp)
    {
        currentAbortSignalHandler = this;

        prevSigIntHandler_ = std::signal(SIGINT, ::signalHandler);
        prevSigTermHandler_ = std::signal(SIGTERM, ::signalHandler);
    }

    ~AbortSignalHandler() {
        std::signal(SIGTERM, prevSigTermHandler_);
        std::signal(SIGINT, prevSigIntHandler_);

        currentAbortSignalHandler = nullptr;
    }

    void signalHandler(int signal) {
        ipp_->cancel();
    }

private:
    Ipp* const ipp_;

    typedef void(*SigHandler)(int);
    SigHandler prevSigIntHandler_;
    SigHandler prevSigTermHandler_;
};

void
signalHandler(int signal) {
    currentAbortSignalHandler->signalHandler(signal);
}

} // namespace

extern "C" {

struct PyIpp {
    PyObject_HEAD

    Ipp ipp;
};

static PyObject*
ippNew(PyTypeObject* type, PyObject* args, PyObject* kwds) {
    auto const self(reinterpret_cast<PyIpp*>(type->tp_alloc(type, 0)));
    if (self != nullptr) {
        // In-place construct the Ipp instance.
        new (&self->ipp) Ipp();
    }
    return (PyObject *) self;
}

static void
ippDealloc(PyIpp* self) {
    // Destruct the ipp instance.
    self->ipp.~Ipp();

    Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyObject*
ippLoadPwalns(PyIpp* self, PyObject* args) {
    // Reads the pwalns from the given file.
    char const* fileName;
    if (!PyArg_ParseTuple(args,"s", &fileName)) {
        return nullptr;
    }

    try {
        self->ipp.loadPwalns(fileName);
    } catch (std::exception const& e) {
        PyErr_SetString(PyExc_RuntimeError, e.what());
        return nullptr;
    }

    Py_RETURN_NONE;
}

static PyObject*
ippLoadGenomeSizes(PyIpp* self, PyObject* args) {
    // Reads the genome sizes from the files in the given directory.
    char const* dirName;
    if (!PyArg_ParseTuple(args,"s", &dirName)) {
        return nullptr;
    }

    try {
        self->ipp.loadGenomeSizes(dirName);
    } catch (std::exception const& e) {
        PyErr_SetString(PyExc_RuntimeError, e.what());
        return nullptr;
    }

    Py_RETURN_NONE;
}

static PyObject*
ippSetHalfLifeDistance(PyIpp* self, PyObject* args) {
    // Sets the half-life distance.
    unsigned halfLifeDistance;
    if (!PyArg_ParseTuple(args,"I", &halfLifeDistance)) {
        return nullptr;
    }

    self->ipp.setHalfLifeDistance(halfLifeDistance);

    Py_RETURN_NONE;
}

static PyObject*
ippProjectCoords(PyIpp* self, PyObject* args) {
    char const* refSpecies;
    char const* qrySpecies;
    PyObject* pyRefCoords; // [ (chrom1, loc1), (chrom2, loc2), ... ]
    unsigned nThreads;
    PyObject* callback;
    int res(PyArg_ParseTuple(args,
                             "ssO!IO",
                             &refSpecies,
                             &qrySpecies,
                             &PyList_Type, &pyRefCoords,
                             &nThreads,
                             &callback));
    if (!res) {
        return nullptr;
    }
    if (!PyCallable_Check(callback)) {
        PyErr_SetString(PyExc_TypeError, "callback parameter must be callable");
        return nullptr;
    }

    // Parse the ref coords from the given list of tuples.
    std::vector<Ipp::Coords> refCoords;
    refCoords.reserve(PyList_Size(pyRefCoords));
    for (unsigned i(0); i < PyList_Size(pyRefCoords); ++i) {
        PyObject* const coordsEntry(PyList_GetItem(pyRefCoords, i));

        char const* refCoordsChromName;
        uint32_t refCoordsLoc;
        if (!PyArg_ParseTuple(coordsEntry,
                              "sI",
                              &refCoordsChromName,
                              &refCoordsLoc)) {
            return nullptr;
        }

        refCoords.emplace_back(self->ipp.chromIdFromName(refCoordsChromName),
                               refCoordsLoc);
    }

    auto const onJobDone = [&](Ipp::Coords const& refCoord,
                               Ipp::CoordProjection const& coordProjection) {
        // Call the given callback from the python script.
        auto const coordsStr = [self](Ipp::Coords const& coords) {
            return format("%s:%u",
                          self->ipp.chromName(coords.chrom).c_str(),
                          coords.loc);
        };
        auto const refAnchorStr = [](Ipp::PwalnEntry const& anchor) {
            return format("%u:%u", anchor.refStart(), anchor.refEnd());
        };
        auto const qryAnchorStr = [](Ipp::PwalnEntry const& anchor) {
            return format("%u:%u", anchor.qryStart(), anchor.qryEnd());
        };

        // Backtrace the shortest path from the reference to the given target
        // species (in reversed order).
        PyObject* const multiShortestPath(PyList_New(0));
        if (coordProjection.multiShortestPath.find(qrySpecies)
            != coordProjection.multiShortestPath.end()) {
            std::string currentSpecies(qrySpecies);
            while (!currentSpecies.empty()) {
                Ipp::ShortestPathEntry const& current(
                    coordProjection.multiShortestPath.at(currentSpecies));
                PyObject* const tuple(Py_BuildValue(
                        "sds(ss)(ss)",
                        currentSpecies.c_str(),
                        current.score,
                        coordsStr(current.coords).c_str(),
                        refAnchorStr(current.anchors.upstream).c_str(),
                        refAnchorStr(current.anchors.downstream).c_str(),
                        qryAnchorStr(current.anchors.upstream).c_str(),
                        qryAnchorStr(current.anchors.downstream).c_str()));
                PyList_Append(multiShortestPath, tuple);
                Py_DECREF(tuple);
                currentSpecies = current.prevSpecies;
            }
        }

        // Reverse the shortest path list to have it in the right order.
        PyList_Reverse(multiShortestPath);

        // Call the callback function.
        PyObject* argList;
        if (coordProjection.direct.has_value()) {
            argList = Py_BuildValue(
                "sds(ss)(ss)O",
                coordsStr(refCoord).c_str(),
                coordProjection.direct->score,
                coordsStr(coordProjection.direct->nextCoords).c_str(),
                refAnchorStr(coordProjection.direct->anchors.upstream).c_str(),
                refAnchorStr(coordProjection.direct->anchors.downstream).c_str(),
                qryAnchorStr(coordProjection.direct->anchors.upstream).c_str(),
                qryAnchorStr(coordProjection.direct->anchors.downstream).c_str(),
                multiShortestPath);
        } else {
            argList = Py_BuildValue("sds(ss)(ss)O",
                                    coordsStr(refCoord).c_str(),
                                    0.0d,
                                    nullptr,
                                    nullptr,
                                    nullptr,
                                    nullptr,
                                    nullptr,
                                    multiShortestPath);
        }
        PyObject* const result(PyObject_CallObject(callback, argList));
        Py_DECREF(argList);
        Py_DECREF(multiShortestPath);

        if (!result) {
            // An error occured.
            throw std::runtime_error("Error in the callback function");
        }
        Py_DECREF(result);
    };

    // Listen for Ctrl-C signals.
    AbortSignalHandler const abortSignalHandler(&self->ipp);

    // Do the coord projection.
    try {
        self->ipp.projectCoords(refSpecies,
                                qrySpecies,
                                refCoords,
                                nThreads,
                                onJobDone);
    } catch (std::exception const& e) {
        if (!PyErr_Occurred()) {
            PyErr_SetString(PyExc_RuntimeError, e.what());
        }
        return nullptr;
    }

    Py_RETURN_NONE;
}

static PyObject*
ippCancel(PyIpp* self, PyObject* args) {
    // Cancel ongoing project_coords() call.
    self->ipp.cancel();

    Py_RETURN_NONE;
}

static PyMethodDef ippMethods[] = {
    {"load_pwalns", (PyCFunction)ippLoadPwalns, METH_VARARGS, "Reads the chromosomes and pwalns from the given file"},
    {"load_genome_sizes", (PyCFunction)ippLoadGenomeSizes, METH_VARARGS, "Reads the genome sizes from the files in the given directory"},
    {"set_half_life_distance", (PyCFunction)ippSetHalfLifeDistance, METH_VARARGS, "Sets the half-life distance"},
    {"project_coords", (PyCFunction)ippProjectCoords, METH_VARARGS, ""},
    {"cancel", (PyCFunction)ippCancel, METH_VARARGS, "Cancel ongoing project_coords() call"},

    {nullptr, nullptr, 0, nullptr} /* Sentinel */
};
static PyType_Slot ippTypeSlots[] = {
    {Py_tp_new, (void*)ippNew},
    {Py_tp_dealloc, (void*)ippDealloc},
    {Py_tp_methods, (void*)ippMethods},

    {0, nullptr}                   /* Sentinel */
};

static PyType_Spec ippTypeSpec = {
    "ipp.Ipp",                     /* name */
    sizeof(PyIpp),                 /* basicsize */
    0,                             /* itemsize */
    Py_TPFLAGS_DEFAULT,            /* flags */
    ippTypeSlots                   /* slots */
};

static struct PyModuleDef ippModule = {
    PyModuleDef_HEAD_INIT,
    "ipp",                         /* m_name */
    nullptr,                       /* m_doc */
    0                              /* m_size */
};

PyMODINIT_FUNC
PyInit_ipp(void) {
    PyObject* const PyIpp_Type(PyType_FromSpec(&ippTypeSpec));
    if (!PyIpp_Type) {
        return nullptr;
    }

    PyObject* const m(PyModule_Create(&ippModule));
    if (!m) {
        return nullptr;
    }

    Py_INCREF(PyIpp_Type);
    if (PyModule_AddObject(m, "Ipp", PyIpp_Type) < 0) {
        Py_DECREF(PyIpp_Type);
        Py_DECREF(m);
        return nullptr;
    }

    return m;
}

} // extern "C"
