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

PyTypeObject* PyIppAnchor_Type(nullptr);
PyTypeObject* PyIppCoords_Type(nullptr);
PyTypeObject* PyIppShortestPathEntry_Type(nullptr);

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
        auto const createPyAnchor = [](Ipp::PwalnEntry const& anchor) {
            PyObject* const pyAnchor(PyStructSequence_New(PyIppAnchor_Type));
            PyStructSequence_SetItem(pyAnchor, 0, PyLong_FromSize_t(anchor.refStart()));
            PyStructSequence_SetItem(pyAnchor, 1, PyLong_FromSize_t(anchor.refEnd()));
            PyStructSequence_SetItem(pyAnchor, 2, PyLong_FromSize_t(anchor.qryStart()));
            PyStructSequence_SetItem(pyAnchor, 3, PyLong_FromSize_t(anchor.qryEnd()));
            return pyAnchor;
        };
        auto const createPyCoords = [self](Ipp::Coords const& coords) {
            PyObject* const pyCoords(PyStructSequence_New(PyIppCoords_Type));
            std::string const chromName(self->ipp.chromName(coords.chrom));
            PyStructSequence_SetItem(pyCoords, 0, PyUnicode_FromString(chromName.c_str()));
            PyStructSequence_SetItem(pyCoords, 1, PyLong_FromSize_t(coords.loc));
            return pyCoords;
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
                PyObject* const pySpe(PyStructSequence_New(PyIppShortestPathEntry_Type));
                PyStructSequence_SetItem(pySpe, 0, PyUnicode_FromString(currentSpecies.c_str()));
                PyStructSequence_SetItem(pySpe, 1, PyFloat_FromDouble(current.score));
                PyStructSequence_SetItem(pySpe, 2, createPyCoords(current.coords));
                PyStructSequence_SetItem(pySpe, 3, createPyAnchor(current.anchors.upstream));
                PyStructSequence_SetItem(pySpe, 4, createPyAnchor(current.anchors.downstream));
                PyList_Append(multiShortestPath, pySpe);
                Py_DECREF(pySpe);
                currentSpecies = current.prevSpecies;
            }
        }

        // Reverse the shortest path list to have it in the right order.
        PyList_Reverse(multiShortestPath);

        // Call the callback function.
        bool const hasDirect(coordProjection.direct.has_value());
        PyObject* const argList(Py_BuildValue(
                "OdOOOO",
                createPyCoords(refCoord),
                hasDirect ? coordProjection.direct->score : 0.0d,
                hasDirect ? createPyCoords(coordProjection.direct->nextCoords) : Py_None,
                hasDirect ? createPyAnchor(coordProjection.direct->anchors.upstream) : Py_None,
                hasDirect ? createPyAnchor(coordProjection.direct->anchors.downstream) : Py_None,
                multiShortestPath));
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

static PyStructSequence_Field ippAnchorFields[] = {
    {"ref_start", nullptr},
    {"ref_end", nullptr},
    {"qry_start", nullptr},
    {"qry_end", nullptr},

    {nullptr, nullptr}             /* Sentinel */
};
static PyStructSequence_Desc ippAnchorTypeDesc = {
    "ipp.Anchor",                  /* name */
    nullptr,                       /* doc */
    ippAnchorFields,               /* fields */
    (sizeof(ippAnchorFields)/sizeof(ippAnchorFields[0]) - 1) /* n_in_sequence */
};

static PyStructSequence_Field ippCoordsFields[] = {
    {"chrom", nullptr},
    {"loc", nullptr},

    {nullptr, nullptr}             /* Sentinel */
};
static PyStructSequence_Desc ippCoordsTypeDesc = {
    "ipp.Coords",                  /* name */
    nullptr,                       /* doc */
    ippCoordsFields,               /* fields */
    (sizeof(ippCoordsFields)/sizeof(ippCoordsFields[0]) - 1) /* n_in_sequence */
};
static PyObject*
ippCoordsToStr(PyObject *self) {
    PyObject* const format(PyUnicode_FromString("%s:%u"));
    PyObject* const ret(PyUnicode_Format(format, self));
    Py_DECREF(format);
    return ret;
}

static PyStructSequence_Field ippShortestPathEntryFields[] = {
    {"species", "species name [string]"},
    {"score", "projection score [float]}"},
    {"coords", "coords [Ipp.Coords]"},
    {"up_anchor", "upstream anchor [Ipp.Anchor]"},
    {"down_anchor", "downstream anchor [Ipp.Anchor]"},

    {nullptr, nullptr}             /* Sentinel */
};
static PyStructSequence_Desc ippShortestPathEntryTypeDesc = {
    "ipp.ShortestPathEntry",       /* name */
    nullptr,                       /* doc */
    ippShortestPathEntryFields,    /* fields */
    (sizeof(ippShortestPathEntryFields)/sizeof(ippShortestPathEntryFields[0]) - 1) /* n_in_sequence */
};

PyMODINIT_FUNC
PyInit_ipp(void) {
    class RefAnchor {
        // Anchor that decrements all registered objects upon destruction.
    public:
        ~RefAnchor() {
            for (PyObject* const obj : objs_) {
                Py_DECREF(obj);
            }
        }
        void add(PyObject* obj) {
            if (obj) {
                objs_.push_back(obj);
            }
        }
        void release() {
            objs_.clear();
        }
    private:
        std::vector<PyObject*> objs_;
    };
    RefAnchor refAnchor;

    // Create the module.
    PyObject* const m(PyModule_Create(&ippModule));
    refAnchor.add(m);
    if (!m) {
        return nullptr;
    }

    // Create the ipp.Ipp type and add it to the module.
    PyObject* const PyIpp_Type(PyType_FromSpec(&ippTypeSpec));
    refAnchor.add(PyIpp_Type);
    if (!PyIpp_Type) {
        return nullptr;
    }
    if (PyModule_AddObject(m, "Ipp", PyIpp_Type) < 0) {
        return nullptr;
    }

    // Create the "ipp.Anchor" named tuple and add it to the module.
    PyIppAnchor_Type = PyStructSequence_NewType(&ippAnchorTypeDesc);
    if (!PyIppAnchor_Type) {
        return nullptr;
    }
    if (PyModule_AddType(m, PyIppAnchor_Type) < 0) {
        return nullptr;
    }

    // Create the "Ipp.Coords" named tuple and add it to the module.
    PyIppCoords_Type = PyStructSequence_NewType(&ippCoordsTypeDesc);
    if (!PyIppCoords_Type) {
        return nullptr;
    }
    PyIppCoords_Type->tp_str = ippCoordsToStr;
    if (PyModule_AddType(m, PyIppCoords_Type) < 0) {
        return nullptr;
    }

    // Create the "Ipp.ShortestPathEntry" named tuple and add it to the module.
    PyIppShortestPathEntry_Type =
        PyStructSequence_NewType(&ippShortestPathEntryTypeDesc);
    if (!PyIppShortestPathEntry_Type) {
        return nullptr;
    }
    if (PyModule_AddType(m, PyIppShortestPathEntry_Type) < 0) {
        return nullptr;
    }

    refAnchor.release();
    return m;
}

} // extern "C"
