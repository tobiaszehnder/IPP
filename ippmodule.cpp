/**
 * Defines the "ipp" Python module and translates the calls to the Ipp class.
 */
#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

#include "ipp.h"

namespace {

template <typename ...Args>
std::string
format(char const* fmt, Args&& ...args) {
    auto const len(std::snprintf(nullptr, 0, fmt, std::forward<Args>(args)...));

    std::string ret(len+1, '\0');
    std::sprintf(ret.data(), fmt, std::forward<Args>(args)...);
    return ret;
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
    unsigned nCores;
    PyObject* callback;
    int res(PyArg_ParseTuple(args,
                             "ssO!IO",
                             &refSpecies,
                             &qrySpecies,
                             &PyList_Type, &pyRefCoords,
                             &nCores,
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
            return format("%u:%u", anchor.refStart, anchor.refEnd);
        };
        auto const qryAnchorStr = [](Ipp::PwalnEntry const& anchor) {
            return format("%u:%u", anchor.qryStart, anchor.qryEnd);
        };

        // Backtrace the shortest path from the reference to the given target
        // species (in reversed order).
        PyObject* const multiShortestPath(PyList_New(0));
        std::string currentSpecies(qrySpecies);
        while (!currentSpecies.empty()) {
            Ipp::ShortestPathEntry const& current(
                coordProjection.multiShortestPath.at(currentSpecies));
            PyObject* const tuple(
                Py_BuildValue("sds(ss)(ss)",
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

    // Do the coord projection.
    try {
        self->ipp.projectCoords(refSpecies,
                                qrySpecies,
                                refCoords,
                                nCores,
                                onJobDone);
    } catch (std::exception const&) {
        return nullptr;
    }

    Py_RETURN_NONE;
}

static PyMethodDef ippMethods[] = {
    {"load_pwalns", (PyCFunction)ippLoadPwalns, METH_VARARGS, "Reads the chromosomes and pwalns from the given file"},
    {"load_genome_sizes", (PyCFunction)ippLoadGenomeSizes, METH_VARARGS, "Reads the genome sizes from the files in the given directory"},
    {"set_half_life_distance", (PyCFunction)ippSetHalfLifeDistance, METH_VARARGS, "Sets the half-life distance"},
    {"project_coords", (PyCFunction)ippProjectCoords, METH_VARARGS, ""},

    {nullptr, nullptr, 0, nullptr}        /* Sentinel */
};
static PyTypeObject PyIpp_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)

    .tp_name = "ipp.Ipp",
    .tp_basicsize = sizeof(PyIpp),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_methods = ippMethods,
    .tp_new = ippNew,
};

static struct PyModuleDef ippModule = {
    PyModuleDef_HEAD_INIT,
    "ipp",   /* name of module */
    nullptr, /* module documentation, may be NULL */
    -1,      /* size of per-interpreter state of the module,
                or -1 if the module keeps state in global variables. */
    nullptr
};

PyMODINIT_FUNC
PyInit_ipp(void) {
    if (PyType_Ready(&PyIpp_Type) < 0) {
        return nullptr;
    }

    PyObject* const m(PyModule_Create(&ippModule));
    if (!m) {
        return nullptr;
    }

    Py_INCREF(&PyIpp_Type);
    if (PyModule_AddObject(m, "Ipp", (PyObject*)&PyIpp_Type) < 0) {
        Py_DECREF(&PyIpp_Type);
        Py_DECREF(m);
        return nullptr;
    }

    return m;
}

} // extern "C"
