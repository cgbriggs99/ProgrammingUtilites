#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "../include/extramath.h"
#include <math.h>

PyObject *pyj0(PyObject *self, PyObject *args) {
  PyObject *z;
  PyArg_ParseTuple(args, "O", z);

  if(PyComplex_CheckExact(z)) {
    double _Complex res = cj0(CMPLX(PyComplex_RealAsDouble(z), PyComplex_ImagAsDouble(z)));
    return PyComplex_FromDoubles(creal(res), cimag(res));
  } else if(PyFloat_CheckExact(z)) {
    return PyFloat_FromDouble(j0(PyFloat_AsDouble(z)));
  } else if(PyLong_CheckExact(z)) {
    return PyFloat_FromDouble(j0(PyLong_AsDouble(z)));
  } else {
    PyErr_SetString(PyExc_TypeError, "j0 can only accept floats, ints, and complexes.");
    return NULL;
  } 
}

static PyMethodDef PythonTemplateMethods[] = {
  {"j0", pyj0, METH_VARARGS, ""},
  {NULL, NULL, 0, NULL}
};

static struct PyModuleDef pythontemplate = {
  PyModuleDef_HEAD_INIT,
  "extramath.bessels",
  "Doc",
  -1,
  PythonTemplateMethods
};

PyMODINIT_FUNC PyInit_{modulename}(void) {
  return PyModule_Create(&pythontemplate);
}
