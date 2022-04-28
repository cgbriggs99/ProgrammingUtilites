
#define PY_SSIZE_T_CLEAN
#include <Python.h>

PyObject *func(PyObject *self, PyObject *args);

static PyMethodDef RootsMethods[] = {
  {"name", func, METH_VARARGS, "doc"},
  {NULL, NULL, 0, NULL}
};

static struct PyModuleDef rootsmodule = {
  PyModuleDef_HEAD_INIT,
  "roots",
  "Contains root finding algorithms.",
  -1,
  RootsMethods
};

PyMODINIT_FUNC PyInit_roots(void) {
  return PyModule_Create(&rootsmodule);
}
