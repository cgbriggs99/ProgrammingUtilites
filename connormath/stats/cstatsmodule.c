
#define PY_SSIZE_T_CLEAN
#include <Python.h>

PyObject *func(PyObject *self, PyObject *args);

static PyMethodDef CStatsMethods[] = {
  {"name", func, METH_VARARGS, "doc"},
  {NULL, NULL, 0, NULL}
};

static struct PyModuleDef cstatsmodule = {
  PyModuleDef_HEAD_INIT,
  "name",
  "Doc",
  -1,
  CStatsMethods
};

PyMODINIT_FUNC PyInit_cstats(void) {
  return PyModule_Create(&cstatsmodule);
}
