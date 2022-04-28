
#define PY_SSIZE_T_CLEAN
#include <Python.h>

PyObject *func(PyObject *self, PyObject *args);

static PyMethodDef IntsMethods[] = {
  {"name", func, METH_VARARGS, "doc"},
  {NULL, NULL, 0, NULL}
};

static struct PyModuleDef intsmodule = {
  PyModuleDef_HEAD_INIT,
  "ints",
  "Various methods for solving integrals.",
  -1,
  IntsMethods
};

PyMODINIT_FUNC PyInit_ints(void) {
  return PyModule_Create(&intsmodule);
}
