
#define PY_SSIZE_T_CLEAN
#include <Python.h>

PyObject *func(PyObject *self, PyObject *args);

static PyMethodDef FuncMethods[] = {
  {"name", func, METH_VARARGS, "doc"},
  {NULL, NULL, 0, NULL}
};

static struct PyModuleDef funcsmodule = {
  PyModuleDef_HEAD_INIT,
  "funcs",
  "Contains several useful functions.",
  -1,
  FuncMethods
};

PyMODINIT_FUNC PyInit_funcs(void) {
  return PyModule_Create(&funcsmodule);
}
