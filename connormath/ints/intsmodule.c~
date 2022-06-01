
#define PY_SSIZE_T_CLEAN
#include <Python.h>

PyObject *func(PyObject *self, PyObject *args);

static PyMethodDef PythonTemplateMethods[] = {
  {"name", func, METH_VARARGS, "doc"},
  {NULL, NULL, 0, NULL}
};

static struct PyModuleDef pythontemplate = {
  PyModuleDef_HEAD_INIT,
  "name",
  "Doc",
  -1,
  PythonTemplateMethods
};

PyMODINIT_FUNC PyInit_{modulename}(void) {
  return PyModule_Create(&pythontemplate);
}
