
#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include "../include/connormath.h"
#include "../include/roots.h"
#include <errno.h>

static double roots_func(double x, void *arg) {
  return (PyFloat_AsDouble(PyObject_Call((PyObject *) arg,
					 PyFloat_FromDouble(x), NULL)));
}

PyObject *py_brent_dekker(PyObject *self, PyObject *args) {
  double x0, x1, conv, tol;
  PyObject *call;

  PyErr_Clear();
  errno = 0;
  if(PyArg_ParseTuple(args, "Odddd", &call, &x0, &x1, &conv, &tol) != 1) {
    return (NULL);
  }
  
  return (PyFloat_FromDouble(brent_dekker(x0, x1, conv,
					  tol, roots_func, call)));
}

PyObject *py_newton(PyObject *self, PyObject *args) {
  double x0, conv;
  PyObject *call, *deriv;

  PyErr_Clear();
  errno = 0;
  if(PyArg_ParseTuple(args, "OOdd", &call, &deriv, &x0, &conv) != 1) {
    return (NULL);
  }

  return (PyFloat_FromDouble(newton(x0, conv, roots_func, call,
				    roots_func, deriv)));
}

PyObject *py_secant(PyObject *self, PyObject *args) {
  double x0, x1, conv;
  PyObject *call;

  PyErr_Clear();
  errno = 0;
  if(PyArg_ParseTuple(args, "Oddd", &call, &x0, &x1, &conv) != 1) {
    return (NULL);
  }

  return (PyFloat_FromDouble(secant(x0, x1, conv,
					  roots_func, call)));
}

PyObject *py_regula_falsi(PyObject *self, PyObject *args) {
  double x0, x1, conv;
  PyObject *call;

  PyErr_Clear();
  errno = 0;
  if(PyArg_ParseTuple(args, "Oddd", &call, &x0, &x1, &conv) != 1) {
    return (NULL);
  }

  return (PyFloat_FromDouble(regula_falsi(x0, x1, conv,
					  roots_func, call)));
}

#define ALLOC_SKIPS 8
PyObject *py_poly_roots(PyObject *self, PyObject *args) {
  PyObject *py_coefs;
  double *coefs, *out, conv;
  int len;

  PyErr_Clear();
  errno = 0;
  if(PyArg_ParseTuple(args, "Od", &py_coefs, &conv) != 1) {
    return (NULL);
  }

  PyObject *iter = PyObject_GetIter(py_coefs);
  if(iter == NULL) {
    return (NULL);
  }

  len = PyObject_LengthHint(py_coefs);
  if(len != -1) {
    coefs = calloc(len, sizeof(double));
  } else {
    coefs = malloc(0, sizeof(double));
  }
  int pos = 0;
  while(1) {
    PyObject *item = PyIter_Next(iter);
    if(item == NULL) {
      break;
    }
    if(pos >= len) {
      coefs = realloc(coefs, (len + ALLOC_SKIP) * sizeof(double));
      len += ALLOC_SKIP;
    }
    coefs[pos] = PyFloat_AsDouble(item);
    pos++;
  }
  if(pos < len) {
    coefs = realloc(coefs, pos * sizeof(double));
    len = pos;
  }
  out = calloc(len, sizeof(double));
  int out_len = poly_roots(coefs, len, out, conv);
  PyObject *outval = PyList_New(out_len);

  for(int i = 0; i < out_len; i++) {
    PyList_SetItem(outval, i, PyFloat_FromDouble(out[i]));
  }
  free(out);
  free(coefs);
  return (outval);
}

static PyMethodDef RootsMethods[] = {
  {"brent_dekker", py_brent_dekker, METH_VARARGS, "Takes arguments: "
   "brent_dekker(func, x0, x1, conv, tolerance"},
  {"newton", py_newton, METH_VARARGS, "Takes arguments: "
   "newton(func, der, x0, conv)"},
  {"secant", py_secant, METH_VARARGS, "Takes arguments: "
   "secant(func, x0, x1, conv)"},
  {"regula_falsi", py_regula_falsi, METH_VARARGS, "Takes arguments: "
   "regula_falsi(func, x0, x1, conv)"},
  {"poly_roots", py_poly_roots, METH_VARARGS, "Takes arguments: "
   "poly_roots(coefs, conv), where coefs starts from the ones place,"
   " then goes to the x place, then x squared, and so on."},
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
