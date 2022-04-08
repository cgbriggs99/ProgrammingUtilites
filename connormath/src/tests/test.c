 
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <dlfcn.h>
#include "test.h"

static test_struct_t *curr_tests;

test_struct_t setup_tests(test_func_t *funcs) {
  test_struct_t out;
  out.funcs = funcs;
  out.succs = 0;
  out.fails = 0;
  out.total = 0;
  out.subsuccs = 0;
  out.subfails = 0;
  out.subtotal = 0;
  return (out);
}

int assert_cont(int cond) {
  if(cond) {
    curr_tests->subsuccs++;
  } else {
    curr_tests->subfails++;
  }
  curr_tests->subtot++;
