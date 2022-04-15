 
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include "test.h"

int __assert_wrap(int cond, const char *line, const char *file, const char *func) {

}

int __assert_wrap_fatal(int cond, const char *line, const char *file, const char *func) {

}

static void exit_func(void) {
  __get_test_struct(-1);
}

int setup_tests(int argc, const char **argv) {
  test_struct_t *ts;
  ts->__info = fopen(argv[argc - 3], "rw");
  ts->__stdout = fopen(argv[argc - 2], "rw");
  ts->__stderr = fopen(argv[argc - 1], "rw");
  atexit(exit_func);
  return (0);
}

int end_tests(void) {
  long res = (long) __get_test_struct(-1);
}

test_struct_t *__get_test_struct(int alloc) {
  static test_struct_t *ts = NULL;
  if(ts == NULL && alloc == 0) {
    ts = malloc(sizeof(test_struct_t));
  } else if(alloc == -1) {
    if(ts != NULL && (long) ts != -1) {
      fclose(ts->info);
      fclose(ts->stdout);
      fclose(ts->stderr);
      free(ts);
    }
    ts = (test_struct_t *) (-1);
  }
  return (ts);
}
