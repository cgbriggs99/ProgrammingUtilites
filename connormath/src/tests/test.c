 
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include "test.h"
#include <string.h>

int __assert_wrap(int cond, int line, const char *file, const char *func,
		  const char *condst) {
  __info_data_head_t data = {line, strlen(file), strlen(func), cond,
			     strlen(condst), 0, errno};
  __output_info(data, file, func, condst);
  return (cond || data.errno);
}

int __assert_wrap_fatal(int cond, int line, const char *file, const char *func,
			const char *condst) {
  __info_data_head_t data = {line, strlen(file), strlen(func), cond,
			     strlen(condst), 1, errno};
  __output_info(data, file, func, condst);
  if(cond || data.errno) {
    exit(cond);
  }
  return (cond);
}

void __output_info(__info_data_head_t info, const char *file, const char *func,
		   const char *condst) {
  fwrite(&info, sizeof(__info_data_head_t), 1, __get_teststruct(0)->__info);
  fwrite(file, sizeof(char), strlen(file), __get_teststruct(0)->__info);
  fwrite(func, sizeof(char), strlen(func), __get_teststruct(0)->__info);
  fwrite(condst, sizeof(char), strlen(condst), __get_test_struct(0)->__info);
}

static void exit_func(void) {
  __get_test_struct(-1);
}

int setup_tests(int argc, const char **argv) {
  test_struct_t *ts;
  ts->__info = fopen(argv[1], "rw");
  ts->__stdout = fopen(argv[2], "rw");
  ts->__stderr = fopen(argv[3], "rw");
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
      fclose(ts->__info);
      fclose(ts->__stdout);
      fclose(ts->__stderr);
      free(ts);
    }
    ts = (test_struct_t *) (-1);
  }
  return (ts);
}
