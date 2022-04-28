
#ifndef __CONNOR_TEST_H__
#define __CONNOR_TEST_H__

typedef struct {
  int passes;
  int fails;
  int runs;
  int warns;
} test_struct_t;

static test_func_t *__test_funcs = NULL;
static test_struct_t *__test_struct = NULL;

// Define here.
static test_struct_t init_tests();

static setup_tests(test_func_t *const tests) {
  __test_struct = malloc(sizeof(__test_struct));
  __test_funcs = tests;
  __test_struct->passes = 0;
  __test_struct->fails = 0;
  __test_struct->runs = 0;
  __test_struct->warns = 0;
  return (__test_struct);
}

#define assert(cond) 


#endif
