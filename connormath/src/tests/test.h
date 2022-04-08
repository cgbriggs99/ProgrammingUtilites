#ifndef __CM_TEST_H__
#define __CM_TEST_H__

typedef int (*test_func_t)(void);

// This will be defined in each file.
typedef struct {
  test_func_t *funcs;
  unsigned int succs;
  unsigned int fails;
  unsigned int total;
  unsigned int subsuccs;
  unsigned int subfails;
  unsigned int subtotal;
} test_struct_t;

extern test_struct_t setup_tests(const test_func_t *funcs);

extern int assert_cont(int cond);

extern int assert_end(int cond);

typedef test_struct_t (*init_tests_func_t)(void);

#endif
