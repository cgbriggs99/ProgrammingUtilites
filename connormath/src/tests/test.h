#ifndef __CM_TEST_H__
#define __CM_TEST_H__

#include <stdio.h>
#include <errno.h>

typedef struct {
  FILE *__info;
  FILE *__stdout;
  FILE *__stderr;
} test_struct_t;

typedef struct {
  uint16_t line_length;
  uint16_t line;
  uint16_t filename_len;
  uint16_t funcname_len;
  int32_t status;
  int32_t errno;
} __info_data_head_t;

extern int setup_tests(int argc, const char **argv);

extern int end_tests(void);

// Redefine functions to be captured.
#define printf(fmt, ...) fprintf(__get_teststruct(0)->__stdout, fmt, __VA_ARGS__)

#define vprintf(fmt, arg) vfprintf(__get_teststruct(0)->__stdout, fmt, arg)

#define puts(str) fputs(str, __get_teststruct(0)->__stdout)

#define putchar(ch) fputc(ch, __get_teststruct(0)->__stdout)

#define putchar_unlocked(ch) fputc_unlocked(ch, __get_teststruct(0)->__stdout)

#define perror(str) fprintf(__get_teststruct(0)->__stderr, "%s: %s\n", str, strerror(errno))

#define stderr (__get_teststruct(0)->__stderr)
#define stdout (__get_teststruct(0)->__stdout)

#define assert(cond) __assert_wrap(cond, __LINE__, __FILE__, __func__)

#define assert_fatal(cond) __assert_wrap_fatal(cond, __LINE__, __FILE__, __func__)

extern int __assert_wrap(int cond, const char *line, const char *file, const char *func);

extern int __assert_wrap_fatal(int cond, const char *line, const char *file, const char *func);

extern test_struct_t *__get_teststruct(int alloc);

extern void __output_info(int status, int line_num, const char *file,
			  const char *func);

#endif
