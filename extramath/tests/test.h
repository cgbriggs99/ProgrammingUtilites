#ifndef __TEST_H__
#define __TEST_H__

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#define ASSERT_WARN(cond, warn_var) errno = 0; if(!cond || errno != 0) { fprintf(stderr, "Warning: assertion failed at line %d in file " __FILE__ ": " #cond "\n", __LINE__); perror(""); warn_var++; }
#define ASSERT_WARN_MSG(cond, warn_var, message, ...) errno = 0; if(!cond || errno != 0) { fprintf(stderr, "Warning: assertion failed at line %d in file " __FILE__ ": " #cond ": " message, __LINE__, ##__VA_ARGS__); perror(""); warn_var++; }

#define ASSERT(cond, err_var) errno = 0; if(!cond || errno != 0) { fprintf(stderr, "Error: assertion failed at line %d in file " __FILE__ ": " #cond "\n", __LINE__); perror(""); err_var++; }
#define ASSERT_MSG(cond, err_var, message, ...) errno = 0; if(!cond || errno != 0) { fprintf(stderr, "Error: assertion failed at line %d in file " __FILE__ ": " #cond ": " message, __LINE__, ##__VA_ARGS__); perror(""); err_var++; }

#define ASSERT_FATAL(cond, warn_var) errno = 0; if(!cond || errno != 0) { fprintf(stderr, "Fatal error: assertion failed at line %d in file " __FILE__ ": " #cond "\n", __LINE__); perror("");return -1;}
#define ASSERT_FATAL_MSG(cond, message, ...) errno = 0; if(!cond || errno != 0) { fprintf(stderr, "Fatal error: assertion failed at line %d in file " __FILE__ ": " #cond ": " message, __LINE__, ##__VA_ARGS__); perror(""); return -1; }

#define CONV 1e-3
#define NEAR(a, b) (fabs(a - b) <= CONV)
#define NEARF(a, b) (fabsf(a - b) <= CONV)
#define NEARL(a, b) (fabsl(a - b) <= CONV)
#define CNEAR(a, b) (cabs(a - b) <= CONV)
#define CNEARF(a, b) (cabsf(a - b) <= CONV)
#define CNEARL(a, b) (cabsl(a - b) <= CONV)

#endif
