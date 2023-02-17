#ifndef __TEST_H__
#define __TEST_H__

#include <stdio.h>
#include <stdlib.h>

#define ASSERT_WARN(cond, warn_var) if(!cond) { fprintf(stderr, "Warning: assertion failed at line %d in file " __FILE__ ": " #cond "\n", __LINE__); warn_var++; }
#define ASSERT_WARN_MSG(cond, warn_var, message, ...) if(!cond) { fprintf(stderr, "Warning: assertion failed at line %d in file " __FILE__ ": " #cond ": " message, __LINE__, ##__VA_ARGS__); warn_var++; }

#define ASSERT(cond, err_var) if(!cond) { fprintf(stderr, "Error: assertion failed at line %d in file " __FILE__ ": " #cond "\n", __LINE__); err_var++; }
#define ASSERT_MSG(cond, err_var, message, ...) if(!cond) { fprintf(stderr, "Error: assertion failed at line %d in file " __FILE__ ": " #cond ": " message, __LINE__, ##__VA_ARGS__); err_var++; }

#define ASSERT_FATAL(cond, warn_var) if(!cond) { fprintf(stderr, "Fatal error: assertion failed at line %d in file " __FILE__ ": " #cond "\n", __LINE__); return -1;}
#define ASSERT_FATAL_MSG(cond, message, ...) if(!cond) { fprintf(stderr, "Fatal error: assertion failed at line %d in file " __FILE__ ": " #cond ": " message, __LINE__, ##__VA_ARGS__); return -1; }

#define CONV 1e-3
#define NEAR(a, b) fabs(a - b) < CONV
#define NEARF(a, b) fabsf(a - b) < CONV
#define NEARL(a, b) fabsl(a - b) < CONV
#define CNEAR(a, b) cabs(a - b) < CONV
#define CNEARF(a, b) cabsf(a - b) < CONV
#define CNEARL(a, b) cabsl(a - b) < CONV

#endif
