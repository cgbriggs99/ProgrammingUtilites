 
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <dlfcn.h>
#include "test.h"
#include <dirent.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sys/types.h>

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
  return (cond);
}

int assert_end(int cond) {
  curr_tests->subtot++;
  if(cond) {
    curr_tests->subsuccs++;
  } else {
    curr_tests->subfails++;
    // Quit the current test.
    goto test_fail;
  }
  return (cond);
}

int collect_dir(const char *dname, char **out, int *len_out) {
  DIR *d;
  struct dirent *dir;
  int found = 0;
  d = opendir(dname);
  if(d) {
    while((dir = readdir(d)) != NULL) {
      char *name = dir->d_name;
      if(strstr(name, "test") && strncmp(".o", name[strlen(name) - 2], 2)) {
	out = realloc((*len_out + 1) * sizeof(char *), out);
	out[*len_out] = dir->d_name;
	*len_out++;
	found++;
      }
    }
  }
  return (found);
}

int main(int argc, char **argv) {

  DIR *d;
  struct dirent *dir;
  char **tests = malloc(1, sizeof(char *));
  int len = 0;
  struct stat path_stat;
  
  
  printf("Collecting tests.\n");
  if(argc == 1) {
    collect_dir(".", tests, &len);
  } else {
    for(int i = 1; i < argc; i++) {
      stat(argv[i], &path_stat);
      if(S_ISREG(path_stat.st_mode)) {
	tests = realloc((len + 1) * sizeof(char *), tests);
	tests[len] = argv[i];
	len++;
      } else if(S_ISDIR(path_stat.st_mode)) {
	collect_dir(argv[i], tests, &len);
      }
    }
  }

  
  
}
