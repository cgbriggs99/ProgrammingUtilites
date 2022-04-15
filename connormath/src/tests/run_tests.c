#include "test.h" 

int collect_dir(const char *dname, char **out, int *len_out) {
  DIR *d;
  struct dirent *dir;
  int found = 0;
  struct stat path_stat;
  d = opendir(dname);
  if(d) {
    while((dir = readdir(d)) != NULL) {
      char *name = dir->d_name;
      stat(name, &path_stat);
      if(strstr(name, "test") && strncmp(".o", &(name[strlen(name) - 2]), 2) == 0 && S_ISREG(path_stat.st_mode)) {
	char *pname = calloc(strlen(dname) + strlen(name) + 2, sizeof(char));
	pname = strncpy(pname, dname, strlen(dname));
	pname[strlen(dname)] = '/';
	pname = strncat(pname, name, strlen(name));
	out = realloc(out, (*len_out + 1) * sizeof(char *));
	out[*len_out] = pname;
	(*len_out)++;
	found++;
      }
    }
  }
  return (found);
}

int main(int argc, char **argv) {
  DIR *d;
  struct dirent *dir;
  char **tests = malloc(1 * sizeof(char *));
  int len = 0;
  struct stat path_stat;
  int *passes, *fails, *totals, *passed, *warns;

  curr_tests = calloc(1, sizeof(test_struct_t));
  
  printf("Collecting tests.\n");
  printf("%d\n", argc);
  if(argc == 1) {
    collect_dir(".", tests, &len);
  } else {
    for(int i = 1; i < argc; i++) {
      stat(argv[i], &path_stat);
      if(S_ISREG(path_stat.st_mode)) {
	tests = realloc(tests, (len + 1) * sizeof(char *));
	tests[len] = argv[i];
	len++;
      } else if(S_ISDIR(path_stat.st_mode)) {
	collect_dir(argv[i], tests, &len);
      }
    }
  }

  if(len == 1) {
    printf("Collected 1 test.\n");
  } else {
    printf("Collected %d tests.\n", len);
  }
  passes = calloc(len, sizeof(int));
  fails = calloc(len, sizeof(int));
  totals = calloc(len, sizeof(int));
  passed = calloc(len, sizeof(int));
  warns = calloc(len, sizeof(int));
  
  for(int i = 0; i < len; i++) {
    dlerror();
    void *handle = dlopen(tests[i], RTLD_NOW);
    char *err = dlerror();
    if(err != NULL) {
      fprintf(stderr, "Failed to open file: %s\n", err);
      continue;
    }
    init_tests_func_t init_func = dlsym(handle, "init_tests");
    if(err != NULL) {
      fprintf(stderr, "Failed to initialize tests: %s\n", err);
    }
    *curr_tests = init_func();
    

    for(int j = 0; curr_tests->funcs[j] != NULL; j++) {
      int ret = curr_tests->funcs[j]();
      if(ret == 0) {
	passed[j] = 1;
	curr_tests->succs++;
	curr_tests->total++;
      } else {
	passed[j] = 0;
	curr_tests->fails++;
	curr_tests->total++;
      }
    }
    passes[i] = curr_tests->succs;
    fails[i] = curr_tests->fails;
    totals[i] = curr_tests->total;
    warns[i] = curr_tests->subfails;
    dlclose(handle);
  }
  printf("Summary:\n");
  int pass = 0, fail = 0, total = 0, warn = 0;
  for(int i = 0; i < len; i++) {
    pass += passes[i];
    fail += fails[i];
    total += totals[i];
    warn += warns[i];
    if(warns[i] == 1) {
      printf("%s: %d passed, %d failed, 1 warning: ", tests[i], passes[i],
	     fails[i]);
    } else {
      printf("%s: %d passed, %d failed, %d warnings: ", tests[i], passes[i],
	     fails[i], warns[i]);
    }
    if(passed[i]) {
      printf("Success\n");
    } else {
      printf("Failed\n");
    }
    
  }
  if(pass == 1) {
    printf("1 pass, ");
  } else {
    printf("%d passes, ", pass);
  }
  if(fail == 1) {
    printf("1 fail, ");
  } else {
    printf("%d fails, ", fail);
  }
  if(warn == 1) {
    printf("1 warning\n");
  } else {
    printf("%d warnings\n", warn);
  }

  
  for(int i = 0; i < len; i++) {
    free(tests[i]);
  }
  free(passes);
  free(fails);
  free(totals);
  free(warns);
  free(passed);
  free(tests);
  return (0);
}
