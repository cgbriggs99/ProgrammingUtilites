#include "test.h"

int get_line_at(FILE *fp, int line, int maxlen, char *out) {
  int currline = 1;
  fseek(fp, 0, SEEK_SET);
  while(!feof(fp) && currline != line) {
    int ch = fgetc(fp);
    if(ch == '\n') {
      currline++;
    }
  }
  if(currline == line) {
    fgets(out, maxlen, fp);
    return (strlen(out));
  } else {
    return (-1);
  }
}

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

int parse_infofile(const char *fname, int *pass, int *fail, int *warn) {
  FILE *fp = fopen(fname, "r");
  __info_data_head_t data;

  while(!feof(fp)) {
    // Get the data packet.
    fread(&data, sizeof(__info_data_head_t), 1, fp);
    if(data.cond == 0) {
      *pass++;
    } else {
      switch(status) {
      case 0:
	*warn++;
	break;
      case 1:
	*fail++;
	break;
      }
    }
    // Skip the error stuff.
    fseek(fp, data.funcname_len + data.filename_len + data.cond_len, SEEK_CUR);
  }
  fclose(fp);
  return (0);
}

int show_errors(const char *fname) {
  FILE *fp = fopen(fname, "r");
  __info_data_head_t data;
  char *filebuf = calloc(BUFSIZ, sizeof(char)),
    *funcbuf = calloc(BUFSIZ, sizeof(char)),
    *condbuf = calloc(BUFSIZ, sizeof(char)),
    *linebuf = calloc(BUFSIZ, sizeof(char));

  while(!feof(fp)) {
    // Read data stuff.
    fread(&data, sizeof(__info_data_head_t), 1, fp);
    if(data.cond != 0 || errno != 0) {
      // Clear memory.
      memset(filebuf, 0, BUFSIZ);
      memset(funcbuf, 0, BUFSIZ);
      memset(condbuf, 0, BUFSIZ);
      memset(linebuf, 0, BUFSIZ);
      // Get debug info.
      fread(filebuf, sizeof(char), data.filename_len, fp);
      fread(funcbuf, sizeof(char), data.funcname_len, fp);
      fread(condbuf, sizeof(char), data.cond_len, fp);
      get_line_at(filebuf, data.line, BUFSIZ, linebuf);
      // Print debug info.
      switch(data.status) {
      case 0:
	printf("Warning at %s:%d in file %s: %s\n%s\n%s = %d\n", funcbuf,
	       data.line, filebuf, strerror(data.errno), linebuf, condbuf,
	       data.cond);
	break;
      case 1:
	printf("Error at %s:%d in file %s: %s\n%s\n%s = %d\n", funcbuf,
	       data.line, filebuf, strerror(data.errno), linebuf, condbuf,
	       data.cond);
	break;
      }
    }
    
  
  
#define INFOFILE ".info"
#define ERRFILE ".error"
#define OUTFILE ".out"

int main(int argc, char **argv) {
  DIR *d;
  struct dirent *dir;
  char **tests = malloc(1 * sizeof(char *));
  int len = 0;
  struct stat path_stat;
  int *passes, *fails, *totals, *passed, *warns;
  char *info = calloc(BUFSIZ, sizeof(char)),
    *out = calloc(BUFSIZ, sizeof(char)),
    *err = calloc(BUFSIZ, sizeof(char));
  char *pass_args[3];

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
    // Clear the buffers.
    memset(out, 0, BUFSIZ);
    memset(info, 0, BUFSIZ);
    memset(err, 0, BUFSIZ);
    // Set the buffers to the test string.
    strcpy(info, tests[i], strlen(tests[i]));
    strcpy(out, tests[i], strlen(tests[i]));
    strcpy(err, tests[i], strlen(tests[i]));
    // Mangle name to get outputs.
    strcat(info, INFOFILE);
    strcat(out, OUTFILE);
    strcat(err, ERRFILE);
    // Set up the arguments to pass.
    pass_args[0] = info;
    pass_args[1] = out;
    pass_args[2] = err;
    // Run the test.
    execv(test[i], pass_args);
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
