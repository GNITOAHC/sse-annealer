#include "util.h"
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

double double_r250 () {
    return (double)rand() / (double)RAND_MAX;
}

/* Free all the allocated memory, expecting the last pointer to be NULL */
void clean_up (void *ptr, ...) {
    va_list args;
    va_start(args, ptr);
    void *cur_ptr = ptr;

    while (cur_ptr != NULL) {
        free(cur_ptr);
        cur_ptr = va_arg(args, void *);
    }

    va_end(args);
    return;
}

int count_file_lines (char *filename) {
    /* count the number of lines in the file called filename */
    FILE *fp  = fopen(filename, "r");
    int ch    = 0;
    int lines = 0;

    if (fp == NULL) return -1;

    while ((ch = fgetc(fp)) != EOF) {
        if (ch == '\n') lines++;
    }
    fclose(fp);
    return lines;
}
