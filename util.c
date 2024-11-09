#include <stdarg.h>
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
