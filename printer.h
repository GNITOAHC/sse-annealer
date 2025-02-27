#pragma once
#include <stdio.h>

/* Print the spins to file according to the original input fptr */
void print_spins(FILE *fptr, const int n, const short *const spins, const double final_eng,
                 const double init_t, const int tau);
