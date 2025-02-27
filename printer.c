#include "defer.h"
#include "printer.h"
#include "stdlib.h"

static char *output_file (const int count, const double init_t, const int tau) {
    char *output = (char *)malloc(sizeof(char) * 50);
    if (output == NULL) {
        perror("malloc");
        return NULL;
    }
    int n;
    n = snprintf(output, 50, "conf_N%d_T%f_tau%d.dat", count, init_t, tau);
    return output;
}

void print_spins (FILE *fptr, const int n, const short *const spins, const double final_eng,
                  const double init_t, const int tau) {
    Deferral;
    /* Check the original input file (for skipping spins that aren't presenting) */
    double values[3] = { 0.0 }; /* 1 to 3 value(s) per line */
    char line[256]   = { 0 };
    short *spins_tmp = (short *)malloc(n * sizeof(short));
    Defer(free(spins_tmp));

    /* Spin values are either 1 or -1 (up or down, 0 means not present) */
    for (int i = 0; i < n; ++i)
        spins_tmp[i] = 0;

    while (fgets(line, 256, fptr)) {
        int count = sscanf(line, "%lf %lf %lf", &values[0], &values[1], &values[2]);

        if (line[0] == '#') continue; /* Skip comments */

        switch (count) {
            case 3:
                spins_tmp[(int)values[0]] = spins[(int)values[0]];
                spins_tmp[(int)values[1]] = spins[(int)values[1]];
                break;
            case 2: spins_tmp[(int)values[0]] = spins[(int)values[0]]; break;
            case 1: break;
            default: fprintf(stderr, "Invalid input file\n"); break;
        }
    }

    char *output = output_file(n, init_t, tau);
    if (output == NULL) {
        perror("malloc");
        exit(0);
    }
    Defer(free(output));
    FILE *fptr_out = fopen(output, "w");
    if (fptr_out == NULL) {
        perror("fopen");
        exit(0);
    }
    fprintf(fptr_out, "%f\n", final_eng);

    for (int i = 0; i < n; ++i) {
        /* if (spins_tmp[i] != 0) { printf("%d %d\n", i, spins_tmp[i]); } */
        if (spins_tmp[i] != 0) { fprintf(fptr_out, "%d %d\n", i, spins_tmp[i]); }
    }
    Return;
}
