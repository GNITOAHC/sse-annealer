#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

#include "argp.h"
#include "globals.h"
#include "mc.h"
#include "util.h"

#define MAX_LINE_LENGTH 256

/* TODO: Add one column & two column support */

void measure_energy(params_t *, constants_t, int);
void input_reader(FILE *source, bond_t **b, lcoeff_t **l, int *c, int *n, int *nb);

/* Deprecated */
void default_setup (bond_t *bonds, short *spins, oper_t *opers, int l, int n, int m) {
    for (int i = 0; i < n; i++) {
        bonds[i].site1 = i;
        bonds[i].site2 = (int)(i / l) * l + (i + 1) % l;
    }
    for (int i = n; i < n * 2; i++) {
        bonds[i].site1 = i % n;
        bonds[i].site2 = (bonds[i].site1 + l) % n;
    }
    for (int i = n * 2; i < n * 3; i++) {
        bonds[i].site1 = i % n;
        bonds[i].site2 = bonds[(i - l) % n].site2;
    }
    for (int i = 0; i < n * 3; i++) {
        bonds[i].val = 1.0;
    }
    for (int i = 0; i < n; i++) {
        spins[i] = (double_r250() < 0.5 ? 1 : -1);
    }
    for (int i = 0; i < m; i++) {
        opers[i].type = IDENT;
    }
}

int main (int argc, char *argv[]) {
    /* Default arguments */
    args_t args = (args_t) {
        .qubo       = 0,    /* Currently not supported */
        .input_file = NULL, /* Currently not supported */
        .init_t     = 0.5,
        .final_t    = 0.000001,
        .init_hx    = 2.0,
        .final_hx   = 0.0,
        .tau        = 4096,
    };
    args_parse(argc, argv, &args);

    double init_t   = args.init_t;
    double final_t  = args.final_t;
    double init_hx  = args.init_hx;
    double final_hx = args.final_hx;
    int tau         = args.tau;

    constants_t constants = (constants_t) {
        .n        = 0,   /* Number of sites */
        .nb       = 0,   /* Number of bonds */
        .m        = 0,   /* Number of operators */
        .tau      = tau, /* Numer of sweeps steps */
        .final_hx = final_hx,
        .init_hx  = init_hx,
        .final_t  = final_t,
        .init_t   = init_t,
    };

    /*
     * ni: Initialize count of not IDENT operator to zero,
     * hx, j, beta, pa1, pa2: Will be set later in for loop,
     * opers, bonds, spins: Will be set later after reading file
     */
    params_t params = params_init();

    FILE *fptr = fopen(args.input_file, "r");
    if (fptr == NULL) {
        perror("fopen");
        exit(0);
    }

    /* Returned pointer should be freed */
    input_reader(fptr, &params.bonds, &params.lcoeffs, &params.constant, &constants.n,
                 &constants.nb);
    fclose(fptr);
    constants.m  = constants.n * 5000;
    params.opers = (oper_t *)malloc(constants.m * sizeof(oper_t));
    params.spins = (short *)malloc(constants.n * sizeof(short));
    for (int i = 0; i < constants.n; i++) {
        params.spins[i] = (double_r250() < 0.5 ? 1 : -1);
    }
    for (int i = 0; i < constants.m; i++) {
        params.opers[i].type = IDENT;
    }
    /* default_setup(params.bonds, params.spins, params.opers, l, n, m); */

    /* exit(0); */

    /* Warm up */
    /* for (int i = 0; i < 100; ++i) */
    /*     mc_sweep(); */

    double jsum = 0.0;
    for (int i = 0; i < constants.nb; ++i)
        jsum += (params.bonds[i].val < 0 ? -params.bonds[i].val : params.bonds[i].val);

    for (int j = 0; j <= tau; ++j) {
        params.hx          = init_hx * (1 - ((double)j / tau)) + final_hx * ((double)j / tau);
        double temperature = init_t * (1 - ((double)j / tau)) + final_t * ((double)j / tau);
        params.beta        = 1.0 / temperature;

        const double hx = params.hx, beta = params.beta;

        const int nb = constants.nb, n = constants.n;

        params.pa1 = (2. * jsum + hx * (double)(n)) * beta;
        params.pa2 = 2 * jsum / (2. * jsum + hx * (double)(n));

        mc_sweep(&params, constants);
        measure_energy(&params, constants, j);
    }

    clean_up(params.spins, params.opers, params.bonds, NULL);
    return 0;
}

void measure_energy (params_t *p, constants_t c, int stp) {
    oper_t *opers = p->opers;
    bond_t *bonds = p->bonds;
    short *spins  = p->spins;

    const int n = c.n, nb = c.nb, m = c.m, tau = c.tau;

    double eng = 0.0;
    for (int b = 0; b < nb; b++) {
        int i = bonds[b].site1;
        int j = bonds[b].site2;
        eng += (spins[i] * spins[j] * bonds[b].val);
    }

    printf("%f\t%f\t%f\n", p->beta, p->hx, eng);
    if (stp == tau) {
        for (int b = 0; b < nb; ++b) {
            int i = bonds[b].site1, j = bonds[b].site2;
            printf("%d(%d)\t%d(%d)\t%f\t%f\t%f\n", spins[i], i, spins[j], j, bonds[b].val,
                   (spins[i] * spins[j] * bonds[b].val), eng);
        }
    }
}

/* Read the input file and set the number of sites and bonds */
void input_reader (FILE *source, bond_t **b, lcoeff_t **l, int *c, int *n, int *nb) {
    int capacity  = 500;
    bond_t *bonds = (bond_t *)malloc(capacity * sizeof(bond_t));

    int l_cap               = 500;
    lcoeff_t *linear_coeffs = (lcoeff_t *)malloc(l_cap * sizeof(lcoeff_t));

    double values[3]           = { 0.0 }; /* 1 to 3 value(s) per line */
    char line[MAX_LINE_LENGTH] = { 0 };
    int k                      = 0; /* Number of bonds */
    int max_site               = 0; /* Number of sites */
    int lc                     = 0; /* Number of linear coefficients */

    for (int i = 0; i < *nb; ++i) {
        bonds[i].site1 = -1;
        bonds[i].site2 = -1;
        bonds[i].val   = 0.0;
    }

    while (fgets(line, MAX_LINE_LENGTH, source)) {
        int count = sscanf(line, "%lf %lf %lf", &values[0], &values[1], &values[2]);
        printf("%d\n", count);
        int i, j;
        double val;
        switch (count) {
            case 3:
                bonds[k].site1 = (int)values[0];
                bonds[k].site2 = (int)values[1];
                bonds[k].val   = values[2];
                ++k;
                break;
            case 2:
                linear_coeffs[k].site = (int)values[0];
                linear_coeffs[k].val  = values[1];
                ++lc;
                break;
            case 1: *c += values[0]; break;
            default: fprintf(stderr, "Invalid input file\n"); break;
        }
        /* int i          = (int)values[0]; */
        /* int j          = (int)values[1]; */
        /* double val     = values[2]; */
        /* bonds[k].site1 = i; */
        /* bonds[k].site2 = j; */
        /* bonds[k].val   = val; */
        /* max_site       = (i > max_site ? i : max_site); */
        /* max_site       = (j > max_site ? j : max_site); */
        /* ++k; */
        if (k == capacity - 2) {
            capacity *= 2;
            bonds = (bond_t *)realloc(bonds, capacity * sizeof(bond_t));
        }
        if (lc == l_cap - 2) {
            l_cap         = l_cap * 2;
            linear_coeffs = (lcoeff_t *)realloc(linear_coeffs, l_cap * sizeof(lcoeff_t));
        }
    }

    for (int i = 0; i < k; ++i) {
        max_site = (bonds[i].site1 > max_site ? bonds[i].site1 : max_site);
        max_site = (bonds[i].site2 > max_site ? bonds[i].site2 : max_site);
    }

    *n  = max_site + 1; /* Number of sites */
    *nb = k;            /* Number of bonds */
    *b  = bonds;
    *l  = linear_coeffs;

    for (int i = 0; i < k; ++i) {
        printf("%d\t\t%d\t%d\t%f\n", i, bonds[i].site1, bonds[i].site2, bonds[i].val);
    }
}
