#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

#include "alltypes.h"
#include "argp.h"
#include "mc.h"
#include "util.h"

#define MAX_LINE_LENGTH 256

/* TODO: Add one column & two column support */

void measure_energy(params_t *, constants_t, int);
void input_reader(FILE *source, bond_t *bonds, int *n, int *nb);

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
    FILE *fptr = fopen("./rec.dat", "r");
    if (fptr == NULL) {
        printf("Cannot open file \n");
        exit(0);
    }
    /* Personal setup (default arguments) */
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

    int l = 9;
    /* ************************** */

    double init_t   = args.init_t;
    double final_t  = args.final_t;
    double init_hx  = args.init_hx;
    double final_hx = args.final_hx;
    int tau         = args.tau;

    int n  = l * l;
    int nb = 3 * n;
    int m  = n * 5000;

    double J              = 1.0;
    constants_t constants = (constants_t) {
        .l        = l,
        .n        = n,   /* Number of sites */
        .nb       = nb,  /* Number of bonds */
        .m        = m,   /* Number of operators */
        .tau      = tau, /* Numer of sweeps steps */
        .final_hx = final_hx,
        .init_hx  = init_hx,
        .final_t  = final_t,
        .init_t   = init_t,
    };

    params_t params = (params_t) {
        .ni    = 0,   /* Initialize count of not IDENT operator to zero */
        .hx    = 0.0, /* Will be set later in for loop */
        .j     = 0.0, /* Will be set later in for loop */
        .beta  = 0.0, /* Will be set later in for loop */
        .pa1   = 0.0, /* Will be set later in for loop */
        .pa2   = 0.0, /* Will be set later in for loop */
        .opers = (oper_t *)malloc(500 * 5000 * sizeof(oper_t)),
        .bonds = (bond_t *)malloc(500 * sizeof(bond_t)),
        .spins = (short *)malloc(n * sizeof(short)),
    };

    input_reader(fptr, params.bonds, &constants.n, &constants.nb);
    fclose(fptr);
    constants.m = constants.n * 5000;
    for (int i = 0; i < constants.n; i++) {
        params.spins[i] = (double_r250() < 0.5 ? 1 : -1);
    }
    for (int i = 0; i < constants.m; i++) {
        params.opers[i].type = IDENT;
    }
    /* default_setup(params.bonds, params.spins, params.opers, l, n, m); */

    /* Warm up */
    /* for (int i = 0; i < 100; ++i) */
    /*     mc_sweep(); */

    double Jsum = 0.0;
    for (int i = 0; i < constants.nb; ++i)
        Jsum += (params.bonds[i].val < 0 ? -params.bonds[i].val : params.bonds[i].val);

    for (int j = 0; j <= tau; ++j) {
        params.hx          = init_hx * (1 - ((double)j / tau)) + final_hx * ((double)j / tau);
        double temperature = init_t * (1 - ((double)j / tau)) + final_t * ((double)j / tau);
        params.beta        = 1.0 / temperature;

        const double hx = params.hx, beta = params.beta;

        const int nb = constants.nb, n = constants.n;

        params.pa1 = (2. * Jsum + hx * (double)(n)) * beta;
        params.pa2 = 2 * Jsum / (2. * Jsum + hx * (double)(n));

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

    const int l = c.l, n = c.n, nb = c.nb, m = c.m, tau = c.tau;

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
void input_reader (FILE *source, bond_t *bonds, int *n, int *nb) {
    double values[3]           = { 0.0 }; // 1 to 3 value(s) per line
    char line[MAX_LINE_LENGTH] = { 0 };
    int k                      = 0;
    int max_site               = 0;

    for (int i = 0; i < *nb; ++i) {
        bonds[i].site1 = -1;
        bonds[i].site2 = -1;
        bonds[i].val   = 0.0;
    }

    while (fgets(line, MAX_LINE_LENGTH, source)) {
        int count      = sscanf(line, "%lf %lf %lf", &values[0], &values[1], &values[2]);
        int i          = (int)values[0];
        int j          = (int)values[1];
        double val     = values[2];
        bonds[k].site1 = i;
        bonds[k].site2 = j;
        bonds[k].val   = val;
        max_site       = (i > max_site ? i : max_site);
        max_site       = (j > max_site ? j : max_site);
        ++k;
    }

    *n  = max_site + 1; /* Number of sites */
    *nb = k;            /* Number of bonds */

    for (int i = 0; i < k; ++i) {
        printf("%d\t\t%d\t%d\t%f\n", i, bonds[i].site1, bonds[i].site2, bonds[i].val);
    }
}
