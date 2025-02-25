#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "argp.h"
#include "globals.h"
#include "mc.h"
#include "printer.h"
#include "util.h"

#define MAX_LINE_LENGTH 256

double measure_energy(params_t *, constants_t, int);
void input_reader(FILE *source, bond_t **b, lcoeff_t **l, int *c, int *n, int *nb, int *);

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

int countlines (char *filename) {
    // count the number of lines in the file called filename
    FILE *fp  = fopen(filename, "r");
    int ch    = 0;
    int lines = 0;

    if (fp == NULL)
        ;
    return 0;

    lines++;
    while ((ch = fgetc(fp)) != EOF) {
        if (ch == '\n') lines++;
    }
    fclose(fp);
    return lines;
}

static double *path_t  = NULL;
static double *path_hx = NULL;

int main (int argc, char *argv[]) {
    args_t args = args_default();
    args_parse(argc, argv, &args);

    int bool_use_path = 0;

    FILE *cfgfptr = fopen("path.txt", "r");
    if (cfgfptr == NULL) {
        /* perror("fopen"); */
        /* exit(0); */
    }
    char cfgline[MAX_LINE_LENGTH] = { 0 };
    double _path_t                = { 0.0 };
    double _path_hx               = { 0.0 };
    int ch = 0, lines = 0;
    while (cfgfptr != NULL && !feof(cfgfptr)) {
        ch = fgetc(cfgfptr);
        if (ch == '\n') { lines++; }
    }
    if (cfgfptr != NULL) fclose(cfgfptr);

    cfgfptr = fopen("path.txt", "r");
    if (lines > 0) args.tau = lines;
    /* printf("Line count: %d\n", lines); */
    int current_line = 0;
    while (cfgfptr != NULL && fgets(cfgline, MAX_LINE_LENGTH, cfgfptr)) {
        int count = sscanf(cfgline, "%lf %lf", &_path_t, &_path_hx);

        if (path_t == NULL && path_hx == NULL) {
            path_t  = (double *)malloc((args.tau + 10) * sizeof(double));
            path_hx = (double *)malloc((args.tau + 10) * sizeof(double));
            printf("Allocated memory\n");
        }

        path_t[current_line]  = _path_t;
        path_hx[current_line] = _path_hx;

        ++current_line;
    }
    if (cfgfptr != NULL) fclose(cfgfptr);

    if (lines > 0) printf("Path: %f\t%f\n", _path_t, _path_hx);

    /* for (int i = 0; i < args.tau; ++i) { */
    /*     printf("%f\t%f\n", path_t[i], path_hx[i]); */
    /* } */
    if (lines > 0) bool_use_path = 1;

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
        .lc_len   = 0, /* Number of linear coefficients */
    };

    /*
     * ni: Initialize count of not IDENT operator to zero,
     * hx, j, beta, pa1, pa2: Will be set later in for loop,
     * opers, bonds, spins: Will be set later after reading file
     */
    params_t params = params_init();

    int lc_len = 0;
    if (args.tri_l != 0) {
        const int l  = args.tri_l;
        const int n  = l * l;
        const int m  = n * 5000;
        params.bonds = (bond_t *)malloc(n * 3 * sizeof(bond_t));
        params.spins = (short *)malloc(n * sizeof(short));
        params.opers = (oper_t *)malloc(m * sizeof(oper_t));
        default_setup(params.bonds, params.spins, params.opers, l, n, m);
        constants.nb = n * 3;
        constants.n  = n;
        constants.m  = m;
    } else {
        FILE *fptr = fopen(args.input_file, "r");
        if (fptr == NULL) {
            perror("fopen");
            exit(0);
        }

        /* Returned pointer should be freed */
        lc_len = 0;
        input_reader(fptr, &params.bonds, &params.lcoeffs, &params.constant, &constants.n,
                     &constants.nb, &lc_len);
        fclose(fptr);
        constants.lc_len = lc_len;
        constants.m      = constants.n * 5000;
        params.opers     = (oper_t *)malloc(constants.m * sizeof(oper_t));
        params.spins     = (short *)malloc(constants.n * sizeof(short));
        for (int i = 0; i < constants.n; i++) {
            params.spins[i] = (double_r250() < 0.5 ? 1 : -1);
        }
        for (int i = 0; i < constants.m; i++) {
            params.opers[i].type = IDENT;
        }
    }
    /* default_setup(params.bonds, params.spins, params.opers, l, n, m); */
    if (args.spin_conf != NULL) {
        FILE *fptr = fopen(args.spin_conf, "r");
        if (fptr == NULL) {
            perror("fopen");
            exit(0);
        }
        char line[MAX_LINE_LENGTH] = { 0 };
        while (fgets(line, MAX_LINE_LENGTH, fptr)) {
            int idx = 0, spin = 0;
            int count = sscanf(line, "%d %d", &idx, &spin);
            if (count != 2) continue;
            params.spins[idx] = spin;
        }
        fclose(fptr);

        params.beta = 1.0 / init_t; /* Initialize beta */
        const double init_eng = measure_energy(&params, constants, 0);
        printf("Initial energy: %f\n", init_eng);
    }

    /* exit(0); */

    /* Warm up */
    /* for (int i = 0; i < 100; ++i) */
    /*     mc_sweep(); */

    double jsum = 0.0;
    for (int i = 0; i < constants.nb; ++i)
        jsum += (params.bonds[i].val < 0 ? -params.bonds[i].val : params.bonds[i].val);

    double lcsum = 0.0;
    for (int i = 0; i < lc_len; ++i) {
        lcsum += (params.lcoeffs[i].val < 0 ? -params.lcoeffs[i].val : params.lcoeffs[i].val);
    }

    for (int j = 0; j < tau; ++j) {
        params.hx          = init_hx * (1 - ((double)j / tau)) + final_hx * ((double)j / tau);
        double temperature = init_t * (1 - ((double)j / tau)) + final_t * ((double)j / tau);

        if (bool_use_path) { /* Use path.txt configuration */
            temperature = path_t[j];
            params.hx   = path_hx[j];
        }

        params.beta = 1.0 / temperature;

        const double hx = params.hx, beta = params.beta;
        const int nb = constants.nb, n = constants.n;

        /* params.pa1 = (2. * jsum + hx * (double)(n)) * beta; */
        /* params.pa2 = 2 * jsum / (2. * jsum + hx * (double)(n)); */
        params.pa1 = (2. * jsum + hx * (double)(n - lc_len)) * beta;
        params.pa2 = 2 * jsum / (2. * jsum + hx * (double)(n - lc_len));

        mc_sweep(&params, constants);
        measure_energy(&params, constants, j);
    }

    if (args.print_conf) {
        /* for (int i = 0; i < constants.n; ++i) { */
        /*     printf("%d\t%d\n", i, params.spins[i]); */
        /* } */
        FILE *fptr = fopen(args.input_file, "r");
        if (fptr == NULL) {
            perror("fopen");
            exit(0);
        }
        const double final_eng = measure_energy(&params, constants, constants.tau);
        print_spins(fptr, constants.n, params.spins, final_eng);
    }

    clean_up(params.spins, params.opers, params.bonds, NULL);
    return 0;
}

/* If the site is present in the lcoeffs, return the value */
static int check_lcoeffs (lcoeff_t *lcoeffs, size_t n, const int site) {
    for (int i = 0; i < n; ++i) {
        if (lcoeffs[i].site == site) return lcoeffs[i].val;
    }
    return -1;
}

double measure_energy (params_t *p, constants_t c, int stp) {
    oper_t *opers = p->opers;
    bond_t *bonds = p->bonds;
    short *spins  = p->spins;

    const int n = c.n, nb = c.nb, m = c.m, tau = c.tau;

    double eng = 0.0;
    for (int b = 0; b < nb; b++) {
        int i = bonds[b].site1;
        int j = bonds[b].site2;
        eng += (spins[i] * spins[j] * bonds[b].val); /* TODO: Add linear coefficients */
        if (check_lcoeffs(p->lcoeffs, c.lc_len, i) != -1) {
            eng += spins[i] * check_lcoeffs(p->lcoeffs, c.lc_len, i);
        }
    }

    /* printf("%f\t%f\t%f\n", p->beta, p->hx, eng); */
    printf("%f\t%f\t%f\n", 1 / p->beta, p->hx, eng);
    if (stp == tau) {
        /* Print the final state */
        for (int b = 0; b < nb; ++b) {
            int i = bonds[b].site1, j = bonds[b].site2;
            /* printf("%d(%d)\t%d(%d)\t%f\t%f\t%f\n", spins[i], i, spins[j], j, bonds[b].val, */
            /*        (spins[i] * spins[j] * bonds[b].val), eng); */
        }
    }

    return eng;
}

/* Read the input file and set the number of sites and bonds */
void input_reader (FILE *source, bond_t **b, lcoeff_t **l, int *c, int *n, int *nb, int *lc_len) {
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
        if (line[0] == '#') continue;
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
                linear_coeffs[lc].site = (int)values[0];
                linear_coeffs[lc].val  = values[1];
                ++lc;
                break;
            case 1: *c += values[0]; break;
            default: fprintf(stderr, "Invalid input file\n"); break;
        }
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

    *n      = max_site + 1; /* Number of sites */
    *nb     = k;            /* Number of bonds */
    *b      = bonds;
    *l      = linear_coeffs;
    *lc_len = lc;

    /* for (int i = 0; i < *lc_len; ++i) { */
    /*     printf("lc_len: %d\t%d\t%f\n", i, linear_coeffs[i].site, linear_coeffs[i].val); */
    /* } */

    /* for (int i = 0; i < k; ++i) { */
    /*     printf("%d\t\t%d\t%d\t%f\n", i, bonds[i].site1, bonds[i].site2, bonds[i].val); */
    /* } */
}
