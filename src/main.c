#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../include/argp.h"
#include "../include/color.h"
#include "globals.h"
#include "mc.h"
#include "printer.h"
#include "util.h"

#define MAX_LINE_LENGTH 256

double measure_energy(params_t *, constants_t, int);
void input_reader(FILE *source, bond_t **b, lcoeff_t **l, int *c, int *n, int *nb, int *);

option_t options[] = {
    {           "help",       no_argument,    0, 'h',      0,             "Display this help and exit" },
    {        "version",       no_argument,    0, 'v',      0,   "Display version information and exit" },
    {           "qubo",       no_argument,    0, 'q',      0,                        "Use QUBO format" },
    {           "file", required_argument,    0, 'f', "FILE",                      "Input hamiltonian" },
    {      "spin-conf", required_argument,    0, 'c', "FILE",             "Initialize spins from file" },
    {      "path-conf", required_argument,    0, 'a', "FILE",              "Initialize path from file" },
    {         "init-t", required_argument,    0, 't', "TEMP",                    "Initial temperature" },
    {        "final-t", required_argument,    0, 'T', "TEMP",                      "Final temperature" },
    {        "init-hx", required_argument,    0, 'x',   "HX",                          "Initial field" },
    {       "final-hx", required_argument,    0, 'X',   "HX",                            "Final field" },
    {            "tau", required_argument,    0, 's',  "TAU",                                    "Tau" },
    {     "print-conf",       no_argument,    0, 'p',      0, "Print final spin configuration to file" },
    { "print-progress",       no_argument,    0, 'P',      0,           "Print the annealing progress" },
    {             NULL,                 0, NULL,   0,   NULL,                                     NULL },
};

typedef struct {
    short qubo;
    char *file, *spin_conf, *path_conf;
    double init_t, final_t, init_hx, final_hx;
    int tau;
    int print_conf, print_progress;
} args_t;

/* load_path_conf: Load path configuration from file, return {path_t, path_hx, tau} */
static double **load_path_conf (char *filename) {
    /* Get the number of lines in the file */
    int *lines = (int *)malloc(sizeof(int)), tau = 0;
    if ((*lines = count_file_lines(filename)) > 0) tau = *lines;
    else {
        free(lines);
        return NULL;
    }

    /* Read the actual file as path configuration */
    double *path_t = (double *)malloc((tau + 10) * sizeof(double));
    if (path_t == NULL) {
        free(lines);
        return NULL;
    }
    double *path_hx = (double *)malloc((tau + 10) * sizeof(double));
    if (path_hx == NULL) {
        clean_up(path_t, lines, NULL);
        return NULL;
    }

    FILE *fptr = fopen(filename, "r");

    int current_line = 0, count = 0;
    char cfgline[MAX_LINE_LENGTH] = { 0 };
    double _path_t = 0.0, _path_hx = 0.0;

    while (fptr != NULL && fgets(cfgline, MAX_LINE_LENGTH, fptr)) {
        count = sscanf(cfgline, "%lf %lf", &_path_t, &_path_hx);

        path_t[current_line]  = _path_t;
        path_hx[current_line] = _path_hx;
        ++current_line;
    }
    fclose(fptr);

    double **path = (double *[]) { path_t, path_hx, (double *)lines };
    return path;
}

int main (int argc, char *argv[]) {
    int opt = 0, opt_index = 0;
    int valid = 1;

    args_t args = {
        .qubo           = 0,
        .file           = NULL,
        .spin_conf      = NULL,
        .path_conf      = NULL,
        .init_t         = 0.5,
        .final_t        = 0.000001,
        .init_hx        = 2.0,
        .final_hx       = 0.0,
        .tau            = 2048,
        .print_conf     = 0,
        .print_progress = 0,
    };

    struct option *long_options = long_opts(options);
    char *short_options         = short_opts(options);
    char *help_message          = help_mes(options);
    char *version_message       = "0.1.0";

    while ((opt = getopt_long(argc, argv, short_options, long_options, &opt_index)) != -1) {
        switch (opt) {
            case 'h': printf("%s\n", help_message); exit(1);
            case 'v': printf("Version %s\n", version_message); exit(1);
            case 'q': args.qubo = 1; break;
            case 'f': args.file = optarg; break;
            case 'c': args.spin_conf = optarg; break;
            case 'a': args.path_conf = optarg; break;
            case 't': args.init_t = atof(optarg); break;
            case 'T': args.final_t = atof(optarg); break;
            case 'x': args.init_hx = atof(optarg); break;
            case 'X': args.final_hx = atof(optarg); break;
            case 's': args.tau = atoi(optarg); break;
            case 'p': args.print_conf = 1; break;
            case 'P': args.print_progress = 1; break;
            case '?': valid = 0; break;
            default: printf("opt = %d\n", opt); break;
        }
    }
    if (!valid) {
        fprintf(stderr, RED "Invalid option\n");
        exit(1);
    }
    free_options(long_options, short_options, help_message);

    int bool_use_path = 0;
    double *path_t = NULL, *path_hx = NULL;
    if (args.path_conf != NULL) {
        double **path = load_path_conf(args.path_conf);
        if (path != NULL) {
            path_t = path[0], path_hx = path[1];
            args.tau = *((int *)path[2]);
            free(path[2]);
            bool_use_path = 1;
        }
    }
    if (bool_use_path)
        printf("Using path configuration: %s, path_t[%d] = %f, path_hx[%d] = %f\n", args.path_conf,
               args.tau, path_t[args.tau], args.tau, path_hx[args.tau]);

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
    if (args.file == NULL) {
        fprintf(stderr, RED "No input file provided\n");
        exit(1);
    }
    FILE *fptr = fopen(args.file, "r");
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
        printf("Initial energy: %f\n", measure_energy(&params, constants, 0));
    }

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
    clean_up(path_t, path_hx, NULL);

    if (args.print_conf) {
        /* for (int i = 0; i < constants.n; ++i) { */
        /*     printf("%d\t%d\n", i, params.spins[i]); */
        /* } */
        FILE *fptr = fopen(args.file, "r");
        if (fptr == NULL) {
            perror("fopen");
            exit(0);
        }
        const double final_eng = measure_energy(&params, constants, constants.tau);
        print_spins(fptr, constants.n, params.spins, final_eng, constants.init_t, constants.tau,
                    constants.init_hx);
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
