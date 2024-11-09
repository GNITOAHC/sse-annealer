#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

#include "alltypes.h"
#include "argp.h"
#include "mc.h"
#include "util.h"

#define MAX_LINE_LENGTH 256

void measure_energy(params_t *, constants_t, int);

int main (int argc, char *argv[]) {
    /* Personal setup (default arguments) */
    args_t args = (args_t) {
        .qubo       = 0,    /* Currently not supported */
        .input_file = NULL, /* Currently not supported */
        .init_t     = 0.5,
        .final_t    = 0.000001,
        .init_hx    = 2.0,
        .final_hx   = 0.0,
        .tau        = 1024,
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
        .n        = n,
        .nb       = nb,
        .m        = m,
        .tau      = tau,
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
        .opers = (oper_t *)malloc(m * sizeof(oper_t)),
        .bonds = (bond_t *)malloc(nb * sizeof(bond_t)),
        .spins = (short *)malloc(n * sizeof(short)),
    };
    for (int i = 0; i < n; i++) {
        params.bonds[i].site1 = i;
        params.bonds[i].site2 = (int)(i / l) * l + (i + 1) % l;
    }
    for (int i = n; i < n * 2; i++) {
        params.bonds[i].site1 = i % n;
        params.bonds[i].site2 = (params.bonds[i].site1 + l) % n;
    }
    for (int i = n * 2; i < n * 3; i++) {
        params.bonds[i].site1 = i % n;
        params.bonds[i].site2 = params.bonds[(i - l) % n].site2;
    }

    for (int i = 0; i < n; i++) {
        params.spins[i] = (double_r250() < 0.5 ? 1 : -1);
    }
    for (int i = 0; i < m; i++) {
        params.opers[i].type = IDENT;
    }

    /* Warm up */
    /* for (int i = 0; i < 100; ++i) */
    /*     mc_sweep(); */

    for (int j = 0; j <= tau; ++j) {
        params.hx          = init_hx * (1 - ((double)j / tau)) + final_hx * ((double)j / tau);
        double temperature = init_t * (1 - ((double)j / tau)) + final_t * ((double)j / tau);
        params.beta        = 1.0 / temperature;

        const double hx = params.hx, beta = params.beta;

        /* TODO: nb * J -> Jsum = abs(bonds[i].val) */
        params.pa1 = (2. * J * (double)(nb) + hx * (double)(n)) * beta;
        params.pa2 = 2 * J * (double)(nb) / (2. * J * (double)(nb) + hx * (double)(n));

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
        // cout << b << "\t" << i << "\t" << j << endl;
        eng += (spins[i] * spins[j]);
    }

    printf("%f\t%f\t%f\n", p->beta, p->hx, eng);
    /* if (stp == tau) final_eng = eng; */
}
