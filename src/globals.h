#pragma once

typedef enum { IDENT, HX, H, JJ, LC } oper_type_t;
typedef struct {
    oper_type_t type;
    int site, site1, site2; /* JJ uses site1 and site2, else uses site */
} oper_t;

typedef struct {
    int site1, site2;
    double val;
} bond_t;

/* Linear coefficients */
typedef struct {
    int site;
    double val;
} lcoeff_t;

typedef struct {
    int ni; /* NONE IDENT operator count */
    double hx, j, beta;
    double pa1, pa2;
    oper_t *opers;
    bond_t *bonds;
    lcoeff_t *lcoeffs;
    int constant;
    short *spins;
} params_t;
params_t params_init();

typedef struct {
    int n, nb, m, tau;
    double final_hx, init_hx, final_t, init_t;
    int lc_len;
} constants_t;
