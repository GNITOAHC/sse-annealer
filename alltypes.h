#pragma once

typedef enum { IDENT, HX, H, JJ } oper_type_t;
typedef struct {
    oper_type_t type;
    int site, site1, site2; /* JJ uses site1 and site2, else uses site */
} oper_t;

typedef struct {
    int site1, site2;
    double val;
} bond_t;

typedef struct {
    int ni;
    double hx, j, beta;
    double pa1, pa2;
    oper_t *opers;
    bond_t *bonds;
    short *spins;
} params_t;

typedef struct {
    int l, n, nb, m, tau;
    double final_hx, init_hx, final_t, init_t;
} constants_t;
