#include "mc.h"
#include "util.h"
#include <stdio.h>
#include <stdlib.h>

void mc_sweep (params_t *p, constants_t c) {
    const int n = c.n, nb = c.nb, m = c.m, tau = c.tau;
    int *vrtx = (int *)malloc(m * 4 * sizeof(int));
    int *stck = (int *)malloc(m * 4 * sizeof(int));
    int *frst = (int *)malloc(n * sizeof(int));
    int *last = (int *)malloc(n * sizeof(int));

    diagonal_update(p, c);
    vertices_link(p, c, frst, last, vrtx);
    cluster_update(p, c, frst, vrtx, stck);
    clean_up(vrtx, stck, frst, last, NULL);
}

/* Return the chosen index of the chosen bond */
static int random_bond (bond_t *bonds, size_t nb) {
    double total_weight = 0.0;
    for (size_t i = 0; i < nb; ++i)
        total_weight += (bonds[i].val < 0.0 ? -bonds[i].val : bonds[i].val);
    double random_value      = ((double)rand() / RAND_MAX) * total_weight;
    double cumulative_weight = 0.0;
    for (int i = 0; i < nb; i++) {
        cumulative_weight += (bonds[i].val < 0.0 ? -bonds[i].val : bonds[i].val);
        if (random_value <= cumulative_weight) { return i; }
    }
    return nb - 1; /* Prevent rounding errors */
}

static int random_h (lcoeff_t *lcoeffs, size_t n) {
    double total_weight = 0.0;
    for (size_t i = 0; i < n; ++i)
        total_weight += (lcoeffs[i].val < 0.0 ? -lcoeffs[i].val : lcoeffs[i].val);
    double random_value      = ((double)rand() / RAND_MAX) * total_weight;
    double cumulative_weight = 0.0;
    for (int i = 0; i < n; i++) {
        cumulative_weight += (lcoeffs[i].val < 0.0 ? -lcoeffs[i].val : lcoeffs[i].val);
        if (random_value <= cumulative_weight) { return i; }
    }
    return n - 1; /* Prevent rounding errors */
}

/* If the site is present in the lcoeffs, return the value */
static int check_lcoeffs (lcoeff_t *lcoeffs, size_t n, const int site) {
    for (int i = 0; i < n; ++i) {
        if (lcoeffs[i].site == site) return lcoeffs[i].val;
    }
    return -1;
}

/* Update the diagnonal component in the Matrix */
void diagonal_update (params_t *p, constants_t c) {
    oper_t *opers    = p->opers;
    bond_t *bonds    = p->bonds;
    short *spins     = p->spins;
    int *ni          = &p->ni;
    const double pa1 = p->pa1, pa2 = p->pa2;
    const double hx = p->hx;

    const int lc_len  = c.lc_len;
    lcoeff_t *lcoeffs = p->lcoeffs;

    const int n = c.n, nb = c.nb, m = c.m, tau = c.tau;

    for (int p = 0; p < m; p++) {
        switch (opers[p].type) {
            case IDENT:
                if (double_r250() * (double)(m - *ni) > pa1) break; /* Do not insert operator */
                /* if (check_lcoeffs(lcoeffs, lc_len, p) > hx) break; */

                if (double_r250() > pa2) { /* Insert H operator */

                    /* Only insert H operator if the site's coefficient is less than hx */
                    /* if (check_lcoeffs(lcoeffs, lc_len, p) != -1) break; */

                    opers[p].type = H;
                    opers[p].site = (int)(double_r250() * n);
                    *ni += 1;
                    break;
                }
                /* Insert JJ operator */
                int a  = random_bond(bonds, nb);
                int o1 = bonds[a].site1, o2 = bonds[a].site2;

                short bonding_sign = 0;
                for (int i = 0; i < nb; ++i) {
                    if (bonds[i].site1 == o1 && bonds[i].site2 == o2) {
                        bonding_sign = bonds[i].val > 0 ? 1 : -1;
                        break;
                    }
                    if (bonds[i].site1 == o2 && bonds[i].site2 == o1) {
                        bonding_sign = bonds[i].val > 0 ? 1 : -1;
                        break;
                    }
                }
                if (spins[o1] == -(spins[o2] * bonding_sign)) { /* anti or not */
                    opers[p].type  = JJ;
                    opers[p].site1 = o1;
                    opers[p].site2 = o2;
                    *ni += 1;
                }
                break;
            case H:
            case JJ:
                if (double_r250() * pa1 < (double)(m - *ni + 1)) { /* remove the operator */
                    opers[p].type = IDENT;
                    *ni -= 1;
                    break;
                }
                int k = 1;
                while (k == 1) {
                    if (double_r250() < pa2) {
                        int b  = random_bond(bonds, nb);
                        int o1 = bonds[b].site1, o2 = bonds[b].site2;

                        /* TODO: Check anti along with the bond between two sites */
                        short bonding_sign = 0;
                        for (int i = 0; i < nb; ++i) {
                            if (bonds[i].site1 == o1 && bonds[i].site2 == o2) {
                                bonding_sign = bonds[i].val > 0 ? 1 : -1;
                                break;
                            }
                            if (bonds[i].site1 == o2 && bonds[i].site2 == o1) {
                                bonding_sign = bonds[i].val > 0 ? 1 : -1;
                                break;
                            }
                        }
                        if (spins[o1] == -(spins[o2] * bonding_sign)) {
                            opers[p].type  = JJ;
                            opers[p].site1 = o1;
                            opers[p].site2 = o2;
                            k              = -1;
                        }
                    } else {
                        opers[p].type = H;
                        opers[p].site = (int)(double_r250() * n);
                        k             = -1;
                    }
                }
                break;
            case HX: spins[opers[p].site] *= -1; break;
            default: break;
        }
    }
}

void vertices_link (params_t *p, constants_t c, int *frst, int *last, int *vrtx) {
    oper_t *opers = p->opers;
    bond_t *bonds = p->bonds;
    short *spins  = p->spins;

    const int n = c.n, nb = c.nb, m = c.m, tau = c.tau;

    int v0, s1, s2, v1, v2;

    for (int i = 0; i < n; ++i)
        last[i] = frst[i] = -1;

    for (int p = 0; p < m; ++p) { /* p => operator serial number */
        v0 = 4 * p;

        switch (opers[p].type) {
            case IDENT:
                for (int i = 0; i < 4; i++)
                    vrtx[v0 + i] = -1;
                break;
            case H:
            case HX:
                s1 = opers[p].site;
                v1 = last[s1];

                if (v1 != -1) {
                    vrtx[v1] = v0;
                    vrtx[v0] = v1;
                } else frst[s1] = v0;
                last[s1]     = v0 + 2;
                vrtx[v0 + 1] = -1;
                vrtx[v0 + 3] = -1;
                break;
            case JJ:
                s1 = opers[p].site1;
                s2 = opers[p].site2;
                v1 = last[s1];
                v2 = last[s2];

                if (v1 != -1) {
                    vrtx[v1] = v0;
                    vrtx[v0] = v1;
                } else frst[s1] = v0;

                if (v2 != -1) {
                    vrtx[v2]     = v0 + 1;
                    vrtx[v0 + 1] = v2;
                } else frst[s2] = v0 + 1;
                last[s1] = v0 + 2;
                last[s2] = v0 + 3;
                break;
            default: break;
        }
    }
    for (int i = 0; i < n; i++) {
        if (frst[i] != -1) {
            vrtx[frst[i]] = last[i];
            vrtx[last[i]] = frst[i];
        }
    }
}

void cluster_update (params_t *p, constants_t c, int *frst, int *vrtx, int *stck) {
    oper_t *opers = p->opers;
    bond_t *bonds = p->bonds;
    short *spins  = p->spins;

    const int n = c.n, nb = c.nb, m = c.m, tau = c.tau;

    int f, s, j, k;
    for (int p = 0; p < m; p++) {
        int i = 4 * p;
        if (vrtx[i] < 0) continue;

        if (double_r250() < 0.5) f = -1;
        else f = -2;
        if (opers[p].type == JJ) {
            s = -1;
            for (int k = i; k < i + 4; k++) {
                if (vrtx[k]) {
                    s += 1;
                    stck[s] = vrtx[k];
                    vrtx[k] = f;
                }
            }
        } else {
            s       = 0;
            stck[0] = vrtx[i];
            vrtx[i] = f;

            if (f == -2) opers[p].type = (opers[p].type == H) ? HX : H;
        }
        while (s > -1) {
            k       = stck[s];
            vrtx[k] = f;
            j       = k / 4;
            s -= 1;
            if (opers[j].type == JJ) {
                for (int k = 4 * j; k < 4 * j + 4; k++) {
                    if (vrtx[k] > -1) {
                        s++;
                        stck[s] = vrtx[k];
                        vrtx[k] = f;
                    }
                }
            } else if (f == -2) {
                opers[j].type = (opers[j].type == H) ? HX : H;
            }
        }
    }
    for (int i = 0; i < n; i++) {
        k = frst[i];
        if (k == -1) {
            if (double_r250() < 0.5) spins[i] *= -1;
        } else if (vrtx[k] == -2) spins[i] *= -1;
    }
}
