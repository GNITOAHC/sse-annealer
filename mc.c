#include "mc.h"
#include "util.h"
#include <stdio.h>
#include <stdlib.h>

void mc_sweep (params_t *p, constants_t c) {
    const int l = c.l, n = c.n, nb = c.nb, m = c.m, tau = c.tau;
    int *vrtx = (int *)malloc(m * 4 * sizeof(int));
    int *stck = (int *)malloc(m * 4 * sizeof(int));
    int *frst = (int *)malloc(n * sizeof(int));
    int *last = (int *)malloc(n * sizeof(int));

    diagonal_update(p, c);                  // 對角更新
    vertices_link(p, c, frst, last, vrtx);  // 頂點連結
    cluster_update(p, c, frst, vrtx, stck); // 叢集更新
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

// Update the diagnonal component in the Matrix
void diagonal_update (params_t *p, constants_t c) {
    oper_t *opers    = p->opers;
    bond_t *bonds    = p->bonds;
    short *spins     = p->spins;
    int *ni          = &p->ni;
    const double pa1 = p->pa1, pa2 = p->pa2;

    const int l = c.l, n = c.n, nb = c.nb, m = c.m, tau = c.tau;

    for (int p = 0; p < m; p++) {
        if (opers[p].type == IDENT) { // I
            if (double_r250() * (double)(m - *ni) < pa1) {
                if (double_r250() < pa2) {
                    int a  = random_bond(bonds, nb);
                    int o1 = bonds[a].site1, o2 = bonds[a].site2;

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
                    if (spins[o1] == -(spins[o2] * bonding_sign)) { // anti or not
                        opers[p].type  = JJ;
                        opers[p].site1 = o1;
                        opers[p].site2 = o2;
                        *ni += 1;
                    }
                } else {
                    opers[p].type = H;
                    opers[p].site = (int)(double_r250() * n);
                    *ni += 1;
                }
            }
        } else if (opers[p].type == H || opers[p].type == JJ) { // h or JJ
            if (double_r250() * pa1 < (double)(m - *ni + 1)) {  // remove the operator
                opers[p].type = IDENT;
                *ni -= 1;
            } // One then goes to the next operator in the list
            else {
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
            }
        } else spins[opers[p].site] *= -1;
    }
}

void vertices_link (params_t *p, constants_t c, int *frst, int *last, int *vrtx) {
    oper_t *opers = p->opers;
    bond_t *bonds = p->bonds;
    short *spins  = p->spins;

    const int l = c.l, n = c.n, nb = c.nb, m = c.m, tau = c.tau;

    int v0, s1, s2, v1, v2;
    // frst => 進入, last => 出來
    // example:
    //   the s-th operator :  oper[0][s]=-2,oper[1][s]=5
    //   means the s-th operator is hx and in site 5

    for (int i = 0; i < n; ++i)
        last[i] = frst[i] = -1;

    for (int p = 0; p < m; ++p) { // p => operator 編號
        v0 = 4 * p;

        if (opers[p].type == IDENT) { // oper[1][p] == -1 => Identity
            for (int i = 0; i < 4; i++)
                vrtx[v0 + i] = -1;
        } else {
            if (opers[p].type == H || opers[p].type == HX) {
                s1 = opers[p].site;
                v1 = last[s1];

                if (v1 != -1) {
                    vrtx[v1] = v0;
                    vrtx[v0] = v1;
                } else frst[s1] = v0;
                last[s1]     = v0 + 2;
                vrtx[v0 + 1] = -1;
                vrtx[v0 + 3] = -1;
            } else {
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
            }
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

    const int l = c.l, n = c.n, nb = c.nb, m = c.m, tau = c.tau;

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
            j       = k / 4; // might
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
