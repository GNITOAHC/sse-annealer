#include "globals.h"
#include <stdlib.h>

params_t params_init () {
    return (params_t) {
        .ni    = 0,
        .hx    = 0.0,
        .j     = 0.0,
        .beta  = 0.0,
        .pa1   = 0.0,
        .pa2   = 0.0,
        .opers = NULL,
        .bonds = NULL,
        .spins = NULL,
    };
}
