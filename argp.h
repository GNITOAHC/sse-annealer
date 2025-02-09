#pragma once

#include <getopt.h>

typedef struct {
    /* Flags, 1 for true, default 0 */
    short qubo;
    /* Input file path */
    char *input_file;
    /* MC Sweep constants */
    double init_t, final_t, init_hx, final_hx;
    int tau;
    /* Print final configuration  */
    int print_conf;
    /* Triangular setup if used */
    int tri_l;
} args_t;

static struct option long_options[] = {
    {       "help",       no_argument, 0, 'h' },
    {    "version",       no_argument, 0, 'v' },
    {       "qubo",       no_argument, 0, 'q' },
    {       "file", required_argument, 0, 'f' },
    {     "init-t", required_argument, 0, 't' },
    {    "final-t", required_argument, 0, 'T' },
    {    "init-hx", required_argument, 0, 'x' },
    {   "final-hx", required_argument, 0, 'X' },
    {        "tau", required_argument, 0, 's' },
    { "print-conf",       no_argument, 0, 'p' },
    {      "tri-l", required_argument, 0, 'l' },
    {            0,                 0, 0,   0 }
};
static char short_options[] = "hvqf:t:T:x:X:s:pl:";

static char version_message[] = "SSE v0.1\n";
static char help_message[]    = "Usage: %s [OPTION]...\n"
                                "  -h, --help       Display this help and exit\n"
                                "  -v, --version    Display version information and exit\n"
                                "  -q, --qubo       Use QUBO format\n"
                                "  -f, --file       Input file path\n"
                                "  -t, --init-t     Initial temperature\n"
                                "  -T, --final-t    Final temperature\n"
                                "  -x, --init-hx    Initial Hamiltonian\n"
                                "  -X, --final-hx   Final Hamiltonian\n"
                                "  -s, --tau        Tau\n"
                                "  -p, --print-conf Print final configuration\n"
                                "  -l. --tri-l      Default triangular setup\n";

args_t args_default();
void args_parse(int argc, char **argv, args_t *d);
void args_print(args_t args);
