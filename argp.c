#include "argp.h"

#include <stdio.h>
#include <stdlib.h>

void args_parse (int argc, char **argv, args_t *d) {
    int opt = 0, opt_index = 0;
    int valid = 1;

    /* Use default values */
    args_t *args = d;

    while ((opt = getopt_long(argc, argv, short_options, long_options, &opt_index)) != -1) {
        switch (opt) {
            case 'h':
                printf("%s\n", help_message);
                exit(1);
                break;
            case 'v': printf("Version %s\n", version_message); break;
            case 'q': args->qubo = 1; break;
            case 'f': args->input_file = optarg; break;
            case 't': args->init_t = atof(optarg); break;
            case 'T': args->final_t = atof(optarg); break;
            case 'x': args->init_hx = atof(optarg); break;
            case 'X': args->final_hx = atof(optarg); break;
            case 's': args->tau = atoi(optarg); break;
            case '?': valid = 0; break;
            default: printf("opt = %d\n", opt); break;
        }
    }

    if (!valid) {
        printf("Error: Invalid option\n");
        exit(1);
    }
    return;
}

void args_print (args_t args) {
    printf("QUBO: %d\n", args.qubo);
    printf("Input file: %s\n", args.input_file);
    printf("Initial temperature: %f\n", args.init_t);
    printf("Final temperature: %f\n", args.final_t);
    printf("Initial field: %f\n", args.init_hx);
    printf("Final field: %f\n", args.final_hx);
    printf("Tau: %d\n", args.tau);
}
