#pragma once

#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    /* name of long option */
    const char *name;
    /*
     * one of no_argument, required_argument, and optional_argument:
     * whether option takes an argument
     */
    int has_arg;
    /* if not NULL, set *flag to val when option found */
    int *flag;
    /* if flag not NULL, value to set *flag to; else return value */
    int val;
    /* argument */
    const char *arg;
    /* description for usage */
    const char *desc;
} option_t;

static inline char *help_mes (option_t *options) {
    int count          = 0;
    size_t buffer_size = 256; /* Initial buffer size */
    char *message      = (char *)malloc(buffer_size);

    if (message == NULL) { return NULL; } /* Handle memory allocation failure */
    message[0] = '\0';                    /* Start with an empty string */

    char option_indication[20] = "Options:\n"; /* Buffer for the option indication */
    strcat(message, option_indication);

    /* Build the help message */
    while (options[count].name != NULL) {
        char option_line[128]; /* Buffer for each option line */
        char arg_string[128];  /* Buffer for argument name */
        char short_option[8];  /* Buffer for short option */

        if (options[count].arg != NULL) {
            snprintf(arg_string, sizeof(arg_string), " <%s>", options[count].arg);
        } else {
            arg_string[0] = '\0';
        }

        if (options[count].val == 0) snprintf(short_option, sizeof(short_option), "  ");
        else snprintf(short_option, sizeof(short_option), "  -%c, ", options[count].val);

        snprintf(option_line, sizeof(option_line), "%s--%s%s\n      %s\n\n", short_option,
                 options[count].name, options[count].has_arg != no_argument ? arg_string : "",
                 options[count].desc ? options[count].desc : "");

        /* Check if we need to expand the buffer */
        if (strlen(message) + strlen(option_line) + 1 > buffer_size) {
            buffer_size *= 2;
            char *new_message = (char *)realloc(message, buffer_size);
            if (new_message == NULL) {
                free(message);
                return NULL; /* Handle memory allocation failure */
            }
            message = new_message;
        }

        strcat(message, option_line);
        count++;
    }

    return message;
}

static inline struct option *long_opts (option_t *options) {
    int count = 0;
    while (options[count].name != NULL) {
        ++count;
    }

    /* Allocate memory for long options */
    struct option *long_options = (struct option *)malloc(sizeof(struct option) * (count + 1));
    if (long_options == NULL) { return NULL; }
    /* Convert option_t to option */
    for (int i = 0; i < count; ++i) {
        long_options[i].name    = options[i].name;
        long_options[i].has_arg = options[i].has_arg;
        long_options[i].flag    = options[i].flag;
        long_options[i].val     = options[i].val;
    }

    /* The last element of the array must be filled with zeros */
    long_options[count].name    = NULL;
    long_options[count].has_arg = 0;
    long_options[count].flag    = NULL;
    long_options[count].val     = 0;

    return long_options;
}

static inline char *short_opts (option_t *options) {
    int count  = 0;
    int length = 0;

    /* Calculate the length of the short options string */
    while (options[count].name != NULL) {
        if (options[count].val != 0) {
            ++length;
            if (options[count].has_arg == required_argument) {
                ++length;
            } else if (options[count].has_arg == optional_argument) {
                length += 2;
            }
        }
        ++count;
    }

    /* Allocate memory for the short options string */
    char *short_opts = (char *)malloc(length + 1);
    if (short_opts == NULL) {
        return NULL; // Handle memory allocation failure
    }

    /* Build the short options string */
    int index = 0;
    for (int i = 0; i < count; ++i) {
        if (options[i].val != 0) {
            short_opts[index++] = (char)options[i].val;
            if (options[i].has_arg == required_argument) {
                short_opts[index++] = ':';
            } else if (options[i].has_arg == optional_argument) {
                short_opts[index++] = ':';
                short_opts[index++] = ':';
            }
        }
    }
    short_opts[index] = '\0'; // Null-terminate the string

    return short_opts;
}

static inline void free_options (struct option *long_opts, char *short_opts, char *help_message) {
    free(long_opts);
    free(short_opts);
    free(help_message);
}
