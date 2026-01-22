#/*******************************************************************************
#  binary_playground.c
#
#  Small demonstration of scalar, vector and struct summation functions.
#
#  Provides three functions (declared in `binary_playground.h`):
#    - `sum_scalars(float a, float b, float c)` returns a+b+c
#    - `sum_vectors(float *v1, const float *v2, const float *v3, size_t len)`
#         computes element-wise v1[i] = v1[i] + v2[i] + v3[i]
#    - `sum_structs(const MyStruct *s1, const MyStruct *s2)`
#         returns a `MyStruct` with fields summed component-wise
#
#  The `main` implements a small CLI to exercise these functions. Example
#  usages:
#
#    Build binary:
#      make
#
#    Run scalar mode (specify three floats):
#      ./binary_playground --mode scalar --a 1.2 --b 3.4 --c 5.6
#
#    Run vector example (uses built-in example vectors):
#      ./binary_playground --mode vector
#
#    Run struct example (uses two example structs):
#      ./binary_playground --mode struct
#
#    Build shared library:
#      make lib
#       
*******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include "binary_playground.h"

/*
 * sum_scalars
 * ---------------
 * Compute the sum of three scalar doubles and return the result.
 */
double sum_scalars(double a, double b, double c) {
    return a + b + c;
}

/*
 * sum_vectors
 * ---------------
 * In-place element-wise sum of three vectors. The result is written into
 * the first vector (`v1`). All arrays are assumed to have at least `len`
 * elements.
 */
void sum_vectors(double *v1, const double *v2, const double *v3, size_t len) {
    for (size_t i = 0; i < len; ++i) {
        v1[i] = v1[i] + v2[i] + v3[i];
    }
}

/*
 * sum_structs
 * ---------------
 * Add corresponding numeric fields of two `MyStruct` instances and return
 * a new `MyStruct` with the aggregated values.
 */
MyStruct sum_structs(const MyStruct *s1, const MyStruct *s2) {
    MyStruct out;
    out.x = s1->x + s2->x;
    out.y = s1->y + s2->y;
    out.z = s1->z + s2->z;
    return out;
}

/*
 * print_usage
 * ---------------
 * Print a short summary of command-line options.
 */
static void print_usage(const char *prog) {
    printf("Usage: %s [--mode scalar|vector|struct] [--use-lib] [--a <float> --b <float> --c <float>]\n", prog);
}

int main(int argc, char **argv) {
    /* Default configuration */
    const char *mode = "scalar"; /* mode can be: scalar, vector, struct */
    /* Default scalar values now set to 1,2,3 as requested */
    double a = 1.0, b = 2.0, c = 3.0; /* scalar inputs for scalar mode */

    /* Command-line options parsed with getopt_long */
    static struct option long_options[] = {
        {"mode", required_argument, 0, 'm'},
        {"a", required_argument, 0, 'a'},
        {"b", required_argument, 0, 'b'},
        {"c", required_argument, 0, 'c'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    int opt;
    while ((opt = getopt_long(argc, argv, "m:a:b:c:h", long_options, NULL)) != -1) {
        switch (opt) {
            case 'm': mode = optarg; break;
            case 'a': a = strtof(optarg, NULL); break;
            case 'b': b = strtof(optarg, NULL); break;
            case 'c': c = strtof(optarg, NULL); break;
            case 'h': print_usage(argv[0]); return 0;
            default: print_usage(argv[0]); return 1;
        }
    }

    /* NOTE: dynamic loading support removed — program always calls the
     * local functions directly. The `Makefile` still offers a `lib` target
     * if you want to build a shared library for other uses, but this
     * executable will not attempt to dlopen it.
     */

    /* Dispatch by mode */
    if (strcmp(mode, "scalar") == 0) {
        /* Direct local call */
        float res = sum_scalars(a, b, c);
        printf("sum_scalars(%f, %f, %f) = %f\n", a, b, c, res);

    } else if (strcmp(mode, "vector") == 0) {
        /* Example fixed-length vectors for demonstration */
        size_t len = 3;
        double v1[3] = {1.0, 2.0, 3.0};
        double v2[3] = {0.5, 0.5, 0.5};
        double v3[3] = {0.1, 0.2, 0.3};

        /* Direct local call */
        sum_vectors(v1, v2, v3, len);

        printf("sum_vectors result: [");
        for (size_t i = 0; i < len; ++i) {
            printf("%f%s", v1[i], (i+1 < len) ? ", " : "");
        }
        printf("]\n");

    } else if (strcmp(mode, "struct") == 0) {
        /* Example structs for demonstration */
        MyStruct s1 = { .x = 1, .y = 2.5f, .z = 3.25 };
        MyStruct s2 = { .x = 4, .y = 1.5f, .z = 0.75 };
        MyStruct out;

        /* Direct local call */
        out = sum_structs(&s1, &s2);

        printf("sum_structs: x=%d, y=%f, z=%f\n", out.x, out.y, out.z);

    } else {
        fprintf(stderr, "Unknown mode '%s'\n", mode);
        print_usage(argv[0]);
        return 1;
    }
    return 0;
}
