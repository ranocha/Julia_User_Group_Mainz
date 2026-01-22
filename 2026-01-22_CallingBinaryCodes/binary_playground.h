#ifndef BINARY_PLAYGROUND_H
#define BINARY_PLAYGROUND_H

#include <stddef.h>

typedef struct {
    int x;
    double y;
    double z;
} MyStruct;

/* Sum three scalar doubles and return the result */
double sum_scalars(double a, double b, double c);

/* In-place sum of three vectors: result stored into v1. `len` is the length. */
void sum_vectors(double *v1, const double *v2, const double *v3, size_t len);

/* Sum fields of two structs and return the resulting struct */
MyStruct sum_structs(const MyStruct *s1, const MyStruct *s2);

#endif /* BINARY_PLAYGROUND_H */
