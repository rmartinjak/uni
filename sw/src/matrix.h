#ifndef MATRIX_H
#define MATRIX_H
#include <stddef.h>

void **matrixalloc(size_t sz, int n, int m);
void matrixfree(void **matrix, int n);

#endif
