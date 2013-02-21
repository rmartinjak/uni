#include "matrix.h"
#include <stdlib.h>
#include <stddef.h>

void *matrixalloc(size_t sz, int n, int m)
{
    void **mat;
    int i;

    mat = malloc(n * sizeof(void *));
    if (mat != NULL) {
        mat[0] = malloc(n * m * sz);
        if (mat[0] == NULL) {
            free(mat);
            return NULL;
        }

        for (i = 1; i < n; i++)
            mat[i] = (char*)mat[i-1] + (m*sz);
    }
    return mat;
}

void matrixfree(void *mat)
{
    void **m = mat;
    free(*m);
    free(mat);
}
