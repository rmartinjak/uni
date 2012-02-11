#include "matrix.h"
#include <stdlib.h>
#include <stddef.h>

void **matrixalloc(size_t sz, int n, int m) {
	int i, k;
	void **mat;

	mat = malloc(n * sizeof(void *));
	if (!mat) {
		return NULL;
	}
	for (i=0; i<n; ++i) {
		mat[i] = malloc(m * sz);

		if (!mat[i]) {
			for (k=0; k<i; ++k)
				free(mat[k]);
			free(mat);
			return NULL;
		}
	}
	return mat;
}

void matrixfree(void **mat, int n) {
	int i;
	for (i=0; i<n; ++i) {
		free(mat[i]);
	}
	free(mat);
}
