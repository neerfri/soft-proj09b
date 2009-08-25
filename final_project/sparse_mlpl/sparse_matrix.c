#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "sparse_matrix.h"

sparse_matrix_arr* allocate_sparse_matrix_arr(int n, int numNnz) {
	sparse_matrix_arr *matrix;
	matrix = (sparse_matrix_arr *)malloc(sizeof(sparse_matrix_arr));
	if (matrix == NULL) {
		MEMORY_ALLOCATION_FAILURE_AT("allocate_sparse_matrix_arr: matrix");
		return NULL;
	}
	matrix->n = n;
	matrix->colind = (int *)malloc(numNnz*sizeof(int));
	if (matrix->colind == NULL) {
		MEMORY_ALLOCATION_FAILURE_AT("allocate_sparse_matrix_arr: matrix->colind");
		free(matrix);
		return(NULL);
	}
	matrix->values = (elem *)calloc(numNnz,sizeof(elem));
	if (matrix->colind == NULL) {
		MEMORY_ALLOCATION_FAILURE_AT("allocate_sparse_matrix_arr: matrix->values");
		free(matrix->colind);
		free(matrix);
		return(NULL);
	}
	matrix->rowptr = (int *)malloc((n+1)*sizeof(int));
	if (matrix->rowptr == NULL) {
		MEMORY_ALLOCATION_FAILURE_AT("allocate_sparse_matrix_arr: matrix->rowptr");
		free(matrix->colind);
		free(matrix->values);
		free(matrix);
		return NULL;
	}
	matrix->rowptr[matrix->n] = numNnz;
	return matrix;
}

void free_sparse_matrix_arr(sparse_matrix_arr* matrix) {
	free(matrix->colind);
	free(matrix->rowptr);
	free(matrix->values);
	free(matrix);
}

void  mult_sparse_arr(const sparse_matrix_arr *A, const elem* v, elem* result) {
	int i, j, val_index;
	/* this scans the matrix lines */
	for(i=0; i<A->n; i++) {
		/* now scan the values in that line */
		for(val_index = A->rowptr[i]; val_index<A->rowptr[i+1]; val_index++) {
			j = A->colind[val_index];
			result[i] += A->values[val_index]*v[j];
		}
	}
}

elem *allocate_vector(int n) {
	int i;
	elem *vector;
	if ((vector = malloc(n*sizeof(elem))) == NULL) {
		MEMORY_ALLOCATION_FAILURE_AT("allocate_vector: vector");
		return NULL;
	}
	for(i=0; i< n; i++) {
		vector[i] = 0;
	}
	return vector;
}
