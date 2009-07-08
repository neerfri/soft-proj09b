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
	matrix->values = (elem *)malloc(numNnz*sizeof(elem));
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
	return matrix;
}

void free_sparse_matrix_arr(sparse_matrix_arr* matrix) {
	free(matrix->colind);
	free(matrix->rowptr);
	free(matrix->values);
	free(matrix);
}

void  mult_sparse_arr(const sparse_matrix_arr *A, const elem* v, elem* result) {
	int j, val_index;
	val_index = 0;
	for(j=0; j < A->n; j++) {
		result[j] = 0;
		for(; val_index<A->rowptr[j+1]; val_index++) {
			result[j] = result[j] + (A->values[val_index]*v[A->colind[val_index]]);
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

elem *allocate_and_read_vector(int n) {
	int i;
	double scanf_receptor;
	elem *vector;
	if ((vector = allocate_vector(n)) == NULL) {
		return NULL;
	}
	for(i=0; i<n; i++) {
		if (scanf("%lf", &scanf_receptor) < 1) {
			fprintf(stderr, "Error reading vector value in pos: %d", i);
			free(vector);
			return NULL;
		}
		vector[i] = scanf_receptor;
	}
	return vector;
}

sparse_matrix_arr *allocate_and_read_matrix(FILE *fp) {
	int n=0, numNnz=0;
	int i, j, val_index, pos;
	double scanf_receptor;
	elem *values;
	sparse_matrix_arr* matrix;
	if (fscanf(fp, "%d", &n) < 1) {
		fprintf(stderr, "Error reading matrix size");
		return NULL;
	}
#ifdef DEBUG
	printf("matrix size: %d\n", n);
#endif
	values = malloc(n*n*sizeof(elem));
	if (values == NULL) {
		MEMORY_ALLOCATION_FAILURE_AT("read_matrix: values");
		return NULL;
	}
	for(j=0; j<n; j++) {
		for(i=0; i<n; i++) {
			pos = i + (j*n);
			if (fscanf(fp, "%lf", &scanf_receptor) < 1) {
				fprintf(stderr, "Error reading matrix values at pos: %d\n", pos);
				free(values);
				return NULL;
			}
			values[pos] = scanf_receptor;
			if (values[pos] != 0)
				numNnz++;
		}
	}
		
	if ((matrix = allocate_sparse_matrix_arr(n, numNnz)) == NULL)
		return NULL;
	val_index = 0;
	matrix->rowptr[0] = 0;
	for(j=0; j<n; j++) {
		for(i=0; i<n; i++) {
			pos = i + (j*n);
			if (values[pos] != 0) {
				matrix->values[val_index] = values[pos];
				matrix->colind[val_index] = i;
#ifdef DEBUG
				printf("Found new non-0 val: %lf, colind: %d, val_index:%d\n", values[pos], matrix->colind[val_index], val_index);
#endif
				val_index++;
			}
		}
		matrix->rowptr[j+1] = val_index;
	}
	free(values);
	return matrix;
}
