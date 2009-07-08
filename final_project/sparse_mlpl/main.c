#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "sparse_matrix.h"

void print_elem_vector(elem *vector, int n) {
	int i;
	for(i=0; i < n;i++) {
		printf("%f ", vector[i]);
	}
}

void print_int_vector(int *vector, int n) {
	int i;
	for(i=0; i < n;i++) {
		printf("%d ", vector[i]);
	}
}

void print_sparse_matrix_data(sparse_matrix_arr *matrix) {
	int n = matrix->rowptr[matrix->n];
	print_elem_vector(matrix->values, n);
	printf("\n");
	print_int_vector(matrix->colind, n);
	printf("\n");
	print_int_vector(matrix->rowptr, n+1);
}

int main(void) {
	elem *vector, *result;
	sparse_matrix_arr* matrix;
	if ((matrix = allocate_and_read_matrix(stdin)) == NULL) {
		exit(EXIT_FAILURE);
	}
#ifdef DEBUG
	print_sparse_matrix_data(matrix);
#endif
	if ((vector = allocate_and_read_vector(matrix->n)) == NULL) {
		free(matrix);
		exit(EXIT_FAILURE);
	}
#ifdef DEBUG
	print_elem_vector(vector, matrix->n);
#endif
	if ((result = allocate_vector(matrix->n)) == NULL) {
		free_sparse_matrix_arr(matrix);
		free(vector);
		exit(EXIT_FAILURE);
	}
	mult_sparse_arr(matrix, vector, result);
	print_elem_vector(result, matrix->n);
	return EXIT_SUCCESS;
}