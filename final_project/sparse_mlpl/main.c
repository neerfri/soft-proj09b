#include "default_includes.h"

int main(void) {
	elem *vector, *result;
	sparse_matrix_arr* matrix;
	if ((matrix = allocate_and_read_matrix(stdin)) == NULL) {
		exit(EXIT_FAILURE);
	}
	if ((vector = allocate_and_read_vector(stdin, matrix->n)) == NULL) {
		free_sparse_matrix_arr(matrix);
		exit(EXIT_FAILURE);
	}
	if ((result = allocate_vector(matrix->n)) == NULL) {
		free(vector);
		free_sparse_matrix_arr(matrix);
		exit(EXIT_FAILURE);
	}
	mult_sparse_arr(matrix, vector, result);
	print_elem_vector(result, matrix->n);
	free(result);
	free(vector);
	free_sparse_matrix_arr(matrix);
	return EXIT_SUCCESS;
}
