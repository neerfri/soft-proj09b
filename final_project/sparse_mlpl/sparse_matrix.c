/*
sparse_matrix_arr* allocate_sparse_matrix_arr(int n, int numNnz);
void free_sparse_matrix_arr(sparse_matrix_arr* matrix);
void  mult_sparse_arr(const sparse_matrix_arr *A, const elem* v, elem* result);
*/

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
	matrix->rowptr = (int *)malloc((numNnz+1)*sizeof(int));
	if (matrix->rowptr == NULL) {
		MEMORY_ALLOCATION_FAILURE_AT("allocate_sparse_matrix_arr: matrix->rowptr");
		free(matrix->colind);
		free(matrix->values);
		free(matrix);
		return NULL;
	}
	return matrix;
}

int read_values(sparse_matrix_arr *matrix) {
	
	return 0;
}

sparse_matrix_arr *read_matrix() {
	int n=0, numNnz=0;
	int i, j, val_index, pos;
	double scanf_receptor;
	elem *values;
	sparse_matrix_arr* matrix;
	if (scanf("%d", &n) < 1) {
		fprintf(stderr, "Error reading matrix size");
		exit(EXIT_FAILURE);
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
			if (scanf("%lf", &scanf_receptor) < 1) {
				fprintf(stderr, "Error reading matrix values at pos: %d\n", pos);
				/*TODO:free memory */
				return NULL;
			}
			values[pos] = scanf_receptor;
			if (values[pos] != 0)
				numNnz++;
		}
	}
	/*
	for(i=0; i<n; i++) {
		for(j=0; j<n; j++) {
			printf("%f ", (float)(values[i + (j*n)]));
		}
		printf("\n");
	}*/
		
	if ((matrix = allocate_sparse_matrix_arr(n, numNnz)) == NULL)
		return NULL;
	val_index = 0;
	matrix->rowptr[0] = 0;
	for(j=0; j<n; j++) {
		matrix->rowptr[j] = val_index;
		for(i=0; i<n; i++) {
			pos = i + (j*n);
			if (values[i + (j*n)] != 0) {
				matrix->values[val_index] = *(values+(i + (j*n)));
				matrix->colind[val_index] = i;
#ifdef DEBUG
				printf("Found new non-0 val: %d, colind: %d, val_index:%d\n", (int)(*(values+(i + (j*n)))), matrix->colind[val_index], val_index);
#endif
				val_index++;
			}
		}
	}
	matrix->rowptr[j] = val_index;
	return matrix;
}

void print_sparse_matrix_data(sparse_matrix_arr *matrix) {
	int i;
	for(i=0; matrix->rowptr[i] < matrix->n;i++) {
		printf("%f ", matrix->values[i]);
	}
	printf("\n");
	for(i=0; matrix->rowptr[i] < matrix->n;i++) {
		printf("%d ", matrix->colind[i]);
	}
	printf("\n");
	for(i=0; matrix->rowptr[i] < matrix->n;i++) {
		printf("%d ", matrix->rowptr[i]);
	}
	printf("%d", matrix->rowptr[i++]);
	printf("\n");
}

int main(void) {
	sparse_matrix_arr* matrix = read_matrix();
	if (matrix == NULL) {
		exit(EXIT_FAILURE);
	}
	print_sparse_matrix_data(matrix);
	return EXIT_SUCCESS;
}
