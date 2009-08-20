#ifndef __FINAL_PROJECT_SHARED_H
#define __FINAL_PROJECT_SHARED_H

#include "sparse_matrix.h"
#include "data_structures.h"

#define MEMORY_ALLOCATION_FAILURE_AT(func) fprintf(stderr, "MEMORY ALLOCATION FAILURE AT '%s'. Aborting\n", func);
#define IS_POSITIVE(X) ((X) > 0.00001)
/* Type definitions */

/* Holds a square matrix including it's size and values.
 * The values of the matrix are of type elem
 */
typedef struct {
	int n;			/* size */
	elem *values;
} square_matrix;

/************************************************/

void free_square_matrix(square_matrix *matrix);

void print_elem_vector(elem *vector, int n);
void print_int_vector(int *vector, int n);
void print_sparse_matrix(sparse_matrix_arr *matrix);

int degree_of_vertice(int i, sparse_matrix_arr *matrix);
square_matrix *allocate_square_matrix(int n);
void print_square_matrix(square_matrix *mat);
square_matrix *calculate_modularity_matrix(sparse_matrix_arr *adj_matrix, int_vector *vgroup);


eigen_pair *calculate_leading_eigen_pair(square_matrix *mod_mat, double precision);

two_division *divide_network_in_two(square_matrix *mod_mat, eigen_pair *leading_eigen_pair);

#endif /* __FINAL_PROJECT_SHARED_H */
