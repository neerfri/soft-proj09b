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
void print_modularity_matrix(mod_matrix *mod_mat);
square_matrix *calculate_modularity_matrix(sparse_matrix_arr *adj_matrix, int_vector *vgroup);
elem_vector *calculate_degree_of_vertices(sparse_matrix_arr *adj_matrix, int_vector *vertices_group);
elem_vector *calculate_F_g_array(mod_matrix *mod_mat);

void vec_normalize(elem_vector *vec);
elem_vector *sparse_mat_vec_multiply(sparse_matrix_arr *A, elem_vector *v);

eigen_pair *calculate_leading_eigen_pair(mod_matrix *Bijtag, double precision);

two_division *divide_network_in_two(mod_matrix *mod_mat, eigen_pair *leading_eigen_pair, int use_improve);
/*int improve_network_division(square_matrix *mod_mat, two_division *division);*/
int improve_network_division(mod_matrix *mod_mat, two_division *division);

mod_matrix *allocate_partial_modularity_matrix(sparse_matrix_arr *adj_matrix, int_vector *vertices_group);
void print_clusters(n_division *division);

/* Runs Algorithm 2 and returns a new devision, given the Adjacency matrix
	If use_imporve works, runs algorithm 4*/
two_division *algorithm2(sparse_matrix_arr *adj_matrix, int_vector *vgroup, double precision, int use_improve);

/* Implementing Algorithm 3 */
n_division* algorithm3 (sparse_matrix_arr *adj_matrix, double precision, int use_improve);

#endif /* __FINAL_PROJECT_SHARED_H */
