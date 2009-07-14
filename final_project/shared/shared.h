#ifndef __FINAL_PROJECT_SHARED_H
#define __FINAL_PROJECT_SHARED_H

#include "sparse_matrix.h"

#define MEMORY_ALLOCATION_FAILURE_AT(func) fprintf(stderr, "MEMORY ALLOCATION FAILURE AT '%s'. Aborting\n", func);

/* Type definitions */

/* Integer array with stored length
 * This structure is used for holding vertices groups.
 **/
typedef struct {
  int     count;    /* size                  */
  int*    vertices;
} int_vector;

/* Holds a square matrix including it's size and values.
 * The values of the matrix are of type elem
 */
typedef struct {
	int n;			/* size */
	elem *values;
} square_matrix;

/************************************************/

void print_elem_vector(elem *vector, int n);
void print_int_vector(int *vector, int n);
void print_sparse_matrix(sparse_matrix_arr *matrix);

sparse_matrix_arr *read_adjacency_matrix(const char* file);
int_vector *read_vertices_group_file(const char* file, int max_count);
int read_n_vertices_group(FILE *fp, int_vector *vertices, int n);

int degree_of_vertice(int i, sparse_matrix_arr *matrix);
square_matrix *allocate_square_matrix(int n);
void print_square_matrix(square_matrix *mat);
square_matrix *calculate_modularity_matrix(sparse_matrix_arr *adj_matrix, int_vector *vgroup);

#endif /* __FINAL_PROJECT_SHARED_H */
