#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "sparse_matrix.h"

#ifndef MEMORY_ALLOCATION_FAILURE_AT
#define MEMORY_ALLOCATION_FAILURE_AT(func) fprintf(stderr, "MEMORY ALLOCATION FAILURE AT '%s'. Aborting\n", func);
#endif

typedef struct
{
  int     count;    /* size                  */
  int*    vertices;
} vertices_group;

typedef struct
{
	int n;			/* size */
	double *values;
} square_double_matrix;

/********************* Temporary Helper methods ************************/
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
/*********************************************************************/

int read_n_vertices_group(FILE *fp, vertices_group *vertices, int n) {
	int scanf_receptor, i, *values;
	if ((values = calloc(n, sizeof(int))) == NULL) {
		MEMORY_ALLOCATION_FAILURE_AT("read_vertices_group: values");
		return -1;
	}
	for(i=0; i < n; i++) {
		if (feof(fp) || fscanf(fp, "%d", &scanf_receptor) < 1) {
			break;
		} else {
			values[i] = scanf_receptor;
		}
	}
	vertices->count = i;
	if ((vertices->vertices = calloc(vertices->count, sizeof(int))) == NULL) {
		MEMORY_ALLOCATION_FAILURE_AT("read_vertices_group: vertices->vertices");
		free(values);
		return -1;
	}
	for(i=0; i<vertices->count; i++) {
		vertices->vertices[i] = values[i];
	}
	free(values);
	return vertices->count;
}

int degree_of_vertice(int i, sparse_matrix_arr *matrix) {
	int j, degree=0;
	for(j=0; j<matrix->rowptr[matrix->n]; j++) {
		if (matrix->colind[j] == i)
			degree++;
	}
	return degree;
}


square_double_matrix *allocate_square_double_matrix(int n) {
	square_double_matrix *matrix;
	if ((matrix = malloc(sizeof(square_double_matrix))) == NULL) {
		MEMORY_ALLOCATION_FAILURE_AT("allocate_square_double_matrix: matrix");
		return NULL;
	}
	if ((matrix->values = calloc(n*n, sizeof(double))) == NULL) {
		MEMORY_ALLOCATION_FAILURE_AT("allocate_square_double_matrix: matrix->values");
		free(matrix);
		return NULL;
	}
	matrix->n = n;
	return matrix;
}

square_double_matrix *calculate_modularity_matrix(sparse_matrix_arr *adj_matrix, vertices_group *vgroup) {
	square_double_matrix *mod_mat;
	int i,j, pos, K_i, K_j, mat_val_index=0;
	double A_i_j;
	double *F_g;
	if ((mod_mat = allocate_square_double_matrix(vgroup->count)) == NULL) {
		return NULL;
	}
	if ((F_g = calloc(vgroup->count, sizeof(double))) == NULL) {
		free(mod_mat);
		return NULL;
	}
	for(i=0; i<vgroup->count; i++) {
		F_g[i] = 0;
		mat_val_index = adj_matrix->rowptr[vgroup->vertices[i]];
		for(j=0; j<vgroup->count; j++) {
			if (j <= adj_matrix->rowptr[vgroup->vertices[i]+1]) {
				/* we are still in values for column j row i; */
				/* first we need to advance the sparse matrix pointer */
				while(mat_val_index < adj_matrix->rowptr[vgroup->vertices[i]+1] &&
						adj_matrix->colind[mat_val_index] < vgroup->vertices[j]) {
					mat_val_index++;
				}
				/* if its in adj_matrix values, take it, otherwise it's 0*/
				/* first condition makes sure we don't overflow from the array boundries */
				if (mat_val_index < adj_matrix->rowptr[adj_matrix->n] &&
						adj_matrix->colind[mat_val_index] == vgroup->vertices[j]) {
					A_i_j = adj_matrix->values[mat_val_index];
					mat_val_index++;
				} else {
					A_i_j = 0;
				}
			} else {
				/* there are no more values for this row, meaning value is 0 */
				A_i_j = 0;
			}
			pos = (i*mod_mat->n) + j;
			/* next few lines claculates: B_i_j = A_i_j - (K_i*K_j)/M */
			K_i = degree_of_vertice(vgroup->vertices[i], adj_matrix);
			K_j = degree_of_vertice(vgroup->vertices[j], adj_matrix);
			mod_mat->values[pos] = A_i_j - (((double)(K_i*K_j))/adj_matrix->rowptr[adj_matrix->n]);
			/* next line calculates: F_g[i] = Sum(Over all j in vgroup)[B[g]_i_j]*/
			F_g[i] = F_g[i] +  mod_mat->values[pos];
		}
	}
	for(i=0; i<vgroup->count; i++) {
		pos = (i*mod_mat->n) + i;
		mod_mat->values[pos] = mod_mat->values[pos] - F_g[i];
	}
	free(F_g);
	return mod_mat;
}

void print_square_double_matrix(square_double_matrix *mat) {
	int i,j;
	printf("%d ", mat->n);
	for(i=0; i< mat->n; i++) {
		printf("\r\n");
		for(j=0; j<mat->n; j++) {
			printf("%f ", mat->values[(i*mat->n)+j]);
		}
	}
}


int main(int argc, char **argv) {
	FILE *fp;
	sparse_matrix_arr* adj_matrix;
	vertices_group *vgroup;
	square_double_matrix *modularity_matrix;
	if (argc < 3) {
		fprintf(stderr, "Invalid Arguments, Aborting.\n");
		fprintf(stderr, "Usage: %s <adjacency-mat-file> <group-file>\n", argv[0]);
		return EXIT_FAILURE;
	}
#ifdef DEBUG
	printf("opening file: %s\n", argv[1]);
#endif
	if ((fp = fopen(argv[1], "r")) == NULL) {
		/*File open error ! abort */
		fprintf(stderr, "Could not open adjacency matrix file: '%s'. Aborting.\n", argv[1]);
		return EXIT_FAILURE;
	}
#ifdef DEBUG
	printf("file opened successfully.\n");
#endif
	if ((adj_matrix = allocate_and_read_matrix(fp)) == NULL) {
		fclose(fp);
		exit(EXIT_FAILURE);
	}
	fclose(fp);
	if ((fp = fopen(argv[2], "r")) == NULL) {
		/*File open error ! abort */
		fprintf(stderr, "Could not open vertices file: '%s'. Aborting.\n", argv[2]);
		free(adj_matrix);
		return EXIT_FAILURE;
	}
	if ((vgroup = malloc(sizeof(vertices_group))) == NULL) {
		MEMORY_ALLOCATION_FAILURE_AT("main: vgroup");
		free(adj_matrix);
		return EXIT_FAILURE;
	}
	read_n_vertices_group(fp, vgroup, adj_matrix->n);
	fclose(fp);
	modularity_matrix = calculate_modularity_matrix(adj_matrix, vgroup);
	print_square_double_matrix(modularity_matrix);
	return EXIT_SUCCESS;
}
