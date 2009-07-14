#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "shared.h"

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

void print_sparse_matrix(sparse_matrix_arr *matrix) {
	int i, j, mat_val_index;
	elem A_i_j;
	for(i=0; i<matrix->n; i++) {
		mat_val_index = matrix->rowptr[i];
		for(j=0; j<matrix->n; j++) {
			if (mat_val_index < matrix->rowptr[i+1]) {
				/* we are still in values for column j row i; */
				/* first we need to advance the sparse matrix pointer */
				while(mat_val_index < matrix->rowptr[i+1] &&
					matrix->colind[mat_val_index] < j) {
					mat_val_index++;
				}
				/* if its in matrix values, take it, otherwise it's 0*/
				/* first condition makes sure we don't overflow from the array boundaries */
				if (mat_val_index < matrix->rowptr[matrix->n] &&
						matrix->colind[mat_val_index] == j) {
					A_i_j = matrix->values[mat_val_index];
					mat_val_index++;
				} else {
					A_i_j = 0;
				}
			} else {
				/* there are no more values for this row, meaning value is 0 */
				A_i_j = 0;
			}
			printf("%d ", (int)A_i_j);
		}
		printf("\n");
	}
}

void free_int_vector(int_vector *vector) {
	free(vector->vertices);
	free(vector);
}

void free_square_matrix(square_matrix *matrix) {
	free(matrix->values);
	free(matrix);
}

sparse_matrix_arr *read_adjacency_matrix(const char* file) {
	FILE *fp;
	sparse_matrix_arr *adj_matrix;
	if ((fp = fopen(file, "r")) == NULL) {
		/*File open error ! abort */
		fprintf(stderr, "Could not open adjacency matrix file: '%s'. Aborting.\n", file);
		return NULL;
	}
	if ((adj_matrix = allocate_and_read_matrix(fp)) == NULL) {
		MEMORY_ALLOCATION_FAILURE_AT("read_adjacency_matrix: adj_matrix");
		fclose(fp);
		return NULL;
	}
	fclose(fp);
	return adj_matrix;
}

int_vector *read_vertices_group_file(const char* file, int max_count) {
	FILE *fp;
	int_vector *vgroup;
	if ((fp = fopen(file, "r")) == NULL) {
		/*File open error ! abort */
		fprintf(stderr, "Could not open vertices file: '%s'. Aborting.\n", file);
		return NULL;
	}
	if ((vgroup = malloc(sizeof(int_vector))) == NULL) {
		MEMORY_ALLOCATION_FAILURE_AT("read_vertices_group_file: vgroup");
		fclose(fp);
		return NULL;
	}
	if (read_n_vertices_group(fp, vgroup, max_count) < 1) {
		/*Error reading values from file or no values found*/
		free(vgroup);
		fclose(fp);
		return NULL;
	}
	fclose(fp);
	return vgroup;
}
int read_n_vertices_group(FILE *fp, int_vector *vertices, int n) {
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

square_matrix *allocate_square_matrix(int n) {
	square_matrix *matrix;
	if ((matrix = malloc(sizeof(square_matrix))) == NULL) {
		MEMORY_ALLOCATION_FAILURE_AT("allocate_square_double_matrix: matrix");
		return NULL;
	}
	if ((matrix->values = calloc(n*n, sizeof(elem))) == NULL) {
		MEMORY_ALLOCATION_FAILURE_AT("allocate_square_double_matrix: matrix->values");
		free(matrix);
		return NULL;
	}
	matrix->n = n;
	return matrix;
}

void print_square_matrix(square_matrix *mat) {
	int i,j;
	printf("%d ", mat->n);
	for(i=0; i< mat->n; i++) {
		printf("\r\n");
		for(j=0; j<mat->n; j++) {
			printf("%f ", mat->values[(i*mat->n)+j]);
		}
	}
}

square_matrix *calculate_modularity_matrix(sparse_matrix_arr *adj_matrix, int_vector *vgroup) {
	square_matrix *mod_mat;
	int i,j, pos, mat_val_index=0;
	elem A_i_j;
	elem *F_g;
	int *K;
	if ((mod_mat = allocate_square_matrix(vgroup->count)) == NULL) {
		MEMORY_ALLOCATION_FAILURE_AT("calculate_modularity_matrix: mod_mat");
		return NULL;
	}
	if ((F_g = calloc(vgroup->count, sizeof(elem))) == NULL) {
		MEMORY_ALLOCATION_FAILURE_AT("calculate_modularity_matrix: F_g");
		free(mod_mat);
		return NULL;
	}
	if ((K = calloc(vgroup->count, sizeof(int))) == NULL) {
		MEMORY_ALLOCATION_FAILURE_AT("calculate_modularity_matrix: K");
		free(mod_mat);
		free(F_g);
		return NULL;
	}
	/*Calculate and store degree of vertices*/
	for(i=0; i<vgroup->count; i++) {
		K[i] = degree_of_vertice(vgroup->vertices[i], adj_matrix);
	}
	for(i=0; i<vgroup->count; i++) {
		F_g[i] = 0;
		mat_val_index = adj_matrix->rowptr[vgroup->vertices[i]];
		for(j=0; j<vgroup->count; j++) {
			if (mat_val_index < adj_matrix->rowptr[vgroup->vertices[i]+1]) {
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
			mod_mat->values[pos] = A_i_j - (((elem)(K[i]*K[j]))/adj_matrix->rowptr[adj_matrix->n]);
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

elem_vector *allocate_elem_vector(int length) {
	elem_vector *vector;
	if((vector = malloc(sizeof(elem_vector))) == NULL) {
		/* Error allocating  */
		/*TODO: handle error*/
		return NULL;
	}
	vector->n = length;
	if ((vector->values = calloc(length, sizeof(elem))) == NULL) {
		/* Error allocating  */
		/*TODO: handle error*/
		return NULL;
	}
	return vector;
}

void free_elem_vector(elem_vector *vector) {
	free(vector->values);
	free(vector);
}

elem_vector *mat_vec_multiply(square_matrix *mat, elem_vector *vec) {
	int i, j;
	elem_vector *result;
	if ((result = allocate_elem_vector(vec->n)) == NULL) {
		/* Error allocating  */
		/*TODO: handle error*/
		return NULL;
	}
	for(i=0; i<mat->n; i++) {
		vec->values[i] = 0;
		for(j=0; j<mat->n; j++) {
			result->values[i] = result->values[i] + (mat->values[(i*mat->n) + j]*vec->values[j]);
		}
	}
	return result;
}

eigen_pair *calculate_leading_eigen_pair(square_matrix *mod_mat) {
	elem_vector *X, *X_next; /* represent X[0], x[1]... from the algorithm */
	int i;
	if ((X = allocate_elem_vector(mod_mat->n)) == NULL) {
		/* Error allocating  */
		/*TODO: handle error*/
		return NULL;
	}
	X->values[0] = 1;
	for (i=1; i<X->n; i++) {
		X->values[i] = 0;
	}
	X_next = mat_vec_multiply(mod_mat, X);
	print_elem_vector(X_next->values, X_next->n);
	return NULL;
}
