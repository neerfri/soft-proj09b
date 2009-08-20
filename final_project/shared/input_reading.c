/*
 * input_reading.c
 *
 *  Created on: Aug 19, 2009
 *      Author: neer
 */

#include "input_reading.h"

elem *allocate_and_read_vector(FILE *fp, int n) {
	int i;
	double scanf_receptor;
	elem *vector;
	if ((vector = allocate_vector(n)) == NULL) {
		return NULL;
	}
	for(i=0; i<n; i++) {
		if (fscanf(fp, "%lf", &scanf_receptor) < 1) {
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
	matrix_cell_link head, tail, temp_cell;
	sparse_matrix_arr* matrix;
	if (fscanf(fp, "%d", &n) < 1) {
		fprintf(stderr, "Error reading matrix size");
		return NULL;
	}
	head = NULL; tail = NULL;
	for(j=0; j<n; j++) {
		for(i=0; i<n; i++) {
			pos = i + (j*n);
			if (fscanf(fp, "%lf", &scanf_receptor) < 1) {
				fprintf(stderr, "Error reading matrix values at pos: %d\n", pos);
				/*free(values);*/
				return NULL;
			}
			if (scanf_receptor != 0) {
				numNnz++;
				if ((temp_cell = new_matrix_cell(scanf_receptor, j, i)) == NULL) {
					/* Error in adding value to list, need to free memory and abort */
					/*TODO: free memory */
				}
				if (tail == NULL) {
					tail = temp_cell;
					head = temp_cell;
				}
				else {
					tail->next = temp_cell;
					tail = tail->next;
				}
			}

		}
	}

	if ((matrix = allocate_sparse_matrix_arr(n, numNnz)) == NULL)
		return NULL;
	val_index = 0;
	matrix->rowptr[0] = 0;
	temp_cell = head;
	for(j=0; j<n; j++) {
		while(temp_cell != NULL && temp_cell->row == j) {
			matrix->values[val_index] = temp_cell->value;
			matrix->colind[val_index] = temp_cell->col;
			temp_cell = temp_cell->next;
			val_index++;
		}
		matrix->rowptr[j+1] = val_index;
	}
	free_matrix_cell_list(head);
	return matrix;
}

sparse_matrix_arr *read_adjacency_matrix_file(const char* file) {
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
	int_list_link head, ptr;
	int scanf_receptor, i;
	head = NULL;
	for(i=0; i < n; i++) {
		if (feof(fp) || fscanf(fp, "%d", &scanf_receptor) < 1) {
			break;
		} else {
			if (head == NULL) {
				/* this is the first value, initialize head pointer*/
				if ((head = malloc(sizeof(int_list_element))) == NULL) {
					MEMORY_ALLOCATION_FAILURE_AT("read_vertices_group: head");
					return -1;
				}
				ptr = head;
			} else {
				/* this is not the first value, add an element to the list */
				if ((ptr->next = malloc(sizeof(int_list_element))) == NULL) {
					MEMORY_ALLOCATION_FAILURE_AT("read_vertices_group: ptr->next");
					free_int_list(head);
					return -1;
				}
				ptr = ptr->next;
			}
			ptr->value = scanf_receptor;
			ptr->next = NULL;
		}
	}
	vertices->n = i;
	if ((vertices->values = calloc(vertices->n, sizeof(int))) == NULL) {
		MEMORY_ALLOCATION_FAILURE_AT("read_vertices_group: vertices->values");
		free_int_list(head);
		return -1;
	}
	ptr = head;
	for(i=0; i<vertices->n; i++) {
		vertices->values[i] = ptr->value;
		ptr = ptr->next;
	}
	free_int_list(head);
	return vertices->n;
}


