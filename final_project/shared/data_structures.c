/*
 * data_structures.c
 *
 *  Created on: Aug 20, 2009
 *      Author: neer
 */

#include "data_structures.h"

matrix_cell_link new_matrix_cell(elem value, int row, int col) {
	matrix_cell_element *new_cell;
	if ((new_cell = malloc(sizeof(matrix_cell_element))) == NULL) {
		MEMORY_ALLOCATION_FAILURE_AT("new_matrix_cell: new_cell");
		return NULL;
	}
	new_cell->value = value;
	new_cell->row = row;
	new_cell->col = col;
	new_cell->next = NULL;
	return new_cell;
}

void free_matrix_cell_list(matrix_cell_link head) {
	matrix_cell_link ptr;
	while(head!=NULL) {
		ptr = head;
		head = head->next;
		free(ptr);
	}
}

sparse_matrix_arr *create_sparse_matrix_from_list(matrix_cell_link head, int n, int numNnz) {
	sparse_matrix_arr *matrix;
	int mat_index, j;
	matrix_cell_link cell_ptr;
	if ((matrix = allocate_sparse_matrix_arr(n, numNnz)) == NULL)
		return NULL;
	mat_index = 0;
	matrix->rowptr[0] = 0;
	cell_ptr = head;
	for(j=0; j<n; j++) {
		while(cell_ptr != NULL && cell_ptr->row == j) {
			matrix->values[mat_index] = cell_ptr->value;
			matrix->colind[mat_index] = cell_ptr->col;
			cell_ptr = cell_ptr->next;
			mat_index++;
		}
		matrix->rowptr[j+1] = mat_index;
	}
	return matrix;
}


void free_int_list(int_list_link head) {
	int_list_link temp;
	while(head!=NULL) {
		temp = head;
		head = head->next;
		free(temp);
	}
}

int_vector *allocate_int_vector(int n) {
	int_vector *result;
	if ((result = malloc(sizeof(int_vector))) == NULL) {
		MEMORY_ALLOCATION_FAILURE_AT("allocate_int_vector: result");
		return NULL;
	}
	result->n = n;
	if ((result->values = calloc(n, sizeof(int))) == NULL) {
		MEMORY_ALLOCATION_FAILURE_AT("allocate_int_vector: result->values");
		return NULL;
	}
	return result;
}

void free_int_vector(int_vector *vector) {
	free(vector->values);
	free(vector);
}

elem_vector *allocate_elem_vector(int length) {
	elem_vector *vector;
	if((vector = malloc(sizeof(elem_vector))) == NULL) {
		/* Error allocating  */
		MEMORY_ALLOCATION_FAILURE_AT("allocate_elem_vector: vector");
		/*TODO: handle error*/
		return NULL;
	}
	vector->n = length;
	if ((vector->values = calloc(length, sizeof(elem))) == NULL) {
		/* Error allocating  */
		MEMORY_ALLOCATION_FAILURE_AT("allocate_elem_vector: vector->values");
		/*TODO: handle error*/
		return NULL;
	}
	return vector;
}

void free_elem_vector(elem_vector *vector) {
	free(vector->values);
	free(vector);
}

void free_mod_matrix(mod_matrix *mat);

sparse_matrix_arr *get_partial_sparse_matrix(sparse_matrix_arr *mat, int_vector *group) {
	int i,j, mat_index;
	int numNnz=0; /* number of values found for the new matrix */
	matrix_cell_link head=NULL, tail=NULL, temp_cell=NULL;
	int_vector *reverse_group; /*this holds a reverse mapping to group->values[i], negative in case non existent*/
	sparse_matrix_arr *result;
	if ((reverse_group = allocate_int_vector(mat->n)) == NULL) {
		return NULL;
	}
	/* initialize reverse_vgroup */
	for(i =0, j=0; i < mat->n; i++) {
		if(group->values[j] == i) {
			reverse_group->values[i] = j;
			j++;
		} else
			reverse_group->values[i] = -1;
	}
	/*the following loop creates a list of cells needed for the new matrix */
	for(i=0; i<group->n; i++) {
		mat_index = mat->rowptr[group->values[i]];
		while(mat_index < mat->rowptr[group->values[i]+1]) {
			j = reverse_group->values[mat->colind[mat_index]]; /* the column in terms of new mat*/
			if(j >= 0) {
				/* current cell should be in the partial matrix */
				if ((temp_cell = new_matrix_cell(mat->values[mat_index], i, j)) == NULL) {
					free_int_vector(reverse_group);
					free_matrix_cell_list(head);
					return NULL;
				}
				if (tail == NULL) {
					tail = temp_cell;
					head = temp_cell;
				}
				else {
					tail->next = temp_cell;
					tail = tail->next;
				}
				numNnz++;
			}
			mat_index++;
		}
	}
	result = create_sparse_matrix_from_list(head, group->n, numNnz);
	free_int_vector(reverse_group);
	free_matrix_cell_list(head);
	return result;
}

eigen_pair *allocate_eigen_pair(elem value, elem_vector *vector) {
	eigen_pair *result;
	if ((result = malloc(sizeof(eigen_pair))) == NULL) {
		/* Error allocating  */
		MEMORY_ALLOCATION_FAILURE_AT("allocate_eigen_pair: result");
		/*TODO: handle error*/
		return NULL;
	}
	result->value = value;
	result->vector = vector;
	return result;
}

void free_eigen_pair(eigen_pair *pair) {
	free_elem_vector(pair->vector);
	free(pair);
}



two_division *allocate_two_division(int n) {
	two_division *result;
	if ((result = malloc(sizeof(two_division))) == NULL) {
		MEMORY_ALLOCATION_FAILURE_AT("allocate_two_division: result");
		return NULL;
	}
	if ((result->division = allocate_int_vector(n)) == NULL) {
		free(result);
		return NULL;
	}
	return result;
}

void free_two_division(two_division *division) {
	free_int_vector(division->division);
	free(division);
}
