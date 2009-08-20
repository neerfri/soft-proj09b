/*
 * data_structures.c
 *
 *  Created on: Aug 20, 2009
 *      Author: neer
 */

#include "data_structures.h"

void free_matrix_cell_list(matrix_cell_link head) {
	matrix_cell_link ptr;
	while(head!=NULL) {
		ptr = head;
		head = head->next;
		free(head);
	}
}


void free_int_list(int_list_link head) {
	int_list_link ptr;
	while(head!=NULL) {
		ptr = head;
		head = head->next;
		free(head);
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
