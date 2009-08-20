/*
 * data_structures.h
 *
 *  Created on: Aug 20, 2009
 *      Author: neer
 */

#ifndef DATA_STRUCTURES_H_
#define DATA_STRUCTURES_H_

#include "sparse_matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

/*
 * A structure to hold a single matrix value with it's position in the matrix
 * this structure is used in the sparse matrix construction process
 */
struct matrix_cell_struct
{
  int    row; /* pointers to where rows begin in the values array. */
  int    col; /* column indices */
  elem   value;
  struct matrix_cell_struct *next;
};
typedef struct matrix_cell_struct matrix_cell_element;
typedef matrix_cell_element *matrix_cell_link;

matrix_cell_link new_matrix_cell(elem value, int row, int col);
void free_matrix_cell_list(matrix_cell_link head);



/*
 * This structure holds a list of integers
 * it is used to read groups of indexes when the group size is unknown
 */
struct int_list_element_struct
{
	int value;
	struct int_list_element_struct *next;
};
typedef struct int_list_element_struct int_list_element;
typedef int_list_element *int_list_link;

void free_int_list(int_list_link head);


/* Integer array with stored length
 * This structure is used for holding vertices groups.
 **/
typedef struct {
	int	n;    /* size                  */
	int *values;
} int_vector;

int_vector *allocate_int_vector(int n);
void free_int_vector(int_vector *vector);


/* Holds a vector of elements of type elem
 */
typedef struct {
	int n;			/* vector's length */
	elem *values;
} elem_vector;

elem_vector *allocate_elem_vector(int length);
void free_elem_vector(elem_vector *vector);


/*
 * Holds eigen-value and eigen-vector pair
 */
typedef struct {
	elem value; /* the eigen value */
	elem_vector *vector; /* the eigen vector */
} eigen_pair;

eigen_pair *allocate_eigen_pair(elem value, elem_vector *vector);
void free_eigen_pair(eigen_pair *pair);



typedef struct {
	elem quality;
	int_vector *division;
} two_division;

two_division *allocate_two_division(int n);
void free_two_division(two_division *division);
#endif /* DATA_STRUCTURES_H_ */
