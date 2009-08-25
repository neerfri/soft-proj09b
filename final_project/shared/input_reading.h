/*
 * input_reading.h
 *
 *  Created on: Aug 19, 2009
 *      Author: neer
 */

#ifndef INPUT_READING_H_
#define INPUT_READING_H_

#include "default_includes.h"

/* read adjacency matrix from a file with path 'file'
 */
sparse_matrix_arr *read_adjacency_matrix_file(const char* file);
/* read vertices groups from a file with path 'file'
 */
int_vector *read_vertices_group_file(const char* file, int max_count);

/* reads a matrix from fp FILE pointer */
sparse_matrix_arr *allocate_and_read_matrix(FILE *fp);

/* reads a double valued vector of length n from file fp
 */
elem *allocate_and_read_vector(FILE *fp, int n);

/* reads up to 'n' vertex indices from fp FILE pointer
 * vertex indices should be presented in an increasing order
 */
int read_n_vertices_group(FILE *fp, int_vector *vertices, int n);


#endif /* INPUT_READING_H_ */
