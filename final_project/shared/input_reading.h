/*
 * input_reading.h
 *
 *  Created on: Aug 19, 2009
 *      Author: neer
 */

#ifndef INPUT_READING_H_
#define INPUT_READING_H_

#include "default_includes.h"

sparse_matrix_arr *read_adjacency_matrix_file(const char* file);
int_vector *read_vertices_group_file(const char* file, int max_count);

sparse_matrix_arr *allocate_and_read_matrix(FILE *fp);
elem *allocate_and_read_vector(FILE *fp, int n);
int read_n_vertices_group(FILE *fp, int_vector *vertices, int n);


#endif /* INPUT_READING_H_ */
