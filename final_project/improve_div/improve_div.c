#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "shared.h"

int main(int argc, char **argv) {
	sparse_matrix_arr* adj_matrix;
	int_vector *vgroup;
	square_matrix *modularity_matrix;
	eigen_pair *leading_eigen_pair;
	two_division *division;
	int i;
	if (argc < 4) {
		fprintf(stderr, "Invalid Arguments, Aborting.\n");
		fprintf(stderr, "Usage: %s <adjacency-mat-file> <group-file> <precision>\n", argv[0]);
		return EXIT_FAILURE;
	}
	if ((adj_matrix = read_adjacency_matrix(argv[1])) == NULL) {
		/*Problem reading adjacency matrix data */
		return EXIT_FAILURE;
	}
	if ((vgroup = read_vertices_group_file(argv[2], adj_matrix->n)) == NULL) {
		/*Problem reading adjacency matrix data */
		free_sparse_matrix_arr(adj_matrix);
		return EXIT_FAILURE;
	}
	if ((modularity_matrix = calculate_modularity_matrix(adj_matrix, vgroup)) == NULL) {
		/* Failed calculating modularity matrix */
		free_sparse_matrix_arr(adj_matrix);
		free_int_vector(vgroup);
		return EXIT_FAILURE;
	}
	if ((leading_eigen_pair = calculate_leading_eigen_pair(modularity_matrix)) == NULL) {
		/* Failed calculating leading eigen pair */
		free_sparse_matrix_arr(adj_matrix);
		free_int_vector(vgroup);
		free_square_matrix(modularity_matrix);
		return EXIT_FAILURE;
	}
	if ((division = divide_network_in_two(modularity_matrix, leading_eigen_pair)) == NULL) {
		/* Failed calculating partition */
		free_sparse_matrix_arr(adj_matrix);
		free_int_vector(vgroup);
		free_square_matrix(modularity_matrix);
		free_eigen_pair(leading_eigen_pair);
		return EXIT_FAILURE;
	}
	division->division->vertices[4] = division->division->vertices[4] * -1;
	improve_network_division(modularity_matrix, division);
	printf("%f\n", division->quality);
	for(i=0; i<division->division->n; i++) {
		if (division->division->vertices[i] == division->division->vertices[0]) {
			printf("%d ", i);
		}
	}
	printf("\n");
	for(i=0; i<division->division->n; i++) {
		if (division->division->vertices[i] != division->division->vertices[0]) {
			printf("%d ", vgroup->vertices[i]);
		}
	}
	/*
	printf("%f\n", leading_eigen_pair->value);
	print_elem_vector(leading_eigen_pair->vector->values, leading_eigen_pair->vector->n);
	*/
	printf("\n");
	return EXIT_SUCCESS;
}
