#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "shared.h"

int main(int argc, char **argv) {
	sparse_matrix_arr* adj_matrix;
	int_vector *vgroup;
	square_matrix *modularity_matrix;
	eigen_pair *leading_eigen_pair;
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
	printf("%f\n", leading_eigen_pair->value);
	print_elem_vector(leading_eigen_pair->vector->values, leading_eigen_pair->vector->n);
	printf("\n");
	return EXIT_SUCCESS;
}
