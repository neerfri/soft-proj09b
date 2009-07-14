/*TODO:
 * validate the vertices group (values should be in range of matrix size and in increasing order)
 * catch when calculate_modularity_matrix returns NULL
 * */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "sparse_matrix.h"
#include "shared.h"

int main(int argc, char **argv) {
	sparse_matrix_arr* adj_matrix;
	int_vector *vgroup;
	square_matrix *modularity_matrix;
	if (argc < 3) {
		fprintf(stderr, "Invalid Arguments, Aborting.\n");
		fprintf(stderr, "Usage: %s <adjacency-mat-file> <group-file>\n", argv[0]);
		return EXIT_FAILURE;
	}
	if ((adj_matrix = read_adjacency_matrix(argv[1])) == NULL) {
		/*Problem reading adjacency matrix data */
		return EXIT_FAILURE;
	}
	if ((vgroup = read_vertices_group_file(argv[2], adj_matrix->n)) == NULL) {
		/*Problem reading adjacency matrix data */
		free(adj_matrix);
		return EXIT_FAILURE;
	}
	modularity_matrix = calculate_modularity_matrix(adj_matrix, vgroup);
	print_square_matrix(modularity_matrix);
	return EXIT_SUCCESS;
}
