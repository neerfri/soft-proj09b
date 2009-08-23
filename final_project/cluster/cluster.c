#include "default_includes.h"

int main(int argc, char **argv) {
	sparse_matrix_arr* adj_matrix;
	double precision;
	two_division *final_division;
	if (argc < 3) {
		fprintf(stderr, "Invalid Arguments, Aborting.\n");
		fprintf(stderr, "Usage: %s <adjacency-mat-file> <group-file> <precision>\n", argv[0]);
		return EXIT_FAILURE;
	}
	if (sscanf(argv[2], "%lf", &precision) < 1) {
		fprintf(stderr, "Invalid precision, Aborting.\n");
		fprintf(stderr, "Usage: %s <adjacency-mat-file> <group-file> <precision>\n", argv[0]);
		return EXIT_FAILURE;
	}
	if ((adj_matrix = read_adjacency_matrix_file(argv[1])) == NULL) {
		/*Problem reading adjacency matrix data */
		return EXIT_FAILURE;
	}
	if ((final_division =  algorithm3(adj_matrix, precision, 1)) == NULL) {
		free_sparse_matrix_arr(adj_matrix);
		return EXIT_FAILURE;
	}
	print_clusters(final_division);
	return EXIT_SUCCESS;
}
