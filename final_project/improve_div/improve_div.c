#include "default_includes.h"

int main(int argc, char **argv) {
	sparse_matrix_arr* adj_matrix;
	int_vector *vgroup;
	mod_matrix *Bijtag;
	eigen_pair *leading_eigen_pair;
	two_division *division;
	double precision;
	int i;
	if (argc < 4) {
		fprintf(stderr, "Invalid Arguments, Aborting.\n");
		fprintf(stderr, "Usage: %s <adjacency-mat-file> <group-file> <precision>\n", argv[0]);
		return EXIT_FAILURE;
	}
	if (sscanf(argv[3], "%lf", &precision) < 1) {
		fprintf(stderr, "Invalid precision, Aborting.\n");
		fprintf(stderr, "Usage: %s <adjacency-mat-file> <group-file> <precision>\n", argv[0]);
		return EXIT_FAILURE;
	}
	if ((adj_matrix = read_adjacency_matrix_file(argv[1])) == NULL) {
		/*Problem reading adjacency matrix data */
		return EXIT_FAILURE;
	}
	if ((vgroup = read_vertices_group_file(argv[2], adj_matrix->n)) == NULL) {
		/*Problem reading adjacency matrix data */
		free_sparse_matrix_arr(adj_matrix);
		return EXIT_FAILURE;
	}
	if ((Bijtag = allocate_partial_modularity_matrix(adj_matrix, vgroup)) == NULL) {
		free_int_vector(vgroup);
		free_sparse_matrix_arr(adj_matrix);
		return EXIT_FAILURE;
	}
	free_sparse_matrix_arr(adj_matrix);
	if ((leading_eigen_pair = calculate_leading_eigen_pair(Bijtag, precision)) == NULL) {
		/* Failed calculating leading eigen pair */
		free_sparse_matrix_arr(adj_matrix);
		free_int_vector(vgroup);
		free_mod_matrix(Bijtag);
		return EXIT_FAILURE;
	}
	if ((division = divide_network_in_two(Bijtag, leading_eigen_pair, 1)) == NULL) {
		/* Failed calculating partition */
		free_sparse_matrix_arr(adj_matrix);
		free_int_vector(vgroup);
		free_mod_matrix(Bijtag);
		free_eigen_pair(leading_eigen_pair);
		return EXIT_FAILURE;
	}
	free_eigen_pair(leading_eigen_pair);
	printf("%f\n", division->quality);
	for(i=0; i<division->s_vector->n; i++) {
		if (division->s_vector->values[i] == division->s_vector->values[0]) {
			printf("%d ", vgroup->values[i]);
		}
	}
	printf("\n");
	for(i=0; i<division->s_vector->n; i++) {
		if (division->s_vector->values[i] != division->s_vector->values[0]) {
			printf("%d ", vgroup->values[i]);
		}
	}
	printf("\n");
	free_int_vector(vgroup);
	free_two_division(division);
	free_mod_matrix(Bijtag);
	return EXIT_SUCCESS;
}
