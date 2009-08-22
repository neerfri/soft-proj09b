#include "default_includes.h"

int main(int argc, char **argv) {
	sparse_matrix_arr* adj_matrix;
	int_vector *vgroup;
	elem *F_g;
	int *K;
	eigen_pair *leading_eigen_pair;
	double precision;
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
	if ((K = calculate_degree_of_vertices(adj_matrix, vgroup)) == NULL) {
		free_sparse_matrix_arr(adj_matrix);
		free_int_vector(vgroup);
		return EXIT_FAILURE;
	}
	if ((F_g = calculate_F_g_array(adj_matrix, vgroup, K)) == NULL) {
		free_sparse_matrix_arr(adj_matrix);
		free_int_vector(vgroup);
		free(K);
		return EXIT_FAILURE;
	}
	leading_eigen_pair = calculate_leading_eigen_pair(adj_matrix, vgroup, K, F_g, precision);
#ifdef FALSE_DEFINITION
	if ((modularity_matrix = calculate_modularity_matrix(adj_matrix, vgroup)) == NULL) {
		/* Failed calculating modularity matrix */
		free_sparse_matrix_arr(adj_matrix);
		free_int_vector(vgroup);
		return EXIT_FAILURE;
	}
	if ((leading_eigen_pair = calculate_leading_eigen_pair(modularity_matrix, precision)) == NULL) {
		/* Failed calculating leading eigen pair */
		free_sparse_matrix_arr(adj_matrix);
		free_int_vector(vgroup);
		free_square_matrix(modularity_matrix);
		return EXIT_FAILURE;
	}
#endif
	printf("%f\n", leading_eigen_pair->value);
	print_elem_vector(leading_eigen_pair->vector->values, leading_eigen_pair->vector->n);
	printf("\n");
	return EXIT_SUCCESS;
}
