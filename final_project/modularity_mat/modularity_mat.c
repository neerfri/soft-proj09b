/*TODO:
 * validate the vertices group (values should be in range of matrix size and in increasing order)
 * catch when calculate_modularity_matrix returns NULL
 * */

#include "default_includes.h"

int main(int argc, char **argv) {
	sparse_matrix_arr* adj_matrix;
	int_vector *vgroup;
	mod_matrix *Bijtag;
	if (argc < 3) {
		fprintf(stderr, "Invalid Arguments, Aborting.\n");
		fprintf(stderr, "Usage: %s <adjacency-mat-file> <group-file>\n", argv[0]);
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
	printf("%d \r\n", vgroup->n);
	print_modularity_matrix(Bijtag);
	free_sparse_matrix_arr(adj_matrix);
	free_int_vector(vgroup);
	return EXIT_SUCCESS;
}
