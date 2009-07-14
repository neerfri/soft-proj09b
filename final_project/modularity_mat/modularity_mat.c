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
	FILE *fp;
	sparse_matrix_arr* adj_matrix;
	int_vector *vgroup;
	square_matrix *modularity_matrix;
	if (argc < 3) {
		fprintf(stderr, "Invalid Arguments, Aborting.\n");
		fprintf(stderr, "Usage: %s <adjacency-mat-file> <group-file>\n", argv[0]);
		return EXIT_FAILURE;
	}
	if ((fp = fopen(argv[1], "r")) == NULL) {
		/*File open error ! abort */
		fprintf(stderr, "Could not open adjacency matrix file: '%s'. Aborting.\n", argv[1]);
		return EXIT_FAILURE;
	}
	if ((adj_matrix = allocate_and_read_matrix(fp)) == NULL) {
		fclose(fp);
		exit(EXIT_FAILURE);
	}
	fclose(fp);
	if ((fp = fopen(argv[2], "r")) == NULL) {
		/*File open error ! abort */
		fprintf(stderr, "Could not open vertices file: '%s'. Aborting.\n", argv[2]);
		free(adj_matrix);
		return EXIT_FAILURE;
	}
	if ((vgroup = malloc(sizeof(int_vector))) == NULL) {
		MEMORY_ALLOCATION_FAILURE_AT("main: vgroup");
		free(adj_matrix);
		return EXIT_FAILURE;
	}
	read_n_vertices_group(fp, vgroup, adj_matrix->n);
	fclose(fp);
	modularity_matrix = calculate_modularity_matrix(adj_matrix, vgroup);
	print_square_matrix(modularity_matrix);
	return EXIT_SUCCESS;
}
