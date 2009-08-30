#include "default_includes.h"

int main(int argc, char **argv) {
	sparse_matrix_arr* adj_matrix;
	int_vector *vgroup, *subgroup;
	elem_vector *s_vector;
	mod_matrix *Bijtag;
	two_division *division;
	int i,j;
	if (argc < 4) {
		fprintf(stderr, "Invalid Arguments, Aborting.\n");
		fprintf(stderr, "Usage: %s <adjacency-mat-file> <group-file> <first-subgroup-file>\n", argv[0]);
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
	if ((subgroup = read_vertices_group_file(argv[3], adj_matrix->n)) == NULL) {
		/*Problem reading adjacency matrix data */
		free_int_vector(vgroup);
		free_sparse_matrix_arr(adj_matrix);
		return EXIT_FAILURE;
	}
	if ((Bijtag = allocate_partial_modularity_matrix(adj_matrix, vgroup)) == NULL) {
		free_int_vector(vgroup);
		free_sparse_matrix_arr(adj_matrix);
		return EXIT_FAILURE;
	}
	free_sparse_matrix_arr(adj_matrix);
	if ((s_vector = allocate_elem_vector(vgroup->n)) == NULL) {
		free_mod_matrix(Bijtag);
		free_int_vector(vgroup);
		return EXIT_FAILURE;
	}
	if ((division = allocate_two_division(s_vector)) == NULL) {
		free_elem_vector(s_vector);
		free_mod_matrix(Bijtag);
		free_int_vector(vgroup);
		return EXIT_FAILURE;
	}

	/*align subgroup indices to be a subgroup of vgroup indices*/
	for(i=0, j=0; i<subgroup->n; i++) {
		while(j<vgroup->n && vgroup->values[j]!=subgroup->values[i])
			j++;
		if (j==vgroup->n) {
			fprintf(stderr, "vertex %d not found in main group file.\n", subgroup->values[i]);
		}
		if (vgroup->values[j]==subgroup->values[i]) {
			subgroup->values[i] = j;
		}
	}

	for(i=0; i<s_vector->n; i++) {
		s_vector->values[i] = 1;
	}
	for(i=0; i<subgroup->n; i++) {
		s_vector->values[subgroup->values[i]] = -1;
	}
	division->quality = 0;
	if(improve_network_division(Bijtag, division) == 0) {
		free_int_vector(vgroup);
		free_int_vector(subgroup);
		free_two_division(division);
		free_mod_matrix(Bijtag);
		return EXIT_FAILURE;
	}
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
	free_int_vector(subgroup);
	free_int_vector(vgroup);
	free_two_division(division);
	free_mod_matrix(Bijtag);
	return EXIT_SUCCESS;
}
