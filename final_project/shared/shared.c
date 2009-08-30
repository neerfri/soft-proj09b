#include "default_includes.h"

void print_elem_vector(elem *vector, int n) {
	int i;
	for(i=0; i < n;i++) {
		printf("%f ", vector[i]);
	}
}

void print_int_vector(int *vector, int n) {
	int i;
	for(i=0; i < n;i++) {
		printf("%d ", vector[i]);
	}
}

int degree_of_vertice(int i, sparse_matrix_arr *matrix) {
	return matrix->rowptr[i+1]-matrix->rowptr[i];
}


/*
 * Computes the value at row i column j of the modularity
 * matrix (non generalized one!!!) for a given vertices group
 */
elem modularity_matrix_cell\
(mod_matrix *mod_mat, int i, int j) {
	elem A_i_j, B_i_j;
	/* An index to scan the sparse matrix values*/
	int mat_val_index;

	/* find the corresponding A_i_j element in the sparse matrix */
	/*jump to the beginning of the row: */
	mat_val_index = mod_mat->A_g->rowptr[i];
	/*search for the item with colind>=j in that row*/
	while(mat_val_index < mod_mat->A_g->rowptr[i+1] &&
			mod_mat->A_g->colind[mat_val_index] < j) {
		mat_val_index++;
	}
	if (mat_val_index < mod_mat->A_g->rowptr[i+1] && \
			mod_mat->A_g->colind[mat_val_index] == j) {
		/* we found the cell value */
		A_i_j = mod_mat->A_g->values[mat_val_index];
	} else {
		/*cell is not in sparse matrix, hence 0*/
		A_i_j = 0;
	}
	/* next lines calculates: B_i_j = A_i_j - (K_i*K_j)/M */
	B_i_j = A_i_j - ((mod_mat->K->values[i]*mod_mat->K->values[j])/mod_mat->total_degree);
	return B_i_j;
}

/*
 * Computes the value at row i column j of the generalized modularity
 * matrix for a given vertices group
 */
elem generalized_modularity_matrix_cell(mod_matrix *mod_mat, int i, int j) {
	elem B_i_j = modularity_matrix_cell(mod_mat, i, j);
	if (i == j)
		return B_i_j - mod_mat->f_g->values[i];
	else
		return B_i_j;
}

void print_modularity_matrix(mod_matrix *mod_mat) {
	int i,j;
	for(i=0; i<mod_mat->A_g->n; i++) {
		for(j=0; j<mod_mat->A_g->n; j++) {
			printf("%f ", generalized_modularity_matrix_cell(mod_mat, i, j));
		}
		printf("\n");
	}
}

/* calculates the first norm of mod_mat
 */
elem calculate_matrix_first_norm(mod_matrix *mod_mat) {
	int new_j,i;
	elem result=0, accumulator;
	elem val_at_i_j;
	for(i=0; i<mod_mat->A_g->n; i++) {
		accumulator = 0;
		for(new_j=0; new_j<mod_mat->A_g->n; new_j++) {
			val_at_i_j = generalized_modularity_matrix_cell(mod_mat, new_j, i);
			accumulator = accumulator + fabs(val_at_i_j);
		}
		if (accumulator > result)
			result = accumulator;
	}
	return result;
}

/* Calculates an array storing the degree of all vertices in 'vertices_group'
 * from the adjacency matrix  'adj_matrix'
 */
elem_vector *calculate_degree_of_vertices\
(sparse_matrix_arr *adj_matrix, int_vector *vertices_group) {
	elem_vector *result;
	int i;
	if ((result = allocate_elem_vector(vertices_group->n)) == NULL) {
		return NULL;
	}
	/*Calculate and store degree of vertices*/
	for(i=0; i<vertices_group->n; i++) {
		result->values[i] = degree_of_vertice(vertices_group->values[i], adj_matrix);
	}
	return result;
}

elem_vector *calculate_F_g_array(mod_matrix *mod_mat) {
	elem_vector *result;
	elem A_i_j;
	int i,j;
	int mat_index; /* An index to scan the sparse matrix values*/
	if ((result = allocate_elem_vector(mod_mat->A_g->n)) == NULL) {
		return NULL;
	}
	for(i=0; i<mod_mat->A_g->n; i++) {
		result->values[i] = 0;
		mat_index = mod_mat->A_g->rowptr[i];
		for(j=0; j<mod_mat->A_g->n; j++) {
			if (mat_index < mod_mat->A_g->rowptr[i+1]) {
				/* we are still in values for row i; */
				/* first we need to advance the sparse matrix pointer to at least column j*/
				while(mat_index < mod_mat->A_g->rowptr[i+1] &&
						mod_mat->A_g->colind[mat_index] < j) {
					mat_index++;
				}
				/* if its in adj_matrix values, take it, otherwise it's 0*/
				/* first condition makes sure we don't overflow from the array boundaries */
				if (mat_index < mod_mat->A_g->rowptr[mod_mat->A_g->n] &&
						mod_mat->A_g->colind[mat_index] == j) {
					A_i_j = mod_mat->A_g->values[mat_index];
					mat_index++;
				} else {
					A_i_j = 0;
				}
			} else {
				/* there are no more values for this row, meaning value is 0 */
				A_i_j = 0;
			}
			/* next two lines calculates: F_g[i] = Sum(Over all j in vgroup)[B[g]_i_j]*/
			result->values[i] += A_i_j;
			result->values[i] -=((mod_mat->K->values[i]*mod_mat->K->values[j])/
					mod_mat->total_degree);
		}
	}
	return result;
}

elem vec_vec_multiply(elem_vector *vec1, elem_vector *vec2) {
	elem result = 0;
	int i;
	for(i=0; i<vec1->n; i++) {
		result = result + (vec1->values[i]*vec2->values[i]);
	}
	return result;
}

/* calculate a vector norm: |vec|
 */
elem vec_norm(elem_vector *vec) {
	return (elem)sqrt(vec_vec_multiply(vec, vec));
}

elem_vector *elem_vec_multiply(elem scalar, elem_vector *vec) {
	elem_vector *result;
	int i;
	if ((result = allocate_elem_vector(vec->n)) == NULL) {
		/* Error allocating  */
		return NULL;
	}
	for(i=0; i<vec->n; i++) {
		result->values[i] = scalar * vec->values[i];
	}
	return result;
}

void vec_normalize(elem_vector *vec) {
	elem factor;
	int i;
	factor = vec_norm(vec);
	for(i=0; i<vec->n; i++) {
		vec->values[i] = vec->values[i]/factor;
	}
}

elem_vector *sparse_mat_vec_multiply(sparse_matrix_arr *A, elem_vector *v) {
	int i, j, val_index;
	elem_vector *result;
	if ((result = allocate_elem_vector(v->n)) == NULL) {
		return NULL;
	}
	/* this scans the matrix lines */
	for(i=0; i<A->n; i++) {
		/* now scan the values in that line */
		for(val_index = A->rowptr[i]; val_index<A->rowptr[i+1]; val_index++) {
			j = A->colind[val_index];
			result->values[i] += A->values[val_index]*v->values[j];
		}
	}
	return result;
}

elem_vector *pure_mod_mat_vec_multiply(mod_matrix *Bijtag, elem_vector *vec) {
	elem_vector *result;
	elem_vector *kgx_M_kg, *Agx;
	int i;
	if ((result = allocate_elem_vector(Bijtag->A_g->n)) == NULL) {
		return NULL;
	}
	/* ( (K[g]*X)/M )*K[g]  */
	if ((kgx_M_kg = elem_vec_multiply(vec_vec_multiply(Bijtag->K, vec)/Bijtag->total_degree,
			Bijtag->K)) == NULL) {
		free(result);
		return NULL;
	}
	if((Agx = sparse_mat_vec_multiply(Bijtag->A_g, vec)) == NULL) {
		free_elem_vector(kgx_M_kg);
		free(result);
		return NULL;
	}
	for (i = 0; i<Bijtag->A_g->n; i++){
		result->values[i] = Agx->values[i]-kgx_M_kg->values[i];
	}
	free_elem_vector(kgx_M_kg);
	free_elem_vector(Agx);
	return result;
}

elem_vector *gen_mod_mat_vec_multiply(mod_matrix *Bijtag, elem_vector *vec) {
	int i;
	elem_vector *result;
	if((result = pure_mod_mat_vec_multiply(Bijtag,vec)) == NULL) {
		return NULL;
	}
	for(i = 0; i<Bijtag->A_g->n; i++){
		result->values[i] -= Bijtag->f_g->values[i]*vec->values[i];
	}
	return result;
}

elem_vector *shifted_gen_mod_mat_vec_multiply(mod_matrix *Bijtag, elem_vector *vec) {
	elem_vector *result;
	int i;
	if((result = gen_mod_mat_vec_multiply(Bijtag, vec)) == NULL) {
		return NULL;
	}
	for(i = 0; i<Bijtag->A_g->n; i++){
		result->values[i]+=Bijtag->norm_1*vec->values[i];
	}
	return result;
}

eigen_pair *calculate_leading_eigen_pair(mod_matrix *Bijtag, double precision) {
	elem_vector *X=NULL, *X_next; /* represent X[0], x[1]... from the algorithm */
	elem eigen_value, convergence_meter;
	elem_vector *tmp_vector;
	eigen_pair *result;
	int i;
	if ((X_next = allocate_elem_vector(Bijtag->A_g->n)) == NULL) {
		/* Error allocating  */
		return NULL;
	}
	X_next->values[0] = 1;
	for (i=1; i<X_next->n; i++) {
		X_next->values[i] = 0;
	}
	if ((tmp_vector = allocate_elem_vector(Bijtag->A_g->n)) == NULL) {
		free_elem_vector(X_next);
		return NULL;
	}
	do {
		if (X!=NULL)
			free_elem_vector(X);
		X = X_next;
		vec_normalize(X);
		if((X_next = shifted_gen_mod_mat_vec_multiply(Bijtag,X)) == NULL) {
			free_elem_vector(X);
			return NULL;
		}

		eigen_value = vec_vec_multiply(X_next, X)/vec_vec_multiply(X,X);

		/*tmp_vector is used to calculate eq 19: |Ax - jx| */
		for (i=0;i<Bijtag->A_g->n;i++) {
			tmp_vector->values[i] = X_next->values[i];
			tmp_vector->values[i] -= eigen_value*X->values[i];
		}
		convergence_meter = vec_norm(tmp_vector)/vec_norm(X);
	} while(convergence_meter>precision);
	eigen_value -=Bijtag->norm_1;
	free_elem_vector(tmp_vector);
	if((result = allocate_eigen_pair(eigen_value, X)) == NULL) {
		free_elem_vector(X);
		free_elem_vector(X_next);
		return NULL;
	}
	free_elem_vector(X_next);
	return result;
}

/* get the s vector corresponding to 'vector'
 * meaning s[i] = sign(vector[i])
 */
elem_vector *get_s_vector_for(elem_vector *vector) {
	elem_vector *s;
	int i;
	if ((s = allocate_elem_vector(vector->n)) == NULL) {
		return NULL;
	}
	for(i=0; i<s->n; i++) {
		if (IS_POSITIVE(vector->values[i])) {
			s->values[i] = 1;
		} else {
			s->values[i] = -1;
		}
	}
	return s;
}

/* calculates 0.5*s*B[g]*s
 */
int calculate_modularity_score(mod_matrix *mod_mat,elem_vector *s, elem *result) {
	elem_vector *product;
	if ((product = gen_mod_mat_vec_multiply(mod_mat, s)) == NULL) {
		return 0;
	}
	*result = 0.5*vec_vec_multiply(s, product);
	free_elem_vector(product);
	return 1;
}

two_division *divide_network_in_two\
(mod_matrix *mod_mat, eigen_pair *leading_eigen_pair, int use_improve) {
	elem_vector *s;
	two_division *result;
	int i;
	if (IS_POSITIVE(leading_eigen_pair->value)) {
		if ((s = get_s_vector_for(leading_eigen_pair->vector)) == NULL) {
			return NULL;
		}
		if ((result = allocate_two_division(s)) == NULL) {
			free_elem_vector(s);
			return NULL;
		}
		if (use_improve != 0) {
			if(improve_network_division(mod_mat, result) == 0) {
				free_two_division(result); /*also frees s !!! */
				return NULL;
			}
		}
		if (calculate_modularity_score(mod_mat, s, &result->quality) == 0) {
			free_two_division(result); /*also frees s !!! */
			return NULL;
		}
		if (IS_POSITIVE(result->quality)) {
			return result;
		}
		free_two_division(result); /*also frees s !!! */
	}
	/*return the trivial division*/
	if ((s = allocate_elem_vector(mod_mat->A_g->n)) == NULL) {
		return NULL;
	}
	if ((result = allocate_two_division(s)) == NULL) {
		free_elem_vector(s);
		return NULL;
	}
	for(i=0; i<result->s_vector->n; i++) {
		result->s_vector->values[i] = 1;
	}
	return result;
}

/*
 * init's the division score vector
*/
int init_division_score(mod_matrix *mod_mat, elem_vector *s_vector, elem *score) {
	int i;
	elem_vector *x_vector;
	if ((x_vector = pure_mod_mat_vec_multiply(mod_mat,s_vector)) == NULL) {
		return 0;
	}
	for(i=0; i<mod_mat->A_g->n; i++) {
		score[i] = -2 * s_vector->values[i]*x_vector->values[i];
		score[i] -= 2 * mod_mat->K->values[i]*mod_mat->K->values[i]/mod_mat->total_degree;
	}
	free_elem_vector(x_vector);
	return 1;
}

void update_division_score\
(mod_matrix *mod_mat, elem_vector *s_vector, int max_index, elem *score) {
	int i;
	for (i=0;i<mod_mat->A_g->n;i++) {
		if (i==max_index)
			score[i]= -1*score[i];
		else
			score[i]-= 4 * s_vector->values[i]*s_vector->values[max_index]*
						modularity_matrix_cell(mod_mat, i, max_index);
	}

}


int improve_network_division(mod_matrix *mod_mat, two_division *division) {
	int *unmoved;
	int *indices;
	int i, j, i_tag, j_tag;
	elem delta_Q;
	elem *score, *improve;
	if ((score = calloc(mod_mat->A_g->n, sizeof(elem))) == NULL) {
		MEMORY_ALLOCATION_FAILURE_AT("improve_network_division: score");
		return 0;
	}
	if ((improve = calloc(mod_mat->A_g->n, sizeof(elem))) == NULL) {
		MEMORY_ALLOCATION_FAILURE_AT("improve_network_division: improve");
		free(score);
		return 0;
	}
	if ((indices = calloc(mod_mat->A_g->n, sizeof(int))) == NULL) {
		MEMORY_ALLOCATION_FAILURE_AT("improve_network_division: indices");
		free(score);
		free(improve);
		return 0;
	}
	if ((unmoved = calloc(mod_mat->A_g->n, sizeof(int))) == NULL) {
		MEMORY_ALLOCATION_FAILURE_AT("improve_network_division: unmoved");
		free(indices);
		free(score);
		free(improve);
		return 0;
	}
	do {
		/* Initialize values in unmoved */
		for(i=0; i<mod_mat->A_g->n; i++) {
			unmoved[i] = 1;
		}

		/* initialize score values */
		if (init_division_score(mod_mat, division->s_vector, score) == 0) {
			free(unmoved);
			free(indices);
			free(score);
			free(improve);
			return 0;
		}



		/* trying to find an improvement of the partition defined by s */
		for(i=0; i<mod_mat->A_g->n; i++) {

			/* find j (ptr_to_j_tag->value) with maximal score[j] */
			for(j_tag = -1, j=0; j<mod_mat->A_g->n; j++) {
				if (unmoved[j]) {
					if (j_tag==-1 || score[j_tag] < score[j]) {
						j_tag = j;
					}
				}
			}


			division->s_vector->values[j_tag] *= -1; /* S[j_tag] = -S[j_tag] */
			indices[i] = j_tag;
			improve[i] = score[j_tag];
			if (i)
				improve[i] += improve[i-1];
			unmoved[j_tag] = 0;

			/* update scores */
			update_division_score(mod_mat, division->s_vector, j_tag, score);
		}
		/* Find the maximum improvement of s and update s accordingly */
		for(i=0, i_tag=0; i<mod_mat->A_g->n; i++)
			if(improve[i_tag]<improve[i])
				i_tag = i;

		for(i=mod_mat->A_g->n-1; i>i_tag; i--)
			division->s_vector->values[indices[i]] *= -1;

		delta_Q = (i_tag==mod_mat->A_g->n-1) ? 0 : improve[i_tag];
		division->quality += delta_Q;

	} while(IS_POSITIVE(delta_Q));
	free(unmoved);
	free(indices);
	free(improve);
	free(score);
	return(1);
}

mod_matrix *allocate_partial_modularity_matrix\
(sparse_matrix_arr *adj_matrix, int_vector *vertices_group) {
	mod_matrix *mod_mat;
	if ((mod_mat = malloc(sizeof(mod_matrix))) == NULL) {
		MEMORY_ALLOCATION_FAILURE_AT("allocate_partial_modularity_matrix: mod_mat");
		return NULL;
	}
	if ((mod_mat->A_g = get_partial_sparse_matrix(adj_matrix, vertices_group)) == NULL) {
		free(mod_mat);
		return NULL;
	}
	if ((mod_mat->K = calculate_degree_of_vertices(adj_matrix, vertices_group)) == NULL) {
		free_sparse_matrix_arr(mod_mat->A_g);
		free(mod_mat);
		return NULL;
	}
	mod_mat->total_degree = adj_matrix->rowptr[adj_matrix->n];
	if ((mod_mat->f_g = calculate_F_g_array(mod_mat)) == NULL) {
		free(mod_mat->K);
		free_sparse_matrix_arr(mod_mat->A_g);
		free(mod_mat);
		return NULL;
	}
	mod_mat->norm_1 = calculate_matrix_first_norm(mod_mat);
	return mod_mat;
}



/* Runs Algorithm 2 and returns a new devision, given the Adjacency matrix
	If use_imporve is on, runs algorithm 4*/
two_division *algorithm2\
(sparse_matrix_arr *adj_matrix, int_vector *vgroup, double precision, int use_improve){
	mod_matrix *Bijtag;
	eigen_pair *leading_eigen_pair;
	two_division *division;

	if ((Bijtag = allocate_partial_modularity_matrix(adj_matrix, vgroup)) == NULL) {
		return NULL;
	}
	if ((leading_eigen_pair = calculate_leading_eigen_pair(Bijtag, precision)) == NULL) {
		/* Failed calculating leading eigen pair */
		free_mod_matrix(Bijtag);
		return NULL;
	}
	if ((division = divide_network_in_two(Bijtag, leading_eigen_pair, 1)) == NULL) {
		/* Failed calculating partition */
		free_mod_matrix(Bijtag);
		free_eigen_pair(leading_eigen_pair);
		return NULL;
	}

	free_mod_matrix(Bijtag);
	free_eigen_pair(leading_eigen_pair);
	return division;
}

/* Implementing Algorithm 3 */

n_division* algorithm3 (sparse_matrix_arr *adj_matrix, double precision, int use_improve){
	int n_groups = 0; /* counts number of groups in P */
	int i, j;
	int n = adj_matrix->n; /* number of vertices */
	int_vector *is_divisable; /* indicated whether a group is divisable */
	int_vector *vgroup;
	int next_divisable_group;
	int group_size;
	two_division* division;
	int n_div1; /* Marks the size of the +1 group in a two-division */

	int_vector* p_groups;
	n_division* res;
	
	/* Allocate and Init res */
	if ((p_groups = allocate_int_vector(n)) == NULL){
		/* allocation failed */ 
		return NULL;
	}
	if ((res = allocate_n_division(p_groups)) == NULL) {
		/* allocation failed */ 
		return NULL;
	}
	res->quality = 0.0;
	for (i = 0; i < n; i++){
		res->p_groups->values[i] = -1;
	}

	/* Allocate and Init is_divisable */
	if ((is_divisable = allocate_int_vector(n+1)) == NULL){
		/* allocation failed */ 
		free_n_division(res);
		return NULL;
	}

	/* First we want to remove singletons */
	for (i = 0; i < n; ++i){
		if (adj_matrix->rowptr[i] == adj_matrix->rowptr[i+1]){
			res->p_groups->values[i] = n_groups;
			is_divisable->values[n_groups] = 0;
			++n_groups;
		}
	}

	if (n_groups == n) {
		next_divisable_group = -1;
	}
	else {
		/* Now fill in the "last" group which is the group with every node
		which isn't a singleton */
		for (i = 0; i < n; i++) {
			if (res->p_groups->values[i] == -1)
				res->p_groups->values[i] = n_groups;
		}
		++n_groups;

		is_divisable->values[n_groups-1] = 1;
		next_divisable_group = n_groups-1;
	}
	
	while (next_divisable_group >= 0){
		group_size = 0;
		for (i = 0; i < n; i ++){
			if (res->p_groups->values[i] == next_divisable_group){
				group_size++;
			}
		}
		if ((vgroup = allocate_int_vector(group_size)) == NULL){
			/* allocation failed */
			free_n_division(res);
			free_int_vector(is_divisable);
			return NULL;
		}

		/* filling vgroup with its relevant vertices */
		j = 0;
		for (i = 0; j < group_size; ++i){
			if (res->p_groups->values[i] == next_divisable_group){
				vgroup->values[j++] = i; 
			}
		}

		if ((division = algorithm2(adj_matrix, vgroup, precision, use_improve)) == NULL) {
			/* Error in algorithm2! */
			free_n_division(res);
			free_int_vector(is_divisable);
			free_int_vector(vgroup);
			return NULL;
		}
		else {
			/* Now we count how big is the +1 group */
			n_div1 = 0;
			for (i = 0; i < group_size; i++) {
				n_div1 += (division->s_vector->values[i] > 0) ? 1 : 0;
			}
			if ((n_div1 == 0) || (n_div1 == group_size)) {
				/* No division found... */
				is_divisable->values[next_divisable_group] = 0;
			}
			else {
				/* Woohoo! New division! 
				The division is a two-division, i.e. +1 / -1, and so we need
				to screen the new groups using vgroup and the division we got */
				res->quality += division->quality;
				for (i = 0; i < group_size; i++) {
					if (division->s_vector->values[i] > 0) {
						res->p_groups->values[vgroup->values[i]] = n_groups;
					}
				}
				n_groups++;

				is_divisable->values[n_groups-1] = (n_div1 > 1) ? 1 : 0;
				is_divisable->values[next_divisable_group] = 
					((group_size-n_div1) > 1) ? 1 : 0;
			}

			free_two_division(division);
		}

		/* Free the vgroup, there is no use of it anymore */
		free_int_vector(vgroup);
		
		/* Check if there are more divisable groups */
		while ((next_divisable_group < n_groups) && 
			(!is_divisable->values[next_divisable_group])) {
				next_divisable_group++;
		}
		if (next_divisable_group == n_groups) next_divisable_group = -1;
	}
	
	free_int_vector(is_divisable);
	
	return res;
}

void print_clusters(n_division *division) {
	int i;
	int last_marker;
	int *was_printed;
	if((was_printed = calloc(division->p_groups->n, sizeof(int))) == NULL) {
		return;
	}
	for(last_marker=0;last_marker<division->p_groups->n; last_marker++) {
		if (!was_printed[last_marker]) {
			for(i=last_marker; i<division->p_groups->n; i++) {
				if(division->p_groups->values[i] == division->p_groups->values[last_marker]) {
					/* this value belongs to the group we now print... */
					printf("%d ", i);
					was_printed[i] = 1;
				}
			}
			printf("\n");
		}
	}
	free(was_printed);
}
