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

void print_sparse_matrix(sparse_matrix_arr *matrix) {
	int i, j, mat_val_index;
	elem A_i_j;
	for(i=0; i<matrix->n; i++) {
		mat_val_index = matrix->rowptr[i];
		for(j=0; j<matrix->n; j++) {
			if (mat_val_index < matrix->rowptr[i+1]) {
				/* we are still in values for column j row i; */
				/* first we need to advance the sparse matrix pointer */
				while(mat_val_index < matrix->rowptr[i+1] &&
					matrix->colind[mat_val_index] < j) {
					mat_val_index++;
				}
				/* if its in matrix values, take it, otherwise it's 0*/
				/* first condition makes sure we don't overflow from the array boundaries */
				if (mat_val_index < matrix->rowptr[matrix->n] &&
						matrix->colind[mat_val_index] == j) {
					A_i_j = matrix->values[mat_val_index];
					mat_val_index++;
				} else {
					A_i_j = 0;
				}
			} else {
				/* there are no more values for this row, meaning value is 0 */
				A_i_j = 0;
			}
			printf("%d ", (int)A_i_j);
		}
		printf("\n");
	}
}

void free_square_matrix(square_matrix *matrix) {
	free(matrix->values);
	free(matrix);
}

int degree_of_vertice(int i, sparse_matrix_arr *matrix) {
	int j, degree=0;
	for(j=0; j<matrix->rowptr[matrix->n]; j++) {
		if (matrix->colind[j] == i)
			degree++;
	}
	return degree;
}

square_matrix *allocate_square_matrix(int n) {
	square_matrix *matrix;
	if ((matrix = malloc(sizeof(square_matrix))) == NULL) {
		MEMORY_ALLOCATION_FAILURE_AT("allocate_square_double_matrix: matrix");
		return NULL;
	}
	if ((matrix->values = calloc(n*n, sizeof(elem))) == NULL) {
		MEMORY_ALLOCATION_FAILURE_AT("allocate_square_double_matrix: matrix->values");
		free(matrix);
		return NULL;
	}
	matrix->n = n;
	return matrix;
}

void print_square_matrix(square_matrix *mat) {
	int i,j;
	printf("%d ", mat->n);
	for(i=0; i< mat->n; i++) {
		printf("\r\n");
		for(j=0; j<mat->n; j++) {
			printf("%f ", mat->values[(i*mat->n)+j]);
		}
	}
}

/*
 * Computes the value at row i column j of the modularity matrix (non generalized one!!!) for a given vertices group
 * parameter K is an array containing the degrees of vertices in the group
 */
elem modularity_matrix_cell(mod_matrix *mod_mat, int i, int j) {
	elem A_i_j, B_i_j;
	int mat_val_index; /* An index to scan the sparse matrix values*/

	/* find the corresponding A_i_j element in the sparse matrix */
	/*jump to the beginning of the row: */
	mat_val_index = mod_mat->A_g->rowptr[i];
	/*search for the item with colind>=j in that row*/
	while(mat_val_index < mod_mat->A_g->rowptr[i+1] &&
			mod_mat->A_g->colind[mat_val_index] < j) {
		mat_val_index++;
	}
	if (mod_mat->A_g->colind[mat_val_index] == j) {
		/* we found the cell value */
		A_i_j = mod_mat->A_g->values[mat_val_index];
	} else {
		/*cell is not in sparse matrix, hence 0*/
		A_i_j = 0;
	}
	/* next lines calculates: B_i_j = A_i_j - (K_i*K_j)/M */
	B_i_j = A_i_j - (((elem)(mod_mat->K->values[i]*mod_mat->K->values[j]))/mod_mat->total_degree);
	return B_i_j;
}

/*
 * Computes the value at row i column j of the generalized modularity matrix for a given vertices group
 * parameter K is an array containing the degrees of vertices in the group
 * parameter F_g should be calculated from calculate_F_g_array
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
		if (i<(mod_mat->A_g->n-1))
			printf("\r\n");
	}
}

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

/*
 * Calculates an array storing the degree of all vertices in 'vertices_group'
 * from the adjacency matrix  'adj_matrix'
 */
elem_vector *calculate_degree_of_vertices(sparse_matrix_arr *adj_matrix, int_vector *vertices_group) {
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
			/* next line calculates: F_g[i] = Sum(Over all j in vgroup)[B[g]_i_j]*/
			result->values[i] += A_i_j - (((elem)(mod_mat->K->values[i]*mod_mat->K->values[j]))/mod_mat->total_degree);
		}
	}
	return result;
}


elem_vector *mat_vec_multiply(square_matrix *mat, elem_vector *vec) {
	int i, j;
	elem_vector *result;
	if ((result = allocate_elem_vector(vec->n)) == NULL) {
		/* Error allocating  */
		/*TODO: handle error*/
		return NULL;
	}
	for(i=0; i<mat->n; i++) {
		result->values[i] = 0;
		for(j=0; j<mat->n; j++) {
			result->values[i] = result->values[i] + (mat->values[(i*mat->n) + j]*vec->values[j]);
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

elem_vector *vec_substraction(elem_vector *vec1, elem_vector *vec2) {
	elem_vector *result;
	int i;
	if ((result = allocate_elem_vector(vec1->n)) == NULL) {
		/* Error allocating  */
		/*TODO: handle error*/
		return NULL;
	}
	for(i=0; i<vec1->n; i++) {
		result->values[i] = vec1->values[i] - vec2->values[i];
	}
	return result;
}

/*
 * this will calculate vec1-vec2 and store it in vec1
 */
void vec_substraction_inplace(elem_vector *vec1, elem_vector *vec2) {
	int i;
	for(i=0; i<vec1->n; i++) {
		vec1->values[i] -= vec2->values[i];
	}
}

/*
 * this will calculate vec1+vec2 and store it in vec1
 */
void vec_addition_inplace(elem_vector *vec1, elem_vector *vec2) {
	int i;
	for(i=0; i<vec1->n; i++) {
		vec1->values[i] += vec2->values[i];
	}
}

elem vec_norm(elem_vector *vec) {
	return (elem)sqrt(vec_vec_multiply(vec, vec));
}

elem_vector *elem_vec_multiply(elem scalar, elem_vector *vec) {
	elem_vector *result;
	int i;
	if ((result = allocate_elem_vector(vec->n)) == NULL) {
		/* Error allocating  */
		/*TODO: handle error*/
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
	int j, val_index;
	elem_vector *result;
	if ((result = allocate_elem_vector(v->n)) == NULL) {
		return NULL;
	}
	for(val_index=0, j=0; val_index<A->rowptr[A->n]; val_index++) {
		result->values[j] = result->values[j] + (A->values[val_index]*v->values[A->colind[val_index]]);
		while (val_index+1>=A->rowptr[j+1])
			j++;
	}
	return result;
}

/*
 * calculates the product: A[g]*vector
 * where A[g] is the sub matrix of A corresponding to
 * the rows and columns in 'vgroup' values
 */
elem_vector *sparse_mat_vec_multiply_by_group(sparse_matrix_arr *A, int_vector *vgroup, elem_vector *vector) {
	int i, j, val_index;
	elem_vector *result;
	int_vector *reverse_vgroup; /*this holds a reverese mapping to vgroup->values[i], negative in case non existent*/
	if ((result = allocate_elem_vector(vgroup->n)) == NULL) {
		return NULL;
	}
	if ((reverse_vgroup = allocate_int_vector(A->n)) == NULL) {
		return NULL;
	}
	/* initialize reverse_vgroup */
	for(i =0, j=0; i < A->n; i++) {
		if(vgroup->values[j] == i) {
			reverse_vgroup->values[i] = j;
			j++;
		} else
			reverse_vgroup->values[i] = -1;
	}
	/*scan the sparse matrix*/
	for(i=0; i < vgroup->n; i++) {
		result->values[i] = 0;
		val_index = A->rowptr[vgroup->values[i]];
		for(;val_index<A->rowptr[vgroup->values[i]+1]; val_index++) {
			if(reverse_vgroup->values[A->colind[val_index]] >= 0) {
				/*the current column show be included in multiplication */
				result->values[i] = result->values[i] + (A->values[val_index]*vector->values[j]);
			}
		}
	}
	free_int_vector(reverse_vgroup);
	return result;
}

sparse_matrix_arr *elem_array_to_diagonal_matrix(elem *arr, int n) {
	sparse_matrix_arr *matrix;
	int i;
	if ((matrix = allocate_sparse_matrix_arr(n,n)) == NULL) {
		return NULL;
	}
	for(i=0; i<n; i++) {
		matrix->values[i] = arr[i];
		matrix->rowptr[i] = i;
		matrix->colind[i] = i;
	}
	return matrix;
}

elem_vector *int_array_to_elem_vector_by_group(int *arr, int_vector *vertices_group) {
	elem_vector *result;
	int i;
	if ((result = allocate_elem_vector(vertices_group->n)) == NULL) {
		return NULL;
	}
	for(i=0; i<result->n; i++)
		result->values[i] = arr[vertices_group->values[i]];
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
	if ((kgx_M_kg = elem_vec_multiply(vec_vec_multiply(Bijtag->K, vec)/Bijtag->total_degree, Bijtag->K)) == NULL) {
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
		/*TODO: handle error*/
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
	if((result = allocate_eigen_pair(eigen_value, X)) == NULL) {
		/*#TODO: free shit... */
		return NULL;
	}
	return result;
}

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

elem calculate_modularity_score(mod_matrix *mod_mat,elem_vector *s) {
	elem_vector *product;
	elem result;
	if ((product = gen_mod_mat_vec_multiply(mod_mat, s)) == NULL) {
		return -1;
	}
	result = 0.5*vec_vec_multiply(s, product);
	free_elem_vector(product);
	return result;
}

two_division *divide_network_in_two(mod_matrix *mod_mat, eigen_pair *leading_eigen_pair, int use_improve) {
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
		result->quality = calculate_modularity_score(mod_mat, s);
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

int improve_network_division(mod_matrix *mod_mat, two_division *division) {
	int_list_link unmoved = NULL, ptr, ptr_to_j_tag;
	int *indices;
	int i, i_tag;
	elem Q_0, delta_Q;
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
	for(i=0; i<mod_mat->A_g->n; i++) {
		if (unmoved == NULL) {
			/* this is the first value, initialize unmoved pointer*/
			if ((unmoved = malloc(sizeof(int_list_element))) == NULL) {
				MEMORY_ALLOCATION_FAILURE_AT("imrove_network_division: unmoved");
				return 0;
			}
			ptr = unmoved;
		} else {
			/* this is not the first value, add an element to the list */
			if ((ptr->next = malloc(sizeof(int_list_element))) == NULL) {
				MEMORY_ALLOCATION_FAILURE_AT("imrove_network_division: ptr->next");
				free_int_list(unmoved);
				return 0;
			}
			ptr = ptr->next;
		}
		ptr->value = i;
		ptr->next = NULL;
	}
	do {
		for(i=0; i<mod_mat->A_g->n; i++) {
			Q_0 = calculate_modularity_score(mod_mat, division->s_vector);
			for(ptr = unmoved; ptr!=NULL; ptr=ptr->next) {
				division->s_vector->values[ptr->value] *= -1; /* s[k] = -s[k] */
				score[ptr->value] = calculate_modularity_score(mod_mat, division->s_vector) - Q_0;
				division->s_vector->values[ptr->value] *= -1; /* s[k] = -s[k] */
			}

			/* find j (ptr_to_j_tag->value) with maximal score[j] */
			for(ptr_to_j_tag = unmoved, ptr = unmoved; ptr!=NULL; ptr=ptr->next)
				if (score[ptr_to_j_tag->value] < score[ptr->value])
					ptr_to_j_tag = ptr;

			division->s_vector->values[ptr_to_j_tag->value] *= -1; /* S[j_tag] = -S[j_tag] */
			indices[i] = ptr_to_j_tag->value;
			if (i==0)
				improve[i] = score[ptr_to_j_tag->value];
			else
				improve[i] = improve[i-1] + score[ptr_to_j_tag->value];

			/* remove j_tag from unmoved list */
			if(unmoved == ptr_to_j_tag) {
				unmoved = unmoved->next;
				free(ptr_to_j_tag);
			} else {
				for(ptr=unmoved; ptr->next != ptr_to_j_tag; ptr = ptr->next);
				ptr->next = ptr->next->next;
				free(ptr_to_j_tag);
			}
		}
		/* Find the maximum improvement of s and update s accordingly */
		for(i=0, i_tag=0; i<mod_mat->A_g->n; i++)
			if(improve[i_tag]<improve[i])
				i_tag = i;
		for(i=mod_mat->A_g->n-1; i>i_tag; i--)
			division->s_vector->values[indices[i]] *= -1;

		delta_Q = (i_tag==mod_mat->A_g->n-1) ? 0 : improve[i_tag];
	} while(IS_POSITIVE(delta_Q));
	return(1);
}

#ifdef FALSE_DEFINITION

int improve_network_division(square_matrix *mod_mat, two_division *s_vector) {
	int *unmoved, *indices;
	int i,k,j;
	elem Q_0, delta_Q;
	elem *score, *improve;
	elem_vector *s;
	if ((s = allocate_elem_vector(mod_mat->n)) == NULL) {
		return -1;
	}
	if ((score = calloc(mod_mat->n, sizeof(elem))) == NULL) {
		MEMORY_ALLOCATION_FAILURE_AT("improve_network_division: score");
		free_elem_vector(s);
		return -1;
	}
	if ((improve = calloc(mod_mat->n, sizeof(elem))) == NULL) {
		MEMORY_ALLOCATION_FAILURE_AT("improve_network_division: score");
		free_elem_vector(s);
		free(score);
		return -1;
	}
	/* initialize S Vector*/
	for(i=0; i<mod_mat->n; i++) {
		s->values[i] = s_vector->s_vector->values[i];
	}
	if ((unmoved = calloc(mod_mat->n, sizeof(int))) == NULL) {
		MEMORY_ALLOCATION_FAILURE_AT("improve_network_division: unmoved");
		free_elem_vector(s);
		free(score);
		free(improve);
		return 1;
	}
	if ((indices = calloc(mod_mat->n, sizeof(int))) == NULL) {
		MEMORY_ALLOCATION_FAILURE_AT("improve_network_division: indices");
		free_elem_vector(s);
		free(score);
		free(improve);
		free(unmoved);
		return 1;
	}
	do {
		for(i=0; i<mod_mat->n; i++) {
				unmoved[i] = 1; /*this marks all vertices as unmoved */
		}
		for(i=0; i<mod_mat->n; i++) {
			Q_0 = calculate_modularity_score(mod_mat, s);
			for(k=0; k<mod_mat->n; k++) {
				if(unmoved[k] == 1) {
					s->values[k] = -1 * s->values[k];
					score[k] = (0.5 * calculate_modularity_score(mod_mat, s)) - Q_0;
					s->values[k] = -1 * s->values[k];
				}
			}
			/* find vertex j with a maximal score */
			j=0;
			for(k=0; k<mod_mat->n; k++) {
				if (unmoved[k] == 1) {
					if (unmoved[j] == 0) {
						j = k; /*The initial chosen j was illegal, fix it. this can only happen once*/
					} else if(score[j] < score[k]){
						j = k;
					}
				}
			}
			/* change vertex j's group */
			s->values[j] = -1 * s->values[j];
			indices[i] = j;
			if (i==0) {
				improve[i] = score[j];
			} else {
				improve[i] = improve[i-1] + score[j];
			}
			unmoved[j] = 0;
		}
		/* find the maximum improvement of s and update s accordingly */
		i = 0;
		for(k=1; k<mod_mat->n; k++) {
			if (improve[i] < improve[k])
				i = k;
		}
		for(k=mod_mat->n-1; k>i+1; k--) {
			j = indices[k];
			s->values[j] = -1 * s->values[j];
		}
		if(i == mod_mat->n-1) {
			delta_Q = 0;
		} else {
			delta_Q = improve[i];
		}
		printf("%f\n", delta_Q);
	} while (IS_POSITIVE(delta_Q));
}

#endif

mod_matrix *allocate_partial_modularity_matrix(sparse_matrix_arr *adj_matrix, int_vector *vertices_group) {
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

two_division *algorithm3(sparse_matrix_arr *adj_matrix, double precision, int use_improve) {
	two_division *final_div;
	elem_vector *groups;
	int i;
	groups = allocate_elem_vector(adj_matrix->n);
	for(i=0; i<adj_matrix->n; i++) {
		groups->values[i] = i%3;
	}
	final_div = allocate_two_division(groups);
	final_div->quality = 123;
	return final_div;
}

void print_clusters(two_division *division) {
	int i;
	int last_marker;
	int *was_printed;
	if((was_printed = calloc(division->s_vector->n, sizeof(int))) == NULL) {
		return;
	}
	for(last_marker=0;last_marker<division->s_vector->n; last_marker++) {
		if (!was_printed[last_marker]) {
			for(i=last_marker; i<division->s_vector->n; i++) {
				if(division->s_vector->values[i] == division->s_vector->values[last_marker]) {
					/* this value belongs to the group we now print... */
					printf("%d ", i);
					was_printed[i] = 1;
				}
			}
			printf("\n");
		}
	}
}
