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
 * Calculates an array storing the degree of all vertices in 'vertices_group'
 * from the adjacency matrix  'adj_matrix'
 */
int *calculate_degree_of_vertices(sparse_matrix_arr *adj_matrix, int_vector *vertices_group) {
	int *degrees_array;
	int i;
	if ((degrees_array = calloc(vertices_group->n, sizeof(int))) == NULL) {
		MEMORY_ALLOCATION_FAILURE_AT("calculate_degree_of_vertices: degrees_array");
		return NULL;
	}
	/*Calculate and store degree of vertices*/
	for(i=0; i<vertices_group->n; i++) {
		degrees_array[i] = degree_of_vertice(vertices_group->values[i], adj_matrix);
	}
	return degrees_array;
}

elem *calculate_F_g_array(sparse_matrix_arr *adj_matrix, int_vector *vertices_group, int *K) {
	elem *F_g;
	elem A_i_j;
	int i,j;
	int mat_val_index; /* An index to scan the sparse matrix values*/
	if ((F_g = calloc(vertices_group->n, sizeof(elem))) == NULL) {
		MEMORY_ALLOCATION_FAILURE_AT("calculate_F_g_array: F_g");
		return NULL;
	}
	for(i=0; i<vertices_group->n; i++) {
		F_g[i] = 0;
		mat_val_index = adj_matrix->rowptr[vertices_group->values[i]];
		for(j=0; j<vertices_group->n; j++) {
			if (mat_val_index < adj_matrix->rowptr[vertices_group->values[i]+1]) {
				/* we are still in values for row i; */
				/* first we need to advance the sparse matrix pointer to at least column j*/
				while(mat_val_index < adj_matrix->rowptr[vertices_group->values[i]+1] &&
						adj_matrix->colind[mat_val_index] < vertices_group->values[j]) {
					mat_val_index++;
				}
				/* if its in adj_matrix values, take it, otherwise it's 0*/
				/* first condition makes sure we don't overflow from the array boundaries */
				if (mat_val_index < adj_matrix->rowptr[adj_matrix->n] &&
						adj_matrix->colind[mat_val_index] == vertices_group->values[j]) {
					A_i_j = adj_matrix->values[mat_val_index];
					mat_val_index++;
				} else {
					A_i_j = 0;
				}
			} else {
				/* there are no more values for this row, meaning value is 0 */
				A_i_j = 0;
			}
			/* next line calculates: F_g[i] = Sum(Over all j in vgroup)[B[g]_i_j]*/
			F_g[i] = F_g[i] +  A_i_j - (((elem)(K[i]*K[j]))/adj_matrix->rowptr[adj_matrix->n]);
		}
	}
	return F_g;
}

/*
 * Computes the value at row i column j of the modularity matrix (non generalized one!!!) for a given vertices group
 * parameter K is an array containing the degrees of vertices in the group
 */
elem modularity_matrix_cell(sparse_matrix_arr *adj_matrix, int_vector *vertices_group, int *K, int i, int j) {
	elem A_i_j;
	int mat_val_index; /* An index to scan the sparse matrix values*/

	/* find the corresponding A_i_j element in the sparse matrix */
	mat_val_index = adj_matrix->rowptr[vertices_group->values[i]];
	if (mat_val_index < adj_matrix->rowptr[vertices_group->values[i]+1]) {
		/* we are still in values for row i; */
		/* first we need to advance the sparse matrix pointer to at least column j*/
		/*TODO: this can be improved using binary search*/
		while(mat_val_index < adj_matrix->rowptr[vertices_group->values[i]+1] &&
				adj_matrix->colind[mat_val_index] < vertices_group->values[j]) {
			mat_val_index++;
		}
		/* if its in adj_matrix values, take it, otherwise it's 0*/
		/* first condition makes sure we don't overflow from the array boundaries */
		if (mat_val_index < adj_matrix->rowptr[adj_matrix->n] &&
				adj_matrix->colind[mat_val_index] == vertices_group->values[j]) {
			A_i_j = adj_matrix->values[mat_val_index];
			mat_val_index++;
		} else {
			A_i_j = 0;
		}
	} else {
		/* there are no more values for this row, meaning value is 0 */
		A_i_j = 0;
	}
	/* next few lines calculates: B_i_j = A_i_j - (K_i*K_j)/M */
	return(A_i_j - (((elem)(K[i]*K[j]))/adj_matrix->rowptr[adj_matrix->n]));
}

/*
 * Computes the value at row i column j of the generalized modularity matrix for a given vertices group
 * parameter K is an array containing the degrees of vertices in the group
 * parameter F_g should be calculated from calculate_F_g_array
 */
elem generalized_modularity_matrix_cell(sparse_matrix_arr *adj_matrix, int_vector *vertices_group, int *K, elem *F_g, int i, int j) {
	elem B_i_j = modularity_matrix_cell(adj_matrix, vertices_group, K, i, j);
	if (i == j)
		return B_i_j - F_g[i];
	else
		return B_i_j;
}

void print_modularity_matrix(sparse_matrix_arr *adj_matrix, int_vector *vertices_group) {
	int i,j;
	elem *F_g;
	int *K;
	if ((K = calculate_degree_of_vertices(adj_matrix, vertices_group)) == NULL) {
		return;
	}
	if ((F_g = calculate_F_g_array(adj_matrix, vertices_group, K)) == NULL) {
		free(K);
		return;
	}
	for(i=0; i<vertices_group->n; i++) {
		for(j=0; j<vertices_group->n; j++) {
			printf("%f ", generalized_modularity_matrix_cell(adj_matrix, vertices_group, K, F_g, i, j));
		}
		if (i<vertices_group->n-1)
			printf("\r\n");
	}
	free(F_g);
	free(K);
}

elem calculate_matrix_first_norm(sparse_matrix_arr *adj_matrix, int_vector *vertices_group, int *K, elem *F_g) {
	int i,j;
	elem result=0, accumulator;
	elem val_at_i_j;
	for(j=0; j<vertices_group->n; j++) {
		accumulator = 0;
		for(i=0; i<vertices_group->n; i++) {
			val_at_i_j = generalized_modularity_matrix_cell(adj_matrix, vertices_group, K, F_g, i, j);
			accumulator = accumulator + fabs(val_at_i_j);
		}
		if (accumulator > result)
			result = accumulator;
	}
	return result;
}

square_matrix *calculate_modularity_matrix(sparse_matrix_arr *adj_matrix, int_vector *vgroup) {
	square_matrix *mod_mat;
	int i,j, pos, mat_val_index=0;
	elem A_i_j;
	elem *F_g;
	int *K;
	if ((mod_mat = allocate_square_matrix(vgroup->n)) == NULL) {
		MEMORY_ALLOCATION_FAILURE_AT("calculate_modularity_matrix: mod_mat");
		return NULL;
	}
	if ((F_g = calloc(vgroup->n, sizeof(elem))) == NULL) {
		MEMORY_ALLOCATION_FAILURE_AT("calculate_modularity_matrix: F_g");
		free(mod_mat);
		return NULL;
	}
	if ((K = calculate_degree_of_vertices(adj_matrix, vgroup)) == NULL) {
		free(mod_mat);
		free(F_g);
		return NULL;
	}
	for(i=0; i<vgroup->n; i++) {
		F_g[i] = 0;
		mat_val_index = adj_matrix->rowptr[vgroup->values[i]];
		for(j=0; j<vgroup->n; j++) {
			if (mat_val_index < adj_matrix->rowptr[vgroup->values[i]+1]) {
				/* we are still in values for column j row i; */
				/* first we need to advance the sparse matrix pointer */
				while(mat_val_index < adj_matrix->rowptr[vgroup->values[i]+1] &&
						adj_matrix->colind[mat_val_index] < vgroup->values[j]) {
					mat_val_index++;
				}
				/* if its in adj_matrix values, take it, otherwise it's 0*/
				/* first condition makes sure we don't overflow from the array boundries */
				if (mat_val_index < adj_matrix->rowptr[adj_matrix->n] &&
						adj_matrix->colind[mat_val_index] == vgroup->values[j]) {
					A_i_j = adj_matrix->values[mat_val_index];
					mat_val_index++;
				} else {
					A_i_j = 0;
				}
			} else {
				/* there are no more values for this row, meaning value is 0 */
				A_i_j = 0;
			}
			pos = (i*mod_mat->n) + j;
			/* next few lines calculates: B_i_j = A_i_j - (K_i*K_j)/M */
			mod_mat->values[pos] = A_i_j - (((elem)(K[i]*K[j]))/adj_matrix->rowptr[adj_matrix->n]);
			/* next line calculates: F_g[i] = Sum(Over all j in vgroup)[B[g]_i_j]*/
			F_g[i] = F_g[i] +  mod_mat->values[pos];
		}
	}
	for(i=0; i<vgroup->n; i++) {
		pos = (i*mod_mat->n) + i;
		mod_mat->values[pos] = mod_mat->values[pos] - F_g[i];
	}
	free(F_g);
	free(K);
	return mod_mat;
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

elem_vector *sparse_mat_vec_multiply(const sparse_matrix_arr *A, elem_vector *v) {
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

eigen_pair *calculate_dominant_eigen_pair(sparse_matrix_arr *adj_matrix, int_vector *vertices_group, int *K, elem *F_g, double precision) {
	elem_vector *X, *X_next, *temp, *temp2; /* represent X[0], x[1]... from the algorithm */
	elem_vector *K_g;
	elem first_norm_of_mod_mat; /* ||B^[g]|| */
	sparse_matrix_arr *D_g;
	elem eigen_value, convergence_meter, M;
	eigen_pair *result;
	int i;
	M = adj_matrix->rowptr[adj_matrix->n];
	first_norm_of_mod_mat = calculate_matrix_first_norm(adj_matrix, vertices_group, K, F_g);
	if ((D_g = elem_array_to_diagonal_matrix(F_g, vertices_group->n)) == NULL) {
		return NULL;
	}
	if ((K_g = int_array_to_elem_vector_by_group(K, vertices_group)) == NULL) {
		free_sparse_matrix_arr(D_g);
		return NULL;
	}
	X = NULL;
	if ((X_next = allocate_elem_vector(vertices_group->n)) == NULL) {
		/* Error allocating  */
		/*TODO: handle error*/
		return NULL;
	}
	X_next->values[0] = 1;
	for (i=1; i<X_next->n; i++) {
		X_next->values[i] = 0;
	}
	do {
		/*TODO: add if statements to all sub methods that allocate memory */
		vec_normalize(X_next);
		X = X_next;
		X_next = sparse_mat_vec_multiply_by_group(adj_matrix, vertices_group, X);
		temp = elem_vec_multiply(vec_vec_multiply(K_g, X)/M, K_g); /* ( (K[g]*X)/M )*K[g]  */
		vec_substraction_inplace(X_next, temp);
		free_elem_vector(temp);
		temp = sparse_mat_vec_multiply(D_g, X); /* D[g]*X */
		vec_substraction_inplace(X_next, temp);
		free_elem_vector(temp);
		temp = elem_vec_multiply(first_norm_of_mod_mat, X);
		vec_addition_inplace(X_next, temp);
		free_elem_vector(temp);

		eigen_value = vec_vec_multiply(X_next, X)/vec_vec_multiply(X,X);
		/* Check convergence as described in eq. 19 */
		temp = elem_vec_multiply(-1*eigen_value, X);
		vec_addition_inplace(temp, X_next); /* Notice that this is actually subtraction... */
		convergence_meter = vec_norm(temp)/vec_norm(X);
		free_elem_vector(temp);

	} while(convergence_meter > precision);
	vec_normalize(X_next);
	if ((result = allocate_eigen_pair(eigen_value, X_next)) == NULL) {
		/* Error creating eigen pair  */
		free_elem_vector(X);
		free_elem_vector(X_next);
		return NULL;
	}
	free_elem_vector(X);
	return result;
}

eigen_pair *calculate_leading_eigen_pair(sparse_matrix_arr *adj_matrix, int_vector *vertices_group, int *K, elem *F_g, double precision) {
	eigen_pair *result;
	elem first_norm_of_mod_mat;
	if ((result = calculate_dominant_eigen_pair(adj_matrix, vertices_group, K, F_g, precision)) == NULL) {
		return NULL;
	}
	first_norm_of_mod_mat = calculate_matrix_first_norm(adj_matrix, vertices_group, K, F_g);
	result->value -= first_norm_of_mod_mat;
	return result;
}

#ifdef FALSE_DEFINITION
eigen_pair *calculate_leading_eigen_pair(square_matrix *mod_mat, double precision) {
	elem_vector *X, *X_next; /* represent X[0], x[1]... from the algorithm */
	elem first_norm_of_mod_mat; /* ||B^[g]|| */
	int i;
	eigen_pair *result;
	elem eigen_value, convergence_meter;
	X = NULL;
	if ((X_next = allocate_elem_vector(mod_mat->n)) == NULL) {
		/* Error allocating  */
		/*TODO: handle error*/
		return NULL;
	}
	X_next->values[0] = 1;
	for (i=1; i<X_next->n; i++) {
		X_next->values[i] = 0;
	}
	do {
		if (X!=NULL)
			free_elem_vector(X); /* this is not needed in first iteration... */
		vec_normalize(X_next);
		X = X_next;
		X_next = mat_vec_multiply(mod_mat, X);
		eigen_value = vec_vec_multiply(X_next, X)/vec_vec_multiply(X,X);
		/* Check convergence as described in eq. 19 */
		convergence_meter = vec_norm(vec_substraction(mat_vec_multiply(mod_mat, X_next),
				elem_vec_multiply(eigen_value, X_next)))/vec_norm(X_next);
	} while(convergence_meter > precision);
	vec_normalize(X_next);
	if ((result = allocate_eigen_pair(eigen_value, X_next)) == NULL) {
		/* Error creating eigen pair  */
		free_elem_vector(X);
		free_elem_vector(X_next);
		return NULL;
	}
	free_elem_vector(X);
	return result;
}


eigen_pair *shift_and_calculate_leading_eigen_pair(square_matrix *mod_mat, double precision) {
	eigen_pair *result;
	square_matrix *shifted_mod_mat;
	elem A_1_norm=0, A_1_norm_summer;
	int pos;
	int i, j;
	if ((shifted_mod_mat = allocate_square_matrix(mod_mat->n)) == NULL) {
		MEMORY_ALLOCATION_FAILURE_AT("shift_and_calculate_leading_eigen_pair: shifted_mod_mat");
		return NULL;
	}
	for(j=0; j<mod_mat->n; j++) {
		for(i=0; i<mod_mat->n; i++) {
			pos = (i*mod_mat->n) + j;
			A_1_norm_summer = A_1_norm_summer + fabs(mod_mat->values[pos]);
		}
		if (A_1_norm < A_1_norm_summer) {
			A_1_norm = A_1_norm_summer;
		}
	}
	for(i=0; i<mod_mat->n; i++) {
		for(j=0; j<mod_mat->n; j++) {
			pos = (i*mod_mat->n) + j;
			if(i==j) {
				shifted_mod_mat->values[pos] = mod_mat->values[pos] + A_1_norm;
			} else {
				shifted_mod_mat->values[pos] = mod_mat->values[pos];
			}
		}
	}
	result = calculate_leading_eigen_pair(shifted_mod_mat, precision);
	free_square_matrix(shifted_mod_mat);
	return result;
}

#endif

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

elem calculate_improvement_in_modularity(square_matrix *mod_mat,elem_vector *s) {
	elem_vector *product;
	elem result;
	if ((product = mat_vec_multiply(mod_mat, s)) == NULL) {
		return -1;
	}
	result = vec_vec_multiply(s, product);
	free_elem_vector(product);
	return result;
}

two_division *divide_network_in_two(square_matrix *mod_mat, eigen_pair *leading_eigen_pair) {
	elem_vector *s;
	two_division *result;
	int i;
	printf("leading_eigen_pair->value: %f\n", leading_eigen_pair->value);
	if (IS_POSITIVE(leading_eigen_pair->value)) {
		if ((s = get_s_vector_for(leading_eigen_pair->vector)) == NULL) {
			return NULL;
		}
		if ((result = allocate_two_division(mod_mat->n)) == NULL) {
			free_elem_vector(s);
			return NULL;
		}
		result->quality = calculate_improvement_in_modularity(mod_mat, s);
		if (IS_POSITIVE(result->quality)) {
			/*convert result to int_vector */
			for(i=0; i<s->n; i++) {
				result->division->values[i] = s->values[i];
			}
			free_elem_vector(s);
			return result;
		}
		free_elem_vector(s);
		free_two_division(result);
		return NULL;
	}
	return NULL;
}

int improve_network_division(square_matrix *mod_mat, two_division *division) {
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
		s->values[i] = division->division->values[i];
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
			Q_0 = 0.5 * calculate_improvement_in_modularity(mod_mat, s);
			for(k=0; k<mod_mat->n; k++) {
				if(unmoved[k] == 1) {
					s->values[k] = -1 * s->values[k];
					score[k] = (0.5 * calculate_improvement_in_modularity(mod_mat, s)) - Q_0;
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
