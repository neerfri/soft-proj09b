#ifndef __SPARSE_MATRIX_H
#define __SPARSE_MATRIX_H

typedef double elem;

#define MEMORY_ALLOCATION_FAILURE_AT(func) fprintf(stderr, "MEMORY ALLOCATION FAILURE AT '%s'. Aborting\n", func); 

typedef struct
{
  int     n;    /* size                  */
  int*    rowptr; /* pointers to where rows begin in the values array. */
  int*    colind; /* column indices */
  elem*   values;
} sparse_matrix_arr;


sparse_matrix_arr* allocate_sparse_matrix_arr(int n, int numNnz);
void free_sparse_matrix_arr(sparse_matrix_arr* matrix);
void  mult_sparse_arr(const sparse_matrix_arr *A, const elem* v, elem* result);

/*sparse_matrix_arr *allocate_and_read_matrix(FILE *fp);*/
elem *allocate_vector(int n);
/*elem *allocate_and_read_vector(int n);*/

#endif /* __SPARSE_MATRIX_H */
