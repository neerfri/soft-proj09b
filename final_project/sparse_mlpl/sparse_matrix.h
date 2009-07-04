#ifndef __SPARSE_MATRIX_H
#define __SPARSE_MATRIX_H

typedef double elem;

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

#endif /* __SPARSE_MATRIX_H */
