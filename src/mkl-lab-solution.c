/*
 ********************************************************************************
 *   Copyright(C) 2016 Filipe Oliveira, Sergio Caldas. Universidade do Minho
 *   All Rights Reserved.
 *
 ********************************************************************************
 *   Content : Simple MKL Sparse Matrix Multiply C example
 *
 ********************************************************************************/

#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include "mkl_spblas.h"

int main(int argc, char* argv[])
{
  clock_t start, stop;
  start = clock();

// Define sparse-matrix M
int mi[5] = {0, 2, 5, 8, 10};
int mj[10] = {0, 1, 0, 1, 2, 1, 2, 3, 2, 3};
double mv[10] = {2.0f, 1.0f, 1.0f, 2.0f, 1.0f, 1.0f, 2.0f, 1.0f, 1.0f, 2.0f};
sparse_matrix_t M;

// Define sparse-matrix N
int ni[5] = {0, 1, 2, 3, 4};
int nj[4] = {0, 1, 2, 3};
double nv[4] = {3.0f, 2.0f, 1.0f, -1.0f};
sparse_matrix_t N;

// create csr matrix
mkl_sparse_d_create_csr(&M, SPARSE_INDEX_BASE_ZERO, 4, 4, mi, mi+1, mj, mv);
mkl_sparse_d_create_csr(&N, SPARSE_INDEX_BASE_ZERO, 4, 4, ni, ni+1, nj, nv);

// do addition
sparse_matrix_t C;
mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE ,M, N, &C);

// free memory
mkl_sparse_destroy(M);
mkl_sparse_destroy(N);
mkl_sparse_destroy(C);


  stop = clock();

  printf("Dgemm_multiply(). Elapsed time = %g seconds\n",
      ((double)(stop - start)) / CLOCKS_PER_SEC);
  return 0;
}

