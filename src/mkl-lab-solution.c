/*
 ********************************************************************************
 *   Copyright(C) 2004-2011 Intel Corporation. All Rights Reserved.
 *   
 *   The source code, information  and  material ("Material") contained herein is
 *   owned  by Intel Corporation or its suppliers or licensors, and title to such
 *   Material remains  with Intel Corporation  or its suppliers or licensors. The
 *   Material  contains proprietary information  of  Intel or  its  suppliers and
 *   licensors. The  Material is protected by worldwide copyright laws and treaty
 *   provisions. No  part  of  the  Material  may  be  used,  copied, reproduced,
 *   modified, published, uploaded, posted, transmitted, distributed or disclosed
 *   in any way  without Intel's  prior  express written  permission. No  license
 *   under  any patent, copyright  or  other intellectual property rights  in the
 *   Material  is  granted  to  or  conferred  upon  you,  either  expressly,  by
 *   implication, inducement,  estoppel or  otherwise.  Any  license  under  such
 *   intellectual  property  rights must  be express  and  approved  by  Intel in
 *   writing.
 *   
 *   *Third Party trademarks are the property of their respective owners.
 *   
 *   Unless otherwise  agreed  by Intel  in writing, you may not remove  or alter
 *   this  notice or  any other notice embedded  in Materials by Intel or Intel's
 *   suppliers or licensors in any way.
 *
 ********************************************************************************
 *   Content : Simple MKL Matrix Multiply C example
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

