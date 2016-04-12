#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

typedef struct {
  int i, j;
  double aij;
} entry;

static int compare_row(const void *A, const void *B) {
  entry *a, *b;
  a = (entry *) A;
  b = (entry *) B;
  if (a->i < b->i) return(-1);
  if (a->i > b->i) return(1);
  if (a->j < b->j) return(-1);
  if (a->j > b->j) return(1);
  /* this should never happen since duplicate entries are not allowed */
  assert(1);
  return(0);
}

static int compare_col(const void *A, const void *B) {
  entry *a, *b;
  a = (entry *) A;
  b = (entry *) B;
  if (a->j < b->j) return(-1);
  if (a->j > b->j) return(1);
  if (a->i < b->i) return(-1);
  if (a->i > b->i) return(1);
  /* this should never happen since duplicate entries are not allowed */
  assert(1);
  return(0);
}

/*converts the given matrix in i/j/Aij format to CSR format. The results are in
 *   bp/bj/bx vectors. Note: Duplicate entries are not allowed.
 *     flag = 0, CSR format
 *       flag = 1, CSC format
 *       */
void convert_format(int ndim, int nnz, int *irow, int *icol, double *a, int* bp, int* bj, double* bx, int flag) {
  int i, j, k, tmp;
  entry *A;
  A =  entry [nnz];
  assert(A);
  /* copy data to sort array */
  for (i=0; i<nnz; i++) {
    A[i].i = irow[i];
    A[i].j = icol[i];
    A[i].aij = a[i];
  }
  /* sort the array */
  if (flag == 0) {
    qsort(A, nnz, sizeof(entry), compare_row);
  } else {
    qsort(A, nnz, sizeof(entry), compare_col);
  }

  printf("After Sort:\n");
  for (i=0; i<nnz; i++) {
    printf("%d\t%d\t%g\n",A[i].i, A[i].j,A[i].aij);
  }


  /* re-do the row array */
  k = bp[0] = j = 0;
  for (i=1; i<nnz; i++) {
    if (flag == 0) {
      tmp = A[i].i;
    } else {
      tmp = A[i].j;
    }
    if (tmp != k) {
      k = tmp;
      bp[++j] = i;
    }
  }
  bp[++j] = nnz;
  for (i=0; i<nnz; i++) {
    if (flag == 0) {
      bj[i] = A[i].j;
    } else {
      bj[i] = A[i].i;
    }
    bx[i] = A[i].aij;
  }
  delete [] A;
}

#if TEST
void main() {
  int i, info, nnz, ndim;
  void    *S;
  int *irow, *icol;
  double* a;

  int ROW[12] = {
    0, 1, 4, 1, 2, 4, 0, 2, 0, 3, 3, 4
  };
  int COL[12] = {
    0, 0, 0, 1, 1, 1, 2, 2, 3, 3, 4, 4
  };
  double AIJ[12] = {
    19, 12, 12, 21, 12, 12, 21, 16, 21, 5, 21, 18
  };
  nnz = 12;
  ndim = 5;
    int ROW[32] = {
      0, 1, 2, 3, 4, 5, 6, 7, 8, 2, 7, 5, 8, 7, 8, 5, 6, 7, 7, 8, 0, 0, 0, 1, 1, 1, 2, 3, 4, 4, 5, 5
    };
  int COL[32] = {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 0, 0, 1, 1, 2, 2, 3, 3, 4, 6, 7, 2, 6, 7, 2, 4, 8, 6, 7, 5, 8, 6, 7 
  };
  double AIJ[32] = {
    14, 14, 16, 14, 14, 16, 16, 71, 16, -1, -3, -2, -1, -4, -2, -1, -1, -3, -4, -4, -5, -1, -6, -1, -3, -1, -3, -3, -1, -1, -2, -4
  };
  nnz = 32;
  ndim = 9;
  irow = new int [ndim+1];
  icol = new int [nnz];
  a = new double [nnz];

  for (i=0; i<nnz; i++) {
    printf("%d\t%d\t%g\n",ROW[i],COL[i],AIJ[i]);
  }
  convert_format(ndim,nnz,ROW,COL,AIJ,irow,icol,a, 1);
  printf("converted format:\n");
  for (i=0; i<=ndim; i++) {
    printf("%d\n",irow[i]+1);
  }
  for (i=0; i<nnz; i++) {
    printf("%d\t%d\t%g\n",i,icol[i]+1,a[i]);
  }
  delete [] irow;
  delete [] icol;
  delete [] a;
  exit(-1);
}
#endif
