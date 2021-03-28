#include "setoper.h"
#include "cdd.h"

#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

#include <Eigen/Dense>

int main(int argc, char *argv[]) {
  Eigen::MatrixXd A_eigen(4, 3);
  A_eigen << 7, -3, 0,
             7,  0, -3,
             1,  1,  0,
             1,  0,  1;
  // Eigen::MatrixXd A_eigen(1, 4);
  // A_eigen << 0,  0,  0, 1;

  dd_PolyhedraPtr poly;
  dd_MatrixPtr A, G;
  dd_rowrange m;
  dd_colrange d;
  dd_ErrorType err;

  dd_set_global_constants();  /* First, this must be called to use cddlib. */

  m = A_eigen.rows();
  d = A_eigen.cols();
  A = dd_CreateMatrix(m, d);
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < d; ++j) {
      dd_set_d(A->matrix[i][j], A_eigen(i, j));
      printf("[%d][%d]: %f ", i, j, A_eigen(i, j));
    }
    printf("\n");
  }
  /* 7 - 3 x1          >= 0
     7         - 3x2   >= 0
     1 +   x1          >= 0
     1         +  x2   >= 0
  */
  A->representation=dd_Inequality;
  poly=dd_DDMatrix2Poly(A, &err);  /* compute the second (generator) representation */
  if (err!=dd_NoError) {
    dd_WriteErrorMessages(stdout,err);
    dd_free_global_constants();
    return -1;
  }
  printf("\nInput is H-representation:\n");
  G=dd_CopyGenerators(poly);
  dd_WriteMatrix(stdout,A);  printf("\n");
  dd_WriteMatrix(stdout,G);

  int lin_num = set_card(G->linset);
  printf("rowsize: %ld, colsize: %ld, lin num: %ld\n", G->rowsize, G->colsize, set_card(G->linset));
  std::vector<int> lin_elements;
  if (lin_num > 0){
    long elem;
    printf("linset: ");
    for (elem=1;elem<=G->linset[0];elem++) {
      if (set_member(elem,G->linset)) {
        lin_elements.push_back(int(elem));
        printf("%ld ",elem);
      }
    }
    printf("\n");
  }

  int Grow = G->rowsize;
  int Gcol = G->colsize;
  Eigen::MatrixXd G_eigen(Grow + lin_num, Gcol);
  for (int i = 0; i < Grow; ++i)
    for (int j = 0; j < Gcol; ++j)
      G_eigen(i,j) = dd_get_d(G->matrix[i][j]);
  for (int i = 0; i < lin_num; ++i) {
    G_eigen.middleRows<1>(Grow + i) = -G_eigen.middleRows<1>(lin_elements[i]-1);
  }
  std::cout << "G_eigen:\n" << G_eigen << std::endl;
  dd_FreeMatrix(A);
  dd_FreeMatrix(G);

  /* Add inequalities:
     7 +  x1   -3x2   = 0
     7 - 3x1   + x2   >= 0
  */
  // m=2;
  // B=dd_CreateMatrix(m,d);
  // dd_set_d(B->matrix[0][0],7.0); dd_set_d(B->matrix[0][1], 1.0); dd_set_d(B->matrix[0][2],-3.0);
  // dd_set_d(B->matrix[1][0],7.0); dd_set_d(B->matrix[1][1],-3.0); dd_set_d(B->matrix[1][2], 1.0);
  // set_addelem(B->linset,1); /* setting the first to be equality */

/* Above dd_set_d is used instead of dd_set_si.  This might be useful when your input is float. Yet,
   0.33333 won't be converted to 1/3 when -DGMPRATIONAL flag is used.  Better alternative might be
   dd_set_si2 function to assign a rational number, e.g.
   dd_set_si2(B->matrix[0][0],1,3).  Use these three assignment functions according to your need.
   These functions are defined in cddmp.h and cddmp.c.
*/

  // dd_DDInputAppend(&poly,B, &err); /* append the two inequalities and compute the generators */
  // if (err!=dd_NoError) goto _L99;
  // A=dd_CopyInequalities(poly);  /* get the inequalities (=input). */
  // G=dd_CopyGenerators(poly);  /* get the generators (=output). */
  // printf("\nNew H-representation with added inequalities:\n");
  // dd_WriteMatrix(stdout,A);  printf("\n");
  // dd_WriteMatrix(stdout,G);

  // dd_FreeMatrix(A);
  // dd_FreeMatrix(B);
  // dd_FreeMatrix(G);
  dd_FreePolyhedra(poly);

  if (err!=dd_NoError) dd_WriteErrorMessages(stdout,err);
  dd_free_global_constants();  /* At the end, this must be called. */
  return 0;
}