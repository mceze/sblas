/*
 *  MxV.c
 *  sblas
 *
 *  Created by Marco Ceze on 1/26/11.
 *  Copyright 2011 University of Michigan. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "sblas.h"

int main()
{
  int ierr, i;
  double s,e, qsi;
  sblas_smat *M, *A, *P1=NULL, *P2=NULL, *AM = NULL, *P =  NULL;
  sblas_svec *b=NULL, *x = NULL, *V = NULL;
  
  //build Block-Jacobi preconditioner with Jacobian only
//  call(sblas_readsmat("eddy_test/R_U.txt", &M));
//  call(sblas_bjac(M, &M2, 40));
//  sblas_destroysmat(M);
  
  //R_U
  call(sblas_readsmat("dRdU_p_lA0.txt", &A));
  call(sblas_readsmat("P.txt", &P));
//  //M
//  call(sblas_readsmat("eddy_test/matrix_bkp/naca_n2/M.txt", &M));
//  
//  //build Block-Jacobi preconditioner for R_U
//  call(sblas_bjac(A, &P1, 40));
//
//  //build Block-Jacobi preconditioner for M
//  call(sblas_bjac(M, &P2, 40));
//  
//  //combine matrices
//  call(sblas_smpm(1.0, A, False, 5000.0, M, False, &AM));
//  //combine preconditioners
//  //call(sblas_smpm(1.0, P1, False, 5000.0, P2, False, &P));
//  call(sblas_bjac(AM, &P, 40));
  
  
//  //no preconditioner
//  call(sblas_createsmat(&M2, A->m, A->n));
//  for (i = 0; i < M2->m; i++) {
//    call(sblas_smatentry(M2, i, i, 1.0));
//  }
  
  call(sblas_readsvec("R.txt", &b));
  
  call(sblas_createsvec(&x, b->m));
  call(sblas_zerovec(x));
  
  s = clock();
  call(sblas_qmr(A, b, x, 1e-15, 200, &P));
  
  e = clock();
  printf("time: %1.3e\n",(e-s)/CLOCKS_PER_SEC);
  
  //write vector out
  call(sblas_writesvecascii("x.txt", x));
  
  sblas_destroysmat(A);
  sblas_destroysmat(M);
  sblas_destroysmat(AM);
  sblas_destroysmat(P1);
  sblas_destroysmat(P2);
  sblas_destroysvec(b);
  sblas_destroysvec(x);
  
  
  return sb_OK;
}