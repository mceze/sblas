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
  int ierr;
  double s,e;
  sblas_smat *A, *H;
  sblas_svec *V, *B=NULL, *x;
  
  call(sblas_readsmat("eddy_test/R_U_p_M.txt", &H));
  
  
  call(sblas_readsvec("eddy_test/R.txt", &B));
  
  
//  call(sblas_bjac(H, &A, 40));
//  
//  call(sblas_writesmatascii("eddy_test/BJac.txt", A, True));
//  
//  sblas_destroysmat(A);
  
  
  
  
//  call(sblas_readsvec("cg/x.txt", &x));
//  
  
//  s = clock();
//  //Matrix vector product
//  call(sblas_smxv(1.0, A, True, V, &B));
//  
//  e = clock();
//  printf("time: %1.3e\n",(e-s)/CLOCKS_PER_SEC);
  
  
//  call(sblas_smxm(1.0, A, True, A, False, &H));
//  
//  
  call(sblas_cpvec(B, &x));
  
  
  call(sblas_zerovec(x));
  
  
  s = clock();
  call(sblas_qmr(H, B, x, 1e-10, 1000));
  
  e = clock();
  printf("time: %1.3e\n",(e-s)/CLOCKS_PER_SEC);
  
  //write vector out
  call(sblas_writesvecascii("x.txt", x));
  
  
  //test in-place vector increment
  call(sblas_svpv(2.3, B, 1.3, x, &V));
  
  
  //test in-place vector increment
  call(sblas_svpv(2.3, B, 1.3, x, &B));
  
  
  sblas_svec *test;
  call(sblas_svpv(1.0, B, -1.0, V, &test));
  
  
  //sblas_destroysmat(A);
  sblas_destroysmat(H);
  sblas_destroysvec(V);
  sblas_destroysvec(test);
  sblas_destroysvec(B);
  sblas_destroysvec(x);
  
  
  return sb_OK;
}