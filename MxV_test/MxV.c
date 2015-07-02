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
  
  ierr = sblas_error(sblas_readsmat("cg/H.txt", &H));
  if (ierr != sb_OK) return ierr;
  
  ierr = sblas_error(sblas_readsvec("cg/b.txt", &B));
  if (ierr != sb_OK) return ierr;
  
  ierr = sblas_error(sblas_readsvec("cg/x.txt", &x));
  if (ierr != sb_OK) return ierr;
  
//  s = clock();
//  //Matrix vector product
//  ierr = sblas_error(sblas_smxv(1.0, A, True, V, &B));
//  if (ierr != sb_OK) return ierr;
//  e = clock();
//  printf("time: %1.3e\n",(e-s)/CLOCKS_PER_SEC);
  
  
//  ierr = sblas_error(sblas_smxm(1.0, A, True, A, False, &H));
//  if (ierr != sb_OK) return ierr;
//  
//  ierr = sblas_error(sblas_cpvec(B, &x));
//  if (ierr != sb_OK) return ierr;
  
  s = clock();
  ierr = sblas_error(sblas_qmr(H, B, x, 1e-13, 1000));
  if (ierr != sb_OK) return ierr;
  e = clock();
  printf("time: %1.3e\n",(e-s)/CLOCKS_PER_SEC);
  
  //write vector out
  ierr = sblas_error(sblas_writesvecascii("x.txt", x));
  if (ierr != sb_OK) return ierr;
  
  //test in-place vector increment
  ierr = sblas_error(sblas_svpv(2.3, B, 1.3, x, &V));
  if (ierr != sb_OK) return ierr;
  
  //test in-place vector increment
  ierr = sblas_error(sblas_svpv(2.3, B, 1.3, x, &B));
  if (ierr != sb_OK) return ierr;
  
  sblas_svec *test;
  ierr = sblas_error(sblas_svpv(1.0, B, -1.0, V, &test));
  if (ierr != sb_OK) return ierr;
  
  //sblas_destroysmat(A);
  sblas_destroysmat(H);
  sblas_destroysvec(V);
  sblas_destroysvec(test);
  sblas_destroysvec(B);
  sblas_destroysvec(x);
  
  
  return sb_OK;
}