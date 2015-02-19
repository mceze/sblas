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
  sblas_smat *A;
  sblas_svec *V, *B=NULL;
  
  ierr = sblas_error(sblas_readsmat("A.txt", &A));
  if (ierr != OK) return ierr;
  
  ierr = sblas_error(sblas_readsvec("R.txt", &V));
  if (ierr != OK) return ierr;
  
  s = clock();
  //Matrix vector product
  ierr = sblas_error(sblas_smxv(1.0, A, True, V, &B));
  if (ierr != OK) return ierr;
  e = clock();
  printf("time: %1.3e\n",(e-s)/CLOCKS_PER_SEC);
  
  //write vector out
  ierr = sblas_error(sblas_writesvecascii("b.txt", B));
  if (ierr != OK) return ierr;
  
  
  sblas_destroysmat(A);
  sblas_destroysvec(V);
  sblas_destroysvec(B);
  
  return OK;
}