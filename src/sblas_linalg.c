/*
 *  sblas_linalg.c
 *  sblas
 *
 *  Created by Marco Ceze on 2/2/11.
 *  Copyright 2011 University of Michigan. All rights reserved.
 *
 */

#include "sblas_linalg.h"
#include "sblas_struct.h"
#include "sblas_def.h"
#include "sblas_enum.h"
#include "sblas_io.h"
#include "sblas_aux.h"
#include "sblas_utils.h"

/* function: sblas_svdv */
/* sparse dot product between 2 vectors*/
int sblas_svdv(double alpha, sblas_svec *Va,
               sblas_svec *Vb, double *value)
{
  int top_i, bot_i, z, begin, end, i, indr;
  int rank;
  sblas_svec *Vl, *Vr;
  
  if (Va == NULL || Vb == NULL)
    return sb_INPUT_ERROR;
  
  if (Va->m != Vb->m)
    return sb_INCOMPATIBLE;
  
  //use the "sparsest" as reference
  if (Va->nZ <= Vb->nZ){
    Vl = Va;
    Vr = Vb;
  }
  else{
    Vl = Vb;
    Vr = Va;
  }
  
  value[0] = 0.0;
  //optimized for identical vectors
  if (Vl == Vr){
    for (z = 0; z < Vl->nZ; z++){
      value[0] += alpha*Vl->val[z]*Vr->val[z];
    }
  }
  else{
    z = 0;
    i = 0;
    begin = 0;
    end = Vr->nZ-1;
    top_i = 0;
    bot_i = Vl->nZ-1;
    while (z <= Vl->nZ){
      //check overlap
      if (Vl->index[top_i] > Vr->index[end]||
          Vl->index[bot_i] < Vr->index[begin])
        return sb_OK;    
      //look at top
      if (begin > end)
        return sb_OK;
      indr = sblas_bsearch(Vl->index[top_i], Vr->index, 
                           begin, end, &rank);
      if (indr >= 0){
        begin = rank+1; 
        value[0] += alpha*Vl->val[top_i]*Vr->val[indr];
      }
      else if (indr != sb_NOT_FOUND)
        return sblas_error(indr);
      z++;
      //look at bottom
      if (bot_i != top_i){
        if (begin > end)
          return sb_OK;
        indr = sblas_bsearch(Vl->index[bot_i], Vr->index, 
                             begin, end, &rank);
        if (indr >= 0){
          end = rank-1; 
          value[0] += alpha*Vl->val[bot_i]*Vr->val[indr];
        }
        else if (indr != sb_NOT_FOUND)
          return sblas_error(indr);
        z++;
      }
      i++;
      top_i = i;
      bot_i = Vl->nZ-1-i;
    }
  }
  
  return sb_OK;
}

/* function: sblas_smxv */
/* sparse matrix-vector product: Axb=c*/
int sblas_smxv(double alpha, sblas_smat *A, 
               enum sblas_bool TrA, 
               sblas_svec *b, sblas_svec **pc)
{
  int ierr, i, m, n;
  double val;
  sblas_svec **Row;
  
  if (A == NULL || b == NULL)
    return sb_INPUT_ERROR;
  
  if (TrA){
    m = A->n;
    n = A->m;
    Row = A->Col;
  }
  else {
    m = A->m;
    n = A->n;
    Row = A->Row;
  }
  
  if (n != b->m)
    return sb_INCOMPATIBLE;
  //create rhs
  ierr = sblas_error(sblas_createsvec(pc, m));
  if (ierr != sb_OK) return ierr;
  
  for (i = 0; i < m; i++){
    ierr = sblas_error(sblas_svdv(alpha, Row[i], b, &val));
    if (ierr != sb_OK) return ierr;
    
    ierr = sblas_error(sblas_svecentry(pc[0], i+1, val));
    if (ierr != sb_OK) return ierr;
  }
  
  return sb_OK;
}

/* function: sblas_smxm */
/* sparse matrix-matrix product: op(A)xop(B)=C*/
int sblas_smxm(double alpha, sblas_smat *A, 
               enum sblas_bool TrA, sblas_smat *B, 
               enum sblas_bool TrB, sblas_smat **pC)
{
  int ierr, mA, nA, mB, nB, i, j;
  double val;
  sblas_svec **Row, **Col;
  
  if (A == NULL || B == NULL)
    return sb_INPUT_ERROR;
  
  if (TrA){
    mA = A->n;
    nA = A->m;
    Row = A->Col;
  }
  else {
    mA = A->m;
    nA = A->n;
    Row = A->Row;
  }
  
  if (TrB){
    mB = B->n;
    nB = B->m;
    Col = B->Row;
  }
  else {
    mB = B->m;
    nB = B->n;
    Col = B->Col;
  }
  
  if (nA != mB)
    return sb_INCOMPATIBLE;
  
  ierr = sblas_error(sblas_createsmat(pC, mA, nB));
  if (ierr != sb_OK) return ierr;
  
  for (i = 0; i < mA; i++){
    for (j = 0; j < nB; j++){
      ierr = sblas_error(sblas_svdv(alpha, Row[i], 
                                    Col[j], &val));
      if (ierr != sb_OK) return ierr;
      
      ierr = sblas_error(sblas_smatentry(pC[0], i+1, 
                                         j+1, val));
      if (ierr != sb_OK) return ierr;
    }
  }
  
  return sb_OK;
}

/* function: sblas_svpv */
/* sparse sum between 2 vectors*/
int sblas_svpv(double a, sblas_svec *Va, double b,
               sblas_svec *Vb, sblas_svec **pVc)
{
  int ierr, z;
  double al, ar;
  sblas_svec *Vl, *Vr;
  
  if (Va == NULL || Vb == NULL)
    return sb_INPUT_ERROR;
  
  if (Va->m != Vb->m)
    return sb_INCOMPATIBLE;
  
  if (Va->nZ <= Vb->nZ){
    Vl = Va;
    al = a;
    Vr = Vb;
    ar = b;    
  }
  else{
    Vl = Vb;
    al = b;
    Vr = Va;
    ar = a;
  }
  
  ierr = sblas_error(sblas_createsvec(pVc, Va->m));
  if (ierr != sb_OK) return ierr;
  
  for (z = 0; z < Vl->nZ; z++){
    ierr = sblas_error(sblas_svecentry(pVc[0], 
                                       Vl->index[z], 
                                       al*Vl->val[z]));
    if (ierr != sb_OK) return ierr;
    
    ierr = sblas_error(sblas_svecentry(pVc[0], 
                                       Vr->index[z], 
                                       ar*Vr->val[z]));
    if (ierr != sb_OK) return ierr;
  }
  
  for (z = Vl->nZ; z < Vr->nZ; z++){
    ierr = sblas_error(sblas_svecentry(pVc[0], 
                                       Vr->index[z], 
                                       ar*Vr->val[z]));
    if (ierr != sb_OK) return ierr;
  }
  
  return sb_OK;
}

/* function: sblas_smpm */
/* sparse sum between 2 matrices*/
int sblas_smpm(double a, sblas_smat *A, 
               enum sblas_bool TrA,double b,
               sblas_smat *B, enum sblas_bool TrB,
               sblas_smat **pC)
{
  int ierr, i;
  sblas_smat *Ma, *Mb, *C;
  
  if (A == NULL || B == NULL)
    return sb_INPUT_ERROR;
  
  if ((Ma = malloc(sizeof(sblas_smat))) == NULL)
    return sb_MEMORY_ERROR;
  if ((Mb = malloc(sizeof(sblas_smat))) == NULL)
    return sb_MEMORY_ERROR;
  
  if (TrA){
    Ma->m = A->n;
    Ma->n = A->m;
    Ma->nZ = A->nZ;
    Ma->Row = A->Col;
    Ma->Col = A->Row;
  }
  else {
    Ma->m = A->m;
    Ma->n = A->n;
    Ma->nZ = A->nZ;
    Ma->Row = A->Row;
    Ma->Col = A->Col;
  }
  if (TrB){
    Mb->m = B->n;
    Mb->n = B->m;
    Mb->nZ = B->nZ;
    Mb->Row = B->Col;
    Mb->Col = B->Row;
  }
  else {
    Mb->m = B->m;
    Mb->n = B->n;
    Mb->nZ = B->nZ;
    Mb->Row = B->Row;
    Mb->Col = B->Col;
  }
  
  if (Ma->m != Mb->m || 
      Ma->n != Mb->n)
    return sb_INCOMPATIBLE;
  
  if ((C = malloc(sizeof(sblas_smat))) == NULL)
    return sb_MEMORY_ERROR;
  C->m = Ma->m;
  C->n = Ma->n;
  C->nZ = 0;
  if ((C->Row = malloc(C->m*sizeof(sblas_svec))) == NULL)
    return sb_MEMORY_ERROR;
  if ((C->Col = malloc(C->n*sizeof(sblas_svec))) == NULL)
    return sb_MEMORY_ERROR;
  
  //add rows
  for (i = 0; i < Ma->m; i++){
    ierr = sblas_error(sblas_svpv(a, Ma->Row[i], 
                                  b, Mb->Row[i], 
                                  &C->Row[i]));
    if (ierr != sb_OK) return ierr;
    C->nZ += C->Row[i]->nZ;
  }
  //add cols
  for (i = 0; i < Ma->n; i++){
    ierr = sblas_error(sblas_svpv(a, Ma->Col[i], 
                                  b, Mb->Col[i], 
                                  &C->Col[i]));
    if (ierr != sb_OK) return ierr;
  }
  
  pC[0] = C;
  
  //clean up
  Ma->Row = NULL;
  Ma->Col = NULL;
  Mb->Row = NULL;
  Mb->Col = NULL;
  free(Ma);
  free(Mb);
  
  return sb_OK;
}
