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
    end = max(Vr->nZ-1,0);
    top_i = 0;
    bot_i = max(Vl->nZ-1,0);
    while (z < Vl->nZ){
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
      bot_i = max(Vl->nZ-1-i,0);
    }
  }
  
  return sb_OK;
}

/* function: sblas_smxv */
/* sparse matrix-vector product: A*b=c
 if Alloc = True, a new vector c is created */
int sblas_smxv(double alpha, sblas_smat *A,
               enum sblas_bool TrA,sblas_svec *b,
               sblas_svec **pc, enum sblas_bool Alloc)
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
    return sblas_error(sb_INCOMPATIBLE);
  if (Alloc){
    //create rhs
    call(sblas_createsvec(pc, m));
    
  }
  else if (pc[0]->m != b->m)
    return sblas_error(sb_INCOMPATIBLE);
  
  for (i = 0; i < m; i++){
    call(sblas_svdv(alpha, Row[i], b, &val));
    
    
    call(sblas_svecentry(pc[0], i, val));
    
  }
  
  return sb_OK;
}

/* function: sblas_smxm */
/* sparse matrix-matrix product: op(A)*op(B)=C*/
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
  
  call(sblas_createsmat(pC, mA, nB));
  
  
  for (i = 0; i < mA; i++){
    for (j = 0; j < nB; j++){
      call(sblas_svdv(alpha, Row[i], 
                                    Col[j], &val));
      
      
      call(sblas_smatentry(pC[0], i,
                                         j, val));
      
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
  int inplace = 0;
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
  
  //check if we are operating in place
  if (pVc[0] == Vl)
    inplace = -1;
  if (pVc[0] == Vr)
    inplace = 1;
  
  if (inplace == 0){
    call(sblas_createsvec(pVc, Va->m));
    
    //sum small array region
    for (z = 0; z < Vl->nZ; z++){
      call(sblas_svecentry(pVc[0],
                                         Vl->index[z],
                                         al*Vl->val[z]));
      
      
      call(sblas_svecentry(pVc[0],
                                         Vr->index[z],
                                         ar*Vr->val[z]));
      
    }
    //sum remainder part
    for (z = Vl->nZ; z < Vr->nZ; z++){
      call(sblas_svecentry(pVc[0],
                                         Vr->index[z],
                                         ar*Vr->val[z]));
      
    }
  }
  else if (inplace == -1){//left vector is updated
    //multiply first
    for (z = 0; z < Vl->nZ; z++){
      Vl->val[z] *= al;
    }
    //add second vector
    for (z = 0; z < Vr->nZ; z++){
      call(sblas_svecentry(pVc[0],
                                         Vr->index[z],
                                         ar*Vr->val[z]));
      
    }
  }
  else if (inplace == 1){
    //multiply first
    for (z = 0; z < Vr->nZ; z++){
      Vr->val[z] *= ar;
    }
    //add second vector
    for (z = 0; z < Vl->nZ; z++){
      call(sblas_svecentry(pVc[0],
                                         Vl->index[z],
                                         al*Vl->val[z]));
      
    }
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
    call(sblas_svpv(a, Ma->Row[i], 
                                  b, Mb->Row[i], 
                                  &C->Row[i]));
    
    C->nZ += C->Row[i]->nZ;
  }
  //add cols
  for (i = 0; i < Ma->n; i++){
    call(sblas_svpv(a, Ma->Col[i], 
                                  b, Mb->Col[i], 
                                  &C->Col[i]));
    
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

/* function: sblas_cpvec */
/* set b = a*/
int sblas_cpvec(sblas_svec *a, sblas_svec **pb)
{
  int ierr;
  sblas_svec *b;
  
  call(sblas_createsvec(&b, a->m));
  
  
  b->nZ = a->nZ;
  b->nZprealloc = a->nZprealloc;
  
  if ((b->index = malloc(b->nZprealloc*sizeof(int))) == NULL)
    return sb_MEMORY_ERROR;
  memcpy(b->index+0, a->index+0, a->nZprealloc*sizeof(int));
  if ((b->val = malloc(b->nZprealloc*sizeof(double))) == NULL)
    return sb_MEMORY_ERROR;
  memcpy(b->val+0, a->val+0, a->nZprealloc*sizeof(double));
  
  (*pb) = b;
  
  return sb_OK;
}

/* function: sblas_zerovec */
/* set Va(i) = 0.0*/
int sblas_zerovec(sblas_svec *Va)
{
  memset(Va->val, 0.0, Va->nZprealloc*sizeof(double));
  Va->nZ = 0;
  
  return sb_OK;
}

/* function: sblas_scalevec */
/* set Va *= a */
int sblas_scalevec(sblas_svec *Va, double a)
{
  int i;
  if (fabs(a) <= MEPS) return sblas_error(sblas_zerovec(Va));
  for (i = 0;i < Va->nZ;i++)Va->val[i] *= a;
  
  return sb_OK;
}

/* function: sblas_svadd */
/* sets Vb = b*Vb + a*Va */
int sblas_svadd(double a, sblas_svec *Va, double b, sblas_svec *Vb)
{
  int ierr, i, idx;
  double val;
  
  if (Va->m != Vb->m) return sblas_error(sb_INCOMPATIBLE);
  
  if (b == 0.0){
    call(sblas_zerovec(Vb));
  }
  else{
    call(sblas_scalevec(Vb, b));
  }
  if (a != 0.0)
    for (i = 0; i < Va->nZ; i++) {
      idx = Va->index[i];
      val = a*Va->val[i];
      call(sblas_svecentry(Vb, idx, val));
      
    }
  
  return sb_OK;
}

/* function: sblas_sv2norm */
/* vector 2-norm */
double sblas_sv2norm(sblas_svec *Va)
{
  double norm = 0.0;
  int i;
  
  for (i = 0; i < Va->nZ; i++)
    norm += Va->val[i]*Va->val[i];
  norm = sqrt(norm);
  
  return norm;
}

/* function: sblas_lu */
/* decomposes full matrix as A = L*U*/
int sblas_lu(int n, double *a, double *l, double *u)
{
  
  int i, j, m;
  
  for(i = 0; i < n; i++) for(j = 0; j < n; j++) l[i*n+j] = u[i*n+j] = 0.0;
  
  for(i = 0; i < n; i++) l[i*n+i] = 1.0;
  
  for(m = 0; m < n; m++){
    for(i = m; i < n; i++){
      u[m*n+i] = a[m*n+i];
      for(j = 0; j < m; j++){
        u[m*n+i] -= l[m*n+j]*u[j*n+i];
      }
    }
    for(i = m+1; i < n; i++){
      l[i*n+m] = a[i*n+m];
      for (j = 0; j < m; j++){
        l[i*n+m] -= l[i*n+j]*u[j*n+m];
      }
      if (fabs(u[m*n+m]) < MEPS)
        return sblas_error(sb_BREAKDOWN);
      l[i*n+m] = l[i*n+m]/u[m*n+m];
    }
  }
  
  return sb_OK;
}

/* function: sblas_luinv */
/* inverts L and U factors*/
int sblas_luinv(int n, double *l, double *u, double *linv, double *uinv)
{
  int i, j, k;
  double sum;
  //Computing inverse of U
  for(j = n-1; j >= 0; j--){
    for(i = n-1; i >= 0; i--){
      if(i == j) uinv[i*n+j] = 1.0/u[i*n+j];
      else {
        sum = 0.0;
        for(k = n-1; k > i; k--) sum += u[i*n+k]*uinv[k*n+j];
        uinv[i*n+j] = -sum/u[i*n+i];
      }
    }
  }
  
  //Computing inverse of L
  for(j = 0; j < n; j++){
    for(i = 0; i < n; i++){
      if(i == j) linv[i*n+j] = 1.0/l[i*n+j];
      else {
        sum = 0.0;
        for(k = 0; k < i; k++) sum += l[i*n+k]*linv[k*n+j];
        linv[i*n+j] = -sum/l[i*n+i];
      }
    }
  }
  
  return sb_OK;
}

/* function: sblas_mxm */
/* product of 2 full matrices*/
int sblas_mxm(int m, int n, int l, double *a, double *b, double *c)
{
  int i, j, k;
  
  for (i = 0; i < m; i++) {
    for (k = 0; k < l; k++) {
      c[i*m+k] = 0.0;
      for (j = 0; j < n; j++) {
        c[i*m+k] += a[i*m+j]*b[j*n+k];
      }
    }
  }
  
  return sb_OK;
}

/* function: sblas_svechad */
/* Hadamard (entrywise) product of 2 sparse vectors */
int sblas_svechad(sblas_svec *a, sblas_svec *b,
                  enum sblas_bool invflag, sblas_svec **pc)
{
  int ierr, i, ai;
  double v;
  
  if (a->m != b->m) return sblas_error(sb_INCOMPATIBLE);
  
  call(sblas_createsvec(pc, a->m));
  
  for (i = 0; i < a->nZ; i++) {
    ai = a->index[i];
    call(sblas_svec_getentry(b, ai, &v));
    if(invflag)
      v = 1.0/v;
    call(sblas_svecentry((*pc), ai, a->val[i]*v));
  }
  
  return sb_OK;
}

