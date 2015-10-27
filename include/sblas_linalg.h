/*
 *  sblas_linalg.h
 *  sblas
 *
 *  Created by Marco Ceze on 2/2/11.
 *  Copyright 2011 University of Michigan. All rights reserved.
 *
 */

#ifndef _sblas_linalg_h
#define _sblas_linalg_h 1

#include "sblas_def.h"
#include "sblas_struct.h"
#include "sblas_enum.h"

/* function: sblas_svdv */
/* sparse dot product between 2 vectors*/
int sblas_svdv(double alpha, sblas_svec *Va,
              sblas_svec *Vb, double *value);

/* function: sblas_smxv */
/* sparse matrix-vector product: A*b=c
 if Alloc = True, a new vector c is created */
int sblas_smxv(double alpha, sblas_smat *A,
               enum sblas_bool TrA,sblas_svec *b,
               sblas_svec **pc, enum sblas_bool Alloc);

/* function: sblas_smxm */
/* sparse matrix-matrix product: op(A)xop(B)=C*/
int sblas_smxm(double alpha, sblas_smat *A, 
               enum sblas_bool TrA, sblas_smat *B, 
               enum sblas_bool TrB, sblas_smat **pC);

/* function: sblas_svpv */
/* sparse sum between 2 vectors*/
int sblas_svpv(double a, sblas_svec *Va, double b,
               sblas_svec *Vb, sblas_svec **pVc);

/* function: sblas_smpm */
/* sparse sum between 2 matrices*/
int sblas_smpm(double a, sblas_smat *A, 
               enum sblas_bool TrA,double b,
               sblas_smat *B, enum sblas_bool TrB,
               sblas_smat **pC);

/* function: sblas_cpvec */
/* set b = a*/
int sblas_cpvec(sblas_svec *a, sblas_svec **pb);

/* function: sblas_zerovec */
/* set Va(i) = 0.0*/
int sblas_zerovec(sblas_svec *Va);

/* function: sblas_scalevec */
/* set Va *= a */
int sblas_scalevec(sblas_svec *Va, double a);

/* function: sblas_svadd */
/* sets Vb = b*Vb + a*Va */
int sblas_svadd(double a, sblas_svec *Va, double b,
                sblas_svec *Vb);

/* function: sblas_sv2norm */
/* vector 2-norm */
double sblas_sv2norm(sblas_svec *Va);

/* function: sblas_lu */
/* decomposes full matrix as A = L*U*/
int sblas_lu(int n, double *a, double *l, double *u);

/* function: sblas_luinv */
/* inverts L and U factors*/
int sblas_luinv(int n, double *l, double *u, double *linv, double *uinv);

/* function: sblas_mxm */
/* product of 2 full matrices*/
int sblas_mxm(int m, int n, int l, double *a, double *b, double *c);

/* function: sblas_svechad */
/* Hadamard (entrywise) product of 2 sparse vectors */
int sblas_svechad(sblas_svec *a, sblas_svec *b,
                  enum sblas_bool invflag, sblas_svec **pc);

#endif


