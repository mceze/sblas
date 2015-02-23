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
/* sparse matrix-vector product: op(A)xb=c*/
int sblas_smxv(double alpha, sblas_smat *A, 
               enum sblas_bool TrA, 
               sblas_svec *b, sblas_svec **pc);

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
/* set a(i) = 0.0*/
int sblas_zerovec(sblas_svec *a);

#endif