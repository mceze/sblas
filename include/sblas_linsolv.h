//
//  sblas_linsolv.h
//  sblas
//
//  Created by Marco Ceze on 2/22/15.
//
//

#ifndef __sblas__sblas_linsolv__
#define __sblas__sblas_linsolv__

#include "sblas.h"

/* function: sblas_bjac */
/* buils a block-Jacobi approximation for Ainv.
 Blocks are at most mxm */
int sblas_bjac(sblas_smat *A, sblas_smat **pM, int m);

/* function: sblas_cg */
/* sparse matrix-vector product: Axb=c*/
int sblas_cg(sblas_smat *A, sblas_svec *b,
                   sblas_svec *x, float const tol, int const niter);

/* function: sblas_qmr_la */
/* Solves A*x = b, using the Quasi-Minimal Residual method*/
int sblas_qmr(sblas_smat *A, sblas_svec *b,
              sblas_svec *x, float const tol,
              int const niter, sblas_smat **pM2);

#endif /* defined(__sblas__sblas_linsolv__) */
