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

/* function: sblas_conjgrad */
/* sparse matrix-vector product: Axb=c*/
int sblas_conjgrad(sblas_smat *A, sblas_svec *b,
                   sblas_svec *x, float const tol, int const niter);

#endif /* defined(__sblas__sblas_linsolv__) */
