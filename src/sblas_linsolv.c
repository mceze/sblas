//
//  sblas_linsolv.c
//  sblas
//
//  Created by Marco Ceze on 2/22/15.
//
//

#include "sblas.h"

/* function: sblas_conjgrad */
/* sparse matrix-vector product: Axb=c*/
int sblas_conjgrad(sblas_smat *A, sblas_svec *b,
                   sblas_svec *x, float const tol, int const niter)
{
  int ierr, it;
  double alpha, beta, den, num, r0;
  sblas_svec *p, *r, *t, *Ap;
  
  //create search direction and residual vectors
  ierr = sblas_error(sblas_createsvec(&p,b->m));
  if (ierr != OK) return ierr;
  
  ierr = sblas_error(sblas_smxv(1.0, A, False, x, &t));
  if (ierr != OK) return ierr;
  //r = b-A*x
  ierr = sblas_error(sblas_svpv(1.0, b, -1.0, t, &r));
  if (ierr != OK) return ierr;
  
  ierr = sblas_error(sblas_cpvec(r, &p));
  if (ierr != OK) return ierr;
  
  ierr = sblas_error(sblas_svdv(1.0, r, r, &num));
  if (ierr != OK) return ierr;
  
  r0 = num;
  printf("|R0| = %1.5e\n",r0);
  
  for (it = 0; it < niter; it++) {
    ierr = sblas_error(sblas_smxv(1.0, A, False, p, &Ap));
    if (ierr != OK) return ierr;
    ierr = sblas_error(sblas_svdv(1.0, p, Ap, &den));
    if (ierr != OK) return ierr;
    //alpha = r^T*r/(p^T*A*p)
    alpha = num/den;
    //x = x+alpha*p
    sblas_destroysvec(t);
    ierr = sblas_error(sblas_svpv(1.0, x, alpha, p, &t));
    if (ierr != OK) return ierr;
    sblas_destroysvec(x);
    ierr = sblas_error(sblas_cpvec(t, &x));
    if (ierr != OK) return ierr;
    //r = r - alpha*A*p
    sblas_destroysvec(t);
    ierr = sblas_error(sblas_cpvec(r, &t));
    if (ierr != OK) return ierr;
    sblas_destroysvec(r);
    ierr = sblas_error(sblas_svpv(1.0, t, -alpha, Ap, &r));
    if (ierr != OK) return ierr;
    den = num;
    ierr = sblas_error(sblas_svdv(1.0, r, r, &num));
    if (ierr != OK) return ierr;
    //check convergence
    if (sqrt(num) < tol) break;
    beta = num/den;
    printf("Iteration %d beta = %1.5e\n",it, beta);
    //p = r + A*p
    sblas_destroysvec(t);
    ierr = sblas_error(sblas_svpv(1.0, r, beta, p, &t));
    if (ierr != OK) return ierr;
    sblas_destroysvec(p);
    ierr = sblas_error(sblas_cpvec(t, &p));
    if (ierr != OK) return ierr;
    
    sblas_destroysvec(Ap);
  }
  
  sblas_destroysvec(p);
  sblas_destroysvec(t);
  sblas_destroysvec(r);
  sblas_destroysvec(Ap);
  
  
  return OK;
}