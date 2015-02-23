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
  ierr = sblas_error(sblas_cpvec(x, &t));
  if (ierr != sb_OK) return ierr;
  ierr = sblas_error(sblas_zerovec(t));
  if (ierr != sb_OK) return ierr;
  ierr = sblas_error(sblas_cpvec(x, &r));
  if (ierr != sb_OK) return ierr;
  ierr = sblas_error(sblas_zerovec(r));
  if (ierr != sb_OK) return ierr;
  ierr = sblas_error(sblas_cpvec(x, &Ap));
  if (ierr != sb_OK) return ierr;
  ierr = sblas_error(sblas_zerovec(Ap));
  if (ierr != sb_OK) return ierr;
  
  ierr = sblas_error(sblas_smxv(1.0, A, False, x, &t));
  if (ierr != sb_OK) return ierr;
  //r = b-A*x
  ierr = sblas_error(sblas_svpv(1.0, b, -1.0, t, &r));
  if (ierr != sb_OK) return ierr;
  
  ierr = sblas_error(sblas_cpvec(r, &p));
  if (ierr != sb_OK) return ierr;
  
  ierr = sblas_error(sblas_svdv(1.0, r, r, &num));
  if (ierr != sb_OK) return ierr;
  
  r0 = sqrt(num);
  printf("|R0| = %1.5e\n",r0);
  
  for (it = 0; it < niter; it++) {
    ierr = sblas_error(sblas_smxv(1.0, A, False, p, &Ap));
    if (ierr != sb_OK) return ierr;
    ierr = sblas_error(sblas_svdv(1.0, p, Ap, &den));
    if (ierr != sb_OK) return ierr;
    //alpha = r^T*r/(p^T*A*p)
    alpha = num/den;
    //x = x+alpha*p
    ierr = sblas_error(sblas_zerovec(t));
    if (ierr != sb_OK) return ierr;
    ierr = sblas_error(sblas_svpv(1.0, x, alpha, p, &t));
    if (ierr != sb_OK) return ierr;
    ierr = sblas_error(sblas_zerovec(x));
    if (ierr != sb_OK) return ierr;
    ierr = sblas_error(sblas_cpvec(t, &x));
    if (ierr != sb_OK) return ierr;
    //r = r - alpha*A*p
    sblas_destroysvec(t);
    ierr = sblas_error(sblas_cpvec(r, &t));
    if (ierr != sb_OK) return ierr;
    ierr = sblas_error(sblas_zerovec(r));
    if (ierr != sb_OK) return ierr;
    ierr = sblas_error(sblas_svpv(1.0, t, -alpha, Ap, &r));
    if (ierr != sb_OK) return ierr;
    den = num;
    ierr = sblas_error(sblas_svdv(1.0, r, r, &num));
    if (ierr != sb_OK) return ierr;
    //check convergence
    if (sqrt(num) < tol) break;
    beta = num/den;
    printf("Iteration %d |R|/|R0| = %1.5e\n",it, sqrt(num)/r0);
    //p = r + A*p
    ierr = sblas_error(sblas_zerovec(t));
    if (ierr != sb_OK) return ierr;
    ierr = sblas_error(sblas_svpv(1.0, r, beta, p, &t));
    if (ierr != sb_OK) return ierr;
    ierr = sblas_error(sblas_zerovec(p));
    if (ierr != sb_OK) return ierr;
    ierr = sblas_error(sblas_cpvec(t, &p));
    if (ierr != sb_OK) return ierr;
    
    ierr = sblas_error(sblas_zerovec(Ap));
    if (ierr != sb_OK) return ierr;
  }
  
  ierr = sblas_error(sblas_writesvecascii("x.txt", x));
  if (ierr != sb_OK) return ierr;
  
  sblas_destroysvec(p);
  sblas_destroysvec(t);
  sblas_destroysvec(r);
  sblas_destroysvec(Ap);
  
  
  return sb_OK;
}