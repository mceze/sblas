//
//  sblas_linsolv.c
//  sblas
//
//  Created by Marco Ceze on 2/22/15.
//
//

#include "sblas.h"

/* function: sblas_conjgrad */
/* Solves A*x = b, with A symmetric, positive definite, 
 using the conjugate gradient method*/
int sblas_conjgrad(sblas_smat *A, sblas_svec *b,
                   sblas_svec *x, float const tol,
                   int const niter)
{
  int ierr, it, i;
  double alpha, beta, den, num, r0, val;
  sblas_svec *p, *r, *t, *Ap, *bm;
  sblas_smat *M, *Am;
  
  //create Jacobi preconditioner
  ierr = sblas_error(sblas_createsmat(&M, A->m, A->n));
  if (ierr != sb_OK) return ierr;
  for (i = 0; i < A->m; i++) {
    ierr = sblas_error(sblas_smat_getentry(A, i, i, &val));
    if (ierr != sb_OK) return ierr;
    ierr = sblas_error(sblas_smatentry(M, i, i, 1.0/val));
    if (ierr != sb_OK) return ierr;
  }
  
  //naive preconditioner application
  ierr = sblas_error(sblas_smxm(1.0, M, False, A, False, &Am));
  if (ierr != sb_OK) return ierr;
  
  ierr = sblas_error(sblas_smxv(1.0, M, False, b, &bm, True));
  if (ierr != sb_OK) return ierr;
  
  //create search direction and residual vectors
  ierr = sblas_error(sblas_cpvec(x, &t));
  if (ierr != sb_OK) return ierr;
  
  ierr = sblas_error(sblas_zerovec(t));
  if (ierr != sb_OK) return ierr;
  
  ierr = sblas_error(sblas_cpvec(x, &Ap));
  if (ierr != sb_OK) return ierr;
  
  ierr = sblas_error(sblas_zerovec(Ap));
  if (ierr != sb_OK) return ierr;

  //r = b-A*x
  ierr = sblas_error(sblas_zerovec(t));
  if (ierr != sb_OK) return ierr;
  
  ierr = sblas_error(sblas_smxv(1.0, Am, False, x, &t, False));
  if (ierr != sb_OK) return ierr;
  
  ierr = sblas_error(sblas_svpv(1.0, bm, -1.0, t, &r));
  if (ierr != sb_OK) return ierr;
  
  //p = r
  ierr = sblas_error(sblas_cpvec(r, &p));
  if (ierr != sb_OK) return ierr;
  
  ierr = sblas_error(sblas_svdv(1.0, r, r, &num));
  if (ierr != sb_OK) return ierr;
  
  r0 = sqrt(num);
  printf("|R0| = %1.5e\n",r0);
  
  for (it = 0; it < niter; it++) {
    ierr = sblas_error(sblas_smxv(1.0, Am, False, p, &Ap, False));
    if (ierr != sb_OK) return ierr;
    ierr = sblas_error(sblas_svdv(1.0, p, Ap, &den));
    if (ierr != sb_OK) return ierr;
    
    //alpha = r^T*r/(p^T*A*p)
    alpha = num/den;
    
    //x = x+alpha*p
    ierr = sblas_error(sblas_svadd(alpha, p, 1.0, x));
    if (ierr != sb_OK) return ierr;
    
    //r = r - alpha*A*p
    ierr = sblas_error(sblas_svadd(-alpha, Ap, 1.0, r));
    if (ierr != sb_OK) return ierr;
    
    
    den = num;
    ierr = sblas_error(sblas_svdv(1.0, r, r, &num));
    if (ierr != sb_OK) return ierr;
    //check convergence
    if (sqrt(num) < tol) break;
    beta = num/den;
    
    printf("Iteration %d |R|/|R0| = %1.5e\n",it, sqrt(num)/r0);
    //p = r + beta*p
    ierr = sblas_error(sblas_svadd(1.0, r, beta, p));
    if (ierr != sb_OK) return ierr;
    
    ierr = sblas_error(sblas_zerovec(Ap));
    if (ierr != sb_OK) return ierr;
  }
  
  
  sblas_destroysvec(p);
  sblas_destroysvec(t);
  sblas_destroysvec(r);
  sblas_destroysvec(Ap);
  sblas_destroysvec(bm);
  sblas_destroysmat(Am);
  sblas_destroysmat(M);
  
  return sb_OK;
}

/* function: sblas_qmr_la */
/* Solves A*x = b, using the Quasi-Minimal Residual method*/
int sblas_qmr(sblas_smat *A, sblas_svec *b,
              sblas_svec *x, float const tol,
              int const niter)
{
  int ierr, it;
  double rho, qsi, gamma, eta, delta, eps, beta, theta;
  double temp, rnorm0,rnorm;
  sblas_svec *r, *v, *w, *p=NULL, *q=NULL, *ptil, *s, *d;
  
  //allocate residual and compute its initial value
  ierr = sblas_error(sblas_smxv(-1.0, A, False, x, &r, True));
  if (ierr != sb_OK) return ierr;
  //add "b" in place
  ierr = sblas_error(sblas_svpv(1.0, r, 1.0, b, &r));
  if (ierr != sb_OK) return ierr;
  
  ierr = sblas_error(sblas_svdv(1.0, r, r, &rnorm0));
  if (ierr != sb_OK) return ierr;
  rnorm0 = sqrt(rnorm0);
  printf("|R0| = %1.5e\n",rnorm0);
  //vtil = r
  ierr = sblas_error(sblas_cpvec(r, &v));
  if (ierr != sb_OK) return ierr;
  //rho = |v|
  ierr = sblas_error(sblas_svdv(1.0, v, v, &rho));
  if (ierr != sb_OK) return ierr;
  rho = sqrt(rho);
  
  //wtil = r
  ierr = sblas_error(sblas_cpvec(r, &w));
  if (ierr != sb_OK) return ierr;
  //psi = |w|
  ierr = sblas_error(sblas_svdv(1.0, w, w, &qsi));
  if (ierr != sb_OK) return ierr;
  qsi = sqrt(qsi);

  ierr = sblas_error(sblas_cpvec(v, &ptil));
  if (ierr != sb_OK) return ierr;
  
  gamma = 1.0; eta = -1.0;
  
  for (it = 0; it < niter; it++) {
    //check for breakdown on Lanczos
    if (rho < MEPS || qsi < MEPS)//(uncurable breakdown)
      return sblas_error(sb_BREAKDOWN);
    //v = v/|v|
    ierr = sblas_error(sblas_scalevec(v, 1.0/rho));
    if (ierr != sb_OK) return ierr;
    //w = w/|w|
    ierr = sblas_error(sblas_scalevec(w, 1.0/qsi));
    if (ierr != sb_OK) return ierr;
    //delta = w^t*v
    ierr = sblas_error(sblas_svdv(1.0, w, v, &delta));
    if (ierr != sb_OK) return ierr;
    //check for breakdown on Lanczos
    if (delta < MEPS)//curable
      return sblas_error(sb_BREAKDOWN);
    //compute p and q
    if (it == 0) {
      ierr = sblas_error(sblas_cpvec(v, &p));
      if (ierr != sb_OK) return ierr;
      ierr = sblas_error(sblas_cpvec(w, &q));
      if (ierr != sb_OK) return ierr;
    }
    else {
      //p = v-((qsi*delta)/eps)*p
      ierr = sblas_error(sblas_svpv(-(qsi*delta)/eps,p, 1.0, v, &p));
      if (ierr != sb_OK) return ierr;
      //q = w-((rho*delta)/eps)*q
      ierr = sblas_error(sblas_svpv(-(rho*delta)/eps,q, 1.0, w, &q));
      if (ierr != sb_OK) return ierr;
    }
    //ptil = A*p
    ierr = sblas_error(sblas_zerovec(ptil));
    if (ierr != sb_OK) return ierr;
    ierr = sblas_error(sblas_smxv(1.0, A, False, p, &ptil, False));
    if (ierr != sb_OK) return ierr;
    //eps = q^T*ptil
    ierr = sblas_error(sblas_svdv(1.0, q, ptil, &eps));
    if (ierr != sb_OK) return ierr;
    //check for breakdown
    if (eps < MEPS)
      return sblas_error(sb_BREAKDOWN);
    beta = eps/delta;
    //check for breakdown
    if (beta < MEPS)
      return sblas_error(sb_BREAKDOWN);
    //v = ptil-beta*v
    ierr = sblas_error(sblas_svpv(-beta,v, 1.0, ptil, &v));
    if (ierr != sb_OK) return ierr;
    //rho = |vtil|
    ierr = sblas_error(sblas_svdv(1.0, v, v, &rho));
    if (ierr != sb_OK) return ierr;
    rho = sqrt(rho);
    //w = A^T*q-beta*w
    ierr = sblas_error(sblas_scalevec(w, -beta));
    if (ierr != sb_OK) return ierr;
    ierr = sblas_error(sblas_smxv(1.0, A, True, q, &w, False));
    if (ierr != sb_OK) return ierr;
    //psi = |w|
    ierr = sblas_error(sblas_svdv(1.0, w, w, &qsi));
    if (ierr != sb_OK) return ierr;
    qsi = sqrt(qsi);
    //theta = rho/(gamma*|beta|)
    theta = rho/(gamma*fabs(beta));
    //gamma = 1/sqrt(1+theta^2)
    temp = gamma;
    gamma = 1.0/sqrt(1.0+theta*theta);
    if (isnan(gamma))
      return sblas_error(sb_BREAKDOWN);
    //eta = -eta*rho*gamma^2/(beta*gamma^2)
    eta = -eta*rho*gamma*gamma/(beta*temp*temp);
    if (it == 0){
      //d = eta*p
      ierr = sblas_error(sblas_cpvec(p, &d));
      if (ierr != sb_OK) return ierr;
      ierr = sblas_error(sblas_scalevec(d, eta));
      if (ierr != sb_OK) return ierr;
      //s = eta*ptil
      ierr = sblas_error(sblas_cpvec(ptil, &s));
      if (ierr != sb_OK) return ierr;
      ierr = sblas_error(sblas_scalevec(s, eta));
      if (ierr != sb_OK) return ierr;
    }
    else {
      //d = eta*p + (theta*gamma)^2*d
      ierr = sblas_error(sblas_svpv((theta*gamma)*(theta*gamma),d, eta, p,&d));
      if (ierr != sb_OK) return ierr;
      //s = eta*ptil + (theta*gamma)^2*s
      ierr = sblas_error(sblas_svpv((theta*gamma)*(theta*gamma),s, eta, ptil,&s));
      if (ierr != sb_OK) return ierr;
    }
    //x = x+d
    ierr = sblas_error(sblas_svpv(1.0,x, 1.0, d,&x));
    if (ierr != sb_OK) return ierr;
    //r = r-s
    ierr = sblas_error(sblas_svpv(1.0,r, -1.0, s,&r));
    if (ierr != sb_OK) return ierr;
    ierr = sblas_error(sblas_svdv(1.0, r, r, &rnorm));
    if (ierr != sb_OK) return ierr;
    rnorm0 = sqrt(rnorm);
    printf("Iteration %d |R|/|R0| = %1.5e\n",it, rnorm/rnorm0);
  }
  
  
  
  
  
  return sb_OK;
}