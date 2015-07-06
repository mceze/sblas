//
//  sblas_linsolv.c
//  sblas
//
//  Created by Marco Ceze on 2/22/15.
//
//

#include "sblas.h"

/* function: sblas_cg */
/* Solves A*x = b, with A symmetric, positive definite, 
 using the conjugate gradient method*/
int sblas_cg(sblas_smat *A, sblas_svec *b,
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
    ierr = sblas_error(sblas_smatentry(M, i, i, 1.0));
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
  double theta_prev, gamma_prev, eta_prev, rho_next, qsi_next;
  double rnorm0,rnorm;
  sblas_svec *r, *v, *w, *p=NULL, *q=NULL, *ptil, *s, *d;
  
  /********************/
  // r = b - A*x
  /********************/
  //allocate residual and compute its initial value
  call(sblas_smxv(-1.0, A, False, x, &r, True));
  //add "b" in place
  call(sblas_svpv(1.0, r, 1.0, b, &r));
  rnorm0 = sblas_sv2norm(r);
  printf("|R0| = %1.5e\n",rnorm0);
  
  /********************/
  //vtil = r
  /********************/
  call(sblas_cpvec(r, &v));
  rho = sblas_sv2norm(v);
  
  /********************/
  //wtil = r
  /********************/
  call(sblas_cpvec(r, &w));
  qsi = sblas_sv2norm(w);
  
  gamma_prev = 1.0; eta_prev = -1.0;
  
  for (it = 0; it < niter; it++) {
    //check for breakdown on Lanczos
    if (rho < MEPS || qsi < MEPS)//(uncurable breakdown)
      return sblas_error(sb_BREAKDOWN);
    /********************/
    //v = v/|v|
    /********************/
    call(sblas_scalevec(v, 1.0/rho));
    
    /********************/
    //w = w/|w|
    /********************/
    call(sblas_scalevec(w, 1.0/qsi));

    /********************/
    //delta = w^t*v
    /********************/
    call(sblas_svdv(1.0, w, v, &delta));
    
    //check for breakdown on Lanczos
    if (fabs(delta) < MEPS)//curable
      return sblas_error(sb_BREAKDOWN);
    //compute p and q
    if (it == 0) {
      //this creates both p and q
      call(sblas_cpvec(v, &p));
      call(sblas_cpvec(w, &q));
    }
    else {
      //p = v-((qsi*delta)/eps)*p
      call(sblas_svpv(-(qsi*delta)/eps,p, 1.0, v, &p));
      
      //q = w-((rho*delta)/eps)*q
      call(sblas_svpv(-(rho*delta)/eps,q, 1.0, w, &q));
      
      call(sblas_zerovec(ptil));
    }
    //ptil = A*p
    call(sblas_smxv(1.0, A, False, p, &ptil, (it==0)?True:False));
    //eps = q^T*ptil
    call(sblas_svdv(1.0, q, ptil, &eps));
    //check for breakdown
    if (fabs(eps) < MEPS)
      return sblas_error(sb_BREAKDOWN);
    beta = eps/delta;
    //check for breakdown
    if (fabs(beta) < MEPS)
      return sblas_error(sb_BREAKDOWN);
    //v = ptil-beta*v
    call(sblas_svpv(-beta,v, 1.0, ptil, &v));
    
    //rho = |vtil|
    rho_next = sblas_sv2norm(v);

    //w = A^T*q-beta*w
    call(sblas_scalevec(w, -beta));
    call(sblas_smxv(1.0, A, True, q, &w, False));
    //qsi = |w|
    qsi_next = sblas_sv2norm(w);

    //theta = rho/(gamma*|beta|)
    theta = rho_next/(gamma_prev*fabs(beta));
    //gamma = 1/sqrt(1+theta^2)
    gamma = 1.0/sqrt(1.0+theta*theta);
    if (fabs(gamma) < MEPS)
      return sblas_error(sb_BREAKDOWN);
    //eta = -eta*rho*gamma^2/(beta*gamma^2)
    eta = -eta_prev*rho*gamma*gamma/(beta*gamma_prev*gamma_prev);
    
    if (it == 0){
      //d = eta*p
      call(sblas_cpvec(p, &d));
      call(sblas_scalevec(d, eta));
      
      //s = eta*ptil
      call(sblas_cpvec(ptil, &s));
      call(sblas_scalevec(s, eta));
    }
    else {
      //d = eta*p + (theta*gamma)^2*d
      call(sblas_svpv((theta_prev*gamma)*(theta_prev*gamma),d, eta, p,&d));
      //s = eta*ptil + (theta*gamma)^2*s
      call(sblas_svpv((theta_prev*gamma)*(theta_prev*gamma),s, eta, ptil,&s));
    }
    //x = x+d
    call(sblas_svpv(1.0,x, 1.0, d,&x));
    //r = r-s
    call(sblas_svpv(1.0,r, -1.0, s,&r));
    
    theta_prev = theta;
    gamma_prev  = gamma;
    eta_prev  = eta;
    rho       = rho_next;
    qsi       = qsi_next;
    
    rnorm = sblas_sv2norm(r);
    printf("Iteration %d |R|/|R0| = %1.5e\n",it, rnorm/rnorm0);
    if (rnorm < tol) break;
  }
  
  sblas_destroysvec(r);
  sblas_destroysvec(v);
  sblas_destroysvec(w);
  sblas_destroysvec(p);
  sblas_destroysvec(q);
  sblas_destroysvec(ptil);
  sblas_destroysvec(s);
  sblas_destroysvec(d);
  
  return sb_OK;
}