//
//  sblas_linsolv.c
//  sblas
//
//  Created by Marco Ceze on 2/22/15.
//
//

#include "sblas.h"

/* function: sblas_bjac */
/* buils a block-Jacobi approximation for Ainv.
 Blocks are at most mxm */
int sblas_bjac(sblas_smat *A, sblas_smat **pM, int m)
{
  int ierr, t, ib, shift, i,j, n;
  int nb;
  double *Ab, *L, *U, *Linv, *Uinv, *Abinv;
  
  if (A->m != A->n) return sblas_error(sb_INPUT_ERROR);
  
  //memory for blocks
  if ((Ab = malloc(m*m*sizeof(double))) == NULL)
    return sblas_error(sb_MEMORY_ERROR);
  if ((L = malloc(m*m*sizeof(double))) == NULL)
    return sblas_error(sb_MEMORY_ERROR);
  if ((U = malloc(m*m*sizeof(double))) == NULL)
    return sblas_error(sb_MEMORY_ERROR);
  if ((Linv = malloc(m*m*sizeof(double))) == NULL)
    return sblas_error(sb_MEMORY_ERROR);
  if ((Uinv = malloc(m*m*sizeof(double))) == NULL)
    return sblas_error(sb_MEMORY_ERROR);
  if ((Abinv = malloc(m*m*sizeof(double))) == NULL)
    return sblas_error(sb_MEMORY_ERROR);
  
  //memory for preconditioner
  call(sblas_createsmat(pM, A->m, A->n));
  
  //number of blocks
  nb = floor(A->m/m);
  //remainder block size
  t = A->m-nb*m;
  
  //loop over mxm blocks and invert them
  for (ib = 0; ib < nb+1; ib++) {
    shift = ib*m;
    if (ib < nb){
      n = m;
    }
    else {
      n = t;
      if (t == 0)
        break;
    }
    //build local block
    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
        call(sblas_smat_getentry(A, shift+i, shift+j, &Ab[i*n+j]));
      }
    }
    //compute Ab = L*U
    call(sblas_lu(n, Ab, L, U));
    //compute  Linv and Uinv
    call(sblas_luinv(n, L, U, Linv, Uinv));
    //compute Ab^(-1) = Uinv^(-1)*Linv^(-1)
    call(sblas_mxm(n, n, n, Uinv, Linv, Abinv));
    //put inverse in sparse matrix
    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
        call(sblas_smatentry((*pM), shift+i, shift+j, Abinv[i*n+j]));
      }
    }
  }
  
  free(Ab);
  free(Abinv);
  free(L);
  free(U);
  free(Linv);
  free(Uinv);
  
  
  return sb_OK;
}

/* function: sblas_cg */
/* Solves A*x = b, with A symmetric, positive definite, 
 using the conjugate gradient method*/
int sblas_cg(sblas_smat *A, sblas_svec *b,
                   sblas_svec *x, float const tol,
                   int const niter)
{
  int ierr, it;
  double alpha, beta, den, num, r0;
  sblas_svec *p, *r, *t, *Ap;
  
  //create search direction and residual vectors
  call(sblas_cpvec(x, &t));
  
  
  call(sblas_zerovec(t));
  
  
  call(sblas_cpvec(x, &Ap));
  
  
  call(sblas_zerovec(Ap));
  

  //r = b-A*x
  call(sblas_zerovec(t));
  
  
  call(sblas_smxv(1.0, A, False, x, &t, False));
  
  
  call(sblas_svpv(1.0, b, -1.0, t, &r));
  
  
  //p = r
  call(sblas_cpvec(r, &p));
  
  
  call(sblas_svdv(1.0, r, r, &num));
  
  
  r0 = sqrt(num);
  printf("|R0| = %1.5e\n",r0);
  
  for (it = 0; it < niter; it++) {
    call(sblas_smxv(1.0, A, False, p, &Ap, False));
    
    call(sblas_svdv(1.0, p, Ap, &den));
    
    
    //alpha = r^T*r/(p^T*A*p)
    alpha = num/den;
    
    //x = x+alpha*p
    call(sblas_svadd(alpha, p, 1.0, x));
    
    
    //r = r - alpha*A*p
    call(sblas_svadd(-alpha, Ap, 1.0, r));
    
    
    
    den = num;
    call(sblas_svdv(1.0, r, r, &num));
    
    //check convergence
    if (sqrt(num) < tol) break;
    beta = num/den;
    
    printf("Iteration %d |R|/|R0| = %1.5e\n",it, sqrt(num)/r0);
    //p = r + beta*p
    call(sblas_svadd(1.0, r, beta, p));
    
    
    call(sblas_zerovec(Ap));
    
  }
  
  
  sblas_destroysvec(p);
  sblas_destroysvec(t);
  sblas_destroysvec(r);
  sblas_destroysvec(Ap);
    
  return sb_OK;
}

/* function: sblas_qmr_la */
/* Solves A*x = b, using the Quasi-Minimal Residual method*/
int sblas_qmr(sblas_smat *A, sblas_svec *b,
              sblas_svec *x, float const tol,
              int const niter, sblas_smat **pM2)
{
  int ierr, it, m = b->m;
  double rho, qsi, gamma, eta, delta, eps, beta, theta;
  double theta_prev, gamma_prev, eta_prev, rho_next, qsi_next;
  double rnorm0,rnorm, rnormtrue;
  sblas_svec *r, *d, *s, *rtest;
  sblas_svec *v, *w, *y, *z, *p, *q;
  sblas_svec *vtil, *wtil, *ytil, *ztil, *ptil;
  sblas_smat *M2;
  
  //allocate memory
  call(sblas_createsvec(&r, m));
  call(sblas_createsvec(&d, m));
  call(sblas_createsvec(&s, m));
  call(sblas_createsvec(&v, m));
  call(sblas_createsvec(&w, m));
  call(sblas_createsvec(&y, m));
  call(sblas_createsvec(&z, m));
  call(sblas_createsvec(&p, m));
  call(sblas_createsvec(&q, m));
  call(sblas_createsvec(&vtil, m));
  call(sblas_createsvec(&wtil, m));
  call(sblas_createsvec(&ytil, m));
  call(sblas_createsvec(&ztil, m));
  call(sblas_createsvec(&ptil, m));
  call(sblas_createsvec(&rtest, m));
  
  //create preconditioner
  if ((*pM2) == NULL){
    call(sblas_bjac(A, &M2, 40));
    (*pM2) = M2;
  }
  else {
    M2 = (*pM2);
  }
  
  //r = -A*x
  //calculate residual and compute its initial value
  call(sblas_zerovec(r));
  call(sblas_smxv(-1.0, A, False, x, &r, False));
  //r += b
  call(sblas_svadd(1.0, b, 1.0, r));
  rnorm0 = sblas_sv2norm(r);
  printf("|R0| = %1.5e\n",rnorm0);
  
  //vtil = r
  call(sblas_svadd(1.0, r, 0.0, vtil));
  //solve M1*y = vtil (for now we set M1 = I)
  call(sblas_svadd(1.0, vtil, 0.0, y));
  rho = sblas_sv2norm(y);
  
  //wtil = r
  call(sblas_svadd(1.0, r, 0.0, wtil));
  //solve M2^T*z = wtil
  call(sblas_zerovec(z));
  call(sblas_smxv(1.0, M2, True, wtil, &z, False));
  qsi = sblas_sv2norm(z);
  
  gamma_prev = 1.0; eta_prev = -1.0;
  
  for (it = 0; it < niter; it++) {
    //check for breakdown on Lanczos
    if (rho < MEPS || qsi < MEPS)//(uncurable breakdown)
      return sblas_error(sb_BREAKDOWN);

    //v = vtil/rho
    call(sblas_svadd(1.0/rho, vtil, 0.0, v));
    //y = y/rho
    call(sblas_scalevec(y, 1.0/rho));
    
    //w = wtil/qsi
    call(sblas_svadd(1.0/qsi, wtil, 0.0, w));
    //z = z/qsi
    call(sblas_scalevec(z, 1.0/qsi));

    //delta = z^t*y
    call(sblas_svdv(1.0, z, y, &delta));
    
    //check for breakdown on Lanczos
    if (fabs(delta) < MEPS)//curable
      return sblas_error(sb_BREAKDOWN);
    
    //solve M2*ytil = y
    call(sblas_zerovec(ytil));
    call(sblas_smxv(1.0, M2, False, y, &ytil, False));
    //solve M1^t*ztil = z (for now we set M1 = I)
    call(sblas_svadd(1.0, z, 0.0, ztil));
    
    //compute p and q
    if (it == 0) {
      //p = ytil
      call(sblas_svadd(1.0, ytil, 0.0, p));
      //q = ztil
      call(sblas_svadd(1.0, ztil, 0.0, q));
    }
    else {
      //p = ytil-((qsi*delta)/eps)*p
      call(sblas_svadd(1.0, ytil, -(qsi*delta)/eps, p));
      //q = ztil-((rho*delta)/eps)*q
      call(sblas_svadd(1.0, ztil, -(rho*delta)/eps, q));
    }
    //ptil = A*p
    call(sblas_zerovec(ptil));
    call(sblas_smxv(1.0, A, False, p, &ptil, False));
    //eps = q^T*ptil
    call(sblas_svdv(1.0, q, ptil, &eps));
    //check for breakdown
    if (fabs(eps) < MEPS)
      return sblas_error(sb_BREAKDOWN);
    beta = eps/delta;
    //check for breakdown
    if (fabs(beta) < MEPS)
      return sblas_error(sb_BREAKDOWN);
    //vtil = ptil-beta*v
    call(sblas_svadd(1.0, ptil,0.0, vtil));
    call(sblas_svadd(-beta,v,1.0, vtil));
    //solve M1*y = vtil (for now we set M1 = I)
    call(sblas_svadd(1.0, vtil, 0.0, y));
    //rho_next = |y|
    rho_next = sblas_sv2norm(y);

    //wtil = A^T*q-beta*w
    call(sblas_zerovec(wtil));
    call(sblas_smxv(1.0, A, True, q, &wtil, False));
    call(sblas_svadd(-beta, w,1.0, wtil));
    //solve M2^t*z = wtil (for now we set M2 = I)
    call(sblas_zerovec(z));
    call(sblas_smxv(1.0, M2, True, wtil, &z, False));
    //qsi_next = |z|
    qsi_next = sblas_sv2norm(z);
    //theta = rho_next/(gamma_prev*|beta|)
    theta = rho_next/(gamma_prev*fabs(beta));
    //gamma = 1/sqrt(1+theta^2)
    gamma = 1.0/sqrt(1.0+theta*theta);
    if (fabs(gamma) < MEPS)
      return sblas_error(sb_BREAKDOWN);
    //eta = -eta_prev*rho*gamma^2/(beta*gamma_prev^2)
    eta = -eta_prev*rho*gamma*gamma/(beta*gamma_prev*gamma_prev);
    
    if (it == 0){
      //d = eta*p
      call(sblas_svadd(eta, p, 0.0, d));
      //s = eta*ptil
      call(sblas_svadd(eta, ptil, 0.0, s));
    }
    else {
      //d = eta*p + (theta_prev*gamma)^2*d
      call(sblas_svadd(eta, p, (theta_prev*gamma)*(theta_prev*gamma), d));
      //s = eta*ptil + (theta_prev*gamma)^2*s
      call(sblas_svadd(eta, ptil, (theta_prev*gamma)*(theta_prev*gamma), s));
    }
    //x = x+d
    call(sblas_svadd(1.0, d, 1.0, x));
    //r = r-s
    call(sblas_svadd(-1.0, s, 1.0, r));
    
    rnorm = sblas_sv2norm(r);
    
    //TEMPORARY
    call(sblas_zerovec(rtest));
    call(sblas_smxv(-1.0, A, False, x, &rtest, False));
    //r += b
    call(sblas_svadd(1.0, b, 1.0, rtest));
    
    rnormtrue = sblas_sv2norm(rtest);
    printf("Iteration %d |R| = %1.5e/%1.5e \n",it, rnorm, rnormtrue);
    //printf("%1.5e\n",rnorm/rnorm0);
    if (rnorm/rnorm0 < tol) break;
    
    theta_prev  = theta;
    gamma_prev  = gamma;
    eta_prev    = eta;
    rho         = rho_next;
    qsi         = qsi_next;
  }
  
  //release memory
  sblas_destroysvec(r);
  sblas_destroysvec(d);
  sblas_destroysvec(s);
  sblas_destroysvec(v);
  sblas_destroysvec(w);
  sblas_destroysvec(y);
  sblas_destroysvec(z);
  sblas_destroysvec(p);
  sblas_destroysvec(vtil);
  sblas_destroysvec(wtil);
  sblas_destroysvec(ytil);
  sblas_destroysvec(ztil);
  sblas_destroysvec(rtest);
  sblas_destroysmat(M2);
  
  return sb_OK;
}