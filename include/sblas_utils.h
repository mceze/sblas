/*
 *  sblas_utils.h
 *  sblas
 *
 *  Created by Marco Ceze on 2/1/11.
 *  Copyright 2011 University of Michigan. All rights reserved.
 *
 */

#ifndef _sblas_utils_h
#define _sblas_utils_h 1

#include "sblas_def.h"
#include "sblas_struct.h"
#include "sblas_enum.h"

/* function: sblas_initsvec */
/* initializes a sblas_svec structure */
int sblas_initsvec(sblas_svec *V, int m);

/* function: sblas_createsvec */
/* creates and initializes a sblas_svec structure */
int sblas_createsvec(sblas_svec **pV, int m);

/* function: sblas_initsmat */
/* initializes a sblas_smat structure */
int sblas_initsmat(sblas_smat *M, int m, int n);

/* function: sblas_createsmat */
/* creates and initializes a sblas_smat structure */
int sblas_createsmat(sblas_smat **pM, int m, int n);

/* function: sblas_destroysvec */
/* destroys svec structure*/
void sblas_destroysvec(sblas_svec *V);

/* function: sblas_destroysmat */
/* destroys smat structure*/
void sblas_destroysmat(sblas_smat *M);

/* function: sblas_svec2dvec */
/* converts a sparse vector into a dense array */
int sblas_svec2dvec(sblas_svec *V, double **pdV);

/* function: sblas_svecentry */
/* add an entry to svec */
int sblas_svecentry(sblas_svec *V, int index, 
                    double value);

/* function: sblas_smatentry */
/* add an entry to smat */
int sblas_smatentry(sblas_smat *M, int row,  
                    int col, double value);

/* function: sblas_svec_getentry */
/* get an entry from vec */
int sblas_svec_getentry(sblas_svec *V, int k, double *value);

/* function: sblas_smat_getentry */
/* get an entry from mat */
int sblas_smat_getentry(sblas_smat *M, int i, int j, double *value);

#endif