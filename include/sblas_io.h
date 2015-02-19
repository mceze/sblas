/*
 *  sblas_io.h
 *  sblas
 *
 *  Created by Marco Ceze on 2/1/11.
 *  Copyright 2011 University of Michigan. All rights reserved.
 *
 */

#ifndef _sblas_io_h
#define _sblas_io_h 1

#include "sblas_def.h"
#include "sblas_struct.h"

/* function:  sblas_ER*/
/* reports the file, line and call  where the error happened*/
int sblas_er( char *file, int line, char *call, int ierr);

/* function: sblas_readsvec */
/* reads ascii  file in sparse matrix format,
 ex: row collumn value\n */
int sblas_readsvec(char *filename, sblas_svec **pV);

/* function: sblas_readsmat */
/* reads ascii  file in sparse matrix format,
 ex: row collumn value\n */
int sblas_readsmat(char *filename, sblas_smat **pM);

/* function: sblas_writesvecascii */
/* writes out a sparse vector in HYPRE format */
int sblas_writesvecascii(char *filename, sblas_svec *V);

/* function: sblas_writesmatascii */
/* writes out a sparse matrix in HYPRE format */
int sblas_writesmatascii(char *filename, sblas_smat *M);

#endif