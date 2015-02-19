/*
 *  sblas_struct.h
 *  sblas
 *
 *  Created by Marco Ceze on 2/1/11.
 *  Copyright 2011 University of Michigan. All rights reserved.
 *
 */

#ifndef _sblas_struct_h
#define _sblas_struct_h 1

#include "sblas_def.h"

/* structure: sblas_svec */
/* sparse vector structure */
typedef struct
{
  int m;                   /* dimension */
  int nZ;                  /* number of non-zero entries */
  int nZprealloc;          /* preallocation size */
  int *index;              /* coordinates of the entries */
  double *val;             /* values */
}
sblas_svec;

/* structure: sblas_smat */
/* sparse matrix structure stored in csr and csc 
 form so operations become faster */
typedef struct
{
  int m;                   /* number of rows */
  int n;                   /* number of collumns */
  int nZ;                  /* number of non-zero entries */
  sblas_svec **Row;        /* one sparse vector per row */
  sblas_svec **Col;        /* one sparse vector per col */
}
sblas_smat;

#endif // end ifndef _sblas_struct_h