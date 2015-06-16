/*
 *  sblas_utils.c
 *  sblas
 *
 *  Created by Marco Ceze on 2/1/11.
 *  Copyright 2011 University of Michigan. All rights reserved.
 *
 */

#include "sblas_def.h"
#include "sblas_utils.h"
#include "sblas_aux.h"
#include "sblas_io.h"

/* function: sblas_initsvec */
/* initializes a sblas_svec structure */
int sblas_initsvec(sblas_svec *V, int m)
{
  if (V == NULL) return sb_INPUT_ERROR;
  
  V->m = m;
  V->nZ = 0;
  V->nZprealloc = 10;
  V->index = malloc(V->nZprealloc*sizeof(int));
  memset(V->index, 0, V->nZprealloc*sizeof(int));
  V->val = malloc(V->nZprealloc*sizeof(double));
  memset(V->val, 0, V->nZprealloc*sizeof(double));
  
  return sb_OK;
}

/* function: sblas_createsvec */
/* creates and initializes a sblas_svec structure */
int sblas_createsvec(sblas_svec **pV, int m)
{
  int ierr;
  
  //allocate and initialize the vector
  if (((*pV) = malloc(sizeof(sblas_svec))) == NULL)
    return sb_MEMORY_ERROR;
  
  ierr = sblas_error(sblas_initsvec((*pV), m));
  if (ierr != sb_OK) return ierr;
  
  return sb_OK;
}

/* function: sblas_initsmat */
/* initializes a sblas_smat structure */
int sblas_initsmat(sblas_smat *M, int m, int n)
{
  int ierr, i;
  if (M == NULL) return sb_INPUT_ERROR;
  
  M->m = m;
  M->n = n;
  M->nZ = 0;
  //init rows and collumns
  if ((M->Row = malloc(M->m*sizeof(sblas_svec))) == NULL)
    return sb_MEMORY_ERROR;
  for (i = 0; i < M->m; i++) {
    ierr = sblas_error(sblas_createsvec(&(M->Row[i]), n));
    if (ierr != sb_OK) return ierr;
  }
  if ((M->Col = malloc(M->n*sizeof(sblas_svec))) == NULL)
    return sb_MEMORY_ERROR;
  for (i = 0; i < M->n; i++) {
    ierr = sblas_error(sblas_createsvec(&(M->Col[i]), m));
    if (ierr != sb_OK) return ierr;
  }
  
  return sb_OK;
}

/* function: sblas_createsmat */
/* creates and initializes a sblas_smat structure */
int sblas_createsmat(sblas_smat **pM, int m, int n)
{
  int ierr;
  
  //allocate and initialize the matrix
  if (((*pM) = malloc(sizeof(sblas_smat))) == NULL)
    return sb_MEMORY_ERROR;
  
  ierr = sblas_error(sblas_initsmat((*pM), m, n));
  if (ierr != sb_OK) return ierr;
  
  return sb_OK;
}

/* function: sblas_destroysvec */
/* destroys svec structure*/
void sblas_destroysvec(sblas_svec *V)
{
  if (V != NULL){
    free(V->index);
    free(V->val);
    free(V);
  }
}

/* function: sblas_destroysmat */
/* destroys smat structure*/
void sblas_destroysmat(sblas_smat *M)
{
  int i;
  
  if (M != NULL){
    for (i = 0; i < M->m; i++)
      sblas_destroysvec(M->Row[i]);
    free(M->Row);
    for (i = 0; i < M->n; i++)
      sblas_destroysvec(M->Col[i]);
    free(M->Col);
    
    free(M);
  }
}

/* function: sblas_svec2dvec */
/* converts a sparse vector into a dense array */
int sblas_svec2dvec(sblas_svec *V, double **pdV)
{
  int i, z;
  
  if (V == NULL)
    return sb_INPUT_ERROR;
  
  //allocate array
  if (((*pdV) = malloc(V->m*sizeof(double))) == NULL)
    return sb_MEMORY_ERROR;
  
  z = 0;
  for (i = 0; i < V->m; i++){
    if (i+1 == V->index[z]){
      (*pdV)[V->index[z]] = V->val[z];
      z++;
    }
    else
      (*pdV)[i] = 0.0;
  }
  
  return sb_OK;
}

/* function: sblas_svecentry */
/* add an entry to svec */
int sblas_svecentry(sblas_svec *V, int index, 
                    double value)
{
  int z, dest = 0, src = 0, movesize, rank;
  
  if (V == NULL)
    return sblas_error(sb_INPUT_ERROR);
  
  //return immediately if value is machine zero
  if (fabs(value) <= MEPS)
    return sb_OK;
  
  //change size if necessary
  if (index >= V->m)
    V->m = index+1;
  
  /* adding entry for the first time? */
  if (V->nZ == 0){
    V->nZ++;
    if (V->nZ > V->nZprealloc){
      V->nZprealloc = min(2*(V->nZ+10),V->m);
      if((V->index = realloc(V->index, V->nZprealloc
                             *sizeof(int))) == NULL)
        return sblas_error(sb_MEMORY_ERROR);
      if((V->val = realloc(V->val, V->nZprealloc
                           *sizeof(double))) == NULL)
        return sblas_error(sb_MEMORY_ERROR);
    }
    V->index[0] = index;
    V->val[0] = value;
  }
  else{
    /* is it a new entry */
    if ((z = sblas_bsearch(index, V->index, 0, 
                           V->nZ-1, &rank)) == sb_NOT_FOUND){
      /* make sure that there is enough space in the 
       arrays for the new entry */
      if (V->nZ+1 > V->nZprealloc){
        V->nZprealloc = V->nZ+1;
        if((V->index = realloc(V->index, V->nZprealloc
                               *sizeof(int))) == NULL)
          return sblas_error(sb_MEMORY_ERROR);
        if((V->val = realloc(V->val, V->nZprealloc
                             *sizeof(double))) == NULL)
          return sblas_error(sb_MEMORY_ERROR);
      }
      //make sure to keep the crescent order
      if (rank == V->nZ){
        movesize = 0;
      }
      else if (rank == -1){
        //move all on position forward
        rank = 0;
        src = 0;
        dest = 1;
        movesize = V->nZ;
      }
      else {
        src = rank;
        dest = rank+1;
        movesize = V->nZ-rank;
      }
      
      if (movesize > 0){
        if (memmove(V->index+dest,V->index+src,
                    movesize*sizeof(int)) == NULL)
          return sblas_error(sb_MEMORY_ERROR);
        if (memmove(V->val+dest,V->val+src,
                    movesize*sizeof(double)) == NULL)
          return sblas_error(sb_MEMORY_ERROR);
      }
      V->index[rank] = index;
      V->val[rank] = value;
      V->nZ++;
    }
    else {//not a new entry
      V->val[z] += value;
      //check for cancelation
      if (fabs(V->val[z]) <= MEPS){
        src = z+1;
        dest = z;
        movesize = V->nZ-src;
        if (memmove(V->index+dest,V->index+src,
                    movesize*sizeof(int)) == NULL)
          return sblas_error(sb_MEMORY_ERROR);
        if (memmove(V->val+dest,V->val+src,
                    movesize*sizeof(double)) == NULL)
          return sblas_error(sb_MEMORY_ERROR);
        V->nZ--;        
      }
    }
  }
  
  return sb_OK;
}

/* function: sblas_smatentry */
/* add an entry to smat */
int sblas_smatentry(sblas_smat *M, int row,  
                    int col, double value)
{
  int ierr;
  
  if (M == NULL)
    return sblas_error(sb_INPUT_ERROR);
  
  //return immediately if value is machine zero
  if (fabs(value) <= MEPS)
    return sb_OK;
  
  //check range
  if (row > M->m || col > M->n)
    return sblas_error(sb_OUT_OF_BOUNDS);
  
  //store nz before and then add the new nz for the row
  M->nZ -= M->Row[row]->nZ;
  
  ierr = sblas_error(sblas_svecentry(M->Row[row], col, value));
  if (ierr != sb_OK) return ierr;
  
  M->nZ += M->Row[row]->nZ;
  
  ierr = sblas_error(sblas_svecentry(M->Col[col], row, value));
  if (ierr != sb_OK) return ierr;
  
  return sb_OK;
}

