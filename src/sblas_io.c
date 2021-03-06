/*
 *  sblas_io.c
 *  sblas
 *
 *  Created by Marco Ceze on 2/1/11.
 *  Copyright 2011 University of Michigan. All rights reserved.
 *
 */

#include "sblas_io.h"
#include "sblas_utils.h"

/* function:  sblas_er*/
/* reports the file, line and call  where 
 the error happened*/
int sblas_er( char *file, int line, char *call, int ierr)
{
  if (ierr == sb_OK) return sb_OK;
  printf("**************Error: %d**************",ierr);
  printf("\nFile: %s\nLine: %d\nCall: %s\n",
         file, line, call);
  fflush(stdout);
  return ierr;
}

/* function: sblas_readsvec */
/* reads ascii  file in sparse matrix format,
 ex: row collumn value\n */
int sblas_readsvec(char *filename, sblas_svec **pV)
{
  int ierr, ref, row, col, i, nI;
  char buf[MAX_LINE_LENGTH];
  double value;
  FILE *fid;
  sblas_svec *V;
  
  //scroll through file and count the number of entries
  if ((fid = fopen(filename, "r")) == NULL)
    return sblas_error(sb_READWRITE_ERROR);
  
  i = 0;
  nI = 0;
  while (!feof(fid)) {
    if(fgets(buf, MAX_LINE_LENGTH, fid) != NULL)
      if (strncmp(buf, "\%", 1) != 0){
        sscanf(buf, "%d %d %lf",&row, &col, &value);
        /* ref stores the collumn value of the vector
         that is assumed to be a collum*/
        if (i == 0)
          ref = col;
        else if (col != ref)//input is not a vector
          return sblas_error(sb_INPUT_ERROR);
        //number of rows in the vector is the maximum index
        if (row > nI)
          nI = row;
        i++;
      }
  }
  
  //create vector
  call(sblas_createsvec(&V,nI));
  
  
  rewind(fid);
  while (!feof(fid)) {
    if(fgets(buf, MAX_LINE_LENGTH, fid) != NULL)
      if (strncmp(buf, "\%", 1) != 0){
        //store ordered
        sscanf(buf, "%d %d %lf",&row, &col, &value);
        
        call(sblas_svecentry(V, row-1, value));
        
      }
  }
  
  pV[0] = V;
  
  fclose(fid);
  
  return sb_OK;
  
}

/* function: sblas_readsmat */
/* reads ascii  file in sparse matrix format,
 ex: row collumn value\n */
int sblas_readsmat(char *filename, sblas_smat **pM)
{
  int ierr, row, col, nrow, ncol;
  char buf[MAX_LINE_LENGTH];
  double value;
  FILE *fid;
  sblas_smat *M;
  
  //scroll through file and count the number of entries
  if ((fid = fopen(filename, "r")) == NULL)
    return sblas_error(sb_READWRITE_ERROR);
  
  ncol = nrow = 0;
  while (!feof(fid)) {
    if(fgets(buf, MAX_LINE_LENGTH, fid) != NULL)
      if (strncmp(buf, "\%", 1) != 0){
        sscanf(buf, "%d %d %lf",&row, &col, &value);
        
        if (row > nrow)
          nrow = row;
        if (col > ncol)
          ncol = col;
      }
  }
  
  //create matrix
  call(sblas_createsmat(&M, nrow, ncol));
  
  
  rewind(fid);
  while (!feof(fid)) {
    if(fgets(buf, MAX_LINE_LENGTH, fid) != NULL)
      if (strncmp(buf, "\%", 1) != 0){
        sscanf(buf, "%d %d %lf",&row, &col, &value);
        
        call(sblas_smatentry(M, row-1,
                                           col-1, value));
        
      }
  }
  
  pM[0] = M;
  
  fclose(fid);
  
  return sb_OK;
  
}

/* function: sblas_writesvecascii */
/* writes out a sparse vector in ijv format */
int sblas_writesvecascii(char *filename, sblas_svec *V)
{
  int ierr, i;
  FILE *fid;
  
  if ((fid = fopen(filename, "w")) == NULL)
    return sblas_error(sb_READWRITE_ERROR);
  //print header
  ierr = fprintf(fid,"%% %d %d\n", 1, V->m);
  if (ierr < sb_OK) return sblas_error(sb_READWRITE_ERROR);
  
  for (i = 0; i < V->nZ; i++){
    ierr = fprintf(fid,"%d 1 %1.15e\n", V->index[i]+1,V->val[i]);
    if (ierr < sb_OK) return sblas_error(sb_READWRITE_ERROR);
  }
  ierr = fclose(fid);
  if (ierr != sb_OK) return sblas_error(sb_READWRITE_ERROR);
  
  return sb_OK;
}


/* function: sblas_writesmatascii */
/* writes out a sparse matrix in ijv format */
int sblas_writesmatascii(char *filename, sblas_smat *M,
                         enum sblas_bool ByRow)
{
  int ierr, i, j;
  FILE *fid;
  
  if ((fid = fopen(filename, "w")) == NULL)
    return sblas_error(sb_READWRITE_ERROR);
  //print header
  ierr = fprintf(fid,"%% %d %d\n", M->m, M->n);
  if (ierr < sb_OK) return sblas_error(sb_READWRITE_ERROR);
  
  if (ByRow){
    for (i = 0; i < M->m; i++)
      for (j = 0; j < M->Row[i]->nZ; j++){
        ierr = fprintf(fid,"%d %d %1.15e\n",
                       i+1, M->Row[i]->index[j]+1,
                       M->Row[i]->val[j]);
        if (ierr < sb_OK) return sblas_error(sb_READWRITE_ERROR);
      }
  }
  else {
    for (i = 0; i < M->n; i++)
      for (j = 0; j < M->Col[i]->nZ; j++){
        ierr = fprintf(fid,"%d %d %1.15e\n",
                       M->Col[i]->index[j]+1, i+1,
                       M->Col[i]->val[j]);
        if (ierr < sb_OK) return sblas_error(sb_READWRITE_ERROR);
      }
  }
  ierr = fclose(fid);
  if (ierr != sb_OK) return sblas_error(sb_READWRITE_ERROR);
  
  return sb_OK;
}

























