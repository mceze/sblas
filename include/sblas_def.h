/*
 *  sblas_def.h
 *  sblas
 *
 *  Created by Marco Ceze on 2/1/11.
 *  Copyright 2011 University of Michigan. All rights reserved.
 *
 */

#ifndef _sblas_def_h
#define _sblas_def_h 1

/* standard headers used */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <math.h>

/* definitions */
#define MEPS 1e-16 // machine precision
#define MAX_LINE_LENGTH 400 //maximum line length
#define MAX_STR_LENGTH 100 //maximum string length

/* simple math macros */
#define max(a,b)  (((a) > (b)) ? (a) : (b))
#define min(a,b)  (((a) < (b)) ? (a) : (b))
#define sign(a) (((a) < 0.0) ? -1 : 1)
#define swap(a,b, t)  {t = a; a = b; b = t;}
#define comp(a,b)  (((a) < (b)) ? -1 : 1)

/* Error Macro: used to report error occurrences */
#define sblas_error(X) (sblas_er( __FILE__, __LINE__, #X, (X)))
#define call(X) do{ierr = sblas_error(X); if(ierr != sb_OK) return ierr;} while(0)

/* Error codes */
#define sb_OK               0
#define sb_READWRITE_ERROR -1
#define sb_INCONSISTENCY   -2
#define sb_INPUT_ERROR     -3
#define sb_INCOMPATIBLE    -4
#define sb_MEMORY_ERROR    -5
#define sb_NOT_FOUND       -6
#define sb_OUT_OF_BOUNDS   -7
#define sb_BREAKDOWN       -8



#endif // end ifndef _sblas_def_h