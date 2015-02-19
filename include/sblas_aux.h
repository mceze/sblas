/*
 *  sblas_aux.h
 *  sblas
 *
 *  Created by Marco Ceze on 2/1/11.
 *  Copyright 2011 University of Michigan. All rights reserved.
 *
 */

#ifndef _sblas_aux_h
#define _sblas_aux_h 1

#include "sblas_def.h"

/* function: sblas_bsearch */
/* finds the position between "begin" and "end" of 
 integer "target" in incresing ordered "set" */
int sblas_bsearch(const int target, const int *set, 
                     int begin, int end, int *rank);

#endif