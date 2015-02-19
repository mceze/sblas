/*
 *  sblas_aux.c
 *  sblas
 *
 *  Created by Marco Ceze on 2/1/11.
 *  Copyright 2011 University of Michigan. All rights reserved.
 *
 */

#include "sblas_aux.h"

/* function: sblas_bsearch */
/* finds the position between "begin" and "end" of 
 integer "target" in incresing ordered "set" */
int sblas_bsearch(const int target, const int *set, 
                 int begin, int end, int *rank)
{
  int index, center;
  int size;
  
  if (end < begin)
    return INPUT_ERROR;
  
  size = end-begin+1;
  if (size == 1){
    if (target == set[begin]){
      if (rank != NULL)
        (*rank) = begin;
      return begin;
    }
    else{
      if (rank != NULL){
        if (target < set[begin])
          (*rank) = begin-1;
        else
          (*rank) = end+1;
      }
      return NOT_FOUND;
    }
  }
  if (target == set[begin]){
    index = begin;
    if (rank != NULL)
      (*rank) = begin;
  }
  
  if (target == set[end]){
    index = end;
    if (rank != NULL)
      (*rank) = end;
  }
  
  //if within range
  if (target > set[end]){
    if (rank != NULL)
      (*rank) = end+1;
    index = NOT_FOUND;
  }
  else if (target < set[begin]){
    if (rank != NULL)
      (*rank) = begin-1;
    index = NOT_FOUND;
  }
  else{
    while (size > 1){
      center = begin + size/2-1;
      if (target > set[center]){
        //in top portion of set
        begin = center+1;
      }
      else {
        end = center;
      }
      size = end-begin+1;
    }
    if (target == set[begin]){
      index = begin;
      if (rank != NULL)
        (*rank) = begin;
    }
    else{
      index = NOT_FOUND;
      if (target > set[begin]){
        if (rank != NULL)
          (*rank) = begin+1;
      }
      else{
        if (rank != NULL)
          (*rank) = begin;
      }
    }
  }
  
  return index;
}

