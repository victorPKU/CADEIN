/***************************************************************************
 *                             CADEIN
 *  Program for cadmium binding protein design   *
 *         Author: Dr. Zhang Changsheng 
 * Copyright (c) 2021 Molecular Design Lab at Peking University  *
 *
 *   contact: changshengzhang@pku.edu.cn
 *
 *  Free for academical purpose and please reference to:
 *  Tailoring Escherichia coli Chemotactic Sensing towards Cadmium 
 *  by Computational Redesign of Ribose-binding Protein. Hengyi Li, 
 *  Changsheng Zhang, Xi Chen and Luhua Lai. mSystems. 2021. *
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "mem.h"

void *BStd_calloc(char *name, unsigned long long nelem,unsigned elsize)
{
  void *p;
  
  p=NULL;
  if ((nelem==0)||(elsize==0))
    p=NULL;
  else
    {
      if ((p=calloc((size_t)nelem,(size_t)elsize))==NULL){
		  printf(" Memery Alloc Error for %s\n", name);
		 exit(0); 
	  }
     }
  return p;
}

void *BStd_realloc(char *name, void *ptr,unsigned long long size)
{
  void *p;
  
  p=NULL;
  if (size==0)
    p=NULL;
  else
    {
      if (ptr==NULL) 
	p=malloc((size_t)size); 
      else 
	p=realloc(ptr,(size_t)size);
      if (p==NULL){ 
		 printf(" Memery Alloc Error for %s\n", name);
		 exit(0); 
	  }
    }
  return p;
}

void BStd_free(void *ptr)
{
  if (ptr != NULL)
    free(ptr);
}
