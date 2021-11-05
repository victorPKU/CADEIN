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

#define BStd_new(ptr,nelem) (ptr)=BStd_calloc(#ptr,(nelem),sizeof(*(ptr)))
#define BStd_renew(ptr,nelem) (ptr)=BStd_realloc(#ptr,(ptr),(nelem)*sizeof(*(ptr)))

void *BStd_calloc(char *name, unsigned long long nelem,unsigned elsize);
void *BStd_realloc(char *name, void *ptr,unsigned long long size);
void  BStd_free( void *ptr);
