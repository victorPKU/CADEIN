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

#ifndef _protein_h
#define _protein_h

#include "datatype.h"

void read_protein(PROTEIN *p, char *fn);
void clean_protein(PROTEIN *p);
int get_atid(PROCHAIN *p,int ir,char *name);
void get_reslink(PROTEIN *p);

#endif

