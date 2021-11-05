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

#ifndef _mutation_h
#define _mutation_h

#include "datatype.h"

#define COLLISIONDIS  0.78

void muta_protein(PROCHAIN *p,RESBACKRUB *backrubs,int mutN[4],MUTATION *mutation[4], int *res_mutindex[4]);
void graft(float main0[][3], float maint[][3],int siden, float side0[][3], float side[][3]);
void get_rotamer_x(float chi[4], float CAx[3], int siden, float side[][3], int muttype);

#endif

