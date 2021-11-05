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

#include "datatype.h"

int combine_coordination(COORGEN *coorlib,int coorpairN,COORPAIR *coorpairs,int *coorpairsim,int *coorpairsimindex,COORCOMBINE *combine_coor,MUTATION *mut[4]);
int arrange_combine(int combineN, COORCOMBINE *combine,int *index, int *resindex);

