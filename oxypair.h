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

int gen_coorpairs(PROCHAIN *p,RESBACKRUB *backrubs,MUTATION *m[4], COORGEN *coorlib, int *res_coorindex, COORPAIR *coorpairs,int *coorpairindex);
void simple_coorlink(int resN, int coorN, COORGEN *coorlib,int *coorpairN,COORPAIR *coorpairs,int *coorpairindex,int *coorpairsim,int *coorpairsimindex);


