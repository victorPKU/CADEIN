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

#ifndef _backrub_h
#define _backrub_h

#include "datatype.h"

void backrub_protein(PROCHAIN *p, RESBACKRUB *backrubs);
int resi_backrub(float CApx[3], float CAnx[3],float mainatx[7][3], float backrub_mainatx[7][3], float theta);

#endif

