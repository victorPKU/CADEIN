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

void record_result(RECORD *r,PROTEIN *p,COORGEN *coorlib,MUTATION *m[4],COORCOMBINE *combine,int *index, int resultN, int *resindex);
void write_result(RECORD *rec,char *fn);
void write_temp_result(RECORD *rec,char *fn);
int read_result(RECORD *rec,char *fn);

