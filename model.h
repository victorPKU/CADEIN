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


#ifndef _MODEL_H_
#define _MODEL_H_


typedef struct{
	char resnum[6];
	char atname[8];
	float x[3];
}MODELATM;

typedef struct{
	int atn;
	MODELATM at[100];
	float cadmium[3];
}MODEL;


#endif /* MODEL_H_ */
