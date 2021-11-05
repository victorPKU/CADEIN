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

#ifndef _rotamer_h
#define _rotamer_h

typedef struct{
	float phi;
	float psi;
	float chi[4];
	float sig[4];
	float p;
}bbdepROT;
typedef struct{
	char res[4];
	int chiN;
	int rotN;
	int phipsi_index[37][37][2];
	bbdepROT rotas[80000];
}ROTLIB;

ROTLIB  rotamer_lib[4];
int MAXROTNUM[4];

void get_phipsi_index(ROTLIB  *lib,float phi,float psi,int *start, int *end);
void read_rotamer_lib(char *fn);
void get_maxrotnum();

#endif

