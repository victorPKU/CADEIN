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

/****************************************************************************************
 *                               AutoMatch  
 * Potential active site searching for enzyme design and target binding protein design  *
 *                      Written by Zhang Changsheng 
 * Copyright (c) 2011 Molecular design lab in Peking university, Beijing, China  *
 *
 *   contact: Author: cszhang@mdl.ipc.pku.edu.cn,
 *            Head of the lab: professor Lai :lhlai@pku.edu.cn  *
 *   http://mdl.ipc.pku.edu.cn
 *
 *  Free for academical usage and please reference to:
 *  Zhang changsheng, Lai luhua. AutoMatch:Potential active site searching for enzyme  
 *    design and target binding protein design. Journal of computational chemistry. 2011. *
 *
 * rotamer.c : read rotamer library file, get rotamers for a particular residue with angles phi,psi.
 ******************************************************************************************/
#include "rotamer.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

void set_phipsi_start(ROTLIB  *lib,float phi,float psi,int start)
{
	int iphi,ipsi;
	
	iphi=(int)((phi+185)/10)%37;
	ipsi=(int)((psi+185)/10)%37;
	lib->phipsi_index[iphi][ipsi][0]=start;
}
void set_phipsi_end(ROTLIB  *lib,float phi,float psi,int end)
{
	int iphi,ipsi;
	
	iphi=(int)((phi+185)/10)%37;
	ipsi=(int)((psi+185)/10)%37;
	lib->phipsi_index[iphi][ipsi][1]=end;
}
void get_phipsi_index(ROTLIB  *lib,float phi,float psi,int *start, int *end)
{
	int iphi,ipsi;
	
	iphi=(int)((phi+185)/10)%37;
	ipsi=(int)((psi+185)/10)%37;
	(*start)= lib->phipsi_index[iphi][ipsi][0];
	(*end)=lib->phipsi_index[iphi][ipsi][1];
}

int get_chiN(char *line)
{
	if(line[3]=='0') return 1;
	if(line[6]=='0') return 2;
	if(line[9]=='0') return 3;
	return 4;
}
int which_res(char *line)
{
	if(!strncmp("GLU",line,3))	return 0;
	if(!strncmp("ASP",line,3))	return 1;
	if(!strncmp("HIS",line,3))	return 2;
	if(!strncmp("CYS",line,3))	return 3;

	return -1;
}

#define RLIBFILELINE  120

void read_rotamer_lib(char *fn)
{
	FILE *pf;
	char line[RLIBFILELINE];
	char pre_line[RLIBFILELINE]="";
	float prob;
	int i=0;
	int res=4;

	if((pf=fopen(fn,"r"))==NULL){
		printf("ERROR: Can not open rotamer library file %s!\n",fn);
		exit(0);
	}
	for(i=0;i<4;i++){
		rotamer_lib[i].chiN=0;
		rotamer_lib[i].rotN=0;
	}
	while(fgets(line,RLIBFILELINE,pf)){
		if(line[0]=='#') continue;
		if(strncmp(line,pre_line,3)){
			if(res<4)
				set_phipsi_end(rotamer_lib+res,rotamer_lib[res].rotas[rotamer_lib[res].rotN-1].phi,rotamer_lib[res].rotas[rotamer_lib[res].rotN-1].psi,rotamer_lib[res].rotN-1);
			res=which_res(line);
			strncpy(rotamer_lib[res].res,line,3);
			rotamer_lib[res].chiN=get_chiN(line+25);
			rotamer_lib[res].rotN=0;
		}
		rotamer_lib[res].rotas[rotamer_lib[res].rotN].phi=atof(line+5);
		rotamer_lib[res].rotas[rotamer_lib[res].rotN].psi=atof(line+10);
		prob=atof(line+35);
		if(prob<1e-6) prob=-6.0;
		else prob=log10(prob);
		rotamer_lib[res].rotas[rotamer_lib[res].rotN].p=prob;
		rotamer_lib[res].rotas[rotamer_lib[res].rotN].chi[0]=atof(line+45);
		rotamer_lib[res].rotas[rotamer_lib[res].rotN].chi[1]=atof(line+53);
		rotamer_lib[res].rotas[rotamer_lib[res].rotN].chi[2]=atof(line+61);
		rotamer_lib[res].rotas[rotamer_lib[res].rotN].chi[3]=atof(line+69);
		if(strncmp(line,pre_line,14)){
			set_phipsi_start(rotamer_lib+res,rotamer_lib[res].rotas[rotamer_lib[res].rotN].phi,rotamer_lib[res].rotas[rotamer_lib[res].rotN].psi,rotamer_lib[res].rotN);
			if(rotamer_lib[res].rotN!=0)
				set_phipsi_end(rotamer_lib+res,rotamer_lib[res].rotas[rotamer_lib[res].rotN-1].phi,rotamer_lib[res].rotas[rotamer_lib[res].rotN-1].psi,rotamer_lib[res].rotN-1);
		}
		rotamer_lib[res].rotN++;
		strcpy(pre_line,line);
	}
	set_phipsi_end(rotamer_lib+res,rotamer_lib[res].rotas[rotamer_lib[res].rotN-1].phi,rotamer_lib[res].rotas[rotamer_lib[res].rotN-1].psi,rotamer_lib[res].rotN-1);

	fclose(pf);
}
void get_maxrotnum()
{
	int r,i,j;
	int n;

	for(r=0;r<4;r++){
		MAXROTNUM[r]=0;
		for(i=0;i<37;i++)
			for(j=0;j<37;j++){
				n=rotamer_lib[r].phipsi_index[i][j][1]-rotamer_lib[r].phipsi_index[i][j][0]+1;
				if(MAXROTNUM[r]<n)
					MAXROTNUM[r]=n;
			}
	}
}
