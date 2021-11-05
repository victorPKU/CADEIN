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

#include "protein.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "geometry.h"
#include "mem.h"

void getxyz(char line[], float xyz[3])
{
	int i;
	for(i=0;i<3;i++)
		xyz[i]=atof(line+8*i);
}

float getrad(char line[])
{
		if(line[0]=='C'){
			if(line[1]==' ')
				return 1.80;
			else if(line[1]=='G'&&(!strncmp(line+4,"ASP",3)||!strncmp(line+4,"ASN",3)||!strncmp(line+4,"HIS",3)||!strncmp(line+4,"PHE",3)||!strncmp(line+4,"TYR",3)||!strncmp(line+4,"TRP",3)))
				return 1.80;
			else if(line[1]=='D'&&(!strncmp(line+4,"GLU",3)||!strncmp(line+4,"GLN",3)))
				return 1.80;
			else if(line[1]=='Z'&&(!strncmp(line+4,"ARG",3)))
				return 1.80;
			else
				return 1.90;
		}
		else if(line[0]=='N'){
			return 1.70;
		}
		else if(line[0]=='O'){
			return 1.60;
		}
		else if(line[0]=='S'){
			return 1.85;
		}
		else
			return 1.85;	
}

#define LINELEN 81
int IsAA(char *line)
{
	if(!strncmp(line, "GLY",3)) return 1;
	if(!strncmp(line, "ALA",3)) return 1;
	if(!strncmp(line, "SER",3)) return 1;
	if(!strncmp(line, "PRO",3)) return 1;
	if(!strncmp(line, "VAL",3)) return 1;
	if(!strncmp(line, "THR",3)) return 1;
	if(!strncmp(line, "CYS",3)) return 1;
	if(!strncmp(line, "ILE",3)) return 1;
	if(!strncmp(line, "LEU",3)) return 1;
	if(!strncmp(line, "ASP",3)) return 1;
	if(!strncmp(line, "ASN",3)) return 1;
	if(!strncmp(line, "GLU",3)) return 1;
	if(!strncmp(line, "GLN",3)) return 1;
	if(!strncmp(line, "LYS",3)) return 1;
	if(!strncmp(line, "HIS",3)) return 1;
	if(!strncmp(line, "MET",3)) return 1;
	if(!strncmp(line, "PHE",3)) return 1;
	if(!strncmp(line, "ARG",3)) return 1;
	if(!strncmp(line, "TYR",3)) return 1;
	if(!strncmp(line, "TRP",3)) return 1;
	if(!strncmp(line, "MSE",3)) return 1;

	return 0;
}
int get_pro_size(FILE *pf, PROCHAIN *p)
{
	int ia,ir;
	char conform=' ';
	char resname[10];
	char line[LINELEN];

	ia=ir=0;

	while(fgets(line,LINELEN,pf)){
		if(!strncmp(line,"ENDMDL",5)) break;
		if(strncmp(line,"ATOM",4)&&strncmp(line,"HETATM",4)) continue;
		if(line[13]=='H'||line[12]=='H') continue;
		if(IsAA(line+17)==0) continue;
		if(line[16]!=' '){
			if(conform==' ')
				conform=line[16];
			else if(line[16]!=conform)
				continue;
		}
		if(ir==0||strncmp(line+17,resname,9)){
			strncpy(resname,line+17,9);
			ir++;
		}
		ia++;
	}

	p->atmn=ia;
	p->resn=ir;

	return ia;
}
void read_structure(FILE *pf, PROCHAIN *p,int count)
{
	int ia,ir;
	char conform=' ';
	char line[LINELEN];

	ia=ir=0;

	while(fgets(line,LINELEN,pf)){
		if(!strncmp(line,"ENDMDL",5)) break;
		if(strncmp(line,"ATOM",4)&&strncmp(line,"HETATM",4)) continue;
		if(line[13]=='H'||line[12]=='H') continue;
		if(IsAA(line+17)==0) continue;
		if(line[16]!=' '){
			if(conform==' ')
				conform=line[16];
			else if(line[16]!=conform)
				continue;
		}
		if(ir==0||strncmp(line+17,p->r[ir-1].name,9)){
			strncpy(p->r[ir].name,line+17,9);
			p->r[ir].start=ia;
			p->r[ir].end=ia;
			p->r[ir].sstart=ia;
			p->r[ir].send=ia;
			p->r[ir].mutability=0;
			ir++;
		}
		getxyz(line+30,p->a[ia].x);
		p->a[ia].vdw_rad=getrad(line+13);
		strncpy(p->a[ia].name,line+13,14);
		
		if(!strncmp(p->a[ia].name,"CA ",3)){
                   if(line[70]=='1') p->r[ir-1].mutability=1;
		}
		p->a[ia].at_res=ir-1;
		if(!strncmp(p->a[ia].name,"N  ",3)||!strncmp(p->a[ia].name,"CA ",3)||!strncmp(p->a[ia].name,"C  ",3)||!strncmp(p->a[ia].name,"O  ",3))
			p->r[ir-1].sstart++;
		if(strncmp(p->a[ia].name,"OXT",3))
			p->r[ir-1].send++;
		p->r[ir-1].end++;
		ia++;
		if(ia==count) break;
	}
	p->atmn=ia;
	p->resn=ir;
}
void read_protein_model(FILE *pf, PROCHAIN *p)
{
	fpos_t ptr;
	int count;

	fgetpos(pf,&ptr);
	count=get_pro_size(pf,p);
	BStd_new(p->a, p->atmn);
	BStd_new(p->r, p->resn);
	fsetpos(pf,&ptr);
	read_structure(pf,p,count);
}
void read_protein(PROTEIN *p, char *fn)
{
	FILE *pf;
	if((pf=fopen(fn,"r"))==NULL){
		printf("Can not open protein structure file %s\n",fn);
		exit(0);
	}

	read_protein_model(pf,&(p->pro));
	strcpy(p->pdbname, fn);

	fclose(pf);
}
void clean_protein(PROTEIN *p)
{
	BStd_free(p->pro.a);
	BStd_free(p->pro.r);
}
int get_atid(PROCHAIN *p,int ir,char *name)
{
	int ia;

	for(ia=p->r[ir].start;ia<p->r[ir].end;ia++){
		if(!strncmp(p->a[ia].name, name,3))
			return ia;
	}
	return -1;
}
#define BLDIS2 4.0*4.0

void get_reslink(PROTEIN *p)
{
	int ir;
	PROCHAIN *pro;
	
	pro=&(p->pro);
	for(ir=0;ir<pro->resn;ir++){
		if(pro->r[ir].end-pro->r[ir].start<4){
			pro->r[ir].pre=-1;
			pro->r[ir].next=-1;
			if(ir!=0) pro->r[ir-1].next=-1;
			if(ir!=pro->resn-1){
				pro->r[ir+1].pre=-1;
				ir++;
			}
			continue;
		}
		if(ir==0){
			pro->r[ir].pre=-1;
			continue;
		}
		if(distance2(pro->a[pro->r[ir].start].x,pro->a[pro->r[ir-1].start].x)<=BLDIS2){
			pro->r[ir].pre=ir-1;
			pro->r[ir-1].next=ir;
		}
		else{
			pro->r[ir].pre=-1;
			pro->r[ir-1].next=-1;
		}
		if(ir==pro->resn-1){
			pro->r[ir].next=-1;
		}
	}
}
