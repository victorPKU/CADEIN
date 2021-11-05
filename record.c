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

#include "record.h"
#include "mutation.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mem.h"

void record_cadmium(RECORD *r,float cadmium[3], int maxscore)
{
	int i;

	for(i=0;i<3;i++)
		r->cadmium[i]=cadmium[i];
	r->coor_score=maxscore;
}
void record_coordintation(RECORD *r,PROCHAIN *p, COORCOMBINE *combine, COORGEN *coorlib,MUTATION *m[4])
{
	int ioxy;
	int restyp,coortyp;
	int source;

	r->res_n=combine->coorresN;
	r->coorn=combine->coorN;
	
	for(ioxy=0;ioxy<combine->coorresN;ioxy++){
		restyp=coorlib[combine->coor[ioxy]].muttype;
		coortyp=coorlib[combine->coor[ioxy]].coorat;
		source=coorlib[combine->coor[ioxy]].source;
		r->res[ioxy]=coorlib[combine->coor[ioxy]].res;
		strcpy(r->res_name[ioxy], p->r[r->res[ioxy]].name);
		r->res_type[ioxy]=restyp;
		r->coor_at[ioxy]=coortyp;
		if(restyp!=-1){
			r->res_rub[ioxy]=m[restyp][source].backrub;
			r->rotamer[ioxy]=m[restyp][source].rota;
			r->rotshift[ioxy]=m[restyp][source].rotshift;
		}
		else{
			r->res_rub[ioxy]=-1;
			r->rotamer[ioxy]=-1;
			r->rotshift[ioxy]=-1;
		}
	}
}
int record_combine(RECORD *rec,int score)
{
	struct LIST *p,*q;
	static int ren=0;
	int index;

	p=q=NULL;
	if(ren>=MAXRECORD && score<=rec[order->i].coor_score)
		return -1;
	if(order==NULL) {
		BStd_new(order,1);
		order->n=NULL;
		p=order;
	}
	else{
		p=order;
		while(p!=NULL){
			if(score>rec[p->i].coor_score) {q=p;p=p->n;}
			else if(p!=order){
				BStd_new(q->n,1);
				q->n->n=p;
				p=q->n;
				break;
			}
			else{
				BStd_new(order,1);
				order->n=p;
				p=order;
				break;
			}
		}
		if(p==NULL){
			BStd_new(q->n,1);
			q->n->n=NULL;
			p=q->n;
		}
	}
	if(ren<MAXRECORD){
		index=ren;
		p->i=ren;
		ren++;
	}
	else{
		index=order->i;
		p->i=order->i;
		p=order->n;
		BStd_free(order);
		order=p;		
	}

	return index;
}
void record_result(RECORD *r,PROTEIN *p,COORGEN *coorlib,MUTATION *m[4],COORCOMBINE *combine,int *index, int resultN, int *resindex)
{
	int score,maxscore;
	int max_index;
	int resi,i;
	int record_index;

	for(resi=0;resi<resultN;resi++){
		maxscore=0;max_index=index[resindex[resi]];
		for(i=resindex[resi];i<resindex[resi+1];i++){
			score=combine[index[i]].combine_score;
			if(score>maxscore){
				maxscore=score;
				max_index=index[i];
			}
		}
		/*printf("%5d %5d %5d %5d\n",resi,maxscore,max_index,resultN);*/
		record_index=record_combine(r,maxscore);
		/*printf("record_combine %5d\n",record_index);*/
		if(record_index!=-1){
			strcpy(r[record_index].pdbname,p->pdbname);
			/*printf("%s %s\n",r->pdbname,p->pdbname);*/
			record_cadmium(r+record_index,combine[max_index].Cd,maxscore);
			/*printf("record_cadmium\n");*/
			record_coordintation(r+record_index,&(p->pro),combine+max_index, coorlib,m);
			/*printf("record_coordintation\n");*/
		}	
	}
}
int arrange_record_0(int result_index[])
{
	struct LIST *p,*q;
	int i=0;

	for(p=order;p!=NULL;){
		result_index[i]=p->i;
		q=p->n;
		p=q;
		i++;
	}
	return i;
}

int arrange_record(int result_index[])
{
	struct LIST *p,*q;
	int i=0;

	for(p=order;p!=NULL;){
		result_index[i]=p->i;
		q=p->n;
		BStd_free(p);
		p=q;
		i++;
	}
	return i;
}

char type2char(int type)
{
	if(type==-1) return 'M';
	if(type==0) return 'E';
	if(type==1) return 'D';
	if(type==2) return 'H';
	if(type==3) return 'C';
	return ' ';
}
char *coordinate2char(int type)
{
	if(type==-1) return "O  ";
	if(type==MUT2GLU) return "OE ";
	if(type==MUT2GLUO1) return "OE1";
	if(type==MUT2GLUO2) return "OE2";
	if(type==MUT2ASP) return "OD ";
	if(type==MUT2ASPO1) return "OD1";
	if(type==MUT2ASPO2) return "OD2";
	if(type==MUT2HISN1) return "ND1";
	if(type==MUT2HISN2) return "NE2";
	if(type==MUT2CYS) return "SG ";

	return "  ";
}

void write_result(RECORD *rec,char *fn)
{
	FILE *rf;
	int ii,i,j;
	int index[MAXRECORD];
	int n;

	n=arrange_record(index);
	if((rf=fopen(fn,"w"))==NULL)
		printf("Can't open %s",fn);

	for(ii=0;ii<n;ii++){
		i=index[n-ii-1];
		/*printf("pdbname %d,%s\n",i,rec[i].pdbname);*/
		fprintf(rf,"%-6d %8.3f%8.3f%8.3f --> %3d %s\n",ii+1,
			rec[i].cadmium[0],rec[i].cadmium[1],rec[i].cadmium[2],
			rec[i].coor_score,rec[i].pdbname);
		for(j=0;j<rec[i].res_n;j++){
			if(rec[i].res_type[j]==-1)
				fprintf(rf,"    %-4d   %c(%s) [%s]\n", j+1,type2char(rec[i].res_type[j]),coordinate2char(rec[i].coor_at[j]),rec[i].res_name[j]);
			else{
				fprintf(rf,"    %-4d   %c(%s) [%s] %4d %4d %4d\n", j+1,type2char(rec[i].res_type[j]),coordinate2char(rec[i].coor_at[j]),rec[i].res_name[j],rec[i].res_rub[j]%BACKRUBAN,rec[i].rotamer[j],rec[i].rotshift[j]);
				/*printf("    %-4d   %c(%s) [%s] %4d %4d \n", j+1,type2char(rec[i].res_type[j]),coordinate2char(rec[i].coordinate_atom[j]),rec[i].res_name[j],rec[i].res_rub[j]%5,rec[i].rotamer[j]);*/
			}
		}
	}
	
	fclose(rf);
}
void write_temp_result(RECORD *rec,char *fn)
{
	FILE *rf;
	int ii,i,j;
	int index[MAXRECORD];
	int n;

	n=arrange_record_0(index);
	if((rf=fopen(fn,"w"))==NULL)
		printf("Can't open %s",fn);

	for(ii=0;ii<n;ii++){
		i=index[n-ii-1];
		fprintf(rf,"%-6d %8.3f%8.3f%8.3f --> %3d %s\n",ii+1,
			rec[i].cadmium[0],rec[i].cadmium[1],rec[i].cadmium[2],
			rec[i].coor_score,rec[i].pdbname);
		for(j=0;j<rec[i].res_n;j++){
			if(rec[i].res_type[j]==-1)
				fprintf(rf,"    %-4d   %c(%s) [%s]\n", j+1,type2char(rec[i].res_type[j]),coordinate2char(rec[i].coor_at[j]),rec[i].res_name[j]);
			else{
				fprintf(rf,"    %-4d   %c(%s) [%s] %4d %4d %4d \n", j+1,type2char(rec[i].res_type[j]),coordinate2char(rec[i].coor_at[j]),rec[i].res_name[j],rec[i].res_rub[j]%BACKRUBAN,rec[i].rotamer[j],rec[i].rotshift[j]);
				/*printf("    %-4d   %c(%s) [%s] %4d %4d \n", j+1,type2char(rec[i].res_type[j]),coordinate2char(rec[i].coordinate_atom[j]),rec[i].res_name[j],rec[i].res_rub[j]%5,rec[i].rotamer[j]);*/
			}
		}
	}
	
	fclose(rf);
}

#define RECORDLINELEN  200
int read_result(RECORD *rec,char *fn)
{
	FILE *rf;
	int c=-1;
	int i;
	int ci;
	char line[RECORDLINELEN];

	if((rf=fopen(fn,"r"))==NULL){
		printf("Can't open %s",fn);
		exit(0);
	}
	while(fgets(line,RECORDLINELEN,rf)){
		if(line[0]!=' '){
			if(c>=MAXRECORD-1) break;
			c++;
			for(i=0;i<3;i++)
				rec[c].cadmium[i]=atof(line+7+8*i);
			rec[c].res_n=0;
			ci=40;
			while(line[ci]!='\n'&&line[ci]!='\0'){
				rec[c].pdbname[ci-40]=line[ci]; ci++;
			}
			rec[c].pdbname[ci-40]='\0';
		}
		else{
			if(line[11]=='M'){
				rec[c].res_type[rec[c].res_n]=-1;
			}
			else if(line[11]=='E'){
				rec[c].res_type[rec[c].res_n]=MUT2GLU;
			}
			else if(line[11]=='D'){
				rec[c].res_type[rec[c].res_n]=MUT2ASP;
			}
			else if(line[11]=='H'){
				rec[c].res_type[rec[c].res_n]=MUT2HIS;
			}
			else if(line[11]=='C'){
				rec[c].res_type[rec[c].res_n]=MUT2CYS;
			}
			strncpy(rec[c].res_name[rec[c].res_n],line+19,9);
			rec[c].res_name[rec[c].res_n][9]='\0';
			rec[c].res_rub[rec[c].res_n]=atoi(line+31);
			rec[c].rotamer[rec[c].res_n]=atoi(line+34);
			rec[c].rotshift[rec[c].res_n]=atoi(line+43);
			rec[c].res_n++;
		}
	}

	fclose(rf);
	
	return c+1;
}

