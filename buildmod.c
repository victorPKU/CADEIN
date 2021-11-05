/***************************************************************************
 *                             CADEIN
 *  Program for cadmium binding protein design   *
 *         Author: Dr. Zhang, Changsheng 
 * Copyright (c) 2021 Molecular Design Lab at Peking University  *
 *
 *   contact: changshengzhang@pku.edu.cn
 *
 *  Free for academical purpose and please reference to:
 *  Tailoring Escherichia coli Chemotactic Sensing towards Cadmium 
 *  by Computational Redesign of Ribose-binding Protein. Hengyi Li, 
 *  Changsheng Zhang, Xi Chen and Luhua Lai. mSystems. 2021. *
 ***************************************************************************/


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "record.h"
#include "mutation.h"
#include "backrub.h"
#include "protein.h"
#include "model.h"
#include "geometry.h"
#include "rotamer.h"
#include "mem.h"

int get_resindex(PROCHAIN *pro, char *name)
{
	int ir;

	for(ir=0;ir<pro->resn;ir++){
		if(!strcmp(pro->r[ir].name,name))
			return ir;
	}

	return -1;
}

void get_native_res(PROCHAIN *pro,int res,MODEL *mod)
{
	int ia;
	int i;

	for(ia=pro->r[res].start;ia<=pro->r[res].end;ia++){
		strncpy(mod->at[mod->atn].atname,pro->a[ia].name,7);
		mod->at[mod->atn].atname[7]='\0';
		strncpy(mod->at[mod->atn].resnum,pro->a[ia].name+8,5);
		mod->at[mod->atn].resnum[5]='\0';
		for(i=0;i<3;i++)
			mod->at[mod->atn].x[i]=pro->a[ia].x[i];
		mod->atn++;
	}
}
void cpy2model(float atmx[4][3],float sideatx[5][3], int muttype,char *resnum, MODEL *mod)
{
	int ia;
	int i;

	for(ia=0;ia<4;ia++){
		for(i=0;i<3;i++)
			mod->at[mod->atn].x[i]=atmx[ia][i];
		if(ia==0){
			mod->at[mod->atn].atname[0]='N';
			mod->at[mod->atn].atname[1]=' ';
			mod->at[mod->atn].atname[2]=' ';
		}
		else if(ia==1){
			mod->at[mod->atn].atname[0]='C';
			mod->at[mod->atn].atname[1]='A';
			mod->at[mod->atn].atname[2]=' ';
		}
		else if(ia==2){
			mod->at[mod->atn].atname[0]='C';
			mod->at[mod->atn].atname[1]=' ';
			mod->at[mod->atn].atname[2]=' ';
		}
		else if(ia==3){
			mod->at[mod->atn].atname[0]='O';
			mod->at[mod->atn].atname[1]=' ';
			mod->at[mod->atn].atname[2]=' ';
		}
		mod->at[mod->atn].atname[3]=' ';
		if(muttype==MUT2GLU){
			mod->at[mod->atn].atname[4]='G';
			mod->at[mod->atn].atname[5]='L';
			mod->at[mod->atn].atname[6]='U';
		}
		else if(muttype==MUT2ASP){
			mod->at[mod->atn].atname[4]='A';
			mod->at[mod->atn].atname[5]='S';
			mod->at[mod->atn].atname[6]='P';
		}
		else if(muttype==MUT2HIS){
			mod->at[mod->atn].atname[4]='H';
			mod->at[mod->atn].atname[5]='I';
			mod->at[mod->atn].atname[6]='S';
		}
		else if(muttype==MUT2CYS){
			mod->at[mod->atn].atname[4]='C';
			mod->at[mod->atn].atname[5]='Y';
			mod->at[mod->atn].atname[6]='S';
		}
		mod->at[mod->atn].atname[7]='\0';
		strcpy(mod->at[mod->atn].resnum,resnum);
		mod->at[mod->atn].resnum[5]='\0';
		mod->atn++;
	}

	for(ia=0;ia<EDHC_sidean[muttype];ia++){
		for(i=0;i<3;i++)
			mod->at[mod->atn].x[i]=sideatx[ia][i];
		mod->at[mod->atn].atname[3]=' ';
		if(ia==0){
			mod->at[mod->atn].atname[0]='C';
			mod->at[mod->atn].atname[1]='B';
			mod->at[mod->atn].atname[2]=' ';
		}
		else if(ia==1&&muttype!=MUT2CYS){
			mod->at[mod->atn].atname[0]='C';
			mod->at[mod->atn].atname[1]='G';
			mod->at[mod->atn].atname[2]=' ';
		}
		if(muttype==MUT2GLU){
			if(ia==2){
				mod->at[mod->atn].atname[0]='C';
				mod->at[mod->atn].atname[1]='D';
				mod->at[mod->atn].atname[2]=' ';
			}
			else if(ia==3){
				mod->at[mod->atn].atname[0]='O';
				mod->at[mod->atn].atname[1]='E';
				mod->at[mod->atn].atname[2]='1';
			}
			else if(ia==4){
				mod->at[mod->atn].atname[0]='O';
				mod->at[mod->atn].atname[1]='E';
				mod->at[mod->atn].atname[2]='2';
			}
			mod->at[mod->atn].atname[4]='G';
			mod->at[mod->atn].atname[5]='L';
			mod->at[mod->atn].atname[6]='U';
		}
		else if(muttype==MUT2ASP){
			if(ia==2){
				mod->at[mod->atn].atname[0]='O';
				mod->at[mod->atn].atname[1]='D';
				mod->at[mod->atn].atname[2]='1';
			}
			else if(ia==3){
				mod->at[mod->atn].atname[0]='O';
				mod->at[mod->atn].atname[1]='D';
				mod->at[mod->atn].atname[2]='2';
			}
			mod->at[mod->atn].atname[4]='A';
			mod->at[mod->atn].atname[5]='S';
			mod->at[mod->atn].atname[6]='P';
		}
		else if(muttype==MUT2HIS){
			if(ia==2){
				mod->at[mod->atn].atname[0]='C';
				mod->at[mod->atn].atname[1]='D';
				mod->at[mod->atn].atname[2]='2';
			}
			else if(ia==3){
				mod->at[mod->atn].atname[0]='N';
				mod->at[mod->atn].atname[1]='D';
				mod->at[mod->atn].atname[2]='1';
			}
			else if(ia==4){
				mod->at[mod->atn].atname[0]='C';
				mod->at[mod->atn].atname[1]='E';
				mod->at[mod->atn].atname[2]='1';
			}
			else if(ia==5){
				mod->at[mod->atn].atname[0]='N';
				mod->at[mod->atn].atname[1]='E';
				mod->at[mod->atn].atname[2]='2';
			}
			mod->at[mod->atn].atname[4]='H';
			mod->at[mod->atn].atname[5]='I';
			mod->at[mod->atn].atname[6]='S';
		}
		else if(muttype==MUT2CYS){
			if(ia==1){
				mod->at[mod->atn].atname[0]='S';
				mod->at[mod->atn].atname[1]='G';
				mod->at[mod->atn].atname[2]=' ';
			}
			mod->at[mod->atn].atname[4]='C';
			mod->at[mod->atn].atname[5]='Y';
			mod->at[mod->atn].atname[6]='S';
		}
		mod->at[mod->atn].atname[7]='\0';
		strcpy(mod->at[mod->atn].resnum,resnum);
		mod->atn++;
	}
}
void build_model(RECORD *rec,PROCHAIN *pro,MODEL *mod)
{
	int ia,i,j;
	int res;
	int rubi;
	int muttype;
	float CApx[3],CAnx[3];
	float mainatx[7][3];
	float backrub_mainatx[7][3];
	float atmx[4][3];
	float sideatx[5][3];
	char resnum[6];

	mod->atn=0;
	for(i=0;i<rec->res_n;i++){
		res=get_resindex(pro,rec->res_name[i]);
		if(rec->res_type[i]==-1||rec->rotamer[i]==-1){
			get_native_res(pro,res,mod);
			continue;
		}
		rubi=rec->res_rub[i];
		if(rubi!=0){
			cpx(pro->a[pro->r[pro->r[res].pre].start+1].x,CApx);
			cpx(pro->a[pro->r[pro->r[res].next].start+1].x,CAnx);
			cpx(pro->a[pro->r[pro->r[res].pre].start+2].x,mainatx[0]);
			cpx(pro->a[pro->r[pro->r[res].pre].start+3].x,mainatx[1]);
			cpx(pro->a[pro->r[res].start+0].x,mainatx[2]);
			cpx(pro->a[pro->r[res].start+1].x,mainatx[3]);
			cpx(pro->a[pro->r[res].start+2].x,mainatx[4]);
			cpx(pro->a[pro->r[res].start+3].x,mainatx[5]);
			cpx(pro->a[pro->r[pro->r[res].next].start].x,mainatx[6]);
			resi_backrub(CApx,CAnx,mainatx,backrub_mainatx,rub_angle_plan[rubi]);
			cpx(backrub_mainatx[2],atmx[0]);
			cpx(backrub_mainatx[3],atmx[1]);
			cpx(backrub_mainatx[4],atmx[2]);
			cpx(backrub_mainatx[5],atmx[3]);
		}
		else{
			for(ia=pro->r[res].start,j=0;ia<pro->r[res].sstart;ia++,j++){
				cpx(pro->a[ia].x,atmx[j]);
			}
		}
		muttype=rec->res_type[i];
		graft(atmx, main_chain_template, EDHC_sidean[muttype], EDHC_template[muttype], sideatx);
		get_rotamer_x(rotamer_lib[muttype].rotas[rec->rotamer[i]].chi, atmx[1], EDHC_sidean[muttype], sideatx, muttype);
		strncpy(resnum,pro->r[res].name+4,5);
		resnum[5]='\0';
		cpy2model(atmx,sideatx, muttype,resnum,mod);
	}
	for(i=0;i<3;i++)
		mod->cadmium[i]=rec->cadmium[i];
}
void write_model(MODEL *mod,char *mfn)
{
	FILE *rf;
	int ia;

	if((rf=fopen(mfn,"w"))==NULL)
		printf("Can't open %s",mfn);

	for(ia=0;ia<mod->atn;ia++){
		fprintf(rf,"ATOM%7d  %s %s    %8.3f%8.3f%8.3f\n",ia+1,mod->at[ia].atname,mod->at[ia].resnum,mod->at[ia].x[0],mod->at[ia].x[1],mod->at[ia].x[2]);
	}
	fprintf(rf,"HETATM%5d CD    CD     1    %8.3f%8.3f%8.3f \n",ia+1,mod->cadmium[0],mod->cadmium[1],mod->cadmium[2]);

	fclose(rf);
}
void get_mfn(int num,char *pdbname,char *dir, char *mfn)
{
	int i;
	int m,n;

	strcpy(mfn,dir);
	m=strlen(dir);
	mfn[m]='/';
	n=strlen(pdbname);
	for(i=0;i<n;i++){
		if(pdbname[i]=='.') break;
		mfn[m+1+i]=pdbname[i];
	}
	n=m+1+i;
	mfn[n]='0'+(num)/1000;
	mfn[n+1]='0'+(num)%1000/100;
	mfn[n+2]='0'+(num)%100/10;
	mfn[n+3]='0'+(num)%10;
	mfn[n+4]='.';
	mfn[n+5]='p';
	mfn[n+6]='d';
	mfn[n+7]='b';
	mfn[n+8]='\0';
}

int main(int argc, char *argv[])
{
	int resultN;
	RECORD rec[MAXRECORD];
	int irec;
	MODEL mod;
	char mfn[100];
	PROTEIN cadein;

	/*printf("read_rotamer_lib\n");*/
	read_rotamer_lib(argv[1]);
	/*printf("read_result\n");*/
	resultN=read_result(rec,argv[2]);
	for(irec=0;irec<resultN;irec++){
		printf("read_protein %s\n",rec[irec].pdbname);
		read_protein(&cadein, rec[irec].pdbname);
		get_reslink(&cadein);
		printf("build_model %5d\n",irec+1);
		build_model(rec+irec,&(cadein.pro), &mod);
		get_mfn(irec+1,rec[irec].pdbname,argv[3],mfn);
		/*printf("  write_model%s\n",mfn);*/
		write_model(&mod,mfn);
		/*printf("clean_protein\n");*/
		clean_protein(&cadein);
	}

  return EXIT_SUCCESS;
}

