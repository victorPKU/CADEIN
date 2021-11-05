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



#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "record.h"
#include "combine.h"
#include "oxypair.h"
#include "oxylib.h"
#include "mutation.h"
#include "backrub.h"
#include "protein.h"
#include "rotamer.h"
#include "mem.h"


#define LIBFILELEN  100
int read_pdb_lib(char lib[][100], char *fn)
{
	FILE *libf;
	char line[LIBFILELEN];
	int c=0;
	int i; 

	if((libf=fopen(fn,"r"))==NULL){
		printf("Can not open  PDB lib file %s.\n", fn);
		exit(0);
	}
	while(fgets(line,LIBFILELEN,libf)){
		if(line[0]=='%')
			continue;
		if(c>=100000) break;
		i=0;
		while(i<LIBFILELEN-1&&line[i]!='\n'&&line[i]!='\0'){
			lib[c][i]=line[i]; i++;
		}
		lib[c][i]='\0';
		c++;
	}
	fclose(libf);

	return c;
}
void clean_libs_0(int *res_mutindex[4])
{
	int i;

	for(i=0;i<4;i++){
		BStd_free(res_mutindex[i]);
	}
}
void clean_libs_1(int *res_coorindex,RESBACKRUB *backrubs)
{
	BStd_free(res_coorindex);
	BStd_free(backrubs);
}
void clean_libs_2(int *coorpairindex)
{
	BStd_free(coorpairindex);
}

void clean_libs_3(COORPAIR *coorpairs,int *coorpairsim,int *coorpairsimindex)
{
	BStd_free(coorpairs);
	BStd_free(coorpairsim);
	BStd_free(coorpairsimindex);
}
void clean_libs_4(MUTATION *m[4],COORGEN *coorlib,COORCOMBINE *combine,int *index,int *resindex)
{
	int i;

	for(i=0;i<4;i++){
		BStd_free(m[i]);
	}
	BStd_free(coorlib);
	BStd_free(combine);
	BStd_free(index);
	BStd_free(resindex);
}

void clean_libs_5(MUTATION *m[4],int *res_mutindex[4],COORGEN *coorlib,int *res_coorindex,RESBACKRUB *backrubs)
{
	int i;

	for(i=0;i<4;i++){
		BStd_free(m[i]);
		BStd_free(res_mutindex[i]);
	}
	BStd_free(coorlib);
	BStd_free(res_coorindex);
	BStd_free(backrubs);
}

int main(int argc, char *argv[])
{
	char pdb_lib[5000][100];
	PROTEIN cadein;
	int mutN[4];
	MUTATION *m[4];
	int *res_mutindex[4];
	int coorN;
	COORGEN *coorlib=NULL;
	int *res_coorindex=NULL;
	RESBACKRUB *backrubs=NULL;
	int coorpairN;
	COORPAIR *coorpairs=NULL;
	int *coorpairindex=NULL;
      int *coorpairsim=NULL;
      int *coorpairsimindex=NULL;
	int combineN;
	COORCOMBINE *combine=NULL;
	int *index=NULL;
	int *resindex=NULL;
	int resultN,allresultN=0;
	RECORD rec[MAXRECORD];
	unsigned long long pN=1;
	int pdb_n,ipdb;

	read_rotamer_lib(argv[1]);
	/*pdb_n=read_pdb_lib(pdb_lib, argv[2]);
	printf("Scaffold number %d\n",pdb_n);
	for(ipdb=0;ipdb<pdb_n;ipdb++){
		strcpy(cadein.pdbname, pdb_lib[ipdb]);
		printf("%s\n",cadein.pdbname);
		read_protein(&cadein, argv[3],pdb_lib[ipdb]);*/
		read_protein(&cadein, argv[2]);
		resultN=0;
		get_reslink(&cadein);
		BStd_new(backrubs, cadein.pro.resn*BACKRUBAN);
		backrub_protein(&(cadein.pro), backrubs);
		muta_protein(&(cadein.pro),backrubs,mutN,m, res_mutindex);
		BStd_new(coorlib,cadein.pro.resn+mutN[0]*3+mutN[1]*3+mutN[2]*2+mutN[3]);
		BStd_new(res_coorindex,cadein.pro.resn*5+1);
		coorN=gen_COOR_lib(&(cadein.pro),m, res_mutindex, coorlib, res_coorindex);
		printf("coorlib %d,%d\n",cadein.pro.resn+mutN[0]*3+mutN[1]*3+mutN[2]*2+mutN[3],coorN);
		/*printf("coorlib 26000: %d,%d, %d,%d\n",coorlib[26000].muttype,coorlib[26000].res,coorlib[26000].coorat,coorlib[26000].source);*/
		BStd_renew(coorlib,coorN);
		if(coorN!=0){
			clean_libs_0(res_mutindex);
			BStd_new(coorpairindex,(cadein.pro.resn*5)*(cadein.pro.resn*5-1));
			BStd_new(coorpairs,pN*coorN*100);
			coorpairN=gen_coorpairs(&(cadein.pro),backrubs,m,coorlib,res_coorindex,coorpairs, coorpairindex);
			clean_libs_1(res_coorindex,backrubs);
			printf("coorpairlib %d,%d\n",coorpairN,coorN*100);
			/*printf("coorlib 26000: %d,%d, %d,%d\n",coorlib[26000].muttype,coorlib[26000].res,coorlib[26000].coorat,coorlib[26000].source);*/
			BStd_renew(coorpairs,coorpairN);
			BStd_new(coorpairsim,coorpairN);
			BStd_new(coorpairsimindex,(cadein.pro.resn*5)*(cadein.pro.resn*5-1)/2+1);
			printf("coorpairsimlib %d,%d\n",coorpairN,(cadein.pro.resn*5)*(cadein.pro.resn*5-1)/2+1);
			simple_coorlink(cadein.pro.resn, coorN, coorlib, &coorpairN,coorpairs,coorpairindex,coorpairsim,coorpairsimindex);
			printf("combinelib %d\n",coorpairN);
			/*printf("coorlib 26000: %d,%d, %d,%d\n",coorlib[26000].muttype,coorlib[26000].res,coorlib[26000].coorat,coorlib[26000].source);*/

			clean_libs_2(coorpairindex);
			BStd_new(combine,coorpairN);
			combineN=combine_coordination(coorlib,coorpairN,coorpairs,coorpairsim,coorpairsimindex,combine,m);
			/*printf("combine 41247: %d %d %d %d\n",combine[41247].coor[0],combine[41247].coor[1],combine[41247].coor[2],combine[41247].coor[3]);*/
			printf("combine end %5d\n",combineN);
			clean_libs_3(coorpairs, coorpairsim, coorpairsimindex);
			if(combineN!=0){
				BStd_new(index,combineN);
				BStd_new(resindex,combineN+1);
				resultN=arrange_combine(combineN, combine,index, resindex);
				/*printf("combine 41247: %d %d %d %d\n",combine[41247].coor[0],combine[41247].coor[1],combine[41247].coor[2],combine[41247].coor[3]);*/
				BStd_renew(resindex,resultN+1);
				allresultN+=resultN;
				record_result(rec,&cadein,coorlib,m, combine,index, resultN,resindex);
				printf("record resultN %5d\n",resultN);
			}
			clean_libs_4(m,coorlib,combine,index,resindex);
		}
		else clean_libs_5(m,res_mutindex,coorlib,res_coorindex,backrubs);
 		clean_protein(&cadein);
		/*if(resultN!=0){
			write_temp_result(rec,argv[4]);
		}	
	}*/
	if(allresultN!=0)
		write_result(rec,argv[3]);	

	return EXIT_SUCCESS;
}
