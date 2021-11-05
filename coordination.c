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

#include "geometry.h"
#include "coordination.h"
#include <stdio.h>
#include <math.h>

#define FAILSCORE -20
#define COBOXYBEST 16
#define COBOXYPASS 8
#define NITROGENBEST 16
#define NITROGENPASS 8
#define SULFURBEST 24
#define SULFURPASS 10
#define OXYBEST 10
#define OXYPASS 4
#define MOXYBEST 16
#define MOXYPASS 6

int Cd2OxyScore(COORGEN *coorlib,float Cd[3])
{
	float Cddis1,Cddis2,Cddis,Cddisdiff;
	float Cddihedral;

	Cddis1=distance(Cd,coorlib->coor);
	Cddis2=distance(Cd,coorlib->thirdatom);
	if(Cddis1>Cddis2){
		Cddis=Cddis1;Cddis1=Cddis2;Cddis2=Cddis;
	}
	Cddisdiff=Cddis2-Cddis1;
	if(Cddis1<2.10||Cddis2>2.70||Cddisdiff>0.4) return FAILSCORE;
	Cddihedral=fabs(TODEG(cal_dih(coorlib->link, coorlib->coor, coorlib->thirdatom, Cd)));

	if(Cddihedral<135) return FAILSCORE;

	if(Cddis1>2.20&&Cddis2<2.55&&Cddisdiff<0.25&&Cddihedral>145) return COBOXYBEST;
	return COBOXYPASS;
}
int CdNScore(COORGEN *coorlib,float Cd[3])
{
	float Cddis1,Cddis2,Cddis,Cddisdiff;
	float Cddihedral;

	Cddis=distance(Cd,coorlib->coor);
	Cddis1=distance(Cd,coorlib->link);
	Cddis2=distance(Cd,coorlib->thirdatom);
	Cddisdiff=fabs(Cddis2-Cddis1);
	if(Cddis<2.05||Cddis>2.60||Cddisdiff>0.5) return FAILSCORE;
	Cddihedral=fabs(TODEG(cal_dih(coorlib->link, coorlib->coor, coorlib->thirdatom, Cd)));

	if(Cddihedral<135) return FAILSCORE;

	if(Cddis1>2.10&&Cddis2<2.45&&Cddisdiff<0.35&&Cddihedral>145) return NITROGENBEST;
	return NITROGENPASS;
}
int CdSScore(COORGEN *coorlib,float Cd[3])
{
	float Cddis;
	float Cdangle;

	Cddis=distance(Cd,coorlib->coor);
	if(Cddis<2.18||Cddis>2.70) return FAILSCORE;
	Cdangle=fabs(TODEG(cal_angle(coorlib->link, coorlib->coor, Cd)));

	if(Cdangle<95&&Cdangle>165) return FAILSCORE;

	if(Cddis>2.25&&Cddis<2.55&&Cdangle>105&&Cdangle<150) return SULFURBEST;
	return SULFURPASS;
}
int Cd1OxyScore(COORGEN *coorlib,float Cd[3])
{
	float Cddis;
	float Cdangle;
	float Cddihedral;

	Cddis=distance(Cd,coorlib->coor);
	if(Cddis<2.10||Cddis>2.70) return FAILSCORE;
	Cdangle=fabs(TODEG(cal_angle(coorlib->link, coorlib->coor, Cd)));
	if(Cdangle<95) return FAILSCORE;
	Cddihedral=fabs(TODEG(cal_dih(coorlib->link, coorlib->coor, coorlib->thirdatom, Cd)));
	if(Cddihedral<120) return FAILSCORE;

	if(Cddis>2.20&&Cddis<2.55&&Cdangle>105&&Cddihedral>140) return OXYBEST;
	return OXYPASS;
}
int CdmOScore(COORGEN *coorlib,float Cd[3])
{
	float Cddis;
	float Cddihedral;

	Cddis=distance(Cd,coorlib->coor);
	if(Cddis<2.10||Cddis>2.70) return FAILSCORE;
	Cddihedral=fabs(TODEG(cal_dih(coorlib->thirdatom, coorlib->link, coorlib->coor, Cd)));
	if(Cddihedral<130) return FAILSCORE;

	if(Cddis>2.20&&Cddis<2.55&&Cddihedral>145) return MOXYBEST;
	return MOXYPASS;
}

#define COORMINDIS1    2.68*2.68
#define COORMINDIS2    2.80*2.80

int CollisionScore0(float coorA[3],float coorB[3])
{
	float dis;

	dis=distance2(coorA,coorB);
	
	if(dis<COORMINDIS1) return -6;
	if(dis<COORMINDIS2) return -2;

	return 0;
}

#define COLLISIONDIS1    2.9*2.9
#define COLLISIONDIS2    3.0*3.0

int CollisionScore1(float A[3],float B[3])
{
	float dis;

	dis=distance2(A,B);
	
	if(dis<COLLISIONDIS1) return -6;
	if(dis<COLLISIONDIS2) return -2;

	return 0;
}

int CollisionScore(COORGEN *coorA,COORGEN *coorB)
{
	int score=0;

	score+=CollisionScore0(coorA->coor,coorB->coor);
	score+=CollisionScore1(coorA->coor,coorB->thirdatom);
	score+=CollisionScore1(coorA->coor,coorB->link);
	score+=CollisionScore1(coorB->coor,coorA->thirdatom);
	score+=CollisionScore1(coorB->coor,coorA->link);
	score+=CollisionScore1(coorA->link,coorB->thirdatom);
	score+=CollisionScore1(coorA->link,coorB->link);
	score+=CollisionScore1(coorB->link,coorA->thirdatom);
	score+=CollisionScore1(coorB->thirdatom,coorA->thirdatom);

	return score;
}

void get_CenterCd(int n, COORGEN *coorlib, int index[6],int typ[6],float Cd[3],char cadno[][MAXSAVEDCAD])
{
	int ioxy,joxy,i,j,k;
	float w[6],wsum=0.0;
	int m[6],ii[6],c,f,N;
	float cadx[MAXSAVEDCAD*MAXSAVEDCAD*MAXSAVEDCAD*MAXSAVEDCAD*MAXSAVEDCAD][6][3];
	float sumdis,minsumdis;
	int mini;

	for(i=0;i<3;i++) Cd[i]=0.0;
	for(ioxy=0;ioxy<n;ioxy++){
		m[ioxy]=1;
		for(k=1;k<MAXSAVEDCAD;k++) if(cadno[ioxy][k]!=cadno[ioxy][0]) m[ioxy]++;
		ii[ioxy]=0;
	}
	c=0;f=n-1;
	while(!(ii[0]>=m[0])){
		for(j=0;j<n;j++){
			if(cadno[j][ii[j]]==0){
				for(i=0;i<3;i++){
					cadx[c][j][i]=coorlib[index[j]].Cd[i];
				}
			}
			else{
				for(i=0;i<3;i++){
					cadx[c][j][i]=coorlib[index[j]].Cds[cadno[j][ii[j]]-1][i];
				}
			}
		}
		c++;
		if(ii[f]<m[f]){ ii[f]++;}
		while(ii[f]>=m[f]&&f>0){
			f--;ii[f]++;
		}
		for(k=f+1;k<n;k++)ii[k]=0;
		f=n-1;
	}
	N=c;
	for(ioxy=0;ioxy<n;ioxy++){
		w[ioxy]=1.0;
		if(typ[ioxy]==DOUBLEOXYGEN||typ[ioxy]==SULFUR) w[ioxy]=2.0;
		else if(typ[ioxy]==NITROGEN||typ[ioxy]==MAINCHN) w[ioxy]=1.5;
		wsum+=w[ioxy];
	}
	if(N==1){
		for(ioxy=0;ioxy<n;ioxy++){
			for(i=0;i<3;i++)
				Cd[i]+=cadx[0][ioxy][i]*w[ioxy];
		}
	}
	else{
		minsumdis=10000.0;
		mini=-1;
		for(c=0;c<N;c++){
			sumdis=0.0;
			for(ioxy=0;ioxy<n;ioxy++){
				for(joxy=0;joxy<n;joxy++){
					sumdis+=distance2(cadx[c][ioxy],cadx[c][joxy]);
				}
			}
			if(sumdis<minsumdis){
				minsumdis=sumdis;
				mini=c;
			}
		}
		for(ioxy=0;ioxy<n;ioxy++){
			for(i=0;i<3;i++)
				Cd[i]+=cadx[mini][ioxy][i]*w[ioxy];
		}
	}
	for(i=0;i<3;i++) Cd[i]/=wsum;
}

int floatroundint(float score)
{
	int intscore;

	if(score>=0.0){
		intscore=(int)score;
		if(score-intscore<0.5)
			return intscore;
		else return intscore+1;
	}	
	else{
		intscore=(int)(-score);
		if(score+intscore>-0.5)
			return -intscore;
		else return -intscore-1;
	}
}

int do_coordination(COORCOMBINE *combine, COORGEN *coorlib,MUTATION *m[4],char cadno[][MAXSAVEDCAD])
{
	int ioxy,joxy;
	int Cdscore=0.0;
	int Collscore=0.0;
	float coor_score0=0.0;
	int coor_score=0;
	float Cd0[3],Cd[3];
	float max_Cd[3];
	int x,y,z;
	int maxscore=0;
	max_Cd[0]=0.0;
	max_Cd[1]=0.0;
	max_Cd[2]=0.0;

	Collscore=0;

	for(ioxy=0;ioxy<combine->coorresN;ioxy++){
		for(joxy=ioxy+1;joxy<combine->coorresN;joxy++){
			Collscore+=CollisionScore(coorlib+combine->coor[ioxy],coorlib+combine->coor[joxy]);
		}
	}

	if(Collscore<-8) return Collscore;

	get_CenterCd(combine->coorresN, coorlib, combine->coor, combine->coortyp,Cd0,cadno);

	for(ioxy=0;ioxy<combine->coorresN;ioxy++){
		if(coorlib[combine->coor[ioxy]].muttype!=-1)
			coor_score0+=m[coorlib[combine->coor[ioxy]].muttype][coorlib[combine->coor[ioxy]].source].mut_score*2.4;
		else
			coor_score0+=8;
	}

	for(x=-4;x<=4;x++){ for(y=-4;y<=4;y++){ for(z=-4;z<=4;z++){
		Cd[0]=Cd0[0]+x*0.1;
		Cd[1]=Cd0[1]+y*0.1;
		Cd[2]=Cd0[2]+z*0.1;
		coor_score=0;
	for(ioxy=0;ioxy<combine->coorresN;ioxy++){
		if(combine->coortyp[ioxy]==MAINCHN)
			Cdscore=CdmOScore(coorlib+combine->coor[ioxy],Cd);
		else if(combine->coortyp[ioxy]==DOUBLEOXYGEN)
			Cdscore=Cd2OxyScore(coorlib+combine->coor[ioxy],Cd);
		else if(combine->coortyp[ioxy]==SINGLEOXYGEN)
			Cdscore=Cd1OxyScore(coorlib+combine->coor[ioxy],Cd);
		else if(combine->coortyp[ioxy]==NITROGEN)
			Cdscore=CdNScore(coorlib+combine->coor[ioxy],Cd);
		else if(combine->coortyp[ioxy]==SULFUR)
			Cdscore=CdSScore(coorlib+combine->coor[ioxy],Cd);
		coor_score+=Cdscore;
	}
	if(maxscore<coor_score){
		maxscore=coor_score;
		max_Cd[0]=Cd[0];
		max_Cd[1]=Cd[1];
		max_Cd[2]=Cd[2];
	}
	}}}
	combine->Cd[0]=max_Cd[0];
	combine->Cd[1]=max_Cd[1];
	combine->Cd[2]=max_Cd[2];

	coor_score=maxscore+floatroundint(coor_score0);

	return coor_score+Collscore;
}

int do_coordination2(COORCOMBINE *combine, COORGEN *coorlib,MUTATION *m[4],char cadno[][MAXSAVEDCAD])
{
	int ioxy,joxy;
	int Cdscore=0;
	int Collscore=0;
	float coor_score0=0.0;
	int coor_score=0;
	float Cd[3];

	Collscore=0;

	for(ioxy=0;ioxy<combine->coorresN;ioxy++){
		for(joxy=ioxy+1;joxy<combine->coorresN;joxy++){
			Collscore+=CollisionScore(coorlib+combine->coor[ioxy],coorlib+combine->coor[joxy]);
		}
	}

	if(Collscore<-8) return Collscore;

	get_CenterCd(combine->coorresN, coorlib, combine->coor, combine->coortyp,Cd,cadno);

	for(ioxy=0;ioxy<combine->coorresN;ioxy++){
		if(coorlib[combine->coor[ioxy]].muttype!=-1)
			coor_score0+=m[coorlib[combine->coor[ioxy]].muttype][coorlib[combine->coor[ioxy]].source].mut_score*2.4;
		else
			coor_score0+=8;
	}

	for(ioxy=0;ioxy<combine->coorresN;ioxy++){
		if(combine->coortyp[ioxy]==MAINCHN)
			Cdscore=CdmOScore(coorlib+combine->coor[ioxy],Cd);
		else if(combine->coortyp[ioxy]==DOUBLEOXYGEN)
			Cdscore=Cd2OxyScore(coorlib+combine->coor[ioxy],Cd);
		else if(combine->coortyp[ioxy]==SINGLEOXYGEN)
			Cdscore=Cd1OxyScore(coorlib+combine->coor[ioxy],Cd);
		else if(combine->coortyp[ioxy]==NITROGEN)
			Cdscore=CdNScore(coorlib+combine->coor[ioxy],Cd);
		else if(combine->coortyp[ioxy]==SULFUR)
			Cdscore=CdSScore(coorlib+combine->coor[ioxy],Cd);
		if(Cdscore<0) return Cdscore;
		coor_score+=Cdscore;
	}

	coor_score=coor_score+floatroundint(coor_score0);

	return coor_score+Collscore;
}
