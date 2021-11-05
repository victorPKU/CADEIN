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

#include "backrub.h"
#include "geometry.h"
#include "mem.h"

#define ROT2LAMBDA  0.7
#define STATAO  111.0
#define TAODEV   5.0
#define TAO0DEV   2.0

int resi_backrub(float CApx[3], float CAnx[3],float mainatx[7][3], float backrub_mainatx[7][3], float theta)
{
	float alpha1,alpha2;
	float tao0,tao;

	rot_axis(7,mainatx,backrub_mainatx,CApx,CAnx,theta);
	alpha1=TODEG(cal_dih(mainatx[1],CApx,backrub_mainatx[3],backrub_mainatx[1]));
	rot_axis(3,backrub_mainatx,backrub_mainatx,CApx,backrub_mainatx[3],alpha1*ROT2LAMBDA);
	alpha2=TODEG(cal_dih(mainatx[6],backrub_mainatx[3],CAnx,backrub_mainatx[6]));
	rot_axis(3,backrub_mainatx+4,backrub_mainatx+4,backrub_mainatx[3],CAnx,alpha2*ROT2LAMBDA);

	tao0=TODEG(cal_angle(mainatx[2],mainatx[3],mainatx[4]));
	tao=TODEG(cal_angle(backrub_mainatx[2],backrub_mainatx[3],backrub_mainatx[4]));

	if(tao>STATAO+TAODEV&&tao>tao0+TAO0DEV)
		return 0;
	else if(tao<STATAO-TAODEV&&tao<tao0-TAO0DEV)
		return 0;
	else
		return 1;
}

void cal_rub_phipsi(float mainatx[7][3], float phipsi[2])
{
	phipsi[0]=TODEG(cal_dih(mainatx[0],mainatx[2],mainatx[3],mainatx[4]));
	phipsi[1]=TODEG(cal_dih(mainatx[2],mainatx[3],mainatx[4],mainatx[6]));
}
#define MAXOXYDIS2   0.3*0.3
#define MAXOXYANGLE2   0.4*0.4
#define MAXNITRODIS2   0.6*0.6

int CO_compatible(float co1[][3], float co2[][3])
{
	float oxydis;
	float angle;
	float v1[3],v2[3];
	int i;

	oxydis=distance2(co1[1],co2[1]);

	if(oxydis>MAXOXYDIS2) return 0;

	for(i=0;i<3;i++){
		v1[i]=co1[1][i]-co1[0][i];
		v2[i]=co2[1][i]-co2[0][i];
	}
	angle=distance2(v1,v2);
	if(angle>MAXOXYANGLE2) return 0;

	return 1;
}
int N_compatible(float n1[3], float n2[3])
{
	float nitrodis;
	
	nitrodis=distance2(n1,n2);

	if(nitrodis>MAXNITRODIS2) return 0;

	return 1;
}
void record_rub(int ir, int rubi, float mainatx[7][3],float mainatx0[7][3],RESBACKRUB *backrub)
{
	int ia,i;

	backrub->at_res=ir;
	backrub->backrub=rubi;
	backrub->exist=1;
	for(ia=0;ia<4;ia++)
		for(i=0;i<3;i++)
			backrub->mx[ia][i]=mainatx[ia+2][i];
	for(i=0;i<3;i++){
		backrub->COx[0][i]=mainatx[0][i];
		backrub->COx[1][i]=mainatx[1][i];
		backrub->Nx[i]=mainatx[6][i];
	}
	cal_rub_phipsi(mainatx, backrub->phipsi);

	backrub->native_compatible[0]=CO_compatible(mainatx+4, mainatx0+4);
	backrub->native_compatible[1]=CO_compatible(mainatx, mainatx0);
}
void record_noexit(int ir,int rubi,RESBACKRUB *backrub)
{
	backrub->at_res=ir;
	backrub->backrub=rubi;
	backrub->exist=0;
}

void cal_uncompatible(RESBACKRUB *backrubs, int ir,RESBACKRUB *backrub0)
{
	int rubi;

	backrub0->pre_uncompatible[0]=0;
	for(rubi=0;rubi<BACKRUBAN;rubi++){
		if(backrubs[(ir-1)*BACKRUBAN+rubi].exist==1){
			if(CO_compatible(backrub0->COx, backrubs[(ir-1)*BACKRUBAN+rubi].mx+2)==0){
				backrub0->pre_uncompatible[0]++;
				backrub0->pre_uncompatible[backrub0->pre_uncompatible[0]]=(ir-1)*BACKRUBAN+rubi;
			}
			else if(N_compatible(backrub0->mx[0], backrubs[(ir-1)*BACKRUBAN+rubi].Nx)==0){
				backrub0->pre_uncompatible[0]++;
				backrub0->pre_uncompatible[backrub0->pre_uncompatible[0]]=(ir-1)*BACKRUBAN+rubi;
			}
		}
	}
}
void resi_cpy2rub(int ir, int rubi, float mainatx[7][3],RESBACKRUB *backrub)
{
	int ia,i;

	backrub->at_res=ir;
	backrub->backrub=rubi;
	backrub->exist=1;
	for(ia=0;ia<4;ia++)
		for(i=0;i<3;i++)
			backrub->mx[ia][i]=mainatx[ia+2][i];

	backrub->native_compatible[0]=1;
	backrub->native_compatible[1]=1;
}

void backrub_protein(PROCHAIN *p, RESBACKRUB *backrubs)
{
	int ir;
	int canbr;
	float CApx[3],CAnx[3];
	float mainatx[7][3];
	float backrub_mainatx[7][3];
	int rubi;

	for(ir=0;ir<p->resn;ir++){
		if(p->r[ir].mutability==0){
			for(rubi=0;rubi<BACKRUBAN;rubi++){
				record_noexit(ir,rubi,backrubs+ir*BACKRUBAN+rubi);
			}
		} 
		else if(p->r[ir].end-p->r[ir].start>=3){
			cpx(p->a[p->r[ir].start+0].x,mainatx[2]);
			cpx(p->a[p->r[ir].start+1].x,mainatx[3]);
			cpx(p->a[p->r[ir].start+2].x,mainatx[4]);
			cpx(p->a[p->r[ir].start+3].x,mainatx[5]);
			if(p->r[ir].pre!=-1&&p->r[ir].next!=-1){
				cpx(p->a[p->r[p->r[ir].pre].start+1].x,CApx);	
				cpx(p->a[p->r[p->r[ir].next].start+1].x,CAnx);
				cpx(p->a[p->r[p->r[ir].pre].start+2].x,mainatx[0]);
				cpx(p->a[p->r[p->r[ir].pre].start+3].x,mainatx[1]);
				cpx(p->a[p->r[p->r[ir].next].start].x,mainatx[6]);
			
				for(rubi=0;rubi<BACKRUBAN;rubi++){
					canbr=resi_backrub(CApx,CAnx,mainatx,backrub_mainatx,rub_angle_plan[rubi]);
					if(canbr==1){
						record_rub(ir,rubi,backrub_mainatx,mainatx,backrubs+ir*BACKRUBAN+rubi);
						cal_uncompatible(backrubs,ir,backrubs+ir*BACKRUBAN+rubi);
					}
					else{
						record_noexit(ir,rubi,backrubs+ir*BACKRUBAN+rubi);
					}
				}
			}
			else {
				for(rubi=0;rubi<BACKRUBAN;rubi++){
					if(rubi==BACKRUBAN/2){
						resi_cpy2rub(ir,rubi,mainatx,backrubs+ir*BACKRUBAN+rubi);
						if(p->r[ir].pre!=-1){
							cpx(p->a[p->r[p->r[ir].pre].start+2].x,mainatx[0]);
							backrubs[ir*BACKRUBAN+rubi].phipsi[0]=TODEG(cal_dih(mainatx[0],mainatx[2],mainatx[3],mainatx[4]));
						}
						else backrubs[ir*BACKRUBAN+rubi].phipsi[0]=-120.0;
						if(p->r[ir].next!=-1){
							cpx(p->a[p->r[p->r[ir].next].start].x,mainatx[6]);
							backrubs[ir*BACKRUBAN+rubi].phipsi[1]=TODEG(cal_dih(mainatx[2],mainatx[3],mainatx[4],mainatx[6]));
						}
						else backrubs[ir*BACKRUBAN+rubi].phipsi[1]=+120.0;
					}
					else{
						record_noexit(ir,rubi,backrubs+ir*BACKRUBAN+rubi);
					}
				}
			}
		}
		else{
			for(rubi=0;rubi<BACKRUBAN;rubi++)
				record_noexit(ir,rubi,backrubs+ir*BACKRUBAN+rubi);
		}
	}
}

