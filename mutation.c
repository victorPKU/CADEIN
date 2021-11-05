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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mutation.h"
#include "mem.h"
#include "geometry.h"
#include "rotamer.h"
#include "superpose.h"

int cal_res_collision(PROCHAIN *p, int ir,int side_atn,float sideatx[][3])
{
	int count=0;
	int ia,ja;
	float dis2;
	float nearest;

	for(ia=0;ia<side_atn-3;ia++){
		for(ja=0;ja<p->atmn;ja++){
			if(p->a[ja].at_res==ir)	
				continue;
			dis2=distance2(sideatx[ia],p->a[ja].x);
			nearest=(p->a[ja].vdw_rad+1.9)*COLLISIONDIS;
			if(dis2<nearest*nearest*0.9*0.9*0.9*0.9) count+=4;
			else if(dis2<nearest*nearest*0.9*0.9) count+=2;
			else if(dis2<nearest*nearest)
				count++;
		}
	}
	for(;ia<side_atn-2;ia++){
		for(ja=0;ja<p->atmn;ja++){
			if(p->a[ja].at_res==ir)	
				continue;
			dis2=distance2(sideatx[ia],p->a[ja].x);
			nearest=(p->a[ja].vdw_rad+1.8)*COLLISIONDIS;
			if(dis2<nearest*nearest*0.9*0.9*0.9*0.9) count+=4;
			else if(dis2<nearest*nearest*0.9*0.9) count+=2;
			else if(dis2<nearest*nearest)
				count++;
		}
	}
	for(;ia<side_atn;ia++){
		for(ja=0;ja<p->atmn;ja++){
			if(p->a[ja].at_res==ir)	
				continue;
			dis2=distance2(sideatx[ia],p->a[ja].x);
			nearest=(p->a[ja].vdw_rad+1.65)*COLLISIONDIS;
			if(dis2<nearest*nearest*0.9*0.9*0.9*0.9) count+=4;
			else if(dis2<nearest*nearest*0.9*0.9) count+=2;
			else if(dis2<nearest*nearest)
				count++;
		}
	}
	return count;
}
void graft(float main0[][3], float maint[][3],int siden, float side0[][3], float side[][3])
{
	float R[3][3];
	float main0c[3][3];
	float maintc[3][3];
	float center1[3];
	float center2[3];

	get_center(3,main0,center1);
	get_center(3,maint,center2);
	irtranslate(3,main0,main0c,center1);
	irtranslate(3,maint,maintc,center2);
	calc_fit_R(3,main0c,maintc,R);
	irtranslate(siden,side0,side,center2);
	do_rot(siden,side,R);
	translate(siden,side,side,center1);
}

void get_rotamer_x(float chi[4], float CAx[3], int siden, float side[][3], int muttype)
{
	rot_axis(siden-1,side+1, side+1, CAx, side[0],EDHC_chi[muttype][0]-chi[0]);
	if(muttype!=MUT2CYS){
		rot_axis(siden-2,side+2, side+2, side[0], side[1],EDHC_chi[muttype][1]-chi[1]);
	}
	if(muttype==MUT2GLU){
		rot_axis(siden-3,side+3, side+3, side[1], side[2],EDHC_chi[muttype][2]-chi[2]);
	}
}

void side_assemble(float mainatx[][3], float sideatx[][3],int muttype,bbdepROT *rota,int shift[4])
{
	float chis[4];
	int i;

	for(i=0;i<4;i++) chis[i]=rota->chi[i]+rota->sig[i]* mut_shift_plan[shift[i]];
	graft(mainatx, main_chain_template, EDHC_sidean[muttype], EDHC_template[muttype], sideatx);
	get_rotamer_x(chis, mainatx[1], EDHC_sidean[muttype], sideatx, muttype);
}
#define MAXCOLLIDE 5.0
int mutate(PROCHAIN *p, int ir, MUTATION *m, float mx[][3], float phipsi[2],int muttype)
{
	int rotai;
	int rota_start,rota_end;
	int mutn=0;
	int collide;
	int shiftn;
	int shifti;
	int shift[4];

	get_phipsi_index(rotamer_lib+muttype,phipsi[0],phipsi[1],&rota_start,&rota_end);
	shiftn=pow(MSHIFTN,EDHC_chiN[muttype]);
	for(rotai=rota_start;rotai<=rota_end;rotai++){
		for(shifti=0;shifti<shiftn;shifti++){
			shift[3]=shifti%MSHIFTN;
			shift[2]=((shifti-shift[3])/MSHIFTN)%MSHIFTN;
			shift[1]=((shifti-shift[3]-shift[2]*MSHIFTN)/(MSHIFTN*MSHIFTN))%MSHIFTN;
			shift[0]=0;
			side_assemble(mx,m->sx, muttype,rotamer_lib[muttype].rotas+rotai,shift);
			collide=cal_res_collision(p,ir,EDHC_sidean[muttype],m->sx);
			/*if(collide<0)
				printf("%8d\n",collide);*/
			if(collide<=MAXCOLLIDE){
				m->rota=rotai;
				m->rotshift=shifti;
				m->satn=EDHC_sidean[muttype];
				m->mut_score=collide;
				m++;
				mutn++;
			}
		}
	}
	return mutn;
}

void get_native(PROCHAIN *p, MUTATION *m, int ir)
{
	int ia,i;
	int start;

	m->at_res=ir;
	m->backrub=ir*BACKRUBAN+BACKRUBAN/2;
	m->rota=-1;
	if(!strncmp(p->r[ir].name,"GLU",3))
		m->muttype=MUT2GLU;
	else if(!strncmp(p->r[ir].name,"ASP",3))
		m->muttype=MUT2ASP;
	else if(!strncmp(p->r[ir].name,"HIS",3))
		m->muttype=MUT2HIS;
	else if(!strncmp(p->r[ir].name,"CYS",3))
		m->muttype=MUT2CYS;
	if(!strncmp(p->r[ir].name,"ASP",3))
		m->satn=4;
	else if(!strncmp(p->r[ir].name,"GLU",3))
		m->satn=5;
	else if(!strncmp(p->r[ir].name,"HIS",3))
		m->satn=6;
	else if(!strncmp(p->r[ir].name,"CYS",3))
		m->satn=2;
	start=p->r[ir].sstart;
	for(ia=0;ia<m->satn;ia++){
		for(i=0;i<3;i++)
			m->sx[ia][i]=p->a[start+ia].x[i];
	}
	m->isnative=2;
	m->mut_score=3.0;
}

float Score_mutation(int muttype,int rotshift[4], float collide, int is_native,float rot_possibility, int rubi)
{
	float rubscore=0.0;
	float collidescore=0.0;
	float mutscore=0.0;
	int i,shifti;
	float shiftscore=0.0;

	if(rubi==BACKRUBAN/2) 
		rubscore=0.0;
	else if(rubi>BACKRUBAN/2) 
		rubscore=-(rubi-BACKRUBAN/2)*2.0;
	else if(rubi<BACKRUBAN/2) 
		rubscore=-(BACKRUBAN/2-rubi)*2.0;
	/*else rubscore=-4.0;*/
	for(i=0;i<4;i++){
		shifti=rotshift[i];
		if(shifti==MSHIFTN/2)
			shiftscore+=0.0;
		else if(shifti>MSHIFTN/2)
			shiftscore+=-(shifti-MSHIFTN/2)*0.3;
		else if(shifti<MSHIFTN/2)
			shiftscore+=-(MSHIFTN/2-shifti)*0.3;
	}

	collidescore=collide*(-1.5);
	if(muttype==MUT2ASP||muttype==MUT2GLU) mutscore=0.8;
	else if(muttype==MUT2HIS) mutscore=1.0;
	else if(muttype==MUT2CYS) mutscore=1.5;

	if(mutscore+rot_possibility*0.9+collidescore+rubscore+is_native*2+shiftscore>10)
		printf("%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n",mutscore*1.0,rot_possibility*0.9,collidescore*1.0,rubscore*1.0,is_native*2.0,shiftscore*1.0);

	return mutscore+rot_possibility*0.9+collidescore+rubscore+is_native*2+shiftscore;
}

#define MINMUTSCORE -6.0
int cp_mut_message(PROCHAIN *p, MUTATION *mi, MUTATION *mo,int ir,int muttype,int rubi)
{
	float collide;
	int shifti;
	int rotshift[4];
	int i,j;
	int isnative=0;
	float mutscore;
	
	if(!strncmp(p->r[ir].name,"GLU",3)){
		if(muttype==MUT2GLU)
			isnative=1;
		else
			isnative=0;
	}
	else if(!strncmp(p->r[ir].name,"ASP",3)){
		if(muttype==MUT2ASP)
			isnative=1;
		else
			isnative=0;
	}
	else if(!strncmp(p->r[ir].name,"HIS",3)){
		if(muttype==MUT2HIS)
			isnative=1;
		else
			isnative=0;
	}
	else if(!strncmp(p->r[ir].name,"CYS",3)){
		if(muttype==MUT2CYS)
			isnative=1;
		else
			isnative=0;
	}
	else
		isnative=0;

	collide=mi->mut_score;
	if(collide<0)
		printf("%8.3f\n",collide);
	shifti=mi->rotshift;
	rotshift[3]=shifti%MSHIFTN;
	rotshift[2]=((shifti-rotshift[3])/MSHIFTN)%MSHIFTN;
	rotshift[1]=((shifti-rotshift[3]-rotshift[2]*MSHIFTN)/(MSHIFTN*MSHIFTN))%MSHIFTN;
	rotshift[0]=0;
	mutscore=Score_mutation(muttype,rotshift,collide, isnative,rotamer_lib[muttype].rotas[mi->rota].p, rubi);

	if(mutscore<MINMUTSCORE+isnative)	return 0;

	mo->at_res=ir;
	mo->backrub=ir*BACKRUBAN+rubi;
	mo->muttype=muttype;
	mo->rota=mi->rota;
	mo->rotshift=mi->rotshift;
	mo->isnative=isnative;
	mo->mut_score=mutscore;
	mo->satn=mi->satn;
	for(i=0;i<mi->satn;i++){
		for(j=0;j<3;j++){
			mo->sx[i][j]=mi->sx[i][j];
		}
	}

	return 1;

}
void res_mut(PROCHAIN *p,RESBACKRUB *backrubs,int mutN[4],MUTATION *m[4], int ir)
{
	int typei;
	int rubi;
	int mn,imut;
	int n0;

	for(rubi=0;rubi<BACKRUBAN;rubi++){
		if(backrubs[rubi].exist==1){
			if(rubi==BACKRUBAN/2){
				if(!strncmp(p->r[ir].name,"GLU",3)&&(p->r[ir].send-p->r[ir].sstart==5)){
					get_native(p,m[0]+mutN[0],ir);
					mutN[0]++;
				}
				if(!strncmp(p->r[ir].name,"ASP",3)&&(p->r[ir].send-p->r[ir].sstart==4)){
					get_native(p,m[1]+mutN[1],ir);
					mutN[1]++;
				}
				if(!strncmp(p->r[ir].name,"HIS",3)&&(p->r[ir].send-p->r[ir].sstart==6)){
					get_native(p,m[2]+mutN[2],ir);
					if(ir==45||ir==116)
						printf("%3d %7d %3d\n",ir,mutN[2],m[2][mutN[2]].rota);
					mutN[2]++;
				}
				if(!strncmp(p->r[ir].name,"CYS",3)&&(p->r[ir].send-p->r[ir].sstart==2)){
					get_native(p,m[3]+mutN[3],ir);
					if(ir==111)
						printf("%3d %7d %3d\n",ir,mutN[3],m[3][mutN[3]].rota);
					mutN[3]++;
				}
			}
			for(typei=0;typei<4;typei++){
				mn=mutate(p, ir,m[typei]+mutN[typei], backrubs[rubi].mx, backrubs[rubi].phipsi,typei);
				n0=mutN[typei];
				for(imut=0;imut<mn;imut++){
					mutN[typei]+=cp_mut_message(p,m[typei]+n0+imut,m[typei]+mutN[typei],ir,typei,rubi);
				}
			}
		}
	}
}
void muta_protein(PROCHAIN *p,RESBACKRUB *backrubs,int mutN[4],MUTATION *mutation[4], int *res_mutindex[4])
{
	int ir,im;
	int shiftn;

	get_maxrotnum();
	for(im=0;im<4;im++){
		shiftn=pow(MSHIFTN,EDHC_chiN[im]);
		BStd_new(mutation[im], p->resn*MAXROTNUM[im]*BACKRUBAN*shiftn );
		BStd_new(res_mutindex[im], p->resn+1);
		mutN[im]=0;
		res_mutindex[im][0]=0;
		/*printf("%10d",p->resn*MAXROTNUM[0]*BACKRUBAN*shiftn);*/
	}
	/*printf("\n");*/
	for(ir=0;ir<p->resn;ir++){
		if(p->r[ir].mutability==1) res_mut(p, backrubs+ir*BACKRUBAN,mutN,mutation,ir);
		for(im=0;im<4;im++)
			res_mutindex[im][ir+1]=mutN[im];
		/*printf("%3d %6d %6d %6d %6d\n",ir+1,mutN[0],mutN[1],mutN[2],mutN[3]);*/
	}
	for(im=0;im<4;im++)
		BStd_renew(mutation[im], mutN[im]);
}
