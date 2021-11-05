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

#include "oxypair.h"
#include "geometry.h"
#include "mem.h"
#include <stdio.h>

#define MINCDDIS2   1.6*1.6
#define MINCOORDIS2   2.6*2.6
#define MINLINKDIS2  2.8*2.8

float assigncadno(int *c,char cadno[MAXSAVEDCAD][2],char a,char b,float Cdd[MAXSAVEDCAD],float Cddis,float Cddismx)
{
	int i,n;
	int maxi;

	n=(*c);
	if(n==0){
		cadno[0][0]=a;
		cadno[0][1]=b;
		Cdd[0]=Cddis;
		(*c)++;
		return Cddis;
	}
	if(n<MAXSAVEDCAD){
		cadno[n][0]=a;
		cadno[n][1]=b;
		Cdd[n]=Cddis;
		(*c)++;
		if(Cddis>Cddismx)
			return Cddis;
		else return Cddismx;
	}
	if(Cddis>Cddismx){
		return Cddismx;
	}
	maxi=0;
	for(i=1;i<n;i++){
		if(Cdd[i]>Cdd[maxi]){
			maxi=i;
		}
	}
	cadno[maxi][0]=a;
	cadno[maxi][1]=b;
	Cdd[maxi]=Cddis;
	maxi=0;
	for(i=1;i<n;i++){
		if(Cdd[i]>Cdd[maxi]){
			maxi=i;
		}
	}
	return Cdd[maxi];
}
int Is_pair(COORGEN *coor1,COORGEN *coor2,char cadno[MAXSAVEDCAD][2])
{
	float Cddis,Cdd[MAXSAVEDCAD],Cddismx=0.0;
	float coor_dis;
	float link_dis;
	float min_dis;
	int ia,ja,ic;
	int c=0;
	
	for(ic=0;ic<MAXSAVEDCAD;ic++){
		cadno[ic][0]=0;
		cadno[ic][1]=0;
	}

	Cddis=distance2(coor1->Cd,coor2->Cd);
	if(Cddis<MINCDDIS2)	Cddismx=assigncadno(&c,cadno,0,0,Cdd,Cddis,MINCDDIS2);
	if(coor1->type==MAINCHN){
		Cddis=distance2(coor1->Cds[0],coor2->Cd);
		if(Cddis<MINCDDIS2)	Cddismx=assigncadno(&c,cadno,1,0,Cdd,Cddis,Cddismx);
		Cddis=distance2(coor1->Cds[1],coor2->Cd);
		if(Cddis<MINCDDIS2)	Cddismx=assigncadno(&c,cadno,2,0,Cdd,Cddis,Cddismx);
	}
	if(coor2->type==MAINCHN){
		Cddis=distance2(coor1->Cd,coor2->Cds[0]);
		if(Cddis<MINCDDIS2)	Cddismx=assigncadno(&c,cadno,0,1,Cdd,Cddis,Cddismx);
		Cddis=distance2(coor1->Cd,coor2->Cds[1]);
		if(Cddis<MINCDDIS2)	Cddismx=assigncadno(&c,cadno,0,2,Cdd,Cddis,Cddismx);
	}
	if(coor2->type==MAINCHN&&coor1->type==MAINCHN){
		Cddis=distance2(coor1->Cds[0],coor2->Cds[0]);
		if(Cddis<MINCDDIS2)	Cddismx=assigncadno(&c,cadno,1,1,Cdd,Cddis,Cddismx);
		Cddis=distance2(coor1->Cds[0],coor2->Cds[1]);
		if(Cddis<MINCDDIS2)	Cddismx=assigncadno(&c,cadno,1,2,Cdd,Cddis,Cddismx);
		Cddis=distance2(coor1->Cds[1],coor2->Cds[0]);
		if(Cddis<MINCDDIS2)	Cddismx=assigncadno(&c,cadno,2,1,Cdd,Cddis,Cddismx);
		Cddis=distance2(coor1->Cds[1],coor2->Cds[1]);
		if(Cddis<MINCDDIS2)	Cddismx=assigncadno(&c,cadno,2,2,Cdd,Cddis,Cddismx);
	}

	if(coor1->type==SULFUR){
		for(ia=0;ia<6;ia++){
			Cddis=distance2(coor1->Cds[ia],coor2->Cd);
			if(Cddis<MINCDDIS2)	Cddismx=assigncadno(&c,cadno,ia+1,0,Cdd,Cddis,Cddismx);
		}
	}
	if(coor2->type==SULFUR){
		for(ia=0;ia<6;ia++){
			Cddis=distance2(coor1->Cd,coor2->Cds[ia]);
			if(Cddis<MINCDDIS2)	Cddismx=assigncadno(&c,cadno,0,ia+1,Cdd,Cddis,Cddismx);
		}
	}
	if(coor2->type==SULFUR&&coor1->type==SULFUR){
		for(ia=0;ia<6;ia++){
			for(ja=0;ja<6;ja++){
				Cddis=distance2(coor1->Cds[ia],coor2->Cds[ja]);
				if(Cddis<MINCDDIS2)	Cddismx=assigncadno(&c,cadno,ia+1,ja+1,Cdd,Cddis,Cddismx);
			}
		}
	}
	if(coor2->type==SULFUR&&coor1->type==MAINCHN){
		for(ia=0;ia<2;ia++){
			for(ja=0;ja<6;ja++){
				Cddis=distance2(coor1->Cds[ia],coor2->Cds[ja]);
				if(Cddis<MINCDDIS2)	Cddismx=assigncadno(&c,cadno,ia+1,ja+1,Cdd,Cddis,Cddismx);
			}
		}
	}
	if(coor2->type==MAINCHN&&coor1->type==SULFUR){
		for(ia=0;ia<6;ia++){
			for(ja=0;ja<2;ja++){
				Cddis=distance2(coor1->Cds[ia],coor2->Cds[ja]);
				if(Cddis<MINCDDIS2)	Cddismx=assigncadno(&c,cadno,ia+1,ja+1,Cdd,Cddis,Cddismx);
			}
		}
	}

	if(c==0) return 0;

	coor_dis=distance2(coor1->coor,coor2->coor);
	if(coor_dis<MINCOORDIS2) return 0;
	/*carbon_dis=distance2(oxy1->carbon,oxy2->carbon);
	if(carbon_dis<MINCARBONDIS2) return 0;*/

	min_dis=distance2(coor1->link,coor2->coor);
	if(min_dis<MINLINKDIS2) return 0;
	link_dis=distance2(coor1->thirdatom,coor2->coor);
	if(link_dis<MINLINKDIS2) return 0;
	link_dis=distance2(coor2->link,coor1->coor);
	if(link_dis<MINLINKDIS2) return 0;
	link_dis=distance2(coor2->thirdatom,coor1->coor);
	if(link_dis<MINLINKDIS2) return 0;
	link_dis=distance2(coor2->link,coor1->link);
	if(link_dis<MINLINKDIS2) return 0;
	link_dis=distance2(coor1->thirdatom,coor2->link);
	if(link_dis<MINLINKDIS2) return 0;
	link_dis=distance2(coor2->thirdatom,coor1->link);
	if(link_dis<MINLINKDIS2) return 0;
	link_dis=distance2(coor2->thirdatom,coor1->thirdatom);
	if(link_dis<MINLINKDIS2) return 0;

	return c;
}

void coor_linkres(COORGEN *coor, int r)
{
	int i;
	int lr;
	int typ,typ0;

	for(i=0;i<coor->linkresN;i++){
		lr=coor->linkres[i];
		if(lr/5==r/5){
			typ0=lr%5;
			typ=r%5;
			if(typ0!=SULFUR){
				if(typ==SULFUR)	coor->linkres[i]=r;
				else if(typ0!=DOUBLEOXYGEN){
					if(typ==DOUBLEOXYGEN) coor->linkres[i]=r;
					else if(typ0!=NITROGEN&&typ==NITROGEN){
						coor->linkres[i]=r;
					}
				}
			}
			break;
		}
	}
	if(i>=coor->linkresN){
		coor->linkres[i]=r;
		(coor->linkresN)++;
		/*if((oxy->linkresN)>=30)
			printf("stop\n");*/
	}
}

int gen_coorpairs(PROCHAIN *p,RESBACKRUB *backrubs,MUTATION *m[4], COORGEN *coorlib, int *res_coorindex, COORPAIR *coorpairs,int *coorpairindex)
{
	int ir,jr;
	int ires,jres;
	int ityp,jtyp;
	int ic;
	int oi,oj;
	int count=0;
	int indexindex;
	int isp;
	int neighbor;
	int rubi=-1;
	int rubj=-1;
	int ucompi;
	char cadno[MAXSAVEDCAD][2];

	/*oxypairindex[0]=0; */
	for(ir=0;ir<p->resn*5;ir++){
		/*printf("%3d %6d\n",ir+1,count);*/
		for(jr=0;jr<ir;jr++){
			indexindex=ir*(ir-1)/2+jr;
			coorpairindex[indexindex*2]=count;
			ires=ir/5;jres=jr/5;
			ityp=ir%5;jtyp=jr%5;
			if((ires==jres)&&jtyp!=MAINCHN){
				coorpairindex[indexindex*2+1]=count;
				continue;
			}
			neighbor=0;
			if(ires==jres+1&&ityp!=MAINCHN&&jtyp!=MAINCHN) neighbor=1;
			else if(ires==jres&&jtyp==MAINCHN) neighbor=2;
			else if(ires==jres+1&&ityp!=MAINCHN&&jtyp==MAINCHN) neighbor=3;

			for(oi=res_coorindex[ir];oi<res_coorindex[ir+1];oi++){
				if(neighbor!=0){
					if(coorlib[oi].muttype==-1) rubi=-1;
					else rubi=m[coorlib[oi].muttype][coorlib[oi].source].backrub;
					if(neighbor==2&&backrubs[rubi].native_compatible[0]==0)
						continue;
					if(neighbor==3&&backrubs[rubi].native_compatible[1]==0)
						continue;
				}
				for(oj=res_coorindex[jr];oj<res_coorindex[jr+1];oj++){
					if(neighbor==1){
						rubj=m[coorlib[oj].muttype][coorlib[oj].source].backrub;
						for(ucompi=0;ucompi<backrubs[rubi].pre_uncompatible[0];ucompi++){
							if(rubj==backrubs[rubi].pre_uncompatible[ucompi+1])
								break;
						}
						if(ucompi<backrubs[rubi].pre_uncompatible[0])
							continue;
					}
					isp=Is_pair(coorlib+oi,coorlib+oj,cadno);
					if(isp!=0){
						coorpairs[count].A_index=oi;
						coorpairs[count].B_index=oj;
						for(ic=0;ic<isp;ic++){
							coorpairs[count].A_cadNo[ic]=cadno[ic][0];
							coorpairs[count].B_cadNo[ic]=cadno[ic][1];
						}
						for(ic=isp;ic<MAXSAVEDCAD;ic++){
							coorpairs[count].A_cadNo[ic]=cadno[0][0];
							coorpairs[count].B_cadNo[ic]=cadno[0][1];
						}
						coor_linkres(coorlib+oi,jr);
						coor_linkres(coorlib+oj,ir);
						count++;
					}
				}
			}
			coorpairindex[indexindex*2+1]=count;
		}
	}

	return count;
}

void simple_coorlink(int resN, int coorN, COORGEN *coorlib, int *coorpairN,COORPAIR *coorpairs,int *coorpairindex,int *coorpairsim,int *coorpairsimindex)
{
	int io;
	int links;
	int ilink;
	int ir,jr;
	int ipair,jpair;
	int oA,oB;
	int indexindex;
	int typ;

	for(io=0;io<(coorN);io++){
		if(coorlib[io].type==DOUBLEOXYGEN||coorlib[io].type==NITROGEN||coorlib[io].type==SULFUR)
			links=2;
		else 
			links=1;
		for(ilink=0;ilink<coorlib[io].linkresN;ilink++){
			typ=(coorlib[io].linkres[ilink])%5;
			if(typ==DOUBLEOXYGEN||typ==NITROGEN||typ==SULFUR)
				links+=2;
			else
				links++;
		}
		/*printf("%3d: %3d,%3d-->%3d\n",io+1,oxylib[io].res,oxylib[io].muttype,links);*/
		if(links<6) coorlib[io].linkresN=0;
	}

	for(ipair=0;ipair<(*coorpairN);ipair++){
		oA=coorpairs[ipair].A_index;
		oB=coorpairs[ipair].B_index;
		if(coorlib[oA].linkresN==0||coorlib[oB].linkresN==0){
			coorpairs[ipair].A_index=-1;
			coorpairs[ipair].B_index=-1;
		}
	}

	*coorpairN=0;
	coorpairsimindex[0]=0;
      jpair=0;
	for(ir=0;ir<resN*5;ir++){
		/*printf("%3d %6d\n",ir+1,(*coorpairN));*/
		for(jr=0;jr<ir;jr++){
			indexindex=ir*(ir-1)/2+jr;
			/*if(coorpairindex[indexindex*2+1]>coorpairindex[indexindex*2])
			printf("%d,%d--%d,%d \n",ir,jr,coorpairindex[indexindex*2],coorpairindex[indexindex*2+1]);*/
			for(ipair=coorpairindex[indexindex*2];ipair<coorpairindex[indexindex*2+1];ipair++){
				if(coorpairs[ipair].A_index!=-1&&coorpairs[ipair].B_index!=-1){
	             		coorpairsim[jpair]=ipair;
	                		jpair++;
				}
			}
            	coorpairsimindex[indexindex+1]=jpair;
		}
	}
	(*coorpairN)=jpair;
}

