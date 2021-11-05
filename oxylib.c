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

#include "oxylib.h"
#include "mem.h"
#include <stdio.h>
#include <string.h>
#include "geometry.h"
#include <math.h>

void get_Cd_doubleoxy(float root[3],float x1[3],float x2[3], float cadmium[3])
{
	float v1[3],v2[3],v3[3];
	int i;
	float len3;

	direction( x1, root, v1);
	direction( x2, root, v2);
	for(i=0;i<3;i++){
		v3[i]=v1[i]+v2[i];
	}
	len3=get_len(v3);
	for(i=0;i<3;i++){
		cadmium[i]=2.80*(v3[i])/len3+root[i];
	}
}
void get_Cd_singleoxy(float root[3],float x1[3],float x2[3], float cadmium[3])
{
	float v1[3],v2[3],v3[3],v4[3];
	int i;

	direction( x1, root, v1);
	direction( x2, root, v2);
	oprod(v1, v2, v3);
	rotvex(v1,v4,v3,60.0);

	for(i=0;i<3;i++){
		cadmium[i]=2.30*v4[i]+x1[i];
	}
}
void get_Cd_mainchain(float root[3],float x1[3],float x2[3], float cadmium[3], float cadmiumsub[][3])
{
	float v[3],v1[3],v2[3],v3[3],v4[3];
	float len;
	int i;

    for(i=0;i<3;i++){
            v[i]=x1[i]-root[i];
    }
    len=get_len(v);
    for(i=0;i<3;i++){
    	cadmium[i]=2.30*(v[i])/len+root[i];
    }

    direction( x1, root, v1);
	direction( x2, root, v2);
	oprod(v1, v2, v3);
	rotvex(v1,v4,v3,60.0);
	for(i=0;i<3;i++){
		cadmiumsub[0][i]=2.30*v4[i]+x1[i];
	}
	rotvex(v1,v4,v3,-60.0);
	for(i=0;i<3;i++){
		cadmiumsub[1][i]=2.30*v4[i]+x1[i];
	}
}
void get_Cd_nitrogen(float root[3],float x1[3],float x2[3], float cadmium[3])
{
	float v1[3],v2[3],v3[3];
	int i;
	float len3;

	direction( x1, root, v1);
	direction( x1, x2, v2);
	for(i=0;i<3;i++){
		v3[i]=v1[i]+v2[i];
	}
	len3=get_len(v3);
	for(i=0;i<3;i++){
		cadmium[i]=2.23*(v3[i])/len3+x1[i];
	}
}
int getmin(float x[3])
{
	float min=1.1;
	int i,mini;

	for(i=0;i<3;i++){
		if(x[i]>0&&x[i]<min){
			min=x[i];
			mini=i;
		}
		else if(x[i]<0&&(-x[i])<min){
			min=-x[i];
			mini=i;
		}
	}
	return mini;
}
void get_Cd_sulfur(float root[3],float x[3], float cadmium[3], float cadmiumsub[][3])
{
    float v[3],v1[3],v2[3],v3[3],v4[3];
    /*float len;*/
    int i,j,k;
    int mini;
    float sum;

   /* for(i=0;i<3;i++){
            v[i]=x[i]-root[i];
    }
    len=get_len(v);
    for(i=0;i<3;i++){
    	cadmium[i]=2.40*(v[i])/len+root[i];
    }*/

    direction( x, root, v);
    mini=getmin(v);
    sum=0.0;
    for(i=0;i<3;i++){
    	if(i!=mini)	sum+=v[i]*v[i];
    }
    for(i=0;i<3;i++){
    	if(i==mini) v1[i]=0.0;
    	else{
    		for(j=0;j<3;j++){
    			if(j==mini) continue;
    			if(j!=i) v1[i]=sqrt(v[j]*v[j]/sum);
    		}
    	}
    }

	oprod(v, v1, v2);
	rotvex(v,v3,v2,60.0);
	for(i=0;i<3;i++){
		cadmium[i]=2.40*v3[i]+x[i];
	}
	for(k=0;k<6;k++){
		rotvex(v3,v4,v,360.0/7*(k+1));
		for(i=0;i<3;i++){
			cadmiumsub[k][i]=2.40*v4[i]+x[i];
		}
	}
}


void get_coor(COORGEN *coor, int res,int type, int source, float x[3], float root[3], float x2[3])
{
	int i;
	
	coor->res=res;
	coor->coorat=type;
	coor->source=source;

	if(type==-1) coor->type=MAINCHN;
	else if(type==MUT2GLU||type==MUT2ASP) coor->type=DOUBLEOXYGEN;
	else if(type==MUT2GLUO1||type==MUT2GLUO2||type==MUT2ASPO1||type==MUT2ASPO2) coor->type=SINGLEOXYGEN;
	else if(type==MUT2HISN1||type==MUT2HISN2) coor->type=NITROGEN;
	else if(type==MUT2CYS) coor->type=SULFUR;

	if(type==-1) coor->muttype=-1;
	else if(type==MUT2GLU||type==MUT2GLUO1||type==MUT2GLUO2) coor->muttype=MUT2GLU;
	else if(type==MUT2ASP||type==MUT2ASPO1||type==MUT2ASPO2) coor->muttype=MUT2ASP;
	else if(type==MUT2HISN1||type==MUT2HISN2) coor->muttype=MUT2HIS;
	else if(type==MUT2CYS) coor->muttype=MUT2CYS;

	for(i=0;i<3;i++){
		coor->coor[i]=x[i];
		coor->link[i]=root[i];
		coor->thirdatom[i]=x2[i];
	}
	coor->linkresN=0;
	if(coor->type==MAINCHN){
		get_Cd_mainchain(root,x,x2,coor->Cd,coor->Cds);
	}
	else if(coor->type==DOUBLEOXYGEN){
		get_Cd_doubleoxy(root,x,x2,coor->Cd);
	}
	else if(coor->type==SINGLEOXYGEN){
		get_Cd_singleoxy(root,x,x2,coor->Cd);
	}
	else if(coor->type==NITROGEN){
		get_Cd_nitrogen(root,x,x2,coor->Cd);
	}
	else if(coor->type==SULFUR){
		get_Cd_sulfur(root,x,coor->Cd,coor->Cds);
	}
}


int gen_COOR_lib(PROCHAIN *p,MUTATION *m[4], int *res_mutindex[4], COORGEN *coorlib, int *res_coorindex)
{
	int ir,ia,im;
	int count=0;
	int start,end,t;
      float *Ox,*Cx,*CAx;

	res_coorindex[0]=0;
	for(ir=0;ir<p->resn;ir++){
		Ox=Cx=CAx=NULL;
            for(ia=p->r[ir].start;ia<p->r[ir].sstart;ia++){
                        if(!strncmp(p->a[ia].name, "O  ",3)) Ox=p->a[ia].x;
                        else if(!strncmp(p->a[ia].name, "C  ",3)) Cx=p->a[ia].x;
                        else if(!strncmp(p->a[ia].name, "CA ",3)) CAx=p->a[ia].x;
                }
             if(Ox==NULL||Cx==NULL||CAx==NULL){
                        for(t=MAINCHN;t<=SULFUR;t++) res_coorindex[ir*5+t+1]=count;
                        continue;
                }
		get_coor(coorlib+count,ir,-1,-1,Ox, Cx, CAx);
		count++;
		res_coorindex[ir*5+MAINCHN+1]=count;

		start=res_mutindex[MUT2GLU][ir];
		end=res_mutindex[MUT2GLU][ir+1];
		for(im=start;im<end;im++){
			get_coor(coorlib+count,ir,MUT2GLU,im,m[MUT2GLU][im].sx[3], m[MUT2GLU][im].sx[2], m[MUT2GLU][im].sx[4]);
			count++;
		}
		start=res_mutindex[MUT2ASP][ir];
		end=res_mutindex[MUT2ASP][ir+1];
		for(im=start;im<end;im++){
			get_coor(coorlib+count,ir,MUT2ASP,im,m[MUT2ASP][im].sx[2], m[MUT2ASP][im].sx[1], m[MUT2ASP][im].sx[3]);
			count++;
		}
		res_coorindex[ir*5+DOUBLEOXYGEN+1]=count;
             if(p->r[ir].mutability==0){
                        for(t=MAINCHN+1;t<=SULFUR;t++) res_coorindex[ir*5+t+1]=count;
                        continue;
                }

		start=res_mutindex[MUT2GLU][ir];
		end=res_mutindex[MUT2GLU][ir+1];
		for(im=start;im<end;im++){
			get_coor(coorlib+count,ir,MUT2GLUO1,im,m[MUT2GLU][im].sx[3], m[MUT2GLU][im].sx[2], m[MUT2GLU][im].sx[4]);
			count++;
			get_coor(coorlib+count,ir,MUT2GLUO2,im,m[MUT2GLU][im].sx[4], m[MUT2GLU][im].sx[2], m[MUT2GLU][im].sx[3]);
			count++;
		}
		start=res_mutindex[MUT2ASP][ir];
		end=res_mutindex[MUT2ASP][ir+1];
		for(im=start;im<end;im++){
			get_coor(coorlib+count,ir,MUT2ASPO1,im,m[MUT2ASP][im].sx[2], m[MUT2ASP][im].sx[1], m[MUT2ASP][im].sx[3]);
			count++;
			get_coor(coorlib+count,ir,MUT2ASPO2,im,m[MUT2ASP][im].sx[3], m[MUT2ASP][im].sx[1], m[MUT2ASP][im].sx[2]);
			count++;
		}
		res_coorindex[ir*5+SINGLEOXYGEN+1]=count;


		start=res_mutindex[MUT2HIS][ir];
		end=res_mutindex[MUT2HIS][ir+1];
		for(im=start;im<end;im++){
			get_coor(coorlib+count,ir,MUT2HISN1,im,m[MUT2HIS][im].sx[3], m[MUT2HIS][im].sx[1], m[MUT2HIS][im].sx[4]);
			count++;
			get_coor(coorlib+count,ir,MUT2HISN2,im,m[MUT2HIS][im].sx[5], m[MUT2HIS][im].sx[2], m[MUT2HIS][im].sx[4]);
			count++;
		}
		res_coorindex[ir*5+NITROGEN+1]=count;

		start=res_mutindex[MUT2CYS][ir];
		end=res_mutindex[MUT2CYS][ir+1];
		for(im=start;im<end;im++){
			get_coor(coorlib+count,ir,MUT2CYS,im,m[MUT2CYS][im].sx[1], m[MUT2CYS][im].sx[0], m[MUT2CYS][im].sx[0]);
			count++;
		}
		res_coorindex[ir*5+SULFUR+1]=count;
		/*printf("S   %6d\n",count);*/
	}
	
	return count;
}


