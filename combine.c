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
#include <math.h>
#include "combine.h"
#include "coordination.h"
#include "mem.h"

int cadnoequ(char xa[MAXSAVEDCAD],char xb[MAXSAVEDCAD],char *ya,char *yb)
{
	int i,j;

	for(i=0;i<MAXSAVEDCAD;i++){
		for(j=0;j<MAXSAVEDCAD;j++){
			if((xa[i]==ya[j])&&(xb[i]==yb[j])) return 1;
		}
	}
	return 0;
}
int cadnoequ2(char xa[MAXSAVEDCAD],char *ya)
{
	int i,j;

	for(i=0;i<MAXSAVEDCAD;i++){
		for(j=0;j<MAXSAVEDCAD;j++){
			if(xa[i]==ya[j]) return 1;
		}
	}
	return 0;
}

int bin_search_coor(COORPAIR *coorpairs,int *coorpairsim,int low,int high,int a,int b,char xa[MAXSAVEDCAD],char xb[MAXSAVEDCAD])
{
	float k;
	int mid;
	float t;
	int i;
	char* ya,*yb;

	if(high-low<64){
		for(i=low;i<=high;i++){
			if(coorpairs[coorpairsim[i]].A_index==a&&coorpairs[coorpairsim[i]].B_index==b){
				ya=coorpairs[coorpairsim[i]].A_cadNo;
				yb=coorpairs[coorpairsim[i]].B_cadNo;
				if(cadnoequ(xa,xb,ya,yb)==1)
					return 1;
			}
		}
	}
	else{
		k=1.0*a*(a-1)/2+b;
		while(low<=high){
			mid=(low+high)/2;
			t=coorpairs[coorpairsim[mid]].A_index;
			t=t*(t-1)/2+coorpairs[coorpairsim[mid]].B_index;
			if(t==k){
				ya=coorpairs[coorpairsim[mid]].A_cadNo;
				yb=coorpairs[coorpairsim[mid]].B_cadNo;
				if(cadnoequ(xa,xb,ya,yb)==1)
					return 1;
			}
			if(k>t) low=mid+1;
			else high=mid-1;
		}
	}
	return -1;
}
void bin_search_many_coor(COORPAIR *coorpairs,int *coorpairsim,int a,int low,int high,int *start,int *end)
{
	int mid;

	while(low<=high){
		mid=(low+high)/2;
		if(a==coorpairs[coorpairsim[mid]].A_index) break;
		if(a>coorpairs[coorpairsim[mid]].A_index) low=mid+1;
		else high=mid-1;
	}
	if(low>high) (*start)=-1;
	else{
		(*start)=mid;
		while(coorpairs[coorpairsim[(*start)]].A_index==a){
			(*start)--;
			if((*start)<low) break;
		}
		(*start)++;
		(*end)=mid;
		while(coorpairs[coorpairsim[(*end)]].A_index==a){
			(*end)++;
			if((*end)>high) break;
		}
	}
}

int search_t0_coor(COORPAIR *coorpairs,int *coorpairsim,int *coorpairsimindex,int a,int rmax,int *t0match)
{
	int ir,i;
	int c=0;
	int start,end;

	for(ir=0;ir<rmax;ir++){
		bin_search_many_coor(coorpairs,coorpairsim,a,coorpairsimindex[rmax*(rmax-1)/2+ir],coorpairsimindex[rmax*(rmax-1)/2+ir+1]-1,&start,&end);
		if(start!=-1){
			for(i=start;i<end;i++){
				t0match[c]=i;
				c++;
			}
		}
	}

	return c;
}

#define MINCOORSCORE 38

int combine_coordination(COORGEN *coorlib,int coorpairN,COORPAIR *coorpairs,int *coorpairsim,int *coorpairsimindex,COORCOMBINE *combine_coor,MUTATION *mut[4])
{
	int m,n;
	int t0[6],t[6],d[6],r[6];
	char x[6][MAXSAVEDCAD];
	int zero;
	int f;
	int res_pair;
	int *t0match[7];
	int t0matchN[7];
	int count=0;
	int j,xi;
	int score;
	/*int resN;*/
	int coorN;
	char* y;
	int pass,passi;
	/*int combinescore=0;*/

	for(n=0;n<6;n++)for(xi=0;xi<MAXSAVEDCAD;xi++){x[n][xi]=0;}
	/*resN=p->resn;*/
	n=0;j=0;zero=0;
	d[0]=d[1]=d[2]=d[3]=d[4]=d[5]=0;
	r[0]=r[1]=r[2]=r[3]=r[4]=r[5]=0;
	t0match[0]=t0match[1]=t0match[2]=t0match[3]=t0match[4]=t0match[5]=NULL;
	t0[0]=t0[1]=t0[2]=t0[3]=t0[4]=t0[5]=-1;

	while(d[0]<coorpairN-1){
		t[0]=coorpairs[coorpairsim[d[0]]].A_index;
		r[0]=(coorlib[t[0]].res)*5+coorlib[t[0]].type;
		for(xi=0;xi<MAXSAVEDCAD;xi++)x[0][xi]=coorpairs[coorpairsim[d[0]]].A_cadNo[xi];
		n++;
		t[n]=coorpairs[coorpairsim[d[n-1]]].B_index;
		r[n]=(coorlib[t[n]].res)*5+coorlib[t[n]].type;
		for(xi=0;xi<MAXSAVEDCAD;xi++)x[n][xi]=coorpairs[coorpairsim[d[n-1]]].B_cadNo[xi];
		f=1;
		do{
			zero=1;
			if(f==1){
				for(m=1;m<n;m++){
					res_pair=r[m]*(r[m]-1)/2+r[n];
					if(bin_search_coor(coorpairs,coorpairsim,coorpairsimindex[res_pair],coorpairsimindex[res_pair+1]-1,t[m],t[n],x[m],x[n])==-1){
						zero=1;
						break;
					}
				}
				if(m>=n){
					n++;
					coorN=0;
					/*combinescore=1;*/
					score=MINCOORSCORE+1;
					if(n>=4){
					/*if(n==3){*/
						/*if(t[0]==50690&&t[1]==49531&&t[2]==28062&&t[3]==26000)
							printf("findfindfindfind%3d %6d %6d %6d %6d\n",count, t[0],t[1],t[2],t[3]);*/
						coorN=0;
						combine_coor[count].coorresN=n;
						for(j=0;j<n;j++){
							if(coorlib[t[j]].type==DOUBLEOXYGEN||coorlib[t[j]].type==NITROGEN||coorlib[t[j]].type==SULFUR)
								coorN+=2;
							else coorN++;
						}
						/*if(coorN>=6&&coorN<=8){*/
						if(coorN>=5&&coorN<=8){
							for(j=0;j<n;j++){
								combine_coor[count].coor[j]=t[j];
								combine_coor[count].coorres[j]=coorlib[t[j]].res;
								combine_coor[count].coorat[j]=coorlib[t[j]].coorat;
								combine_coor[count].coortyp[j]=coorlib[t[j]].type;
							}
							combine_coor[count].coorN=coorN;
							score=do_coordination(combine_coor+count,coorlib,mut,x);
							if(score>=MINCOORSCORE){
								/*printf("%3d (%d %6d) %6d %6d %6d %6d\n",count, coorN, score,t[0],t[1],t[2],t[3]);*/
								combine_coor[count].combine_score=score;
								count++;
							}
						}
					/*else if(n>=3){
						combine_coor[count].coorresN=n;
						for(j=0;j<n;j++){
							combine_coor[count].coor[j]=t[j];
							combine_coor[count].coorres[j]=coorlib[t[j]].res;
							combine_coor[count].coorat[j]=coorlib[t[j]].coorat;
							combine_coor[count].coortyp[j]=coorlib[t[j]].type;
						}
						combinescore=do_coordination2(combine_coor+count,coorlib,mut,x);
					}
					if(n==6||coorN>=8||combinescore<0||score<MINCOORSCORE){*/
						n--;
						zero=1;
					}
					else
						zero=0;
				}
			}
			if(zero==0){
				if(t[0]!=t0[n]){
					BStd_renew(t0match[n],coorpairsimindex[r[0]*(r[0]-1)/2+r[0]]-coorpairsimindex[r[0]*(r[0]-1)/2]);
					t0matchN[n]=search_t0_coor(coorpairs,coorpairsim,coorpairsimindex,t[0],r[0],t0match[n]);
					t0[n]=t[0];
				}
				if(t0matchN[n]!=0){
					pass=0; passi=0;
					while(pass==0&&passi<t0matchN[n]){
						y=coorpairs[coorpairsim[t0match[n][passi]]].A_cadNo;
						pass=cadnoequ2(x[0],y);
						if(pass==0)	passi++;
					}
					if(pass==0) f=0;
					else{
						t[n]=coorpairs[coorpairsim[t0match[n][passi]]].B_index;
						r[n]=(coorlib[t[n]].res)*5+coorlib[t[n]].type;
						for(xi=0;xi<MAXSAVEDCAD;xi++) x[n][xi]=coorpairs[coorpairsim[t0match[n][passi]]].B_cadNo[xi];
						d[n-1]=passi;
						f=1;
					}
				}
				else f=0;
			}
			else{
				d[n-1]++;
				if(d[n-1]>=t0matchN[n]) f=0;
				else{
					pass=0; passi=d[n-1];
					while(pass==0&&passi<t0matchN[n]){
						y=coorpairs[coorpairsim[t0match[n][passi]]].A_cadNo;
						pass=cadnoequ2(x[0],y);
						if(pass==0)	passi++;
					}
					if(pass==0) f=0;
					else{
						t[n]=coorpairs[coorpairsim[t0match[n][passi]]].B_index;
						r[n]=(coorlib[t[n]].res)*5+coorlib[t[n]].type;
						for(xi=0;xi<MAXSAVEDCAD;xi++) x[n][xi]=coorpairs[coorpairsim[t0match[n][passi]]].B_cadNo[xi];
						f=1;
					}
				}
			}
			if(f==0) {
				n--;
				zero=1;
			}
		}while(!(f==0&&n==1));
		d[0]++;
		n--;
	}

	for(j=0;j<4;j++){
		BStd_free(t0match[j]);
	}

	return count;
}
int compare_combine(COORCOMBINE *combineA,COORCOMBINE *combineB)
{
	int i;
	
	if(combineA->coorresN<combineB->coorresN)
		return 1;
	if(combineA->coorresN>combineB->coorresN)
		return 2;
	for(i=0;i<combineA->coorresN;i++){
		if(combineA->coorres[i]<combineB->coorres[i]) return 1;
		if(combineA->coorres[i]>combineB->coorres[i]) return 2;
	}
	for(i=0;i<combineA->coorresN;i++){
		if(combineA->coorat[i]<combineB->coorat[i]) return 1;
		if(combineA->coorat[i]>combineB->coorat[i]) return 2;
	}
	return 0;
}

void shellsort_lable(int combineN, COORCOMBINE *combine,int *index)
{
	int i,j,flag=1;
	int n=combineN;
	int gap=sqrt(n)*2;
	int tmp;

	for(i=0;i<combineN;i++) index[i]=i;
	while(gap>1){
		gap=gap/2;
		do{
			flag=0;
			for(i=0;i<n-gap;i++){
				j=i+gap;
				if(compare_combine(combine+index[i],combine+index[j])==1){
					tmp=index[i];
					index[i]=index[j];
					index[j]=tmp;
					flag=1;
				}
			}
		}while(flag!=0);
	}
}

int combine_res_index(int combineN, COORCOMBINE *combine,int *index, int *resindex)
{
	int i;
	int count=0;

	resindex[0]=0;

	for(i=1;i<combineN;i++){
		while(i<combineN&&compare_combine(combine+index[resindex[count]],combine+index[i])==0){
			i++;
		}
		resindex[count+1]=i;
		count++;
	}

	return count;
}

int arrange_combine(int combineN, COORCOMBINE *combine,int *index, int *resindex)
{
	int rescombineN=0;

	shellsort_lable(combineN, combine,index);
	rescombineN=combine_res_index(combineN, combine,index, resindex);
	
	return rescombineN;
}

