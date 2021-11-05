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

#ifndef _datatype_h
#define _datatype_h


typedef struct{
	int at_res;
	char name[15];
	float vdw_rad;
	float x[3];
}ATM;

typedef struct{
	char name[10];
	int mutability;
	int start;
	int end;
	int sstart;
	int send;
	int pre;
	int next;
}RES;

typedef struct{
	int atmn;
	int resn;
	ATM *a;
	RES *r;
}PROCHAIN;

typedef struct{
	char pdbname[100];
	/*int model_n;*/
	PROCHAIN pro;
}PROTEIN;

#define BACKRUBAN   3

typedef struct{
	int at_res;
	int backrub;
	int exist;
	float mx[4][3];
	float COx[2][3];
	float Nx[3];
	float phipsi[2];
	int native_compatible[2];
	int pre_uncompatible[BACKRUBAN];
} RESBACKRUB;

extern float rub_angle_plan[BACKRUBAN];

#define MUT2GLU 0
#define MUT2ASP 1
#define MUT2HIS 2
#define MUT2CYS 3

#define MSHIFTN 1

extern float mut_shift_plan[MSHIFTN];

typedef struct{
	int at_res;
	int backrub;
	int rota;
	int rotshift;
	int muttype;
	int satn;
	int isnative;
	float sx[6][3];
	float mut_score;
} MUTATION;

extern float main_chain_template[3][3];
extern float EDHC_template[4][6][3];
extern int EDHC_sidean[4];
extern int EDHC_chiN[4];
extern float EDHC_chi[4][3];

#define MAINCHN 0
#define DOUBLEOXYGEN 1
#define SINGLEOXYGEN 2
#define NITROGEN 3
#define SULFUR 4

#define MUT2GLUO1 4
#define MUT2ASPO1 5
#define MUT2GLUO2 6
#define MUT2ASPO2 7
#define MUT2HISN1 8
#define MUT2HISN2 9

typedef struct{
	int res;
	int muttype;
	int source;
	int type;
	int coorat;
	float coor[3];
	float link[3];
	float thirdatom[3];
	float Cd[3];
	float Cds[6][3];
	int linkresN;
	int linkres[30];
}COORGEN;

#define MAXSAVEDCAD 4

typedef struct{
	int A_index;
	int B_index;
	char A_cadNo[MAXSAVEDCAD];
	char B_cadNo[MAXSAVEDCAD];
}COORPAIR;

typedef struct{
	int coorN;
	int coorresN;
	int coor[6];
	int coorres[6];
	int coorat[6];
	int coortyp[6];
	int combine_score;
	float Cd[3];
}COORCOMBINE;

typedef struct{
	char pdbname[100];
	float cadmium[3];
	int res_n;
	int coorn;
	int res[6];
	char res_name[6][10];
	int res_type[6];
	int coor_at[6];
	int res_rub[6];
	int rotamer[6];
	int rotshift[6];
	int coor_score;
}RECORD;
struct LIST{
	int i;
	struct LIST *n;
};
struct LIST *order;

#define MAXRECORD 10000

#endif

