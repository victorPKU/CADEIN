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

#include "datatype.h"

/*float rub_angle_plan[BACKRUBAN]={ -8.0,-4.0,0.0,4.0, 8.0};
float rub_angle_plan[BACKRUBAN]={0.0};
float rub_angle_plan[BACKRUBAN]={ -9.0, -6.0,-3.0,0.0,3.0, 6.0, 9.0};*/
float rub_angle_plan[BACKRUBAN]={ -6.0,0.0,6.0};

/*float mut_shift_plan[MSHIFTN]={ -0.7, -0.35, 0.0, 0.35, 0.7};
float mut_shift_plan[MSHIFTN]={ -0.65,  0.0, 0.65};*/
float mut_shift_plan[MSHIFTN]={   0.0};

float main_chain_template[3][3]={
	{ 1.449, 0.000, 0.000}, /*N*/
	{ 0.000, 0.000, 0.000}, /*CA*/
	{-0.523, 1.429, 0.000}, /*C*/
};

float EDHC_template[4][6][3]={
	{
	{-0.514,-0.773, 1.211}, /*GLU*/
	{-2.037,-0.851, 1.156}, 
	{-2.583,-1.712, 2.287}, 
	{-1.749,-2.214, 3.071}, 
	{-3.824,-1.851, 2.346}
	},
	{
	{-0.514,-0.773, 1.211}, /*ASP*/
	{-2.012,-1.024, 1.105},
	{-2.591,-0.584, 0.088},
	{-2.550,-1.652, 2.043}
	},
	{
	{-0.514,-0.773, 1.211}, /*HIS*/
	{-2.005,-0.786, 1.265},
	{-2.895,-1.480, 0.512},
	{-2.709,-0.010, 2.159},
	{-3.989,-0.264, 1.971},
	{-4.150,-1.140, 0.975}
	},
	{
	{-0.514,-0.773, 1.211}, /*CYS*/
	{-2.245,-0.493, 1.643}
	}
};
int EDHC_sidean[4]={5,4,6,2};
int EDHC_chiN[4]={3,2,2,1};

float EDHC_chi[4][3]={
	{-176.20, 175.38,   0.04},
	{-169.10,   0.02},
	{179.25,-105.16},
	{163.26}
};


