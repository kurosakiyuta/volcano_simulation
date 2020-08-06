#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "solubilityH2O.h"
#include "viscosity.h"
#include "/Users/DELL/source/C++/���l�v�ZC++/���l�v�ZC++/nrutil.h"
#include "rk4c.h"
#include "density.h"
#include "velocity.h"
#include "porosity.h"
#include "speed of sound.h"
#include "Exsolvedgas.h"
#include "wt to mol.h"
#include "Mach number.h"
#include "rk4(crystallization).h"
#include "xeq.h"
#include "relativeviscosity.h"
#include "density(liquid).h"
#include "xc.h"
#include "equilibrium.h"
#include "meltcomposition.h"
#include "wttomol.h"
#include "MPF(bimodal).h"
#include "MPF(mono).h"

#pragma warning(disable : 4996)


#define n 1 /*�v�f���i���������j*/
#define nc 7 /*�v�f���i���������j*/
#define g 9.80665

#define SiO2 0
#define TiO2 1
#define Al2O3 2
#define Fe2O3 3
#define FeO 4
#define MnO 5
#define MgO 6
#define CaO 7
#define Na2O 8
#define K2O 9
#define P2O5 10
#define F2O 0
#define H2O 11

#define R 461.88888
#define Mw 0.01802
#define cpl 2650  //�Β��΂̖��x high albite Rock forming minerals 
#define cpy 3585  //�P�΂̖��x Enstatite-ferrosilite �̒��Ԓn
#define co 3700  //�I���r���̖��x
#define cb 3200  //�����̖��x�i�K��
#define cq 2196 //wiki
#define cs 3600
#define cm 2400 //melt density
#define ccl 3585//�C���̕K�v��
#define cor 3585


/************/

#define Qi 35.999518* pow(10,6)  //����,�̐ϗ���kg/m~2s
#define z0 6000
#define Pii 150 * pow(10,6)

//******�ϐ��萔******/
#define cm 2400 //melt density
#define T 1010+273.65  //Kelvin
#define G pow(10,-0)
#define ��q 0.8
#define h  0.01
#define r 30
#define rp 4

double tpl[7] = { 62.838700967511066 , -102.95890679690591 , 92.72184036973631 , -26.844087952274144 };
double tpy[7] = { 16.8719986417709 , -16.89469887330756 , 19.36556465691195 , -6.331621152958288 };
double to[7] = { 0.7754808809083791 , -68.57050654653938 , 3426.3309616672377 , -75287.0696997242 };
double tq[7] = {};
double ts[7] = { 11.585166136451374 , -5.659273484872328 , 3.3929146667837213 , -0.5537835441752051 };
double tcl[7] = {};
double tor[7] = {};

double tw[2] = { -9.721707911529391 , 0.6004911319169789 };

double comp_melt_ini[11] = { 57.0211,	0.8541,	17.4296,	1.1859,	3.3842,	0.3558,	1.3805,	5.4095,	5.9228,	3.2678,	0.8105 };

//                          0=SiO2    1=TiO2	2=Al2O3	3=Fe2O3 4=FeO   5=MnO   6=MgO	7=CaO  	8=Na2O  9=K2O   10= P2O5
double comp_pl[11] = {};
double comp_py[11] = {};
double comp_o[11] = {};
double comp_q[11] = {};
double comp_s[11] = {};
double comp_cl[11] = {};
double comp_or[11] = {};

double c[8] = { cpl,cpy,co,cq,cs,cm,ccl,cor };

/************************************�������̒�`**********************************/
/*********************************�C�A�̂Ȃ��^��������*****************************/
void ryuutai(double z,double b[], double P[], double dPdz[],double watersuminitial,double vis,double bi[]) {
	double Cw, c_bulk, v, visco,S,M,xc[nc];
	crystalwt(b, xc, c);
	v = velosity(P[0], Qi, T, tw,watersuminitial, c, b, r);/*���x��*/
	c_bulk = density(P[0], T, tw, watersuminitial, c,b);/*���x������*/
	S = crystalviscosity_shear(b,bi, rp, nc, v,r);/*�S���̌����ˑ���*/
	visco = S * vis;/*�S���W��[Pas]  [mol%]*/
	M = machnumber(P[0], Qi, T, tw, watersuminitial, c, b,xc, r);/*�}�b�n��*/

	/*���C*/double Ff =( 8 * visco*v )/ (c_bulk*r*r);

	dPdz[0] = (c_bulk * (g + Ff) )/ (M*M-1);
}
/******************************�C�A���̉^��������**********************************/
void kihoryu(double z,double b[], double P[], double dPdz[], double watersuminitial,double vis,double bi[]) {//�C�A���̉^��������
	double Cw, c_bulk, v, visco, S, M, xc[nc];
	crystalwt(b, xc, c);
	v = velosity(P[0], Qi, T, tw, watersuminitial, c, b, r);/*���x��*/
	Cw = solubilityw(P[0], tw);/*�n��x��*/
	c_bulk = density(P[0], T, tw, watersuminitial, c, b);/*���x������*/
	S = crystalviscosity_shear(b,bi, rp, nc, v, r);/*�S���̌����ˑ���*/
	visco = S * vis;/*�S���W��[Pas]  [mol%]*/
	M = machnumber(P[0], Qi, T, tw, watersuminitial, c, b, xc, r);/*�}�b�n��*/

	/*���C*/double Ff = (8 * visco*v )/ (c_bulk*r*r);

	dPdz[0] = (c_bulk*(g + Ff)) / (M*M - 1);
}
/*********************************�������̉^��������*******************************/
void hunmryu(double z,double b[], double P[], double dPdz[], double watersuminitial,double vis,double bi[]) {
	double Cw, c_bulk, v,M, xc[nc];
	crystalwt(b, xc, c);
	v = velosity(P[0], Qi, T, tw, watersuminitial, c, b, r);/*���x��*/
	M = machnumber(P[0], Qi, T, tw, watersuminitial, c, b, xc, r);/*�}�b�n��*/
	c_bulk = density(P[0], T, tw, watersuminitial, c, b);

	/*���C*/double ff = (0.0025*v*v)/r;

	dPdz[0] = (c_bulk * (g + ff)) / (M*M - 1);
}

/******************************************������****************************/
void crystallization(double z,double P[], double b[], double dbdz[], double watersuminitial) {
	double beq[nc], xeq[nc];
	double	v = velosity(P[0], Qi, T, tw, watersuminitial, c, b, r);/*���x��*/
	crystalequilibrium(P[0], beq, xeq, tpl, tpy, to, tq, ts, tcl, tor, c);


//	if (P[0] / pow(10,8) > 1.5) { beq[0] = 0; };


	dbdz[0] = (G / v)*pow(b[0], 2 / 3)*(beq[0] - b[0]);
	dbdz[1] = (G / v)*pow(b[1], 2 / 3)*(beq[1] - b[1]);
	dbdz[2] = (G / v)*pow(b[2], 2 / 3)*(beq[2] - b[2]);
	dbdz[3] = (G / v)*pow(b[3], 2 / 3)*(beq[3] - b[3]);
	dbdz[4] = (G / v)*pow(b[4], 2 / 3)*(beq[4] - b[4]);
	dbdz[5] = (G / v)*pow(b[5], 2 / 3)*(beq[5] - b[5]);
	dbdz[6] = (G / v)*pow(b[6], 2 / 3)*(beq[6] - b[6]);

	if ( beq[1] - b[1] < 0 ) {
	dbdz[1] = (1 / v)*pow(b[1], 2 / 3)*(beq[1] - b[1]);
	}

/*	if (dbdz[0] < 0) { dbdz[0] = 0; };
	if (dbdz[1] < 0) { dbdz[1] = 0; };
	if (dbdz[2] < 0) { dbdz[2] = 0; };
	if (dbdz[3] < 0) { dbdz[3] = 0; };
	if (dbdz[4] < 0) { dbdz[4] = 0; };
	if (dbdz[5] < 0) { dbdz[5] = 0; };
	if (dbdz[6] < 0) { dbdz[6] = 0; };*/

//	double beq = bi + (1 - bi)*(C0 + C1 * log(P[0] / 1000000) + C2 * log(P[0] / 1000000)*log(P[0] / 1000000));
//	dbdz[0] = (G / v)*pow(b[0], 2 / 3)*(beq - b[0]);
	//�������Ɋւ��Ẵ����Q�N�b�^���s���Ƃ��AP���x�N�g���ő�����Ă邩��A�C���Ɖt���𕪂��čl����Ƃ���P��萔�Ƃ��ē��ꂽ����������������Ȃ�
}

int main(void) {
	double Cw,Ci, v,bt,vis,S,��c;

	int  i = 0;
	int a = (1 / h);
	double z;
	double P[n], dPdz[n], Pout[n], b[nc],bi[nc], dbdz[nc], bout[nc], xc[nc], xeq[nc], beq[nc], xci[nc], comp_melt[11], comp_melt_mol[12];

	/*******��������********/
	z = 0.0;
	P[0] = Pii;
	crystalequilibrium(P[0], beq, xeq, tpl, tpy, to, tq, ts, tcl, tor, c);
	bi[0] = beq[0];
	bi[1] = beq[1];
	bi[2] = beq[2];
	bi[3] = beq[3];
	bi[4] = beq[4];
	bi[5] = beq[5];
	bi[6] = beq[6];
	b[0] = bi[0];
	b[1] = bi[1];
	b[2] = bi[2];
	b[3] = bi[3];
	b[4] = bi[4];
	b[5] = bi[5];
	b[6] = bi[6];
	printf("initial xeq=%f\n", xeq[5]);

	��c = MPF_bimo(rp, b,bi);
	crystalwt(b, xc, c);
	double xm = (1 - xc[0] - xc[1] - xc[2] - xc[3] - xc[4] - xc[5] - xc[6]);

	Ci = solubilityw(P[0], tw);//�������Ă��邶�傤����
	if (solubilityw(P[0], tw) > 3) { Ci = 3; }
	double xgasi = 0;//�}�O�}���܂�ŋC�A���Ȃ�����
	double xwi = Ci / 100;
	double watersuminitial = (1 - xgasi)*xm*xwi + xgasi;//�����n��ʂ�O�a�����āA�C�A���Ȃ��Ɖ��肵���Ƃ�

	for (int i = 0; i < nc; i++) { xci[i] = xc[i]; };
	composiotion_melt(comp_melt_ini, comp_melt, comp_pl, comp_py, comp_o, comp_q, comp_s, comp_cl, comp_or, xc, xci, P[0], tw,watersuminitial);
	wttomol(comp_melt, comp_melt_mol, P[0], tw, xc, watersuminitial);
	vis = viscosity(comp_melt_mol, T);
	double visi = vis;

	/***********************/
	printf("beq[0]=%f,beq[1]=%f,beq[2]=%f,beq[3]=%f,beq[4]=%f,beq[5]=%f\n", beq[0], beq[1], beq[2], beq[3], beq[4], beq[5]);

	FILE *date, *command, *gp;
	const char *date_file, *bat_file;
	char run_command[64];

	/*------------�f�[�^�t�@�C���̍쐬--------------*/

	date_file = "date.txt";
	date = fopen(date_file, "w");

	printf("   z        P1           \n");

	/**********�C�A�̂Ȃ�����**************/
		while (solubilityw(P[0],tw) > Ci && z <= z0) {
			ryuutai(z,b, P, dPdz, watersuminitial,vis,bi);
			crystallization(z, P, b, dbdz, watersuminitial);
			rk4(P, dPdz, n, z, h, Pout,b, watersuminitial,vis, bi, ryuutai);
			rk4c(b, dbdz, nc, n, z, h, P,  bout,watersuminitial, crystallization);

			P[0] = Pout[0];
			b[0] = bout[0];
			b[1] = bout[1];
			b[2] = bout[2];
			b[3] = bout[3];
			b[4] = bout[4];
			b[5] = bout[5];
			b[6] = bout[6];
			bi[0] = b[0];
			bi[1] = b[1];
			bi[2] = b[2];
			bi[3] = b[3];
			bi[4] = b[4];
			bi[5] = b[5];
			bi[6] = b[6];
			z += h;
			crystalwt(b, xc, c);
			crystalequilibrium(P[0], beq, xeq, tpl, tpy, to, tq, ts, tcl, tor, c);
			composiotion_melt(comp_melt_ini, comp_melt, comp_pl, comp_py, comp_o, comp_q, comp_s, comp_cl, comp_or, xc, xci, P[0], tw,watersuminitial);
			wttomol(comp_melt, comp_melt_mol, P[0], tw,xc, watersuminitial);
			vis = viscosity(comp_melt_mol, T);
			bt = b[0] + b[1] + b[2] + b[3] + b[4] + b[5] + b[6];
			i++;
//					printf("z=%f, P=%f,  ��c=%f ,vis=%f S=%f, q=%f bt=%f, bpl=%f\n", z, P[0] / pow(10, 6), ��c, viscosity(comp_melt_mol, T),crystalviscosity_shear(b,bi, rp, nc, velosity(P[0], Qi, T, tw, watersuminitial, c, b, r), r), porosity(P[0], T, c, tw, watersuminitial, xc), velosity(P[0], Qi, T, tw, watersuminitial, c, b, r),bt, b[0]);

			if ((i % a) == 0) {
				��c = MPF_bimo(rp, b, bi);
				Cw = solubilityw(P[0], tw);
				v = velosity(P[0], Qi, T, tw, watersuminitial, c, b, r);
				double M = machnumber(P[0], Qi, T, tw, watersuminitial, c, b, xc, r);
				double S = crystalviscosity_shear(b, bi,rp, nc, v, r);
				double c_bulk = density(P[0], T, tw, watersuminitial, c, b);
				double q = porosity(P[0], T, c, tw, watersuminitial, xc);
				ryuutai(z, b, P, dPdz, watersuminitial, vis, bi);
				double dudz = -(M*M*dPdz[0]) / (c_bulk*velosity(P[0], Qi, T, tw, watersuminitial, c, b, r));
				/*���C*/double Ff = (8 * vis*S*v) / (c_bulk*r*r);
				double nn= 1 - 0.2*(rp)*pow((bt / ��c), 4);
				fprintf(date, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", z, P[0] / 1000000, beq[0], beq[1], beq[2], beq[3], beq[4], beq[5], beq[6], b[0], b[1], b[2], b[3], b[4], b[5], b[6], bt, ��c, v, Cw, q,nn, log10(S*vis), dudz, vis / visi, S, c_bulk, comp_melt[0], comp_melt[1], comp_melt[2], comp_melt[3], comp_melt[4], comp_melt[5], comp_melt[6], comp_melt[7], comp_melt[8], comp_melt[9], comp_melt[10], dPdz[0] / pow(10, 6), c_bulk*v*dudz, c_bulk*Ff, c_bulk*g, M); //         �f�[�^�̏�������
			}
		}
	printf("%f�C�A�����܂�܂����I\n",z);
	/*****************�C�A��****************/
	bt = b[0] + b[1] + b[2] + b[3] + b[4] + b[5] + b[6];

	
	while (porosity(P[0], T, c, tw, watersuminitial, xc) < ��q && z <= z0 && MPF_bimo(rp, b, bi) - 0.0 > bt) {
		kihoryu(z, b, P, dPdz, watersuminitial, vis,bi);
		crystallization(z, P, b, dbdz, watersuminitial);
		rk4(P, dPdz, n, z, h, Pout, b, watersuminitial, vis, bi, kihoryu);
		rk4c(b, dbdz, nc, n, z, h, P, bout, watersuminitial, crystallization);
		P[0] = Pout[0];
		b[0] = bout[0];
		b[1] = bout[1];
		b[2] = bout[2];
		b[3] = bout[3];
		b[4] = bout[4];
		b[5] = bout[5];
		b[6] = bout[6];
		z += h;


		crystalwt(b, xc, c);
		crystalequilibrium(P[0], beq, xeq, tpl, tpy, to, tq, ts, tcl, tor, c);
		composiotion_melt(comp_melt_ini, comp_melt, comp_pl, comp_py, comp_o, comp_q, comp_s, comp_cl, comp_or, xc, xci, P[0], tw, watersuminitial);
		wttomol(comp_melt, comp_melt_mol, P[0], tw, xc, watersuminitial);
		vis = viscosity(comp_melt_mol, T);
		bt = b[0] + b[1] + b[2] + b[3] + b[4] + b[5] + b[6];
//		printf("z=%f, P=%f SiO2wtpl=%f ,vis=%f\n", z, P[0] / pow(10, 6), b[0], viscosity(comp_melt_mol, T) / visi);
//		printf("z=%f, P=%f,  ��c=%f ,vis=%f S=%f, q=%f bt=%f, bpl=%f\n", z, P[0] / pow(10, 6), ��c, viscosity(comp_melt_mol, T), crystalviscosity_shear(b, bi, rp, nc, velosity(P[0], Qi, T, tw, watersuminitial, c, b, r), r), porosity(P[0], T, c, tw, watersuminitial, xc), velosity(P[0], Qi, T, tw, watersuminitial, c, b, r), bt, b[0]);


		i++;
		if ((i % a) == 0) {
			Cw = solubilityw(P[0], tw);
			v = velosity(P[0], Qi, T, tw, watersuminitial, c, b, r);
			double M = machnumber(P[0], Qi, T, tw, watersuminitial, c, b, xc, r);
			double 	S = crystalviscosity_shear(b, bi,rp, nc, v, r);
			double c_bulk = density(P[0], T, tw, watersuminitial, c, b);
			double q = porosity(P[0], T, c, tw, watersuminitial, xc);
			kihoryu(z, b, P, dPdz, watersuminitial, vis, bi);
			double dudz = -(M*M*dPdz[0]) / (c_bulk*velosity(P[0], Qi, T, tw, watersuminitial, c, b, r));
			/*���C*/double Ff = (8 * vis*S*v) / (c_bulk*r*r);
			double nn = 1 - 0.2*(rp)*pow((bt / ��c), 4);
			��c = MPF_bimo(rp, b, bi);
			fprintf(date, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", z, P[0] / 1000000, beq[0], beq[1], beq[2], beq[3], beq[4], beq[5], beq[6], b[0], b[1], b[2], b[3], b[4], b[5], b[6], bt, ��c, v, Cw, q,nn, log10(S*vis), dudz, vis / visi, S, c_bulk, comp_melt[0], comp_melt[1], comp_melt[2], comp_melt[3], comp_melt[4], comp_melt[5], comp_melt[6], comp_melt[7], comp_melt[8], comp_melt[9], comp_melt[10], dPdz[0] / pow(10, 6), c_bulk*v*dudz, c_bulk*Ff, c_bulk*g, M); //         �f�[�^�̏�������
		}
	}
	double bt_fragmentation = b[0] + b[1] + b[2] + b[3] + b[4]+b[5] + b[6];
	double bf_fragmentation = ((b[0] + b[1] + b[2] + b[3] + b[4] + b[5] + b[6] - (bi[0] + bi[1] + bi[2] + bi[3] + bi[4] + bi[5] + bi[6]))) / bt;

	printf("%.1f %.2f�j�� bf_fragmentation=%f MP=%f \n",z,P[0] / 1000000,bf_fragmentation*100, MPF_bimo(rp, b, bi));
	/***************������******************/
	while (machnumber(P[0], Qi, T, tw, watersuminitial, c, b, xc, r) <= 1 && z <= z0) {
		hunmryu(z, b, P, dPdz, watersuminitial, vis, bi);
		crystallization(z, P, b, dbdz, watersuminitial);
		rk4(P, dPdz, n, z, h, Pout, b, watersuminitial, vis, bi, hunmryu);
		rk4c(b, dbdz, nc, n, z, h, P, bout, watersuminitial, crystallization);
		P[0] = Pout[0];

		z += h;

		crystalwt(b, xc, c);
		crystalequilibrium(P[0], beq, xeq, tpl, tpy, to, tq, ts, tcl, tor, c);
		composiotion_melt(comp_melt_ini, comp_melt, comp_pl, comp_py, comp_o, comp_q, comp_s, comp_cl, comp_or, xc, xci, P[0], tw, watersuminitial);
		wttomol(comp_melt, comp_melt_mol, P[0], tw, xc, watersuminitial);
		vis = viscosity(comp_melt_mol, T);
		bt = b[0] + b[1] + b[2] + b[3] + b[4] + b[5] + b[6];
//		printf("z=%f, P=%f SiO2wt=%f ,vis=%f\n", z, P[0] / pow(10, 6), comp_melt[0], viscosity(comp_melt_mol, T) / visi);
	//	printf("z=%f, P=%f,  ��c=%f ,vis=%f, bt=%f, bpl=%f\n", z, P[0] / pow(10, 6), ��c, viscosity(comp_melt_mol, T)*crystalviscosity_shear(b,bi, rp, nc, velosity(P[0], Qi, T, tw, watersuminitial, c, b, r), r),bt, b[0]);

		i++;
		if ((i % a) == 0) {
			Cw = solubilityw(P[0], tw);
			v = velosity(P[0], Qi, T, tw, watersuminitial, c, b, r);
			vis = 0;
			double M = machnumber(P[0], Qi, T, tw, watersuminitial, c, b, xc, r);
			double S = 0;
			double c_bulk = density(P[0], T, tw, watersuminitial, c, b);
			double q = porosity(P[0], T, c, tw, watersuminitial, xc);
			hunmryu(z, b, P, dPdz, watersuminitial, vis, bi);
			double dudz = -(M*M*dPdz[0]) / (c_bulk*velosity(P[0], Qi, T, tw, watersuminitial, c, b, r));
			double _ = 0;
			fprintf(date, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", z, P[0] / 1000000, beq[0], beq[1], beq[2], beq[3], beq[4], beq[5], beq[6], b[0], b[1], b[2], b[3], b[4], b[5], b[6], bt,��c, v, Cw, q,_, _, dudz,_, _, c_bulk, comp_melt[0], comp_melt[1], comp_melt[2], comp_melt[3], comp_melt[4], comp_melt[5], comp_melt[6], comp_melt[7], comp_melt[8], comp_melt[9], comp_melt[10], dPdz[0] / pow(10, 6),_,_,_, M); //         �f�[�^�̏�������
		}
	}
	bt = b[0] + b[1] + b[2] + b[3] + b[4] + b[5] + b[6];
	printf("z=%.0f P[0]=%f bt=%f q=%f M=%.2f",z,P[0]/1000000,bt, porosity(P[0], T, c, tw, watersuminitial,xc), machnumber(P[0], Qi, T, tw, watersuminitial, c ,b,xc, r));

	fclose(date);

	/*-------------------�o�b�`�t�@�C���̍쐬--------------------*/
	bat_file = "command.gp";
	command = fopen(bat_file, "w");

	fprintf(command, "set yrange  0:5000]\n"); // �͈͂̎w��(�ȗ���)
	fprintf(command, "set xrange [0:0.5]\n");
	fprintf(command, "plot '%s' using 3:1 with lines,'' using 4:1 with lines,'' using 5:1 with lines,'' using 6:1 with lines,'' using 7:1 with lines,'' using 8:1 with lines,'' using 9:1 with lines\n", date_file);
	fprintf(command, "pause -1 'Hit any key to close plot window.'");

	fclose(command);
/**/
	/*---------------------gnuplot�̎��s----------------------*/
	sprintf(run_command, "C:/PROGRA~1/gnuplot/bin/wgnuplot.exe \"%s\"", bat_file);
	system(run_command);

	return 0;
}
