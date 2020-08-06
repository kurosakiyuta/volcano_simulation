#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "solubilityH2O.h"
#include "viscosity.h"
#include "/Users/DELL/source/C++/数値計算C++/数値計算C++/nrutil.h"
#include "rk4c.h"
#include "density.h"
#include "velocity.h"
#include "porosity.h"
#include "speed of sound.h"
#include "Exsolvedgas.h"
#include "wt to mol.h"
#include "Mach number.h"
#include "rk4(crystallization).h"
#include "βeq.h"
#include "crystalviscosity.h"
#include "specific entropy(gas etc).h"
#include "specific entropy(melt etc).h"
#include "specific internal energy(gas etc).h"
#include "specific internal energy(melt etc).h"
#include "density(melt etc).h"
#include "density(liquid).h"
#include "density(gas).h"
#include  "entropy(liquid).h"
#include "rk4(temperature).h"
#include "energy(Liquid).h"
#include "rk4(Solubility).h"

#pragma warning(disable : 4996)


#define n 1 /*要素数（方程式数）*/
#define nc 3 /*要素数（方程式数）*/
#define nT 1 /*要素数（方程式数）*/
#define nd 1
#define g 9.80665

#define SiO2 48.89      //Metrich et al (2001) ST133s
#define TiO2 0.98
#define Al2O3 18.05
#define FeO 3.86 +4.90       
#define MnO  0.16
#define P2O5 0.45
#define MgO 5.97
#define CaO 10.78
#define Na2O 2.6
#define F2O 0
#define K2O 2.07
#define SiO2mol 47.55
#define TiO2mol 0.72
#define Al2O3mol 10.33
#define FeOmol 4.21         /*mol%JR1*/
#define MnOmol 0.12
#define P2O5mol 0.18
#define MgOmol 8.64
#define CaOmol 11.23
#define Na2Omol 2.76
#define F2Omol 0
#define K2Omol 1.29
#define summolfruction 1713.704616 

#define H 0.001
#define C 0.012
#define O 0.016
#define F 0.019        /*原子の分子量kg/mol*/
#define Na 0.023
#define Mg 0.0243
#define Al 0.027
#define Si 0.028
#define RIN 0.031
#define K 0.039
#define Ca 0.04
#define Ti 0.047867
#define Mn 0.055
#define Fe 0.055845

#define R 461.88888
#define Mw 0.01802
//#define cpl 2650  //斜長石の密度 high albite Rock forming minerals 
//#define cpy 3585  //輝石の密度 Enstatite-ferrosilite の中間地
//#define co 3700  //オリビンの密度
//#define cb 3200  //結晶の密度（適当
/******変数定数******/
#define r 0.8 //radius
//#define cl 2400 //melt density
#define Ci 4.0 //水の初期溶解量

#define Tii 1050+273.65  //Kelvin
#define A -4.55       //粘性係数の定数
#define AI 0.02+0.02-0.09 //cation mol fruction　Na+K-Al
/************/
#define Qi 70000  //多分,体積流量kg/m~2s
#define z0 10000

#define C0 0.842
#define C1 -0.0954
#define C2 -0.0116 /*結晶化のパラメータ*/
#define bpli 0.0
#define bpyi 0.20
#define boi 0.05
#define bi 0.1

#define G pow(10,1)
#define Gs 0.5
#define φx 0.35
#define φq 0.8
//                      0=melt,            1=disH2O,        2= pl,           3=py,        4=  olivine ,         5=exs H2O
double P0[6] = { 2.5*pow(10,8) ,	  pow(10,8)  ,	2.5*pow(10,8) ,	2.5*pow(10,8) ,	2.5*pow(10,8) };
double C0k[6] = { 1500          ,            407 ,            2000,          2000 ,          2000 };
double c0k[6] = { 2650          ,           1000 ,            2800,          3250 ,          3300 };
double Y[6] = { 2.09          ,           1.11 ,             3.4,           3.4 ,           3.4 };
double cv[6] = { 707          ,           3637 ,            360 ,           345 ,           362 ,                 1571};
double ek[6] = { 0          ,              0 ,              0 ,             0 ,             0 ,                      0};
double sk[6] = { 0          ,            1235,              0 ,             0 ,             0 };
double T0[6] = { 298.65      ,          298.65,         298.65 ,        298.65 ,        298.65 };

/************************************方程式の定義**********************************/
/*********************************気泡のない運動方程式*****************************/
void ryuutai(double z,double xd[],double xc[], double T[], double P[], double dPdz[]) {
	double Cw, c, v, vis, Cwmol,S;
	double X = exsolvegas(P[0], T[0], Ci, xd);
	double cm = densitymelt(P[0], T[0], cv[0], c0k[0], C0k[0], Y[0], P0[0]);
	double cd = densitymelt(P[0], T[0], cv[1], c0k[1], C0k[1], Y[1], P0[1]);
	double cpl = densitymelt(P[0], T[0], cv[2], c0k[2], C0k[2], Y[2], P0[2]);
	double cpy = densitymelt(P[0], T[0], cv[3], c0k[3], C0k[3], Y[3], P0[3]);
	double co = densitymelt(P[0], T[0], cv[4], c0k[4], C0k[4], Y[4], P0[4]);
	double cl = densityliquid(P[0], xc, xd, cpl, cpy, co, cm, cd);
	double cg = densitygas(P[0], T[0]);
	/*速度式*/v = velosity(P[0], Qi, T[0], xd, Ci, cl, cg, r);
	/*溶解度則*/Cw = solubilityw(P[0], T[0], AI);
	/*密度方程式*/c = density(P[0], T[0], xd, Ci, cl, cg);
	/*溶解した水(mol%)*/Cwmol = (Cw * 100 / (H * 2 + O)) / (summolfruction);
	/*粘性の結晶依存性*/S = crystalviscosity(xc, cpl, cpy, co, cl);
	/*粘性係数[Pas]  [mol%]*/vis = S * viscosity(SiO2mol, TiO2mol, Al2O3mol, FeOmol, MnOmol, P2O5mol, MgOmol, CaOmol, Na2Omol, Cwmol, F2Omol, K2Omol, T[0], A);

	/*摩擦*/double Ff = 8 * vis*v / (c*v*v);



	dPdz[0] = -c * (g + Ff);

}
/******************************気泡流の運動方程式**********************************/
void kihoryu(double z,double xd[],double xc[],double T[], double P[], double dPdz[]) {//気泡流の運動方程式
	double Cw, c, v, vis,M,Cwmol,S,X;
	/*揮発した揮発性成分【wt％】*/X = exsolvegas(P[0], T[0], Ci, xd);
	double cm = densitymelt(P[0], T[0], cv[0], c0k[0], C0k[0], Y[0], P0[0]);
	double cd = densitymelt(P[0], T[0], cv[1], c0k[1], C0k[1], Y[1], P0[1]);
	double cpl = densitymelt(P[0], T[0], cv[2], c0k[2], C0k[2], Y[2], P0[2]);
	double cpy = densitymelt(P[0], T[0], cv[3], c0k[3], C0k[3], Y[3], P0[3]);
	double co = densitymelt(P[0], T[0], cv[4], c0k[4], C0k[4], Y[4], P0[4]);
	double cl = densityliquid(P[0], xc, xd, cpl, cpy, co, cm, cd);
	double cg = densitygas(P[0], T[0]);
	/*速度式*/v = velosity(P[0], Qi,T[0], xd, Ci, cl, cg, r);
	/*溶解度則*/Cw = solubilityw(P[0], T[0], AI);
	/*密度方程式*/c = density(P[0], T[0], xd, Ci, cl, cg);
	/*溶解した水(mol%)*/Cwmol = (Cw * 100 / (H * 2 + O)) / (summolfruction);
	/*粘性の結晶依存性*/S = crystalviscosity(xc, cpl, cpy, co, cl);
	/*粘性係数[Pas]  [mol%]*/vis = S*viscosity(SiO2mol, TiO2mol, Al2O3mol, FeOmol, MnOmol, P2O5mol, MgOmol, CaOmol, Na2Omol, Cwmol, F2Omol, K2Omol, T[0], A);
	/*マッハ数*/M = machnumber(P[0], Qi, T[0], xd, Ci, cl, cg, r);

	/*摩擦*/double Ff = 8 * vis*v / (c*v*v);

	dPdz[0] = c*(g+Ff)/(M*M-1) ;

}
/*********************************噴霧流の運動方程式*******************************/
void hunmryu(double z,double xd[],double xc[],double T[], double P[], double dPdz[]) {
	double Cw, c, v, M,X;
	/*揮発した揮発性成分【wt％】*/X = exsolvegas(P[0], T[0], Ci, xd);
	double cm = densitymelt(P[0], T[0], cv[0], c0k[0], C0k[0], Y[0], P0[0]);
	double cd = densitymelt(P[0], T[0], cv[1], c0k[1], C0k[1], Y[1], P0[1]);
	double cpl = densitymelt(P[0], T[0],cv[2], c0k[2], C0k[2], Y[2], P0[2]);
	double cpy = densitymelt(P[0], T[0], cv[3], c0k[3], C0k[3], Y[3], P0[3]);
	double co = densitymelt(P[0], T[0], cv[4], c0k[4], C0k[4], Y[4], P0[4]);
	double cl = densityliquid(P[0], xc, xd, cpl, cpy, co, cm, cd);
	double cg = densitygas(P[0], T[0]);
	/*溶解度則*/Cw = solubilityw(P[0], T[0], AI);
	/*密度方程式*/c = density(P[0], T[0], xd, Ci, cl, cg);
	/*速度式*/v = velosity(P[0], Qi, T[0], xd, Ci, cl, cg, r);
	/*マッハ数*/M = machnumber(P[0], Qi, T[0], xd, Ci, cl, cg, r);
	/*摩擦*/double ff = (0.0025*v*v)/r;

	dPdz[0] = c * (g + ff) / (M*M - 1);
}
/******************************************結晶化**********************************/
void crystallization(double z,double P[],double xd[], double T[],double xc[], double dxcdz[]) {
	double X = exsolvegas(P[0], T[0], Ci, xd);
	double cm = densitymelt(P[0], T[0], cv[0], c0k[0], C0k[0], Y[0], P0[0]);
	double cd = densitymelt(P[0], T[0], cv[1], c0k[1], C0k[1], Y[1], P0[1]);
	double cpl = densitymelt(P[0], T[0], cv[2], c0k[2], C0k[2], Y[2], P0[2]);
	double cpy = densitymelt(P[0], T[0], cv[3], c0k[3], C0k[3], Y[3], P0[3]);
	double co = densitymelt(P[0], T[0], cv[4], c0k[4], C0k[4], Y[4], P0[4]);
	double cl = densityliquid(P[0], xc, xd, cpl, cpy, co, cm, cd);
	double cg = densitygas(P[0], T[0]);
	double Cw = solubilityw(P[0], T[0], AI);
	double tpl[10] = {/*1*/-3.0938*pow(10,-9),/*2*/-8.0678*pow(10,-6),/*3*/2.2949*pow(10,-3),/*4*/-6.9366*pow(10,-8),/*5*/-6.6318*pow(10,-4),/*6*/-1.0765*pow(10,-5),/*7*/1.1705*pow(10,-4),/*8*/1.4834*pow(10,-2),/*9*/5.3169*pow(10,-1),/*10*/-6.1889*pow(10,-0) };//MELTsからのパラメータ
	double tpy[10] ={/*1*/-3.2977*pow(10,-9),/*2*/-1.8604*pow(10,-5),/*3*/-6.0187*pow(10,-3),/*4*/4.7402*pow(10,-7),/*5*/-9.4442*pow(10,-4),/*6*/1.1645*pow(10,-5),/*7*/-4.7521*pow(10,-4),/*8*/3.9398*pow(10,-2),/*9*/9.8480*pow(10,-1),/*10*/-2.0651*pow(10,1) };//MELTsからのパラメータ
	double to[10] = {/*1*/-1.5580*pow(10,-9),/*2*/2.8401*pow(10,-6),/*3*/6.4915*pow(10,-3),/*4*/-3.9540*pow(10,-8),/*5*/2.7014*pow(10,-4),/*6*/-2.2924*pow(10,-6),/*7*/5.0093*pow(10,-5),/*8*/-7.3693*pow(10,-3),/*9*/-3.3874*pow(10,-1),/*10*/4.7447*pow(10,0) };//MELTsからのパラメータ
	double xcpleq = xcequilibrium(T[0],  P[0], Cw, Ci, cpl, tpl, AI, cl,xd[0]);
	double xcpyeq = xcequilibrium(T[0], P[0], Cw, Ci, cpy, tpy, AI, cl, xd[0]);
	double xcoeq = xcequilibrium(T[0], P[0], Cw, Ci, co, to, AI, cl, xd[0]);
	double v = velosity(P[0], Qi, T[0], xd, Ci, cl, cg, r);
	double  q = porosity(P[0], T[0], cl, xd, Ci);
	double xdeq = Cw/100;
	if (solubilityw(P[0], T[0], AI) > Ci) { xdeq = Ci/100; };

	double  φg = (1 / Gs)*(xd[0] - xdeq)*(1 - q)*cl*(1 - xc[0] - xc[1] - xc[2]);

	dxcdz[0] = -(xc[0] - xcpleq) / (v*G) + (xc[0] / ((1 - q)*cl*v))*φg;
	dxcdz[1] = -(xc[1] - xcpyeq)/ (v*G) + (xc[1] / ((1 - q)*cl*v))*φg;
	dxcdz[2] = -(1 / (v*G))*(xc[2] - xcoeq) + (xc[2] / ((1 - q)*cl*v))*φg;

//	printf("xpy=%f dxpydz=%f xceq=%f\n", xc[1], dxcdz[1],xcpyeq);
//	double beq = bi + (1 - bi)*(C0 + C1 * log(P[0] / 1000000) + C2 * log(P[0] / 1000000)*log(P[0] / 1000000));
//	dbdz[0] = (G / v)*pow(b[0], 2 / 3)*(beq - b[0]);
	//結晶化に関してのルンゲクッタを行うとき、Pをベクトルで代入してるから、気相と液相を分けて考えるときにPを定数として入れた方がいいかもしれない
}
/**********************************solubility**************************************/
void solubility(double z, double P[], double T[], double xc[], double xd[], double dxddz[]) {
	double X = exsolvegas(P[0], T[0], Ci, xd);
	double cm = densitymelt(P[0], T[0], cv[0], c0k[0], C0k[0], Y[0], P0[0]);
	double cd = densitymelt(P[0], T[0], cv[1], c0k[1], C0k[1], Y[1], P0[1]);
	double cpl = densitymelt(P[0], T[0], cv[2], c0k[2], C0k[2], Y[2], P0[2]);
	double cpy = densitymelt(P[0], T[0], cv[3], c0k[3], C0k[3], Y[3], P0[3]);
	double co = densitymelt(P[0], T[0], cv[4], c0k[4], C0k[4], Y[4], P0[4]);
	double cl = densityliquid(P[0], xc, xd, cpl, cpy, co, cm, cd);
	double cg = densitygas(P[0], T[0]);
	double  q = porosity(P[0], T[0], cl, xd, Ci);
	double v = velosity(P[0], Qi, T[0], xd, Ci, cl, cg, r);

	double xdeq= solubilityw(P[0], T[0], AI)/100;
	if (solubilityw(P[0], T[0], AI) > Ci) { xdeq = Ci/100; };

	double  φg = (1 / Gs)*(xd[0] - xdeq)*(1 - q)*cl*(1 - xc[0] - xc[1] - xc[2]);

	dxddz[0] = -((1 - xd[0])*φg) / (cl*(1 - q)*v);
}
/**********************************温度変化****************************************/
void Temperature (double z,double xc[],double xd[],double dxcdz[],double dxddz[],double P[] ,double dPdz[], double T[],double dTdz[] ){
	double sm, sd, scpl, scpy, sco, sl, sg, el, eg, em, ed, ecpl, ecpy, eco,X,cm,cd,cpl,cpy,co,cl,cg;
	/*揮発した揮発性成分【wt％】*/X = exsolvegas(P[0], T[0], Ci, xd);
	double xm = 1 - xd[0] - xc[0] - xc[1] - xc[2];
	cm = densitymelt(P[0], T[0], cv[0],c0k[0],C0k[0],Y[0],P0[0]);
	cd = densitymelt(P[0], T[0], cv[1], c0k[1], C0k[1], Y[1], P0[1]);
	cpl = densitymelt(P[0], T[0], cv[2], c0k[2], C0k[2], Y[2], P0[2]);
	cpy = densitymelt(P[0], T[0], cv[3], c0k[3], C0k[3], Y[3], P0[3]);
	co = densitymelt(P[0], T[0], cv[4], c0k[4], C0k[4], Y[4], P0[4]);
	cl = densityliquid(P[0], xc, xd, cpl, cpy, co, cm, cd);

	cg = densitygas(P[0], T[0]);

	sm = entropymelt(P[0], T[0], sk[0], cv[0], c0k[0], T0[0], Y[0], cm);
	sd = entropymelt(P[0], T[0], sk[1], cv[1], c0k[1], T0[1], Y[1], cd);
	scpl = entropymelt(P[0], T[0], sk[2], cv[2], c0k[2], T0[2], Y[2], cpl);
	scpy = entropymelt(P[0], T[0], sk[3], cv[3], c0k[3], T0[3], Y[3], cpy);
	sco = entropymelt(P[0], T[0], sk[4], cv[4], c0k[4], T0[4], Y[4], co);
	sl = entropyliquid(sm,scpl,scpy, sco, sd,xc,xd, cpl, cpy,co,cl);

	sg = entropygas(P[0], T[0], cg);

	em = energymelt(P[0], T[0], ek[0], cv[0], c0k[0], C0k[0], Y[0], P0[0], cm);
	ed = energymelt(P[0], T[0], ek[1], cv[1], c0k[1], C0k[1], Y[1], P0[1], cd);
	ecpl = energymelt(P[0], T[0], ek[2], cv[2], c0k[2], C0k[2], Y[2], P0[2], cpl);
	ecpy = energymelt(P[0], T[0], ek[3], cv[3], c0k[3], C0k[3], Y[3], P0[3], cpy);
	eco = energymelt(P[0], T[0], ek[4], cv[4], c0k[4], C0k[4], Y[4], P0[4], co);
	el = energyliquid(xd, xc, em, ed, ecpl, ecpy, eco);

	eg = energygas(T[0], ek[5], cv[5]);

	double/*溶解度則*/Cw = solubilityw(P[0], T[0], AI);
	double/*密度方程式*/c = density(P[0], T[0], xd, Ci, cl, cg);
	double/*溶解した水(mol%)*/Cwmol = (Cw * 100 / (H * 2 + O)) / (summolfruction);
	double/*粘性の結晶依存性*/S = crystalviscosity(xc, cpl, cpy, co, cl);
	double/*粘性係数[Pas]  [mol%]*/vis = S * viscosity(SiO2mol, TiO2mol, Al2O3mol, FeOmol, MnOmol, P2O5mol, MgOmol, CaOmol, Na2Omol, Cwmol, F2Omol, K2Omol, T[0], A);
	double/*速度式*/v = velosity(P[0], Qi, T[0], xd, Ci, cl, cg, r);
	double  q = porosity(P[0], T[0], cl, xd, Ci);
	double M = machnumber(P[0], Qi, T[0], xd, Ci, cl, cg, r);

	double tpl[10] = {/*1*/-3.0938*pow(10,-9),/*2*/-8.0678*pow(10,-6),/*3*/2.2949*pow(10,-3),/*4*/-6.9366*pow(10,-8),/*5*/-6.6318*pow(10,-4),/*6*/-1.0765*pow(10,-5),/*7*/1.1705*pow(10,-4),/*8*/1.4834*pow(10,-2),/*9*/5.3169*pow(10,-1),/*10*/-6.1889*pow(10,-0) };//MELTsからのパラメータ
	double tpy[10] = {/*1*/-3.2977*pow(10,-9),/*2*/-1.8604*pow(10,-5),/*3*/-6.0187*pow(10,-3),/*4*/4.7402*pow(10,-7),/*5*/-9.4442*pow(10,-4),/*6*/1.1645*pow(10,-5),/*7*/-4.7521*pow(10,-4),/*8*/3.9398*pow(10,-2),/*9*/9.8480*pow(10,-1),/*10*/-2.0651*pow(10,1) };//MELTsからのパラメータ
	double to[10] = {/*1*/-1.5580*pow(10,-9),/*2*/2.8401*pow(10,-6),/*3*/6.4915*pow(10,-3),/*4*/-3.9540*pow(10,-8),/*5*/2.7014*pow(10,-4),/*6*/-2.2924*pow(10,-6),/*7*/5.0093*pow(10,-5),/*8*/-7.3693*pow(10,-3),/*9*/-3.3874*pow(10,-1),/*10*/4.7447*pow(10,0) };//MELTsからのパラメータ
	double xcpleq = xcequilibrium(T[0], P[0], Cw, Ci, cpl, tpl, AI, cl, xd[0]);
	double xcpyeq = xcequilibrium(T[0], P[0], Cw, Ci, cpy, tpy, AI, cl, xd[0]);
	double xcoeq = xcequilibrium(T[0], P[0], Cw, Ci, co, to, AI, cl, xd[0]);

	double xdeq = Cw/100;
	if (solubilityw(P[0], T[0], AI) > Ci) { xdeq = Ci/100; };

	double  φg = (1 / Gs)*(xd[0] - xdeq)*(1 - q)*cl*(1 - xc[0] - xc[1] - xc[2]);

	double ζ[5];
	ζ[0] = (c0k[0] * C0k[0] * C0k[0] - Y[0] * P0[0]) / Y[0];//melt
	ζ[1] = (c0k[1] * C0k[1] * C0k[1] - Y[1] * P0[1]) / Y[1];//dissolved gas
	ζ[2] = (c0k[2] * C0k[2] * C0k[2] - Y[2] * P0[2]) / Y[2];//plgioclase
	ζ[3] = (c0k[3] * C0k[3] * C0k[3] - Y[3] * P0[3]) / Y[3];//pyroxene
	ζ[4] = (c0k[4] * C0k[4] * C0k[4] - Y[4] * P0[4]) / Y[4];//oliven

	dxddz[0] = -((1 - xd[0])*φg) / (cl*(1 - q)*v);
	dxcdz[0] = -(1 / (v*G))*(xc[0] - xcpleq) + (xc[0] /( (1 - q)*cl*v))*φg;
	dxcdz[1] = -(1 / (v*G))*(xc[1] - xcpyeq) + (xc[1] / ((1 - q)*cl*v))*φg;
	dxcdz[2] = -(1 / (v*G))*(xc[2] - xcoeq) + (xc[2] / ((1 - q)*cl*v))*φg;
	double Aa = xm * (cv[0] + ζ[0] / (cm*T[0])) + xd[0] * (cv[1] + ζ[1] / (cd*T[0])) + xc[0] * (cv[2] + ζ[2] / (cpl*T[0])) + xc[1] * (cv[3] + ζ[3] / (cpy*T[0])) + xc[2] * (cv[4] + ζ[4] / (co*T[0]));
	double N = (1 - q)*cl*Aa + q * cg*cv[5];

	double B = ζ[0] * xm / (cm*cm*cv[0] * (Y[0] - 1)*T[0]) + ζ[1] * xd[0] / (cd*cd*cv[1] * (Y[1] - 1)*T[0]) + ζ[2] * xc[0] / (cpl*cpl*cv[2] * (Y[2] - 1)*T[0]) + ζ[3] * xc[1] / (cpy*cpy*cv[3] * (Y[3] - 1)*T[0]) + ζ[4] * xc[2] / (co*co*cv[4] * (Y[4] - 1)*T[0]);
	double M1 = φg * (el - eg);
	double M2 = (1 - q)*cl*(( - dxddz[0] - dxcdz[0] - dxcdz[1] - dxcdz[2])*em + dxddz[0] * ed + dxcdz[0] * ecpl + dxcdz[1] * ecpy + dxcdz[2] * eco);
	double M3 = ((1 - q)*cl*B + (P[0]*M*M) / (c*v))*dPdz[0];

	dTdz[0] = (M1 - M2 + M3) / N;

//	printf("N=%.0f M1=%.0f M2=%.0f M3=%.0f dxddz=%.2f  dxcpldz=%.2f dxcpydz=%.2f dxcodz=%.2f  φ=%.1f　xdeq=%f  dTdz=%f dPdz=%f\n", N, M1, M2, M3, dxddz[0],dxcdz[0],dxcdz[1],dxcdz[2],φg,xdeq,dTdz[0],dPdz[0]);
//	printf("xcpleq=%f xcpyeq=%f xcoeq=%f \n", xcpleq, xcpyeq, xcoeq);

}
int main(void) {
	double Cw, v,X,cpl,cpy,co,cl,cg,cd,cm,c;
	double h, z;
	double P[n], dPdz[n], Pout[n], xc[nc], dxcdz[nc], xcout[nc],dTdz[nT],T[nT],Tout[nT],xd[nd],dxddz[nd],xdout[nd],b[nc];
	/*******初期条件********/
	z = 0.0;
	P[0] = 250 * 1000000;
	T[0] = Tii;
	Cw = solubilityw(P[0], T[0], AI);

	 cm = densitymelt(P[0], T[0], cv[0], c0k[0], C0k[0], Y[0], P0[0]);
	 cd = densitymelt(P[0], T[0], cv[1], c0k[1], C0k[1], Y[1], P0[1]);
	 cpl = densitymelt(P[0], T[0], cv[2], c0k[2], C0k[2], Y[2], P0[2]);
	 cpy = densitymelt(P[0], T[0], cv[3], c0k[3], C0k[3], Y[3], P0[3]);
	 co = densitymelt(P[0], T[0], cv[4], c0k[4], C0k[4], Y[4], P0[4]);

	cl = densityliquid(P[0], xc, xd, cpl, cpy, co, cm, cd);
	 cg = densitygas(P[0], T[0]);

	xd[0] = Ci/100;
	if (Cw < Ci) { xd[0] = Cw/100; };

	double tpl[10] = {/*1*/-3.0938*pow(10,-9),/*2*/-8.0678*pow(10,-6),/*3*/2.2949*pow(10,-3),/*4*/-6.9366*pow(10,-8),/*5*/-6.6318*pow(10,-4),/*6*/-1.0765*pow(10,-5),/*7*/1.1705*pow(10,-4),/*8*/1.4834*pow(10,-2),/*9*/5.3169*pow(10,-1),/*10*/-6.1889*pow(10,-0) };//MELTsからのパラメータ
	double tpy[10] = {/*1*/-3.2977*pow(10,-9),/*2*/-1.8604*pow(10,-5),/*3*/-6.0187*pow(10,-3),/*4*/4.7402*pow(10,-7),/*5*/-9.4442*pow(10,-4),/*6*/1.1645*pow(10,-5),/*7*/-4.7521*pow(10,-4),/*8*/3.9398*pow(10,-2),/*9*/9.8480*pow(10,-1),/*10*/-2.0651*pow(10,1) };//MELTsからのパラメータ
	double to[10]  = {/*1*/-1.5580*pow(10,-9),/*2*/2.8401*pow(10,-6),/*3*/6.4915*pow(10,-3),/*4*/-3.9540*pow(10,-8),/*5*/2.7014*pow(10,-4),/*6*/-2.2924*pow(10,-6),/*7*/5.0093*pow(10,-5),/*8*/-7.3693*pow(10,-3),/*9*/-3.3874*pow(10,-1),/*10*/4.7447*pow(10,0) };//MELTsからのパラメータ
	double xcpleq = xcequilibrium(T[0], P[0], Cw, Ci, cpl, tpl, AI, cl, xd[0]);
	double xcpyeq = xcequilibrium(T[0], P[0], Cw, Ci, cpy, tpy, AI, cl, xd[0]);
	double xcoeq = xcequilibrium(T[0], P[0], Cw, Ci, co, to, AI, cl, xd[0]);

	xc[0] = xcpleq;
	xc[1] = xcpyeq;
	xc[2] = xcoeq;

	/***********************/



	double em = energymelt(P[0], T[0], ek[0], cv[0], c0k[0], C0k[0], Y[0], P0[0], cm);
	double ed = energymelt(P[0], T[0], ek[1], cv[1], c0k[1], C0k[1], Y[1], P0[1], cd);
	double ecpl = energymelt(P[0], T[0], ek[2], cv[2], c0k[2], C0k[2], Y[2], P0[2], cpl);
	double ecpy = energymelt(P[0], T[0], ek[3], cv[3], c0k[3], C0k[3], Y[3], P0[3], cpy);
	double eco = energymelt(P[0], T[0], ek[4], cv[4], c0k[4], C0k[4], Y[4], P0[4], co);
	double el = energyliquid(xd, xc, em, ed, ecpl, ecpy, eco);
	double eg = energygas(T[0], ek[5], cv[5]);

	 v = velosity(P[0], Qi, T[0], xd, Ci, cl, cg, r);
	double sm = entropymelt(P[0], T[0], sk[0], cv[0], c0k[0], T0[0], Y[0], cm);
	double sd = entropymelt(P[0], T[0], sk[1], cv[1], c0k[1], T0[1], Y[1], cd);
	double scpl = entropymelt(P[0], T[0], sk[2], cv[2], c0k[2], T0[2], Y[2], cpl);
	double scpy = entropymelt(P[0], T[0], sk[3], cv[3], c0k[3], T0[3], Y[3], cpy);
	double sco = entropymelt(P[0], T[0], sk[4], cv[4], c0k[4], T0[4], Y[4], co);
	double sl = entropyliquid(sm, scpl, scpy, sco, sd, xc, xd, cpl, cpy, co, cl);
	double sg = entropygas(P[0], T[0], cg);
	/*溶解度則*/Cw = solubilityw(P[0], T[0], AI);
	/*密度方程式*/c = density(P[0], T[0], xd, Ci, cl, cg);
	double/*溶解した水(mol%)*/Cwmol = (Cw * 100 / (H * 2 + O)) / (summolfruction);
	double/*粘性の結晶依存性*/S = crystalviscosity(xd, cpl, cpy, co, cl);
	double/*粘性係数[Pas]  [mol%]*/vis = S * viscosity(SiO2mol, TiO2mol, Al2O3mol, FeOmol, MnOmol, P2O5mol, MgOmol, CaOmol, Na2Omol, Cwmol, F2Omol, K2Omol, T[0], A);
	/*速度式*/v = velosity(P[0], Qi, T[0], xd, Ci, cl, cg, r);
	double  q = porosity(P[0], T[0], cl, xd, Ci);
	double xdeq = solubilityw(P[0], T[0], AI) / 100;
	if (solubilityw(P[0], T[0], AI) > Ci) { xdeq = Ci / 100; };

	double  φg = (1 / Gs)*(xd[0] - xdeq)*(1 - q)*cl*(1 - xc[0] - xc[1] - xc[2]);
	dxcdz[1] = -(1 / (v*G))*(xc[1] - xcpyeq) + (xc[1] / ((1 - q)*cl*v))*φg;
//	printf("%.2f %.1f T=%.0f cl=%.0f sl=%.0f xd=%.2f xpl=%.2f xpy=%.3f xo=%.3f  cm=%.4f cd=%.4f cpl=%.4f cpy=%.4f co=%.4f em=%f epl=%fel=%f eg=%f \n", z, P[0] / 1000000, T[0],cl,sl,xd[0],xc[0],xc[1],xc[2],cm,cd,cpl,cpy,co,em,ecpl,el,eg);
//	printf("%.2f %.1f T=%.0f dtdz=%f cl=%.0f xc[1]=%f xcpyeq=%f dxcdz[1]=%f dxcdz[左]=%f dxcdz[右]=%f  \n", z, P[0] / 1000000, T[0], dTdz[0], cl, xc[1], xcpyeq, dxcdz[1], -(xc[1] - xcpyeq) / (v*G), (xc[1] / ((1 - q)*cl*v))*φg);
	FILE *date, *command, *gp;
	const char *date_file, *bat_file;
	char run_command[64];

	/*------------データファイルの作成--------------*/

	date_file = "date.txt";
	date = fopen(date_file, "w");

	h = 0.05;

	printf("   z        P1           \n");

	/**********気泡のない流れ**************/
//	X = exsolvegas(P[0], T[0], Ci, AI);
//	cm = densitymelt(P[0], T[0], cv[0], c0k[0], C0k[0], Y[0], P0[0]);
//	cd = densitymelt(P[0], T[0], cv[1], c0k[1], C0k[1], Y[1], P0[1]);
//	cpl = densitymelt(P[0], T[0], cv[2], c0k[2], C0k[2], Y[2], P0[2]);
//	cpy = densitymelt(P[0], T[0], cv[3], c0k[3], C0k[3], Y[3], P0[3]);
//	co = densitymelt(P[0], T[0], cv[4], c0k[4], C0k[4], Y[4], P0[4]);
//	cl = densityliquid(P[0], b[0], b[1], b[2], X, Ci, cpl, cpy, co, cm, cd);
//	cg = densitygas(P[0], T[0]);
//	c = density(P[0], T[0], AI, Ci, cl, cg);
//	double q = porosity(P[0], T[0], cl, AI, Ci);
//	printf("%.2f %.1f  c=%f cm=%f ad=%f, cpl=%f, cpy=%f, co=%f, cg=%f,cl=%f q=%f\n", z, P[0] / 1000000, c, cm, cd, cpl, cpy, co, cg,cl, q);


	if (solubilityw(P[0], T[0], AI) >= Ci) {/***溶解度が初期溶解量を超えていれば、気泡は生まれない***/
		Cw = Ci;
		for (int i = 0; solubilityw(P[0], T[0], AI) > Ci - 0.01&& z <= z0; i++) {
			crystallization(z, P,xd,T, xc, dxcdz);
			rk4c(xc, dxcdz, nc, n, z, h, P, xcout,T, xd,crystallization);
			ryuutai(z,xd,xc,T, P, dPdz);
			Temperature(z, xc,xd,dxcdz,dxddz, P,dPdz, T, dTdz);
			solubility(z, P, T, xc, xd, dxddz);
			rk4(P, dPdz, n, z, h, Pout, xc,T,xd,xc, ryuutai);
			rk4T( T,dTdz,nT,z, h,Tout,xc, P, dPdz,xc,xd,dxcdz,dxddz, Temperature);
			rk4S(xd, dxddz, nd, z, h, xdout, xc, T, xd, xc, P, solubility);
			P[0] = Pout[0];
			xc[0] = xcout[0];
			xc[1] = xcout[1];
			xc[2] = xcout[2];
			xd[0] = xdout[0];
			T[0] = Tout[0];
			z += h;
			double X = exsolvegas(P[0], T[0], Ci, xd);
			 cm = densitymelt(P[0], T[0], cv[0], c0k[0], C0k[0], Y[0], P0[0]);
			 cd = densitymelt(P[0], T[0], cv[1], c0k[1], C0k[1], Y[1], P0[1]);
			 cpl = densitymelt(P[0], T[0], cv[2], c0k[2], C0k[2], Y[2], P0[2]);
			 cpy = densitymelt(P[0], T[0], cv[3], c0k[3], C0k[3], Y[3], P0[3]);
			 co = densitymelt(P[0], T[0], cv[4], c0k[4], C0k[4], Y[4], P0[4]);
			 cl = densityliquid(P[0], xc, xd, cpl, cpy, co, cm, cd);
			 cg = densitygas(P[0], T[0]);


			double Cw = solubilityw(P[0], T[0], AI);
			double tpl[10] = {/*1*/-3.0938*pow(10,-9),/*2*/-8.0678*pow(10,-6),/*3*/2.2949*pow(10,-3),/*4*/-6.9366*pow(10,-8),/*5*/-6.6318*pow(10,-4),/*6*/-1.0765*pow(10,-5),/*7*/1.1705*pow(10,-4),/*8*/1.4834*pow(10,-2),/*9*/5.3169*pow(10,-1),/*10*/-6.1889*pow(10,-0) };//MELTsからのパラメータ
			double tpy[10] = {/*1*/-3.2977*pow(10,-9),/*2*/-1.8604*pow(10,-5),/*3*/-6.0187*pow(10,-3),/*4*/4.7402*pow(10,-7),/*5*/-9.4442*pow(10,-4),/*6*/1.1645*pow(10,-5),/*7*/-4.7521*pow(10,-4),/*8*/3.9398*pow(10,-2),/*9*/9.8480*pow(10,-1),/*10*/-2.0651*pow(10,1) };//MELTsからのパラメータ
			double to[10] = {/*1*/-1.5580*pow(10,-9),/*2*/2.8401*pow(10,-6),/*3*/6.4915*pow(10,-3),/*4*/-3.9540*pow(10,-8),/*5*/2.7014*pow(10,-4),/*6*/-2.2924*pow(10,-6),/*7*/5.0093*pow(10,-5),/*8*/-7.3693*pow(10,-3),/*9*/-3.3874*pow(10,-1),/*10*/4.7447*pow(10,0) };//MELTsからのパラメータ
			double xcpleq = xcequilibrium(T[0], P[0], Cw, Ci, cpl, tpl, AI, cl, xd[0]);
			double xcpyeq = xcequilibrium(T[0], P[0], Cw, Ci, cpy, tpy, AI, cl, xd[0]);
			double xcoeq = xcequilibrium(T[0], P[0], Cw, Ci, co, to, AI, cl, xd[0]);
			double v = velosity(P[0], Qi, T[0], xd, Ci, cl, cg, r);
			double  q = porosity(P[0], T[0], cl, xd, Ci);

			fprintf(date, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", z, P[0] / 1000000, T[0], xc[0], xcpleq, xc[1], xcpyeq, xc[2], xcoeq, q); //         データの書き込み
			double sm = entropymelt(P[0], T[0], sk[0], cv[0], c0k[0], T0[0], Y[0], cm);
			double sd = entropymelt(P[0], T[0], sk[1], cv[1], c0k[1], T0[1], Y[1], cd);
			double scpl = entropymelt(P[0], T[0], sk[2], cv[2], c0k[2], T0[2], Y[2], cpl);
			double scpy = entropymelt(P[0], T[0], sk[3], cv[3], c0k[3], T0[3], Y[3], cpy);
			double sco = entropymelt(P[0], T[0], sk[4], cv[4], c0k[4], T0[4], Y[4], co);
			double sl = entropyliquid(sm, scpl, scpy, sco, sd, xc,xd, cpl, cpy, co, cl);
			double sg = entropygas(P[0], T[0], cg);
			/*溶解度則*/Cw = solubilityw(P[0], T[0], AI);
			double/*密度方程式*/c = density(P[0], T[0], xd, Ci, cl, cg);
			double/*溶解した水(mol%)*/Cwmol = (Cw * 100 / (H * 2 + O)) / (summolfruction);
			double/*粘性の結晶依存性*/S = crystalviscosity(xd, cpl, cpy, co, cl);
			double/*粘性係数[Pas]  [mol%]*/vis = S * viscosity(SiO2mol, TiO2mol, Al2O3mol, FeOmol, MnOmol, P2O5mol, MgOmol, CaOmol, Na2Omol, Cwmol, F2Omol, K2Omol, T[0], A);
			/*速度式*/v = velosity(P[0], Qi, T[0], xd, Ci, cl, cg, r);
			double M = machnumber(P[0], Qi, T[0], xd, Ci, cl, cg, r);
			double xdeq = solubilityw(P[0], T[0], AI) / 100;
			if (solubilityw(P[0], T[0], AI) > Ci) { xdeq = Ci / 100; };

			double  φg = (1 / Gs)*(xd[0] - xdeq)*(1 - q)*cl*(1 - xc[0] - xc[1] - xc[2]);

			double xm = 1 - xd[0] - xc[0] - xc[1] - xc[2];
			double ζ[5];
			ζ[0] = (c0k[0] * C0k[0] * C0k[0] - Y[0] * P0[0]) / Y[0];//melt
			ζ[1] = (c0k[1] * C0k[1] * C0k[1] - Y[1] * P0[1]) / Y[1];//dissolved gas
			ζ[2] = (c0k[2] * C0k[2] * C0k[2] - Y[2] * P0[2]) / Y[2];//plgioclase
			ζ[3] = (c0k[3] * C0k[3] * C0k[3] - Y[3] * P0[3]) / Y[3];//pyroxene
			ζ[4] = (c0k[4] * C0k[4] * C0k[4] - Y[4] * P0[4]) / Y[4];//oliven

			dxddz[0] = -((1 - xd[0])*φg) / (cl*(1 - q)*v);
			dxcdz[0] = -(1 / (v*G))*(xc[0] - xcpleq) + (xc[0] /( (1 - q)*cl*v))*φg;
			dxcdz[1] = -(1 / (v*G))*(xc[1] - xcpyeq) + (xc[1] / ((1 - q)*cl*v))*φg;
			dxcdz[2] = -(1 / (v*G))*(xc[2] - xcoeq) + (xc[2] / ((1 - q)*cl*v))*φg;
			double Aa = xm * (cv[0] + ζ[0] / (cm*T[0])) + xd[0] * (cv[1] + ζ[1] / (cd*T[0])) + xc[0] * (cv[2] + ζ[2] / (cpl*T[0])) + xc[1] * (cv[3] + ζ[3] / (cpy*T[0])) + xc[2] * (cv[4] + ζ[4] / (co*T[0]));
			double N = (1 - q)*cl*Aa + q * cg*cv[5];

			double B = ζ[0] * xm / (cm*cm*cv[0] * (Y[0] - 1)*T[0]) + ζ[1] * xd[0] / (cd*cd*cv[1] * (Y[1] - 1)*T[0]) + ζ[2] * xc[0] / (cpl*cpl*cv[2] * (Y[2] - 1)*T[0]) + ζ[3] * xc[1] / (cpy*cpy*cv[3] * (Y[3] - 1)*T[0]) + ζ[4] * xc[2] / (co*co*cv[4] * (Y[4] - 1)*T[0]);
			double M1 = φg * (el - eg);
			double M2 = (1 - q)*cl*((-dxddz[0] - dxcdz[0] - dxcdz[1] - dxcdz[2])*em + dxddz[0] * ed + dxcdz[0] * ecpl + dxcdz[1] * ecpy + dxcdz[2] * eco);
			double M3 = ((1 - q)*cl*B + (P[0] * M*M) / (c*v))*dPdz[0];


//			printf("%.2f %.1f q=%.3f T=%.0f dtdz=%f sum=%f M1=%.0f M2=%.0f M3=%.0f xd=%f,xcpl=%f,xcpy=%f,xco=%f dPdz=%.0f\n", z, P[0] / 1000000, q,T[0],dTdz[0], (M1 - M2 + M3) / N,M1,M2,M3,xd[0],xc[0],xc[1],xc[2],dPdz[0]);
//			printf("%.2f %.1f q=%.3f T=%.0f dtdz=%f cl=%.0f xc[1]=%f xcpyeq=%f dxcdz[1]=%f dxcdz[左]=%f dxcdz[右]=%f  \n", z, P[0] / 1000000,q, T[0], dTdz[0], cl,xc[1],xcpyeq, dxcdz[0], -(xc[1] - xcpyeq) / (v*G) , (xc[1] / ((1 - q)*cl*v))*φg);
		}
	}
	printf("%f  T=%f Cw=%f 気泡が生まれました！\n",z,T[0], solubilityw(P[0], T[0], AI));


	/*****************気泡流****************/

	X = exsolvegas(P[0], T[0], Ci, xd);
	cm = densitymelt(P[0], T[0], cv[0], c0k[0], C0k[0], Y[0], P0[0]);
	cd = densitymelt(P[0], T[0], cv[1], c0k[1], C0k[1], Y[1], P0[1]);
	cpl = densitymelt(P[0], T[0], cv[2], c0k[2], C0k[2], Y[2], P0[2]);
	cpy = densitymelt(P[0], T[0], cv[3], c0k[3], C0k[3], Y[3], P0[3]);
	co = densitymelt(P[0], T[0], cv[4], c0k[4], C0k[4], Y[4], P0[4]);
	cl = densityliquid(P[0],xc,xd, cpl, cpy, co, cm, cd);
	cg = densitygas(P[0], T[0]);

	for (int i = 0; porosity(P[0], T[0], cl, xd, Ci) < φq && z <= z0 && (xc[0] + xc[1] + xc[2]) < φx; i++) {
//	for (int i = 0; z <= z0 ; i++) {
			kihoryu(z, xd, xc, T, P, dPdz);
		crystallization(z, P, xd, T, xc, dxcdz);
		solubility(z, P, T, xc, xd, dxddz);
		Temperature(z, xc, xd, dxcdz, dxddz, P, dPdz, T, dTdz);
		rk4(P, dPdz, n, z, h, Pout, xc, T, xd, xc, kihoryu);
		rk4c(xc, dxcdz, nc, n, z, h, P, xcout, T, xd, crystallization);
		rk4T(T, dTdz, nT, z, h, Tout, xc, P, dPdz, xc, xd, dxcdz, dxddz, Temperature);
		rk4S(xd, dxddz, nd, z, h, xdout, xc, T, xd, xc, P, solubility);
		P[0] = Pout[0];
		xc[0] = xcout[0];
		xc[1] = xcout[1];
		xc[2] = xcout[2];
		xd[0] = xdout[0];
		T[0] = Tout[0];
		z += h;
		double X = exsolvegas(P[0], T[0], Ci, xd);
		cm = densitymelt(P[0], T[0], cv[0], c0k[0], C0k[0], Y[0], P0[0]);
		cd = densitymelt(P[0], T[0], cv[1], c0k[1], C0k[1], Y[1], P0[1]);
		cpl = densitymelt(P[0], T[0], cv[2], c0k[2], C0k[2], Y[2], P0[2]);
		cpy = densitymelt(P[0], T[0], cv[3], c0k[3], C0k[3], Y[3], P0[3]);
		co = densitymelt(P[0], T[0], cv[4], c0k[4], C0k[4], Y[4], P0[4]);
		cl = densityliquid(P[0], xc, xd, cpl, cpy, co, cm, cd);
		cg = densitygas(P[0], T[0]);

		double Cw = solubilityw(P[0], T[0], AI);
		double tpl[10] = {/*1*/-3.0938*pow(10,-9),/*2*/-8.0678*pow(10,-6),/*3*/2.2949*pow(10,-3),/*4*/-6.9366*pow(10,-8),/*5*/-6.6318*pow(10,-4),/*6*/-1.0765*pow(10,-5),/*7*/1.1705*pow(10,-4),/*8*/1.4834*pow(10,-2),/*9*/5.3169*pow(10,-1),/*10*/-6.1889*pow(10,-0) };//MELTsからのパラメータ
		double tpy[10] = {/*1*/-3.2977*pow(10,-9),/*2*/-1.8604*pow(10,-5),/*3*/-6.0187*pow(10,-3),/*4*/4.7402*pow(10,-7),/*5*/-9.4442*pow(10,-4),/*6*/1.1645*pow(10,-5),/*7*/-4.7521*pow(10,-4),/*8*/3.9398*pow(10,-2),/*9*/9.8480*pow(10,-1),/*10*/-2.0651*pow(10,1) };//MELTsからのパラメータ
		double to[10] = {/*1*/-1.5580*pow(10,-9),/*2*/2.8401*pow(10,-6),/*3*/6.4915*pow(10,-3),/*4*/-3.9540*pow(10,-8),/*5*/2.7014*pow(10,-4),/*6*/-2.2924*pow(10,-6),/*7*/5.0093*pow(10,-5),/*8*/-7.3693*pow(10,-3),/*9*/-3.3874*pow(10,-1),/*10*/4.7447*pow(10,0) };//MELTsからのパラメータ
		double xcpleq = xcequilibrium(T[0], P[0], Cw, Ci, cpl, tpl, AI, cl, xd[0]);
		double xcpyeq = xcequilibrium(T[0], P[0], Cw, Ci, cpy, tpy, AI, cl, xd[0]);
		double xcoeq = xcequilibrium(T[0], P[0], Cw, Ci, co, to, AI, cl, xd[0]);
		double v = velosity(P[0], Qi, T[0], xd, Ci, cl, cg, r);
		double  q = porosity(P[0], T[0], cl, xd, Ci);
		fprintf(date, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", z, P[0] / 1000000, T[0], xc[0], xcpleq, xc[1], xcpyeq, xc[2], xcoeq,q); //         データの書き込み
		double sm = entropymelt(P[0], T[0], sk[0], cv[0], c0k[0], T0[0], Y[0], cm);
		double sd = entropymelt(P[0], T[0], sk[1], cv[1], c0k[1], T0[1], Y[1], cd);
		double scpl = entropymelt(P[0], T[0], sk[2], cv[2], c0k[2], T0[2], Y[2], cpl);
		double scpy = entropymelt(P[0], T[0], sk[3], cv[3], c0k[3], T0[3], Y[3], cpy);
		double sco = entropymelt(P[0], T[0], sk[4], cv[4], c0k[4], T0[4], Y[4], co);
		double sl = entropyliquid(sm, scpl, scpy, sco, sd, xc, xd, cpl, cpy, co, cl);
		double sg = entropygas(P[0], T[0], cg);
		/*溶解度則*/Cw = solubilityw(P[0], T[0], AI);
		double/*密度方程式*/c = density(P[0], T[0], xd, Ci, cl, cg);
		double/*溶解した水(mol%)*/Cwmol = (Cw * 100 / (H * 2 + O)) / (summolfruction);
		double/*粘性の結晶依存性*/S = crystalviscosity(xd, cpl, cpy, co, cl);
		double/*粘性係数[Pas]  [mol%]*/vis = S * viscosity(SiO2mol, TiO2mol, Al2O3mol, FeOmol, MnOmol, P2O5mol, MgOmol, CaOmol, Na2Omol, Cwmol, F2Omol, K2Omol, T[0], A);
		/*速度式*/v = velosity(P[0], Qi, T[0], xd, Ci, cl, cg, r);
		double M = machnumber(P[0], Qi, T[0], xd, Ci, cl, cg, r);


		double xdeq = solubilityw(P[0], T[0], AI) / 100;
		if (solubilityw(P[0], T[0], AI) > Ci) { xdeq = Ci / 100; };

		double  φg = (1 / Gs)*(xd[0] - xdeq)*(1 - q)*cl*(1 - xc[0] - xc[1] - xc[2]);

		double xm = 1 - xd[0] - xc[0] - xc[1] - xc[2];
		double ζ[5];
		ζ[0] = (c0k[0] * C0k[0] * C0k[0] - Y[0] * P0[0]) / Y[0];//melt
		ζ[1] = (c0k[1] * C0k[1] * C0k[1] - Y[1] * P0[1]) / Y[1];//dissolved gas
		ζ[2] = (c0k[2] * C0k[2] * C0k[2] - Y[2] * P0[2]) / Y[2];//plgioclase
		ζ[3] = (c0k[3] * C0k[3] * C0k[3] - Y[3] * P0[3]) / Y[3];//pyroxene
		ζ[4] = (c0k[4] * C0k[4] * C0k[4] - Y[4] * P0[4]) / Y[4];//oliven

		dxddz[0] = -((1 - xd[0])*φg) / (cl*(1 - q)*v);
		dxcdz[0] = -(1 / (v*G))*(xc[0] - xcpleq) + (xc[0] / ((1 - q)*cl*v))*φg;
		dxcdz[1] = -(1 / (v*G))*(xc[1] - xcpyeq) + (xc[1] / ((1 - q)*cl*v))*φg;
		dxcdz[2] = -(1 / (v*G))*(xc[2] - xcoeq) + (xc[2] / ((1 - q)*cl*v))*φg;
		double Aa = xm * (cv[0] + ζ[0] / (cm*T[0])) + xd[0] * (cv[1] + ζ[1] / (cd*T[0])) + xc[0] * (cv[2] + ζ[2] / (cpl*T[0])) + xc[1] * (cv[3] + ζ[3] / (cpy*T[0])) + xc[2] * (cv[4] + ζ[4] / (co*T[0]));
		double N = (1 - q)*cl*Aa + q * cg*cv[5];

		double B = ζ[0] * xm / (cm*cm*cv[0] * (Y[0] - 1)*T[0]) + ζ[1] * xd[0] / (cd*cd*cv[1] * (Y[1] - 1)*T[0]) + ζ[2] * xc[0] / (cpl*cpl*cv[2] * (Y[2] - 1)*T[0]) + ζ[3] * xc[1] / (cpy*cpy*cv[3] * (Y[3] - 1)*T[0]) + ζ[4] * xc[2] / (co*co*cv[4] * (Y[4] - 1)*T[0]);
		double M1 = φg * (el - eg);
		double M2 = (1 - q)*cl*((-dxddz[0] - dxcdz[0] - dxcdz[1] - dxcdz[2])*em + dxddz[0] * ed + dxcdz[0] * ecpl + dxcdz[1] * ecpy + dxcdz[2] * eco);
		double M3 = ((1 - q)*cl*B + (P[0] * M*M) / (c*v))*dPdz[0];

//		printf("%.2f %.1f q=%.3f T=%.0f dtdz=%f sum=%f M1=%.0f M2=%.0f M3=%.0f xd=%f,xcpl=%f,xcpy=%f,xco=%f dPdz=%.0f\n", z, P[0] / 1000000, q, T[0], dTdz[0], (M1 - M2 + M3) / N,M1, M2, M3, xd[0], xc[0], xc[1], xc[2], dPdz[0]);
//		printf("%.2f %.1f q=%.3f T=%.0f dtdz=%f cl=%.0f xc[1]=%f xcpyeq=%f dxcdz[1]=%f dxcdz[左]=%f dxcdz[右]=%f  \n", z, P[0] / 1000000,q, T[0], dTdz[0], cl, xc[1], xcpyeq, dxcdz[0], -(xc[1] - xcpyeq) / (v*G), (xc[1] / ((1 - q)*cl*v))*φg);

	}

	printf("%.1f %.2f T=%f q=%f破砕破砕破砕破砕破砕破砕破砕破砕破砕\n\n\n\n",z,P[0] / 1000000,T[0],porosity(P[0], T[0], cl, xd, Ci));
	v = velosity(P[0], Qi, T[0], xd, Ci, cl, cg, r);

	
	/***************噴霧流******************/


	for (int i = 0; soundspeed(P[0], cl, T[0], Ci, xd, cg)-velosity(P[0], Qi, T[0], xd, Ci, cl, cg, r) > 0 && z <= z0; i++) {
		hunmryu(z, xd, xc, T, P, dPdz);
		crystallization(z, P, xd, T, xc, dxcdz);
		Temperature(z, xc, xd, dxcdz, dxddz, P, dPdz, T, dTdz);
		solubility(z, P, T, xc, xd, dxddz);
		rk4(P, dPdz, n, z, h, Pout, xc, T, xd, xc, hunmryu);
		rk4c(xc, dxcdz, nc, n, z, h, P, xcout, T, xd, crystallization);
		rk4S(xd, dxddz, nd, z, h, xdout, xc, T, xd, xc, P, solubility);
		rk4T(T, dTdz, nT, z, h, Tout, xc, P, dPdz, xc, xd, dxcdz, dxddz, Temperature);

		P[0] = Pout[0];
		xc[0] = xcout[0];
		xc[1] = xcout[1];
		xc[2] = xcout[2];
		xd[0] = xdout[0];
		T[0] = Tout[0];
		z += h;
		double X = exsolvegas(P[0], T[0], Ci, xd);
		cm = densitymelt(P[0], T[0], cv[0], c0k[0], C0k[0], Y[0], P0[0]);
		cd = densitymelt(P[0], T[0], cv[1], c0k[1], C0k[1], Y[1], P0[1]);
		cpl = densitymelt(P[0], T[0], cv[2], c0k[2], C0k[2], Y[2], P0[2]);
		cpy = densitymelt(P[0], T[0], cv[3], c0k[3], C0k[3], Y[3], P0[3]);
		co = densitymelt(P[0], T[0], cv[4], c0k[4], C0k[4], Y[4], P0[4]);
		cl = densityliquid(P[0], xc, xd, cpl, cpy, co, cm, cd);
		cg = densitygas(P[0], T[0]);

		double Cw = solubilityw(P[0], T[0], AI);
		double tpl[10] = {/*1*/-3.0938*pow(10,-9),/*2*/-8.0678*pow(10,-6),/*3*/2.2949*pow(10,-3),/*4*/-6.9366*pow(10,-8),/*5*/-6.6318*pow(10,-4),/*6*/-1.0765*pow(10,-5),/*7*/1.1705*pow(10,-4),/*8*/1.4834*pow(10,-2),/*9*/5.3169*pow(10,-1),/*10*/-6.1889*pow(10,-0) };//MELTsからのパラメータ
		double tpy[10] = {/*1*/-3.2977*pow(10,-9),/*2*/-1.8604*pow(10,-5),/*3*/-6.0187*pow(10,-3),/*4*/4.7402*pow(10,-7),/*5*/-9.4442*pow(10,-4),/*6*/1.1645*pow(10,-5),/*7*/-4.7521*pow(10,-4),/*8*/3.9398*pow(10,-2),/*9*/9.8480*pow(10,-1),/*10*/-2.0651*pow(10,1) };//MELTsからのパラメータ
		double to[10] = {/*1*/-1.5580*pow(10,-9),/*2*/2.8401*pow(10,-6),/*3*/6.4915*pow(10,-3),/*4*/-3.9540*pow(10,-8),/*5*/2.7014*pow(10,-4),/*6*/-2.2924*pow(10,-6),/*7*/5.0093*pow(10,-5),/*8*/-7.3693*pow(10,-3),/*9*/-3.3874*pow(10,-1),/*10*/4.7447*pow(10,0) };//MELTsからのパラメータ
		double xcpleq = xcequilibrium(T[0], P[0], Cw, Ci, cpl, tpl, AI, cl, xd[0]);
		double xcpyeq = xcequilibrium(T[0], P[0], Cw, Ci, cpy, tpy, AI, cl, xd[0]);
		double xcoeq = xcequilibrium(T[0], P[0], Cw, Ci, co, to, AI, cl, xd[0]);
		double v = velosity(P[0], Qi, T[0], xd, Ci, cl, cg, r);
		double  q = porosity(P[0], T[0], cl, xd, Ci);
		fprintf(date, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", z, P[0] / 1000000, T[0], xc[0],xcpleq, xc[1],xcpyeq, xc[2],xcoeq, q); //         データの書き込み
		double sm = entropymelt(P[0], T[0], sk[0], cv[0], c0k[0], T0[0], Y[0], cm);
		double sd = entropymelt(P[0], T[0], sk[1], cv[1], c0k[1], T0[1], Y[1], cd);
		double scpl = entropymelt(P[0], T[0], sk[2], cv[2], c0k[2], T0[2], Y[2], cpl);
		double scpy = entropymelt(P[0], T[0], sk[3], cv[3], c0k[3], T0[3], Y[3], cpy);
		double sco = entropymelt(P[0], T[0], sk[4], cv[4], c0k[4], T0[4], Y[4], co);
		double sl = entropyliquid(sm, scpl, scpy, sco, sd, xc, xd, cpl, cpy, co, cl);
		double sg = entropygas(P[0], T[0], cg);
		/*溶解度則*/Cw = solubilityw(P[0], T[0], AI);
		double/*密度方程式*/c = density(P[0], T[0], xd, Ci, cl, cg);
		double/*溶解した水(mol%)*/Cwmol = (Cw * 100 / (H * 2 + O)) / (summolfruction);
		double/*粘性の結晶依存性*/S = crystalviscosity(xd, cpl, cpy, co, cl);
		double/*粘性係数[Pas]  [mol%]*/vis = S * viscosity(SiO2mol, TiO2mol, Al2O3mol, FeOmol, MnOmol, P2O5mol, MgOmol, CaOmol, Na2Omol, Cwmol, F2Omol, K2Omol, T[0], A);
		/*速度式*/v = velosity(P[0], Qi, T[0], xd, Ci, cl, cg, r);
		double M = machnumber(P[0], Qi, T[0], xd, Ci, cl, cg, r);
		double xdeq = solubilityw(P[0], T[0], AI) / 100;
		if (solubilityw(P[0], T[0], AI) > Ci) { xdeq = Ci / 100; };
		double  φg = (1 / Gs)*(xd[0] - xdeq)*(1 - q)*cl*(1 - xc[0] - xc[1] - xc[2]);
		double xm = 1 - xd[0] - xc[0] - xc[1] - xc[2];
		double ζ[5];


		ζ[0] = (c0k[0] * C0k[0] * C0k[0] - Y[0] * P0[0]) / Y[0];//melt
		ζ[1] = (c0k[1] * C0k[1] * C0k[1] - Y[1] * P0[1]) / Y[1];//dissolved gas
		ζ[2] = (c0k[2] * C0k[2] * C0k[2] - Y[2] * P0[2]) / Y[2];//plgioclase
		ζ[3] = (c0k[3] * C0k[3] * C0k[3] - Y[3] * P0[3]) / Y[3];//pyroxene
		ζ[4] = (c0k[4] * C0k[4] * C0k[4] - Y[4] * P0[4]) / Y[4];//oliven

		dxddz[0] = -((1 - xd[0])*φg) / (cl*(1 - q)*v);
		dxcdz[0] = -(1 / (v*G))*(xc[0] - xcpleq) + (xc[0] / ((1 - q)*cl*v))*φg;
		dxcdz[1] = -(1 / (v*G))*(xc[1] - xcpyeq) + (xc[1] / ((1 - q)*cl*v))*φg;
		dxcdz[2] = -(1 / (v*G))*(xc[2] - xcoeq) + (xc[2] / ((1 - q)*cl*v))*φg;
		double Aa = xm * (cv[0] + ζ[0] / (cm*T[0])) + xd[0] * (cv[1] + ζ[1] / (cd*T[0])) + xc[0] * (cv[2] + ζ[2] / (cpl*T[0])) + xc[1] * (cv[3] + ζ[3] / (cpy*T[0])) + xc[2] * (cv[4] + ζ[4] / (co*T[0]));
		double N = (1 - q)*cl*v*Aa + q * v*cg*cv[5];

		double B = ζ[0] * xm / (cm*cm*cv[0] * (Y[0] - 1)*T[0]) + ζ[1] * xd[0] / (cd*cd*cv[1] * (Y[1] - 1)*T[0]) + ζ[2] * xc[0] / (cpl*cpl*cv[2] * (Y[2] - 1)*T[0]) + ζ[3] * xc[1] / (cpy*cpy*cv[3] * (Y[3] - 1)*T[0]) + ζ[4] * xc[2] / (co*co*cv[4] * (Y[4] - 1)*T[0]);
		double M1 = φg * (el - eg);
		double M2 = (1 - q)*cl*v*((-dxddz[0] - dxcdz[0] - dxcdz[1] - dxcdz[2])*em + dxddz[0] * ed + dxcdz[0] * ecpl + dxcdz[1] * ecpy + dxcdz[2] * eco);
		double M3 = ((1 - q)*cl*B*v + (P[0] * M*M) / (c*v))*dPdz[0];


//		printf("%.2f %.1f T=%.0f dtdz=%f cl=%.0f N=%.0f M1=%.0f M2=%.0f M3=%.0f dxddz=%f dxcdz[0]=%f dxcdz[1]=%f dxcdz[2]=%f  \n", z, P[0] / 1000000, T[0],dTdz[0],cl,N,M1,M2,M3, dxddz [0],dxcdz[0] ,dxcdz[1], dxcdz[2]);
//		printf("%.2f %.1f q=%.3f T=%.0f dtdz=%f M1=%.0f M2=%.0f M3=%.0f xd=%f,xcpl=%f,xcpy=%f,xco=%f dPdz=%.0f\n", z, P[0] / 1000000, q, T[0], dTdz[0], M1, M2, M3, xd[0], xc[0], xc[1], xc[2], dPdz[0]);
	}

	printf("%.1f P=%fチョークチョークチョークチョークチョークチョークチョークチョークチョークチョークチョークチョークチョークチョークチョーク\n\n\n\n",z,P[0] / 1000000);
	fclose(date);

	/*-------------------バッチファイルの作成--------------------*/
	bat_file = "command.gp";
	command = fopen(bat_file, "w");

	fprintf(command, "set yrange [-6.5:10000]\n"); // 範囲の指定(省略可)
	fprintf(command, "set xrange [1280:1350]\n");
//	fprintf(command, "set xrange [0:1]\n");
	fprintf(command, "plot '%s' using 3:1 with lines\n", date_file);
//	fprintf(command, "plot '%s' using 3:1 with lines,'' using 4:1 with lines,'' using 5:1 with lines,'' using 6:1 with lines,'' using 7:1 with lines,'' using 8:1 with lines,'' using 9:1 with lines,'' using 10:1 with lines\n", date_file);
	fprintf(command, "pause -1 'Hit any key to close plot window.'");
	fclose(command);
/**/
	/*---------------------gnuplotの実行---------,'' using 5:1 with lines,'' using 6:1 with lines-------------*/
	sprintf(run_command, "C:/PROGRA~1/gnuplot/bin/wgnuplot.exe \"%s\"", bat_file);
	system(run_command);

	return 0;
}
