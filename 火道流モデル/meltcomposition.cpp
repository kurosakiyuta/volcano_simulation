#include <stdio.h>
#include <math.h>
#include "meltcomposition.h"
#include "solubilityH2O.h"
#include "Exsolvedgas.h"

void composiotion_melt(double comp_melt_ini[], double comp_melt[], double comp_pl[], double comp_py[], double comp_o[], double comp_q[], double comp_s[],double comp_cl[], double comp_or[], double xc[],double xci[],double P, double tw[],double watersuminitial) {
// 0=SiO2  1=TiO2	2=Al2O3	3=Fe2O3 4=FeO   5=MnO   6=MgO	7=CaO  	8=Na2O  9=K2O   10= P2O5
	double xm = 1 - xc[0] - xc[1] - xc[2] - xc[3] - xc[4] - xc[5] - xc[6];
	double xmi = 1 - xci[0] - xci[1] - xci[2] - xci[3] - xci[4] - xci[5] - xci[6];

//初期状態における組成iの全体の中の割合
	double comp_ini[11];
	for (int i = 0; i < 11; i++) {
		comp_ini[i] = comp_melt_ini[i] * xmi
			+ comp_pl[i] * xci[0]
			+ comp_py[i] * xci[1]
			+ comp_o[i] * xci[2]
			+ comp_q[i] * xci[3]
			+ comp_s[i] * xci[4]
			+ comp_cl[i] * xci[5]
			+ comp_or[i] * xci[6];
		//	printf("comp_ini[%d]=%f\n",i, comp_ini[i]);
	}

//成分iの結晶中の成分量
	double comp_crystal[11];
	for (int i = 0; i < 11; i++) {
		comp_crystal[i] = xc[0] * comp_pl[i] + xc[1] * comp_py[i] + xc[2] * comp_o[i] + xc[3] * comp_q[i] + xc[4] * comp_s[i] + xc[5] * comp_cl[i] + xc[6] * comp_or[i];
	}
//成分iのメルト中の成分量
	for (int i = 0; i < 11; i++) {
		comp_melt[i] = (comp_ini[i] - comp_crystal[i]) / xm;
//		printf("comp_melt[%d]=%f\n", i, comp_melt[i]);

	}


//水を除く成分量の正規化
	double comp_melt_sum = 0;
	for (int i = 0; i < 11; i++) {
		comp_melt_sum += comp_melt[i];
	}
	for (int i = 0; i < 11; i++) {
		comp_melt[i] = comp_melt[i] / comp_melt_sum;

	}

//溶解量は固定であるから、それを考えてメルト全体で再度正規化する。
//メルトには水が入っているので、それを入れる。そして、全体の正規化した値を得る

	double n = exsolvegas(P, watersuminitial, tw, xc);
	double xw = solubilityw(P, tw) / 100;
	if (xw*xm*(1-n)+n > watersuminitial) { xw = (watersuminitial-n)/(xm*(1-n)); };//watersuminitialは、全岩中であり、溶解量はメルト中。メルト割合＊液体割合＝全体の割合＊
	for (int i = 0; i < 11; i++) {
		comp_melt[i] = comp_melt[i] / (1 + xw);
//		printf("comp_meltfinal[%d]=%f\n", i, comp_melt[i]);

	}

//正規化されたメルトの組成が破壊的に求められた。水の溶解量を変化させないようにしたので、メルトの溶解量はそのまま。

}

