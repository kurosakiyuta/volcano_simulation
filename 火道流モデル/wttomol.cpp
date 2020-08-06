#include <stdio.h>
#include <math.h>
#include "solubilityH2O.h"
#include "Exsolvedgas.h"


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
void wttomol(double comp_melt[],double comp_melt_mol[],double P, double tw[],double xc[], double watersuminitial) {
	/*Cwは、溶存水のwt％・・H2Oは揮発性成分の水ｗｔ％*/
	double xm = 1 - xc[0] - xc[1] - xc[2] - xc[3] - xc[4] - xc[5] - xc[6];
	double n = exsolvegas(P, watersuminitial, tw, xc);
	double xw = solubilityw(P, tw) / 100;
	if (xw*xm*(1 - n) + n > watersuminitial) { xw = (watersuminitial - n) / (xm*(1 - n)); };//watersuminitialは、全岩中であり、溶解量はメルト中。メルト割合＊液体割合＝全体の割合＊

	double summolfruction =
		comp_melt[SiO2]/ (Si + 2 * O)
		+ comp_melt[TiO2] / (Ti + 2 * O) 
		+ comp_melt[Al2O3] / (2 * Al + 3 * O) 
		+ comp_melt[Fe2O3] / (2 * Fe + 3 * O) 
		+ comp_melt[FeO] / (Fe + O) 
		+ comp_melt[MnO] / (Mn + O) 
		+ comp_melt[P2O5] / (2 * P + 5 * O) 
		+ comp_melt[MgO] / (Mg + O) 
		+ comp_melt[CaO] / (Ca + O) 
		+ comp_melt[Na2O] / (Na + 2 * O) 
		+ F2O / (2 * F + O) 
		+ comp_melt[K2O] / (2 * K + O) 
		+ xw / (H * 2 + O);



	comp_melt_mol[SiO2] = 100*(comp_melt[SiO2] / (Si + 2 * O)) / summolfruction;
	comp_melt_mol[TiO2] = 100 * (comp_melt[TiO2] / (Ti + 2 * O)) / summolfruction;
	comp_melt_mol[Al2O3]= 100 * (comp_melt[Al2O3] / (2 * Al + 3 * O)) / summolfruction;
	comp_melt_mol[FeO] = 100 * ((comp_melt[FeO] / (Fe + O)) / summolfruction);
	comp_melt_mol[Fe2O3]= 100 * ((comp_melt[Fe2O3] / (2 * Fe + 3 * O)) / summolfruction);
	comp_melt_mol[MnO] = 100 * (comp_melt[MnO] / (Mn + O)) / summolfruction;
	comp_melt_mol[MgO] = 100 * (comp_melt[MgO] / (Mg + O)) / summolfruction;
	comp_melt_mol[CaO] = 100 * (comp_melt[CaO] / (Ca + O)) / summolfruction;
	comp_melt_mol[Na2O] = 100 * (comp_melt[Na2O] / (Na + 2 * O)) / summolfruction;
	comp_melt_mol[K2O] = 100 * (comp_melt[K2O] / (2 * K + O)) / summolfruction;
	comp_melt_mol[P2O5] = 100 * (comp_melt[P2O5] / (2 * P + 5 * O)) / summolfruction;
	comp_melt_mol[H2O] = 100 * (xw / (H * 2 + O)) / summolfruction;

//	for (int i = 0; i < 12; i++) {
//		printf("comp_melt_mo[%d]=%f\n", i, comp_melt_mol[i]);
//	}

}
