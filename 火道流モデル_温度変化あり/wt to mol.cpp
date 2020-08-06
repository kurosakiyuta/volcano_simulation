#include <stdio.h>/*各成分のwt％をmol%に変換*/
#include <math.h>
#include "wt to mol.h"

#define H 0.001
#define C 0.012
#define O 0.016
#define F 0.019        /*原子の分子量kg/mol*/
#define Na 0.023
#define Mg 0.0243
#define Al 0.027
#define Si 0.028
#define P 0.031
#define K 0.039
#define Ca 0.04
#define Ti 0.047867
#define Mn 0.055
#define Fe 0.055845
void wttomol(double SiO2,double TiO2,double Al2O3,double FeO,double MnO,double P2O5,double MgO,double CaO,double Na2O,double F2O,double K2O,double Cw,double Ci) {

	/*Cwは、溶存水のwt％・・H2Oは揮発性成分の水ｗｔ％*/
	double SiO2mol, TiO2mol, Al2O3mol, FeOmol, MnOmol, P2O5mol, MgOmol, CaOmol, Na2Omol, F2Omol, K2Omol, Cwmol, H2Omol,H2O;
	H2O = Ci - Cw;

	double summolfruction = SiO2 / (Si + 2 * O) + TiO2 / (Ti + 2 * O) + Al2O3 / (2 * Al + 3 * O) + FeO / (Fe + O) +
		MnO / (Mn + O) + P2O5 / (2 * P + 5 * O) + MgO / (Mg + O) + CaO / (Ca + O) + Na2O / (Na + 2 * O) +
		F2O / (2 * F + O) + K2O / (2 * K + O) + Cw / (H * 2 + O) + H2O / (H * 2 + O);

	SiO2mol = (SiO2 / (Si + 2 * O)) / summolfruction;
	TiO2mol = (TiO2 / (Ti + 2 * O)) / summolfruction;
	Al2O3mol = (Al2O3 / (2 * Al + 3 * O)) / summolfruction;
	FeOmol = (FeO / (Fe + O)) / summolfruction;
	MnOmol = (MnO / (Mn + O)) / summolfruction;
	P2O5mol = (P2O5 / (2 * P + 5 * O)) / summolfruction;
	MgOmol = (MgO / (Mg + O)) / summolfruction;
	CaOmol = (CaO / (Ca + O)) / summolfruction;
	Na2Omol = (Na2O / (Na + 2 * O)) / summolfruction;
	F2Omol = (F2O / (2 * F + O)) / summolfruction;
	K2Omol = (K2O / (2 * K + O)) / summolfruction;
	Cwmol = (Cw / (H * 2 + O)) / summolfruction;
	H2Omol = (H2O / (H * 2 + O)) / summolfruction;
}