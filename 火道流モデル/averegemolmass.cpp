#include <stdio.h>
#include <math.h>
#include "averagemolmass.h"


#define H 0.001
#define C 0.012
#define O 0.016
#define F 0.019          /*mol%確認が必要*/
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



/*水分子以外の平均モル質量*/

double averagemolmass(double  SiO2, double TiO2, double Al2O3, double FeO, double MnO, double P2O5, double MgO, double CaO, double Na2O, double F2O, double K2O) {

	double aveM = (SiO2 * (Si + 2 * O) + TiO2 * (Ti + 2 * O) + Al2O3 * (2 * Al + 3 * O) + FeO * (Fe + O) +
		MnO * (Mn + O) + P2O5 * (2 * P + 5 * O) + MgO * (Mg + O) + CaO * (Ca + O) + Na2O * (Na + 2 * O) +
		F2O * (2 * F + O) + K2O * (2 * K + O))/100;

	return aveM;
}