#include <stdio.h>
#include <math.h>
#include "solubilityH2O.h"
#include "viscosity.h"
#include "/Users/DELL/source/C++/数値計算C++/数値計算C++/nrutil.h"
#include "/Users/DELL/source/C++/数値計算C++/数値計算C++/rk4.h"
#include "density.h"
#include "velocity.h"
#include "porosity.h"
#include "dcdp.h"
#include"Exsolvedgas.h"
#include "averagemolmass.h"
#include "Liquid Water Mol.h"

double watermol(double P, double T, double Ci, double AI, double  SiO2, double TiO2, double Al2O3,
	double FeO, double MnO, double P2O5, double MgO, double CaO, double Na2O, double F2O, double K2O) {


	double X, Ce, Mave, Mw;

	Mw = 0.018;
	X = exsolvegas(P, T, Ci, AI);
	Ce = Ci - X;
	Mave = averagemolmass(SiO2, TiO2, Al2O3, FeO, MnO, P2O5, MgO, CaO, Na2O, F2O, K2O);
	double Wmol = 100*Ce / ((Mw*(100 - X - Ce) / Mave) + Ce + X);
	return Wmol;
}