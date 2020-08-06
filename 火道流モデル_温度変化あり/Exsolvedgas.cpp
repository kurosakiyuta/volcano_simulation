#include <stdio.h>
#include <math.h>
#include "solubilityH2O.h"
#include"Exsolvedgas.h"

double exsolvegas(double P, double T, double Ci, double xd[]) {//wt%
	double X = 0;
	X = Ci - xd[0]*100;
	

	return X;
}