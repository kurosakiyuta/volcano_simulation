#include <stdio.h>
#include <math.h>
#include "solubilityH2O.h"
#include "porosity.h"
#include"Exsolvedgas.h"

double porosity(double P, double T, double cl, double xd[], double Ci) {
	double R, Mw,j,X;
	R = 461.88888;
	Mw= 0.01802;
	X = exsolvegas(P, T, Ci, xd);
	double n = X / 100;
	double A = P / (n*R*T);
	double B = cl / (1 - n);
	j = B / (A + B);

	return j;
}