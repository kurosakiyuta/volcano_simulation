#include <stdio.h>
#include <math.h>
#include "solubilityH2O.h"
#include "porosity.h"
#include "Exsolvedgas.h"
#include "density(liquid).h"

double porosity(double P, double T, double c[], double tw[], double watersuminitial,double xc[]) {
	double R, Mw,j,n;
	R = 461.88888;
	Mw= 0.01802;
	double cl = densityliquidwt(xc, c);
	n = exsolvegas(P, watersuminitial, tw,xc );
	double A = P / (n*R*T);
	double B = cl / (1 - n);
	j = B / (A + B);



	return j;
}