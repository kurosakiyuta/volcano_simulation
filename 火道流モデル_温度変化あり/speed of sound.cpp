#include <stdio.h>
#include <math.h>
#include "solubilityH2O.h"
#include "density.h"
#include "speed of sound.h"
#include "Exsolvedgas.h"
#include "porosity.h"
double soundspeed(double P, double cl, double T, double Ci, double xd[], double cg) {
	double X = exsolvegas(P, T, Ci, xd);
	double R = 461.88888;
	double n = X / 100;
	double c = density(P, T, xd, Ci, cl,cg);
	double q = porosity(P, T, cl, xd, Ci);

	double cp = P/(c*q) ;

	double a = sqrt(cp);
	return a ;
}