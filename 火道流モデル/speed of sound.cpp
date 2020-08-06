#include <stdio.h>
#include <math.h>
#include "solubilityH2O.h"
#include "density.h"
#include "speed of sound.h"
#include "Exsolvedgas.h"
#include "porosity.h"
#include "density(liquid).h"

double soundspeed(double P, double T, double watersuminitial, double tw[], double c[],double b[],double xc[]) {
	double n = exsolvegas(P, watersuminitial,tw, xc);
	double R = 461.88888;
	double Mw = 0.01802;
	double c_bulk = density(P, T,tw, watersuminitial, c, b);
	double cl = densityliquidvol(b, c);
	double q = porosity(P, T, c, tw, watersuminitial,xc);

	double cp = P/(c_bulk*q) ;

	double a = sqrt(cp);
	return a ;

}