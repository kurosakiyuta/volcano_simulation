#include <stdio.h>
#include <math.h>
#include "solubilityH2O.h"
#include "Exsolvedgas.h"

double exsolvegas(double P, double watersuminitial, double tw[], double xc[]) {
	double x, xm;

	xm = (1 - xc[0] - xc[1] - xc[2] - xc[3] - xc[4] - xc[5] - xc[6]);
	double Cw = solubilityw(P, tw);
	double xw = Cw / 100;
	x = 0;

	if ((watersuminitial - xm * xw) >= 0) { x = (watersuminitial - xm * xw) / (1 - xm * xw); };

	return x;
}

/*	double s = 4 * pow(10, -6);
	double m = 0.5;
	double n0 = Ci;

	if (n0 > s*pow(P, m)) {
		X = (100 * (n0 - s * pow(P, m))) / (1 - s * pow(P, m));
	}*/
