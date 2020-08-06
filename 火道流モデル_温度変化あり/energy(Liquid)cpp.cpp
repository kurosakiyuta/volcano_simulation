#include <stdio.h>
#include <math.h>
#include "energy(Liquid).h"

double energyliquid(double xd[], double xc[],double em,double ed, double ecpl, double ecpy, double eco) {
	double xm = 1 - xd[0] - xc[0] - xc[1] - xc[2];
	double el = (xm)*em + xd[0] * ed + xc[0] * ecpl + xc[1] * ecpy + xc[2] * eco;
	return el;
}