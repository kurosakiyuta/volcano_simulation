#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "solubilityH2O.h"
#include "É¿eq.h"

double xcequilibrium(double T, double P,double Cw, double Ci,double cj, double t[],double AI,double cl,double xd) {

	double Pbar = P / 100000;
	double Td = T - 273.65;
	double d = xd * 100;

	double xceq = (t[0]*Pbar*Pbar + t[1]*Td*Td + t[2]*d*d + t[3]*Pbar*Td + t[4]*Td*d + t[5]*d*Pbar + t[6]*Pbar + t[7]*Td + t[8]*d + t[9]);
	if (xceq < 0) { xceq = 0; };
	return xceq;
}