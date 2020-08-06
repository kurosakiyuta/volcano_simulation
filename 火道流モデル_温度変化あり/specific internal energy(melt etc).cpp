#include <stdio.h>
#include <math.h>
#include "specific internal energy(melt etc).h"

double energymelt(double P, double T, double ek, double cv, double c0, double C0, double Y, double P0, double ck) {
	//ek=e_, co= Éœ0, Y=É¡,

	double e = ek + cv * T + (c0*C0*C0 - Y * P0) / (Y*ck);
	return e;
}