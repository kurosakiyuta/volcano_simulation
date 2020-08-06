#include <stdio.h>
#include <math.h>
#include "density(liquid).h"

double densityliquidwt(double xc[] ,double c[]) {

	double xm = 1 - xc[0] - xc[1] - xc[2] - xc[3] - xc[4] - xc[5] - xc[6];
	double clre = xm /c[5] + xc[0] / c[0] + xc[1] / c[1] + xc[2] / c[2] + xc[3] / c[3] + xc[4] / c[4] + xc[5] / c[6] + xc[6] / c[7];
	double cl = 1 / clre;
	return cl;
}

double densityliquidvol(double b[], double c[]) {

	double am = 1 - b[0] - b[1] - b[2] - b[3] - b[4] - b[5] - b[6];
	double cl = am * c[5] + b[0] * c[0] + b[1] * c[1] + b[2] * c[2] + b[3] * c[3] + b[4] * c[4] + b[5] * c[5] + b[6] * c[6];
	return cl;
}


