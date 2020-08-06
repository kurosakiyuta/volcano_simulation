#include <stdio.h>
#include <math.h>
#include "density(melt etc).h"

double densitymelt(double P, double T, double cv, double c0, double C0, double Y, double P0) {

	double A = (c0*C0*C0 - Y * P0) / (Y);
	double c = (P + A) / (cv*(Y - 1)*T);

	return c;
}