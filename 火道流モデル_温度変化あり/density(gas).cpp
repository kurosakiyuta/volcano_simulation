#include <stdio.h>
#include <math.h>
#include "density(gas).h"

double densitygas(double P, double T) {

	double R = 461.88888;
	double cg = P / (R*T);
	return cg;
}