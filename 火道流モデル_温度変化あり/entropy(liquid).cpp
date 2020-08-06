#include <stdio.h>
#include <math.h>
#include "entropy(liquid).h"

double entropyliquid(double sm, double scpl, double scpy, double sco,double sd,double xc[],double xd[], double cpl, double cpy, double co, double cl) {
	double xm = 1 - xd[0] - xc[0] - xc[1] - xc[2];
	double sl = (xm)*sm + xc[0] * scpl + xc[1] * scpy + xc[2] * sco + xd[0] * sd;
	return sl;
}