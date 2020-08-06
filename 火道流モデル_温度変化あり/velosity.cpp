#include <stdio.h>
#include <math.h>
#include "density.h"
#include "velocity.h"

double velosity(double P, double Qi, double T, double xd[], double Ci, double cl,double cg,double r) {

	double c, v;
	c = density(P, T, xd, Ci,cl,cg);
	double A = 3.14159*r*r;
	v = Qi / (c*A);
	return v;
}