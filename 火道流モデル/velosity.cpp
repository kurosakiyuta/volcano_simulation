#include <stdio.h>
#include <math.h>
#include "solubilityH2O.h"
#include "density.h"
#include "velocity.h"

double velosity(double P, double Qi, double T, double tw[], double watersuminitial, double c[], double b[],double r) {

	double c_bulk, v,A;
	A = 3.14159*r*r;

	c_bulk = density(P, T, tw, watersuminitial, c, b);
	v = Qi / (c_bulk*A);

	return v;
}