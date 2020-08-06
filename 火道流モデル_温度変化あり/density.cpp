#include <stdio.h>
#include <math.h>
#include "solubilityH2O.h"
#include "density.h"
#include "Exsolvedgas.h"
#include "porosity.h"


double density(double P, double T, double xd[], double Ci, double cl,double cg) {


	double q = porosity(P, T, cl, xd, Ci);
	
	double c = cl * (1 - q) + q * cg;
	return c;
}
