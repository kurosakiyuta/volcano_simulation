#include <stdio.h>
#include <math.h>
#include "solubilityH2O.h"
#include "density.h"
#include "Exsolvedgas.h"
#include "porosity.h"
#include "density(liquid).h"

double density(double P, double T, double tw[] ,  double watersuminitial, double c[],double b[]) {
	int i;
	double c_bulk, X,Mw,R;
	Mw = 0.01802;
	R = 461.88888;
	double cl = densityliquidvol(b, c);
	double q = porosity(P, T, c,tw, watersuminitial,b);
	
	c_bulk = (q*P) / (R*T) + (1 - q)*cl;

	return c_bulk;
}
