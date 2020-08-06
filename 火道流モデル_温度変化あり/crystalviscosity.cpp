#include <stdio.h>
#include <math.h>
#include "crystalviscosity.h"

double crystalviscosity(double xc[],double cpl,double cpy, double co,double cl) {//costa et al 2009, vona et al (2011)ƒpƒ‰ƒƒ^Œˆ’è
	double a, B, c, d, e, g, F,x[3];//g=ƒÓ
	B = 2.8;    //B
	c = 0.39;   //ƒÓ*
	d = 0.03;   //ƒÄ
	e = 0.84;   //ƒÁ
	a = 2 - e;

	g = (xc[0]+xc[1]+xc[2]) / c;
	F = (1 - d)*erf(sqrt(3.1416)*g*(1 + pow(g, a)) / (2 * (1 - d)));

	double S = (1 + pow(g, a)) / pow(1 - F, B*c);
	return S;
}