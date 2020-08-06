#include <stdio.h>
#include <math.h>
#include "relativeviscosity.h"
#include "MPF(bimodal).h"
double crystalviscosity_shear(double b[],double bi[],double rp,double nc,double v,double r) {//costa et al 2009, vona et al (2011)パラメタ決定
	double bsum=0;
	double φ = MPF_bimo(rp, b, bi);
	double shear = v / r;

	for (int i=0; i < nc; i++) {
		bsum += b[i];
	}

	double ff = (bsum / φ);
	if (ff > 0.8) { ff = 0.8; };
	if (rp > 10) { rp = 10; };
	double n = 1 - 0.2*(rp)*pow(ff, 4);
//	printf("n=%f\n", n);
	double R = pow(shear, n - 1);
	double S = pow(1 - (bsum / φ),-2);

	return S * R;
}

double crystalviscosity(double b[], double φ, double nc) {//costa et al 2009, vona et al (2011)パラメタ決定
	double bsum = 0;

	for (int i = 0; i < nc; i++) {
		bsum += b[i];
	}

	double S = pow(1 - (bsum / φ), -2);
	return S;
}