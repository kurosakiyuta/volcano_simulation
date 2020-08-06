#include <stdio.h>
#include <math.h>
#include "solubilityH2O.h"

/**********************************************/
/*深さと温度の変化に伴うマントル中の水の溶解度*/
/********つまり、マントル中の水のwt%***********/
/**********************************************/

double solubilityw(double P, double tw[]) {
	double Cw, lnCw;

	lnCw = tw[0]+ tw[1] * (log(P));
	Cw = exp(lnCw);

	return Cw;//[wt%]
}



/*	double SP;
	double MP = P * pow(10,-6);

	SP = sqrt(MP);

	Cw = (-0.231 + 651.1 / T)*SP + (0.03424 - 32.57 / T + 0.2447*AI)*MP;//MPa

/*	double sigma = 5.0025*pow(10, -7);
	double e = 0.6;
	Cw =100* sigma * pow(P, e);*/
