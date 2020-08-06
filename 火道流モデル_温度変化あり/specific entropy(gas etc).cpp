#include <stdio.h>
#include <math.h>
#include "specific entropy(gas etc).h"

double entropygas(double P, double T,  double ck) {
	//ek=e_, co= Éœ0, Y=É¡,
	double R = 461.88888;
	double c0 = pow(10,8)/(R*(373));
	double T0 = 373;
	double Y = 1.324;
	double sk = 188.825;
	double cv = 1571;


	double s = cv * log((T / T0)*pow((c0 / ck), Y - 1))+sk;
	return s;
}