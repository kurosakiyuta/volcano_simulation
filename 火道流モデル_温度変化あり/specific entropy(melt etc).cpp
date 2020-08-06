#include <stdio.h>
#include <math.h>
#include "specific entropy(melt etc).h"

double entropymelt(double P, double T, double sk, double cv, double c0, double T0, double Y, double ck) {
	//ek=e_, co= Éœ0, Y=É¡,

	double s = sk + cv * log((T / T0)*pow((c0 / ck), Y - 1));
	return s;
}