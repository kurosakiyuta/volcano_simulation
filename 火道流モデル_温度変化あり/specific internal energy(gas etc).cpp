#include <stdio.h>
#include <math.h>
#include "specific internal energy(gas etc).h"

double energygas(double T, double ek, double cv) {
	//ek=e_, co= ��0, Y=��,

	double e = cv * T + ek;
	return e;
}