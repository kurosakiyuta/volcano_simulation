#include <stdio.h>
#include <math.h>
#include "Mach number.h"
#include "speed of sound.h"
#include "velocity.h"

double machnumber(double P, double Qi, double T, double xd[], double Ci, double cl,double cg,double r) {

	double a = soundspeed(P, cl, T, Ci, xd, cg);
	double M = velosity(P, Qi, T, xd, Ci,cl, cg,r)/a;
	
	return M;
}