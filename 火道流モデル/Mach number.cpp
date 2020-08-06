#include <stdio.h>
#include <math.h>
#include "Mach number.h"
#include "speed of sound.h"
#include "velocity.h"

double machnumber(double P, double Qi, double T, double tw[], double watersuminitial, double c[], double b[],double xc[], double r) {

	double a = soundspeed(P, T, watersuminitial, tw, c, b,xc);
	double M = velosity(P, Qi, T, tw, watersuminitial,c, b, r)/a;
	return M;
}