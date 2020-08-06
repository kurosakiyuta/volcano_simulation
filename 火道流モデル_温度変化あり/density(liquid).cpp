#include <stdio.h>
#include <math.h>
#include "density(liquid).h"

double densityliquid( double P,double xc[],double xd[], double cpl, double cpy, double co, double cm, double cd) {
	//double xd = (1-X/100)*(Ci - X)/100;
	//double A = 1 - xd + xd / cd;
	//double cl = (cm *(1 - bpl - bpy - bo) + cpl * bpl + cpy * bpy + co * bo)/A;
	double xm = 1 - xd[0] - xc[0] - xc[1] - xc[2];
	double a = (xm / cm) + (xd[0] / cd) + (xc[0] / cpl) + (xc[1] / cpy) + (xc[2] / co);
	double cl = 1 / a;
	return cl;

}