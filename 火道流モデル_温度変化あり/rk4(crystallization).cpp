#include <stdio.h>
#include "/Users/DELL/source/C++/数値計算C++/数値計算C++/nrutil.h"
#include "rk4(crystallization).h"


void rk4c(double y[], double dydx[], int nc,int n, double x, double h,double P[], double yout[],double T[],double xd[],
	void(*devris)(double,double[],double[], double[] ,double[], double[])) { //devris=(z, P,xd,T, xc, dxcdz)
	int i;
	double xh, hh, h6, *dym, *dyt, *yt;

	dyt = vector(0, nc - 1);
	dym = vector(0, nc - 1);
	yt = vector(0, nc - 1);

	hh = h * 0.5;
	h6 = h / 6;
	xh = x + hh;


	for (i = 0; i <= nc-1; i++) {
		yt[i] = y[i] + hh * dydx[i];
	}

	(*devris)(xh,P,xd,T, yt, dyt);

	for (i = 0; i <= nc-1; i++) {
		yt[i] = y[i] + hh * dyt[i];
	}

	(*devris)(xh, P,xd,T, yt, dym);


	for (i = 0; i <= nc-1; i++) {
		yt[i] = y[i] + h * dym[i];
		dym[i] += dyt[i];
	}

	(*devris)(x + h, P,xd, T,yt, dyt);



	for (i = 0; i <= nc-1; i++)
	{
		yout[i] = y[i] + h6 * (dydx[i] + dyt[i] + 2.0*dym[i]);
	}

	free_vector(yt, 0, n - 1);
	free_vector(dyt, 0, n - 1);
	free_vector(dym, 0, n - 1);
}




