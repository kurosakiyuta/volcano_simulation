#include <stdio.h>
#include "/Users/DELL/source/C++/数値計算C++/数値計算C++/nrutil.h"
#include "rk4(temperature).h"


void rk4T(double y[], double dydx[], int n, double x, double h, double yout[], double b[],double P[],double dPdz[],double xc[],double xd[],double dxcdz[], double dxddz[],
	void(*devris)(double, double[], double[], double[],double[] ,double[],double[],double[],double[] )) { //devris=(z, xc,xd,dxcdz,dxddz, P,dPdz, T, dTdz)
	int i;
	double xh, hh, h6, *dym, *dyt, *yt;

	dyt = vector(0, n - 1);
	dym = vector(0, n - 1);
	yt = vector(0, n - 1);


	hh = h * 0.5;
	h6 = h / 6;
	xh = x + hh;


	for (i = 0; i <= n - 1; i++) {
		yt[i] = y[i] + hh * dydx[i];
	}



	(*devris)(xh, xc,xd,dxcdz,dxddz,P,dPdz, yt, dyt);

	for (i = 0; i <= n - 1; i++) {
		yt[i] = y[i] + hh * dyt[i];
	}

	(*devris)(xh, xc, xd, dxcdz, dxddz,P,dPdz, yt, dym);



	for (i = 0; i <= n - 1; i++) {
		yt[i] = y[i] + h * dym[i];
		dym[i] += dyt[i];
	}

	(*devris)(x + h, xc, xd, dxcdz, dxddz, P,dPdz, yt, dyt);



	for (i = 0; i <= n - 1; i++)
	{
		yout[i] = y[i] + h6 * (dydx[i] + dyt[i] + 2.0*dym[i]);
	}

	free_vector(yt, 0, n - 1);
	free_vector(dyt, 0, n - 1);
	free_vector(dym, 0, n - 1);

}

