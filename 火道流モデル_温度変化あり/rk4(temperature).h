#pragma once
void rk4T(double y[], double dydx[], int n, double x, double h, double yout[], double b[], double P[], double dPdz[], double xc[], double xd[], double dxcdz[], double dxddz[],
	void(*devris)(double, double[], double[], double[], double[], double[], double[], double[], double[]));