#pragma once
void rk4S(double y[], double dydx[], int n, double x, double h, double yout[], double b[], double T[], double xd[], double xc[], double P[],
	void(*devris)(double, double[], double[], double[], double[], double[]));