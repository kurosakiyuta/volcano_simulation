#pragma once

void rk4(double y[], double dydx[], int n, double x, double h, double yout[], double b[], double T[], double xd[], double xc[],
	void(*devris)(double, double[], double[], double[], double[], double[]));