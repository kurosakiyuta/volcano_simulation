#pragma once

void rk4c(double y[], double dydx[], int nc, int n, double x, double h, double P[], double yout[], double T[], double xd[],
	void(*devris)(double, double[], double[], double[], double[], double[]));