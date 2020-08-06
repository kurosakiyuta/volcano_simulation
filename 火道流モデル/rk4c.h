#pragma once

void rk4(double y[], double dydx[], int n, double x, double h, double yout[],double b[],double watersuminitial,double vis,double bi[],
	void(*devris)(double,double [], double[], double[],double,double, double []));