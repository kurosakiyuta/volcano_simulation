#include <stdio.h>
#include <math.h>


double MPF_mono(double rp) {

	double b = 1.08;//ariance of the logÅ]Gaussian function
	double f = 0.656;//the maximum packing fraction for particles with rp = 1

	double a;
	double aa = log10(rp);
	a = exp(-(aa*aa / (2 * b*b)));


	double maxf = f * a;
	return maxf;
}