#include <stdio.h>
#include <math.h>
#include "MPF(mono).h"


double MPF_bimo(double rp,double b[],double bi[]) {
	double MPF_mono_small = MPF_mono(rp);
	double MPF_mono_large = MPF_mono(1);
	double λ = 3;
	double b1 = 2;
	double b2 = 0.5;
	double bt = b[0] + b[1] + b[2] + b[3] + b[4] + b[5] + b[6];
	//何をマイクロライトとして考えるかは注意が必要
	double bf = ((b[0] + b[1] + b[2] + b[3] + b[4] + b[5] + b[6] - (bi[0] + bi[1] + bi[2] + bi[3] + bi[4] + bi[5] + bi[6]))) / bt;
	if (bt == 0) { bf = 0; };
	if (bf < 0) { bf = 0; };

	double bbb = 1 - pow((1 - (1 / λ)), 1.79);
	double bb = pow(bbb, b1);

	double aa= 1 - pow((1 - (1 / λ)), 1.13);
	double a = pow(aa, b2);

	double MPF_small;
	MPF_small=MPF_mono_small/(1-(1-bf)*(1-MPF_mono_small+bb*(MPF_mono_small-1)));

	double MPF_large = MPF_mono_large / (1 - bf * (1 - a));

	double MPF = MPF_large;

	if (MPF_large > MPF_small) {

		MPF = MPF_small;

	}

	return MPF;
}