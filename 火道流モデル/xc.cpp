#include <stdio.h>
#include <math.h>
#include "density(liquid).h"
#include "xc.h"

void crystalwt(double b[],double xc[] ,double c[]) {

	double cl = densityliquidvol(b, c);
	xc[0] = (c[0] / cl)*b[0];
	xc[1] = (c[1] / cl)*b[1];
	xc[2] = (c[2] / cl)*b[2];
	xc[3] = (c[3] / cl)*b[3];
	xc[4] = (c[4] / cl)*b[4];
	xc[5] = (c[5] / cl)*b[5];
	xc[6] = (c[6] / cl)*b[6];

}

void crystalwttovol(double b[], double xc[], double c[]){
	double cl = densityliquidwt(xc, c);
	printf("cl=%f", cl);
	printf("xc0=%f\n", xc[1]);
	b[0] = (cl / c[0])*xc[0];
	b[1] = (cl / c[0])*xc[1];
	b[2] = (cl / c[0])*xc[2];
	b[3] = (cl / c[0])*xc[3];
	b[4] = (cl / c[0])*xc[4];
	b[5] = (cl / c[0])*xc[5];
	b[6] = (cl / c[0])*xc[6];
}