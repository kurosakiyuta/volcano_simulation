#include <stdio.h>
#include <math.h>
#include "xeq.h"
#include "density(liquid).h"
#include "equilibrium.h"

void crystalequilibrium(double P, double beq[], double xeq[], double tpl[], double tpy[], double to[], double tq[], double ts[], double tcl[], double tor[], double c[]) {

	xeq[0] = xequilibrium(P, c[0], tpl);
	xeq[1] = xequilibrium(P, c[1], tpy);
	xeq[2] = xequilibrium(P, c[2], to);
	xeq[3] = xequilibrium(P, c[3], tq);
	xeq[4] = xequilibrium(P, c[4], ts);
	xeq[5] = xequilibrium(P, c[4], tcl);
	xeq[6] = xequilibrium(P, c[4], tor);


	double cleq = densityliquidwt(xeq, c);

	beq[0] = (cleq / c[0])*xeq[0];
	beq[1] = (cleq / c[1])*xeq[1];
	beq[2] = (cleq / c[2])*xeq[2];
	beq[3] = (cleq / c[3])*xeq[3];
	beq[4] = (cleq / c[4])*xeq[4];
	beq[5] = (cleq / c[5])*xeq[5];
	beq[6] = (cleq / c[6])*xeq[6];

	//if (beq[1] > 0.29) { beq[1] = 0.29; };


	if (beq[0] < 0) { beq[0] = 0; };
	if (beq[1] < 0) { beq[1] = 0; };
	if (beq[2] < 0) { beq[2] = 0; };
	if (beq[3] < 0) { beq[3] = 0; };
	if (beq[4] < 0) { beq[4] = 0; };
	if (beq[5] < 0) { beq[5] = 0; };
	if (beq[6] < 0) { beq[6] = 0; };

}

