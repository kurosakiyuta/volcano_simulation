#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "solubilityH2O.h"
#include "xeq.h"

double xequilibrium(double P,double cj,double t[]) {
	double Pbar = P / pow(10,8);

	double xeq = (t[0]
		+ t[1] * Pbar
		+ t[2] * Pbar*Pbar
		+ t[3] * Pbar*Pbar*Pbar
		+ t[4] * Pbar*Pbar*Pbar*Pbar
		+ t[5] * Pbar*Pbar*Pbar*Pbar*Pbar
		+ t[6] * Pbar*Pbar*Pbar*Pbar*Pbar*Pbar
		)/100;


	if (xeq < 0) { xeq = 0; };
//	printf("xeq=%f\n", xeq);
	return xeq;
}



/*		t[0] * Pbar*Pbar
		+ t[1] * Td*Td
		+ t[2] * Cw*Cw
		+ t[3] * Pbar*Td
		+ t[4] * Td*Cw
		+ t[5] * Cw*Pbar
		+ t[6] * Pbar
		+ t[7] * Td
		+ t[8] * Cw
		+ t[9]
		+ t[10] * Pbar * Pbar*Pbar
		+ t[11] * Cw*Cw*Cw
		+ t[12] * Td*Td*Td
		+ t[13] * Pbar*Pbar*Cw
		+ t[14] * Pbar*Cw*Cw
		+ t[15] * Pbar*Cw*Td
		+ t[16] * Pbar*Td*Td
		+ t[17] * Cw*Cw*Td
		+ t[18] * Cw*Td*Td
		+ t[13] * Pbar*Pbar*Pbar*Pbar
		+ t[14] * Cw*Cw*Cw*Cw
		+ t[15] * Td*Td*Td*Td*/