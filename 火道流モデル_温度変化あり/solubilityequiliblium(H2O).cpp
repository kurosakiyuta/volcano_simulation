#include <stdio.h>
#include <math.h>
#include "solubilityH2O.h"

/**********************************************/
/*�[���Ɖ��x�̕ω��ɔ����}���g�����̐��̗n��x*/
/********�܂�A�}���g�����̐���wt%***********/
/**********************************************/

double solubilityw(double P, double T, double AI) {
	double Cw;
	double SP;
	double MP = P / 1000000;

	SP = sqrt(MP);

	Cw = (-0.231 + 651.1 / T)*SP + (0.03424 - 32.57 / T + 0.2447*AI)*MP;//MPa

	return Cw;//[wt%]
}


