#include <stdio.h>
#include <math.h>
#include "viscosity.h"
/**********************************************/
/************粘性係数の組成依存性**************/
/**mainで定義すべき変数は、各組成比とAとTとCw**/
/**********************************************/

double viscosity(double SiO2, double TiO2, double Al2O3, double FeO, double MnO, double P2O5, double MgO,
	double CaO, double Na2O, double Cw, double F2O, double K2O, double T, double A) {

	double V = Cw + F2O;
	double TA = TiO2 + Al2O3;
	double FM = FeO + MnO + MgO;
	double NK = Na2O + K2O;
	double b[10] = { 159.6,  -173.3,  72.1 ,  75.7  , /*5*/  -39.0 ,  -84.1  ,  141.5  ,  -2.43  ,  -0.91  ,   17.6 };
	double c[7] =  { 2.75,  15.7  ,  8.3  ,  10.2  ,  -12.3  ,  -99.5  ,  0.3 };
	double Mb[10] = { /*0*/SiO2 + TiO2,/*1*/Al2O3, /*2*/FeO + MnO + P2O5,/*3*/ MgO, /*4*/CaO,/*5*/ Na2O + V, /*6*/V + log(1 + Cw),
						/*7*/(SiO2+TiO2)*(FM),/*8*/(SiO2 + TA + P2O5)*(NK+Cw), /*9*/ Al2O3*NK };
	double Mc[7] = {/*0*/SiO2, /*1*/TA,/*2*/FM,/*3*/CaO,/*4*/NK,/*5*/log(1 + V),/*6*/(Al2O3 + FM + CaO - P2O5)*(NK + V) };

	double B = (b[0] * Mb[0]) + (b[1] * Mb[1]) + (b[2] * Mb[2]) + (b[3] * Mb[3]) + (b[4] * Mb[4]) + (b[5] * Mb[5]) + (b[6] * Mb[6]) +
				(b[7] * Mb[7]) + (b[8] * Mb[8]) + (b[9] * Mb[9]);

	double C = (c[0] * Mc[0]) + (c[1] * Mc[1]) + (c[2] * Mc[2]) + (c[3] * Mc[3]) + (c[4] * Mc[4]) + (c[5] * Mc[5])
					+ (c[6] * Mc[6]);

	double vis = pow(10,A + B / (T - C));
	/*float BB[10];
	for (int i = 0; i <= 9; i++) {
		BB[i] = b[i] * Mb[i];
		printf("BB[%d}%f\n",i, BB[i]);
	}

	printf("%f\n%f\n", B, C);*/
	return vis;
}

/*int main() {
	
	float SiO2=62.38;
	float TiO2 = 0.41;
	float Al2O3 = 11.79;
	float FeO=0.03;
	float MnO = 0.02;
	float P2O5 = 0.05;
	float MgO = 4.8;
	float CaO=9.73;
	float Na2O = 3.41;

	float F2O = 0;
	float K2O = 0.59;
	float T = 1273;
	float A = -4.55;
			float Cw;
		Cw = 0;
	for (int i = 0; Cw < 7; i++) {

		float a = viscosity(SiO2, TiO2, Al2O3, FeO, MnO, P2O5, MgO, CaO, Na2O, Cw, F2O, K2O, T, A);
		printf("%f      %f\n", Cw,a);
		Cw += 0.1;
	}


}*/