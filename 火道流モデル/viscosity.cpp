#include <stdio.h>
#include <math.h>
#include "viscosity.h"
/**********************************************/
/************îSê´åWêîÇÃëgê¨àÀë∂ê´**************/
/**mainÇ≈íËã`Ç∑Ç◊Ç´ïœêîÇÕÅAäeëgê¨î‰Ç∆AÇ∆TÇ∆Cw**/
/**********************************************/
#define SiO2 0
#define TiO2 1
#define Al2O3 2
#define Fe2O3 3
#define FeO 4
#define MnO 5
#define MgO 6
#define CaO 7
#define Na2O 8
#define K2O 9
#define P2O5 10
#define F2O 0
#define H2O 11

double viscosity(double comp_melt_mol[], double T) {
	double A = -4.55;//îSê´åWêîÇÃíËêî
	double V = comp_melt_mol[H2O] + F2O;
	double TA = comp_melt_mol[TiO2] + comp_melt_mol[Al2O3];
	double FM = comp_melt_mol[FeO] + comp_melt_mol[MnO] + comp_melt_mol[MgO];
	double NK = comp_melt_mol[Na2O] + comp_melt_mol[K2O];
	double b[10] = { 159.6,  -173.3,  72.1 ,  75.7  , /*5*/  -39.0 ,  -84.1  ,  141.5  ,  -2.43  ,  -0.91  ,   17.6 };
	double c[7] =  { 2.75,  15.7  ,  8.3  ,  10.2  ,  -12.3  ,  -99.5  ,  0.3 };

	double Mb[10] = { /*0*/comp_melt_mol[SiO2] + comp_melt_mol[TiO2],/*1*/comp_melt_mol[Al2O3], /*2*/comp_melt_mol[FeO] + comp_melt_mol[MnO] + comp_melt_mol[P2O5],
		/*3*/ comp_melt_mol[MgO],/*4*/comp_melt_mol[CaO],/*5*/ comp_melt_mol[Na2O] + V, /*6*/V + log(1 + comp_melt_mol[H2O]),
						/*7*/(comp_melt_mol[SiO2]+ comp_melt_mol[TiO2])*(FM),/*8*/(comp_melt_mol[SiO2] + TA + comp_melt_mol[P2O5])*(NK+ comp_melt_mol[H2O]),
										/*9*/ comp_melt_mol[Al2O3]*NK };

	double Mc[7] = {/*0*/comp_melt_mol[SiO2], /*1*/TA,/*2*/FM,/*3*/comp_melt_mol[CaO],/*4*/NK,
		/*5*/log(1 + V),/*6*/(comp_melt_mol[Al2O3] + FM + comp_melt_mol[CaO] - comp_melt_mol[P2O5])*(NK + V) };

	double B = (b[0] * Mb[0]) + (b[1] * Mb[1]) + (b[2] * Mb[2]) + (b[3] * Mb[3]) + (b[4] * Mb[4]) + (b[5] * Mb[5]) + (b[6] * Mb[6]) +
				(b[7] * Mb[7]) + (b[8] * Mb[8]) + (b[9] * Mb[9]);

	double C = (c[0] * Mc[0]) + (c[1] * Mc[1]) + (c[2] * Mc[2]) + (c[3] * Mc[3]) + (c[4] * Mc[4]) + (c[5] * Mc[5])
					+ (c[6] * Mc[6]);

	double vis = pow(10,A + B / (T - C));	
	return vis;

}
	/*float BB[10];
	for (int i = 0; i <= 9; i++) {
		BB[i] = b[i] * Mb[i];
		printf("BB[%d}%f\n",i, BB[i]);
	}

	printf("%f\n%f\n", B, C);*/

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