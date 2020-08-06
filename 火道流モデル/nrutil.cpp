#include <stdio.h>
#include <stddef.h>
#include<stdlib.h>
#include "/Users/DELL/source/C++/数値計算C++/数値計算C++/nrutil.h"
#define NR_END 1
#define FREE_AGE char*

void nrerror(const char error_text[])/* Numerical Recipes standard error handler */
{
	fprintf(stderr, "Nemerical Recipes run-time error...\n");
	fprintf(stderr, "%s\n",error_text);
	fprintf(stderr, "...now exitng to system...\n");
	exit(1);
}/***************************************不確定要素*******************************************************/
/***********************************//*これを使うのは配列の要素数を関数ファイル内で定義せずともできるため*/
/* nlからnhまでのベクトルの割り当て*//*配列の要素数を変数のまま宣言できる*********************************/
/***********************************//*なお、その要素の幅は、0~n-1であることに注意************************/
/***********************************//*ニューメリカルの教科書では、１～ｎで行っている*********************/
/*********************************************************************************************************/
	double *vector(long nl, long nh)
	{
		double *v;
        /*ニューメリカルのベクトル返還とは変えてるから気を付けて！！*/
		v = (double *)malloc((size_t) ((nh - nl + 1 /* + NR_END*/) * sizeof(double)));
		if(!v) nrerror("allocation failure in vector()");
		return v - nl /* + NR_END*/ ;
	}

	void free_vector(double *v, long nl, long nh)/* bector()で割り当てられたベクトルの解放　*/
	{

		free((FREE_AGE)(v + nl));
	}

