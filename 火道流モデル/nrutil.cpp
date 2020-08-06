#include <stdio.h>
#include <stddef.h>
#include<stdlib.h>
#include "/Users/DELL/source/C++/���l�v�ZC++/���l�v�ZC++/nrutil.h"
#define NR_END 1
#define FREE_AGE char*

void nrerror(const char error_text[])/* Numerical Recipes standard error handler */
{
	fprintf(stderr, "Nemerical Recipes run-time error...\n");
	fprintf(stderr, "%s\n",error_text);
	fprintf(stderr, "...now exitng to system...\n");
	exit(1);
}/***************************************�s�m��v�f*******************************************************/
/***********************************//*������g���͔̂z��̗v�f�����֐��t�@�C�����Œ�`�����Ƃ��ł��邽��*/
/* nl����nh�܂ł̃x�N�g���̊��蓖��*//*�z��̗v�f����ϐ��̂܂ܐ錾�ł���*********************************/
/***********************************//*�Ȃ��A���̗v�f�̕��́A0~n-1�ł��邱�Ƃɒ���************************/
/***********************************//*�j���[�����J���̋��ȏ��ł́A�P�`���ōs���Ă���*********************/
/*********************************************************************************************************/
	double *vector(long nl, long nh)
	{
		double *v;
        /*�j���[�����J���̃x�N�g���Ԋ҂Ƃ͕ς��Ă邩��C��t���āI�I*/
		v = (double *)malloc((size_t) ((nh - nl + 1 /* + NR_END*/) * sizeof(double)));
		if(!v) nrerror("allocation failure in vector()");
		return v - nl /* + NR_END*/ ;
	}

	void free_vector(double *v, long nl, long nh)/* bector()�Ŋ��蓖�Ă�ꂽ�x�N�g���̉���@*/
	{

		free((FREE_AGE)(v + nl));
	}

