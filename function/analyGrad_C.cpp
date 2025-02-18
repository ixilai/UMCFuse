#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "mex.h"
#include <stdio.h>

#define  ML    80
#define  NDP  722  
int W, H;
int *rpos;
double *AS;
double *SA;
void analyGrad(double *sI, int *ndg, int wi);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

	const mxArray *sI_m = prhs[0], *ndg_m = prhs[1], *rpos_m = prhs[2];
	int wi = (int)mxGetScalar(prhs[3]);
	int rh = ML, rw = NDP;
	W = mxGetDimensions(sI_m)[1];
	H = mxGetDimensions(sI_m)[0];

	double *sI = (double *)malloc(sizeof(double)*W*H);
	int *ndg = (int *)malloc(sizeof(int)*W*H);
	rpos = (int *)malloc(sizeof(int)*rw*rh);
	double *sI_mp = (double *)mxGetPr(sI_m); 
	double *ndg_mp = (double *)mxGetPr(ndg_m), *rpos_mp = (double *)mxGetPr(rpos_m);

	int i, j;
	for(j=0; j<H; j++)
	{
		for(i=0; i<W; i++)
		{
			*(sI+j*W+i) = *(sI_mp+i*H+j);     *(ndg+j*W+i) = (int)*(ndg_mp+i*H+j);
		}
	}

	for(j=0; j<rh; j++)
	{
		for(i=0; i<rw; i++)
		{
				*(rpos+j*rw+i) = (int)*(rpos_mp+i*rh+j);
		}
	}

	AS = (double *)malloc(sizeof(double)*W*H);
	SA = (double *)malloc(sizeof(double)*W*H);

	analyGrad(sI, ndg, wi);


	plhs[0] = mxCreateDoubleMatrix((mwSize)H, (mwSize)W, mxREAL);
	plhs[1] = mxCreateDoubleMatrix((mwSize)H, (mwSize)W, mxREAL);
	double *AS_mp = (double *)mxGetPr(plhs[0]), *SA_mp = (double *)mxGetPr(plhs[1]);

	for(j=0; j<H; j++)
	{
		for(i=0; i<W; i++)
		{
			*(AS_mp + i*H + j) = *(AS + j*W + i);  *(SA_mp + i*H + j) = *(SA + j*W + i);
		}
	}

	free(sI); free(ndg); free(rpos); free(AS); free(SA);

}


void analyGrad(double *sI, int *ndg, int wi)
{
	int rw = 2*wi+2, rh = NDP; 
	int rhh = rh/2;
	int *rp = (int *)malloc(sizeof(int)*rw*rh);
	memset(rp, 0, sizeof(int)*rw*rh);
	int i, j;
	//for(j=0; j<722; j++)
	//{
	//	for(i=0; i<wi; i++)
	//	{
	//		if (j<361)
	//			*(rp + j*rw + i) = -(*(rpos + (wi-1-i)*722 + j));
	//		else
	//			*(rp + j*rw + i) = -(*(rpos + (wi-1-i)*722 + j))*W;
	//	}
	//}
	for(j=0; j<rh; j++)
	{
		for(i=0; i<wi+1; i++)
		{
			if (j<rhh)
				*(rp + j*rw + wi+1 + i) = *(rpos + i*rh + j);
			else
				*(rp + j*rw + wi+1 + i) = (*(rpos + i*rh + j))*W;
		}
	}
	int wid = 2*wi;
	for(j=0; j<rh; j++)
	{
		for(i=0; i<wi; i++)
		{
				*(rp + j*rw + i) = -(*(rp + j*rw + wid - i));
		}
	}

	//FILE *fp = fopen("D:\\t1.txt", "w");
	//for(j=0; j<rh; j++)
	//{
	//	fprintf(fp, "%d:  ", j);
	//	for(i=0; i<rw; i++)
	//	{
	//			fprintf(fp, "%d, ", *(rp+j*rw+i)); 
	//	}
	//	fprintf(fp, "\n");
	//}
	//fclose(fp);

//	clock_t m_begin = clock();  // time measurement;

	int  op = 0;
	int  rwhh = rw*rhh;
	double *pAS = AS, *pSA = SA;
	for(j=0; j<H; j++)
	{
		for(i=0; i<W; i++)
		{
			int dg = *ndg++; 
			int *rpdx = rp + dg*rw;
			int *rpdy = rpdx + rwhh;
			int *rpdx1 = rpdx, *rpdy1 = rpdy;
			double v1 = -1;
			double v2, vp;
			double sa = 0;
			double as = 0;
			if ( j<wi+1 || j>H-wi-2 || i<wi+1 || i>W-wi-2 )
			{
				for(int k=0; k<rw; k++)
				{
					int pdx = *rpdx++;
					int pdy = *rpdy++;
					int xd = pdx;
					if (xd+i<0 || xd+i>W-1) continue;
					else
					{
						int yd = pdy/W;
						if (yd+j<0 || yd+j>H-1) continue;
					}

					if(v1<0) {v1 = *(sI + op + pdy + pdx ); v2 = v1; vp = v1; }
					else
					{
						v2 = *(sI + op + pdy + pdx); 
						sa += fabs(v2-vp);
						vp = v2;
					}                
				}
				if (v1>=0)  as = fabs(v2-v1);
			}
			else
			{
				for(int k=0; k<rw; k++)
				{
					int pdx = *rpdx1++;
					int pdy = *rpdy1++;
					if(v1<0) {v1 = *(sI + op + pdy + pdx); v2 = v1; vp = v1; }
					else
					{
						v2 = *(sI + op + pdy + pdx); 
						sa += fabs(v2-vp);
						vp = v2;
					}
				}
				as = fabs(v2-v1);
			}
			op++;
			*pAS++ = as;  *pSA++ = sa; 
		}
	}

//	mexPrintf("Elapsed time is %f seconds.\n", double(clock()-m_begin)/CLOCKS_PER_SEC);

	free(rp);
}