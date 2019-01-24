#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "function.h"
#include "common.h"


//================alloc and free matrix
float ***Creat3dArray(int m, int n, int k)
{
	float ***tt=(float ***)malloc(sizeof(float **)*m);
	for (int i=0;i<m;i++)
	{
		tt[i]=(float **)malloc(sizeof(float *)*n);
		for (int j=0;j<n;j++)
		{
			tt[i][j]=(float*)malloc(sizeof(float)*k);
		}
	}
	return tt;
}
void free3dArray(float ***tt, int m, int n, int k)
{
	if (tt!=NULL)
	{
		for (int i=0;i<m;i++)
		{
			for (int j=0;j<n;j++)
			{
				free((tt[i][j]));
			}
			free(tt[i]);
		}
		free(tt);
		tt=NULL;
	}
}
float **Creat2dArray(int m, int n)
{
	float **tt=(float**)malloc(sizeof(float *)*m);
	for (int i=0;i<m;i++)
	{
		tt[i]=(float*)malloc(sizeof(float)*n);
	}
	return tt;
}
void free2dArray(float **tt, int m, int n)
{
	if (tt!=NULL)
	{
		for (int i=0;i<m;i++)
		{
			free(tt[i]);
		}
		free(tt);
		tt=NULL;
	}
}
int **Creat2dArray_int(int m, int n)
{
	int **tt=(int**)malloc(sizeof(int *)*m);
	for (int i=0;i<m;i++)
	{
		tt[i]=(int*)malloc(sizeof(int)*n);
	}
	return tt;
}
void free2dArray_int(int **tt, int m, int n)
{
	if (tt!=NULL)
	{
		for (int i=0;i<m;i++)
		{
			free(tt[i]);
		}
		free(tt);
		tt=NULL;
	}
}
//================initial matrix(one or two dimension)
void fmemset1(float *p, int len)
{
	for(int j=0;j<len;j++)
		p[j]=0.0;
}
void fmemset2(float **p, int nz, int nx)
{
	for(int iz=0;iz<nz;iz++)
		for(int ix=0;ix<nx;ix++)
			p[iz][ix]=0.0;
}
void fmemset1v(float *p, int len, float v)
{
	for(int j=0;j<len;j++)
		p[j]=v;
}

void fmemset2v(float **p, int nz, int nx, float v)
{
	for(int iz=0;iz<nz;iz++)
		for(int ix=0;ix<nx;ix++)
			p[iz][ix]=v;
}
void fmemset1vp(float *p, int len, float *vp)
{
	for(int j=0;j<len;j++)
		p[j]=vp[j];
}
void fmemset2vp(float **p, int nz, int nx, float **vp)
{
	for(int iz=0;iz<nz;iz++)
		for(int ix=0;ix<nx;ix++)
			p[iz][ix]=vp[iz][ix];
}
//================max or min value of matrix
float absMaxval2_AB(float **A, float **B, int nz, int nx)
{
	float data = fabs(A[0][0]*B[0][0]);
	for (int iz=0; iz<nz; iz++)
	{
		for (int ix=0; ix<nx; ix++)
		{
			data = MAX(data,fabs(sqrt(A[iz][ix])*sqrt(B[iz][ix])));
		}
	}
	return data;
}
float Maxval1(float *v, int n)
{
	float a=v[0];
	for (int i=0;i<n;i++)
	{
		a=MAX(a,v[i]);
	}
	return a;
}
float Minval1(float *v, int n)
{
	float a=v[0];
	for (int i=0;i<n;i++)
	{
		a=MIN(a,v[i]);
	}
	return a;
}
float Maxval2(float **v, int nz, int nx)
{
	float data = v[0][0];
	for (int iz=0; iz<nz; iz++)
	{
		for (int ix=0; ix<nx; ix++)
		{
			data = MAX(data,v[iz][ix]);
		}
	}
	return data;
}
float Minval2(float **v, int nz, int nx)
{
	float data = v[0][0];
	for (int iz=0; iz<nz; iz++)
	{
		for (int ix=0; ix<nx; ix++)
		{
			data = MIN(data,v[iz][ix]);
		}
	}
	return data;
}
float absMaxval1(float *v, int n)
{
	float a=fabs(v[0]);
	for (int i=0;i<n;i++)
	{
		a=MAX(a,fabs(v[i]));
	}
	return a;
}
float absMinval1(float *v, int n)
{
	float a=fabs(v[0]);
	for (int i=0;i<n;i++)
	{
		a=MIN(a,fabs(v[i]));
	}
	return a;
}
float absMaxval2(float **v, int nz, int nx)
{
	float data = fabs(v[0][0]);
	for (int iz=0; iz<nz; iz++)
	{
		for (int ix=0; ix<nx; ix++)
		{
			data = MAX(data,fabs(v[iz][ix]));
		}
	}
	return data;
}
float absMinval2(float **v, int nz, int nx)
{
	float data = fabs(v[0][0]);
	for (int iz=0; iz<nz; iz++)
	{
		for (int ix=0; ix<nx; ix++)
		{
			data = MIN(data,fabs(v[iz][ix]));
		}
	}
	return data;
}
//================sum
float sum1(float *data, int n)
{
	float s=0;
	for (int i=0;i<n;i++)
	{
		s=s+data[i];
	}
	return s;
}
float sum2(float **data, int nx, int nz)
{
	float s=0;
	for (int i=0;i<nz;i++)
	{
		for (int j=0;j<nx;j++)
			s=s+data[i][j];
	}
	return s;
}
float sum1abs(float *data, int n)
{
	float s=0;
	for (int i=0;i<n;i++)
	{
		s=s+fabs(data[i]);
	}
	return s;
}
float sum2abs(float **data, int nx, int nz)
{
	float s=0;
	for (int i=0;i<nz;i++)
	{
		for (int j=0;j<nx;j++)
			s=s+fabs(data[i][j]);
	}
	return s;
}
void arrayabs(float *data1, float *data, int n)
{
	for (int i=0;i<n;i++)
	{
		data1[i]=fabs(data[i]);
	}
}
//============================file operation
int exist(char *buff)
{
	FILE *fp;
	if ((fp=fopen(buff,"r")) == NULL)
		return 0;
	else{
		fclose(fp);
		return 1;}
}
long filesize(FILE *fp)
{
	long length;
	fseek(fp,0,2);
	length = ftell(fp);
	fseek(fp,0,0);
	return length;
}
