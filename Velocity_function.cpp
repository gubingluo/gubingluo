#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>
#include "function.h"
#include "common.h"


void pmlvelsmooth1d(float *vp, int nx, int nz, int npml)
{
	int ix,iz;
	float **vvp;
	vvp=Creat2dArray(nz,nx);
	for (ix=0;ix<npml;ix++)
	{
		for (iz=1;iz<nz-1;iz++)
		{
			vp[iz*nx + ix]=vp[iz*nx + 2*npml-ix-1];
			vp[iz*nx + nx-npml+ix]=vp[iz*nx + nx-npml-ix-1];					
		}
	}
	for (iz=0;iz<npml;iz++)
	{
		for (ix=1;ix<nx-1;ix++)
		{
			vp[iz*nx + ix]=vp[(2*npml-iz-1)*nx + ix];
			vp[(nz-npml+iz)*nx + ix]=vp[(nz-npml-iz-1)*nx + ix];			
		}
	}
	float Error=1.0;
	int Num=0;
	float MaxIter=1000;
	for (ix=0;ix<nx;ix++)
	{
		for (iz=0;iz<nz;iz++)
		{
			vvp[iz][ix]=vp[iz*nx + ix];
		}
	}
	while (Num <= MaxIter && Error >0.0005)
	{
		for (ix=1;ix<npml;ix++)
		{
			for (iz=1;iz<nz-1;iz++)
			{
				vvp[iz][ix]=0.25*(vvp[iz][ix+1]+vvp[iz][ix-1]+
					vvp[iz+1][ix]+vvp[iz-1][ix]);
				
				vvp[iz][ix+nx-npml-1]=0.25*(vvp[iz][ix+nx-npml]+vvp[iz][ix+nx-npml-2]+
					vvp[iz+1][ix+nx-npml-1]+vvp[iz-1][ix+nx-npml-1]);
			}
		}
		for (ix=1;ix<nx-1;ix++)
		{
			for (iz=1;iz<npml;iz++)
			{
				vvp[iz][ix]=0.25*(vvp[iz][ix+1]+vvp[iz][ix-1]+
					vvp[iz+1][ix]+vvp[iz-1][ix]);
				
				vvp[iz+nz-npml-1][ix]=0.25*(vvp[iz+nz-npml-1][ix+1]+vvp[iz+nz-npml-1][ix-1]+
					vvp[iz+nz-npml][ix]+vvp[iz+nz-npml-2][ix]);
			}
		}
		for (iz=1;iz<nz-1;iz++)
		{
			vvp[iz][0]=0.25*(vvp[iz][0]+vvp[iz][1]+vvp[iz+1][0]+vvp[iz-1][0]);
			
			vvp[iz][nx-1]=0.25*(vvp[iz][nx-1]+vvp[iz][nx-2]+vvp[iz+1][nx-1]+vvp[iz-1][nx-1]);
		}
		for (ix=1;ix<nx-1;ix++)
		{
			vvp[0][ix]=0.25*(vvp[0][ix]+vvp[1][ix]+vvp[0][ix-1]+vvp[0][ix+1]);
			
			vvp[nz-1][ix]=0.25*(vvp[nz-1][ix]+vvp[nz-2][ix]+vvp[nz-1][ix-1]+vvp[nz-1][ix+1]);
		}
		vvp[0][0]=0.25*(vvp[0][0]+vvp[1][0]+vvp[0][1]+vvp[1][1]);
		
		vvp[nz-1][0]=0.25*(vvp[nz-1][0]+vvp[nz-2][0]+vvp[nz-1][1]+vvp[nz-2][1]);
		
		vvp[0][nx-1]=0.25*(vvp[0][nx-1]+vvp[0][nx-2]+vvp[1][nx-1]+vvp[1][nx-2]);
		
		vvp[nz-1][nx-1]=0.25*(vvp[nz-1][nx-1]+vvp[nz-2][nx-1]+vvp[nz-1][nx-2]+vvp[nz-2][nx-2]);
		
		float maxvel=-1.0e+5;
		for (iz=0;iz<nz;iz++)
		{
			for (ix=0;ix<nx;ix++)
			{
				maxvel=MAX(fabs(vp[iz*nx + ix]-vvp[iz][ix])/fabs(vp[iz*nx + ix]),maxvel);
			}
		}
		Error=maxvel;
		for (iz=0;iz<nz;iz++)
		{
			for (ix=0;ix<nx;ix++)
			{
				vp[iz*nx + ix]=vvp[iz][ix];
			}
		}
		Num=Num+1;
	}
	free2dArray(vvp,nz,nx);
}
void pmlvelsmooth2d(float **vp, int nx, int nz, int npml)
{
	int ix,iz;
	float **vvp;
	vvp=Creat2dArray(nz,nx);
	for (ix=0;ix<npml;ix++)
	{
		for (iz=1;iz<nz-1;iz++)
		{
			vp[iz][ix]=vp[iz][2*npml-ix-1];
			vp[iz][nx-npml+ix]=vp[iz][nx-npml-ix-1];					
		}
	}
	for (iz=0;iz<npml;iz++)
	{
		for (ix=1;ix<nx-1;ix++)
		{
			vp[iz][ix]=vp[2*npml-iz-1][ix];
			vp[nz-npml+iz][ix]=vp[nz-npml-iz-1][ix];			
		}
	}
	float Error=1.0;
	int Num=0;
	float MaxIter=1000;
	for (ix=0;ix<nx;ix++)
	{
		for (iz=0;iz<nz;iz++)
		{
			vvp[iz][ix]=vp[iz][ix];
		}
	}
	while (Num <= MaxIter && Error >0.0005)
	{
		for (ix=1;ix<npml;ix++)
		{
			for (iz=1;iz<nz-1;iz++)
			{
				vvp[iz][ix]=0.25*(vvp[iz][ix+1]+vvp[iz][ix-1]+
					vvp[iz+1][ix]+vvp[iz-1][ix]);
				
				vvp[iz][ix+nx-npml-1]=0.25*(vvp[iz][ix+nx-npml]+vvp[iz][ix+nx-npml-2]+
					vvp[iz+1][ix+nx-npml-1]+vvp[iz-1][ix+nx-npml-1]);
			}
		}
		for (ix=1;ix<nx-1;ix++)
		{
			for (iz=1;iz<npml;iz++)
			{
				vvp[iz][ix]=0.25*(vvp[iz][ix+1]+vvp[iz][ix-1]+
					vvp[iz+1][ix]+vvp[iz-1][ix]);
				
				vvp[iz+nz-npml-1][ix]=0.25*(vvp[iz+nz-npml-1][ix+1]+vvp[iz+nz-npml-1][ix-1]+
					vvp[iz+nz-npml][ix]+vvp[iz+nz-npml-2][ix]);
			}
		}
		for (iz=1;iz<nz-1;iz++)
		{
			vvp[iz][0]=0.25*(vvp[iz][0]+vvp[iz][1]+vvp[iz+1][0]+vvp[iz-1][0]);
			
			vvp[iz][nx-1]=0.25*(vvp[iz][nx-1]+vvp[iz][nx-2]+vvp[iz+1][nx-1]+vvp[iz-1][nx-1]);
		}
		for (ix=1;ix<nx-1;ix++)
		{
			vvp[0][ix]=0.25*(vvp[0][ix]+vvp[1][ix]+vvp[0][ix-1]+vvp[0][ix+1]);
			
			vvp[nz-1][ix]=0.25*(vvp[nz-1][ix]+vvp[nz-2][ix]+vvp[nz-1][ix-1]+vvp[nz-1][ix+1]);
		}
		vvp[0][0]=0.25*(vvp[0][0]+vvp[1][0]+vvp[0][1]+vvp[1][1]);
		
		vvp[nz-1][0]=0.25*(vvp[nz-1][0]+vvp[nz-2][0]+vvp[nz-1][1]+vvp[nz-2][1]);
		
		vvp[0][nx-1]=0.25*(vvp[0][nx-1]+vvp[0][nx-2]+vvp[1][nx-1]+vvp[1][nx-2]);
		
		vvp[nz-1][nx-1]=0.25*(vvp[nz-1][nx-1]+vvp[nz-2][nx-1]+vvp[nz-1][nx-2]+vvp[nz-2][nx-2]);
		
		float maxvel=-1.0e+5;
		for (iz=0;iz<nz;iz++)
		{
			for (ix=0;ix<nx;ix++)
			{
				maxvel=MAX(fabs(vp[iz][ix]-vvp[iz][ix])/fabs(vp[iz][ix]),maxvel);
			}
		}
		Error=maxvel;
		for (iz=0;iz<nz;iz++)
		{
			for (ix=0;ix<nx;ix++)
			{
				vp[iz][ix]=vvp[iz][ix];
			}
		}
		Num=Num+1;
	}
	free2dArray(vvp,nz,nx);
}
void velsmooth1d(float *vp,int n1,int n2,int nsp)
{
	// nsp 窗口大小
	// n1 行 nz
	// n2 列 nx
	// a 速度
	// b 扩展速度
	int n1e,n2e,i1,i2,i11,i22;
	double PI=3.141592653;
	float **a,**b;
	double a1,b1,dist1,dist2;
	n1e=n1+2*nsp; //数值维度扩展
	n2e=n2+2*nsp; //数组维度扩展
	// 开辟空间  b[n1e][n2e]
	a=(float**)calloc(n1,sizeof(float*));
    for(i1=0;i1<n1;i1++)
	{
		a[i1]=(float*)calloc(n2,sizeof(float));
	}
	b=(float**)calloc(n1e,sizeof(float*));
    for(i1=0;i1<n1e;i1++)
	{
		b[i1]=(float*)calloc(n2e,sizeof(float));
	}
	for(i1=0;i1<n1;i1++)
	{
		for(i2=0;i2<n2;i2++)
		{
			a[i1][i2]=vp[i1*n2+i2];
		}
	}
	//中间
	for(i1=0;i1<n1;i1++)
	{
		for(i2=0;i2<n2;i2++)
		{
			b[i1+nsp][i2+nsp]=a[i1][i2];
		}
	}
	//左边-右边
	for(i1=0;i1<n1;i1++)
	{
		for(i2=0;i2<nsp;i2++)
		{
			b[i1+nsp][i2]=a[i1][0];
			b[i1+nsp][i2+n2+nsp]=a[i1][n2-1];
		}
	}
	//上边-下边
	for(i1=0;i1<nsp;i1++)
	{
		for(i2=0;i2<n2;i2++)
		{
			b[i1][i2+nsp]=a[0][i2];
			b[nsp+n1+i1][i2+nsp]=a[n1-1][i2];
		}
	}
	//左上角-右上角-左下角-右下角
	for(i1=0;i1<nsp;i1++)
	{
		for(i2=0;i2<nsp;i2++)
		{
			b[i1][i2]=a[0][0];
			b[i1][nsp+n2+i2]=a[0][n2-1];
			b[i1+nsp+n1][i2]=a[n1-1][0];
			b[i1+nsp+n1][i2+nsp+n2]=a[n1-1][n2-1];
		}
	}
	//
	for(i1=nsp;i1<n1+nsp;i1++)
	{
		for(i2=nsp;i2<n2+nsp;i2++)
		{
			a1=0;
			for(i11=i1-nsp;i11<=i1+nsp;i11++)
			{
				for(i22=i2-nsp;i22<=i2+nsp;i22++)
				{
					dist1=i11-i1;
					dist2=i22-i2;
					b1=exp(-(dist1*dist1+dist2*dist2)/(2.0*nsp/3.0*nsp/3.0));
					a1+=b1*b[i11][i22];
				}
			}
			a[i1-nsp][i2-nsp]=a1;
		}
	}
	a1=0;
	for(i11=0;i11<=2.0*nsp;i11++)
	{
		for(i22=0;i22<=2.0*nsp;i22++)
		{
			dist1=i11-nsp;
			dist2=i22-nsp;
			b1=exp(-(dist1*dist1+dist2*dist2)/(2.0*nsp/3.0*nsp/3.0));
			a1+=b1;
		}
	}
	for(i1=0;i1<n1;i1++)
	{
		for(i2=0;i2<n2;i2++)
		{
			a[i1][i2]/=a1;
		}
	}
	for(i1=0;i1<n1;i1++)
	{
		for(i2=0;i2<n2;i2++)
		{
			vp[i1*n2+i2] = a[i1][i2];
		}
	}
	free2dArray(a,n1,n2);
	free2dArray(b,n1e,n2e);
}
void velsmooth2d(float **vp,int n1,int n2,int nsp)
{
	// nsp 窗口大小
	// n1 行 nz
	// n2 列 nx
	// a 速度
	// b 扩展速度
	int n1e,n2e,i1,i2,i11,i22;
	double PI=3.141592653;
	float **a,**b;
	double a1,b1,dist1,dist2;
	n1e=n1+2*nsp; //数值维度扩展
	n2e=n2+2*nsp; //数组维度扩展
	// 开辟空间  b[n1e][n2e]
	a=(float**)calloc(n1,sizeof(float*));
    for(i1=0;i1<n1;i1++)
	{
		a[i1]=(float*)calloc(n2,sizeof(float));
	}
	b=(float**)calloc(n1e,sizeof(float*));
    for(i1=0;i1<n1e;i1++)
	{
		b[i1]=(float*)calloc(n2e,sizeof(float));
	}
	for(i1=0;i1<n1;i1++)
	{
		for(i2=0;i2<n2;i2++)
		{
			a[i1][i2]=vp[i1][i2];
		}
	}
	//中间
	for(i1=0;i1<n1;i1++)
	{
		for(i2=0;i2<n2;i2++)
		{
			b[i1+nsp][i2+nsp]=a[i1][i2];
		}
	}
	//左边-右边
	for(i1=0;i1<n1;i1++)
	{
		for(i2=0;i2<nsp;i2++)
		{
			b[i1+nsp][i2]=a[i1][0];
			b[i1+nsp][i2+n2+nsp]=a[i1][n2-1];
		}
	}
	//上边-下边
	for(i1=0;i1<nsp;i1++)
	{
		for(i2=0;i2<n2;i2++)
		{
			b[i1][i2+nsp]=a[0][i2];
			b[nsp+n1+i1][i2+nsp]=a[n1-1][i2];
		}
	}
	//左上角-右上角-左下角-右下角
	for(i1=0;i1<nsp;i1++)
	{
		for(i2=0;i2<nsp;i2++)
		{
			b[i1][i2]=a[0][0];
			b[i1][nsp+n2+i2]=a[0][n2-1];
			b[i1+nsp+n1][i2]=a[n1-1][0];
			b[i1+nsp+n1][i2+nsp+n2]=a[n1-1][n2-1];
		}
	}
	//
	for(i1=nsp;i1<n1+nsp;i1++)
	{
		for(i2=nsp;i2<n2+nsp;i2++)
		{
			a1=0;
			for(i11=i1-nsp;i11<=i1+nsp;i11++)
			{
				for(i22=i2-nsp;i22<=i2+nsp;i22++)
				{
					dist1=i11-i1;
					dist2=i22-i2;
					b1=exp(-(dist1*dist1+dist2*dist2)/(2.0*nsp/3.0*nsp/3.0));
					a1+=b1*b[i11][i22];
				}
			}
			a[i1-nsp][i2-nsp]=a1;
		}
	}
	a1=0;
	for(i11=0;i11<=2.0*nsp;i11++)
	{
		for(i22=0;i22<=2.0*nsp;i22++)
		{
			dist1=i11-nsp;
			dist2=i22-nsp;
			b1=exp(-(dist1*dist1+dist2*dist2)/(2.0*nsp/3.0*nsp/3.0));
			a1+=b1;
		}
	}
	for(i1=0;i1<n1;i1++)
	{
		for(i2=0;i2<n2;i2++)
		{
			a[i1][i2]/=a1;
		}
	}
	for(i1=0;i1<n1;i1++)
	{
		for(i2=0;i2<n2;i2++)
		{
			vp[i1][i2] = a[i1][i2];
		}
	}
	free2dArray(a,n1,n2);
	free2dArray(b,n1e,n2e);
}

static void convolve_cwp_s(int,int,float*,int,int,float*,int,int,float*);

void convolve_cwp (int lx, int ifx, float *x,
	   int ly, int ify, float *y,
	   int lz, int ifz, float *z)
/*****************************************************************************
Compute z = x convolved with y; i.e.,

           ifx+lx-1
    z[i] =   sum    x[j]*y[i-j]  ;  i = ifz,...,ifz+lz-1
            j=ifx
******************************************************************************
Input:
lx		length of x array
ifx		sample index of first x
x		array[lx] to be convolved with y
ly		length of y array
ify		sample index of first y
y		array[ly] with which x is to be convolved
lz		length of z array
ifz		sample index of first z

Output:
z		array[lz] containing x convolved with y
******************************************************************************
Notes:
The x samples are contained in x[0], x[1], ..., x[lx-1]; likewise for
the y and z samples.  The sample indices of the first x, y, and z values
determine the location of the origin for each array.  For example, if
z is to be a weighted average of the nearest 5 samples of y, one might
use 
	...
	x[0] = x[1] = x[2] = x[3] = x[4] = 1.0/5.0;
	conv(5,-2,x,lx,0,y,ly,0,z);
	...
In this example, the filter x is symmetric, with index of first sample = -2.

This function is optimized for architectures that can simultaneously perform
a multiply, add, and one load from memory; e.g., the IBM RISC System/6000.
Because, for each value of i, it accumulates the convolution sum z[i] in a
scalar, this function is not likely to be optimal for vector architectures.
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 11/23/91
*****************************************************************************/
#ifdef SIMPLE_CONV
{
	int ilx=ifx+lx-1,ily=ify+ly-1,ilz=ifz+lz-1,i,j,jlow,jhigh;
	float sum;
	
	x -= ifx;  y -= ify;  z -= ifz;
	for (i=ifz; i<=ilz; ++i) {
		jlow = i-ily;  if (jlow<ifx) jlow = ifx;
		jhigh = i-ify;  if (jhigh>ilx) jhigh = ilx;
		for (j=jlow,sum=0.0; j<=jhigh; ++j)
			sum += x[j]*y[i-j];
		z[i] = sum;
	}
}
#else
{
	int ilx=ifx+lx-1,ily=ify+ly-1,ilz=ifz+lz-1,
		i,j,ilow,ihigh,jlow,jhigh;
	float sa,sb,xa,xb,ya,yb,*t;

	/* if x is longer than y, swap x and y */
	if (lx>ly) {
		i = ifx;  ifx = ify;  ify = i;
		i = ilx;  ilx = ily;  ily = i;
		i = lx;  lx = ly;  ly = i;
		t = x;  x = y;  y = t;
	}
	
	/* handle short x with special code */
	if (lx>=1 && lx<=30) {
		convolve_cwp_s(lx,ifx,x,ly,ify,y,lz,ifz,z);
		return;
	}
	
	/* adjust pointers for indices of first samples */
	x -= ifx;
	y -= ify;
	z -= ifz;
		
	/* OFF LEFT:  i < ify+ifx */
	
	/* zero output for all i */
	ilow = ifz;
	ihigh = ify+ifx-1;  if (ihigh>ilz) ihigh = ilz;
	for (i=ilow; i<=ihigh; ++i)
		z[i] = 0.0;

	/* ROLLING ON:  ify+ifx <= i < ify+ilx */
	
	/* if necessary, do one i so that number of j in overlap is odd */
	if (i<ify+ilx && i<=ilz) {
		jlow = ifx;
		jhigh = i-ify;
		if ((jhigh-jlow)%2) {
			sa = 0.0;
			for (j=jlow; j<=jhigh; ++j)
				sa += x[j]*y[i-j];
			z[i++] = sa;
		}
	}
	
	/* loop over pairs of i and j */
	ilow = i;
	ihigh = ilx+ify-1;  if (ihigh>ilz) ihigh = ilz;
	jlow = ifx;
	jhigh = ilow-ify;
	for (i=ilow; i<ihigh; i+=2,jhigh+=2) {
		sa = sb = 0.0;
		xb = x[jhigh+1];
		yb = 0.0;
		for (j=jhigh; j>=jlow; j-=2) {
			sa += xb*yb;
			ya = y[i-j];
			sb += xb*ya;
			xa = x[j];
			sa += xa*ya;
			yb = y[i+1-j];
			sb += xa*yb;
			xb = x[j-1];
		}
		z[i] = sa;
		z[i+1] = sb;
	}
	
	/* if number of i is odd */
	if (i==ihigh) {
		jlow = ifx;
		jhigh = i-ify;
		sa = 0.0;
		for (j=jlow; j<=jhigh; ++j)
			sa += x[j]*y[i-j];
		z[i++] = sa;
	}
	
	/* MIDDLE:  ify+ilx <= i <= ily+ifx */
	
	/* determine limits for i and j */
	ilow = i;
	ihigh = ily+ifx;  if (ihigh>ilz) ihigh = ilz;
	jlow = ifx;
	jhigh = ilx;
	
	/* if number of j is even, do j in pairs with no leftover */
	if ((jhigh-jlow)%2) {
		for (i=ilow; i<ihigh; i+=2) {
			sa = sb = 0.0;
			yb = y[i+1-jlow];
			xa = x[jlow];
			for (j=jlow; j<jhigh; j+=2) {
				sb += xa*yb;
				ya = y[i-j];
				sa += xa*ya;
				xb = x[j+1];
				sb += xb*ya;
				yb = y[i-1-j];
				sa += xb*yb;
				xa = x[j+2];
			}
			z[i] = sa;
			z[i+1] = sb;
		}
	
	/* else, number of j is odd, so do j in pairs with leftover */
	} else {
		for (i=ilow; i<ihigh; i+=2) {
			sa = sb = 0.0;
			yb = y[i+1-jlow];
			xa = x[jlow];
			for (j=jlow; j<jhigh; j+=2) {
				sb += xa*yb;
				ya = y[i-j];
				sa += xa*ya;
				xb = x[j+1];
				sb += xb*ya;
				yb = y[i-1-j];
				sa += xb*yb;
				xa = x[j+2];
			}
			z[i] = sa+x[jhigh]*y[i-jhigh];
			z[i+1] = sb+x[jhigh]*y[i+1-jhigh];
		}
	}
	
	/* if number of i is odd */
	if (i==ihigh) {
		sa = 0.0;
		for (j=jlow; j<=jhigh; ++j)
			sa += x[j]*y[i-j];
		z[i++] = sa;
	}

	/* ROLLING OFF:  ily+ifx < i <= ily+ilx */
	
	/* if necessary, do one i so that number of j in overlap is even */
	if (i<=ily+ilx && i<=ilz) {
		jlow = i-ily;
		jhigh = ilx;
		if (!((jhigh-jlow)%2)) {
			sa = 0.0;
			for (j=jlow; j<=jhigh; ++j)
				sa += x[j]*y[i-j];
			z[i++] = sa;
		}
	}
	
	/* number of j is now even, so loop over both i and j in pairs */
	ilow = i;
	ihigh = ily+ilx;  if (ihigh>ilz) ihigh = ilz;
	jlow = ilow-ily;
	jhigh = ilx-2; /* Dave's new patch */
        for (i=ilow; i<ihigh; i+=2,jlow+=2) {
                sa = sb = 0.0;
                xa = x[jlow];
                yb = 0.0;
                for (j=jlow; j<jhigh; j+=2) {
                        sb += xa*yb;
                        ya = y[i-j];
                        sa += xa*ya;
                        xb = x[j+1];
                        sb += xb*ya;
                        yb = y[i-1-j];
                        sa += xb*yb;
                        xa = x[j+2];
                }
                sb += xa*yb;
                ya = y[i-j];
                sa += xa*ya;
                xb = x[j+1];
                sb += xb*ya;
                yb = y[i-1-j];
                sa += xb*yb;
                z[i] = sa;
                z[i+1] = sb;
        }
	
	/* if number of i is odd */
	if (i==ihigh) {
		jlow = i-ily;
		jhigh = ilx;
		sa = 0.0;
		for (j=jlow; j<=jhigh; ++j)
			sa += x[j]*y[i-j];
		z[i++] = sa;
	}
	
	/* OFF RIGHT:  ily+ilx < i */
	
	/* zero output for all i */
	ilow = i;
	ihigh = ilz;
	for (i=ilow; i<=ihigh; ++i)
		z[i] = 0.0;
}

/* internal function optimized for short x */
static void convolve_cwp_s (int lx, int ifx, float *x,
	   int ly, int ify, float *y, 
	   int lz, int ifz, float *z)
{
	int ilx=ifx+lx-1,ily=ify+ly-1,ilz=ifz+lz-1,
		i,j,ilow,ihigh,jlow,jhigh;
	float x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,
		x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,
		x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,
		ya,yb,z0,z1,sum;
	
	x -= ifx;
	y -= ify;
	z -= ifz;
	
		
	/* OFF LEFT:  i < ifx+ify */
	ilow = ifz;
	ihigh = ify+ifx-1;  if (ihigh>ilz) ihigh = ilz;
	for (i=ilow; i<=ihigh; ++i)
		z[i] = 0.0;
	
	/* ROLLING ON:  ify+ifx <= i < ify+ilx */
	ilow = ify+ifx;  if (ilow<ifz) ilow = ifz;
	ihigh = ify+ilx-1;  if (ihigh>ilz) ihigh = ilz;
	jlow = ifx;
	jhigh = ilow-ify;
	for (i=ilow; i<=ihigh; ++i,++jhigh) {
		for (j=jlow,sum=0.0; j<=jhigh; ++j)
			sum += x[j]*y[i-j];
		z[i] = sum;
	}
	
	/* MIDDLE:  ify+ilx <= i <= ily+ifx */
	ilow = ify+ilx;  if (ilow<ifz) ilow = ifz;
	ihigh = ily+ifx;  if (ihigh>ilz) ihigh = ilz;
	if (lx==1) {
		x0 = x[ifx];
		for (i=ilow; i<=ihigh-1; i+=2) {
			ya = y[i+1-ifx];  z1 = x0*ya;
			yb = y[i-ifx];  z0 = x0*yb;
			z[i+1] = z1;
			z[i] = z0;
		}
	} else if (lx==2) {
		x0 = x[ifx];
		x1 = x[ifx+1];
		for (i=ilow; i<=ihigh-1; i+=2) {
			ya = y[i+1-ifx];  z1 = x0*ya;
			yb = y[i-ifx];  z0 = x0*yb;  z1 += x1*yb;
			ya = y[i-ifx-1];  z0 += x1*ya;
			z[i+1] = z1;
			z[i] = z0;
		}
	} else if (lx==3) {
		x0 = x[ifx];
		x1 = x[ifx+1];
		x2 = x[ifx+2];
		for (i=ilow; i<=ihigh-1; i+=2) {
			ya = y[i+1-ifx];  z1 = x0*ya;
			yb = y[i-ifx];  z0 = x0*yb;  z1 += x1*yb;
			ya = y[i-ifx-1];  z0 += x1*ya;  z1 += x2*ya;
			yb = y[i-ifx-2];  z0 += x2*yb;
			z[i+1] = z1;
			z[i] = z0;
		}
	} else if (lx==4) {
		x0 = x[ifx];
		x1 = x[ifx+1];
		x2 = x[ifx+2];
		x3 = x[ifx+3];
		for (i=ilow; i<=ihigh-1; i+=2) {
			ya = y[i+1-ifx];  z1 = x0*ya;
			yb = y[i-ifx];  z0 = x0*yb;  z1 += x1*yb;
			ya = y[i-ifx-1];  z0 += x1*ya;  z1 += x2*ya;
			yb = y[i-ifx-2];  z0 += x2*yb;  z1 += x3*yb;
			ya = y[i-ifx-3];  z0 += x3*ya;
			z[i+1] = z1;
			z[i] = z0;
		}
	} else if (lx==5) {
		x0 = x[ifx];
		x1 = x[ifx+1];
		x2 = x[ifx+2];
		x3 = x[ifx+3];
		x4 = x[ifx+4];
		for (i=ilow; i<=ihigh-1; i+=2) {
			ya = y[i+1-ifx];  z1 = x0*ya;
			yb = y[i-ifx];  z0 = x0*yb;  z1 += x1*yb;
			ya = y[i-ifx-1];  z0 += x1*ya;  z1 += x2*ya;
			yb = y[i-ifx-2];  z0 += x2*yb;  z1 += x3*yb;
			ya = y[i-ifx-3];  z0 += x3*ya;  z1 += x4*ya;
			yb = y[i-ifx-4];  z0 += x4*yb;
			z[i+1] = z1;
			z[i] = z0;
		}
	} else if (lx==6) {
		x0 = x[ifx];
		x1 = x[ifx+1];
		x2 = x[ifx+2];
		x3 = x[ifx+3];
		x4 = x[ifx+4];
		x5 = x[ifx+5];
		for (i=ilow; i<=ihigh-1; i+=2) {
			ya = y[i+1-ifx];  z1 = x0*ya;
			yb = y[i-ifx];  z0 = x0*yb;  z1 += x1*yb;
			ya = y[i-ifx-1];  z0 += x1*ya;  z1 += x2*ya;
			yb = y[i-ifx-2];  z0 += x2*yb;  z1 += x3*yb;
			ya = y[i-ifx-3];  z0 += x3*ya;  z1 += x4*ya;
			yb = y[i-ifx-4];  z0 += x4*yb;  z1 += x5*yb;
			ya = y[i-ifx-5];  z0 += x5*ya;
			z[i+1] = z1;
			z[i] = z0;
		}
	} else if (lx==7) {
		x0 = x[ifx];
		x1 = x[ifx+1];
		x2 = x[ifx+2];
		x3 = x[ifx+3];
		x4 = x[ifx+4];
		x5 = x[ifx+5];
		x6 = x[ifx+6];
		for (i=ilow; i<=ihigh-1; i+=2) {
			ya = y[i+1-ifx];  z1 = x0*ya;
			yb = y[i-ifx];  z0 = x0*yb;  z1 += x1*yb;
			ya = y[i-ifx-1];  z0 += x1*ya;  z1 += x2*ya;
			yb = y[i-ifx-2];  z0 += x2*yb;  z1 += x3*yb;
			ya = y[i-ifx-3];  z0 += x3*ya;  z1 += x4*ya;
			yb = y[i-ifx-4];  z0 += x4*yb;  z1 += x5*yb;
			ya = y[i-ifx-5];  z0 += x5*ya;  z1 += x6*ya;
			yb = y[i-ifx-6];  z0 += x6*yb;
			z[i+1] = z1;
			z[i] = z0;
		}
	} else if (lx==8) {
		x0 = x[ifx];
		x1 = x[ifx+1];
		x2 = x[ifx+2];
		x3 = x[ifx+3];
		x4 = x[ifx+4];
		x5 = x[ifx+5];
		x6 = x[ifx+6];
		x7 = x[ifx+7];
		for (i=ilow; i<=ihigh-1; i+=2) {
			ya = y[i+1-ifx];  z1 = x0*ya;
			yb = y[i-ifx];  z0 = x0*yb;  z1 += x1*yb;
			ya = y[i-ifx-1];  z0 += x1*ya;  z1 += x2*ya;
			yb = y[i-ifx-2];  z0 += x2*yb;  z1 += x3*yb;
			ya = y[i-ifx-3];  z0 += x3*ya;  z1 += x4*ya;
			yb = y[i-ifx-4];  z0 += x4*yb;  z1 += x5*yb;
			ya = y[i-ifx-5];  z0 += x5*ya;  z1 += x6*ya;
			yb = y[i-ifx-6];  z0 += x6*yb;  z1 += x7*yb;
			ya = y[i-ifx-7];  z0 += x7*ya;
			z[i+1] = z1;
			z[i] = z0;
		}
	} else if (lx==9) {
		x0 = x[ifx];
		x1 = x[ifx+1];
		x2 = x[ifx+2];
		x3 = x[ifx+3];
		x4 = x[ifx+4];
		x5 = x[ifx+5];
		x6 = x[ifx+6];
		x7 = x[ifx+7];
		x8 = x[ifx+8];
		for (i=ilow; i<=ihigh-1; i+=2) {
			ya = y[i+1-ifx];  z1 = x0*ya;
			yb = y[i-ifx];  z0 = x0*yb;  z1 += x1*yb;
			ya = y[i-ifx-1];  z0 += x1*ya;  z1 += x2*ya;
			yb = y[i-ifx-2];  z0 += x2*yb;  z1 += x3*yb;
			ya = y[i-ifx-3];  z0 += x3*ya;  z1 += x4*ya;
			yb = y[i-ifx-4];  z0 += x4*yb;  z1 += x5*yb;
			ya = y[i-ifx-5];  z0 += x5*ya;  z1 += x6*ya;
			yb = y[i-ifx-6];  z0 += x6*yb;  z1 += x7*yb;
			ya = y[i-ifx-7];  z0 += x7*ya;  z1 += x8*ya;
			yb = y[i-ifx-8];  z0 += x8*yb;
			z[i+1] = z1;
			z[i] = z0;
		}
	} else if (lx==10) {
		x0 = x[ifx];
		x1 = x[ifx+1];
		x2 = x[ifx+2];
		x3 = x[ifx+3];
		x4 = x[ifx+4];
		x5 = x[ifx+5];
		x6 = x[ifx+6];
		x7 = x[ifx+7];
		x8 = x[ifx+8];
		x9 = x[ifx+9];
		for (i=ilow; i<=ihigh-1; i+=2) {
			ya = y[i+1-ifx];  z1 = x0*ya;
			yb = y[i-ifx];  z0 = x0*yb;  z1 += x1*yb;
			ya = y[i-ifx-1];  z0 += x1*ya;  z1 += x2*ya;
			yb = y[i-ifx-2];  z0 += x2*yb;  z1 += x3*yb;
			ya = y[i-ifx-3];  z0 += x3*ya;  z1 += x4*ya;
			yb = y[i-ifx-4];  z0 += x4*yb;  z1 += x5*yb;
			ya = y[i-ifx-5];  z0 += x5*ya;  z1 += x6*ya;
			yb = y[i-ifx-6];  z0 += x6*yb;  z1 += x7*yb;
			ya = y[i-ifx-7];  z0 += x7*ya;  z1 += x8*ya;
			yb = y[i-ifx-8];  z0 += x8*yb;  z1 += x9*yb;
			ya = y[i-ifx-9];  z0 += x9*ya;
			z[i+1] = z1;
			z[i] = z0;
		}
	} else if (lx==11) {
		x0 = x[ifx];
		x1 = x[ifx+1];
		x2 = x[ifx+2];
		x3 = x[ifx+3];
		x4 = x[ifx+4];
		x5 = x[ifx+5];
		x6 = x[ifx+6];
		x7 = x[ifx+7];
		x8 = x[ifx+8];
		x9 = x[ifx+9];
		x10 = x[ifx+10];
		for (i=ilow; i<=ihigh-1; i+=2) {
			ya = y[i+1-ifx];  z1 = x0*ya;
			yb = y[i-ifx];  z0 = x0*yb;  z1 += x1*yb;
			ya = y[i-ifx-1];  z0 += x1*ya;  z1 += x2*ya;
			yb = y[i-ifx-2];  z0 += x2*yb;  z1 += x3*yb;
			ya = y[i-ifx-3];  z0 += x3*ya;  z1 += x4*ya;
			yb = y[i-ifx-4];  z0 += x4*yb;  z1 += x5*yb;
			ya = y[i-ifx-5];  z0 += x5*ya;  z1 += x6*ya;
			yb = y[i-ifx-6];  z0 += x6*yb;  z1 += x7*yb;
			ya = y[i-ifx-7];  z0 += x7*ya;  z1 += x8*ya;
			yb = y[i-ifx-8];  z0 += x8*yb;  z1 += x9*yb;
			ya = y[i-ifx-9];  z0 += x9*ya;  z1 += x10*ya;
			yb = y[i-ifx-10];  z0 += x10*yb;
			z[i+1] = z1;
			z[i] = z0;
		}
	} else if (lx==12) {
		x0 = x[ifx];
		x1 = x[ifx+1];
		x2 = x[ifx+2];
		x3 = x[ifx+3];
		x4 = x[ifx+4];
		x5 = x[ifx+5];
		x6 = x[ifx+6];
		x7 = x[ifx+7];
		x8 = x[ifx+8];
		x9 = x[ifx+9];
		x10 = x[ifx+10];
		x11 = x[ifx+11];
		for (i=ilow; i<=ihigh-1; i+=2) {
			ya = y[i+1-ifx];  z1 = x0*ya;
			yb = y[i-ifx];  z0 = x0*yb;  z1 += x1*yb;
			ya = y[i-ifx-1];  z0 += x1*ya;  z1 += x2*ya;
			yb = y[i-ifx-2];  z0 += x2*yb;  z1 += x3*yb;
			ya = y[i-ifx-3];  z0 += x3*ya;  z1 += x4*ya;
			yb = y[i-ifx-4];  z0 += x4*yb;  z1 += x5*yb;
			ya = y[i-ifx-5];  z0 += x5*ya;  z1 += x6*ya;
			yb = y[i-ifx-6];  z0 += x6*yb;  z1 += x7*yb;
			ya = y[i-ifx-7];  z0 += x7*ya;  z1 += x8*ya;
			yb = y[i-ifx-8];  z0 += x8*yb;  z1 += x9*yb;
			ya = y[i-ifx-9];  z0 += x9*ya;  z1 += x10*ya;
			yb = y[i-ifx-10];  z0 += x10*yb;  z1 += x11*yb;
			ya = y[i-ifx-11];  z0 += x11*ya;
			z[i+1] = z1;
			z[i] = z0;
		}
	} else if (lx==13) {
		x0 = x[ifx];
		x1 = x[ifx+1];
		x2 = x[ifx+2];
		x3 = x[ifx+3];
		x4 = x[ifx+4];
		x5 = x[ifx+5];
		x6 = x[ifx+6];
		x7 = x[ifx+7];
		x8 = x[ifx+8];
		x9 = x[ifx+9];
		x10 = x[ifx+10];
		x11 = x[ifx+11];
		x12 = x[ifx+12];
		for (i=ilow; i<=ihigh-1; i+=2) {
			ya = y[i+1-ifx];  z1 = x0*ya;
			yb = y[i-ifx];  z0 = x0*yb;  z1 += x1*yb;
			ya = y[i-ifx-1];  z0 += x1*ya;  z1 += x2*ya;
			yb = y[i-ifx-2];  z0 += x2*yb;  z1 += x3*yb;
			ya = y[i-ifx-3];  z0 += x3*ya;  z1 += x4*ya;
			yb = y[i-ifx-4];  z0 += x4*yb;  z1 += x5*yb;
			ya = y[i-ifx-5];  z0 += x5*ya;  z1 += x6*ya;
			yb = y[i-ifx-6];  z0 += x6*yb;  z1 += x7*yb;
			ya = y[i-ifx-7];  z0 += x7*ya;  z1 += x8*ya;
			yb = y[i-ifx-8];  z0 += x8*yb;  z1 += x9*yb;
			ya = y[i-ifx-9];  z0 += x9*ya;  z1 += x10*ya;
			yb = y[i-ifx-10];  z0 += x10*yb;  z1 += x11*yb;
			ya = y[i-ifx-11];  z0 += x11*ya;  z1 += x12*ya;
			yb = y[i-ifx-12];  z0 += x12*yb;
			z[i+1] = z1;
			z[i] = z0;
		}
	} else if (lx==14) {
		x0 = x[ifx];
		x1 = x[ifx+1];
		x2 = x[ifx+2];
		x3 = x[ifx+3];
		x4 = x[ifx+4];
		x5 = x[ifx+5];
		x6 = x[ifx+6];
		x7 = x[ifx+7];
		x8 = x[ifx+8];
		x9 = x[ifx+9];
		x10 = x[ifx+10];
		x11 = x[ifx+11];
		x12 = x[ifx+12];
		x13 = x[ifx+13];
		for (i=ilow; i<=ihigh-1; i+=2) {
			ya = y[i+1-ifx];  z1 = x0*ya;
			yb = y[i-ifx];  z0 = x0*yb;  z1 += x1*yb;
			ya = y[i-ifx-1];  z0 += x1*ya;  z1 += x2*ya;
			yb = y[i-ifx-2];  z0 += x2*yb;  z1 += x3*yb;
			ya = y[i-ifx-3];  z0 += x3*ya;  z1 += x4*ya;
			yb = y[i-ifx-4];  z0 += x4*yb;  z1 += x5*yb;
			ya = y[i-ifx-5];  z0 += x5*ya;  z1 += x6*ya;
			yb = y[i-ifx-6];  z0 += x6*yb;  z1 += x7*yb;
			ya = y[i-ifx-7];  z0 += x7*ya;  z1 += x8*ya;
			yb = y[i-ifx-8];  z0 += x8*yb;  z1 += x9*yb;
			ya = y[i-ifx-9];  z0 += x9*ya;  z1 += x10*ya;
			yb = y[i-ifx-10];  z0 += x10*yb;  z1 += x11*yb;
			ya = y[i-ifx-11];  z0 += x11*ya;  z1 += x12*ya;
			yb = y[i-ifx-12];  z0 += x12*yb;  z1 += x13*yb;
			ya = y[i-ifx-13];  z0 += x13*ya;
			z[i+1] = z1;
			z[i] = z0;
		}
	} else if (lx==15) {
		x0 = x[ifx];
		x1 = x[ifx+1];
		x2 = x[ifx+2];
		x3 = x[ifx+3];
		x4 = x[ifx+4];
		x5 = x[ifx+5];
		x6 = x[ifx+6];
		x7 = x[ifx+7];
		x8 = x[ifx+8];
		x9 = x[ifx+9];
		x10 = x[ifx+10];
		x11 = x[ifx+11];
		x12 = x[ifx+12];
		x13 = x[ifx+13];
		x14 = x[ifx+14];
		for (i=ilow; i<=ihigh-1; i+=2) {
			ya = y[i+1-ifx];  z1 = x0*ya;
			yb = y[i-ifx];  z0 = x0*yb;  z1 += x1*yb;
			ya = y[i-ifx-1];  z0 += x1*ya;  z1 += x2*ya;
			yb = y[i-ifx-2];  z0 += x2*yb;  z1 += x3*yb;
			ya = y[i-ifx-3];  z0 += x3*ya;  z1 += x4*ya;
			yb = y[i-ifx-4];  z0 += x4*yb;  z1 += x5*yb;
			ya = y[i-ifx-5];  z0 += x5*ya;  z1 += x6*ya;
			yb = y[i-ifx-6];  z0 += x6*yb;  z1 += x7*yb;
			ya = y[i-ifx-7];  z0 += x7*ya;  z1 += x8*ya;
			yb = y[i-ifx-8];  z0 += x8*yb;  z1 += x9*yb;
			ya = y[i-ifx-9];  z0 += x9*ya;  z1 += x10*ya;
			yb = y[i-ifx-10];  z0 += x10*yb;  z1 += x11*yb;
			ya = y[i-ifx-11];  z0 += x11*ya;  z1 += x12*ya;
			yb = y[i-ifx-12];  z0 += x12*yb;  z1 += x13*yb;
			ya = y[i-ifx-13];  z0 += x13*ya;  z1 += x14*ya;
			yb = y[i-ifx-14];  z0 += x14*yb;
			z[i+1] = z1;
			z[i] = z0;
		}
	} else if (lx==16) {
		x0 = x[ifx];
		x1 = x[ifx+1];
		x2 = x[ifx+2];
		x3 = x[ifx+3];
		x4 = x[ifx+4];
		x5 = x[ifx+5];
		x6 = x[ifx+6];
		x7 = x[ifx+7];
		x8 = x[ifx+8];
		x9 = x[ifx+9];
		x10 = x[ifx+10];
		x11 = x[ifx+11];
		x12 = x[ifx+12];
		x13 = x[ifx+13];
		x14 = x[ifx+14];
		x15 = x[ifx+15];
		for (i=ilow; i<=ihigh-1; i+=2) {
			ya = y[i+1-ifx];  z1 = x0*ya;
			yb = y[i-ifx];  z0 = x0*yb;  z1 += x1*yb;
			ya = y[i-ifx-1];  z0 += x1*ya;  z1 += x2*ya;
			yb = y[i-ifx-2];  z0 += x2*yb;  z1 += x3*yb;
			ya = y[i-ifx-3];  z0 += x3*ya;  z1 += x4*ya;
			yb = y[i-ifx-4];  z0 += x4*yb;  z1 += x5*yb;
			ya = y[i-ifx-5];  z0 += x5*ya;  z1 += x6*ya;
			yb = y[i-ifx-6];  z0 += x6*yb;  z1 += x7*yb;
			ya = y[i-ifx-7];  z0 += x7*ya;  z1 += x8*ya;
			yb = y[i-ifx-8];  z0 += x8*yb;  z1 += x9*yb;
			ya = y[i-ifx-9];  z0 += x9*ya;  z1 += x10*ya;
			yb = y[i-ifx-10];  z0 += x10*yb;  z1 += x11*yb;
			ya = y[i-ifx-11];  z0 += x11*ya;  z1 += x12*ya;
			yb = y[i-ifx-12];  z0 += x12*yb;  z1 += x13*yb;
			ya = y[i-ifx-13];  z0 += x13*ya;  z1 += x14*ya;
			yb = y[i-ifx-14];  z0 += x14*yb;  z1 += x15*yb;
			ya = y[i-ifx-15];  z0 += x15*ya;
			z[i+1] = z1;
			z[i] = z0;
		}
	} else if (lx==17) {
		x0 = x[ifx];
		x1 = x[ifx+1];
		x2 = x[ifx+2];
		x3 = x[ifx+3];
		x4 = x[ifx+4];
		x5 = x[ifx+5];
		x6 = x[ifx+6];
		x7 = x[ifx+7];
		x8 = x[ifx+8];
		x9 = x[ifx+9];
		x10 = x[ifx+10];
		x11 = x[ifx+11];
		x12 = x[ifx+12];
		x13 = x[ifx+13];
		x14 = x[ifx+14];
		x15 = x[ifx+15];
		x16 = x[ifx+16];
		for (i=ilow; i<=ihigh-1; i+=2) {
			ya = y[i+1-ifx];  z1 = x0*ya;
			yb = y[i-ifx];  z0 = x0*yb;  z1 += x1*yb;
			ya = y[i-ifx-1];  z0 += x1*ya;  z1 += x2*ya;
			yb = y[i-ifx-2];  z0 += x2*yb;  z1 += x3*yb;
			ya = y[i-ifx-3];  z0 += x3*ya;  z1 += x4*ya;
			yb = y[i-ifx-4];  z0 += x4*yb;  z1 += x5*yb;
			ya = y[i-ifx-5];  z0 += x5*ya;  z1 += x6*ya;
			yb = y[i-ifx-6];  z0 += x6*yb;  z1 += x7*yb;
			ya = y[i-ifx-7];  z0 += x7*ya;  z1 += x8*ya;
			yb = y[i-ifx-8];  z0 += x8*yb;  z1 += x9*yb;
			ya = y[i-ifx-9];  z0 += x9*ya;  z1 += x10*ya;
			yb = y[i-ifx-10];  z0 += x10*yb;  z1 += x11*yb;
			ya = y[i-ifx-11];  z0 += x11*ya;  z1 += x12*ya;
			yb = y[i-ifx-12];  z0 += x12*yb;  z1 += x13*yb;
			ya = y[i-ifx-13];  z0 += x13*ya;  z1 += x14*ya;
			yb = y[i-ifx-14];  z0 += x14*yb;  z1 += x15*yb;
			ya = y[i-ifx-15];  z0 += x15*ya;  z1 += x16*ya;
			yb = y[i-ifx-16];  z0 += x16*yb;
			z[i+1] = z1;
			z[i] = z0;
		}
	} else if (lx==18) {
		x0 = x[ifx];
		x1 = x[ifx+1];
		x2 = x[ifx+2];
		x3 = x[ifx+3];
		x4 = x[ifx+4];
		x5 = x[ifx+5];
		x6 = x[ifx+6];
		x7 = x[ifx+7];
		x8 = x[ifx+8];
		x9 = x[ifx+9];
		x10 = x[ifx+10];
		x11 = x[ifx+11];
		x12 = x[ifx+12];
		x13 = x[ifx+13];
		x14 = x[ifx+14];
		x15 = x[ifx+15];
		x16 = x[ifx+16];
		x17 = x[ifx+17];
		for (i=ilow; i<=ihigh-1; i+=2) {
			ya = y[i+1-ifx];  z1 = x0*ya;
			yb = y[i-ifx];  z0 = x0*yb;  z1 += x1*yb;
			ya = y[i-ifx-1];  z0 += x1*ya;  z1 += x2*ya;
			yb = y[i-ifx-2];  z0 += x2*yb;  z1 += x3*yb;
			ya = y[i-ifx-3];  z0 += x3*ya;  z1 += x4*ya;
			yb = y[i-ifx-4];  z0 += x4*yb;  z1 += x5*yb;
			ya = y[i-ifx-5];  z0 += x5*ya;  z1 += x6*ya;
			yb = y[i-ifx-6];  z0 += x6*yb;  z1 += x7*yb;
			ya = y[i-ifx-7];  z0 += x7*ya;  z1 += x8*ya;
			yb = y[i-ifx-8];  z0 += x8*yb;  z1 += x9*yb;
			ya = y[i-ifx-9];  z0 += x9*ya;  z1 += x10*ya;
			yb = y[i-ifx-10];  z0 += x10*yb;  z1 += x11*yb;
			ya = y[i-ifx-11];  z0 += x11*ya;  z1 += x12*ya;
			yb = y[i-ifx-12];  z0 += x12*yb;  z1 += x13*yb;
			ya = y[i-ifx-13];  z0 += x13*ya;  z1 += x14*ya;
			yb = y[i-ifx-14];  z0 += x14*yb;  z1 += x15*yb;
			ya = y[i-ifx-15];  z0 += x15*ya;  z1 += x16*ya;
			yb = y[i-ifx-16];  z0 += x16*yb;  z1 += x17*yb;
			ya = y[i-ifx-17];  z0 += x17*ya;
			z[i+1] = z1;
			z[i] = z0;
		}
	} else if (lx==19) {
		x0 = x[ifx];
		x1 = x[ifx+1];
		x2 = x[ifx+2];
		x3 = x[ifx+3];
		x4 = x[ifx+4];
		x5 = x[ifx+5];
		x6 = x[ifx+6];
		x7 = x[ifx+7];
		x8 = x[ifx+8];
		x9 = x[ifx+9];
		x10 = x[ifx+10];
		x11 = x[ifx+11];
		x12 = x[ifx+12];
		x13 = x[ifx+13];
		x14 = x[ifx+14];
		x15 = x[ifx+15];
		x16 = x[ifx+16];
		x17 = x[ifx+17];
		x18 = x[ifx+18];
		for (i=ilow; i<=ihigh-1; i+=2) {
			ya = y[i+1-ifx];  z1 = x0*ya;
			yb = y[i-ifx];  z0 = x0*yb;  z1 += x1*yb;
			ya = y[i-ifx-1];  z0 += x1*ya;  z1 += x2*ya;
			yb = y[i-ifx-2];  z0 += x2*yb;  z1 += x3*yb;
			ya = y[i-ifx-3];  z0 += x3*ya;  z1 += x4*ya;
			yb = y[i-ifx-4];  z0 += x4*yb;  z1 += x5*yb;
			ya = y[i-ifx-5];  z0 += x5*ya;  z1 += x6*ya;
			yb = y[i-ifx-6];  z0 += x6*yb;  z1 += x7*yb;
			ya = y[i-ifx-7];  z0 += x7*ya;  z1 += x8*ya;
			yb = y[i-ifx-8];  z0 += x8*yb;  z1 += x9*yb;
			ya = y[i-ifx-9];  z0 += x9*ya;  z1 += x10*ya;
			yb = y[i-ifx-10];  z0 += x10*yb;  z1 += x11*yb;
			ya = y[i-ifx-11];  z0 += x11*ya;  z1 += x12*ya;
			yb = y[i-ifx-12];  z0 += x12*yb;  z1 += x13*yb;
			ya = y[i-ifx-13];  z0 += x13*ya;  z1 += x14*ya;
			yb = y[i-ifx-14];  z0 += x14*yb;  z1 += x15*yb;
			ya = y[i-ifx-15];  z0 += x15*ya;  z1 += x16*ya;
			yb = y[i-ifx-16];  z0 += x16*yb;  z1 += x17*yb;
			ya = y[i-ifx-17];  z0 += x17*ya;  z1 += x18*ya;
			yb = y[i-ifx-18];  z0 += x18*yb;
			z[i+1] = z1;
			z[i] = z0;
		}
	} else if (lx==20) {
		x0 = x[ifx];
		x1 = x[ifx+1];
		x2 = x[ifx+2];
		x3 = x[ifx+3];
		x4 = x[ifx+4];
		x5 = x[ifx+5];
		x6 = x[ifx+6];
		x7 = x[ifx+7];
		x8 = x[ifx+8];
		x9 = x[ifx+9];
		x10 = x[ifx+10];
		x11 = x[ifx+11];
		x12 = x[ifx+12];
		x13 = x[ifx+13];
		x14 = x[ifx+14];
		x15 = x[ifx+15];
		x16 = x[ifx+16];
		x17 = x[ifx+17];
		x18 = x[ifx+18];
		x19 = x[ifx+19];
		for (i=ilow; i<=ihigh-1; i+=2) {
			ya = y[i+1-ifx];  z1 = x0*ya;
			yb = y[i-ifx];  z0 = x0*yb;  z1 += x1*yb;
			ya = y[i-ifx-1];  z0 += x1*ya;  z1 += x2*ya;
			yb = y[i-ifx-2];  z0 += x2*yb;  z1 += x3*yb;
			ya = y[i-ifx-3];  z0 += x3*ya;  z1 += x4*ya;
			yb = y[i-ifx-4];  z0 += x4*yb;  z1 += x5*yb;
			ya = y[i-ifx-5];  z0 += x5*ya;  z1 += x6*ya;
			yb = y[i-ifx-6];  z0 += x6*yb;  z1 += x7*yb;
			ya = y[i-ifx-7];  z0 += x7*ya;  z1 += x8*ya;
			yb = y[i-ifx-8];  z0 += x8*yb;  z1 += x9*yb;
			ya = y[i-ifx-9];  z0 += x9*ya;  z1 += x10*ya;
			yb = y[i-ifx-10];  z0 += x10*yb;  z1 += x11*yb;
			ya = y[i-ifx-11];  z0 += x11*ya;  z1 += x12*ya;
			yb = y[i-ifx-12];  z0 += x12*yb;  z1 += x13*yb;
			ya = y[i-ifx-13];  z0 += x13*ya;  z1 += x14*ya;
			yb = y[i-ifx-14];  z0 += x14*yb;  z1 += x15*yb;
			ya = y[i-ifx-15];  z0 += x15*ya;  z1 += x16*ya;
			yb = y[i-ifx-16];  z0 += x16*yb;  z1 += x17*yb;
			ya = y[i-ifx-17];  z0 += x17*ya;  z1 += x18*ya;
			yb = y[i-ifx-18];  z0 += x18*yb;  z1 += x19*yb;
			ya = y[i-ifx-19];  z0 += x19*ya;
			z[i+1] = z1;
			z[i] = z0;
		}
	} else if (lx==21) {
		x0 = x[ifx];
		x1 = x[ifx+1];
		x2 = x[ifx+2];
		x3 = x[ifx+3];
		x4 = x[ifx+4];
		x5 = x[ifx+5];
		x6 = x[ifx+6];
		x7 = x[ifx+7];
		x8 = x[ifx+8];
		x9 = x[ifx+9];
		x10 = x[ifx+10];
		x11 = x[ifx+11];
		x12 = x[ifx+12];
		x13 = x[ifx+13];
		x14 = x[ifx+14];
		x15 = x[ifx+15];
		x16 = x[ifx+16];
		x17 = x[ifx+17];
		x18 = x[ifx+18];
		x19 = x[ifx+19];
		x20 = x[ifx+20];
		for (i=ilow; i<=ihigh-1; i+=2) {
			ya = y[i+1-ifx];  z1 = x0*ya;
			yb = y[i-ifx];  z0 = x0*yb;  z1 += x1*yb;
			ya = y[i-ifx-1];  z0 += x1*ya;  z1 += x2*ya;
			yb = y[i-ifx-2];  z0 += x2*yb;  z1 += x3*yb;
			ya = y[i-ifx-3];  z0 += x3*ya;  z1 += x4*ya;
			yb = y[i-ifx-4];  z0 += x4*yb;  z1 += x5*yb;
			ya = y[i-ifx-5];  z0 += x5*ya;  z1 += x6*ya;
			yb = y[i-ifx-6];  z0 += x6*yb;  z1 += x7*yb;
			ya = y[i-ifx-7];  z0 += x7*ya;  z1 += x8*ya;
			yb = y[i-ifx-8];  z0 += x8*yb;  z1 += x9*yb;
			ya = y[i-ifx-9];  z0 += x9*ya;  z1 += x10*ya;
			yb = y[i-ifx-10];  z0 += x10*yb;  z1 += x11*yb;
			ya = y[i-ifx-11];  z0 += x11*ya;  z1 += x12*ya;
			yb = y[i-ifx-12];  z0 += x12*yb;  z1 += x13*yb;
			ya = y[i-ifx-13];  z0 += x13*ya;  z1 += x14*ya;
			yb = y[i-ifx-14];  z0 += x14*yb;  z1 += x15*yb;
			ya = y[i-ifx-15];  z0 += x15*ya;  z1 += x16*ya;
			yb = y[i-ifx-16];  z0 += x16*yb;  z1 += x17*yb;
			ya = y[i-ifx-17];  z0 += x17*ya;  z1 += x18*ya;
			yb = y[i-ifx-18];  z0 += x18*yb;  z1 += x19*yb;
			ya = y[i-ifx-19];  z0 += x19*ya;  z1 += x20*ya;
			yb = y[i-ifx-20];  z0 += x20*yb;
			z[i+1] = z1;
			z[i] = z0;
		}
	} else if (lx==22) {
		x0 = x[ifx];
		x1 = x[ifx+1];
		x2 = x[ifx+2];
		x3 = x[ifx+3];
		x4 = x[ifx+4];
		x5 = x[ifx+5];
		x6 = x[ifx+6];
		x7 = x[ifx+7];
		x8 = x[ifx+8];
		x9 = x[ifx+9];
		x10 = x[ifx+10];
		x11 = x[ifx+11];
		x12 = x[ifx+12];
		x13 = x[ifx+13];
		x14 = x[ifx+14];
		x15 = x[ifx+15];
		x16 = x[ifx+16];
		x17 = x[ifx+17];
		x18 = x[ifx+18];
		x19 = x[ifx+19];
		x20 = x[ifx+20];
		x21 = x[ifx+21];
		for (i=ilow; i<=ihigh-1; i+=2) {
			ya = y[i+1-ifx];  z1 = x0*ya;
			yb = y[i-ifx];  z0 = x0*yb;  z1 += x1*yb;
			ya = y[i-ifx-1];  z0 += x1*ya;  z1 += x2*ya;
			yb = y[i-ifx-2];  z0 += x2*yb;  z1 += x3*yb;
			ya = y[i-ifx-3];  z0 += x3*ya;  z1 += x4*ya;
			yb = y[i-ifx-4];  z0 += x4*yb;  z1 += x5*yb;
			ya = y[i-ifx-5];  z0 += x5*ya;  z1 += x6*ya;
			yb = y[i-ifx-6];  z0 += x6*yb;  z1 += x7*yb;
			ya = y[i-ifx-7];  z0 += x7*ya;  z1 += x8*ya;
			yb = y[i-ifx-8];  z0 += x8*yb;  z1 += x9*yb;
			ya = y[i-ifx-9];  z0 += x9*ya;  z1 += x10*ya;
			yb = y[i-ifx-10];  z0 += x10*yb;  z1 += x11*yb;
			ya = y[i-ifx-11];  z0 += x11*ya;  z1 += x12*ya;
			yb = y[i-ifx-12];  z0 += x12*yb;  z1 += x13*yb;
			ya = y[i-ifx-13];  z0 += x13*ya;  z1 += x14*ya;
			yb = y[i-ifx-14];  z0 += x14*yb;  z1 += x15*yb;
			ya = y[i-ifx-15];  z0 += x15*ya;  z1 += x16*ya;
			yb = y[i-ifx-16];  z0 += x16*yb;  z1 += x17*yb;
			ya = y[i-ifx-17];  z0 += x17*ya;  z1 += x18*ya;
			yb = y[i-ifx-18];  z0 += x18*yb;  z1 += x19*yb;
			ya = y[i-ifx-19];  z0 += x19*ya;  z1 += x20*ya;
			yb = y[i-ifx-20];  z0 += x20*yb;  z1 += x21*yb;
			ya = y[i-ifx-21];  z0 += x21*ya;
			z[i+1] = z1;
			z[i] = z0;
		}
	} else if (lx==23) {
		x0 = x[ifx];
		x1 = x[ifx+1];
		x2 = x[ifx+2];
		x3 = x[ifx+3];
		x4 = x[ifx+4];
		x5 = x[ifx+5];
		x6 = x[ifx+6];
		x7 = x[ifx+7];
		x8 = x[ifx+8];
		x9 = x[ifx+9];
		x10 = x[ifx+10];
		x11 = x[ifx+11];
		x12 = x[ifx+12];
		x13 = x[ifx+13];
		x14 = x[ifx+14];
		x15 = x[ifx+15];
		x16 = x[ifx+16];
		x17 = x[ifx+17];
		x18 = x[ifx+18];
		x19 = x[ifx+19];
		x20 = x[ifx+20];
		x21 = x[ifx+21];
		x22 = x[ifx+22];
		for (i=ilow; i<=ihigh-1; i+=2) {
			ya = y[i+1-ifx];  z1 = x0*ya;
			yb = y[i-ifx];  z0 = x0*yb;  z1 += x1*yb;
			ya = y[i-ifx-1];  z0 += x1*ya;  z1 += x2*ya;
			yb = y[i-ifx-2];  z0 += x2*yb;  z1 += x3*yb;
			ya = y[i-ifx-3];  z0 += x3*ya;  z1 += x4*ya;
			yb = y[i-ifx-4];  z0 += x4*yb;  z1 += x5*yb;
			ya = y[i-ifx-5];  z0 += x5*ya;  z1 += x6*ya;
			yb = y[i-ifx-6];  z0 += x6*yb;  z1 += x7*yb;
			ya = y[i-ifx-7];  z0 += x7*ya;  z1 += x8*ya;
			yb = y[i-ifx-8];  z0 += x8*yb;  z1 += x9*yb;
			ya = y[i-ifx-9];  z0 += x9*ya;  z1 += x10*ya;
			yb = y[i-ifx-10];  z0 += x10*yb;  z1 += x11*yb;
			ya = y[i-ifx-11];  z0 += x11*ya;  z1 += x12*ya;
			yb = y[i-ifx-12];  z0 += x12*yb;  z1 += x13*yb;
			ya = y[i-ifx-13];  z0 += x13*ya;  z1 += x14*ya;
			yb = y[i-ifx-14];  z0 += x14*yb;  z1 += x15*yb;
			ya = y[i-ifx-15];  z0 += x15*ya;  z1 += x16*ya;
			yb = y[i-ifx-16];  z0 += x16*yb;  z1 += x17*yb;
			ya = y[i-ifx-17];  z0 += x17*ya;  z1 += x18*ya;
			yb = y[i-ifx-18];  z0 += x18*yb;  z1 += x19*yb;
			ya = y[i-ifx-19];  z0 += x19*ya;  z1 += x20*ya;
			yb = y[i-ifx-20];  z0 += x20*yb;  z1 += x21*yb;
			ya = y[i-ifx-21];  z0 += x21*ya;  z1 += x22*ya;
			yb = y[i-ifx-22];  z0 += x22*yb;
			z[i+1] = z1;
			z[i] = z0;
		}
	} else if (lx==24) {
		x0 = x[ifx];
		x1 = x[ifx+1];
		x2 = x[ifx+2];
		x3 = x[ifx+3];
		x4 = x[ifx+4];
		x5 = x[ifx+5];
		x6 = x[ifx+6];
		x7 = x[ifx+7];
		x8 = x[ifx+8];
		x9 = x[ifx+9];
		x10 = x[ifx+10];
		x11 = x[ifx+11];
		x12 = x[ifx+12];
		x13 = x[ifx+13];
		x14 = x[ifx+14];
		x15 = x[ifx+15];
		x16 = x[ifx+16];
		x17 = x[ifx+17];
		x18 = x[ifx+18];
		x19 = x[ifx+19];
		x20 = x[ifx+20];
		x21 = x[ifx+21];
		x22 = x[ifx+22];
		x23 = x[ifx+23];
		for (i=ilow; i<=ihigh-1; i+=2) {
			ya = y[i+1-ifx];  z1 = x0*ya;
			yb = y[i-ifx];  z0 = x0*yb;  z1 += x1*yb;
			ya = y[i-ifx-1];  z0 += x1*ya;  z1 += x2*ya;
			yb = y[i-ifx-2];  z0 += x2*yb;  z1 += x3*yb;
			ya = y[i-ifx-3];  z0 += x3*ya;  z1 += x4*ya;
			yb = y[i-ifx-4];  z0 += x4*yb;  z1 += x5*yb;
			ya = y[i-ifx-5];  z0 += x5*ya;  z1 += x6*ya;
			yb = y[i-ifx-6];  z0 += x6*yb;  z1 += x7*yb;
			ya = y[i-ifx-7];  z0 += x7*ya;  z1 += x8*ya;
			yb = y[i-ifx-8];  z0 += x8*yb;  z1 += x9*yb;
			ya = y[i-ifx-9];  z0 += x9*ya;  z1 += x10*ya;
			yb = y[i-ifx-10];  z0 += x10*yb;  z1 += x11*yb;
			ya = y[i-ifx-11];  z0 += x11*ya;  z1 += x12*ya;
			yb = y[i-ifx-12];  z0 += x12*yb;  z1 += x13*yb;
			ya = y[i-ifx-13];  z0 += x13*ya;  z1 += x14*ya;
			yb = y[i-ifx-14];  z0 += x14*yb;  z1 += x15*yb;
			ya = y[i-ifx-15];  z0 += x15*ya;  z1 += x16*ya;
			yb = y[i-ifx-16];  z0 += x16*yb;  z1 += x17*yb;
			ya = y[i-ifx-17];  z0 += x17*ya;  z1 += x18*ya;
			yb = y[i-ifx-18];  z0 += x18*yb;  z1 += x19*yb;
			ya = y[i-ifx-19];  z0 += x19*ya;  z1 += x20*ya;
			yb = y[i-ifx-20];  z0 += x20*yb;  z1 += x21*yb;
			ya = y[i-ifx-21];  z0 += x21*ya;  z1 += x22*ya;
			yb = y[i-ifx-22];  z0 += x22*yb;  z1 += x23*yb;
			ya = y[i-ifx-23];  z0 += x23*ya;
			z[i+1] = z1;
			z[i] = z0;
		}
	} else if (lx==25) {
		x0 = x[ifx];
		x1 = x[ifx+1];
		x2 = x[ifx+2];
		x3 = x[ifx+3];
		x4 = x[ifx+4];
		x5 = x[ifx+5];
		x6 = x[ifx+6];
		x7 = x[ifx+7];
		x8 = x[ifx+8];
		x9 = x[ifx+9];
		x10 = x[ifx+10];
		x11 = x[ifx+11];
		x12 = x[ifx+12];
		x13 = x[ifx+13];
		x14 = x[ifx+14];
		x15 = x[ifx+15];
		x16 = x[ifx+16];
		x17 = x[ifx+17];
		x18 = x[ifx+18];
		x19 = x[ifx+19];
		x20 = x[ifx+20];
		x21 = x[ifx+21];
		x22 = x[ifx+22];
		x23 = x[ifx+23];
		x24 = x[ifx+24];
		for (i=ilow; i<=ihigh-1; i+=2) {
			ya = y[i+1-ifx];  z1 = x0*ya;
			yb = y[i-ifx];  z0 = x0*yb;  z1 += x1*yb;
			ya = y[i-ifx-1];  z0 += x1*ya;  z1 += x2*ya;
			yb = y[i-ifx-2];  z0 += x2*yb;  z1 += x3*yb;
			ya = y[i-ifx-3];  z0 += x3*ya;  z1 += x4*ya;
			yb = y[i-ifx-4];  z0 += x4*yb;  z1 += x5*yb;
			ya = y[i-ifx-5];  z0 += x5*ya;  z1 += x6*ya;
			yb = y[i-ifx-6];  z0 += x6*yb;  z1 += x7*yb;
			ya = y[i-ifx-7];  z0 += x7*ya;  z1 += x8*ya;
			yb = y[i-ifx-8];  z0 += x8*yb;  z1 += x9*yb;
			ya = y[i-ifx-9];  z0 += x9*ya;  z1 += x10*ya;
			yb = y[i-ifx-10];  z0 += x10*yb;  z1 += x11*yb;
			ya = y[i-ifx-11];  z0 += x11*ya;  z1 += x12*ya;
			yb = y[i-ifx-12];  z0 += x12*yb;  z1 += x13*yb;
			ya = y[i-ifx-13];  z0 += x13*ya;  z1 += x14*ya;
			yb = y[i-ifx-14];  z0 += x14*yb;  z1 += x15*yb;
			ya = y[i-ifx-15];  z0 += x15*ya;  z1 += x16*ya;
			yb = y[i-ifx-16];  z0 += x16*yb;  z1 += x17*yb;
			ya = y[i-ifx-17];  z0 += x17*ya;  z1 += x18*ya;
			yb = y[i-ifx-18];  z0 += x18*yb;  z1 += x19*yb;
			ya = y[i-ifx-19];  z0 += x19*ya;  z1 += x20*ya;
			yb = y[i-ifx-20];  z0 += x20*yb;  z1 += x21*yb;
			ya = y[i-ifx-21];  z0 += x21*ya;  z1 += x22*ya;
			yb = y[i-ifx-22];  z0 += x22*yb;  z1 += x23*yb;
			ya = y[i-ifx-23];  z0 += x23*ya;  z1 += x24*ya;
			yb = y[i-ifx-24];  z0 += x24*yb;
			z[i+1] = z1;
			z[i] = z0;
		}
	} else if (lx==26) {
		x0 = x[ifx];
		x1 = x[ifx+1];
		x2 = x[ifx+2];
		x3 = x[ifx+3];
		x4 = x[ifx+4];
		x5 = x[ifx+5];
		x6 = x[ifx+6];
		x7 = x[ifx+7];
		x8 = x[ifx+8];
		x9 = x[ifx+9];
		x10 = x[ifx+10];
		x11 = x[ifx+11];
		x12 = x[ifx+12];
		x13 = x[ifx+13];
		x14 = x[ifx+14];
		x15 = x[ifx+15];
		x16 = x[ifx+16];
		x17 = x[ifx+17];
		x18 = x[ifx+18];
		x19 = x[ifx+19];
		x20 = x[ifx+20];
		x21 = x[ifx+21];
		x22 = x[ifx+22];
		x23 = x[ifx+23];
		x24 = x[ifx+24];
		x25 = x[ifx+25];
		for (i=ilow; i<=ihigh-1; i+=2) {
			ya = y[i+1-ifx];  z1 = x0*ya;
			yb = y[i-ifx];  z0 = x0*yb;  z1 += x1*yb;
			ya = y[i-ifx-1];  z0 += x1*ya;  z1 += x2*ya;
			yb = y[i-ifx-2];  z0 += x2*yb;  z1 += x3*yb;
			ya = y[i-ifx-3];  z0 += x3*ya;  z1 += x4*ya;
			yb = y[i-ifx-4];  z0 += x4*yb;  z1 += x5*yb;
			ya = y[i-ifx-5];  z0 += x5*ya;  z1 += x6*ya;
			yb = y[i-ifx-6];  z0 += x6*yb;  z1 += x7*yb;
			ya = y[i-ifx-7];  z0 += x7*ya;  z1 += x8*ya;
			yb = y[i-ifx-8];  z0 += x8*yb;  z1 += x9*yb;
			ya = y[i-ifx-9];  z0 += x9*ya;  z1 += x10*ya;
			yb = y[i-ifx-10];  z0 += x10*yb;  z1 += x11*yb;
			ya = y[i-ifx-11];  z0 += x11*ya;  z1 += x12*ya;
			yb = y[i-ifx-12];  z0 += x12*yb;  z1 += x13*yb;
			ya = y[i-ifx-13];  z0 += x13*ya;  z1 += x14*ya;
			yb = y[i-ifx-14];  z0 += x14*yb;  z1 += x15*yb;
			ya = y[i-ifx-15];  z0 += x15*ya;  z1 += x16*ya;
			yb = y[i-ifx-16];  z0 += x16*yb;  z1 += x17*yb;
			ya = y[i-ifx-17];  z0 += x17*ya;  z1 += x18*ya;
			yb = y[i-ifx-18];  z0 += x18*yb;  z1 += x19*yb;
			ya = y[i-ifx-19];  z0 += x19*ya;  z1 += x20*ya;
			yb = y[i-ifx-20];  z0 += x20*yb;  z1 += x21*yb;
			ya = y[i-ifx-21];  z0 += x21*ya;  z1 += x22*ya;
			yb = y[i-ifx-22];  z0 += x22*yb;  z1 += x23*yb;
			ya = y[i-ifx-23];  z0 += x23*ya;  z1 += x24*ya;
			yb = y[i-ifx-24];  z0 += x24*yb;  z1 += x25*yb;
			ya = y[i-ifx-25];  z0 += x25*ya;
			z[i+1] = z1;
			z[i] = z0;
		}
	} else if (lx==27) {
		x0 = x[ifx];
		x1 = x[ifx+1];
		x2 = x[ifx+2];
		x3 = x[ifx+3];
		x4 = x[ifx+4];
		x5 = x[ifx+5];
		x6 = x[ifx+6];
		x7 = x[ifx+7];
		x8 = x[ifx+8];
		x9 = x[ifx+9];
		x10 = x[ifx+10];
		x11 = x[ifx+11];
		x12 = x[ifx+12];
		x13 = x[ifx+13];
		x14 = x[ifx+14];
		x15 = x[ifx+15];
		x16 = x[ifx+16];
		x17 = x[ifx+17];
		x18 = x[ifx+18];
		x19 = x[ifx+19];
		x20 = x[ifx+20];
		x21 = x[ifx+21];
		x22 = x[ifx+22];
		x23 = x[ifx+23];
		x24 = x[ifx+24];
		x25 = x[ifx+25];
		x26 = x[ifx+26];
		for (i=ilow; i<=ihigh-1; i+=2) {
			ya = y[i+1-ifx];  z1 = x0*ya;
			yb = y[i-ifx];  z0 = x0*yb;  z1 += x1*yb;
			ya = y[i-ifx-1];  z0 += x1*ya;  z1 += x2*ya;
			yb = y[i-ifx-2];  z0 += x2*yb;  z1 += x3*yb;
			ya = y[i-ifx-3];  z0 += x3*ya;  z1 += x4*ya;
			yb = y[i-ifx-4];  z0 += x4*yb;  z1 += x5*yb;
			ya = y[i-ifx-5];  z0 += x5*ya;  z1 += x6*ya;
			yb = y[i-ifx-6];  z0 += x6*yb;  z1 += x7*yb;
			ya = y[i-ifx-7];  z0 += x7*ya;  z1 += x8*ya;
			yb = y[i-ifx-8];  z0 += x8*yb;  z1 += x9*yb;
			ya = y[i-ifx-9];  z0 += x9*ya;  z1 += x10*ya;
			yb = y[i-ifx-10];  z0 += x10*yb;  z1 += x11*yb;
			ya = y[i-ifx-11];  z0 += x11*ya;  z1 += x12*ya;
			yb = y[i-ifx-12];  z0 += x12*yb;  z1 += x13*yb;
			ya = y[i-ifx-13];  z0 += x13*ya;  z1 += x14*ya;
			yb = y[i-ifx-14];  z0 += x14*yb;  z1 += x15*yb;
			ya = y[i-ifx-15];  z0 += x15*ya;  z1 += x16*ya;
			yb = y[i-ifx-16];  z0 += x16*yb;  z1 += x17*yb;
			ya = y[i-ifx-17];  z0 += x17*ya;  z1 += x18*ya;
			yb = y[i-ifx-18];  z0 += x18*yb;  z1 += x19*yb;
			ya = y[i-ifx-19];  z0 += x19*ya;  z1 += x20*ya;
			yb = y[i-ifx-20];  z0 += x20*yb;  z1 += x21*yb;
			ya = y[i-ifx-21];  z0 += x21*ya;  z1 += x22*ya;
			yb = y[i-ifx-22];  z0 += x22*yb;  z1 += x23*yb;
			ya = y[i-ifx-23];  z0 += x23*ya;  z1 += x24*ya;
			yb = y[i-ifx-24];  z0 += x24*yb;  z1 += x25*yb;
			ya = y[i-ifx-25];  z0 += x25*ya;  z1 += x26*ya;
			yb = y[i-ifx-26];  z0 += x26*yb;
			z[i+1] = z1;
			z[i] = z0;
		}
	} else if (lx==28) {
		x0 = x[ifx];
		x1 = x[ifx+1];
		x2 = x[ifx+2];
		x3 = x[ifx+3];
		x4 = x[ifx+4];
		x5 = x[ifx+5];
		x6 = x[ifx+6];
		x7 = x[ifx+7];
		x8 = x[ifx+8];
		x9 = x[ifx+9];
		x10 = x[ifx+10];
		x11 = x[ifx+11];
		x12 = x[ifx+12];
		x13 = x[ifx+13];
		x14 = x[ifx+14];
		x15 = x[ifx+15];
		x16 = x[ifx+16];
		x17 = x[ifx+17];
		x18 = x[ifx+18];
		x19 = x[ifx+19];
		x20 = x[ifx+20];
		x21 = x[ifx+21];
		x22 = x[ifx+22];
		x23 = x[ifx+23];
		x24 = x[ifx+24];
		x25 = x[ifx+25];
		x26 = x[ifx+26];
		x27 = x[ifx+27];
		for (i=ilow; i<=ihigh-1; i+=2) {
			ya = y[i+1-ifx];  z1 = x0*ya;
			yb = y[i-ifx];  z0 = x0*yb;  z1 += x1*yb;
			ya = y[i-ifx-1];  z0 += x1*ya;  z1 += x2*ya;
			yb = y[i-ifx-2];  z0 += x2*yb;  z1 += x3*yb;
			ya = y[i-ifx-3];  z0 += x3*ya;  z1 += x4*ya;
			yb = y[i-ifx-4];  z0 += x4*yb;  z1 += x5*yb;
			ya = y[i-ifx-5];  z0 += x5*ya;  z1 += x6*ya;
			yb = y[i-ifx-6];  z0 += x6*yb;  z1 += x7*yb;
			ya = y[i-ifx-7];  z0 += x7*ya;  z1 += x8*ya;
			yb = y[i-ifx-8];  z0 += x8*yb;  z1 += x9*yb;
			ya = y[i-ifx-9];  z0 += x9*ya;  z1 += x10*ya;
			yb = y[i-ifx-10];  z0 += x10*yb;  z1 += x11*yb;
			ya = y[i-ifx-11];  z0 += x11*ya;  z1 += x12*ya;
			yb = y[i-ifx-12];  z0 += x12*yb;  z1 += x13*yb;
			ya = y[i-ifx-13];  z0 += x13*ya;  z1 += x14*ya;
			yb = y[i-ifx-14];  z0 += x14*yb;  z1 += x15*yb;
			ya = y[i-ifx-15];  z0 += x15*ya;  z1 += x16*ya;
			yb = y[i-ifx-16];  z0 += x16*yb;  z1 += x17*yb;
			ya = y[i-ifx-17];  z0 += x17*ya;  z1 += x18*ya;
			yb = y[i-ifx-18];  z0 += x18*yb;  z1 += x19*yb;
			ya = y[i-ifx-19];  z0 += x19*ya;  z1 += x20*ya;
			yb = y[i-ifx-20];  z0 += x20*yb;  z1 += x21*yb;
			ya = y[i-ifx-21];  z0 += x21*ya;  z1 += x22*ya;
			yb = y[i-ifx-22];  z0 += x22*yb;  z1 += x23*yb;
			ya = y[i-ifx-23];  z0 += x23*ya;  z1 += x24*ya;
			yb = y[i-ifx-24];  z0 += x24*yb;  z1 += x25*yb;
			ya = y[i-ifx-25];  z0 += x25*ya;  z1 += x26*ya;
			yb = y[i-ifx-26];  z0 += x26*yb;  z1 += x27*yb;
			ya = y[i-ifx-27];  z0 += x27*ya;
			z[i+1] = z1;
			z[i] = z0;
		}
	} else if (lx==29) {
		x0 = x[ifx];
		x1 = x[ifx+1];
		x2 = x[ifx+2];
		x3 = x[ifx+3];
		x4 = x[ifx+4];
		x5 = x[ifx+5];
		x6 = x[ifx+6];
		x7 = x[ifx+7];
		x8 = x[ifx+8];
		x9 = x[ifx+9];
		x10 = x[ifx+10];
		x11 = x[ifx+11];
		x12 = x[ifx+12];
		x13 = x[ifx+13];
		x14 = x[ifx+14];
		x15 = x[ifx+15];
		x16 = x[ifx+16];
		x17 = x[ifx+17];
		x18 = x[ifx+18];
		x19 = x[ifx+19];
		x20 = x[ifx+20];
		x21 = x[ifx+21];
		x22 = x[ifx+22];
		x23 = x[ifx+23];
		x24 = x[ifx+24];
		x25 = x[ifx+25];
		x26 = x[ifx+26];
		x27 = x[ifx+27];
		x28 = x[ifx+28];
		for (i=ilow; i<=ihigh-1; i+=2) {
			ya = y[i+1-ifx];  z1 = x0*ya;
			yb = y[i-ifx];  z0 = x0*yb;  z1 += x1*yb;
			ya = y[i-ifx-1];  z0 += x1*ya;  z1 += x2*ya;
			yb = y[i-ifx-2];  z0 += x2*yb;  z1 += x3*yb;
			ya = y[i-ifx-3];  z0 += x3*ya;  z1 += x4*ya;
			yb = y[i-ifx-4];  z0 += x4*yb;  z1 += x5*yb;
			ya = y[i-ifx-5];  z0 += x5*ya;  z1 += x6*ya;
			yb = y[i-ifx-6];  z0 += x6*yb;  z1 += x7*yb;
			ya = y[i-ifx-7];  z0 += x7*ya;  z1 += x8*ya;
			yb = y[i-ifx-8];  z0 += x8*yb;  z1 += x9*yb;
			ya = y[i-ifx-9];  z0 += x9*ya;  z1 += x10*ya;
			yb = y[i-ifx-10];  z0 += x10*yb;  z1 += x11*yb;
			ya = y[i-ifx-11];  z0 += x11*ya;  z1 += x12*ya;
			yb = y[i-ifx-12];  z0 += x12*yb;  z1 += x13*yb;
			ya = y[i-ifx-13];  z0 += x13*ya;  z1 += x14*ya;
			yb = y[i-ifx-14];  z0 += x14*yb;  z1 += x15*yb;
			ya = y[i-ifx-15];  z0 += x15*ya;  z1 += x16*ya;
			yb = y[i-ifx-16];  z0 += x16*yb;  z1 += x17*yb;
			ya = y[i-ifx-17];  z0 += x17*ya;  z1 += x18*ya;
			yb = y[i-ifx-18];  z0 += x18*yb;  z1 += x19*yb;
			ya = y[i-ifx-19];  z0 += x19*ya;  z1 += x20*ya;
			yb = y[i-ifx-20];  z0 += x20*yb;  z1 += x21*yb;
			ya = y[i-ifx-21];  z0 += x21*ya;  z1 += x22*ya;
			yb = y[i-ifx-22];  z0 += x22*yb;  z1 += x23*yb;
			ya = y[i-ifx-23];  z0 += x23*ya;  z1 += x24*ya;
			yb = y[i-ifx-24];  z0 += x24*yb;  z1 += x25*yb;
			ya = y[i-ifx-25];  z0 += x25*ya;  z1 += x26*ya;
			yb = y[i-ifx-26];  z0 += x26*yb;  z1 += x27*yb;
			ya = y[i-ifx-27];  z0 += x27*ya;  z1 += x28*ya;
			yb = y[i-ifx-28];  z0 += x28*yb;
			z[i+1] = z1;
			z[i] = z0;
		}
	} else if (lx==30) {
		x0 = x[ifx];
		x1 = x[ifx+1];
		x2 = x[ifx+2];
		x3 = x[ifx+3];
		x4 = x[ifx+4];
		x5 = x[ifx+5];
		x6 = x[ifx+6];
		x7 = x[ifx+7];
		x8 = x[ifx+8];
		x9 = x[ifx+9];
		x10 = x[ifx+10];
		x11 = x[ifx+11];
		x12 = x[ifx+12];
		x13 = x[ifx+13];
		x14 = x[ifx+14];
		x15 = x[ifx+15];
		x16 = x[ifx+16];
		x17 = x[ifx+17];
		x18 = x[ifx+18];
		x19 = x[ifx+19];
		x20 = x[ifx+20];
		x21 = x[ifx+21];
		x22 = x[ifx+22];
		x23 = x[ifx+23];
		x24 = x[ifx+24];
		x25 = x[ifx+25];
		x26 = x[ifx+26];
		x27 = x[ifx+27];
		x28 = x[ifx+28];
		x29 = x[ifx+29];
		for (i=ilow; i<=ihigh-1; i+=2) {
			ya = y[i+1-ifx];  z1 = x0*ya;
			yb = y[i-ifx];  z0 = x0*yb;  z1 += x1*yb;
			ya = y[i-ifx-1];  z0 += x1*ya;  z1 += x2*ya;
			yb = y[i-ifx-2];  z0 += x2*yb;  z1 += x3*yb;
			ya = y[i-ifx-3];  z0 += x3*ya;  z1 += x4*ya;
			yb = y[i-ifx-4];  z0 += x4*yb;  z1 += x5*yb;
			ya = y[i-ifx-5];  z0 += x5*ya;  z1 += x6*ya;
			yb = y[i-ifx-6];  z0 += x6*yb;  z1 += x7*yb;
			ya = y[i-ifx-7];  z0 += x7*ya;  z1 += x8*ya;
			yb = y[i-ifx-8];  z0 += x8*yb;  z1 += x9*yb;
			ya = y[i-ifx-9];  z0 += x9*ya;  z1 += x10*ya;
			yb = y[i-ifx-10];  z0 += x10*yb;  z1 += x11*yb;
			ya = y[i-ifx-11];  z0 += x11*ya;  z1 += x12*ya;
			yb = y[i-ifx-12];  z0 += x12*yb;  z1 += x13*yb;
			ya = y[i-ifx-13];  z0 += x13*ya;  z1 += x14*ya;
			yb = y[i-ifx-14];  z0 += x14*yb;  z1 += x15*yb;
			ya = y[i-ifx-15];  z0 += x15*ya;  z1 += x16*ya;
			yb = y[i-ifx-16];  z0 += x16*yb;  z1 += x17*yb;
			ya = y[i-ifx-17];  z0 += x17*ya;  z1 += x18*ya;
			yb = y[i-ifx-18];  z0 += x18*yb;  z1 += x19*yb;
			ya = y[i-ifx-19];  z0 += x19*ya;  z1 += x20*ya;
			yb = y[i-ifx-20];  z0 += x20*yb;  z1 += x21*yb;
			ya = y[i-ifx-21];  z0 += x21*ya;  z1 += x22*ya;
			yb = y[i-ifx-22];  z0 += x22*yb;  z1 += x23*yb;
			ya = y[i-ifx-23];  z0 += x23*ya;  z1 += x24*ya;
			yb = y[i-ifx-24];  z0 += x24*yb;  z1 += x25*yb;
			ya = y[i-ifx-25];  z0 += x25*ya;  z1 += x26*ya;
			yb = y[i-ifx-26];  z0 += x26*yb;  z1 += x27*yb;
			ya = y[i-ifx-27];  z0 += x27*ya;  z1 += x28*ya;
			yb = y[i-ifx-28];  z0 += x28*yb;  z1 += x29*yb;
			ya = y[i-ifx-29];  z0 += x29*ya;
			z[i+1] = z1;
			z[i] = z0;
		}
	}
	if (ihigh>=ilow && (ihigh-ilow)%2==0) {
		ilow = ihigh;
		jlow = ifx;
		jhigh = ilx;
		for (i=ilow; i<=ihigh; ++i) {
			for (j=jlow,sum=0.0; j<=jhigh; ++j)
				sum += x[j]*y[i-j];
			z[i] = sum;
		}
	}
	
	/* ROLLING OFF:  ily+ifx < i <= ily+ilx */
	ilow = ily+ifx+1;  if (ilow<ifz) ilow = ifz;
	ihigh = ily+ilx;  if (ihigh>ilz) ihigh = ilz;
	jlow = ilow-ily;
	jhigh = ilx;
	for (i=ilow; i<=ihigh; ++i,++jlow) {
		for (j=jlow,sum=0.0; j<=jhigh; ++j)
			sum += x[j]*y[i-j];
		z[i] = sum;
	}
	
	/* OFF RIGHT:  ily+ilx < i */
	ilow = ily+ilx+1;  if (ilow<ifz) ilow = ifz;
	ihigh = ilz;
	for (i=ilow; i<=ihigh; ++i)
		z[i] = 0.0;
}
#endif

void gaussian1d_smoothing (int ns, int nsr, float *data)
/******************************************************************************
Input:
ns		number of samples in the input data 
nsr		width (in samples) of the gaussian for which 
		amplitude > 0.5*max amplitude
data		1-D array[ns] of data to smooth

Output:
data		1-D array[ns] of smoothed data
******************************************************************************/
{
	int is;				/* loop counter */
	float sum=0.0;
	float fcut;
	float r;
	float fcutr=1.0/nsr;
	static int n;
	static int mean;
	static float fcutl;
	static float s[401];		/* smoothing filter array */
	float *temp;			/* temporary array */

	/* allocate space */
	temp=(float *)malloc(ns*sizeof(float));

	/* save input fcut */
	fcut=fcutr;

	/* don't smooth if nsr equal to zero */
	if (nsr==0 || ns<=1) return;

	/* if halfwidth more than 100 samples, truncate */
	if (nsr>100) fcut=1.0/100;

	/* initialize smoothing function if not the same as the last one used */
	if (fcut !=fcutl) {
		fcutl=fcut;

		/* set span of 3, at width of 1.5*exp(-PI*1.5**2)=1/1174 */
		n=3.0/fcut+0.5;
		n=2*n/2+1;		/* make it odd for symmetry */

		/* mean is the index of the zero in the smoothing wavelet */
		mean=n/2;

		/* s(n) is the smoothing gaussian */
		for (is=1; is<=n; is++) {
			r=is-mean-1;
			r= -r*r*fcut*fcut*pi;
			s[is-1]=exp(r);
		}
			
		/* normalize to unit area, will preserve DC frequency at full 
		amplitude. Frequency at fcut will be half amplitude */
		for (is=0; is<n; is++) sum +=s[is];
		for (is=0; is<n; is++) s[is] /=sum;
	}

	/* convolve by gaussian into buffer */
	if (1.01/fcutr>(float)ns) {

		/* replace drastic smoothing by averaging */
		sum=0.0;
		for (is=0; is<ns; is++) sum +=data[is];
		sum /=ns;
		for (is=0; is<ns; is++) data[is]=sum;	

	} else {

		/* convolve with gaussian */
		convolve_cwp (n, -mean, s, ns, -mean, data, ns, -mean, temp);

		/* copy filtered data back to output array */
		for (is=0; is<ns; is++) data[is]=temp[is];
	}

	/* free allocated space */
	free(temp);
}
void interpvel1(float *vp0, float *vp, int nx0, int nz0, int nx, int nz, int lx, int lz)
{
	int ix,iz,il;
	float *vel=(float *)malloc(nx*nz0*sizeof(float));
	// interpolate model along x-direction
	#pragma omp parallel for default(shared) private(it,ix,il)
	for (iz=0;iz<nz0;iz++)
	{
		for (ix=0;ix<nx0-1;ix++)
		{
			for (il=0;il<lx;il++)
			{
				vel[iz*nx+ix*lx+il] = ((lx-il)*vp0[iz*nx0+ix] + il*vp0[iz*nx0+ix+1])/lx;
			}
		}
		vel[iz*nx+(nx0-1)*lx] = vp0[iz*nx0+nx0-1];
	}
	// interpolate model along z-direction
	#pragma omp parallel for default(shared) private(it,ix,il)
	for (ix=0;ix<nx;ix++)
	{
		for (iz=0;iz<nz0-1;iz++)
		{
			for (il=0;il<lz;il++)
			{
				vp[(iz*lz+il)*nx+ix] = ((lz-il)*vel[iz*nx+ix] + il*vel[(iz+1)*nx+ix])/lz;
			}
		}
		vp[(nz0-1)*lz*nx+ix] = vel[(nz0-1)*nx+ix];
	}
	free(vel);
}
void interpvel2(float **vp0, float **vp, int nx0, int nz0, int nx, int nz, int lx, int lz)
{
	int ix,iz,il;
	float **vel;
	vel = Creat2dArray(nz0,nx);
	// interpolate model along x-direction
	#pragma omp parallel for default(shared) private(it,ix,il)
	for (iz=0;iz<nz0;iz++)
	{
		for (ix=0;ix<nx0-1;ix++)
		{
			for (il=0;il<lx;il++)
			{
				vel[iz][ix*lx+il] = ((lx-il)*vp0[iz][ix] + il*vp0[iz][ix+1])/lx;
			}
		}
		vel[iz][(nx0-1)*lx] = vp0[iz][nx0-1];
	}
	// interpolate model along z-direction
	#pragma omp parallel for default(shared) private(it,ix,il)
	for (ix=0;ix<nx;ix++)
	{
		for (iz=0;iz<nz0-1;iz++)
		{
			for (il=0;il<lz;il++)
			{
				vp[iz*lz+il][ix] = ((lz-il)*vel[iz][ix] + il*vel[iz+1][ix])/lz;
			}
		}
		vp[(nz0-1)*lz][ix] = vel[nz0-1][ix];
	}
	free2dArray(vel,nz0,nx);
}
void extendvel1(float *vp, float *vp0, int nx, int nz, int npml)
{
	int ix,iz;
	int nxpad,nzpad;
	nxpad = nx + 2*npml;
	nzpad = nz + 2*npml;
	for (iz=0;iz<nz;iz++){
		for (ix=0;ix<nx;ix++)
			vp[(iz+npml)*nxpad + ix+npml] = vp0[iz*nx + ix];}
	// extend vel to pml layer along x-direction
	for (iz=npml;iz<nzpad-npml;iz++){
		for (ix=0;ix<npml;ix++){
			vp[iz*nxpad+ix] = vp[iz*nxpad+npml];
			vp[iz*nxpad+nxpad-1-ix] = vp[iz*nxpad+nxpad-npml-1];}}
	// extend vel to pml layer along z-direction
	for (ix=0;ix<nxpad;ix++){
		for (iz=0;iz<npml;iz++){
			vp[iz*nxpad+ix] = vp[npml*nxpad+ix];
			vp[(nzpad-1-iz)*nxpad+ix] = vp[(nzpad-1-npml)*nxpad+ix];}}
}
void extendvel2(float **vp, float **vp0, int nx, int nz, int npml)
{
	int ix,iz;
	int nxpad,nzpad;
	nxpad = nx + 2*npml;
	nzpad = nz + 2*npml;
	for (iz=0;iz<nz;iz++){
		for (ix=0;ix<nx;ix++)
			vp[iz+npml][ix+npml] = vp0[iz][ix];}
	// extend vel to pml layer along x-direction
	for (iz=npml;iz<nzpad-npml;iz++){
		for (ix=0;ix<npml;ix++){
			vp[iz][ix] = vp[iz][npml];
			vp[iz][nxpad-1-ix] = vp[iz][nxpad-npml-1];}}
	// extend vel to pml layer along z-direction
	for (ix=0;ix<nxpad;ix++){
		for (iz=0;iz<npml;iz++){
			vp[iz][ix] = vp[npml][ix];
			vp[nzpad-1-iz][ix] = vp[nzpad-1-npml][ix];}}
}
void extractvel1(float *temp, float *vp, int nx, int nz, int nx1, int nx2)
{
	int ix,iz;
	int nxlength = nx2-nx1+1;
	for (iz=0;iz<nz;iz++)
	{
		for (ix=nx1;ix<=nx2;ix++)
		{
			temp[iz*nxlength+ix-nx1]=vp[iz*nx+ix]*vp[iz*nx+ix];
		}
	}
}
void extractrho1(float *temp, float *rho, int nx, int nz, int nx1, int nx2)
{
	int ix,iz;
	int nxlength = nx2-nx1+1;
	for (iz=0;iz<nz;iz++)
	{
		for (ix=nx1;ix<=nx2;ix++)
		{
			temp[iz*nxlength+ix-nx1]=rho[iz*nx+ix];
		}
	}
}
void extractvel2(float **temp, float *vp, int nx, int nz, int nx1, int nx2)
{
	int ix,iz;
	for (iz=0;iz<nz;iz++)
	{
		for (ix=nx1;ix<=nx2;ix++)
		{
			temp[iz][ix-nx1]=vp[iz*nx+ix]*vp[iz*nx+ix];
		}
	}
}
void extractrho2(float **temp, float *rho, int nx, int nz, int nx1, int nx2)
{
	int ix,iz;
	for (iz=0;iz<nz;iz++)
	{
		for (ix=nx1;ix<=nx2;ix++)
		{
			temp[iz][ix-nx1]=rho[iz*nx+ix];
		}
	}
}
