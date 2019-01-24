#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include "function.h"
#include "common.h"


void scaleimage1d(float *illum, float *grad, int n)
{
	int i;
	float maxillum = 0.0001*absMaxval1(illum, n);
	#pragma omp parallel for default(shared) private(i)
	for (i=0; i<n; i++)
	{
		grad[i] /= illum[i] + maxillum;
	}
}
void scaleimage2d(float **illum, float **grad, int nz, int nx)
{
	int ix,iz;
	float maxillum = 0.0001*absMaxval2(illum,nz,nx);
	#pragma omp parallel for default(shared) private(ix,iz)
	for (iz=0; iz<nz; iz++)
	{
		for (ix=0; ix<nx; ix++)
		{
			grad[iz][ix] /= illum[iz][ix] + maxillum;
		}
	}
}
void illumsmooth1d(float *Illum, int nx, int nz, int nsr)
{
	for (int iz=0;iz<nz;iz++)
		gaussian1d_smoothing(nx, nsr, &Illum[iz*nx]);
}
void illumsmooth2d(float **Illum, int nx, int nz, int nsr)
{
	for (int iz=0;iz<nz;iz++)
		gaussian1d_smoothing(nx, nsr, Illum[iz]);
}
void laplacefd1d(int nx, int nz, int dx, int dz, float *image)
{
	int ix,iz;

	float **image1,temp1,temp2;
	image1=Creat2dArray(nz,nx);

	for (ix = 0; ix < nx; ix++)
	{
		for (iz = 0; iz < nz; iz++)
		{
			image1[iz][ix] = image[iz*nx+ix];
		}
	}
	for (ix = 1; ix < nx-1; ix++)
	{
		for (iz = 1; iz < nz-1; iz++)
		{
			temp1 = (image1[iz][ix+1] + image1[iz][ix-1] - 2.0*image1[iz][ix])/(dx*dx);
			temp2 = (image1[iz+1][ix] + image1[iz-1][ix] - 2.0*image1[iz][ix])/(dz*dz);

			image[iz*nx+ix] = (temp1 + temp2);
		}
	}
	for (ix = 0; ix < nx; ix++)
	{
		image[ix] = 0.0;
		image[(nz-1)*nx+ix] = 0.0;
	}
	for (iz = 0; iz < nz; iz++)
	{
		image[iz*nx] = 0.0;
		image[iz*nx+nx-1] = 0.0;
	}
	free2dArray(image1,nz,nx);
}
void laplacefd2d(int nx, int nz, int dx, int dz, float **image)
{
	int ix,iz;

	float **image1,temp1,temp2;
	image1=Creat2dArray(nz,nx);

	for (ix = 0; ix < nx; ix++)
	{
		for (iz = 0; iz < nz; iz++)
		{
			image1[iz][ix] = image[iz][ix];
		}
	}
	for (ix = 1; ix < nx-1; ix++)
	{
		for (iz = 1; iz < nz-1; iz++)
		{
			temp1 = (image1[iz][ix+1] + image1[iz][ix-1] - 2.0*image1[iz][ix])/(dx*dx);
			temp2 = (image1[iz+1][ix] + image1[iz-1][ix] - 2.0*image1[iz][ix])/(dz*dz);

			image[iz][ix] = (temp1 + temp2);
		}
	}
	for (ix = 0; ix < nx; ix++)
	{
		image[0][ix] = 0.0;
		image[nz-1][ix] = 0.0;
	}
	for (iz = 0; iz < nz; iz++)
	{
		image[iz][0] = 0.0;
		image[iz][nx-1] = 0.0;
	}
	free2dArray(image1,nz,nx);
}
void taperh1d(int nx, int nz, int nwinlen, float gd, float *g)
{
	int j,k;
	float *wz=(float *)malloc(nz*sizeof(float));
	
	fmemset1v(wz, nz, 1.0);
	for (j=0;j<nwinlen;j++)
	{
		wz[j]=exp(-0.5*gd*gd*(j-nwinlen)*(j-nwinlen)/(4.0*nwinlen*nwinlen))                                  ;
	}

	for(int ix=0;ix<nx;ix++)
	{
		for (int iz=0;iz<nz;iz++)
		{
			g[iz*nx+ix]*=wz[iz];
		}
	}
	free(wz);
}
void taperh2d(int nx, int nz, int nwinlen, float gd, float **g)
{
	int j,k;
	float *wz=(float *)malloc(nz*sizeof(float));
	
	fmemset1v(wz, nz, 1.0);
	for (j=0;j<nwinlen;j++)
	{
		wz[j]=exp(-0.5*gd*gd*(j-nwinlen)*(j-nwinlen)/(4.0*nwinlen*nwinlen));
	}

	for(int ix=0;ix<nx;ix++)
	{
		for (int iz=0;iz<nz;iz++)
		{
			g[iz][ix] *= wz[iz];
		}
	}
	free(wz);
}

