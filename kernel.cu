#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <cuda.h>
#include <cuda_runtime.h>

#include "common.h"



__global__ void cuda_ricker_wavelet(float *d_wavelet, float fdom, float dt, int nt)
// generate ricker wavelet with time deley
// grid:block ((nt+511)/512,512)
{
	int it = threadIdx.x + blockDim.x*blockIdx.x;
	float temp = pi*fdom*fabs(it*dt - 1.0/fdom);
	temp *=temp;
	if (it < nt){
          d_wavelet[it] = (1.0 - 2.0*temp)*expf(-temp);}
}
__global__ void cuda_source(float *d_source, int nsx, int nsz, int nxlength, int nz, float amp, float alp, float dx2, float dz2)
// calculate source: source = amp*exp(-alp*alp*dist^2)
// grid:block <(nz+blocksizez-1)/blocksizez,(nz+blocksizex-1)/blocksizex),(blocksizez,blocksizex)>
{
	int idz = threadIdx.x + blockDim.x * blockIdx.x;
	int idx = threadIdx.y + blockDim.y * blockIdx.y;
	int id = idz*nxlength + idx;
	float dist;

	if (idx < nxlength && idz < nz){
          dist = (idx - nsx)*(idx - nsx)*dx2 +  (idz - nsz)*(idz - nsz)*dz2;
          d_source[id] = amp*expf(-alp*alp*dist);}
}
__global__ void cuda_add_source(float *d_p, float *d_source, float *d_wlt, float dt, int add, int nxlength, int nz, int it)
// add = 1 , add the source; add = 0, substract the source
// grid:block <(nz+blocksizez-1)/blocksizez,(nz+blocksizex-1)/blocksizex),(blocksizez,blocksizex)>
{
	int idz = threadIdx.x + blockDim.x*blockIdx.x;
	int idx = threadIdx.y + blockDim.y*blockIdx.y;
	int id = idz*nxlength + idx;

	if (idx < nxlength && idz < nz)
		if (add == 1)     d_p[id] += dt*d_source[id]*d_wlt[it];
		else              d_p[id] -= dt*d_source[id]*d_wlt[it];
}
__global__ void cuda_record1(float *d_p, float *d_seis, int npml, int nxlength)
// extract seismogram at time it
// grid:block ((nxlength-2*npml+511)/512,512)
{
	int idx = threadIdx.x + blockDim.x * blockIdx.x;
	int nxlocal = nxlength;
	if (idx < nxlocal)      d_seis[idx] = d_p[npml*(nxlength+2*npml) + idx+npml];
}
__global__ void cuda_record2(float *d_p, float *d_seis, int npml, int nxlength, int noffset)
// extract seismogram at time it
// grid:block ((NX+511)/512,512)
{
	int idx = threadIdx.x + blockDim.x * blockIdx.x;
	int nxlocal = nxlength - noffset;
	if (idx < nxlocal)      d_seis[idx] = d_p[npml*(nxlength+2*npml) + idx+npml+noffset];
}
__global__ void cuda_record3(float *d_p, float *d_seis, int npml, int nxlength, int noffset)
// extract seismogram at time it
// grid:block ((NX+511)/512,512)
{
	int idx = threadIdx.x + blockDim.x * blockIdx.x;
	int nxlocal = nxlength - noffset;
	if (idx < nxlocal)      d_seis[idx] = d_p[npml*(nxlength+2*npml) + idx+npml];
}
__global__ void cuda_insert_record1(float *d_p, float *d_seis, int npml, int nxlength, float dt)
// insert seismogram at time it
// grid:block ((nxlength-2*npml+511)/512,512)
{
	int idx = threadIdx.x + blockDim.x * blockIdx.x;
	int nxlocal = nxlength;
	if (idx < nxlocal)      d_p[npml*(nxlength+2*npml)+idx+npml] = dt*d_seis[idx];//d_p[npml*nxlength + idx+npml] += dt*d_seis[idx];
}
__global__ void cuda_insert_record2(float *d_p, float *d_seis, int npml, int nxlength, int noffset, float dt)
// insert seismogram at time it
// grid:block ((NX+511)/512,512)
{
	int idx = threadIdx.x + blockDim.x * blockIdx.x;
	int nxlocal = nxlength-noffset;
	if (idx < nxlocal)     d_p[npml*(nxlength+2*npml) + idx+npml+noffset]  = dt*d_seis[idx];//d_p[npml*nxlength + idx+npml+noffset]  += dt*d_seis[idx];
}
__global__ void cuda_insert_record3(float *d_p, float *d_seis, int npml, int nxlength, int noffset, float dt)
// insert seismogram at time it
// grid:block ((NX+511)/512,512)
{
	int idx = threadIdx.x + blockDim.x * blockIdx.x;
	int nxlocal = nxlength-noffset;
	if (idx < nxlocal)     d_p[npml*(nxlength+2*npml) + idx+npml]  = dt*d_seis[idx];//d_p[npml*nxlength + idx+npml+noffset]  += dt*d_seis[idx];
}
__global__ void cuda_mute1(float *d_seis, float *d_vp, int nsx, int nsz, int nt, int npml, int nxlength, int nw, int tlength, float fdom, float dx2, float dz2, float _dt)
// < mute direct wave >
// grid:block <(nt+blocksizez-1)/blocksizez,(nxlength-2*npml+blocksizex-1)/blocksizex),(blocksizez,blocksizex)>
{
	int idt = threadIdx.x + blockDim.x * blockIdx.x;
	int idx = threadIdx.y + blockDim.y * blockIdx.y;
	int nxlocal = nxlength;	
	int id = idt*nxlocal + idx;
	int cut,temp;
	
      float vtemp,dist;

	if (idt < nt && idx < nxlocal)
	{
		vtemp = 0.482*(1.0/sqrtf(d_vp[npml*(nxlength+2*npml) + idx + npml]) + 1.0/sqrtf(d_vp[nsz*(nxlength+2*npml) + nsx]));
		dist = (idx + npml - nsx)*(idx + npml - nsx)*dx2 + (npml - nsz)*(npml - nsz)*dz2;
	
		dist = sqrtf(dist);
		temp = tlength<nt?tlength:nt;
		cut = ((int)(dist*vtemp*_dt)+nw)<nt?((int)(dist*vtemp*_dt)+nw):nt;
		cut = cut>temp?cut:temp;
		cut = cut<nt?cut:nt;

		if (idt <= cut && idx < nxlocal)    d_seis[id] = 0.0;
	}
}
__global__ void cuda_mute2(float *d_seis, float *d_vp, int nsx, int nsz, int nt, int npml, int nxlength, int noffset, int nw, int tlength, float fdom, float dx2, float dz2, float _dt)
// < mute direct wave >
// grid:block <(nt+blocksizez-1)/blocksizez,(NX+blocksizex-1)/blocksizex),(blocksizez,blocksizex)>
{
	int idt = threadIdx.x + blockDim.x * blockIdx.x;
	int idx = threadIdx.y + blockDim.y * blockIdx.y;
	int nxlocal = nxlength-noffset;	
	int id = idt*nxlocal + idx;
	int cut,temp;
	
      float vtemp,dist;

	if (idt < nt && idx < nxlocal)
	{
		vtemp = 0.482*(1.0/sqrtf(d_vp[npml*(nxlength+2*npml) + idx + npml + noffset]) + 1.0/sqrtf(d_vp[nsz*(nxlength+2*npml) + nsx]));
		dist = (idx + npml + noffset - nsx)*(idx + npml + noffset - nsx)*dx2 + (npml - nsz)*(npml - nsz)*dz2;
	
		dist = sqrtf(dist);
		temp = tlength<nt?tlength:nt;
		cut = ((int)(dist*vtemp*_dt)+nw)<nt?((int)(dist*vtemp*_dt)+nw):nt;
		cut = cut>temp?cut:temp;
		cut = cut<nt?cut:nt;

		if (idt <= cut && idx < nxlocal)    d_seis[id] = 0.0;
	}
}
__global__ void cuda_mute3(float *d_seis, float *d_vp, int nsx, int nsz, int nt, int npml, int nxlength, int noffset, int nw, int tlength, float fdom, float dx2, float dz2, float _dt)
// < mute direct wave >
// grid:block <(nt+blocksizez-1)/blocksizez,(NX+blocksizex-1)/blocksizex),(blocksizez,blocksizex)>
{
	int idt = threadIdx.x + blockDim.x * blockIdx.x;
	int idx = threadIdx.y + blockDim.y * blockIdx.y;
	int nxlocal = nxlength-noffset;	
	int id = idt*nxlocal + idx;
	int cut,temp;
	
      float vtemp,dist;

	if (idt < nt && idx < nxlocal)
	{
		vtemp = 0.482*(1.0/sqrtf(d_vp[npml*(nxlength+2*npml) + idx + npml]) + 1.0/sqrtf(d_vp[nsz*(nxlength+2*npml) + nsx]));
		dist = (idx + npml - nsx)*(idx + npml - nsx)*dx2 + (npml - nsz)*(npml - nsz)*dz2;
	
		dist = sqrtf(dist);
		temp = tlength<nt?tlength:nt;
		cut = ((int)(dist*vtemp*_dt)+nw)<nt?((int)(dist*vtemp*_dt)+nw):nt;
		cut = cut>temp?cut:temp;
		cut = cut<nt?cut:nt;

		if (idt <= cut && idx < nxlocal)    d_seis[id] = 0.0;
	}
}
__global__ void cuda_pmlCoeffpx(float *d_ddx, float vpmax, float dx, int npml, int nxlength)
// < pml absorb coefficients >
// grid:block ((nxlength+511)/512,512)
{
	int idx = threadIdx.x + blockDim.x * blockIdx.x;
	
	float R = 1.0e-10;
	float d = npml*dx;
	float d0 = 1.5*vpmax*logf(1/R)/d;
	float temp;
	
	if(idx < npml){
		temp = (float)(npml - idx);
		temp = temp/npml;
		temp = temp*temp;
		d_ddx[idx] = d0*temp;}
	if (idx >= npml && idx < nxlength-npml)   
		d_ddx[idx] = 0.0;
	if (idx >= nxlength-npml && idx < nxlength){                       
		temp = (float)(idx-nxlength+npml+1);
		temp = temp/npml;
		temp = temp*temp;
		d_ddx[idx] = d0*temp;}
}
__global__ void cuda_pmlCoeffpz(float *d_ddz, float vpmax, float dz, int npml, int nz)
// < pml absorb coefficients >
// grid:block ((nz+511)/512,512)
{
	int idx = threadIdx.x + blockDim.x * blockIdx.x;
	
	float R = 1.0e-10;
	float d = npml*dz;
	float d0 = 1.5*vpmax*logf(1/R)/d;
	float temp;

	if(idx < npml){
		temp = (float)(npml - idx);
		temp = temp/npml;
		temp = temp*temp;
		d_ddz[idx] = d0*temp;}
	if (idx >= npml && idx <nz-npml)   
		d_ddz[idx] = 0.0;
	if (idx >= nz-npml && idx < nz){                   
		temp = (float)(idx-nz+npml+1);
		temp = temp/npml;
		temp = temp*temp;
		d_ddz[idx] = d0*temp;}
}
__global__ void cuda_pmlCoeffvx(float *d_ddxVx, float vpmax, float dx, int npml, int nxlength)
// < pml absorb coefficients >
// grid:block ((nxlength-1+511)/512,512)   call para:nxlength -> nxlength-1
{
	int idx = threadIdx.x + blockDim.x * blockIdx.x;
	
	float R = 1.0e-10;
	float d = npml*dx;
	float d0 = 1.5*vpmax*logf(1/R)/d;
	float temp;
	
	if(idx<npml){                                 
		temp = (float)(npml - idx) - 0.5;
		temp = temp/npml;
		temp = temp*temp;
		d_ddxVx[idx] = d0*temp;}
	if (idx >= npml && idx < nxlength-npml)                
		d_ddxVx[idx] = 0.0;
	if (idx >= nxlength-npml && idx < nxlength){                   
		temp = (float)(idx-nxlength+npml+1) - 0.5;
		temp = temp/npml;
		temp = temp*temp;
		d_ddxVx[idx] = d0*temp;}
}
__global__ void cuda_pmlCoeffvz(float *d_ddzVz, float vpmax, float dz, int npml, int nz)
// < pml absorb coefficients >
// grid:block ((nz-1+511)/512,512)   call para:nz -> nz-1 
{
	int idx = threadIdx.x + blockDim.x * blockIdx.x;
	
	float R = 1.0e-10;
	float d = npml*dz;
	float d0 = 1.5*vpmax*logf(1/R)/d;
	float temp;

	if(idx < npml){
		temp = (float)(npml - idx) - 0.5;
		temp = temp/npml;
		temp = temp*temp;
		d_ddzVz[idx] = d0*temp;}
	if (idx >= npml && idx < nz-npml)                 
		d_ddzVz[idx] = 0.0;
	if (idx >= nz-npml && idx < nz){                      
		temp = (float)(idx-nz+npml+1) - 0.5;
		temp = temp/npml;
		temp = temp*temp;
		d_ddzVz[idx] = d0*temp;}
}
//__global__ void cuda_init_norder(int **d_norder, int nxlength, int nz, int N)
__global__ void cuda_norder(int *d_norder, int nxlength, int nz)
// < set P wavefield difference order of each grid >
// grid:block <(nz+blocksizez-1)/blocksizez,(nxlength-2*npml+blocksizex-1)/blocksizex),(blocksizez,blocksizex)>
{
	int idz = threadIdx.x + blockDim.x * blockIdx.x;
	int idx = threadIdx.y + blockDim.y * blockIdx.y;
	int id = idz*nxlength + idx;

	if (idz >= 1 && idz <N/2)
	{
		if (idx >= 1 && idx < N/2)                  		  d_norder[id] = 2*(idx<idz?idx:idz);
		if (idx >= N/2 && idx < nxlength-N/2)  	  	  	  d_norder[id] = 2*idz;
		if (idx >= nxlength-N/2 && idx < nxlength-1)       	  d_norder[id] = 2*(idz<(nxlength-idx-1)?idz:(nxlength-idx-1));
	}
	if (idz >= N/2 && idz < nz-N/2)
	{
		if (idx >= 1 && idx < N/2)                  		  d_norder[id] = 2*idx;
		if (idx >= N/2 && idx < nxlength-N/2)  	  	  	  d_norder[id] = N;
		if (idx >= nxlength-N/2 && idx < nxlength-1)       	  d_norder[id] = 2*(nxlength-idx-1);
	}
	if (idz >= nz-N/2 && idz < nz-1)
	{
		if (idx >= 1 && idx < N/2)                  	  	  d_norder[id] = 2*(idx<(nz-idz-1)?idx:(nz-idz-1));
		if (idx >= N/2 && idx < nxlength-N/2 )       	  	  d_norder[id] = 2*(nz-idz-1);
		if (idx >= nxlength-N/2 && idx < nxlength-1)       	  d_norder[id] = 2*((nz-idz-1)<(nxlength-idx-1)?(nz-idz-1):(nxlength-idx-1));
	}
}
//__global__ void cuda_init_norderx(int *d_norderx, int nxlength, int N)
__global__ void cuda_norderx(int *d_norderx, int nxlength)
// < set Vx wavefield difference coefficients >
// grid:block ((nxlength+511)/512,512)    call para:nxlength -> nxlength-1 
{
	int idx = threadIdx.x + blockDim.x * blockIdx.x;
	if (idx >=0 && idx < N/2)                       		d_norderx[idx] = 2*(idx+1);
	if (idx >= N/2 && idx < nxlength-N/2)    				d_norderx[idx] = N;
	if (idx >= nxlength-N/2 && idx < nxlength)         		d_norderx[idx] = 2*(nxlength-idx);
}
//__global__ void cuda_init_norderz(int *d_norderz, int nz, int N)
__global__ void cuda_norderz(int *d_norderz, int nz)
// < set Vz wavefield difference coefficients >
// grid:block ((nz+511)/512,512)    para:nz -> nz-1 
{
	int idx = threadIdx.x + blockDim.x * blockIdx.x;
	if (idx >=0 && idx < N/2)                       		d_norderz[idx] = 2*(idx+1);
	if (idx >= N/2 && idx < nz-N/2)    					d_norderz[idx] = N;
	if (idx >= nz-N/2 && idx < nz)         				d_norderz[idx] = 2*(nz-idx);
}
__global__ void cuda_backward_vx(float *d_p, float *d_vx, float *d_rho, float *d_diffcoef, float _dtx, int npml, int nxlength, int nz)
// < Vx: inner update>   iz = npml: nz-npml-1  ix = npml: nxlength - npml -2
// grid:block <(nz+blocksizez-1)/blocksizez,(nxlength-2*npml-1+blocksizex-1)/blocksizex),(blocksizez,blocksizex)>
{
	int idz = threadIdx.x + blockDim.x * blockIdx.x + npml;
	int idx = threadIdx.y + blockDim.y * blockIdx.y + npml;
	int idzl = threadIdx.x;	
	int idxl = threadIdx.y;

	int radius = N/2;
	int tidx  =  idxl + radius - 1;

      float diffx;
      int i;

	__shared__ float pp[Block_Sizez][Block_Sizex + N - 1];

	if (idz < nz-npml)
	{
		if (idxl <  radius - 1)                 pp[idzl][idxl]         = d_p[idz*nxlength + idx - radius + 1];
		if (idxl >= blockDim.y - radius)        pp[idzl][idxl + N - 1] = d_p[idz*nxlength + idx + radius];


		pp[idzl][tidx] = d_p[idz*nxlength + idx];
		__syncthreads();

		if (idx < nxlength-npml-1)
		{
			diffx = 0.0;
			for (i = 1; i <= radius; i++)      
				diffx += d_diffcoef[(radius-1)*radius + i - 1]*(pp[idzl][tidx + i] - pp[idzl][tidx - i + 1]);
			diffx = diffx*_dtx/d_rho[idz*nxlength + idx];
			d_vx[idz*(nxlength - 1) + idx] += diffx; 
		}
	}
}
__global__ void cuda_backward_vz(float *d_p, float *d_vz, float *d_rho, float *d_diffcoef, float _dtz, int npml, int nxlength, int nz)
// < Vz: inner update> iz = npml: nz-npml-2  ix = npml: nxlength-npml
// grid:block <(nz+blocksizez-1)/blocksizez,(nxlength-2*npml-1+blocksizex-1)/blocksizex),(blocksizez,blocksizex)>
{
	int idz = threadIdx.x + blockDim.x * blockIdx.x + npml;
	int idx = threadIdx.y + blockDim.y * blockIdx.y + npml;
	int idzl = threadIdx.x;	
	int idxl = threadIdx.y;

	int radius = N/2;
	int tidz  =  idzl + radius - 1;

      float diffz;
      int i;

	__shared__ float pp[Block_Sizez + N - 1][Block_Sizex];

	if (idx < nxlength-npml)
	{
		if (idzl < radius - 1)                  pp[idzl][idxl]         = d_p[(idz - radius + 1)*nxlength + idx];
		if (idzl >= blockDim.x - radius)        pp[idzl + N - 1][idxl] = d_p[(idz + radius)*nxlength + idx];

		pp[tidz][idxl] = d_p[idz*nxlength + idx];
		__syncthreads();

		if (idz < nz-npml-1)
		{
			diffz = 0.0;
			for (i = 1; i <= radius; i++)      
				diffz += d_diffcoef[(radius-1)*radius + i - 1]*(pp[tidz + i][idxl] - pp[tidz - i + 1][idxl]);
			diffz = diffz*_dtz/d_rho[idz*nxlength + idx];
			d_vz[idz*nxlength + idx] += diffz; 
		}
	}
}
__global__ void cuda_backward_p(float *d_p, float *d_vx, float *d_vz, float *d_rho, float *d_vp, float *d_diffcoef, float _dtx, float _dtz, int npml, int nxlength, int nz)
// < P: inner update>
// grid:block <(nz-2*npml+blocksizez-1)/blocksizez,(nxlength-2*npml+blocksizex-1)/blocksizex),(blocksizez,blocksizex)>
{
	int idz = threadIdx.x + blockDim.x * blockIdx.x + npml;
	int idx = threadIdx.y + blockDim.y * blockIdx.y + npml;
	int idzl = threadIdx.x;
	int idxl = threadIdx.y;
      int radius = N/2;
	int tidz  =  idzl + radius;
	int tidx  =  idxl + radius;

	__shared__ float vvx[Block_Sizez][Block_Sizex + N - 1];
	__shared__ float vvz[Block_Sizez + N - 1][Block_Sizex];

       float diffx,diffz;
       int i;	

	if (idxl < radius)                     vvx[idzl][idxl]         = d_vx[idz*(nxlength - 1) + idx - radius];
	if (idxl >= blockDim.y - radius + 1)   vvx[idzl][idxl + N - 1] = d_vx[idz*(nxlength - 1) + idx + radius - 1];
	if (idzl < radius)                     vvz[idzl][idxl]         = d_vz[(idz - radius)*nxlength + idx];
	if (idzl >= blockDim.x - radius + 1)   vvz[idzl + N - 1][idxl] = d_vz[(idz + radius - 1)*nxlength + idx];

	vvx[idzl][tidx] = d_vx[idz*(nxlength - 1) + idx];
	vvz[tidz][idxl] = d_vz[idz*nxlength + idx];
	__syncthreads();

	if (idz < nz-npml && idx < nxlength - npml)
	{
            diffx = 0.0;
            diffz = 0.0;
		for (i = 1; i <= radius; i++){      
		    diffx += d_diffcoef[(radius-1)*radius + i - 1]*(vvx[idzl][tidx + i - 1] - vvx[idzl][tidx - i]);
                diffz += d_diffcoef[(radius-1)*radius + i - 1]*(vvz[tidz + i - 1][idxl] - vvz[tidz - i][idxl]); }

		diffx = diffx*_dtx*d_rho[idz*nxlength + idx]*d_vp[idz*nxlength + idx];
		diffz = diffz*_dtz*d_rho[idz*nxlength + idx]*d_vp[idz*nxlength + idx];

		d_p[idz*nxlength + idx] += diffx + diffz; 
	}
}
__global__ void cuda_forward_vx(float *d_p, float *d_vx, float *d_rho, float *d_diffcoef, float _dtx, int npml, int nxlength, int nz)
// < Vx: inner update>   iz = 0: nz-1  ix = npml: nxlength - npml -2
// grid:block <(nz+blocksizez-1)/blocksizez,(nxlength-2*npml-1+blocksizex-1)/blocksizex),(blocksizez,blocksizex)>
{
	int idz = threadIdx.x + blockDim.x * blockIdx.x;
	int idx = threadIdx.y + blockDim.y * blockIdx.y + npml;
	int idzl = threadIdx.x;	
	int idxl = threadIdx.y;

	int radius = N/2;
	int tidx  =  idxl + radius - 1;

      float diffx;
      int i;

	__shared__ float pp[Block_Sizez][Block_Sizex + N - 1];

	if (idz < nz)
	{
		if (idxl <  radius - 1)                 pp[idzl][idxl]         = d_p[idz*nxlength + idx - radius + 1];
		if (idxl >= blockDim.y - radius)        pp[idzl][idxl + N - 1] = d_p[idz*nxlength + idx + radius];


		pp[idzl][tidx] = d_p[idz*nxlength + idx];
		__syncthreads();

		if (idx < nxlength-npml-1)
		{
			diffx = 0.0;
			for (i = 1; i <= radius; i++)      
				diffx += d_diffcoef[(radius-1)*radius + i - 1]*(pp[idzl][tidx + i] - pp[idzl][tidx - i + 1]);
			diffx = diffx*_dtx/d_rho[idz*nxlength + idx];
			d_vx[idz*(nxlength - 1) + idx] -= diffx; 
		}
	}
}
__global__ void cuda_forward_vz(float *d_p, float *d_vz, float *d_rho, float *d_diffcoef, float _dtz, int npml, int nxlength, int nz)
// < Vz: inner update> iz = npml: nz-npml-2  ix = 0: nxlength
// grid:block <(nz+blocksizez-1)/blocksizez,(nxlength-2*npml-1+blocksizex-1)/blocksizex),(blocksizez,blocksizex)>
{
	int idz = threadIdx.x + blockDim.x * blockIdx.x + npml;
	int idx = threadIdx.y + blockDim.y * blockIdx.y;
	int idzl = threadIdx.x;	
	int idxl = threadIdx.y;

	int radius = N/2;
	int tidz  =  idzl + radius - 1;

      float diffz;
      int i;

	__shared__ float pp[Block_Sizez + N - 1][Block_Sizex];

	if (idx < nxlength)
	{
		if (idzl < radius - 1)                  pp[idzl][idxl]         = d_p[(idz - radius + 1)*nxlength + idx];
		if (idzl >= blockDim.x - radius)        pp[idzl + N - 1][idxl] = d_p[(idz + radius)*nxlength + idx];

		pp[tidz][idxl] = d_p[idz*nxlength + idx];
		__syncthreads();

		if (idz < nz-npml-1)
		{
			diffz = 0.0;
			for (i = 1; i <= radius; i++)      
				diffz += d_diffcoef[(radius-1)*radius + i - 1]*(pp[tidz + i][idxl] - pp[tidz - i + 1][idxl]);
			diffz = diffz*_dtz/d_rho[idz*nxlength + idx];
			d_vz[idz*nxlength + idx] -= diffz; 
		}
	}
}
__global__ void cuda_forward_p(float *d_p, float *d_vx, float *d_vz, float *d_rho, float *d_vp, float *d_diffcoef, float _dtx, float _dtz, int npml, int nxlength, int nz)
// < P: inner update>
// grid:block <(nz-2*npml+blocksizez-1)/blocksizez,(nxlength-2*npml+blocksizex-1)/blocksizex),(blocksizez,blocksizex)>
{
	int idz = threadIdx.x + blockDim.x * blockIdx.x + npml;
	int idx = threadIdx.y + blockDim.y * blockIdx.y + npml;
	int idzl = threadIdx.x;
	int idxl = threadIdx.y;
      int radius = N/2;
	int tidz  =  idzl + radius;
	int tidx  =  idxl + radius;

	__shared__ float vvx[Block_Sizez][Block_Sizex + N - 1];
	__shared__ float vvz[Block_Sizez + N - 1][Block_Sizex];

       float diffx,diffz;
       int i;	

	if (idxl < radius)                     vvx[idzl][idxl]         = d_vx[idz*(nxlength - 1) + idx - radius];
	if (idxl >= blockDim.y - radius + 1)   vvx[idzl][idxl + N - 1] = d_vx[idz*(nxlength - 1) + idx + radius - 1];
	if (idzl < radius)                     vvz[idzl][idxl]         = d_vz[(idz - radius)*nxlength + idx];
	if (idzl >= blockDim.x - radius + 1)   vvz[idzl + N - 1][idxl] = d_vz[(idz + radius - 1)*nxlength + idx];

	vvx[idzl][tidx] = d_vx[idz*(nxlength - 1) + idx];
	vvz[tidz][idxl] = d_vz[idz*nxlength + idx];
	__syncthreads();

	if (idz < nz-npml && idx < nxlength - npml)
	{
            diffx = 0.0;
            diffz = 0.0;
		for (i = 1; i <= radius; i++){      
		    diffx += d_diffcoef[(radius-1)*radius + i - 1]*(vvx[idzl][tidx + i - 1] - vvx[idzl][tidx - i]);
                diffz += d_diffcoef[(radius-1)*radius + i - 1]*(vvz[tidz + i - 1][idxl] - vvz[tidz - i][idxl]); }

		diffx = diffx*_dtx*d_rho[idz*nxlength + idx]*d_vp[idz*nxlength + idx];
		diffz = diffz*_dtz*d_rho[idz*nxlength + idx]*d_vp[idz*nxlength + idx];

		d_p[idz*nxlength + idx] -= diffx + diffz; 
	}
}
__global__ void cuda_pml_vxlr(float *d_p, float *d_vx, float *d_rho, float *d_diffcoef, float *d_ddxVx, float _dtx, float dt, int npml, int nxlength, int nz, int *d_norderx)
// < Vx: pml left and right >
// grid:block <(nz+blocksizez-1)/blocksizez,2),(blocksizez,blocksizex)>
{
	int idz = threadIdx.x + blockDim.x * blockIdx.x;
	int idx = threadIdx.y + blockDim.y * blockIdx.y;
	int idzl = threadIdx.x;
	int idxl = threadIdx.y;
	int tidx;
	int radius = N/2;
	float diffx,c1,c2;

      int radiuslx1, radiuslx2, radiuslx;
      int radiusx;

	__shared__ float pp[Block_Sizez][Block_Sizex + N/2];

      int i;
	
	if (idz < nz)
	{

		if (blockIdx.y == 0)
		{
			//left
			idx = idx;

			radiuslx1 = 1;                        // left  order
			radiuslx2 = radius;                   // right order 
			radiuslx  = radiuslx1 + radiuslx2;    // 
			

			tidx = idxl + radiuslx1 - 1;

			if (idxl <  radiuslx1 - 1)                pp[idzl][idxl]                = d_p[idz*nxlength + idx - radiuslx1 + 1];
			if (idxl >= blockDim.y - radiuslx2)       pp[idzl][idxl + radiuslx - 1] = d_p[idz*nxlength + idx + radiuslx2];


			pp[idzl][tidx] = d_p[idz*nxlength + idx];
			__syncthreads();

			if (idx < npml)
			{
				diffx = 0.0;
				radiusx = d_norderx[idx]/2;           // local order - x

				for (i = 1; i <= radiusx; i++)      
					diffx += d_diffcoef[(radiusx-1)*radius + i - 1]*(pp[idzl][tidx + i] - pp[idzl][tidx - i + 1]);

				diffx = diffx*_dtx/d_rho[idz*nxlength + idx];
				c1 = 1.0 - 0.5*dt*d_ddxVx[idx];
				c2 = 2.0 - c1;
				d_vx[idz*(nxlength - 1) + idx] = (c1*d_vx[idz*(nxlength - 1) + idx] - diffx)/c2; 
			}
		}
		if (blockIdx.y == 1)
		{
			//right
			idx = nxlength -2*npml + idx - 1;

			radiuslx1 = radius;                        // left  order
			radiuslx2 = 1;                             // right order 
			radiuslx  = radiuslx1 + radiuslx2;         // 
			

			tidx = idxl + radiuslx1 - 1;

			if (idxl <  radiuslx1 - 1)                pp[idzl][idxl]                = d_p[idz*nxlength + idx - radiuslx1 + 1];
			if (idxl >= blockDim.y - radiuslx2)       pp[idzl][idxl + radiuslx - 1] = d_p[idz*nxlength + idx + radiuslx2];

			pp[idzl][tidx] = d_p[idz*nxlength + idx];
			__syncthreads();

			if (idx >= nxlength - npml - 1 && idx < nxlength - 1)
			{
				diffx = 0.0;
				radiusx = d_norderx[idx]/2;           // local order - x
				//vx
				for (i = 1; i <= radiusx; i++)      
					diffx += d_diffcoef[(radiusx-1)*radius + i - 1]*(pp[idzl][tidx + i] - pp[idzl][tidx - i + 1]);

				diffx = diffx*_dtx/d_rho[idz*nxlength + idx];
				c1 = 1.0 - 0.5*dt*d_ddxVx[idx];
				c2 = 2.0 - c1;
				d_vx[idz*(nxlength - 1) + idx] = (c1*d_vx[idz*(nxlength - 1) + idx] - diffx)/c2; 
			}
		}
	}
}
__global__ void cuda_pml_vztb(float *d_p, float *d_vz, float *d_rho, float *d_diffcoef, float *d_ddzVz, float _dtz, float dt, int npml, int nxlength, int nz, int *d_norderz)
// <Vx Vz: pml top and bottom >
// grid:block <2,(nxlength+blocksizez-1)/blocksizez),(blocksizez,blocksizex)>
{
	int idz = threadIdx.x + blockDim.x * blockIdx.x;
	int idx = threadIdx.y + blockDim.y * blockIdx.y;
	int idzl = threadIdx.x;
	int idxl = threadIdx.y;
	int tidz;
	int radius = N/2;
	float diffz,c1,c2;

      int radiuslz1, radiuslz2, radiuslz;
      int radiusz;

	__shared__ float pp[Block_Sizez + N/2][Block_Sizex];

      int i;

	
	if (idx < nxlength)
	{
		if (blockIdx.x == 0)
		{
			//top
			idz = idz;
			
			radiuslz1 = 1;
			radiuslz2 = radius;
			radiuslz  = radiuslz1 + radiuslz2;

			tidz = idzl + radiuslz1 - 1;

			if (idzl < radiuslz1 - 1)                 pp[idzl][idxl]                = d_p[(idz - radiuslz1 + 1)*nxlength + idx];
			if (idzl >= blockDim.x - radiuslz2)       pp[idzl + radiuslz - 1][idxl] = d_p[(idz + radiuslz2)*nxlength + idx];

			pp[tidz][idxl] = d_p[idz*nxlength + idx];
			__syncthreads();

			if (idz < npml)
			{
				diffz = 0.0;
				radiusz = d_norderz[idz]/2;		  // local order - z
				for (i = 1; i <= radiusz; i++)      
					diffz += d_diffcoef[(radiusz-1)*radius + i - 1]*(pp[tidz + i][idxl] - pp[tidz - i + 1][idxl]);

				diffz = diffz*_dtz/d_rho[idz*nxlength + idx];
				c1 = 1.0 - 0.5*dt*d_ddzVz[idz];
				c2 = 2.0 - c1;
				d_vz[idz*nxlength + idx] = (c1*d_vz[idz*nxlength + idx] - diffz)/c2; 
			}
		}
		if (blockIdx.x == 1)
		{
			//bottom
			idz = nz -2*npml + idz - 1;

			radiuslz1 = radius;
			radiuslz2 = 1;
			radiuslz  = radiuslz1 + radiuslz2;

			tidz = idzl + radiuslz1 - 1;

			if (idzl < radiuslz1 - 1)                 pp[idzl][idxl]                = d_p[(idz - radiuslz1 + 1)*nxlength + idx];
			if (idzl >= blockDim.x - radiuslz2)       pp[idzl + radiuslz - 1][idxl] = d_p[(idz + radiuslz2)*nxlength + idx];

			pp[tidz][idxl] = d_p[idz*nxlength + idx];
			__syncthreads();

			if (idz >= nz - npml - 1 && idz < nz - 1)
			{
				diffz = 0.0;
				radiusz = d_norderz[idz]/2;		  // local order - z
				for (i = 1; i <= radiusz; i++)      
					diffz += d_diffcoef[(radiusz-1)*radius + i - 1]*(pp[tidz + i][idxl] - pp[tidz - i + 1][idxl]);
				diffz = diffz*_dtz/d_rho[idz*nxlength + idx];
				c1 = 1.0 - 0.5*dt*d_ddzVz[idz];
				c2 = 2.0 - c1;
				d_vz[idz*nxlength + idx] = (c1*d_vz[idz*nxlength + idx] - diffz)/c2; 
			}
		}
	}
}
__global__ void cuda_pml_plr(float *d_p, float *d_vx, float *d_vz, float *d_pl1, float *d_pl2, float *d_pr1, float *d_pr2, float *d_rho, float *d_vp, float *d_diffcoef, float *d_ddx, float *d_ddz, float _dtx, float _dtz, float dt, int npml, int nxlength, int nz, int *d_norder)
// < P: pml left and right >
// grid:block <(nz+blocksizez-1)/blocksizez,2),(blocksizez,blocksizex)>
// blockIdx.y == 0 : pl1 dvx/dx pl2 dvz/dz
// blockIdx.y == 1 : pr1 dvx/dx pr2 dvz/dz
{
	int idz = threadIdx.x + blockDim.x * blockIdx.x + npml;
	int idx = threadIdx.y + blockDim.y * blockIdx.y;
	int idzl = threadIdx.x;
	int idxl = threadIdx.y;
	int tidx,tidz;
	int radius = N/2;
	float diffx,diffz,c1x,c2x,c1z,c2z;

	__shared__ float vvx[Block_Sizez][Block_Sizex + N/2];  // left: + radius   right: + radius - 1
	__shared__ float vvz[Block_Sizez + N - 1][Block_Sizex];  

      int radiuslx1, radiuslx2, radiuslx,
          radiuslz1, radiuslz2, radiuslz;
      int radiusl;
      int i;


	if (blockIdx.y == 0)
	{
		
		//left			
		radiuslx1 = 0;
		radiuslx2 = radius;
		radiuslx  = radiuslx1 + radiuslx2;

		radiuslz1 = radius;
		radiuslz2 = radius;
		radiuslz  = radiuslz1 + radiuslz2;

		tidx = idxl + radiuslx1;
		tidz = idzl + radiuslz1;

		if (idxl < radiuslx1)                     vvx[idzl][idxl]                = d_vx[idz*(nxlength - 1) + idx - radiuslx1];
		if (idxl >= blockDim.y - radiuslx2 + 1)   vvx[idzl][idxl + radiuslx - 1] = d_vx[idz*(nxlength - 1) + idx + radiuslx2 - 1];
		vvx[idzl][tidx] = d_vx[idz*(nxlength - 1) + idx];

		if (idzl < radiuslz1)                     vvz[idzl][idxl]                = d_vz[(idz - radiuslz1)*nxlength + idx];
		if (idzl >= blockDim.x - radiuslz2 + 1)   vvz[idzl + radiuslz - 1][idxl] = d_vz[(idz + radiuslz2 -1)*nxlength + idx];
		vvz[tidz][idxl] = d_vz[idz*nxlength + idx];		
	
		__syncthreads();

		if (idx > 0 && idz < nz - npml)
		{
			radiusl = d_norder[idz*nxlength + idx]/2;
			diffx = 0.0;
			diffz = 0.0;
			for (i = 1; i <= radiusl; i++){      diffx += d_diffcoef[(radiusl-1)*radius + i - 1]*(vvx[idzl][tidx + i - 1] - vvx[idzl][tidx - i]);
								   	 diffz += d_diffcoef[(radiusl-1)*radius + i - 1]*(vvz[tidz + i - 1][idxl] - vvz[tidz - i][idxl]); }

			diffx = diffx*_dtx*d_rho[idz*nxlength + idx]*d_vp[idz*nxlength + idx];
			diffz = diffz*_dtz*d_rho[idz*nxlength + idx]*d_vp[idz*nxlength + idx];

		
			c1x = 1.0 - 0.5*dt*d_ddx[idx];
			c2x = 2.0 - c1x;

			c1z = 1.0;
			c2z = 1.0;

			d_pl1[idz*npml + idx] = (c1x*d_pl1[idz*npml + idx] - diffx)/c2x; 
			d_pl2[idz*npml + idx] = (c1z*d_pl2[idz*npml + idx] - diffz)/c2z; 

			d_p[idz*nxlength + idx] = d_pl1[idz*npml + idx] + d_pl2[idz*npml + idx];
		}
	}
	if (blockIdx.y == 1)
	{
		//right
		idx = nxlength -2*npml + idx;
			
		radiuslx1 = radius;
		radiuslx2 = 0;
		radiuslx  = radiuslx1 + radiuslx2;

		radiuslz1 = radius;
		radiuslz2 = radius;
		radiuslz  = radiuslz1 + radiuslz2;

		tidx = idxl + radiuslx1;
		tidz = idzl + radiuslz1;

		if (idxl < radiuslx1)                     vvx[idzl][idxl]                = d_vx[idz*(nxlength - 1) + idx - radiuslx1];
		if (idxl >= blockDim.y - radiuslx2 + 1)   vvx[idzl][idxl + radiuslx - 1] = d_vx[idz*(nxlength - 1) + idx + radiuslx2 - 1];
		vvx[idzl][tidx] = d_vx[idz*(nxlength - 1) + idx];

		if (idzl < radiuslz1)                     vvz[idzl][idxl]                = d_vz[(idz - radiuslz1)*nxlength + idx];
		if (idzl >= blockDim.x - radiuslz2 + 1)   vvz[idzl + radiuslz - 1][idxl] = d_vz[(idz + radiuslz2 -1)*nxlength + idx];
		vvz[tidz][idxl] = d_vz[idz*nxlength + idx];		
	
		__syncthreads();

		if (idx < nxlength-1 && idz < nz - npml)
		{
			radiusl = d_norder[idz*nxlength + idx]/2;
			diffx = 0.0;
			diffz = 0.0;
			for (i = 1; i <= radiusl; i++){      diffx += d_diffcoef[(radiusl-1)*radius + i - 1]*(vvx[idzl][tidx + i - 1] - vvx[idzl][tidx - i]);
								   	 diffz += d_diffcoef[(radiusl-1)*radius + i - 1]*(vvz[tidz + i - 1][idxl] - vvz[tidz - i][idxl]); }

			diffx = diffx*_dtx*d_rho[idz*nxlength + idx]*d_vp[idz*nxlength + idx];
			diffz = diffz*_dtz*d_rho[idz*nxlength + idx]*d_vp[idz*nxlength + idx];

		
			c1x = 1.0 - 0.5*dt*d_ddx[idx];
			c2x = 2.0 - c1x;

			c1z = 1.0;
			c2z = 1.0;

//			d_pr1[idz*npml + idx - nxlength +2*npml] = (c1x*d_pr1[idz*npml + idx - nxlength +2*npml] - diffx)/c2x; 
//			d_pr2[idz*npml + idx - nxlength +2*npml] = (c1z*d_pr2[idz*npml + idx - nxlength +2*npml] - diffz)/c2z; 

//			d_p[idz*nxlength + idx] = d_pr1[idz*npml+ idx - nxlength +2*npml] + d_pr2[idz*npml + idx - nxlength +2*npml];

			d_pr1[idz*npml + idx - nxlength + npml] = (c1x*d_pr1[idz*npml + idx - nxlength + npml] - diffx)/c2x; 
			d_pr2[idz*npml + idx - nxlength + npml] = (c1z*d_pr2[idz*npml + idx - nxlength + npml] - diffz)/c2z; 

			d_p[idz*nxlength + idx] = d_pr1[idz*npml+ idx - nxlength + npml] + d_pr2[idz*npml + idx - nxlength + npml];
		}
	}
}
__global__ void cuda_pml_ptb(float *d_p, float *d_vx, float *d_vz, float *d_pt1, float *d_pt2, float *d_pb1, float *d_pb2, float *d_rho, float *d_vp, float *d_diffcoef, float *d_ddz, float _dtx, float _dtz, float dt, int npml, int nxlength, int nz, int *d_norder)
// < P: pml left and right >
// grid:block <(2,nxlength-2*npml+blocksizez-1)/blocksizez),(blocksizez,blocksizex)>
// blockIdx.y == 0 : pt1 dvx/dx pt2 dvz/dz
// blockIdx.y == 1 : pb1 dvx/dx pb2 dvz/dz
{
	int idz = threadIdx.x + blockDim.x * blockIdx.x;
	int idx = threadIdx.y + blockDim.y * blockIdx.y + npml;
	int idzl = threadIdx.x;
	int idxl = threadIdx.y;
	int tidx,tidz;
	int radius = N/2;
	float diffx,diffz,c1x,c2x,c1z,c2z;

	__shared__ float vvx[Block_Sizez][Block_Sizex + N - 1];  // left: + radius - 1  right: + radius
	__shared__ float vvz[Block_Sizez + N/2][Block_Sizex];  

      int radiuslx1, radiuslx2, radiuslx,
          radiuslz1, radiuslz2, radiuslz;
      int radiusl;
      int i;

	if (blockIdx.x == 0)
	{
		//top
		radiuslx1 = radius;
		radiuslx2 = radius;
		radiuslx  = radiuslx1 + radiuslx2;

		radiuslz1 = 0;
		radiuslz2 = radius;
		radiuslz  = radiuslz1 + radiuslz2;

		tidx = idxl + radiuslx1;
		tidz = idzl + radiuslz1;

		if (idxl < radiuslx1)                     vvx[idzl][idxl]                = d_vx[idz*(nxlength - 1) + idx - radiuslx1];
		if (idxl >= blockDim.y - radiuslx2 + 1)   vvx[idzl][idxl + radiuslx - 1] = d_vx[idz*(nxlength - 1) + idx + radiuslx2 - 1];
		vvx[idzl][tidx] = d_vx[idz*(nxlength - 1) + idx];

		if (idzl < radiuslz1)                     vvz[idzl][idxl]                = d_vz[(idz - radiuslz1)*nxlength + idx];
		if (idzl >= blockDim.x - radiuslz2 + 1)   vvz[idzl + radiuslz - 1][idxl] = d_vz[(idz + radiuslz2 -1)*nxlength + idx];
		vvz[tidz][idxl] = d_vz[idz*nxlength + idx];		
	
		__syncthreads();

		if (idz > 0 && idx < nxlength - npml)
		{
			radiusl = d_norder[idz*nxlength + idx]/2;
			diffx = 0.0;
			diffz = 0.0;
			for (i = 1; i <= radiusl; i++){      diffx += d_diffcoef[(radiusl-1)*radius + i - 1]*(vvx[idzl][tidx + i - 1] - vvx[idzl][tidx - i]);
								   	 diffz += d_diffcoef[(radiusl-1)*radius + i - 1]*(vvz[tidz + i - 1][idxl] - vvz[tidz - i][idxl]); }

			diffx = diffx*_dtx*d_rho[idz*nxlength + idx]*d_vp[idz*nxlength + idx];
			diffz = diffz*_dtz*d_rho[idz*nxlength + idx]*d_vp[idz*nxlength + idx];

		
			c1x = 1.0;
			c2x = 1.0;

			c1z = 1.0 - 0.5*dt*d_ddz[idz];
			c2z = 2.0 - c1z;

			d_pt1[idz*(nxlength - 2*npml) + idx - npml] = (c1x*d_pt1[idz*(nxlength - 2*npml) + idx - npml] - diffx)/c2x; 
			d_pt2[idz*(nxlength - 2*npml) + idx - npml] = (c1z*d_pt2[idz*(nxlength - 2*npml) + idx - npml] - diffz)/c2z; 

			d_p[idz*nxlength + idx] = d_pt1[idz*(nxlength - 2*npml) + idx - npml] + d_pt2[idz*(nxlength - 2*npml) + idx - npml];
		}
	}
	if (blockIdx.x == 1)
	{
		//bottom
		idz = nz - 2*npml + idz;

		radiuslx1 = radius;
		radiuslx2 = radius;
		radiuslx  = radiuslx1 + radiuslx2;

		radiuslz1 = radius;
		radiuslz2 = 0;
		radiuslz  = radiuslz1 + radiuslz2;

		tidx = idxl + radiuslx1;
		tidz = idzl + radiuslz1;

		if (idxl < radiuslx1)                     vvx[idzl][idxl]                = d_vx[idz*(nxlength - 1) + idx - radiuslx1];
		if (idxl >= blockDim.y - radiuslx2 + 1)   vvx[idzl][idxl + radiuslx - 1] = d_vx[idz*(nxlength - 1) + idx + radiuslx2 - 1];
		vvx[idzl][tidx] = d_vx[idz*(nxlength - 1) + idx];

		if (idzl < radiuslz1)                     vvz[idzl][idxl]                = d_vz[(idz - radiuslz1)*nxlength + idx];
		if (idzl >= blockDim.x - radiuslz2 + 1)   vvz[idzl + radiuslz - 1][idxl] = d_vz[(idz + radiuslz2 -1)*nxlength + idx];
		vvz[tidz][idxl] = d_vz[idz*nxlength + idx];		
	
		__syncthreads();

		if (idz < nz-1 && idx < nxlength - npml)
		{
			radiusl = d_norder[idz*nxlength + idx]/2;
			diffx = 0.0;
			diffz = 0.0;
			for (i = 1; i <= radiusl; i++){      diffx += d_diffcoef[(radiusl-1)*radius + i - 1]*(vvx[idzl][tidx + i - 1] - vvx[idzl][tidx - i]);
								   	 diffz += d_diffcoef[(radiusl-1)*radius + i - 1]*(vvz[tidz + i - 1][idxl] - vvz[tidz - i][idxl]); }

			diffx = diffx*_dtx*d_rho[idz*nxlength + idx]*d_vp[idz*nxlength + idx];
			diffz = diffz*_dtz*d_rho[idz*nxlength + idx]*d_vp[idz*nxlength + idx];

		
			c1x = 1.0;
			c2x = 1.0;

			c1z = 1.0 - 0.5*dt*d_ddz[idz];
			c2z = 2.0 - c1z;


			d_pb1[(idz - nz + npml)*(nxlength - 2*npml) + idx - npml] = (c1x*d_pb1[(idz - nz + npml)*(nxlength - 2*npml) + idx - npml] - diffx)/c2x; 
			d_pb2[(idz - nz + npml)*(nxlength - 2*npml) + idx - npml] = (c1z*d_pb2[(idz - nz + npml)*(nxlength - 2*npml) + idx - npml] - diffz)/c2z;

			d_p[idz*nxlength + idx] = d_pb1[(idz - nz + npml)*(nxlength - 2*npml) + idx - npml] + d_pb2[(idz - nz + npml)*(nxlength - 2*npml) + idx - npml];
		}
	}
}
__global__ void cuda_pml_pconner(float *d_p, float *d_vx, float *d_vz, float *d_pl1, float *d_pl2, float *d_pr1, float *d_pr2, float *d_rho, float *d_vp, float *d_diffcoef, float *d_ddx, float *d_ddz, float _dtx, float _dtz, float dt, int npml, int nxlength, int nz, int *d_norder)
// < P: pml left and right >
// grid:block <(nz+blocksizez-1)/blocksizez,2),(blocksizez,blocksizex)>
// blockIdx.y == 0 : pl1 dvx/dx pl2 dvz/dz
// blockIdx.y == 1 : pr1 dvx/dx pr2 dvz/dz
{
	int idz = threadIdx.x + blockDim.x * blockIdx.x;
	int idx = threadIdx.y + blockDim.y * blockIdx.y;
	int idzl = threadIdx.x;
	int idxl = threadIdx.y;
	int tidx,tidz;
	int radius = N/2;
	float diffx,diffz,c1x,c2x,c1z,c2z;

	__shared__ float vvx[Block_Sizez][Block_Sizex + N/2];  // left: + radius   right: + radius - 1
	__shared__ float vvz[Block_Sizez + N/2][Block_Sizex];  

      int radiuslx1, radiuslx2, radiuslx,
          radiuslz1, radiuslz2, radiuslz;
      int radiusl;
      int i;
	if (blockIdx.y == 0 && blockIdx.x == 0)  // left -up conner
	{
		
		//left			
		radiuslx1 = 0;
		radiuslx2 = radius;
		radiuslx  = radiuslx1 + radiuslx2;

		radiuslz1 = 0;
		radiuslz2 = radius;
		radiuslz  = radiuslz1 + radiuslz2;

		tidx = idxl + radiuslx1;
		tidz = idzl + radiuslz1;

		if (idxl < radiuslx1)                     vvx[idzl][idxl]                = d_vx[idz*(nxlength - 1) + idx - radiuslx1];
		if (idxl >= blockDim.y - radiuslx2 + 1)   vvx[idzl][idxl + radiuslx - 1] = d_vx[idz*(nxlength - 1) + idx + radiuslx2 - 1];
		vvx[idzl][tidx] = d_vx[idz*(nxlength - 1) + idx];

		if (idzl < radiuslz1)                     vvz[idzl][idxl]                = d_vz[(idz - radiuslz1)*nxlength + idx];
		if (idzl >= blockDim.x - radiuslz2 + 1)   vvz[idzl + radiuslz - 1][idxl] = d_vz[(idz + radiuslz2 -1)*nxlength + idx];
		vvz[tidz][idxl] = d_vz[idz*nxlength + idx];		
	
		__syncthreads();

		if (idx > 0 && idz > 0)
		{
			radiusl = d_norder[idz*nxlength + idx]/2;
			diffx = 0.0;
			diffz = 0.0;
			for (i = 1; i <= radiusl; i++){      diffx += d_diffcoef[(radiusl-1)*radius + i - 1]*(vvx[idzl][tidx + i - 1] - vvx[idzl][tidx - i]);
								   	 diffz += d_diffcoef[(radiusl-1)*radius + i - 1]*(vvz[tidz + i - 1][idxl] - vvz[tidz - i][idxl]); }

			diffx = diffx*_dtx*d_rho[idz*nxlength + idx]*d_vp[idz*nxlength + idx];
			diffz = diffz*_dtz*d_rho[idz*nxlength + idx]*d_vp[idz*nxlength + idx];

		
			c1x = 1.0 - 0.5*dt*d_ddx[idx];
			c2x = 2.0 - c1x;

			c1z = 1.0 - 0.5*dt*d_ddz[idz];
			c2z = 2.0 - c1z;

			d_pl1[idz*npml + idx] = (c1x*d_pl1[idz*npml + idx] - diffx)/c2x; 
			d_pl2[idz*npml + idx] = (c1z*d_pl2[idz*npml + idx] - diffz)/c2z; 

			d_p[idz*nxlength + idx] = d_pl1[idz*npml + idx] + d_pl2[idz*npml + idx];
		}
	}
	if (blockIdx.y == 1 && blockIdx.x == 0)  // right - up conner
	{
		//right
		idx = nxlength -2*npml + idx;
			
		radiuslx1 = radius;
		radiuslx2 = 0;
		radiuslx  = radiuslx1 + radiuslx2;

		radiuslz1 = 0;
		radiuslz2 = radius;
		radiuslz  = radiuslz1 + radiuslz2;

		tidx = idxl + radiuslx1;
		tidz = idzl + radiuslz1;

		if (idxl < radiuslx1)                     vvx[idzl][idxl]                = d_vx[idz*(nxlength - 1) + idx - radiuslx1];
		if (idxl >= blockDim.y - radiuslx2 + 1)   vvx[idzl][idxl + radiuslx - 1] = d_vx[idz*(nxlength - 1) + idx + radiuslx2 - 1];
		vvx[idzl][tidx] = d_vx[idz*(nxlength - 1) + idx];

		if (idzl < radiuslz1)                     vvz[idzl][idxl]                = d_vz[(idz - radiuslz1)*nxlength + idx];
		if (idzl >= blockDim.x - radiuslz2 + 1)   vvz[idzl + radiuslz - 1][idxl] = d_vz[(idz + radiuslz2 -1)*nxlength + idx];
		vvz[tidz][idxl] = d_vz[idz*nxlength + idx];		
	
		__syncthreads();

		if (idx < nxlength-1 && idz > 0)
		{
			radiusl = d_norder[idz*nxlength + idx]/2;
			diffx = 0.0;
			diffz = 0.0;
			for (i = 1; i <= radiusl; i++){      diffx += d_diffcoef[(radiusl-1)*radius + i - 1]*(vvx[idzl][tidx + i - 1] - vvx[idzl][tidx - i]);
								   	 diffz += d_diffcoef[(radiusl-1)*radius + i - 1]*(vvz[tidz + i - 1][idxl] - vvz[tidz - i][idxl]); }

			diffx = diffx*_dtx*d_rho[idz*nxlength + idx]*d_vp[idz*nxlength + idx];
			diffz = diffz*_dtz*d_rho[idz*nxlength + idx]*d_vp[idz*nxlength + idx];

		
			c1x = 1.0 - 0.5*dt*d_ddx[idx];
			c2x = 2.0 - c1x;

			c1z = 1.0 - 0.5*dt*d_ddz[idz];
			c2z = 2.0 - c1z;

			d_pr1[idz*npml + idx - nxlength + npml] = (c1x*d_pr1[idz*npml + idx - nxlength + npml] - diffx)/c2x; 
			d_pr2[idz*npml + idx - nxlength + npml] = (c1z*d_pr2[idz*npml + idx - nxlength + npml] - diffz)/c2z; 

			d_p[idz*nxlength + idx] = d_pr1[idz*npml+ idx - nxlength + npml] + d_pr2[idz*npml + idx - nxlength + npml];
		}
	}
	if (blockIdx.y == 0 && blockIdx.x == 1)  // left -bottom conner
	{
		
		//left
		idz = nz -2*npml + idz;			
		radiuslx1 = 0;
		radiuslx2 = radius;
		radiuslx  = radiuslx1 + radiuslx2;

		radiuslz1 = radius;
		radiuslz2 = 0;
		radiuslz  = radiuslz1 + radiuslz2;

		tidx = idxl + radiuslx1;
		tidz = idzl + radiuslz1;

		if (idxl < radiuslx1)                     vvx[idzl][idxl]                = d_vx[idz*(nxlength - 1) + idx - radiuslx1];
		if (idxl >= blockDim.y - radiuslx2 + 1)   vvx[idzl][idxl + radiuslx - 1] = d_vx[idz*(nxlength - 1) + idx + radiuslx2 - 1];
		vvx[idzl][tidx] = d_vx[idz*(nxlength - 1) + idx];

		if (idzl < radiuslz1)                     vvz[idzl][idxl]                = d_vz[(idz - radiuslz1)*nxlength + idx];
		if (idzl >= blockDim.x - radiuslz2 + 1)   vvz[idzl + radiuslz - 1][idxl] = d_vz[(idz + radiuslz2 -1)*nxlength + idx];
		vvz[tidz][idxl] = d_vz[idz*nxlength + idx];		
	
		__syncthreads();

		if (idx > 0 && idz < nz-1)
		{
			radiusl = d_norder[idz*nxlength + idx]/2;
			diffx = 0.0;
			diffz = 0.0;
			for (i = 1; i <= radiusl; i++){      diffx += d_diffcoef[(radiusl-1)*radius + i - 1]*(vvx[idzl][tidx + i - 1] - vvx[idzl][tidx - i]);
								   	 diffz += d_diffcoef[(radiusl-1)*radius + i - 1]*(vvz[tidz + i - 1][idxl] - vvz[tidz - i][idxl]); }

			diffx = diffx*_dtx*d_rho[idz*nxlength + idx]*d_vp[idz*nxlength + idx];
			diffz = diffz*_dtz*d_rho[idz*nxlength + idx]*d_vp[idz*nxlength + idx];

		
			c1x = 1.0 - 0.5*dt*d_ddx[idx];
			c2x = 2.0 - c1x;

			c1z = 1.0 - 0.5*dt*d_ddz[idz];
			c2z = 2.0 - c1z;

			d_pl1[idz*npml + idx] = (c1x*d_pl1[idz*npml + idx] - diffx)/c2x; 
			d_pl2[idz*npml + idx] = (c1z*d_pl2[idz*npml + idx] - diffz)/c2z; 

			d_p[idz*nxlength + idx] = d_pl1[idz*npml + idx] + d_pl2[idz*npml + idx];
		}
	}
	if (blockIdx.y == 1 && blockIdx.x == 1)  // right - bottom conner
	{
		//right
		idx = nxlength -2*npml + idx;
		idz = nz -2*npml + idz;	

		radiuslx1 = radius;
		radiuslx2 = 0;
		radiuslx  = radiuslx1 + radiuslx2;

		radiuslz1 = radius;
		radiuslz2 = 0;
		radiuslz  = radiuslz1 + radiuslz2;

		tidx = idxl + radiuslx1;
		tidz = idzl + radiuslz1;

		if (idxl < radiuslx1)                     vvx[idzl][idxl]                = d_vx[idz*(nxlength - 1) + idx - radiuslx1];
		if (idxl >= blockDim.y - radiuslx2 + 1)   vvx[idzl][idxl + radiuslx - 1] = d_vx[idz*(nxlength - 1) + idx + radiuslx2 - 1];
		vvx[idzl][tidx] = d_vx[idz*(nxlength - 1) + idx];

		if (idzl < radiuslz1)                     vvz[idzl][idxl]                = d_vz[(idz - radiuslz1)*nxlength + idx];
		if (idzl >= blockDim.x - radiuslz2 + 1)   vvz[idzl + radiuslz - 1][idxl] = d_vz[(idz + radiuslz2 -1)*nxlength + idx];
		vvz[tidz][idxl] = d_vz[idz*nxlength + idx];		
	
		__syncthreads();

		if (idx < nxlength-1 && idz < nz-1)
		{
			radiusl = d_norder[idz*nxlength + idx]/2;
			diffx = 0.0;
			diffz = 0.0;
			for (i = 1; i <= radiusl; i++){      diffx += d_diffcoef[(radiusl-1)*radius + i - 1]*(vvx[idzl][tidx + i - 1] - vvx[idzl][tidx - i]);
								   	 diffz += d_diffcoef[(radiusl-1)*radius + i - 1]*(vvz[tidz + i - 1][idxl] - vvz[tidz - i][idxl]); }

			diffx = diffx*_dtx*d_rho[idz*nxlength + idx]*d_vp[idz*nxlength + idx];
			diffz = diffz*_dtz*d_rho[idz*nxlength + idx]*d_vp[idz*nxlength + idx];

		
			c1x = 1.0 - 0.5*dt*d_ddx[idx];
			c2x = 2.0 - c1x;

			c1z = 1.0 - 0.5*dt*d_ddz[idz];
			c2z = 2.0 - c1z;

			d_pr1[idz*npml + idx - nxlength + npml] = (c1x*d_pr1[idz*npml + idx - nxlength + npml] - diffx)/c2x; 
			d_pr2[idz*npml + idx - nxlength + npml] = (c1z*d_pr2[idz*npml + idx - nxlength + npml] - diffz)/c2z; 

			d_p[idz*nxlength + idx] = d_pr1[idz*npml+ idx - nxlength + npml] + d_pr2[idz*npml + idx - nxlength + npml];
		}
	}
}
__global__ void save_d_vxpml (float *d_vx, float *d_vxspmllr, int nxlength, int nz, int npml)
//
// <((nz-2*npml)+blocksizez-1)/blocksizez,2>  <blocksizez, N/2>
{
	int idz = threadIdx.x + blockDim.x * blockIdx.x;
	int idx = threadIdx.y + blockDim.y * blockIdx.y;

	int radius = N/2;
	int id;
	int nl = radius*(nz-2*npml);
	if (idz < nz-2*npml)
	{
		if (blockIdx.y == 0)  //left 
		{	
			id = idz*radius + idx;
			if (idx < radius)
				d_vxspmllr[id] = d_vx[(idz+npml)*(nxlength - 1) + idx - radius + npml];
		}
		if (blockIdx.y == 1) //right
		{
			id = idz*radius + idx - radius;
			if (idx < N)
				d_vxspmllr[id + nl] = d_vx[(idz+npml)*(nxlength - 1) + idx - radius + nxlength - npml - 1];
		}
	}

}
__global__ void save_d_vzpml (float *d_vz, float *d_vzspmltb, int nxlength, int nz, int npml)
//
// <2,((nxlength-2*npml)+blocksizex-1)/blocksizex>  <N/2,blocksizex >
{
	int idz = threadIdx.x + blockDim.x * blockIdx.x;
	int idx = threadIdx.y + blockDim.y * blockIdx.y;

	int radius = N/2;
	int id;
	int nl = radius*(nxlength-2*npml);
	if (idx < nxlength-2*npml)
	{
		if (blockIdx.x == 0)  //top
		{	
			id = idz*(nxlength-2*npml) + idx;
			if (idz < radius)
				d_vzspmltb[id] = d_vz[(idz - radius + npml)*nxlength + idx + npml];
		}
		if (blockIdx.x == 1) //bottom
		{
			id = (idz-radius)*(nxlength-2*npml) + idx;
			if (idz < N)
				d_vzspmltb[id + nl] = d_vz[(idz - radius + nz -npml -1)*nxlength + idx + npml];
		}
	}

}
__global__ void save_d_ppmllr (float *d_p, float *d_pspmllr, int nxlength, int nz, int npml)
//
// <((nz-2*npml)+blocksizez-1)/blocksizez,2>  <blocksizez, N/2>
{
	int idz = threadIdx.x + blockDim.x * blockIdx.x;
	int idx = threadIdx.y + blockDim.y * blockIdx.y;

	int radius = N/2;
	int id;
	int nl = radius*(nz-2*npml);
	if (idz < nz-2*npml)
	{
		if (blockIdx.y == 0)  //left
		{	
			id = idz*radius + idx;
			if (idx < radius)
				d_pspmllr[id] = d_p[(idz+npml)*nxlength + idx - radius + npml];
		}
		if (blockIdx.y == 1) //right
		{
			id = idz*radius + idx - radius;
			if (idx < N)
				d_pspmllr[id + nl] = d_p[(idz+npml)*nxlength + idx - radius + nxlength - npml];
		}
	}

}
__global__ void save_d_ppmltb (float *d_p, float *d_pspmltb, int nxlength, int nz, int npml)
//
// <2,((nxlength-2*npml)+blocksizex-1)/blocksizex>  <N/2,blocksizex >
{
	int idz = threadIdx.x + blockDim.x * blockIdx.x;
	int idx = threadIdx.y + blockDim.y * blockIdx.y;

	int radius = N/2;
	int id;
	int nl = radius*(nxlength-2*npml);
	if (idx < nxlength-2*npml)
	{
		if (blockIdx.x == 0)  //top
		{	
			id = idz*(nxlength-2*npml) + idx;
			if (idz < radius)
				d_pspmltb[id] = d_p[(idz - radius + npml)*nxlength + idx + npml];
		}
		if (blockIdx.x == 1) //bottom
		{
			id = (idz-radius)*(nxlength-2*npml) + idx;
			if (idz < N)
				d_pspmltb[id + nl] = d_p[(idz - radius + nz -npml)*nxlength + idx + npml];
		}
	}

}
__global__ void read_d_vxpml (float *d_vx, float *d_vxspmllr, int nxlength, int nz, int npml)
//
// <((nz-2*npml)+blocksizez-1)/blocksizez,2>  <blocksizez, N/2>
{
	int idz = threadIdx.x + blockDim.x * blockIdx.x;
	int idx = threadIdx.y + blockDim.y * blockIdx.y;

	int radius = N/2;
	int id;
	int nl = radius*(nz-2*npml);
	if (idz < nz-2*npml)
	{
		if (blockIdx.y == 0)  //left
		{	
			id = idz*radius + idx;
			if (idx < radius)
				d_vx[(idz+npml)*(nxlength - 1) + idx - radius + npml] = d_vxspmllr[id];
		}
		if (blockIdx.y == 1) //right
		{
			id = idz*radius + idx - radius;
			if (idx < N)
				d_vx[(idz+npml)*(nxlength - 1) + idx - radius + nxlength - npml - 1] = d_vxspmllr[id + nl];
		}
	}

}
__global__ void read_d_vzpml (float *d_vz, float *d_vzspmltb, int nxlength, int nz, int npml)
//
// <2,((nxlength-2*npml)+blocksizex-1)/blocksizex>  <N/2,blocksizex >
{
	int idz = threadIdx.x + blockDim.x * blockIdx.x;
	int idx = threadIdx.y + blockDim.y * blockIdx.y;

	int radius = N/2;
	int id;
	int nl = radius*(nxlength-2*npml);
	if (idx < nxlength-2*npml)
	{
		if (blockIdx.x == 0)  //top
		{	
			id = idz*(nxlength-2*npml) + idx;
			if (idz < radius)
				d_vz[(idz - radius + npml)*nxlength + idx + npml] = d_vzspmltb[id];
		}
		if (blockIdx.x == 1) //bottom
		{
			id = (idz-radius)*(nxlength-2*npml) + idx;
			if (idz < N)
				d_vz[(idz - radius + nz -npml -1)*nxlength + idx + npml] = d_vzspmltb[id + nl];
		}
	}

}
__global__ void read_d_ppmllr (float *d_p, float *d_pspmllr, int nxlength, int nz, int npml)
//
// <((nz-2*npml)+blocksizez-1)/blocksizez,2>  <blocksizez, N/2>
{
	int idz = threadIdx.x + blockDim.x * blockIdx.x;
	int idx = threadIdx.y + blockDim.y * blockIdx.y;

	int radius = N/2;
	int id;
	int nl = radius*(nz-2*npml);
	if (idz < nz-2*npml)
	{
		if (blockIdx.y == 0)  //left
		{	
			id = idz*radius + idx;
			if (idx < radius)
				d_p[(idz+npml)*nxlength + idx - radius + npml] = d_pspmllr[id];
		}
		if (blockIdx.y == 1) //right
		{
			id = idz*radius + idx - radius;
			if (idx < N)
				d_p[(idz+npml)*nxlength + idx - radius + nxlength - npml] = d_pspmllr[id + nl];
		}
	}

}
__global__ void read_d_ppmltb (float *d_p, float *d_pspmltb, int nxlength, int nz, int npml)
//
// <2,((nxlength-2*npml)+blocksizex-1)/blocksizex>  <N/2,blocksizex >
{
	int idz = threadIdx.x + blockDim.x * blockIdx.x;
	int idx = threadIdx.y + blockDim.y * blockIdx.y;

	int radius = N/2;
	int id;
	int nl = radius*(nxlength-2*npml);
	if (idx < nxlength-2*npml)
	{
		if (blockIdx.x == 0)  //top
		{	
			id = idz*(nxlength-2*npml) + idx;
			if (idz < radius)
				d_p[(idz - radius + npml)*nxlength + idx + npml] = d_pspmltb[id];
		}
		if (blockIdx.x == 1) //bottom
		{
			id = (idz-radius)*(nxlength-2*npml) + idx;
			if (idz < N)
				d_p[(idz - radius + nz -npml)*nxlength + idx + npml] = d_pspmltb[id + nl];
		}
	}
}
__global__ void cuda_cross_coorelation(float *d_ps,float *d_p,float *d_imagesg, float *d_imagess, int nxpad, int nzpad, int npml)
{
	int idz = threadIdx.x + blockDim.x * blockIdx.x + npml;
	int idx = threadIdx.y + blockDim.y * blockIdx.y + npml;

	int idp = idz*nxpad+idx;
	int idimage = (idz-npml)*(nxpad-2*npml)+idx-npml;


	if (idz < nzpad-npml && idx < nxpad - npml)
	{
		d_imagesg[idimage] += d_ps[idp]*d_p[idp];
		d_imagess[idimage] += d_ps[idp]*d_ps[idp];
	}	

}
__global__ void cuda_wavefield_decomposition(float *d_ps, float *d_vxs, float *d_vzs, float *d_p,  float *d_vx,  float *d_vz,
							   float *d_g2ud, float *d_g2du, float *d_g2lr, float *d_g2rl, int nxpad, int nzpad, int npml)
{
	int idz = threadIdx.x + blockDim.x * blockIdx.x + npml;
	int idx = threadIdx.y + blockDim.y * blockIdx.y + npml;
	int id = idz*nxpad+idx;
	int idpart = (idz-npml)*(nxpad - 2*npml) + idx - npml;

	float psu,psd,psl,psr,pu,pd,pl,pr;
	float vectxs,vectzs,
		vectxr,vectzr;

	if (idz < nzpad-npml && idx < nxpad - npml)
	{
		vectxs = - d_ps[id]*d_vxs[id];
		vectzs = - d_ps[id]*d_vzs[id];
		vectxr = - d_p[id]*d_vx[id];
		vectzr = - d_p[id]*d_vz[id];

		psu = (vectzs < 0)?d_ps[id]:0.0;
		psd = (vectzs > 0)?d_ps[id]:0.0;
		psl = (vectxs < 0)?d_ps[id]:0.0;
		psr = (vectxs > 0)?d_ps[id]:0.0;	

		pu = (vectzr < 0)?d_p[id]:0.0;
		pd = (vectzr > 0)?d_p[id]:0.0;
		pl = (vectxr < 0)?d_p[id]:0.0;
		pr = (vectxr > 0)?d_p[id]:0.0;

		d_g2ud[idpart] += psu*pd;
		d_g2du[idpart] += psd*pu;
		d_g2lr[idpart] += psl*pr;
		d_g2rl[idpart] += psr*pl;		
	}

}
__global__ void cuda_poynting(float *d_ps, float *d_vxs, float *d_vzs, float *d_p,  float *d_vx,  float *d_vz, float *d_vp, float *d_g31, float *d_g32, int nxpad, int nzpad, int npml)
{
	int idz = threadIdx.x + blockDim.x * blockIdx.x + npml;
	int idx = threadIdx.y + blockDim.y * blockIdx.y + npml;
	int id = idz*nxpad+idx;

	float ksx, ksz,
		krx, krz;
	float wtheta;
	int iw;
	float dw = pi/180;


	if (idz < nzpad-npml && idx < nxpad - npml)
	{
		ksx = d_ps[id]*d_vxs[id];
		ksz = d_ps[id]*d_vz[id];

		krx = d_p[id]*d_vx[id];
		krz = d_p[id]*d_vz[id];

		wtheta = (ksx*krx + ksz*krz)/((sqrtf((ksx*ksx+ksz*ksz)*(krx*krx+krz*krz)) > 0)?sqrtf((ksx*ksx+ksz*ksz)*(krx*krx+krz*krz)):eps);

		if (wtheta >= -1.0 && wtheta <= 1.0){
			wtheta = 0.5*acosf(wtheta);
			wtheta = pi/2.0 - wtheta;
			if (wtheta > 0 && wtheta <= pi/3.0)
			{
				iw = (int)(wtheta/dw);
				d_g31[iw*(nzpad-2*npml)*(nxpad-2*npml) + (idz-npml)*(nxpad-2*npml) + idx-npml] += d_ps[id]*d_p[id]*cosf(wtheta);
				d_g32[iw*(nzpad-2*npml)*(nxpad-2*npml) + (idz-npml)*(nxpad-2*npml) + idx-npml] += d_ps[id]*d_p[id]*16.0*pi*sqrtf(d_vp[id])/sinf(wtheta);
			}
		}	
	}

}
__global__ void cuda_energynorm(float *d_ps, float *d_ps_pre, float *d_p,  float *d_p_pre,  float *d_vp, float *d_g4, float dx, float dz, float dt, int nxpad, int nzpad, int npml)
{
	int idz = threadIdx.x + blockDim.x * blockIdx.x + npml;
	int idx = threadIdx.y + blockDim.y * blockIdx.y + npml;
	int idzl = threadIdx.x;
	int idxl = threadIdx.y;
	int tidz  =  idzl + 1;
	int tidx  =  idxl + 1;
	int id = idz*nxpad+idx;
	int idpart = (idz-npml)*(nxpad-2*npml)+idx-npml;
	int radius = 1;
	float psx,psz,px,pz;
	__shared__ float p[Block_Sizez + 2][Block_Sizex + 2];
	__shared__ float ps[Block_Sizez + 2][Block_Sizex + 2];	


	// source wavefield
	if (idxl < radius)                     ps[tidz][idxl]          = d_ps[idz*nxpad + idx - radius];
	if (idxl >= blockDim.y - radius)       ps[tidz][idxl + 2]      = d_ps[idz*nxpad + idx + radius];
	if (idzl < radius)                     ps[idzl][tidx]          = d_ps[(idz - radius)*nxpad + idx];
	if (idzl >= blockDim.y - radius)       ps[idzl + 2][tidx]      = d_ps[(idz + radius)*nxpad + idx];
	// receiver wavefield
	if (idxl < radius)                     p[tidz][idxl]          = d_p[idz*nxpad + idx - radius];
	if (idxl >= blockDim.y - radius)       p[tidz][idxl + 2]      = d_p[idz*nxpad + idx + radius];
	if (idzl < radius)                     p[idzl][tidx]          = d_p[(idz - radius)*nxpad + idx];
	if (idzl >= blockDim.y - radius)       p[idzl + 2][tidx]      = d_p[(idz + radius)*nxpad + idx];

	ps[tidz][tidx]  = d_ps[id];
	p[tidz][tidx]   = d_p[id];

	__syncthreads();


	if (idz < nzpad-npml && idx < nxpad - npml)
	{
		psx = 0.5*(ps[tidz][tidx + 1] - ps[tidz][tidx - 1])/dx;
		psz = 0.5*(ps[tidz + 1][tidx] - ps[tidz - 1][tidx])/dz;
		px = 0.5*(p[tidz][tidx + 1] - p[tidz][tidx - 1])/dx;
		pz = 0.5*(p[tidz + 1][tidx] - p[tidz - 1][tidx])/dz;

		d_g4[idpart] += -(d_ps[id]-d_ps_pre[id])*(d_p[id]-d_p_pre[id])/(dt*dt*d_vp[id])
									  + psx*px + psz*pz;
	}

}
__global__ void cuda_stack_udlr(float *d_g2ud, float *d_g2du, float *d_g2lr, float *d_g2rl, float *d_g2, int nx, int nz)
{
	int idz = threadIdx.x + blockDim.x * blockIdx.x;
	int idx = threadIdx.y + blockDim.y * blockIdx.y;
	int id = idz*nx+idx;

	if (idz < nz && idx < nx){
		d_g2[id] = d_g2ud[id] + d_g2du[id] + d_g2lr[id] + d_g2rl[id];
	}

}
__global__ void cuda_stack_theta(float *d_g31, float *d_g3, int nx, int nz)
{
	int idz = threadIdx.x + blockDim.x * blockIdx.x;
	int idx = threadIdx.y + blockDim.y * blockIdx.y;
	int id = idz*nx+idx;

	if (idz < nz && idx < nx)
		for (int iw=0;iw<=60;iw++)
			d_g3[id] += d_g31[iw*nx*nz+id];

}
__global__ void cuda_applyic(float *d_image, float *d_imagesg,  int nx, int nxlength, int nx1, int nz)
{
	int idz = threadIdx.x + blockDim.x * blockIdx.x;
	int idx = threadIdx.y + blockDim.y * blockIdx.y;
	int idpart = idz*nxlength + idx;
	int idfull = idz*nx + idx + nx1;

	if (idz < nz && idx < nxlength)
	{
		d_image[idfull] += d_imagesg[idpart];
	}	
}
__global__ void cuda_applyics(float *d_image, float *d_imagesg, float *d_imagess, int nx, int nxlength, int nx1, int nz)
{
	int idz = threadIdx.x + blockDim.x * blockIdx.x;
	int idx = threadIdx.y + blockDim.y * blockIdx.y;
	int idpart = idz*nxlength + idx;
	int idfull = idz*nx + idx + nx1;

	if (idz < nz && idx < nxlength)
	{
		d_image[idfull] += d_imagesg[idpart]/(d_imagess[idpart] + eps);
	}	
}
__global__ void cuda_applyimagingconditionmulti(float *d_Image, float *d_Imagesg, float *d_Imagess, int nx, int nz)
{
	int idz = threadIdx.x + blockDim.x * blockIdx.x;
	int idx = threadIdx.y + blockDim.y * blockIdx.y;
	int idfull = idz*nx + idx;

	if (idz < nz && idx < nx)
	{	
		d_Image[idfull] = d_Imagesg[idfull]/(d_Imagess[idfull] + eps);
	}	
}
__global__ void cuda_taper_calculate(float *d_taper, int nxlength, int nz, int nwin, float alpha)
{
	int idz = threadIdx.x + blockDim.x * blockIdx.x;
	int idx = threadIdx.y + blockDim.y * blockIdx.y;
	int id = idz*nxlength + idx;
	if (idz < nz){
	if (idx < nwin)
		d_taper[id] = expf(-alpha*alpha*(nwin-idx)*(nwin-idx));
	if (idx >= nxlength - nwin && idx < nxlength)
		d_taper[id] = expf(-alpha*alpha*(idx - nxlength + nwin + 1)*(idx - nxlength + nwin + 1));
	if (idx >= nwin && idx < nxlength - nwin)
		d_taper[id] = 1.0;
	}
}
__global__ void cuda_taper_image(float *d_g, float *d_taper, int nxlength, int nz)
{
	int idz = threadIdx.x + blockDim.x * blockIdx.x;
	int idx = threadIdx.y + blockDim.y * blockIdx.y;
	int id = idz*nxlength + idx;

	if (idz < nz && idx < nxlength)
	{
		d_g[id] *= d_taper[id];
	}	
}
void check_gpu_error(const char *msg)
// < check GPU errors >
{
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err){
		printf("Cuda error: %s: %s",msg,cudaGetErrorString(err));
		exit(0);
	}
}
