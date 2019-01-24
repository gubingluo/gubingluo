#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "function.h"
#include "sgyhead.h"
#include "common.h"
#include "kernel.cuh"

// shotprofileId[5][nshot]
// 1 行: sx             shotprofileId[0][ishot]
// 2 行: ntrace         shotprofileId[1][ishot]
// 3 行: itraceAll      shotprofileId[2][ishot]
// 4 行: offset         shotprofileId[3][ishot]

extern "C" void rtm_real(int idevice, int nshot, int startshot, int dshot, int medium_flag, int flag_layout, int sgyflag, int endianflag, int ismth,
	                   int nx, int nz, int nt, float dx, float dz, float dt, int nx0, int nz0, int nt0, float dx0, float dz0, float dt0,
		             int npml, int lx, int lz, float tsnap, float fdom, float spz0, float amp, float alp, float direct, float direct0,  
		             float *vp0, float *rho0, float *diffcoef, float *wavelet, char filename_shotgather[40], int **shotprofileId, 
				 float *image1,  float *image2s, float *image3s, float *image4s, float *image5s, 
		 		 float *image2m, float *image3m, float *image4m, float *image5m, float *illum)
{
// common parameter
	int nxpad,nzpad;
	int ntsnap = (int)(tsnap/dt);
	int nw = (int)(direct/(fdom*dt)),
	    tlength = (int)(direct0/dt);
	int it,ix,ishot;
	long long offset_trace;
	FILE *fp;	
// wavefield extrapolation parameter
	float mstimer;
	char buffrecord[40];
	float vpmax;
	float spx,spz;

	static float dx2,dz2,_dt,_dtx,_dtz;
	static int nsx,nsz,nxz;

	static dim3 dimGrid,dimBlock,
		      dimGridp,dimGridvx,dimGridvz,dimGridvxb,dimGridvzb,
		      dimGridvxlr,dimGridvztb,
		      dimGridplr,dimGridptb,
			dimGridpcooner,
		      dimGridrecord,
			dimGridpmllr,dimGridpmltb,dimBlockpmllr,dimBlockpmltb,
			dimGridfull;

	// variables on host
	float *vp,*rho,*temp,*record,*tr;
	float *vxspmllr,*vzspmltb,*pspmllr,*pspmltb;
	float *p;

	// variables on device
	float *d_wavelet,*d_diffcoef;
	float *d_source,*d_record,*d_vp,*d_rho;
	float *d_p,*d_vx,*d_vz,*d_p_pre,
		*d_ps,*d_vxs,*d_vzs,*d_ps_pre;
	float *d_pl1,*d_pl2,*d_pr1,*d_pr2,
	      *d_pt1,*d_pt2,*d_pb1,*d_pb2;
	float *d_ddx,*d_ddz,*d_ddxVx,*d_ddzVz;
	int *d_norder,*d_norderx,*d_norderz;

	float *d_vxspmllr,*d_vzspmltb,*d_pspmllr,*d_pspmltb;
	// single shot image
	float *d_g1,                                 // cross-coorelation
		*d_g2,*d_g2ud,*d_g2du,*d_g2lr,*d_g2rl, // wavefiled-decomposition
		*d_g3,*d_g3_true,*d_g31,*d_g32,                   // poynting vector   d_g32 true amplitude
		*d_g4;                                 // energy norm

	float *d_image1,
		*d_image2s,*d_image3s,*d_image4s,*d_image5s,
		*d_image2m,*d_image3m,*d_image4m,*d_image5m,
		*d_imagetrue,*d_illum; 
	float *imagetrue,*d_Illum;        
	float *d_taper;
	int nwin=50;
	float alpha=0.06;

	float **seisobs,**seiscal;

      _dtx = dt/dx;
	_dtz = dt/dz;
	_dt  = 1.0/dt;
	dx2 = dx*dx;
	dz2 = dz*dz;

	cudaSetDevice(idevice);
	check_gpu_error("Failed to initialize device");
	
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	cudaMalloc(&d_wavelet,       nt*sizeof(float));
	cudaMemset(d_wavelet,   0,   nt*sizeof(float));
	cuda_ricker_wavelet<<<(nt+511)/512,512>>>(d_wavelet, fdom, dt, nt);

	cudaMalloc(&d_image1,       nx*nz*sizeof(float));
	cudaMalloc(&d_image2s,      nx*nz*sizeof(float));	
	cudaMalloc(&d_image3s,      nx*nz*sizeof(float));	
	cudaMalloc(&d_image4s,      nx*nz*sizeof(float));	
	cudaMalloc(&d_image5s,      nx*nz*sizeof(float));	
	cudaMalloc(&d_image2m,      nx*nz*sizeof(float));	
	cudaMalloc(&d_image3m,      nx*nz*sizeof(float));	
	cudaMalloc(&d_image4m,      nx*nz*sizeof(float));	
	cudaMalloc(&d_image5m,      nx*nz*sizeof(float));	
	cudaMalloc(&d_Illum,        nx*nz*sizeof(float));
	cudaMalloc(&d_imagetrue,    nx*nz*sizeof(float));	
	cudaMemset(d_image1,    0,  nx*nz*sizeof(float));
	cudaMemset(d_image2s,   0,  nx*nz*sizeof(float));
	cudaMemset(d_image3s,   0,  nx*nz*sizeof(float));
	cudaMemset(d_image4s,   0,  nx*nz*sizeof(float));
	cudaMemset(d_image5s,   0,  nx*nz*sizeof(float));
	cudaMemset(d_image2m,   0,  nx*nz*sizeof(float));
	cudaMemset(d_image3m,   0,  nx*nz*sizeof(float));
	cudaMemset(d_image4m,   0,  nx*nz*sizeof(float));
	cudaMemset(d_image5m,   0,  nx*nz*sizeof(float));
	cudaMemset(d_Illum,     0,  nx*nz*sizeof(float));
	cudaMemset(d_imagetrue, 0,  nx*nz*sizeof(float));

	imagetrue=(float *)malloc(nx*nz*sizeof(float));
	memset(imagetrue, 0, nx*nz*sizeof(float));
	for (ishot = startshot; ishot <= nshot; ishot=ishot+dshot)
	{
		cudaEventRecord(start);
		//==============================================================================
		//单边排列（排列位于右侧，最小偏移距为offset）
		int noffset,nx1,nx2,nxlength;
		int noffset0,nx10,nx20,nxlength0;
		float offsetmax,offset;

		offsetmax = (float)(shotprofileId[3][ishot-1] + (shotprofileId[1][ishot-1] - 1)*dx0);
		offset    = (float)(shotprofileId[3][ishot-1]);
		spx       = (float)(shotprofileId[0][ishot-1]);
		spz       = spz0;
		// given
		noffset0 = (int)(offset/dx0);
		nx10 = (int)(spx/dx0);
		nx20 = (int)((spx + offsetmax)/dx0);
		nxlength0 = nx20 - nx10 + 1;
		// optimal
		noffset = (int)(offset/dx);
		nx1 = (int)(spx/dx);
		nx2 = (int)((spx + offsetmax)/dx);
		nxlength = nx2 - nx1 + 1;
		// optimal extend
		spx = spx - nx1*dx;   // local
		spx = spx + npml*dx;
		spz = spz + npml*dz;
		nxpad = nxlength + 2*npml;
		nzpad = nz + 2*npml;
			
		nsx = (int)(spx/dx);				
		nsz = (int)(spz/dz);
		nxz = nxpad*nzpad;
		// block and thread
		dimBlock = dim3(Block_Sizez, Block_Sizex);
		dimBlockpmllr = dim3(Block_Sizez,N/2);
		dimBlockpmltb = dim3(N/2,Block_Sizex);
		dimGrid  = dim3((nzpad+Block_Sizez-1)/Block_Sizez,           (nxpad+Block_Sizex-1)/Block_Sizex); 
            dimGridp  = dim3((nzpad-2*npml+Block_Sizez-1)/Block_Sizez,   (nxpad-2*npml+Block_Sizex-1)/Block_Sizex);
		dimGridvx = dim3((nzpad+Block_Sizez-1)/Block_Sizez,          (nxpad-2*npml-1+Block_Sizex-1)/Block_Sizex);
		dimGridvz = dim3((nzpad-2*npml-1+Block_Sizez-1)/Block_Sizez, (nxpad+Block_Sizex-1)/Block_Sizex);
		dimGridvxb = dim3((nzpad-2*npml+Block_Sizez-1)/Block_Sizez,  (nxpad-2*npml-1+Block_Sizex-1)/Block_Sizex);
		dimGridvzb = dim3((nzpad-2*npml-1+Block_Sizez-1)/Block_Sizez,(nxpad-2*npml+Block_Sizex-1)/Block_Sizex);

		dimGridvxlr = dim3((nzpad+Block_Sizez-1)/Block_Sizez,2);
		dimGridvztb = dim3(2,(nxpad+Block_Sizex-1)/Block_Sizex);
		dimGridplr  = dim3((nzpad-2*npml+Block_Sizez-1)/Block_Sizez,2);
		dimGridptb  = dim3(2,(nxpad-2*npml+Block_Sizex-1)/Block_Sizex);
		dimGridpcooner = dim3(2,2);

		dimGridrecord = dim3((nt+Block_Sizez-1)/Block_Sizez,     (nxpad-2*npml-noffset+Block_Sizex-1)/Block_Sizex);
		dimGridpmllr  = dim3((nzpad-2*npml+Block_Sizez-1)/Block_Sizez,2);
		dimGridpmltb  = dim3(2,(nxpad-2*npml+Block_Sizex-1)/Block_Sizex);
		//模拟参数换算与准备
		record = (float *)malloc(nt*(nxlength-noffset)*sizeof(float));
		temp = (float *)malloc(nz*nxlength*sizeof(float));
		vp = (float *)malloc(nzpad*nxpad*sizeof(float));
		rho= (float *)malloc(nzpad*nxpad*sizeof(float));			
		p  = (float *)malloc(nzpad*nxpad*sizeof(float));
		vxspmllr = (float *)malloc(N*(nzpad-2*npml)*nt*sizeof(float));
		vzspmltb = (float *)malloc(N*(nxpad-2*npml)*nt*sizeof(float));
		pspmllr = (float *)malloc(N*(nzpad-2*npml)*nt*sizeof(float));
		pspmltb = (float *)malloc(N*(nxpad-2*npml)*nt*sizeof(float));
		tr = (float *)malloc(nt0*sizeof(float));
						
		memset(record, 0, nt*(nxlength-noffset)*sizeof(float));
		memset(temp,   0, nz*nxlength*sizeof(float));
		memset(vp,     0, nxpad*nzpad*sizeof(float));
		memset(rho,    0, nxpad*nzpad*sizeof(float));
		memset(p,      0, nxpad*nzpad*sizeof(float));
		memset(vxspmllr,0, N*(nzpad-2*npml)*nt*sizeof(float));
		memset(vzspmltb,0, N*(nxpad-2*npml)*nt*sizeof(float));
		memset(pspmllr, 0, N*(nzpad-2*npml)*nt*sizeof(float));
		memset(pspmltb, 0, N*(nxpad-2*npml)*nt*sizeof(float));
		//===============================================================================
		extractvel1(temp, vp0,  nx, nz, nx1, nx2);
		extendvel1(vp,  temp, nxlength, nz, npml);	
		extractrho1(temp, rho0, nx, nz, nx1, nx2);					
		extendvel1(rho, temp, nxlength, nz, npml);
		free(temp);
		// pml layers smooth
		if (medium_flag){
			pmlvelsmooth1d(vp,  nxpad, nzpad, npml);
				pmlvelsmooth1d(rho, nxpad, nzpad, npml);}		
		vpmax = sqrtf(Maxval1(vp, nzpad*nxpad));

		// alloc device memory
		cudaMalloc(&d_diffcoef,   (N/2)*(N/2)*sizeof(float));
		cudaMalloc(&d_record,     (nxlength-noffset)*nt*sizeof(float));
		cudaMalloc(&d_source,     nxz*sizeof(float));
		cudaMalloc(&d_vp,         nxz*sizeof(float));
		cudaMalloc(&d_rho,        nxz*sizeof(float));

		cudaMalloc(&d_p,          nxz*sizeof(float));
		cudaMalloc(&d_p_pre,      nxz*sizeof(float));
		cudaMalloc(&d_vx,         (nxpad-1)*nzpad*sizeof(float));
		cudaMalloc(&d_vz,         nxpad*(nzpad-1)*sizeof(float));
		cudaMalloc(&d_ps,         nxz*sizeof(float));
		cudaMalloc(&d_ps_pre,     nxz*sizeof(float));
		cudaMalloc(&d_vxs,        (nxpad-1)*nzpad*sizeof(float));
		cudaMalloc(&d_vzs,        nxpad*(nzpad-1)*sizeof(float));
		
		cudaMalloc(&d_pl1,        npml*nzpad*sizeof(float));
		cudaMalloc(&d_pl2,        npml*nzpad*sizeof(float));
		cudaMalloc(&d_pr1,        npml*nzpad*sizeof(float));
		cudaMalloc(&d_pr2,        npml*nzpad*sizeof(float));
		cudaMalloc(&d_pt1,        npml*(nxpad-2*npml)*sizeof(float));
		cudaMalloc(&d_pt2,        npml*(nxpad-2*npml)*sizeof(float));
		cudaMalloc(&d_pb1,        npml*(nxpad-2*npml)*sizeof(float));
		cudaMalloc(&d_pb2,        npml*(nxpad-2*npml)*sizeof(float));

		cudaMalloc(&d_ddx,        nxpad*sizeof(float));
		cudaMalloc(&d_ddz,        nzpad*sizeof(float));
		cudaMalloc(&d_ddxVx,      (nxpad-1)*sizeof(float));
		cudaMalloc(&d_ddzVz,      (nzpad-1)*sizeof(float));
		cudaMalloc(&d_norder,     nxz*sizeof(int));
		cudaMalloc(&d_norderx,    (nxpad-1)*sizeof(int));
		cudaMalloc(&d_norderz,    (nzpad-1)*sizeof(int));

		cudaMalloc(&d_vxspmllr,   N*(nzpad-2*npml)*sizeof(float));
		cudaMalloc(&d_vzspmltb,   N*(nxpad-2*npml)*sizeof(float));
		cudaMalloc(&d_pspmllr,    N*(nzpad-2*npml)*sizeof(float));
		cudaMalloc(&d_pspmltb,    N*(nxpad-2*npml)*sizeof(float));	
			
		cudaMalloc(&d_g1,          (nxpad-2*npml)*(nzpad-2*npml)*sizeof(float));	
		cudaMalloc(&d_g2,          (nxpad-2*npml)*(nzpad-2*npml)*sizeof(float));		
		cudaMalloc(&d_g3,          (nxpad-2*npml)*(nzpad-2*npml)*sizeof(float));	
		cudaMalloc(&d_g3_true,     (nxpad-2*npml)*(nzpad-2*npml)*sizeof(float));	
		cudaMalloc(&d_g4,          (nxpad-2*npml)*(nzpad-2*npml)*sizeof(float));
		cudaMalloc(&d_g2ud,        (nxpad-2*npml)*(nzpad-2*npml)*sizeof(float));	
		cudaMalloc(&d_g2du,        (nxpad-2*npml)*(nzpad-2*npml)*sizeof(float));
		cudaMalloc(&d_g2lr,        (nxpad-2*npml)*(nzpad-2*npml)*sizeof(float));
		cudaMalloc(&d_g2rl,        (nxpad-2*npml)*(nzpad-2*npml)*sizeof(float));
		cudaMalloc(&d_g31,         (nxpad-2*npml)*(nzpad-2*npml)*61*sizeof(float));	
		cudaMalloc(&d_g32,         (nxpad-2*npml)*(nzpad-2*npml)*61*sizeof(float));
		cudaMalloc(&d_illum,       (nxpad-2*npml)*(nzpad-2*npml)*sizeof(float));
		cudaMalloc(&d_taper,       (nxpad-2*npml)*(nzpad-2*npml)*sizeof(float));	

		cudaMemcpy(d_diffcoef, diffcoef, (N/2)*(N/2)*sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(d_vp,       vp,       nxz*sizeof(float),         cudaMemcpyHostToDevice);
		cudaMemcpy(d_rho,      rho,      nxz*sizeof(float),         cudaMemcpyHostToDevice);
		
		cudaMemset(d_record, 0, (nxlength-noffset)*nt*sizeof(float));
		cudaMemset(d_p,      0, nxpad*nzpad*sizeof(float));
		cudaMemset(d_vx,     0, (nxpad-1)*nzpad*sizeof(float));
		cudaMemset(d_vz,     0, nxpad*(nzpad-1)*sizeof(float));
		cudaMemset(d_p_pre,  0, nxpad*nzpad*sizeof(float));
		cudaMemset(d_ps,     0, nxpad*nzpad*sizeof(float));
		cudaMemset(d_vxs,    0, (nxpad-1)*nzpad*sizeof(float));
		cudaMemset(d_vzs,    0, nxpad*(nzpad-1)*sizeof(float));
		cudaMemset(d_ps_pre, 0, nxpad*nzpad*sizeof(float));		

		cudaMemset(d_pl1,    0, npml*nzpad*sizeof(float));
		cudaMemset(d_pl2,    0, npml*nzpad*sizeof(float));
		cudaMemset(d_pr1,    0, npml*nzpad*sizeof(float));
		cudaMemset(d_pr2,    0, npml*nzpad*sizeof(float));
		cudaMemset(d_pt1,    0, npml*(nxpad-2*npml)*sizeof(float));
		cudaMemset(d_pt2,    0, npml*(nxpad-2*npml)*sizeof(float));
		cudaMemset(d_pb1,    0, npml*(nxpad-2*npml)*sizeof(float));
		cudaMemset(d_pb2,    0, npml*(nxpad-2*npml)*sizeof(float));
		cudaMemset(d_ddx,    0, nxpad*sizeof(float));
		cudaMemset(d_ddz,    0, nzpad*sizeof(float));
		cudaMemset(d_ddxVx,  0, (nxpad-1)*sizeof(float));
		cudaMemset(d_ddzVz,  0, (nzpad-1)*sizeof(float));
		cudaMemset(d_norder, 0, nxpad*nzpad*sizeof(int));
		cudaMemset(d_norderx,0, (nxpad-1)*sizeof(int));
		cudaMemset(d_norderz,0, (nzpad-1)*sizeof(int));
		cudaMemset(d_vxspmllr,0, N*(nzpad-2*npml)*sizeof(float));
		cudaMemset(d_vzspmltb,0, N*(nxpad-2*npml)*sizeof(float));
		cudaMemset(d_pspmllr, 0, N*(nzpad-2*npml)*sizeof(float));
		cudaMemset(d_pspmltb, 0, N*(nxpad-2*npml)*sizeof(float));

		cudaMemset(d_g1,      0, (nxpad-2*npml)*(nzpad-2*npml)*sizeof(float));	
		cudaMemset(d_g2,      0, (nxpad-2*npml)*(nzpad-2*npml)*sizeof(float));	
		cudaMemset(d_g3,      0, (nxpad-2*npml)*(nzpad-2*npml)*sizeof(float));
		cudaMemset(d_g3_true, 0, (nxpad-2*npml)*(nzpad-2*npml)*sizeof(float));
		cudaMemset(d_g4,      0, (nxpad-2*npml)*(nzpad-2*npml)*sizeof(float));
		cudaMemset(d_g2ud,    0, (nxpad-2*npml)*(nzpad-2*npml)*sizeof(float));
		cudaMemset(d_g2du,    0, (nxpad-2*npml)*(nzpad-2*npml)*sizeof(float));
		cudaMemset(d_g2lr,    0, (nxpad-2*npml)*(nzpad-2*npml)*sizeof(float));
		cudaMemset(d_g2rl,    0, (nxpad-2*npml)*(nzpad-2*npml)*sizeof(float));
		cudaMemset(d_g31,     0, (nxpad-2*npml)*(nzpad-2*npml)*61*sizeof(float));
		cudaMemset(d_g32,     0, (nxpad-2*npml)*(nzpad-2*npml)*61*sizeof(float));
		cudaMemset(d_illum,   0, (nxpad-2*npml)*(nzpad-2*npml)*sizeof(float));
		cudaMemset(d_taper,   0, (nxpad-2*npml)*(nzpad-2*npml)*sizeof(float));

		cuda_source<<<dimGrid,dimBlock>>>(d_source, nsx, nsz, nxpad, nzpad, amp, alp, dx2, dz2);
		cuda_pmlCoeffpx<<<(nxpad+127)/128,128>>>(d_ddx, vpmax, dx, npml, nxpad);
		cuda_pmlCoeffpz<<<(nzpad+127)/128,128>>>(d_ddz, vpmax, dz, npml, nzpad);
		cuda_pmlCoeffvx<<<(nxpad+126)/128,128>>>(d_ddxVx, vpmax, dx, npml, nxpad-1);
		cuda_pmlCoeffvz<<<(nzpad+126)/128,128>>>(d_ddzVz, vpmax, dz, npml, nzpad-1);
		cuda_norder<<<dimGrid,dimBlock>>>       (d_norder, nxpad, nzpad);
		cuda_norderx<<<(nxpad+126)/128,128>>>   (d_norderx, nxpad-1);
		cuda_norderz<<<(nzpad+126)/128,128>>>   (d_norderz, nzpad-1);	
		while (2*nwin >= nxlength)
			nwin -= 3;
		cuda_taper_calculate<<<dimGridp,dimBlock>>>(d_taper, nxpad-2*npml, nzpad-2*npml, nwin, alpha);
		printf("N0. %d shot source wavefield calculating......\n",ishot);			
		// calculate source wavefiled to save pml layer
		for (it=0; it<nt; it++)
		{
			if (it%ntsnap == 0){
				cudaMemcpy(p, d_p, nxz*sizeof(float), cudaMemcpyDeviceToHost);
				printf("source-the current shot: %d\ttime: %f s; wavefield: %.5e\n",ishot,it*dt, absMaxval1(p, nxz));}
			cuda_forward_vx<<<dimGridvx,dimBlock>>>(d_p, d_vx, d_rho, d_diffcoef, _dtx, npml, nxpad, nzpad);
			cuda_forward_vz<<<dimGridvz,dimBlock>>>(d_p, d_vz, d_rho, d_diffcoef, _dtz, npml, nxpad, nzpad);
			cuda_pml_vxlr<<<dimGridvxlr,dimBlock>>>(d_p, d_vx, d_rho, d_diffcoef, d_ddxVx, _dtx, dt, npml, nxpad, nzpad, d_norderx);
			cuda_pml_vztb<<<dimGridvztb,dimBlock>>>(d_p, d_vz, d_rho, d_diffcoef, d_ddzVz, _dtz, dt, npml, nxpad, nzpad, d_norderz);
			cuda_forward_p<<<dimGridp,dimBlock>>>(d_p, d_vx, d_vz, d_rho, d_vp, d_diffcoef, _dtx, _dtz, npml, nxpad, nzpad);
			cuda_pml_plr<<<dimGridplr,dimBlock>>>(d_p, d_vx, d_vz, d_pl1, d_pl2, d_pr1, d_pr2, d_rho, d_vp, d_diffcoef, d_ddx, d_ddz,_dtx, _dtz, dt, npml, nxpad, nzpad, d_norder);
			cuda_pml_ptb<<<dimGridptb,dimBlock>>>(d_p, d_vx, d_vz, d_pt1, d_pt2, d_pb1, d_pb2, d_rho, d_vp, d_diffcoef, d_ddz, _dtx, _dtz, dt, npml, nxpad, nzpad, d_norder);
			cuda_pml_pconner<<<dimGridpcooner,dimBlock>>>(d_p, d_vx, d_vz, d_pl1, d_pl2, d_pr1, d_pr2, d_rho, d_vp, d_diffcoef, d_ddx, d_ddz,_dtx, _dtz, dt, npml, nxpad, nzpad, d_norder);
			cuda_add_source<<<dimGrid,dimBlock>>>(d_p, d_source, d_wavelet, dt, 1, nxpad, nzpad, it);

			save_d_vxpml<<<dimGridpmllr,dimBlockpmllr>>>(d_vx, d_vxspmllr, nxpad, nzpad, npml);
			save_d_vzpml<<<dimGridpmltb,dimBlockpmltb>>>(d_vz, d_vzspmltb, nxpad, nzpad, npml);
			save_d_ppmllr<<<dimGridpmllr,dimBlockpmllr>>>(d_p, d_pspmllr,  nxpad, nzpad, npml);
			save_d_ppmltb<<<dimGridpmltb,dimBlockpmltb>>>(d_p, d_pspmltb,  nxpad, nzpad, npml);
			cudaMemcpy(&vxspmllr[it*N*(nzpad-2*npml)],       d_vxspmllr, N*(nzpad-2*npml)*sizeof(float),       cudaMemcpyDeviceToHost);
			cudaMemcpy(&vzspmltb[it*N*(nxpad-2*npml)],       d_vzspmltb, N*(nxpad-2*npml)*sizeof(float),       cudaMemcpyDeviceToHost);	
			cudaMemcpy(&pspmllr[it*N*(nzpad-2*npml)],        d_pspmllr,  N*(nzpad-2*npml)*sizeof(float),       cudaMemcpyDeviceToHost);
			cudaMemcpy(&pspmltb[it*N*(nxpad-2*npml)],        d_pspmltb,  N*(nxpad-2*npml)*sizeof(float),       cudaMemcpyDeviceToHost);
		}

		// initial source wavefiled
		// save last snap used to reconstruction source wavefield
		cudaMemcpy(d_ps, d_p,  nxz*sizeof(float), cudaMemcpyDeviceToDevice);
		// initial receiver wavefield
		cudaMemset(d_p,      0, nxpad*nzpad*sizeof(float));
		cudaMemset(d_vx,     0, (nxpad-1)*nzpad*sizeof(float));
		cudaMemset(d_vz,     0, nxpad*(nzpad-1)*sizeof(float));
		cudaMemset(d_pl1,    0, npml*nzpad*sizeof(float));
		cudaMemset(d_pl2,    0, npml*nzpad*sizeof(float));
		cudaMemset(d_pr1,    0, npml*nzpad*sizeof(float));
		cudaMemset(d_pr2,    0, npml*nzpad*sizeof(float));
		cudaMemset(d_pt1,    0, npml*(nxpad-2*npml)*sizeof(float));
		cudaMemset(d_pt2,    0, npml*(nxpad-2*npml)*sizeof(float));
		cudaMemset(d_pb1,    0, npml*(nxpad-2*npml)*sizeof(float));
		cudaMemset(d_pb2,    0, npml*(nxpad-2*npml)*sizeof(float));
		cudaMemset(d_vxspmllr,0, N*(nzpad-2*npml)*sizeof(float));
		cudaMemset(d_vzspmltb,0, N*(nxpad-2*npml)*sizeof(float));
		cudaMemset(d_pspmllr, 0, N*(nzpad-2*npml)*sizeof(float));
		cudaMemset(d_pspmltb, 0, N*(nxpad-2*npml)*sizeof(float));		


		if (ismth)
		{
			velsmooth1d(vp,  nzpad, nxpad, ismth);
			velsmooth1d(rho, nzpad, nxpad, ismth);
			cudaMemcpy(d_vp,       vp,       nxz*sizeof(float),         cudaMemcpyHostToDevice);
			cudaMemcpy(d_rho,      rho,      nxz*sizeof(float),         cudaMemcpyHostToDevice);
		}
		//===============================================================================
		// prepare seismic profile
		seisobs=Creat2dArray(nt0,nxlength0-noffset0);
		seiscal=Creat2dArray(nt, nxlength-noffset);
		sprintf(buffrecord,"./seisReal/%s",filename_shotgather);		
		offset_trace = (shotprofileId[2][ishot-1] - 1)*(240 + nt0*sizeof(float)) + 3600*sgyflag;
		fp = fopen(buffrecord,"rb");
		fseek(fp,offset_trace,0);
		for (ix=0; ix<shotprofileId[1][ishot-1]; ix++)
		{		
			fseek(fp,240L,1);
			fread(tr,sizeof(float),nt0,fp);
			float_to_float(tr, nt0, endianflag);
			for (it=0; it<nt0; it++)
				seisobs[it][ix] = tr[it];

		}
		fclose(fp);
		Interpseis2d(seiscal,seisobs,nxlength-noffset,nxlength0-noffset0,nt,nt0,dt,dt0);
		for (it=0;it<nt;it++)
			for (ix=0;ix<nxlength-noffset;ix++)
				record[it*(nxlength-noffset)+ix] = seiscal[it][ix];

		free2dArray(seiscal, nt,  nxlength-noffset);
		free2dArray(seisobs, nt0, nxlength0-noffset0);
		cudaMemcpy(d_record,  record,    (nxlength-noffset)*nt*sizeof(float), cudaMemcpyHostToDevice);
		cuda_mute2<<<dimGridrecord,dimBlock>>>(d_record, d_vp, nsx, nsz, nt, npml, nxlength, noffset, nw, tlength, fdom, dx2, dz2, _dt);		
		// implement RTM				
		// insert seismic record for the last time slice
		cuda_insert_record2<<<(nxlength-noffset+127)/128,128>>>(d_p, &d_record[(nt-1)*(nxlength - noffset)], npml, nxlength, noffset, dt);
		// backforward record
		printf("source wavefiled prepared over...\nBegin backward......\n");
		for (it=nt-2; it>=0; it--)
		{
			if (it%ntsnap == 0){
				cudaMemcpy(p, d_p, nxz*sizeof(float), cudaMemcpyDeviceToHost);
				printf("backward-the current shot: %d\ttime: %f s; wavefield: %.5e\n",ishot, it*dt, absMaxval1(p, nxz));}	

			cudaMemcpy(d_ps_pre, d_ps, nxz*sizeof(float), cudaMemcpyDeviceToDevice);	
			cudaMemcpy(d_p_pre,  d_p,  nxz*sizeof(float), cudaMemcpyDeviceToDevice);
			
			// source wavefield 1: read vx vz pml; 2: calculate inner vx vz; 3: read p pml; 4; calculate inner p
			// 1  --  2
			cudaMemcpy(d_vxspmllr, &vxspmllr[(it+1)*N*(nzpad-2*npml)],       N*(nzpad-2*npml)*sizeof(float),       cudaMemcpyHostToDevice);
			cudaMemcpy(d_vzspmltb, &vzspmltb[(it+1)*N*(nxpad-2*npml)],       N*(nxpad-2*npml)*sizeof(float),       cudaMemcpyHostToDevice);				
			read_d_vxpml<<<dimGridpmllr,dimBlockpmllr>>>(d_vxs, d_vxspmllr, nxpad, nzpad, npml);
			read_d_vzpml<<<dimGridpmltb,dimBlockpmltb>>>(d_vzs, d_vzspmltb, nxpad, nzpad, npml);					
			cuda_backward_vx<<<dimGridvxb,dimBlock>>>(d_ps, d_vxs, d_rho, d_diffcoef, _dtx, npml, nxpad, nzpad);
			cuda_backward_vz<<<dimGridvzb,dimBlock>>>(d_ps, d_vzs, d_rho, d_diffcoef, _dtz, npml, nxpad, nzpad);
			// 3  --  4
			cudaMemcpy(d_pspmllr, &pspmllr[it*N*(nzpad-2*npml)],        N*(nzpad-2*npml)*sizeof(float),       cudaMemcpyHostToDevice);
			cudaMemcpy(d_pspmltb, &pspmltb[it*N*(nxpad-2*npml)],        N*(nxpad-2*npml)*sizeof(float),       cudaMemcpyHostToDevice);
			read_d_ppmllr<<<dimGridpmllr,dimBlockpmllr>>>(d_ps, d_pspmllr, nxpad, nzpad, npml);
			read_d_ppmltb<<<dimGridpmltb,dimBlockpmltb>>>(d_ps, d_pspmltb, nxpad, nzpad, npml);	
			cuda_backward_p<<<dimGridp,dimBlock>>>(d_ps, d_vxs, d_vzs, d_rho, d_vp, d_diffcoef, _dtx, _dtz, npml, nxpad, nzpad);
			// insert source 
			cuda_add_source<<<dimGrid,dimBlock>>>(d_ps, d_source, d_wavelet, dt, 2, nxpad, nzpad, it);
			// receiver wavefield
			cuda_forward_vx<<<dimGridvx,dimBlock>>>(d_p, d_vx, d_rho, d_diffcoef, _dtx, npml, nxpad, nzpad);
			cuda_forward_vz<<<dimGridvz,dimBlock>>>(d_p, d_vz, d_rho, d_diffcoef, _dtz, npml, nxpad, nzpad);
			cuda_pml_vxlr<<<dimGridvxlr,dimBlock>>>(d_p, d_vx, d_rho, d_diffcoef, d_ddxVx, _dtx, dt, npml, nxpad, nzpad, d_norderx);
			cuda_pml_vztb<<<dimGridvztb,dimBlock>>>(d_p, d_vz, d_rho, d_diffcoef, d_ddzVz, _dtz, dt, npml, nxpad, nzpad, d_norderz);

			cuda_forward_p<<<dimGridp,dimBlock>>>(d_p, d_vx, d_vz, d_rho, d_vp, d_diffcoef, _dtx, _dtz, npml, nxpad, nzpad);
			cuda_pml_plr<<<dimGridplr,dimBlock>>>  (d_p, d_vx, d_vz, d_pl1, d_pl2, d_pr1, d_pr2, d_rho, d_vp, d_diffcoef, d_ddx, d_ddz,_dtx, _dtz, dt, npml, nxpad, nzpad, d_norder);
			cuda_pml_ptb<<<dimGridptb,dimBlock>>>  (d_p, d_vx, d_vz, d_pt1, d_pt2, d_pb1, d_pb2, d_rho, d_vp, d_diffcoef, d_ddz, _dtx, _dtz, dt, npml, nxpad, nzpad, d_norder);
			cuda_pml_pconner<<<dimGridpcooner,dimBlock>>>(d_p, d_vx, d_vz, d_pl1, d_pl2, d_pr1, d_pr2, d_rho, d_vp, d_diffcoef, d_ddx, d_ddz,_dtx, _dtz, dt, npml, nxpad, nzpad, d_norder);
			// insert source 
			cuda_insert_record2<<<(nxlength-noffset+127)/128,128>>>(d_p, &d_record[it*(nxlength-noffset)], npml, nxlength, noffset, dt);
			// imaging condition:
			cuda_cross_coorelation<<<dimGridp,dimBlock>>>(d_ps,d_p,d_g1,d_illum,nxpad,nzpad,npml);
			cuda_wavefield_decomposition<<<dimGridp,dimBlock>>>(d_ps,d_vxs,d_vzs,d_p,d_vx,d_vz,d_g2ud,d_g2du,d_g2lr,d_g2rl,nxpad,nzpad,npml);
			cuda_poynting<<<dimGridp,dimBlock>>>(d_ps,d_vxs,d_vzs,d_p,d_vx,d_vz,d_vp,d_g31,d_g32,nxpad,nzpad,npml);
			cuda_energynorm<<<dimGridp,dimBlock>>>(d_ps,d_ps_pre,d_p,d_p_pre,d_vp,d_g4,dx,dz,dt,nxpad,nzpad,npml);
		}
		// abtain g2 and g3
		cuda_stack_udlr<<<dimGridp,dimBlock>>>(d_g2ud,d_g2du,d_g2lr,d_g2rl,d_g2,nxlength,nz);
		cuda_stack_theta<<<dimGridp,dimBlock>>>(d_g31,d_g3,nxlength,nz);
		cuda_stack_theta<<<dimGridp,dimBlock>>>(d_g32,d_g3_true,nxlength,nz);
		// taper image
		cuda_taper_image<<<dimGridp,dimBlock>>>(d_g1, d_taper, nxlength, nz);
		cuda_taper_image<<<dimGridp,dimBlock>>>(d_g2, d_taper, nxlength, nz);
		cuda_taper_image<<<dimGridp,dimBlock>>>(d_g3, d_taper, nxlength, nz);
		cuda_taper_image<<<dimGridp,dimBlock>>>(d_g4, d_taper, nxlength, nz);
		cuda_taper_image<<<dimGridp,dimBlock>>>(d_g3_true, d_taper, nxlength, nz);
		cuda_taper_image<<<dimGridp,dimBlock>>>(d_illum, d_taper, nxlength, nz);
		// single-shot normalized			
		cuda_applyics<<<dimGridp,dimBlock>>>(d_image2s,d_g1,d_illum,nx,nxlength,nx1,nz);
		cuda_applyics<<<dimGridp,dimBlock>>>(d_image3s,d_g2,d_illum,nx,nxlength,nx1,nz);
		cuda_applyics<<<dimGridp,dimBlock>>>(d_image4s,d_g3,d_illum,nx,nxlength,nx1,nz);	
		cuda_applyics<<<dimGridp,dimBlock>>>(d_image5s,d_g4,d_illum,nx,nxlength,nx1,nz);	
		// multi-shot normalized
		cuda_applyic<<<dimGridp,dimBlock>>> (d_image1,  d_g1,nx,nxlength,nx1,nz);
		cuda_applyic<<<dimGridp,dimBlock>>> (d_image2m, d_g1,nx,nxlength,nx1,nz);
		cuda_applyic<<<dimGridp,dimBlock>>> (d_image3m, d_g2,nx,nxlength,nx1,nz);
		cuda_applyic<<<dimGridp,dimBlock>>> (d_image4m, d_g3,nx,nxlength,nx1,nz);
		cuda_applyic<<<dimGridp,dimBlock>>> (d_image5m, d_g4,nx,nxlength,nx1,nz);
		cuda_applyic<<<dimGridp,dimBlock>>> (d_imagetrue, d_g3_true,nx,nxlength,nx1,nz);
		cuda_applyic<<<dimGridp,dimBlock>>> (d_Illum,   d_illum,nx,nxlength,nx1,nz);
		// output temp image
		if ((ishot-1)%50 == 0)
		{		
			// single-shot normalized	
			// Image2s
			cudaMemcpy(image2s, d_image2s, nx*nz*sizeof(float), cudaMemcpyDeviceToHost);
			sprintf(buffrecord,"./output/Image2stemp.su");
			Output1d(image2s, nz, nx, dx, buffrecord, 1);
			// Image3s
			cudaMemcpy(image3s, d_image3s, nx*nz*sizeof(float), cudaMemcpyDeviceToHost);
			sprintf(buffrecord,"./output/Image3stemp.su");
			Output1d(image3s, nz, nx, dx, buffrecord, 1);
			// Image4s
			cudaMemcpy(image4s, d_image4s, nx*nz*sizeof(float), cudaMemcpyDeviceToHost);
			sprintf(buffrecord,"./output/Image4stemp.su");
			Output1d(image4s, nz, nx, dx, buffrecord, 1);
			// Image5s
			cudaMemcpy(image5s, d_image5s, nx*nz*sizeof(float), cudaMemcpyDeviceToHost);
			sprintf(buffrecord,"./output/Image5stemp.su");
			Output1d(image5s, nz, nx, dx, buffrecord, 1);		
			// multishot normlized
			// Image1
			cudaMemcpy(image1, d_image1, nx*nz*sizeof(float), cudaMemcpyDeviceToHost);
			sprintf(buffrecord,"./output/Image1temp.su");
			Output1d(image1, nz, nx, dx, buffrecord, 1);
			// Illum
			cudaMemcpy(illum,  d_Illum,  nx*nz*sizeof(float), cudaMemcpyDeviceToHost);
			sprintf(buffrecord,"./output/Illumtemp.su");
			Output1d(illum, nz, nx, dx, buffrecord, 1);
			// Image2m
			cudaMemcpy(image2m, d_image2m, nx*nz*sizeof(float), cudaMemcpyDeviceToHost);
			sprintf(buffrecord,"./output/Image2mtemp.su");
			Output1d(image2m, nz, nx, dx, buffrecord, 1);
			// Image3m
			cudaMemcpy(image3m, d_image3m, nx*nz*sizeof(float), cudaMemcpyDeviceToHost);
			sprintf(buffrecord,"./output/Image3mtemp.su");
			Output1d(image3m, nz, nx, dx, buffrecord, 1);
			// Image4m
			cudaMemcpy(image4m, d_image4m, nx*nz*sizeof(float), cudaMemcpyDeviceToHost);
			sprintf(buffrecord,"./output/Image4mtemp.su");
			Output1d(image4m, nz, nx, dx, buffrecord, 1);
			// Image5m
			cudaMemcpy(image5m, d_image5m, nx*nz*sizeof(float), cudaMemcpyDeviceToHost);
			sprintf(buffrecord,"./output/Image5mtemp.su");
			Output1d(image5m, nz, nx, dx, buffrecord, 1);
			// Imagetrue
			cudaMemcpy(imagetrue, d_imagetrue, nx*nz*sizeof(float), cudaMemcpyDeviceToHost);
			sprintf(buffrecord,"./output/Imagetruetemp.su");
			Output1d(imagetrue, nz, nx, dx, buffrecord, 1);	
		}
		free(record);
		free(vp);
		free(rho);   		   
		free(vxspmllr);
		free(vzspmltb);
		free(pspmllr);
		free(pspmltb);
		free(p);
		free(tr);

		cudaFree(d_diffcoef);cudaFree(d_record);cudaFree(d_source);
		cudaFree(d_vp);cudaFree(d_rho);
		cudaFree(d_p); cudaFree(d_vx); cudaFree(d_vz); cudaFree(d_p_pre);
		cudaFree(d_ps);cudaFree(d_vxs);cudaFree(d_vzs);cudaFree(d_ps_pre);
		cudaFree(d_pl1);cudaFree(d_pl2);cudaFree(d_pr1);cudaFree(d_pr2);
		cudaFree(d_pt1);cudaFree(d_pt2);cudaFree(d_pb1);cudaFree(d_pb2);

		cudaFree(d_ddx);cudaFree(d_ddz);cudaFree(d_ddxVx);cudaFree(d_ddzVz);
		cudaFree(d_norder);cudaFree(d_norderx);cudaFree(d_norderz);

		cudaFree(d_vxspmllr);cudaFree(d_vzspmltb);cudaFree(d_pspmllr);cudaFree(d_pspmltb);	

		cudaFree(d_g1);  cudaFree(d_g2);  cudaFree(d_g3);  cudaFree(d_g4);	
		cudaFree(d_g2ud);cudaFree(d_g2du);cudaFree(d_g2lr);cudaFree(d_g2rl);	
		cudaFree(d_g31); cudaFree(d_g32); cudaFree(d_g3_true);		
		cudaFree(d_illum);
		cudaFree(d_taper);

		cudaEventRecord(stop); 
		cudaEventSynchronize(stop); 
	   	cudaEventElapsedTime(&mstimer, start, stop); 

		printf("%d shot finished: %g (s)\n",ishot, mstimer*1.e-3); 
	}
	cudaMemcpy(imagetrue, d_imagetrue, nx*nz*sizeof(float), cudaMemcpyDeviceToHost);
	sprintf(buffrecord,"./output/%dImagetrue.su",idevice);
	Output1d(imagetrue, nz, nx, dx, buffrecord, 1);	
	
	cudaMemcpy(image2s, d_image2s, nx*nz*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(image3s, d_image3s, nx*nz*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(image4s, d_image4s, nx*nz*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(image5s, d_image5s, nx*nz*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(illum,   d_Illum,   nx*nz*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(image1,  d_image1,  nx*nz*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(image2m, d_image2m, nx*nz*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(image3m, d_image3m, nx*nz*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(image4m, d_image4m, nx*nz*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(image5m, d_image5m, nx*nz*sizeof(float), cudaMemcpyDeviceToHost);

	cudaFree(d_image1);
	cudaFree(d_image2s);cudaFree(d_image3s);cudaFree(d_image4s);cudaFree(d_image5s);
	cudaFree(d_image2m);cudaFree(d_image3m);cudaFree(d_image4m);cudaFree(d_image5m);
	cudaFree(d_Illum);
	cudaFree(d_imagetrue);
	free(imagetrue);
}
