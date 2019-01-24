#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include "fftw3.h"
#include "function.h"
#include "common.h"

// sinc interpolate functions
double modZeroBessel(double x)
{
	int i;
	double x_2 = x/2;
	double num = 1;
	double fact = 1;
	double result = 1;

	for (i=1 ; i<20 ; i++) {
		num *= x_2 * x_2;
		fact *= i;
		result += num / (fact * fact);
	}
	return result;
}
void Kaiser(float *wkaiser, int windowLength, double beta)
{
	double m_2 = (double)(windowLength-1)/2.0;
	double denom = modZeroBessel(beta);					// Denominator of Kaiser function

	if (wkaiser == NULL) {
		wkaiser = (float *) malloc(windowLength * sizeof(float));
		if (wkaiser == NULL) {
			fprintf(stderr, "Could not allocate memory for window\n");
		}
	}
	int n;
	for (n=0; n<windowLength; n++)
	{
		double val = ((n) - m_2) / m_2;
		val = 1 - (val * val);
		wkaiser[n] = modZeroBessel(beta * sqrt(val)) / denom;
	}
}
float Sinc(float x)
{
	float sinc;
	if (x == 0.0)
		sinc = cos(x);
	else
		sinc = sin(pi*x)/(pi*x);
	return sinc;
}
float sum_ABC(float *A, float *B, float *C, int length)
{
	float sum=0.0;
	for (int il=0; il<length; il++)
		sum = sum + A[il]*B[il]*C[il];
	return sum;
}
void Sincinterp(float **seisin, float **seisout, int nx, int nt, int nt0, float dt, float dt0)
{
	// windowed sinc interpolation
	// seisin:  input record  size [nt0][nx]
	// seisout: output record size [nt][nx]
	int ix,it,ik,halfnw,
		it0,itstart,itend,lenit,inwstart,inwend;
	float nw,beta,t,xtemp;
	nw = 21;
	beta = 2.0;
	halfnw = (nw-1)/2;
	float *tn=(float *)malloc(nt0*sizeof(float));
	float *wkaiser=(float *)malloc(nw*sizeof(float));
	float *tr,*wpart,*sinc;
	for (it=0;it<nt0;it++)
		tn[it] = it*dt0;

	Kaiser(wkaiser, nw, beta);

	for (ix=0;ix<nx;ix++)
	{
		for (it=0;it<nt;it++)
		{
			it0 = (int)((it+1)*dt/dt0);
			itstart = MAX(0,    it0-halfnw);
			itend   = MIN(nt0-1,it0+halfnw);

			lenit = itend - itstart + 1;

			inwstart = (itstart > halfnw)*0 + (itstart <= halfnw)*(nw - lenit);
			inwend   = (itstart <= nt0 - nw)*(nw-1) + (itstart > nt0 - nw)*(nt0 - itstart);
			
			t = it*dt;

			tr   = (float *)malloc(lenit*sizeof(float));
			sinc = (float *)malloc(lenit*sizeof(float));
			wpart= (float *)malloc(lenit*sizeof(float));

			for (ik=0; ik<lenit; ik++){
				tr[ik] = seisin[ik+itstart][ix];
				wpart[ik] = wkaiser[ik+inwstart];
				
				
				xtemp = (t-tn[ik+itstart])/dt0;

				sinc[ik] = Sinc(xtemp);}

			seisout[it][ix] = sum_ABC(tr, sinc, wpart, lenit);

			free(tr);
			free(sinc);
			free(wpart);
		}		
	}
	free(tn);
	free(wkaiser);
}
// GFKI interpolate functions
float fftw_maxfabs(fftw_complex *in, int length)
{
	float result=0.0;
	int il;
	for (il=0;il<length;il++)
		result = MAX(result, fftw_fabs(in[il]));
	return result;
}
void GeneralFKinterp(float **seisin, float **seisout, int nx, int nx0, int nt0, int lx)
{
	// General F-K interpolation
	// seisin:  input record  size [nt0][nx0]
	// seisout: output record size [nt0][nx]
	int it,it1,ix,il;
	int nf = 1 + (int)(nt0/2);
	int Nx = lx*nx0,
		Nt = lx*nt0,
		nx_nt = nx0*nt0,
		Nx_nt = Nx*nt0,
		Nx_Nt = Nx*Nt,
		Nx_nf = Nx*nf;
	float zmax;
	fftw_complex *B,*K,*C,*S,*Z,*H,*G1,*G;
	fftw_plan p1,p2,p3;
	B = (fftw_complex *)fftw_malloc(nx_nt*sizeof(fftw_complex));
	K = (fftw_complex *)fftw_malloc(nx_nt*sizeof(fftw_complex));
	C = (fftw_complex *)fftw_malloc(Nx_nt*sizeof(fftw_complex));
	S = (fftw_complex *)fftw_malloc(Nx_Nt*sizeof(fftw_complex));
	Z = (fftw_complex *)fftw_malloc(Nx_Nt*sizeof(fftw_complex));
	H = (fftw_complex *)fftw_malloc(Nx_nf*sizeof(fftw_complex));
	G1 = (fftw_complex *)fftw_malloc(Nx_nf*sizeof(fftw_complex));
	G = (fftw_complex *)fftw_malloc(Nx_nt*sizeof(fftw_complex));
	// step 1)
	#pragma omp parallel for default(shared) private(it,ix)
	for (ix=0;ix<nx0;ix++)
		for (it=0;it<nt0;it++){
			B[ix*nt0+it][0] = seisin[it][ix];
			B[ix*nt0+it][1] = 0.0;}

	p1 = fftw_plan_dft_2d(nx0,nt0,B,K,FFTW_FORWARD,FFTW_ESTIMATE);
	fftw_execute(p1);

	#pragma omp parallel for default(shared) private(it,ix,il)
	for (il=0; il<lx; il++)
		for (ix=0; ix<nx0; ix++)
			for (it=0; it<nt0; it++){
				C[(ix+il*nx0)*nt0+it][0] = K[ix*nt0+it][0];
				C[(ix+il*nx0)*nt0+it][1] = K[ix*nt0+it][1];}
	// step 2)
	memset(S,0,Nx_Nt*sizeof(fftw_complex));
	memset(Z,0,Nx_Nt*sizeof(fftw_complex));
	#pragma omp parallel for default(shared) private(it,ix)
	for (ix=0;ix<nx0;ix++)
		for (it=0;it<nt0;it++){
			S[ix*Nt+it][0] = B[ix*nt0+it][0];
			S[ix*Nt+it][1] = B[ix*nt0+it][1];}
	p2 = fftw_plan_dft_2d(Nx,Nt,S,S,FFTW_FORWARD,FFTW_ESTIMATE);
	fftw_execute(p2);
	#pragma omp parallel for default(shared) private(it,ix,il)
	for (il=0;il<lx;il++)
		for (ix=0; ix<nx0; ix++)
			for (it=0; it<Nt; it++){
				Z[ix*Nt+it][0] = Z[ix*Nt+it][0] + S[(ix+il*nx0)*Nt+it][0];
				Z[ix*Nt+it][1] = Z[ix*Nt+it][1] + S[(ix+il*nx0)*Nt+it][1];}

	#pragma omp parallel for default(shared) private(it,ix)
	for (ix=0; ix<nx0; ix++)
		for (it=0; it<Nt; it++){
			Z[ix*Nt+it][0] = Z[ix*Nt+it][0]/lx;
			Z[ix*Nt+it][1] = Z[ix*Nt+it][1]/lx;}
	#pragma omp parallel for default(shared) private(it,ix,il)
	for (il=1;il<lx;il++)
		for (ix=0; ix<nx0; ix++)
			for (it=0; it<Nt; it++){
				Z[(ix+il*nx0)*Nt+it][0] = Z[ix*Nt+it][0];
				Z[(ix+il*nx0)*Nt+it][1] = Z[ix*Nt+it][1];}
	zmax = fftw_maxfabs(Z, Nx_Nt);
	#pragma omp parallel for default(shared) private(it,ix)
	for (ix=0; ix<Nx; ix++)
		for (it=0; it<Nt; it++){
			Z[ix*Nt+it][0] = Z[ix*Nt+it][0] + zmax*0.00001;}
	#pragma omp parallel for default(shared) private(it,ix)
	for (ix=0; ix<Nx; ix++)
		for (it=0; it<nf; it++){
			H[ix*nf+it][0] = fftw_divi_real(S[ix*Nt+it],Z[ix*Nt+it]);
			H[ix*nf+it][1] = fftw_divi_imag(S[ix*Nt+it],Z[ix*Nt+it]);
			if (fftw_fabs(H[ix*nf+it]) > lx*lx){
				H[ix*nf+it][0] = (float)(lx);
				H[ix*nf+it][1] = 0.0;}}
	// step 3)	
	for (it=0; it<nf; it++){
		#pragma omp parallel for default(shared) private(ix)
		for (ix=0; ix<Nx; ix++){
			G1[it*Nx+ix][0] = fftw_multip_real(H[ix*nf+it],C[ix*nt0+it]);
			G1[it*Nx+ix][1] = fftw_multip_imag(H[ix*nf+it],C[ix*nt0+it]);}
		p3 = fftw_plan_dft_1d(Nx,&G1[it*Nx],&G1[it*Nx],FFTW_BACKWARD,FFTW_ESTIMATE);
		fftw_execute(p3);}
	#pragma omp parallel for default(shared) private(it,ix)
	for (ix=0; ix<Nx; ix++)
		for (it=0; it<nf; it++){
			G[ix*nt0+it][0] = G1[it*Nx+ix][0];
			G[ix*nt0+it][1] = G1[it*Nx+ix][1];}
	for (ix=0; ix<Nx; ix++){
		for (it=nf,it1=(int)((nt0+1)/2); it<nt0; it++,it1--){
			G[ix*nt0+it][0] =  G1[it1*Nx+ix][0];
			G[ix*nt0+it][1] = -G1[it1*Nx+ix][1];}
		p3 = fftw_plan_dft_1d(nt0,&G[ix*nt0],&G[ix*nt0],FFTW_BACKWARD,FFTW_ESTIMATE);
		fftw_execute(p3);}

	for (ix=0;ix<nx;ix++)
		for (it=0;it<nt0;it++)
			seisout[it][ix] = fftw_real(G[ix*nt0+it]);

	fftw_execute(p1); 
	fftw_execute(p2); 
	fftw_execute(p3); 
	fftw_free(B);
	fftw_free(K);
	fftw_free(C);
	fftw_free(S);
	fftw_free(Z);
	fftw_free(H);
	fftw_free(G1);
	fftw_free(G);
}
// interpolate seismic profile
void Interpseis2d(float **seiscal, float **seisobs, int nx, int nx0, int nt, int nt0, float dt, float dt0)
{
	float **seistemp;
	int lx = (int)(nx-1)/(nx0-1);
	
	if (nx > nx0 && nt > nt0){
		// x- and t- interpolate
		seistemp = Creat2dArray(nt0,nx);
		GeneralFKinterp(seisobs, seistemp, nx, nx0, nt0, lx);
		Sincinterp(seistemp, seiscal, nx, nt, nt0, dt, dt0);
		free2dArray(seistemp,nt0,nx);}
	else if (nx > nx0 && nt <= nt0){
		// x-interpolate
		GeneralFKinterp(seisobs, seiscal, nx, nx0, nt0, lx);}
	else if (nx <= nx0 && nt > nt0){
		// t-interpolate
		Sincinterp(seisobs, seiscal, nx, nt, nt0, dt, dt0);}
	else{
		seiscal = seisobs;}
}
