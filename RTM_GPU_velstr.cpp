//////////////////////////////////////////////////////////////////////////
//      2D finite difference time domain acoustic wave 
// two-order displacement wave equation forward simulation
//   multi-shots for least square reverse time migration
//////////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <time.h>

#include "function.h"
#include "sgyhead.h"
#include "common.h"

void get_constant(float *vp, float fdom,
		      float dx0, float dz0, int nx0, int nz0, int nt0, float dt0,
			      float *dx, float *dz, int *nx, int *nz, int *nt, float *dt,  int *lx, int *lz, float *vpmin, float *vpmax);
void scanshotgather(char filename_shotgather[40], int **shotprofileId, int sgyflag, int *endianflag);


extern "C" void rtm_model(int idevice, int nshot, int startshot, int dshot, int mode, int medium_flag, int flag_layout, int sgyflag, int endianflag, int ismth,
	                    int nx, int nz, int nt, float dx, float dz, float dt, int npml,  float direct, float direct0, 
	                    float tsnap, float fdom, float spx0, float spz0, float dspx, float offset, float offsetmax, float Amp, float alp,
	                    float *vp0, float *rho0, float *diffcoef, float *wavelet, 
				  float *image1,  float *image2s, float *image3s, float *image4s, float *image5s, 
		 		  float *image2m, float *image3m, float *image4m, float *image5m, float *Illum);

extern "C" void rtm_real(int idevice, int nshot, int startshot, int dshot, int medium_flag, int flag_layout, int sgyflag, int endianflag, int ismth,
	                   int nx, int nz, int nt, float dx, float dz, float dt, int nx0, int nz0, int nt0, float dx0, float dz0, float dt0,
		             int npml, int lx, int lz, float tsnap, float fdom, float spz0, float Amp, float alp, int nw, int tlength, 
		             float *vp, float *rho, float *diffcoef, float *wavelet, char filename_shotgather[40], int **shotprofileId,
				 float *image1,  float *image2s, float *image3s, float *image4s, float *image5s, 
		 		 float *image2m, float *image3m, float *image4m, float *image5m, float *Illum);

int main(int argc, char* argv[])
{
	float dt,dx,dz,xmax,zmax;
	float fdom,Amp,alp,spx0,spz0,dspx;
	float dt0,dx0,dz0;
	float vpmax,vpmin;
	float offset,offsetmax;
	float tsnap;
	float direct,direct0;

	int nt0,nx0,nz0;
	int lx,lz;
	int nt,nx,nz,startshot,nshot,dshot,nshotall,npml,ismth,idevice,
	    medium_flag,flag_layout,dataflag,mode,sgyflag,endianflag,rhoflag;

	char buff[100],vp_file[40],rho_file[40];
	char cmd[40],cmd1[100];
	char filename_shotgather[40];
	int i,ishot,ix,iz;
	clock_t start,stop;

	FILE *fp=NULL;	

	if ((fp=fopen("./input/cmd.in","r"))==NULL)
	{
		printf("cmd.in is not exist !\n");
		exit(0);
	}
	fgets(buff,100,fp);
	fscanf(fp,"%s", cmd);
	fclose(fp);

	sprintf(cmd1,"./input/%s",cmd);

	if ((fp=fopen(cmd1,"r"))==NULL)
	{
		printf("%s is not exist !\n",cmd);
		exit(0);
	}
	fgets(buff,100,fp);
	fscanf(fp,"%f", &dt0);fscanf(fp,"%d",&nt0);fscanf(fp,"%f",&tsnap);
	fgets(buff,100,fp);fgets(buff,100,fp);
	fscanf(fp,"%f", &xmax);fscanf(fp,"%f", &zmax);fscanf(fp,"%f", &dx0);fscanf(fp,"%f", &dz0);			
	fgets(buff,100,fp);fgets(buff,100,fp);
	fscanf(fp,"%f", &fdom);fscanf(fp,"%f", &Amp);fscanf(fp,"%f", &alp); 			
	fgets(buff,100,fp);fgets(buff,100,fp);
	fscanf(fp,"%f", &spx0);fscanf(fp,"%f", &spz0);fscanf(fp,"%f", &dspx);
	fgets(buff,100,fp);fgets(buff,100,fp);
	fscanf(fp,"%f", &offset); fscanf(fp,"%f", &offsetmax);				 
	fgets(buff,100,fp);fgets(buff,100,fp);
	fscanf(fp,"%d", &startshot);fscanf(fp,"%d", &nshot);fscanf(fp,"%d", &dshot);fscanf(fp,"%d", &nshotall);
	fgets(buff,100,fp);fgets(buff,100,fp);
	fscanf(fp,"%d", &npml);
	fgets(buff,100,fp);fgets(buff,100,fp);
	fscanf(fp,"%d", &medium_flag);fscanf(fp,"%d", &ismth);
	fgets(buff,100,fp);fgets(buff,100,fp);
	fscanf(fp,"%f", &direct);fscanf(fp,"%f", &direct0); 
	fgets(buff,100,fp);fgets(buff,100,fp);
	fscanf(fp,"%d", &idevice); 
	fgets(buff,100,fp);fgets(buff,100,fp);
	fscanf(fp,"%d", &flag_layout);
	fgets(buff,100,fp);fgets(buff,100,fp);
	fscanf(fp,"%d", &mode);
	fgets(buff,100,fp);fgets(buff,100,fp);
	fscanf(fp,"%d", &dataflag);
	fgets(buff,100,fp);fgets(buff,100,fp);
	fscanf(fp,"%d", &sgyflag);
	fgets(buff,100,fp);fgets(buff,100,fp);
	fscanf(fp,"%d", &endianflag);
	fgets(buff,100,fp);fgets(buff,100,fp);
	fscanf(fp,"%d", &rhoflag);
	fgets(buff,100,fp);fgets(buff,100,fp);
	fscanf(fp,"%s", vp_file);
	if (dataflag == 1)
	{
		fgets(buff,100,fp);fgets(buff,100,fp);
		fscanf(fp,"%s", filename_shotgather);
	}
	if (rhoflag == 1)
	{
		fgets(buff,100,fp);fgets(buff,100,fp);
		fscanf(fp,"%s", rho_file);
	}
	fclose(fp);
	nx0=(int)(xmax/dx0)+1;
	nz0=(int)(zmax/dz0)+1;
	if (spx0<0 || spz0 <0 || spx0+(nshot-1)*dspx>(nx0-1)*dx0 ||spx0+(nshot-1)*dspx<0)
	{
		printf("the shot position out of the model !\nplease check the parameter !\n");
		exit(0);
	}
	// optimal
	float *vp0=(float *)malloc((nx0*nz0)*sizeof(float));
	ReadVelMod(vp_file,nx0,nz0,vp0);
	get_constant(vp0, fdom, dx0, dz0, nx0, nz0, nt0, dt0,
			 &dx, &dz,  &nx, &nz, &nt, &dt, &lx, &lz, &vpmin, &vpmax); 

	// interpolate model (velocity)
	float *vp=(float *)malloc((nx*nz)*sizeof(float)),
		*rho=(float *)malloc((nx*nz)*sizeof(float));

	interpvel1(vp0, vp, nx0, nz0, nx, nz, lx, lz);
	if (rhoflag){
		ReadVelMod(rho_file,nx0,nz0,vp0);
		int rhomax;
		rhomax = Maxval1(vp0, nx0*nz0);
		if (rhomax < 10.0)
			for (int ix=0;ix<nx0*nz0;ix++)
				vp0[ix] *= 1000.0;
		interpvel1(vp0, rho, nx0, nz0, nx, nz, lx, lz);}
	else
		fmemset1v(rho, nz*nx, 2000.0);
	free(vp0);
	
	//======================================================================================//
	//===================main program=======================================================//
	//======================================================================================//	

	printf("=====================================================\n");
	printf("=====================================================\n");
	printf("Prestack Reverse Time Migration by Bingluo Gu from UPC.\n");
	printf("            version 4.0  date 2015.11.25               \n");
	printf("=====================================================\n");
	printf("=====================================================\n");
	printf("================Required parameters:=================\n");
	printf(" velfile=   %s\n",vp_file);
	printf(" nshot  =   %d\t",nshot);printf(" fshot  =   %d\t",startshot);printf(" dshot  =   %d\n",dshot);
	printf(" nz     =   %d\t",nz0);  printf(" nx     =   %d\t",nx0);      printf(" nt     =   %d\n",nt0);
	printf(" dz     =   %.1f\t",dz0);printf(" dx     =   %.1f\t",dx0);    printf(" dt     =   %.4f\n",dt0);
	printf(" t      =   %f\n",nt0*dt0);
	printf(" vmin   =   %.2f\n",vpmin);
	printf(" vmax   =   %.2f\n",vpmax);
	printf(" fdom      =  %f  dominant frequency\n",fdom);
	printf(" npml      =  %2d  pml absorb layers number\n",npml);
	printf(" Norder    =  %2d  order of finite difference\n",N);
	printf(" dataflag  =  %2d  0: model data:                 1: real data\n",dataflag);
	printf(" mode      =  %2d  0: generate seismic profile;   1: seismic profile exist\n",mode);
	printf(" format:      %2d  0: su;                         1: sgy\n",sgyflag);
	printf(" Endian:      %2d  0: little edition;             1: big edition\n",endianflag);
	printf(" layout:      %2d  1: Bilateral;  2: full array;  3: Unilateral(r);  4: Unilateral(l)\n",flag_layout);
	printf("=================Optional parameters:================\n");
	printf(" nz     =   %d\t",nz);  printf(" nx   =   %d\t",nx);  printf(" nt   =   %d\n",nt);
	printf(" dz     =   %.1f\t",dz);printf(" dx   =   %.1f\t",dx);printf(" dt   =   %.4f\n",dt);
	printf(" lz     =   %d\t",lz);  printf(" lx   =   %d\n",lx);
	printf(" t      =   %f\n",nt*dt);
	printf("=======================================================\n");
	// prepare difference coefficients
	float *diffcoef=(float *)malloc((N/2)*(N/2)*sizeof(float));
	memset(diffcoef,0,(N/2)*(N/2)*sizeof(float));
	for (ix=0;ix<N/2;ix++)
	{
		int N1 = 2*(ix+1);
		float *diff2temp = (float *)malloc((N1/2)*sizeof(float));
		Diff_coeffVelStrLiuyang(diff2temp, N1);

		for (iz=0;iz<N1/2;iz++)
		{
			diffcoef[ix*N/2 + iz] = diff2temp[iz];
		}
	}
	// prepare wavelet
	float *wavelet=(float *)malloc(nt*sizeof(float));
	Ricker(wavelet, fdom, nt, dt, 1);
	// main function used to migration
	// alloc image
	float *image2m=(float *)malloc(nx*nz*sizeof(float)),
		*image3m=(float *)malloc(nx*nz*sizeof(float)),
		*image4m=(float *)malloc(nx*nz*sizeof(float)),
		*image5m=(float *)malloc(nx*nz*sizeof(float)),
		*image1 =(float *)malloc(nx*nz*sizeof(float)),
		*image2s=(float *)malloc(nx*nz*sizeof(float)),
		*image3s=(float *)malloc(nx*nz*sizeof(float)),
		*image4s=(float *)malloc(nx*nz*sizeof(float)),
		*image5s=(float *)malloc(nx*nz*sizeof(float)),
		*Illum =(float *)malloc(nx*nz*sizeof(float));
	int **shotprofileId;
	memset(image1, 0,nx*nz*sizeof(float));
	memset(image2s,0,nx*nz*sizeof(float));
	memset(image3s,0,nx*nz*sizeof(float));
	memset(image4s,0,nx*nz*sizeof(float));
	memset(image5s,0,nx*nz*sizeof(float));
	memset(image2m,0,nx*nz*sizeof(float));
	memset(image3m,0,nx*nz*sizeof(float));
	memset(image4m,0,nx*nz*sizeof(float));
	memset(image5m,0,nx*nz*sizeof(float));
	memset(Illum, 0,nx*nz*sizeof(float));

	shotprofileId = Creat2dArray_int(4,nshotall);

	if (!dataflag)
	{
		// use model data to implement RTM
		rtm_model(idevice, nshot, startshot, dshot, mode, medium_flag, flag_layout, sgyflag, endianflag, ismth,
			nx, nz, nt, dx, dz, dt, npml, direct, direct0, 
			tsnap, fdom, spx0, spz0, dspx, offset, offsetmax, Amp, alp,
			vp, rho, diffcoef, wavelet, image1, image2s, image3s, image4s, image5s, image2m, image3m, image4m, image5m, Illum);
	}
	else 
	{
		// scan shot gather to get starttrace number and nx of each shot
		sprintf(buff,"shotprofileId.dat");
		if (!exist(buff)){
			start = clock();
			scanshotgather(filename_shotgather, shotprofileId, sgyflag, &endianflag);
			fp = fopen("shotprofileId.dat","w");
			for (ix=0;ix<nshotall;ix++)
				fprintf(fp,"%d %d %d %d\n",shotprofileId[0][ix],shotprofileId[1][ix],shotprofileId[2][ix],shotprofileId[3][ix]);
			fclose(fp);
			stop = clock();
			printf("The time-elapse of scan shot gather: %f ms\n\n",(double)(stop-start)/CLOCKS_PER_SEC);}
		else{
			fp = fopen("shotprofileId.dat","r");
			for (ix=0;ix<nshotall;ix++)
				fscanf(fp,"%d %d %d %d",&shotprofileId[0][ix],&shotprofileId[1][ix],&shotprofileId[2][ix],&shotprofileId[3][ix]);
			fclose(fp);}			
		// use real data to implement RTM
		rtm_real(idevice, nshot, startshot, dshot, medium_flag, flag_layout, sgyflag, endianflag,ismth,
			nx, nz, nt, dx, dz, dt, nx0, nz0, nt0, dx0, dz0, dt0, 
			npml, lx, lz, tsnap, fdom, spz0, Amp, alp,direct, direct0, 
			vp, rho, diffcoef, wavelet, filename_shotgather, shotprofileId, 
			image1, image2s, image3s, image4s, image5s, image2m, image3m, image4m, image5m, Illum);
	}
	// output image
	char buffimage[40];	
	//===========================================================================================================
	// multi-shot normalized
	sprintf(buffimage,"./output/%dIllum.su",idevice);
	Output1d(Illum,  nz, nx, dz, buffimage, 1);
	// Image2: normalized cross-coorelation
	sprintf(buffimage,"./output/%dImage2m.su",idevice);
//	scaleimage1d(Illum, image2m, nz*nx);
	Output1d(image2m,  nz, nx, dz, buffimage, 1);
	// Image3: wavefield decomposed
	sprintf(buffimage,"./output/%dImage3m.su",idevice);
//	scaleimage1d(Illum, image3m, nz*nx);
	Output1d(image3m,  nz, nx, dz, buffimage, 1);
	// Image4: poynting vector
	sprintf(buffimage,"./output/%dImage4m.su",idevice);
//	scaleimage1d(Illum, image4m, nz*nx);
	Output1d(image4m,  nz, nx, dz, buffimage, 1);
	// Image5: energy norm
	sprintf(buffimage,"./output/%dImage5m.su",idevice);
//	scaleimage1d(Illum, image5m, nz*nx);
	Output1d(image5m,  nz, nx, dz, buffimage, 1);
	//=============================================================================================================
	// single-shot normalized
	// Image1: cross-coorelation
	sprintf(buffimage,"./output/%dImage1.su",idevice);	
	Output1d(image1,  nz, nx, dz, buffimage, 1);
	// Image2: normalized cross-coorelation
	sprintf(buffimage,"./output/%dImage2s.su",idevice);
	Output1d(image2s,  nz, nx, dz, buffimage, 1);
	// Image3: wavefield decomposed
	sprintf(buffimage,"./output/%dImage3s.su",idevice);
	Output1d(image3s,  nz, nx, dz, buffimage, 1);
	// Image4: poynting vector
	sprintf(buffimage,"./output/%dImage4s.su",idevice);
	Output1d(image4s,  nz, nx, dz, buffimage, 1);
	// Image5: energy norm
	sprintf(buffimage,"./output/%dImage5s.su",idevice);
	Output1d(image5s,  nz, nx, dz, buffimage, 1);

	free(vp);
	free(rho);
	free(wavelet);
	free(image1);
	free(image2s);
	free(image3s);
	free(image4s);
	free(image5s);
	free(image2m);
	free(image3m);
	free(image4m);
	free(image5m);
	free(Illum);
	free(diffcoef);
	free2dArray_int(shotprofileId, 4, nshot);

	printf("The program calculate over !\n");
	return 0;
}

void get_constant(float *vp, float fdom, 
			float dx0, float dz0, int nx0, int nz0, int nt0, float dt0,
			float *dx, float *dz, int *nx, int *nz, int *nt, float *dt,  int *lx, int *lz, float *vpmin, float *vpmax)            
{
    	int i,j;
	float hmin;
	//==================initial new parameter======================//
	*dx = dx0;
	*dz = dz0;
	*nx = nx0;
	*nz = nz0;
	*nt = nt0;
	*dt = dt0;
	*lx = 1;
	*lz = 1;
	*vpmax = Maxval1(vp, (*nz)*(*nx));
	*vpmin = Minval1(vp, (*nz)*(*nx));
	/*************** determine mininum spatial sampling interval */
	       
	/****** determine time sampling interval to ensure stability====*/
	float dx_max=(*vpmin)/fdom/7.25;  //dx_max=(*vpmin)/fdom*0.1;
    	float dz_max=dx_max;  

	//==================calculate new parameter===================//
    	while(dx_max < *dx){
		*dx = *dx/2;
		*lx = (*lx)*2;
		*nx = (*nx - 1)*2 + 1;}
    	while(dz_max < *dz){
		*dz = *dz/2;
		*lz = (*lz)*2;
		*nz = (*nz - 1)*2 + 1;}  
	//============ determine time sampling interval to ensure stability ===//
	hmin = ((*dx)<(*dz))?(*dx):(*dz);  
	float *diffcoef = (float *)malloc(N/2*sizeof(float));
	Diff_coeffVelStrLiuyang(diffcoef, N);
    	float dt_max=hmin/(sqrtf(2.0)*sum1abs(diffcoef, N/2)*(*vpmax));
	free(diffcoef);
	float t;
	t = nt0*dt0;
	if (dt_max < dt0){
		*dt = (int)(dt_max*10000.0)*0.0001;
		*nt=(int)(t/(*dt));}
}
void scanshotgather(char filename_shotgather[40], int **shotprofileId, int sgyflag, int *endianflag)
{
	// shotprofileId[5][nshot]
	// 1 行: sx             shotprofileId[0][ishot]
	// 2 行: ntrace         shotprofileId[1][ishot]
	// 3 行: itraceAll      shotprofileId[2][ishot]
	// 4 行: offset         shotprofileId[3][ishot]

	char buff[100];
	short int dt,nt;
	int iskip;
	sgyheadkey headerkey;

	FILE *fp = NULL,*fp1 = NULL;
	sprintf(buff,"./seisReal/%s",filename_shotgather);

	// check Endian
	*endianflag = 0;
	fp=fopen(buff,"rb"); if(fp == NULL) printf("please check the file or file folder name\n");

	if (sgyflag)
		fseek(fp,3600L,1);
	fread_sgyhead(&headerkey, fp, *endianflag);
	fclose(fp);

	dt = headerkey.Dt;
	nt = headerkey.Nt;

	if (dt <= 0 || nt <= 0){
		swap_sgyheadkey(&headerkey);
		
		dt = headerkey.Dt;
		nt = headerkey.Nt;
		if (dt > 0 && nt > 0){
			*endianflag = 1;}}
	else{
		*endianflag = 0;}

	printf("nt %d\n",nt);
	printf("dt %d\n",dt);
	// initialize out data
	int ishot, 
	    itraceAll, itraceAll_old, 
	    ntrace;	
	long length,lengthcur;

	iskip = nt*sizeof(float);

	fp=fopen(buff,"rb");
	if (sgyflag)
		fseek(fp,3600L,1);
	fread_sgyhead(&headerkey, fp, *endianflag);
	fclose(fp);

	ishot = 0;
	itraceAll = 1;
	itraceAll_old = 1;
	ntrace = itraceAll - itraceAll_old;	

	shotprofileId[0][0] = headerkey.sx;
	shotprofileId[1][0] = ntrace;
	shotprofileId[2][0] = itraceAll;
	shotprofileId[3][0] = headerkey.Offset;

	fp=fopen(buff,"rb");
	length = filesize(fp);
	lengthcur = sgyflag*3600;
	if (sgyflag)
		fseek(fp,3600L,1);
	do
	{				
		fread_sgyhead(&headerkey, fp, *endianflag);
		fseek(fp,iskip,1);
		lengthcur = lengthcur + iskip + 240;
		if (shotprofileId[0][ishot] != headerkey.sx)
		{
			ishot = ishot + 1;
			ntrace = itraceAll - itraceAll_old;

			shotprofileId[0][ishot] = headerkey.sx;
			shotprofileId[1][ishot - 1] = ntrace;
			shotprofileId[2][ishot] = itraceAll;
			shotprofileId[3][ishot] = headerkey.Offset;

			itraceAll_old = itraceAll;
		}
		itraceAll = itraceAll + 1;
	}while(!feof(fp) && headerkey.Nt > 0 && headerkey.Dt > 0 && lengthcur <= length);
	ntrace = itraceAll - itraceAll_old - 1;
	shotprofileId[1][ishot] = ntrace;
	fclose(fp);	
}
