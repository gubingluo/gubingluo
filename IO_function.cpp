#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "function.h"
#include "sgyhead.h"
#include "common.h"


void ReadVelMod(char velfile[40], int nx ,int nz ,float *vp)
{
	int ix,iz;
	float **vel;
	char cmd[100];
	sprintf(cmd,"./input/%s",velfile);
	
	vel=Creat2dArray(nx,nz);
	FILE *fp=NULL;
	
	fp=fopen(cmd,"rb");
	if (fp==NULL)
	{
		printf("The file %s open failed !\n",velfile);
	}
	for (ix=0;ix<nx;ix++)
	{
		fseek(fp,240L,1);
		fread(vel[ix],sizeof(float),nz,fp);
	}
	fclose(fp);

	for (iz=0;iz<nz;iz++)
		for (ix=0;ix<nx;ix++)
			vp[iz*nx+ix] = vel[ix][iz];
	free2dArray(vel,nx,nz);
}
// I/O functions
void Outseisfile1d_single(float *seiscal, sgyheadkey headerkey, int outflag, int endianflag, char buff[40])
{
	// outflag:    0 - seiscal[nt][nx]; 1 - seiscal[nx][nt]
	// sgyflag:    0 - su;              1: sgy
	// endianflag: 0 - little endian;   1: big endian
	int it,ix;
	int iTraceAll,iTrace,offset,gx;
	short int nt;
	int nx;
	float dx;
	sgyhead header;
	FILE *fp = NULL;

	nt = headerkey.Nt;
	nx = headerkey.Ntrace;

	float *temp = (float *)malloc(nt*sizeof(float ));
	fp=fopen(buff,"ab+");
	if (fp==NULL){printf("The file %s open failed !\n",buff);}

	iTraceAll = headerkey.iTraceAll;
	iTrace = headerkey.iTrace;
	offset = headerkey.Offset;
	gx = headerkey.gx;
	dx = headerkey.Dx;

	for (ix=0;ix<nx;ix++){
		init_sgyhead(&header);
		headerkey.iTraceAll = iTraceAll + ix;
		headerkey.iTrace = iTrace + ix;
		headerkey.Offset = offset + (int)(ix*dx);
		headerkey.gx     = gx     + (int)(ix*dx);	
		headerkey.Dx     = dx;
		headerkey.Ntrace = nx;	

		fwrite_sgyhead(&header, headerkey, fp, endianflag);
		if (!outflag){
			for (it=0; it<nt; it++)
				temp[it] = seiscal[it*nx+ix];}
		else{
			for (it=0; it<nt; it++)
				temp[it] = seiscal[ix*nt+it];}
		float_to_float(temp, nt, endianflag);
		fwrite(temp,sizeof(float),nt,fp);}
	fclose(fp);
	free(temp);
}
void Outseisfile2d_single(float **seiscal, sgyheadkey headerkey, int outflag, int endianflag, char buff[40])
{
	// outflag:    0 - seiscal[nt][nx]; 1 - seiscal[nx][nt]
	// sgyflag:    0 - su;              1: sgy
	// endianflag: 0 - little endian;   1: big endian
	int it,ix;
	int iTraceAll,iTrace,offset,gx;
	short int nt;
	int nx;
	float dx;
	sgyhead header;
	FILE *fp = NULL;

	nt = headerkey.Nt;
	nx = headerkey.Ntrace;

	float *temp = (float *)malloc(nt*sizeof(float ));
	fp=fopen(buff,"ab+");
	if (fp==NULL){printf("The file %s open failed !\n",buff);}

	iTraceAll = headerkey.iTraceAll;
	iTrace = headerkey.iTrace;
	offset = headerkey.Offset;
	gx = headerkey.gx;
	dx = headerkey.Dx;

	for (ix=0;ix<nx;ix++){
		init_sgyhead(&header);
		headerkey.iTraceAll = iTraceAll + ix;
		headerkey.iTrace = iTrace + ix;
		headerkey.Offset = offset + (int)(ix*dx);
		headerkey.gx     = gx     + (int)(ix*dx);	
		headerkey.Dx     = dx;
		headerkey.Ntrace = nx;	

		fwrite_sgyhead(&header, headerkey, fp, endianflag);
		if (!outflag){
			for (it=0; it<nt; it++)
				temp[it] = seiscal[it][ix];}
		else{
			for (it=0; it<nt; it++)
				temp[it] = seiscal[ix][it];}
		float_to_float(temp, nt, endianflag);
		fwrite(temp,sizeof(float),nt,fp);}
	fclose(fp);
	free(temp);
}
void Outseisfile1d(float *seiscal, sgyheadkey headerkey, int outflag, int sgyflag, int endianflag, char buff[40])
{
	// outflag:    0 - seiscal[nt][nx]; 1 - seiscal[nx][nt]
	// sgyflag:    0 - su;              1: sgy
	// endianflag: 0 - little endian;   1: big endian
	int it,ix;
	int iTraceAll,iTrace,offset,gx;
	short int nt;
	int nx;
	float dx;
	sgyhead header;
	FILE *fp = NULL;

	nt = headerkey.Nt;
	nx = headerkey.Ntrace;

	float *temp = (float *)malloc(nt*sizeof(float ));
	fp=fopen(buff,"wb");
	if (fp==NULL){printf("The file %s open failed !\n",buff);}
	if (sgyflag){
		short int volume[1800];
		memset(volume,0,1800*sizeof(short int));
		short_to_short(volume,1800,endianflag);
		fwrite(volume, sizeof(short int), 1800, fp);}

	iTraceAll = headerkey.iTraceAll;
	iTrace = headerkey.iTrace;
	offset = headerkey.Offset;
	gx = headerkey.gx;
	dx = headerkey.Dx;

	for (ix=0;ix<nx;ix++){
		init_sgyhead(&header);
		headerkey.iTraceAll = iTraceAll + ix;
		headerkey.iTrace = iTrace + ix;
		headerkey.Offset = offset + (int)(ix*dx);
		headerkey.gx     = gx     + (int)(ix*dx);	
		headerkey.Dx     = dx;
		headerkey.Ntrace = nx;	

		fwrite_sgyhead(&header, headerkey, fp, endianflag);
		if (!outflag){
			for (it=0; it<nt; it++)
				temp[it] = seiscal[it*nx+ix];}
		else{
			for (it=0; it<nt; it++)
				temp[it] = seiscal[ix*nt+it];}
		float_to_float(temp, nt, endianflag);
		fwrite(temp,sizeof(float),nt,fp);}
	fclose(fp);
	free(temp);
}
void Outseisfile2d(float **seiscal, sgyheadkey headerkey, int outflag, int sgyflag, int endianflag, char buff[40])
{
	// outflag:    0 - seiscal[nt][nx]; 1 - seiscal[nx][nt]
	// sgyflag:    0 - su;              1: sgy
	// endianflag: 0 - little endian;   1: big endian
	int it,ix;
	int iTraceAll,iTrace,offset,gx;
	short int nt;
	int nx;
	float dx;
	sgyhead header;
	FILE *fp = NULL;

	nt = headerkey.Nt;
	nx = headerkey.Ntrace;

	float *temp = (float *)malloc(nt*sizeof(float ));
	fp=fopen(buff,"wb");
	if (fp==NULL){printf("The file %s open failed !\n",buff);}
	if (sgyflag){
		short int volume[1800];
		memset(volume,0,1800*sizeof(short int));
		short_to_short(volume,1800,endianflag);
		fwrite(volume, sizeof(short int), 1800, fp);}

	iTraceAll = headerkey.iTraceAll;
	iTrace = headerkey.iTrace;
	offset = headerkey.Offset;
	gx = headerkey.gx;
	dx = headerkey.Dx;

	for (ix=0;ix<nx;ix++){
		init_sgyhead(&header);
		headerkey.iTraceAll = iTraceAll + ix;
		headerkey.iTrace = iTrace + ix;
		headerkey.Offset = offset + (int)(ix*dx);
		headerkey.gx     = gx     + (int)(ix*dx);	
		headerkey.Dx     = dx;
		headerkey.Ntrace = nx;	

		fwrite_sgyhead(&header, headerkey, fp, endianflag);
		if (!outflag){
			for (it=0; it<nt; it++)
				temp[it] = seiscal[it][ix];}
		else{
			for (it=0; it<nt; it++)
				temp[it] = seiscal[ix][it];}
		float_to_float(temp, nt, endianflag);
		fwrite(temp,sizeof(float),nt,fp);}
	fclose(fp);
	free(temp);
}
void get_parameter(char buff[40], int sgyflag, int *ntrace, short int *nsample, short int *dt, int *endianflag)
{
	sgyheadkey headerkey;
	short int volume[1800];
	int iskip;
	FILE *fp=NULL;
	// assume the file is little endian su file
	*endianflag = 0;

	fp=fopen(buff,"rb");
	if (sgyflag)
		fread(volume,sizeof(short int), 1800, fp);
	fread_sgyhead(&headerkey, fp, *endianflag);
	fclose(fp);

	*dt = headerkey.Dt;
	*nsample = headerkey.Nt;

	if (*dt <= 0 || *nsample <= 0){
		swap_sgyheadkey(&headerkey);
		
		*dt = headerkey.Dt;
		*nsample = headerkey.Nt;
		if (*dt > 0 && *nsample > 0){
			*endianflag = 1;}}
	else{
		*endianflag = 0;}
	fp=fopen(buff,"rb");	
	if (sgyflag){
		fseek(fp,3600L,1);}
	fread_sgyhead(&headerkey, fp, *endianflag);
	fclose(fp);
	*dt = headerkey.Dt;
	*nsample = headerkey.Nt;	
	*ntrace = 0;
	float *temp =(float *)malloc((*nsample)*sizeof(float));
	fp=fopen(buff,"rb");
	if (sgyflag)
		fseek(fp,3600L,1);

	while(!feof(fp)){
		fseek(fp,240L,1);
		fread(temp, sizeof(float), *nsample, fp);
		*ntrace = *ntrace + 1;}
	fclose(fp);
	free(temp);
	*ntrace = *ntrace - 1;
	
}
void Inseisfile1d(float *seisobs, sgyheadkey *headerkey, int inflag, int sgyflag, char buff[40])
{
	// inflag:     0 - seisabs[nt][nx]; 1 - seisabs[nx][nt]
	// sgyflag:    0 - su;              1 - sgy
	// endianflag: 0 - little endian;   1 - big endian
	int it,ix;
	int nx;
	short nt;
	short dt;
	int endianflag;
	FILE *fp = NULL;

	get_parameter(buff,sgyflag,&nx,&nt,&dt,&endianflag);

	float *temp =(float *)malloc(nt*sizeof(float));
	fp=fopen(buff,"rb");
	if (sgyflag)
		fseek(fp,3600L,1);
	fread_sgyhead(headerkey, fp, endianflag);
	fclose(fp);
	
	headerkey[0].Dt = dt;
	headerkey[0].Nt = nt;
	headerkey[0].Ntrace = nx;

	fp=fopen(buff,"rb");
	if (sgyflag)
		fseek(fp,3600L,1);
	for (ix=0;ix<nx;ix++){		
		fseek(fp,240L,1);
		fread(temp, sizeof(float), nt, fp);
		float_to_float(temp, nt, endianflag);
		if (!inflag)
			for (it=0; it<nt; it++)
				seisobs[it*nx+ix] = temp[it];
		else
			for (it=0; it<nt; it++)
				seisobs[ix*nt+it] = temp[it];}
	fclose(fp);
	free(temp);
}
void Inseisfile2d(float **seisobs, sgyheadkey *headerkey, int inflag, int sgyflag, char buff[40])
{
	// inflag:     0 - seisabs[nt][nx]; 1 - seisabs[nx][nt]
	// sgyflag:    0 - su;              1 - sgy
	// endianflag: 0 - little endian;   1 - big endian
	int it,ix;
	int nx;
	short nt;
	short dt;
	int endianflag;
	FILE *fp = NULL;

	get_parameter(buff,sgyflag,&nx,&nt,&dt,&endianflag);

	float *temp =(float *)malloc(nt*sizeof(float));
	fp=fopen(buff,"rb");
	if (sgyflag)
		fseek(fp,3600L,1);
	fread_sgyhead(headerkey, fp, endianflag);
	fclose(fp);
	
	headerkey[0].Dt = dt;
	headerkey[0].Nt = nt;
	headerkey[0].Ntrace = nx;

	fp=fopen(buff,"rb");
	if (sgyflag)
		fseek(fp,3600L,1);
	for (ix=0;ix<nx;ix++){		
		fseek(fp,240L,1);
		fread(temp, sizeof(float), nt, fp);
		float_to_float(temp, nt, endianflag);
		if (!inflag)
			for (it=0; it<nt; it++)
				seisobs[it][ix] = temp[it];
		else
			for (it=0; it<nt; it++)
				seisobs[ix][it] = temp[it];}
	fclose(fp);
	free(temp);
}
void Output1d(float *record, int nt, int nx, float dt, char buff[40], int Out_flag)
{
	int it,ix;
	FILE *fp=NULL;
	if (Out_flag==1)
	{
		float *temp =(float *)malloc(nt*sizeof(float));
		fp=fopen(buff,"wb");
		if (fp==NULL)
		{
			printf("The file %s open failed !\n",buff);
		}
		short int header[120];
		for (it=0;it<120;it++)
		{
			header[it]=0;
		}
		header[57]=(short int)(nt);
		header[58]=(short int)(dt*1000000.0);         // dt
		header[104]=(short int)(nx);
		for (ix=0;ix<nx;ix++)
		{
			header[0]=ix+1;
			fwrite(header,2,120,fp);

			for (it=0; it<nt; it++)
				temp[it] = record[it*nx+ix];

			fwrite(temp,sizeof(float),nt,fp);
		}
		fclose(fp);
		free(temp);
	}
	else
	{
		fp=fopen(buff,"wb");
		if (fp==NULL)
		{
			printf("The file %s open failed !\n",buff);
		}
		short int header[120];
		for (it=0;it<120;it++)
		{
			header[it]=0;
		}
		header[57]=(short int)(nt);
		header[58]=(short int)(dt*1000000.0);         // dt
		header[104]=(short int)(nx);
		for (ix=0;ix<nx;ix++)
		{
			header[0]=ix+1;
			fwrite(header,2,120,fp);
			fwrite(&record[ix*nt],sizeof(float),nt,fp);
		}
		fclose(fp);
	} 
}
void Input1d(float *record, int nt, int nx, char buff[40],int In_flag)
{
	int ix,it;
	FILE *fp=NULL;
	if (In_flag==1)
	{
		float *temp =(float *)malloc(nt*sizeof(float));
		fp=fopen(buff,"rb");
		if (fp==NULL)
		{
			printf("The file %s open failed !\n",buff);
		}
		for (ix=0;ix<nx;ix++)
		{
			fseek(fp,240L,1);
			fread(temp,sizeof(float),nt,fp);
			for (it=0; it <nt; it++)
				record[it*nx+ix] = temp[it];
		}
		fclose(fp);
		free(temp);
	} 
	else
	{
		fp=fopen(buff,"rb");
		if (fp==NULL)
		{
			printf("The file %s open failed !\n",buff);
		}
		for (ix=0;ix<nx;ix++)
		{
			fseek(fp,240L,1);
			fread(&record[ix*nt],sizeof(float),nt,fp);
		}
		fclose(fp);
	}	
}

void Output2d(float **record, int nt, int nx, float dt, char buff[40], int Out_flag)
{
	int it,ix;
	FILE *fp=NULL;
	if (Out_flag==1)
	{
		float *temp =(float *)malloc(nt*sizeof(float));
		fp=fopen(buff,"wb");
		if (fp==NULL)
		{
			printf("The file %s open failed !\n",buff);
		}
		short int header[120];
		for (it=0;it<120;it++)
		{
			header[it]=0;
		}
		header[57]=(short int)(nt);
		header[58]=(short int)(dt*1000000.0);         // dt
		header[104]=(short int)(nx);
		for (ix=0;ix<nx;ix++)
		{
			header[0]=ix+1;
			fwrite(header,2,120,fp);

			for (it=0; it<nt; it++)
				temp[it] = record[it][ix];

			fwrite(temp,sizeof(float),nt,fp);
		}
		fclose(fp);
		free(temp);
	}
	else
	{
		fp=fopen(buff,"wb");
		if (fp==NULL)
		{
			printf("The file %s open failed !\n",buff);
		}
		short int header[120];
		for (it=0;it<120;it++)
		{
			header[it]=0;
		}
		header[57]=(short int)(nt);
		header[58]=(short int)(dt*1000000.0);         // dt
		header[104]=(short int)(nx);
		for (ix=0;ix<nx;ix++)
		{
			header[0]=ix+1;
			fwrite(header,2,120,fp);
			fwrite(&record[ix],sizeof(float),nt,fp);
		}
		fclose(fp);
	} 
}
void Input2d(float **record, int nt, int nx, char buff[40],int In_flag)
{
	int ix,it;
	FILE *fp=NULL;
	if (In_flag==1)
	{
		float *temp =(float *)malloc(nt*sizeof(float));
		fp=fopen(buff,"rb");
		if (fp==NULL)
		{
			printf("The file %s open failed !\n",buff);
		}
		for (ix=0;ix<nx;ix++)
		{
			fseek(fp,240L,1);
			fread(temp,sizeof(float),nt,fp);
			for (it=0; it <nt; it++)
				record[it][ix] = temp[it];
		}
		fclose(fp);
		free(temp);
	} 
	else
	{
		fp=fopen(buff,"rb");
		if (fp==NULL)
		{
			printf("The file %s open failed !\n",buff);
		}
		for (ix=0;ix<nx;ix++)
		{
			fseek(fp,240L,1);
			fread(&record[ix],sizeof(float),nt,fp);
		}
		fclose(fp);
	}	
}
