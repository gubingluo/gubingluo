#ifndef FUNCTION_H
#define FUNCTION_H

//=============================================================//
//=======================function.cpp==========================//
//=============================================================//
//================alloc and free matrix
float ***Creat3dArray(int m, int n, int k);
void free3dArray(float ***tt, int m, int n, int k);
float **Creat2dArray(int m, int n);
void free2dArray(float **tt, int m, int n);
int **Creat2dArray_int(int m, int n);
void free2dArray_int(int **tt, int m, int n);
//================initial matrix(one or two dimension)
void fmemset1(float *p, int len);
void fmemset2(float **p, int nz, int nx);
void fmemset3(float ***p, int nz, int nx, int ny);
void fmemset1v(float *p, int len, float v);
void fmemset2v(float **p, int nz, int nx, float v);
void fmemset1vp(float *p, int len, float *vp);
void fmemset2vp(float **p, int nz, int nx, float **vp);
//================max or min value of matrix
float absMaxval2_AB(float **A, float **B, int nz, int nx);
float Maxval1(float *v, int n);
float Minval1(float *v, int n);
float Maxval2(float **v, int nz, int nx);
float Minval2(float **v, int nz, int nx);
float absMaxval1(float *v, int n);
float absMinval1(float *v, int n);
float absMaxval2(float **v, int nz, int nx);
float absMinval2(float **v, int nz, int nx);
//================sum
float sum1(float *data, int n);
float sum2(float **data, int nx, int nz);
float sum1abs(float *data, int n);
float sum2abs(float **data, int nx, int nz);
void arrayabs(float *data1, float *data, int n);
//============================file operation
int exist(char *buff);
long filesize(FILE *fp);
//=============================================================//
//=======================velocity_function.cpp=================//
//=============================================================//
void pmlvelsmooth1d(float *vp, int nx, int nz, int npml);
void pmlvelsmooth2d(float **vp, int nx, int nz, int npml);
void velsmooth1d(float *vp,int n1,int n2,int nsp);
void velsmooth2d(float **vp,int n1,int n2,int nsp);
void gaussian1d_smoothing (int ns, int nsr, float *data);
void interpvel1(float *vp0,  float *vp,  int nx0, int nz0, int nx, int nz, int lx, int lz);
void interpvel2(float **vp0, float **vp, int nx0, int nz0, int nx, int nz, int lx, int lz);
void extendvel1(float *vp0, float *vp, int nx, int nz, int npml);
void extendvel2(float **vp0, float **vp, int nx, int nz, int npml);
void extractvel1(float *temp, float *vp, int nx, int nz, int nx1, int nx2);
void extractrho1(float *temp, float *rho, int nx, int nz, int nx1, int nx2);
void extractvel2(float **temp, float *vp, int nx, int nz, int nx1, int nx2);
void extractrho2(float **temp, float *rho, int nx, int nz, int nx1, int nx2);

//=============================================================//
//=======================IO_function.cpp=======================//
//=============================================================//
void ReadVelMod(char velfile[40], int nx ,int nz ,float *vp);
void Output1d(float *record,  int nt, int nx, float dt, char buff[40], int Out_flag);
void Output2d(float **record, int nt, int nx, float dt, char buff[40], int Out_flag);
void Input1d(float *record,  int nt, int nx, char buff[40],int In_flag);
void Input2d(float **record, int nt, int nx, char buff[40],int In_flag);
//=============================================================//
//=======================Interp_function.cpp===================//
//=============================================================//
double modZeroBessel(double x);
void Kaiser(double *wkaiser, int windowLength, double beta);
float Sinc(float x);
float sum_ABC(float *A, float *B, float *C, int length);
void Sincinterp(float **seisin, float **seisout, int nx, int nt, int nt0, float dt, float dt0);
void GeneralFKinterp(float **seisin, float **seisout, int nx, int nx0, int nt0, int lx);
void Interpseis2d(float **seiscal, float **seisobs, int nx, int nx0, int nt, int nt0, float dt, float dt0);
//=============================================================//
//=======================FD_function.cpp===================//
//=============================================================//
void Ricker(float *wavelet, float fdom, int nt, float dt, int Flag);
void Diff_coeffVelStr(float *data, int N_order);
void Diff_coeffVelStrLiuyang(float *data, int N_order);

//=============================================================//
//=======================image_function.cpp====================//
//=============================================================//
void scaleillum1d(float *illum, float *grad, int n);
void scaleillum2d(float **illum, float **grad, int nz, int nx, int npml);
void scaleimage1d(float *illum, float *grad, int n);
void scaleimage2d(float **illum, float **grad, int nz, int nx);
void illumsmooth1d(float *Illum, int nx, int nz, int nsr);
void illumsmooth2d(float **Illum, int nx, int nz, int nsr);
void laplacefd1d(int nx, int nz, int dx, int dz, float *image);
void laplacefd2d(int nx, int nz, int dx, int dz, float **image);
void taperh1d(int nx, int nz, int nwinlen, float gd, float *g);
void taperh2d(int nx, int nz, int nwinlen, float gd, float **g);

#endif
