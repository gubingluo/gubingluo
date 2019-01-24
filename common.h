#ifndef Esc
#define Esc   0x1b
#endif
#ifndef Enter
#define Enter 0x0d
#endif

//#define MPICH_SKIP_MPICXX

#ifndef pi
#define pi 3.1415926535897932
#endif

#ifndef eps
#define eps 1.0e-6
#endif

#ifndef MAX
#define	MAX(x,y) ((x) > (y) ? (x) : (y))
#endif
#ifndef MIN
#define	MIN(x,y) ((x) < (y) ? (x) : (y))
#endif

#ifndef fftw_fabs
#define	fftw_fabs(x) (sqrtf(x[0]*x[0] + x[1]*x[1]))
#endif

#ifndef fftw_real
#define	fftw_real(x) (x[0])
#endif

#ifndef fftw_multip_real
#define fftw_multip_real(x,y) (x[0]*y[0] - x[1]*y[1])
#endif

#ifndef fftw_multip_imag
#define fftw_multip_imag(x,y) (x[0]*y[1] + x[1]*y[0])
#endif

#ifndef fftw_divi_real
#define fftw_divi_real(x,y) ((x[0]*y[0] + x[1]*y[1])/(y[0]*y[0] + y[1]*y[1]))
#endif

#ifndef fftw_divi_imag
#define fftw_divi_imag(x,y) ((x[1]*y[0] - x[0]*y[1])/(y[0]*y[0] + y[1]*y[1]))
#endif


#ifndef sign
#define sign(x) (x>=0)?1:-1
#endif

#define Block_Sizex 32
#define Block_Sizez 32

#define N 20
