#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "function.h"
#include "common.h"

void Ricker(float *wavelet, float fdom, int nt, float dt, int Flag)
{
	
	float t, t0,i1;
	int i;
	for (i=0;i<nt;i++)
	{
		wavelet[i]=0.0;
	}
	if (Flag==1)         // Ricker
	{
		t0=1.0/fdom;
		for (i=0;i<nt;i++)
		{
			i1=(float)(i);
			t=pi*fdom*(i1*dt-t0);
			wavelet[i]=(1.0-2.0*powf(t,2.0))*expf(-powf(t,2.0));
		}
	}
	else if (Flag==2)    // integral(Ricker)
	{
		for (i=0;i<nt;i++)
		{
			i1=(float)(i);
			t=pi*fdom*(i1*dt-1.0/fdom);
			wavelet[i]=-2.0*t*expf(-powf(t,2.0));
		}
	}
}

void Diff_coeffVelStr(float *data, int N_order)
{
	int N1=N_order/2;
	int m,i;
	if (N1 == 1)
	{
		data[0] = 1.0;
	}
	else if (N1 == 2)
	{
		data[0] = 9.0/8.0;
		data[1] = -1.0/24.0;
	}
	else
	{
		for (m=1; m<=N1; m++)
		{
			float a,b;
			a = 1.0;
			b = 1.0;
			for (i=1; i<=N1; i++)
			{
				if (i == m)
					continue;
				else
				{
					a = a*(2.0*i-1.0)*(2.0*i-1.0);
					b = b*fabs((2.0*m-1.0)*(2.0*m-1.0)-(2.0*i-1.0)*(2.0*i-1.0));
				}
			}
			a = a*powf(-1.0,m+1.0);
			b = b*(2.0*m-1.0);
			data[m-1] = a/b;
		}
	}
}
void Diff_coeffVelStrLiuyang(float *data, int N_order)
{
	int N1=N_order/2;
	int m,i;

	if (N1 > 10)
		N1 = 10;
	if (N1 == 1)
	{
		data[0] = 1.0;
	}
	else if (N1 == 2)
	{
		data[0] = 0.1129042e+1;
		data[1] = -0.4301412e-1;
	}
	else if (N1 == 3)
	{
		data[0] = 0.1185991e+1;
		data[1] = -0.7249965e-1;
		data[2] = 0.6301572e-2;
	}
	else if (N1 == 4)
	{
		data[0] = 0.1217990e+1;
		data[1] = -0.9382142e-1;
		data[2] = 0.1507536e-1;
		data[3] = -0.1700324e-2;
	}
	else if (N1 == 5)
	{
		data[0] = 0.1236607e+1;
		data[1] = -0.1082265e+0;
		data[2] = 0.2343440e-1;
		data[3] = -0.5033546e-2;
		data[4] = 0.6817483e-3;
	}
	else if (N1 == 6)
	{
		data[0] = 0.1247662e+1;
		data[1] = -0.1175538e+0;
		data[2] = 0.2997970e-1;
		data[3] = -0.8719078e-2;
		data[4] = 0.2215897e-2;
		data[5] = -0.3462075e-3;
	}
	else if (N1 == 7)
	{
		data[0] = 0.1254799E+1;
		data[1] = -0.1238928E+0;
		data[2] = 0.3494371E-1;
		data[3] = -0.1208897E-1;
		data[4] = 0.4132531E-2;
		data[5] = -0.1197110E-2;
		data[6] = 0.2122227E-3;
	}
	else if (N1 == 8)
	{
		data[0] = 0.1259312E+1;
		data[1] = -0.1280347E+0;
		data[2] = 0.3841945E-1;
		data[3] = -0.1473229E-1;
		data[4] = 0.5924913E-2;
		data[5] = -0.2248618E-2;
		data[6] = 0.7179226E-3;
		data[7] = -0.1400855E-3;
	}
	else if (N1 == 9)
	{
		data[0] = 0.1262502E+1;
		data[1] = -0.1310244E+0;
		data[2] = 0.4103928E-1;
		data[3] = -0.1686807E-1;
		data[4] = 0.7530520E-2;
		data[5] = -0.3345071E-2;
		data[6] = 0.1380367E-2;
		data[7] = -0.4808410E-3;
		data[8] = 0.1023759E-3;
	}
	else if (N1 == 10)
	{
		data[0] = 0.1264748E+1;
		data[1] = -0.1331606E+0;
		data[2] = 0.4296909E-1;
		data[3] = -0.1851897E-1;
		data[4] = 0.8861071E-2;
		data[5] = -0.4347073E-2;
		data[6] = 0.2076101E-2;
		data[7] = -0.9164925E-3;
		data[8] = 0.3437446E-3;
		data[9] = -0.7874250E-4;
	}
}

