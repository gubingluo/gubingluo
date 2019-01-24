#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sgyhead.h"
#include "common.h"


void swap_short_2(short *tni2)  
/************************************************************************** 
swap_short_2        swap a short integer 
***************************************************************************/  
{  
 *tni2=(((*tni2>>8)&0xff) | ((*tni2&0xff)<<8));    
}  
  
void swap_u_short_2(unsigned short *tni2)  
/************************************************************************** 
swap_u_short_2      swap an unsigned short integer 
***************************************************************************/  
{  
 *tni2=(((*tni2>>8)&0xff) | ((*tni2&0xff)<<8));    
}  
  
void swap_int_4(int *tni4)  
/************************************************************************** 
swap_int_4      swap a 4 byte integer 
***************************************************************************/  
{  
 *tni4=(((*tni4>>24)&0xff) | ((*tni4&0xff)<<24) |  
        ((*tni4>>8)&0xff00) | ((*tni4&0xff00)<<8));    
}  
  
void swap_u_int_4(unsigned int *tni4)  
/************************************************************************** 
swap_u_int_4        swap an unsigned integer 
***************************************************************************/  
{  
 *tni4=(((*tni4>>24)&0xff) | ((*tni4&0xff)<<24) |  
        ((*tni4>>8)&0xff00) | ((*tni4&0xff00)<<8));    
}  
  
void swap_long_4(long *tni4)  
/************************************************************************** 
swap_long_4     swap a long integer 
***************************************************************************/  
{  
 *tni4=(((*tni4>>24)&0xff) | ((*tni4&0xff)<<24) |  
        ((*tni4>>8)&0xff00) | ((*tni4&0xff00)<<8));    
}  
  
void swap_u_long_4(unsigned long *tni4)  
/************************************************************************** 
swap_u_long_4       swap an unsigned long integer 
***************************************************************************/  
{  
 *tni4=(((*tni4>>24)&0xff) | ((*tni4&0xff)<<24) |  
        ((*tni4>>8)&0xff00) | ((*tni4&0xff00)<<8));    
}  
  
void swap_float_4(float *tnf4)  
/************************************************************************** 
swap_float_4        swap a float 
***************************************************************************/  
{  
 int *tni4=(int *)tnf4;  
 *tni4=(((*tni4>>24)&0xff) | ((*tni4&0xff)<<24) |  
        ((*tni4>>8)&0xff00) | ((*tni4&0xff00)<<8));    
}

void float_to_float(float *trace, int nt, int endianflag)
{	
	if (endianflag){
	for (int it=0; it<nt; it++){
		swap_float_4(&trace[it]);}}
}
void int_to_int(int *trace, int nt, int endianflag)
{
	if (endianflag){
	for (int it=0; it<nt; it++)
		swap_int_4(&trace[it]);}
}
void short_to_short(short int *trace, int nt, int endianflag)
{
	if (endianflag){
	for (int it=0; it<nt; it++)
		swap_short_2(&trace[it]);}
}

void init_sgyhead(sgyhead *header)
{
	header[0].tracl = 0;	                                          //1-4             1 itraceAll                   read and write
	header[0].tracr = 0;	                                          //5-8
	header[0].fldr = 0;	                                          //9-12
	header[0].tracf = 0;                                            	//13-16           2 itrace at each shot         read and write
	header[0].ep = 0;	                                          	//17-20
	header[0].cdp = 0;	                                          //21-24
	header[0].cdpt = 0;	                                          //25-28
	header[0].trid = 0;	                                          //29-30
	header[0].nvs = 0;	                                          //31-32
	header[0].nhs = 0;	                                          //33-34
	header[0].duse = 0;	                                          //35-36
	header[0].offset = 0;	                                          //37-40           3  offset for each trace at each shot (change with gx)  read and write
	header[0].gelev = 0;	                                          //41-44
	header[0].selev = 0;	                                          //45-48
	header[0].sdepth = 0;	                                          //49-52
	header[0].gdel = 0;	                                          //53-56
	header[0].sdel = 0;	                                          //57-60
	header[0].swdep = 0;	                                          //61-64
	header[0].gwdep = 0;	                                          //65-68               
	header[0].scalel = 0;                                         	//69-70
	header[0].scalco = 0;                                         	//71-72
	header[0].sx = 0;                                              	//73-76           4  sx     real position at whole model        read and write
	header[0].sy = 0;                                              	//77-80
	header[0].gx = 0;	                                          	//81-84           5  gx     real position at whole model        read and write
	header[0].gy = 0;                                              	//85-88 
	header[0].counit = 0;                                         	//89-90
	header[0].wevel = 0;                                          	//91-92
	header[0].swevel = 0;	                                    	//93-94
	header[0].sut = 0;	                                          //95-96
	header[0].gut = 0;                                            	//97-98
	header[0].sstat = 0;                                          	//99-100
	header[0].gstat = 0;                                          	//101-102
	header[0].tstat = 0;                                          	//103-104
	header[0].laga = 0;	                                          //105-106
	header[0].lagb = 0;	                                          //107-108
	header[0].delrt = 0;                                          	//109-110
	header[0].muts = 0;	                                          //111-112
	header[0].mute = 0;	                                          //113-114
	header[0].ns = 0;                                    			//115-116         6   nt    number of sample                    read and write
	header[0].dt = 0;                                    			//117-118         7   dt    frequency of sample (us)            read and write
	header[0].gain = 0;	                                          //119-120
	header[0].igc = 0;                                            	//121-122
	header[0].igi = 0;                                            	//123-124
	header[0].corr = 0;                                           	//125-126
	header[0].sfs = 0;                                            	//127-128
	header[0].sfe = 0;                                            	//129-130
	header[0].slen = 0;	                                          //131-132
	header[0].styp = 0;                                           	//133-134
	header[0].stas = 0;	                                          //135-136
	header[0].stae = 0;	                                          //137-138
	header[0].tatyp = 0;                                          	//139-140
	header[0].afilf = 0;	                                    	//141-142
	header[0].afils = 0;	                                    	//143-144
	header[0].nofilf = 0;                                         	//145-146
	header[0].nofils = 0;                                         	//147-148
	header[0].lcf = 0;	                                          //149-150
	header[0].hcf = 0;                                            	//151-152
	header[0].lcs = 0;                                            	//153-154
	header[0].hcs = 0;	                                          //155-156
	header[0].year = 0;	                                          //157-158
	header[0].day = 0;	                                          //159-160
	header[0].hour = 0;	                                          //161-162
	header[0].minute = 0;                                         	//163-164
	header[0].sec = 0;	                                          //165-166
	header[0].timbas = 0;                                         	//167-168
	header[0].trwf = 0;                                           	//169-170
	header[0].grnors = 0;                                         	//171-172
	header[0].grnofr = 0;                                         	//173-174
	header[0].grnlof = 0;                                         	//175-176
	header[0].gaps = 0;	                                          //177-178
	header[0].otrav = 0;                                          	//179-180
	/* cwp local assignments */
	header[0].d1 = 0.0;                                             	//181-184
	header[0].f1 = 0.0;                                             	//185-188
	header[0].d2 = 0.0;	/* sample spacing between traces  byte# *///189-192 	      8   dx        write but not read		
	header[0].f2 = 0.0;	                                          //193-196
	header[0].ungpow = 0.0;                                         	//197-200
	header[0].unscale = 0.0;                                        	//201-204
	header[0].ntr = 0;                                              	//205-208         9   nx        write but not read
	header[0].mark = 0;	                                          //209-210
      header[0].shortpad = 0;                                       	//211-212
	header[0].unass[0] = 0;                                      	//213-240
	header[0].unass[1] = 0;                                      	//213-240
	header[0].unass[2] = 0;                                      	//213-240
	header[0].unass[3] = 0;                                      	//213-240
	header[0].unass[4] = 0;                                      	//213-240
	header[0].unass[5] = 0;                                      	//213-240
	header[0].unass[6] = 0;                                      	//213-240
	header[0].unass[7] = 0;                                      	//213-240
	header[0].unass[8] = 0;                                      	//213-240
	header[0].unass[9] = 0;                                      	//213-240
	header[0].unass[10] = 0;                                      	//213-240
	header[0].unass[11] = 0;                                      	//213-240
	header[0].unass[12] = 0;                                      	//213-240
	header[0].unass[13] = 0;                                      	//213-240

}
void swap_sgyhead(sgyhead *header)
{
	swap_int_4(&header[0].tracl);	                                          	//1-4             1 itraceAll                   read and write
	swap_int_4(&header[0].tracr);	                                          	//5-8
	swap_int_4(&header[0].fldr);	                                          	//9-12
	swap_int_4(&header[0].tracf);                                            	//13-16           2 itrace at each shot         read and write
	swap_int_4(&header[0].ep);	                                          	//17-20
	swap_int_4(&header[0].cdp);	                                          	//21-24
	swap_int_4(&header[0].cdpt);	                                         		//25-28
	swap_short_2(&header[0].trid);	                                          //29-30
	swap_short_2(&header[0].nvs);	                                          	//31-32
	swap_short_2(&header[0].nhs);	                                          	//33-34
	swap_short_2(&header[0].duse);	                                          //35-36
	swap_int_4(&header[0].offset);	                                          //37-40           3  offset for each trace at each shot (change with gx)  read and write
	swap_int_4(&header[0].gelev);	                                          	//41-44
	swap_int_4(&header[0].selev);	                                          	//45-48
	swap_int_4(&header[0].sdepth);	                                          //49-52
	swap_int_4(&header[0].gdel);	                                          	//53-56
	swap_int_4(&header[0].sdel);	                                          	//57-60
	swap_int_4(&header[0].swdep);	                                          	//61-64
	swap_int_4(&header[0].gwdep);	                                          	//65-68               
	swap_short_2(&header[0].scalel);                                         	//69-70
	swap_short_2(&header[0].scalco);                                         	//71-72
	swap_int_4(&header[0].sx);                                             		//73-76           4  sx     real position at whole model        read and write
	swap_int_4(&header[0].sy);                                              	//77-80
	swap_int_4(&header[0].gx);	                                          	//81-84           5  gx     real position at whole model        read and write
	swap_int_4(&header[0].gy);                                              	//85-88 
	swap_short_2(&header[0].counit);                                         	//89-90
	swap_short_2(&header[0].wevel);                                          	//91-92
	swap_short_2(&header[0].swevel);	                                    	//93-94
	swap_short_2(&header[0].sut);	                                          	//95-96
	swap_short_2(&header[0].gut);                                            	//97-98
	swap_short_2(&header[0].sstat);                                          	//99-100
	swap_short_2(&header[0].gstat);                                          	//101-102
	swap_short_2(&header[0].tstat);                                          	//103-104
	swap_short_2(&header[0].laga);	                                          //105-106
	swap_short_2(&header[0].lagb);	                                          //107-108
	swap_short_2(&header[0].delrt);                                          	//109-110
	swap_short_2(&header[0].muts);	                                          //111-112
	swap_short_2(&header[0].mute);	                                          //113-114
	swap_short_2(&header[0].ns);                                    		//115-116         6   nt    number of sample                    read and write
	swap_short_2(&header[0].dt);                                    		//117-118         7   dt    frequency of sample (us)            read and write
	swap_short_2(&header[0].gain);	                                          //119-120
	swap_short_2(&header[0].igc);                                            	//121-122
	swap_short_2(&header[0].igi);                                            	//123-124
	swap_short_2(&header[0].corr);                                           	//125-126
	swap_short_2(&header[0].sfs);                                            	//127-128
	swap_short_2(&header[0].sfe);                                            	//129-130
	swap_short_2(&header[0].slen);	                                          //131-132
	swap_short_2(&header[0].styp);                                           	//133-134
	swap_short_2(&header[0].stas);	                                          //135-136
	swap_short_2(&header[0].stae);	                                          //137-138
	swap_short_2(&header[0].tatyp);                                          	//139-140
	swap_short_2(&header[0].afilf);	                                    	//141-142
	swap_short_2(&header[0].afils);	                                    	//143-144
	swap_short_2(&header[0].nofilf);                                         	//145-146
	swap_short_2(&header[0].nofils);                                         	//147-148
	swap_short_2(&header[0].lcf);	                                          //149-150
	swap_short_2(&header[0].hcf);                                            	//151-152
	swap_short_2(&header[0].lcs);                                            	//153-154
	swap_short_2(&header[0].hcs);	                                          //155-156
	swap_short_2(&header[0].year);	                                          //157-158
	swap_short_2(&header[0].day);	                                          //159-160
	swap_short_2(&header[0].hour);	                                          //161-162
	swap_short_2(&header[0].minute);                                         	//163-164
	swap_short_2(&header[0].sec);	                                          //165-166
	swap_short_2(&header[0].timbas);                                         	//167-168
	swap_short_2(&header[0].trwf);                                           	//169-170
	swap_short_2(&header[0].grnors);                                         	//171-172
	swap_short_2(&header[0].grnofr);                                         	//173-174
	swap_short_2(&header[0].grnlof);                                         	//175-176
	swap_short_2(&header[0].gaps);	                                          //177-178
	swap_short_2(&header[0].otrav);                                          	//179-180
	/* cwp local assignments */
	swap_float_4(&header[0].d1);                                           //181-184
	swap_float_4(&header[0].f1);                                           //185-188
	swap_float_4(&header[0].d2);								//189-192 	      8   dx        write but not read		
	swap_float_4(&header[0].f2);	                                         	//193-196
	swap_float_4(&header[0].ungpow);                                       //197-200
	swap_float_4(&header[0].unscale);                                      //201-204
	swap_int_4(&header[0].ntr);                                              	//205-208         9   nx        write but not read
	swap_short_2(&header[0].mark);	                                          //209-210
      swap_short_2(&header[0].shortpad);                                       	//211-212
	swap_short_2(&header[0].unass[0]);                                      	//213-240
	swap_short_2(&header[0].unass[1]);                                      	//213-240
	swap_short_2(&header[0].unass[2]);                                      	//213-240
	swap_short_2(&header[0].unass[3]);                                      	//213-240
	swap_short_2(&header[0].unass[4]);                                      	//213-240
	swap_short_2(&header[0].unass[5]);                                      	//213-240
	swap_short_2(&header[0].unass[6]);                                      	//213-240
	swap_short_2(&header[0].unass[7]);                                      	//213-240
	swap_short_2(&header[0].unass[8]);                                      	//213-240
	swap_short_2(&header[0].unass[9]);                                      	//213-240
	swap_short_2(&header[0].unass[10]);                                      	//213-240
	swap_short_2(&header[0].unass[11]);                                      	//213-240
	swap_short_2(&header[0].unass[12]);                                      	//213-240
	swap_short_2(&header[0].unass[13]);                                      	//213-240

}
void swap_sgyheadkey(sgyheadkey *headerkey)
{
	swap_int_4(&headerkey[0].iTraceAll);
	swap_int_4(&headerkey[0].iTrace);
	swap_int_4(&headerkey[0].Offset);
	swap_int_4(&headerkey[0].sx);
	swap_int_4(&headerkey[0].gx);
	swap_short_2(&headerkey[0].Nt);
	swap_short_2(&headerkey[0].Dt);
	swap_float_4(&headerkey[0].Dx);
	swap_int_4(&headerkey[0].Ntrace);
}
void fread_sgyhead(sgyheadkey *headerkey, FILE *fp, int Endian)
{
	if (!Endian)
	{
		fread(&headerkey[0].iTraceAll,sizeof(int),1,fp);        //4
		fseek(fp,8L,1);                                         //12
		fread(&headerkey[0].iTrace,sizeof(int),1,fp);           //16
		fseek(fp,20L,1);                                        //36
		fread(&headerkey[0].Offset,sizeof(int),1,fp);           //40
		fseek(fp,32L,1);		                                //72
		fread(&headerkey[0].sx,sizeof(int),1,fp);               //76
		fseek(fp,4L,1);                                         //80
		fread(&headerkey[0].gx,sizeof(int),1,fp);               //84
		fseek(fp,30L,1);                                        //114
		fread(&headerkey[0].Nt,sizeof(short int),1,fp);         //116
		fread(&headerkey[0].Dt,sizeof(short int),1,fp);         //118
		fseek(fp,70L,1);                                        //188
		fread(&headerkey[0].Dx,sizeof(float),1,fp);             //192
		fseek(fp,12L,1);                                        //204
		fread(&headerkey[0].Ntrace,sizeof(int),1,fp);           //208
		fseek(fp,32L,1);                                        //240
	}
	else
	{
		fread(&headerkey[0].iTraceAll,sizeof(int),1,fp);
		swap_int_4(&headerkey[0].iTraceAll);
		fseek(fp,8L,1);
		fread(&headerkey[0].iTrace,sizeof(int),1,fp);
		swap_int_4(&headerkey[0].iTrace);
		fseek(fp,20L,1);
		fread(&headerkey[0].Offset,sizeof(int),1,fp);
		swap_int_4(&headerkey[0].Offset);
		fseek(fp,32L,1);
		fread(&headerkey[0].sx,sizeof(int),1,fp);
		swap_int_4(&headerkey[0].sx);
		fseek(fp,4L,1);
		fread(&headerkey[0].gx,sizeof(int),1,fp);
		swap_int_4(&headerkey[0].gx);
		fseek(fp,30L,1);
		fread(&headerkey[0].Nt,sizeof(short int),1,fp);
		swap_short_2(&headerkey[0].Nt);
		fread(&headerkey[0].Dt,sizeof(short int),1,fp);
		swap_short_2(&headerkey[0].Dt);
		fseek(fp,70L,1);
		fread(&headerkey[0].Dx,sizeof(float),1,fp);
		swap_float_4(&headerkey[0].Dx);
		fseek(fp,12L,1);
		fread(&headerkey[0].Ntrace,sizeof(int),1,fp);
		swap_int_4(&headerkey[0].Ntrace);
		fseek(fp,32L,1);
	}
}
void fwrite_sgyhead(sgyhead *header, sgyheadkey headerkey, FILE *fp, int Endian)
{
	if (!Endian)
	{
		header[0].tracl       =   headerkey.iTraceAll;
		header[0].tracf	    =   headerkey.iTrace;
		header[0].offset      =   headerkey.Offset;		
		header[0].sx          =   headerkey.sx;
		header[0].gx          =   headerkey.gx;
		header[0].ns          =   headerkey.Nt;
		header[0].dt          =   headerkey.Dt;
		header[0].d2          =   headerkey.Dx;
		header[0].ntr         =   headerkey.Ntrace;

		fwrite(header,sizeof(sgyhead),1,fp);
	}
	else
	{
		header[0].tracl       =   headerkey.iTraceAll;
		header[0].tracf	    =   headerkey.iTrace;
		header[0].offset      =   headerkey.Offset;		
		header[0].sx          =   headerkey.sx;
		header[0].gx          =   headerkey.gx;
		header[0].ns          =   headerkey.Nt;
		header[0].dt          =   headerkey.Dt;
		header[0].d2          =   headerkey.Dx;
		header[0].ntr         =   headerkey.Ntrace;

		swap_sgyhead(header);

		fwrite(header,sizeof(sgyhead),1,fp);
	}
}


