#ifndef SGYHEAD_H
#define SGYHEAD_H

typedef struct SgyHead
{
	int tracl;	                                          //1-4             1 itraceAll                   read and write
	int tracr;	                                          //5-8
	int fldr;	                                          //9-12
	int tracf;                                            //13-16           2 itrace at each shot         read and write
	int ep;	                                          //17-20
	int cdp;	                                          //21-24
	int cdpt;	                                          //25-28
	short trid;	                                          //29-30
	short nvs;	                                          //31-32
	short nhs;	                                          //33-34
	short duse;	                                          //35-36
	int offset;	                                          //37-40           3  offset for each trace at each shot (change with gx)  read and write
	int gelev;	                                          //41-44
	int selev;	                                          //45-48
	int sdepth;	                                          //49-52
	int gdel;	                                          //53-56
	int sdel;	                                          //57-60
	int swdep;	                                          //61-64
	int gwdep;	                                          //65-68               
	short scalel;                                         //69-70
	short scalco;                                         //71-72
	int  sx;                                              //73-76           4  sx     real position at whole model        read and write
	int  sy;                                              //77-80
	int  gx;	                                          //81-84           5  gx     real position at whole model        read and write
	int  gy;                                              //85-88 
	short counit;                                         //89-90
	short wevel;                                          //91-92
	short swevel;	                                    //93-94
	short sut;	                                          //95-96
	short gut;                                            //97-98
	short sstat;                                          //99-100
	short gstat;                                          //101-102
	short tstat;                                          //103-104
	short laga;	                                          //105-106
	short lagb;	                                          //107-108
	short delrt;                                          //109-110
	short muts;	                                          //111-112
	short mute;	                                          //113-114
	short ns;                                    //115-116         6   nt    number of sample                    read and write
	short dt;                                    //117-118         7   dt    frequency of sample (us)            read and write
	short gain;	                                          //119-120
	short igc;                                            //121-122
	short igi;                                            //123-124
	short corr;                                           //125-126
	short sfs;                                            //127-128
	short sfe;                                            //129-130
	short slen;	                                          //131-132
	short styp;                                           //133-134
	short stas;	                                          //135-136
	short stae;	                                          //137-138
	short tatyp;                                          //139-140
	short afilf;	                                    //141-142
	short afils;	                                    //143-144
	short nofilf;                                         //145-146
	short nofils;                                         //147-148
	short lcf;	                                          //149-150
	short hcf;                                            //151-152
	short lcs;                                            //153-154
	short hcs;	                                          //155-156
	short year;	                                          //157-158
	short day;	                                          //159-160
	short hour;	                                          //161-162
	short minute;                                         //163-164
	short sec;	                                          //165-166
	short timbas;                                         //167-168
	short trwf;                                           //169-170
	short grnors;                                         //171-172
	short grnofr;                                         //173-174
	short grnlof;                                         //175-176
	short gaps;	                                          //177-178
	short otrav;                                          //179-180
	/* cwp local assignments */
	float d1;                                             //181-184
	float f1;                                             //185-188
	float d2;	/* sample spacing between traces  byte# *///189-192 	      8   dx        write but not read		
	float f2;	                                          //193-196
	float ungpow;                                         //197-200
	float unscale;                                        //201-204
	int ntr;                                              //205-208         9   nx        write but not read
	short mark;	                                          //209-210
      short shortpad;                                       //211-212
	short unass[14];                                      //213-240

}sgyhead;


typedef struct SgyHeadKey
{
	int iTraceAll;                                //0            1   TraceSequenceLine
	int iTrace;						    //12	       2   TraceNumber
	int Offset;                                   //36           3   Offset £¨for each trace£©
	int sx;                                       //72           4   sx
	int gx;                                       //80           5   gx
	short Nt;                            //114          6   nt
	short Dt;                            //116          7   dt£¨us£©
	float Dx;                                     //188          8   dx  
	int Ntrace;				          	    //208          9   nx
}sgyheadkey;

//=============================================================//
//=======================swapbytes_function.cpp==================//
//=============================================================//
void swap_short_2(short *tni2);
void swap_u_short_2(unsigned short *tni2); 
void swap_int_4(int *tni4);
void swap_u_int_4(unsigned int *tni4); 
void swap_long_4(long *tni4); 
void swap_u_long_4(unsigned long *tni4);
void swap_float_4(float *tnf4);
//=============================================================//
//=======================sgyhead_function.cpp==================//
//=============================================================//
void init_sgyhead(sgyhead *header);
void swap_sgyhead(sgyhead *header);
void swap_sgyheadkey(sgyheadkey *headerkey);
void fread_sgyhead(sgyheadkey *headerkey, FILE *fp, int Endian);
void fwrite_sgyhead(sgyhead *header, sgyheadkey headerkey, FILE *fp, int Endian);
void float_to_float(float *trace, int nt, int endianflang);
void int_to_int(int *trace, int nt, int endianflag);
void short_to_short(short int *trace, int nt, int endianflag);
//=============================================================//
//=======================IO_function.cpp=======================//
//=============================================================//  
void Outseisfile2d_single(float **seiscal, sgyheadkey headerkey, int outflag, int endianflag, char buff[40]);
void Outseisfile2d(float **seiscal, sgyheadkey headerkey, int outflag, int sgyflag, int endianflag, char buff[40]);
void Inseisfile2d(float **seisobs, sgyheadkey *headerkey, int inflag, int sgyflag, char buff[40]);

void Outseisfile1d_single(float *seiscal, sgyheadkey headerkey, int outflag, int endianflag, char buff[40]);
void Outseisfile1d(float *seiscal, sgyheadkey headerkey, int outflag, int sgyflag, int endianflag, char buff[40]);
void Inseisfile1d(float *seisobs, sgyheadkey *headerkey, int inflag, int sgyflag, char buff[40]);

#endif
