/*
 *  PROGRAM NAME:  multifario
 *
 *  (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
 *  CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
 *
 *  Please refer to the LICENSE file in the top directory
 *
 *      author: Mike Henderson mhender@watson.ibm.com
 */

static char *id="@(#) $Id: Wunderlich.c,v 1.2 2007/06/08 21:01:20 mhender Exp $";

#include <MFAtlas.h>
#include <MFNRegion.h>
#include <MFNVector.h>
#include <math.h>
#include <MFFortran.h>
#include <MFMultifariosMethod.h>

int DrawWunderlich(MFNVector,double*,void*,MFErrorHandler);

static void F(int *n,double *u,int *m,double *f, void *data, MFErrorHandler e)
 {
  static double pi = 3.141592653589;
  static double l1 = 9;
  static double ct1 = 1;
  static double st1 = 0;
  static double l2 = 3;
  static double l3 = 3;
  static double l4 = 9;
  static double l5 = 6;
  static double l6 = 6;
  static double l7 = 9;
  static double l8 = 3;
  static double l9 = 3;
  static double l10 = 9;
  static double l11 = 6;
  static double l12 = 6;
  
  double st2;
  double ct2;
  double st3;
  double ct3;
  double st4;
  double ct4;
  double st5;
  double ct5;
  double st6;
  double ct6;
  double st7;
  double ct7;
  double st8;
  double ct8;
  double st9;
  double ct9;
  double st10;
  double ct10;
  double st11;
  double ct11;
  double st12;
  double ct12;
  
  st2=u[0];
  ct2=u[1];
  st3=u[2];
  ct3=u[3];
  st4=u[4];
  ct4=u[5];
  st5=u[6];
  ct5=u[7];
  st6=u[8];
  ct6=u[9];
  st7=u[10];
  ct7=u[11];
  st8=u[12];
  ct8=u[13];
  st9=u[14];
  ct9=u[15];
  st10=u[16];
  ct10=u[17];
  st11=u[18];
  ct11=u[19];
  st12=u[20];
  ct12=u[21];
  
  f[0]=l2*ct2 + l11*ct11 + l10*ct10 - l12*ct12 - l3*ct3 - l1*ct1;
  f[1]=l2*st2 + l11*st11 + l10*st10 - l12*st12 - l3*st3 - l1*st1;
  f[2]=l6*ct6 + l9*ct9 - l12*ct12 - l3*ct3;
  f[3]=l6*st6 + l9*st9 - l12*st12 - l3*st3;
  f[4]=l5*ct5 + l8*ct8 - l11*ct11 - l2*ct2;
  f[5]=l5*st5 + l8*st8 - l11*st11 - l2*st2;
  f[6]=l5*ct5 + l7*ct7 - l6*ct6 - l1*ct1;
  f[7]=l5*st5 + l7*st7 - l6*st6 - l1*st1;
  f[8]=l2*ct2 + l4*ct4 - l3*ct3 - l1*ct1;
  f[9]=l2*st2 + l4*st4 - l3*st3 - l1*st1;
  
  f[10]=st2*st2 + ct2*ct2 - 1;
  f[11]=st3*st3 + ct3*ct3 - 1;
  f[12]=st4*st4 + ct4*ct4 - 1;
  f[13]=st5*st5 + ct5*ct5 - 1;
  f[14]=st6*st6 + ct6*ct6 - 1;
  f[15]=st7*st7 + ct7*ct7 - 1;
  f[16]=st8*st8 + ct8*ct8 - 1;
  f[17]=st9*st9 + ct9*ct9 - 1;
  f[18]=st10*st10 + ct10*ct10 - 1;
  f[19]=st11*st11 + ct11*ct11 - 1;
  f[20]=st12*st12 + ct12*ct12 - 1;

  return;
 }

static void dF(int *n,double *u,int *m,double *df, void *data, MFErrorHandler e)
 {
  static double pi = 3.141592653589;
  static double l1 = 9;
  static double ct1 = 1;
  static double st1 = 0;
  static double l2 = 3;
  static double l3 = 3;
  static double l4 = 9;
  static double l5 = 6;
  static double l6 = 6;
  static double l7 = 9;
  static double l8 = 3;
  static double l9 = 3;
  static double l10 = 9;
  static double l11 = 6;
  static double l12 = 6;
  
  double st2;
  double ct2;
  double st3;
  double ct3;
  double st4;
  double ct4;
  double st5;
  double ct5;
  double st6;
  double ct6;
  double st7;
  double ct7;
  double st8;
  double ct8;
  double st9;
  double ct9;
  double st10;
  double ct10;
  double st11;
  double ct11;
  double st12;
  double ct12;

  int i;

  for(i=0;i<20*21;i++)df[i]=0.;
  
  st2=u[0];
  ct2=u[1];
  st3=u[2];
  ct3=u[3];
  st4=u[4];
  ct4=u[5];
  st5=u[6];
  ct5=u[7];
  st6=u[8];
  ct6=u[9];
  st7=u[10];
  ct7=u[11];
  st8=u[12];
  ct8=u[13];
  st9=u[14];
  ct9=u[15];
  st10=u[16];
  ct10=u[17];
  st11=u[18];
  ct11=u[19];
  st12=u[20];
  ct12=u[21];
  
  df[0+20* 1]=l2;
  df[0+20*19]=l11;
  df[0+20*17]=l10;
  df[0+20*21]=-l12;
  df[0+20* 3]=-l3;

  df[1+20* 0]=l2;
  df[1+20*18]=l11;
  df[1+20*16]=l10;
  df[1+20*20]=-l12;
  df[1+20* 2]=-l3;

  df[2+20* 9]=l6;
  df[2+20*15]=l9;
  df[2+20*21]=-l12;
  df[2+20* 3]=l3;

  df[3+20* 8]=l6;
  df[3+20*14]=l9;
  df[3+20*20]=-l12;
  df[3+20* 2]=-l3;

  df[4+20* 7]=l5;
  df[4+20*13]=l8;
  df[4+20*19]=-l11;
  df[4+20* 1]=-l2;

  df[5+20* 6]=l5;
  df[5+20*12]=l8;
  df[5+20*18]=-l11;
  df[5+20* 0]=-l2;

  df[6+20* 7]=l5;
  df[6+20*11]=l7;
  df[6+20* 9]=-l6;

  df[7+20* 6]=l5;
  df[7+20*10]=l7;
  df[7+20* 8]=-l6;

  df[8+20* 1]=l2;
  df[8+20* 5]=l4;
  df[8+20* 3]=-l3;

  df[9+20* 0]=l2;
  df[9+20* 4]=l4;
  df[9+20* 2]=-l3;

  df[10+20* 0]=2*st2;
  df[10+20* 1]=2*ct2;

  df[11+20* 2]=2*st3;
  df[11+20* 3]=2*ct3;

  df[12+20* 4]=2*st4;
  df[12+20* 5]=2*ct4;

  df[13+20* 6]=2*st5;
  df[13+20* 7]=2*ct5;

  df[14+20* 8]=2*st6;
  df[14+20* 9]=2*ct6;

  df[15+20*10]=2*st7;
  df[15+20*11]=2*ct7;

  df[16+20*12]=2*st8;
  df[16+20*13]=2*ct8;

  df[17+20*14]=2*st9;
  df[17+20*15]=2*ct9;

  df[18+20*16]=2*st10;
  df[18+20*17]=2*ct10;

  df[19+20*18]=2*st11;
  df[19+20*19]=2*ct11;

  df[20+20*20]=2*st12;
  df[20+20*21]=2*ct12;

  return;
 }

static void ddF(int *n,double *u,int *m, double *ddf, void *data, MFErrorHandler e)
 {
  static double pi = 3.141592653589;
  static double l1 = 9;
  static double ct1 = 1;
  static double st1 = 0;
  static double l2 = 3;
  static double l3 = 3;
  static double l4 = 9;
  static double l5 = 6;
  static double l6 = 6;
  static double l7 = 9;
  static double l8 = 3;
  static double l9 = 3;
  static double l10 = 9;
  static double l11 = 6;
  static double l12 = 6;
  
  double st2;
  double ct2;
  double st3;
  double ct3;
  double st4;
  double ct4;
  double st5;
  double ct5;
  double st6;
  double ct6;
  double st7;
  double ct7;
  double st8;
  double ct8;
  double st9;
  double ct9;
  double st10;
  double ct10;
  double st11;
  double ct11;
  double st12;
  double ct12;

  int i;

  for(i=0;i<20*21*21;i++)ddf[i]=0.;
  
  st2=u[0];
  ct2=u[1];
  st3=u[2];
  ct3=u[3];
  st4=u[4];
  ct4=u[5];
  st5=u[6];
  ct5=u[7];
  st6=u[8];
  ct6=u[9];
  st7=u[10];
  ct7=u[11];
  st8=u[12];
  ct8=u[13];
  st9=u[14];
  ct9=u[15];
  st10=u[16];
  ct10=u[17];
  st11=u[18];
  ct11=u[19];
  st12=u[20];
  ct12=u[21];
  
  ddf[10+20*( 0+20* 0)]=2;
  ddf[10+20*( 1+20* 1)]=2;

  ddf[11+20*( 2+20* 2)]=2;
  ddf[11+20*( 3+20* 3)]=2;

  ddf[12+20*( 4+20* 4)]=2;
  ddf[12+20*( 5+20* 5)]=2;

  ddf[13+20*( 6+20* 6)]=2;
  ddf[13+20*( 7+20* 7)]=2;

  ddf[14+20*( 8+20* 8)]=2;
  ddf[14+20*( 9+20* 9)]=2;

  ddf[15+20*(10+20*10)]=2;
  ddf[15+20*(11+20*11)]=2;

  ddf[16+20*(12+20*12)]=2;
  ddf[16+20*(13+20*13)]=2;

  ddf[17+20*(14+20*14)]=2;
  ddf[17+20*(15+20*15)]=2;

  ddf[18+20*(16+20*16)]=2;
  ddf[18+20*(17+20*17)]=2;

  ddf[19+20*(18+20*18)]=2;
  ddf[19+20*(19+20*19)]=2;

  ddf[20+20*(20+20*20)]=2;
  ddf[20+20*(21+20*21)]=2;

  return;
 }

int main(int argc, char *argv[])
 {
  MFImplicitMF M;
  int i,j;
  int n,k;
  MFNRegion Omega;
  MFAtlas S;
  MFNVector u0;
  MFContinuationMethod H;
  MFErrorHandler e;

  e=MFCreateErrorHandler();

  n=21;
  k=1;

  M=MFIMFCreateAlgebraicSubroutine(n,k,F,dF,ddF,NULL,e);
  MFIMFSetProjectForDraw(M,DrawWunderlich,e);

  Omega=MFNRegionCreateHyperCube(n,1.1,e);

  u0=MFIMFVectorFactory(M,e);
  MFNVSetC(u0, 0,0.656064,e);
  MFNVSetC(u0, 1,0.726026,e);
  MFNVSetC(u0, 2,0.687668,e);
  MFNVSetC(u0, 3,0.754705,e);
  MFNVSetC(u0, 4,0.656066,e);
  MFNVSetC(u0, 5,0.726027,e);
  MFNVSetC(u0, 6,0.687666,e);
  MFNVSetC(u0, 7,0.754703,e);
  MFNVSetC(u0, 8,-0.000568,e);
  MFNVSetC(u0, 9,0.000568,e);
  MFNVSetC(u0,10,1.000000,e);
  MFNVSetC(u0,11,1.000000,e);
  MFNVSetC(u0,12,0.660063,e);
  MFNVSetC(u0,13,0.723646,e);
  MFNVSetC(u0,14,0.690172,e);
  MFNVSetC(u0,15,0.751211,e);
  MFNVSetC(u0,16,0.660105,e);
  MFNVSetC(u0,17,0.723645,e);
  MFNVSetC(u0,18,0.690173,e);
  MFNVSetC(u0,19,0.751173,e);
  MFNVSetC(u0,20,-0.000941,e);
  MFNVSetC(u0,21,0.000943,e);

  H=MFCreateMultifariosMethod(e);
  MFMultifarioSetRealParameter(H,"epsilon",.01,e);
  MFMultifarioSetIntegerParameter(H,"maxCharts",100,e);
  MFMultifarioSetIntegerParameter(H,"verbose",1,e);
  MFMultifarioSetIntegerParameter(H,"page",1,e);
  MFMultifarioSetIntegerParameter(H,"dumpToPlotFile",1,e);
  MFMultifarioSetIntegerParameter(H,"dumpToCenterFile",0,e);
  MFMultifarioSetFilename(H,"SphereSub",e);

  S=MFComputeAtlas(H,M,Omega,u0,e);

  MFCloseAtlas(H,S,e);
  printf("Done computating Atlas\n");fflush(stdout);

  MFFreeAtlas(S,e);
  MFFreeImplicitMF(M,e);
  MFFreeNRegion(Omega,e);
  MFFreeNVector(u0,e);
  MFFreeContinuationMethod(H,e);

  MFFreeErrorHandler(e);

  return 0;
 }

int DrawWunderlich(MFNVector v,double *u,void *d, MFErrorHandler e)
 {
  if(u==(double*)NULL||v==(MFNVector)NULL) return 3;

  u[0]=MFNV_C(v,0,e);
  u[1]=MFNV_C(v,1,e);
  u[2]=MFNV_C(v,2,e);

  return 0;
 }
