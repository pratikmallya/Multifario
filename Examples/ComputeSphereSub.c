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

static char *id="@(#) $Id: ComputeSphereSub.c,v 1.4 2011/07/21 17:43:45 mhender Exp $";

#include <MFAtlas.h>
#include <MFNRegion.h>
#include <MFNVector.h>
#include <math.h>
#include <MFFortran.h>
#include <MFMultifariosMethod.h>

int DrawSphere(MFNVector,double*,void*,MFErrorHandler);

static void F(int *n,double *z,int *m,double *f, void *data, MFErrorHandler e)
 {
  f[0]= z[0]*z[0]+z[1]*z[1]+z[2]*z[2]-1.;

  return;
 }

static void dF(int *n,double *z,int *m,double *df, void *data, MFErrorHandler e)
 {
  df[0+1* 0]= 2*z[0];
  df[0+1* 1]= 2*z[1];
  df[0+1* 2]= 2*z[2];

  return;
 }

static void ddF(int *n,double *z,int *m, double *ddf, void *data, MFErrorHandler e)
 {
  ddf[0+1*(0+3*0)]= 2;
  ddf[0+1*(1+3*0)]= 0;
  ddf[0+1*(2+3*0)]= 0;

  ddf[0+1*(0+3*1)]= 0;
  ddf[0+1*(1+3*1)]= 2;
  ddf[0+1*(2+3*1)]= 0;

  ddf[0+1*(0+3*2)]= 0;
  ddf[0+1*(1+3*2)]= 0;
  ddf[0+1*(2+3*2)]= 2;

  return;
 }

int main(int argc, char *argv[])
 {
  MFImplicitMF M;
  int i,j;
  int n,k;
  MFNRegion Omega;
  MFAtlas S;
  MFNVector u0[2];
  MFContinuationMethod H;
  MFErrorHandler e;

  e=MFCreateErrorHandler();

  n=3;
  k=2;

  M=MFIMFCreateAlgebraicSubroutine(n,k,F,dF,ddF,NULL,e);
  MFIMFSetProjectForDraw(M,DrawSphere,e);

  Omega=MFNRegionCreateHyperCube(n,1.1,e);

  u0[0]=MFIMFVectorFactory(M,e);
  for(j=0;j<n;j++)MFNVSetC(u0[0],j,0.,e);
  MFNVSetC(u0[0],1,1.,e);

  u0[1]=MFIMFVectorFactory(M,e);
  for(j=0;j<n;j++)MFNVSetC(u0[1],j,0.,e);
  MFNVSetC(u0[1],2,1.,e);

  H=MFCreateMultifariosMethod(e);
  MFMultifarioSetRealParameter(H,"epsilon",.01,e);
  MFMultifarioSetIntegerParameter(H,"maxCharts",100,e);
  MFMultifarioSetIntegerParameter(H,"verbose",1,e);
  MFMultifarioSetIntegerParameter(H,"page",1,e);
  MFMultifarioSetIntegerParameter(H,"dumpToPlotFile",1,e);
  MFMultifarioSetIntegerParameter(H,"dumpToCenterFile",0,e);
  MFMultifarioSetFilename(H,"SphereSub",e);

  S=MFComputeAtlasMultiple(H,M,Omega,2,u0,e);

  MFCloseAtlas(H,S,e);
  printf("Done computating Atlas\n");fflush(stdout);

  MFFreeAtlas(S,e);
  MFFreeImplicitMF(M,e);
  MFFreeNRegion(Omega,e);
  MFFreeNVector(u0[0],e);
  MFFreeNVector(u0[1],e);
  MFFreeContinuationMethod(H,e);

  MFFreeErrorHandler(e);

  return 0;
 }

int DrawSphere(MFNVector v,double *u,void *d, MFErrorHandler e)
 {
  if(u==(double*)NULL||v==(MFNVector)NULL) return 3;

  u[0]=MFNV_C(v,0,e);
  u[1]=MFNV_C(v,1,e);
  u[2]=MFNV_C(v,2,e);

  return 0;
 }
