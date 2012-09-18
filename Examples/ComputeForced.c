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

static char *id="@(#) $Id: ComputeForced.c,v 1.4 2011/07/21 17:43:45 mhender Exp $";

#include <MFAtlas.h>
#include <MFNRegion.h>
#include <MFNVector.h>
#include <MFForced.h>
#include <math.h>
#include <MFMultifariosMethod.h>

#define MFPI 3.14159265358979323846264338327950288

int main(int argc, char *argv[])
 {
  MFImplicitMF M;
  int i;
  int n;
  MFNRegion Omega;
  MFAtlas S;
  MFNVector u0;
  FILE *fid;
  MFNVector ug;
  MFNVector t[2];
  MFNKMatrix Tan;
  double *tanData;
  double tau;
  MFErrorHandler eh;

  int nt=600;

  double p=1.;
  int q=4;
  double e,B;
  MFContinuationMethod H;

  eh=MFCreateErrorHandler();

  M=MFIMFCreateForcedOscillator(nt,q,eh);
  n=nt+2;

  H=MFCreateMultifariosMethod(eh);
  MFMultifarioSetRealParameter(H,"epsilon",1.,eh);            /* Max distance from TS to M */
  MFMultifarioSetIntegerParameter(H,"maxCharts",-1,eh);       /* -1 means as many as needed */
  MFMultifarioSetIntegerParameter(H,"verbose",1,eh);          /* Write info to stdout */
  MFMultifarioSetIntegerParameter(H,"page",1,eh);             /* Page out non-interior polyhedra */
  MFMultifarioSetIntegerParameter(H,"dumpToPlotFile",1,eh);   /* Write polyhedra to a plotfile */
  MFMultifarioSetIntegerParameter(H,"dumpToCenterFile",0,eh); /* Write points to a file */

#define FOURONE
#ifdef FOURONE
  MFMultifarioSetFilename(H,"ForcedOscillatorFourOneA",eh);
  p=1.;
  q=4;
  e=11.4;
  B=.05*20;
#endif
#ifdef THREEEONE
  MFMultifarioSetFilename(H,"ForcedOscillatorThreeOne",eh);
  p=1.;
  q=3;
  e=7.5;
  B=.04*20;
#endif

  ug=MFIMFVectorFactory(M,eh);

  for(i=0;i<nt;i++)
   {
    tau=q*2*MFPI*i/nt;
    MFNVSetC(ug,i,2*sin(p*tau/q),eh);
   }
  MFNVSetC(ug,nt+0,e,eh);
  MFNVSetC(ug,nt+1,B,eh);
  ForcedWrite(ug,eh);

  tanData=(double*)malloc(n*2*sizeof(double));
  if(tanData==NULL)
   {
    printf("Couldn't allocate the tangent\n");return 8;
   }
  for(i=0;i<nt;i++)
    tanData[i+n*0]=0.;
  tanData[nt+0+n*0]=1.;
  tanData[nt+1+n*0]=0.;

  for(i=0;i<nt;i++)
    tanData[i+n*1]=0.;
  tanData[nt+0+n*1]=0.;
  tanData[nt+1+n*1]=1.;

  Tan=MFCreateNKMatrixWithData(n,2,tanData,eh);

  u0=MFIMFVectorFactory(M,eh);
  if(!MFIMFProject(M,ug,Tan,u0,eh))
   {
    printf("Initial project failed\n");return 8;
   }
  MFFreeNVector(ug,eh);
  MFFreeNKMatrix(Tan,eh);

/*Omega=MFNRegionCreateForcedOscillator(n,q,0.,20.,0.,.5*20,eh);*/
  Omega=MFNRegionCreateForcedOscillator(n,q,0.,30.,0.,.8*20,eh);

  S=MFComputeAtlas(H,M,Omega,u0,eh);
  MFCloseAtlas(H,S,eh);
  printf("Done writing Atlas\n");fflush(stdout);

  MFFreeAtlas(S,eh);
  MFFreeImplicitMF(M,eh);
  MFFreeNRegion(Omega,eh);
  MFFreeNVector(u0,eh);
  ForcedCloser(eh);
  MFFreeContinuationMethod(H,eh);
  MFFreeErrorHandler(eh);

  return 0;
 }
