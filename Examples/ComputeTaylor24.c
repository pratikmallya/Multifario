/* 
    %W%
    %D% %T%
   
    PROGRAM NAME:  multifario

    (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
    CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.

    Please refer to the LICENSE file in the top directory

*/
/*      author: Mike Henderson mhender@watson.ibm.com */

#include <MFAtlas.h>
#include <MFNRegion.h>
#include <MFNVector.h>
#include <MFPrint.h>
#include <math.h>
#include <MFMultifariosMethod.h>

int MFTaylor24ProjectToDraw(MFNVector,double*,void*,MFErrorHandler);

void findSolution(double,double,double,double,double,double,double*,double*);
int main(int argc, char *argv[])
 {
  MFImplicitMF M;
  MFNRegion Omega;
  int n;
  MFAtlas S;
  MFNVector u0;
  MFNVector ll,ur;
  MFContinuationMethod H;
  MFErrorHandler e;

  e=MFCreateErrorHandler();

/* Model for Taylor (2,4) mode interaction */

  printf("Expression is %s\n","[x*(2*x*x+y*y-r+s)+x*y,y*(x*x+2*y*y-r-s)-x*x]");fflush(stdout);
  M=MFIMFCreateAlgebraicExpressionWithRadius("[x,y,r,s]","[x*(2*x*x+y*y-r+s)+x*y,y*(x*x+2*y*y-r-s)-x*x]",.05,e);
  MFIMFSetProjectForDraw(M,MFTaylor24ProjectToDraw,e);

  n=MFIMF_N(M,e);
  ll=MFIMFVectorFactory(M,e);
  ur=MFIMFVectorFactory(M,e);
  MFNVSetC(ll,0,-5.,e);MFNVSetC(ur,0,5.,e);
  MFNVSetC(ll,1,-5.,e);MFNVSetC(ur,1,5.,e);
  MFNVSetC(ll,2,-1.,e);MFNVSetC(ur,2,1.,e);
  MFNVSetC(ll,3,-1.,e);MFNVSetC(ur,3,1.,e); 
  Omega=MFNRegionCreateHyperCubeByCorners(n,ll,ur,e);
  MFFreeNVector(ll,e);
  MFFreeNVector(ur,e);
  
  H=MFCreateMultifariosMethod(e);
  MFMultifarioSetRealParameter(H,"epsilon",.01,e);
  MFMultifarioSetIntegerParameter(H,"maxCharts",-1,e);
  MFMultifarioSetIntegerParameter(H,"verbose",1,e);
  MFMultifarioSetIntegerParameter(H,"page",1,e);
  MFMultifarioSetIntegerParameter(H,"dumpToPlotFile",1,e);
  MFMultifarioSetIntegerParameter(H,"dumpToCenterFile",0,e);
  MFMultifarioSetFilename(H,"Circle",e);

  u0=MFIMFVectorFactory(M,e);

/* Trivial */

  MFNVSetC(u0,0,0.,e);
  MFNVSetC(u0,1,0.,e);
  MFNVSetC(u0,2,-.4,e);
  MFNVSetC(u0,3,0.,e);
  MFMultifarioSetFilename(H,"Taylor24",e);

  S=MFComputeAtlas(H,M,Omega,u0,e);

  MFCloseAtlas(H,S,e);
  MFFreeAtlas(S,e);
  MFFreeImplicitMF(M,e);
  MFFreeNRegion(Omega,e);
  MFFreeNVector(u0,e);
  MFFreeContinuationMethod(H,e);
  MFFreeErrorHandler(e);

  return(0);
 }

int MFTaylor24ProjectToDraw(MFNVector u, double *x, void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFTaylor24ProjectToDraw"};

  if(x==(double*)NULL)return 3;

  x[0]=MFNV_C(u,2,e);
  x[1]=MFNV_C(u,3,e);
  x[2]=.5*(MFNV_C(u,0,e)+MFNV_C(u,1,e));

  return 0;
 }
