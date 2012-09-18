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

void findSolution(double,double,double,double,double,double,double*,double*);
int main(int argc, char *argv[])
 {
  MFImplicitMF M;
  int i,j,n;
  int n0;
  MFNRegion Omega;
  MFAtlas S;
  MFNVector u0;
  MFNVector ll,ur;
  MFContinuationMethod H;
  MFErrorHandler e;

  e=MFCreateErrorHandler();

/* Complex Cusp Catastrophe */

/*   u*(u*u-r)-s = 0  */
/*   (x+i*y)*((x+i*y)*(x+i*y)-r)-s = 0  */
/*   (x+i*y)*((x*x-y*y-r)+i*2*x*y)-s = 0  */
/*   x*((x*x-y*y-r)+i*2*x*y)+i*y*((x*x-y*y-r)+i*2*x*y)-s = 0  */

/*   x*(x*x-y*y-r)-2*x*y*y - s = 0 */
/*   2*x*x*y+y*(x*x-y*y-r) = 0  */

  M=MFIMFCreateAlgebraicExpressionWithRadius("[x,y,r,s]","[x*(x*x-y*y-r)-2*x*y*y-s,2*x*x*y+y*(x*x-2*y*y-r)]",.05,e);

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
  MFMultifarioSetRealParameter(H,"Epsilon",.01,e);
  MFMultifarioSetRealParameter(H,"DotMin",.9,e);
  MFMultifarioSetIntegerParameter(H,"maxCharts",100000,e);
  MFMultifarioSetIntegerParameter(H,"verbose",1,e);
  MFMultifarioSetIntegerParameter(H,"page",0,e);
  MFMultifarioSetIntegerParameter(H,"branchSwitch",200,e);
  MFMultifarioSetIntegerParameter(H,"dumpToPlotFile",1,e);   /* Write polyhedra to a plotfile */
  MFMultifarioSetIntegerParameter(H,"dumpToCenterFile",0,e); /* Write points to a file */
  MFMultifarioSetIntegerParameter(H,"dumpToRestartFile",0,e);

  u0=MFIMFVectorFactory(M,e);

/* Trivial */

  MFNVSetC(u0,0,-.4,e);
  MFNVSetC(u0,1,0.,e);
  MFNVSetC(u0,2,0.,e);
  MFNVSetC(u0,3,-.4*.4*.4,e);
  MFMultifarioSetFilename(H,"Cusp",e);

  S=MFComputeAtlas(H,M,Omega,u0,e);

  MFCloseAtlas(H,S,e);
  MFFreeAtlas(S,e);
  MFFreeContinuationMethod(H,e);
  MFFreeImplicitMF(M,e);
  MFFreeNRegion(Omega,e);
  MFFreeNVector(u0,e);
  MFFreeErrorHandler(e);

  return 0;
 }
