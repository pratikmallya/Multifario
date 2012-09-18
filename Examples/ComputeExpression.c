/* 
    @(#)ComputeExpression.c	1.7
    03/02/25 10:33:35
   
    PROGRAM NAME:  multifario

    (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
    CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.

    Please refer to the LICENSE file in the top directory
*/
/*      author: Mike Henderson mhender@watson.ibm.com */

#include <MFAtlas.h>
#include <MFNRegion.h>
#include <MFNVector.h>
#include <math.h>
#include <MFMultifariosMethod.h>

int main(int argc, char *argv[])
 {
  MFImplicitMF M;
  int n;
  MFNRegion Omega;
  MFAtlas S;
  MFNVector u0;
  FILE *fid;
  MFContinuationMethod H;
  MFErrorHandler e;

  e=MFCreateErrorHandler();

  M=MFIMFCreateAlgebraicExpression("[x,y,z]","[x**2+y**2-1+.1*sin(6*3.1415926*x),z-x*y]",e);
  n=MFIMF_N(M,e);
  Omega=MFNRegionCreateHyperCube(n,1.1,e);

  u0=MFIMFVectorFactory(M,e);
  MFNVSetC(u0,0, 0.,e);
  MFNVSetC(u0,1, 1.,e);
  MFNVSetC(u0,2, 0.,e);

  H=MFCreateMultifariosMethod(e);
  MFMultifarioSetRealParameter(H,"epsilon",.001,e);
  MFMultifarioSetIntegerParameter(H,"maxCharts",3000,e);
  MFMultifarioSetIntegerParameter(H,"verbose",1,e);
  MFMultifarioSetIntegerParameter(H,"page",0,e);
  MFMultifarioSetIntegerParameter(H,"dumpToPlotFile",1,e);
  MFMultifarioSetIntegerParameter(H,"dumpToCenterFile",0,e);
  MFMultifarioSetFilename(H,"Expression",e);

  S=MFComputeAtlas(H,M,Omega,u0,e);

  MFCloseAtlas(H,S,e);
  printf("Done computing Atlas\n");fflush(stdout);

  MFFreeAtlas(S,e);
  MFFreeContinuationMethod(H,e);
  MFFreeImplicitMF(M,e);
  MFFreeNRegion(Omega,e);
  MFFreeNVector(u0,e);
  MFFreeErrorHandler(e);

  return(0);
 }
