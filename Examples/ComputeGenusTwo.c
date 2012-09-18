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

static char *id="@(#) $Id: ComputeGenusTwo.c,v 1.3 2011/07/21 17:43:45 mhender Exp $";

#include <MFAtlas.h>
#include <MFErrorHandler.h>
#include <MFMultifariosMethod.h>

int main(int argc, char *argv[])
 {
  MFImplicitMF M;
  int n;
  MFNRegion Omega;
  MFAtlas S;
  MFNVector u0,u;
  MFNKMatrix Phi;
  FILE *fid;
  MFContinuationMethod H;
  MFErrorHandler e;

  e=MFCreateErrorHandler();

  M=MFIMFCreateAlgebraicExpression("[x,y,z]",
  "[((sqrt((x-.8)**2+y**2)-.8)**2+z**2 -.5**2)*((sqrt((x+.8)**2+y**2)-.8)**2+z**2 -.5**2)-.05]",e);
  n=MFIMF_N(M,e);
  Omega=MFNRegionCreateHyperCube(n,2.4,e);

/* This bit starts with an approximate point on M, and calls project to get a point on M */

  u=MFIMFVectorFactory(M,e);
  MFNVSetC(u,0, 1.2,e);
  MFNVSetC(u,1, 0.,e);
  MFNVSetC(u,2, 0.,e);
  u0=MFIMFVectorFactory(M,e);
  Phi=MFIMFMatrixFactory(M,e);

  MFNKMSetC(Phi,0,0,0.,e);
  MFNKMSetC(Phi,1,0,1.,e);
  MFNKMSetC(Phi,2,0,0.,e);
  MFNKMSetC(Phi,0,1,0.,e);
  MFNKMSetC(Phi,1,1,0.,e);
  MFNKMSetC(Phi,2,1,1.,e);

  if(!MFIMFProject(M,u,Phi,u0,e))return 8;

  H=MFCreateMultifariosMethod(e);
  MFMultifarioSetRealParameter(H,"epsilon",.01,e);
  MFMultifarioSetIntegerParameter(H,"maxCharts",-1,e);
  MFMultifarioSetIntegerParameter(H,"verbose",1,e);
  MFMultifarioSetIntegerParameter(H,"page",1,e);
  MFMultifarioSetIntegerParameter(H,"dumpToPlotFile",1,e);
  MFMultifarioSetIntegerParameter(H,"dumpToCenterFile",0,e);
  MFMultifarioSetFilename(H,"GenusTwo",e);

  S=MFComputeAtlas(H,M,Omega,u0,e);

  MFCloseAtlas(H,S,e);
  printf("Done computing Atlas\n");fflush(stdout);

  MFFreeAtlas(S,e);
  MFFreeImplicitMF(M,e);
  MFFreeContinuationMethod(H,e);
  MFFreeNRegion(Omega,e);
  MFFreeNVector(u0,e);

  MFFreeErrorHandler(e);

  return(0);
 }
