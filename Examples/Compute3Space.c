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

static char *id="@(#) $Id: Compute3Space.c,v 1.3 2011/07/21 17:43:45 mhender Exp $";

#include <MFAtlas.h>
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

  M=MFIMFCreateNSpaceWithRadius(3,.1,e);
  n=MFIMF_N(M,e);
  Omega=MFNRegionCreateCube(-1.,-1.,-1.,1.,1.,1.,e);

  u0=MFIMFVectorFactory(M,e);
  MFNVSetC(u0,0, 0.,e);
  MFNVSetC(u0,1, 0.,e);
  MFNVSetC(u0,2, 0.,e);

  H=MFCreateMultifariosMethod(e);
  MFMultifarioSetRealParameter(H,"epsilon",.05,e);
  MFMultifarioSetIntegerParameter(H,"maxCharts",-1,e);
  MFMultifarioSetIntegerParameter(H,"verbose",1,e);
  MFMultifarioSetIntegerParameter(H,"page",1,e);
  MFMultifarioSetIntegerParameter(H,"dumpToPlotFile",1,e);
  MFMultifarioSetIntegerParameter(H,"dumpToCenterFile",0,e);
  MFMultifarioSetFilename(H,"3Space",e);

  S=MFComputeAtlas(H,M,Omega,u0,e);

  MFCloseAtlas(H,S,e);
  printf("Done computing Atlas\n");fflush(stdout);

  MFFreeAtlas(S,e);
  MFFreeImplicitMF(M,e);
  MFFreeContinuationMethod(H,e);
  MFFreeNRegion(Omega,e);
  MFFreeNVector(u0,e);

  MFFreeErrorHandler(e);

  return 0;
 }
