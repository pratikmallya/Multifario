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

static char *id="@(#) $Id: ComputePlaneClip.c,v 1.3 2011/07/21 17:43:45 mhender Exp $";

#include <MFAtlas.h>
#include <MFMultifariosMethod.h>

double XGreaterThanMinusOne(MFNVector u, MFErrorHandler e){return -MFNV_C(u,0,e)-1.;}
double XLessThanOne(MFNVector u, MFErrorHandler e){return MFNV_C(u,0,e)-1.;}
double YGreaterThanMinusOne(MFNVector u, MFErrorHandler e){return -MFNV_C(u,1,e)-1.;}
double YLessThanOne(MFNVector u, MFErrorHandler e){return MFNV_C(u,1,e)-1.;}

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

  M=MFIMFCreatePlane(e);
  n=MFIMF_N(M,e);
  Omega=MFNRegionCreateRectangle(-1.,-1.,1.,1.,e);

  u0=MFIMFVectorFactory(M,e);
  MFNVSetC(u0,0, 0.,e);
  MFNVSetC(u0,1, 0.,e);

  H=MFCreateMultifariosMethod(e);
  MFMultifarioSetRealParameter(H,"epsilon",.1,e);
  MFMultifarioSetIntegerParameter(H,"maxCharts",-1,e);
  MFMultifarioSetIntegerParameter(H,"verbose",1,e);
  MFMultifarioSetIntegerParameter(H,"page",1,e);
  MFMultifarioSetIntegerParameter(H,"dumpToPlotFile",1,e);
  MFMultifarioSetIntegerParameter(H,"dumpToCenterFile",0,e);
  MFMultifarioSetFilename(H,"PlaneClip",e);

  MFMultifarioAddClipF(H,XGreaterThanMinusOne,e);
  MFMultifarioAddClipF(H,XLessThanOne,e);
  MFMultifarioAddClipF(H,YGreaterThanMinusOne,e);
  MFMultifarioAddClipF(H,YLessThanOne,e);

  S=MFComputeAtlas(H,M,Omega,u0,e);

  MFCloseAtlas(H,S,e);
  printf("Done computing Atlas\n");fflush(stdout);

  MFFreeAtlas(S,e);
  MFFreeImplicitMF(M,e);
  MFFreeNRegion(Omega,e);
  MFFreeNVector(u0,e);
  MFFreeContinuationMethod(H,e);

  MFFreeErrorHandler(e);

  return 0;
 }
