
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

static char *id="@(#) $Id: ComputeSphere.c,v 1.3 2011/07/21 17:43:45 mhender Exp $";
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
#include <MFMultifariosMethod.h>

int main(int argc, char *argv[])
 {
  MFImplicitMF M;
  MFNRegion Omega;
  MFAtlas A;
  MFNVector u0;
  MFContinuationMethod H;
  MFErrorHandler e;

  e=MFCreateErrorHandler();

  M=MFIMFCreateAlgebraicExpression("[x,y,z]",
                                   "[x*x+y*y+z*z-1.]",e);
  MFIMFSetR(M,.1,e);
  Omega=MFNRegionCreateHyperCube(3,1.1,e);

  u0=MFIMFVectorFactory(M,e);
  MFNVSetC(u0,0, 0.,e);
  MFNVSetC(u0,1, 1.,e);
  MFNVSetC(u0,2, 0.,e);


  H=MFCreateMultifariosMethod(e);
  MFMultifarioSetRealParameter(H,"epsilon",.01,e);            /* Max distance from TS to M */
  MFMultifarioSetIntegerParameter(H,"maxCharts",-1,e);       /* -1 means as many as needed */
  MFMultifarioSetIntegerParameter(H,"verbose",1,e);          /* Write info to stdout */
  MFMultifarioSetIntegerParameter(H,"page",1,e);             /* Page out non-interior polyhedra */
  MFMultifarioSetIntegerParameter(H,"dumpToPlotFile",1,e);   /* Write polyhedra to a plotfile */
  MFMultifarioSetIntegerParameter(H,"dumpToCenterFile",0,e); /* Write points to a file */
  MFMultifarioSetFilename(H,"Sphere",e);

  A=MFComputeAtlas(H,M,Omega,u0,e);

  MFCloseAtlas(H,A,e);
  printf("Done computing Atlas\n");fflush(stdout);

  MFFreeAtlas(A,e);
  MFFreeImplicitMF(M,e);
  MFFreeNRegion(Omega,e);
  MFFreeNVector(u0,e);
  MFFreeContinuationMethod(H,e);

  MFFreeErrorHandler(e);

  return 0;
 }
