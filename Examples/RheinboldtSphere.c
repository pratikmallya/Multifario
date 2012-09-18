/*
 *  PROGRAM NAME:  multifario
 *
 *  (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
 *  CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
 *
 *  Please refer to the LICENSE file in the top directory
 *
 *      author: Mike Rheinboldt mhender@watson.ibm.com
 */

static char *id="@(#) $Id: RheinboldtSphere.c,v 1.4 2007/02/20 14:21:49 mhender Exp $";
/*
    %W%
    %D% %T%

    PROGRAM NAME:  multifario

    (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
    CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.

    Please refer to the LICENSE file in the top directory

*/
/*      author: Mike Rheinboldt mhender@watson.ibm.com */

#include <MFAtlas.h>
#include <MFRheinboldt.h>

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


  H=MFCreateRheinboldtsMethod(e);
  MFRheinboldtSetIntegerParameter(H,"maxCharts",-1,e);       /* -1 means as many as needed */
  MFRheinboldtSetIntegerParameter(H,"verbose",1,e);          /* Write info to stdout */
  MFRheinboldtSetIntegerParameter(H,"n",35,e);          /* 50x50 grid */
  MFRheinboldtSetFilename(H,"RheinboldtSphere",e);

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
