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

static char *id="@(#) $Id: PeriodFourResonance.c,v 1.2 2011/07/21 17:43:45 mhender Exp $";
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
  MFNVector u;
  MFNVector u0;
  MFNKMatrix Phi;
  MFNVector du0;
  MFNVector du1;
  MFContinuationMethod H;
  MFErrorHandler e;
  double x0,y0;
  double x1,y1;
  double a,b;
  double ll[10]={-4.,-12.,-1000.,-1000.,-1000.,-1000.,-1000.,-1000.,-1000.,-1000.};
  double ur[10]={ 4., 8., 1000., 1000., 1000., 1000., 1000., 1000., 1000., 1000.};
  MFNVector LL,UR;
  int i;

  e=MFCreateErrorHandler();

  M=MFIMFCreateAlgebraicExpression("[a,b,x0,y0,x1,y1,x2,y2,x3,y3]", "[x1-a*x0-y0,y1-b-x0*x0, x2-a*x1-y1,y2-b-x1*x1, x3-a*x2-y2,y3-b-x2*x2, x0-a*x3-y3,y0-b-x3*x3]",e);
  MFIMFSetR(M,.3,e);


  UR=MFCreateNVectorWithData(10,ur,e);
  LL=MFCreateNVectorWithData(10,ll,e);

  Omega=MFNRegionCreateHyperCubeByCorners(10,LL,UR,e);

  u =MFIMFVectorFactory(M,e);
  u0=MFIMFVectorFactory(M,e);

  a=0.001193;
  b=-.766831;
  MFNVSetC(u,0,a,e);
  MFNVSetC(u,1,b,e);

  x0=-.499487;
  y0=-.366997;
  MFNVSetC(u,2,x0,e);
  MFNVSetC(u,3,y0,e);
  printf(" (x0,y0)=(%f,%f)\n",x0,y0);fflush(stdout);

  x1=a*x0+y0;
  y1=b+x0*x0;
  x0=x1;
  y0=y1;
  MFNVSetC(u,4,x0,e);
  MFNVSetC(u,5,y0,e);
  printf(" (x1,y1)=(%f,%f)\n",x0,y0);fflush(stdout);

  x1=a*x0+y0;
  y1=b+x0*x0;
  x0=x1;
  y0=y1;
  MFNVSetC(u,6,x0,e);
  MFNVSetC(u,7,y0,e);
  printf(" (x2,y2)=(%f,%f)\n",x0,y0);fflush(stdout);

  x1=a*x0+y0;
  y1=b+x0*x0;
  x0=x1;
  y0=y1;
  MFNVSetC(u,8,x0,e);
  MFNVSetC(u,9,y0,e);
  printf(" (x3,y3)=(%f,%f)\n",x0,y0);fflush(stdout);

  Phi=MFIMFMatrixFactory(M,e);
  du0=MFMColumn(Phi,0,e);
  du1=MFMColumn(Phi,1,e);
  for(i=0;i<10;i++){MFNVSetC(du0,i,0.,e);MFNVSetC(du1,i,0.,e);}
  MFNVSetC(du0,0,1.,e);
  MFNVSetC(du1,1,1.,e);

  MFFreeNVector(du0,e);
  MFFreeNVector(du1,e);

  MFIMFProject(M,u,Phi,u0,e);

  H=MFCreateMultifariosMethod(e);
  MFMultifarioSetRealParameter(H,"epsilon",.1,e);            /* Max distance from TS to M */
  MFMultifarioSetIntegerParameter(H,"maxCharts",-1,e);       /* -1 means as many as needed */
  MFMultifarioSetIntegerParameter(H,"verbose",1,e);          /* Write info to stdout */
  MFMultifarioSetIntegerParameter(H,"page",1,e);             /* Page out non-interior polyhedra */
  MFMultifarioSetIntegerParameter(H,"dumpToPlotFile",1,e);   /* Write polyhedra to a plotfile */
  MFMultifarioSetIntegerParameter(H,"dumpToCenterFile",0,e); /* Write points to a file */
  MFMultifarioSetFilename(H,"PeriodFourResonance",e);

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
