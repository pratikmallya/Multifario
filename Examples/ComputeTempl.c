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

static char *id="@(#) $Id: ComputeTempl.c,v 1.2 2007/06/08 21:01:20 mhender Exp $";

#include <MFAtlas.h>
#include <MFMultifariosMethod.h>

MFImplicitMF MFIMFCreateTempl(double stiff, MFErrorHandler e);

int main(int argc, char *argv[])
 {
  MFImplicitMF M;
  int n;
  MFNRegion Omega;
  MFAtlas S;
  MFNVector u0;
  MFNVector minn,maxx;
  FILE *fid;
  /*double Rxy=.8;
    double Rz=.5;*/
  double stiff=15.;
  MFContinuationMethod H;
  MFErrorHandler e;

  /*printf("Templ, xy R=%lf, z R=%lf (surface area=%lf)\n",Rxy,Rz,4*3.1415926*3.1415926*Rxy*Rz);fflush(stdout);*/
  printf("Templ, stiffness=%lf\n",stiff);fflush(stdout);

  e=MFCreateErrorHandler();

  /*M=MFIMFCreateTempl(0.,0.,0.,Rxy,Rz,e);*/
  M=MFIMFCreateTempl(stiff,e);
  n=MFIMF_N(M,e);
  printf("Templ, n=%lf why is this 0 ? \n",n);fflush(stdout);
  /*Omega=MFNRegionCreateCube(-5.,-5.,-5.,5.,25.,5.,e);*/
  minn=MFIMFVectorFactory(M,e);
  maxx=MFIMFVectorFactory(M,e);
  /*changed some of these lines*/
  MFNVSetC(minn,0,-5.,e);
  MFNVSetC(minn,1,-5.,e);
  MFNVSetC(minn,2,-5.,e);
  /*MFNVSetC(minn,3,-5.,e);*/
  MFNVSetC(maxx,0,5.,e);
  MFNVSetC(maxx,1,25.,e);
  MFNVSetC(maxx,2,5.,e);
  /*MFNVSetC(maxx,3,5.,e);*/
  Omega=MFNRegionCreateHyperCubeByCorners(n,minn,maxx,e);
  /*Omega=MFNRegionCreateRectangle(-5.,-5.,5.,30.,e);*/
  printf("Vector dim is %lf why is this wrong? \n",MFNV_NC(maxx,e));fflush(stdout);
  printf("Vector comp 2 is %lf this is correct \n",MFNV_C(maxx,2,e));fflush(stdout);

  u0=MFIMFVectorFactory(M,e);
  MFNVSetC(u0,0, 0.,e);
  MFNVSetC(u0,1, 0.,e);
  MFNVSetC(u0,2, 0.,e);
  /*MFNVSetC(u0,3, 0.,e);*/

  H=MFCreateMultifariosMethod(e);
  MFMultifarioSetRealParameter(H,"epsilon",.01,e);
  MFMultifarioSetIntegerParameter(H,"maxCharts",-1,e);
  MFMultifarioSetIntegerParameter(H,"verbose",1,e);
  MFMultifarioSetIntegerParameter(H,"page",1,e);
  MFMultifarioSetIntegerParameter(H,"dumpToPlotFile",1,e);
  MFMultifarioSetIntegerParameter(H,"dumpToCenterFile",0,e);
  MFMultifarioSetFilename(H,"Templ",e);

  S=MFComputeAtlas(H,M,Omega,u0,e);

  MFCloseAtlas(H,S,e);
  printf("Done computing Atlas\n");fflush(stdout);

  MFFreeAtlas(S,e);
  MFFreeImplicitMF(M,e);
  MFFreeNRegion(Omega,e);
  MFFreeNVector(u0,e);
  MFFreeNVector(minn,e);
  MFFreeNVector(maxx,e);
  MFFreeContinuationMethod(H,e);

  MFFreeErrorHandler(e);

  return 0;
 }
