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

static char *id="@(#) $Id: LorenzU0.c,v 1.4 2010/05/21 12:22:04 mhender Exp $";

#include <math.h>
#include <IMF.h>
#include <IMFFlow.h>
#include <IMFThreeDFlows.h>

int main(int argc, char *argv[])
 {
  IMFFlow L;
  MFAtlas A;
  MFNRegion Omega;
  MFNVector u0;
  MFNVector ll;
  MFNVector ur;
  char name[1024];
  int i;

  double eps=.5;
  double dt=.01;
/*double tmax=130.;*/
  double tmax=1000.;
  int maxInterp=10;
/*int maxInterp=50000;*/
/*int maxCharts=900000;*/
/*int maxCharts=300;*/
  int maxCharts=1000;
  MFErrorHandler e;

  e=MFCreateErrorHandler();

/*maxInterp=1;*/
  tmax=120.;
  if(argc>1)sscanf(argv[1],"%lf",&tmax);

  L=IMFCreateStandardLorenzFlow(e);

/*ll=MFIMFVectorFactory(L,e);*/
  ll=MFCreateNVector(3,e);
  MFNVSetC(ll,0,-70.,e);
  MFNVSetC(ll,1,-70.,e);
  MFNVSetC(ll,2,-30.,e);

/*ur=MFIMFVectorFactory(L,e);*/
  ur=MFCreateNVector(3,e);
  MFNVSetC(ur,0,70.,e);
  MFNVSetC(ur,1,70.,e);
  MFNVSetC(ur,2,100.,e);

/* For IMFSurvey */

  MFNVSetC(ll,0,-700.,e);
  MFNVSetC(ll,1,-700.,e);
  MFNVSetC(ll,2,-300.,e);
  MFNVSetC(ur,0,700.,e);
  MFNVSetC(ur,1,700.,e);
  MFNVSetC(ur,2,1000.,e);

/* end IMFSurvey */

  Omega=MFNRegionCreateHyperCubeByCorners(3,ll,ur,e);
  MFFreeNVector(ll,e);
  MFFreeNVector(ur,e);

/*u0=MFIMFVectorFactory(L,e);*/
  u0=MFCreateNVector(3,e);
  MFNVSetC(u0,0,0.,e);
  MFNVSetC(u0,1,0.,e);
  MFNVSetC(u0,2,0.,e);

  i=tmax;
  sprintf(name,"LorenzU0%3.3d",i);
  A=IMFComputeStableInvariantManifold(L,name,u0,(MFKVector)NULL,Omega,eps,dt,tmax,maxInterp,maxCharts,2.,e);

  MFFreeAtlas(A,e);
  MFFreeNVector(u0,e);
  MFFreeNRegion(Omega,e);
  IMFFreeFlow(L,e);

  MFFreeErrorHandler(e);

  return 0;
 }
