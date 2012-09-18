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

static char *id="@(#) $Id: AUTINT.c,v 1.6 2007/07/06 12:38:35 mhender Exp $";

#include <multifarioConfig.h>
#include <MFAUTO.h>

int MFNVectorGetNRefs(MFNVector);

/* AUTINT, int - Braatu's equation with boundary condition and integral constraint */

int f(integer ndim, const doublereal *u, const integer *icp,
         const doublereal *par, integer ijac,
         doublereal *f, doublereal *dfdu, doublereal *dfdp)
 {
  doublereal e;

  e=exp(u[0]);

  f[0]=u[1];
  f[1]=-par[0]*e;

  if(ijac==0)return;

  dfdu[0+2*0]=0.;
  dfdu[1+2*0]=-par[0]*e;

  dfdu[0+2*1]=1.;
  dfdu[1+2*1]=0.;

  dfdp[0+1*0]=0.;
  dfdp[1+1*0]=-e;

  return 0;
 }

int a(integer ndim, const doublereal *par, const integer *icp, integer nbc,
         const doublereal *u0, const doublereal *u1, integer ijac,
         doublereal *fb, doublereal *dbc)
 {
  fb[0]=u0[0]-u1[0]-par[1];

  if(ijac==0)return;

  dbc[0+1*0]=1.;
  dbc[0+1*1]=0.;

  dbc[0+1*2]=-1.;
  dbc[0+1*3]= 0.;

  dbc[0+1*4]= 0.;
  dbc[0+1*5]=-1.;
  dbc[0+1*6]= 0.;

  return 0;
 }

int ic(integer ndim, const doublereal *par, const integer *icp, integer nint,
         const doublereal *u, const doublereal *uold, const doublereal *udot,
         const doublereal *upold, integer ijac,
         doublereal *fi, doublereal *dint)
 {
  fi[0]=u[0]-par[2];

  if(ijac==0)return 0;

  dint[0+1*0]=1.;
  dint[0+1*1]=0.;

  dint[0+1*2]= 0.;
  dint[0+1*3]= 0.;
  dint[0+1*4]=-1.;

  return 0;
 }

int stpnt(integer ndim, doublereal t, doublereal *u, doublereal *par)
 {
  par[0]=0.;
  par[1]=0.;
  par[2]=0.;

  u[0]=0.;
  u[1]=0.;

  return 0;
 }

int plt(integer ndim, const void *u, doublereal *par)
 {
/* If algebraically defined parameters are included this does the projection */
  return 0;
 }

int main(int argc,char *argv[])
{
  MFAUTOTPBVP tpbvp;
  MFNSpace space;
  MFNRegion Omega;
  MFImplicitMF M;
  MFContinuationMethod A;
  MFNVector u0;
  MFNKMatrix Phi0;
  MFAtlas S;
  MFErrorHandler e;

  doublereal rl0[1]={0.};
  doublereal rl1[1]={20.};
  doublereal a0=0.;
  doublereal a1=100.;
  doublereal ds=0.1;

  integer npar=3;
  integer nicp=1;
  integer icp[1]={0};
  doublereal par[3];

  doublereal thu[2]={1.,1.};
  doublereal thl[3]={1.,1.,1.};

  integer k=1;
  integer jac=1;
  integer ntst=5,ncol=4,ndim=2;
  integer nbc=1,nic=1;

  e=MFCreateErrorHandler();

  tpbvp=MFCreateAUTOTPBVP(k,ndim,f,jac,nbc,a,nic,ic,npar,nicp,icp,ntst,ncol,plt,e);

  space=MFCreateAUTONSpace(tpbvp,thu,thl,e);
  M=MFCreateAUTOBV(tpbvp,space,e);
  MFIMFSetR(M,ds,e);

  MFAUTOBVSetRealParameter(M,"dlmin",0.0001,e);
  MFAUTOBVSetRealParameter(M,"dsmax",2.,e);
  MFAUTOBVSetRealParameter(M,"ds0",0.01,e);

  MFAUTOAddUserZero(M,0,1.,e);
  MFAUTOAddUserZero(M,0,3.,e);
  MFAUTODetectLimitPoints(M,e);

  Omega=MFNRegionCreateAUTO(space,k,icp,rl0,rl1,a0,a1,e);

  A=MFCreateAUTOsMethod(e);

  par[0]=0.;
  par[1]=0.;
  par[2]=0.;
  MFAUTOGetStartPoint(M,tpbvp,stpnt,par,&u0,&Phi0,e);

  S=MFComputeAtlasWithTangent(A,M,Omega,u0,Phi0,e);
  MFCloseAtlas(A,S,e);

  MFFreeNRegion(Omega,e);
  MFFreeNSpace(space,e);
  MFFreeNVector(u0,e);
  if(Phi0!=(MFNKMatrix)NULL)MFFreeNKMatrix(Phi0,e);
  MFFreeImplicitMF(M,e);
  MFFreeContinuationMethod(A,e);
  MFFreeAtlas(S,e);

  MFFreeErrorHandler(e);

  return 0;
} 
