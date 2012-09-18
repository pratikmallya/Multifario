/*
    %W%
    %D% %T%

*/
/*      author: Mike Henderson mhender@watson.ibm.com */
/*      date:   July 18, 2002 Modified AUTO's main program       */

#include <multifarioConfig.h>
#include <MFAUTO.h>
#include <MFMultifariosMethod.h>
int MFNVectorGetNRefs(MFNVector);

/* AUTBVP */

int f (integer ndim, const doublereal *u, const integer *icp,
          const doublereal *par, integer ijac,
          doublereal *f, doublereal *dfdu, doublereal *dfdp) {
  doublereal d1;
  doublereal pi;

  pi = atan(1.) * 4;
  f[0] = u[1];
  d1 = par[0] * pi;

  f[1] = -(d1 * d1) * u[0] + u[0] * u[0];

  if(ijac)
   {
    dfdu[0+2*0]=0.;            dfdu[0+2*1]=1.;
    dfdu[1+2*0]=-d1*d1+2*u[0]; dfdu[1+2*1]=0.;
    dfdp[0]=0.;                dfdp[1]=-2*d1*u[0];
   }

  return 0;
}
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
int stpnt (integer ndim, doublereal t,
           doublereal *u, doublereal *par) {
  par[0] = 0.;
  u[0] = (doublereal)0.;
  u[1] = (doublereal)0.;

  return 0;
}
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
int a (integer ndim, const doublereal *par, const integer *icp,
          integer nbc, const doublereal *u0, const doublereal *u1, integer ijac,
          doublereal *fb, doublereal *dfb)
 {
  fb[0] = u0[0];
  fb[1] = u1[0];

  if(ijac)
   {
/* u0[0] u0[1] u1[0] u1[1] par[0] */
    dfb[0+2*0]=1.;  dfb[0+2*1]=0.; dfb[0+2*2]=0.;  dfb[0+2*3]=0.; dfb[0+2*4]=0.;
    dfb[1+2*0]=0.;  dfb[1+2*1]=0.; dfb[1+2*2]=1.;  dfb[1+2*3]=0.; dfb[1+2*4]=0.;
   }

  return 0;
}
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
int ic(integer ndim, const doublereal *par, const integer *icp,
          integer nint, const doublereal *u, const doublereal *uold,
          const doublereal *udot, const doublereal *upold, integer ijac,
          doublereal *fi, doublereal *dint) {
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
  MFContinuationMethod H;
  MFContinuationMethod A;
  MFNVector u0;
  MFNKMatrix Phi0;
  MFAtlas S;

  double epsilon=1.e-6;
  int through=0;
  doublereal rl0[1]={0.};
  doublereal rl1[1]={5.};
  doublereal a0=-50.;
  doublereal a1=50.;
  doublereal ds=.2;

  integer npar=1;
  integer nicp=1;
  integer icp[1]={0};
  doublereal par[1];

  doublereal thu[2]={1.,1.};
  doublereal thl[1]={1.};

  integer k=1;
  integer jac=1;
  integer ntst=5,ncol=4,ndim=2;
  integer nbc=2,nic=0;

  MFErrorHandler e;

  e=MFCreateErrorHandler();

  tpbvp=MFCreateAUTOTPBVP(k,ndim,f,jac,nbc,a,nic,ic,npar,nicp,icp,ntst,ncol,plt,e);

  space=MFCreateAUTONSpace(tpbvp,thu,thl,e);
  M=MFCreateAUTOBV(tpbvp,space,e);
  MFIMFSetR(M,ds,e);
  Omega=MFNRegionCreateAUTO(space,npar,icp,rl0,rl1,a0,a1,e);

  A=MFCreateAUTOsMethod(e);

  H=MFCreateMultifariosMethod(e);
  MFMultifarioSetRealParameter(H,"epsilon",epsilon,e);
  MFMultifarioSetIntegerParameter(H,"branchSwitch",through,e);
  MFMultifarioSetIntegerParameter(H,"maxCharts",100,e);
  MFMultifarioSetIntegerParameter(H,"verbose",1,e);
  MFMultifarioSetIntegerParameter(H,"page",1,e);
  MFMultifarioSetIntegerParameter(H,"dumpToPlotFile",1,e);
  MFMultifarioSetIntegerParameter(H,"dumpToCenterFile",1,e);
  MFMultifarioSetFilename(H,"AUTBVP",e);

/* Add par */
  par[0]=0.;
printf("get start point\n");fflush(stdout);
  MFAUTOGetStartPoint(M,tpbvp,stpnt,par,&u0,&Phi0,e);

printf("compute atlas\n");fflush(stdout);
  if(Phi0!=(MFNKMatrix)NULL)S=MFComputeAtlasWithTangent(H,M,Omega,u0,Phi0,e);
    else S=MFComputeAtlas(H,M,Omega,u0,e);

  MFCloseAtlas(H,S,e);
  MFFreeNRegion(Omega,e);
  MFFreeNSpace(space,e);
  MFFreeNVector(u0,e);
  MFFreeImplicitMF(M,e);
  MFFreeContinuationMethod(H,e);
  MFFreeContinuationMethod(A,e);
  MFFreeErrorHandler(e);

  return 0;
} 
