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
#include <MFNRegion.h>
#include <MFNVector.h>
#include <MFPrint.h>
#include <math.h>
#include <MFMultifariosMethod.h>

int MFTaylorAProjectToDraw(MFNVector,double*,void*,MFErrorHandler);

void findSolution(double,double,double,double,double,double,double*,double*,MFErrorHandler);

int main(int argc, char *argv[])
 {
  MFImplicitMF M;
  int i,j,n;
  int n0;
  MFNRegion Omega;
  MFAtlas S;
  MFNVector u0[20];
  MFNVector ug;
  MFNVector ll,ur;
  MFNKMatrix Tan;
  MFErrorHandler e;
  FILE *fid=(FILE*)NULL;
  char string[1024]="";
  char f1[1024]="";
  char b1[1024]="";
  char f2[1024]="";
  char b2[1024]="";

  double r0=78.53836;
  double g0=2.881799;
  double f1r=0.00117;
  double f1g=-0.0137;
  double f1rr=-0.00000854;
  double f1rg=-0.000407;
  double f1gg=0.00212;
  double a1=3.67;
  double b10=-0.0975;
  double b1r=-0.00392;
  double b1g=0.0543;
  double f2r=0.000681;
  double f2g=0.00955;
  double f2rr=-0.00000521;
  double f2rg=-0.000216;
  double f2gg=-0.00985;
  double a2=1.19;
  double b20=0.0331;
  double b2r=0.000476;
  double b2g=-0.00546;

  double x1=.1;
  double x2=-.1;
  double R=1.;
  double G=0.;
  double LR,Lg,Lx;
  double LEq1,LEq2;
  double vf1,vf2,vb1,vb2;
  MFNVector basis[2]={(MFNVector)NULL,(MFNVector)NULL};
  MFContinuationMethod H;

  double c0,c1,c2,c3,roots[3];
  int    nroots;

  e=MFCreateErrorHandler();

/* Type A mode interaction */

/* Scale R by 10 */

/*L=1./10.;
  f1r=L*f1r;
  f1rr=L*L*f1r;
  f1rg=L*f1rg;
  b1r=L*b1r;
  f2r=L*f2r;
  f2rr=L*L*f2r;
  f2rg=L*f2rg;
  b2r=L*b2r;*/

/* Scale R and g 10 */

  Lx=1./18.;
  LR=2.;
  Lg=.2;
  LR=.2*LR;
  Lg=.2*Lg;
  Lx=.1*Lx;
  LEq1=pow(Lx,3.);       LEq1=1.;
  LEq2=pow(Lx,3.);       LEq2=1.;

  printf("Lx=%lf, LR=%lf, Lg=%lf,    LEq1=%lf, LEq2=%lf\n",Lx,LR,Lg,LEq1,LEq2);fflush(stdout);

  f1r =LR   *f1r /Lx/Lx;
  f1g =Lg   *f1g /Lx/Lx;
  f1rr=LR*LR*f1rr/Lx/Lx;
  f1rg=LR*Lg*f1rg/Lx/Lx;
  f1gg=Lg*Lg*f1gg/Lx/Lx;

  b1r =LR*b1r/Lx;
  b1g =Lg*b1g/Lx;

  f2r =LR   *f2r /Lx/Lx;
  f2g =Lg   *f2g /Lx/Lx;
  f2rr=LR*LR*f2rr/Lx/Lx;
  f2rg=LR*Lg*f2rg/Lx/Lx;
  f2gg=Lg*Lg*f2gg/Lx/Lx;

  b2r =LR*b2r/Lx;
  b2g =Lg*b2g/Lx;

  sprintf(f1,"(%lf)*R+(%lf)*G+(%lf)*R*R+(%lf)*R*G+(%lf)*G*G",f1r,f1g,f1rr,f1rg,f1gg);
  sprintf(b1,"(%lf)+(%lf)*R+(%lf)*G",b10,b1r,b1g);
  sprintf(f2,"(%lf)*R+(%lf)*G+(%lf)*R*R+(%lf)*R*G+(%lf)*G*G",f2r,f2g,f2rr,f2rg,f2gg);
  sprintf(b2,"(%lf)+(%lf)*R+(%lf)*G",b20,b2r,b2g);

  sprintf(string,"[%lf*x1*(x1*x1+(%lf)*x2*x2-(%s)+(%s)*x2),%lf*(x2*((%lf)*x1*x1+x2*x2-(%s))+(%s)*x1*x1)]",LEq1,a1,f1,b1,LEq2,a2,f2,b2);

  printf("Expression is %s\n",string);fflush(stdout);
  M=MFIMFCreateAlgebraicExpressionWithRadius("[x1,x2,R,G]",string,.2,e); /* was .2 */
  MFIMFSetProjectForDraw(M,MFTaylorAProjectToDraw,e);

  n=MFIMF_N(M,e);
  ll=MFIMFVectorFactory(M,e);
  ur=MFIMFVectorFactory(M,e);
  MFNVSetC(ll,0,-20.,e);MFNVSetC(ur,0,20.,e);
  MFNVSetC(ll,1,-20.,e);MFNVSetC(ur,1,20.,e);
  MFNVSetC(ll,2,-5. ,e);MFNVSetC(ur,2,5.,e);
  MFNVSetC(ll,3,-5. ,e);MFNVSetC(ur,3,5.,e); 
  Omega=MFNRegionCreateHyperCubeByCorners(n,ll,ur,e);
  MFFreeNVector(ll,e);
  MFFreeNVector(ur,e);
  

  H=MFCreateMultifariosMethod(e);
  MFMultifarioSetRealParameter(H,"epsilon",.01,e);
  MFMultifarioSetIntegerParameter(H,"maxCharts",-1,e);
  MFMultifarioSetIntegerParameter(H,"verbose",1,e);
  MFMultifarioSetIntegerParameter(H,"page",1,e);
  MFMultifarioSetIntegerParameter(H,"dumpToPlotFile",1,e);
  MFMultifarioSetIntegerParameter(H,"dumpToCenterFile",0,e);
  MFMultifarioSetFilename(H,"Circle",e);
  MFMultifarioSetIntegerParameter(H,"branchSwitch",100,e);
  MFMultifarioSetIntegerParameter(H,"dumpToPlotFile",1,e);
  MFMultifarioSetIntegerParameter(H,"dumpToRestartFile",1,e);
  MFMultifarioSetIntegerParameter(H,"dumpToRestartFileEvery",10,e);
  MFMultifarioSetRealParameter(H,"dotmin",.9,e);

  for(i=0;i<20;i++)u0[i]=MFIMFVectorFactory(M,e);

/* Trivial */

  n0=0;
  MFNVSetC(u0[n0],0,0.,e);
  MFNVSetC(u0[n0],1,0.,e);
  MFNVSetC(u0[n0],2,-4.5,e);
  MFNVSetC(u0[n0],3,4.5,e);
  n0++;
  MFMultifarioSetFilename(H,"TaylorA",e);

  printf("There are %d initial solutions\n",n0);fflush(stdout);
  n0=1.;
  S=MFComputeAtlasMultiple(H,M,Omega,n0,u0,e);

  MFCloseAtlas(H,S,e);
  MFFreeAtlas(S,e);
  MFFreeImplicitMF(M,e);
  MFFreeNRegion(Omega,e);
  for(i=0;i<20;i++)MFFreeNVector(u0[i],e);
  MFFreeContinuationMethod(H,e);
  MFFreeErrorHandler(e);

  return(0);
 }

int MFTaylorAProjectToDraw(MFNVector u, double *x, void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFTaylorAProjectToDraw"};

  if(x==(double*)NULL)return 3;

  x[0]=MFNV_C(u,3,e);
  x[1]=MFNV_C(u,2,e);
  x[2]=.5*(MFNV_C(u,0,e)+MFNV_C(u,1,e));

  return 0;
 }
