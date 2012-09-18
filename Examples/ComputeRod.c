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

static char *id="@(#) $Id: ComputeRod.c,v 1.6 2011/07/21 17:43:45 mhender Exp $";

#include <MFAtlas.h>
#include <MFNRegion.h>
#include <MFNVector.h>
#include <MFDX.h>
#include <MFFortran.h>
#include <MFTPBVP.h>
#include <math.h>
#include <MFMultifariosMethod.h>
static char ComputeRodErrorHandlerMsg[256]="";

#define round(x) ((x)>=0?(long)((x)+0.5):(long)((x)-0.5))

int MFRodProjectToDraw(MFNVector,double*,void*,MFErrorHandler);

double MzPositive(MFNVector,MFErrorHandler);
double MzNegative(MFNVector,MFErrorHandler);
double MzSmallPositive(MFNVector,MFErrorHandler);
double MzSmallNegative(MFNVector,MFErrorHandler);
double Mx0Positive(MFNVector,MFErrorHandler);
double Mx0Negative(MFNVector,MFErrorHandler);
double FPositive(MFNVector,MFErrorHandler);
double DNegative(MFNVector,MFErrorHandler);
double DPositive(MFNVector,MFErrorHandler);
double FLessThanSixteen(MFNVector,MFErrorHandler);
double DSmallNegative(MFNVector,MFErrorHandler);
double DSmallPositive(MFNVector,MFErrorHandler);
double DSmallerNegative(MFNVector,MFErrorHandler);
double DSmallerPositive(MFNVector,MFErrorHandler);

double K(double);
double Ki(double);

double MFTPBVPTestTangent(MFImplicitMF,MFNVector,MFNKMatrix,MFErrorHandler);
MFNKMatrix ComputeTangent(MFImplicitMF,int,int,int,int,MFNVector,MFErrorHandler);
int MFTPBVPAnalyze(MFImplicitMF,MFNVector,MFErrorHandler);
void MFTPBVPNSpaceSetPeriodicParameter(MFNSpace,int,double,MFErrorHandler);
int MFRodStopN(MFImplicitMF,MFNVector,MFNKMatrix,MFNVector,MFNKMatrix,void*,MFErrorHandler);
int MFRodStopD(MFImplicitMF,MFNVector,MFNKMatrix,MFNVector,MFNKMatrix,void*,MFErrorHandler);
int MFRodProjectToDraw(MFNVector,double*,void*,MFErrorHandler);
MFNVector GetInitialPoint(MFImplicitMF,double,double,double,double,double,MFErrorHandler);
MFNVector GetInitialPointMx0Mz(MFImplicitMF,double,double,double,double,MFErrorHandler);
MFNVector GetInitialPointMx0f(MFImplicitMF,double,double,double,double,MFErrorHandler);
MFNVector GetInitialPointI(MFImplicitMF,double,double,double,double,MFErrorHandler);
void flipf(int,int,MFNVector,MFErrorHandler);
void flipMz(int,int,MFNVector,MFErrorHandler);
void flipMx0(int,int,MFNVector,MFErrorHandler);
void flipMx0AndMz(int,int,MFNVector,MFErrorHandler);
void moveToPrinciple(int,int,MFNVector,MFErrorHandler);
void moveToPositive(int,int,MFNVector,MFErrorHandler);

void MFTPBVPEvaluateIntegralConstraints(MFImplicitMF,MFNVector,MFNVector,double*, MFErrorHandler);
void MFTPBVPEvaluateBoundaryConditions(MFImplicitMF,MFNVector,MFNVector,double*,MFErrorHandler);

int GetNullVector(int,int,int,int,double*,double*,double*,MFErrorHandler);

struct MFTPBVPData
 {
  int nx;
  int nu;
  int np;
  int nbc;
  int nic;
  int k;
  MFTPBVPFFUNCTION f;
  MFTPBVPFFUNCTION fu;
  MFTPBVPFFUNCTION fl;
  MFTPBVPAFUNCTION a;
  MFTPBVPAFUNCTION au;
  MFTPBVPAFUNCTION al;
  MFTPBVPLFUNCTION l;
  MFTPBVPLFUNCTION lu;
  MFTPBVPLFUNCTION ll;
  MFTPBVPMFUNCTION m;
  MFTPBVPMFUNCTION ml;
  MFNSpace space;
 };

#define PI 3.14159265358979323846264338327950288

/*
u'=f(u,t,l);                     nu differential equations
a(u(0),u(1),l)=0;                nbc=nu equations
int l(u(t),t,l) dt + m(l) = 0;   nc equations
*/

/* mz  a - p[0] */
/* mx0 c - p[1] */
/* f   t - p[2] */
/*     theta0 - p[3] */

void f(double r, int nu, double *u, int np, double *p, double *u0, double *l0, double *f, MFErrorHandler e)
 {
  f[0]=PI*u[1];
  f[1]=PI*(p[2]*u[0]*u[5]-p[0]*u[3]);
  f[2]=PI*u[3];
  f[3]=PI*(p[2]*u[2]*u[5]-p[1]*u[5]+p[0]*u[1]);
  f[4]=PI*u[5];
  f[5]=PI*(-p[2]*u[2]*u[3]-p[2]*u[0]*u[1]+p[1]*u[3]);

  return;
 }

void fu(double r, int nu, double *u, int np, double *p, double *u0, double *l0, double *fu, MFErrorHandler e)
 {
  int i;

  for(i=0;i<nu*nu;i++)fu[i]=0.;

  fu[0+nu*1]=PI;

  fu[1+nu*0]=PI*p[2]*u[5];
  fu[1+nu*3]=-PI*p[0];
  fu[1+nu*5]=PI*p[2]*u[0];

  fu[2+nu*3]=PI;

  fu[3+nu*1]=PI*p[0];
  fu[3+nu*2]=PI*p[2]*u[5];
  fu[3+nu*5]=PI*(p[2]*u[2]-p[1]);

  fu[4+nu*5]=PI;

  fu[5+nu*0]=-PI*p[2]*u[1];
  fu[5+nu*1]=-PI*p[2]*u[0];
  fu[5+nu*2]=-PI*p[2]*u[3];
  fu[5+nu*3]=PI*(-p[2]*u[2]+p[1]);

  return;
 }

void fl(double r, int nu, double *u, int np, double *p, double *u0, double *l0, double *fl, MFErrorHandler e)
 {
  int i;

  for(i=0;i<nu*np;i++)fl[i]=0.;

  fl[1+nu*0]=-PI*u[3];
  fl[1+nu*2]=PI*u[0]*u[5];
  fl[3+nu*0]=PI*u[1];
  fl[3+nu*1]=-PI*u[5];
  fl[3+nu*2]=PI*u[2]*u[5];
  fl[5+nu*1]=PI*u[3];
  fl[5+nu*2]=-PI*(u[2]*u[3]+u[0]*u[1]);

  return;
 }

void a(int nbc, int nu, double *uL, double *uR, int np, double *p, double *u0L, double *u0R,double *l0, double *a, MFErrorHandler e)
 {
  a[0]=uL[0];
  a[1]=uL[1]-sin(p[3]);
  a[2]=uL[2];
  a[3]=uL[3];
  a[4]=uL[4];
  a[5]=uL[5]-cos(p[3]);
  a[6]=uR[3];
  a[7]=uR[0]*uR[5]-uR[4]*uR[1];

  return;
 }

void au(int nbc, int nu, double *uL, double *uR, int np, double *p, double *u0L, double *u0R,double *l0, double *au, MFErrorHandler e)
 {
  int i;

  for(i=0;i<nbc*2*nu;i++)au[i]=0.;

  au[0+nbc*0]=1.;
  au[1+nbc*1]=1.;
  au[2+nbc*2]=1.;
  au[3+nbc*3]=1.;
  au[4+nbc*4]=1.;
  au[5+nbc*5]=1.;
  au[6+nbc*(6+3)]=1.;
  au[7+nbc*(6+0)]=uR[5];
  au[7+nbc*(6+1)]=-uR[4];
  au[7+nbc*(6+4)]=-uR[1];
  au[7+nbc*(6+5)]=uR[0];

  return;
 }

void al(int nbc, int nu, double *uL, double *uR, int np, double *p, double *u0L, double *u0R,double *l0, double *al, MFErrorHandler e)
 {
  int i;

  for(i=0;i<nbc*np;i++)al[i]=0.;

  al[1+nbc*3]=-cos(p[3]);
  al[5+nbc*3]= sin(p[3]);

  return;
 }

void l(int nic, double r, int nu, double *u, int np, double *p,double *u0,double *l0, double *l, MFErrorHandler e)
 {
  int i;

  for(i=0;i<nic;i++)l[i]=0.;

  if(nic>0)
    l[0]=u[1]*u[1]+u[3]*u[3]+u[5]*u[5];
  if(nic>1)
    l[1]=u[0]*u[0]+u[2]*(u[2]-p[1])+u[5];

  return;
 }

void lu(int nic, double r, int nu, double *u, int np, double *p,double *u0,double *l0, double *lu, MFErrorHandler e)
 {
  int i;

  for(i=0;i<nu*nic;i++)lu[i]=0.;

  if(nic>0)
   {
    lu[0+nic*1]=2*u[1];
    lu[0+nic*3]=2*u[3];
    lu[0+nic*5]=2*u[5];
   }
  if(nic>1)
   {
    lu[1+nic*0]=2*u[0];
    lu[1+nic*2]=2*u[2]-p[1];
    lu[1+nic*5]=1.;
   }

  return;
 }

void ll(int nic, double r, int nu, double *u, int np, double *p,double *u0,double *l0, double *ll, MFErrorHandler e)
 {
  int i;

  for(i=0;i<np*nic;i++)ll[i]=0.;
  if(nic>0)ll[1+nic*8]=-u[2];

  return;
 }

void m(int nic, int np, double *p, double *p0, double *m, MFErrorHandler e)
 {
  int i;

  for(i=0;i<nic;i++)m[i]=0.;
  if(nic>0)m[0]=-pow(p[2],4);
  if(nic>1)m[1]=-p[3];

  return;
 }

void ml(int nic, int np, double *p, double *p0, double *ml, MFErrorHandler e)
 {
  int i;

  for(i=0;i<np*nic;i++)ml[i]=0.;

  if(nic>0)ml[0+nic*2]=-4.*pow(p[2],3);
  if(nic>1)ml[1+nic*3]=-1.;

  return;
 }

#define U0I {u0[n0]=GetInitialPointI(M,Mx0,Mz,F,theta0,e);if(u0[n0]!=NULL){printf("u0[%d] is (%lf,%lf,%lf,%lf)\n",n0,MFNV_C(u0[n0],nx*nu+0,e),MFNV_C(u0[n0],nx*nu+1,e),MFNV_C(u0[n0],nx*nu+2,e),MFNV_C(u0[n0],nx*nu+3,e));n0++;}}

#define U0Mx0Mz {u0[n0]=GetInitialPointMx0Mz(M,Mx0,Mz,F,theta0,e);if(u0[n0]!=NULL){printf("u0[%d] is (%lf,%lf,%lf,%lf)\n",n0,MFNV_C(u0[n0],nx*nu+0,e),MFNV_C(u0[n0],nx*nu+1,e),MFNV_C(u0[n0],nx*nu+2,e),MFNV_C(u0[n0],nx*nu+3,e));n0++;}}

#define FLIPMz {u0[n0]=MFCloneNVector(u0[n0-1],e);flipMz(nx,nu,u0[n0],e);n0++;}
#define FLIPMx0 {u0[n0]=MFCloneNVector(u0[n0-1],e);flipMx0(nx,nu,u0[n0],e);n0++;}

int main(int argc, char *argv[])
 {
  MFImplicitMF M;
  int i,j,n;
  MFNRegion Omega;
  MFAtlas S;
  MFNVector u0[500];
  int n0;
  FILE *fid;
  MFNVector ug;
  MFNKMatrix Tan;

  double C,A,t;
  double Mx0,Mz,F,theta0;
  double r,dr,x1,y1,x2,y2,x3,y3;
  double xy0[6],p0[4],p1[4];
  double error,e0,e1,corr;
  double dA,dC,dy1,dy3,a00,a01,a10,a11,de0,de1;
  int itimes;
  MFNVector tug;
  double eps;

  int nx=150;
  int np=4;
  int nu=6;
  int nbc=8;
  int nic=0;
  int k;
  MFContinuationMethod H;
  char name[80];
  MFErrorHandler e;

  e=MFCreateErrorHandler();

/* ( a,-c, t, th) <-> ( a, c, t,2PI-th) */
/* (-a, c, t, th) <-> ( a, c, t,2PI-th) */
/* ( a, c,-t, th) <-> ( a, c, t,th-PI) (x,y,z,dx,dy,dz)-> (-x,-y,-z,-dx,-dy,-dz) */

  k=nu+np-nbc-nic;
  
  M=MFIMFCreateTPBVP(k,nx,nu,np,f,fu,fl,nbc,a,au,al,nic,l,lu,ll,m,ml,e);
  MFTPBVPNSpaceSetPeriodicParameter(MFIMFNSpace(M,e),3,2*PI,e);
  MFIMFSetStop(M,MFRodStopN,e);
  MFIMFSetR(M,.3,e);
  MFIMFSetProjectForDraw(M,MFRodProjectToDraw,e);
  n=MFIMF_N(M,e);
  printf("TPBVP calculated k to be %d\n",MFIMF_K(M,e));

  p0[0]= -0.05;   /* Mz .1  */
  p0[1]= -10.5;   /* Mx0.1  */
  p0[2]= -.1;   /* f  */
  p0[3]= -1*PI;  /* th */

  p1[0]= 10.5;  /* a  Mz */
  p1[1]=  0.05;   /* c  Mx0 */
  p1[2]= 16.;   /* 16 f  */
  p1[3]= 1*PI;  /* th */
  Omega=MFNRegionCreateTPBVP(nx,nu,np,p0,p1,-200.,200.,e);

/* 2- Mz<0, Mx0<0 */
  n0=0;
  Mx0=-1.5483;Mz= 1.9552;F=4.1227;theta0=-1.4082;U0Mx0Mz
  Mx0=-0.2279;Mz= 1.2260;F=4.5294;theta0=-0.9955;U0Mx0Mz
  printf("Number of converged initial points = %d\n",n0);fflush(stdout);

  H=MFCreateMultifariosMethod(e);
  MFMultifarioSetRealParameter(H,"epsilon",.1*60,e); /* This looks strange, but it's because of the 6 eqs, on intervals of length PI */
  MFMultifarioSetIntegerParameter(H,"maxCharts",-1,e);
  MFMultifarioSetIntegerParameter(H,"verbose",1,e);
  MFMultifarioSetIntegerParameter(H,"page",1,e);
  MFMultifarioSetIntegerParameter(H,"pageEvery",100,e);
  MFMultifarioSetIntegerParameter(H,"branchSwitch",0,e);
  MFMultifarioSetIntegerParameter(H,"dumpToPlotFile",1,e);
  MFMultifarioSetIntegerParameter(H,"dumpToCenterFile",0,e);
  MFMultifarioSetFilename(H,"Rod",e);

  MFMultifarioAddClipF(H,MzPositive,e);
  MFMultifarioAddClipF(H,Mx0Negative,e);
  MFMultifarioAddClipF(H,FPositive,e);
  MFMultifarioAddClipF(H,DNegative,e);
  MFMultifarioAddClipF(H,FLessThanSixteen,e);

  S=MFComputeAtlasMultiple(H,M,Omega,n0,u0,e);
  MFFlushAtlas(H,S,e);
  MFFreeNRegion(Omega,e);
  for(i=0;i<n0;i++)MFFreeNVector(u0[i],e);
  p0[0]= -10.5;   /* Mz .1  */
  p0[1]=  0.05;   /* Mx0.1  */
  p0[2]= -.1;   /* f  */
  p0[2]=  .05;   /* f  */
  p0[3]= -4*PI;  /* th */

  p1[0]= -.05;  /* a  Mz */
  p1[1]= 10.5;   /* c  Mx0 */
  p1[2]= 16.;   /* 16 f  */
  p1[3]= 4*PI;  /* th */
  Omega=MFNRegionCreateTPBVP(nx,nu,np,p0,p1,-200.,200.,e);

/* 2- Mz<0, Mx0>0 */
  n0=0;
  Mx0=1.4595;Mz=-0.2071;F=3.5751;theta0=-1.0814;U0Mx0Mz
  Mx0=1.3491;Mz=-0.0881;F=3.3506;theta0=-1.1738;U0Mx0Mz
/* 4- Mz<0, Mx0>0 */
  Mx0=1.4294;Mz=-0.2155;F=3.9047;theta0=-0.9269;U0Mx0Mz
  printf("Number of converged initial points = %d\n",n0,e);fflush(stdout);

  MFMultifariosMethodClearClipF(H,e);
  MFMultifarioAddClipF(H,MzNegative,e);
  MFMultifarioAddClipF(H,Mx0Positive,e);
  MFMultifarioAddClipF(H,DNegative,e);
  MFMultifarioAddClipF(H,FPositive,e);
  MFMultifarioAddClipF(H,FLessThanSixteen,e);

  MFExtendAtlasMultiple(S,H,M,Omega,n0,u0,e);
  MFFlushAtlas(H,S,e);
  MFFreeNRegion(Omega,e);
  for(i=0;i<n0;i++)MFFreeNVector(u0[i],e);

  p0[0]= -.15;   /* Mz .1  */
  p0[1]= -0.01;   /* Mx0.1  */
  p0[2]= -.15;   /* f  */
  p0[2]=  .05;   /* f  */
  p0[3]= -4*PI;  /* th */

  p1[0]=  .15;  /* a  Mz */
  p1[1]= 10.5;   /* c  Mx0 */
  p1[2]= 16.;   /* 16 f  */
  p1[3]= 4*PI;  /* th */
  Omega=MFNRegionCreateTPBVP(nx,nu,np,p0,p1,-200.,200.,e);
/* 2- Mz=0, Mx0>0 */
  n0=0;
  Mx0=1.6306;Mz=-0.0109;F=1.8838;theta0=-0.4894;U0Mx0Mz

  printf("Number of converged initial points = %d\n",n0);fflush(stdout);

  MFMultifariosMethodClearClipF(H,e);
  MFMultifarioAddClipF(H,MzSmallPositive,e);
  MFMultifarioAddClipF(H,MzSmallNegative,e);
  MFMultifarioAddClipF(H,Mx0Positive,e);
  MFMultifarioAddClipF(H,FPositive,e);
  MFMultifarioAddClipF(H,DNegative,e);
  MFMultifarioAddClipF(H,FLessThanSixteen,e);

  MFExtendAtlasMultiple(S,H,M,Omega,n0,u0,e);
  MFFlushAtlas(H,S,e);
  MFFreeNRegion(Omega,e);
  for(i=0;i<n0;i++)MFFreeNVector(u0[i],e);

  p0[0]= -10.5;   /* Mz .1  */
  p0[1]= -0.05;   /* Mx0.1  */
  p0[2]= -.1;   /* f  */
  p0[2]=  .05;   /* f  */
  p0[3]= -4*PI;  /* th */

  p1[0]= -.05;  /* a  Mz */
  p1[1]= 10.5;   /* c  Mx0 */
  p1[2]= 16.;   /* 16 f  */
  p1[3]= 4*PI;  /* th */
  Omega=MFNRegionCreateTPBVP(nx,nu,np,p0,p1,-200.,200.,e);

/* 2- Mz<0, Mx0>0 */
  n0=0;
  Mx0=-1.5436;Mz= 2.1364;F= 0.5304;theta0=-.3170;U0Mx0Mz
  flipMx0AndMz(nx,nu,u0[n0-1],e);

  printf("Number of converged initial points = %d\n",n0);fflush(stdout);

  MFMultifariosMethodClearClipF(H,e);
  MFMultifarioAddClipF(H,MzNegative,e);
  MFMultifarioAddClipF(H,Mx0Positive,e);
  MFMultifarioAddClipF(H,DNegative,e);
  MFMultifarioAddClipF(H,FPositive,e);
  MFMultifarioAddClipF(H,FLessThanSixteen,e);

  MFExtendAtlasMultiple(S,H,M,Omega,n0,u0,e);
  MFFlushAtlas(H,S,e);
  MFFreeNRegion(Omega,e);
  for(i=0;i<n0;i++)MFFreeNVector(u0[i],e);
  p0[0]=  0.05;   /* Mz .1  */
  p0[1]= -0.1;   /* Mx0.1  */
  p0[2]= -.1;   /* f  */
  p0[2]=  .05;   /* f  */
  p0[3]= -2*PI;  /* th */

  p1[0]=  10.5;   /* a  Mz */
  p1[1]=  10.5;   /* c  Mx0 */
  p1[2]= 16.;   /* 16 f  */
  p1[3]= 2*PI;  /* th */
  Omega=MFNRegionCreateTPBVP(nx,nu,np,p0,p1,-200.,200.,e);

/* 2- Buckled Ring */
  n0=0;

/* 4-  */
  Mx0= 0.5276;Mz= 0.8750;F= 3.1051;theta0=-1.3618;U0Mx0Mz
/* 4-  */
  Mx0= 0.2372;Mz= 0.9332;F= 4.0547;theta0=-1.0424;U0Mx0Mz
/* 4-  (new this trip)*/
  Mx0= 0.7051;Mz= 0.5268;F= 3.7947;theta0=-1.0477;U0Mx0Mz 
/* 2-  */
  Mx0= 0.8811;Mz= 0.3322;F= 2.7780;theta0=-1.6469;U0I

  printf("Number of converged initial points %d\n",n0);fflush(stdout);

  MFMultifariosMethodClearClipF(H,e);
  MFMultifarioAddClipF(H,MzPositive,e);
  MFMultifarioAddClipF(H,Mx0Positive,e);
  MFMultifarioAddClipF(H,FPositive,e);
  MFMultifarioAddClipF(H,DNegative,e);
  MFMultifarioAddClipF(H,FLessThanSixteen,e);

  MFExtendAtlasMultiple(S,H,M,Omega,n0,u0,e);

  printf("Done Computing, now trying to write it\n");fflush(stdout);
  MFCloseAtlas(H,S,e);

  MFFreeAtlas(S,e);
  MFFreeImplicitMF(M,e);
  MFFreeContinuationMethod(H,e);
  MFFreeNRegion(Omega,e);
  for(i=0;i<n0;i++)MFFreeNVector(u0[i],e);

  MFFreeErrorHandler(e);

  printf("Normal termination of ComputeRod\n");fflush(stdout);
  return 0;
 }

MFNKMatrix ComputeTangent(MFImplicitMF M, int nx, int nu, int np, int nic, MFNVector ug, MFErrorHandler e)
 {
  MFNVector dug[4];
  MFNVector tug;
  MFNKMatrix Tan0;
  static double *u0=NULL;
  static double *p0=NULL;
  static double *TMat=NULL;
  static double *bc=NULL;
  static double *dbc=NULL;
  static double *ic=NULL;
  static double *dic=NULL;
  int i,j,n;
  double eps;

  u0=(double*)realloc((void*)u0,nu*sizeof(double));
  p0=(double*)realloc((void*)p0,np*sizeof(double));
  TMat=(double*)realloc((void*)TMat,(nu+np)*(nu+np)*sizeof(double));
  bc=(double*)realloc((void*)bc,nu*sizeof(double));
  dbc=(double*)realloc((void*)dbc,nu*sizeof(double));
  ic=(double*)realloc((void*)ic,nic*sizeof(double));
  dic=(double*)realloc((void*)dic,nic*sizeof(double));

  printf("In ComputeTangent\n");

  for(j=0;j<nu;j++)u0[j]=0.;
  for(j=0;j<np;j++)p0[j]=0.;
  eps=1.e-3;

  tug=MFTPBVPIntegrateForTangent(M,ug,u0,p0,e);
  MFTPBVPEvaluateBoundaryConditions(M,tug,tug,bc,e);
  MFTPBVPEvaluateIntegralConstraints(M,tug,tug,ic,e);
  MFFreeNVector(tug,e);

  for(i=0;i<nu+np;i++)
   {
    for(j=0;j<nu;j++)u0[j]=0.;
    for(j=0;j<np;j++)p0[j]=0.;

    if(i<nu)u0[i]=eps;
      else p0[i-nu]=eps;

    tug=MFTPBVPIntegrateForTangent(M,ug,u0,p0,e);
    MFTPBVPEvaluateBoundaryConditions(M,tug,tug,dbc,e);
    MFTPBVPEvaluateIntegralConstraints(M,tug,tug,dic,e);
    MFFreeNVector(tug,e);

    for(j=0;j<nu;j++)
      TMat[j+(nu+np)*i]=(dbc[j]-bc[j])/eps;
    for(j=0;j<nic;j++)
      TMat[nu+j+(nu+np)*i]=(dic[j]-ic[j])/eps;
    for(j=0;j<np-nic;j++)
      TMat[nu+nic+j+(nu+np)*i]=0.;

   }

  for(i=0;i<nu+np;i++)
   {
    if(i==0)printf("T=[%10.3le",TMat[i+(nu+np)*0]);
     else printf("  [%10.3le",TMat[i+(nu+np)*0]);
    for(j=1;j<nu+np;j++)
      printf(" %10.3le",TMat[i+(nu+np)*j]);
    printf("]\n");
   }

/* for each null vector of a */
  printf("Find tangent\n");
  n=1;
  for(j=0;j<n;j++)
   {
    n=GetNullVector(j,nu,np,j,TMat,u0,p0,e);
    if(j==0){printf("There were %d zero eigenvalues\n",n);fflush(stdout);}
    printf("%d[",j);
    for(i=0;i<nu;i++)
     printf(" %10.3le",u0[i]);
    for(i=0;i<np;i++)
     printf(" %10.3le",p0[i]);
    printf("]\n");
 
    dug[j]=MFTPBVPIntegrateForTangent(M,ug,u0,p0,e);
   }
  Tan0=MFCreateNKMatrix(MFIMF_K(M,e),dug,e);
/*MFGramSchmidt(Tan0);*/
/*printf("error in Phi %le\n",MFTPBVPTestTangent(M,ug,Tan0));fflush(stdout);*/

  return Tan0;
 }

int GetNullVector(int time, int nu, int np, int ev,double *A,double *x,double *p, MFErrorHandler e)
 {
  static char RoutineName[]={"GetNullVector"};
  static int i,j,l;
  static double mag;
  static int info;
  static char jobvl;
  static char jobvr;
  static int ldA;
  static double *wr=NULL;
  static double *wi=NULL;
  static int ldvl;
  static double *vl=NULL;
  static int ldvr;
  static double *vr=NULL;
  static double *work=NULL;
  static int lwork;

#ifdef HAVE_LIBLAPACK

  ldA=nu+np;
  wr=(double*)realloc((void*)wr,ldA*sizeof(double));
  wi=(double*)realloc((void*)wi,ldA*sizeof(double));
  ldvl=ldA;
  vl=(double*)realloc((void*)vl,ldvl*ldvl*sizeof(double));
  ldvr=ldA;
  vr=(double*)realloc((void*)vr,ldvr*ldvr*sizeof(double));
  lwork=4*ldA;
  work=(double*)realloc((void*)work,lwork*sizeof(double));
  info=0;
  jobvl='N';
  jobvr='V';
  if(time==0)
     Call_DGEEV(&jobvl,&jobvr,&ldA,A,&ldA,wr,wi,vl,&ldvl,vr,&ldvr,work,&lwork,&info);

  j=0;
  for(i=0;i<nu+np;i++)
   {
    mag=sqrt(wr[i]*wr[i]+wi[i]*wi[i]);
    if(mag<1.e-10)
     {
      if(j==ev)
       {
        for(l=0;l<nu;l++)x[l]=vr[l+(nu+np)*i];
        for(l=0;l<np;l++)p[l]=vr[nu+l+(nu+np)*i];
       }
      j++;
     }
   }

  return j;
#else
  sprintf(ComputeRodErrorHandlerMsg,"The ComputeRod example requires dgeev from LAPACK.");
  MFSetError(e,12,RoutineName,ComputeRodErrorHandlerMsg,__LINE__,__FILE__);
  return -1;
#endif
 }

MFNVector GetInitialPoint(MFImplicitMF M,double C,double A, double t, double theta0,double dA, MFErrorHandler e)
 {
  MFNVector u0,ug,ut[2];
  MFNKMatrix Tan;
  double *r0,dr;
  int i,n,nx,nu,np;
  int rc;

  printf("GetInitialPoint\n");fflush(stdout);

  n=MFIMF_N(M,e);
  nx=MFTPBVPGetNX(M,e);
  nu=MFTPBVPGetNU(M,e);
  np=MFTPBVPGetNP(M,e);
  r0=(double*)malloc(nx*sizeof(double));

  dr=1./(nx-2);
  for(i=0;i<nx;i++)r0[i]=(i-.5)*dr;

/* Assumes theta0= 0 or PI, C=0 */

  ug=MFIMFVectorFactory(M,e);
  ut[0]=MFIMFVectorFactory(M,e);
  ut[1]=MFIMFVectorFactory(M,e);
  for(i=0;i<nx;i++)
   {
    MFNVSetC(ug,nu*i+0,0.,e);
    MFNVSetC(ug,nu*i+1,0.,e);
    MFNVSetC(ug,nu*i+2,0.,e);
    MFNVSetC(ug,nu*i+3,0.,e);
    MFNVSetC(ug,nu*i+4,cos(theta0)*r0[i]*PI,e);
    MFNVSetC(ug,nu*i+5,cos(theta0),e);

    MFNVSetC(ug,nu*nx+np+i,r0[i],e);

    MFNVSetC(ut[0],nu*i+0,0.,e);
    MFNVSetC(ut[0],nu*i+1,0.,e);
    MFNVSetC(ut[0],nu*i+2,0.,e);
    MFNVSetC(ut[0],nu*i+3,0.,e);
    MFNVSetC(ut[0],nu*i+4,0.,e);
    MFNVSetC(ut[0],nu*i+5,0.,e);

    MFNVSetC(ut[0],nu*nx+np+i,r0[i],e);

    MFNVSetC(ut[1],nu*i+0,0.,e);
    MFNVSetC(ut[1],nu*i+1,0.,e);
    MFNVSetC(ut[1],nu*i+2,0.,e);
    MFNVSetC(ut[1],nu*i+3,0.,e);
    MFNVSetC(ut[1],nu*i+4,0.,e);
    MFNVSetC(ut[1],nu*i+5,0.,e);
    MFNVSetC(ut[1],nu*nx+np+i,r0[i],e);
   }
  MFNVSetC(ug,nu*nx+0,A,e);
  MFNVSetC(ug,nu*nx+1,C,e);
  MFNVSetC(ug,nu*nx+2,t,e);
  MFNVSetC(ug,nu*nx+3,theta0,e);

  MFNVSetC(ut[0],nu*nx+0,0.,e);
  MFNVSetC(ut[0],nu*nx+1,dA,e);
  MFNVSetC(ut[0],nu*nx+2,1.,e);
  MFNVSetC(ut[0],nu*nx+3,0.,e);

  MFNVSetC(ut[1],nu*nx+0,0.,e);
  MFNVSetC(ut[1],nu*nx+1,1.,e);
  MFNVSetC(ut[1],nu*nx+2,0.,e);
  MFNVSetC(ut[1],nu*nx+3,0.,e);
/* Try */
  MFNVSetC(ug,nu*nx+0,A+.1,e);
  MFNVSetC(ut[0],nu*nx+0,0.,e);
  MFNVSetC(ut[0],nu*nx+1,0.,e);
  MFNVSetC(ut[0],nu*nx+2,1.,e);
  MFNVSetC(ut[0],nu*nx+3,0.,e);
  MFNVSetC(ut[1],nu*nx+0,1.,e);
  MFNVSetC(ut[1],nu*nx+1,0.,e);
  MFNVSetC(ut[1],nu*nx+2,0.,e);
  MFNVSetC(ut[1],nu*nx+3,0.,e);
/* end Try */
  Tan=MFCreateNKMatrix(2,ut,e);
  Tan=MFIMFTangentSpace(M,ug,e);
  printf("ug=(%le,%le,%le,%le)\n",MFNV_C(ug,nu*nx+0,e),MFNV_C(ug,nu*nx+1,e),MFNV_C(ug,nu*nx+2,e),MFNV_C(ug,nu*nx+3,e));fflush(stdout);
  MFFreeNVector(ut[0],e);
  MFFreeNVector(ut[1],e);
  MFFreeNKMatrix(Tan,e);
  free(r0);
  return ug;

  u0=MFIMFVectorFactory(M,e);
  rc=MFIMFProject(M,ug,Tan,u0,e);
  printf("u0=(%le,%le,%le,%le)\n",MFNV_C(u0,nu*nx+0,e),MFNV_C(u0,nu*nx+1,e),MFNV_C(u0,nu*nx+2,e),MFNV_C(u0,nu*nx+3,e));fflush(stdout);

  MFFreeNVector(ug,e);

  MFFreeNVector(ut[0],e);
  MFFreeNVector(ut[1],e);
  MFFreeNKMatrix(Tan,e);
  free(r0);
  if(rc){printf("Project Successful\n");fflush(stdout);}
  if(rc)return u0;
  printf("Project Failed\n");fflush(stdout);
  MFFreeNVector(u0,e);
  return NULL;
 }

MFNVector GetInitialPointMx0Mz(MFImplicitMF M,double Mx0,double Mz, double f, double theta0, MFErrorHandler e)
 {
  double r,dr;
  double xy0[6],p0[4],p1[4];
  double *r0;
  double error,e0,e1,corr;
  double dMz,dMx0,dy1,dy3,a00,a01,a10,a11,de0,de1;
  int itimes;
  MFNVector ug,tug;
  MFNVector u0;
  MFNKMatrix Tan;
  int i,n,nx,nu,np;
  double xnorm;
  int verbose=0;

  if(verbose){printf("GetInitialPointMx0Mz\n");fflush(stdout);}

  n=MFIMF_N(M,e);
  nx=MFTPBVPGetNX(M,e);
  nu=MFTPBVPGetNU(M,e);
  np=MFTPBVPGetNP(M,e);

  r0=(double*)malloc(nx*sizeof(double));
  dr=1./(nx-2);
  for(i=0;i<nx;i++)r0[i]=(i-.5)*dr;

  error=1.;
  itimes=0;

  xy0[0]=0.;
  xy0[1]=sin(theta0);
  xy0[2]=0.;
  xy0[3]=0.;
  xy0[4]=0.;
  xy0[5]=cos(theta0);
  p0[0]=Mz;
  p0[1]=Mx0;
  p0[2]=f;
  p0[3]=theta0;
  if(verbose){printf("Iter (%lf,%lf,%lf,%lf)\n",p0[1],p0[0],p0[2],p0[3]);fflush(stdout);}

  ug=MFTPBVPIntegrateForInitialSolution(M,xy0,p0,r0,e);

  e0=(MFNV_C(ug,nu*(nx-1)+3,e)+MFNV_C(ug,nu*(nx-2)+3,e))/2;
  e1=(MFNV_C(ug,nu*(nx-1)  ,e)+MFNV_C(ug,nu*(nx-2)  ,e))/2*
     (MFNV_C(ug,nu*(nx-1)+5,e)+MFNV_C(ug,nu*(nx-2)+5,e))/2
    -(MFNV_C(ug,nu*(nx-1)+4,e)+MFNV_C(ug,nu*(nx-2)+4,e))/2*
     (MFNV_C(ug,nu*(nx-1)+1,e)+MFNV_C(ug,nu*(nx-2)+1,e))/2;

  error=sqrt(e0*e0+e1*e1);

/*printf(" Iteration to find u0:\n");*/
  while(error>1.e-10 && itimes<200)
   {
    if(verbose){printf("%d error=%le Mx0(C)=%lf, Mz(A)=%lf, f(t)=%lf, th0=%lf\n",itimes,error,p0[0],p0[1],p0[2],p0[3]);}

    dMz=1.e-5;
    dMx0=1.e-5;

    xy0[0]=0.;
    xy0[1]=sin(theta0);
    xy0[2]=0.;
    xy0[3]=0.;
    xy0[4]=0.;
    xy0[5]=cos(theta0);
    p0[0]=Mz+dMz;
    p0[1]=Mx0;
    p0[2]=f;
    p0[3]=theta0;

    tug=MFTPBVPIntegrateForInitialSolution(M,xy0,p0,r0,e);

    de0=(MFNV_C(tug,nu*(nx-1)+3,e)+MFNV_C(tug,nu*(nx-2)+3,e))/2;
    de1=(MFNV_C(tug,nu*(nx-1)  ,e)+MFNV_C(tug,nu*(nx-2)  ,e))/2*
         (MFNV_C(tug,nu*(nx-1)+5,e)+MFNV_C(tug,nu*(nx-2)+5,e))/2
   -(MFNV_C(tug,nu*(nx-1)+4,e)+MFNV_C(tug,nu*(nx-2)+4,e))/2*
         (MFNV_C(tug,nu*(nx-1)+1,e)+MFNV_C(tug,nu*(nx-2)+1,e))/2;
    MFFreeNVector(tug,e);
    a00=(de0-e0)/dMz;
    a10=(de1-e1)/dMz;

    xy0[0]=0.;
    xy0[1]=sin(theta0);
    xy0[2]=0.;
    xy0[3]=0.;
    xy0[4]=0.;
    xy0[5]=cos(theta0);
    p0[0]=Mz;
    p0[1]=Mx0+dMx0;
    p0[2]=f;
    p0[3]=theta0;

    tug=MFTPBVPIntegrateForInitialSolution(M,xy0,p0,r0,e);
    de0=(MFNV_C(tug,nu*(nx-1)+3,e)+MFNV_C(tug,nu*(nx-2)+3,e))/2;
    de1=(MFNV_C(tug,nu*(nx-1)  ,e)+MFNV_C(tug,nu*(nx-2)  ,e))/2*
        (MFNV_C(tug,nu*(nx-1)+5,e)+MFNV_C(tug,nu*(nx-2)+5,e))/2
       -(MFNV_C(tug,nu*(nx-1)+4,e)+MFNV_C(tug,nu*(nx-2)+4,e))/2*
        (MFNV_C(tug,nu*(nx-1)+1,e)+MFNV_C(tug,nu*(nx-2)+1,e))/2;
    MFFreeNVector(tug,e);
    a01=(de0-e0)/dMx0;
    a11=(de1-e1)/dMx0;

    dMz =-( a11*e0-a01*e1)/(a00*a11-a01*a10);
    dMx0=-(-a10*e0+a00*e1)/(a00*a11-a01*a10);
    corr=sqrt(dMz*dMz+dMx0*dMx0);
    xnorm=sqrt(Mz*Mz+Mx0*Mx0);
    if(verbose){printf("      Delta=(%lf,%lf) ||=%lf X=(%lf,%lf),||=%lf\n",dMx0,dMz,corr,Mx0,Mz,sqrt(Mz*Mz+Mx0*Mx0));fflush(stdout);}
/*  if(corr>.1){dMz=.1*dMz/corr;dMx0=.1*dMx0/corr;}*/
    if(0&&corr>1e-2*xnorm){dMz=1.e-2*xnorm*dMz/corr;dMx0=1.e-2*xnorm*dMx0/corr;}
    Mz=Mz+dMz;
    Mx0=Mx0+dMx0;
    itimes++;

    xy0[0]=0.;
    xy0[1]=sin(theta0);
    xy0[2]=0.;
    xy0[3]=0.;
    xy0[4]=0.;
    xy0[5]=cos(theta0);
    p0[0]=Mz;
    p0[1]=Mx0;
    p0[2]=f;
    p0[3]=theta0;

    MFFreeNVector(ug,e);
    ug=MFTPBVPIntegrateForInitialSolution(M,xy0,p0,r0,e);

    e0=(MFNV_C(ug,nu*(nx-1)+3,e)+MFNV_C(ug,nu*(nx-2)+3,e))/2;
    e1=(MFNV_C(ug,nu*(nx-1)  ,e)+MFNV_C(ug,nu*(nx-2)  ,e))/2*
       (MFNV_C(ug,nu*(nx-1)+5,e)+MFNV_C(ug,nu*(nx-2)+5,e))/2
      -(MFNV_C(ug,nu*(nx-1)+4,e)+MFNV_C(ug,nu*(nx-2)+4,e))/2*
       (MFNV_C(ug,nu*(nx-1)+1,e)+MFNV_C(ug,nu*(nx-2)+1,e))/2;

    error=sqrt(e0*e0+e1*e1);
   }
  if(error<1.e-8)
   {
    if(verbose){printf("  got(%lf,%lf,%lf,%lf)\n",p0[1],p0[0],p0[2],p0[3]);fflush(stdout);}
    Tan=MFIMFTangentSpace(M,ug,e);
    u0=MFIMFVectorFactory(M,e);
    MFIMFProject(M,ug,Tan,u0,e);
    MFFreeNKMatrix(Tan,e);
   }else{
    u0=NULL;
    if(verbose){printf("  failed\n");fflush(stdout);}
   }
  MFFreeNVector(ug,e);
  free(r0);
  return u0;
 }

void flipf(int nx,int nu,MFNVector u, MFErrorHandler e)
 {
  int i;
  double th0;

  for(i=0;i<nx;i++)
   {
    MFNVSetC(u,nu*i+0,-MFNV_C(u,nu*i+0,e),e);
    MFNVSetC(u,nu*i+1,-MFNV_C(u,nu*i+1,e),e);
    MFNVSetC(u,nu*i+2,-MFNV_C(u,nu*i+2,e),e);
    MFNVSetC(u,nu*i+3,-MFNV_C(u,nu*i+3,e),e);
    MFNVSetC(u,nu*i+4,-MFNV_C(u,nu*i+4,e),e);
    MFNVSetC(u,nu*i+5,-MFNV_C(u,nu*i+5,e),e);
   }
  MFNVSetC(u,nx*nu+2,-MFNV_C(u,nx*nu+2,e),e);
  th0=MFNV_C(u,nx*nu+3,e)+3.1415926;
  if(th0>2*3.1415926)th0=th0-2*3.1415926;
  MFNVSetC(u,nx*nu+3,th0,e);

  return;
 }

void flipMx0AndMz(int nx,int nu,MFNVector u, MFErrorHandler e)
 {
  int i;
  double th0;

  for(i=0;i<nx;i++)
   {
    MFNVSetC(u,nu*i+2,-MFNV_C(u,nu*i+2,e),e);
    MFNVSetC(u,nu*i+3,-MFNV_C(u,nu*i+3,e),e);
   }
  MFNVSetC(u,nx*nu+0,-MFNV_C(u,nx*nu+0,e),e);
  MFNVSetC(u,nx*nu+1,-MFNV_C(u,nx*nu+1,e),e);

  return;
 }

void flipMz(int nx,int nu,MFNVector u, MFErrorHandler e)
 {
  int i;
  double th0;

  for(i=0;i<nx;i++)
   {
    MFNVSetC(u,nu*i+0,-MFNV_C(u,nu*i+0,e),e);
    MFNVSetC(u,nu*i+1,-MFNV_C(u,nu*i+1,e),e);
   }
  MFNVSetC(u,nx*nu+0,-MFNV_C(u,nx*nu+0,e),e);
  th0=-MFNV_C(u,nx*nu+3,e);
  while(th0>2*3.1415926)th0=th0-2*3.1415926;
  while(th0<0.)th0=th0+2*3.1415926;
  MFNVSetC(u,nx*nu+3,th0,e);

  return;
 }

void flipMx0(int nx,int nu,MFNVector u, MFErrorHandler e)
 {
  flipMx0AndMz(nx,nu,u,e);
  flipMz(nx,nu,u,e);
  return;
 }

void moveToPrinciple(int nx,int nu,MFNVector u, MFErrorHandler e)
 {
  printf("moveToPrinciple, was (%lf,%lf,%lf,%lf)\n",MFNV_C(u,nx*nu+1,e),MFNV_C(u,nx*nu+0,e),MFNV_C(u,nx*nu+2,e),MFNV_C(u,nx*nu+3,e));
  while(MFNV_C(u,nx*nu+3,e)<-PI){MFNVSetC(u,nx*nu+3,MFNV_C(u,nx*nu+3,e)+2*PI,e);}
  while(MFNV_C(u,nx*nu+3,e)> PI){MFNVSetC(u,nx*nu+3,MFNV_C(u,nx*nu+3,e)-2*PI,e);}

  if(MFNV_C(u,nx*nu+3,e)>PI/2)flipf(nx,nu,u,e);
  while(MFNV_C(u,nx*nu+3,e)<-PI){MFNVSetC(u,nx*nu+3,MFNV_C(u,nx*nu+3,e)+2*PI,e);}
  while(MFNV_C(u,nx*nu+3,e)> PI){MFNVSetC(u,nx*nu+3,MFNV_C(u,nx*nu+3,e)-2*PI,e);}
  if(MFNV_C(u,nx*nu+3,e)<-PI/2)flipf(nx,nu,u,e);
  while(MFNV_C(u,nx*nu+3,e)<-PI){MFNVSetC(u,nx*nu+3,MFNV_C(u,nx*nu+3,e)+2*PI,e);}
  while(MFNV_C(u,nx*nu+3,e)> PI){MFNVSetC(u,nx*nu+3,MFNV_C(u,nx*nu+3,e)-2*PI,e);}
  while(MFNV_C(u,nx*nu+3,e)< 0.)flipMz(nx,nu,u,e);
  printf("       result is     (%lf,%lf,%lf,%lf)\n",MFNV_C(u,nx*nu+1,e),MFNV_C(u,nx*nu+0,e),MFNV_C(u,nx*nu+2,e),MFNV_C(u,nx*nu+3,e));

  return;
 }

MFNVector GetInitialPointI(MFImplicitMF M,double Mx0,double Mz, double f, double theta0, MFErrorHandler e)
 {
  double r,dr;
  double xy0[6],p0[4],p1[4];
  double *r0;
  MFNVector u;
  MFNVector ug;
  MFNVector du[2];
  MFNKMatrix Tan;
  int i,n,nx,nu,np;

  printf("GetInitialPointI\n");fflush(stdout);

  n=MFIMF_N(M,e);
  nx=MFTPBVPGetNX(M,e);
  nu=MFTPBVPGetNU(M,e);
  np=MFTPBVPGetNP(M,e);

  r0=(double*)malloc(nx*sizeof(double));
  dr=1./(nx-2);
  for(i=0;i<nx;i++)r0[i]=(i-.5)*dr;

  xy0[0]=0.;
  xy0[1]=sin(theta0);
  xy0[2]=0.;
  xy0[3]=0.;
  xy0[4]=0.;
  xy0[5]=cos(theta0);
  p0[0]=Mz;
  p0[1]=Mx0;
  p0[2]=f;
  p0[3]=theta0;
  printf("     (%lf,%lf,%lf,%lf)\n",p0[1],p0[0],p0[2],p0[3]);fflush(stdout);

  ug=MFTPBVPIntegrateForInitialSolution(M,xy0,p0,r0,e);
  Tan=MFIMFTangentSpace(M,ug,e);
  u=MFIMFVectorFactory(M,e);
/*
  du[0]=MFIMFVectorFactory(M,e);
  du[1]=MFIMFVectorFactory(M,e);
  for(i=0;i<nx*nu;i++)
   {
    MFNVSetC(du[0],i,0.,e);
    MFNVSetC(du[1],i,0.,e);
   }
  MFNVSetC(du[0],nx*nu+0,0.,e);MFNVSetC(du[1],nx*nu+0,0.,e);
  MFNVSetC(du[0],nx*nu+1,0.,e);MFNVSetC(du[1],nx*nu+1,0.,e);
  MFNVSetC(du[0],nx*nu+2,1.,e);MFNVSetC(du[1],nx*nu+2,0.,e);
  MFNVSetC(du[0],nx*nu+3,0.,e);MFNVSetC(du[1],nx*nu+3,1.,e);
  Tan=MFCreateNKMatrix(2,du,e);
*/

  if(!MFIMFProject(M,ug,Tan,u,e))
   {
    printf("     (%lf,%lf,%lf,%lf) failed\n",p0[1],p0[0],p0[2],p0[3]);fflush(stdout);
    MFFreeNVector(u,e);
    u=NULL;
   }else{
    printf("  got(%lf,%lf,%lf,%lf)\n",MFNV_C(u,nx*nu+1,e),MFNV_C(u,nx*nu+0,e),MFNV_C(u,nx*nu+2,e),MFNV_C(u,nx*nu+3,e),e);fflush(stdout);
   }

  free(r0);
  MFFreeNVector(ug,e);
/*
  MFFreeNVector(du[0],e);
  MFFreeNVector(du[1],e);*/

  MFFreeNKMatrix(Tan,e);

  return u;
 }

void moveToPositive(int nx,int nu,MFNVector u, MFErrorHandler e)
 {
  printf("moveToPositive, was (%lf,%lf,%lf,%lf)\n",MFNV_C(u,nx*nu+1,e),MFNV_C(u,nx*nu+0,e),MFNV_C(u,nx*nu+2,e),MFNV_C(u,nx*nu+3,e));

  if(MFNV_C(u,nx*nu+0,e)<0)flipMz(nx,nu,u,e);
  if(MFNV_C(u,nx*nu+1,e)<0)flipMx0(nx,nu,u,e);
  if(MFNV_C(u,nx*nu+2,e)<0)flipf(nx,nu,u,e);

  while(MFNV_C(u,nx*nu+3,e)<0){MFNVSetC(u,nx*nu+3,MFNV_C(u,nx*nu+3,e)+2*PI,e);}
  while(MFNV_C(u,nx*nu+3,e)>2*PI){MFNVSetC(u,nx*nu+3,MFNV_C(u,nx*nu+3,e)-2*PI,e);}

  printf("       result is     (%lf,%lf,%lf,%lf)\n",MFNV_C(u,nx*nu+1,e),MFNV_C(u,nx*nu+0,e),MFNV_C(u,nx*nu+2,e),MFNV_C(u,nx*nu+3,e));

  return;
 }

MFNVector GetInitialPointMx0f(MFImplicitMF M,double Mz,double Mx0, double f, double theta0, MFErrorHandler e)
 {
  double r,dr;
  double xy0[6],p0[4],p1[4];
  double *r0;
  double error,e0,e1,corr;
  double df,dMx0,dy1,dy3,a00,a01,a10,a11,de0,de1;
  int itimes;
  MFNVector ug,tug;
  MFNVector u0;
  MFNKMatrix Tan;
  int i,n,nx,nu,np;
  double xnorm;
  printf("GetInitialPointMx0f\n");fflush(stdout);

  n=MFIMF_N(M,e);
  nx=MFTPBVPGetNX(M,e);
  nu=MFTPBVPGetNU(M,e);
  np=MFTPBVPGetNP(M,e);

  r0=(double*)malloc(nx*sizeof(double));
  dr=1./(nx-2);
  for(i=0;i<nx;i++)r0[i]=(i-.5)*dr;

  error=1.;
  itimes=0;

  xy0[0]=0.;
  xy0[1]=sin(theta0);
  xy0[2]=0.;
  xy0[3]=0.;
  xy0[4]=0.;
  xy0[5]=cos(theta0);
  p0[0]=Mz;
  p0[1]=Mx0;
  p0[2]=f;
  p0[3]=theta0;
  printf("Iter (%lf,%lf,%lf,%lf)\n",p0[1],p0[0],p0[2],p0[3]);fflush(stdout);

  ug=MFTPBVPIntegrateForInitialSolution(M,xy0,p0,r0,e);

  e0=(MFNV_C(ug,nu*(nx-1)+3,e)+MFNV_C(ug,nu*(nx-2)+3,e))/2;
  e1=(MFNV_C(ug,nu*(nx-1)  ,e)+MFNV_C(ug,nu*(nx-2)  ,e))/2*
     (MFNV_C(ug,nu*(nx-1)+5,e)+MFNV_C(ug,nu*(nx-2)+5,e))/2
    -(MFNV_C(ug,nu*(nx-1)+4,e)+MFNV_C(ug,nu*(nx-2)+4,e))/2*
     (MFNV_C(ug,nu*(nx-1)+1,e)+MFNV_C(ug,nu*(nx-2)+1,e))/2;

  error=sqrt(e0*e0+e1*e1);

/*printf(" Iteration to find u0:\n");*/
  while(error>1.e-10 && itimes<200)
   {
    printf("%d error=%le Mx0(C)=%lf, Mz(A)=%lf, f(t)=%lf, th0=%lf\n",itimes,error,p0[0],p0[1],p0[2],p0[3]);

    df=1.e-5;
    dMx0=1.e-5;

    xy0[0]=0.;
    xy0[1]=sin(theta0);
    xy0[2]=0.;
    xy0[3]=0.;
    xy0[4]=0.;
    xy0[5]=cos(theta0);
    p0[0]=Mz;
    p0[1]=Mx0;
    p0[2]=f+df;
    p0[3]=theta0;

    tug=MFTPBVPIntegrateForInitialSolution(M,xy0,p0,r0,e);

    de0=(MFNV_C(tug,nu*(nx-1)+3,e)+MFNV_C(tug,nu*(nx-2)+3,e))/2;
    de1=(MFNV_C(tug,nu*(nx-1)  ,e)+MFNV_C(tug,nu*(nx-2)  ,e))/2*
        (MFNV_C(tug,nu*(nx-1)+5,e)+MFNV_C(tug,nu*(nx-2)+5,e))/2
       -(MFNV_C(tug,nu*(nx-1)+4,e)+MFNV_C(tug,nu*(nx-2)+4,e))/2*
        (MFNV_C(tug,nu*(nx-1)+1,e)+MFNV_C(tug,nu*(nx-2)+1,e))/2;
    MFFreeNVector(tug,e);
    a00=(de0-e0)/df;
    a10=(de1-e1)/df;

    xy0[0]=0.;
    xy0[1]=sin(theta0);
    xy0[2]=0.;
    xy0[3]=0.;
    xy0[4]=0.;
    xy0[5]=cos(theta0);
    p0[0]=Mz;
    p0[1]=Mx0+dMx0;
    p0[2]=f;
    p0[3]=theta0;

    tug=MFTPBVPIntegrateForInitialSolution(M,xy0,p0,r0,e);
    de0=(MFNV_C(tug,nu*(nx-1)+3,e)+MFNV_C(tug,nu*(nx-2)+3,e))/2;
    de1=(MFNV_C(tug,nu*(nx-1)  ,e)+MFNV_C(tug,nu*(nx-2)  ,e))/2*
        (MFNV_C(tug,nu*(nx-1)+5,e)+MFNV_C(tug,nu*(nx-2)+5,e))/2
       -(MFNV_C(tug,nu*(nx-1)+4,e)+MFNV_C(tug,nu*(nx-2)+4,e))/2*
        (MFNV_C(tug,nu*(nx-1)+1,e)+MFNV_C(tug,nu*(nx-2)+1,e))/2;
    MFFreeNVector(tug,e);
    a01=(de0-e0)/dMx0;
    a11=(de1-e1)/dMx0;

    df  =-( a11*e0-a01*e1)/(a00*a11-a01*a10);
    dMx0=-(-a10*e0+a00*e1)/(a00*a11-a01*a10);
    corr=sqrt(df*df+dMx0*dMx0);
    xnorm=sqrt(f*f+Mx0*Mx0);
    printf("      Delta=(%lf,%lf) ||=%lf X=(%lf,%lf),||=%lf\n",dMx0,df,corr,Mx0,f,sqrt(Mx0*Mx0+f*f));fflush(stdout);
/*  if(corr>.1){df =.1*df /corr;dMx0=.1*dMx0/corr;}*/
/*  if(corr>1e-2*xnorm){df =1.e-2*xnorm*df /corr;dMx0=1.e-2*xnorm*dMx0/corr;}*/
    f=f+df;
    Mx0=Mx0+dMx0;
    itimes++;

    xy0[0]=0.;
    xy0[1]=sin(theta0);
    xy0[2]=0.;
    xy0[3]=0.;
    xy0[4]=0.;
    xy0[5]=cos(theta0);
    p0[0]=Mz;
    p0[1]=Mx0;
    p0[2]=f;
    p0[3]=theta0;

    MFFreeNVector(ug,e);
    ug=MFTPBVPIntegrateForInitialSolution(M,xy0,p0,r0,e);

    e0=(MFNV_C(ug,nu*(nx-1)+3,e)+MFNV_C(ug,nu*(nx-2)+3,e))/2;
    e1=(MFNV_C(ug,nu*(nx-1)  ,e)+MFNV_C(ug,nu*(nx-2)  ,e))/2*
       (MFNV_C(ug,nu*(nx-1)+5,e)+MFNV_C(ug,nu*(nx-2)+5,e))/2
      -(MFNV_C(ug,nu*(nx-1)+4,e)+MFNV_C(ug,nu*(nx-2)+4,e))/2*
       (MFNV_C(ug,nu*(nx-1)+1,e)+MFNV_C(ug,nu*(nx-2)+1,e))/2;

    error=sqrt(e0*e0+e1*e1);
   }
  if(error<1.e-8)
   { 
    printf("  got(%lf,%lf,%lf,%lf)\n",p0[1],p0[0],p0[2],p0[3]);fflush(stdout);
    Tan=MFIMFTangentSpace(M,ug,e);
    u0=MFIMFVectorFactory(M,e);
    MFIMFProject(M,ug,Tan,u0,e);
    MFFreeNKMatrix(Tan,e);
   }else{
    printf("  failed\n");fflush(stdout);
    u0=NULL;
   }
  MFFreeNVector(ug,e);
  free(r0);

  return u0;
 }

int MFRodStopN(MFImplicitMF M,MFNVector U0,MFNKMatrix Phi0,MFNVector U1,MFNKMatrix Phi1,void *data, MFErrorHandler e)
 {
  double u0,um,up;
  double theta0,Mx0,Mz,f;
  double delta,Omega,T,gnu,m;
  double x1,y1,x2,y2,x3,y3,del0,del1;
  int s00,s01,s02;
  int s10,s11,s12;
  int N0,N1;
  int nx=150;
  int np=4;
  int nu=6;

  Mz    =MFNV_C(U0,nx*nu+0,e);
  Mx0   =MFNV_C(U0,nx*nu+1,e);
  f     =MFNV_C(U0,nx*nu+2,e);
  theta0=MFNV_C(U0,nx*nu+3,e);

  s00=1;
  if(theta0<PI)s00=-1;
  s02=1;
  if(theta0<-PI)s02=-1;
  s01=1;
  if(Mz<0)s01=-1;

  u0=cos(theta0);
  delta=(Mz*Mz+Mx0*Mx0)*(Mz*Mz+Mx0*Mx0)+8*f*(2*f+u0*(Mx0*Mx0-Mz*Mz))-16*Mz*Mx0*f*sin(theta0);
  up=(Mz*Mz+Mx0*Mx0+sqrt(delta))/4/f;
  um=(Mz*Mz+Mx0*Mx0-sqrt(delta))/4/f;

  gnu=up-u0;
  if(um<u0)gnu=up-um;
  m=fabs(u0-um)/gnu;
  if(m>0)T=2*sqrt(2.)*K(sqrt(m))/sqrt(fabs(f*gnu));
   else  T=2*sqrt(2.)*Ki(sqrt(-m))/sqrt(fabs(f*gnu));
  Omega=2*PI/T;

  N0=round(Omega);
  if(um<u0)N0=N0+1;

  x1=(MFNV_C(U0,nu*(nx-2)+0,e)+MFNV_C(U0,nu*(nx-1)+0,e))/2;
  y1=(MFNV_C(U0,nu*(nx-2)+1,e)+MFNV_C(U0,nu*(nx-1)+1,e))/2;
  x2=(MFNV_C(U0,nu*(nx-2)+2,e)+MFNV_C(U0,nu*(nx-1)+2,e))/2;
  y2=(MFNV_C(U0,nu*(nx-2)+3,e)+MFNV_C(U0,nu*(nx-1)+3,e))/2;
  x3=(MFNV_C(U0,nu*(nx-2)+4,e)+MFNV_C(U0,nu*(nx-1)+4,e))/2;
  y3=(MFNV_C(U0,nu*(nx-2)+5,e)+MFNV_C(U0,nu*(nx-1)+5,e))/2;
  del0=1.-(x1*y1+x2*y2+x3*y3)/3.1415926;
  if(del0<1.)N0=-N0;

  Mz    =MFNV_C(U1,nx*nu+0,e);
  Mx0   =MFNV_C(U1,nx*nu+1,e);
  f     =MFNV_C(U1,nx*nu+2,e);
  theta0=MFNV_C(U1,nx*nu+3,e);

  s10=1;
  if(theta0<PI)s10=-1;
  s12=1;
  if(theta0<-PI)s12=-1;
  s11=1;
  if(Mz<0)s11=-1;

  u0=cos(theta0);
  delta=(Mz*Mz+Mx0*Mx0)*(Mz*Mz+Mx0*Mx0)+8*f*(2*f+u0*(Mx0*Mx0-Mz*Mz))-16*Mz*Mx0*f*sin(theta0);
  up=(Mz*Mz+Mx0*Mx0+sqrt(delta))/4/f;
  um=(Mz*Mz+Mx0*Mx0-sqrt(delta))/4/f;

  gnu=up-u0;
  if(um<u0)gnu=up-um;
  m=fabs(u0-um)/gnu;
  if(m>0)T=2*sqrt(2.)*K(sqrt(m))/sqrt(fabs(f*gnu));
   else  T=2*sqrt(2.)*Ki(sqrt(-m))/sqrt(fabs(f*gnu));
  Omega=2*PI/T;

  N1=round(Omega);
  if(um<u0)N1=N1+1;

  x1=(MFNV_C(U1,nu*(nx-2)+0,e)+MFNV_C(U1,nu*(nx-1)+0,e))/2;
  y1=(MFNV_C(U1,nu*(nx-2)+1,e)+MFNV_C(U1,nu*(nx-1)+1,e))/2;
  x2=(MFNV_C(U1,nu*(nx-2)+2,e)+MFNV_C(U1,nu*(nx-1)+2,e))/2;
  y2=(MFNV_C(U1,nu*(nx-2)+3,e)+MFNV_C(U1,nu*(nx-1)+3,e))/2;
  x3=(MFNV_C(U1,nu*(nx-2)+4,e)+MFNV_C(U1,nu*(nx-1)+4,e))/2;
  y3=(MFNV_C(U1,nu*(nx-2)+5,e)+MFNV_C(U1,nu*(nx-1)+5,e))/2;
  del1=1.-(x1*y1+x2*y2+x3*y3)/3.1415926;
  if(del1<1.)N1=-N1;

/*printf("StopN. N0=%d, del0=%lf, s00=%d, s01=%d\n",N0,del0,s00,s01);
  printf("       N1=%d, del1=%lf, s10=%d, s11=%d\n",N1,del1,s10,s11);
  printf("  result=%d\n",N0!=N1 || del0>1 || del1>1 || s00*s10<0&&s01*s11<0);*/

  return N0!=N1 || del0>1 || del1>1 || s00*s10<0&&s01*s11<0
                                    || s00*s10<0&&s02*s12<0;
 }


int MFRodStopD(MFImplicitMF M, MFNVector U0,MFNKMatrix Phi0,MFNVector U1,MFNKMatrix Phi1,void *data, MFErrorHandler e)
 {
  double x1,y1,x2,y2,x3,y3;
  double d0,d1;
  int nx=150;
  int np=4;
  int nu=6;

  x1=(MFNV_C(U0,nu*(nx-2)+0,e)+MFNV_C(U0,nu*(nx-1)+0,e))/2;
  y1=(MFNV_C(U0,nu*(nx-2)+1,e)+MFNV_C(U0,nu*(nx-1)+1,e))/2;
  x2=(MFNV_C(U0,nu*(nx-2)+2,e)+MFNV_C(U0,nu*(nx-1)+2,e))/2;
  y2=(MFNV_C(U0,nu*(nx-2)+3,e)+MFNV_C(U0,nu*(nx-1)+3,e))/2;
  x3=(MFNV_C(U0,nu*(nx-2)+4,e)+MFNV_C(U0,nu*(nx-1)+4,e))/2;
  y3=(MFNV_C(U0,nu*(nx-2)+5,e)+MFNV_C(U0,nu*(nx-1)+5,e))/2;
  d0=x1*y1+x2*y2+x3*y3;

  x1=(MFNV_C(U1,nu*(nx-2)+0,e)+MFNV_C(U1,nu*(nx-1)+0,e))/2;
  y1=(MFNV_C(U1,nu*(nx-2)+1,e)+MFNV_C(U1,nu*(nx-1)+1,e))/2;
  x2=(MFNV_C(U1,nu*(nx-2)+2,e)+MFNV_C(U1,nu*(nx-1)+2,e))/2;
  y2=(MFNV_C(U1,nu*(nx-2)+3,e)+MFNV_C(U1,nu*(nx-1)+3,e))/2;
  x3=(MFNV_C(U1,nu*(nx-2)+4,e)+MFNV_C(U1,nu*(nx-1)+4,e))/2;
  y3=(MFNV_C(U1,nu*(nx-2)+5,e)+MFNV_C(U1,nu*(nx-1)+5,e))/2;
  d1=x1*y1+x2*y2+x3*y3;

  return fabs(d0)>.05 || fabs(d1)>.05;
 }

int MFRodProjectToDraw(MFNVector u, double *x, void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFRodProjectToDraw"};
  struct MFTPBVPData *data;
  int i,j;
  int n,nu,nx,np;
  double min,max;
  double a,c,t;
  double x1,y1,x2,y2,x3,y3,kp0,kp;
  double mx0,mz,f;
  double theta0,m3;
  double u0,up,um,delta;
  double T,Omega,R,Wr;
  double dist,dth0;
  double d1x,d1y,d1z;
  double gnu;
  double m;
  int verbose=0;
  double Delta,del;
  int N;
  char name[1024];
  FILE *fid;

  data=(struct MFTPBVPData*)d;

  if(x==NULL)return 9;

  n=MFNV_NC(u,e);

  nu=data->nu;
  nx=data->nx;
  np=data->np;

  if(verbose){printf("%s, n=%d, nu=%d, nx=%d, np=%d, comp=%d\n",RoutineName,n,nu,nx,np,nx*nu+np+nx+1);fflush(stdout);}
  min=0;max=0;
  a=MFNV_C(u,nx*nu+0,e);
  c=MFNV_C(u,nx*nu+1,e);
  t=MFNV_C(u,nx*nu+2,e);

  y1=(MFNV_C(u,1,e)+MFNV_C(u,nu+1,e))/2;
  y3=(MFNV_C(u,5,e)+MFNV_C(u,nu+5,e))/2;

  mz=a;
  mx0=c;
  f=t;
  theta0=MFNV_C(u,nx*nu+3,e);
  m3=mx0*y1+mz*y3;

  u0=cos(theta0);
  delta=(mz*mz+mx0*mx0)*(mz*mz+mx0*mx0)+8*f*(2*f+u0*(mx0*mx0-mz*mz))-16*mz*mx0*f*sin(theta0);
  if(fabs(f)>1.e-5)
   {
    up=(mz*mz+mx0*mx0+sqrt(delta))/4/f;
    um=(mz*mz+mx0*mx0-sqrt(delta))/4/f;
   }else{
    up=1000.;
    um=0.;
   }

  gnu=up-u0;
  if(um<u0)gnu=up-um;
  m=fabs(u0-um)/gnu;
  if(verbose)
   {
    if(m>0){printf("   before K(%lf)\n",sqrt(m));fflush(stdout);}
     else {printf("   before Ki(%lf)\n",sqrt(-m));fflush(stdout);}
   }
  if(fabs(m)>1.e-5)
   {
    if(m>0)T=2*sqrt(2.)*K(sqrt(m))/sqrt(fabs(f*gnu));
     else  T=2*sqrt(2.)*Ki(sqrt(-m))/sqrt(fabs(f*gnu));
   }else T=PI/2.;
  if(verbose){printf("   after  K()\n");fflush(stdout);}
  Omega=2*PI/T;

  N=round(Omega);
  if(um<u0)N=N+1;

  if(N>1.e5)N=0;
  if(N<-1.e5)N=0;
  if(N!=N)N=0;

  x1=(MFNV_C(u,nu*(nx-2)+0,e)+MFNV_C(u,nu*(nx-1)+0,e))/2;
  y1=(MFNV_C(u,nu*(nx-2)+1,e)+MFNV_C(u,nu*(nx-1)+1,e))/2;
  x2=(MFNV_C(u,nu*(nx-2)+2,e)+MFNV_C(u,nu*(nx-1)+2,e))/2;
  y2=(MFNV_C(u,nu*(nx-2)+3,e)+MFNV_C(u,nu*(nx-1)+3,e))/2;
  x3=(MFNV_C(u,nu*(nx-2)+4,e)+MFNV_C(u,nu*(nx-1)+4,e))/2;
  y3=(MFNV_C(u,nu*(nx-2)+5,e)+MFNV_C(u,nu*(nx-1)+5,e))/2;
  del=1.-(x1*y1+x2*y2+x3*y3)/3.1415926;
  if(del<1.)N=-N;

  R=0.;
  Wr=0.;

  d1x=0.;
  d1y=1.;
  d1z=0.;
  x1=(MFNV_C(u,0,e)+MFNV_C(u,nu+0,e))/2;
  y1=(MFNV_C(u,1,e)+MFNV_C(u,nu+1,e))/2;
  x2=(MFNV_C(u,2,e)+MFNV_C(u,nu+2,e))/2;
  y2=(MFNV_C(u,3,e)+MFNV_C(u,nu+3,e))/2;
  x3=(MFNV_C(u,4,e)+MFNV_C(u,nu+4,e))/2;
  y3=(MFNV_C(u,5,e)+MFNV_C(u,nu+5,e))/2;
  kp0=t*x1*y1+(t*x2-c)*y2;
  for(i=1;i<nx-1;i++)
   {
    x1=(MFNV_C(u,nu*i+0,e)+MFNV_C(u,nu*(i+1)+0,e))/2;
    y1=(MFNV_C(u,nu*i+1,e)+MFNV_C(u,nu*(i+1)+1,e))/2;
    x2=(MFNV_C(u,nu*i+2,e)+MFNV_C(u,nu*(i+1)+2,e))/2;
    y2=(MFNV_C(u,nu*i+3,e)+MFNV_C(u,nu*(i+1)+3,e))/2;
    x3=(MFNV_C(u,nu*i+4,e)+MFNV_C(u,nu*(i+1)+4,e))/2;
    y3=(MFNV_C(u,nu*i+5,e)+MFNV_C(u,nu*(i+1)+5,e))/2;

    kp=t*x1*y1+(t*x2-c)*y2;
    if(kp*kp0<=0.)
     {
      if(i==1)
       {
        if(kp>0.)min+=1;
        if(kp<0.)max+=1;
       }else{
        if(kp>0.)min+=2;
        if(kp<0.)max+=2;
       }
     }
    kp0=kp;
   }

  x[0]=mx0;
  x[1]=mz;
  x[2]=f;
  x[3]=theta0;
  x[4]=x1*y1+x2*y2+x3*y3;
  x[5]=up;
  x[6]=um;
  x[7]=N;
  x[8]=Omega;

  for(i=0;i<9;i++)if(x[i]!=x[i])x[i]=0.;

  if(verbose){printf("done %s\n",RoutineName);fflush(stdout);}
  return 0;
 }

double MzPositive(MFNVector u, MFErrorHandler e)
 {
  double x1,y1,x2,y2,x3,y3;
  double d,Mx0,Mz,f;
  int nx=150;
  int np=4;
  int nu=6;

  Mz    =MFNV_C(u,nx*nu+0,e);
  Mx0   =MFNV_C(u,nx*nu+1,e);
  f     =MFNV_C(u,nx*nu+2,e);
  x1=(MFNV_C(u,nu*(nx-2)+0,e)+MFNV_C(u,nu*(nx-1)+0,e))/2;
  y1=(MFNV_C(u,nu*(nx-2)+1,e)+MFNV_C(u,nu*(nx-1)+1,e))/2;
  x2=(MFNV_C(u,nu*(nx-2)+2,e)+MFNV_C(u,nu*(nx-1)+2,e))/2;
  y2=(MFNV_C(u,nu*(nx-2)+3,e)+MFNV_C(u,nu*(nx-1)+3,e))/2;
  x3=(MFNV_C(u,nu*(nx-2)+4,e)+MFNV_C(u,nu*(nx-1)+4,e))/2;
  y3=(MFNV_C(u,nu*(nx-2)+5,e)+MFNV_C(u,nu*(nx-1)+5,e))/2;
  d=x1*y1+x2*y2+x3*y3-1;

  return -Mz;
 }

double Mx0Positive(MFNVector u, MFErrorHandler e)
 {
  double x1,y1,x2,y2,x3,y3;
  double d,Mx0,Mz,f;
  int nx=150;
  int np=4;
  int nu=6;

  Mz    =MFNV_C(u,nx*nu+0,e);
  Mx0   =MFNV_C(u,nx*nu+1,e);
  f     =MFNV_C(u,nx*nu+2,e);
  x1=(MFNV_C(u,nu*(nx-2)+0,e)+MFNV_C(u,nu*(nx-1)+0,e))/2;
  y1=(MFNV_C(u,nu*(nx-2)+1,e)+MFNV_C(u,nu*(nx-1)+1,e))/2;
  x2=(MFNV_C(u,nu*(nx-2)+2,e)+MFNV_C(u,nu*(nx-1)+2,e))/2;
  y2=(MFNV_C(u,nu*(nx-2)+3,e)+MFNV_C(u,nu*(nx-1)+3,e))/2;
  x3=(MFNV_C(u,nu*(nx-2)+4,e)+MFNV_C(u,nu*(nx-1)+4,e))/2;
  y3=(MFNV_C(u,nu*(nx-2)+5,e)+MFNV_C(u,nu*(nx-1)+5,e))/2;
  d=x1*y1+x2*y2+x3*y3-1;

  return -Mx0;
 }

double Mx0Negative(MFNVector u, MFErrorHandler e)
 {
  double x1,y1,x2,y2,x3,y3;
  double d,Mx0,Mz,f;
  int nx=150;
  int np=4;
  int nu=6;

  Mz    =MFNV_C(u,nx*nu+0,e);
  Mx0   =MFNV_C(u,nx*nu+1,e);
  f     =MFNV_C(u,nx*nu+2,e);
  x1=(MFNV_C(u,nu*(nx-2)+0,e)+MFNV_C(u,nu*(nx-1)+0,e))/2;
  y1=(MFNV_C(u,nu*(nx-2)+1,e)+MFNV_C(u,nu*(nx-1)+1,e))/2;
  x2=(MFNV_C(u,nu*(nx-2)+2,e)+MFNV_C(u,nu*(nx-1)+2,e))/2;
  y2=(MFNV_C(u,nu*(nx-2)+3,e)+MFNV_C(u,nu*(nx-1)+3,e))/2;
  x3=(MFNV_C(u,nu*(nx-2)+4,e)+MFNV_C(u,nu*(nx-1)+4,e))/2;
  y3=(MFNV_C(u,nu*(nx-2)+5,e)+MFNV_C(u,nu*(nx-1)+5,e))/2;
  d=x1*y1+x2*y2+x3*y3-1;

  return Mx0;
 }

double FPositive(MFNVector u, MFErrorHandler e)
 {
  double x1,y1,x2,y2,x3,y3;
  double d,Mx0,Mz,f;
  int nx=150;
  int np=4;
  int nu=6;

  Mz    =MFNV_C(u,nx*nu+0,e);
  Mx0   =MFNV_C(u,nx*nu+1,e);
  f     =MFNV_C(u,nx*nu+2,e);
  x1=(MFNV_C(u,nu*(nx-2)+0,e)+MFNV_C(u,nu*(nx-1)+0,e))/2;
  y1=(MFNV_C(u,nu*(nx-2)+1,e)+MFNV_C(u,nu*(nx-1)+1,e))/2;
  x2=(MFNV_C(u,nu*(nx-2)+2,e)+MFNV_C(u,nu*(nx-1)+2,e))/2;
  y2=(MFNV_C(u,nu*(nx-2)+3,e)+MFNV_C(u,nu*(nx-1)+3,e))/2;
  x3=(MFNV_C(u,nu*(nx-2)+4,e)+MFNV_C(u,nu*(nx-1)+4,e))/2;
  y3=(MFNV_C(u,nu*(nx-2)+5,e)+MFNV_C(u,nu*(nx-1)+5,e))/2;
  d=x1*y1+x2*y2+x3*y3-1;

  return -f;
 }

double DNegative(MFNVector u, MFErrorHandler e)
 {
  double x1,y1,x2,y2,x3,y3;
  double d,Mx0,Mz,f;
  int nx=150;
  int np=4;
  int nu=6;

  Mz    =MFNV_C(u,nx*nu+0,e);
  Mx0   =MFNV_C(u,nx*nu+1,e);
  f     =MFNV_C(u,nx*nu+2,e);
  x1=(MFNV_C(u,nu*(nx-2)+0,e)+MFNV_C(u,nu*(nx-1)+0,e))/2;
  y1=(MFNV_C(u,nu*(nx-2)+1,e)+MFNV_C(u,nu*(nx-1)+1,e))/2;
  x2=(MFNV_C(u,nu*(nx-2)+2,e)+MFNV_C(u,nu*(nx-1)+2,e))/2;
  y2=(MFNV_C(u,nu*(nx-2)+3,e)+MFNV_C(u,nu*(nx-1)+3,e))/2;
  x3=(MFNV_C(u,nu*(nx-2)+4,e)+MFNV_C(u,nu*(nx-1)+4,e))/2;
  y3=(MFNV_C(u,nu*(nx-2)+5,e)+MFNV_C(u,nu*(nx-1)+5,e))/2;
  d=1-(x1*y1+x2*y2+x3*y3)/3.1415926;

  return d-1;
 }

double DPositive(MFNVector u, MFErrorHandler e)
 {
  double x1,y1,x2,y2,x3,y3;
  double d,Mx0,Mz,f;
  int nx=150;
  int np=4;
  int nu=6;

  Mz    =MFNV_C(u,nx*nu+0,e);
  Mx0   =MFNV_C(u,nx*nu+1,e);
  f     =MFNV_C(u,nx*nu+2,e);
  x1=(MFNV_C(u,nu*(nx-2)+0,e)+MFNV_C(u,nu*(nx-1)+0,e))/2;
  y1=(MFNV_C(u,nu*(nx-2)+1,e)+MFNV_C(u,nu*(nx-1)+1,e))/2;
  x2=(MFNV_C(u,nu*(nx-2)+2,e)+MFNV_C(u,nu*(nx-1)+2,e))/2;
  y2=(MFNV_C(u,nu*(nx-2)+3,e)+MFNV_C(u,nu*(nx-1)+3,e))/2;
  x3=(MFNV_C(u,nu*(nx-2)+4,e)+MFNV_C(u,nu*(nx-1)+4,e))/2;
  y3=(MFNV_C(u,nu*(nx-2)+5,e)+MFNV_C(u,nu*(nx-1)+5,e))/2;
  d=1-(x1*y1+x2*y2+x3*y3)/3.1415926;

  return 1-d;
 }

double FLessThanSixteen(MFNVector u, MFErrorHandler e)
 {
  double x1,y1,x2,y2,x3,y3;
  double d,Mx0,Mz,f;
  int nx=150;
  int np=4;
  int nu=6;

  Mz    =MFNV_C(u,nx*nu+0,e);
  Mx0   =MFNV_C(u,nx*nu+1,e);
  f     =MFNV_C(u,nx*nu+2,e);
  x1=(MFNV_C(u,nu*(nx-2)+0,e)+MFNV_C(u,nu*(nx-1)+0,e))/2;
  y1=(MFNV_C(u,nu*(nx-2)+1,e)+MFNV_C(u,nu*(nx-1)+1,e))/2;
  x2=(MFNV_C(u,nu*(nx-2)+2,e)+MFNV_C(u,nu*(nx-1)+2,e))/2;
  y2=(MFNV_C(u,nu*(nx-2)+3,e)+MFNV_C(u,nu*(nx-1)+3,e))/2;
  x3=(MFNV_C(u,nu*(nx-2)+4,e)+MFNV_C(u,nu*(nx-1)+4,e))/2;
  y3=(MFNV_C(u,nu*(nx-2)+5,e)+MFNV_C(u,nu*(nx-1)+5,e))/2;
  d=x1*y1+x2*y2+x3*y3-1;

  return f-16.;
 }

double MzNegative(MFNVector u, MFErrorHandler e)
 {
  double x1,y1,x2,y2,x3,y3;
  double d,Mx0,Mz,f;
  int nx=150;
  int np=4;
  int nu=6;

  Mz    =MFNV_C(u,nx*nu+0,e);
  Mx0   =MFNV_C(u,nx*nu+1,e);
  f     =MFNV_C(u,nx*nu+2,e);
  x1=(MFNV_C(u,nu*(nx-2)+0,e)+MFNV_C(u,nu*(nx-1)+0,e))/2;
  y1=(MFNV_C(u,nu*(nx-2)+1,e)+MFNV_C(u,nu*(nx-1)+1,e))/2;
  x2=(MFNV_C(u,nu*(nx-2)+2,e)+MFNV_C(u,nu*(nx-1)+2,e))/2;
  y2=(MFNV_C(u,nu*(nx-2)+3,e)+MFNV_C(u,nu*(nx-1)+3,e))/2;
  x3=(MFNV_C(u,nu*(nx-2)+4,e)+MFNV_C(u,nu*(nx-1)+4,e))/2;
  y3=(MFNV_C(u,nu*(nx-2)+5,e)+MFNV_C(u,nu*(nx-1)+5,e))/2;
  d=x1*y1+x2*y2+x3*y3-1;

  return Mz;
 }

double MzSmallPositive(MFNVector u, MFErrorHandler e)
 {
  double x1,y1,x2,y2,x3,y3;
  double d,Mx0,Mz,f;
  int nx=150;
  int np=4;
  int nu=6;

  Mz    =MFNV_C(u,nx*nu+0,e);
  Mx0   =MFNV_C(u,nx*nu+1,e);
  f     =MFNV_C(u,nx*nu+2,e);
  x1=(MFNV_C(u,nu*(nx-2)+0,e)+MFNV_C(u,nu*(nx-1)+0,e))/2;
  y1=(MFNV_C(u,nu*(nx-2)+1,e)+MFNV_C(u,nu*(nx-1)+1,e))/2;
  x2=(MFNV_C(u,nu*(nx-2)+2,e)+MFNV_C(u,nu*(nx-1)+2,e))/2;
  y2=(MFNV_C(u,nu*(nx-2)+3,e)+MFNV_C(u,nu*(nx-1)+3,e))/2;
  x3=(MFNV_C(u,nu*(nx-2)+4,e)+MFNV_C(u,nu*(nx-1)+4,e))/2;
  y3=(MFNV_C(u,nu*(nx-2)+5,e)+MFNV_C(u,nu*(nx-1)+5,e))/2;
  d=x1*y1+x2*y2+x3*y3-1;

  return Mz-.05;
 }

double MzSmallNegative(MFNVector u, MFErrorHandler e)
 {
  double x1,y1,x2,y2,x3,y3;
  double d,Mx0,Mz,f;
  int nx=150;
  int np=4;
  int nu=6;

  Mz    =MFNV_C(u,nx*nu+0,e);
  Mx0   =MFNV_C(u,nx*nu+1,e);
  f     =MFNV_C(u,nx*nu+2,e);
  x1=(MFNV_C(u,nu*(nx-2)+0,e)+MFNV_C(u,nu*(nx-1)+0,e))/2;
  y1=(MFNV_C(u,nu*(nx-2)+1,e)+MFNV_C(u,nu*(nx-1)+1,e))/2;
  x2=(MFNV_C(u,nu*(nx-2)+2,e)+MFNV_C(u,nu*(nx-1)+2,e))/2;
  y2=(MFNV_C(u,nu*(nx-2)+3,e)+MFNV_C(u,nu*(nx-1)+3,e))/2;
  x3=(MFNV_C(u,nu*(nx-2)+4,e)+MFNV_C(u,nu*(nx-1)+4,e))/2;
  y3=(MFNV_C(u,nu*(nx-2)+5,e)+MFNV_C(u,nu*(nx-1)+5,e))/2;
  d=x1*y1+x2*y2+x3*y3-1;

  return -Mz-.05;
 }

double DSmallPositive(MFNVector u, MFErrorHandler e)
 {
  double x1,y1,x2,y2,x3,y3;
  double d,Mx0,Mz,f;
  int nx=150;
  int np=4;
  int nu=6;

  Mz    =MFNV_C(u,nx*nu+0,e);
  Mx0   =MFNV_C(u,nx*nu+1,e);
  f     =MFNV_C(u,nx*nu+2,e);
  x1=(MFNV_C(u,nu*(nx-2)+0,e)+MFNV_C(u,nu*(nx-1)+0,e))/2;
  y1=(MFNV_C(u,nu*(nx-2)+1,e)+MFNV_C(u,nu*(nx-1)+1,e))/2;
  x2=(MFNV_C(u,nu*(nx-2)+2,e)+MFNV_C(u,nu*(nx-1)+2,e))/2;
  y2=(MFNV_C(u,nu*(nx-2)+3,e)+MFNV_C(u,nu*(nx-1)+3,e))/2;
  x3=(MFNV_C(u,nu*(nx-2)+4,e)+MFNV_C(u,nu*(nx-1)+4,e))/2;
  y3=(MFNV_C(u,nu*(nx-2)+5,e)+MFNV_C(u,nu*(nx-1)+5,e))/2;
  d=1-(x1*y1+x2*y2+x3*y3)/3.1415926;

  return d-1-.1;
 }

double DSmallNegative(MFNVector u, MFErrorHandler e)
 {
  double x1,y1,x2,y2,x3,y3;
  double d,Mx0,Mz,f;
  int nx=150;
  int np=4;
  int nu=6;

  Mz    =MFNV_C(u,nx*nu+0,e);
  Mx0   =MFNV_C(u,nx*nu+1,e);
  f     =MFNV_C(u,nx*nu+2,e);
  x1=(MFNV_C(u,nu*(nx-2)+0,e)+MFNV_C(u,nu*(nx-1)+0,e))/2;
  y1=(MFNV_C(u,nu*(nx-2)+1,e)+MFNV_C(u,nu*(nx-1)+1,e))/2;
  x2=(MFNV_C(u,nu*(nx-2)+2,e)+MFNV_C(u,nu*(nx-1)+2,e))/2;
  y2=(MFNV_C(u,nu*(nx-2)+3,e)+MFNV_C(u,nu*(nx-1)+3,e))/2;
  x3=(MFNV_C(u,nu*(nx-2)+4,e)+MFNV_C(u,nu*(nx-1)+4,e))/2;
  y3=(MFNV_C(u,nu*(nx-2)+5,e)+MFNV_C(u,nu*(nx-1)+5,e))/2;
  d=1-(x1*y1+x2*y2+x3*y3)/3.1415926;

  return 1-d-.1;
 }

double DSmallerPositive(MFNVector u, MFErrorHandler e)
 {
  double x1,y1,x2,y2,x3,y3;
  double d,Mx0,Mz,f;
  int nx=150;
  int np=4;
  int nu=6;

  Mz    =MFNV_C(u,nx*nu+0,e);
  Mx0   =MFNV_C(u,nx*nu+1,e);
  f     =MFNV_C(u,nx*nu+2,e);
  x1=(MFNV_C(u,nu*(nx-2)+0,e)+MFNV_C(u,nu*(nx-1)+0,e))/2;
  y1=(MFNV_C(u,nu*(nx-2)+1,e)+MFNV_C(u,nu*(nx-1)+1,e))/2;
  x2=(MFNV_C(u,nu*(nx-2)+2,e)+MFNV_C(u,nu*(nx-1)+2,e))/2;
  y2=(MFNV_C(u,nu*(nx-2)+3,e)+MFNV_C(u,nu*(nx-1)+3,e))/2;
  x3=(MFNV_C(u,nu*(nx-2)+4,e)+MFNV_C(u,nu*(nx-1)+4,e))/2;
  y3=(MFNV_C(u,nu*(nx-2)+5,e)+MFNV_C(u,nu*(nx-1)+5,e))/2;
  d=1-(x1*y1+x2*y2+x3*y3)/3.1415926;

  return d-1-.05;
 }

double DSmallerNegative(MFNVector u, MFErrorHandler e)
 {
  double x1,y1,x2,y2,x3,y3;
  double d,Mx0,Mz,f;
  int nx=150;
  int np=4;
  int nu=6;

  Mz    =MFNV_C(u,nx*nu+0,e);
  Mx0   =MFNV_C(u,nx*nu+1,e);
  f     =MFNV_C(u,nx*nu+2,e);
  x1=(MFNV_C(u,nu*(nx-2)+0,e)+MFNV_C(u,nu*(nx-1)+0,e))/2;
  y1=(MFNV_C(u,nu*(nx-2)+1,e)+MFNV_C(u,nu*(nx-1)+1,e))/2;
  x2=(MFNV_C(u,nu*(nx-2)+2,e)+MFNV_C(u,nu*(nx-1)+2,e))/2;
  y2=(MFNV_C(u,nu*(nx-2)+3,e)+MFNV_C(u,nu*(nx-1)+3,e))/2;
  x3=(MFNV_C(u,nu*(nx-2)+4,e)+MFNV_C(u,nu*(nx-1)+4,e))/2;
  y3=(MFNV_C(u,nu*(nx-2)+5,e)+MFNV_C(u,nu*(nx-1)+5,e))/2;
  d=1-(x1*y1+x2*y2+x3*y3)/3.1415926;

  return 1-d-.05;
 }

double K(double p)
 {
/* This from Spanier and Oldham "An Atlas of Functions" page 61:8 */
  double g,a,t,e,s,k;

  if(p<.067)return .5*PI*(16.-5*p*p)/(16-9*p*p);
  if(1.-p<1.e-15)return 1.e20;

  g=sqrt(1-p*p);
  a=1.;
  t=1.;
  e=1.+g*g;
  s=g*a;
  while(t*(a*a-s)>6e-12)
   {
    s=g*a;
    a=(a+g)/2;
    g=sqrt(s);
    t=2*t;
    e-=t*(a*a-s);
   }
  k=.5*PI/g;

  return k;
 }

double Ki(double p)
 {
  double x;
  double result;

  x=p/sqrt(1+p*p);
  result=K(x)/sqrt(1.+x*x);

  return result;
 }
