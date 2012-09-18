/* 
    %W%
    %D% %T%
   
    PROGRAM NAME:  multifario

    (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
    CORPORATION 12/1/2003.  ALL RIGHTS RESERVED.

    Please refer to the LICENSE file in the top directory

*/
/*      author: Mike Henderson mhender@watson.ibm.com */

#include <MFAtlas.h>
#include <MFTPBVP.h>
#include <math.h>
#include <MFMultifariosMethod.h>

int MFDomokosProjectToDraw(MFNVector,double*,void*,MFErrorHandler);
void MFTPBVPSetStability(MFImplicitMF,MFNVector,MFNKMatrix,void*,MFErrorHandler);
int MFStopTPBVP(MFImplicitMF,MFNVector,MFNKMatrix,MFNVector,MFNKMatrix,void*,MFErrorHandler);

#define PI 3.14159265358979323846264338327950288

/*
u'=f(u,t,l);                     nu differential equations
a(u(0),u(1),l)=0;                nbc=nu equations
int l(u(t),t,l) dt + m(l) = 0;   nc equations
*/

/*   M=createTPBVP("[a,M,y]","[H,V]","[M,-H*sin(a)+V*cos(a),sin(a)]","[a(0),y(0),a(1),y(1)]"); */

void f(double r, int nu, double *u, int np, double *p, double *u0, double *l0, double *f,MFErrorHandler e)
 {
/*  (alpha,M,y;H,V)  */

  f[0]=2*PI*u[1];
  f[1]=-p[0]*sin(2*PI*u[0])+p[1]*cos(2*PI*u[0]);
  f[2]=sin(2*PI*u[0]);

  return;
 }

void fu(double r, int nu, double *u, int np, double *p, double *u0, double *l0, double *fu,MFErrorHandler e)
 {
  int i;

  for(i=0;i<nu*nu;i++)fu[i]=0.;

  fu[0+nu*1]=2.*PI;
  fu[1+nu*0]=-2*PI*p[0]*cos(2*PI*u[0])-2*PI*p[1]*sin(2*PI*u[0]);
  fu[2+nu*0]=2*PI*cos(2*PI*u[0]);

  return;
 }

void fl(double r, int nu, double *u, int np, double *p, double *u0, double *l0, double *fl,MFErrorHandler e)
 {
  int i;

  for(i=0;i<nu*np;i++)fl[i]=0.;

  fl[1+nu*0]=-sin(2*PI*u[0]);
  fl[1+nu*1]= cos(2*PI*u[0]);

  return;
 }

void a(int nbc, int nu, double *uL, double *uR, int np, double *p, double *u0L, double *u0R,double *l0, double *a,MFErrorHandler e)
 {
/*  (alpha,M,y;H,V)  */

  a[0]=uL[0];
  a[1]=uR[0];
  a[2]=uL[2];
  a[3]=uR[2];

  return;
 }

void au(int nbc, int nu, double *uL, double *uR, int np, double *p, double *u0L, double *u0R,double *l0, double *au,MFErrorHandler e)
 {
  int i;

  for(i=0;i<nbc*2*nu;i++)au[i]=0.;

  au[0+nbc*(0+0*nu)]=1.;
  au[1+nbc*(0+1*nu)]=1.;
  au[2+nbc*(2+0*nu)]=1.;
  au[3+nbc*(2+1*nu)]=1.;

  return;
 }

void al(int nbc, int nu, double *uL, double *uR, int np, double *p, double *u0L, double *u0R,double *l0, double *al,MFErrorHandler e)
 {
  int i;

  for(i=0;i<nbc*np;i++)al[i]=0.;

  return;
 }

void l(int nic, double r, int nu, double *u, int np, double *p,double *u0,double *l0, double *l,MFErrorHandler e)
 {
  int i;

  for(i=0;i<nic;i++)l[i]=0.;

  return;
 }

void lu(int nic, double r, int nu, double *u, int np, double *p,double *u0,double *l0, double *lu,MFErrorHandler e)
 {
  int i;

  for(i=0;i<nu*nic;i++)lu[i]=0.;

  return;
 }

void ll(int nic, double r, int nu, double *u, int np, double *p,double *u0,double *l0, double *ll,MFErrorHandler e)
 {
  int i;

  for(i=0;i<np*nic;i++)ll[i]=0.;

  return;
 }

void m(int nic, int np, double *p, double *p0, double *m,MFErrorHandler e)
 {
  int i;

  for(i=0;i<nic;i++)m[i]=0.;

  return;
 }

void ml(int nic, int np, double *p, double *p0, double *ml,MFErrorHandler e)
 {
  int i;

  for(i=0;i<np*nic;i++)ml[i]=0.;

  return;
 }

#define NX 100

int main(int argc, char *argv[])
 {
  int i,j,n;
  MFImplicitMF M;
  MFNRegion Omega;
  MFAtlas S;
  MFNVector ug;
  MFNVector u0;
  MFNKMatrix Tan;
  MFContinuationMethod H;
  MFErrorHandler e;
  double p0[2];
  double p1[2];
  double *r0;
  double dr;
  double xy0[3];

  int nx=NX; /* 150 */
  int np=2;
  int nu=3;
  int nbc=4;
  int nic=0;
  int k;

  e=MFCreateErrorHandler();

  k=nu+np-nbc-nic;
  
  M=MFIMFCreateTPBVP(k,nx,nu,np,f,fu,fl,nbc,a,au,al,nic,l,lu,ll,m,ml,e);
  MFIMFSetR(M,6./(4*PI*PI),e); /* (because the trivial branch has no scale! */
  MFIMFSetR(M,.06/(4*PI*PI),e); /* (because the trivial branch has no scale! */
  n=MFIMF_N(M,e);
  MFIMFSetProjectForDraw(M,MFDomokosProjectToDraw,e);
  MFIMFSetSetStability(M,MFTPBVPSetStability,e);
  MFIMFSetStop(M,MFStopTPBVP,e);

  p0[0]= -600./(4*PI*PI); p1[0]=600./(4*PI*PI);   /* H */
  p0[1]= -20.; p1[1]=20.;   /* V */
  Omega=MFNRegionCreateTPBVP(nx,nu,np,p0,p1,-200.,200.,e);

  r0=(double*)malloc((nx+1)*sizeof(double));
  dr=1./(nx-2);
  for(i=0;i<nx+1;i++)r0[i]=(i-.5)*dr;

  xy0[0]=0.;
  xy0[1]=0.;
  xy0[2]=0.;
  p0[0]=1.15;
  p0[1]=0.;

  printf("  parameters   (%lf,%lf)\n",p0[0],p0[1]);fflush(stdout);

  printf("Construct the initial guess \n");fflush(stdout);
  ug=MFTPBVPIntegrateForInitialSolution(M,xy0,p0,r0,e);

  Tan=MFIMFTangentSpace(M,ug,e);
  if(Tan==(MFNKMatrix)NULL)
   {
    fprintf(stderr,"Failed to find tangent at guess\n");
    return 8;
   }

  u0=MFIMFVectorFactory(M,e);
  if(!MFIMFProject(M,ug,Tan,u0,e))
   {
    fprintf(stderr,"Initial guess failed to project\n");
    return 8;
   }
  MFFreeNKMatrix(Tan,e);
  MFFreeNVector(ug,e);

  H=MFCreateMultifariosMethod(e);
  MFMultifarioSetRealParameter(H,"epsilon",.03,e);
  MFMultifarioSetRealParameter(H,"dotmin",.9,e);
  MFMultifarioSetIntegerParameter(H,"maxCharts",-1,e);
  MFMultifarioSetIntegerParameter(H,"verbose",1,e);
  MFMultifarioSetIntegerParameter(H,"page",1,e);
  MFMultifarioSetIntegerParameter(H,"branchSwitch",10,e);
  MFMultifarioSetFilename(H,"Domokos",e);

  printf("ComputeAtlas\n");fflush(stdout);
  S=MFComputeAtlas(H,M,Omega,u0,e);
  MFCloseAtlas(H,S,e);

  MFFreeAtlas(S,e);
  MFFreeContinuationMethod(H,e);
  MFFreeImplicitMF(M,e);
  MFFreeNRegion(Omega,e);
  MFFreeNVector(u0,e);
  MFFreeErrorHandler(e);

  return(0);
 }

int MFDomokosProjectToDraw(MFNVector u, double *x, void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFDomokosProjectToDraw"};
  int n;
  int nx=NX;
  int np=2;
  int nu=3;

  if(x==(double*)NULL)return 3;

  n=MFNV_NC(u,e);
  x[0]=MFNV_C(u,nx*nu,e);
  x[2]=MFNV_C(u,nx*nu+1,e);
  x[1]=MFNV_C(u,1,e);

  return 3;
 }
