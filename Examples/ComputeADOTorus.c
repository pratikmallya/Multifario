/*
 *  PROGRAM NAME:  multifario
 *
 *  (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
 *  CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
 *
 *  Please refer to the LICENSE file in the top directory
 *
 *      author: Mike Henderson mhender@watson.ibm.com
 *      date:   November 11, 1997
 *              February 2, 1999   Ported to C
 *              April 29, 2002     Made into a base class
 */

static char *id="@(#) $Id: ComputeADOTorus.c,v 1.2 2011/07/21 17:43:45 mhender Exp $";

static char MFNSpaceErrorMsg[256]="";

#include <multifarioConfig.h>
#include <math.h>
#include <stdio.h>

#include <MFAtlas.h>
#include <MFAtlasFriends.h>
#include <MFChart.h>
#include <MFErrorHandler.h>
#include <MFPrint.h>
#include <MFBinaryTree.h>
#include <MFMultifariosMethod.h>

#include <IMF.h>
#include <IMFFlow.h>
#include <IMFFlat.h>
#include <IMFIntegrateFat.h>
#include <IMFInterpolation.h>
#include <IMFExpansion.h>
#include <IMFExpansionSpace.h>
#include <IMFExpansionPt.h>

#ifdef __cplusplus
 extern "C" {
#endif

static char IMFIntegrateFatErrorMsg[256];
#define round(x) ((x)>=0?(long)((x)+0.5):(long)((x)-0.5))

/* #define DOINTERPANIM */
#define DOTRAJ

extern int IMFNInterper;
extern FILE *IMFInterper;
extern int IMFNInterpee;
extern FILE *IMFInterpee;
extern int IMFNInterpT;
extern FILE *IMFInterpT;
extern int IMFNCircle;
extern FILE *IMFCircle;
extern int IMFNTraj;
extern FILE *IMFTraj;
FILE *IMFInterp;
int IMFNInterp;
FILE *IMFBC;
int IMFNBC;

void         MFChartSetPaged(MFChart,MFErrorHandler);
void         MFAtlasPageOutChartsNotNearBoundary(MFAtlas,int,int,char*,MFErrorHandler);
void         MFAtlasSetNearRtn(MFAtlas,int (*)(MFAtlas,MFChart,MFChart,MFErrorHandler),MFErrorHandler);
MFChart      MFAtlasChart(MFAtlas,int,MFErrorHandler);
MFBinaryTree MFAtlasGetBB(MFAtlas,MFErrorHandler);

MFNVector IMFGetInterpolationPointOnList(MFAtlas,IMFFlow,MFKVector,MFAtlas,double,MFNRegion,int,int*,MFErrorHandler);
IMFFlow   IMFCreateADOTorusFlow(MFErrorHandler);
MFNVector IMFGetInterpolationPointOnList(MFAtlas,IMFFlow,MFKVector,MFAtlas,double,MFNRegion,int,int*,MFErrorHandler);
MFAtlas   IMFComputeCovering(IMFFlow,int,MFNVector,MFKVector,MFNRegion,double,double,double,int,int,double,double,MFErrorHandler);

int     CIMFIsNear(MFAtlas,MFChart,MFChart,MFErrorHandler);
MFChart CIMFExtendAtlasAlongFatTraj(MFAtlas I, IMFFlow F, char *name, MFNVector u0, MFKVector p0, double dt, double ti, double tf,double epsilon,MFNRegion Omega, int maxSteps, int maxCharts,double Rmax,int nSkip, int skipCheck, MFErrorHandler e);
MFAtlas CIMFComputeInvariantManifold(IMFFlow F,MFKVector p0,char *name, MFAtlas c,MFNRegion Omega, double eps, double dt, double tmax, int maxInterp, int maxCharts, double Rmax,MFErrorHandler e);
void CIMFGetBoundaryConditions(MFAtlas I,MFChart chart,int nT, MFChart *T, MFErrorHandler e);

int ADOTorusProject(MFNVector u, double *x, void *d, MFErrorHandler e);
int ADOTorusProjectBB(MFNVector u, double *x, void *d, MFErrorHandler e);

static int MEHKellerInt(IMFFlow,double*,MFKVector,double,double,int,MFErrorHandler);
static void ADOTorus(double *u, double *p, double *F,  void *data, MFErrorHandler e);
static void dADOTorus(double *u, double *p, double *dF,  void *data, MFErrorHandler e);
static void dADOTorusdp(double *u, double *p, double *dF,  void *data, MFErrorHandler e);
static void ddADOTorus(double *u, double *p, double *dF,  void *data, MFErrorHandler e);

/*-------------------------------------------------------------------------------------------------------*/

int main(int argc, char *argv[])
 {
  MFAtlas A;

  IMFFlow F;
  MFNVector u0;
  MFKVector p0;
  MFNRegion Omega;
  double eps;
  double dt;
  double tmax;
  int maxInterp;
  int maxCharts;
  double R0;
  double Rmax;
  MFErrorHandler e;
  MFNVector ll,ur;
  int torusNumber;

  e=MFCreateErrorHandler();
  F=IMFCreateADOTorusFlow(e);
  if(argc<2)torusNumber=10;
   else     torusNumber=atoi(argv[1]);

  ll=MFCreateNVector(4,e);
  ur=MFCreateNVector(4,e);
  MFNVSetC(ll,0,-2.,e);MFNVSetC(ur,0,2.,e);
  MFNVSetC(ll,1,-2.,e);MFNVSetC(ur,1,2.,e);
  MFNVSetC(ll,2,-2.,e);MFNVSetC(ur,2,2.,e);
  MFNVSetC(ll,3,-2.,e);MFNVSetC(ur,3,2.,e);
  Omega=MFNRegionCreateHyperCubeByCorners(4,ll,ur,e);
  MFFreeNVector(ll,e);
  MFFreeNVector(ur,e);

  u0=MFFlowNVectorFactory(F,e);
  p0=MFFlowKVectorFactory(F,e);

  MFNVSetC(u0,0,1.,e);
  MFNVSetC(u0,1,0.,e);
  MFNVSetC(u0,2,1.,e);
  MFNVSetC(u0,3,0.,e);
  MFKVSetC(p0,0,0.25359,e);

  tmax=24.;
  maxInterp=100000;
  maxCharts=152000;  /* 005 */
  maxCharts=101900;  /* 001 */

  R0=0.04;
  Rmax=.04;
  eps=2.e-2;
  dt=1.e-4;
  A=IMFComputeCovering(F,torusNumber,u0,p0,Omega,eps,dt, tmax, maxInterp, maxCharts, R0, Rmax, e);
  if(A!=NULL){printf("A has %d charts\n",MFAtlasNumberOfCharts(A,e));fflush(stdout);}
  
  if(A!=NULL)MFFreeAtlas(A,e);
  IMFFreeFlow(F,e);
  MFFreeKVector(p0,e);
  MFFreeNRegion(Omega,e);
  MFFreeErrorHandler(e);

  return 0;
 }

MFAtlas IMFComputeCovering(IMFFlow F,int nTorus, MFNVector v0, MFKVector p0, MFNRegion Omega, double eps, double dt, double tmax, int maxInterp, int maxCharts, double R0, double Rmax, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFComputeCovering"};
  MFAtlas I;
  MFAtlas c;
  MFContinuationMethod H;
  MFImplicitMF M;
  FILE *fTorus;
  int nc,nx,ny;
  double lambda;
  int i,j,n,m;
  double **x;
  double a,b;
  MFNVector u0;
  char fileName[]="/home/mhender/tori/ado/ado010.data";

/*  4   70   70    .25359 */
/*  1    1   1.00002    .00000   1.00002    .00000   -.13765   -.13765 */

  sprintf(fileName,"/home/mhender/tori/ado/ado%3.3d.data",nTorus);
  fTorus=fopen(fileName,"r");

  if(fTorus==NULL)
   {
    printf("Could not open %s\n",fileName);fflush(stdout);
    return NULL;
   }
  fscanf(fTorus,"%d %d %d %lf\n",&nc,&nx,&ny,&lambda);

  MFKVSetC(p0,0,lambda,e);
  printf("parameter value is %lf\n",lambda);fflush(stdout);

  printf("Read:\n");fflush(stdout);
  x=(double**)malloc(nc*sizeof(double*));
  for(j=0;j<nc;j++)x[j]=(double*)malloc(nx*sizeof(double));

  for(i=0;i<nx;i++)
   {
    fscanf(fTorus,"%d %d ",&n,&m);
    for(j=0;j<nc;j++)fscanf(fTorus,"%lf",&(x[j][i]));
    fscanf(fTorus,"%lf %lf\n",&a,&b);
   }
  fclose(fTorus);

/* Fit circular spline to the computed data */

  M=MFIMFCreatePeriodicSpline(nx,nc,x,e);
  MFIMFSetProjectForDraw(M,ADOTorusProject,e);
  MFIMFSetProjectForBB(M,ADOTorusProjectBB,e);
  MFIMFSetR(M,R0,e);

  u0=MFIMFVectorFactory(M,e);
  for(j=0;j<nc;j++)MFNVSetC(u0,j,x[j][0],e);

  H=MFCreateMultifariosMethod(e);
  MFMultifarioSetRealParameter(H,"epsilon",eps,e);
  MFMultifarioSetIntegerParameter(H,"maxCharts",-1,e);
  MFMultifarioSetRealParameter(H,"maxR",Rmax,e);
  MFMultifarioSetIntegerParameter(H,"verbose",1,e);
  MFMultifarioSetIntegerParameter(H,"page",0,e);
  MFMultifarioSetIntegerParameter(H,"dumpToPlotFile",1,e);
  MFMultifarioSetIntegerParameter(H,"dumpToCenterFile",0,e);
  sprintf(fileName,"ADO%3.3dInitialSection",nTorus);
  MFMultifarioSetFilename(H,fileName,e);

  c=MFComputeAtlas(H,M,Omega,u0,e);
  MFCloseAtlas(H,c,e);

  sprintf(fileName,"ADO%3.3dTorus",nTorus);
  I=CIMFComputeInvariantManifold(F,p0,fileName,c,Omega,eps,dt,tmax,maxInterp,maxCharts,Rmax,e);

  MFFreeNVector(u0,e);
  MFFreeAtlas(c,e);
  MFFreeImplicitMF(M,e);
  MFFreeContinuationMethod(H,e);

  return I;
 }

IMFFlow IMFCreateADOTorusFlow(MFErrorHandler e)
 {
  static char RoutineName[]={"IMFCreateADOTorusFlow"};
  return IMFCreateFlow(4,1,ADOTorus,dADOTorus,dADOTorusdp,ddADOTorus,NULL,NULL,NULL,e);
 }

static void ADOTorus(double *u, double *p, double *F,  void *data, MFErrorHandler e)
 {
  static char RoutineName[]={"ADOTorus"};
  static double a1=1.;
  static double a2=1.;
  static double b1=.55;
  static double b2=.55;

  F[0] = a1*u[0]+b1*u[1]-(u[0]*u[0]+u[1]*u[1])*u[0]-p[0]*(u[0]+u[1]-u[2]-u[3]);
  F[1] =-b1*u[0]+a1*u[1]-(u[0]*u[0]+u[1]*u[1])*u[1]-p[0]*(u[0]+u[1]-u[2]-u[3]);
  F[2] = a2*u[2]+b2*u[3]-(u[2]*u[2]+u[3]*u[3])*u[2]+p[0]*(u[0]+u[1]-u[2]-u[3]);
  F[3] =-b2*u[2]+a2*u[3]-(u[2]*u[2]+u[3]*u[3])*u[3]+p[0]*(u[0]+u[1]-u[2]-u[3]);

  return;
 }

static void dADOTorus(double *u, double *p, double *dF,  void *data, MFErrorHandler e)
 {
  static char RoutineName[]={"dADOTorus"};
  static double a1=1.;
  static double a2=1.;
  static double b1=.55;
  static double b2=.55;

  dF[0+4*0]= a1-(3*u[0]*u[0]+u[1]*u[1])-p[0];
  dF[0+4*1]= b1-2*u[1]*u[0]-p[0];
  dF[0+4*2]= p[0];
  dF[0+4*3]= p[0];

  dF[1+4*0]=-b1-2*u[0]*u[1]-p[0];
  dF[1+4*1]=a1-(u[0]*u[0]+3*u[1]*u[1])-p[0];
  dF[1+4*2]=p[0];
  dF[1+4*3]=p[0];

  dF[2+4*0]= p[0];
  dF[2+4*1]= p[0];
  dF[2+4*2]= a2-(3*u[2]*u[2]+u[3]*u[3])-p[0];
  dF[2+4*3]=b2-2*u[3]*u[2];

  dF[3+4*0]=p[0];
  dF[3+4*1]=p[0];
  dF[3+4*2]=-b2-2*u[2]*u[3]-p[0];
  dF[3+4*3]=a2-(u[2]*u[2]+3*u[3]*u[3])-p[0];

  return;
 }

static void dADOTorusdp(double *u, double *p, double *dp,  void *data, MFErrorHandler e)
 {
  static char RoutineName[]={"dADOTorusdp"};
  static double a1=1.;
  static double a2=1.;
  static double b1=.55;
  static double b2=.55;

  dp[0]=-(u[0]+u[1]-u[2]-u[3]);
  dp[1]=-(u[0]+u[1]-u[2]-u[3]);
  dp[2]=  u[0]+u[1]-u[2]-u[3];
  dp[3]=  u[0]+u[1]-u[2]-u[3];

  return;
 }

IMFFlow IMFCreateADOFlow(MFErrorHandler e)
 {
  return IMFCreateFlow(4,1,ADOTorus,dADOTorus,dADOTorusdp,ddADOTorus,NULL,NULL,NULL,e);
 }

static void ddADOTorus(double *u, double *p, double *ddF,  void *data, MFErrorHandler e)
 {
  static char RoutineName[]={"dADOTorus"};
  static double a1=1.;
  static double a2=1.;
  static double b1=.55;
  static double b2=.55;

  ddF[0+4*(0+4*0)]= -6*u[0];
  ddF[0+4*(1+4*0)]= -2*u[1];
  ddF[0+4*(2+4*0)]=  0.;
  ddF[0+4*(3+4*0)]=  0.;

  ddF[1+4*(0+4*0)]=-2*u[1];
  ddF[1+4*(1+4*0)]=-2*u[0];
  ddF[1+4*(2+4*0)]= 0.;
  ddF[1+4*(3+4*0)]= 0.;

  ddF[2+4*(0+4*0)]= 0.;
  ddF[2+4*(1+4*0)]= 0.;
  ddF[2+4*(2+4*0)]= 0.;
  ddF[2+4*(3+4*0)]= 0.;

  ddF[3+4*(0+4*0)]= 0.;
  ddF[3+4*(1+4*0)]= 0.;
  ddF[3+4*(2+4*0)]= 0.;
  ddF[3+4*(3+4*0)]= 0.;

/* ----------------------------------- */

  ddF[0+4*(0+4*1)]= -2*u[1];
  ddF[0+4*(1+4*1)]= -2*u[0];
  ddF[0+4*(2+4*1)]=  0.;
  ddF[0+4*(3+4*1)]=  0.;

  ddF[1+4*(0+4*1)]=-2*u[0];
  ddF[1+4*(1+4*1)]=-6*u[1];
  ddF[1+4*(2+4*1)]=  0.;
  ddF[1+4*(3+4*1)]=  0.;

  ddF[2+4*(0+4*1)]=  0.;
  ddF[2+4*(1+4*1)]=  0.;
  ddF[2+4*(2+4*1)]=  0.;
  ddF[2+4*(3+4*1)]=  0.;

  ddF[3+4*(0+4*1)]=  0.;
  ddF[3+4*(1+4*1)]=  0.;
  ddF[3+4*(2+4*1)]=  0.;
  ddF[3+4*(3+4*1)]=  0.;

/* ----------------------------------- */

  ddF[0+4*(0+4*2)]=  0.;
  ddF[0+4*(1+4*2)]=  0.;
  ddF[0+4*(2+4*2)]=  0.;
  ddF[0+4*(3+4*2)]=  0.;

  ddF[1+4*(0+4*2)]=  0.;
  ddF[1+4*(1+4*2)]=  0.;
  ddF[1+4*(2+4*2)]=  0.;
  ddF[1+4*(3+4*2)]=  0.;

  ddF[2+4*(0+4*2)]=  0.;
  ddF[2+4*(1+4*2)]=  0.;
  ddF[2+4*(2+4*2)]= -6*u[2];
  ddF[2+4*(3+4*2)]= -2*u[3];

  ddF[3+4*(0+4*2)]=  0.;
  ddF[3+4*(1+4*2)]=  0.;
  ddF[3+4*(2+4*2)]=-2*u[3];
  ddF[3+4*(3+4*2)]=-2*u[2];

/* ----------------------------------- */

  ddF[0+4*(0+4*3)]=   0.;
  ddF[0+4*(1+4*3)]=   0.;
  ddF[0+4*(2+4*3)]=   0.;
  ddF[0+4*(3+4*3)]=   0.;

  ddF[1+4*(0+4*3)]=  0.;
  ddF[1+4*(1+4*3)]=  0.;
  ddF[1+4*(2+4*3)]=  0.;
  ddF[1+4*(3+4*3)]=  0.;

  ddF[2+4*(0+4*3)]=  0.;
  ddF[2+4*(1+4*3)]=  0.;
  ddF[2+4*(2+4*3)]= -2*u[3];
  ddF[2+4*(3+4*3)]= -2*u[2];

  ddF[3+4*(0+4*3)]=  0.;
  ddF[3+4*(1+4*3)]=  0.;
  ddF[3+4*(2+4*3)]= -2*u[2];
  ddF[3+4*(3+4*3)]= -6*u[3];

  return;
 }

static void ddADOTorusdp(double *u, double *p, double *dFdp,  void *data, MFErrorHandler e)
 {
  static char RoutineName[]={"dADOTorus"};
  static double a1=1.;
  static double a2=1.;
  static double b1=.55;
  static double b2=.55;

  dFdp[0+4*0]= -1.;
  dFdp[0+4*1]= -1.;
  dFdp[0+4*2]=  1.;
  dFdp[0+4*3]=  1.;

  dFdp[1+4*0]= -1.;
  dFdp[1+4*1]= -1.;
  dFdp[1+4*2]=  1.;
  dFdp[1+4*3]=  1.;

  dFdp[2+4*0]=  1.;
  dFdp[2+4*1]=  1.;
  dFdp[2+4*2]= -1.;
  dFdp[2+4*3]= -1.;

  dFdp[3+4*0]=  1.;
  dFdp[3+4*1]=  1.;
  dFdp[3+4*2]= -1.;
  dFdp[3+4*3]= -1.;

  return;
 }

static void ddADOTorusddp(double *u, double *p, double *ddp,  void *data, MFErrorHandler e)
 {
  static char RoutineName[]={"dADOTorusdp"};
  static double a1=1.;
  static double a2=1.;
  static double b1=.55;
  static double b2=.55;

  ddp[0]= 0.;
  ddp[1]= 0.;
  ddp[2]= 0.;
  ddp[3]= 0.;

  return;
 }

int ADOTorusProject(MFNVector u, double *x, void *d, MFErrorHandler e)
 {
  double r1,r2;
  double cost1,sint1;
  double cost2,sint2;
  double t1,t2;

  if(x==NULL)return 3;

  r1=sqrt(MFNV_C(u,0,e)*MFNV_C(u,0,e)+MFNV_C(u,1,e)*MFNV_C(u,1,e));
  r2=sqrt(MFNV_C(u,2,e)*MFNV_C(u,2,e)+MFNV_C(u,3,e)*MFNV_C(u,3,e));

  cost1=MFNV_C(u,0,e)/r1;
  sint1=MFNV_C(u,1,e)/r1;
  cost2=MFNV_C(u,2,e)/r2;
  sint2=MFNV_C(u,3,e)/r2;

  t1=atan2(sint1,cost1);
  t2=atan2(sint2,cost2);

  r1=2.;

  x[0]=(r1+r2*cost2)*cost1;
  x[1]=(r1+r2*cost2)*sint1;
  x[2]=r2*sint2;

  return 0;
 }

int ADOTorusProjectBB(MFNVector u, double *x, void *d, MFErrorHandler e)
 {
  int i,n;

  n=MFNV_NC(u,e);
  if(x==NULL)return n;

  for(i=0;i<n;i++)x[i]=MFNV_C(u,i,e);
  return 0;
 }

/*! \fn MFChart CIMFExtendAtlasAlongFatTraj(MFAtlas I, IMFFlow F, char *name, MFNVector u0, MFKVector p0, double dt, double ti, double tf,double epsilon,MFNRegion Omega, int maxSteps, int maxCharts,double Rmax,int nSkip,int skipCheck, MFErrorHandler e);
 *  \brief Integrates a fat trajectory and adds the points to an atlas.
 *
 * \param I The atlas to which new charts will be added.
 * \param F The flow.
 * \param name The name of the problem (used to page results to disk).
 * \param u0 The initial point on the fat trajectory.
 * \param p0 The flow parameters.
 * \param dt The initial time step.
 * \param ti The initial time.
 * \param tf The upper limit on time.
 * \param epsilon The maximum allowed size of the second order terms within a ball. Determines the stepsize.
 * \param Omega A region to limit the computation.
 * \param maxSteps The maximum number of steps to take along the fat trajectory.
 * \param maxCharts The maximum number of charts allowed in the atlas.
 * \param Rmax The largest radius of a ball.
 * \param nSkip The number of steps to make without adding points to the atlas.
 * \param skipCheck Used to skip the initial check for closeness to other points. Interpolated points can be close.
 * \param e An error handler.
 * \returns chart The last chart (in atlas I) on the fat trajectory.
 */
MFChart CIMFExtendAtlasAlongFatTraj(MFAtlas I, IMFFlow F, char *name, MFNVector u0, MFKVector p0, double dt, double ti, double tf,double epsilon,MFNRegion Omega, int maxSteps, int maxCharts,double Rmax,int nSkip, int skipCheck, MFErrorHandler e)
 {
  static char RoutineName[]={"CIMFExtendAtlasAlongFatTraj"};

  double t,tout;
  MFNKMatrix TS;
  double R;
  int i,j;
  int p,ni;
  double d;
  int nsteps;
  int neq,n;
  IMFExpansion E;
  MFNVector u;
  MFNVector f,g;
  double *fp;
  MFKVector s;
  double *y0;
  double *y;
  MFListOfCharts L;
  MFBinaryTree BTree;
  MFNSpace space;
  int k,prev;
  IMFFlow Fp;
  int nearChart;
  double dMin;
  double Rf;
  MFChart c;
  int bail;
  int verbose=0;
  double z[3];
  double *BBCenter=NULL;
  int BBDimension;
  MFChart thisChart;
  MFChart prevChart;
  MFChart firstChart;
  MFChart chart;
  MFNVector half;
  double tol1;
  double tol2;
  int nMin;
  double ratio;
  double rnear;

  int nTrajectory;

  firstChart=NULL;

  nTrajectory=0;

  if(MFAtlasNumberOfCharts(I,e)>0)
   {
    BBDimension=MFIMFProjectToBB(MFAtlasMF(I,e),MFAtlasChartCenter(I,0,e),NULL,e);
    BBCenter=(double*)malloc(BBDimension*sizeof(double));
    printf("   BBDimension=%d\n",BBDimension);fflush(stdout);
   }

#ifdef DOTRAJ
  if(IMFTraj==NULL)
   {
    IMFTraj=fopen("IMFTraj.dx","w");
    IMFNTraj=0;
   }
#endif

  E=IMFExpansionNVGetE(u0,e);
#ifdef MFALLOWVERBOSE
  if(verbose){printf("Initial expansion\n");IMFPrintExpansion(stdout,E,e);fflush(stdout);}
#endif
  neq=IMFExpansionDataLn(E,e);
  space=MFIMFNSpace(MFAtlasMF(I,e),e);
  BTree=MFAtlasGetBB(I,e);

  tol1=0.5; tol2=0.7;nMin=4;

  half=MFIMFVectorFactory(MFAtlasMF(I,e),e);

/* find closest point to keep trajectories from getting too close */

  R=Rmax;
  if(MFAtlasNumberOfCharts(I,e)>0 && !skipCheck)
   {
    MFIMFProjectToBB(MFAtlasMF(I,e),u0,BBCenter,e);
    L=MFCreateListOfNearbyCharts(BTree,BBCenter,R,e);
    ni=MFNumberOfIntersectingCharts(L,e);
    nearChart=-1;
    dMin=2*R;
    if(ni==0)dMin=0;
    for(p=0;p<ni;p++)
     {
      j=MFIntersectingChart(L,p,e);
      if(MFChartPaged(MFAtlasChart(I,j,e),e))continue;
      d=MFNSpaceDistance(space,MFAtlasCenterOfChart(I,j,e),u0,e);
      if(d<tol1*MFAtlasChartRadius(I,j,e)&&(nearChart==-1||d<dMin)) /* @@ */
       {
        dMin=d;
        nearChart=j;
        ratio=d/MFAtlasChartRadius(I,j,e);
       }
     }
    printf("   Check that initial point is away from others. dMin=%lf\n",dMin);fflush(stdout);
    if(nearChart==-1){printf("    point is OK, continue.\n");fflush(stdout);}
     else            {printf("    point is too close to chart %d (ratio is %lf, tol is %lf), go on to next.\n",j,ratio,tol1);fflush(stdout);}
    MFFreeListOfIntersectingCharts(L,e);
    if(nearChart!=-1){if(BBCenter!=NULL)free(BBCenter);MFFreeNVector(half,e);return NULL;}
   }else{
    if(skipCheck){printf("   Initial point is an interpolation point. Skip check if point is away from others.\n");fflush(stdout);}
      else       {printf("   There are no points yet (%d), skip check if point is away from others.\n",MFAtlasNumberOfCharts(I,e));fflush(stdout);}
    nearChart=-1;
   }
  if(nearChart!=-1)
   {
    if(BBCenter!=NULL)free(BBCenter);
    MFFreeNVector(half,e);
    return NULL;
   }

  t=ti;

  n=IMFFlowNU(F,e);
  k=IMFExpansionK(E,e);

  printf("CreateFatFlow\n");fflush(stdout);
  Fp=IMFCreateFatFlow(F,k,e);

  fp=(double*)malloc(n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(fp==NULL)
   {
    sprintf(IMFIntegrateFatErrorMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(double));
    MFSetError(e,12,RoutineName,IMFIntegrateFatErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    if(BBCenter!=NULL)free(BBCenter);
    MFFreeNVector(half,e);
    return NULL;
   }
#endif

  f=MFCreateWrappedNVector(n,fp,e);
  g=MFCreateNVector(n,e);
  s=MFCreateKVector(IMFExpansionK(E,e),e);

  nsteps=0;
  u=u0;
  MFRefNVector(u,e);
  prevChart=NULL;
  while(1)
   {
    bail=0;
    if(t>tf){printf("Ending integration, final time reached (%lf>%lf)\n",t,tf);fflush(stdout);bail=1;}

    if(nsteps>=maxSteps){printf("Ending integration, too many steps\n");fflush(stdout);goto FreeAndReturn;}

    if(MFAtlasNumberOfCharts(I,e)>=maxCharts){printf("Ending integration, too many charts\n");fflush(stdout);goto FreeAndReturn;}

    if(!MFNRegionInterior(Omega,u,e)){printf("Ending integration, point is outside Omega. ");MFPrintNVector(stdout,u,e);printf("\n");fflush(stdout);bail=1;}

    IMFEvaluateFlow(F,u,p0,fp,e);
    d=sqrt(MFNSpaceInner(space,f,f,e));
    if(d<epsilon){printf("Ending integration, |F| (%le<%le) too small (near fixed point)\n",d,epsilon);fflush(stdout);bail=1;}

    E=IMFExpansionNVGetE(u,e);
    TS=IMFExpansionTS(E,e);

    MFNSpaceScale(space,1./d,f,f,e);
    MFMVMulT(space,TS,f,s,e);
    MFMVMul(space,TS,s,g,e);
    d=MFNSpaceDistance(space,f,g,e);
    if(0&&d>epsilon){printf("Ending integration, F too far from tangent space, d=%l2\n",d);fflush(stdout);bail=1;}

    R=IMFExpansionR(E,epsilon,e);
    Rf=IMFFlowR(F,epsilon,u,p0,TS,e);
    if(Rf<R)R=Rf;
    if(R>Rmax||R!=R)R=Rmax;
    c=MFCreateChart(MFAtlasMF(I,e),u,TS,R,e);
    if(0&&R<epsilon){printf("Ending integration, R too small\n");fflush(stdout);bail=1;}

/* find closest point to keep trajectories from getting too close */

    nearChart=-1;
    if(MFAtlasNumberOfCharts(I,e)>0&&nTrajectory>nMin)
     {
      MFIMFProjectToBB(MFAtlasMF(I,e),MFChartCenter(c,e),BBCenter,e);
      L=MFCreateListOfNearbyCharts(BTree,BBCenter,R,e);
      ni=MFNumberOfIntersectingCharts(L,e);
      nearChart=-1;
      dMin=2*R;
      if(ni==0)dMin=0;
      for(p=0;p<ni;p++)
       {
        j=MFIntersectingChart(L,p,e);
        if(MFChartPaged(MFAtlasChart(I,j,e),e))continue;
        if(!IMFIsNear(I,MFAtlasChart(I,j,e),c,e))continue;
        d=MFNSpaceDistance(space,MFAtlasCenterOfChart(I,j,e),u,e);
        if(d<tol2*MFAtlasChartRadius(I,j,e)&&(nearChart==-1||d<dMin)) /* @@ */
         {
          dMin=d;
          ratio=d/MFAtlasChartRadius(I,j,e);
          nearChart=j;
          rnear=j;
         }
       }
/* (this+next)/2 */
      for(p=0;p<ni;p++)
       {
        j=MFIntersectingChart(L,p,e);
        if(MFChartPaged(MFAtlasChart(I,j,e),e))continue;
        if(!IMFIsNear(I,MFAtlasChart(I,j,e),c,e))continue;
        chart=MFAtlasChart(I,j,e);
        chart=MFChartGetNextChart(chart,e);
        if(chart!=NULL && !MFChartPaged(chart,e))
         {
          MFNSpaceAdd(space,MFChartCenter(chart,e),MFAtlasCenterOfChart(I,j,e),half,e);
          MFNSpaceScale(space,.5,half,half,e);
          d=MFNSpaceDistance(space,half,u,e);
          if(d<tol2*MFAtlasChartRadius(I,j,e)&&(nearChart==-1||d<dMin)) /* @@ */
           {
            dMin=d;
            nearChart=j;
            ratio=d/MFAtlasChartRadius(I,j,e);
            rnear=j+.5;
           }
         }
/* (this+last)/2 */
        chart=MFAtlasChart(I,j,e);
        chart=MFChartGetPrevChart(chart,e);
        if(chart!=NULL && !MFChartPaged(chart,e))
         {
          MFNSpaceAdd(space,MFChartCenter(chart,e),MFAtlasCenterOfChart(I,j,e),half,e);
          MFNSpaceScale(space,.5,half,half,e);
          d=MFNSpaceDistance(space,half,u,e);
          if(d<tol2*MFAtlasChartRadius(I,j,e)&&(nearChart==-1||d<dMin)) /* @@ */
           {
            dMin=d;
            nearChart=j;
            ratio=d/MFAtlasChartRadius(I,j,e);
            rnear=j-.5;
           }
         }
       }
/*    printf("   Check that this point is away from others.\n");fflush(stdout);
      if(nearChart==-1){printf("    point is OK, continue.\n");fflush(stdout);}
       else            {printf("    point is too close to chart %.1lf. ratio is %lf, tol is %lf\n",rnear,ratio,tol2);fflush(stdout);}
 */
      MFFreeListOfIntersectingCharts(L,e);
      MFFreeChart(c,e);
     }else{
/*
      if(!(nTrajectory<nMin)){printf("   Trajectory is too short, skip check if point is away from others.\n");fflush(stdout);}
        else       {printf("   There are no points yet (%d), skip check if point is away from others.\n",MFAtlasNumberOfCharts(I,e));fflush(stdout);}
 */
      nearChart=-1;
     }

/* Passed checks, add chart */

/* NV Types: 
           0     normal
           1     fixed point or on c
           2     interpolation point
           3     skipped point, not on manifold - don't plot
*/

    prev=-1;
    if(  ((nearChart==-1||IMFExpansionNVGetType(u,e)==1)||IMFExpansionNVGetType(u,e)==2) && nsteps>=nSkip )
     {
      printf("%d) Point %d on fat traj., R=%lf, t=%lf ",MFAtlasNumberOfCharts(I,e),nsteps,R,t);MFPrintNVector(stdout,u,e);fflush(stdout);
       {
        MFNVector col;

        col=MFMColumn(TS,0,e);
        printf("    phi_0  ");MFPrintNVector(stdout,col,e);printf("\n");fflush(stdout);
        MFFreeNVector(col,e);
        col=MFMColumn(TS,1,e);
        printf("    phi_1  ");MFPrintNVector(stdout,col,e);printf("\n");fflush(stdout);
        MFFreeNVector(col,e);
       }
      prev=MFAtlasAddChartWithAll(I,u,TS,R,e);
      if(prev>-1)
       {
        thisChart=MFAtlasChart(I,prev,e);
        MFChartSetPrevChart(thisChart,prevChart,e);
        if(prevChart!=NULL)MFChartSetNextChart(prevChart,thisChart,e);
         else              if(firstChart==NULL)firstChart=thisChart;
        prevChart=thisChart;
        nTrajectory++;
       }
      if(BTree==NULL)BTree=MFAtlasGetBB(I,e);
      if(BBCenter==NULL)
       {
        BBDimension=MFIMFProjectToBB(MFAtlasMF(I,e),MFAtlasChartCenter(I,0,e),NULL,e);
        BBCenter=(double*)malloc(BBDimension*sizeof(double));
       }
      if(0&&MFAtlasNumberOfCharts(I,e)%1000==0){printf("Paging\n");MFAtlasPageOutChartsNotNearBoundary(I,1,0,name,e);}
      if(prev>0 && (IMFExpansionNVGetType(u,e)==1||IMFExpansionNVGetType(u,e)==2))MFChartSetSingular(MFAtlasChart(I,prev,e),e);
       else if(prev>1)MFChartSetNonSingular(MFAtlasChart(I,prev,e),e);
#ifdef DOTRAJ
         {
          int m;
          m=MFIMFProjectToDraw(MFAtlasMF(I,e),NULL,NULL,e);
          MFIMFProjectToDraw(MFAtlasMF(I,e),u0,z,e);
          for(i=0;i<m;i++)fprintf(IMFTraj," %lf",z[i]);fprintf(IMFTraj,"\n");fflush(IMFTraj);IMFNTraj++;
          MFIMFProjectToDraw(MFAtlasMF(I,e),u,z,e);
          for(i=0;i<m;i++)fprintf(IMFTraj," %lf",z[i]);fprintf(IMFTraj,"\n");fflush(IMFTraj);IMFNTraj++;
         }
#endif

     }else{
      if(nearChart!=-1){printf("Ending integration, came too close to another trajectory\n");fflush(stdout);goto FreeAndReturn;}

      printf("*) Point %d on fat traj. ",nsteps);MFPrintNVector(stdout,u,e);printf(", R=%lf, t=%lf\n",R,t);fflush(stdout);
      if(nsteps<nSkip){printf("    is being skipped, %d<%d\n",nsteps,nSkip);fflush(stdout);}
      prev=MFAtlasAddChartToList(I,MFCreateChart(MFAtlasMF(I,e),u,TS,R,e),e);

      thisChart=MFAtlasChart(I,prev,e);
      MFChartSetPrevChart(thisChart,prevChart,e);
      if(prevChart!=NULL)MFChartSetNextChart(prevChart,thisChart,e);
      prevChart=thisChart;
  
      nTrajectory++;

#ifdef DOTRAJ
       {
        int m;
        m=MFIMFProjectToDraw(MFAtlasMF(I,e),NULL,NULL,e);
        MFIMFProjectToDraw(MFAtlasMF(I,e),u0,z,e);
        for(i=0;i<m;i++)fprintf(IMFTraj," %lf",z[i]);fprintf(IMFTraj,"\n");fflush(IMFTraj);IMFNTraj++;
        MFIMFProjectToDraw(MFAtlasMF(I,e),u,z,e);
        for(i=0;i<m;i++)fprintf(IMFTraj," %lf",z[i]);fprintf(IMFTraj,"\n");fflush(IMFTraj);IMFNTraj++;
       }
#endif
      MFChartSetSingular(MFAtlasChart(I,prev,e),e);
      IMFExpansionNVSetType(u,3,e);

     }

    if(bail){if(prev>-1)MFChartSetSingular(MFAtlasChart(I,prev,e),e);goto FreeAndReturn;}

    nsteps++;

/* Now find the next */

    y0=IMFExpansionData(E,e);
    y=(double*)malloc(IMFExpansionDataLn(E,e)*sizeof(double));

#ifndef MFNOSAFETYNET
    if(y==NULL)
     {
      sprintf(IMFIntegrateFatErrorMsg,"Out of memory, trying to allocate %d bytes",IMFExpansionDataLn(E,e)*sizeof(double));
      MFSetError(e,12,RoutineName,IMFIntegrateFatErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      if(BBCenter!=NULL)free(BBCenter);
      MFFreeNVector(half,e);
      return NULL;
     }
#endif

    for(i=0;i<IMFExpansionDataLn(E,e);i++)y[i]=y0[i];
#ifdef DOTRAJ
    if(0&&nsteps>nSkip)
     {
      int m;
      MFNVector u;
      m=MFIMFProjectToDraw(MFAtlasMF(I,e),NULL,NULL,e);
      u=MFCreateWrappedNVector(IMFExpansionDataLn(E,e),y,e);
      MFIMFProjectToDraw(MFAtlasMF(I,e),u,z,e);
      for(i=0;i<m;i++)fprintf(IMFTraj," %lf",z[i]);fprintf(IMFTraj,"\n");fflush(IMFTraj);IMFNTraj++;
      MFFreeNVector(u,e);
     }
#endif

    d=0;
    tout=t+R;
    while(t+1.e-7<=tout && t<tf && d<1.03*R)
     {
      if(!MEHKellerInt(Fp,y,p0,t,tout,2,e)){
#ifdef DOTRAJ
  if(0){
                int m;
              MFNVector u;
          m=MFIMFProjectToDraw(MFAtlasMF(I,e),NULL,NULL,e);
          u=MFCreateWrappedNVector(IMFExpansionDataLn(E,e),y,e);
          MFIMFProjectToDraw(MFAtlasMF(I,e),u,z,e);
          for(i=0;i<m;i++)fprintf(IMFTraj," %lf",z[i]);fprintf(IMFTraj,"\n");fflush(IMFTraj);IMFNTraj++;
          MFFreeNVector(u,e);
         }
#endif
        printf("Ending integration, Integrator failed to converge\n");fflush(stdout);goto FreeAndReturn;}
/* was ,5); */
    
      for(i=0;i<n;i++)
       {
        if(y[i]!=y[i]){printf("Ending integration, Integrator returned a NaN\n");fflush(stdout);goto FreeAndReturn;}
       }
      d=0.;for(i=0;i<n;i++)d+=pow(y[i]-y0[i],2);d=sqrt(d);
      t=tout;
      tout=t+.1*R;
     }
#ifdef DOTRAJ
     if(0){
      int m;
      MFNVector u;
     m=MFIMFProjectToDraw(MFAtlasMF(I,e),NULL,NULL,e);
    u=MFCreateWrappedNVector(IMFExpansionDataLn(E,e),y,e);
          MFIMFProjectToDraw(MFAtlasMF(I,e),u,z,e);
          for(i=0;i<m;i++)fprintf(IMFTraj," %lf",z[i]);fprintf(IMFTraj,"\n");fflush(IMFTraj);IMFNTraj++;
          MFFreeNVector(u,e);
         }
#endif
    E=IMFCreateExpansion(n,k,e);
    IMFExpansionSetDerivatives(E,y,y+n,y+n*k,NULL,e);
    free(y);

    u0=u;
    u=IMFCreateExpansionNVector(E,t,IMFExpansionNVGetSigma(u0,e),prev,0,e);
    IMFExpansionNVSetS0(u,IMFExpansionNVGetS0(u0,e),e);
    IMFExpansionNVSetChart0(u,IMFExpansionNVGetChart0(u0,e),e);
    MFFreeNVector(u0,e);
    IMFFreeExpansion(E,e);
   }

FreeAndReturn:
  MFFreeNVector(u,e);
  MFFreeNVector(f,e);
  free(fp);
  IMFFreeFlow(Fp,e);
  MFFreeNVector(g,e);
  MFFreeKVector(s,e);

/*printf("done %s\n",RoutineName);fflush(stdout);*/
  if(BBCenter!=NULL)free(BBCenter);
  MFFreeNVector(half,e);
  return firstChart;
 }

/*! \fn MFAtlas CIMFComputeInvariantManifold(IMFFlow F,MFKVector p0,char *name, MFAtlas c,MFNRegion Omega, double eps, double dt, double tmax, int maxInterp, int maxCharts, double Rmax, MFErrorHandler e);
 * \brief Computes an atlas of charts that cover the image of a manifold under a flow. A streamsurface.
 *
 * \param L The flow.
 * \param name A name to used for paging files and so forth.
 * \param u0 The fixed point.
 * \param p0 The parameters of the flow for the fixed point.
 * \param c The manifold of initial conditions.
 * \param Omega A region to bound the computation.
 * \param eps The tolerance on the size of the quadratic terms over a balls. Controls the stepsize.
 * \param dt The initial timestep to use along a fat trajectory.
 * \param tmax The upper limit on trajectory length.
 * \param maxInterp The upper limit on the number of interpolation performed.
 * \param maxCharts The upper limit on the number of charts in the atlas.
 * \param RMax An upper limit to impose on the radius of the balls along the fat trajectories.
 * \param e An MFErrorHandler to handle exceptions and errors.
 * \returns A new Atlas
 */
MFAtlas CIMFComputeInvariantManifold(IMFFlow F,MFKVector p0,char *name, MFAtlas c,MFNRegion Omega, double eps, double dt, double tmax, int maxInterp, int maxCharts, double Rmax,MFErrorHandler e)
 {
  static char RoutineName[]={"CIMFComputeInvariantManifold"};
  double s0[2],s1[2];
  double *ddu;
  double *du;
  MFKVector s;
  MFNVector ug,v;
  MFNKMatrix Phi;
  MFNKMatrix Psi;
  int i;
  int ii,jj;
  IMFExpansion cE,U;
  MFContinuationMethod H;
  MFAtlas A;
  MFAtlas I;
  MFImplicitMF M;
  double R,r;
  char sname[1024];
  FILE *fid;
  MFNKMatrix TS;
  MFNVector ui;
  double t;
  int chart;
  MFNVector sigma;
  MFKVector zero;
  double *u0;
  int n,k;
  int verbose=0;
  MFChart *Trajectory=NULL;
  int nTrajectories;
  int mTrajectories=10000;

  if(verbose){printf("in %s\n",RoutineName);fflush(stdout);}

  Trajectory=(MFChart*)realloc(Trajectory,mTrajectories*sizeof(MFChart));

  n=MFAtlasN(c,e);
  k=MFAtlasK(c,e);
  M=IMFCreateFlat(n,k+1,e);
  MFIMFSetProjectForDraw(M,MFIMFGetProjectForDraw(MFAtlasMF(c,e),e),e);
  MFIMFSetProjectForSave(M,MFIMFGetProjectForSave(MFAtlasMF(c,e),e),e);
  MFIMFSetProjectForBB  (M,MFIMFGetProjectForBB  (MFAtlasMF(c,e),e),e);
  MFIMFSetR(M,Rmax/2,e);
  I=MFCreateAtlas(M,e);
  printf("%s, Create Atlas I=0x%8.8x\n",RoutineName,I);fflush(stdout);

  H=MFCreateMultifariosMethod(e);
  MFMultifarioSetIntegerParameter(H,"page",0,e);
  MFMultifarioSetIntegerParameter(H,"maxCharts",-1,e);
  MFMultifarioSetIntegerParameter(H,"dumpToPlotFile",1,e);
  MFMultifarioSetFilename(H,name,e);

  u0=(double*)malloc(n*sizeof(double));

  IMFBC=fopen("IMFBC.dx","w");
  IMFNBC=0;
  nTrajectories=0;

#ifndef MFNOSAFETYNET
  if(u0==NULL)
   {
    sprintf(MFNSpaceErrorMsg,"Out of memory trying to allocate %d bytes\n",n*k*sizeof(double));
    MFSetError(e,12,RoutineName,MFNSpaceErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  du=(double*)malloc(n*k*sizeof(double));

#ifndef MFNOSAFETYNET
  if(du==NULL)
   {
    sprintf(MFNSpaceErrorMsg,"Out of memory trying to allocate %d bytes\n",n*k*sizeof(double));
    MFSetError(e,12,RoutineName,MFNSpaceErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  ddu=(double*)malloc(n*k*k*sizeof(double));

#ifndef MFNOSAFETYNET
  if(ddu==NULL)
   {
    sprintf(MFNSpaceErrorMsg,"Out of memory trying to allocate %d bytes\n",n*sizeof(double));
    MFSetError(e,12,RoutineName,MFNSpaceErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  zero=MFCreateKVector(k-1,e);

  s0[0]=-r; s1[0]= r;
  s0[1]=-r; s1[1]= r;
  for(i=0;i<n*k*k;i++)ddu[i]=0.;
  for(i=0;i<MFAtlasNumberOfCharts(c,e);i++)
   {
    if(verbose){printf("   chart %d on the initial manifold\n",i);fflush(stdout);}
    if(MFNV_CStar(MFAtlasChartCenter(c,i,e),e)!=NULL)
     {
      for(ii=0;ii<n;ii++)
        u0[ii]=MFNV_CStar(MFAtlasChartCenter(c,i,e),e)[ii];
     }else{
      for(ii=0;ii<n;ii++)
        u0[ii]=MFNV_C(MFAtlasChartCenter(c,i,e),ii,e);
     }
    Phi=MFAtlasChartTangentSpace(c,i,e);

    if(MFNKM_CStar(Phi,e)!=NULL)
     {
      for(ii=0;ii<n*k;ii++)
        du[ii]=MFNKM_CStar(Phi,e)[ii];
     }else{
      for(ii=0;ii<n;ii++)
       for(jj=0;jj<k;jj++)
        du[ii+n*jj]=MFNKMGetC(Phi,ii,jj,e);
     }

    cE=IMFCreateExpansion(n,k,e);

/* Assumes it's a straight line  */
/* Assumes Phi is a dense matrix */


    IMFExpansionSetDerivatives(cE,u0,du,ddu,NULL,e);
    U=IMFInflateExpansionWithFlow(cE,F,p0,e);

    ui=IMFCreateExpansionNVector(U,0.,MFAtlasCenterOfChart(c,i,e),-1,1,e);
    IMFExpansionNVSetChart0(ui,i,e);
    IMFExpansionNVSetS0(ui,zero,e);
    printf("*) Point on c ");MFPrintNVector(stdout,ui,e);fflush(stdout);
    printf("   Extend along fat traj\n");fflush(stdout);
    if(!(nTrajectories<mTrajectories))
     {
      mTrajectories+=10000;
      Trajectory=(MFChart*)realloc(Trajectory,mTrajectories*sizeof(MFChart));
     }
    Trajectory[nTrajectories]=CIMFExtendAtlasAlongFatTraj(I,F,name,ui,p0,dt,0.,tmax,eps,Omega,maxCharts,maxCharts,Rmax,0,0,e);
    if(Trajectory[nTrajectories]!=NULL)nTrajectories++;
    printf("   done fat traj\n\n");fflush(stdout);


    if(0&&MFAtlasNumberOfCharts(I,e)%1000==0)MFAtlasPageOutChartsNotNearBoundary(I,1,0,name,e);
    MFFreeNVector(ui,e);

    IMFFreeExpansion(cE,e);
    IMFFreeExpansion(U,e);
   }

  printf("*** Done covering c\n");fflush(stdout);
#ifdef DOINTERPOLATIONPT
  sprintf(sname,"%s.intpt",name);
  IMFInterp=fopen(sname,"w");
  IMFNInterp=0;
#endif

  i=0;
  while(i<maxInterp&& (ui=IMFGetInterpolationPoint(I,F,p0,A,tmax,NULL,e))!=NULL && MFAtlasNumberOfCharts(I,e)<maxCharts)
   {
    printf("\n*** Interpolated Point %d",i);fflush(stdout);
    MFPrintNVector(stdout,ui,e);
    printf("\n");fflush(stdout);

#ifdef DOINTERPOLATIONPT
    for(j=0;j<3;j++)fprintf(IMFInterp," %lf",IMFExpansionU(IMFExpansionNVGetE(ui,e))[j],e);
    fprintf(IMFInterp,"\n");fflush(IMFInterp);
    IMFNInterp++;
#endif
    if(nTrajectories>=mTrajectories)
     {
      mTrajectories+=10000;
      Trajectory=(MFChart*)realloc(Trajectory,mTrajectories*sizeof(MFChart));
     }
    Trajectory[nTrajectories]=CIMFExtendAtlasAlongFatTraj(I,F,name,ui,p0,dt,0.,tmax,eps,Omega,maxCharts,maxCharts,Rmax,0,1,e);
    if(Trajectory[nTrajectories]!=NULL)nTrajectories++;

    if(0&&MFAtlasNumberOfCharts(I,e)%1000==0)MFAtlasPageOutChartsNotNearBoundary(I,1,0,name,e);
    MFFreeNVector(ui,e);
    i++;
   }

  if(i<maxInterp)
    {printf("*** No more interpolation points, total was %d\n",i);fflush(stdout);}
   else
    {printf("*** Max number of interpolation points reached, %d\n",i);fflush(stdout);}

  MFFlushAtlas(H,I,e);
  MFCloseAtlas(H,I,e);

  MFFreeContinuationMethod(H,e);
  MFFreeImplicitMF(M,e);

  free(u0);
  free(du);
  free(ddu);

  {
   MFChart chart;
   MFChart nextChart;
   double tf;
   double t;

   for(i=0;i<nTrajectories;i++)
    {
     chart=Trajectory[i];
     while((nextChart=MFChartGetNextChart(chart,e))!=NULL)chart=nextChart;
     tf=IMFExpansionNVGetT(MFChartCenter(chart,e),e);
     chart=Trajectory[i];
     while((nextChart=MFChartGetNextChart(chart,e))!=NULL)
      {
       t=IMFExpansionNVGetT(MFChartCenter(chart,e),e);
       IMFExpansionNVSetT(MFChartCenter(chart,e),t/tf,e);
       chart=nextChart;
      }
     t=IMFExpansionNVGetT(MFChartCenter(chart,e),e);
     IMFExpansionNVSetT(MFChartCenter(chart,e),t/tf,e);
    }
  }

  for(i=0;i<nTrajectories;i++)
   {
    printf("Trajectory %d",i);
    CIMFGetBoundaryConditions(I,Trajectory[i],nTrajectories,Trajectory,e);
   }

#ifdef DOINTERPOLATIONPT
  if(IMFInterp!=NULL)
   {
    fclose(IMFInterp);

    sprintf(sname,"i%s.dx",name);
    IMFInterp=fopen(sname,"w");
    sprintf(sname,"%s.intpt",name);
    if(IMFNInterp>0)
     {
      fprintf(IMFInterp,"object \"Vertices\" class array type float rank 1 shape 3 items %d data file %s.dx\n",sname,IMFNInterp);
     }else{
      fprintf(IMFInterp,"object \"Vertices\" class array type float rank 1 shape 3 items 1 data follows\n");
      fprintf(IMFInterp,"  0. 0. 0.\n");
     }
    fprintf(IMFInterper,"object \"interpolated points\" class field\n");
    fprintf(IMFInterper,"component \"positions\" value \"Vertices\"\n");
    fclose(IMFInterper);
   }
#endif

#ifdef DOINTERPANIM

  if(IMFInterper!=NULL)
   {
    fclose(IMFInterper);
    IMFInterper=fopen("gInterp1.dx","w");
    if(IMFNInterper>0)
     {
      fprintf(IMFInterper,"object \"Vertices\" class array type float rank 1 shape 3 items %d data file Interp1.dx\n",IMFNInterper);
     }else{
      fprintf(IMFInterper,"object \"Vertices\" class array type float rank 1 shape 3 items 1 data follows\n");
      fprintf(IMFInterper,"  0. 0. 0.\n");
     }
    fprintf(IMFInterper,"object \"interpolated points\" class field\n");
    fprintf(IMFInterper,"component \"positions\" value \"Vertices\"\n");
    fclose(IMFInterper);
   }

  if(IMFInterpee!=NULL)
   {
    fclose(IMFInterpee);
    IMFInterpee=fopen("gInterp2.dx","w");
    if(IMFNInterpee>0)
     {
      fprintf(IMFInterpee,"object \"Vertices\" class array type float rank 1 shape 3 items %d data file Interp2.dx\n",IMFNInterpee);
      fprintf(IMFInterpee,"object \"Lines\" class array type int rank 1 shape 2 items %d data follows\n",IMFNInterpee/2);
      for(i=0;i<IMFNInterpee;i+=2)
      fprintf(IMFInterpee,"   %d %d\n",i,i+1);
     }else{
      fprintf(IMFInterpee,"object \"Vertices\" class array type float rank 1 shape 3 items 2 data follows\n");
      fprintf(IMFInterpee,"  0. 0. 0.\n");
      fprintf(IMFInterpee,"  0. 0. 0.001\n");
      fprintf(IMFInterpee,"object \"Lines\" class array type int rank 1 shape 2 items 1 data follows\n");
      fprintf(IMFInterpee,"   %d %d\n",0,1);
     }
    fprintf(IMFInterpee,"attribute \"ref\" string \"positions\"\n");
    fprintf(IMFInterpee,"attribute \"element type\" string \"lines\"\n");

    fprintf(IMFInterpee,"object \"interpolation points\" class field\n");
    fprintf(IMFInterpee,"component \"positions\" value \"Vertices\"\n");
    fprintf(IMFInterpee,"component \"connections\" value \"Lines\"\n");
    fclose(IMFInterpee);
   }

  if(IMFInterpT!=NULL)
   {
    fclose(IMFInterpT);
    IMFInterpT=fopen("gInterp3.dx","w");
    if(IMFNInterpT>0)
     {
      fprintf(IMFInterpT,"object \"Vertices\" class array type float rank 1 shape 3 items %d data file Interp3.dx\n",IMFNInterpT
);
      fprintf(IMFInterpT,"object \"Lines\" class array type int rank 1 shape 2 items %d data follows\n",IMFNInterpT-1);
      for(i=0;i<IMFNInterpT-1;i++)
      fprintf(IMFInterpT,"   %d %d\n",i,i+1);
     }else{
      fprintf(IMFInterpT,"object \"Vertices\" class array type float rank 1 shape 3 items 2 data follows\n");
      fprintf(IMFInterpT,"  0. 0. 0.\n");
      fprintf(IMFInterpT,"  0. 0. 0.001\n");
      fprintf(IMFInterpT,"object \"Lines\" class array type int rank 1 shape 2 items 1 data follows\n");
      fprintf(IMFInterpT,"   %d %d\n",0,1);
     }
    fprintf(IMFInterpT,"attribute \"ref\" string \"positions\"\n");
    fprintf(IMFInterpT,"attribute \"element type\" string \"lines\"\n");

    fprintf(IMFInterpT,"object \"interpolation points\" class field\n");
    fprintf(IMFInterpT,"component \"positions\" value \"Vertices\"\n");
    fprintf(IMFInterpT,"component \"connections\" value \"Lines\"\n");
    fclose(IMFInterpT);
   }

  if(IMFCircle!=NULL)
   {
    fclose(IMFCircle);
    IMFCircle=fopen("gIMFCircle.dx","w");
    fprintf(IMFCircle,"object \"Vertices\" class array type float rank 1 shape 3 items %d data file IMFCircle.dx\n",IMFNCircle);
    fprintf(IMFCircle,"object \"Lines\" class array type int rank 1 shape 2 items %d data follows\n",IMFNCircle/2);
    for(i=0;i<IMFNCircle;i+=2)
    fprintf(IMFCircle,"   %d %d\n",i,i+1);
    fprintf(IMFCircle,"attribute \"ref\" string \"positions\"\n");
    fprintf(IMFCircle,"attribute \"element type\" string \"lines\"\n");

    fprintf(IMFCircle,"object \"interpolation points\" class field\n");
    fprintf(IMFCircle,"component \"positions\" value \"Vertices\"\n");
    fprintf(IMFCircle,"component \"connections\" value \"Lines\"\n");
    fclose(IMFCircle);
   }
#endif

#ifdef DOTRAJ
  if(IMFTraj!=NULL)
   {
    printf("dump trajectory file\n");fflush(stdout);
    fclose(IMFTraj);
    IMFTraj=fopen("gIMFTraj.dx","w");
    fprintf(IMFTraj,"object \"Vertices\" class array type float rank 1 shape 3 items %d data file IMFTraj.dx\n",IMFNTraj);
    fprintf(IMFTraj,"object \"Lines\" class array type int rank 1 shape 2 items %d data follows\n",IMFNTraj/2);
    for(i=0;i<IMFNTraj;i+=2)
    fprintf(IMFTraj,"   %d %d\n",i,i+1);
    fprintf(IMFTraj,"attribute \"ref\" string \"positions\"\n");
    fprintf(IMFTraj,"attribute \"element type\" string \"lines\"\n");

    fprintf(IMFTraj,"object \"trajectories\" class field\n");
    fprintf(IMFTraj,"component \"positions\" value \"Vertices\"\n");
    fprintf(IMFTraj,"component \"connections\" value \"Lines\"\n");
    fclose(IMFTraj);
   }
#endif

  if(IMFBC!=NULL)
   {
    fclose(IMFBC);
    IMFBC=fopen("gIMFBC.dx","w");
    fprintf(IMFBC,"object \"Vertices\" class array type float rank 1 shape 3 items %d data file IMFBC.dx\n",IMFNBC);
    fprintf(IMFBC,"object \"Lines\" class array type int rank 1 shape 2 items %d data follows\n",IMFNBC/2);
    for(i=0;i<IMFNBC;i+=2)
    fprintf(IMFBC,"   %d %d\n",i,i+1);
    fprintf(IMFBC,"attribute \"ref\" string \"positions\"\n");
    fprintf(IMFBC,"attribute \"element type\" string \"lines\"\n");

    fprintf(IMFBC,"object \"trajectories\" class field\n");
    fprintf(IMFBC,"component \"positions\" value \"Vertices\"\n");
    fprintf(IMFBC,"component \"connections\" value \"Lines\"\n");
    fclose(IMFBC);
   }

  printf("done %s\n",RoutineName);fflush(stdout);


  return I;
 }

static int MEHKellerInt(IMFFlow F,double *y,MFKVector p0,double ti,double tf,int nsteps, MFErrorHandler e)
 {
  static char RoutineName[]={"MEHKellerInt"};
  double t,dt;
  static int meq=-1;
  static double *y0=NULL;
  static double *ym=NULL;
  static double *y1=NULL;
  static double *f=NULL;
  static double *A=NULL;
  static double *b=NULL;
  double fNorm;
  int i,j,itimes;
  double error;
  int verbose=0;
  int neq;
  MFNVector vy;
  MFNVector vym;
  MFNVector vy0;

  neq=IMFFlowNU(F,e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s, neq=%d integrate from t=%lf to %lf\n",RoutineName,neq,ti,tf);fflush(stdout);}
#endif

  if(meq<neq)
   {
    y0=(double*)realloc((void*)y0,neq*sizeof(double));

#ifndef MFNOSAFETYNET
    if(y0==NULL)
     {
      sprintf(IMFIntegrateFatErrorMsg,"Out of memory, trying to allocate %d bytes",neq*sizeof(double));
      MFSetError(e,12,RoutineName,IMFIntegrateFatErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return 0;
     }
#endif

    ym=(double*)realloc((void*)ym,neq*sizeof(double));

#ifndef MFNOSAFETYNET
    if(ym==NULL)
     {
      sprintf(IMFIntegrateFatErrorMsg,"Out of memory, trying to allocate %d bytes",neq*sizeof(double));
      MFSetError(e,12,RoutineName,IMFIntegrateFatErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return 0;
     }
#endif

    y1=(double*)realloc((void*)y1,neq*sizeof(double));

#ifndef MFNOSAFETYNET
    if(y1==NULL)
     {
      sprintf(IMFIntegrateFatErrorMsg,"Out of memory, trying to allocate %d bytes",neq*sizeof(double));
      MFSetError(e,12,RoutineName,IMFIntegrateFatErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return 0;
     }
#endif

    f =(double*)realloc((void*)f ,neq*sizeof(double));

#ifndef MFNOSAFETYNET
    if(f==NULL)
     {
      sprintf(IMFIntegrateFatErrorMsg,"Out of memory, trying to allocate %d bytes",neq*sizeof(double));
      MFSetError(e,12,RoutineName,IMFIntegrateFatErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return 0;
     }
#endif

    A =(double*)realloc((void*)A ,neq*neq*sizeof(double));

#ifndef MFNOSAFETYNET
    if(A==NULL)
     {
      sprintf(IMFIntegrateFatErrorMsg,"Out of memory, trying to allocate %d bytes",neq*neq*sizeof(double));
      MFSetError(e,12,RoutineName,IMFIntegrateFatErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return 0;
     }
#endif

    b =(double*)realloc((void*)b ,neq*sizeof(double));

#ifndef MFNOSAFETYNET
    if(b==NULL)
     {
      sprintf(IMFIntegrateFatErrorMsg,"Out of memory, trying to allocate %d bytes",neq*sizeof(double));
      MFSetError(e,12,RoutineName,IMFIntegrateFatErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return 0;
     }
#endif

    meq=neq;
   }

  vy=MFCreateWrappedNVector(neq,y,e);
  vym=MFCreateWrappedNVector(neq,ym,e);
  vy0=MFCreateWrappedNVector(neq,y0,e);

  IMFEvaluateFlow(F,vy,p0,f,e);
  fNorm=0.;for(i=0;i<neq;i++)fNorm+=f[i]*f[i]; fNorm=sqrt(fNorm);
  if(fNorm>10.)fNorm=10.;


  nsteps=round(10*(tf-ti)*fNorm+2);

  t=ti;
  dt=(tf-ti)/nsteps;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("   nsteps=%d, fNorm=%lf\n",nsteps,fNorm);fflush(stdout);}
#endif

  for(j=0;j<neq;j++)y0[j]=y[j];
  for(i=0;i<nsteps;i++)
   {

#ifdef MFALLOWVERBOSE
    if(verbose){printf("Keller Box for IVP, step %d/%d, t=%lf\n",i,nsteps,t);fflush(stdout);  }
#endif

    IMFEvaluateFlow(F,vy0,p0,f,e);

    for(j=0;j<neq;j++)y1[j]=y0[j]+dt*f[j];

    error=1.;
    itimes=0;
    while(error>1.e-5)
     {
      for(j=0;j<neq;j++)ym[j]=.5*(y0[j]+y1[j]);
      IMFEvaluateFlow(F,vym,p0,f,e);
      IMFEvaluateDerivativeOfFlow(F,vym,p0,A,e);
      for(j=0;j<neq*neq;j++)A[j]=-.5*dt*A[j];
      error=0.;
      for(j=0;j<neq;j++)
       {
        b[j]=-(y1[j]-y0[j]-dt*f[j]);
        error+=b[j]*b[j];
        A[j+neq*j]=A[j+neq*j]+1.;
       }

#ifdef MFALLOWVERBOSE
      if(verbose){printf("Iteration %d, error=%le\n",itimes,error);fflush(stdout);  }
#endif

      MFSolveFull(neq,A,b,e);
      for(j=0;j<neq;j++)y1[j]+=b[j];
      itimes++;
      if(itimes>100)return 0;
     }
    t+=dt;
    for(j=0;j<neq;j++)y0[j]=y1[j];
   }
  for(j=0;j<neq;j++)y[j]=y1[j];

  MFFreeNVector(vy,e);
  MFFreeNVector(vym,e);
  MFFreeNVector(vy0,e);
  return 1;
 }

int CIMFIsNear(MFAtlas A,MFChart c0,MFChart c1,MFErrorHandler e)
 {
  static char RoutineName[]={"IsNear"};
  MFNVector u0,u1;
  double t0,t1;

  u0=MFChartCenter(c0,e);
  u1=MFChartCenter(c1,e);

/*t0=IMFExpansionNVGetT(u0,e);
  t1=IMFExpansionNVGetT(u1,e);
  printf("IMFIsNeart, nCharts=%d, t0=%lf, t1=%lf\n",MFAtlasNumberOfCharts(A,e),t0,t1);fflush(stdout);*/
/*  if(fabs(t0-t1)>30. && t0>0. && t1>0. )return 0;   Not for CompactInvariant */

  return 1;
 }

void CIMFGetBoundaryConditions(MFAtlas I,MFChart chart,int nT, MFChart *T, MFErrorHandler e)
 {
  static char RoutineName[]={"CIMFGetBoundaryConditions"};
  int ni;
  int i,j,l;
  int m; 
  double z0[1000];
  double z [1000];
  MFChart nextChart;
  MFNVector u0;
  double R;
  MFPolytope P;
  int fi;
  int traj;
  int ii;
  double time;
  int jj;

  R=MFChartRadius(chart,e);
  u0=MFChartCenter(chart,e);
  R=0.04;

  traj=-1;for(ii=0;ii<nT;ii++)if(chart==T[ii])traj=ii;
  jj=MFChartGetPositionInAtlas(chart,e);
  time=IMFExpansionNVGetT(MFChartCenter(chart,e),e);
  printf("	u_{%d}(%lf)=%d",traj,time,jj);fflush(stdout);

  while((nextChart=MFChartGetNextChart(chart,e))!=NULL)chart=nextChart;

  jj=MFChartGetPositionInAtlas(chart,e);
  time=IMFExpansionNVGetT(MFChartCenter(chart,e),e);
  printf("	u_{%d}(%lf)=%d",traj,time,jj);fflush(stdout);

  u0=MFChartCenter(chart,e);
  m=MFIMFProjectToDraw(MFAtlasMF(I,e),NULL,NULL,e);
    MFIMFProjectToDraw(MFAtlasMF(I,e),u0,z0,e);

  P=MFChartPolytope(chart,e);
  ni=MFPolytopeNumberOfFaces(P,e);
  for(i=0;i<ni;i++)
   {
    fi=MFPolytopeFaceIndex(P,i,e);
    if(fi>3)
     {
      j=MFAtlasHalfSpaceLeftChart(I,fi,e);
      if(MFAtlasChart(I,j,e)==chart)j=MFAtlasHalfSpaceRightChart(I,fi,e);
  
      if(j>-1)
       {
        traj=-1;
        time=IMFExpansionNVGetT(MFAtlasChartCenter(I,j,e),e);
        for(ii=0;ii<nT;ii++)
         {
          nextChart=T[ii];
          if(MFChartGetPositionInAtlas(nextChart,e)==j)traj=ii;
          while((nextChart=MFChartGetNextChart(nextChart,e))!=NULL)
            if(MFChartGetPositionInAtlas(nextChart,e)==j)traj=ii;
         }

        printf("	u_{%d}(%lf)=%d",traj,IMFExpansionNVGetT(MFAtlasChartCenter(I,j,e),e),j);fflush(stdout);

        for(l=0;l<m;l++)fprintf(IMFBC," %lf",z0[l]);fprintf(IMFBC,"\n");fflush(IMFBC);IMFNBC++;
        MFIMFProjectToDraw(MFAtlasMF(I,e),MFAtlasChartCenter(I,j,e),z,e);
        for(l=0;l<m;l++)fprintf(IMFBC," %lf",z[l]);fprintf(IMFBC,"\n");fflush(IMFBC);IMFNBC++;
       }
     }
   }
  printf("\n");fflush(stdout);

  return;
 }

#ifdef __cplusplus
 }
#endif
