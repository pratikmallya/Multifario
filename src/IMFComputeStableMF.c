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

static char *id="@(#) $Id: IMFComputeStableMF.c,v 1.7 2011/07/21 17:42:46 mhender Exp $";

static char MFNSpaceErrorMsg[256]="";

#include <IMF.h>
#include <IMFExpansion.h>
#include <IMFFlow.h>
#include <IMFFixedPt.h>
#include <math.h>
#include <IMFSphereOnExpansion.h>
#include <IMFExpansionSpace.h>
#include <IMFExpansionPt.h>
#include <IMFIntegrateFat.h>
#include <IMFInterpolation.h>
#include <IMFFlat.h>
#include <stdio.h>
#include <MFAtlas.h>
#include <MFMultifariosMethod.h>
#include <MFErrorHandler.h>
#include <MFPrint.h>

#ifdef __cplusplus
 extern "C" {
#endif

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

/*
#define DOINTERPANIM
*/
#define DOTRAJ

void MFAtlasPageOutChartsNotNearBoundary(MFAtlas,int,int,char*,MFErrorHandler);
void MFAtlasSetNearRtn(MFAtlas,int (*)(MFAtlas,MFChart,MFChart,MFErrorHandler),MFErrorHandler);
int IMFIsNear(MFAtlas,MFChart,MFChart,MFErrorHandler);
MFChart MFAtlasChart(MFAtlas,int,MFErrorHandler);

MFNVector IMFGetInterpolationPointOnList(MFAtlas,IMFFlow,MFKVector,MFAtlas,double,MFNRegion,int,int*,MFErrorHandler);

/*! \fn MFAtlas IMFComputeStableInvariantManifold(IMFFlow L,char *name, MFNVector u0, MFKVector p0, MFNRegion Omega, double eps, double dt, double tmax, int maxInterp, int maxCharts, double R0,MFErrorHandler e);
 * \brief Computes an atlas of charts that cover the stable manifold of a hyerbolic fixed point.
 *
 * \param L The flow.
 * \param name A name to used for paging files and so forth.
 * \param u0 The fixed point.
 * \param p0 The parameters of the flow for the fixed point.
 * \param Omega A region to bound the computation.
 * \param eps The tolerance on the size of the quadratic terms over a balls. Controls the stepsize.
 * \param dt The initial timestep to use along a fat trajectory.
 * \param tmax The upper limit on trajectory length.
 * \param maxInterp The upper limit on the number of interpolation performed.
 * \param maxCharts The upper limit on the number of charts in the atlas.
 * \param R0 The radius for the initial ball about the fixedpoint that serves as initial surface for the manifold. 
 * \param e An MFErrorHandler to handle exceptions and errors.
 * \returns A new Atlas
 */
MFAtlas IMFComputeStableInvariantManifold(IMFFlow L,char *name, MFNVector u0, MFKVector p0, MFNRegion Omega, double eps, double dt, double tmax, int maxInterp, int maxCharts, double R0, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFComputeStableInvariantManifold"};
  IMFFlow F;
  MFAtlas A;

  F=IMFCreateBackwardFlow(L,e);
  A=IMFComputeUnstableInvariantManifold(F,name,u0,p0,Omega,eps,dt,tmax,maxInterp,maxCharts,R0,e);
  IMFFreeFlow(F,e);
  return A;
 }

/*! \fn MFAtlas IMFComputeUnstableInvariantManifold(IMFFlow F,char *name, MFNVector u0, MFKVector p0, MFNRegion Omega, double eps, double dt, double tmax, int maxInterp, int maxCharts, double R0, MFErrorHandler e);
 * \brief Computes an atlas of charts that cover the unstable manifold of a hyerbolic fixed point.
 *
 * \param L The flow.
 * \param name A name to used for paging files and so forth.
 * \param u0 The fixed point.
 * \param p0 The parameters of the flow for the fixed point.
 * \param Omega A region to bound the computation.
 * \param eps The tolerance on the size of the quadratic terms over a balls. Controls the stepsize.
 * \param dt The initial timestep to use along a fat trajectory.
 * \param tmax The upper limit on trajectory length.
 * \param maxInterp The upper limit on the number of interpolation performed.
 * \param maxCharts The upper limit on the number of charts in the atlas.
 * \param R0 The radius for the initial ball about the fixedpoint that serves as initial surface for the manifold. 
 * \param e An MFErrorHandler to handle exceptions and errors.
 * \returns A new Atlas
 */
MFAtlas IMFComputeUnstableInvariantManifold(IMFFlow F,char *name, MFNVector u0, MFKVector p0, MFNRegion Omega, double eps, double dt, double tmax, int maxInterp, int maxCharts, double R0, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFComputeUnstableInvariantManifold"};
  IMFExpansion u,a;
  double s0[2],s1[2];
  MFKVector s;
  MFNVector ug,v;
  MFNVector ustar;
  MFNKMatrix Phi0;
  MFNKMatrix Phi;
  MFNKMatrix Psi;
  MFImplicitMF c;
  int i,j;
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
  MFNVector ut;
  double t;
  int chart;
  MFNVector sigma;
  MFKVector zero;
  int n,k;
  int verbose=0;

  Phi=IMFGetBasisForUnstableInvariantSubspace(F,u0,p0,e);
  Psi=IMFGetBasisForOrthogonalComplement(Phi,e);

  n=MFNKMatrixN(Phi,e);
  k=MFNKMatrixK(Phi,e);
  u=IMFCreateExpansion(n,k,e);
  a=IMFCreateExpansion(k,k,e);
  IMFFindExpansionNearFixedPt(u0,p0,Phi,Psi,F,u,a,e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("done IMFFindExpansionNearFixedPt\n");fflush(stdout);printf("Expansion at fixed point:\n");IMFPrintExpansion(stdout,u,e);fflush(stdout);}
#endif

  zero=MFCreateKVector(k-1,e);

  R=IMFExpansionR(u,eps,e);if(R!=R||R>R0)R=R0;
  r=.5*R;
/* For U1!    R=R0; r=.05; */
  c=IMFCreateSphereOnExpansion(u,F,p0,eps,R,r,e);

  s =MFCreateKVector(2,e);
  MFKVSetC(s,0,R,e);
  MFKVSetC(s,1,0.,e);

  ug=MFCreateNVector(n,e);
  IMFEvaluateExpansion(u,s,ug,e);

  ut=MFCreateNVector(n+k,e);
  for(i=0;i<n;i++)MFNVSetC(ut,i,MFNV_C(ug,i,e),e);
  for(i=0;i<k;i++)MFNVSetC(ut,n+i,MFKV_C(s,i,e),e);
  MFFreeNVector(ug,e);

  ustar=MFCreateNVector(n+k,e);

  printf("Cover the initial surface\n");
  Phi0=MFIMFTangentSpace(c,ut,e);
  if(MFIMFProject(c,ut,Phi0,ustar,e))
   {

    H=MFCreateMultifariosMethod(e);
    MFMultifarioSetRealParameter(H,"epsilon",.8,e);
    MFMultifarioSetIntegerParameter(H,"maxCharts",-1,e);
    MFMultifarioSetIntegerParameter(H,"verbose",1,e);
    MFMultifarioSetIntegerParameter(H,"page",0,e);
    MFMultifarioSetIntegerParameter(H,"dumpToPlotFile",0,e);
    MFMultifarioSetIntegerParameter(H,"dumpToCenterFile",0,e);
    MFMultifarioSetFilename(H,"SphereOnU0",e);

    A=MFComputeAtlas(H,c,Omega,ustar,e);
    MFFreeNVector(ut,e);
    MFCloseAtlas(H,A,e);
    MFFreeContinuationMethod(H,e);
    printf("Done computing Atlas for initial values\n");fflush(stdout);
   }else{
    printf("Projection of initial point failed\n");fflush(stdout);
    return NULL;
   }

  M=IMFCreateFlat(n,k,e);
  I=MFCreateAtlas(M,e);
  MFAtlasSetNearRtn(I,IMFIsNear,e);

  H=MFCreateMultifariosMethod(e);
  MFMultifarioSetIntegerParameter(H,"page",0,e);
  MFMultifarioSetIntegerParameter(H,"dumpToPlotFile",1,e);
  MFMultifarioSetFilename(H,name,e);

  TS=IMFExpansionTS(u,e);
  ui=IMFCreateExpansionNVector(u,-100.,u0,-1,1,e);
  MFAtlasAddChartWithAll(I,ui,TS,R,e);
  printf("0) Fixed Point. ");MFPrintNVector(stdout,ui,e);printf(", R=%lf\n",R);fflush(stdout);
  MFFreeNVector(ui,e);
  MFFreeNKMatrix(TS,e);

  s0[0]=-r; s1[0]= r;
  s0[1]=-r; s1[1]= r;
  for(i=0;i<MFAtlasNumberOfCharts(A,e);i++)
   {
    cE=IMFSphereOnExpansionGetLocal(A,MFAtlasCenterOfChart(A,i,e),e);
    U=IMFInflateExpansionWithFlow(cE,F,p0,e);

    ui=IMFCreateExpansionNVector(U,0.,MFAtlasCenterOfChart(A,i,e),-1,1,e);
    IMFExpansionNVSetChart0(ui,i,e);
    IMFExpansionNVSetS0(ui,zero,e);
    printf("*) Point on c");fflush(stdout);
    IMFExtendAtlasAlongFatTraj(I,F,name,ui,p0,dt,0.,tmax,eps,Omega,maxCharts,maxCharts,1.,0,e);
    if(0&&MFAtlasNumberOfCharts(I,e)%1000==0){printf("Paging\n");MFAtlasPageOutChartsNotNearBoundary(I,1,0,name,e);}
    MFFreeNVector(ui,e);

    IMFFreeExpansion(cE,e);
    IMFFreeExpansion(U,e);
    if(MFAtlasNumberOfCharts(I,e)>maxCharts)
     {
      printf("Maximum number of charts exceeded (%d>%d).\n",MFAtlasNumberOfCharts(I,e),maxCharts,e);fflush(stdout);
      goto DoneCover;
     }
   }

  printf("*** Done covering c\n");fflush(stdout);
#ifdef DOINTERPOLATIONPT
  sprintf(sname,"%s.intpt",name);
  IMFInterp=fopen(sname,"w");
  IMFNInterp=0;
#endif

  i=0;
  printf("*** Begin interpolation\n");fflush(stdout);
  while(i<maxInterp&& (ui=IMFGetInterpolationPoint(I,F,p0,A,tmax,NULL,e))!=NULL && MFAtlasNumberOfCharts(I,e)<maxCharts)
   {
    printf("*** Interpolated Point %d\n",i);fflush(stdout);
#ifdef DOINTERPOLATIONPT
    for(j=0;j<3;j++)fprintf(IMFInterp," %lf",IMFExpansionU(IMFExpansionNVGetE(ui,e))[j],e);
    fprintf(IMFInterp,"\n");fflush(IMFInterp);
    IMFNInterp++;
#endif
    if(i<maxInterp-1)IMFExtendAtlasAlongFatTraj(I,F,name,ui,p0,dt,IMFExpansionNVGetT(ui,e),tmax,eps,Omega,maxCharts,maxCharts,1.,0,e);
    if(0&&MFAtlasNumberOfCharts(I,e)%1000==0){printf("Paging\n");MFAtlasPageOutChartsNotNearBoundary(I,1,0,name,e);}
    MFFreeNVector(ui,e);
    i++;
    if(MFAtlasNumberOfCharts(I,e)>maxCharts)
     {
      printf("Maximum number of charts exceeded (%d>%d).\n",MFAtlasNumberOfCharts(I,e),maxCharts);fflush(stdout);
      goto DoneCover;
     }
   }

  if(i<maxInterp)
    {printf("*** No more interpolation points, total was %d\n",i);fflush(stdout);}
   else
    {printf("*** Max number of interpolation points reached, %d\n",i);fflush(stdout);}

DoneCover:
  printf("*** Done cover\n");fflush(stdout);

  MFFreeImplicitMF(c,e);
  MFFreeKVector(s,e);

  MFCloseAtlas(H,I,e);
  MFFreeContinuationMethod(H,e);

  MFFreeNKMatrix(Phi0,e);

  MFFreeNVector(ustar,e);

  MFFreeNKMatrix(Phi,e);

  MFFreeNKMatrix(Psi,e);

  IMFFreeExpansion(u,e);

  IMFFreeExpansion(a,e);

  MFFreeImplicitMF(M,e);

  MFFreeAtlas(A,e);

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

  printf("done %s\n",RoutineName);fflush(stdout);

  return I;
 }

/*! \fn MFAtlas IMFComputeInvariantManifold(IMFFlow F,MFKVector p0,char *name, MFAtlas c,MFNRegion Omega, double eps, double dt, double tmax, int maxInterp, int maxCharts, double Rmax, MFErrorHandler e);
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
MFAtlas IMFComputeInvariantManifold(IMFFlow F,MFKVector p0,char *name, MFAtlas c,MFNRegion Omega, double eps, double dt, double tmax, int maxInterp, int maxCharts, double Rmax,MFErrorHandler e)
 {
  static char RoutineName[]={"IMFComputeInvariantManifold"};
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

  if(verbose){printf("in %s\n",RoutineName);fflush(stdout);}

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
    if(0&&i==0) /* Why??? */
     {
      TS=IMFExpansionTS(cE,e);
      R=MFAtlasChartRadius(c,i,e);
      MFAtlasAddChartWithAll(I,ui,TS,R,e);
     }
    IMFExpansionNVSetChart0(ui,i,e);
    IMFExpansionNVSetS0(ui,zero,e);
    printf("*) Point on c\n");fflush(stdout);
    MFPrintNVector(stdout,ui,e);
    printf("   Extend along fat traj\n");fflush(stdout);
    IMFExtendAtlasAlongFatTraj(I,F,name,ui,p0,dt,0.,tmax,eps,Omega,maxCharts,maxCharts,Rmax,0,e);
    printf("   done fat traj\n");fflush(stdout);
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
    printf("*** Interpolated Point %d\n",i);fflush(stdout);
#ifdef DOINTERPOLATIONPT
    for(j=0;j<3;j++)fprintf(IMFInterp," %lf",IMFExpansionU(IMFExpansionNVGetE(ui,e))[j],e);
    fprintf(IMFInterp,"\n");fflush(IMFInterp);
    IMFNInterp++;
#endif
    IMFExtendAtlasAlongFatTraj(I,F,name,ui,p0,dt,IMFExpansionNVGetT(ui,e),tmax,eps,Omega,maxCharts,maxCharts,Rmax,0,e);
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

  printf("done %s\n",RoutineName);fflush(stdout);


  return I;
 }

int IMFIsNear(MFAtlas A,MFChart c0,MFChart c1,MFErrorHandler e)
 {
  static char RoutineName[]={"IsNear"};
  MFNVector u0,u1;
  double t0,t1;

  u0=MFChartCenter(c0,e);
  u1=MFChartCenter(c1,e);

  t0=IMFExpansionNVGetT(u0,e);
  t1=IMFExpansionNVGetT(u1,e);
/*printf("IMFIsNeart, nCharts=%d, t0=%lf, t1=%lf\n",MFAtlasNumberOfCharts(A,e),t0,t1);fflush(stdout);*/
  if(fabs(t0-t1)>30. && t0>0. && t1>0. )return 0;

  return 1;
 }

MFAtlas IMFComputeStableInvariantManifold2(IMFFlow L,char *name, MFNVector u0,MFKVector p0,MFNRegion Omega, double eps, double dt, double tmax, int maxInterp, int maxCharts, double R0, MFErrorHandler e)
 {
  static char RoutineName[]={"ComputeStableInvariantManifold2"};
  IMFFlow F;
  MFAtlas A;

  F=IMFCreateBackwardFlow(L,e);
  A=IMFComputeUnstableInvariantManifold2(F,name,u0,p0,Omega,eps,dt,tmax,maxInterp,maxCharts,R0,e);
  IMFFreeFlow(F,e);
  return A;
 }

MFAtlas IMFComputeUnstableInvariantManifold2(IMFFlow F,char *name, MFNVector u0,MFKVector p0,MFNRegion Omega, double eps, double dt, double tmax, int maxInterp, int maxCharts, double R0,MFErrorHandler e)
 {
  static char RoutineName[]={"ComputeUnstableInvariantManifold2"};
  IMFExpansion u,a;
  double s0[2],s1[2];
  MFKVector s;
  MFNVector ug,v;
  MFNVector ustar;
  MFNKMatrix Phi0;
  MFNKMatrix Phi;
  MFNKMatrix Psi;
  MFImplicitMF c;
  int i,j,l;
  IMFExpansion cE,U;
  MFContinuationMethod H;
  MFAtlas A;
  MFAtlas I;
  MFImplicitMF M;
  double R,r;
  double Rf,Rmax;
  char sname[1024];
  FILE *fid;
  MFNKMatrix TS;
  MFNVector ui;
  MFNVector ut;
  double t;
  int chart;
  MFNVector sigma;
  int n,k;
  int nSkip;
  int iring,maxRings;
  double epsilon=1.;

  int nCharts;
  int *chartList=NULL;
  int mCharts;
  MFChart *cList=NULL;
  MFChart nchart;
  double dist,d;
  MFKVector zero;

  Rmax=1.;

  Phi=IMFGetBasisForUnstableInvariantSubspace(F,u0,p0,e);
  Psi=IMFGetBasisForOrthogonalComplement(Phi,e);

  n=MFNKMatrixN(Phi,e);
  k=MFNKMatrixK(Phi,e);
  printf("n=%d, k=%d\n",n,k);fflush(stdout);
  u=IMFCreateExpansion(n,k,e);
  a=IMFCreateExpansion(k,k,e);
  IMFFindExpansionNearFixedPt(u0,p0,Phi,Psi,F,u,a,e);
  zero=MFCreateKVector(k-1,e);

  R=IMFExpansionR(u,eps,e);if(R!=R||R>R0)R=R0;
  r=.5*R;
  c=IMFCreateSphereOnExpansion(u,F,p0,eps,R,r,e);

  s =MFCreateKVector(2,e);
  MFKVSetC(s,0,R,e);
  MFKVSetC(s,1,0.,e);

  ug=MFCreateNVector(n,e);
  IMFEvaluateExpansion(u,s,ug,e);

  ut=MFCreateNVector(n+k,e);
  for(i=0;i<n;i++)MFNVSetC(ut,i,MFNV_C(ug,i,e),e);
  for(i=0;i<k;i++)MFNVSetC(ut,n+i,MFKV_C(s,i,e),e);
  MFFreeNVector(ug,e);

  ustar=MFCreateNVector(n+k,e);

  Phi0=MFIMFTangentSpace(c,ut,e);
  if(MFIMFProject(c,ut,Phi0,ustar,e))
   {

    H=MFCreateMultifariosMethod(e);
    MFMultifarioSetRealParameter(H,"epsilon",.8,e);
    MFMultifarioSetIntegerParameter(H,"maxCharts",-1,e);
    MFMultifarioSetIntegerParameter(H,"verbose",0,e);
    MFMultifarioSetIntegerParameter(H,"page",0,e);
    MFMultifarioSetIntegerParameter(H,"dumpToPlotFile",0,e);
    MFMultifarioSetIntegerParameter(H,"dumpToCenterFile",0,e);
    MFMultifarioSetFilename(H,"SphereOnU0",e);

    A=MFComputeAtlas(H,c,Omega,ustar,e);
    MFFreeNVector(ut,e);
    MFCloseAtlas(H,A,e);
    MFFreeContinuationMethod(H,e);
    printf("Done computing Atlas\n");fflush(stdout);
   }else{
    printf("Projection of initial point failed\n");fflush(stdout);
    return NULL;
   }

  M=IMFCreateFlat(n,k,e);
  I=MFCreateAtlas(M,e);
/*MFAtlasSetNearRtn(I,IMFIsNear,e);*/

  H=MFCreateMultifariosMethod(e);
  MFMultifarioSetIntegerParameter(H,"page",0,e);
  MFMultifarioSetIntegerParameter(H,"dumpToPlotFile",1,e);
  MFMultifarioSetFilename(H,name,e);

  TS=IMFExpansionTS(u,e);
  ui=IMFCreateExpansionNVector(u,-100.,u0,-1,1,e);
  chart=MFAtlasAddChartWithAll(I,ui,TS,R,e);
  printf("%d) Fixed Point. ",chart);MFPrintNVector(stdout,ui,e);printf(", R=%lf\n",R);fflush(stdout);
  MFFreeNVector(ui,e);
  MFFreeNKMatrix(TS,e);

/* --------------------------------------------------------------------- */

  s0[0]=-r; s1[0]= r;
  s0[1]=-r; s1[1]= r;
  for(i=0;i<MFAtlasNumberOfCharts(A,e);i++)
   {
    cE=IMFSphereOnExpansionGetLocal(A,MFAtlasCenterOfChart(A,i,e),e);
    U=IMFInflateExpansionWithFlow(cE,F,p0,e);

    ui=IMFCreateExpansionNVector(U,0.,MFAtlasCenterOfChart(A,i,e),-1,1,e);
    IMFExpansionNVSetChart0(ui,i,e);
    IMFExpansionNVSetS0(ui,zero,e);

    TS=IMFExpansionTS(U,e);
    R=IMFExpansionR(U,.5*epsilon,e);
    Rf=IMFFlowR(F,.5*epsilon,ui,p0,TS,e);
    if(Rf<R)R=Rf;
    if(R>Rmax||R!=R)R=Rmax;
    chart=MFAtlasAddChartWithAll(I,ui,TS,R,e);
    printf("%d) Point on c ",chart);MFPrintNVector(stdout,ui,e);printf(", R=%lf\n",R);fflush(stdout);

    MFFreeNKMatrix(TS,e);
    MFFreeNVector(ui,e);
    IMFFreeExpansion(cE,e);
    IMFFreeExpansion(U,e);
    if(MFAtlasNumberOfCharts(I,e)>maxCharts)
     {
      printf("Maximum number of charts exceeded (%d>%d).\n",MFAtlasNumberOfCharts(I,e),maxCharts);fflush(stdout);
      goto DoneCover;
     }
   }

  printf("*** Done covering c\n");fflush(stdout);

/* --------------------------------------------------------------------- */

  i=0;
  iring=0;
  maxRings=1;
  while(iring<maxRings&&i<maxInterp&&MFAtlasNumberOfCharts(I,e)<maxCharts&&MFAtlasNumberOfChartsWithBoundary(I,e)>0)
   {
    nCharts=MFAtlasNumberOfChartsWithBoundary(I,e);
    chartList=(int*)realloc((void*)chartList,nCharts*sizeof(int));

#ifndef MFNOSAFETYNET
    if(chartList==NULL)
     {
      sprintf(MFNSpaceErrorMsg,"Out of memory trying to allocate %d bytes\n",nCharts*sizeof(int));
      MFSetError(e,12,RoutineName,MFNSpaceErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return NULL;
     }
#endif

    mCharts=0;
    cList=(MFChart*)realloc((void*)cList,nCharts*sizeof(MFChart));

#ifndef MFNOSAFETYNET
    if(cList==NULL)
     {
      sprintf(MFNSpaceErrorMsg,"Out of memory trying to allocate %d bytes\n",nCharts*sizeof(MFChart));
      MFSetError(e,12,RoutineName,MFNSpaceErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return NULL;
     }
#endif

    for(j=0;j<nCharts;j++)
     chartList[j]=MFAtlasChartWithBoundary(I,j,e);

    printf("Extend ring %d\n",iring);fflush(stdout);
    for(j=0;j<nCharts;j++)
     {
      ui=MFAtlasChartCenter(I,chartList[j],e);
      nSkip=1;
      printf("-------Front point %d/%d\n",j,nCharts);fflush(stdout);
      nchart=IMFStepAlongFatTraj(I,F,ui,p0,MFAtlasChartRadius(I,chartList[j],e),dt,eps,Omega,Rmax,e);

      if(nchart!=NULL)
       {
        cList[mCharts]=nchart;
        mCharts++;
       }

      printf("-------\n");fflush(stdout);
     }
    iring++;

    for(j=0;j<mCharts;j++)
     {
      for(l=0;l<j;l++)
       {
        d=MFNSpaceDistance(MFIMFNSpace(MFAtlasMF(I,e),e),MFChartCenter(cList[j],e),MFChartCenter(cList[l],e),e)/MFChartRadius(cList[j],e);
        if(l==0||d<dist)dist=d;
       }
      if(1||j==0||dist>1.)
        MFAtlasAddChartWithAll(I,MFChartCenter(cList[j],e),MFChartTangentSpace(cList[j],e),MFChartRadius(cList[j],e),e);
     }

    for(j=0;j<mCharts;j++)MFFreeChart(cList[j],e);

    for(j=0;j<nCharts;j++)
     {
      if(nCharts>1&&!MFChartHasBoundary(MFAtlasChart(I,chartList[j],e),e) )
       {
        if(j!=nCharts-1)chartList[j]=chartList[nCharts-1];
        nCharts--;
        if(j>0)j--;
       }
     }

    while(0&&i<maxInterp&& (ui=IMFGetInterpolationPointOnList(I,F,p0,A,tmax,NULL,nCharts,chartList,e))!=NULL && MFAtlasNumberOfCharts(I,e)<maxCharts)
     {
      printf("*** Interpolated Point %d\n",i);fflush(stdout);

      U=IMFExpansionNVGetE(ui,e);
      TS=IMFExpansionTS(U,e);
      R=IMFExpansionR(U,epsilon,e);
      Rf=IMFFlowR(F,epsilon,ui,p0,TS,e);
      if(Rf<R)R=Rf;
      if(R>Rmax||R!=R)R=Rmax;
      MFAtlasAddChartWithAll(I,ui,TS,R,e);

      MFFreeNKMatrix(TS,e);
      IMFFreeExpansion(cE,e);
      IMFFreeExpansion(U,e);
      MFFreeNVector(ui,e);
      if(MFAtlasNumberOfCharts(I,e)>maxCharts)
       {
        printf("Maximum number of charts exceeded (%d>%d).\n",MFAtlasNumberOfCharts(I,e),maxCharts);fflush(stdout);
        goto DoneCover;
       }
      i++;

      for(j=0;j<nCharts;j++)
       {
        if(nCharts>1&&!MFChartHasBoundary(MFAtlasChart(I,chartList[j],e),e) )
         {
          if(j!=nCharts-1)chartList[j]=chartList[nCharts-1];
          nCharts--;
          if(j>0)j--;
         }
       }
     }
   }

  if(i<maxInterp)
    {printf("*** No more boundary, total was %d interpolated, %d charts\n",i,MFAtlasNumberOfCharts(I,e));fflush(stdout);}
   else
    {printf("*** Max number of interpolation points reached, %d\n",i);fflush(stdout);}

/* --------------------------------------------------------------------- */

DoneCover:

  MFFreeImplicitMF(c,e);
  MFFreeKVector(s,e);

  printf("Flushing invariant manifold\n");
  MFCloseAtlas(H,I,e);
  MFFreeContinuationMethod(H,e);

  MFFreeNKMatrix(Phi0,e);

  MFFreeNVector(ustar,e);
  MFFreeNKMatrix(Phi,e);
  MFFreeNKMatrix(Psi,e);

  IMFFreeExpansion(u,e);
  IMFFreeExpansion(a,e);
  MFFreeImplicitMF(M,e);
  MFFreeAtlas(A,e);
  free(chartList);
  free(cList);
  free(zero);

  return I;
 }

#ifdef __cplusplus
}
#endif
