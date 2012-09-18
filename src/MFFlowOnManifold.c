/*
 *  PROGRAM NAME:  multifario
 *
 *  (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
 *  CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
 *
 *  Please refer to the LICENSE file in the top directory
 *
 *      author: Mike Henderson mhender@watson.ibm.com
 *      date:   May 20, 2010       Adapted from IMFComputeStableMF
 */

static char *id="@(#) $Id: MFFlowOnManifold.c,v 1.2 2011/07/21 17:42:46 mhender Exp $";

static char MFFlowOnManifoldErrorMsg[256]="";

#include <multifarioConfig.h>
#include <IMF.h>
#include <IMFFlow.h>
#include <math.h>
#include <stdio.h>
#include <MFAtlas.h>
#include <MFChart.h>
#include <MFErrorHandler.h>
#include <MFPrint.h>

#ifdef __cplusplus
 extern "C" {
#endif

FILE *IMFTraj=NULL;
int   IMFNTraj=0;
FILE *IMFTrajSingle=NULL;
int   IMFNTrajSingle=0;

FILE *IMFBVP=NULL;
FILE *IMFENDS=NULL;
FILE *IMFENDPTS=NULL;
#define WRITEBVP
static int fpos=0;

/*void MFAtlasPageOutChartsNotNearBoundary(MFAtlas,int,int,char*,MFErrorHandler);*/
void MFAtlasWriteOutChartsNotPaged(MFAtlas,int,int,char*,MFErrorHandler);

void MFAtlasSetNearRtn(MFAtlas,int (*)(MFAtlas,MFChart,MFChart,MFErrorHandler),MFErrorHandler);
MFChart MFAtlasChart(MFAtlas,int,MFErrorHandler);
int MFPolytopeVertexNumberOfIndices(MFPolytope,int,MFErrorHandler);

static void MFNullVectorSmall(int,double*,double*,MFErrorHandler);
static void MFSolveFullSmall(int,double*,double*,MFErrorHandler);
double *MFKV_CStar(MFKVector,MFErrorHandler);
double *MFNV_CStar(MFNVector,MFErrorHandler);
int MFChartChangedFlag(MFChart,MFErrorHandler);
void MFAtlasRemoveChartFromBoundaryList(MFAtlas,int,MFErrorHandler);
void MFChartResetChangedFlag(MFChart,MFErrorHandler);
int MFPolytopeIntersectIndexSets(MFPolytope,int,int,int*,MFErrorHandler);
static void MFSolveFullSmall(int,double*,double*,MFErrorHandler);
int  MFSolveFull(int,double*,double*,MFErrorHandler);

MFBinaryTree MFAtlasGetBB(MFAtlas,MFErrorHandler);
MFChart MFAtlasChart(MFAtlas,int,MFErrorHandler);

static int MEHKellerInt(IMFFlow,double*,MFKVector,double,double,int,MFErrorHandler);
int MFIsNear(MFAtlas,MFChart,MFChart,MFErrorHandler);
void MFChartSetPaged(MFChart,MFErrorHandler);

#define round(x) ((x)>=0?(long)((x)+0.5):(long)((x)-0.5))
int MFAtlasHalfSpaceLeftChart(MFAtlas,int,MFErrorHandler);
int MFAtlasHalfSpaceRightChart(MFAtlas,int,MFErrorHandler);

void MFNullVectorSmall(int k,double *A,double *phi, MFErrorHandler e);
void MFSolveFullSmall(int k,double *A,double *b, MFErrorHandler e);

MFNVector MFGetInterpolationPoint(MFAtlas A,IMFFlow F, MFKVector p0, MFNRegion Omega, MFChart *, MFErrorHandler e);
void MFExtendAtlasAlongTrajOnManifold(MFAtlas A, IMFFlow F, MFChart chart, char *name, MFNVector u0, MFKVector p0, double dt, double ti, double tf,double epsilon,MFNRegion Omega, int maxSteps, int maxCharts,double Rmax,double Dmin, int nSkip, MFErrorHandler e);
MFNVector MFInterpolateOnDualFace(MFAtlas A, IMFFlow F, MFKVector p0, MFNVector z, int chart, int *inter, MFErrorHandler e);

#define kmax 10

/*! \fn MFAtlas MFCoverFlowOnManifold(MFImplicitMF M, IMFFlow L,char *name, int nInitial, MFNVector *u0, MFKVector p0, MFNRegion Omega, double eps, double dt, double tmax, double dmin, int maxInterp, int maxCharts, MFErrorHandler e);
 * \brief Computes an atlas of charts that cover the unstable manifold of a hyerbolic fixed point.
 *
 * \param M The manifold.
 * \param L The flow.
 * \param name A name to used for paging files and so forth.
 * \param nInitial The number of initial points.
 * \param u0 An array of initial points on M.
 * \param p0 The parameters of the flow for the fixed point.
 * \param Omega A region to bound the computation.
 * \param eps The tolerance on the size of the quadratic terms over a balls. Controls the stepsize.
 * \param dt The initial timestep to use along a fat trajectory.
 * \param tmax The upper limit on trajectory length.
 * \param dmin The fractional minimum distance between trajs.
 * \param maxInterp The upper limit on the number of interpolation performed.
 * \param maxCharts The upper limit on the number of charts in the atlas.
 * \param e An MFErrorHandler to handle exceptions and errors.
 * \returns A new Atlas
 */
MFAtlas MFCoverFlowOnManifold(MFImplicitMF M, IMFFlow F,char *name, int nInitial, MFNVector *u0, MFKVector p0, MFNRegion Omega, double eps, double dt, double tmax, double dmin, int maxInterp, int maxCharts, MFErrorHandler e)
 {
  static char RoutineName[]={"MFCoverFlowOnManifold"};
  int i,j;
  MFAtlas A;
  MFNVector ui;
  int n,k;
  int verbose=0;
  MFChart chart;
  char filename[1024];

  printf("*** %s\n",RoutineName);fflush(stdout);

  if(IMFTraj==NULL)
   {
    sprintf(filename,"%sTraj.dx",name);
    IMFTraj=fopen(filename,"w");
    IMFNTraj=0;
   }

#if 1
  if(IMFTrajSingle==NULL)
   {
    sprintf(filename,"%sTrajSingle.dx",name);
    IMFTrajSingle=fopen(filename,"w");
    IMFNTrajSingle=0;
   }
#endif

#ifdef WRITEBVP
  if(IMFBVP==NULL)
   {
    sprintf(filename,"%sBVP.dx",name);
    IMFBVP=fopen(filename,"w");
   }
#endif

  if(IMFENDPTS==NULL)
   {
    sprintf(filename,"%sENDPTS.dx",name);
    IMFENDPTS=fopen(filename,"w");
   }

  A=MFCreateAtlas(M,e);

  for(i=0;i<nInitial;i++)
    MFExtendAtlasAlongTrajOnManifold(A,F,(MFChart)NULL,name,u0[i],p0,dt,0.,tmax,eps,Omega,maxCharts,maxCharts,1.,dmin,0,e);

  printf("*** Done trajectories from initial points\n");fflush(stdout);

  i=0;
  printf("*** Begin interpolation\n");fflush(stdout);
  while(i<maxInterp&& (ui=MFGetInterpolationPoint(A,F,p0,Omega,&chart,e))!=NULL && (maxCharts<0 || MFAtlasNumberOfCharts(A,e)<maxCharts) )
   {
    printf("*** Interpolated Point %d, ui (0x%8.8x), index=%d\n",i,ui,MFNVGetIndex(ui,e));fflush(stdout);

    if(i<maxInterp-1)
     {
      MFExtendAtlasAlongTrajOnManifold(A,F,chart,name,ui,p0,dt,0.,tmax,eps,Omega,maxCharts,maxCharts,1.,dmin,0,e);
     }

    MFFreeNVector(ui,e);

    i++;
    if(maxCharts>0 && MFAtlasNumberOfCharts(A,e)>maxCharts)
     {
      printf("Maximum number of charts exceeded (%d>%d).\n",MFAtlasNumberOfCharts(A,e),maxCharts);fflush(stdout);
      goto DoneCover;
     }
   }

  if(i<maxInterp)
    {printf("*** No more interpolation points, total was %d\n",i);fflush(stdout);}
   else
    {printf("*** Max number of interpolation points reached, %d\n",i);fflush(stdout);}

DoneCover:
  printf("*** Done cover\n");fflush(stdout);
  MFAtlasWriteOutChartsNotPaged(A,1,0,name,e);

  printf("done %s\n",RoutineName);fflush(stdout);

  if(IMFTraj!=NULL)
   {
    printf("dump trajectory file\n");fflush(stdout);
    fclose(IMFTraj);
    sprintf(filename,"g%sTraj.dx",name);
    IMFTraj=fopen(filename,"w");
    sprintf(filename,"%sTraj.dx",name);
    fprintf(IMFTraj,"object \"Vertices\" class array type float rank 1 shape 3 items %d data file %sTraj.dx\n",IMFNTraj,name);
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

#if 1
   {
    int end;
    int i,j,n;
    int adj;
    int hs;
    int nEnds;
    MFPolytope P;
    int ndim;
    double *x;

#if 1
    fflush(IMFTrajSingle);
    fclose(IMFTrajSingle);
#endif

    fflush(IMFENDPTS);
    fclose(IMFENDPTS);

    sprintf(filename,"%sENDS.dx",name);
    IMFENDS=fopen(filename,"w");

    sprintf(filename,"%sENDPTS.dx",name);
    IMFENDPTS=fopen(filename,"r");

    nEnds=0;
    while(!feof(IMFENDPTS))
     {
      fscanf(IMFENDPTS,"%d\n",&end);
      P=MFChartPolytope(MFAtlasChart(A,end,e),e);
      printf("   chart %d, P=0x%8.8x\n",end,P);fflush(stdout);
      MFPrintPolytope(stdout,P,e);fflush(stdout);
      for(i=0;i<MFPolytopeNumberOfFaces(P,e);i++)
       {
        hs=MFPolytopeFaceIndex(P,i,e);
        adj=MFAtlasHalfSpaceLeftChart(A,hs,e);
        if(adj==end) adj=MFAtlasHalfSpaceRightChart(A,hs,e);
        printf("      face %d, hs=%d, left=%d right=%d\n",i,hs,MFAtlasLeftPolytope(A,hs,e),MFAtlasRightPolytope(A,hs,e));fflush(stdout);
        fprintf(IMFENDS,"   %d %d\n",end,adj);fflush(IMFENDS);
        nEnds++;
       }
     }
    fflush(IMFENDPTS);
    fclose(IMFENDPTS);
    fflush(IMFENDS);
    fclose(IMFENDS);

    sprintf(filename,"g%sENDS.dx",name);
    IMFENDS=fopen(filename,"w");
    fprintf(IMFENDS,"object \"Vertices\" class array type float rank 1 shape 3 items %d data file %sTrajTrajSingledx\n",IMFNTrajSingle,name);
    fprintf(IMFENDS,"object \"Lines\" class array type int rank 1 shape 2 items %d data file %sENDS.dx\n",nEnds,name);
    fprintf(IMFENDS,"attribute \"ref\" string \"positions\"\n");
    fprintf(IMFENDS,"attribute \"element type\" string \"lines\"\n");
    fprintf(IMFENDS,"\n");

    fprintf(IMFENDS,"object \"ends\" class field\n");
    fprintf(IMFENDS,"component \"positions\" file %s.dx,\"Vertices\"\n",name);
    fprintf(IMFENDS,"component \"connections\" value \"Lines\"\n");
    printf("Close gIMFENDS\n");fflush(stdout);
    fclose(IMFENDS);
   }
#endif

  return A;
 }

static int trajectoryNumber=0;

/*! \fn void MFExtendAtlasAlongTrajOnManifold(MFAtlas A, IMFFlow F, MFChart chart, char *name, MFNVector u0, MFKVector p0, double dt, double ti, double tf,double epsilon,MFNRegion Omega, int maxSteps, int maxCharts,double Rmax,double Dmin, int nSkip, MFErrorHandler e)
 *  \brief Integrates a trajectory on a manifold and adds the points to an atlas.
 *
 * \param A The atlas to which new charts will be added.
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
 * \param Dmin The closest (relative to radius approach allowed.
 * \param nSkip The number of steps to make without adding points to the atlas.
 */
void MFExtendAtlasAlongTrajOnManifold(MFAtlas A, IMFFlow F, MFChart chart, char *name, MFNVector u0, MFKVector p0, double dt, double ti, double tf,double epsilon,MFNRegion Omega, int maxSteps, int maxCharts,double Rmax,double Dmin, int nSkip, MFErrorHandler e)
 {
  static char RoutineName[]={"MFExtendAtlasAlongTrajOnManifold"};

  double t,tout;
  MFNKMatrix TS;
  double R;
  int i,j;
  int p,ni;
  double d;
  int nsteps;
  int n;
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
  int end;
  int nearChart;
  double dMin;
  double Rf;
  MFChart c;
  int bail;
  int verbose=1;
  double z[3];
  double *BBCenter=NULL;
  int BBDimension;
  MFImplicitMF M;

  printf("*** %s\n",RoutineName);fflush(stdout);

#ifdef WRITEBVP
   {
    fprintf(IMFBVP,"Trajectory %d\n",trajectoryNumber);
   }
#endif
  trajectoryNumber++;

  if(MFAtlasNumberOfCharts(A,e)>0)
   {
    BBDimension=MFIMFProjectToBB(MFAtlasMF(A,e),MFAtlasChartCenter(A,0,e),NULL,e);
    BBCenter=(double*)malloc(BBDimension*sizeof(double));
   }

  M    =MFAtlasMF(A,e);
  space=MFIMFNSpace(M,e);
  BTree=MFAtlasGetBB(A,e);

/* find closest point to keep trajectories from getting too close */

  if(MFAtlasNumberOfCharts(A,e)>0)
   {
    MFIMFProjectToBB(MFAtlasMF(A,e),u0,BBCenter,e);
    L=MFCreateListOfNearbyCharts(BTree,BBCenter,R,e);
    ni=MFNumberOfIntersectingCharts(L,e);
    nearChart=-1;
    dMin=2*R;
    for(p=0;p<ni;p++)
     {
      j=MFIntersectingChart(L,p,e);
      if(MFChartPaged(MFAtlasChart(A,j,e),e))continue;
      d=MFNSpaceDistance(space,MFAtlasCenterOfChart(A,j,e),u0,e);
      if(d<Dmin*MFAtlasChartRadius(A,j,e)&&(nearChart==-1||d<dMin))
       {
        dMin=d;
        nearChart=j;
       }
     }

    if(verbose)
     {
      printf("   Check that initial point is away from others. dMin=%lf, nearChart=%d\n",dMin,nearChart);fflush(stdout);
      if(nearChart==-1){printf("    point on c is OK, continue.\n");fflush(stdout);}
       else            {printf("    point on c is too close to another, go on to next.\n");fflush(stdout);}
     }
    MFFreeListOfIntersectingCharts(L,e);
   }else{
    if(verbose){printf("   There are no other points, skip check if point is away from others.\n");fflush(stdout);}
    nearChart=-1;
   }
  if(nearChart!=-1)
   {
    printf("    point on c is too close to another. dMin=%lf, nearChart=%d\n",dMin,j);fflush(stdout);
    printf("done IMFExtendAtlasAlongFatTraj\n\n");fflush(stdout);
    if(BBCenter!=NULL)free(BBCenter);
    if(chart!=(MFChart)NULL)MFChartSetSingular(chart,e);
    return;
   }

  t=ti;

  n=IMFFlowNU(F,e);
  k=MFAtlasK(A,e);

  fp=(double*)malloc(n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(fp==NULL)
   {
    sprintf(MFFlowOnManifoldErrorMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(double));
    MFSetError(e,12,RoutineName,MFFlowOnManifoldErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    if(BBCenter!=NULL)free(BBCenter);
    return;
   }
#endif

  f=MFCreateWrappedNVector(n,fp,e);
  g=MFCreateNVector(n,e);
  s=MFCreateKVector(k,e);

/*{
   int m;
   m=MFIMFProjectToDraw(MFAtlasMF(A,e),NULL,NULL,e);
   MFIMFProjectToDraw(MFAtlasMF(A,e),u0,z,e);
   for(i=0;i<m;i++)fprintf(IMFTrajSingle," %lf",z[i]);fprintf(IMFTrajSingle,"\n");fflush(IMFTrajSingle);IMFNTrajSingle++;
  }*/

  nsteps=0;
  u=u0;
  MFRefNVector(u0,e);

  while(1)
   {
    if(verbose)printf(" step %d\n",nsteps);fflush(stdout);

    bail=0;
    if(t>tf){printf("Ending integration, final time reached (%lf>%lf)\n",t,tf);fflush(stdout);bail=1;}

    if(maxSteps>0 && nsteps>=maxSteps){printf("Ending integration, too many steps\n");fflush(stdout);goto FreeAndReturn;}

    if(maxCharts>0 && MFAtlasNumberOfCharts(A,e)>=maxCharts){printf("Ending integration, too many charts\n");fflush(stdout);goto FreeAndReturn;}

    if(!MFNRegionInterior(Omega,u,e)){printf("Ending integration, point is outside Omega. ");MFPrintNVector(stdout,u,e);printf("\n");fflush(stdout);bail=1;}

    IMFEvaluateFlow(F,u,p0,fp,e);
    if(1){printf("u: ");MFPrintNVector(stdout,u,e);printf("\n");fflush(stdout);}
    if(1){printf("f: ");MFPrintNVector(stdout,f,e);printf("\n");fflush(stdout);}
    if(1){printf("|u|: %lf\n",sqrt( MFNV_C(u,0,e)*MFNV_C(u,0,e) + MFNV_C(u,1,e)*MFNV_C(u,1,e) + MFNV_C(u,2,e)*MFNV_C(u,2,e)));fflush(stdout);}
    if(1){printf("u.f: %lf\n",MFNV_C(u,0,e)*MFNV_C(f,0,e) + MFNV_C(u,1,e)*MFNV_C(f,1,e) + MFNV_C(u,2,e)*MFNV_C(f,2,e));fflush(stdout);}
    if(1){printf("fp: (%lf,%lf,%lf)\n",fp[0],fp[1],fp[2]);fflush(stdout);}

    d=sqrt(MFNSpaceInner(space,f,f,e));
    if(0)printf("|F(u,l)|=%le\n",d);fflush(stdout);
    if(d<1.e-5){printf("Ending integration, |F| (%le<%le) too small (near fixed point)\n",d,epsilon);fflush(stdout);bail=1;}

    TS=MFIMFTangentSpace(M,u,e);

    MFNSpaceScale(space,1./d,f,f,e);
    printf("TS=0x%8.8x\n",TS);fflush(stdout);
    printf("f =0x%8.8x\n",f );fflush(stdout);
    printf("g =0x%8.8x\n",g );fflush(stdout);
    printf("s =0x%8.8x\n",s );fflush(stdout);
    MFMVMulT(space,TS,f,s,e);
    MFMVMul(space,TS,s,g,e);
    d=MFNSpaceDistance(space,f,g,e);
    if(0)printf("distance from u+F/|F(u,l)| to the TS = %le\n",d);fflush(stdout);

    R=MFIMFScale(M,u,TS,e);
    Rf=IMFFlowR(F,epsilon,u,p0,TS,e);
    if(Rf<R)R=Rf;

    if(Rmax<1.e-7){Rmax=R;printf("R reset to %lf\n",Rmax);fflush(stdout);}

    if(!(Rf==Rf&&R==R))
     {
      sprintf(MFFlowOnManifoldErrorMsg,"Both Rf and R are nans");
      IMFEvaluateFlow(F,u,p0,fp,e);
      MFSetError(e,12,RoutineName,MFFlowOnManifoldErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      if(BBCenter!=NULL)free(BBCenter);
      return;
     }

    if(R>Rmax||R!=R)R=Rmax;

    c=MFCreateChart(MFAtlasMF(A,e),u,TS,R,e);

/* find closest point to keep trajectories from getting too close */

    nearChart=-1;
    if(t>ti && MFAtlasNumberOfCharts(A,e)>0)
     {
      MFIMFProjectToBB(MFAtlasMF(A,e),MFChartCenter(c,e),BBCenter,e);
      L=MFCreateListOfNearbyCharts(BTree,BBCenter,R,e);
      ni=MFNumberOfIntersectingCharts(L,e);
      nearChart=-1;
      dMin=2*R;
      if(ni==0)dMin=0;
      for(p=0;p<ni;p++)
       {
        j=MFIntersectingChart(L,p,e);
        if(MFChartPaged(MFAtlasChart(A,j,e),e))continue;
        if(!MFIsNear(A,MFAtlasChart(A,j,e),c,e))continue;
  
        d=MFNSpaceDistance(space,MFAtlasCenterOfChart(A,j,e),u,e);
        if(d<Dmin*MFAtlasChartRadius(A,j,e)&&(nearChart==-1||d<dMin))
         {
          dMin=d;
          nearChart=j;
         }
       }
      if(1)
       {
        printf("   Check that this point is away from others. dMin=%lf\n",dMin);fflush(stdout);
        if(nearChart==-1){printf("    point is OK, continue.\n");fflush(stdout);}
         else            {printf("    point is too close to another.\n");fflush(stdout);}
       }
      MFFreeListOfIntersectingCharts(L,e);
      MFFreeChart(c,e);
     }else{
      if(0)printf("   There are no other points, skip check if point is away from others.\n");fflush(stdout);
      nearChart=-1;
     }

/* Passed checks, add chart */

    prev=-1;
    if(nearChart==-1)
     {
      printf("%d) Point %d on traj.",MFAtlasNumberOfCharts(A,e),nsteps);MFPrintNVector(stdout,u,e);printf(", R=%lf, t=%lf\n",R,t);fflush(stdout);
      prev=MFAtlasAddChartWithAll(A,u,TS,R,e);
      end=prev;
       {
        int m;
        m=MFIMFProjectToDraw(MFAtlasMF(A,e),NULL,NULL,e);
        MFIMFProjectToDraw(MFAtlasMF(A,e),u0,z,e);
        for(i=0;i<m;i++)fprintf(IMFTraj," %lf",z[i]);fprintf(IMFTraj,"\n");fflush(IMFTraj);IMFNTraj++;
        MFIMFProjectToDraw(MFAtlasMF(A,e),u,z,e);
        for(i=0;i<m;i++)fprintf(IMFTraj," %lf",z[i]);fprintf(IMFTraj,"\n");fflush(IMFTraj);IMFNTraj++;

        for(i=0;i<m;i++)fprintf(IMFTrajSingle," %lf",z[i]);fprintf(IMFTrajSingle,"\n");fflush(IMFTrajSingle);IMFNTrajSingle++;
       }
#ifdef WRITEBVP
       {
        int i,j;
        int n;
        int k;
        MFNVector phi;

        n=MFAtlasN(A,e);
        k=MFAtlasK(A,e);
        fprintf(IMFBVP,"  chart %d ",prev);
        for(i=0;i<n;i++)fprintf(IMFBVP," %lf",MFNV_C(u0,i,e));fprintf(IMFBVP,"\n");fflush(IMFBVP);
        for(i=0;i<k;i++)
         {
          phi=MFMColumn(TS,i,e);
          fprintf(IMFBVP,"  tangent %d ",i);
          for(j=0;j<n;j++)fprintf(IMFBVP," %lf",MFNV_C(phi,j,e));fprintf(IMFBVP,"\n");fflush(IMFBVP);
          MFFreeNVector(phi,e);
         }
       }
#endif
      if(BTree==NULL)BTree=MFAtlasGetBB(A,e);
      if(BBCenter==NULL)
       {
        BBDimension=MFIMFProjectToBB(MFAtlasMF(A,e),MFAtlasChartCenter(A,0,e),NULL,e);
        BBCenter=(double*)malloc(BBDimension*sizeof(double));
       }
     }else{
      printf("*) Point %d on traj. ",nsteps);MFPrintNVector(stdout,u,e);printf(", R=%lf, t=%lf\n",R,t);fflush(stdout);
      printf("    is too near chart %d (distance = %lf)\n",nearChart,dMin);fflush(stdout);

      prev=MFAtlasAddChartToList(A,MFCreateChart(MFAtlasMF(A,e),u,TS,R,e),e);
      int m;
      m=MFIMFProjectToDraw(MFAtlasMF(A,e),NULL,NULL,e);
      MFIMFProjectToDraw(MFAtlasMF(A,e),u,z,e);
      for(i=0;i<m;i++)fprintf(IMFTraj," %lf",z[i]);fprintf(IMFTraj,"\n");fflush(IMFTraj);IMFNTraj++;

      t=2*tf+1;
      MFChartSetPaged(MFAtlasChart(A,prev,e),e);
     }
    printf("  done adding new point\n");fflush(stdout);

    if(bail){if(prev>-1)MFChartSetSingular(MFAtlasChart(A,prev,e),e);goto FreeAndReturn;}

    nsteps++;

/* Now find the next */

    y0=MFNV_CStar(u,e);
    y=(double*)malloc(n*sizeof(double));

#ifndef MFNOSAFETYNET
    if(y==NULL)
     {
      sprintf(MFFlowOnManifoldErrorMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(double));
      MFSetError(e,12,RoutineName,MFFlowOnManifoldErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      if(BBCenter!=NULL)free(BBCenter);
      return;
     }
#endif

    for(i=0;i<n;i++)y[i]=y0[i];

    d=0;
    tout=t+R;
    while(t+1.e-7<=tout && t<tf && d<1.03*R)
     {
      if(!MEHKellerInt(F,y,p0,t,tout,2,e)){
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

    u0=u;
    u=MFCreateNVectorWithData(n,y,e);
    MFFreeNVector(u0,e);
   }

FreeAndReturn:
  MFFreeNVector(u,e);
  MFFreeNVector(f,e);
  free(fp);
  MFFreeNVector(g,e);
  MFFreeKVector(s,e);

  fprintf(IMFENDPTS,"%d\n",end);fflush(IMFENDPTS);
  printf("done IMFExtendAtlasAlongFatTraj\n\n");fflush(stdout);
  if(BBCenter!=NULL)free(BBCenter);
  return;
 }

/*! \fn MFNVector MFGetInterpolationPoint(MFAtlas A,IMFFlow F, MFKVector p0, MFNRegion Omega, MFChart *c, MFErrorHandler e)
 *  \brief Finds a point on the boundary of the current charts in the atlas at which the flow is outward.
 *
 *  \param A The atlas.
 *  \param F The flow.
 *  \param p0 The parameters for the flow.
 *  \param Omega A region which limits the computation.
 *  \returns A new point satisfying the requirements.
 */
MFNVector MFGetInterpolationPoint(MFAtlas A,IMFFlow F, MFKVector p0, MFNRegion Omega, MFChart *c, MFErrorHandler e)
 {
/*
    edges whose index is all >m give points where edge crosses sphere
    clip points against all faces with index >m not in index set of edge
    mark intersection point
    index of edge is ind1,ind2
*/

  static char RoutineName[]={"MFGetInterpolationPoint"};
  int verbose=0;

  MFPolytope P;
  int iBndChart,nBndCharts;
  int chart;
  int doit;

  int i,j,k,l;
  int n,m;
  int n1,n2;
  double d0,d1,t;
  double r;
  int good0;
  int good1;
  int inter[kmax];
  MFNVector result;
  double R;

  int ii,jj;
  int noncube;
  double a[(kmax-1)*kmax],b[kmax-1],aa[(kmax-1)*(kmax-1)],asave[(kmax-1)*kmax];
  int iface[kmax-1];
  double d[kmax];
  double org[kmax],nrm[kmax];
  double oo;
  double *nn;
  double *f;
  MFKVector vf;
  MFNVector z;
  MFNVector vp;
  double *p;
  MFNVector vF;
  double *pF;
  int nv;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s\n",RoutineName);fflush(stdout);}
#endif

  *c=(MFChart)NULL;

  nBndCharts=MFAtlasNumberOfChartsWithBoundary(A,e);
  if(verbose){printf("There are %d charts on the boundary list\n",nBndCharts);fflush(stdout);}

  k=MFAtlasK(A,e);
  if(k>kmax)
   {
    fprintf(stderr,"%s, I statically allocated some things here for k<10! Complain.\n");
    fflush(stderr);
    abort();
   }

  m=1;
  for(i=0;i<k;i++)m=m*2;
  m=m-1;

  vf=MFCreateKVector(k,e);
  f=MFKV_CStar(vf,e);
  z=MFCreateNVector(MFAtlasN(A,e),e);
  vp=MFCreateNVector(MFAtlasN(A,e),e);
  p=MFNV_CStar(vp,e);
  vF=MFCreateNVector(MFAtlasN(A,e),e);
  pF=MFNV_CStar(vF,e);

/* -------------------------------------------------------------------------------- */

/* Clean boundary list */

  for(iBndChart=0;iBndChart<nBndCharts;iBndChart++)
   {
    chart=MFAtlasChartWithBoundary(A,iBndChart,e);
    P=MFChartPolytope(MFAtlasChart(A,chart,e),e);

    d0=0.;
    for(i=0;i<MFPolytopeNumberOfVertices(P,e);i++)
     {
      d1=MFPolytopeRadiusOfVertex(MFChartPolytope(MFAtlasChart(A,chart,e),e),i,e);
      if(d1>d0)d0=d1;
     }

    R=MFChartRadius(MFAtlasChart(A,chart,e),e);
    if(MFChartIsSingular(MFAtlasChart(A,chart,e),e) ||d0<R )
     {
      MFAtlasRemoveChartFromBoundaryList(A,iBndChart,e);
      iBndChart--;
      nBndCharts--;
     }
   }

/* -------------------------------------------------------------------------------- */

  for(iBndChart=0;iBndChart<nBndCharts;iBndChart++)
   {
    chart=MFAtlasChartWithBoundary(A,iBndChart,e);
    P=MFChartPolytope(MFAtlasChart(A,chart,e),e);

#ifdef MFALLOWVERBOSE
    if(verbose){printf("The %dth chart on the boundary list is chart %d\n",iBndChart,chart);fflush(stdout);}
#endif

    nv=MFPolytopeNumberOfVertices(P,e);

    R=MFChartRadius(MFAtlasChart(A,chart,e),e);
    for(i=0;i<nv;i++)
     {
      n1=MFPolytopeVertexNumberOfIndices(P,i,e);
      for(j=i+1;j<nv;j++)
       {
        n2=MFPolytopeVertexNumberOfIndices(P,j,e);
        m=n1;if(n2>m)m=n2;

        if(m>0)
         {
          n=MFPolytopeIntersectIndexSets(P,i,j,inter,e);
         }else{
          n=0;
         }

        noncube=1;for(l=0;l<n;l++)noncube=noncube&&inter[l]>m;

        if(n==k-1&&noncube)
         {

#ifdef MFALLOWVERBOSE
          if(verbose){printf("  Indices of this non-cube edge:\n");fflush(stdout);}
#endif

          for(ii=0;ii<k-1;ii++)
           {
            iface[ii]=-1;
            for(l=0;l<MFPolytopeNumberOfFaces(P,e);l++)
             {
              if(MFPolytopeFaceIndex(P,l,e)==inter[ii])iface[ii]=l;
             }

#ifdef MFALLOWVERBOSE
            if(verbose){printf("  %d (face %d)\n",inter[ii],iface[ii]);fflush(stdout);}
#endif

           }

          for(ii=0;ii<k-1;ii++)
           {
            b[ii]=MFPolytopeFaceOrigin(P,iface[ii],e);
            nn=MFKV_CStar(MFPolytopeFaceNormal(P,iface[ii],e),e);
            d0;for(l=0;l<k;l++)d0+=nn[l]*nn[l];
            for(l=0;l<k;l++)
             {
              a[ii+(k-1)*l]=nn[l]/d0;
              asave[ii+(k-1)*l]=a[ii+(k-1)*l];
             }
            b[ii]=b[ii]/d0;
           }

          for(ii=0;ii<k-1;ii++)
           {
            for(jj=0;jj<k-1;jj++)
             {
              aa[ii+(k-1)*jj]=0;
              for(l=0;l<k;l++)aa[ii+(k-1)*jj]+=a[ii+(k-1)*l]*a[jj+(k-1)*l];
             }
           }

          MFNullVectorSmall(k,a,nrm,e);
          MFSolveFullSmall(k-1,aa,b,e);

          for(ii=0;ii<k;ii++)
           {
            org[ii]=0.;
            for(l=0;l<k-1;l++)org[ii]+=asave[l+(k-1)*ii]*b[l];
           }

#ifdef MFALLOWVERBOSE
          if(verbose)
           {
            printf("  Done finding edge equation org=(%lf",org[0]);fflush(stdout);
            for(ii=1;ii<k;ii++)printf(",%lf",org[ii]);
            printf(")\n");fflush(stdout);
            printf("                             nrm=(%lf",nrm[0]);fflush(stdout);
            for(ii=1;ii<k;ii++)printf(",%lf",nrm[ii]);
            printf(")\n");fflush(stdout);
    
            printf("Check: is org on all the faces?\n");fflush(stdout);
            d0=0.;for(l=0;l<k;l++)d0+=org[l]*nrm[l];
            printf("            nrm.org=%lf\n",d0);fflush(stdout);
            for(ii=0;ii<k-1;ii++)
             {
              oo=MFPolytopeFaceOrigin(P,iface[ii],e);
              nn=MFKV_CStar(MFPolytopeFaceNormal(P,iface[ii],e),e);
              d0=-oo;for(l=0;l<k;l++)d0+=nn[l]*org[l];
              printf("   face %d: org.n-o=%lf\n",ii,d0);fflush(stdout);
              d0=0.;for(l=0;l<k;l++)d0+=nn[l]*nrm[l];
              printf("            nrm.n=%lf\n",d0);fflush(stdout);
             }
           }
#endif

          r=R*R;for(l=0;l<k;l++)r-=org[l]*org[l];
          if(r>0)
           {
            r=sqrt(r);

#ifdef MFALLOWVERBOSE
            if(verbose)
             {
              d0=0.;for(l=0;l<k;l++)d0+=(org[l]+r*nrm[l])*(org[l]+r*nrm[l]);d0=sqrt(d0);
              printf("  This edge crosses the ball. |pt0|=%lf, R=%lf\n",d0,R);fflush(stdout);
              d0=0.;for(l=0;l<k;l++)d0+=(org[l]+r*nrm[l])*(org[l]+r*nrm[l]);d0=sqrt(d0);
              printf("                              |pt1|=%lf, R=%lf\n",d0,R);fflush(stdout);
             }
#endif

            good0=1;
            good1=1;
            for(ii=0;ii<MFPolytopeNumberOfFaces(P,e);ii++)
             {
              doit=1;for(l=0;l<k-1;l++)if(ii==iface[l])doit=0;
              if(doit)
               {
                oo=MFPolytopeFaceOrigin(P,ii,e);
                nn=MFKV_CStar(MFPolytopeFaceNormal(P,ii,e),e);
                d0=-oo;for(l=0;l<k;l++)d0+=(org[l]+r*nrm[l])*nn[l];
                if(d0>0)good0=0;
                d1=-oo;for(l=0;l<k;l++)d1+=(org[l]-r*nrm[l])*nn[l];
                if(d1>0)good1=0;
               }
             }

#ifdef MFALLOWVERBOSE
            if(verbose)
             {
              if(good0){printf("  The +r point is in the polytope.\n");fflush(stdout);}
                  else {printf("  The +r point is not in the polytope.\n");fflush(stdout);}
             }
#endif

            if(good0)
             {
              for(l=0;l<k;l++)f[l]=org[l]+r*nrm[l];
              MFChartPointInTangentSpace(MFAtlasChart(A,chart,e),vf,vp,e);
              IMFEvaluateFlow(F,vp,p0,pF,e);
              MFChartProjectIntoTangentSpace(MFAtlasChart(A,chart,e),vF,vf,e);

#ifdef MFALLOWVERBOSE
              if(verbose){printf("  flow is ");MFPrintNVector(stdout,vF,e);printf("\n");fflush(stdout);}
              if(verbose){printf("  in TS flow is ");MFPrintKVector(stdout,vf,e);printf("\n");fflush(stdout);}
#endif
              d1=0.;for(l=0;l<k;l++)d1+=f[l]*nrm[l];
              t=-r/d1;

#ifdef MFALLOWVERBOSE
              if(verbose){printf("  f.nrm=%lf\n",d1);fflush(stdout);}
              if(verbose){printf("  t=%lf\n",t);fflush(stdout);}
#endif
              if(t<0.&&(Omega==NULL||MFNRegionInterior(Omega,vp,e)))
               {
                good0=1;
                for(ii=0;ii<k-1;ii++)
                 {
                  d0=0.;for(l=0;l<k;l++)d0+=(org[l]+r*nrm[l]+t*f[l])*asave[ii+(k-1)*l];
                  d1=0.;for(l=0;l<k;l++)d1+= org[l]                 *asave[ii+(k-1)*l];
                  if(d0<0 || d0>d1)good0=0;
                 }
                if(1||good0)
                 {
                  for(l=0;l<k;l++)f[l]=org[l]+r*nrm[l]+t*f[l];
                  MFChartPointInTangentSpace(MFAtlasChart(A,chart,e),vf,z,e);
                  result=MFInterpolateOnDualFace(A,F,p0,z,chart,inter,e);
                  if(result!=NULL)
                   {
#ifdef MFALLOWVERBOSE
                    if(verbose){printf("Result is on chart %d\n",chart);fflush(stdout);}
                    if(verbose){printf("Result is ");MFPrintNVector(stdout,result,e);printf("\n");fflush(stdout);}
#endif
                    MFChartResetChangedFlag(MFAtlasChart(A,chart,e),e);
                    goto FreeAndReturn;
                   }
                 }
               }
             }

#ifdef MFALLOWVERBOSE
            if(verbose)
             {
              if(good1){printf("  The -r point is in the polytope.\n");fflush(stdout);}
                  else {printf("  The +r point is not in the polytope.\n");fflush(stdout);}
             }
#endif

            if(good1)
             {
              for(l=0;l<k;l++)f[l]=org[l]-r*nrm[l];
              MFChartPointInTangentSpace(MFAtlasChart(A,chart,e),vf,vp,e);
              IMFEvaluateFlow(F,vp,p0,pF,e);
              MFChartProjectIntoTangentSpace(MFAtlasChart(A,chart,e),vF,vf,e);

#ifdef MFALLOWVERBOSE
              if(verbose){printf("  flow is ");MFPrintNVector(stdout,vF,e);printf("\n");fflush(stdout);}
              if(verbose){printf("  in TS flow is ");MFPrintKVector(stdout,vf,e);printf("\n");fflush(stdout);}
#endif
              d1=0.;for(l=0;l<k;l++)d1+=f[l]*nrm[l];
              t=r/d1;

#ifdef MFALLOWVERBOSE
              if(verbose){printf("  f.nrm=%lf\n",d1);fflush(stdout);}
              if(verbose){printf("  t=%lf\n",t);fflush(stdout);}
#endif

              if(t<0.&&(Omega==NULL||MFNRegionInterior(Omega,vp,e)))
               {
                good1=1;
                for(ii=0;ii<k-1;ii++)
                 {
                  d0=0.;for(l=0;l<k;l++)d0+=(org[l]-r*nrm[l]+t*f[l])*asave[ii+(k-1)*l];
                  d1=0.;for(l=0;l<k;l++)d1+= org[l]                 *asave[ii+(k-1)*l];
                  if(d0<0 || d0>d1)good1=0;
                 }
                if(1||good1)
                 {
                  for(l=0;l<k;l++)f[l]=org[l]-r*nrm[l]+t*f[l];
                  MFChartPointInTangentSpace(MFAtlasChart(A,chart,e),vf,z,e);
                  result=MFInterpolateOnDualFace(A,F,p0,z,chart,inter,e);
                  if(result!=NULL)
                   {
#ifdef MFALLOWVERBOSE
                    if(verbose){printf("Result is on chart %d\n",chart);fflush(stdout);}
                    if(verbose){printf("Result is ");MFPrintNVector(stdout,result,e);printf("\n");fflush(stdout);}
#endif
                    MFChartResetChangedFlag(MFAtlasChart(A,chart,e),e);
                    goto FreeAndReturn;
                   }
                 }
               }
             }
           }
         }
       }
     }
    MFChartResetChangedFlag(MFAtlasChart(A,chart,e),e);
   }

  result=NULL;
  printf("No such points on the boundary!!!\n");

FreeAndReturn:
  MFFreeKVector(vf,e);
  MFFreeNVector(z,e);
  MFFreeNVector(vp,e);
  MFFreeNVector(vF,e);

  *c=MFAtlasChart(A,chart,e);
  return result;
 }

MFNVector MFInterpolateOnDualFace(MFAtlas A, IMFFlow F, MFKVector p0, MFNVector z, int chart, int *inter, MFErrorHandler e)
 {
  static char RoutineName[]={"MFInterpolateOnDualFace"};
  int i,k,n;
  int j;
  MFNVector c[kmax];
  double a[(kmax-1)*(kmax-1)];
  double t[kmax];
  double d;
  MFNSpace space;
  MFNVector diffi;
  MFNVector diffj;
  MFNVector diffz;
  MFNVector result;

  double *u;
  int good;
  int verbose=0;

  k=MFAtlasK(A,e);
  n=MFAtlasN(A,e);
  space=MFIMFNSpace(MFAtlasMF(A,e),e);

  c[0]=MFAtlasChartCenter(A,chart,e);
  diffi=MFCreateNVector(n,e);
  diffj=MFCreateNVector(n,e);
  diffz=MFCreateNVector(n,e);

  for(i=0;i<k-1;i++)
   {
    j=MFAtlasHalfSpaceRightChart(A,inter[i],e);
    if(j==chart)j=MFAtlasHalfSpaceLeftChart(A,inter[i],e);
    c[i+1]=MFAtlasChartCenter(A,j,e);
   }

  MFNSpaceDirection(space,z,c[0],diffz,e);
  for(i=0;i<k-1;i++)
   {
    MFNSpaceDirection(space,c[i+1],c[0],diffi,e);
    a[i+(k-1)*i]=MFNSpaceInner(space,diffi,diffi,e);
    for(j=i+1;j<k-1;j++)
     {
      MFNSpaceDirection(space,c[j+1],c[0],diffj,e);
      a[i+(k-1)*j]=MFNSpaceInner(space,diffi,diffj,e);
      a[i+(k-1)*j]=a[j+(k-1)*i];
     }
    t[i+1]=MFNSpaceInner(space,diffz,diffi,e);
   }

  MFSolveFullSmall(k-1,a,t+1,e);

  t[0]=1.;for(i=1;i<k;i++)t[0]-=t[i];

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("Barycentric coordinates (%lf",t[0]);
    for(i=1;i<k;i++)printf(",%lf",t[i]);
    printf(")\n");
   }
#endif

  good=1;
  for(i=0;i<k;i++)if(t[i]<0.||t[i]>1.)good=0;
  for(i=0;i<k;i++)if(t[i]<.1||t[i]>.9)good=0;  /* EXPERIMENT */
  if(!good)
   {
    MFFreeNVector(diffi,e);
    MFFreeNVector(diffj,e);
    MFFreeNVector(diffz,e);
    return NULL;
   }

  if(k!=2)
   {
    fprintf(stderr,"Interpolation of more than 2d is not yet supported\n");
    fflush(stderr);
    abort();
   }

  u=(double*)malloc(n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(u==NULL)
   {
    sprintf(MFFlowOnManifoldErrorMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(double));
    MFSetError(e,12,RoutineName,MFFlowOnManifoldErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  for(j=0;j<n;j++)u[j]=0.;

  for(i=0;i<k;i++)
   for(j=0;j<n;j++)u[j]+=t[i]*(MFNV_CStar(c[i],e))[j];

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("Barycentric coordinates (%lf",t[0]);
    for(i=1;i<k;i++)printf(",%lf",t[i]);
    printf(")\n");
   }
#endif

  result=MFCreateNVectorWithData(n,u,e);

  MFFreeNVector(diffi,e);
  MFFreeNVector(diffj,e);
  MFFreeNVector(diffz,e);

  free(u);

  return result;
 }

int MEHKellerInt(IMFFlow F,double *y,MFKVector p0,double ti,double tf,int nsteps, MFErrorHandler e)
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
      sprintf(MFFlowOnManifoldErrorMsg,"Out of memory, trying to allocate %d bytes",neq*sizeof(double));
      MFSetError(e,12,RoutineName,MFFlowOnManifoldErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return 0;
     }
#endif

    ym=(double*)realloc((void*)ym,neq*sizeof(double));

#ifndef MFNOSAFETYNET
    if(ym==NULL)
     {
      sprintf(MFFlowOnManifoldErrorMsg,"Out of memory, trying to allocate %d bytes",neq*sizeof(double));
      MFSetError(e,12,RoutineName,MFFlowOnManifoldErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return 0;
     }
#endif

    y1=(double*)realloc((void*)y1,neq*sizeof(double));

#ifndef MFNOSAFETYNET
    if(y1==NULL)
     {
      sprintf(MFFlowOnManifoldErrorMsg,"Out of memory, trying to allocate %d bytes",neq*sizeof(double));
      MFSetError(e,12,RoutineName,MFFlowOnManifoldErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return 0;
     }
#endif

    f =(double*)realloc((void*)f ,neq*sizeof(double));

#ifndef MFNOSAFETYNET
    if(f==NULL)
     {
      sprintf(MFFlowOnManifoldErrorMsg,"Out of memory, trying to allocate %d bytes",neq*sizeof(double));
      MFSetError(e,12,RoutineName,MFFlowOnManifoldErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return 0;
     }
#endif

    A =(double*)realloc((void*)A ,neq*neq*sizeof(double));

#ifndef MFNOSAFETYNET
    if(A==NULL)
     {
      sprintf(MFFlowOnManifoldErrorMsg,"Out of memory, trying to allocate %d bytes",neq*neq*sizeof(double));
      MFSetError(e,12,RoutineName,MFFlowOnManifoldErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return 0;
     }
#endif

    b =(double*)realloc((void*)b ,neq*sizeof(double));

#ifndef MFNOSAFETYNET
    if(b==NULL)
     {
      sprintf(MFFlowOnManifoldErrorMsg,"Out of memory, trying to allocate %d bytes",neq*sizeof(double));
      MFSetError(e,12,RoutineName,MFFlowOnManifoldErrorMsg,__LINE__,__FILE__);
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

void MFNullVectorSmall(int k,double *A,double *phi, MFErrorHandler e)
 {
  static char RoutineName[]={"MFNullVectorSmall"};

#ifdef HAVE_LAPACK

/*Find the null vector of a (k-1)xk matrix */

  double t0,t1;

  if(k==2)
   {
    t1=sqrt(A[0]*A[0]+A[1]*A[1]);
    t0    = A[1]/t1;
    phi[1]=-A[0]/t1;
    phi[0]=t0;
   }else if(k==3)
   {
    t0=sqrt(A[0]*A[0]+A[2]*A[2]+A[4]*A[4]);
    t1=sqrt(A[1]*A[1]+A[3]*A[3]+A[5]*A[5]);
    phi[0]=(A[2]*A[5]-A[3]*A[4])/t0/t1;
    phi[1]=(A[4]*A[1]-A[5]*A[0])/t0/t1;
    phi[2]=(A[0]*A[3]-A[1]*A[2])/t0/t1;
   }else{  
    int i,j;
    char jobvl='N';
    char jobvr='A';
    int lda;
    double *s;
    double *U;
    int ldU;
    double *V;
    int ldV;
    double *work;
    int lwork;
    int info=0;
    int m;
    double t;

    lda=k-1;
    ldU=k-1;
    U=NULL;

    ldV=k;
    V=(double*)malloc(k*(k-1)*sizeof(double));

#ifndef MFNOSAFETYNET
  if(V==NULL)
   {
    sprintf(MFFlowOnManifoldErrorMsg,"Out of memory, trying to allocate %d bytes",k*(k-1)*sizeof(double));
    MFSetError(e,12,RoutineName,MFFlowOnManifoldErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

    s=(double*)malloc(k*sizeof(double));

#ifndef MFNOSAFETYNET
  if(s==NULL)
   {
    sprintf(MFFlowOnManifoldErrorMsg,"Out of memory, trying to allocate %d bytes",k*sizeof(double));
    MFSetError(e,12,RoutineName,MFFlowOnManifoldErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

    i=-1;
    info=0;
    jobvl='N';
    jobvr='A';
    CALLDGESVD(&jobvl,&jobvr,&lda,&k,A,&lda,s,U,&ldU,V,&ldV,&t,&i,&info);
    lwork=round(t);

    work=(double*)malloc(lwork*sizeof(double));

#ifndef MFNOSAFETYNET
  if(work==NULL)
   {
    sprintf(MFFlowOnManifoldErrorMsg,"Out of memory, trying to allocate %d bytes",lwork*sizeof(double));
    MFSetError(e,12,RoutineName,MFFlowOnManifoldErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

    info=0;
    jobvl='N';
    jobvr='A';
    CALLDGESVD(&jobvl,&jobvr,&lda,&k,A,&lda,s,U,&ldU,V,&ldV,work,&lwork,&info);
  
    m=0;
    for(i=0;i<k;i++)
     {
      if(fabs(s[i])<1.e-7)
       {
        for(j=0;j<k;j++)phi[j]=V[i+ldV*j];
        m++;
       }
     }
    free(V);
    free(s);
    free(work);
   }
  return;
#else
  sprintf(MFFlowOnManifoldErrorMsg,"IMFInterpolation requires dgesvd from Lapack");
  MFSetError(e,12,RoutineName,MFFlowOnManifoldErrorMsg,__LINE__,__FILE__);
  return;
#endif
 }

void MFSolveFullSmall(int k,double *A,double *b, MFErrorHandler e)
 {
  static char RoutineName[]={"MFSolveFullSmall"};
/*solves the kxk system of linear eqs.;*/
  double t;

  if(k==1)
   {
    b[0]=b[0]/A[0];
   }else if(k==2)
   {
    t   =(A[3]*b[0]-A[2]*b[1])/(A[0]*A[3]-A[1]*A[2]);
    b[1]=(A[0]*b[1]-A[1]*b[0])/(A[0]*A[3]-A[1]*A[2]);
    b[0]=t;
   }else{  
    MFSolveFull(k,A,b,e);
   }

  return;
 }

int MFIsNear(MFAtlas A,MFChart c0,MFChart c1,MFErrorHandler e)
 {
  static char RoutineName[]={"IsNear"};
  MFNVector u0,u1;
  double t0,t1;

  u0=MFChartCenter(c0,e);
  u1=MFChartCenter(c1,e);

  return 1;
 }

#ifdef __cplusplus
}
#endif
