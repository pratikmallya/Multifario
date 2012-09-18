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

static char *id="@(#) $Id: IMFIntegrateFat.c,v 1.8 2011/07/21 17:42:46 mhender Exp $";

static char IMFIntegrateFatErrorMsg[256];

#include <MFAtlas.h>
#include <MFBinaryTree.h>
#include <MFChart.h>
#include <IMFFlow.h>
#include <IMFExpansion.h>
#include <IMFExpansionSpace.h>
#include <IMFExpansionPt.h>
#include <IMFIntegrateFat.h>
#include <MFPrint.h>
#include <math.h>
#include <stdio.h>
#include <IMF.h>

#ifdef __cplusplus
 extern "C" {
#endif

#define DOTRAJ

#define round(x) ((x)>=0?(long)((x)+0.5):(long)((x)-0.5))

MFBinaryTree MFAtlasGetBB(MFAtlas,MFErrorHandler);
MFChart MFAtlasChart(MFAtlas,int,MFErrorHandler);

FILE *IMFTraj=NULL;
int   IMFNTraj=0;

static int MEHKellerInt(IMFFlow,double*,MFKVector,double,double,int,MFErrorHandler);
int IMFIsNear(MFAtlas,MFChart,MFChart,MFErrorHandler);
void MFChartSetPaged(MFChart,MFErrorHandler);
int ADoTorusProject(MFNVector u, double *x, void *d, MFErrorHandler e);

void IMFExtendAtlasAlongFatTraj(MFAtlas I, IMFFlow F, char *name, MFNVector u0, MFKVector p0, double dt, double ti, double tf,double epsilon,MFNRegion Omega, int maxSteps, int maxCharts,double Rmax,int nSkip, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFExtendAtlasAlongFatTraj"};

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

/* find closest point to keep trajectories from getting too close */

  if(MFAtlasNumberOfCharts(I,e)>0)
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
      if(d<.5*MFAtlasChartRadius(I,j,e)&&(nearChart==-1||d<dMin))
       {
        dMin=d;
        nearChart=j;
       }
     }
    printf("   Check that initial point is away from others. dMin=%lf\n",dMin);fflush(stdout);
    if(nearChart==-1){printf("    point on c is OK, continue.\n");fflush(stdout);}
     else            {printf("    point on c is too close to another, go on to next.\n");fflush(stdout);}
    MFFreeListOfIntersectingCharts(L,e);
    if(nearChart!=-1){if(BBCenter!=NULL)free(BBCenter);return;}
   }else{
    printf("   There are no other points, skip check if point is away from others.\n");fflush(stdout);
    nearChart=-1;
   }
  if(nearChart!=-1)
   {
    if(BBCenter!=NULL)free(BBCenter);
    return;
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
    return;
   }
#endif

  f=MFCreateWrappedNVector(n,fp,e);
  g=MFCreateNVector(n,e);
  s=MFCreateKVector(IMFExpansionK(E,e),e);

  nsteps=0;
  u=u0;
  MFRefNVector(u,e);
  while(1)
   {
    printf(" step %d\n",nsteps);fflush(stdout);

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
    if(MFAtlasNumberOfCharts(I,e)>0)
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
        if(d<.5*MFAtlasChartRadius(I,j,e)&&(nearChart==-1||d<dMin))
         {
          dMin=d;
          nearChart=j;
         }
       }
      printf("   Check that this point is away from others. dMin=%lf\n",dMin);fflush(stdout);
      if(nearChart==-1){printf("    point is OK, continue.\n");fflush(stdout);}
       else            {printf("    point is too close to another.\n");fflush(stdout);}
      MFFreeListOfIntersectingCharts(L,e);
      MFFreeChart(c,e);
     }else{
      printf("   There are no other points, skip check if point is away from others.\n");fflush(stdout);
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
      printf("%d) Point %d on fat traj.",MFAtlasNumberOfCharts(I,e),nsteps);MFPrintNVector(stdout,u,e);printf(", R=%lf, t=%lf\n",R,t);fflush(stdout);
      prev=MFAtlasAddChartWithAll(I,u,TS,R,e);
      if(BTree==NULL)BTree=MFAtlasGetBB(I,e);
      if(BBCenter==NULL)
       {
        BBDimension=MFIMFProjectToBB(MFAtlasMF(I,e),MFAtlasChartCenter(I,0,e),NULL,e);
        BBCenter=(double*)malloc(BBDimension*sizeof(double));
        printf("%s, BBDimension=%d\n",RoutineName,BBDimension);fflush(stdout);
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
      printf("*) Point %d on fat traj. ",nsteps);MFPrintNVector(stdout,u,e);printf(", R=%lf, t=%lf\n",R,t);fflush(stdout);
      if(nsteps<nSkip){printf("    is being skipped, %d<%d\n",nsteps,nSkip);fflush(stdout);}
        else {printf("    is too near chart %d (distance = %lf)\n",nearChart,dMin);fflush(stdout);}
      prev=MFAtlasAddChartToList(I,MFCreateChart(MFAtlasMF(I,e),u,TS,R,e),e);
      t=2*tf+1;
      if(0&&MFAtlasNumberOfCharts(I,e)%1000==0){printf("Paging\n");MFAtlasPageOutChartsNotNearBoundary(I,1,0,name,e);}
      IMFExpansionNVSetType(u,3,e);
      MFChartSetPaged(MFAtlasChart(I,prev,e),e);
     }
    printf("  done adding new point\n");fflush(stdout);

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
      return;
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

/*printf("done IMFExtendAtlasAlongFatTraj\n");fflush(stdout);*/
  if(BBCenter!=NULL)free(BBCenter);
  return;
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

MFChart IMFStepAlongFatTraj(MFAtlas I, IMFFlow F, MFNVector u0, MFKVector p0, double R0, double dt, double epsilon,MFNRegion Omega, double Rmax, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFStepAlongFatTraj"};

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
  double *BBCenter=NULL;
  int BBDimension;

  if(MFAtlasNumberOfCharts(I,e)>0)
   {
    BBDimension=MFIMFProjectToBB(MFAtlasMF(I,e),MFAtlasChartCenter(I,0,e),NULL,e);
    BBCenter=(double*)malloc(BBDimension*sizeof(double));
    printf("%s, BBDimension=%d\n",RoutineName,BBDimension);fflush(stdout);
   }

  E=IMFExpansionNVGetE(u0,e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("Initial expansion\n");IMFPrintExpansion(stdout,E,e);fflush(stdout);}
#endif

  neq=IMFExpansionDataLn(E,e);
  space=MFIMFNSpace(MFAtlasMF(I,e),e);
  BTree=MFAtlasGetBB(I,e);

  t=0.;
  R=R0;

  n=IMFFlowNU(F,e);
  k=IMFExpansionK(E,e);

  Fp=IMFCreateFatFlow(F,k,e);

  fp=(double*)malloc(n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(fp==NULL)
   {
    sprintf(IMFIntegrateFatErrorMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(double));
    MFSetError(e,12,RoutineName,IMFIntegrateFatErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    if(BBCenter!=NULL)free(BBCenter);
    return 0;
   }
#endif

  f=MFCreateWrappedNVector(n,fp,e);
  g=MFCreateNVector(n,e);
  s=MFCreateKVector(IMFExpansionK(E,e),e);

  nsteps=0;
  u=u0;
  MFRefNVector(u,e);
  y=(double*)malloc(IMFExpansionDataLn(E,e)*sizeof(double));

#ifndef MFNOSAFETYNET
  if(y==NULL)
   {
    sprintf(IMFIntegrateFatErrorMsg,"Out of memory, trying to allocate %d bytes",IMFExpansionDataLn(E,e)*sizeof(double));
    MFSetError(e,12,RoutineName,IMFIntegrateFatErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    if(BBCenter!=NULL)free(BBCenter);
    return 0;
   }
#endif

  while(1)
   {
    E=IMFExpansionNVGetE(u,e);

/* find the next */

    y0=IMFExpansionData(E,e);
    for(i=0;i<IMFExpansionDataLn(E,e);i++)y[i]=y0[i];

    d=0;
    tout=t+R;
    while(t+1.e-7<=tout && d<1.03*R)
     {
      if(!MEHKellerInt(Fp,y,p0,t,tout,2,e)){
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
    E=IMFCreateExpansion(n,k,e);
    IMFExpansionSetDerivatives(E,y,y+n,y+n*k,NULL,e);

    u0=u;
    u=IMFCreateExpansionNVector(E,t,IMFExpansionNVGetSigma(u0,e),prev,0,e);
    MFFreeNVector(u0,e);
    IMFFreeExpansion(E,e);
    bail=0;

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
    if(BBCenter!=NULL)
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
        if(MFAtlasChartRadius(I,j,e)&&(nearChart==-1||d<dMin))
         {
          dMin=d;
          nearChart=j;
         }
       }
      printf("   Check that this point is away from others. dMin=%lf\n",dMin);fflush(stdout);
      if(nearChart==-1){printf("    point is OK, continue.\n");fflush(stdout);}
       else            {printf("    point is too close to another.\n");fflush(stdout);}
      MFFreeListOfIntersectingCharts(L,e);
     }else{
      printf("   There are no other points, skip check if point is away from others.\n");fflush(stdout);
      nearChart=-1;
     }

/* Passed checks, add chart */

    prev=-1;
    if((nearChart==-1||IMFExpansionNVGetType(u,e)==1)||IMFExpansionNVGetType(u,e)==2)
     {
      goto FreeAndReturn;
     }else{
      printf("*) Point %d on fat traj. ",nsteps);MFPrintNVector(stdout,u,e);printf(", R=%lf, t=%lf\n",R,t);fflush(stdout);
      printf("    is too near chart %d (distance = %lf)\n",nearChart,dMin);fflush(stdout);
      printf("    is too near chart %d (distance = %lf)\n",nearChart,dMin);fflush(stdout);
      bail=1;
     }

    if(bail){if(prev>-1)MFChartSetSingular(MFAtlasChart(I,prev,e),e);MFFreeChart(c,e);c=NULL;goto FreeAndReturn;}

    MFFreeChart(c,e);
    nsteps++;
   }

FreeAndReturn:
  MFFreeNVector(u,e);
  MFFreeNVector(f,e);
  free(fp);
  IMFFreeFlow(Fp,e);
  MFFreeNVector(g,e);
  MFFreeKVector(s,e);
  free(y);

/*printf("done %s\n",RoutineName);fflush(stdout);*/
  if(BBCenter!=NULL)free(BBCenter);
  return c;
 }

#ifdef __cplusplus
}
#endif
