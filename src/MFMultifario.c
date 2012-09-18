/*
 *  PROGRAM NAME:  multifario
 *
 *  (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
 *  CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
 *
 *  Please refer to the LICENSE file in the top directory
 *
 *      author: Mike Multifario mhender@watson.ibm.com
 *      date:   April 1, 2002
 */

static char *id="@(#) $Id: MFMultifario.c,v 1.9 2011/07/21 17:42:46 mhender Exp $";

static char MFMultifarioErrorMsg[256];

#include <MFAtlas.h>
#include <MFContinuationMethod.h>
#include <MFMultifariosMethod.h>
#include <MFImplicitMF.h>
#include <MFNRegion.h>
#include <MFNVector.h>
#include <MFPrint.h>
#include <stdio.h>
#include <string.h>

#ifdef __cplusplus
 extern "C" {
#endif

static FILE *restartReg=NULL;
static FILE *restartSing=NULL;
static double *x=NULL;

void MFAtlasPageOutChartsNotNearBoundary(MFAtlas,int,int,char*,MFErrorHandler);
void MFAtlasWriteOutChartsNotPaged(MFAtlas,int,int,char*,MFErrorHandler);
int MFNVectorGetNRefs(MFNVector,MFErrorHandler);
MFChart MFAtlasChart(MFAtlas,int,MFErrorHandler);

void MFAtlasPageOutAllCharts(MFAtlas,MFErrorHandler);

static void MFMultifarioExtendAtlas(MFContinuationMethod H,MFAtlas A, MFImplicitMF M, MFNRegion Omega, int m, MFNVector *u0, MFNKMatrix *Phi,MFErrorHandler);
static void MFMultifarioMethodCopyClipF(MFContinuationMethod,MFAtlas,MFErrorHandler);
static void MFMultifarioCloseAtlas(MFContinuationMethod,MFAtlas,MFErrorHandler);
static void MFMultifarioFlushAtlas(MFContinuationMethod,MFAtlas,MFErrorHandler);
static void MFMultifarioFreeParmBlock(void*,MFErrorHandler);

void MFMultifarioProcessBifurctation(MFAtlas S, int chart, MFNVector u, MFNKMatrix Phi,MFErrorHandler);

struct MFMultifarioParmBlockST {
                             int verbose;
                             int maxCharts;
                             double minR;
                             double maxR;
                             double epsilon;
                             double dotmin;
                             int page;
                             int pageEvery;
                             int useBB;
                             int dumpToPlotFile;
                             int dumpToCenterFile;
                             int dumpToRestartFile;
                             int dumpToRestartFileEvery;
                             char *fileName;
                             int checkPoint;
                             int checkPointEvery;
                             int branchSwitch;
                             int nClipF;
                             int mClipF;
                             double (**clipF)(MFNVector,MFErrorHandler);
                            };
typedef struct MFMultifarioParmBlockST *MFMultifarioParmBlock;

MFContinuationMethod MFCreateMultifariosMethod(MFErrorHandler e)
 {
  static char RoutineName[]={"MFCreateMultifariosMethod"};
  MFContinuationMethod result;
  MFMultifarioParmBlock parms;

  parms=(struct MFMultifarioParmBlockST*)malloc(sizeof(struct MFMultifarioParmBlockST));

#ifndef MFNOSAFETYNET
  if(parms==NULL)
   {
    sprintf(MFMultifarioErrorMsg,"Out of space trying to allocate %d bytes.",sizeof(struct MFMultifarioParmBlockST));
    MFSetError(e,12,RoutineName,MFMultifarioErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    MFFreeContinuationMethod(result,e);
    return NULL;
   }
#endif

  parms->verbose=0;
  parms->maxCharts=-1;
  parms->minR=.01;
  parms->maxR=1.;
  parms->epsilon=1.e-7;
  parms->dotmin=.2;
  parms->page=1;
  parms->pageEvery=1000;
  parms->useBB=1;
  parms->dumpToPlotFile=1;
  parms->dumpToCenterFile=1;
  parms->dumpToRestartFile=1;
  parms->dumpToRestartFileEvery=100;
  parms->checkPoint=0;
  parms->checkPointEvery=100;
  parms->branchSwitch=1;
  parms->fileName=(char*)malloc(6*sizeof(char));

#ifndef MFNOSAFETYNET
  if(parms->fileName==NULL)
   {
    sprintf(MFMultifarioErrorMsg,"Out of space trying to allocate %d bytes.",6*sizeof(char));
    MFSetError(e,12,RoutineName,MFMultifarioErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    MFFreeContinuationMethod(result,e);
    return NULL;
   }
#endif

  strcpy(parms->fileName,"Atlas");

  parms->nClipF=0;
  parms->mClipF=0;
  parms->clipF=NULL;

  result=MFCreateContinuationMethodBaseClass("Multifario",e);
  MFContinuationMethodSetParmBlock(result,parms,e);
  MFContinuationMethodSetFreeParmBlock(result,MFMultifarioFreeParmBlock,e);
  MFContinuationMethodSetCloseAtlas(result,MFMultifarioCloseAtlas,e);
  MFContinuationMethodSetExtendAtlasMultipleWithTangents(result,MFMultifarioExtendAtlas,e);
  MFContinuationMethodSetFlushAtlas(result,MFMultifarioFlushAtlas,e);

  return result;
 }

void MFMultifarioFreeParmBlock(void *parmBlock, MFErrorHandler e)
 { 
  MFMultifarioParmBlock parms;

  parms=(MFMultifarioParmBlock)(parmBlock);
  if(parms->fileName!=NULL)free(parms->fileName);
  if(parms->clipF!=NULL)free(parms->clipF);
  free(parms);

  return;
 }

char *MFMultifarioGetFilename(MFContinuationMethod H, MFErrorHandler e)
 {
  static char RoutineName[]={"MFMultifarioGetFileName"};
  MFMultifarioParmBlock parms;

#ifdef MFNOCONFIDENCE
  if(H==NULL)
   {
    sprintf(MFMultifarioErrorMsg,"MFContinuation Method (argument 1) is NULL.");
    MFSetError(e,12,RoutineName,MFMultifarioErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(strcmp(MFContinuationMethodGetType(H,e),"Multifario"))
   {
    sprintf(MFMultifarioErrorMsg,"MFContinuation Method (argument 1) is not Multifario's Method.");
    MFSetError(e,12,RoutineName,MFMultifarioErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  parms=(MFMultifarioParmBlock)(MFContinuationMethodGetParmBlock(H,e));
  return parms->fileName;
 }

void MFMultifarioSetFilename(MFContinuationMethod H,char *name, MFErrorHandler e)
 {
  static char RoutineName[]={"MFMultifarioSetFileName"};
  MFMultifarioParmBlock parms;

#ifdef MFNOCONFIDENCE
  if(H==NULL)
   {
    sprintf(MFMultifarioErrorMsg,"MFContinuation Method (argument 1) is NULL.");
    MFSetError(e,12,RoutineName,MFMultifarioErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(strcmp(MFContinuationMethodGetType(H,e),"Multifario"))
   {
    sprintf(MFMultifarioErrorMsg,"MFContinuation Method (argument 1) is not Multifario's Method.");
    MFSetError(e,12,RoutineName,MFMultifarioErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  parms=(MFMultifarioParmBlock)(MFContinuationMethodGetParmBlock(H,e));
  free(parms->fileName);
  parms->fileName=(char*)malloc((strlen(name)+1)*sizeof(char));

#ifndef MFNOSAFETYNET
  if(parms->fileName==NULL)
   {
    sprintf(MFMultifarioErrorMsg,"Out of space trying to allocate %d bytes.",(strlen(name)+1)*sizeof(char));
    MFSetError(e,12,RoutineName,MFMultifarioErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  strcpy(parms->fileName,name);
  return;
 }

void MFMultifarioAddClipF(MFContinuationMethod H,double (*clipF)(MFNVector,MFErrorHandler), MFErrorHandler e)
 {
  static char RoutineName[]={"MFMultifarioAddClipF"};
  MFMultifarioParmBlock thisMethod;
  int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s\n",RoutineName);fflush(stdout);}
#endif

  thisMethod=(MFMultifarioParmBlock)(MFContinuationMethodGetParmBlock(H,e));

  if(thisMethod->nClipF>=thisMethod->mClipF)
   {
    thisMethod->mClipF+=10;
    thisMethod->clipF=(double (**)(MFNVector,MFErrorHandler))realloc(thisMethod->clipF,(thisMethod->mClipF)*sizeof(double (*)(MFNVector,MFErrorHandler)));

#ifndef MFNOSAFETYNET
  if(thisMethod->clipF==NULL)
   {
    sprintf(MFMultifarioErrorMsg,"Out of space trying to allocate %d bytes.",(thisMethod->mClipF)*sizeof(double (*)(MFNVector)));
    MFSetError(e,12,RoutineName,MFMultifarioErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

   }
  thisMethod->clipF[thisMethod->nClipF]=clipF;
  thisMethod->nClipF++;
 }

void MFMultifarioMethodCopyClipF(MFContinuationMethod H,MFAtlas A, MFErrorHandler e)
 {
  static char RoutineName[]={"MFMultifariosMethodCopyClipF"};
  int i;
  MFMultifarioParmBlock thisMethod;
  int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s\n",RoutineName);fflush(stdout);}
#endif

  thisMethod=(MFMultifarioParmBlock)(MFContinuationMethodGetParmBlock(H,e));

  MFAtlasClearClipF(A,e);
  for(i=0;i<thisMethod->nClipF;i++)MFAtlasAddClipF(A,thisMethod->clipF[i],e);

  return;
 }

void MFMultifariosMethodClearClipF(MFContinuationMethod H, MFErrorHandler e)
 {
  static char RoutineName[]={"MFMultifariosMethodClearClipF"};
  MFMultifarioParmBlock thisMethod;
  int i;
  int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s\n",RoutineName);fflush(stdout);}
#endif

  thisMethod=(MFMultifarioParmBlock)(MFContinuationMethodGetParmBlock(H,e));

  for(i=0;i<thisMethod->nClipF;i++)thisMethod->clipF[i]=NULL;

  thisMethod->nClipF=0;

  return;
 }

void MFMultifarioExtendAtlas(MFContinuationMethod H, MFAtlas S, MFImplicitMF M, MFNRegion Omega, int m, MFNVector *u0, MFNKMatrix *Phi0, MFErrorHandler e)
 {
  static char RoutineName[]={"MFMultifarioExtendAtlas"};
  int i,j;
  int kkk;
  double delta=0.;
  int chart;
  int n,N;
  MFNVector u;
  MFNKMatrix Phi;
  double R,Rs;
  FILE *fid;
  double volume;
  double epsilon;
  double dotmin;
  double minR;
  int kmax;
  int page,pageEvery;
  int checkPoint,checkPointEvery;
  int branchSwitch;
  int dumpToPlotFile;
  int dumpToCenterFile;
  int rc;
  int verbose;
  int switched;
  int dumpToRestartFile;
  int dumpToRestartFileEvery;
  char *name;
  char FileName[4096];
  int bifpt;

/* If M!=AtlasM change it. What if operations between n vectors don't work? */

  if(m==0)return;

  epsilon=MFMultifarioGetRealParameter(H,"epsilon",e);
  dotmin=MFMultifarioGetRealParameter(H,"dotmin",e);
  minR=MFMultifarioGetRealParameter(H,"minR",e);
  if(MFIMFGetRMin(M,e)<minR)minR=MFIMFGetRMin(M,e);                   /* !!!! Added to accomadate AUTO */
  kmax=MFMultifarioGetIntegerParameter(H,"maxCharts",e);
  verbose=MFMultifarioGetIntegerParameter(H,"verbose",e);
  page=MFMultifarioGetIntegerParameter(H,"page",e);
  pageEvery=MFMultifarioGetIntegerParameter(H,"pageEvery",e);
  checkPoint=MFMultifarioGetIntegerParameter(H,"checkPoint",e);
  checkPointEvery=MFMultifarioGetIntegerParameter(H,"checkPointEvery",e);
  branchSwitch=MFMultifarioGetIntegerParameter(H,"branchSwitch",e);
  dumpToPlotFile=MFMultifarioGetIntegerParameter(H,"dumpToPlotFile",e);
  dumpToCenterFile=MFMultifarioGetIntegerParameter(H,"dumpToCenterFile",e);
  name=MFMultifarioGetFilename(H,e);
  switched=0;
  MFMultifarioMethodCopyClipF(H,S,e);
  dumpToRestartFile=MFMultifarioGetIntegerParameter(H,"dumpToRestartFile",e);
  dumpToRestartFileEvery=MFMultifarioGetIntegerParameter(H,"dumpToRestartFileEvery",e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s, epsilon=%lf, kmax=%d, verbose=%d\n",RoutineName,epsilon,kmax,verbose);fflush(stdout);}
#endif

  kkk=0;
  MFAtlasSetEpsilon(S,epsilon,e);
  MFAtlasSetDotMin(S,dotmin,e);
  MFAtlasSetVerbose(S,verbose,e);
  MFAtlasSetRMin(S,minR,e);
  MFIMFSetRMin(M,minR,e);
  for(i=0;i<m;i++)
   {
    if(Phi0==NULL||Phi0[i]==NULL)
     {
      rc=MFAtlasAddChartWithCenter(S,u0[i],e);
      MFIMFSetStability(MFAtlasMF(S,e),u0[i],MFAtlasChartTangentSpace(S,rc,e),e);
      if(rc<0){printf("line %d in file %s, AddChartWithCenter returned rc=%d\n",__LINE__,__FILE__,rc);fflush(stdout);}
     }else{
      R=MFIMFScale(M,u0[i],Phi0[i],e);
      MFIMFSetStability(MFAtlasMF(S,e),u0[i],Phi0[i],e);
      rc=MFAtlasAddChartWithAll(S,u0[i],Phi0[i],R,e);
      if(rc<0){printf("line %d in file %s, AddChartWithAll returned rc=%d\n",__LINE__,__FILE__,rc);fflush(stdout);}
     }

#ifdef MFALLOWVERBOSE
    if(verbose){printf("%d) Add initial points to Atlas ",rc);MFPrintNVector(stdout,u0[i],e);printf(" index %d, bif? %d\n",MFNVGetIndex(u0[i],e),MFNVGetIndex2(u0[i],e));fflush(stdout);}
    if(verbose){printf("     R=%lf\n",MFAtlasChartRadius(S,rc,e));fflush(stdout);}
#endif

    MFChartSetReferenceNumber(MFAtlasChart(S,rc,e),0,e);
#ifdef MFVOLUME
    volume=MFVolumeOfAtlas(S,Omega,e);

#ifdef MFALLOWVERBOSE
    if(verbose){printf("  Volume of Atlas %lf\n",volume);fflush(stdout);}
#endif

#endif
    kkk++;

    if(MFChartReferenceNumber(MFAtlasChart(S,rc,e),e)>dumpToRestartFileEvery)
     {
/*    printf("Reference Number has exceeded the threshold\n");fflush(stdout);*/
      MFChartSetReferenceNumber(MFAtlasChart(S,rc,e),0,e);
      MFUpdateNeighborsReferenceMarks(S,rc,dumpToRestartFileEvery,e);
     }
    if(dumpToRestartFile && MFChartReferenceNumber(MFAtlasChart(S,rc,e),e)==0)
     {
      if(MFChartIsSingular(MFAtlasChart(S,rc,e),e))
       {
        if(restartSing==NULL)
         {
          strcpy(FileName,name);
          strcat(FileName,".singular");
          restartSing=fopen(FileName,"w");
          if(restartSing==NULL)abort();
         }
        N=MFIMFProjectToDraw(MFAtlasMF(S,e),NULL,NULL,e);
        if(x==NULL)
         {
          x=(double*)malloc(N*sizeof(double));

#ifndef MFNOSAFETYNET
          if(x==NULL)
           {
            sprintf(MFMultifarioErrorMsg,"Out of space trying to allocate %d bytes.",N*sizeof(double));
            MFSetError(e,12,RoutineName,MFMultifarioErrorMsg,__LINE__,__FILE__);
            MFErrorHandlerOutOfMemory(e);
            return;
           }
#endif
         }

        MFIMFProjectToDraw(MFAtlasMF(S,e),MFChartCenter(MFAtlasChart(S,rc,e),e),x,e);
        for(j=0;j<N;j++)fprintf(restartSing," %lf",x[j]);
        fprintf(restartSing,"\n");
       }else{
        if(restartReg==NULL)
         {
          strcpy(FileName,name);
          strcat(FileName,".regular");
          restartReg=fopen(FileName,"w");
          if(restartReg==NULL)abort();
         }
        N=MFIMFProjectToDraw(MFAtlasMF(S,e),NULL,NULL,e);
        if(x==NULL)
         {
          x=(double*)malloc(N*sizeof(double));

#ifndef MFNOSAFETYNET
          if(x==NULL)
           {
            sprintf(MFMultifarioErrorMsg,"Out of space trying to allocate %d bytes.",N*sizeof(double));
            MFSetError(e,12,RoutineName,MFMultifarioErrorMsg,__LINE__,__FILE__);
            MFErrorHandlerOutOfMemory(e);
            return;
           }
#endif

         }
        MFIMFProjectToDraw(MFAtlasMF(S,e),MFChartCenter(MFAtlasChart(S,rc,e),e),x,e);
        for(j=0;j<N;j++)fprintf(restartReg," %lf",x[j]);
        fprintf(restartReg,"\n");
       }
     }

    if(page && kkk%pageEvery==0)
     {
      MFAtlasPageOutChartsNotNearBoundary(S,dumpToPlotFile,dumpToCenterFile,name,e);
     }
    if(checkPoint && kkk%checkPointEvery==0)
     {

#ifdef MFALLOWVERBOSE
      if(verbose){printf("Check pointing, step %d\n",kkk);fflush(stdout);}
#endif

      fid=fopen("Checkpoint.atlas","w");
      MFWriteAtlas(fid,S,e);
      fclose(fid);
     }
   }

  n=MFAtlasN(S,e);
  u=MFCloneNVector(u0[0],e);

GoOn:
  rc=1;
  while(rc>0&&(kmax<0||kkk<kmax)&&(chart=MFAtlasPointOnBoundaryInsideRegion(S,Omega,u,&Phi,&delta,e))>-1)
   {

/* PointOnBoundary always returns an interior point or -1. It must also somehow add co-dimension one boundary charts. */

    R=MFIMFScale(MFAtlasMF(S,e),u,Phi,e);
    Rs=2*MFAtlasChartSuggestedRadius(S,chart,e);
    if(0&&R>Rs)R=Rs;

    bifpt=0;
    if(MFNVGetIndex2(u,e)<0)
     {
      bifpt=1;
      MFNVSetIndex2(u,0,e);
     }

    if(!bifpt)MFIMFSetStability(MFAtlasMF(S,e),u,Phi,e);

    if((rc=MFAtlasAddChartWithAll(S,u,Phi,R,e))>0)
     {

#ifdef MFALLOWVERBOSE
      if(verbose){printf("%d) Add boundary point to Atlas ",rc);MFPrintNVector(stdout,u,e);printf(" R=%lf index %d, bif? %d\n",R,MFNVGetIndex(u,e),MFNVGetIndex2(u,e));fflush(stdout);}
#endif

      if(bifpt)
       {
        MFMultifarioProcessBifurctation(S,chart,u,Phi,e);
       }

      if(MFChartReferenceNumber(MFAtlasChart(S,rc,e),e)>dumpToRestartFileEvery)
       {
/*      printf("Reference Number has exceeded the threshold\n");fflush(stdout);*/
        MFChartSetReferenceNumber(MFAtlasChart(S,rc,e),0,e);
        MFUpdateNeighborsReferenceMarks(S,rc,dumpToRestartFileEvery,e);
       }
/*    printf(" Reference Number is %d\n",MFChartReferenceNumber(MFAtlasChart(S,rc,e),e));*/
      if(dumpToRestartFile && MFChartReferenceNumber(MFAtlasChart(S,rc,e),e)==0)
       {
        if(MFChartIsSingular(MFAtlasChart(S,rc,e),e))
         {
          if(restartSing==NULL)
           {
            strcpy(FileName,name);
            strcat(FileName,".singular");
            restartSing=fopen(FileName,"w");
            if(restartSing==NULL)abort();
           }
          N=MFIMFProjectToDraw(MFAtlasMF(S,e),NULL,NULL,e);
          if(x==NULL)
           {
            x=(double*)malloc(N*sizeof(double));

#ifndef MFNOSAFETYNET
            if(x==NULL)
             {
              sprintf(MFMultifarioErrorMsg,"Out of space trying to allocate %d bytes.",N*sizeof(double));
              MFSetError(e,12,RoutineName,MFMultifarioErrorMsg,__LINE__,__FILE__);
              MFErrorHandlerOutOfMemory(e);
              return;
             }
#endif
           }

          MFIMFProjectToDraw(MFAtlasMF(S,e),MFChartCenter(MFAtlasChart(S,rc,e),e),x,e);
          for(j=0;j<N;j++)fprintf(restartSing," %lf",x[j]);
          fprintf(restartSing,"\n");
         }else{
          if(restartReg==NULL)
           {
            strcpy(FileName,name);
            strcat(FileName,".regular");
            restartReg=fopen(FileName,"w");
            if(restartReg==NULL)abort();
           }
          N=MFIMFProjectToDraw(MFAtlasMF(S,e),NULL,NULL,e);
          if(x==NULL)
           {
            x=(double*)malloc(N*sizeof(double));

#ifndef MFNOSAFETYNET
            if(x==NULL)
             {
              sprintf(MFMultifarioErrorMsg,"Out of space trying to allocate %d bytes.",N*sizeof(double));
              MFSetError(e,12,RoutineName,MFMultifarioErrorMsg,__LINE__,__FILE__);
              MFErrorHandlerOutOfMemory(e);
              return;
             }
#endif
           }
          MFIMFProjectToDraw(MFAtlasMF(S,e),MFChartCenter(MFAtlasChart(S,rc,e),e),x,e);
          for(j=0;j<N;j++)fprintf(restartReg," %lf",x[j]);
          fprintf(restartReg,"\n");
         }
       }
#ifdef MFVOLUME
      volume=MFVolumeOfAtlas(S,Omega);

#ifdef MFALLOWVERBOSE
      if(verbose){printf("  Volume of Atlas %lf\n",volume);fflush(stdout);}
#endif

#endif

      kkk++;
      if(page && kkk%pageEvery==0)
       {
        MFAtlasPageOutChartsNotNearBoundary(S,dumpToPlotFile,dumpToCenterFile,name,e);
       }
      if(checkPoint && kkk%checkPointEvery==0)
       {

#ifdef MFALLOWVERBOSE
        if(verbose){printf("Check pointing, step %d\n",kkk);fflush(stdout);}
#endif

        fid=fopen("Checkpoint.atlas","w");
        MFWriteAtlas(fid,S,e);
        fclose(fid);
       }
     }
    MFFreeNVector(u,e);
    MFFreeNKMatrix(Phi,e);
    u=MFCloneNVector(u0[0],e);
   }

#ifdef MFALLOWVERBOSE
  if(verbose&&chart<0){printf("No more points on the boundary, switch=%d, switched=%d\n",branchSwitch,switched);fflush(stdout);}
   else if(verbose){printf("Maximum number of charts added (%d), kmax=%d\n",kkk,kmax);fflush(stdout);}
#endif

/*branchSwitch=0;*/

  if((switched<branchSwitch) && branchSwitch!=0 && chart<0)
   {
    int schart,bchart;
    MFNVector u1,u;
    MFNKMatrix Phi;

    schart=MFAtlasGetSingularChartWithBoundary(S,Omega,e);
    if(schart>0)
     {

#ifdef MFALLOWVERBOSE
      if(1||verbose){printf("Branch switching, singular chart with boundary is %d, index %d, Radius %lf\n",schart,MFNVGetIndex2(MFAtlasChartCenter(S,schart,e),e),MFAtlasChartRadius(S,schart,e));fflush(stdout);}
#endif

/* This is the branch that passes on through. */

      u1=MFAtlasGetPointOnBoundaryChart(S,Omega,schart,1.,e);
      if(u1==NULL)goto GoOn;

      Phi=MFAtlasChartTangentSpace(S,schart,e);
      R=MFIMFScale(MFAtlasMF(S,e),u1,Phi,e);
      MFIMFSetStability(MFAtlasMF(S,e),u1,Phi,e);
      bchart=MFAtlasAddChartWithAll(S,u1,Phi,R,e);

#ifdef MFALLOWVERBOSE
      if(verbose&&bchart>-1){printf("%d) Add bifurcating chart to Atlas ",bchart);MFPrintNVector(stdout,u1,e);printf(" index %d, bif? %d\n",MFNVGetIndex(u1,e),MFNVGetIndex2(u1,e));fflush(stdout);}
#endif

      switched+=1;

/*    if this is not a smooth continuation of a branch (-99), then somehow need to remove the (-98) charts once a point on the bifurcating sheet has been found. But can't remove them all, because they may all not generate the same branch.*/
/* For each, generate a point on either side and interpolate the tangent normal to the singular curve? */
/* How to get the branch structure? Boundaries have no bifurcating branch. Nor do things like Hopf and bifs to different classes */

      goto GoOn;
     }
   }

  if(u!=NULL)MFFreeNVector(u,e);

  return;
 }

void MFMultifarioCloseAtlas(MFContinuationMethod H, MFAtlas S, MFErrorHandler e)
 {
  static char RoutineName[]={"MFCloseAtlas"};
  int page;
  int dumpToPlotFile;
  int dumpToCenterFile;
  char *name;

  page=MFMultifarioGetIntegerParameter(H,"page",e);
  dumpToPlotFile=MFMultifarioGetIntegerParameter(H,"dumpToPlotFile",e);
  dumpToCenterFile=MFMultifarioGetIntegerParameter(H,"dumpToCenterFile",e);
  name=MFMultifarioGetFilename(H,e);

  if(dumpToPlotFile)MFAtlasWriteOutChartsNotPaged(S,dumpToPlotFile,dumpToCenterFile,name,e);
  if(restartSing!=NULL){fclose(restartSing);restartSing=NULL;}
  if(restartReg!=NULL){fclose(restartReg);restartReg=NULL;}
  if(x!=NULL){free(x);x=NULL;}

  MFAtlasClosePlotfile(S,e);
  MFAtlasCloseCenterfile(S,e);

  return;
 }

void MFMultifarioFlushAtlas(MFContinuationMethod H, MFAtlas S, MFErrorHandler e)
 {
  static char RoutineName[]={"MFFlushAtlas"};
  int page;
  int dumpToPlotFile;
  int dumpToCenterFile;
  char *name;

  page=MFMultifarioGetIntegerParameter(H,"page",e);
  dumpToPlotFile=MFMultifarioGetIntegerParameter(H,"dumpToPlotFile",e);
  dumpToCenterFile=MFMultifarioGetIntegerParameter(H,"dumpToCenterFile",e);
  name=MFMultifarioGetFilename(H,e);

  if(page)MFAtlasWriteOutChartsNotPaged(S,dumpToPlotFile,dumpToCenterFile,name,e);

  return;
 }

int MFMultifarioSetIntegerParameter(MFContinuationMethod M, char *parameterName, int value, MFErrorHandler e)
 {
  static char RoutineName[]={"MFMultifarioSetIntegerParameter"};
  struct MFMultifarioParmBlockST *data;

  if(strcmp(MFContinuationMethodGetType(M,e),"Multifario"))
   {
    sprintf(MFMultifarioErrorMsg,"Parameter 1 must be a Multifarios algorithm");
    MFSetError(e,4,RoutineName,MFMultifarioErrorMsg,__LINE__,__FILE__);
    return 0;
   }

  data=(struct MFMultifarioParmBlockST*)MFContinuationMethodGetParmBlock(M,e);

  if(!strcmp(parameterName,"verbose"))data->verbose=value;
   else if(!strcmp(parameterName,"maxCharts"))data->maxCharts=value;
   else if(!strcmp(parameterName,"page"))data->page=value;
   else if(!strcmp(parameterName,"pageEvery"))data->pageEvery=value;
   else if(!strcmp(parameterName,"useBB"))data->useBB=value;
   else if(!strcmp(parameterName,"dumpToPlotFile"))data->dumpToPlotFile=value;
   else if(!strcmp(parameterName,"dumpToCenterFile"))data->dumpToCenterFile=value;
   else if(!strcmp(parameterName,"dumpToRestartFile"))data->dumpToRestartFile=value;
   else if(!strcmp(parameterName,"dumpToRestartFileEvery"))data->dumpToRestartFileEvery=value;
   else if(!strcmp(parameterName,"checkPoint"))data->checkPoint=value;
   else if(!strcmp(parameterName,"checkPointEvery"))data->checkPointEvery=value;
   else if(!strcmp(parameterName,"branchSwitch"))data->branchSwitch=value;
   else return 0;

  return 1;
 }

int MFMultifarioSetRealParameter(MFContinuationMethod M, char *parameterName, double value, MFErrorHandler e)
 {
  static char RoutineName[]={"MFMultifarioSetRealParameter"};
  struct MFMultifarioParmBlockST *data;

  if(strcmp(MFContinuationMethodGetType(M,e),"Multifario"))
   {
    sprintf(MFMultifarioErrorMsg,"Parameter 1 must be an Multifario");
    MFSetError(e,4,RoutineName,MFMultifarioErrorMsg,__LINE__,__FILE__);
    return 0;
   }

  data=(struct MFMultifarioParmBlockST*)MFContinuationMethodGetParmBlock(M,e);

  if(!strcmp(parameterName,"minR"))data->minR=value;
   else if(!strcmp(parameterName,"maxR"))data->maxR=value;
   else if(!strcmp(parameterName,"epsilon"))data->epsilon=value;
   else if(!strcmp(parameterName,"dotmin"))data->dotmin=value;
   else return 0;

  return 1;
 }

int MFMultifarioGetIntegerParameter(MFContinuationMethod M, char *parameterName, MFErrorHandler e)
 {
  static char RoutineName[]={"MFMultifarioGetIntegerParameter"};
  struct MFMultifarioParmBlockST *data;

  if(strcmp(MFContinuationMethodGetType(M,e),"Multifario"))
   {
    sprintf(MFMultifarioErrorMsg,"Parameter 1 must be an Multifario");
    MFSetError(e,4,RoutineName,MFMultifarioErrorMsg,__LINE__,__FILE__);
    return 0;
   }

  data=(struct MFMultifarioParmBlockST*)MFContinuationMethodGetParmBlock(M,e);

  if(!strcmp(parameterName,"verbose"))return data->verbose;
   else if(!strcmp(parameterName,"maxCharts"))return data->maxCharts;
   else if(!strcmp(parameterName,"page"))return data->page;
   else if(!strcmp(parameterName,"pageEvery"))return data->pageEvery;
   else if(!strcmp(parameterName,"useBB"))return data->useBB;
   else if(!strcmp(parameterName,"dumpToPlotFile"))return data->dumpToPlotFile;
   else if(!strcmp(parameterName,"dumpToCenterFile"))return data->dumpToCenterFile;
   else if(!strcmp(parameterName,"dumpToRestartFile"))return data->dumpToRestartFile;
   else if(!strcmp(parameterName,"dumpToRestartFileEvery"))return data->dumpToRestartFileEvery;
   else if(!strcmp(parameterName,"checkPoint"))return data->checkPoint;
   else if(!strcmp(parameterName,"checkPointEvery"))return data->checkPointEvery;
   else if(!strcmp(parameterName,"branchSwitch"))return data->branchSwitch;
   else{
    sprintf(MFMultifarioErrorMsg,"Parameter %s is not a integer valued parameter",parameterName);
    MFSetError(e,4,RoutineName,MFMultifarioErrorMsg,__LINE__,__FILE__);
    return 0;
   }

  return 0;
 }

double MFMultifarioGetRealParameter(MFContinuationMethod M, char *parameterName, MFErrorHandler e)
 {
  static char RoutineName[]={"MFMultifarioGetRealParameter"};
  struct MFMultifarioParmBlockST *data;

  if(strcmp(MFContinuationMethodGetType(M,e),"Multifario"))
   {
    sprintf(MFMultifarioErrorMsg,"Parameter 1 must be an Multifario");
    MFSetError(e,4,RoutineName,MFMultifarioErrorMsg,__LINE__,__FILE__);
    return 0;
   }

  data=(struct MFMultifarioParmBlockST*)MFContinuationMethodGetParmBlock(M,e);

  if(!strcmp(parameterName,"minR"))return data->minR;
   else if(!strcmp(parameterName,"maxR"))return data->maxR;
   else if(!strcmp(parameterName,"epsilon"))return data->epsilon;
   else if(!strcmp(parameterName,"dotmin"))return data->dotmin;
   else{
    sprintf(MFMultifarioErrorMsg,"Parameter %s is not a real valued parameter",parameterName);
    MFSetError(e,4,RoutineName,MFMultifarioErrorMsg,__LINE__,__FILE__);
    return 0.;
   }

  return 0.;
 }

#ifdef __cplusplus
}
#endif

double MFNKMatrixDot(MFNSpace space, MFNKMatrix A, MFNKMatrix B, MFErrorHandler e);

void MFMultifarioProcessBifurctation(MFAtlas A, int chart, MFNVector u, MFNKMatrix Phi, MFErrorHandler e)
 {
  static char RoutineName[]={"MFMultifarioProcessBifurcation"};

/* --------------------- Got back a singular chart -------------------------------------- */

  double R;
  int bchart;
  int i;
  MFNKMatrix Tan;
  MFNVector x,y;
  int dir;
  MFChart C;
  double vnorm;
  int rc;
  int verbose=0;

  MFNVector w;
  MFNVector v;
  MFNVector up;
  MFKVector s;
  MFNKMatrix Psi;

  if(1||verbose){printf("%s      Singular Point Detected and located.\n",RoutineName);fflush(stdout);}

  R=MFChartRadius(MFAtlasChart(A,chart,e),e);

/* v is the normal to the co--dimension 1 singular "boundary manifold" */

/* Phi is the estimated tangent of the bifurcating branch */

  w=MFCloneNVector(u,e);
  rc=MFIMFSingular(MFAtlasMF(A,e),u,Phi,w,e); /* w is the "other" null vector */

#ifdef MFALLOWVERBOSE
  if(verbose){printf("     return code from MFIMFSingular %d\n",rc);fflush(stdout);}
#endif

/* Compute v, the normal to the co--dimension 1 singular "boundary manifold" */

  v=MFCloneNVector(w,e);
  up=MFCloneNVector(w,e);
  s=MFCreateKVector(MFAtlasK(A,e),e);

  MFNSpaceScale(MFAtlasNSpace(A,e),.5*MFChartRadius(MFAtlasChart(A,chart,e),e),w,up,e);
  MFNSpaceAdd(MFAtlasNSpace(A,e),u,up,up,e);
  Psi=MFIMFTangentSpace(MFAtlasMF(A,e),up,e);

/* y is the vector from the center of the original chart to the singular point */

  y=MFCloneNVector(u,e);
  R=1.*MFChartRadius(MFAtlasChart(A,chart,e),e);
  MFNSpaceDirection(MFAtlasNSpace(A,e),u,MFChartCenter(MFAtlasChart(A,chart,e),e),y,e);
  dir=1;
  if(MFNSpaceInner(MFAtlasNSpace(A,e),v,y,e)<0)dir=-1;
  MFFreeNVector(y,e);

/* Do this if it is a bifurcation point (not limit point or special point) */

#ifndef MFNOSAFETYNET
  if(Psi==NULL)
   {
    printf("*** MFIMFTangentSpace failed at predicted point, exempting vertex %d of chart %d\n",v,chart);fflush(stdout);
    return;
   }
#endif

  MFMVMulT(MFAtlasNSpace(A,e),Psi,w,s,e);        /* v= Psi Psi^T phi_k */
  MFMVMul (MFAtlasNSpace(A,e),Psi,s,v,e);

  MFMVMulT(MFAtlasNSpace(A,e),Phi,v,s,e);        /* v= Phi Phi^T Psi Psi^T phi_k */
  MFMVMul (MFAtlasNSpace(A,e),Phi,s,v,e);

  MFFreeKVector(s,e);

  vnorm=sqrt(MFNSpaceInner(MFAtlasNSpace(A,e),v,v,e));

#ifndef MFNOSAFETYNET
  if(fabs(vnorm)<1.e-10)
   {
    sprintf(MFMultifarioErrorMsg,"Norm of v is zero.");
    MFSetError(e,12,RoutineName,MFMultifarioErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  MFNSpaceScale(MFAtlasNSpace(A,e),1./vnorm,v,v,e);

/* done computing v */

/* make the singular curve normal the first basis vector */

  Tan=MFCloneNKMatrix(Phi,e);
  MFGramSchmidtReplace(MFAtlasNSpace(A,e),Tan,v,w,e);

  MFFreeNKMatrix(Tan,e);

/* --------------------------- disks on the biufurcating branch ----------------------------------------------- */

  x=MFCloneNVector(u,e);
  C=MFCreateChart(MFAtlasMF(A,e),x,Psi,R,e);
  MFFreeNVector(x,e);

  s=MFCreateKVector(MFAtlasK(A,e),e);
  for(i=0;i<MFAtlasK(A,e);i++)MFKVSetC(s,i, R/sqrt(MFAtlasK(A,e)),e);
  for(i=0;i<MFAtlasK(A,e);i++)MFKVSetC(s,i, 0.    ,e);
                              MFKVSetC(s,0, 0.8*R ,e);

  x=MFCloneNVector(u,e);
  rc=MFChartEvaluate(C,s,x,e);
  Tan=MFIMFTangentSpaceWithGuess(MFAtlasMF(A,e),x,Psi,e);
  rc=MFAtlasAddChart(A,MFCreateChart(MFAtlasMF(A,e),x,Tan,R,e),e);
  if(rc>-1){printf("%d) Bifurcating branch          ",rc);MFPrintNVector(stdout,x,e);printf(" index %d, bif? %d\n",MFNVGetIndex(x,e),MFNVGetIndex2(x,e));fflush(stdout);}
  MFFreeKVector(s,e);
  MFFreeNVector(x,e);

/*printf("Dot between tan %d and %d is %lf\n",rc,chart,MFNKMatrixDot(MFAtlasNSpace(A,e),Tan,Phi,e));*/

  MFFreeNKMatrix(Tan,e);

  s=MFCreateKVector(MFAtlasK(A,e),e);
  for(i=0;i<MFAtlasK(A,e);i++)MFKVSetC(s,i,-R/sqrt(MFAtlasK(A,e)),e);
  for(i=0;i<MFAtlasK(A,e);i++)MFKVSetC(s,i, 0.    ,e);
                              MFKVSetC(s,0,-0.8*R ,e);

  x=MFCloneNVector(u,e);
  rc=MFChartEvaluate(C,s,x,e);
  Tan=MFIMFTangentSpaceWithGuess(MFAtlasMF(A,e),x,Psi,e);
  rc=MFAtlasAddChart(A,MFCreateChart(MFAtlasMF(A,e),x,Tan,R,e),e);
  if(rc>-1){printf("%d) Bifurcating branch          ",rc);MFPrintNVector(stdout,x,e);printf(" index %d, bif? %d\n",MFNVGetIndex(x,e),MFNVGetIndex2(x,e));fflush(stdout);}
  MFFreeKVector(s,e);
  MFFreeNVector(x,e);
  MFFreeNKMatrix(Tan,e);

/*printf("Dot between tan %d and %d is %lf\n",rc,chart,MFNKMatrixDot(MFAtlasNSpace(A,e),Tan,Phi,e));*/

/*MFFreeChart(C,e);*/

  MFFreeNVector(v,e);
  MFFreeNVector(w,e);
  MFFreeNKMatrix(Psi,e);
  return;
 }

#if 0
/* Create BoundaryChart creates a half disk */

  if(fabs(vnorm) < 1.e-10)
   {
    x=MFCloneNVector(u,e);
    MFNVSetIndex2(x,-99,e);   /* This depends on the type of the point */

    R=1.*MFChartRadius(MFAtlasChart(A,chart,e),e);
    {printf("%d) Add chart to Atlas ",MFAtlasNumberOfCharts(A,e));MFPrintNVector(stdout,u,e);printf(" index %d, bif? %d\n",MFNVGetIndex(u,e),MFNVGetIndex2(u,e));fflush(stdout);}
    rc=MFAtlasAddChart(A,MFCreateChart(MFAtlasMF(A,e),x,Phi,R,e),e);
    MFFreeNVector(x,e);

    MFIMFSetStability(MFAtlasMF(A,e),w,Psi,e);
    MFNVSetIndex2(w,-98,e);
    bchart=MFAtlasAddChartWithAll(A,w,Psi,R,e);
   
    for(i=1;i<MFAtlasK(A,e);i++)MFKVSetC(s,i,0.,e);MFKVSetC(s,0,2.8*R,e);
   
    x=MFCloneNVector(u,e);
    MFNVSetIndex2(x,-98,e);
    rc=MFChartEvaluate(MFAtlasChart(A,bchart,e),s,x,e);
    Tan=MFIMFTangentSpaceWithGuess(MFAtlasMF(A,e),x,Psi,e);
    MFIMFSetStability(MFAtlasMF(A,e),x,Tan,e);
    rc=MFAtlasAddChart(A,MFCreateBoundaryChart(MFAtlasMF(A,e),x,Tan,1, 1,R,e),e);
    rc=MFAtlasAddChart(A,MFCreateBoundaryChart(MFAtlasMF(A,e),x,Tan,1,-1,R,e),e);
   
    MFFreeNVector(x,e);
    MFFreeNKMatrix(Tan,e);
   
    if(rc>-1){printf("%d) Bifurcating branch          ",rc);MFPrintNVector(stdout,x,e);printf(" index %d, bif? %d\n",MFNVGetIndex(x,e),MFNVGetIndex2(x,e));fflush(stdout);}
   
    MFKVSetC(s,0,-1.8*R,e);
    x=MFCloneNVector(u,e);
    MFNVSetIndex2(x,-98,e);
    rc=MFChartEvaluate(MFAtlasChart(A,bchart,e),s,x,e);
    Tan=MFIMFTangentSpaceWithGuess(MFAtlasMF(A,e),x,Psi,e);
    MFIMFSetStability(MFAtlasMF(A,e),x,Tan,e);
    rc=MFAtlasAddChart(A,MFCreateBoundaryChart(MFAtlasMF(A,e),x,Tan,1, 1,R,e),e);
    rc=MFAtlasAddChart(A,MFCreateBoundaryChart(MFAtlasMF(A,e),x,Tan,1,-1,R,e),e);
   
    MFFreeNKMatrix(Tan,e);
   
    if(rc>-1){printf("%d) Bifurcating branch          ",rc);MFPrintNVector(stdout,x,e);printf(" index %d, bif? %d\n",MFNVGetIndex(x,e),MFNVGetIndex2(x,e));fflush(stdout);}
   
   }else{
    if(1||verbose){printf("      Singular Point Detected and located. |vnorm| > 1.e-10\n");fflush(stdout);}

    R=1.*MFChartRadius(MFAtlasChart(A,chart,e),e);

/* V is the normal to the co--dimension 1 singular "boundary manifold" */

    MFNSpaceScale(MFAtlasNSpace(A,e),1./vnorm,v,v,e);
    MFGramSchmidtReplace(MFAtlasNSpace(A,e),Phi,v,v,e);         /* make the singular curve normal the first basis vector */

/* Phi is the estimated tangent of the bifurcating branch */

/* y is the vector from the center of the original chart to the singular point */

    y=MFCloneNVector(u,e);
    R=1.*MFChartRadius(MFAtlasChart(A,chart,e),e);
    MFNSpaceDirection(MFAtlasNSpace(A,e),u,MFChartCenter(MFAtlasChart(A,chart,e),e),y,e);
    dir=1;
    if(MFNSpaceInner(MFAtlasNSpace(A,e),v,y,e)<0)dir=-1;
    MFFreeNVector(y,e);

    MFNVSetIndex2(MFAtlasCenterOfChart(A,chart,e),-1,e);

/* ---------------------- Half disk forward into the current branch ------------------------------------------- */
#if 0
    x=MFCloneNVector(u,e);
    MFNVSetIndex2(x,-99,e);   /* This depends on the type of the point, says not to continue this until switch */

    rc=MFAtlasAddChart(A,MFCreateBoundaryChart(MFAtlasMF(A,e),x,Phi,1,-dir,R,e),e);
    {printf("%d) singular half disk back into branch  ",rc);MFPrintNVector(stdout,x,e);printf(" index %d, bif? %d\n",MFNVGetIndex(x,e),MFNVGetIndex2(x,e));fflush(stdout);}
    MFFreeNVector(x,e);
#endif

#if 0
/* ---------------------- Original code ----------------------------------------------------------------------- */

/* ---------------------- Full disk at the singular point ----------------------------------------------------- */
  
    MFNSpaceScale(MFAtlasNSpace(A,e),1./vnorm,V,V,e);
    MFGramSchmidtReplace(MFAtlasNSpace(A,e),Phi,V,w,e);
   
    R=1.*MFChartRadius(MFAtlasChart(A,chart,e),e);
    x=MFCloneNVector(u,e);
    MFIMFSetStability(MFAtlasMF(A,e),w,Phi,e);
    MFNVSetIndex2(w,-98,e);
    C=MFCreateChart(MFAtlasMF(A,e),w,Phi,R,e);
    rc=MFAtlasAddChart(A,C,e);
    {printf("%d) Add boundary chart to Atlas ",rc);MFPrintNVector(stdout,u,e);printf(" index %d, bif? %d\n",MFNVGetIndex(u,e),MFNVGetIndex2(u,e));fflush(stdout);}
   
    MFFreeChart(C,e);


/* ---------------------- Half disks on the current branch ---------------------------------------------------- */
   
    s=MFCreateKVector(MFAtlasK(A,e),e);
    for(i=1;i<MFAtlasK(A,e);i++)MFKVSetC(s,i,0.,e);MFKVSetC(s,0,1.8*R,e);
   
    x=MFCloneNVector(u,e);
    MFNVSetIndex2(x,-98,e);
    rc=MFChartEvaluate(MFAtlasChart(A,bchart,e),s,x,e);
    Tan=MFIMFTangentSpaceWithGuess(MFAtlasMF(A,e),x,Phi,e);
    MFIMFSetStability(MFAtlasMF(A,e),x,Tan,e);

    rc=MFAtlasAddChart(A,MFCreateBoundaryChart(MFAtlasMF(A,e),x,Tan,1, 1,R,e),e);
    if(rc>-1){printf("%d) Bifurcating branch          ",rc);MFPrintNVector(stdout,x,e);printf(" index %d, bif? %d\n",MFNVGetIndex(x,e),MFNVGetIndex2(x,e));fflush(stdout);}

    rc=MFAtlasAddChart(A,MFCreateBoundaryChart(MFAtlasMF(A,e),x,Tan,1,-1,R,e),e);
    if(rc>-1){printf("%d) Bifurcating branch          ",rc);MFPrintNVector(stdout,x,e);printf(" index %d, bif? %d\n",MFNVGetIndex(x,e),MFNVGetIndex2(x,e));fflush(stdout);}

    MFFreeNVector(x,e);
    MFFreeNKMatrix(Tan,e);
    MFFreeKVector(s,e);

/* ---------------------- Half disks on the biufurcating branch ----------------------------------------------- */
   
   
    s=MFCreateKVector(MFAtlasK(A,e),e);
    for(i=1;i<MFAtlasK(A,e);i++)MFKVSetC(s,i,0.,e);MFKVSetC(s,0,-1.8*R,e);
    x=MFCloneNVector(u,e);
    MFNVSetIndex2(x,-98,e);

    rc=MFChartEvaluate(MFAtlasChart(A,bchart,e),s,x,e);
    Tan=MFIMFTangentSpaceWithGuess(MFAtlasMF(A,e),x,Phi,e);
    MFIMFSetStability(MFAtlasMF(A,e),x,Tan,e);

    rc=MFAtlasAddChart(A,MFCreateBoundaryChart(MFAtlasMF(A,e),x,Tan,1, 1,R,e),e);
    if(rc>-1){printf("%d) Bifurcating branch          ",rc);MFPrintNVector(stdout,x,e);printf(" index %d, bif? %d\n",MFNVGetIndex(x,e),MFNVGetIndex2(x,e));fflush(stdout);}

    rc=MFAtlasAddChart(A,MFCreateBoundaryChart(MFAtlasMF(A,e),x,Tan,1,-1,R,e),e);
    if(rc>-1){printf("%d) Bifurcating branch          ",rc);MFPrintNVector(stdout,x,e);printf(" index %d, bif? %d\n",MFNVGetIndex(x,e),MFNVGetIndex2(x,e));fflush(stdout);}
   
    MFFreeNVector(x,e);
    MFFreeNKMatrix(Tan,e);
    MFFreeKVector(s,e);
   
#endif


#if 0
/* This computes the bifurcating branches. It's been disabled because it should be put into it's own cell. */

  
    MFNSpaceScale(MFAtlasNSpace(A,e),1./vnorm,V,V,e);
    MFGramSchmidtReplace(MFAtlasNSpace(A,e),Phi,V,w,e);
   
    R=1.*MFChartRadius(MFAtlasChart(A,chart,e),e);
    x=MFCloneNVector(u,e);
   
    MFIMFSetStability(MFAtlasMF(A,e),x,Phi,e);
    MFNVSetIndex2(x,-98,e);
    C=MFCreateChart(MFAtlasMF(A,e),x,Phi,R,e);
    MFFreeChart(x,e);

    s=MFCreateKVector(MFAtlasK(A,e),e);
    for(i=1;i<MFAtlasK(A,e);i++)MFKVSetC(s,i,0.,e);MFKVSetC(s,0,1.8*R,e);
   
    x=MFCloneNVector(u,e);
    MFNVSetIndex2(x,-98,e);
    rc=MFChartEvaluate(C,s,x,e);
    Tan=MFIMFTangentSpaceWithGuess(MFAtlasMF(A,e),x,Phi,e);
    MFIMFSetStability(MFAtlasMF(A,e),x,Tan,e);
    rc=MFAtlasAddChart(A,MFCreateChart(MFAtlasMF(A,e),x,Tan,R,e),e);
    if(rc>-1){printf("%d) Bifurcating branch          ",rc);MFPrintNVector(stdout,x,e);printf(" index %d, bif? %d\n",MFNVGetIndex(x,e),MFNVGetIndex2(x,e));fflush(stdout);}

    MFFreeNVector(x,e);
    MFFreeNKMatrix(Tan,e);
   
    MFKVSetC(s,0,-1.8*R,e);
    x=MFCloneNVector(u,e);
    MFNVSetIndex2(x,-98,e);
    rc=MFChartEvaluate(C,s,x,e);
    Tan=MFIMFTangentSpaceWithGuess(MFAtlasMF(A,e),x,Phi,e);
    MFIMFSetStability(MFAtlasMF(A,e),x,Tan,e);
    rc=MFAtlasAddChart(A,MFCreateChart(MFAtlasMF(A,e),x,Tan,R,e),e);
    if(rc>-1){printf("%d) Bifurcating branch          ",rc);MFPrintNVector(stdout,x,e);printf(" index %d, bif? %d\n",MFNVGetIndex(x,e),MFNVGetIndex2(V,e));fflush(stdout);}

    MFFreeNVector(x,e);
    MFFreeNKMatrix(Tan,e);
    MFFreeKVector(s,e);
    MFFreeChart(C,e);
#endif
   

#if 0

/* -----------------------Half Disk on the new branch ------------------------------------------- */

    MFGramSchmidtReplace(MFAtlasNSpace(A,e),Phi,V,w,e);         /* replace the singular curve tangent with the first basis vector */
  
    x=MFCloneNVector(u,e);
    MFNSpaceScale(MFAtlasNSpace(A,e),1.,u,x,e);
    C=MFCreateChart(MFAtlasMF(A,e),x,Phi,.2*R,e);
    MFFreeNVector(x,e); 

    s=MFCreateKVector(MFAtlasK(A,e),e);
    for(i=1;i<MFAtlasK(A,e);i++)MFKVSetC(s,i,0.,e);MFKVSetC(s,0,1.8*R,e);

    x=MFCloneNVector(u,e);
    rc=MFChartEvaluate(C,s,x,e);
    MFNVSetIndex2(x,-98,e);
    MFIMFSetStability(MFAtlasMF(A,e),x,Phi,e);

    Tan=MFIMFTangentSpaceWithGuess(MFAtlasMF(A,e),x,Phi,e);
    MFFreeChart(C,e);

    rc=MFAtlasAddChart(A,MFCreateBoundaryChart(MFAtlasMF(A,e),x,Tan,1,1,R,e),e);
    {printf("%d) singular half disk off of    branch  ",MFAtlasNumberOfCharts(A,e));MFPrintNVector(stdout,x,e);printf(" index %d, bif? %d\n",MFNVGetIndex(x,e),MFNVGetIndex2(x,e));fflush(stdout);}
    {printf("Tan=\n");MFPrintNKMatrix(stdout,Tan,e);printf("line %d routine %s\n",__LINE__,RoutineName);fflush(stdout);}
    MFFreeNVector(x,e);
#endif

/*  MFNVSetIndex2(u,-99,e);*/

#if 0
/* ---------------------- Half disk back    into the current branch ------------------------------------------- */

    x=MFCloneNVector(u,e);
    MFNVSetIndex2(x,-99,e);   /* This depends on the type of the point, says not to continue this until switch */
    MFNVSetIndex(x,index0,e);

    rc=MFAtlasAddChart(A,MFCreateBoundaryChart(MFAtlasMF(A,e),x,Phi,1,dir,R,e),e);
    {printf("%d) singular half disk continuing  ",rc);MFPrintNVector(stdout,x,e);printf(" index %d, bif? %d\n",MFNVGetIndex(x,e),MFNVGetIndex2(x,e));fflush(stdout);}
    if(0){printf("Tan=\n");MFPrintNKMatrix(stdout,Phi,e);printf("line %d routine %s\n",__LINE__,RoutineName);fflush(stdout);}
    MFFreeNVector(x,e);

/* ---------------------- Half disk forward into the current branch ------------------------------------------- */

    x=MFCloneNVector(u,e);
    MFNVSetIndex2(x,-98,e);   /* This says not to continue this until switch */
    MFNVSetIndex(x,index1,e);

    rc=MFAtlasAddChart(A,MFCreateBoundaryChart(MFAtlasMF(A,e),x,Phi,1,-dir,R,e),e);
    {printf("%d) singular half disk continuing  ",rc);MFPrintNVector(stdout,x,e);printf(" index %d, bif? %d\n",MFNVGetIndex(x,e),MFNVGetIndex2(x,e));fflush(stdout);}
    if(0){printf("Tan=\n");MFPrintNKMatrix(stdout,Phi,e);printf("line %d routine %s\n",__LINE__,RoutineName);fflush(stdout);}
    MFFreeNVector(x,e);

/* ---------------------- Singular disk on the bifurcating branch ------------------------------------------- */
   
    x=MFCloneNVector(u,e);
    MFNVSetIndex2(w,-98,e); 

    MFIMFSetStability(MFAtlasMF(A,e),w,Phi,e);
    bchart=MFAtlasAddChart(A,MFCreateChart(MFAtlasMF(A,e),w,Phi,R,e),e);
    if(bchart>-1){printf("%d)  point tan1 ",MFAtlasNumberOfCharts(A,e));MFPrintNVector(stdout,u,e);printf(" index %d, bif? %d\n",MFNVGetIndex(u,e),MFNVGetIndex2(u,e));fflush(stdout);}
    if(bchart>-1){printf("Tan=\n");MFPrintNKMatrix(stdout,Phi,e);printf("line %d routine %s\n",__LINE__,RoutineName);fflush(stdout);}
    MFFreeNVector(x,e);
   
/* ---------------------- Half disk up into the bifurcating branch ------------------------------------------- */

    
    s=MFCreateKVector(MFAtlasK(A,e),e);
    for(i=1;i<MFAtlasK(A,e);i++)MFKVSetC(s,i,0.,e);MFKVSetC(s,0,1.8*R,e);

    x=MFCloneNVector(u,e);
    MFNVSetIndex2(x,-98,e);

    rc=MFChartEvaluate(MFAtlasChart(A,bchart,e),s,x,e);
    Tan=MFIMFTangentSpaceWithGuess(MFAtlasMF(A,e),x,Phi,e);
    MFIMFSetStability(MFAtlasMF(A,e),x,Tan,e);
    rc=MFAtlasAddChart(A,MFCreateChart(MFAtlasMF(A,e),x,Tan,R,e),e);
    if(rc>-1){printf("%d) A step out the continuing branch ",rc);MFPrintNVector(stdout,x,e);printf(" index %d, bif? %d\n",MFNVGetIndex(x,e),MFNVGetIndex2(x,e));fflush(stdout);}
    if(rc>-1){printf("Tan=\n");MFPrintNKMatrix(stdout,Tan,e);printf("line %d routine %s\n",__LINE__,RoutineName);fflush(stdout);}
/*
    rc=MFAtlasAddChart(A,MFCreateBoundaryChart(MFAtlasMF(A,e),x,Tan,1,-1,R,e),e);
    if(rc>-1){printf("%d) half disk bifurcating branch down ",rc);MFPrintNVector(stdout,x,e);printf(" index %d, bif? %d\n",MFNVGetIndex(x,e),MFNVGetIndex2(x,e));fflush(stdout);}
    if(rc>-1){printf("Tan=\n");MFPrintNKMatrix(stdout,Tan,e);printf("line %d routine %s\n",__LINE__,RoutineName);fflush(stdout);}
 */
    MFFreeNVector(x,e);
    MFFreeKVector(s,e);
    MFFreeNKMatrix(Tan,e);
    
/* ---------------------- Half disk up into the bifurcating branch ------------------------------------------- */
  
    s=MFCreateKVector(MFAtlasK(A,e),e);
    for(i=1;i<MFAtlasK(A,e);i++)MFKVSetC(s,i,0.,e);MFKVSetC(s,0,-1.8*R,e);

    x=MFCloneNVector(u,e);
    MFNVSetIndex2(x,-98,e);

    rc=MFChartEvaluate(MFAtlasChart(A,bchart,e),s,x,e);
    Tan=MFIMFTangentSpaceWithGuess(MFAtlasMF(A,e),x,Phi,e);
    MFIMFSetStability(MFAtlasMF(A,e),x,Tan,e);
    rc=MFAtlasAddChart(A,MFCreateBoundaryChart(MFAtlasMF(A,e),x,Tan,1, 1,R,e),e);
    if(rc>-1){printf("%d) half disk bifurcating branch up ",rc);MFPrintNVector(stdout,x,e);printf(" index %d, bif? %d\n",MFNVGetIndex(x,e),MFNVGetIndex2(x,e));fflush(stdout);}
    if(rc>-1){printf("Tan=\n");MFPrintNKMatrix(stdout,Tan,e);printf("line %d routine %s\n",__LINE__,RoutineName);fflush(stdout);}

    MFFreeNVector(x,e);
    MFFreeKVector(s,e);
    MFFreeNKMatrix(Tan,e);

/* ----------------------------------------------------------------------------------------------------------- */
#endif
#endif
