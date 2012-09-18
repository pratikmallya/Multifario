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

static char *id="@(#) $Id: MFAtlas.c,v 1.16 2011/07/21 17:42:46 mhender Exp $";

#include <MFAtlas.h>
#include <MFErrorHandler.h>
#include <MFNVector.h>
#include <MFKVector.h>
#include <MFNSpace.h>
#include <MFNRegion.h>
#include <MFImplicitMF.h>
#include <MFListOfCharts.h>
#include <MFPolytope.h>
#include <MFPrint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>

#ifdef __cplusplus
 extern "C" {
#endif

double MFNKMatrixDot(MFNSpace,MFNKMatrix,MFNKMatrix,MFErrorHandler);
double MFFullProdEV(int,int,double*,MFErrorHandler);
double MFAtlasDet(int,double*,MFErrorHandler);
int MFPolytopeVertexNumberOfIndices(MFPolytope,int,MFErrorHandler);

#ifndef DBL_QNAN
#define DBL_QNAN 1.e200
#endif

#ifndef DBL_INFINITY
#define DBL_INFINITY 1.e200
#endif

void MFSetError(MFErrorHandler,int,char*,char*,int,char*);
static char MFAtlasErrorMsg[256]="";
void MFSetChartPolytope(MFChart,MFPolytope,MFErrorHandler);
void MFSetChartRadius(MFChart,double,MFErrorHandler);
void MFAtlasSetK(MFAtlas,int,MFErrorHandler);
void MFAtlasSetN(MFAtlas,int,MFErrorHandler);
void MFAtlasPageOutChartsNotNearBoundary(MFAtlas,int,int,char*,MFErrorHandler);
void MFAtlasWriteOutChartsNotPaged(MFAtlas,int,int,char*,MFErrorHandler);
void MFWriteChartToPlotFile(FILE*,MFAtlas,MFChart,int,MFErrorHandler);
int *MFPolytopeVertexIndexSet(MFPolytope,int,MFErrorHandler);
int MFPolytopeIntersectIndexSets(MFPolytope,int,int,int*,MFErrorHandler);

void MFIMFSetSpace(MFImplicitMF,MFNSpace,MFErrorHandler);

#define REALLYPAGE
/* #define TEST 1*/

#ifndef NOBINARYTREE
#include <MFBinaryTree.h>
#endif

/* Internal Routines */

void MFAddInitialHyperplanes(MFAtlas,MFErrorHandler);
int MFAtlasAddChartToList(MFAtlas,MFChart,MFErrorHandler);
MFChart MFAtlasChart(MFAtlas,int,MFErrorHandler);
int MFPolytopeSmallestVertexIndex(MFPolytope,MFErrorHandler);
int MFPolytopeLargestVertexIndex(MFPolytope,MFErrorHandler);
int MFPolytopeVertexLargestIndex(MFPolytope,int,MFErrorHandler);
int MFTestPolytope(MFPolytope,MFErrorHandler);

/* Parameters */

double MFEpsilon=0.2;
double MFRMin=0.001;
double MFDotMin=0.2;
double MFVerbose=1;

/* For Testing */

int MFTestBoundaryPoint(MFAtlas,MFNVector,int,MFErrorHandler);
int MFTestAtlas(MFAtlas,MFErrorHandler);      /* No vertex may be in another PT */
int MFChartTotallyInterior(MFChart,MFKVector,double,MFErrorHandler);

/* For the BoundaryList */

#define MFAtlasIncrementToAllocateBoundaryList 10
void MFAtlasAddChartToBoundaryList(MFAtlas,int,MFErrorHandler);
void MFAtlasRemoveChartFromBoundaryList(MFAtlas,int,MFErrorHandler);
void MFAtlasRemoveChartFromBoundaryListByChart(MFAtlas,int,MFErrorHandler);

/* For the ChartList */

#define MFAtlasIncrementToAllocate 10
int MFAddHalfSpaceToAtlas(MFAtlas,int,int,int,MFKVector,double,MFErrorHandler);

struct MFAtlasSt
 {
  int    n;
  int    k;

  MFImplicitMF M;
  MFNSpace NSpace;

  int    nCharts;
  int    mCharts;
  MFChart *chart;

/* List of Half Spaces */

  int    offset;
  int    nHS;
  int    mHS;
  double *onrm;
  int    nGlobalIndex;
  int    *globalIndex;
  MFKVector *nrm;
  int    *leftChart;
  int    *rightChart;

#ifndef NOBINARYTREE
  MFBinaryTree BTree;
#endif

  FILE *plotFileFid;
  double expFactor;
  FILE *pageFileFid;
  FILE *centerFileFid;
  int   pageFileIndex;
  int   nPagedOut;

#ifndef NOBOUNDARYLIST
  int nBnd;
  int mBnd;
  int *bnd;
#endif

/* For post processing */

  int nclipf;
  int mclipf;
  double (**clipf)(MFNVector,MFErrorHandler);
  int nclipIndx;
  int *clipIndx;

/* For additional culling of nearby points */

  int (*isNear)(MFAtlas,MFChart,MFChart,MFErrorHandler);
 };

void MFPrintChartSummary(MFAtlas A, MFErrorHandler e)
 {
  int i;

  printf("Atlas charts\n");
  for(i=0;i<A->nCharts;i++)
   {
    printf(" Chart %d is 0x%8.8x, ",i,A->chart[i]);
    if(!MFChartPaged(A->chart[i],e))
     {
      printf(" Center is 0x%8.8x, TS is 0x%8.8x",MFChartCenter(A->chart[i],e),MFChartTangentSpace(A->chart[i],e));
      printf(" index %d, bif? %d\n",MFNVGetIndex(MFChartCenter(A->chart[i],e),e),MFNVGetIndex2(MFChartCenter(A->chart[i],e),e));
     }else
      printf(" Paged out\n");
   }
  fflush(stdout);
  return;
 }

MFAtlas MFCreateAtlas(MFImplicitMF M, MFErrorHandler e)
 {
  static char RoutineName[]={"MFCreateAtlas"};
  MFAtlas A;
  int i,j,ij;

#ifdef MFNOCONFIDENCE
  if(M==NULL)
   {
    sprintf(MFAtlasErrorMsg,"Manifold (argument 1) is NULL.");
    MFSetError(e,12,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    return NULL;
   }
#endif

  A=(MFAtlas)malloc(sizeof(struct MFAtlasSt));

#ifndef MFNOSAFETYNET
  if(A==NULL)
   {
    sprintf(MFAtlasErrorMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct MFAtlasSt));
    MFSetError(e,12,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  A->M=M;
  MFRefImplicitMF(M,e);

  A->NSpace=MFIMFNSpace(M,e);
  MFRefNSpace(A->NSpace,e);

  A->n=MFIMF_N(M,e);
  A->k=MFIMF_K(M,e);

/* List of Charts */

  A->nCharts=0;
  A->mCharts=0;
  A->chart=NULL;

/* List of Half spaces */

  A->nHS=0;
  A->mHS=0;
  A->nrm=NULL;
  A->onrm=NULL;
  A->leftChart=NULL;
  A->rightChart=NULL;
  A->nGlobalIndex=0;
  A->globalIndex=NULL;
  A->offset=0;

/* Paging */

  A->pageFileFid=NULL;
  A->pageFileIndex=0;
  A->nPagedOut=0;
  A->plotFileFid=NULL;
  A->expFactor=1.;
  A->centerFileFid=NULL;

/* Half spaces defining a hypercube */

  if(A->k>0)MFAddInitialHyperplanes(A,e);

/* BinaryTree (created when first chart is added) */

#ifndef NOBINARYTREE
  A->BTree=NULL;
#endif

/* List of boundary Charts (created when first chart is added) */

  A->nBnd=0;
  A->mBnd=0;
  A->bnd=NULL;

  A->nclipf=0;
  A->mclipf=0;
  A->clipf=NULL;
  A->nclipIndx=-10;
  A->clipIndx=NULL;

  A->isNear=NULL;

  return A;
 }

void MFFreeAtlas(MFAtlas A, MFErrorHandler e)
 {
  static char RoutineName[]={"MFFreeAtlas"};
  int i;
  int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("MFFreeAtlas\n");fflush(stdout);}
#endif

#ifdef MFNOCONFIDENCE
  if(A==NULL)
   {
    sprintf(MFAtlasErrorMsg,"Pointer to Atlas (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

#ifdef MFALLOWVERBOSE
  if(verbose){printf("  Delete M\n");fflush(stdout);}
#endif

  if(A->M!=NULL)MFFreeImplicitMF(A->M,e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("  Delete NSpace\n");fflush(stdout);}
#endif

  if(A->NSpace!=NULL)MFFreeNSpace(A->NSpace,e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("  Delete List of charts\n");fflush(stdout);}
#endif

  if(A->chart!=NULL)
   {
    for(i=0;i<A->nCharts;i++)
     {
      if(A->chart[i]!=NULL)MFFreeChart(A->chart[i],e);
     }
    free(A->chart);
   }

/* List of half spaces */

#ifdef MFALLOWVERBOSE
  if(verbose){printf("  Delete List of onrm\n");fflush(stdout);}
#endif

  if(A->onrm!=NULL)free(A->onrm);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("  Delete List of globalIndex\n");fflush(stdout);}
#endif

  if(A->globalIndex!=NULL)free(A->globalIndex);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("  Delete List of leftChart\n");fflush(stdout);}
#endif
  if(A->leftChart!=NULL)free(A->leftChart);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("  Delete List of rightChart\n");fflush(stdout);}
#endif
  if(A->rightChart!=NULL)free(A->rightChart);
  if(A->nrm!=NULL)
   {

#ifdef MFALLOWVERBOSE
    if(verbose){printf("  Delete List of nrm\n");fflush(stdout);}
#endif

    if(verbose){printf("    length %d\n",A->nHS);fflush(stdout);}
    for(i=0;i<A->nHS;i++)
     {
      if(verbose){printf("    Half space %d\n",i);fflush(stdout);}
      if(A->nrm[i]!=NULL)MFFreeKVector(A->nrm[i],e);
     }
    if(verbose){printf("    Delete A->nrm\n");fflush(stdout);}
    free(A->nrm);
   }

#ifndef NOBINARYTREE

#ifdef MFALLOWVERBOSE
  if(verbose){printf("  Delete BTree\n");fflush(stdout);}
#endif /* ALLOWVERBOSE */

  if(A->BTree!=NULL)MFFreeBinaryTree(A->BTree,e);
#endif  /* NOBINARYTREE */

  if(verbose){printf("  Delete boundary list\n");fflush(stdout);}

  if(A->bnd!=NULL)free(A->bnd);

  if(verbose){printf("  Delete List of clipping indices\n");fflush(stdout);}
  if(A->clipIndx!=NULL)free(A->clipIndx);
  if(verbose){printf("  Delete List of clipping functions\n");fflush(stdout);}
  if(A->clipf!=NULL)free(A->clipf);

  if(verbose){printf("  Delete page file id\n");fflush(stdout);}
  if(A->pageFileFid!=NULL)fclose(A->pageFileFid);

  if(verbose){printf("  Delete A\n");fflush(stdout);}
  free(A);

#ifdef MFALLOWVERBOSE
  if(verbose)printf("done MFFreeAtlas\n");
#endif

  return;
 }

int MFAtlasK(MFAtlas A, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAtlasK"};

#ifdef MFNOCONFIDENCE
  if(A==NULL)
   {
    sprintf(MFAtlasErrorMsg,"Pointer to Atlas (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    return -1;
   }
#endif

  return A->k;
 }

int MFAtlasN(MFAtlas A, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAtlasN"};

#ifdef MFNOCONFIDENCE
  if(A==NULL)
   {
    sprintf(MFAtlasErrorMsg,"Pointer to Atlas (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    return -1;
   }
#endif

  return A->n;
 }

int MFAtlasNumberOfCharts(MFAtlas A, MFErrorHandler e)
 {
  static char RoutineName[]={"MFNumberOfCharts"};

#ifdef MFNOCONFIDENCE
  if(A==NULL)
   {
    sprintf(MFAtlasErrorMsg,"Pointer to Atlas (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    return -1;
   }
#endif

  return A->nCharts;
 }

double MFAtlasChartRadius(MFAtlas A,int chart, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAtlasChartRadius"};

#ifdef MFNOCONFIDENCE
  if(A==NULL)
   {
    sprintf(MFAtlasErrorMsg,"Pointer to Atlas (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    return -1.;
   }

  if(chart<0||chart>=A->nCharts)
   {
    sprintf(MFAtlasErrorMsg,"Chart %d (argument 2) is invalid. Must be positive and less than %d (total number of charts).",chart,A->nCharts);
    MFSetError(e,4,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    return -1.;
   }
#endif

  return MFChartRadius(A->chart[chart],e);
 }

MFNVector MFAtlasCenterOfChart(MFAtlas A,int chart, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAtlasCenterOfChart"};

#ifdef MFNOCONFIDENCE
  if(A==NULL)
   {
    sprintf(MFAtlasErrorMsg,"Pointer to Atlas (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    return NULL;
   }

  if(chart<0||chart>=A->nCharts)
   {
    sprintf(MFAtlasErrorMsg,"Chart %d (argument 2) is invalid. Must be positive and less than %d (total number of charts).",chart,A->nCharts);
    MFSetError(e,4,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    return NULL;
   }
#endif

  return MFChartCenter(A->chart[chart],e);
 }

int MFAtlasIsPointInChart(MFAtlas A,int chart, MFKVector pt, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAtlasIsPointInChart"};

#ifdef MFNOCONFIDENCE
  if(A==NULL)
   {
    sprintf(MFAtlasErrorMsg,"Pointer to Atlas (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    return -1;
   }

  if(chart<0||chart>=A->nCharts)
   {
    sprintf(MFAtlasErrorMsg,"Chart %d (argument 2) is invalid. Must be positive and less than %d (total number of charts).",chart,A->nCharts);
    MFSetError(e,4,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    return -1;
   }

  if(pt==NULL)
   {
    sprintf(MFAtlasErrorMsg,"Point (argument 3) is NULL.");
    MFSetError(e,4,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    return -1;
   }
#endif

  return MFChartInterior(A->chart[chart],pt,e);
 }

void MFAtlasEvaluateChart(MFAtlas A,int chart,MFKVector s,MFNVector u, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAtlasEvaluateChart"};

#ifdef MFNOCONFIDENCE
  if(A==NULL)
   {
    sprintf(MFAtlasErrorMsg,"Pointer to Atlas (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(chart<0||chart>=A->nCharts)
   {
    sprintf(MFAtlasErrorMsg,"Chart %d (argument 2) is invalid. Must be positive and less than %d (total number of charts).",chart,A->nCharts);
    MFSetError(e,4,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(s==NULL)
   {
    sprintf(MFAtlasErrorMsg,"Point in domain (argument 3) is NULL.");
    MFSetError(e,4,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(u==NULL)
   {
    sprintf(MFAtlasErrorMsg,"Point in range (argument 4) is NULL.");
    MFSetError(e,4,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  MFChartEvaluate(A->chart[chart],s,u,e);

  return;
 }

int MFAtlasNumberOfChartsWithBoundary(MFAtlas A, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAtlasNumberOfChartsWithBoundary"};

#ifdef MFNOCONFIDENCE
  if(A==NULL)
   {
    sprintf(MFAtlasErrorMsg,"Pointer to Atlas (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    return -1;
   }
#endif

  return A->nBnd;
 }

int MFAtlasChartWithBoundary(MFAtlas A,int i, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAtlasChartWithBoundary"};

#ifdef MFNOCONFIDENCE
  if(A==NULL)
   {
    sprintf(MFAtlasErrorMsg,"Pointer to Atlas (argument 1) is NULL");
    printf("%s : %s\n",RoutineName,MFAtlasErrorMsg);fflush(stdout);
    MFSetError(e,4,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    return -1;
   }

  if(i<0||i>=A->nBnd)
   {
    sprintf(MFAtlasErrorMsg,"Index of boundary chart %d (argument 2) is invalid. Must be positive and less than %d (number of charts on boundary)",i,A->nBnd);
    printf("%s : %s\n",RoutineName,MFAtlasErrorMsg);fflush(stdout);
    MFSetError(e,4,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    return -1;
   }
#endif

  return A->bnd[i];
 }

int MFAtlasAddChartWithCenter(MFAtlas A,MFNVector u, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAtlasAddChartWithCenter"};
  int result;
  double R;
  MFNKMatrix Phi;

  Phi=MFIMFTangentSpace(A->M,u,e);
  if(Phi==NULL)return 0;
  R=MFIMFScale(MFAtlasMF(A,e),u,Phi,e);
  result=MFAtlasAddChartWithAll(A,u,Phi,R,e);
  MFFreeNKMatrix(Phi,e);

  return result;
 }

/* Add radius of chart approx came from and distance 
   between u and the chart */
 
int MFAtlasAddChartWithApproxTS(MFAtlas A,MFNVector u,MFNKMatrix TS,double delta, double R0, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAtlasAddChartWithApproxTS"};
  double R;
  MFNKMatrix Phi;
  int chart;
  int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s\n",RoutineName);fflush(stdout);}
#endif

#ifdef MFNOCONFIDENCE
  if(A==NULL)
   {
    sprintf(MFAtlasErrorMsg,"Pointer to Atlas (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    return -1;
   }

  if(u==NULL)
   {
    sprintf(MFAtlasErrorMsg,"Pointer to u (argument 2) is NULL");
    MFSetError(e,4,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    return -1;
   }
#endif

#ifdef MFALLOWVERBOSE
  if(verbose){printf("       calculate TS \n");fflush(stdout);}
#endif

  if(TS==NULL)
   {
#ifdef MFALLOWVERBOSE
    if(verbose){printf("         without guess\n");fflush(stdout);}
#endif

    Phi=MFIMFTangentSpace(A->M,u,e);
   }else{
#ifdef MFALLOWVERBOSE
    if(verbose){printf("         with guess\n");fflush(stdout);}
#endif

    Phi=MFIMFTangentSpaceWithGuess(A->M,u,TS,e);
   }
  if(Phi==NULL)return -1;

  R=MFIMFScale(MFAtlasMF(A,e),u,Phi,e);
  chart=MFAtlasAddChartWithAll(A,u,Phi,R,e);
  MFFreeNKMatrix(Phi,e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf(" done %s\n",RoutineName);fflush(stdout);}
#endif

  return chart;
 }

/* --------------------------------------------------------- */

int MFAtlasPointOnBoundaryInsideRegion(MFAtlas A,MFNRegion Omega,MFNVector u, MFNKMatrix *TS, double *dist, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAtlasPointOnBoundaryInsideRegion"};
  int v,l,n;
  int pass;
  int nChartsTried;
  int chart;
  int bchart;
  MFNVector up;
  MFKVector s;
  MFKVector s2;
  MFNVector w;
  MFNKMatrix Phi;
  double dot;
  double rv;
  double R;
  double delta;
  int rc;
  int abandon;
  double t,tl,tr;
  int bisection;
/*double tol=1.e-10;*/
  double tol=1.e-4;
  int verbose=0;

  int index0,index1;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s\n",RoutineName);fflush(stdout);}
#endif

  up=MFCloneNVector(u,e);
  s=MFCreateKVector(A->k,e);
  nChartsTried=0;
  bisection=0;

  *TS=NULL;

  pass=1;
  bchart=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s, there are %d charts on the boundary list\n",RoutineName,A->nBnd);fflush(stdout);}
  if(verbose){printf("bchart<A->nBnd=%d pass<2 = %d\n",bchart<A->nBnd,pass<2);fflush(stdout);}
#endif /* ALLOWVERBOSE */

  while(bchart<A->nBnd && pass<2)
   {
    chart=A->bnd[bchart];

    nChartsTried++;
    if(MFChartIsSingular(A->chart[chart],e))goto nextchart;

tryagain:
    n=MFPolytopeNumberOfVertices(MFChartPolytope(A->chart[chart],e),e);

#ifdef MFALLOWVERBOSE
    if(verbose)
     {
      printf(" chart %d has %d vertices, R=%lf, MinR=%lf\n",chart,n,MFChartRadius(A->chart[chart],e),MFIMFGetRMin(A->M,e));fflush(stdout);
      for(v=0;v<n;v++)
      printf("   Vertex %d has mark %d, index %d, radius %le\n",v,
           MFPolytopeGetVertexMark(MFChartPolytope(A->chart[chart],e),v,e),
           MFNVGetIndex2(MFChartCenter(A->chart[chart],e),e),
           MFPolytopeRadiusOfVertex(MFChartPolytope(A->chart[chart],e),v,e));
     }
#endif

    v=0;
    while(v<n)
     {
      if(MFChartRadius(A->chart[chart],e)<MFIMFGetRMin(A->M,e))
       {
#ifdef MFALLOWVERBOSE
        if(verbose){printf("   Radius of vertex too small %le < %le\n",MFChartRadius(A->chart[chart],e),MFIMFGetRMin(A->M,e));fflush(stdout);}
#endif
        goto nextvertex;
       }

       if(MFPolytopeGetVertexMark(MFChartPolytope(A->chart[chart],e),v,e)!=0)
       {
#ifdef MFALLOWVERBOSE
        if(verbose){printf("   Vertex mark not zero, ignoring the vertex.\n");fflush(stdout);}
#endif
        goto nextvertex;
       }

       if(MFNVGetIndex2(MFChartCenter(A->chart[chart],e),e)==-98)
       {
#ifdef MFALLOWVERBOSE
        if(verbose){printf("   Chart tagged as -98, a type of singular chart. Skipping the chart.\n");fflush(stdout);}
#endif
        goto nextchart;
       }

       if(MFNVGetIndex2(MFChartCenter(A->chart[chart],e),e)==-99)
       {
#ifdef MFALLOWVERBOSE
        if(verbose){printf("   Chart tagged as -99, a type of singular chart. Skipping the chart.\n");fflush(stdout);}
#endif
        goto nextchart;
       }


#ifdef MFALLOWVERBOSE
      if(verbose){printf(" vertex %d \n",v);fflush(stdout);}
#endif
      if(pass==1 || ( pass==0 && MFPolytopeVertexLargestIndex(MFChartPolytope(A->chart[chart],e),v,e) >= A->offset))
       {
        rv=MFPolytopeRadiusOfVertex(MFChartPolytope(A->chart[chart],e),v,e);
        R=MFChartRadius(A->chart[chart],e);

#ifdef MFALLOWVERBOSE
        if(verbose){printf("  Trying vertex %d of chart %d (Radius %le)\n",v,chart,R);fflush(stdout);}
#endif
        if(rv>=R)
         {
#ifdef MFALLOWVERBOSE
          if(verbose){printf("  Exterior vertex. ");fflush(stdout);}
#endif
          t=1.;
bsect:

#ifdef MFALLOWVERBOSE
          if(verbose){printf("t=%le, bisection=%d\n",t,bisection);fflush(stdout);}
#endif

          if(MFIMFGetRMin(A->M,e)>0.&&t*R<MFIMFGetRMin(A->M,e)&&!bisection){printf("R*t too small, aborting vertex.R=%le, t=%le, RMin=%le\n",R,t,MFIMFGetRMin(A->M,e));fflush(stdout);goto nextvertex;}
          if(MFIMFGetRMin(A->M,e)<=0.&&t<0.001&&!bisection){printf("t too small, aborting vertex\n");fflush(stdout);goto nextvertex;}
          MFPolytopeVertex(MFChartPolytope(A->chart[chart],e),v,s,e);
          MFKVScale(t*R/rv,s,e);
          MFChartPointInTangentSpace((A->chart[chart]),s,up,e);

#ifdef MFALLOWVERBOSE
          if(verbose){printf("up(");MFPrintKVector(stdout,s,e);printf(")=");MFPrintNVector(stdout,up,e);printf("\n");fflush(stdout);}
#endif
          if(!MFNRegionInterior(Omega,up,e))
           {

/* Instead invoke bisection, state=looking for boundary crossing */

/* Omega should be a n-dimensional manifold complex embedded in n-space */

/* At a boundary crossing need to get "active" boundaries and the "cone" of tangents */

#ifdef MFALLOWVERBOSE
            if(verbose){printf("     Not in Omega.\n");fflush(stdout);}
#endif
            MFPolytopeSetVertexMark(MFChartPolytope(A->chart[chart],e),v,1,e);
            goto nextvertex;
           }
#ifdef MFALLOWVERBOSE
          if(verbose){printf("     In Omega.\n");fflush(stdout);}
#endif

          rc=MFChartEvaluate(A->chart[chart],s,u,e);
          if(!rc)
           {

#ifdef MFALLOWVERBOSE
            if(verbose)
             {
              printf("%s - ChartEvaluate failed at vertex %d of chart %d, chart R=%lf vertex R %lf\n",RoutineName,chart,v,R,rv);fflush(stdout);
              printf("        rc=%d,  R*t=%lf, RMin=%lf\n",rc,t*R,MFIMFGetRMin(A->M,e));fflush(stdout);
             }
#endif

            if(rc<0)goto nextvertex;

            if(!bisection)
             {
              t=.8*t;
              if(MFIMFGetRMin(A->M,e)>0.&&t*R<MFIMFGetRMin(A->M,e)){printf("R*t too small, aborting vertex\n");fflush(stdout);goto nextvertex;}
              if(MFIMFGetRMin(A->M,e)<=0.&&t<0.001){printf("t too small, aborting vertex\n");fflush(stdout);goto nextvertex;}
              MFSetSuggestedChartRadius(A->chart[chart],t*R,e);
              goto bsect;
             }else{
              tr=t+.1*tol;
              tl=t-.1*tol;
              goto donebsect;
             }
           }

/* Was correction orthogonal? */

#ifdef MFALLOWVERBOSE
          if(verbose)
           {
            s2=MFCreateKVector(A->k,e);
            MFChartProjectIntoTangentSpace(A->chart[chart],up,s2,e);
            printf("\n     s(up) is ");MFPrintKVector(stdout,s2,e);printf(" should be the same as s ");MFPrintKVector(stdout,s,e);printf("\n");fflush(stdout);
            MFChartProjectIntoTangentSpace(A->chart[chart],u,s2,e);
            printf("     s(u ) is ");MFPrintKVector(stdout,s2,e);printf(" should be the same as s ");MFPrintKVector(stdout,s,e);printf("\n\n");fflush(stdout);
            MFFreeKVector(s2,e);
           }
#endif
 
          delta=MFNSpaceDistance(A->NSpace,u,up,e);

#ifdef MFALLOWVERBOSE
          if(verbose)
           {
            printf("     The correction was %lf, epsilon=%lf\n",delta,MFEpsilon);fflush(stdout);
            printf("     u[%d] is ",chart);MFPrintNVector(stdout,MFChartCenter(A->chart[chart],e),e);printf("\n");fflush(stdout);
            printf("     s is ");MFPrintKVector(stdout,s,e);printf("\n");fflush(stdout);
            printf("     up is ");MFPrintNVector(stdout,up,e);printf("\n");fflush(stdout);
            printf("     u  is ");MFPrintNVector(stdout,u,e);printf("\n");fflush(stdout);
/*          printf("     Phi[%d] is 0x%8.8x\n",chart,MFChartTangentSpace(A->chart[chart],e));fflush(stdout);
 *          MFPrintNKMatrix(stdout,MFChartTangentSpace(A->chart[chart],e),e);printf("\n");fflush(stdout);
 */
           }
#endif

          if(delta>MFEpsilon&&!bisection)
           {

#ifdef MFALLOWVERBOSE
            if(verbose){printf("%s - Correction too large at R=%lf, (%lf>%lf) trying again at %lf. Radius of vertex is %lf\n",RoutineName,t*R,delta,MFEpsilon,.8*t*R,rv);fflush(stdout);}
#endif

            t=.8*t;
/*          if(t<.001)goto nextvertex;*/
            MFSetSuggestedChartRadius(A->chart[chart],t*R,e);
            goto bsect;
           }else if(bisection)
           {

#ifdef MFALLOWVERBOSE
            if(verbose){printf("     Bisection: interval [%lf, %lf, %lf ]\n",tl,t,tr);fflush(stdout);}
#endif

            if(*TS!=NULL){MFFreeNKMatrix(*TS,e);*TS=NULL;}
            *TS=MFIMFTangentSpaceWithGuess(A->M,u,MFChartTangentSpace(A->chart[chart],e),e);
            if(*TS==NULL)
             {
              bisection=0;
              MFPolytopeSetVertexMark(MFChartPolytope(A->chart[chart],e),v,2,e);
              printf("*** MFIMFTangentSpace failed at bisection point, exempting vertex %d of chart %d\n",v,chart);fflush(stdout);
              goto nextvertex;
             }
            if(delta<MFEpsilon&&!MFIMFStop(A->M,MFChartCenter(A->chart[chart],e),MFChartTangentSpace(A->chart[chart],e),u,*TS,e))
             {

#ifdef MFALLOWVERBOSE
              if(verbose)printf("                 Keep Right\n");
#endif
              tl=t;
             }else{

#ifdef MFALLOWVERBOSE
              if(verbose)printf("                 Keep Left\n");
#endif
              tr=t;
             }
            t=(tr+tl)/2;

#ifdef MFALLOWVERBOSE
            if(verbose){printf("            new interval [%lf, %lf, %lf ] width %le \n",tl,t,tr,tr-tl);fflush(stdout);}
#endif

donebsect:
            if(tr-tl<tol)
             {
              if(fabs(t)>1.e-6)
               {

#ifdef MFALLOWVERBOSE
                if(verbose){printf("     Interval less than tol %le (%le), declaring success, t in [%lf,%lf]\n",tol,tr-tl,tl,tr);}
#endif
                MFPolytopeSetVertexMark(MFChartPolytope(A->chart[chart],e),v,3,e);

/* If a crossing of Omega's boundary: */

/* else If a crossing of a singular "boundary": */

                MFNVSetIndex(u,MFNVGetIndex(MFChartCenter(A->chart[chart],e),e),e);

/* Do this if it is a bifurcation point (not limit point or special point) */

                if(*TS!=NULL){MFFreeNKMatrix(*TS,e);*TS=NULL;}
                *TS=MFCloneNKMatrix(MFChartTangentSpace(A->chart[chart],e),e);

                w=MFCloneNVector(up,e);
                rc=MFIMFSingular(A->M,u,*TS,w,e); /* w is the "other" null vector */

#ifdef MFALLOWVERBOSE
                if(verbose){printf("     return code from MFIMFSingular %d\n",rc);fflush(stdout);}
#endif

                if(rc==1)
                 {
                  if(1||verbose){printf("%s      Singular Point Detected and located.\n",RoutineName);fflush(stdout);}
                  MFNVSetIndex2(u,-100,e);

                  MFPolytopeSetVertexMark(MFChartPolytope(A->chart[chart],e),v,1,e);
                  MFFreeNVector(up,e);
                  return chart;
#if 0

/* Then this is a "like--type" bifurcation */

                  double vnorm;

/* Compute V, the normal to the co--dimension 1 singular "boundary manifold" */

                  MFNVector V;
                  MFKVector s;
                  MFNKMatrix Psi;

                  V=MFCloneNVector(w,e);
                  s=MFCreateKVector(MFAtlasK(A,e),e);

                  MFNSpaceScale(A->NSpace,.5*MFChartRadius(A->chart[chart],e),w,up,e);
                  MFNSpaceAdd(A->NSpace,u,up,up,e);
                  Psi=MFIMFTangentSpace(A->M,up,e);
                  if(Psi==NULL)
                   {
                    bisection=0;
                    MFPolytopeSetVertexMark(MFChartPolytope(A->chart[chart],e),v,4,e);
                    printf("*** MFIMFTangentSpace failed at predicted point, exempting vertex %d of chart %d\n",v,chart);fflush(stdout);
                    goto nextvertex;
                   }
                  Phi=MFCloneNKMatrix(MFChartTangentSpace(A->chart[chart],e),e);

                  MFMVMulT(A->NSpace,Psi,w,s,e);        /* V= Psi Psi^T phi_k */
                  MFMVMul (A->NSpace,Psi,s,V,e);

                  MFMVMulT(A->NSpace,Phi,V,s,e);        /* V= Phi Phi^T Psi Psi^T phi_k */
                  MFMVMul (A->NSpace,Phi,s,V,e);

/* done computing V */

                  MFFreeNVector(w,e);
                  MFFreeKVector(s,e);

                  w=MFCloneNVector(u,e);
                  vnorm=sqrt(MFNSpaceInner(A->NSpace,V,V,e));

/* Create BoundaryChart creates a half disk */

                  if(fabs(vnorm) < 1.e-10)
                   {
                    double R;
                    int bchart;
                    int i;
                    MFNKMatrix Tan;
                    MFNVector x;
                    x=MFCloneNVector(u,e);
                    MFNVSetIndex2(x,-99,e);   /* This depends on the type of the point */

                    R=1.*MFChartRadius(A->chart[chart],e);
                    {printf("%d) Add chart to Atlas ",A->nCharts);MFPrintNVector(stdout,u,e);printf(" index %d, bif? %d\n",MFNVGetIndex(u,e),MFNVGetIndex2(u,e));fflush(stdout);}
                    rc=MFAtlasAddChart(A,MFCreateChart(A->M,x,Phi,R,e),e);
                    MFFreeNVector(x,e);

                    MFIMFSetStability(A->M,w,Psi,e);
                    MFNVSetIndex2(w,-98,e);
                    bchart=MFAtlasAddChartWithAll(A,w,Psi,R,e);
   
                    for(i=1;i<A->k;i++)MFKVSetC(s,i,0.,e);MFKVSetC(s,0,2.8*R,e);
   
                    x=MFCloneNVector(u,e);
                    MFNVSetIndex2(x,-98,e);
                    rc=MFChartEvaluate(A->chart[bchart],s,x,e);
                    Tan=MFIMFTangentSpaceWithGuess(A->M,x,Psi,e);
                    MFIMFSetStability(A->M,x,Tan,e);
                    rc=MFAtlasAddChart(A,MFCreateBoundaryChart(A->M,x,Tan,1, 1,R,e),e);
                    rc=MFAtlasAddChart(A,MFCreateBoundaryChart(A->M,x,Tan,1,-1,R,e),e);
   
                    MFFreeNVector(x,e);
                    MFFreeNKMatrix(Tan,e);
   
                    if(rc>-1){printf("%d) Bifurcating branch          ",rc);MFPrintNVector(stdout,x,e);printf(" index %d, bif? %d\n",MFNVGetIndex(x,e),MFNVGetIndex2(x,e));fflush(stdout);}
   
                    MFKVSetC(s,0,-1.8*R,e);
                    x=MFCloneNVector(u,e);
                    MFNVSetIndex2(x,-98,e);
                    rc=MFChartEvaluate(A->chart[bchart],s,x,e);
                    Tan=MFIMFTangentSpaceWithGuess(A->M,x,Psi,e);
                    MFIMFSetStability(A->M,x,Tan,e);
                    rc=MFAtlasAddChart(A,MFCreateBoundaryChart(A->M,x,Tan,1, 1,R,e),e);
                    rc=MFAtlasAddChart(A,MFCreateBoundaryChart(A->M,x,Tan,1,-1,R,e),e);
   
                    MFFreeNKMatrix(Tan,e);
   
                    if(rc>-1){printf("%d) Bifurcating branch          ",rc);MFPrintNVector(stdout,x,e);printf(" index %d, bif? %d\n",MFNVGetIndex(x,e),MFNVGetIndex2(x,e));fflush(stdout);}
   
                   }else{


#if 0
                    double R;
                    int bchart;
                    int i;
                    MFNKMatrix Tan;
                    MFNVector x,y;
                    int dir;
                    MFChart C;

                    if(1||verbose){printf("      Singular Point Detected and located. |vnorm| > 1.e-10\n");fflush(stdout);}

                    R=1.*MFChartRadius(A->chart[chart],e);

          /* V is the normal to the co--dimension 1 singular "boundary manifold" */

                    MFNSpaceScale(A->NSpace,1./vnorm,V,V,e);
                    MFGramSchmidtReplace(A->NSpace,Phi,V,V,e);         /* make the singular curve normal the first basis vector */

          /* Phi is the estimated tangent of the bifurcating branch */

          /* y is the vector from the center of the original chart to the singular point */

                    y=MFCloneNVector(u,e);
                    R=1.*MFChartRadius(A->chart[chart],e);
                    MFNSpaceDirection(A->NSpace,u,MFChartCenter(A->chart[chart],e),y,e);
                    dir=1;
                    if(MFNSpaceInner(A->NSpace,V,y,e)<0)dir=-1;
                    MFFreeNVector(y,e);
#endif

/* ---------------------- Half disk forward into the current branch ------------------------------------------- */
#if 0
                    x=MFCloneNVector(u,e);
                    MFNVSetIndex2(x,-99,e);   /* This depends on the type of the point, says not to continue this until switch */

                    rc=MFAtlasAddChart(A,MFCreateBoundaryChart(A->M,x,Phi,1,-dir,R,e),e);
                    {printf("%d) singular half disk back into branch  ",rc);MFPrintNVector(stdout,x,e);printf(" index %d, bif? %d\n",MFNVGetIndex(x,e),MFNVGetIndex2(x,e));fflush(stdout);}
                    MFFreeNVector(x,e);
#endif

#if 0
/* ---------------------- Original code ----------------------------------------------------------------------- */

/* ---------------------- Full disk at the singular point ----------------------------------------------------- */
  
                    MFNSpaceScale(A->NSpace,1./vnorm,V,V,e);
                    MFGramSchmidtReplace(A->NSpace,Phi,V,w,e);
   
                    R=1.*MFChartRadius(A->chart[chart],e);
                    x=MFCloneNVector(u,e);
   
                    MFIMFSetStability(A->M,w,Phi,e);
                    MFNVSetIndex2(w,-98,e);
                    C=MFCreateChart(A->M,w,Phi,R,e);
                    rc=MFAtlasAddChart(A,C,e);
                    {printf("%d) Add boundary chart to Atlas ",rc);MFPrintNVector(stdout,u,e);printf(" index %d, bif? %d\n",MFNVGetIndex(u,e),MFNVGetIndex2(u,e));fflush(stdout);}
   
                    MFFreeChart(C,e);


/* ---------------------- Half disks on the current branch ---------------------------------------------------- */
   
                    s=MFCreateKVector(MFAtlasK(A,e),e);
                    for(i=1;i<A->k;i++)MFKVSetC(s,i,0.,e);MFKVSetC(s,0,1.8*R,e);
   
                    x=MFCloneNVector(u,e);
                    MFNVSetIndex2(x,-98,e);
                    rc=MFChartEvaluate(A->chart[bchart],s,x,e);
                    Tan=MFIMFTangentSpaceWithGuess(A->M,x,Phi,e);
                    MFIMFSetStability(A->M,x,Tan,e);

                    rc=MFAtlasAddChart(A,MFCreateBoundaryChart(A->M,x,Tan,1, 1,R,e),e);
                    if(rc>-1){printf("%d) Bifurcating branch          ",rc);MFPrintNVector(stdout,x,e);printf(" index %d, bif? %d\n",MFNVGetIndex(x,e),MFNVGetIndex2(x,e));fflush(stdout);}

                    rc=MFAtlasAddChart(A,MFCreateBoundaryChart(A->M,x,Tan,1,-1,R,e),e);
                    if(rc>-1){printf("%d) Bifurcating branch          ",rc);MFPrintNVector(stdout,x,e);printf(" index %d, bif? %d\n",MFNVGetIndex(x,e),MFNVGetIndex2(x,e));fflush(stdout);}

                    MFFreeNVector(x,e);
                    MFFreeNKMatrix(Tan,e);
                    MFFreeKVector(s,e);

/* ---------------------- Half disks on the biufurcating branch ----------------------------------------------- */
   
   
                    s=MFCreateKVector(MFAtlasK(A,e),e);
                    for(i=1;i<A->k;i++)MFKVSetC(s,i,0.,e);MFKVSetC(s,0,-1.8*R,e);
                    x=MFCloneNVector(u,e);
                    MFNVSetIndex2(x,-98,e);

                    rc=MFChartEvaluate(A->chart[bchart],s,x,e);
                    Tan=MFIMFTangentSpaceWithGuess(A->M,x,Phi,e);
                    MFIMFSetStability(A->M,x,Tan,e);

                    rc=MFAtlasAddChart(A,MFCreateBoundaryChart(A->M,x,Tan,1, 1,R,e),e);
                    if(rc>-1){printf("%d) Bifurcating branch          ",rc);MFPrintNVector(stdout,x,e);printf(" index %d, bif? %d\n",MFNVGetIndex(x,e),MFNVGetIndex2(x,e));fflush(stdout);}

                    rc=MFAtlasAddChart(A,MFCreateBoundaryChart(A->M,x,Tan,1,-1,R,e),e);
                    if(rc>-1){printf("%d) Bifurcating branch          ",rc);MFPrintNVector(stdout,x,e);printf(" index %d, bif? %d\n",MFNVGetIndex(x,e),MFNVGetIndex2(x,e));fflush(stdout);}
   
                    MFFreeNVector(x,e);
                    MFFreeNKMatrix(Tan,e);
                    MFFreeKVector(s,e);
   
#endif


#if 0
/* This computes the bifurcating branches. It's been disabled because it should be put into it's own cell. */

  
                    MFNSpaceScale(A->NSpace,1./vnorm,V,V,e);
                    MFGramSchmidtReplace(A->NSpace,Phi,V,w,e);
   
                    R=1.*MFChartRadius(A->chart[chart],e);
                    x=MFCloneNVector(u,e);
   
                    MFIMFSetStability(A->M,x,Phi,e);
                    MFNVSetIndex2(x,-98,e);
                    C=MFCreateChart(A->M,x,Phi,R,e);
                    MFFreeChart(x,e);

                    s=MFCreateKVector(MFAtlasK(A,e),e);
                    for(i=1;i<A->k;i++)MFKVSetC(s,i,0.,e);MFKVSetC(s,0,1.8*R,e);
   
                    x=MFCloneNVector(u,e);
                    MFNVSetIndex2(x,-98,e);
                    rc=MFChartEvaluate(C,s,x,e);
                    Tan=MFIMFTangentSpaceWithGuess(A->M,x,Phi,e);
                    MFIMFSetStability(A->M,x,Tan,e);
                    rc=MFAtlasAddChart(A,MFCreateChart(A->M,x,Tan,R,e),e);
                    if(rc>-1){printf("%d) Bifurcating branch          ",rc);MFPrintNVector(stdout,x,e);printf(" index %d, bif? %d\n",MFNVGetIndex(x,e),MFNVGetIndex2(x,e));fflush(stdout);}

                    MFFreeNVector(x,e);
                    MFFreeNKMatrix(Tan,e);
   
                    MFKVSetC(s,0,-1.8*R,e);
                    x=MFCloneNVector(u,e);
                    MFNVSetIndex2(x,-98,e);
                    rc=MFChartEvaluate(C,s,x,e);
                    Tan=MFIMFTangentSpaceWithGuess(A->M,x,Phi,e);
                    MFIMFSetStability(A->M,x,Tan,e);
                    rc=MFAtlasAddChart(A,MFCreateChart(A->M,x,Tan,R,e),e);
                    if(rc>-1){printf("%d) Bifurcating branch          ",rc);MFPrintNVector(stdout,x,e);printf(" index %d, bif? %d\n",MFNVGetIndex(x,e),MFNVGetIndex2(V,e));fflush(stdout);}

                    MFFreeNVector(x,e);
                    MFFreeNKMatrix(Tan,e);
                    MFFreeKVector(s,e);
                    MFFreeChart(C,e);
#endif
   

#if 0

/* -----------------------Half Disk on the new branch ------------------------------------------- */

                    MFGramSchmidtReplace(A->NSpace,Phi,V,w,e);         /* replace the singular curve tangent with the first basis vector */
  
                    x=MFCloneNVector(u,e);
                    MFNSpaceScale(A->NSpace,1.,u,x,e);
                    C=MFCreateChart(A->M,x,Phi,.2*R,e);
                    MFFreeNVector(x,e); 

                    s=MFCreateKVector(MFAtlasK(A,e),e);
                    for(i=1;i<A->k;i++)MFKVSetC(s,i,0.,e);MFKVSetC(s,0,1.8*R,e);

                    x=MFCloneNVector(u,e);
                    rc=MFChartEvaluate(C,s,x,e);
                    MFNVSetIndex2(x,-98,e);
                    MFIMFSetStability(A->M,x,Phi,e);

                    Tan=MFIMFTangentSpaceWithGuess(A->M,x,Phi,e);
                    MFFreeChart(C,e);

                    rc=MFAtlasAddChart(A,MFCreateBoundaryChart(A->M,x,Tan,1,1,R,e),e);
                    {printf("%d) singular half disk off of    branch  ",A->nCharts);MFPrintNVector(stdout,x,e);printf(" index %d, bif? %d\n",MFNVGetIndex(x,e),MFNVGetIndex2(x,e));fflush(stdout);}
                    {printf("Tan=\n");MFPrintNKMatrix(stdout,Tan,e);printf("line %d routine %s\n",__LINE__,RoutineName);fflush(stdout);}
                    MFFreeNVector(x,e);
#endif

/*                  MFNVSetIndex2(u,-99,e);*/

#if 0
/* ---------------------- Half disk back    into the current branch ------------------------------------------- */

                    x=MFCloneNVector(u,e);
                    MFNVSetIndex2(x,-99,e);   /* This depends on the type of the point, says not to continue this until switch */
                    MFNVSetIndex(x,index0,e);

                    rc=MFAtlasAddChart(A,MFCreateBoundaryChart(A->M,x,Phi,1,dir,R,e),e);
                    {printf("%d) singular half disk continuing  ",rc);MFPrintNVector(stdout,x,e);printf(" index %d, bif? %d\n",MFNVGetIndex(x,e),MFNVGetIndex2(x,e));fflush(stdout);}
                    if(0){printf("Tan=\n");MFPrintNKMatrix(stdout,Phi,e);printf("line %d routine %s\n",__LINE__,RoutineName);fflush(stdout);}
                    MFFreeNVector(x,e);

/* ---------------------- Half disk forward into the current branch ------------------------------------------- */

                    x=MFCloneNVector(u,e);
                    MFNVSetIndex2(x,-98,e);   /* This says not to continue this until switch */
                    MFNVSetIndex(x,index1,e);

                    rc=MFAtlasAddChart(A,MFCreateBoundaryChart(A->M,x,Phi,1,-dir,R,e),e);
                    {printf("%d) singular half disk continuing  ",rc);MFPrintNVector(stdout,x,e);printf(" index %d, bif? %d\n",MFNVGetIndex(x,e),MFNVGetIndex2(x,e));fflush(stdout);}
                    if(0){printf("Tan=\n");MFPrintNKMatrix(stdout,Phi,e);printf("line %d routine %s\n",__LINE__,RoutineName);fflush(stdout);}
                    MFFreeNVector(x,e);

/* ---------------------- Singular disk on the bifurcating branch ------------------------------------------- */
   
                    x=MFCloneNVector(u,e);
                    MFNVSetIndex2(w,-98,e); 

                    MFIMFSetStability(A->M,w,Phi,e);
                    bchart=MFAtlasAddChart(A,MFCreateChart(A->M,w,Phi,R,e),e);
                    if(bchart>-1){printf("%d)  point tan1 ",A->nCharts);MFPrintNVector(stdout,u,e);printf(" index %d, bif? %d\n",MFNVGetIndex(u,e),MFNVGetIndex2(u,e));fflush(stdout);}
                    if(bchart>-1){printf("Tan=\n");MFPrintNKMatrix(stdout,Phi,e);printf("line %d routine %s\n",__LINE__,RoutineName);fflush(stdout);}
                    MFFreeNVector(x,e);
   
/* ---------------------- Half disk up into the bifurcating branch ------------------------------------------- */

    
                    s=MFCreateKVector(MFAtlasK(A,e),e);
                    for(i=1;i<A->k;i++)MFKVSetC(s,i,0.,e);MFKVSetC(s,0,1.8*R,e);

                    x=MFCloneNVector(u,e);
                    MFNVSetIndex2(x,-98,e);

                    rc=MFChartEvaluate(A->chart[bchart],s,x,e);
                    Tan=MFIMFTangentSpaceWithGuess(A->M,x,Phi,e);
                    MFIMFSetStability(A->M,x,Tan,e);
                    rc=MFAtlasAddChart(A,MFCreateChart(A->M,x,Tan,R,e),e);
                    if(rc>-1){printf("%d) A step out the continuing branch ",rc);MFPrintNVector(stdout,x,e);printf(" index %d, bif? %d\n",MFNVGetIndex(x,e),MFNVGetIndex2(x,e));fflush(stdout);}
                    if(rc>-1){printf("Tan=\n");MFPrintNKMatrix(stdout,Tan,e);printf("line %d routine %s\n",__LINE__,RoutineName);fflush(stdout);}
/*
                    rc=MFAtlasAddChart(A,MFCreateBoundaryChart(A->M,x,Tan,1,-1,R,e),e);
                    if(rc>-1){printf("%d) half disk bifurcating branch down ",rc);MFPrintNVector(stdout,x,e);printf(" index %d, bif? %d\n",MFNVGetIndex(x,e),MFNVGetIndex2(x,e));fflush(stdout);}
                    if(rc>-1){printf("Tan=\n");MFPrintNKMatrix(stdout,Tan,e);printf("line %d routine %s\n",__LINE__,RoutineName);fflush(stdout);}
 */
                    MFFreeNVector(x,e);
                    MFFreeKVector(s,e);
                    MFFreeNKMatrix(Tan,e);
    
/* ---------------------- Half disk up into the bifurcating branch ------------------------------------------- */
  
                    s=MFCreateKVector(MFAtlasK(A,e),e);
                    for(i=1;i<A->k;i++)MFKVSetC(s,i,0.,e);MFKVSetC(s,0,-1.8*R,e);

                    x=MFCloneNVector(u,e);
                    MFNVSetIndex2(x,-98,e);

                    rc=MFChartEvaluate(A->chart[bchart],s,x,e);
                    Tan=MFIMFTangentSpaceWithGuess(A->M,x,Phi,e);
                    MFIMFSetStability(A->M,x,Tan,e);
                    rc=MFAtlasAddChart(A,MFCreateBoundaryChart(A->M,x,Tan,1, 1,R,e),e);
                    if(rc>-1){printf("%d) half disk bifurcating branch up ",rc);MFPrintNVector(stdout,x,e);printf(" index %d, bif? %d\n",MFNVGetIndex(x,e),MFNVGetIndex2(x,e));fflush(stdout);}
                    if(rc>-1){printf("Tan=\n");MFPrintNKMatrix(stdout,Tan,e);printf("line %d routine %s\n",__LINE__,RoutineName);fflush(stdout);}

                    MFFreeNVector(x,e);
                    MFFreeKVector(s,e);
                    MFFreeNKMatrix(Tan,e);

/* ----------------------------------------------------------------------------------------------------------- */
   
                    if(1||verbose){printf("      done handling Singular Point\n");fflush(stdout);}
   
                   }

                  MFFreeNVector(V,e);
                  MFFreeNKMatrix(Phi,e);
#endif
#endif
                 }else{
                  MFFreeNVector(w,e);
                 }
                if(*TS!=NULL){MFFreeNKMatrix(*TS,e);*TS=NULL;}

                bisection=0;
                if(verbose){printf("Try again\n");fflush(stdout);}
                goto tryagain;
               }else{

#ifdef MFALLOWVERBOSE
                if(verbose){printf("     Same as center (index %d -- type %d), going on to next vertex.\n",MFNVGetIndex(MFChartCenter(A->chart[chart],e),e),MFNVGetIndex2(MFChartCenter(A->chart[chart],e),e));}
#endif

                bisection=0;
                goto nextvertex;
               }
             }
            goto bsect;
           }else{
            int in0,in1;
            int out0,out1;

#ifdef MFALLOWVERBOSE
            if(verbose){printf("     Check for stability change.\n");fflush(stdout);}
#endif

            in0=MFNVGetIndex(MFChartCenter(A->chart[chart],e),e);
            in1=MFNVGetIndex(u,e);
            if(*TS!=NULL){MFFreeNKMatrix(*TS,e);*TS=NULL;}
            *TS=MFIMFTangentSpaceWithGuess(A->M,u,MFChartTangentSpace(A->chart[chart],e),e);
            if(*TS==NULL)return -1;
            MFIMFSetStability(A->M,MFChartCenter(A->chart[chart],e),MFChartTangentSpace(A->chart[chart],e),e);
            MFIMFSetStability(A->M,u,*TS,e);
            out0=MFNVGetIndex(MFChartCenter(A->chart[chart],e),e);
            out1=MFNVGetIndex(u,e);
            if(in0!=out0){printf("   Index of center of chart %d changed from %d to %d\n",chart,in0,out0);fflush(stdout);}
            if(in1!=out1){printf("   Index of new point           changed from %d to %d\n",in1,out1);fflush(stdout);}
            if(MFIMFStop(A->M,MFChartCenter(A->chart[chart],e),MFChartTangentSpace(A->chart[chart],e),u,*TS,e))
             {

#ifdef MFALLOWVERBOSE
              if(verbose){printf("     Stepped over stability change, invoke bisection\n");fflush(stdout);}
#endif
              t=.5;
              tl=0.;
              tr=1.;

/* Invoke bisection, state=looking for singular boundary crossing */

#ifdef MFALLOWVERBOSE
              if(verbose){printf("            interval [%lf, %lf, %lf ] width %le \n",tl,t,tr,tr-tl);fflush(stdout);}
#endif
              index0=MFNVGetIndex(MFChartCenter(A->chart[chart],e),e);
              index1=MFNVGetIndex(u,e);
              bisection=1;
              goto bsect;
             }else{
#ifdef MFALLOWVERBOSE
              if(verbose){printf("     No stability change\n");fflush(stdout);}
#endif

/*          dot=MFNKMatrixDot(A->NSpace,MFChartTangentSpace(A->chart[chart],e),*TS,e);*/
            dot=0.;

              if(0&&dot>.2)
               {
                printf("%s - TangentSpaces too different (dot=%lf), at t=%lf, trying again at %lf.\n",RoutineName,dot,t,t/2);fflush(stdout);

#ifdef MFALLOWVERBOSE
                if(verbose)
#endif
                 {
                  printf("%s - TangentSpaces too different (dot=%lf), at R=%lf, trying again at %lf. Radius of vertex is %lf\n",RoutineName,dot,R,.8*R,rv);fflush(stdout);
                  printf("     Old is \n");MFPrintNKMatrix(stdout,MFChartTangentSpace(A->chart[chart],e),e);printf("\n");fflush(stdout);
                  printf("     New is \n");MFPrintNKMatrix(stdout,Phi,e);printf("\n");fflush(stdout);
                 }
                t=t/2;
                goto bsect;
               }else{

#ifdef MFALLOWVERBOSE
                if(verbose){printf("%s - TangentSpaces OK (dot=%lf), at t=%lf\n",RoutineName,dot,t);fflush(stdout);}
#endif
               }
             }
           }

#ifdef MFALLOWVERBOSE
          if(verbose){printf("      %s - Succeeded at R=%lf, delta=%lf, chart radius=%lf. Radius of vertex is %lf and is on chart %d\n",RoutineName,R,delta,MFChartRadius(A->chart[chart],e),rv,chart);fflush(stdout);}
#endif
          if(dist!=NULL)
            *dist=MFNSpaceDistance(A->NSpace,u,up,e);
#ifdef MFALLOWVERBOSE
          if(verbose){printf("Distance from point to tangent space is %le\n",*dist);fflush(stdout);}
#endif

#ifdef MFTEST
          if(!MFTestBoundaryPoint(A,u,chart,e))
           {
            sprintf(MFAtlasErrorMsg,"Test of Manifold BPt failed.");
            MFSetError(e,12,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
            return -1;
           }
#endif

#ifdef MFALLOWVERBOSE
          if(verbose)
           {
            printf("    found after examining %d charts\n",nChartsTried);fflush(stdout);
            printf("Chart %d is now \n",chart);
            MFPrintChart(stdout,A->chart[chart],e);
            fflush(stdout);
           }
#endif

          if(*TS==NULL)
           {
            t=.8*t;
            goto bsect;
           }

          MFPolytopeSetVertexMark(MFChartPolytope(A->chart[chart],e),v,1,e);
          MFFreeNVector(up,e);
          MFFreeKVector(s,e);
          return chart;
         }else{
#ifdef MFALLOWVERBOSE
          if(verbose){printf("  Interior vertex. ");fflush(stdout);}
#endif
         }
       }
#ifdef MFALLOWVERBOSE
      if(verbose)printf("   Fall through to nextvertex (n=%d)\n",n);
#endif
nextvertex:
      t=1.;
      v++;
#ifdef MFALLOWVERBOSE
      if(verbose){printf("   nextvertex (v=%d, n=%d)\n",v,n);fflush(stdout);}
#endif

#ifdef MFALLOWVERBOSE
      if(verbose&&v<n)printf("     Vertex %d has mark %d, index %d, radius %le/%le\n",v,
             MFPolytopeGetVertexMark(MFChartPolytope(A->chart[chart],e),v,e),
             MFNVGetIndex2(MFChartCenter(A->chart[chart],e),e),
             MFChartRadius(A->chart[chart],e),MFIMFGetRMin(A->M,e),e);
#endif

     }

#ifdef MFALLOWVERBOSE
    if(verbose)
     {
      printf("  Done testing all vertices.\n");
      fflush(stdout);
     }
#endif

nextchart:

#ifdef MFALLOWVERBOSE
    if(verbose)
     {
      printf("  Remove chart (%d) from the boundary list and go on the the next chart.\n",bchart);
      fflush(stdout);
     }
#endif

/* If not "boundary" and any neighbor is "boundary", check all vertices */

/* O.W. */

    if(A->nBnd>0)MFAtlasRemoveChartFromBoundaryList(A,bchart,e);

    if(bchart>=A->nBnd && A->nBnd>0)
     {
      bchart=0;
      chart=-1;
      pass++;
     }
   }
  MFFreeNVector(up,e);
  MFFreeKVector(s,e);

#ifdef MFALLOWVERBOSE
    if(verbose){printf(" done %s\n",RoutineName);fflush(stdout);}
#endif

  return -1;
 }

int MFAddHalfSpaceToAtlas(MFAtlas A,int gI,int left,int right,MFKVector nrm, double onrm, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAddHalfSpaceToAtlas"};
  int result;
  int i,n;
  int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s\n",RoutineName);fflush(stdout);}
#endif

  if(A->nHS>=A->mHS)
   {
    A->mHS+=MFAtlasIncrementToAllocate;

    A->globalIndex=(int*)realloc((void*)A->globalIndex,A->mHS*sizeof(int));
#ifndef MFNOSAFETYNET
    if(A->globalIndex==NULL)
     {
      sprintf(MFAtlasErrorMsg,"Out of memory (1), trying to allocate %d bytes\n",A->mHS*sizeof(int));
      MFSetError(e,12,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return -1;
     }
#endif
    for(i=A->nHS;i<A->mHS;i++)A->globalIndex[i]=-1;

    A->onrm=(double*)realloc((void*)A->onrm,A->mHS*sizeof(double));
#ifndef MFNOSAFETYNET
    if(A->onrm==NULL)
     {
      sprintf(MFAtlasErrorMsg,"Out of memory (2), trying to allocate %d bytes\n",A->mHS*sizeof(double));
      MFSetError(e,12,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return -1;
     }
#endif
    for(i=A->nHS;i<A->mHS;i++)A->onrm[i]=DBL_INFINITY;

    A->nrm=(MFKVector*)realloc((void*)A->nrm,A->mHS*sizeof(MFKVector));
#ifndef MFNOSAFETYNET
    if(A->nrm==NULL)
     {
      sprintf(MFAtlasErrorMsg,"Out of memory (3), trying to allocate %d bytes",A->mHS*sizeof(MFKVector));
      MFSetError(e,12,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return -1;
     }
#endif
    for(i=A->nHS;i<A->mHS;i++)A->nrm[i]=NULL;

    A->leftChart=(int*)realloc((void*)A->leftChart,A->mHS*sizeof(int));
#ifndef MFNOSAFETYNET
    if(A->leftChart==NULL)
     {
      sprintf(MFAtlasErrorMsg,"Out of memory (4), trying to allocate %d bytes",A->mHS*sizeof(int));
      MFSetError(e,12,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return -1;
     }
#endif
    for(i=A->nHS;i<A->mHS;i++)A->leftChart[i]=-1;

    A->rightChart=(int*)realloc((void*)A->rightChart,A->mHS*sizeof(int));
#ifndef MFNOSAFETYNET
    if(A->rightChart==NULL)
     {
      sprintf(MFAtlasErrorMsg,"Out of memory (5), trying to allocate %d bytes",A->mHS*sizeof(int));
      MFSetError(e,12,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return -1;
     }
#endif
    for(i=A->nHS;i<A->mHS;i++)A->rightChart[i]=-1;
   }

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s -- done allocs\n",RoutineName);fflush(stdout);}
#endif

  result=A->nHS;
  A->nHS++;
  A->globalIndex[result]=gI;
  A->onrm[result]=onrm;
  A->nrm[result]=nrm;
  MFRefKVector(nrm,e);
  A->leftChart[result]=left;
  A->rightChart[result]=right;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("done %s\n",RoutineName);fflush(stdout);}
#endif

  return result;
 }

MFNVector MFAtlasChartCenter(MFAtlas A,int chart, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAtlasChartCenter"};
  MFNVector result;

  result=MFChartCenter(A->chart[chart],e);

  return result;
 }

int MFAtlasOffset(MFAtlas A, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAtlasOffset"};
  int result;

  result=A->offset;

  return result;
 }

int MFAtlasGlobalIndex(MFAtlas A,int i, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAtlasGlobalIndex"};
  int result;

  result=A->globalIndex[i];

  return result;
 }

int MFAtlasIsHalfSpaceHyperCube(MFAtlas A,int i, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAtlasIsHalfSpaceHyperCube"};
  int result;

  result=i<A->offset;

  return result;
 }

int MFAtlasLeftPolytope(MFAtlas A,int i, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAtlasLeftPolytope"};
  int result;

  result=A->leftChart[i];

  return result;
 }

int MFAtlasRightPolytope(MFAtlas A,int i, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAtlasRightPolytope"};
  int result;

  result=A->rightChart[i];

  return result;
 }

int MFAtlasIsHalfSpaceInBothPolytopes(MFAtlas A,int hs, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAtlasIsHalfSpaceInBothPolytopes"};
  int inLeft,inRight;
  int i,n;
  int result;
  int ghs;
  int g;
  MFPolytope P;

  ghs=MFAtlasGlobalIndex(A,hs,e);

  inLeft=0;
  P=MFChartPolytope(A->chart[A->leftChart[hs]],e);
  n=MFPolytopeNumberOfFaces(P,e);
  for(i=0;i<n;i++)
   {
    g=MFAtlasGlobalIndex(A,MFPolytopeFaceIndex(P,i,e),e);
    if(g==ghs)inLeft=1;
   }

  inRight=0;
  P=MFChartPolytope(A->chart[A->leftChart[hs]],e);
  n=MFPolytopeNumberOfFaces(P,e);
  for(i=0;i<n;i++)
   {
    g=MFAtlasGlobalIndex(A,MFPolytopeFaceIndex(P,i,e),e);
    if(g==ghs)inRight=1;
   }

  result=inLeft&&inRight;

  return result;
 }

void MFAtlasAddChartToBoundaryList(MFAtlas A,int ball, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAtlasAddChartToBoundaryList"};

  if(A->nBnd>=A->mBnd)
   {
    A->mBnd+=MFAtlasIncrementToAllocateBoundaryList;
    A->bnd=(int*)realloc((void*)(A->bnd),A->mBnd*sizeof(int));
#ifndef MFNOSAFETYNET
    if(A->bnd==NULL)
     {
      sprintf(MFAtlasErrorMsg,"Out of memory, trying to allocate %d bytes",A->mBnd*sizeof(int));
      MFSetError(e,12,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif
   }

  A->bnd[A->nBnd]=ball;
  MFChartSetPositionInBoundaryList(A->chart[ball],A->nBnd,e);
  (A->nBnd)++;

  return;
 }

void MFAtlasRemoveChartFromBoundaryList(MFAtlas A,int ball, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAtlasRemoveChartFromBoundaryList"};
  int i;
  int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s, boundary chart # %d\n",RoutineName,ball);fflush(stdout);}
#endif

  MFChartSetPositionInBoundaryList(A->chart[A->bnd[ball]],-1,e);
  if(ball<A->nBnd-1)
   {
    for(i=ball;i<A->nBnd-1;i++)
     {
      A->bnd[i]=A->bnd[i+1];
      MFChartSetPositionInBoundaryList(A->chart[A->bnd[i]],i,e);
     }
    A->bnd[A->nBnd-1]=-1;
    (A->nBnd)--;
   }else{
    A->bnd[ball]=-1;
    A->nBnd--;
   }

  return;
 }

void MFAddInitialHyperplanes(MFAtlas A, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAddInitialHyperplanes"};
  int i,j,ij;
  MFKVector nrm;
  double onrm;
#ifndef MFSIMPLEX
  double *c;
  int carry;
#else
  MFPolytope simplex;
#endif

#ifndef MFSIMPLEX
  A->offset=1;for(i=0;i<A->k;i++)A->offset*=2;

  c=(double*)malloc(A->k*sizeof(double));
#ifndef MFNOSAFETYNET
  if(c==NULL)
   {
    sprintf(MFAtlasErrorMsg,"Out of memory, trying to allocate %d bytes",A->k*sizeof(double));
    MFSetError(e,12,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif
  for(j=0;j<A->k;j++)c[j]=-1.;

  for(i=0;i<A->offset;i++)
   {
    nrm=MFCreateKVector(A->k,e);
    for(j=0;j<A->k;j++)MFKVSetC(nrm,j,c[j],e);
    onrm=1;
    if(i%2==0)onrm=-1;

    carry=1;
    j=0;
    while(carry&&j<A->k)
     {
      if(c[j]<0)
       {
        c[j]=1;
        carry=0;
       }else{
        c[j]=-1;
        carry=1;
        j++;
       }
     }
    ij=MFAddHalfSpaceToAtlas(A,-1,0,0,nrm,onrm,e);
    MFFreeKVector(nrm,e);
   }
  free(c);
#else
  A->offset=A->k+1;
  simplex=MFCreateSimplexAtOrigin(A->k,1.,e);
  for(i=0;i<A->offset;i++)
   {
    nrm=MFPolytopeFaceNormal(simplex,i,e);
    MFRefKVector(nrm,e);
    onrm=MFPolytopeFaceOrigin(simplex,i,e);
    ij=MFAddHalfSpaceToAtlas(A,-1,0,0,nrm,onrm,e);
    MFFreeKVector(nrm,e);
   }
  MFFreePolytope(simplex,e);
#endif

  return;
 }

int MFAtlasAddChartToList(MFAtlas A,MFChart chart, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAtlasAddChartToList"};
  int i;

  if(!(A->nCharts<A->mCharts))
   {
    A->mCharts+=MFAtlasIncrementToAllocate;

    A->chart=(MFChart*)realloc((void*)A->chart,A->mCharts*sizeof(MFChart));
#ifndef MFNOSAFETYNET
    if(A->chart==NULL)
     {
      sprintf(MFAtlasErrorMsg,"Out of memory, trying to allocate %d bytes",A->mCharts*sizeof(MFChart));
      MFSetError(e,12,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return -1;
     }
#endif
    for(i=A->nCharts;i<A->mCharts;i++)A->chart[i]=NULL;
   }

  A->chart[A->nCharts]=chart;
  MFChartSetPositionInAtlas(chart,A->nCharts,e);
  A->nCharts++;

  return A->nCharts-1;
 }

int MFTestBoundaryPoint(MFAtlas A,MFNVector u, int chart, MFErrorHandler e)
 {
  static char RoutineName[]={"MFTestBoundaryPoint"};
  int i,n;
  double ds,du,dn;
  MFKVector s;
  MFNVector up;
  MFNVector diff;
  int good;
  int verbose=0;

  n=MFAtlasNumberOfCharts(A,e);
  s=MFCreateKVector(A->k,e);
  up=MFCloneNVector(u,e);
  diff=MFCloneNVector(u,e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("Testing point ");MFPrintNVector(stdout,u,e);printf(" on chart %d\n",chart);fflush(stdout);}
#endif

  good=1;
  for(i=0;i<n;i++)
   {
    du=MFNSpaceDistance(A->NSpace,MFChartCenter(A->chart[i],e),u,e);
/*  if(du<.98*MFChartRadius(A->chart[i]),e)*/
     {
      MFChartProjectIntoTangentSpace(A->chart[i],u,s,e);
      ds=MFKVNorm(s,e);

#ifdef MFALLOWVERBOSE
      if(verbose){printf("   %d, center=",i);MFPrintNVector(stdout,MFChartCenter(A->chart[i],e),e);printf(" R=%lf, distance=%lf",MFChartRadius(A->chart[i],e),ds);}
#endif

/*    if(ds<.98*MFChartRadius(A->chart[i]),e)*/
       {
        MFChartPointInTangentSpace(A->chart[i],s,up,e);
        dn=MFNSpaceDistance(A->NSpace,up,u,e);
        if(i!=chart &&
           ds<.98*MFChartRadius(A->chart[i],e) && 
           du<.98*MFChartRadius(A->chart[i],e) &&
           dn<MFEpsilon
          )
         {
          printf("   Test point is on chart %d, ",chart);
          MFPrintNVector(stdout,u,e);
          printf("\n");
          printf("   Failed against chart %d, center=",i,chart);
          MFPrintNVector(stdout,MFChartCenter(A->chart[i],e),e);
          printf(" R=%lf, distance=%lf",MFChartRadius(A->chart[i],e),du,e);
          printf(" tang=%lf",ds);fflush(stdout);
          printf(" perp=%lf\n",dn);fflush(stdout);
          printf("NO GOOD!\n");fflush(stdout);
          if(ds<.5*MFChartRadius(A->chart[i],e))good=0;
          if(du<MFChartRadius(A->chart[i],e))
             printf(" distance in N-space %lf is < radius %lf!\n",du,MFChartRadius(A->chart[i],e));fflush(stdout);
          if(ds<MFChartRadius(A->chart[i],e))
             printf(" distance in tangent space %lf is < radius %lf!\n",ds,MFChartRadius(A->chart[i],e));fflush(stdout);
          if(dn<MFEpsilon)
             printf(" distance in normal space %lf is < epsilon %lf!\n",dn,MFEpsilon);fflush(stdout);
          MFChartProjectIntoTangentSpace(A->chart[chart],u,s,e);
          printf(" In chart %d point is ",chart);MFPrintKVector(stdout,s,e);printf("\n");fflush(stdout);
          MFChartProjectIntoTangentSpace(A->chart[i],u,s,e);
          printf(" In chart %d point is ",i);MFPrintKVector(stdout,s,e);printf("\n");fflush(stdout);
         }
       }
     }
   }
  MFFreeKVector(s,e);
  MFFreeNVector(up,e);
  MFFreeNVector(diff,e);
  return good;
 }

MFChart MFAtlasChart(MFAtlas A, int chart, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAtlasChart"};

#ifdef MFNOCONFIDENCE
  if(A==NULL)
   {
    sprintf(MFAtlasErrorMsg,"Pointer to Atlas (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    return NULL;
   }

  if(chart<0||chart>=A->nCharts)
   {
    sprintf(MFAtlasErrorMsg,"Chart %d (argument 2) is invalid. Must be positive and less than %d (total number of charts).",chart,A->nCharts);
    MFSetError(e,12,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    return NULL;
   }
#endif

  return A->chart[chart];
 }

int MFAtlasNumberOfHalfSpaces(MFAtlas A, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAtlasNumberOfHalfSpaces"};

#ifdef MFNOCONFIDENCE
  if(A==NULL)
   {
    sprintf(MFAtlasErrorMsg,"Pointer to Atlas (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    return -1;
   }
#endif

  return A->nHS;
 }

int MFAtlasHalfSpaceIndex(MFAtlas A,int i, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAtlasHalfSpaceIndex"};

#ifdef MFNOCONFIDENCE
  if(A==NULL)
   {
    sprintf(MFAtlasErrorMsg,"Pointer to Atlas (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    return -1;
   }

  if(i<0||i>=A->nHS)
   {
    sprintf(MFAtlasErrorMsg,"Halfspace %d (argument 2) is invalid. Must be positive and less than %d (total number of Halfspaces).",i,A->nHS);
    MFSetError(e,4,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    return -1;
   }
#endif

  return A->globalIndex[i];
 }

int MFAtlasHalfSpaceLeftChart(MFAtlas A,int i, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAtlasHalfSpaceLeftChart"};

#ifdef MFNOCONFIDENCE
  if(A==NULL)
   {
    sprintf(MFAtlasErrorMsg,"Pointer to Atlas (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    return -1;
   }

  if(i<0||i>=A->nHS)
   {
    sprintf(MFAtlasErrorMsg,"Halfspace %d (argument 2) is invalid. Must be positive and less than %d (total number of Halfspaces).",i,A->nHS);
    MFSetError(e,4,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    return -1;
   }

  if(A->leftChart[i]<0||A->leftChart[i]>=A->nCharts)
   {
    sprintf(MFAtlasErrorMsg,"Warning, Halfspace %d has an invalid left chart (%d). Must be positive and less than %d (total number of charts).",i,A->leftChart[i],A->nCharts);
    MFSetError(e,4,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    return -1;
   }
#endif

  return A->leftChart[i];
 }

int MFAtlasHalfSpaceRightChart(MFAtlas A,int i, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAtlasHalfSpaceRightChart"};

#ifdef MFNOCONFIDENCE
  if(A==NULL)
   {
    sprintf(MFAtlasErrorMsg,"Pointer to Atlas (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    return -1;
   }

  if(i<0||i>=A->nHS)
   {
    sprintf(MFAtlasErrorMsg,"Halfspace %d (argument 2) is invalid. Must be positive and less than %d (total number of Halfspaces).",i,A->nHS);
    MFSetError(e,4,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    return -1;
   }

  if(A->rightChart[i]<0||A->rightChart[i]>=A->nCharts)
   {
    sprintf(MFAtlasErrorMsg,"Warning, Halfspace %d has an invalid right chart (%d). Must be positive and less than %d (total number of charts).",i,A->rightChart[i],A->nCharts);
    MFSetError(e,4,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    return -1;
   }
#endif

  return A->rightChart[i];
 }

MFKVector MFAtlasHalfSpaceNormal(MFAtlas A,int i, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAtlasHalfSpaceNormal"};

#ifdef MFNOCONFIDENCE
  if(A==NULL)
   {
    sprintf(MFAtlasErrorMsg,"Pointer to Atlas (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    return NULL;
   }

  if(i<0||i>=A->nHS)
   {
    sprintf(MFAtlasErrorMsg,"Halfspace %d (argument 2) is invalid. Must be positive and less than %d (total number of Halfspaces).",i,A->nHS);
    MFSetError(e,4,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    return NULL;
   }
#endif

  return A->nrm[i];
 }

double MFAtlasHalfSpaceOrigin(MFAtlas A,int i, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAtlasHalfSpaceOrigin"};

#ifdef MFNOCONFIDENCE
  if(A==NULL)
   {
    sprintf(MFAtlasErrorMsg,"Pointer to Atlas (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    return DBL_QNAN;
   }

  if(i<0||i>=A->nHS)
   {
    sprintf(MFAtlasErrorMsg,"Halfspace %d (argument 2) is invalid. Must be positive and less than %d (total number of Halfspaces).",i,A->nHS);
    MFSetError(e,4,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    return DBL_QNAN;
   }
#endif

  return A->onrm[i];
 }

MFNKMatrix MFAtlasChartTangentSpace(MFAtlas A,int chart, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAtlasChartTangentSpace"};

#ifdef MFNOCONFIDENCE
  if(A==NULL)
   {
    sprintf(MFAtlasErrorMsg,"Pointer to Atlas (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    return NULL;
   }

  if(chart<0||chart>=A->nCharts)
   {
    sprintf(MFAtlasErrorMsg,"Chart %d (argument 2) is invalid. Must be positive and less than %d (total number of charts).",chart,A->nCharts);
    MFSetError(e,4,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    return NULL;
   }
#endif

  return MFChartTangentSpace(A->chart[chart],e);
 }

void MFAtlasSetVerbose(MFAtlas A, int verbose, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAtlasSetVerbose"};

  MFVerbose=verbose;
  return;
 }

void MFAtlasSetDotMin(MFAtlas A, double d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAtlasSetDotMin"};

  MFDotMin=d;
  return;
 }

void MFAtlasSetEpsilon(MFAtlas A, double e, MFErrorHandler err)
 {
  static char RoutineName[]={"MFAtlasSetEpsilon"};

  MFEpsilon=e;
  return;
 }

void MFAtlasSetRMin(MFAtlas A, double r, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAtlasSetRMin"};

  MFRMin=r;
  return;
 }

MFAtlas MFReadAtlas(FILE *fid, MFErrorHandler e)
 {
  static char RoutineName[]={"MFReadAtlas"};
  MFAtlas A;
  int i;
  double tmp;
  int itmp;
  int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s\n");fflush(stdout);}
#endif

  A=(MFAtlas)malloc(sizeof(struct MFAtlasSt));

#ifndef MFNOSAFETYNET
  if(A==NULL)
   {
    sprintf(MFAtlasErrorMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct MFAtlasSt));
    MFSetError(e,12,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  fscanf(fid,"%d\n",&(A->n));

#ifdef MFALLOWVERBOSE
  if(verbose){printf("   n=%d\n",A->n);fflush(stdout);}
#endif

  fscanf(fid,"%d\n",&(A->k));

#ifdef MFALLOWVERBOSE
  if(verbose){printf("   k=%d\n",A->k);fflush(stdout);}
#endif

  A->M=MFReadImplicitMF(fid,e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("   M=%8.8x\n",A->M);fflush(stdout);}
#endif

  A->NSpace=MFReadNSpace(fid,e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("   NSpace=%8.8x\n",A->NSpace);fflush(stdout);}
#endif

  MFIMFSetSpace(A->M,A->NSpace,e);

/* List of Charts */

  fscanf(fid,"%d\n",&(A->nCharts));

#ifdef MFALLOWVERBOSE
  if(verbose){printf("   nCharts=%d\n",A->nCharts);fflush(stdout);}
#endif

  fscanf(fid,"%d\n",&(A->mCharts));

#ifdef MFALLOWVERBOSE
  if(verbose){printf("   mCharts=%d\n",A->mCharts);fflush(stdout);}
#endif

  A->chart=(MFChart*)malloc(A->mCharts*sizeof(MFChart));

#ifndef MFNOSAFETYNET
  if(A->chart==NULL)
   {
    sprintf(MFAtlasErrorMsg,"Out of memory, trying to allocate %d bytes",A->mCharts*sizeof(MFChart));
    MFSetError(e,12,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  for(i=0;i<A->nCharts;i++)
   {
    A->chart[i]=MFReadChart(fid,A,e);
    if(A->chart[i]==NULL){printf("   chart %d is NULL!\n",i);fflush(stdout);}

#ifdef MFALLOWVERBOSE
    if(verbose){printf("   chart[%d]=0x%8.8x\n",i,A->chart[i]);fflush(stdout);}
#endif

   }

/* List of HalfSpaces */

  fscanf(fid,"%d\n",&(A->offset));

#ifdef MFALLOWVERBOSE
  if(verbose){printf("   offset=%d\n",A->offset);fflush(stdout);}
#endif

  fscanf(fid,"%d\n",&(A->nHS));

#ifdef MFALLOWVERBOSE
  if(verbose){printf("   nHS=%d\n",A->nHS);fflush(stdout);}
#endif

  fscanf(fid,"%d\n",&(A->mHS));

#ifdef MFALLOWVERBOSE
  if(verbose){printf("   mHS=%d\n",A->mHS);fflush(stdout);}
#endif

  A->onrm=(double*)malloc(A->mHS*sizeof(double));
#ifndef MFNOSAFETYNET
  if(A->onrm==NULL)
   {
    sprintf(MFAtlasErrorMsg,"Out of memory, trying to allocate %d bytes",A->k*A->mHS*sizeof(double));
    MFSetError(e,12,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

#ifdef MFALLOWVERBOSE
  if(verbose){printf("   read onrm\n");fflush(stdout);}
#endif

  for(i=0;i<A->nHS;i++)fscanf(fid,"%lf\n",&(A->onrm[i]));
  fscanf(fid,"%d\n",&(A->nGlobalIndex));

  A->globalIndex=(int*)malloc(A->mHS*sizeof(int));

#ifndef MFNOSAFETYNET
  if(A->globalIndex==NULL)
   {
    sprintf(MFAtlasErrorMsg,"Out of memory, trying to allocate %d bytes",A->mHS*sizeof(int));
    MFSetError(e,12,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

#ifdef MFALLOWVERBOSE
  if(verbose){printf("   read globalIndex\n");fflush(stdout);}
#endif

  for(i=0;i<A->nHS;i++)fscanf(fid,"%ld\n",&(A->globalIndex[i]));

  A->nrm=(MFKVector*)malloc(A->mHS*sizeof(MFKVector));
#ifndef MFNOSAFETYNET
  if(A->nrm==NULL)
   {
    sprintf(MFAtlasErrorMsg,"Out of memory, trying to allocate %d bytes",A->mHS*sizeof(MFKVector));
    MFSetError(e,12,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

#ifdef MFALLOWVERBOSE
  if(verbose){printf("   read normals\n");fflush(stdout);}
#endif

  for(i=0;i<A->nHS;i++)A->nrm[i]=MFReadKVector(fid,e);

  A->leftChart=(int*)malloc(A->mHS*sizeof(int));

#ifndef MFNOSAFETYNET
  if(A->leftChart==NULL)
   {
    sprintf(MFAtlasErrorMsg,"Out of memory, trying to allocate %d bytes",A->mHS*sizeof(int));
    MFSetError(e,12,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

#ifdef MFALLOWVERBOSE
  if(verbose){printf("   read chart on left\n");fflush(stdout);}
#endif

  for(i=0;i<A->nHS;i++)fscanf(fid,"%ld\n",&(A->leftChart[i]));


  A->rightChart=(int*)malloc(A->mHS*sizeof(int));
#ifndef MFNOSAFETYNET
  if(A->rightChart==NULL)
   {
    sprintf(MFAtlasErrorMsg,"Out of memory, trying to allocate %d bytes",A->mHS*sizeof(int));
    MFSetError(e,12,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

#ifdef MFALLOWVERBOSE
  if(verbose){printf("   read chart on right\n");fflush(stdout);}
#endif

  for(i=0;i<A->nHS;i++)fscanf(fid,"%ld\n",&(A->rightChart[i]));
#ifdef DONTREAD
  for(i=0;i<A->nHS;i++)fscanf(fid,"%ld\n",&itmp);
#endif

/* Binary tree */

#ifdef MFALLOWVERBOSE
  if(verbose){printf("   read binary tree\n");fflush(stdout);}
#endif

  A->BTree=MFReadBinaryTree(fid,e);

/* List of Boundarys */

  fscanf(fid,"%d\n",&(A->nBnd));
  fscanf(fid,"%d\n",&(A->mBnd));

#ifdef MFALLOWVERBOSE
  if(verbose){printf("   nBnd=%d\n",A->nBnd);fflush(stdout);}
#endif

  A->bnd=(int*)malloc(A->mBnd*sizeof(int));

#ifndef MFNOSAFETYNET
  if(A->bnd==NULL)
   {
    sprintf(MFAtlasErrorMsg,"Out of memory, trying to allocate %d bytes",A->mBnd*sizeof(int));
    MFSetError(e,12,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  for(i=0;i<A->nBnd;i++)
   {
    fscanf(fid,"%d\n",&(A->bnd[i]));
    MFChartSetPositionInBoundaryList(A->chart[A->bnd[i]],i,e);
   }

#ifdef MFALLOWVERBOSE
  if(verbose){printf("   done %s\n",RoutineName);fflush(stdout);}
#endif

  return A;
 }

void MFWriteAtlas(FILE *fid,MFAtlas A, MFErrorHandler e)
 {
  static char RoutineName[]={"MFWriteAtlas"};
  int i;
  int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s\n");fflush(stdout);}
#endif

  fprintf(fid,"%d\n",A->n);
  fprintf(fid,"%d\n",A->k);
  MFWriteImplicitMF(fid,A->M,e);
  MFWriteNSpace(fid,A->NSpace,e);

/* List of Charts */

  fprintf(fid,"%d\n",A->nCharts);
  fprintf(fid,"%d\n",A->mCharts);
  for(i=0;i<A->nCharts;i++)
   {
    MFWriteChart(fid,A->chart[i],e);
    MFChartClean(A->chart[i],e);
   }

/* List of HalfSpaces */

  fprintf(fid,"%d\n",A->offset);
  fprintf(fid,"%d\n",A->nHS);
  fprintf(fid,"%d\n",A->mHS);
  for(i=0;i<A->nHS;i++)
    fprintf(fid,"%lf\n",A->onrm[i]);
  fprintf(fid,"%d\n",A->nGlobalIndex);
  for(i=0;i<A->nHS;i++)fprintf(fid,"%ld\n",A->globalIndex[i]);
  for(i=0;i<A->nHS;i++)MFWriteKVector(fid,A->nrm[i],e);
  for(i=0;i<A->nHS;i++)fprintf(fid,"%ld\n",A->leftChart[i]);
  for(i=0;i<A->nHS;i++)fprintf(fid,"%ld\n",A->rightChart[i]);

/* Binary tree */

#ifndef NOBINARYTREE
  MFWriteBinaryTree(fid,A->BTree,e);
#endif

#ifndef NOBOUNDARYLIST
/* List of Boundarys */

  fprintf(fid,"%d\n",A->nBnd);
  fprintf(fid,"%d\n",A->mBnd);
  for(i=0;i<A->nBnd;i++)fprintf(fid,"%d\n",A->bnd[i]);
#endif

  return;
 }

MFImplicitMF MFAtlasMF(MFAtlas A, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAtlasMF"};

#ifdef MFNOCONFIDENCE
  if(A==NULL)
   {
    sprintf(MFAtlasErrorMsg,"Pointer to Atlas (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    return NULL;
   }
#endif

  return A->M;
 }

int MFTestAtlas(MFAtlas A, MFErrorHandler e)
 {
  static char RoutineName[]={"MFTestAtlas"};
  int i,n,j,m,l;
  MFPolytope P;
  MFChart chart;
  MFKVector s;
  MFNVector u;
  int good;
  int verbose=0;

  n=MFAtlasNumberOfCharts(A,e);
  s=MFCreateKVector(A->k,e);
  if(n>0)u=MFCloneNVector(MFChartCenter(A->chart[0],e),e);
   else u=NULL;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("Testing Atlas ");fflush(stdout);}
#endif

/* 1. Is any polytope vertex inside another any other polytope? */

  good=1;
  for(i=0;i<n;i++)
   {
    chart=MFAtlasChart(A,i,e);
    P=MFChartPolytope(chart,e);
    good=good&&MFTestPolytope(P,e);
   }

  MFFreeKVector(s,e);
  if(u!=NULL)MFFreeNVector(u,e);

  return good;
 }

void MFAtlasReduceChartRadius(MFAtlas A,int chart,double R, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAtlasReduceChartRadius"};
  int k;
  MFChart c;
  MFPolytope P,Q;
  double R0;
  double Ri,Rj;
  MFNVector ui,uj;
  MFNKMatrix Phij;

  int nn,f;
  int *ncharts;
  double *no;
  MFKVector *nnrm;
  int *nindx;
  int i,j,m;
  int nHS;
  double *o;
  MFKVector *nrm;
  int qindx;
  double qo;
  MFKVector qnrm;
  MFNVector dir;
  MFKVector projdir;
  int *indx;
  double dist;
  int BBDimension=0;
  double *BBCenter=NULL;
  int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s\n",RoutineName);fflush(stdout);}
#endif

  c=MFAtlasChart(A,chart,e);
  k=MFAtlasK(A,e);
  R0=MFChartRadius(c,e);
  P=MFChartPolytope(c,e);
  dir=MFCloneNVector(MFAtlasChartCenter(A,chart,e),e);
  projdir=MFCreateKVector(A->k,e);

  /* Check that R<R_chart */

#ifdef MFNOCONFIDENCE
  if(R>=R0)
   {
    sprintf(MFAtlasErrorMsg,"New radius (%lf) must be less than old (%lf), chart will not be changed.",R,R0);
    MFSetError(e,4,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  /* Make a list of all charts that border chart */

#ifdef MFALLOWVERBOSE
  if(verbose){printf("Reducing chart %d from R=%lf to %lf\n",chart,R0,R);fflush(stdout);}
#endif

  nn=MFPolytopeNumberOfFaces(P,e);
  ncharts=(int*)malloc(nn*sizeof(int));

#ifndef MFNOSAFETYNET
  if(ncharts==NULL)
   {
    sprintf(MFAtlasErrorMsg,"Out of memory, trying to allocate %d bytes",nn*sizeof(int));
    MFSetError(e,12,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  nindx=(int*)malloc(nn*sizeof(int));

#ifndef MFNOSAFETYNET
  if(nindx==NULL)
   {
    sprintf(MFAtlasErrorMsg,"Out of memory, trying to allocate %d bytes",nn*sizeof(int));
    MFSetError(e,12,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  no=(double*)malloc(nn*sizeof(double));

#ifndef MFNOSAFETYNET
  if(no==NULL)
   {
    sprintf(MFAtlasErrorMsg,"Out of memory, trying to allocate %d bytes",nn*sizeof(int));
    MFSetError(e,12,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  nnrm=(MFKVector*)malloc(nn*sizeof(MFKVector));

#ifndef MFNOSAFETYNET
  if(nnrm==NULL)
   {
    sprintf(MFAtlasErrorMsg,"Out of memory, trying to allocate %d bytes",nn*sizeof(MFKVector));
    MFSetError(e,12,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

#ifdef MFALLOWVERBOSE
  if(verbose){printf(" chart %d adjoins %d other charts\n",chart,nn);fflush(stdout);}
#endif

  nn=0;
  for(i=0;i<MFPolytopeNumberOfFaces(P,e);i++)
   {
    f=MFPolytopeFaceIndex(P,i,e);
    if(!MFAtlasIsHalfSpaceHyperCube(A,f,e))
     {
      nindx[nn]=f;
      if((ncharts[nn]=MFAtlasHalfSpaceRightChart(A,f,e))==chart)
          ncharts[nn]=MFAtlasHalfSpaceLeftChart(A,f,e);
      no[nn]=MFAtlasHalfSpaceOrigin(A,f,e);
      nnrm[nn]=MFAtlasHalfSpaceNormal(A,f,e);

#ifdef MFALLOWVERBOSE
      if(verbose){printf("    %d\n",ncharts[nn]);fflush(stdout);}
#endif

      nn++;
     }
   }

  /* Remove all halfspaces from chart */

#ifdef MFALLOWVERBOSE
  if(verbose){printf(" Remove the polytope of thisAtlas chart %d",chart);fflush(stdout);}
#endif

  if(P!=NULL)MFFreePolytope(P,e);

  /* Set P_chart=C_R */
#ifdef MFALLOWVERBOSE
  if(verbose){printf(" Create a new hypercube for the polytope of thisAtlas chart %d\n",chart);fflush(stdout);}
#endif

#ifndef MFSIMPLEX
  P=MFCreateHyperCubeAtOrigin(k,1.05*R,e);
#else
  P=MFCreateSimplexAtOrigin(k,1.05*R,e);
#endif
  MFSetChartPolytope(c,P,e);
  MFSetChartRadius(c,R,e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf(" For each adjoining chart:\n");fflush(stdout);}
#endif

  for(i=0;i<nn;i++)
   {

#ifdef MFALLOWVERBOSE
    if(verbose){printf("chart %d:\n",ncharts[i]);fflush(stdout);}
#endif

    Q=MFChartPolytope(A->chart[ncharts[i]],e);

  /* Make list of half spaces to be removed from Q */

    nHS=MFPolytopeNumberOfFaces(Q,e);
    o=(double*)malloc(nHS*sizeof(double));

#ifndef MFNOSAFETYNET
    if(o==NULL)
     {
      sprintf(MFAtlasErrorMsg,"Out of memory, trying to allocate %d bytes",nHS*sizeof(double));
      MFSetError(e,12,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

    nrm=(MFKVector*)malloc(nHS*sizeof(MFKVector));

#ifndef MFNOSAFETYNET
    if(nrm==NULL)
     {
      sprintf(MFAtlasErrorMsg,"Out of memory, trying to allocate %d bytes",nHS*sizeof(MFKVector));
      MFSetError(e,12,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

    indx=(int*)malloc(nHS*sizeof(int));

#ifndef MFNOSAFETYNET
    if(indx==NULL)
     {
      sprintf(MFAtlasErrorMsg,"Out of memory, trying to allocate %d bytes",nHS*sizeof(int));
      MFSetError(e,12,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

#ifdef MFALLOWVERBOSE
    if(verbose){printf("     there are %d faces:\n",nHS);fflush(stdout);}
#endif

    m=0;
    qindx=-1;
    for(j=0;j<nHS;j++)
     {
      f=MFPolytopeFaceIndex(Q,j,e);
      if(!MFAtlasIsHalfSpaceHyperCube(A,f,e))
       {
        if(MFAtlasHalfSpaceRightChart(A,f,e)!=chart
                        && MFAtlasHalfSpaceLeftChart(A,f,e)!= chart)
         {
          indx[m]=f;
          o[m]=MFAtlasHalfSpaceOrigin(A,f,e);
          nrm[m]=MFAtlasHalfSpaceNormal(A,f,e);
          m++;
         }else{
          qindx=f;
          qo=MFAtlasHalfSpaceOrigin(A,f,e);
          qnrm=MFAtlasHalfSpaceNormal(A,f,e);
         }
       }
     }
    if(0&&qindx<0){printf("Problem, qindx[%d] is <0\n",i);fflush(stdout);abort();}
    if(Q!=NULL)MFFreePolytope(Q,e);

  /* Set Q=C_R */

#ifdef MFALLOWVERBOSE
   if(verbose){printf("  Create a new Hypercube for chart %d\n",ncharts[i]);fflush(stdout);}
#endif

    R0=MFAtlasChartRadius(A,ncharts[i],e);

#ifndef MFSIMPLEX
    Q=MFCreateHyperCubeAtOrigin(k,1.05*R0,e);
#else
    Q=MFCreateSimplexAtOrigin(k,1.05*R0,e);
#endif

    MFSetChartPolytope(A->chart[ncharts[i]],Q,e);

  /* Subtract each HS from Q, omitting HS_P */

#ifdef MFALLOWVERBOSE
    if(verbose){printf("  Subtract %d halfspaces from chart %d\n",m,ncharts[i]);fflush(stdout);}
#endif

    for(j=0;j<m;j++)
     {

#ifdef MFALLOWVERBOSE
      if(verbose){printf("     %d, %d\n",j,indx[j]);fflush(stdout);}
#endif

      MFSubtractHalfSpaceFromPolytope(Q,indx[j],nrm[j],o[j],e);
     }

#ifdef MFALLOWVERBOSE
    if(verbose){MFPrintPolytope(stdout,Q,e);fflush(stdout);}
#endif

    free(o);
    free(nrm);
    free(indx);

#ifdef MFALLOWVERBOSE
    if(verbose){printf("Done, now subtract reduced/old intersections\n");fflush(stdout);}
#endif

/* Now deal with reduced polytope */

    Ri=MFChartRadius(c,e);
    ui=MFChartCenter(c,e);
    Rj=MFChartRadius(A->chart[ncharts[i]],e);
    uj=MFChartCenter(A->chart[ncharts[i]],e);
    Phij=MFChartTangentSpace(A->chart[ncharts[i]],e);

  /* Subtract new HS_P from Q */

    if(qindx>-1)
     {

#ifdef MFALLOWVERBOSE
      if(verbose){printf("  Subtract updated halfspaces from chart %d\n",ncharts[i]);fflush(stdout);}
#endif

      MFNSpaceDirection(A->NSpace,uj,ui,dir,e);
      MFChartProjectVectorIntoTangentSpace(A->chart[ncharts[i]],dir,qnrm,e);
  printf("Check if two charts overlap Line %d, File %s\n",__LINE__,__FILE__);
      dist=MFNSpaceInner(A->NSpace,dir,dir,e);

      qo=.5*(dist+Rj*Rj-Ri*Ri);
      A->onrm[qindx]=qo;

#ifdef MFALLOWVERBOSE
      if(verbose)
       {
        printf("  qindx=%d\n",qindx);fflush(stdout);
        printf("   direction from %d to %d is ",ncharts[i],chart);MFPrintNVector(stdout,dir,e);printf("\n");fflush(stdout);
        printf("   in TS %d is ",ncharts[i]);MFPrintKVector(stdout,qnrm,e);printf("\n");fflush(stdout);
        printf("   distance in TS is %lf\n",sqrt(dist));fflush(stdout);
       }
#endif
  
      if(1)MFSubtractHalfSpaceFromPolytope(Q,qindx,qnrm,qo,e);

  /* Subtract new HS_Q from P */

#ifdef MFALLOWVERBOSE
      if(verbose){printf("  Subtract updated halfspaces from chart %d\n",chart);fflush(stdout);}
#endif

      MFNSpaceDirection(A->NSpace,ui,uj,dir,e);
      MFChartProjectVectorIntoTangentSpace(c,dir,nnrm[i],e);
  printf("Find halfspace to subtract Line %d, File %s\n",__LINE__,__FILE__);
      dist=MFNSpaceInner(A->NSpace,dir,dir,e);
      no[i]=.5*(dist+Ri*Ri-Rj*Rj);
      A->onrm[nindx[i]]=no[i];

#ifdef MFALLOWVERBOSE
      if(verbose)
       {
        printf("  nindx[i]=%d\n",nindx[i]);fflush(stdout);
        printf("   direction from %d to %d is ",chart,ncharts[i]);MFPrintNVector(stdout,dir,e);printf("\n");fflush(stdout);
        printf("   in TS %d is ",chart);MFPrintKVector(stdout,nnrm[i],e);printf("\n");fflush(stdout);
        printf("   distance in TS is %lf\n",sqrt(dist));fflush(stdout);
       }
#endif

      MFSubtractHalfSpaceFromPolytope(P,nindx[i],nnrm[i],no[i],e);
     }
   }

#ifdef MFALLOWVERBOSE
  if(verbose){printf("  free temporary\n");fflush(stdout);}
#endif

  free(ncharts);
  free(nindx);
  free(no);
  free(nnrm);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("done %s\n",RoutineName);fflush(stdout);}
#endif

  BBDimension=MFIMFProjectToBB(A->M,MFChartCenter(A->chart[chart],e),NULL,e);
  BBCenter=(double*)malloc(BBDimension*sizeof(double));

#ifndef MFNOSAFETYNET
  if(BBCenter==NULL)
   {
    sprintf(MFAtlasErrorMsg,"Out of memory, trying to allocate %d bytes",BBDimension*sizeof(double));
    MFSetError(e,12,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

#ifdef MFALLOWVERBOSE
  if(verbose){printf("       %s, BBDimension=%d\n",RoutineName,BBDimension);fflush(stdout);}
#endif

  MFIMFProjectToBB(A->M,MFChartCenter(A->chart[chart],e),BBCenter,e);

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("      BBCenter=(%lf",BBCenter[0]);
    for(i=1;i<BBDimension;i++)printf(",%lf",BBCenter[i]);
    printf(" )\n");fflush(stdout);
   }
#endif

  MFFreeNVector(dir,e);
  MFFreeKVector(projdir,e);

  MFRecomputeBoundingBoxes(A->BTree,chart,BBCenter,MFChartRadius(A->chart[chart],e),e);

  return;
 }

MFNSpace MFAtlasNSpace(MFAtlas A, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAtlasNSpace"};

  return A->NSpace;
 }

#if 0
double MFNKMatrixDot(MFNSpace space,MFNKMatrix A,MFNKMatrix B, MFErrorHandler e)
 {
  static char RoutineName[]={"MFNKMatrixDot"};
  double *Q;
  double dot;
  MFNVector col;
  MFKVector C;
  int k,n;
  int i,j;
  int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("%s\n",RoutineName);fflush(stdout);
    printf("Phi_i:\n");MFPrintNKMatrix(stdout,A,e);
    printf("Phi_j:\n");MFPrintNKMatrix(stdout,B,e);
   }
#endif

  n=MFNKMatrixN(A,e);
  k=MFNKMatrixK(A,e);
  if(k>3)return 0.;
  Q=(double*)malloc(k*k*sizeof(double));

#ifndef MFNOSAFETYNET
  if(Q==NULL)
   {
    sprintf(MFAtlasErrorMsg,"Out of memory, trying to allocate %d bytes",k*k*sizeof(double));
    MFSetError(e,12,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0.;
   }
#endif

  C=MFCreateKVector(k,e);
  for(i=0;i<k;i++)
   {
    col=MFMColumn(B,i,e);
    MFMVMulT(space,A,col,C,e);
    for(j=0;j<k;j++)Q[i+j*k]=MFKV_C(C,j,e);
    MFFreeNVector(col,e);
   }

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("Phi_i^TPhi_j=\n");
    for(i=0;i<k;i++)
     {
      printf("[");
      for(j=0;j<k;j++)
       {
        if(j>0)printf(",");
        printf("%lf",Q[i+k*j]);
       }
      printf("]\n");
     }
    fflush(stdout);
   }
#endif

  dot=pow(fabs(MFAtlasDet(k,Q,e)),1./k);
  free(Q);
  MFFreeKVector(C,e);
  return dot;
 }
#endif

double MFNKMatrixDot(MFNSpace space,MFNKMatrix A,MFNKMatrix B, MFErrorHandler e)
 {
  double dot;
  double tmpdot;
  MFNVector col;
  MFNVector proj;
  MFKVector tmp;
  int i,n,k;

  n=MFNKMatrixN(B,e);
  k=MFNKMatrixK(B,e);

  col=MFMColumn(B,0,e);
  proj=MFCloneNVector(col,e);
  MFFreeNVector(col,e);

  tmp=MFCreateKVector(k,e);

  dot=1;
  for(i=0;i<k;i++)
   {
    col=MFMColumn(B,i,e);

    MFMVMulT(space,A,col,tmp,e);        /* v= Phi Phi^T Psi Psi^T phi_k */
    MFMVMul (space,A,tmp,proj,e);

    tmpdot=MFNSpaceInner(space,col,proj,e);
    if(tmpdot<dot)dot=tmpdot;
    MFFreeNVector(col,e);
   }
  MFFreeKVector(tmp,e);
  MFFreeNVector(proj,e);

  return dot;
 }

double MFAtlasDet(int n, double *A, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAtlasDet"};

  if(n==1)return A[0];
   else if(n==2)return A[0]*A[3]-A[1]*A[2];
   else if(n==3)return A[0]*(A[4]*A[8]-A[5]*A[7])
                      -A[1]*(A[3]*A[8]-A[5]*A[6])
                      +A[2]*(A[3]*A[7]-A[4]*A[6]);
  
  return MFFullProdEV(n,0,A,e);
 }

int MFAtlasAddChartWithAll(MFAtlas A,MFNVector u,MFNKMatrix Phi,double R, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAtlasAddChartWithAll"};
  return MFAtlasAddChart(A,MFCreateChart(A->M,u,Phi,R,e),e);
 }

int MFAtlasAddChart(MFAtlas A,MFChart C, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAtlasAddChart"};
  int j,jj;
  int ii,n;
  MFKVector nrm;
  double onrm;
  double dist,sdist;
  double ddist;
  double tdist;
  double dot;
  MFNVector dir;
  MFPolytope P;
  int chart;
  int bchart;
  double Rstar;
  MFNVector uj;
  MFNKMatrix Phij;
  double Rj;
  int nInt,ibnd;
  int ni;
  double RMin;
  double R0;
  double RMax;
  double sqrtTwo;
  int BBDimension=0;
  double *BBCenter=NULL;
  int mark;
  int type0,type1;
  int boundary;

  int ij,ji;
  MFListOfCharts L;

  int verbose=0;

  MFNVector u;
  MFNKMatrix Phi;
  double R;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s\n",RoutineName);fflush(stdout);}
#endif

  sqrtTwo=sqrt(2.);

#ifdef MFNOCONFIDENCE
  if(A==NULL)
   {
    sprintf(MFAtlasErrorMsg,"Pointer to Atlas (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    return -1;
   }

  if(C==NULL)
   {
    sprintf(MFAtlasErrorMsg,"Pointer to Chart (argument 2) is NULL");
    MFSetError(e,4,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    return -1;
   }
#endif

  RMin=-1.;
  RMax=-1.;
  Rstar=-1.;

  u=MFChartCenter(C,e);
  Phi=MFChartTangentSpace(C,e);
  R=MFChartRadius(C,e);

#ifdef MFALLOWVERBOSE
  if(verbose)printf("     R=%lf\n",R);
#endif

/* Add it to the list */

/* Defer to end? */

#ifdef MFALLOWVERBOSE
  if(verbose){printf("       add chart to list\n");fflush(stdout);}
#endif

  chart=MFAtlasAddChartToList(A,C,e);

  MFChartSetNearBoundary(A->chart[chart],1,e); /* Only if not paging! @@@ */
  if(MFNVGetIndex2(u,e)==-99||MFNVGetIndex2(u,e)==-98)
   {
    MFChartSetSingular(A->chart[chart],e);
    MFChartSetNearBoundary(A->chart[chart],1,e);
   }

  dir=MFCloneNVector(u,e);

/* Intersections with other charts */

  nInt=0;
  mark=-1;

/* Add new chart to binary tree */

  BBDimension=MFIMFProjectToBB(A->M,u,NULL,e);
  BBCenter=(double*)malloc(BBDimension*sizeof(double));

#ifndef MFNOSAFETYNET
  if(BBCenter==NULL)
   {
    sprintf(MFAtlasErrorMsg,"Out of memory trying to allocate %d bytes",BBDimension*sizeof(double));
    MFSetError(e,12,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return -1;
   }
#endif

  if(A->BTree==NULL)A->BTree=MFCreateBinaryTree(BBDimension,e);

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("       Project chart center onto space of binary tree\n");fflush(stdout);
    printf("       %s, BBDimension=%d\n",RoutineName,BBDimension);fflush(stdout);
   }
#endif

  MFIMFProjectToBB(A->M,u,BBCenter,e);

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("      BBCenter=(%lf",BBCenter[0]);
    for(jj=1;jj<BBDimension;jj++)printf(",%lf",BBCenter[jj]);
    printf(" )\n");fflush(stdout);
   }
#endif

/* Find minimum radius of intersecting charts */
/* Also find minimum index of intersecting charts (add one to get distance from last labeled pt)*/
/*                                                (reset index to zero after saving)*/

#define CHECKNEIGHBORS_NO
#ifdef CHECKNEIGHBORS

  R0=R;
  if(A->nCharts>0)
   {

#ifdef MFALLOWVERBOSE
    if(verbose){printf("checking neighboring radii\n");fflush(stdout);}
#endif

    L=MFCreateListOfNearbyCharts(A->BTree,BBCenter,R,e);
    ni=MFNumberOfIntersectingCharts(L,e);

#ifdef MFALLOWVERBOSE
    if(verbose){printf("  There are %d neighbors\n",ni);fflush(stdout);}
#endif

    RMin=0.;
    RMax=0.;
    for(jj=0;jj<ni;jj++)
     {
      j=MFIntersectingChart(L,jj,e);
      if(jj==0||MFChartSuggestedRadius(A->chart[j],e)<RMin)
        RMin=MFChartSuggestedRadius(A->chart[j],e);
      if(jj==0||MFChartSuggestedRadius(A->chart[j],e)>RMax)
        RMax=MFChartSuggestedRadius(A->chart[j],e);
  
#ifdef MFALLOWVERBOSE
      if(verbose){printf("    neighbor %d, radius=%lf\n",j,MFChartSuggestedRadius(A->chart[j]));fflush(stdout);}
#endif

     }
    MFFreeListOfIntersectingCharts(L,e);
/*  if(ni>0&&R>sqrt(2.)*RMin)R=.8*sqrt(2.)*RMin;*/
   }

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("Done checking neighboring radii, minimum nearby radius=%lf, max. %lf. reduced from %lf to %lf\n",RMin,RMax,R0,R);
    printf("Done checking neighboring radii\n");fflush(stdout);
    printf("       do intersections\n");fflush(stdout);
    printf("       Add chart to binary tree\n");fflush(stdout);
   }
#endif
#endif  /* CHECKNEIGHBORS */

  if(A->nCharts>1)
   {

#ifdef MFALLOWVERBOSE
    if(verbose){printf("       Create list of intersecting charts\n");fflush(stdout);}
#endif

    L=MFCreateListOfIntersectingCharts(A->BTree,chart,BBCenter,MFChartRadius(A->chart[chart],e),e);

#ifdef MFALLOWVERBOSE
    if(verbose){printf("       Loop through list\n");fflush(stdout);}
#endif

    ni=MFNumberOfIntersectingCharts(L,e);

#ifdef MFALLOWVERBOSE
    if(verbose){printf("  There are %d nearby charts out of %d\n",ni,A->nCharts);fflush(stdout);}
#endif

    for(jj=0;jj<ni;jj++)
     {
      j=MFIntersectingChart(L,jj,e);
      if(MFChartPaged(A->chart[j],e))continue;

      if(A->isNear!=NULL && !(A->isNear)(A,A->chart[j],A->chart[chart],e))continue;

#ifdef MFNOCONFIDENCE
      if(!MFChartNearBoundary(A->chart[j],e))
       {
        sprintf(MFAtlasErrorMsg,"Intersecting chart is not near the boundary.");
        MFSetError(e,4,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
       }
#endif
      if(!MFChartNearBoundary(A->chart[j],e))continue;

#if 0

      type0=MFNVGetIndex2(MFChartCenter(A->chart[j],e),e);
      type1=MFNVGetIndex2(MFChartCenter(A->chart[chart],e),e);

      boundary=MFPolytopeLargestRadiusOfVertex(MFChartPolytope(A->chart[chart],e),e)/MFChartRadius(A->chart[chart],e)>1.;

/*    if(  type0<0  && !(type1<0) && boundary)continue;
      if(!(type0<0) &&   type1<0  && boundary)continue;
 */

/* -99 is the chart on the old branch */

      if(type0==-99&&type1==  0)continue;
      if(type0==  0&&type1==-99)continue;

/* -98 is the chart on the new branch */

      if(type0==-98&&type1==  0)continue;
      if(type0==  0&&type1==-98)continue;
   
      if(type0==-99&&type1==-98)continue;
      if(type0==-98&&type1==-99)continue;

/*    type0=MFNVGetIndex(MFChartCenter(A->chart[j],e),e);
      type1=MFNVGetIndex(MFChartCenter(A->chart[chart],e),e);

      if(type0!=type1)continue;*/
#endif

/* Is the norm of Phi_j^T Phi_i within epsilon of 1? */

      Rj=MFChartRadius(A->chart[j],e);
      uj=MFChartCenter(A->chart[j],e);
      Phij=MFChartTangentSpace(A->chart[j],e);

#ifdef MFALLOWVERBOSE
      if(verbose){printf("\n   chart %d (R=%lf) and %d(R=%lf)\n",chart,R,j,Rj);fflush(stdout);}
#endif

      MFNSpaceDirection(A->NSpace,u,uj,dir,e);

#ifdef MFALLOWVERBOSE
      if(verbose){printf("   direction from %d to %d is ",chart,j);MFPrintNVector(stdout,dir,e);printf("\n");fflush(stdout);}
#endif

      nrm=MFCreateKVector(A->k,e);
      MFChartProjectVectorIntoTangentSpace(A->chart[chart],dir,nrm,e);

#ifdef MFALLOWVERBOSE
      if(verbose){printf("   in TS %d is ",chart);MFPrintKVector(stdout,nrm,e);printf("\n");fflush(stdout);}
#endif

      dist=MFKVDot(nrm,nrm,e);
      sdist=sqrt(dist);
      ddist=MFNSpaceDistance(A->NSpace,u,uj,e);

#ifdef MFALLOWVERBOSE
      if(verbose){printf("   distance in n-space is %lf, in TS is %lf\n",ddist,sdist);fflush(stdout);}
#endif

      if(ddist<1.e-10 && !(MFNVGetIndex2(u,e)<0||MFNVGetIndex2(uj,e)<0))
       {
#ifdef MFNOCONFIDENCE
        sprintf(MFAtlasErrorMsg,"Two points the same! %d and %d.",chart,j);
        fprintf(stderr,"Index2(u)=%d, Index2(uj)=%d.\n",MFNVGetIndex2(u,e),MFNVGetIndex2(uj,e));fflush(stderr);
        MFSetError(e,12,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
#endif
        MFFreeKVector(nrm,e);
        return -12;
       }
  
#ifdef MFALLOWVERBOSE
      if(verbose){printf("   check dist > (R+Rj)*(R+Rj) || ddist > R+Rj\n",dist,(R+Rj)*(R+Rj),ddist,R+Rj);
                  printf("   check %lf > %lf || %lf > %lf\n",dist,(R+Rj)*(R+Rj),ddist,R+Rj);fflush(stdout);}
#endif

      if(dist>(R+Rj)*(R+Rj) || ddist>R+Rj)
       {

#ifdef MFALLOWVERBOSE
        if(verbose){printf("  These do not intersect - sphere test\n");fflush(stdout);}
#endif

        MFFreeKVector(nrm,e);
        continue;
       }
  
#ifdef MFUSETANGENTTEST
      if(0&&MFIMFStop(A->M,u,Phi,uj,Phij,e))
       {

#ifdef MFALLOWVERBOSE
        if(verbose){printf("  These do not intersect - index changes\n");fflush(stdout);}
#endif

        MFFreeKVector(nrm,e);
        continue;
       }

#define normalTolerance .01
      if(0&&(tdist=MFNSpaceTangentDistance(A->NSpace,Phi, MFChartTangentSpace(A->chart[j],e),e))>normalTolerance*R)
       {

#ifdef MFALLOWVERBOSE
        if(verbose){printf("  These do not intersect - tangent test %le\n",tdist);fflush(stdout);}
#endif

        MFFreeKVector(nrm,e);
        continue;
       }
#endif /* MFUSETANGENTTEST */

#ifdef MFUSECYLINDERTEST
      if(0&&ddist*ddist-dist>16*MFEpsilon*MFEpsilon)
       {

#ifdef MFALLOWVERBOSE
        if(verbose){printf("  These do not intersect - cylinder test\n");fflush(stdout);}
#endif

        MFFreeKVector(nrm,e);
        continue;
       }
#endif  /* MFUSECYLINDERTEST */

      dot=MFNKMatrixDot(A->NSpace,MFChartTangentSpace(A->chart[j],e),Phi,e);

#ifdef MFALLOWVERBOSE
      if(verbose){printf("Matrix Dot with %d is %lf\n",j,dot);fflush(stdout);}
#endif

      if(dot<MFDotMin) /* Was .9 for Rod */
       {
        MFFreeKVector(nrm,e);

#ifdef MFALLOWVERBOSE
        if(verbose){printf("These do not intersect, (%d,%d) - too different a TS %lf\n",chart,j,dot);fflush(stdout);}
#endif

        continue;   /*  Tangent Space  */
       }

#ifdef MFALLOWVERBOSE
      if(verbose){printf("This pair (%d,%d) has an OK difference in TS %lf\n",chart,j,dot);fflush(stdout);}
#endif

      if(!MFChartNearBoundary(A->chart[j],e))
       {
        printf("   This pair intersect %d and %d, and %d is not near the boundary!\n",chart,j,j);fflush(stdout);
       }

#ifdef MFALLOWVERBOSE
      if(verbose){printf("   These Intersect\n");fflush(stdout);}
#endif

      if(RMin<0.||Rj/sqrtTwo<RMin)RMin=Rj/sqrtTwo;
      if(RMax<0.||Rj*sqrtTwo>RMax)RMax=Rj*sqrtTwo;
      nInt++;
      if(R>sdist+Rj && -R<sdist-Rj)
       {
        printf("%s -- Warning, New chart (%d) covers chart %d\n",RoutineName,chart,j);fflush(stdout);
        printf("      new chart %d is interval [%lf,%lf], chart %d is interval [%lf,%lf]\n",chart,-R,R,j,sdist-Rj,sdist+Rj);fflush(stdout);
        printf("      new chart radius %lf, old chart radius %lf\n",R,Rj);fflush(stdout);
        if(R>sdist+Rj)Rstar=sdist+Rj;
        if(-R<sdist-Rj)Rstar=Rj-sdist;
/**     chart=-1; */
        continue;
       }
      if(sdist-Rj<-R && R<sdist+Rj)
       {
        printf("%s -- Warning, chart %d covers new chart %d\n",RoutineName,j,chart);fflush(stdout);
        printf("      new chart %d is interval [%lf,%lf], chart %d is interval [%lf,%lf]\n",chart,-R,R,j,sdist-Rj,sdist+Rj);fflush(stdout);
        printf("      new chart radius %lf, old chart radius %lf\n",R,Rj);fflush(stdout);
        MFFreeKVector(nrm,e);
/*      A->nCharts--;*/
/*      chart=-1;*/
        continue;
       }

      if(MFChartIsSingular(A->chart[chart],e))
       {
        if(MFChartIsSingular(A->chart[j],e)&&mark==-1||MFChartReferenceNumber(A->chart[j],e)<mark)mark=MFChartReferenceNumber(A->chart[j],e);
       }else{
        if(mark==-1||MFChartReferenceNumber(A->chart[j],e)<mark)mark=MFChartReferenceNumber(A->chart[j],e);
       }
  
      onrm=.5*(dist+R*R-Rj*Rj);

#ifdef MFALLOWVERBOSE
      if(verbose){printf("   onrm_ij=%lf nrm_ij=",onrm);MFPrintKVector(stdout,nrm,e);printf("\n");fflush(stdout);}
#endif

      ij=MFAddHalfSpaceToAtlas(A,A->nGlobalIndex,chart,j,nrm,onrm,e);
      MFSubtractHalfSpaceFromChart(A->chart[chart],ij,A->nrm[ij],A->onrm[ij],e);
  
      MFFreeKVector(nrm,e);
      MFNSpaceDirection(A->NSpace,uj,u,dir,e);

#ifdef MFALLOWVERBOSE
      if(verbose){printf("   direction from %d to %d is ",j,chart);MFPrintNVector(stdout,dir,e);printf("\n");fflush(stdout);}
#endif

      nrm=MFCreateKVector(A->k,e);
      MFChartProjectVectorIntoTangentSpace(A->chart[j],dir,nrm,e);

#ifdef MFALLOWVERBOSE
      if(verbose){printf("   in TS %d is ",j);MFPrintKVector(stdout,nrm,e);printf("\n");fflush(stdout);}
#endif

      dist=MFKVDot(nrm,nrm,e);

#ifdef MFALLOWVERBOSE
      if(verbose){printf("   distance in TS is %lf\n",sqrt(dist));fflush(stdout);}
#endif
  
      onrm=.5*(dist+Rj*Rj-R*R);

#ifdef MFALLOWVERBOSE
      if(verbose){printf("   onrm_ji=%lf nrm_ji=",onrm);MFPrintKVector(stdout,nrm,e);printf("\n");fflush(stdout);}
#endif

      ji=MFAddHalfSpaceToAtlas(A,A->nGlobalIndex,j,chart,nrm,onrm,e);
      MFSubtractHalfSpaceFromChart(A->chart[j],ji,A->nrm[ji],A->onrm[ji],e);
      MFFreeKVector(nrm,e);

      if( (ibnd=MFChartGetPositionInBoundaryList(A->chart[j],e))>-1 )
       {
        P=MFChartPolytope(A->chart[j],e);
        if(MFPolytopeLargestRadiusOfVertex(P,e)<MFChartRadius(A->chart[j],e))
            MFAtlasRemoveChartFromBoundaryList(A,ibnd,e);
       }
  
      A->nGlobalIndex++;
     }

    MFFreeListOfIntersectingCharts(L,e);
   }

#ifdef MFALLOWVERBOSE
  if(0&&MFVerbose)
   {
    printf("   %d actually intersected\n",nInt);fflush(stdout);
   }

  if(verbose)
   {
    if(Rstar>=0.)printf("       Rstar=%lf\n",Rstar);
    printf("       update boundary list\n");fflush(stdout);
   }
#endif

  if(chart==-1)return -1;

  MFBinaryTreeAddChart(A->BTree,chart,BBCenter,MFChartRadius(A->chart[chart],e),e);
  P=MFChartPolytope(A->chart[chart],e);
#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("       add chart to boundary list\n");fflush(stdout);
    printf("         nvertices = %d\n",MFPolytopeNumberOfVertices(P,e));fflush(stdout);
    printf("         Smallest offset %d, Offset=%d\n",MFPolytopeSmallestVertexIndex(P,e),A->offset);fflush(stdout);
    printf("         Largest radius %lf, chart Radius=%d\n",MFPolytopeLargestRadiusOfVertex(P,e),MFChartRadius(A->chart[chart],e));fflush(stdout);
   }
#endif
  if(  MFPolytopeNumberOfVertices(P,e)>0
    && MFPolytopeSmallestVertexIndex(P,e)<A->offset 
    && MFPolytopeLargestRadiusOfVertex(P,e)>MFChartRadius(A->chart[chart],e)) MFAtlasAddChartToBoundaryList(A,chart,e);

  MFFreeNVector(dir,e);
  free(BBCenter);

  MFChartSetReferenceNumber(A->chart[chart],mark+1,e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf(" done %s\n",RoutineName);fflush(stdout);}
#endif

  return chart;
 }

/* --------------------------------------------------------- */

int MFAtlasPointOnBoundaryWOProject(MFAtlas A,MFNRegion Omega,MFKVector s, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAtlasPointOnBoundaryWOProject"};
  int v,l,n;
  int nChartsTried;
  int chart;
  int bchart;
  MFNVector up;
  MFNKMatrix Phi;
  double rv;
  double R;
  int rc;
  int abandon;
  int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s\n",RoutineName);fflush(stdout);}
#endif

  if(MFAtlasNumberOfCharts(A,e)>0)up=MFCloneNVector(MFAtlasChartCenter(A,0,e),e);
   else               up=NULL;
  nChartsTried=0;

#ifndef NOBOUNDARYLIST
  bchart=0;
  {printf("%s, there are %d charts on the boundary list\n",RoutineName,A->nBnd);fflush(stdout);}
  while(bchart<A->nBnd)
   {
    chart=A->bnd[bchart];
#else
  chart=0;
  while(chart<A->nCharts)
   {
#endif
tryagain:
    nChartsTried++;
    n=MFPolytopeNumberOfVertices(MFChartPolytope(A->chart[chart],e),e);
    v=0;
    while(v<n && MFChartRadius(A->chart[chart],e)>MFIMFGetRMin(A->M,e))
     {

#ifdef MFALLOWVERBOSE
      if(verbose){printf("  Trying vertex %d of chart %d\n",v,chart);fflush(stdout);}
#endif

      rv=MFPolytopeRadiusOfVertex(MFChartPolytope(A->chart[chart],e),v,e);
      R=MFChartRadius(A->chart[chart],e);

      if(rv>=R)
       {

#ifdef MFALLOWVERBOSE
        if(verbose){printf("  Exterior!\n");fflush(stdout);}
#endif

        MFPolytopeVertex(MFChartPolytope(A->chart[chart],e),v,s,e);
        MFKVScale(R/rv,s,e);
        MFChartPointInTangentSpace(A->chart[chart],s,up,e);

        if(MFNRegionInterior(Omega,up,e))
         {

#ifdef MFALLOWVERBOSE
          if(verbose){printf("    In Omega!\n");fflush(stdout);}
#endif

          MFFreeNVector(up,e);

#ifdef MFALLOWVERBOSE
          if(verbose)
           {
	    printf("    found after examining %d chart\n",nChartsTried);fflush(stdout);
            printf("Chart %d is now \n",chart);
            MFPrintChart(stdout,A->chart[chart],e);
            fflush(stdout);
           }
#endif

          return chart;
         }else{

#ifdef MFALLOWVERBOSE
          if(verbose){printf("    Not in Omega!\n");fflush(stdout);}
#endif

         }
       }
      v++;
     }
    MFAtlasRemoveChartFromBoundaryList(A,0,e);
    if(bchart>=A->nBnd)
     {
      bchart=0;
      chart=-1;
     }
   }
  if(up!=NULL)MFFreeNVector(up,e);

  return -1;
 }

void MFAtlasSetK(MFAtlas A,int k, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAtlasSetK"};

#ifdef MFNOCONFIDENCE
  if(A==NULL)
   {
    sprintf(MFAtlasErrorMsg,"Pointer to Atlas (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(A->nCharts>0&&k!=A->k)
   {
    sprintf(MFAtlasErrorMsg,"Changing k after Atlas has charts may lead to errors later.");
    MFSetError(e,4,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
   }
#endif

  A->k=k;
  MFAddInitialHyperplanes(A,e);

  return;
 }

void MFAtlasSetN(MFAtlas A,int n, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAtlasSetK"};

#ifdef MFNOCONFIDENCE
  if(A==NULL)
   {
    sprintf(MFAtlasErrorMsg,"Pointer to Atlas (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(A->nCharts>0&&n!=A->n)
   {
    sprintf(MFAtlasErrorMsg,"Changing n after Atlas has charts may lead to errors later.");
    MFSetError(e,4,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
   }
#endif

  A->n=n;
  MFAddInitialHyperplanes(A,e);

  return;
 }

void MFWriteChartCenters(FILE *fid,MFAtlas A, MFErrorHandler e)
 {
  static char RoutineName[]={"MFWriteChartCenters"};
  int i;
  int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s\n");fflush(stdout);}
#endif

  fprintf(fid,"%d\n",A->n);
  fprintf(fid,"%d\n",A->k);

/* List of Charts */

  fprintf(fid,"%d\n",A->nCharts);
  fprintf(fid,"%d\n",A->mCharts);
  for(i=0;i<A->nCharts;i++)
   {
    MFWriteNVector(fid,MFChartCenter(A->chart[i],e),e);
    MFWriteNKMatrix(fid,MFChartTangentSpace(A->chart[i],e),e);
    fprintf(fid,"%lf\n",MFChartRadius(A->chart[i],e));
   }

  return;
 }


int MFReadChartCenters(FILE *fid,MFNVector **c,MFNKMatrix **t,double**R, MFErrorHandler e)
 {
  static char RoutineName[]={"MFReadChartCenters"};
  int i,n,k;
  int nCharts;
  int mCharts;
  int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose){scanf("%s\n");fflush(stdout);}
#endif

  fscanf(fid,"%d\n",&n);
  fscanf(fid,"%d\n",&k);

/* List of Charts */

  fscanf(fid,"%d\n",&nCharts);
  fscanf(fid,"%d\n",&mCharts);
  *c=(MFNVector*)realloc((void*)(*c),nCharts*sizeof(MFNVector));

#ifndef MFNOSAFETYNET
  if(*c==NULL)
   {
    sprintf(MFAtlasErrorMsg,"Out of memory trying to allocate %d bytes",nCharts*sizeof(MFNVector));
    MFSetError(e,12,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return -1;
   }
#endif

  *t=(MFNKMatrix*)realloc((void*)(*t),nCharts*sizeof(MFNKMatrix));

#ifndef MFNOSAFETYNET
  if(*t==NULL)
   {
    sprintf(MFAtlasErrorMsg,"Out of memory trying to allocate %d bytes",nCharts*sizeof(MFNKMatrix));
    MFSetError(e,12,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return -1;
   }
#endif

  *R=(double*)realloc((void*)(*R),nCharts*sizeof(double));

#ifndef MFNOSAFETYNET
  if(*R==NULL)
   {
    sprintf(MFAtlasErrorMsg,"Out of memory trying to allocate %d bytes",nCharts*sizeof(double));
    MFSetError(e,12,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return -1;
   }
#endif

  for(i=0;i<nCharts;i++)
   {
    (*c)[i]=MFReadNVector(fid,e);
    (*t)[i]=MFReadNKMatrix(fid,e);
    fscanf(fid,"%lf\n",&((*R)[i]));
   }

  return nCharts;
 }

void MFAtlasPageOutChartsNotNearBoundary(MFAtlas A, int dumpToPlotFile, int dumpToCenterFile, char *FileName, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAtlasPageOutChartsNotNearBoundary"};
  int i,ii,j,n,l,m;
  int nmoved;
  MFPolytope P,Q;
  int verbose=0;
  char name[10000];

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("PAGING!!! (%s) -- %d charts in Atlas\n",RoutineName,A->nCharts);fflush(stdout);
    MFPrintChartSummary(A,e);
   }
#endif

  for(i=0;i<A->nCharts;i++)
    MFChartSetNearBoundary(A->chart[i],0,e);

  for(i=0;i<A->nCharts;i++)
   {
    if(!MFChartPaged(A->chart[i],e))
     {
      if(MFNVGetIndex2(MFChartCenter(A->chart[i],e),e)==-99||MFNVGetIndex2(MFChartCenter(A->chart[i],e),e)==-98)
       {
        MFChartSetNearBoundary(A->chart[i],1,e);

        P=MFChartPolytope(A->chart[i],e);
        n=MFPolytopeNumberOfFaces(P,e);
        for(j=0;j<n;j++)
         {
          ii=MFPolytopeFaceIndex(P,j,e);
          if(ii>-1 && ii< A->nHS
             && A->rightChart[ii]>-1 && A->rightChart[ii]<A->nCharts)
           {
            MFChartSetNearBoundary(A->chart[A->rightChart[ii]],1,e);
            Q=MFChartPolytope(A->chart[A->rightChart[ii]],e);
            if(Q!=NULL)
             {
              m=MFPolytopeNumberOfFaces(Q,e);
              for(l=0;l<m;l++)
               {
                ii=MFPolytopeFaceIndex(Q,l,e);
                if(ii>-1 && ii< A->nHS && ii >= A->offset
                   && A->rightChart[ii]>-1 && A->rightChart[ii]<A->nCharts)
                 {
                  MFChartSetNearBoundary(A->chart[A->rightChart[ii]],1,e);
                 }
               }
             }
           }
         }

       }
     }
   }

#ifdef MFALLOWVERBOSE
  if(verbose){fprintf(stdout,"PAGING mark boundary charts(%d on boundary)!!!\n",A->nBnd);fflush(stdout);}
#endif

  for(i=0;i<A->nBnd;i++)
   {
    MFChartSetNearBoundary(A->chart[A->bnd[i]],1,e);
    P=MFChartPolytope(A->chart[A->bnd[i]],e);
    n=MFPolytopeNumberOfFaces(P,e);
    for(j=0;j<n;j++)
     {
      ii=MFPolytopeFaceIndex(P,j,e);
      if(ii>-1 && ii< A->nHS 
         && A->rightChart[ii]>-1 && A->rightChart[ii]<A->nCharts)
       {
        MFChartSetNearBoundary(A->chart[A->rightChart[ii]],1,e);
        Q=MFChartPolytope(A->chart[A->rightChart[ii]],e);
        if(Q!=NULL)
         {
          m=MFPolytopeNumberOfFaces(Q,e);
          for(l=0;l<m;l++)
           {
            ii=MFPolytopeFaceIndex(Q,l,e);
            if(ii>-1 && ii< A->nHS && ii >= A->offset
               && A->rightChart[ii]>-1 && A->rightChart[ii]<A->nCharts)
             {
              MFChartSetNearBoundary(A->chart[A->rightChart[ii]],1,e);
             }
           }
         }
       }
     }
   }
#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    m=0;
    for(i=0;i<A->nCharts;i++)if(MFChartNearBoundary(A->chart[i],e))m++;
    fprintf(stdout,"       %d charts marked as near boundary\n",m);fflush(stdout);
   }
#endif

#ifdef MFALLOWVERBOSE
  if(verbose){fprintf(stdout,"PAGING those not marked\n");fflush(stdout);}
#endif

#ifdef REALLYPAGE

#ifdef MFALLOWVERBOSE
  if(verbose){printf("page those not marked\n");fflush(stdout);}
#endif

  nmoved=0;
  for(i=0;i<A->nCharts;i++)
   {
    if(!MFChartNearBoundary(A->chart[i],e))
     {
      if(dumpToPlotFile && A->plotFileFid==NULL)
       {
        strcpy(name,FileName);
        strcat(name,".plotfile");
        A->plotFileFid=fopen(name,"w+b");
        if(A->plotFileFid==NULL)abort();

#ifdef MFALLOWVERBOSE
        if(verbose){fprintf(stdout,"Opened plotfile, 0x%8.8x\n",A->plotFileFid);fflush(stdout);}
#endif

        fprintf(A->plotFileFid,"Dimension of vertices, %d\n",MFIMFProjectToDraw(A->M,MFAtlasChartCenter(A,0,e),NULL,e));
        fprintf(A->plotFileFid,"Dimension of manifold, %d\n",A->k);
       }
      if(dumpToCenterFile && A->centerFileFid==NULL)
       {
        strcpy(name,FileName);
        strcat(name,".centers");
        A->centerFileFid=fopen(name,"w+b");
        if(A->centerFileFid==NULL)abort();

#ifdef MFALLOWVERBOSE
        if(verbose){fprintf(stdout,"Opened centerfile, 0x%8.8x\n",A->centerFileFid);fflush(stdout);}
#endif

        fprintf(A->centerFileFid,"%d\n",A->n);
        fprintf(A->centerFileFid,"%d\n",A->k);
        fprintf(A->centerFileFid,"nCharts goes here\n");
        fprintf(A->centerFileFid,"mCharts goes here\n");
       }

#ifdef MFALLOWVERBOSE
      if(verbose){printf(" chart %d is not marked\n",i);fflush(stdout);}
#endif

      if(!MFChartPaged(A->chart[i],e))
       {
        if(dumpToPlotFile)
         {
          MFNVector col;

#ifdef MFALLOWVERBOSE
          if(verbose){fprintf(stdout,"Paging Chart %d\n",i);fflush(stdout);}
#endif

          MFWriteChartToPlotFile(A->plotFileFid,A,A->chart[i],i,e);
          fflush(A->plotFileFid);
         }
        if(dumpToCenterFile)
         {
          MFWriteNVector(A->centerFileFid,MFChartCenter(A->chart[i],e),e);
          MFWriteNKMatrix(A->centerFileFid,MFChartTangentSpace(A->chart[i],e),e);
          fprintf(A->centerFileFid,"%lf\n",MFChartRadius(A->chart[i],e));
          fflush(A->centerFileFid);
         }
        if(MFChartPageOut(A->chart[i],A->pageFileFid,A->pageFileIndex,e))A->pageFileIndex++;
        nmoved++;
       }
     }
   }
  A->nPagedOut+=nmoved;

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("done PAGING!!! (%s) moved %d out. A total of %d are out, out of %d\n",RoutineName,nmoved,A->nPagedOut,A->nCharts);
    fflush(stdout);
    MFPrintChartSummary(A,e);
   }
#endif

#endif /* REALLYPAGE */

  return;
 }

int MFAtlasIsChartSingular(MFAtlas thisAtlas,int chart, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAtlasIsChartSingular"};

#ifdef MFNOCONFIDENCE
  if(thisAtlas==NULL)
   {
    sprintf(MFAtlasErrorMsg,"Pointer to Atlas (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  return MFChartIsSingular(thisAtlas->chart[chart],e);
 }

int MFAtlasIsChartNearBoundary(MFAtlas thisAtlas,int chart, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAtlasIsChartNearBoundary"};

#ifdef MFNOCONFIDENCE
  if(thisAtlas==NULL)
   {
    sprintf(MFAtlasErrorMsg,"Pointer to Atlas (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  return MFChartNearBoundary(thisAtlas->chart[chart],e);
 }

int MFAtlasNumberOfNeighborsOfChart(MFAtlas thisAtlas,int chart, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAtlasNumberOfNeighborsOfChart"};

#ifdef MFNOCONFIDENCE
  if(thisAtlas==NULL)
   {
    sprintf(MFAtlasErrorMsg,"Pointer to Atlas (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(chart<0 || chart>=thisAtlas->nCharts)
   {
    sprintf(MFAtlasErrorMsg,"Chart number %d (argument 2) is invalid, must be in [0,%d)",chart,thisAtlas->nCharts);
    MFSetError(e,4,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  return MFPolytopeNumberOfFaces(MFChartPolytope(thisAtlas->chart[chart],e),e);
 }

int MFAtlasChartNeighbor(MFAtlas thisAtlas,int chart,int neighbor, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAtlasChartNeighbor"};
  int n;

#ifdef MFNOCONFIDENCE
  if(thisAtlas==NULL)
   {
    sprintf(MFAtlasErrorMsg,"Pointer to Atlas (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(chart<0 || chart>= thisAtlas->nCharts)
   {
    sprintf(MFAtlasErrorMsg,"Chart number %d (argument 2) is invalid, must be in [0,%d)",chart,thisAtlas->nCharts);
    MFSetError(e,4,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    return;
   }

  n=MFChartNumberOfFaces(thisAtlas->chart[chart],e);
  if(neighbor<0 || neighbor>=n)
   {
    sprintf(MFAtlasErrorMsg,"Neighbor number %d (argument 3) is invalid, must be in [0,%d)",neighbor,n);
    MFSetError(e,4,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  n=MFPolytopeFaceIndex(MFChartPolytope(thisAtlas->chart[chart],e),neighbor,e);

  if(thisAtlas->leftChart[n]==chart)return thisAtlas->rightChart[n];
   else return thisAtlas->leftChart[n];
 }

void MFWriteChartToPlotFile(FILE *fid, MFAtlas A, MFChart thisChart, int chart, MFErrorHandler e)
 {
  static char RoutineName[]={"MFWriteChartToPlotFile"};
  int d;
  double *x;
  MFKVector s;
  double *ps;
  MFNVector u;
  int i,j,k,l;
  int n,m,ne,ie;
  int nv;
  int c;
  int *index;
  int *indx;
  int verbose=0;
  MFNKMatrix Phi;
  MFNVector col;
  int bnd,sing;
  double *dv;
  MFPolytope P;

  if(MFChartPaged(thisChart,e))
   {
#ifdef MFNOCONFIDENCE
    sprintf(MFAtlasErrorMsg,"Chart (argument 3) has been paged out");
    MFSetError(e,4,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
#endif
    return;
   }

  k=MFChartK(thisChart,e);
  if(MFPolytopeNumberOfVertices(MFChartPolytope(thisChart,e),e)<=k)return;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s\n",RoutineName);fflush(stdout);}
#endif

  if(fid==NULL){printf("%s - fid is NULL %s(%d)\n",RoutineName,__LINE__,__FILE__);fflush(stdout);abort();}
  if(A==NULL){printf("%s - A is NULL %s(%d)\n",RoutineName,__LINE__,__FILE__);fflush(stdout);abort();}
  if(thisChart==NULL){printf("%s - Chart is NULL %s(%d)\n",RoutineName,__LINE__,__FILE__);fflush(stdout);abort();}

  d=MFIMFProjectToDraw(A->M,MFChartCenter(thisChart,e),NULL,e);
  if(A->M==NULL){printf("%s - A->M is NULL %s(%d)\n",RoutineName,__LINE__,__FILE__);fflush(stdout);abort();}

  s=MFCreateKVector(MFChartK(thisChart,e),e);
  ps=MFKV_CStar(s,e);
  u=MFCloneNVector(MFChartCenter(thisChart,e),e);
  x=(double*)malloc(d*sizeof(double));

#ifndef MFNOSAFETYNET
  if(x==NULL)
   {
    sprintf(MFAtlasErrorMsg,"Out of memory trying to allocate %d bytes",d*sizeof(double));
    MFSetError(e,12,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  for(j=0;j<A->nclipf;j++)
   {
    n=MFPolytopeNumberOfVertices(MFChartPolytope(thisChart,e),e);
    if(n>0)
     {
      dv=(double*)malloc(n*sizeof(double));

#ifndef MFNOSAFETYNET
      if(dv==NULL)
       {
        sprintf(MFAtlasErrorMsg,"Out of memory trying to allocate %d bytes",n*sizeof(double));
        MFSetError(e,12,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
        MFErrorHandlerOutOfMemory(e);
        return;
       }
#endif

      for(i=0;i<n;i++)
       {
        MFPolytopeVertex(MFChartPolytope(thisChart,e),i,s,e);
        MFChartPointInTangentSpace(thisChart,s,u,e);
        dv[i]=A->clipf[j](u,e);
       }
      MFClipPolytope(MFChartPolytope(thisChart,e),A->clipIndx[j],dv,0,e);
      if(dv!=NULL)free(dv);
     }
   }

  if(MFPolytopeNumberOfVertices(MFChartPolytope(thisChart,e),e)<=k)
   {
    MFFreeKVector(s,e);
    MFFreeNVector(u,e);
    free(x);
    return;
   }

  Phi=MFChartTangentSpace(thisChart,e);

#ifdef MFNOCONFIDENCE
  if(Phi==NULL)
   {
    sprintf(MFAtlasErrorMsg,"Tangent space is NULL");
    MFSetError(e,12,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  col=MFMColumn(Phi,0,e);

#ifdef MFNOCONFIDENCE
  if(col==NULL)
   {
    sprintf(MFAtlasErrorMsg,"First tangent vector is NULL");
    MFSetError(e,12,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
	    return;
   }
#endif

  MFFreeNVector(col,e);
  if(k>1)
   {
    col=MFMColumn(Phi,1,e);

#ifdef MFNOCONFIDENCE
  if(col==NULL)
   {
    sprintf(MFAtlasErrorMsg,"Second tangent vector is NULL");
    MFSetError(e,12,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

    MFFreeNVector(col,e);
   }

  indx=(int*)malloc(k*sizeof(int));

#ifndef MFNOSAFETYNET
      if(indx==NULL)
       {
        sprintf(MFAtlasErrorMsg,"Out of memory trying to allocate %d bytes",k*sizeof(int));
        MFSetError(e,12,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
        MFErrorHandlerOutOfMemory(e);
        return;
       }
#endif

  ne=0;
  n=MFPolytopeNumberOfVertices(MFChartPolytope(thisChart,e),e);
  for(i=0;i<n-1;i++)
    for(j=i+1;j<n;j++)
     {
      if(MFPolytopeIntersectIndexSets(MFChartPolytope(thisChart,e),i,j,indx,e)>=k-1)ne++;
     }

  m=MFPolytopeNumberOfFaces(MFChartPolytope(thisChart,e),e);
  bnd=MFPolytopeLargestRadiusOfVertex(MFChartPolytope(thisChart,e),e)>MFChartRadius(thisChart,e);
  sing=MFChartIsSingular(thisChart,e);


/* For Oriol */
  if(0){
   double *y;
   double dist;
   printf("Writing chart to plotfile\n");fflush(stdout);

   y=(double*)malloc(d*sizeof(double));

   printf("Checking the distance to the center from each vertex \n");fflush(stdout);

   printf("Project the center\n");fflush(stdout);
   for(j=0;j<k;j++)ps[i]=0.;
   MFChartPointInTangentSpace(thisChart,s,u,e);
   MFIMFProjectToDraw(A->M,u,y,e);
   for(i=0;i<n;i++)
    {
   printf("Project vertex %d out of %d, dimension=%d\n",i,n,d);fflush(stdout);
     MFPolytopeVertex(MFChartPolytope(thisChart,e),i,s,e);
     MFChartPointInTangentSpace(thisChart,s,u,e);
     MFIMFProjectToDraw(A->M,u,x,e);
     dist=0;
     for(j=0;j<d;j++)dist+=(x[j]-y[j])*(x[j]-y[j]);
     dist=sqrt(dist);
   printf("   dist to center=%lf, radius=%lf, dist/radius=%lf\n",dist,MFChartRadius(thisChart,e),dist/MFChartRadius(thisChart,e));fflush(stdout);
     if(dist>4*MFChartRadius(thisChart,e))
      {
       printf("   dist to center is more than 4*radius, skipping writing this chartf\n");fflush(stdout);
/*     free(y);
 *     return;
 */
      }
    }
   free(y);
  }

  nv=n+1;
  fprintf(fid,"Polyhedron %d, R=%lf, %d vertices, %d edges, %d faces, boundary %d, singular %d\n",chart,MFChartRadius(thisChart,e),nv,ne,m,bnd,sing);
  for(i=0;i<n;i++)
   {
    MFPolytopeVertex(MFChartPolytope(thisChart,e),i,s,e);
    for(j=0;j<k;j++)ps[j]=(A->expFactor)*ps[j]; /* This allows for eliminating gaps between plates */
    MFChartPointInTangentSpace(thisChart,s,u,e);
    MFNVSetIndex(u,MFNVGetIndex(MFChartCenter(thisChart,e),e),e);
    MFNVSetIndex2(u,MFNVGetIndex2(MFChartCenter(thisChart,e),e),e);
    MFIMFProjectToDraw(A->M,u,x,e);
    for(j=0;j<d;j++)if(x[j]!=x[j]){printf("In %s, x[%d] is NaN!\n",RoutineName,j);fflush(stdout);abort();}
    fprintf(fid,"Vertex %d (%lf",i,x[0]);
    for(j=1;j<d;j++)fprintf(fid,",%lf",x[j]);
    fprintf(fid,")");
    m=MFPolytopeVertexNumberOfIndices(MFChartPolytope(thisChart,e),i,e);
    index=MFPolytopeVertexIndexSet(MFChartPolytope(thisChart,e),i,e);
    fprintf(fid,", %d [%d",m,index[0]);
    for(j=1;j<m;j++)fprintf(fid,",%d",index[j]);
    fprintf(fid,"]\n");
   }

/* Center */

  MFIMFProjectToDraw(A->M,MFChartCenter(thisChart,e),x,e);
  fprintf(fid,"Vertex %d (%lf",i,x[0]);
  for(j=1;j<d;j++)fprintf(fid,",%lf",x[j]);
  fprintf(fid,")");
  fprintf(fid,", %d [ ]\n",0);

  MFFreeNVector(u,e);
  MFFreeKVector(s,e);
  free(x);
  fflush(fid);

  ie=0;
  for(i=0;i<n-1;i++)
   {
    for(j=i+1;j<n;j++)
     {
      if((m=MFPolytopeIntersectIndexSets(MFChartPolytope(thisChart,e),i,j,indx,e))>=k-1)
       {
        fprintf(fid,"Edge %d (%d,%d)",ie,i,j);
        fprintf(fid,", %d [",m);
        if(m>0)fprintf(fid,"%d",indx[0]);
        for(l=1;l<m;l++)fprintf(fid,",%d",indx[l]);
        fprintf(fid,"]\n");
        ie++;
       }
     }
   }
  fflush(fid);
  free(indx);

  m=MFPolytopeNumberOfFaces(MFChartPolytope(thisChart,e),e);
  for(i=0;i<m;i++)
   {
    j=MFPolytopeFaceIndex(MFChartPolytope(thisChart,e),i,e);
    if(j<A->offset)c=-1;
    else if(A->leftChart[j]==chart)c=A->rightChart[j];
     else c=A->leftChart[j];
    fprintf(fid,"Face %d neighbor %d\n",j,c);
   }

  fflush(fid);

Return:

#ifdef MFALLOWVERBOSE
  if(verbose){printf("done %s\n",RoutineName);fflush(stdout);}
#endif

  return;
 }

void MFAtlasPageOutAllCharts(MFAtlas A, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAtlasPageOutAllCharts"};
  int i,ii,j,n,l,m;
  int nmoved;
  MFPolytope P,Q;
  int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose){fprintf(stdout,"PAGING!!! (%s) -- %d charts in Atlas\n",RoutineName,A->nCharts);fflush(stdout);}
#endif

#ifdef REALLYPAGE

#ifdef MFALLOWVERBOSE
  if(verbose){printf("       Page all not paged out\n");fflush(stdout);}
#endif

  nmoved=0;
  for(i=0;i<A->nCharts;i++)
   {
    if(A->plotFileFid==NULL)
     {
      A->plotFileFid=fopen("Atlas.plotfile","w+b");
      if(A->plotFileFid==NULL)abort();

#ifdef MFALLOWVERBOSE
      if(verbose){fprintf(stdout,"Opened plotfile, 0x%8.8x\n",A->plotFileFid);fflush(stdout);}
#endif

      fprintf(A->plotFileFid,"Dimension of vertices, %d\n",MFIMFProjectToDraw(A->M,MFAtlasChartCenter(A,0,e),NULL,e));
      fprintf(A->plotFileFid,"Dimension of manifold, %d\n",A->k);
     }
    if(A->centerFileFid==NULL)
     {
      A->centerFileFid=fopen("Atlas.centers","w+b");

#ifdef MFNOCONFIDENCE
      if(A->centerFileFid==NULL)
       {
        sprintf(MFAtlasErrorMsg,"Unable to open the centerfile \"%s\" for output\n","Atlas.centers");
        MFSetError(e,12,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
        return;
       }
#endif

#ifdef MFALLOWVERBOSE
      if(verbose){fprintf(stdout,"Opened centerfile, 0x%8.8x\n",A->centerFileFid);fflush(stdout);}
#endif

      fprintf(A->centerFileFid,"%d\n",A->n);
      fprintf(A->centerFileFid,"%d\n",A->k);
      fprintf(A->centerFileFid,"nCharts goes here\n");
      fprintf(A->centerFileFid,"mCharts goes here\n");
     }
    if(!MFChartPaged(A->chart[i],e))
     {
      MFWriteChartToPlotFile(A->plotFileFid,A,A->chart[i],i,e);
      MFWriteNVector(A->centerFileFid,MFChartCenter(A->chart[i],e),e);
      MFWriteNKMatrix(A->centerFileFid,MFChartTangentSpace(A->chart[i],e),e);
      fprintf(A->centerFileFid,"%lf\n",MFChartRadius(A->chart[i],e));
      if(MFChartPageOut(A->chart[i],A->pageFileFid,A->pageFileIndex,e))A->pageFileIndex++;
      nmoved++;
      fflush(A->plotFileFid);
      fflush(A->centerFileFid);
     }
   }
  A->nPagedOut+=nmoved;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("done PAGING!!! (%s) moved %d out. A total of %d are out, out of %d\n",RoutineName,nmoved,A->nPagedOut,A->nCharts);fflush(stdout);}
#endif

#endif

  return;
 }

void MFAtlasWriteOutChartsNotPaged(MFAtlas A, int dumpToPlotFile, int dumpToCenterFile, char *FileName, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAtlasWriteOutChartsNotPaged"};
  int i,ii,j,n,l,m;
  int nmoved;
  MFPolytope P,Q;
  int verbose=0;
  char name[10000];

#ifdef MFALLOWVERBOSE
  if(verbose){fprintf(stdout,"Dumping charts not paged!!! (%s) -- %d charts in Atlas\n",RoutineName,A->nCharts);fflush(stdout);}
#endif

  nmoved=0;
  for(i=0;i<A->nCharts;i++)
   {
    if(dumpToPlotFile && A->plotFileFid==NULL)
     {
      strcpy(name,FileName);
      strcat(name,".plotfile");
      A->plotFileFid=fopen(name,"w+b");
      if(A->plotFileFid==NULL)abort();

#ifdef MFALLOWVERBOSE
      if(verbose){fprintf(stdout,"Opened plotfile, 0x%8.8x\n",A->plotFileFid);fflush(stdout);}
#endif

      fprintf(A->plotFileFid,"Dimension of vertices, %d\n",MFIMFProjectToDraw(A->M,MFAtlasChartCenter(A,0,e),NULL,e));
      fprintf(A->plotFileFid,"Dimension of manifold, %d\n",A->k);
     }
    if(dumpToCenterFile && A->centerFileFid==NULL)
     {
      strcpy(name,FileName);
      strcat(name,".centers");
      A->centerFileFid=fopen(name,"w+b");
      if(A->centerFileFid==NULL)abort();

#ifdef MFALLOWVERBOSE
      if(verbose){fprintf(stdout,"Opened centerfile, 0x%8.8x\n",A->centerFileFid);fflush(stdout);}
#endif

      fprintf(A->centerFileFid,"%d\n",A->n);
      fprintf(A->centerFileFid,"%d\n",A->k);
      fprintf(A->centerFileFid,"nCharts goes here\n");
      fprintf(A->centerFileFid,"mCharts goes here\n");
     }
    if(!MFChartPaged(A->chart[i],e))
     {
      if(dumpToPlotFile)
       {
        MFWriteChartToPlotFile(A->plotFileFid,A,A->chart[i],i,e);
        fflush(A->plotFileFid);
       }
      if(dumpToCenterFile)
       {
        MFWriteNVector(A->centerFileFid,MFChartCenter(A->chart[i],e),e);
        MFWriteNKMatrix(A->centerFileFid,MFChartTangentSpace(A->chart[i],e),e);
        fprintf(A->centerFileFid,"%lf\n",MFChartRadius(A->chart[i],e));
        fflush(A->centerFileFid);
       }

#ifdef MFALLOWVERBOSE
      if(verbose){printf("Chart %d/%d written to plotfile and center file\n",i,A->nCharts);}
#endif

      nmoved++;
     }else{

#ifdef MFALLOWVERBOSE
      if(verbose){printf("Chart %d/%d has already been written to plotfile and center file\n",i,A->nCharts);}
#endif

     }
   }

#ifdef MFALLOWVERBOSE
  if(verbose){printf("done Dumping!!! (%s) wrote %d out. A total of %d are out, out of %d\n",RoutineName,nmoved,A->nPagedOut,A->nCharts);fflush(stdout);}
#endif

  return;
 }

void MFAtlasReleaseSingularCharts(MFAtlas A, int type0, int type1, MFErrorHandler e)
 {
  int chart,n,v;

  for(chart=0;chart<A->nCharts;chart++)
   {
    if(!MFChartPaged(A->chart[chart],e))
     {
      if(MFNVGetIndex2(MFChartCenter(A->chart[chart],e),e)==type0)
         MFNVSetIndex2(MFChartCenter(A->chart[chart],e),type1,e);
       else if(MFNVGetIndex2(MFChartCenter(A->chart[chart],e),e)==type1)
         MFNVSetIndex2(MFChartCenter(A->chart[chart],e),type0,e);
     }
   }
  return;
 }

int MFAtlasGetSingularChartWithBoundary(MFAtlas A, MFNRegion Omega, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAtlasGetSingularChartWithBoundary"};
  int chart,n,v;
  double R,rv;
  MFKVector s;
  MFNVector up;
  int pass;
  int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s\n",RoutineName);fflush(stdout);}
#endif

#ifdef MFNOCONFIDENCE
  if(A==NULL)
   {
    sprintf(MFAtlasErrorMsg,"MFAtlas (first argument is NULL");
    MFSetError(e,12,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    return -1;
   }
#endif

  s=MFCreateKVector(A->k,e);
  up=NULL;
  for(chart=0;chart<A->nCharts;chart++)
   if(up==NULL&&!MFChartPaged(A->chart[chart],e))up=MFCloneNVector(MFChartCenter(A->chart[chart],e),e);
  if(up==NULL)return -1;

  for(pass=0;pass<2;pass++)
   {
    for(chart=0;chart<A->nCharts;chart++)
     {

#ifdef MFALLOWVERBOSE
      if(verbose){printf("chart %d/%d\n",chart,A->nCharts);fflush(stdout);}
#endif

      if(!MFChartPaged(A->chart[chart],e)
        &&MFChartIsSingular(A->chart[chart],e)&&MFNVGetIndex2(MFChartCenter(A->chart[chart],e),e)==-99
        &&MFChartRadius(A->chart[chart],e)>MFIMFGetRMin(A->M,e))
       {

#ifdef MFALLOWVERBOSE
        if(verbose){printf("  Is Singular\n");fflush(stdout);}
#endif

        n=MFPolytopeNumberOfVertices(MFChartPolytope(A->chart[chart],e),e);
        for(v=0;v<n;v++)
         {
          rv=MFPolytopeRadiusOfVertex(MFChartPolytope(A->chart[chart],e),v,e);
          R=MFChartRadius(A->chart[chart],e);
          if(rv>=R && MFPolytopeGetVertexMark(MFChartPolytope(A->chart[chart],e),v,e)==0)
           {

#ifdef MFALLOWVERBOSE
            if(verbose){printf("  Vertex %d is exterior\n",v);fflush(stdout);}
#endif

            MFPolytopeVertex(MFChartPolytope(A->chart[chart],e),v,s,e);
            MFKVScale(R/rv,s,e);
            MFChartPointInTangentSpace((A->chart[chart]),s,up,e);
            if(MFNRegionInterior(Omega,up,e))
             {
              MFFreeKVector(s,e);
              MFFreeNVector(up,e);

#ifdef MFALLOWVERBOSE
              if(verbose){printf("done %s\n",RoutineName);fflush(stdout);}
#endif

              printf("*** Vertex not in region, exempting vertex %d of chart %d\n",v,chart);fflush(stdout);
              MFPolytopeSetVertexMark(MFChartPolytope(A->chart[chart],e),v,1,e);
              return chart;
             }
           }
         }
#ifdef MFALLOWVERBOSE
        if(verbose){printf("  No exterior vertices\n");fflush(stdout);}
#endif

       }else{

#ifdef MFALLOWVERBOSE
        if(verbose){printf("  Not singular\n");fflush(stdout);}
#endif

       }
     }
    MFAtlasReleaseSingularCharts(A,-98,-99,e);
   }

  MFFreeKVector(s,e);
  MFFreeNVector(up,e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("done %s\n",RoutineName);fflush(stdout);}
#endif

  return -1;
 }

MFNVector MFAtlasGetPointOnBoundaryChart(MFAtlas A, MFNRegion Omega, int chart, double t0, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAtlasGetPointOnBoundaryChart"};
  int n,v;
  int rc;
  double R,rv,t;
  MFKVector s;
  MFNVector up;
  MFNVector u;
  int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s, chart %d\n",RoutineName,chart);fflush(stdout);}
#endif

#ifdef MFNOCONFIDENCE
  if(A==NULL)
   {
    sprintf(MFAtlasErrorMsg,"MFAtlas (first argument is NULL");
    MFSetError(e,12,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    return NULL;
   }

  if(chart<0||chart>=A->nCharts)
   {
    sprintf(MFAtlasErrorMsg,"Chart %d is bad, must be in range [0,%d)",chart,A->nCharts);
    MFSetError(e,12,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    return NULL;
   }

  if(!MFChartIsSingular(A->chart[chart],e))
   {
    sprintf(MFAtlasErrorMsg,"Chart %d is bad, must be a singular chart",chart);
    MFSetError(e,12,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    return NULL;
   }

  if(MFChartRadius(A->chart[chart],e)<MFIMFGetRMin(A->M,e))
   {
    sprintf(MFAtlasErrorMsg,"Chart %d is bad, radius is too small (%lf<%lf)",chart,
                                              MFChartRadius(A->chart[chart],e),MFIMFGetRMin(A->M,e));
    MFSetError(e,12,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    return NULL;
   }
#endif

  s=MFCreateKVector(A->k,e);
  up=MFCloneNVector(MFChartCenter(A->chart[0],e),e);
  u=MFCloneNVector(up,e);

  n=MFPolytopeNumberOfVertices(MFChartPolytope(A->chart[chart],e),e);
  for(v=0;v<n;v++)
   {
    rv=MFPolytopeRadiusOfVertex(MFChartPolytope(A->chart[chart],e),v,e);
    R=MFChartRadius(A->chart[chart],e);

#ifdef MFALLOWVERBOSE
    if(verbose){printf(" vertex %d/%d  radius=%lf/%lf\n",v,n,rv,R);fflush(stdout);}
#endif

    if(rv>=R)
     {

#ifdef MFALLOWVERBOSE
      if(verbose){printf(" exterior\n");fflush(stdout);}
#endif

      MFPolytopeVertex(MFChartPolytope(A->chart[chart],e),v,s,e);
      MFKVScale(R/rv,s,e);
      MFChartPointInTangentSpace(A->chart[chart],s,up,e);
      if(MFNRegionInterior(Omega,up,e))
       {

#ifdef MFALLOWVERBOSE
        if(verbose){printf(" in Omega\n");fflush(stdout);}
#endif

        t=t0;
        while(t>1.e-4)
         {
          rc=MFChartEvaluate(A->chart[chart],s,u,e);

#ifdef MFALLOWVERBOSE
          if(verbose){printf(" rc from MFChartEvaluate=%d, t=%lf, s=",rc,t);MFPrintKVector(stdout,s,e);printf("\n");fflush(stdout);}
#endif
          if(rc)
           {
            MFFreeKVector(s,e);
            MFFreeNVector(up,e);

#ifdef MFALLOWVERBOSE
            if(verbose){printf("done %s\n",RoutineName);fflush(stdout);}
#endif
            return u;
           }else{
            t=t*.8;
            MFKVScale(.8,s,e);
           }
         }
       }else{
#ifdef MFALLOWVERBOSE
        if(verbose){printf(" not in Omega\n");fflush(stdout);}
#endif
       }
     }
   }

  sprintf(MFAtlasErrorMsg,"Couldn't get a point on boundary of chart %d",chart);
  MFSetError(e,4,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
  return NULL;
 }

void MFAtlasRemoveChartFromBoundaryListByChart(MFAtlas A,int chart, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAtlasRemoveChartFromBoundaryListByChart"};
  int i;

  printf("%s, chart %d\n",RoutineName,chart);fflush(stdout);
  for(i=0;i<A->nBnd;i++)
   {
    printf("      bnd[%d]=%d\n",i,A->bnd[i]);fflush(stdout);
    if(A->bnd[i]==chart)
     {
      MFAtlasRemoveChartFromBoundaryList(A,i,e);
      return;
     }
   }
  printf("  chart %d not found on boundary list!\n",chart);fflush(stdout);
  return;
 }

double MFAtlasChartSuggestedRadius(MFAtlas A,int chart, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAtlasChartSuggestedRadius"};

#ifdef MFNOCONFIDENCE
  if(A==NULL)
   {
    sprintf(MFAtlasErrorMsg,"Pointer to Atlas (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    return -1.;
   }

  if(chart<0||chart>=A->nCharts)
   {
    sprintf(MFAtlasErrorMsg,"Chart %d (argument 2) is invalid. Must be positive and less than %d (total number of charts).",chart,A->nCharts);
    MFSetError(e,4,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    return -1.;
   }
#endif

  return MFChartSuggestedRadius(A->chart[chart],e);
 }

void MFAtlasAddClipF(MFAtlas thisAtlas,double (*clipf)(MFNVector,MFErrorHandler), MFErrorHandler e)
 {
  static char RoutineName[]={"MFAtlasAddClipF"};
  int verbose=0;

  if(thisAtlas->nclipf>=thisAtlas->mclipf)
   {
    thisAtlas->mclipf+=10;
    thisAtlas->clipf=(double (**)(MFNVector,MFErrorHandler))realloc(thisAtlas->clipf,(thisAtlas->mclipf)*sizeof(double (*)(MFNVector,MFErrorHandler)));

#ifndef MFNOSAFETYNET
    if(thisAtlas->clipf==NULL)
     {
      sprintf(MFAtlasErrorMsg,"Out of memory trying to allocate %d bytes",(thisAtlas->mclipf)*sizeof(double (*)(MFNVector,MFErrorHandler)));
      MFSetError(e,12,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

    thisAtlas->clipIndx=(int*)realloc(thisAtlas->clipIndx,(thisAtlas->mclipf)*sizeof(int));

#ifndef MFNOSAFETYNET
    if(thisAtlas->clipIndx==NULL)
     {
      sprintf(MFAtlasErrorMsg,"Out of memory trying to allocate %d bytes",(thisAtlas->mclipf)*sizeof(int));
      MFSetError(e,12,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

   }
  thisAtlas->clipf[thisAtlas->nclipf]=clipf;
  thisAtlas->clipIndx[thisAtlas->nclipf]=thisAtlas->nclipIndx;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s, Just added Clipping function %d index %d\n",RoutineName,thisAtlas->nclipf,thisAtlas->nclipIndx);fflush(stdout);}
#endif

  thisAtlas->nclipf++;
  thisAtlas->nclipIndx--;
 }

MFBinaryTree MFAtlasGetBB(MFAtlas A, MFErrorHandler e)
 {
  return A->BTree;
 }

void MFAtlasSetExpFactor(MFAtlas A, double s, MFErrorHandler e)
 {
  A->expFactor=s;
  return;
 }

double MFAtlasGetExpFactor(MFAtlas A, MFErrorHandler e)
 {
  return A->expFactor;
 }

void MFAtlasClearClipF(MFAtlas thisAtlas, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAtlasClearClipF"};
  int i;
  int verbose=0;

  for(i=0;i<thisAtlas->nclipf;i++)thisAtlas->clipf[i]=NULL;
  thisAtlas->nclipf=0;
  return;
 }

void MFAtlasClosePlotfile(MFAtlas thisAtlas, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAtlasClosePlotfile"};

  printf("Close Plotfile\n");fflush(stdout);

  if(thisAtlas->plotFileFid!=NULL)
   {
    fclose(thisAtlas->plotFileFid);
    thisAtlas->plotFileFid=NULL;
   }
  return;
 }

void MFAtlasCloseCenterfile(MFAtlas thisAtlas, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAtlasCloseCenterfile"};

  if(thisAtlas->centerFileFid!=NULL)
   {
    fclose(thisAtlas->centerFileFid);
    thisAtlas->centerFileFid=NULL;
   }
  return;
 }

void MFAtlasSetNearRtn(MFAtlas A,int (*isNear)(MFAtlas,MFChart,MFChart, MFErrorHandler e), MFErrorHandler e)
 {
  A->isNear=isNear;
 }

void MFUpdateNeighborsReferenceMarks(MFAtlas A, int ichart, int maxdepth, MFErrorHandler e)
 {
  static char RoutineName[]={"MFUpdateNeighborsReferenceMarks"};
  int i,n,hs,j,l;
  int minmark,mark;
  int neighbor;
  MFPolytope P;
  MFChart chart;
  MFChart nchart;
  static int verbose=0;
  int nNeighs;
  int mNeighs;
  int *neigh;
  int *bneigh;
  int N;
  int found;
  int depth;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s: chart=%d, maxdepth=%d\n",RoutineName,ichart,maxdepth);fflush(stdout);}
#endif

  nNeighs=1;
  mNeighs=100;

  neigh=(int*)malloc(mNeighs*sizeof(int));

#ifndef MFNOSAFETYNET
    if(neigh==NULL)
     {
      sprintf(MFAtlasErrorMsg,"Out of memory trying to allocate %d bytes",mNeighs*sizeof(int));
      MFSetError(e,12,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

  bneigh=(int*)malloc(mNeighs*sizeof(int));

#ifndef MFNOSAFETYNET
    if(bneigh==NULL)
     {
      sprintf(MFAtlasErrorMsg,"Out of memory trying to allocate %d bytes",mNeighs*sizeof(int));
      MFSetError(e,12,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

  neigh[0]=ichart;
  bneigh[0]=1;

  for(depth=0;depth<maxdepth;depth++)
   {
    N=nNeighs;
    for(i=0;i<N;i++)
     {
      if(bneigh[i])
       {
        chart=MFAtlasChart(A,neigh[i],e);
        P=MFChartPolytope(chart,e);
        if(P!=NULL)
         {
          n=MFPolytopeNumberOfFaces(P,e);
          bneigh[i]=0;
          for(j=0;j<n;j++)
           {
            hs=MFPolytopeFaceIndex(P,j,e);
            neighbor=-1;
            if(hs>=MFAtlasOffset(A,e) && hs< MFAtlasNumberOfHalfSpaces(A,e))neighbor=MFAtlasRightPolytope(A,hs,e);
  
            if(neighbor>-1 && neighbor<MFAtlasNumberOfCharts(A,e))
             {
              found=0;
              for(l=0;l<nNeighs;l++)if(neigh[l]==neighbor)found=1;
              if(!found)
               {
                if(nNeighs+1>mNeighs)
                 {
                  mNeighs+=100;
                  neigh=(int*)realloc((void*)neigh,mNeighs*sizeof(int));

#ifndef MFNOSAFETYNET
                  if(neigh==NULL)
                   {
                    sprintf(MFAtlasErrorMsg,"Out of memory trying to allocate %d bytes",mNeighs*sizeof(int));
                    MFSetError(e,12,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
                    MFErrorHandlerOutOfMemory(e);
                    return;
                   }
#endif

                  bneigh=(int*)realloc((void*)bneigh,mNeighs*sizeof(int));

#ifndef MFNOSAFETYNET
                  if(bneigh==NULL)
                   {
                    sprintf(MFAtlasErrorMsg,"Out of memory trying to allocate %d bytes",mNeighs*sizeof(int));
                    MFSetError(e,12,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
                    MFErrorHandlerOutOfMemory(e);
                    return;
                   }
#endif

                  if(neigh==NULL)abort();
                  if(bneigh==NULL)abort();
                 }
                neigh[nNeighs]=neighbor;
                bneigh[nNeighs]=1;
                nNeighs++;
               }
             }
           }
         }
       }
     }
   }

#ifdef MFALLOWVERBOSE
  if(verbose){printf("  %d neighbors\n",nNeighs);fflush(stdout);}
#endif

  for(i=1;i<nNeighs;i++)
   {

#ifdef MFALLOWVERBOSE
    if(verbose){printf("  %d: chart=%d\n",i,neigh[i]);fflush(stdout);}
#endif

    chart=MFAtlasChart(A,neigh[i],e);
    P=MFChartPolytope(chart,e);
    if(P!=NULL)
     {
      n=MFPolytopeNumberOfFaces(P,e);
      minmark=-1;
      for(j=0;j<n;j++)
       {
        hs=MFPolytopeFaceIndex(P,j,e);
        
        neighbor=-1;
        if(hs>MFAtlasOffset(A,e) && hs< MFAtlasNumberOfHalfSpaces(A,e))neighbor=MFAtlasRightPolytope(A,hs,e);
    
        if(neighbor>-1 && neighbor<MFAtlasNumberOfCharts(A,e))
         {
          nchart=MFAtlasChart(A,neighbor,e);
          mark=MFChartReferenceNumber(nchart,e);
          if((MFChartIsSingular(chart,e)&&MFChartIsSingular(nchart,e))||(!MFChartIsSingular(chart,e)&&!MFChartIsSingular(nchart,e)))
           {
            if(minmark==-1||mark<minmark)minmark=mark;
           }
         }
       }
      if(minmark!=-1)
       {
        mark=MFChartReferenceNumber(chart,e);

#ifdef MFALLOWVERBOSE
        if(verbose){printf("    old mark:%d new mark:%d\n",mark,minmark+1);fflush(stdout);}
#endif

        if(minmark+1<mark)MFChartSetReferenceNumber(chart,minmark+1,e);
       }
     }
   }

  free(neigh);
  free(bneigh);
  return;
 }

#ifdef __cplusplus
}
#endif
