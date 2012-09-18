/* 
 *  PROGRAM NAME:  multifario
 *
 *  (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
 *  CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
 *
 *  Please refer to the LICENSE file in the top directory
 *
 *  author: Mike Henderson mhender@watson.ibm.com
 */

static char *id="@(#) $Id: MFChart.c,v 1.8 2011/07/21 17:42:46 mhender Exp $";

#include <MFAtlas.h>
#include <MFChart.h>
#include <MFKVector.h>
#include <MFNVector.h>
#include <MFNKMatrix.h>
#include <MFImplicitMF.h>
#include <MFPrint.h>
#include <MFErrorHandler.h>
#include <stdlib.h>
#include <string.h>

#define MFREALLYTHROWITAWAY

static char MFChartErrorMsg[256]="";

#ifdef __cplusplus
 extern "C" {
#endif

int MFPolytopeTotallyInterior(MFPolytope,MFKVector,double,MFErrorHandler);
int *MFPolytopeVertexIndexSet(MFPolytope,int,MFErrorHandler);
int MFPolytopeVertexNumberOfIndices(MFPolytope,int,MFErrorHandler);
int MFChartPageIn(MFChart,FILE*,MFErrorHandler);

struct MFChartSt {
                  int k;
                  int n;
                  MFImplicitMF M;

                  MFPolytope P;
                  int        allAlphasPositive;
                  int        singular;

                  double     R;
                  double     SuggestR;
                  MFNVector  u;
                  MFNKMatrix Phi;
                  int        posInBList;
                  int        posInAtlas;

                  FILE       *fid;
                  int        paged;
                  int        nearBoundary;
                  int        changed;
                  long       indexInPageFile;
                  struct MFChartSt *prev;
                  struct MFChartSt *next;

                  int        refNumber;
                  int        nRefs;
                 };

#define DENSE 0
#define AUTO 1
#define LOCA 2

MFChart MFCreateChartWithCubeSize(MFImplicitMF M,MFNVector u,MFNKMatrix TS, double R, double L, MFErrorHandler e)
 {
  static char RoutineName[]={"MFCreateChartWithCubeSize"};
  MFChart thisChart;

#ifdef MFNOCONFIDENCE
  if(M==NULL)
   {
    sprintf(MFChartErrorMsg,"M (argument 1) is NULL");
    MFSetError(e,8,RoutineName,MFChartErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(MFIMF_K(M,e)<1)
   {
    sprintf(MFChartErrorMsg,"Manifold has illegal dimension %d. Must be positive",MFIMF_K(M,e));
    MFSetError(e,4,RoutineName,MFChartErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(strcmp(MFNVGetId(u,e),"LOCA") && MFIMF_N(M,e)<MFIMF_K(M,e))
   {
    sprintf(MFChartErrorMsg,"Manifold has illegal embedding dimension (%d). Must be greater or equal to it\'s dimension %d.",MFIMF_N(M,e),MFIMF_K(M,e));
    MFSetError(e,12,RoutineName,MFChartErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(u==NULL)
   {
    sprintf(MFChartErrorMsg,"u (argument 2) is NULL");
    MFSetError(e,8,RoutineName,MFChartErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  thisChart=(MFChart)(MFChart)malloc(sizeof(struct MFChartSt));

#ifndef MFSAFETYNET
  if(thisChart==NULL)
   {
    sprintf(MFChartErrorMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct MFChartSt));
    MFSetError(e,12,RoutineName,MFChartErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  thisChart->M=M;
  MFRefImplicitMF(M,e);
  thisChart->u=u;
  MFRefNVector(u,e);
  thisChart->n=MFIMF_N(M,e);
  thisChart->k=MFIMF_K(M,e);

  thisChart->Phi=TS;
  MFRefNKMatrix(thisChart->Phi,e);
  thisChart->R=R;
  thisChart->SuggestR=R;
  thisChart->refNumber=0;

  thisChart->P=MFCreateHyperCubeAtOrigin(thisChart->k,L,e);
/*
 *thisChart->P=MFCreateSimplexAtOrigin(thisChart->k,L,e);
 */

#ifdef MFNOCONFIDENCE
  if(TS==NULL)
   {
    sprintf(MFChartErrorMsg,"TangentSpace (argument 3) is NULL ");
    MFSetError(e,12,RoutineName,MFChartErrorMsg,__LINE__,__FILE__);
    return NULL;
   }
#endif

  thisChart->allAlphasPositive=1;
  thisChart->singular=0;
  thisChart->fid=NULL;
  thisChart->paged=0;
  thisChart->nearBoundary=0;
  thisChart->changed=1;
  thisChart->indexInPageFile=-1;
  thisChart->posInBList=-1;
  thisChart->posInAtlas=-1;
  thisChart->next=NULL;
  thisChart->prev=NULL;
  thisChart->paged=0;
  thisChart->nRefs=1;

  if(thisChart==(MFChart)0x0087dd58){printf("%s 0x%8.8x, center 0x%8.8x, nRefs now %d\n",RoutineName,thisChart,thisChart->u,thisChart->nRefs);fflush(stdout);}

  return thisChart;
 }

MFChart MFCreateChart(MFImplicitMF M,MFNVector u,MFNKMatrix TS, double R, MFErrorHandler e)
 {
  static char RoutineName[]={"MFCreateChart"};
  return MFCreateChartWithCubeSize(M,u,TS,R,1.05*R,e);
 }

void MFSubtractHalfSpaceFromChart(MFChart chart,int i,MFKVector n,double o, MFErrorHandler e)
 {
  static char RoutineName[]={"MFSubtractHalfSpaceFromChart"};

#ifdef MFNOCONFIDENCE
  if(chart==NULL)
   {
    sprintf(MFChartErrorMsg,"Pointer to Chart (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFChartErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  MFSubtractHalfSpaceFromPolytope(chart->P,i,n,o,e);
  if(o<0.)chart->allAlphasPositive=0;
  chart->changed=1;
  return;
 }

void MFRefChart(MFChart chart, MFErrorHandler e)
 {
  static char RoutineName[]={"MFRefChart"};

#ifdef MFNOCONFIDENCE
  if(chart==NULL)
   {
    sprintf(MFChartErrorMsg,"Pointer to Chart (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFChartErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  chart->nRefs++;
  if(chart==(MFChart)0x00898f28){printf("%s 0x%8.8x, center 0x%8.8x, nRefs now %d\n",RoutineName,chart,chart->u,chart->nRefs);fflush(stdout);}
  return;
 }

void MFFreeChart(MFChart chart, MFErrorHandler e)
 {
  static char RoutineName[]={"MFFreeChart"};

#ifdef MFNOCONFIDENCE
  if(chart==NULL)
   {
    sprintf(MFChartErrorMsg,"Pointer to Chart (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFChartErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif
  if(chart==NULL)return;

  chart->nRefs--;

  if(chart->nRefs>0)return;
  if(chart==(MFChart)0x0087dd58){printf("%s 0x%8.8x, center 0x%8.8x, nRefs=%d\n",RoutineName,chart,chart->u,chart->nRefs);fflush(stdout);}

  if(chart==(MFChart)0x0087dd58){printf("  free chart->M 0x%8.8x\n",chart->M);fflush(stdout);}
  if(chart->M!=NULL){MFFreeImplicitMF(chart->M,e);chart->M=NULL;};
  if(chart==(MFChart)0x0087dd58){printf("  free chart->u 0x%8.8x\n",chart->u);fflush(stdout);}
  if(chart->u!=NULL){MFFreeNVector(chart->u,e);chart->u=NULL;};
  if(chart==(MFChart)0x0087dd58){printf("  free chart->P 0x%8.8x\n",chart->P);fflush(stdout);}
  if(chart->P!=NULL)
   {
    MFFreePolytope(chart->P,e);
    chart->P=NULL;
   }
  if(chart==(MFChart)0x0087dd58){printf("  free chart->Phi 0x%8.8x\n",chart->Phi);fflush(stdout);}
  if(chart->Phi!=NULL){MFFreeNKMatrix(chart->Phi,e);chart->Phi=NULL;};
  if(chart==(MFChart)0x0087dd58){printf("  free chart 0x%8.8x\n",chart);fflush(stdout);}
  free(chart);
  if(chart==(MFChart)0x0087dd58){printf("  done");fflush(stdout);}

  return;
 }

MFPolytope MFChartPolytope(MFChart chart, MFErrorHandler e)
 {
  static char RoutineName[]={"MFChartPolytope"};

#ifdef MFNOCONFIDENCE
  if(chart==NULL)
   {
    sprintf(MFChartErrorMsg,"Pointer to Chart (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFChartErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  return chart->P;
 }

void MFSetChartPolytope(MFChart chart, MFPolytope P, MFErrorHandler e)
 {
  static char RoutineName[]={"MFSetChartPolytope"};

#ifdef MFNOCONFIDENCE
  if(chart==NULL)
   {
    sprintf(MFChartErrorMsg,"Pointer to Chart (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFChartErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  chart->P=P;

  return;
 }

void MFSetChartRadius(MFChart chart, double R, MFErrorHandler e)
 {
  static char RoutineName[]={"MFSetChartRadius"};
  static int verbose=0;

#ifdef MFNOCONFIDENCE
  if(chart==NULL)
   {
    sprintf(MFChartErrorMsg,"Pointer to Chart (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFChartErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s chart 0x%8.8x, R was %lf, now %lf\n",RoutineName,chart,chart->R,R);fflush(stdout);}
#endif
  chart->R=R;

  return;
 }

MFNVector MFChartCenter(MFChart chart, MFErrorHandler e)
 {
  static char RoutineName[]={"MFChartCenter"};

#ifdef MFNOCONFIDENCE
  if(chart==NULL)
   {
    sprintf(MFChartErrorMsg,"Pointer to Chart (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFChartErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  if(chart->paged)MFChartPageIn(chart,chart->fid,e);

  return chart->u;
 }

MFNKMatrix MFChartTangentSpace(MFChart chart, MFErrorHandler e)
 {
  static char RoutineName[]={"MFChartTangentSpace"};

#ifdef MFNOCONFIDENCE
  if(chart==NULL)
   {
    sprintf(MFChartErrorMsg,"Pointer to Chart (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFChartErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  if(chart->paged)MFChartPageIn(chart,chart->fid,e);

  return chart->Phi;
 }

double MFChartRadius(MFChart chart, MFErrorHandler e)
 {
  static char RoutineName[]={"MFChartRadius"};

#ifdef MFNOCONFIDENCE
  if(chart==NULL)
   {
    sprintf(MFChartErrorMsg,"Pointer to Chart (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFChartErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  return chart->R;
 }

int MFChartEvaluate(MFChart chart,MFKVector s,MFNVector u, MFErrorHandler e)
 {
  static char RoutineName[]={"MFChartEvaluate"};
  MFNVector u0;
  int rc;
  int verbose=0;

#ifdef MFNOCONFIDENCE
  if(chart==NULL)
   {
    sprintf(MFChartErrorMsg,"Pointer to Chart (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFChartErrorMsg,__LINE__,__FILE__);
    return 0;
   }
#endif

  if(chart->paged)MFChartPageIn(chart,chart->fid,e);

/*u0=MFCloneNVector(u);
#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    fprintf(stderr,"0x%8.8x=u0\n",u0);
    fprintf(stderr,"0x%8.8x=u0(0x%8.8x).data\n",MFNV_CStar(u0,e),u0);
    fflush(stderr);
   }
#endif
  if(chart->paged)MFChartPageIn(chart,chart->fid,e);
  MFChartPointInTangentSpace(chart,s,u0,e);
  rc=MFIMFProject(chart->M,u0,chart->Phi,u,e);
  MFFreeNVector(u0);*/

  rc=MFIMFProjectFromCenter(chart->M,chart->u,chart->Phi,s,u,e);
  return rc;
 }

int MFChartInterior(MFChart chart,MFKVector s, MFErrorHandler e)
 {
  static char RoutineName[]={"MFChartInterior"};
  int result;

#ifdef MFNOCONFIDENCE
  if(chart==NULL)
   {
    sprintf(MFChartErrorMsg,"Pointer to Chart (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFChartErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  if(MFKVNorm(s,e)>chart->R)
   { 
    return 0;
   }else{
    result=MFPolytopeInterior(chart->P,s,e);
    return result;
   }
 }

void MFChartProjectIntoTangentSpace(MFChart chart,MFNVector u,MFKVector s, MFErrorHandler e)
 {
  static char RoutineName[]={"MFChartProjectIntoTangentSpace"};
  MFNVector v;

#ifdef MFNOCONFIDENCE
  if(chart==NULL)
   {
    sprintf(MFChartErrorMsg,"Chart (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFChartErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(u==NULL)
   {
    sprintf(MFChartErrorMsg,"u (argument 2) is NULL");
    MFSetError(e,8,RoutineName,MFChartErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(s==NULL)
   {
    sprintf(MFChartErrorMsg,"s (argument 3) is NULL");
    MFSetError(e,8,RoutineName,MFChartErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(chart->n==0)
   {
    sprintf(MFChartErrorMsg,"Chart has zero dimension");
    MFSetError(e,8,RoutineName,MFChartErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  if(chart->paged)MFChartPageIn(chart,chart->fid,e);

#ifdef MFNOCONFIDENCE
  if(chart->u==NULL)
   {
    sprintf(MFChartErrorMsg,"Chart center is NULL");
    MFSetError(e,8,RoutineName,MFChartErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  v=MFCloneNVector(u,e);

  MFNSpaceDirection(MFIMFNSpace(chart->M,e),chart->u,u,v,e);
  MFMVMulT(MFIMFNSpace(chart->M,e),chart->Phi,v,s,e);
  MFFreeNVector(v,e);

  return;
 }

void MFChartProjectVectorIntoTangentSpace(MFChart chart,MFNVector u,MFKVector s, MFErrorHandler e)
 {
  static char RoutineName[]={"MFChartProjectVectorIntoTangentSpace"};

#ifdef MFNOCONFIDENCE
  if(chart==NULL)
   {
    sprintf(MFChartErrorMsg,"Chart (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFChartErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(u==NULL)
   {
    sprintf(MFChartErrorMsg,"u (argument 2) is NULL");
    MFSetError(e,8,RoutineName,MFChartErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(s==NULL)
   {
    sprintf(MFChartErrorMsg,"s (argument 3) is NULL");
    MFSetError(e,8,RoutineName,MFChartErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(chart->n==0)
   {
    sprintf(MFChartErrorMsg,"Chart has zero dimension");
    MFSetError(e,8,RoutineName,MFChartErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  if(chart->paged)MFChartPageIn(chart,chart->fid,e);

#ifdef MFNOCONFIDENCE
  if(chart->u==NULL)
   {
    sprintf(MFChartErrorMsg,"Chart center is NULL");
    MFSetError(e,8,RoutineName,MFChartErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  MFMVMulT(MFIMFNSpace(chart->M,e),chart->Phi,u,s,e);

  return;
 }

void MFChartPointInTangentSpace(MFChart chart,MFKVector s,MFNVector u, MFErrorHandler e)
 {
  static char RoutineName[]={"MFChartPointInTangentSpace"};
  MFNVector v;
  int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s\n",RoutineName);fflush(stdout);}
#endif

#ifdef MFNOCONFIDENCE
  if(chart==NULL)
   {
    sprintf(MFChartErrorMsg,"Pointer to Chart (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFChartErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  if(chart->paged)MFChartPageIn(chart,chart->fid,e);

  v=MFCloneNVector(u,e);

  MFMVMul(MFIMFNSpace(chart->M,e),chart->Phi,s,v,e);
  MFNSpaceAdd(MFIMFNSpace(chart->M,e),v,chart->u,u,e);
#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    MFPrintNVector(stdout,chart->u,e);
    printf("+");
    MFPrintNKMatrix(stdout,chart->Phi,e);
    MFPrintKVector(stdout,s,e);
    printf("=");
    MFPrintNVector(stdout,u,e);
    printf("\n");
    fflush(stdout);
   }
#endif

  MFFreeNVector(v,e);

  return;
 }

int MFChartHasBoundary(MFChart chart, MFErrorHandler e)
 {
  static char RoutineName[]={"MFChartHasBoundary"};
  int vertex,nVertices;

#ifdef MFNOCONFIDENCE
  if(chart==NULL)
   {
    sprintf(MFChartErrorMsg,"Pointer to Chart (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFChartErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  nVertices=MFPolytopeNumberOfVertices(chart->P,e);
  for(vertex=0;vertex<nVertices;vertex++)
   {
    if(MFPolytopeRadiusOfVertex(chart->P,vertex,e)>chart->R)
     {
      return 1;
     }
   }
  return 0;
 }

int MFChartK(MFChart chart, MFErrorHandler e)
 {
  static char RoutineName[]={"MFChartK"};

#ifdef MFNOCONFIDENCE
  if(chart==NULL)
   {
    sprintf(MFChartErrorMsg,"Pointer to Chart (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFChartErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  return chart->k;
 }

int MFChartN(MFChart chart, MFErrorHandler e)
 {
  static char RoutineName[]={"MFChartN"};

#ifdef MFNOCONFIDENCE
  if(chart==NULL)
   {
    sprintf(MFChartErrorMsg,"Pointer to Chart (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFChartErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  return chart->n;
 }

void MFWriteChart(FILE *fid, MFChart chart, MFErrorHandler e)
 {
  static char RoutineName[]={"MFWriteChart"};
  int i;

  if(chart->paged)
   {
/*  MFChartPageIn(chart,chart->fid,e);*/
    fprintf(fid,"%s\n","Chart");
    fprintf(fid,"%d %d %d %d %d %lf %d %d\n",chart->n,chart->k,chart->allAlphasPositive,chart->paged,chart->nearBoundary,chart->R,chart->singular,chart->nRefs,chart->refNumber);
    MFWritePolytope(fid,chart->P,e);
   }else{
    fprintf(fid,"%s\n","Chart",e);
    fprintf(fid,"%d %d %d %d %d %lf %d %d\n",chart->n,chart->k,chart->allAlphasPositive,chart->paged,chart->nearBoundary,chart->R,chart->singular,chart->nRefs,chart->refNumber,e);
    MFWritePolytope(fid,chart->P,e);
    MFWriteNVector(fid,chart->u,e);
    MFWriteNKMatrix(fid,chart->Phi,e);
   }

  return;
 }

MFChart MFReadChart(FILE *fid,MFAtlas A, MFErrorHandler e)
 {
  static char RoutineName[]={"MFReadChart"};
  MFChart chart;
  char tag[100]="";

  fscanf(fid,"%s\n",tag);
#ifdef MFNOCONFIDENCE
  if(strcmp(tag,"Chart"))
   {
    sprintf(MFChartErrorMsg,"Next Object is not a Chart! (%s)\n",RoutineName,tag);
    MFSetError(e,12,RoutineName,MFChartErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  chart=(MFChart)malloc(sizeof(struct MFChartSt));

#ifdef MFNOCONFIDENCE
  if(chart==NULL)
   {
    sprintf(MFChartErrorMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct MFChartSt));
    MFSetError(e,12,RoutineName,MFChartErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif
  fscanf(fid,"%d %d %d %d %d %lf %d %d\n",&(chart->n),&(chart->k),&(chart->allAlphasPositive),&(chart->paged),&(chart->nearBoundary),&(chart->R),&(chart->singular),&(chart->nRefs),&(chart->refNumber));

  chart->M=MFAtlasMF(A,e);
  chart->P=MFReadPolytope(fid,e);
  if(chart->paged)
   {
    chart->u=NULL;
    chart->Phi=NULL;
   }else{
    chart->u=MFReadNVector(fid,e);
    chart->Phi=MFReadNKMatrix(fid,e);
   }
  chart->posInBList=-1;
  chart->posInAtlas=-1;

  return chart;
 }

int MFChartTotallyInterior(MFChart chart,MFKVector s,double eps, MFErrorHandler e)
 {
  static char RoutineName[]={"MFChartTotallyInterior"};
  int result;

#ifdef MFNOCONFIDENCE
  if(chart==NULL)
   {
    sprintf(MFChartErrorMsg,"Pointer to Chart (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFChartErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  if(MFKVNorm(s,e)>chart->R)
   { 
    return 0;
   }else{
    result=MFPolytopeTotallyInterior(chart->P,s,eps,e);
    return result;
   }
 }

int MFChartGetPositionInBoundaryList(MFChart thisChart, MFErrorHandler e)
 {
  return thisChart->posInBList;
 }

void MFChartSetPositionInBoundaryList(MFChart thisChart, int pos, MFErrorHandler e)
 {
  thisChart->posInBList=pos;

  return;
 }

int MFChartGetPositionInAtlas(MFChart thisChart, MFErrorHandler e)
 {
  return thisChart->posInAtlas;
 }

void MFChartSetPositionInAtlas(MFChart thisChart, int pos, MFErrorHandler e)
 {
  thisChart->posInAtlas=pos;

  return;
 }

MFImplicitMF MFChartGetManifold(MFChart thisChart, MFErrorHandler e)
 {
  return thisChart->M;
 }

int MFChartIsSingular(MFChart thisChart, MFErrorHandler e)
 {
  return thisChart->singular>0;
 }

void MFChartSetSingular(MFChart thisChart, MFErrorHandler e)
 {
  thisChart->singular=1;
  return;
 }

void MFChartSetNonSingular(MFChart thisChart, MFErrorHandler e)
 {
  thisChart->singular=0;
  return;
 }

int MFChartReferenceNumber(MFChart thisChart, MFErrorHandler e)
 {
  return thisChart->refNumber;
 }

void MFChartSetReferenceNumber(MFChart thisChart, int n, MFErrorHandler e)
 {
  thisChart->refNumber=n;
  return;
 }

int MFChartPageOut(MFChart thisChart, FILE *fid, int index, MFErrorHandler e)
 {
  static char RoutineName[]={"MFChartPageOut"};
  static int i,j;
  double *c;
  int verbose=0;
  int rc;

#ifdef MFNOCONFIDENCE
  if(thisChart==NULL)
   {
    sprintf(MFChartErrorMsg,"Chart is NULL");
    MFSetError(e,0,RoutineName,MFChartErrorMsg,__LINE__,__FILE__);
    return 0;
   }
#endif

/* Returns 1 if written, 0 if already paged out */

  if(thisChart->paged)
   {
#ifdef MFALLOWVERBOSE
    if(verbose){printf("Already paged out\n");fflush(stdout);}
#endif
    return 0;
   }

#ifdef MFREALLYPAGE
#ifdef MFALLOWVERBOSE
  if(verbose){fprintf(stderr,"PAGING out chart, fseek to end of file 0x%8.8x\n",fid);fflush(stderr);}
#endif

  rc=fseek(fid,(long)0,SEEK_END);

#ifdef MFALLOWVERBOSE
  if(verbose){fprintf(stderr,"  return code from fseek is %d\n",rc);fprintf(stderr,"ftell= %d\n",ftell(fid));fflush(stderr);}
  if(verbose){fprintf(stderr,"  write to position %d\n",index);fprintf(stderr,"ftell= %d\n",ftell(fid));fflush(stderr);}
#endif

  c=MFNV_CStar(thisChart->u);

#ifdef MFALLOWVERBOSE
  if(verbose){fprintf(stderr,"PAGING   u starts (%lf,%lf,%lf...\n",c[0],c[1],c[2]);fflush(stderr);}
#endif
  thisChart->fid=fid;
  thisChart->paged=1;
  thisChart->indexInPageFile=ftell(fid);
  for(i=0;i<thisChart->n;i++)fprintf(fid,"%lf ",c[i]);
  fprintf(fid,"\n");

  c=MFNKM_CStar(thisChart->Phi,e);  /* BAD What if it's not dense? */

#ifdef MFALLOWVERBOSE
  if(verbose){fprintf(stderr,"PAGING   Phi starts (%lf,%lf,%lf...\n",c[0],c[1],c[2]);fflush(stderr);}
#endif

  for(j=0;j<thisChart->k;j++)
  for(i=0;i<thisChart->n;i++)fprintf(fid,"%lf ",c[i+thisChart->n*j]);
  fprintf(fid,"\n");
#endif

#ifdef MFREALLYTHROWITAWAY
#ifdef MFALLOWVERBOSE
  if(verbose){printf("Discarding contents of chart\n");fflush(stdout);}
#endif

  MFFreeNVector(thisChart->u,e);thisChart->u=NULL;
  MFFreeNKMatrix(thisChart->Phi,e);thisChart->Phi=NULL;
  MFFreePolytope(thisChart->P,e);thisChart->P=NULL;
#endif

  thisChart->paged=1;

  return 1;
 }

int MFChartPageIn(MFChart thisChart, FILE *fid, MFErrorHandler e)
 {
  static int i,j;
  double *c;
  long pos;
  int verbose=0;
  static char RoutineName[]={"MFChartPageIn"};

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("%s!!!!!\n",RoutineName);
    fflush(stdout);
   }
#endif

  sprintf(MFChartErrorMsg,"Routine %s is broken, FATAL ERROR",RoutineName);
  MFSetError(e,16,RoutineName,MFChartErrorMsg,__LINE__,__FILE__);
  return -1;

/* Returns 1 if paged in, 0 if not paged out */

  if(!(thisChart->paged))return 0;

#ifdef MFALLOWVERBOSE
  if(verbose){fprintf(stderr,"PAGING IN CHART! %d in file\n",thisChart->indexInPageFile);fflush(stderr);}
#endif

  if(thisChart->u!=NULL && thisChart->Phi!=NULL)
   {
#ifdef MFALLOWVERBOSE
    if(verbose){fprintf(stderr,"PAGING Already in\n");fflush(stderr);}
#endif
    return 1;
   }

  if(thisChart->u!=NULL)MFFreeNVector(thisChart->u,e);
  thisChart->u=MFCreateNVector(thisChart->n,e);

  if(thisChart->Phi!=NULL)MFFreeNKMatrix(thisChart->Phi,e);
  thisChart->Phi=MFCreateNKMatrixWithData(thisChart->n,thisChart->k,NULL,e);

#ifdef MFALLOWVERBOSE
  if(verbose){fprintf(stderr,"PAGING IN CHART! read from file\n");}
#endif

  fseek(fid,thisChart->indexInPageFile,SEEK_SET);

#ifdef MFALLOWVERBOSE
  if(verbose){fprintf(stderr,"ftell= %d\n",ftell(fid));}
#endif

  c=MFNV_CStar(thisChart->u,e);
  for(i=0;i<thisChart->n;i++)fscanf(fid,"%lf ",&(c[i]));
  fscanf(fid,"\n");

#ifdef MFALLOWVERBOSE
  if(verbose){fprintf(stderr,"PAGING   u starts (%lf,%lf,%lf...\n",c[0],c[1],c[2]);fflush(stderr);}
#endif

  c=MFNKM_CStar(thisChart->Phi,e);
  for(j=0;j<thisChart->k;j++)
  for(i=0;i<thisChart->n;i++)fscanf(fid,"%lf ",&(c[i+thisChart->n*j]));
  fscanf(fid,"\n");

#ifdef MFALLOWVERBOSE
  if(verbose){fprintf(stderr,"PAGING   Phi starts (%lf,%lf,%lf...\n",c[0],c[1],c[2]);fflush(stderr);}
#endif

  fseek(fid,0,SEEK_END);

  return 1;
 }

void MFChartClean(MFChart thisChart, MFErrorHandler e)
 {
  if(!(thisChart->paged))return;

  if(thisChart->u!=NULL){MFFreeNVector(thisChart->u,e);thisChart->u=NULL;}
  if(thisChart->Phi!=NULL){MFFreeNKMatrix(thisChart->Phi,e);thisChart->Phi=NULL;}

  return;
 }

int MFChartPaged(MFChart thisChart, MFErrorHandler e)
 {
  static char RoutineName[]={"MFChartPaged"};

#ifdef MFNOCONFIDENCE
  if(thisChart==NULL)
   {
    sprintf(MFChartErrorMsg,"Chart is NULL");
    MFSetError(e,0,RoutineName,MFChartErrorMsg,__LINE__,__FILE__);
    return 0;
   }
#endif

  return thisChart->paged;
 }

int MFChartNearBoundary(MFChart thisChart, MFErrorHandler e)
 {
  return thisChart->nearBoundary;
 }

void MFChartSetNearBoundary(MFChart thisChart,int nearBoundary, MFErrorHandler e)
 {
  thisChart->nearBoundary=nearBoundary;
  return;
 }

int MFChartNumberOfFaces(MFChart thisChart, MFErrorHandler e)
 {
  return MFPolytopeNumberOfFaces(thisChart->P,e);
 }

void MFSetSuggestedChartRadius(MFChart chart, double R, MFErrorHandler e)
 {
  static char RoutineName[]={"MFSetSuggestedChartRadius"};

#ifdef MFNOCONFIDENCE
  if(chart==NULL)
   {
    sprintf(MFChartErrorMsg,"Pointer to Chart (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFChartErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  chart->SuggestR=R;

  return;
 }

double MFChartSuggestedRadius(MFChart chart, MFErrorHandler e)
 {
  static char RoutineName[]={"MFChartSuggestedRadius"};

#ifdef MFNOCONFIDENCE
  if(chart==NULL)
   {
    sprintf(MFChartErrorMsg,"Pointer to Chart (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFChartErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  return chart->SuggestR;
 }

void MFChartResetChangedFlag(MFChart chart, MFErrorHandler e)
 {
  chart->changed=0;
  return;
 }

int MFChartChangedFlag(MFChart chart, MFErrorHandler e)
 {
  return chart->changed;
 }

void MFChartSetPaged(MFChart thisChart, MFErrorHandler e)
 {
  static char RoutineName[]={"MFChartSetPaged"};

#ifdef MFNOCONFIDENCE
  if(thisChart==NULL)
   {
    sprintf(MFChartErrorMsg,"Chart is NULL");
    MFSetError(e,0,RoutineName,MFChartErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  thisChart->paged=1;
  return;
 }

MFChart MFCreateBoundaryChart(MFImplicitMF M,MFNVector u,MFNKMatrix TS, int m, int dir, double R,MFErrorHandler e)
 {
  static char RoutineName[]={"MFCreateBoundaryChart"};
  MFChart thisChart;
  MFKVector n;
  int i,k;

#ifdef MFNOCONFIDENCE
  if(M==NULL)
   {
    sprintf(MFChartErrorMsg,"M (argument 1) is NULL");
    MFSetError(e,8,RoutineName,MFChartErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(MFIMF_K(M,e)<1)
   {
    sprintf(MFChartErrorMsg,"Manifold has illegal dimension %d. Must be positive",MFIMF_K(M,e));
    MFSetError(e,4,RoutineName,MFChartErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(strcmp(MFNVGetId(u,e),"LOCA") && MFIMF_N(M,e)<MFIMF_K(M,e))
   {
    sprintf(MFChartErrorMsg,"Manifold has illegal embedding dimension (%d). Must be greater or equal to it\'s dimension %d.",MFIMF_N(M,e),MFIMF_K(M,e));
    MFSetError(e,12,RoutineName,MFChartErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(u==NULL)
   {
    sprintf(MFChartErrorMsg,"u (argument 2) is NULL");
    MFSetError(e,8,RoutineName,MFChartErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  thisChart=MFCreateChart(M,u,TS,R,e);

  k=MFIMF_K(M,e);
  n=MFCreateKVector(k,e);
  for(i=0;i<k;i++)(MFKV_CStar(n,e))[i]=0.;

  if(m>k-1){printf("m(%d)>k+1(%d)\n",m,k+1);fflush(stdout);return;}
  
  for(i=0;i<m;i++)
   {
    if(i>0)(MFKV_CStar(n,e))[i-1]=0.;
    (MFKV_CStar(n,e))[i]=-dir;
    MFSubtractHalfSpaceFromChart(thisChart,-2*(i+1),n,0.,e);
   }
  MFFreeKVector(n,e);

  return thisChart;
 }

MFChart MFChartGetNextChart(MFChart thisChart, MFErrorHandler e)
 {
  return thisChart->next;
 }

MFChart MFChartGetPrevChart(MFChart thisChart, MFErrorHandler e)
 {
  return thisChart->prev;
 }

void MFChartSetNextChart(MFChart thisChart, MFChart next, MFErrorHandler e)
 {
  thisChart->next=next;
  return;
 }

void MFChartSetPrevChart(MFChart thisChart, MFChart prev, MFErrorHandler e)
 {
  thisChart->prev=prev;
  return;
 }

#ifdef __cplusplus
}
#endif
