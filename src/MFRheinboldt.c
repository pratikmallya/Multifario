/*
 *  PROGRAM NAME:  multifario
 *
 *  (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
 *  CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
 *
 *  Please refer to the LICENSE file in the top directory
 *
 *      author: Mike Henderson mhender@watson.ibm.com
 *      date:   November 21, 2002
 */

static char *id="@(#) $Id: MFRheinboldt.c,v 1.7 2007/11/16 14:27:34 mhender Exp $";

static char MFRheinboldtErrorMsg[256];

#include <MFAtlas.h>
#include <MFContinuationMethod.h>
#include <MFImplicitMF.h>
#include <MFNRegion.h>
#include <MFNVector.h>
#include <MFPrint.h>
#include <MFPolytope.h>
#include <stdio.h>
#include <math.h>
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
int *MFPolytopeVertexIndexSet(MFPolytope,int,MFErrorHandler);
int MFPolytopeIntersectIndexSets(MFPolytope,int,int,int*,MFErrorHandler);
int MFPolytopeVertexNumberOfIndices(MFPolytope,int,MFErrorHandler);

static void MFRheinboldtExtendAtlas(MFContinuationMethod H,MFAtlas A, MFImplicitMF M, MFNRegion Omega, int m, MFNVector *u0, MFNKMatrix *Phi,MFErrorHandler);
static void MFRheinboldtCloseAtlas(MFContinuationMethod,MFAtlas,MFErrorHandler);
static void MFRheinboldtFlushAtlas(MFContinuationMethod,MFAtlas,MFErrorHandler);
static void MFRheinboldtFreeParmBlock(void*,MFErrorHandler);
static void MFRheinboldtWriteChartToPlotFile(FILE*,MFAtlas,MFPolytope,int,int,int*,MFErrorHandler);

struct MFRheinboldtParmBlockST {
                             int verbose;
                             int n;
                             int maxCharts;
                             double delta;
                             char *fileName;
                            };
typedef struct MFRheinboldtParmBlockST *MFRheinboldtParmBlock;

MFContinuationMethod MFCreateRheinboldtsMethod(MFErrorHandler e)
 {
  static char RoutineName[]={"MFCreateRheinboldtsMethod"};
  MFContinuationMethod result;
  MFRheinboldtParmBlock parms;

  parms=(struct MFRheinboldtParmBlockST*)malloc(sizeof(struct MFRheinboldtParmBlockST));

#ifndef MFNOSAFETYNET
  if(parms==NULL)
   {
    sprintf(MFRheinboldtErrorMsg,"Out of space trying to allocate %d bytes.",sizeof(struct MFRheinboldtParmBlockST));
    MFSetError(e,12,RoutineName,MFRheinboldtErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    MFFreeContinuationMethod(result,e);
    return NULL;
   }
#endif

  parms->verbose=0;
  parms->n=10;
  parms->delta=.1;
  parms->fileName=(char*)malloc(6*sizeof(char));

#ifndef MFNOSAFETYNET
  if(parms->fileName==NULL)
   {
    sprintf(MFRheinboldtErrorMsg,"Out of space trying to allocate %d bytes.",6*sizeof(char));
    MFSetError(e,12,RoutineName,MFRheinboldtErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    MFFreeContinuationMethod(result,e);
    return NULL;
   }
#endif

  strcpy(parms->fileName,"Atlas");

  result=MFCreateContinuationMethodBaseClass("Rheinboldt",e);
  MFContinuationMethodSetParmBlock(result,parms,e);
  MFContinuationMethodSetFreeParmBlock(result,MFRheinboldtFreeParmBlock,e);
  MFContinuationMethodSetCloseAtlas(result,MFRheinboldtCloseAtlas,e);
  MFContinuationMethodSetExtendAtlasMultipleWithTangents(result,MFRheinboldtExtendAtlas,e);
  MFContinuationMethodSetFlushAtlas(result,MFRheinboldtFlushAtlas,e);

  return result;
 }

void MFRheinboldtFreeParmBlock(void *parmBlock, MFErrorHandler e)
 { 
  MFRheinboldtParmBlock parms;

  parms=(MFRheinboldtParmBlock)(parmBlock);
  if(parms->fileName!=NULL)free(parms->fileName);
  free(parms);

  return;
 }

char *MFRheinboldtGetFilename(MFContinuationMethod H, MFErrorHandler e)
 {
  static char RoutineName[]={"MFRheinboldtGetFileName"};
  MFRheinboldtParmBlock parms;

#ifdef MFNOCONFIDENCE
  if(H==NULL)
   {
    sprintf(MFRheinboldtErrorMsg,"MFContinuation Method (argument 1) is NULL.");
    MFSetError(e,12,RoutineName,MFRheinboldtErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(strcmp(MFContinuationMethodGetType(H,e),"Rheinboldt"))
   {
    sprintf(MFRheinboldtErrorMsg,"MFContinuation Method (argument 1) is not Rheinboldt's Method.");
    MFSetError(e,12,RoutineName,MFRheinboldtErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  parms=(MFRheinboldtParmBlock)(MFContinuationMethodGetParmBlock(H,e));
  return parms->fileName;
 }

void MFRheinboldtSetFilename(MFContinuationMethod H,char *name, MFErrorHandler e)
 {
  static char RoutineName[]={"MFRheinboldtSetFileName"};
  MFRheinboldtParmBlock parms;

#ifdef MFNOCONFIDENCE
  if(H==NULL)
   {
    sprintf(MFRheinboldtErrorMsg,"MFContinuation Method (argument 1) is NULL.");
    MFSetError(e,12,RoutineName,MFRheinboldtErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(strcmp(MFContinuationMethodGetType(H,e),"Rheinboldt"))
   {
    sprintf(MFRheinboldtErrorMsg,"MFContinuation Method (argument 1) is not Rheinboldt's Method.");
    MFSetError(e,12,RoutineName,MFRheinboldtErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  parms=(MFRheinboldtParmBlock)(MFContinuationMethodGetParmBlock(H,e));
  free(parms->fileName);
  parms->fileName=(char*)malloc((strlen(name)+1)*sizeof(char));

#ifndef MFNOSAFETYNET
  if(parms->fileName==NULL)
   {
    sprintf(MFRheinboldtErrorMsg,"Out of space trying to allocate %d bytes.",(strlen(name)+1)*sizeof(char));
    MFSetError(e,12,RoutineName,MFRheinboldtErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  strcpy(parms->fileName,name);
  return;
 }

void MFRheinboldtExtendAtlas(MFContinuationMethod H, MFAtlas A, MFImplicitMF M, MFNRegion Omega, int m, MFNVector *u0, MFNKMatrix *Phi0, MFErrorHandler e)
 {
  static char RoutineName[]={"MFRheinboldtExtendAtlas"};
  int i,j;
  int verbose;
  int n,k;
  int nx;
  double h;
  int chart;
  MFKVector S;
  double *s;
  MFNVector u;
  MFNVector unew;
  MFNKMatrix Phi;
  char *name;
  int dir;
  int kkk,kmax;
  int rc;
  int *c;

  if(m==0)return;

  h=MFRheinboldtGetRealParameter(H,"delta",e);
  verbose=MFRheinboldtGetIntegerParameter(H,"verbose",e);
  kmax=MFRheinboldtGetIntegerParameter(H,"maxCharts",e);
  nx=MFRheinboldtGetIntegerParameter(H,"n",e);
  name=MFRheinboldtGetFilename(H,e);

  n=MFAtlasN(A,e);
  k=MFAtlasK(A,e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s, kmax=%d, verbose=%d\n",RoutineName,kmax,verbose);fflush(stdout);}
#endif

  u=MFCloneNVector(u0[0],e);
  
  S=MFCreateKVector(k,e);
  s=MFKV_CStar(S,e);

  c=(int*)malloc(k*sizeof(int));

#ifndef MFNOSAFETYNET
  if(c==NULL)
   {
    sprintf(MFRheinboldtErrorMsg,"Out of space trying to allocate %d bytes.",k*sizeof(int));
    MFSetError(e,12,RoutineName,MFRheinboldtErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  for(i=0;i<k;i++)c[i]=0;

  kkk=0;
  MFAtlasSetVerbose(A,verbose,e);
  if(Phi0==NULL||Phi0[i]==NULL)
   {
    Phi=MFIMFTangentSpace(M,u0[0],e);
    rc=MFAtlasAddChartToList(A,MFCreateChart(MFAtlasMF(A,e),u0[0],Phi,h,e),e);
    MFFreeNKMatrix(Phi,e);
   }else{
    rc=MFAtlasAddChartToList(A,MFCreateChart(MFAtlasMF(A,e),u0[0],Phi0[0],h,e),e);
   }
  kkk++;

  c[0]=1;

  rc=1;
  chart=0;
  dir=0;

  while(rc&&(kmax==-1||kkk<kmax))
   {
    for(i=0;i<k;i++)s[i]=0.;
    s[dir]=h;

    unew=MFIMFVectorFactory(M,e);
    MFIMFProjectFromCenter(M,MFAtlasChartCenter(A,chart,e),MFAtlasChartTangentSpace(A,chart,e),S,unew,e);
    Phi=MFIMFTangentSpaceWithGuess(M,unew,MFAtlasChartTangentSpace(A,chart,e),e);

    rc=MFAtlasAddChartToList(A,MFCreateChart(MFAtlasMF(A,e),unew,Phi,h,e),e);
    printf("Chart %3d, [%2d",rc,c[0]);for(i=1;i<k;i++)printf(",%2d",c[i]);printf("], dir=%d, starting from chart %3d\n",dir,chart);fflush(stdout);
    kkk++;

    MFFreeNVector(unew,e);
    MFFreeNKMatrix(Phi,e);

    c[0]++;

    i=0;
    while(i<k && !(c[i]<nx))
     {
      c[i]=0;
      if(i+1<k)c[i+1]++;
       else c[k-1]=nx;
      i++;
     }
    if(c[k-1]>nx-1)return;

    dir=0;
    while(c[dir]==0&&dir<k)dir++;

    chart=0;
    for(i=0;i<k;i++)
     {
      chart=c[k-1-i]+nx*chart;
      if(k-1-i==dir)chart--;
     }
   }

  return;
 }

void MFRheinboldtCloseAtlas(MFContinuationMethod H, MFAtlas A, MFErrorHandler e)
 {
  static char RoutineName[]={"MFRheinboldtCloseAtlas"};
  FILE *fid;
  char *name;
  char *string;
  MFPolytope cube;
  int chart;
  int i,k,nx;
  int rc;
  int *c;

  string=MFRheinboldtGetFilename(H,e);
  name=(char*)malloc((strlen(string)+10)*sizeof(char));

#ifndef MFNOSAFETYNET
  if(name==NULL)
   {
    sprintf(MFRheinboldtErrorMsg,"Out of space trying to allocate %d bytes.",(strlen(string)+10)*sizeof(char));
    MFSetError(e,12,RoutineName,MFRheinboldtErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  strcpy(name,string);
  strcat(name,".plotfile");
  fid=fopen(name,"w");
  free(name);

  k=MFAtlasK(A,e);
  nx=MFRheinboldtGetIntegerParameter(H,"n",e);

  fprintf(fid,"Dimension of vertices, %d\n",MFIMFProjectToDraw(MFAtlasMF(A,e),MFAtlasChartCenter(A,0,e),NULL,e));
  fprintf(fid,"Dimension of manifold, %d\n",k);

  c=(int*)malloc(k*sizeof(int));

#ifndef MFNOSAFETYNET
  if(c==NULL)
   {
    sprintf(MFRheinboldtErrorMsg,"Out of space trying to allocate %d bytes.",k*sizeof(int));
    MFSetError(e,12,RoutineName,MFRheinboldtErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  cube=MFCreateHyperCubeAtOrigin(k,1.,e);
  for(i=0;i<k;i++)c[i]=0;

  rc=1;
  while(rc)
   {
    chart=0;
    for(i=0;i<k;i++)chart=c[k-1-i]+nx*chart;

    MFRheinboldtWriteChartToPlotFile(fid,A,cube,nx,chart,c,e);

    c[0]++;

    i=0;
    while(i<k-1 && !(c[i]<nx-1))
     {
      c[i]=0;
      c[i+1]++;
      i++;
     }
    if(c[k-1]>nx-2)rc=0;
   }

  fclose(fid);
  free(c);
  MFFreePolytope(cube,e);

  return;
 }

void MFRheinboldtFlushAtlas(MFContinuationMethod H, MFAtlas S, MFErrorHandler e)
 {
  static char RoutineName[]={"MFFlushAtlas"};

  return;
 }

/*!    \fn int MFRheinboldtSetIntegerParameter(MFImplicitMF M, char *parameterName, int value)
 *     \brief Allows the user to set integer parameters.
 *
 * Legal integer parameter names and their default values.
 *  <ul>
 *    <li>                   verbose                      default:    0
 *    <li>                   n                            default:    10
 *  </ul>
 *     \param M An Rheinboldt solution manifold.
 *     \param parameterName Which parameter to set.
 *     \param value The new value.
 *     \returns FALSE if the parameter name does not match a parameter.
 */
int MFRheinboldtSetIntegerParameter(MFContinuationMethod M, char *parameterName, int value, MFErrorHandler e)
 {
  static char RoutineName[]={"MFRheinboldtSetIntegerParameter"};
  struct MFRheinboldtParmBlockST *data;

  if(strcmp(MFContinuationMethodGetType(M,e),"Rheinboldt"))
   {
    sprintf(MFRheinboldtErrorMsg,"Parameter 1 must be a Rheinboldts algorithm");
    MFSetError(e,4,RoutineName,MFRheinboldtErrorMsg,__LINE__,__FILE__);
    return 0;
   }

  data=(struct MFRheinboldtParmBlockST*)MFContinuationMethodGetParmBlock(M,e);

  if(!strcmp(parameterName,"verbose"))data->verbose=value;
   else if(!strcmp(parameterName,"n"))data->n=value;
   else if(!strcmp(parameterName,"maxCharts"))data->maxCharts=value;
   else return 0;

  return 1;
 }

/*!    \fn int MFRheinboldtSetRealParameter(MFContinuationMethod M, char *parameterName, double value)
 *     \brief Allows the user to set real valued parameters.
 *
 * Legal real parameter names and their default values.
 *   <ul>
 *     <li>                   delta                        default:    0.1
 *   </ul>
 *    \param M An Rheinboldt solution manifold.
 *    \param parameterName Which parameter to set.
 *    \param value The new value.
 *    \returns FALSE if the parameter name does not match a parameter.
 */
int MFRheinboldtSetRealParameter(MFContinuationMethod M, char *parameterName, double value, MFErrorHandler e)
 {
  static char RoutineName[]={"MFRheinboldtSetRealParameter"};
  struct MFRheinboldtParmBlockST *data;

  if(strcmp(MFContinuationMethodGetType(M,e),"Rheinboldt"))
   {
    sprintf(MFRheinboldtErrorMsg,"Parameter 1 must be an Rheinboldt");
    MFSetError(e,4,RoutineName,MFRheinboldtErrorMsg,__LINE__,__FILE__);
    return 0;
   }

  data=(struct MFRheinboldtParmBlockST*)MFContinuationMethodGetParmBlock(M,e);

  if(!strcmp(parameterName,"delta"))data->delta=value;
   else return 0;

  return 1;
 }

/*!    \fn int MFRheinboldtGetIntegerParameter(MFContinuationMethod M, char *parameterName)
 *     \brief Allows the user to set integer parameters. These are usually read from a file. Instead, default values
 *            are used when the Rheinboldt is created, and the user may change them.
 *
 * Legal integer parameter names.
 *  <ul>
 *    <li>                   verbose                      1 if output wanted, 0 OW
 *    <li>                   n                            The number of mesh points in any one direction.
 *  </ul>
 *
 *     \param M An Rheinboldt solution manifold.
 *     \param parameterName Which parameter value to retreive. A warning is issued if the parameter name does not match a parameter.
 *     \returns The current value of the parameter.
 */
int MFRheinboldtGetIntegerParameter(MFContinuationMethod M, char *parameterName, MFErrorHandler e)
 {
  static char RoutineName[]={"MFRheinboldtGetIntegerParameter"};
  struct MFRheinboldtParmBlockST *data;

  if(strcmp(MFContinuationMethodGetType(M,e),"Rheinboldt"))
   {
    sprintf(MFRheinboldtErrorMsg,"Parameter 1 must be an Rheinboldt");
    MFSetError(e,4,RoutineName,MFRheinboldtErrorMsg,__LINE__,__FILE__);
    return 0;
   }

  data=(struct MFRheinboldtParmBlockST*)MFContinuationMethodGetParmBlock(M,e);

  if(!strcmp(parameterName,"verbose"))return data->verbose;
   else if(!strcmp(parameterName,"n"))return data->n;
   else if(!strcmp(parameterName,"maxCharts"))return data->maxCharts;
   else{
    sprintf(MFRheinboldtErrorMsg,"Parameter %s is not a integer valued parameter",parameterName);
    MFSetError(e,4,RoutineName,MFRheinboldtErrorMsg,__LINE__,__FILE__);
    return 0;
   }

  return 0;
 }

/*!    \fn double MFRheinboldtGetRealParameter(MFContinuationMethod M, char *parameterName)
 *     \brief Allows the user to set real valued parameters. These are usually read from a file. Instead, default values
 *            are used when the Rheinboldt is created, and the user may change them.
 *
 * Legal real parameter names.
 *    <ul>
 *     <li>                   delta     The mesh spacing in any one direction.
 *    </ul>
 *     \param M An Rheinboldt solution manifold.
 *     \param parameterName Which parameter value to retreive. A warning is issued if the parameter name does not match a parameter.
 *     \returns The current value of the parameter.
 */
double MFRheinboldtGetRealParameter(MFContinuationMethod M, char *parameterName, MFErrorHandler e)
 {
  static char RoutineName[]={"MFRheinboldtGetRealParameter"};
  struct MFRheinboldtParmBlockST *data;

#ifdef MFNOCONFIDENCE
  if(strcmp(MFContinuationMethodGetType(M,e),"Rheinboldt"))
   {
    sprintf(MFRheinboldtErrorMsg,"Parameter 1 must be an Rheinboldt");
    MFSetError(e,4,RoutineName,MFRheinboldtErrorMsg,__LINE__,__FILE__);
    return 0;
   }
#endif

  data=(struct MFRheinboldtParmBlockST*)MFContinuationMethodGetParmBlock(M,e);
  if(!strcmp(parameterName,"delta"))return data->delta;

#ifdef MFNOCONFIDENCE
   sprintf(MFRheinboldtErrorMsg,"Parameter %s is not a real valued parameter",parameterName);
   MFSetError(e,4,RoutineName,MFRheinboldtErrorMsg,__LINE__,__FILE__);
   return 0.;
#endif

  return 0.;
 }

void MFRheinboldtWriteChartToPlotFile(FILE *fid, MFAtlas A, MFPolytope P, int nx, int chart, int *c0, MFErrorHandler e)
 {
  static char RoutineName[]={"MFRheinboldtWriteChartToPlotFile"};
  int d;
  double *x;
  MFKVector s;
  double *ps;
  MFNVector u;
  MFChart thisMethod;
  int i,j,k,l;
  int n,m,ne,ie;
  int nv;
  int c;
  int *index;
  int *indx;
  int verbose=0;
  MFImplicitMF M;
  int bnd,sing;
  int chartDelta;

  k=MFAtlasK(A,e);
  if(MFPolytopeNumberOfVertices(P,e)<=k)return;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s\n",RoutineName);fflush(stdout);}
#endif

  if(fid==NULL){printf("%s - fid is NULL %s(%d)\n",RoutineName,__LINE__,__FILE__);fflush(stdout);abort();}
  if(A==NULL){printf("%s - A is NULL %s(%d)\n",RoutineName,__LINE__,__FILE__);fflush(stdout);abort();}
  if(thisMethod==NULL){printf("%s - Chart is NULL %s(%d)\n",RoutineName,__LINE__,__FILE__);fflush(stdout);abort();}

  M=MFAtlasMF(A,e);

  d=MFIMFProjectToDraw(M,MFAtlasChartCenter(A,chart,e),NULL,e);

  s=MFCreateKVector(k,e);
  ps=MFKV_CStar(s,e);
  x=(double*)malloc(d*sizeof(double));

#ifndef MFNOSAFETYNET
  if(x==NULL)
   {
    sprintf(MFRheinboldtErrorMsg,"Out of memory trying to allocate %d bytes",d*sizeof(double));
    MFSetError(e,12,RoutineName,MFRheinboldtErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  if(MFPolytopeNumberOfVertices(P,e)<=k)
   {
    MFFreeKVector(s,e);
    free(x);
    return;
   }

  n=MFPolytopeNumberOfVertices(P,e);

  indx=(int*)malloc(k*sizeof(int));

#ifndef MFNOSAFETYNET
      if(indx==NULL)
       {
        sprintf(MFRheinboldtErrorMsg,"Out of memory trying to allocate %d bytes",k*sizeof(int));
        MFSetError(e,12,RoutineName,MFRheinboldtErrorMsg,__LINE__,__FILE__);
        MFErrorHandlerOutOfMemory(e);
        return;
       }
#endif

  ne=0;
  for(i=0;i<n-1;i++)
    for(j=i+1;j<n;j++)
      if(MFPolytopeIntersectIndexSets(P,i,j,indx,e)>=k-1)ne++;

  m=MFPolytopeNumberOfFaces(P,e);
  bnd=0;
  sing=0;

  thisMethod=MFAtlasChart(A,chart,e);

  nv=n+1;
  fprintf(fid,"Polyhedron %d, R=%lf, %d vertices, %d edges, %d faces, boundary %d, singular %d\n",chart,MFChartRadius(MFAtlasChart(A,chart,e),e),nv,ne,m,bnd,sing);
  for(i=0;i<n;i++)
   {
    MFPolytopeVertex(P,i,s,e);

    chartDelta=0;
    for(j=0;j<k;j++)
     {
      if(ps[k-1-j]>0.)chartDelta=c0[k-1-j]+1+nx*chartDelta;
       else chartDelta=c0[k-1-j]+nx*chartDelta;
     }

    u=MFAtlasChartCenter(A,chartDelta,e);
    MFIMFProjectToDraw(M,u,x,e);
    fprintf(fid,"Vertex %d (%lf",i,x[0]);
    for(j=1;j<d;j++)fprintf(fid,",%lf",x[j]);
    fprintf(fid,")");
    m=MFPolytopeVertexNumberOfIndices(P,i,e);
    index=MFPolytopeVertexIndexSet(P,i,e);
    fprintf(fid,", %d [%d",m,index[0]);
    for(j=1;j<m;j++)fprintf(fid,",%d",index[j]);
    fprintf(fid,"]\n");
   }

/* Center */

  MFIMFProjectToDraw(M,MFAtlasChartCenter(A,chart,e),x,e);
  fprintf(fid,"Vertex %d (%lf",i,x[0]);
  for(j=1;j<d;j++)fprintf(fid,",%lf",x[j]);
  fprintf(fid,")");
  fprintf(fid,", %d [ ]\n",0);

  MFFreeKVector(s,e);
  free(x);
  fflush(fid);

  ie=0;
  for(i=0;i<n-1;i++)
   {
    for(j=i+1;j<n;j++)
     {
      if((m=MFPolytopeIntersectIndexSets(P,i,j,indx,e))>=k-1)
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

  m=MFPolytopeNumberOfFaces(P,e);
  for(i=0;i<m;i++)
   {
    j=MFPolytopeFaceIndex(P,i,e);
    fprintf(fid,"Face %d neighbor %d\n",j,-1);
   }

  fflush(fid);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("done %s\n",RoutineName);fflush(stdout);}
#endif

  return;
 }

#ifdef __cplusplus
}
#endif
