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

static char *id="@(#) $Id: MFPolytope.c,v 1.3 2007/02/13 01:22:34 mhender Exp $";

static char MFPolytopeErrorMsg[256]="";

#include <math.h>
#include <MFKVector.h>
#include <MFPolytope.h>
#include <MFPrint.h>
#include <MFErrorHandler.h>
#include <stdlib.h>
#include <string.h>

#ifdef __cplusplus
 extern "C" {
#endif

/* Internal */

int MFPolytopeVertexLargestIndex(MFPolytope,int,MFErrorHandler);
int MFPolytopeVertexSmallestIndex(MFPolytope,int,MFErrorHandler);
int MFPolytopeVertexNumberOfIndices(MFPolytope,int,MFErrorHandler);
int MFPolytopeIntersectIndexSets(MFPolytope,int,int,int*,MFErrorHandler);
int MFPolytopeAddVertexIndexIntoSet(MFPolytope,int,int,MFErrorHandler);
int MFTestPolytope(MFPolytope,MFErrorHandler);
int MFPolytopeSmallestVertexIndex(MFPolytope,MFErrorHandler);
int *MFPolytopeVertexIndexSet(MFPolytope,int,MFErrorHandler);

int MFSolveFull(int,double*,double*,MFErrorHandler);

static double MFPolytopeVertexLargestVertexRadius(MFPolytope,MFErrorHandler);
static int MFPolytopeLargestVertexIndex(MFPolytope,MFErrorHandler);
static int MFPolytopeUnionOfIndexSets(MFPolytope,int,int,int*,MFErrorHandler);
static void MFPolytopeUpdateFaceList(MFPolytope,MFErrorHandler);
static void MFPolytopeMergeCloseVertices(MFPolytope,double,MFErrorHandler);
static void MFPolytopeMergeVertices(MFPolytope,int,int,MFErrorHandler);
static int MFPolytopeVerticesOnSameEdge(MFPolytope,int,int,MFErrorHandler);
static int MFPolytopeVerticesOnEdgeWithIndexLessThan(MFPolytope,int,int,int,MFErrorHandler);
static int MFTestVertexIndexOrder(int,int*,MFErrorHandler);
static int MFPolytopeTestVertexList(MFPolytope,MFErrorHandler);
static int MFPolytopeAddVertex(MFPolytope,MFErrorHandler);
static void MFRecomputePolytopeFromFaces(MFPolytope,MFErrorHandler);

struct MFPolytopeSt{
                int k;                /* Dimension of the vertices */
                int n;                /* Number of vertices */
                double R;             /* Radius */
                int m;                /* Space for vertices */
                double *v;            /* vertices */
                int *nIndices;        /* Number of indices for each vertex */
                int *mIndices;        /* Space for indices for each vertex */
                int **indices;        /* Indices of each vertex */
                int *mark;            /* Mark on each vertex */

                int nFaces;
                int mFaces;
                int *face;          /* Index of face */
                int *nFaceV;        /* Number of vertices on each face */
                MFKVector *faceN;   /* Normal of face */
                double *faceO;      /* Origin of face */
               };

#define MFPolytopeIncrementToAllocateToVertexList 10
#define MFPolytopeIncrementToAllocateToIndexSet 2
#define MFPolytopeIncrementToAllocateToFaceList 2

MFPolytope MFCreateHyperCubeAtOrigin(int k, double R, MFErrorHandler e)
 {
  static char RoutineName[]={"MFCreateHyperCubeAtOrigin"};
  int i,j,l;
  double *c;
  MFPolytope result;
  int carry;
  int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose)printf("%s\n",RoutineName);fflush(stdout);
#endif

/* Allocate storage */

  result=(MFPolytope)malloc(sizeof(struct MFPolytopeSt));

#ifndef MFNOSAFETYNET
  if(result==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct MFPolytopeSt));
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

#ifdef MFALLOWVERBOSE
  if(verbose)printf("MFCreateHyperCubeAtOrigin 0x%8.8x\n",result);
#endif

  result->k=k;
  
  result->n=1;for(i=0;i<k;i++)result->n*=2;
  result->m=result->n;
  result->v=(double*)malloc(result->n*k*sizeof(double));

#ifndef MFNOSAFETYNET
  if(result->v==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",result->n*k*sizeof(double));
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  for(i=0;i<result->n*k;i++)result->v[i]=0.;

  result->nIndices=(int*)malloc(result->n*sizeof(int));

#ifndef MFNOSAFETYNET
  if(result->nIndices==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",result->n*sizeof(int));
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  result->mIndices=(int*)malloc(result->n*sizeof(int));

#ifndef MFNOSAFETYNET
  if(result->mIndices==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",result->n*sizeof(int));
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  result->indices=(int**)malloc(result->n*sizeof(int*));

#ifndef MFNOSAFETYNET
  if(result->indices==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",result->n*sizeof(int*));
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

#ifdef MFALLOWVERBOSE
  if(verbose)printf("   result(0x%8.8x)->indices=0x%8.8x, %d int*s\n",result,result->indices,result->n);
#endif

  for(i=0;i<result->n;i++)
   {
    result->nIndices[i]=0;
    result->mIndices[i]=k;
    result->indices[i]=(int*)malloc(k*sizeof(int));

#ifndef MFNOSAFETYNET
    if(result->indices[i]==NULL)
     {
      sprintf(MFPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",k*sizeof(int));
      MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return NULL;
     }
#endif

#ifdef MFALLOWVERBOSE
    if(verbose)printf("   result(0x%8.8x)->indices[%d]=0x%8.8x, %d ints\n",result,i,result->indices[i],k);
#endif

   }

  result->mark=(int*)malloc(result->n*sizeof(int));

#ifndef MFNOSAFETYNET
  if(result->mark==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",result->n*sizeof(int));
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  c=(double*)malloc(k*sizeof(double));

#ifndef MFNOSAFETYNET
  if(c==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",k*sizeof(double));
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  for(j=0;j<k;j++)c[j]=-R;

  for(i=0;i<result->n;i++)
   {
    l=i*k;

#ifdef MFALLOWVERBOSE
    if(verbose)printf("vertex %d\n",i);
#endif

    for(j=0;j<k;j++)result->v[j+l]=c[j];

#ifdef MFALLOWVERBOSE
    if(verbose){printf(" v=[");for(j=0;j<k;j++){if(j>0)printf(",");printf("%lf",result->v[j+l]);}printf("]\n");}
#endif

    result->nIndices[i]=k;
    result->mIndices[i]=k;
    result->mark[i]=0;
    for(j=0;j<k;j++)
     {
      if(c[j]<0)(result->indices[i])[j]=2*j;
       else (result->indices[i])[j]=2*j+1;
     }

#ifdef MFALLOWVERBOSE
    if(verbose){printf(" i=[");for(j=0;j<k;j++){if(j>0)printf(",");printf("%d",(result->indices[i])[j]);}printf("]\n");}
#endif

    carry=1;
    j=0;
    while(carry&&j<k)
     {
      if(c[j]<0)
       {
        c[j]=R;
        carry=0;
       }else{
        c[j]=-R;
        carry=1;
        j++;
       }
     }
   }

  if(c!=NULL)free(c);

  result->nFaces=2*k;
  result->mFaces=result->nFaces;
  result->face=(int*)malloc(result->mFaces*sizeof(int));

#ifndef MFNOSAFETYNET
  if(result->face==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",result->mFaces*sizeof(int));
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  result->nFaceV=(int*)malloc(result->mFaces*sizeof(int));

#ifndef MFNOSAFETYNET
  if(result->nFaceV==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",result->mFaces*sizeof(int));
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  result->faceN=(MFKVector*)malloc(result->mFaces*sizeof(MFKVector));

#ifndef MFNOSAFETYNET
  if(result->faceN==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",result->mFaces*sizeof(MFKVector));
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  for(i=0;i<result->mFaces;i++)result->faceN[i]=NULL;

  result->faceO=(double*)malloc(result->mFaces*sizeof(double));

#ifndef MFNOSAFETYNET
  if(result->faceO==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",result->mFaces*sizeof(double));
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  for(i=0;i<2*k;i++)
   {
    result->face[i]=i;
    result->nFaceV[i]=0;
    result->faceN[i]=MFCreateKVector(k,e);
    if(i%2==0)
      MFKVSetC(result->faceN[i],i/2,-1,e);
     else
      MFKVSetC(result->faceN[i],i/2, 1,e);
    result->faceO[i]=R;
   }

  for(i=0;i<result->n;i++)
   {
    for(j=0;j<result->nIndices[i];j++)
     {
      result->nFaceV[(result->indices[i])[j]]++;
     }
   }

  result->R=R*sqrt(1.*k);

  return result;
 }

MFPolytope MFCreateSimplexAtOrigin(int k, double R, MFErrorHandler e)
 {
  static char RoutineName[]={"MFCreateSimplexAtOrigin"};
  int i,j,l;
  MFPolytope result;
  int l0,m;
  double d;
  double h;
  double Delta;
  double Rd;
  double Rdm;
  int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose)printf("%s\n",RoutineName);fflush(stdout);
#endif

/* Allocate storage */

  result=(MFPolytope)malloc(sizeof(struct MFPolytopeSt));

#ifndef MFNOSAFETYNET
  if(result==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct MFPolytopeSt));
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

#ifdef MFALLOWVERBOSE
  if(verbose)printf("MFCreateSimplexAtOrigin 0x%8.8x\n",result);
#endif

  result->k=k;
  
  result->n=k+1;
  result->m=result->n;
  result->v=(double*)malloc(result->n*k*sizeof(double));

#ifndef MFNOSAFETYNET
  if(result->v==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",result->n*k*sizeof(double));
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  for(i=0;i<result->n*k;i++)result->v[i]=0.;

  result->v[0+k*0]=-R;for(j=1;j<k;j++)result->v[j+k*0]=0.;
  result->v[0+k*1]= R;for(j=1;j<k;j++)result->v[j+k*1]=0.;
  Rdm=R;

  for(i=2;i<k+1;i++)
   {
    Delta=sqrt(4*R*R-Rdm*Rdm);
    Rd=2*R*R/Delta;

    for(l=0;l<i;l++)result->v[i-1+k*l]=Rd-Delta;
    for(j=0;j<i-1;j++)result->v[j+k*i]=0.;result->v[i-1+k*i]=Rd;for(j=i;j<k;j++)result->v[j+k*i]=0.;

    Rdm=Rd;
   }

/* Scale to contain ball radius R */

  d=0.;h=0.;
  for(j=0;j<k;j++){d+=-result->v[j]*result->v[j+k];h+=result->v[j]*result->v[j];}
  for(i=0;i<k*(k+1);i++)result->v[i]=result->v[i]*R*sqrt(h)/d;

  result->nIndices=(int*)malloc(result->n*sizeof(int));

#ifndef MFNOSAFETYNET
  if(result->nIndices==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",result->n*sizeof(int));
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  result->mIndices=(int*)malloc(result->n*sizeof(int));

#ifndef MFNOSAFETYNET
  if(result->mIndices==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",result->n*sizeof(int));
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  result->mark=(int*)malloc(result->n*sizeof(int));

#ifndef MFNOSAFETYNET
  if(result->mark==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",result->n*sizeof(int));
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  result->indices=(int**)malloc(result->n*sizeof(int*));

#ifndef MFNOSAFETYNET
  if(result->indices==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",result->n*sizeof(int*));
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

#ifdef MFALLOWVERBOSE
  if(verbose)printf("   result(0x%8.8x)->indices=0x%8.8x, %d int*s\n",result,result->indices,result->n);
#endif

  for(i=0;i<result->n;i++)
   {
    result->nIndices[i]=k;
    result->mIndices[i]=k;
    result->mark[i]=0;
    result->indices[i]=(int*)malloc(k*sizeof(int));

#ifndef MFNOSAFETYNET
    if(result->indices[i]==NULL)
     {
      sprintf(MFPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",k*sizeof(int));
      MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return NULL;
     }
#endif

#ifdef MFALLOWVERBOSE
    if(verbose)printf("   result(0x%8.8x)->indices[%d]=0x%8.8x, %d ints\n",result,i,result->indices[i],k);
#endif

   }

  result->nFaces=k+1;
  result->mFaces=result->nFaces;
  result->face=(int*)malloc(result->mFaces*sizeof(int));

#ifndef MFNOSAFETYNET
  if(result->face==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",result->mFaces*sizeof(int));
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  result->nFaceV=(int*)malloc(result->mFaces*sizeof(int));

#ifndef MFNOSAFETYNET
  if(result->nFaceV==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",result->mFaces*sizeof(int));
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  result->faceN=(MFKVector*)malloc(result->mFaces*sizeof(MFKVector));

#ifndef MFNOSAFETYNET
  if(result->faceN==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",result->mFaces*sizeof(MFKVector));
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  for(i=0;i<result->mFaces;i++)result->faceN[i]=NULL;

  result->faceO=(double*)malloc(result->mFaces*sizeof(double));

#ifndef MFNOSAFETYNET
  if(result->faceO==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",result->mFaces*sizeof(double));
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  d=0;
  for(j=0;j<k;j++){d+=result->v[j]*result->v[j];}
  d=1./sqrt(d);
  for(i=0;i<k+1;i++)
   {
    result->face[i]=i;
    result->nFaceV[i]=0;
    result->faceN[i]=MFCreateKVector(k,e);
    result->faceO[i]=R;
    for(j=0;j<k;j++)
     {
      MFKVSetC(result->faceN[i],j,-result->v[j+k*i]*d,e);
      result->nFaceV[j]=k;
     }
    m=0;
    for(j=0;j<k+1;j++)
     {
      if(j!=i)
       {
        (result->indices[i])[m]=j;
        m++;
       }
     }
   }

  result->R=R;

  return result;
 }

void MFFreePolytope(MFPolytope P, MFErrorHandler e)
 {
  static char RoutineName[]={"MFFreePolytope"};
  int i;

#ifdef MFNOCONFIDENCE
  if(P==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"Polytope (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  if(P->v!=NULL){free(P->v);P->v=NULL;}
  if(P->nIndices!=NULL){free(P->nIndices);P->nIndices=NULL;}
  if(P->mIndices!=NULL){free(P->mIndices);P->mIndices=NULL;}
  if(P->mark!=NULL){free(P->mark);P->mark=NULL;}
  for(i=0;i<P->n;i++)
   {
    if(P->indices[i]!=NULL){free(P->indices[i]);P->indices[i]=NULL;}
   }
  if(P->indices!=NULL){free(P->indices);P->indices=NULL;}
  if(P->face!=NULL){free(P->face);P->face=NULL;}
  if(P->nFaceV!=NULL){free(P->nFaceV);P->nFaceV=NULL;}
  if(P->faceN!=NULL)
   {
    for(i=0;i<P->mFaces;i++) /* was nFaces */
     {
      if(P->faceN[i]!=NULL)
        MFFreeKVector(P->faceN[i],e);
      P->faceN[i]=NULL;
     }
    free(P->faceN);P->faceN=NULL;
   }
  if(P->faceO!=NULL){free(P->faceO);P->faceO=NULL;}
  free(P);

  return;
 }

void MFSubtractHalfSpaceFromPolytope(MFPolytope P,int index, MFKVector nrm, double on, MFErrorHandler e)
 {
  static char RoutineName[]={"MFSubtractHalfSpaceFromPolytope"};
  int i,j,k,l;
  int n;
  int flag;
  double *d;
  int *inter;
  int n1,n2;
  double *vi,*vj,*v;
  double t;
  int nold,m,newV;
  int ln,li,lj;
  int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("MFSubtractHalfSpaceFromPolytope 0x%8.8x\n",P);
    printf("  Polytope\n");
    MFPrintPolytope(stdout,P,e);
    printf("  HalfSpace\n");

    fprintf(stdout,"                     (%d) x.",index);
    fprintf(stdout,"(");
    for(j=0;j<P->k;j++)
     {
      if(j>0)fprintf(stdout,",");
      t=MFKV_C(nrm,j,e);
      fprintf(stdout,"%le",t);
     }
    fprintf(stdout,")-%le",on);
    fprintf(stdout,"<=0   nrm=0x%8.8x\n",nrm);
   }
#endif

#ifdef MFNOCONFIDENCE
  if(P==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"Polytope (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(nrm==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"nrm (argument 3) is NULL");
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    return;
   }

  for(i=0;i<P->nFaces;i++)
   {
    if(P->faceN[i]!=NULL)
     {
      t=fabs(on-P->faceO[i]);
      for(j=0;j<P->k;j++) t+=fabs(MFKV_C(nrm,j,e)-MFKV_C(P->faceN[i],j,e));
      if(t<1.e-7)
       {
        sprintf(MFPolytopeErrorMsg,"Face being subtracted is already present");
        MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
        return;
       }
     }
   }

  if(!MFPolytopeTestVertexList(P,e))
   {
    sprintf(MFPolytopeErrorMsg,"Vertex List is bas.");
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  if(P->n==0)
   {
    return;
   }

  for(i=0;i<P->nFaces;i++)
   if(P->face[i]==index)
#ifndef MFNOSAFETYNET
    {
     sprintf(MFPolytopeErrorMsg,"Subtracting a halfspace index already present",index);
     MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
     return;
    }
#endif

  d=(double*)malloc(P->n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(d==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",P->n*sizeof(double));
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  t=fabs(on);
  for(j=0;j<P->k;j++)if(fabs(MFKV_C(nrm,j,e))>t)t=fabs(MFKV_C(nrm,j,e));

#ifdef MFALLOWVERBOSE
  if(verbose)printf(" inf norm of (on,norm) is %le\n",t);
#endif

  for(i=0;i<P->n;i++)
   {
    d[i]=-on;
    for(j=0;j<P->k;j++)d[i]+=P->v[j+P->k*i]*MFKV_C(nrm,j,e);
    d[i]=d[i]/t;

#ifdef MFALLOWVERBOSE
    if(verbose)printf("vertex %d is on side %le\n",i,d[i]);
#endif

   }

  nold=P->n;
  for(i=0;i<nold;i++)
   {
    n1=P->nIndices[i];
    for(j=i+1;j<nold;j++)
     {

/* This number should be relative to the longest edge of the polytope */

      if(fabs(d[i])>1.e-14 && fabs(d[i])>1.e-14 &&
            (d[i]<0&&d[j]>0 || d[i]>0 && d[j]<0))
       {
#ifdef MFALLOWVERBOSE
        if(verbose)
         {
          printf("  thisP pair (%d,%d) straddles the new plane\n",i,j);
          printf("    %d indices (",i);
          for(l=0;l<P->nIndices[i];l++){if(l>0)printf(",");printf("%d",(P->indices[i])[l]);}
          printf(") d=%le\n",d[i]);
          printf("    %d indices (",j);
          for(l=0;l<P->nIndices[j];l++){if(l>0)printf(",");printf("%d",(P->indices[j])[l]);}
          printf(") d=%le\n",d[j]);
          t=0.;
          for(l=0;l<P->k;l++)t+=(P->v[l+P->k*i]-P->v[l+P->k*j])*(P->v[l+P->k*i]-P->v[l+P->k*j]);
          t=sqrt(t);
          printf("    edge is length %le\n",t);fflush(stdout);
         }
#endif

        n2=P->nIndices[j];
        m=n1;if(n2>m)m=n2;

        if(m>0)
         {
          inter=(int*)malloc(m*sizeof(int));

#ifndef MFNOSAFETYNET
          if(inter==NULL)
           {
            sprintf(MFPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",m*sizeof(int));
            MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
            MFErrorHandlerOutOfMemory(e);
            return;
           }
#endif

          n=MFPolytopeIntersectIndexSets(P,i,j,inter,e);
  
#ifdef MFALLOWVERBOSE
          if(verbose)
           {
            printf("    intersection indices (");
            for(l=0;l<n;l++){if(l>0)printf(",");printf("%d",inter[l]);}
            printf(")\n");
           }
#endif

         }else{
          inter=NULL;
          n=0;
         }

        if(n<(P->k)-1)
         {
#ifdef MFALLOWVERBOSE
          if(verbose)printf("  Not on an edge\n");
#endif

          if(inter!=NULL)free(inter);
          inter=NULL;
         }else{

#ifdef MFALLOWVERBOSE
          if(verbose)printf("  On an edge\n");
#endif

          newV=MFPolytopeAddVertex(P,e);
          t=-d[i]/(d[j]-d[i]);
          li=P->k*i;
          lj=P->k*j;
          ln=P->k*newV;

#ifdef MFALLOWVERBOSE
          if(verbose)printf("  d[%d]=%le, d[%d]=%le, t=%le\n",i,d[i],j,d[j],t);
#endif

          for(l=0;l<P->k;l++)
           {
            P->v[l+ln]=P->v[l+li]+t*(P->v[l+lj]-P->v[l+li]);

#ifdef MFALLOWVERBOSE
/*          if(verbose)printf("  v[%d]_%d=(vi[%d]_%d)-t*(vj[%d]_%d-vi[%d]_%d)=%le-(%le)*((%le)-(%le))=%le=%le\n",newV,l,i,l,j,l,i,l,P->v[l+li],t,P->v[l+lj]-P->v[l+li],P->v[l+li]-t*(P->v[l+lj]-P->v[l+li]),P->v[l+ln]);*/
#endif

           }
          P->nIndices[newV]=n;
          P->mIndices[newV]=m;
          P->mark[newV]=0;
          P->indices[newV]=inter;
          d=(double*)realloc((void*)d,(P->m)*sizeof(double));

#ifdef MFNOCONFIDENCE
          if(d==NULL)
           {
            sprintf(MFPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",P->m*sizeof(double));
            MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
            MFErrorHandlerOutOfMemory(e);
            return;
           }
#endif

          d[newV]=-100.;
          MFPolytopeAddVertexIndexIntoSet(P,newV,index,e);

#ifdef MFALLOWVERBOSE
          if(verbose)
           {
            printf("> adding new vertex P(%0x8.8x)->indices[%d]=[",P,newV);
            for(l=0;l<P->nIndices[newV];l++)
             {
              if(l>0)printf(",");
              printf("%d",(P->indices[newV])[l]);
             }
            printf("]\n");
           }
#endif
         }
       }
     }
   }

  for(i=0;i<P->n;i++)
   {
    if(d[i]>0)
     {

#ifdef MFALLOWVERBOSE
      if(verbose){printf("Vertex %d is on the wrong side, d[%d]=%le, removing it.\n",i,i,d[i]);fflush(stdout);}
#endif

      if(i<P->n-1)
       {
        for(j=0;j<P->k;j++)P->v[j+P->k*i]=P->v[j+P->k*(P->n-1)];

#ifdef MFALLOWVERBOSE
        if(verbose)
         {
          printf("set P->indices[%d]=P->indices[%d]=[",i,P->n-1);
          for(l=0;l<P->nIndices[P->n-1];l++)
           {
            if(l>0)printf(",");
            printf("%d",(P->indices[P->n-1])[l]);
           }
          printf("] d=%le\n",d[P->n-1]);
         }
#endif

        if(P->indices[i]!=NULL)free(P->indices[i]);
        P->indices[i]=P->indices[P->n-1];
        P->nIndices[i]=P->nIndices[P->n-1];
        P->mIndices[i]=P->mIndices[P->n-1];
        P->mark[i]=P->mark[P->n-1];
        d[i]=d[P->n-1];
       }else{
        if(P->indices[i]!=NULL)free(P->indices[i]);
       }
      P->indices[P->n-1]=NULL;
      P->nIndices[P->n-1]=0;
      P->mIndices[P->n-1]=0;
      P->mark[P->n-1]=0;
      d[P->n-1]=0.;
      (P->n)--;
      i--;
     }else if(d[i]==0)
     {

#ifdef MFALLOWVERBOSE
      if(verbose){printf("Vertex %d is on the hyperplane, adding index to it.\n",i);fflush(stdout);}
#endif

      MFPolytopeAddVertexIndexIntoSet(P,i,index,e);
     }else{

#ifdef MFALLOWVERBOSE
      if(verbose)
       {
        printf("Vertex %d is on the right side, d[%d]=%le, keeping it.\n",i,i,d[i]);
        printf("    P->indices[%d]=[",i);
        for(l=0;l<P->nIndices[i];l++)
         {
          if(l>0)printf(",");
          printf("%d",(P->indices[i])[l]);
         }
        printf("] d=%le\n",d[i]);
       }
#endif

     }
   }

  if(d!=NULL)free(d);

/* Now remove redundant half spaces */
/*   For each appearing index, count occurances. */
/*    If it appears fewer than k-1 times, remove it */
/*    else If dim(span(vertices))<k-1 times, remove it */

  if(P->nFaces>=P->mFaces)
   {
    P->mFaces+=MFPolytopeIncrementToAllocateToFaceList;

    P->face=(int*)realloc((void*)P->face,P->mFaces*sizeof(int));

#ifndef MFNOSAFETYNET
    if(P->face==NULL)
     {
      sprintf(MFPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",P->mFaces*sizeof(int));
      MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

    P->nFaceV=(int*)realloc((void*)P->nFaceV,P->mFaces*sizeof(int));

#ifndef MFNOSAFETYNET
    if(P->nFaceV==NULL)
     {
      sprintf(MFPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",P->mFaces*sizeof(int));
      MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

    P->faceN=(MFKVector*)realloc((void*)(P->faceN),P->mFaces*sizeof(MFKVector));

#ifndef MFNOSAFETYNET
    if(P->faceN==NULL)
     {
      sprintf(MFPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",P->mFaces*sizeof(MFKVector));
      MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

    for(i=P->nFaces;i<P->mFaces;i++)P->faceN[i]=NULL;
  
    P->faceO=(double*)realloc((void*)P->faceO,P->mFaces*sizeof(double));

#ifndef MFNOSAFETYNET
    if(P->faceO==NULL)
     {
      sprintf(MFPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",P->mFaces*sizeof(double));
      MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif
   }

  P->face[P->nFaces]=index;
  if(nrm==NULL)abort();

  P->faceN[P->nFaces]=nrm;
  MFRefKVector(nrm,e);
  P->faceO[P->nFaces]=on;
  P->nFaceV[P->nFaces]=0;
  P->nFaces++;

#ifdef MFALLOWVERBOSE
  if(verbose){MFPrintPolytope(stdout,P,e);printf("\n   -- call UpdateFaceList, P has %d vertices \n",P->n);}
#endif

  if(!MFPolytopeTestVertexList(P,e))
   {
    printf("%s line %d, file %s. Polytope vertex list is bad.\n",RoutineName,__LINE__,__FILE__);
    fflush(stdout);
/*  abort();*/
   }

  MFPolytopeUpdateFaceList(P,e);

#ifdef MFALLOWVERBOSE
  if(verbose)printf("   done MFSubtractHalfSpaceFromPolytope\n");
#endif

  return;
 }

/*----------------------------------------------------------------------*/

int MFPolytopeDimension(MFPolytope P, MFErrorHandler e)
 {
  static char RoutineName[]={"MFPolytopeDimension"};

#ifdef MFNOCONFIDENCE
  if(P==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"First argument is NULL");
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    return -1;
   }
#endif

  return P->k;
 }

int MFPolytopeNumberOfVertices(MFPolytope P, MFErrorHandler e)
 {
  static char RoutineName[]={"MFPolytopeNumberOfVertices"};

#ifdef MFNOCONFIDENCE
  if(P==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"First argument is NULL");
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    return -1;
   }
#endif

  return P->n;
 }

void MFPolytopeVertex(MFPolytope P,int i,MFKVector v, MFErrorHandler e)
 {
  static char RoutineName[]={"MFPolytopeVertex"};
  int j;

#ifdef MFNOCONFIDENCE
  if(P==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"First argument is NULL");
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(i<0 || i>P->n)
   {
    sprintf(MFPolytopeErrorMsg,"Invalid vertex number %d (must be in [0,%d) )",i,P->n);
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  for(j=0;j<P->k;j++)
    MFKVSetC(v,j,P->v[j+i*P->k],e);

  return;
 }

int MFPolytopeNumberOfVertexIndices(MFPolytope P,int i, MFErrorHandler e)
 {
  static char RoutineName[]={"MFPolytopeNumberOfVertexIndices"};

#ifdef MFNOCONFIDENCE
  if(P==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"First argument is NULL");
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    return -1;
   }
  if(i<0 || i>=P->n)
   {
    sprintf(MFPolytopeErrorMsg,"Invalid vertex number %d (must be in [0,%d) )",i,P->n);
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    return -1;
   }
#endif

  return P->nIndices[i];
 }

int MFPolytopeVertexIndex(MFPolytope P,int i,int j, MFErrorHandler e)
 {
  static char RoutineName[]={"MFPolytopeVertexIndex"};

#ifdef MFNOCONFIDENCE
  if(P==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"First argument is NULL");
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    return -1;
   }
  if(i<0 || i>=P->n)
   {
    sprintf(MFPolytopeErrorMsg,"Invalid vertex number %d (must be in [0,%d) )",i,P->n);
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    return -1;
   }
  if(P->indices==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"Polytope has NULL vertex index list");
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    return -1;
   }
  if(j<0 || j>=P->nIndices[i])
   {
    sprintf(MFPolytopeErrorMsg,"Invalid index number %d (must be in [0,%d) )",j,P->nIndices[i]);
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    return -1;
   }
#endif

  return (P->indices[i])[j];
 }

int *MFPolytopeVertexIndexSet(MFPolytope P,int i, MFErrorHandler e)
 {
  static char RoutineName[]={"MFPolytopeVertexIndexSet"};
  int j;

#ifdef MFNOCONFIDENCE
  if(P==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"First argument is NULL");
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    return NULL;
   }
  if(i<0 || i>=P->n)
   {
    sprintf(MFPolytopeErrorMsg,"Invalid vertex number %d (must be in [0,%d) )",i,P->n);
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    return NULL;
   }

  if(P->indices==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"Polytope has NULL vertex index list");
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    return NULL;
   }
#endif

  return P->indices[i];
 }

int MFPolytopeVertexNumberOfIndices(MFPolytope P,int i, MFErrorHandler e)
 {
  static char RoutineName[]={"MFPolytopeVertexNumberOfIndices"};
  int j;

#ifdef MFNOCONFIDENCE
  if(P==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"First argument is NULL");
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    return -1;
   }
  if(i<0 || i>=P->n)
   {
    sprintf(MFPolytopeErrorMsg,"Invalid vertex number %d (must be in [0,%d) )",i,P->n);
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    return -1;
   }

  if(P->nIndices==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"Polytope has NULL vertex index length list\n");
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    return -1;
   }
#endif

  return P->nIndices[i];
 }

/*----------------------------------------------------------------------*/

int MFPolytopeAddVertex(MFPolytope P, MFErrorHandler e)
 {
  static char RoutineName[]={"MFPolytopeAddVertex"};
  int i;
  int verbose=0;

#ifdef MFNOCONFIDENCE
  if(P==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"First argument is NULL");
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    return -1;
   }
#endif

#ifdef MFALLOWVERBOSE
  if(verbose)printf("MFPolytopeAddVertex 0x%8.8x\n",P);
#endif

  if(P->n>=P->m)
   {

#ifdef MFALLOWVERBOSE
    if(verbose)printf(" P->n(%d)>=P->m(%d), creating more space\n",P->n,P->m);
#endif

    P->m+=MFPolytopeIncrementToAllocateToVertexList;

    P->v=(double*)realloc((void*)P->v,P->k*P->m*sizeof(double));

#ifndef MFNOSAFETYNET
    if(P->v==NULL)
     {
      sprintf(MFPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",P->m*sizeof(double));
      MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return -1;
     }
#endif

    for(i=P->k*P->n;i<P->k*P->m;i++)P->v[i]=0.;

    P->nIndices=(int*)realloc((void*)P->nIndices,P->m*sizeof(int));

#ifndef MFNOSAFETYNET
    if(P->nIndices==NULL)
     {
      sprintf(MFPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",P->m*sizeof(int));
      MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return -1;
     }
#endif
    for(i=P->n;i<P->m;i++)P->nIndices[i]=0;

    P->mIndices=(int*)realloc((void*)P->mIndices,P->m*sizeof(int));

#ifndef MFNOSAFETYNET
    if(P->mIndices==NULL)
     {
      sprintf(MFPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",P->m*sizeof(int));
      MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return -1;
     }
#endif

    for(i=P->n;i<P->m;i++)P->mIndices[i]=0;

    P->mark=(int*)realloc((void*)P->mark,P->m*sizeof(int));

#ifndef MFNOSAFETYNET
    if(P->mark==NULL)
     {
      sprintf(MFPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",P->m*sizeof(int));
      MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return -1;
     }
#endif

    for(i=P->n;i<P->m;i++)P->mark[i]=0;

    P->indices=(int**)realloc((void*)P->indices,P->m*sizeof(int*));

#ifndef MFNOSAFETYNET
    if(P->indices==NULL)
     {
      sprintf(MFPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",P->m*sizeof(int*));
      MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return -1;
     }
#endif

    for(i=P->n;i<P->m;i++)P->indices[i]=NULL;

#ifdef MFALLOWVERBOSE
    if(verbose)printf(" now space for P->m(%d) vertices.\n");
#endif

   }

#ifdef MFALLOWVERBOSE
  if(verbose)printf(" returning vertex %d.\n",P->n);
#endif

  (P->n)++;
  return P->n-1;
 }

int MFPolytopeAddVertexIndexIntoSet(MFPolytope P,int v,int index, MFErrorHandler e)
 {
  static char RoutineName[]={"MFPolytopeAddVertexIndexIntoSet"};
  int i,j;
  int *tmp;
  int verbose=0;

#ifdef MFNOCONFIDENCE
  if(P==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"First argument is NULL");
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(v<0 || v>=P->n)
   {
    fprintf(stdout,"ERROR: MFPolytopeAddVertexIndexIntoSet - Invalid face number %d (must be in [0,%d) )\n",v,P->n);
    fprintf(stderr,"MFPolytopeAddVertexIndexIntoSet - Invalid face number %d (must be in [0,%d) )\n",v,P->n);
    return;
   }
#endif

  if(P->nIndices[v]>=P->mIndices[v])
   {
    P->mIndices[v]+=MFPolytopeIncrementToAllocateToIndexSet;
    P->indices[v]=(int*)realloc((void*)P->indices[v],P->mIndices[v]*sizeof(int));

#ifndef MFNOSAFETYNET
    if(P->indices[v]==NULL)
     {
      sprintf(MFPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",P->mIndices[v]*sizeof(int));
      MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return -1;
     }
#endif

    for(i=P->nIndices[v];i<P->mIndices[v];i++)(P->indices[v])[i]=-1;
   }

/* Insert [9] into [11,49] is wrong? */

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("Inserting [%d] into [%d",index,(P->indices[v])[0]);
    for(i=1;i<P->nIndices[v];i++)printf(",%d",(P->indices[v])[i]);
    printf("]\n");fflush(stdout);
   }
#endif

  j=0;
  while(j<P->nIndices[v] && (P->indices[v])[j]<index)j++;

#ifdef MFALLOWVERBOSE
  if(verbose)printf("  first index greater than %d is %d, n=%d\n",index,j,P->nIndices[v]);
#endif

  for(i=P->nIndices[v]-1;i>=j;i--)
   {

#ifdef MFALLOWVERBOSE
    if(verbose){printf("(P->indices[v])[%d]=(P->indices[v])[%d]=%d;\n",i+1,i,(P->indices[v])[i]);fflush(stdout);}
#endif

    (P->indices[v])[i+1]=(P->indices[v])[i];
   }

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("  shifted [");
    if(j!=0)printf("%d",(P->indices[v])[0]);
    for(i=1;i<P->nIndices[v]+1;i++){printf(",");if(j!=i)printf("%d",(P->indices[v])[i]);}
    printf("]\n");fflush(stdout);
   }
#endif

  (P->indices[v])[j]=index;
  (P->nIndices[v])++;

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("  final [%d",(P->indices[v])[0]);
    for(i=1;i<P->nIndices[v];i++)printf(",%d",(P->indices[v])[i]);
    printf("]\n");fflush(stdout);
   }
#endif

  return j;
 }

int MFPolytopeIntersectIndexSets(MFPolytope P,int v1,int v2,int *inter, MFErrorHandler e)
 {
  static char RoutineName[]={"MFPolytopeIntersectIndexSets"};

/* Returns the number of entries two index sets have in common  */
/*   Assumes that the index sets are sorted in increasing order */

  int n;
  int n1,n2;
  int *inter1;
  int *inter2;
  int i1,i2;

#ifdef MFNOCONFIDENCE
  if(P==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"First argument is NULL");
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  n1=P->nIndices[v1];
  inter1=P->indices[v1];
  n2=P->nIndices[v2];
  inter2=P->indices[v2];

  n=0;
  i1=0;i2=0;
  while(i1<n1&&i2<n2)
   {
    if(inter1[i1]==inter2[i2])
     {
      inter[n]=inter1[i1];
      n++;
      i1++;i2++;
     }else if(inter1[i1]<inter2[i2])
     {
      i1++;
     }else if(inter1[i1]>inter2[i2])
     {
      i2++;
     }
   }

  return n;
 }

int MFPolytopeNumberOfFaces(MFPolytope P, MFErrorHandler e)
 {
  static char RoutineName[]={"MFPolytopeNumberOfFaces"};

#ifdef MFNOCONFIDENCE
  if(P==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"First argument is NULL");
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  return P->nFaces;
 }

int MFPolytopeFaceIndex(MFPolytope P,int face, MFErrorHandler e)
 {
  static char RoutineName[]={"MFPolytopeFaceIndex"};

#ifdef MFNOCONFIDENCE
  if(P==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"First argument is NULL");
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  return P->face[face];
 }

MFKVector MFPolytopeFaceNormal(MFPolytope P,int face, MFErrorHandler e)
 {
  static char RoutineName[]={"MFPolytopeFaceNormal"};

#ifdef MFNOCONFIDENCE
  if(P==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"First argument is NULL");
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  return P->faceN[face];
 }

double MFPolytopeFaceOrigin(MFPolytope P,int face, MFErrorHandler e)
 {
  static char RoutineName[]={"MFPolytopeFaceOrigin"};

#ifdef MFNOCONFIDENCE
  if(P==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"First argument is NULL");
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  return P->faceO[face];
 }

int MFPolytopeNumberOfVerticesOnFace(MFPolytope P,int face, MFErrorHandler e)
 {
  static char RoutineName[]={"MFPolytopeNumberOfVerticesOnFace"};

#ifdef MFNOCONFIDENCE
  if(P==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"First argument is NULL");
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  return P->nFaceV[face];
 }

double MFPolytopeRadiusOfVertex(MFPolytope P,int vertex, MFErrorHandler e)
 {
  static char RoutineName[]={"MFPolytopeRadiusOfVertex"};
  int i;
  double r;

#ifdef MFNOCONFIDENCE
  if(P==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"First argument is NULL");
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  if(vertex<0 || vertex>=P->n)
   {
    sprintf(MFPolytopeErrorMsg,"Invalid vertex=%d, P->n=%d\n",vertex,P->n);
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    return -1.;
   }

  r=0.;
  for(i=0;i<P->k;i++)
   {
    r+=P->v[i+P->k*vertex]*P->v[i+P->k*vertex];
   }
  r=sqrt(r);
  return r;
 }


void MFPolytopeUpdateFaceList(MFPolytope P, MFErrorHandler e)
 {
  static char RoutineName[]={"MFPolytopeUpdateFaceList"};
  int i,j;
  int index;
  int l,jj;
  int verbose=0;

#ifdef MFNOCONFIDENCE
  if(P==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"First argument is NULL");
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

#ifdef MFALLOWVERBOSE
  if(verbose)printf("MFPolytopeUpdateFaceList 0x%8.8x\n",P);
#endif

#ifdef MFALLOWVERBOSE
  if(verbose)printf(" This polytope has %d vertices and %d faces\n",P->n,P->nFaces);
#endif

  for(i=0;i<P->nFaces;i++)P->nFaceV[i]=0;

  for(i=0;i<P->n;i++)
   {

#ifdef MFALLOWVERBOSE
    if(verbose){printf(" Vertex %d has %d indices\n",i,P->nIndices[i]);fflush(stdout);
                printf("\nP->nIndices[%d]=%d, P->indices[%d]=0x%8.8x\n",i,P->nIndices[i],i,P->indices[i]);fflush(stdout);}
#endif

    for(j=0;j<P->nIndices[i];j++)
     {

#ifdef MFALLOWVERBOSE
      if(verbose){printf("   Index %d is %d\n",j,(P->indices[i])[j]);fflush(stdout);}
#endif

      index=(P->indices[i])[j];
      jj=-1;
      for(l=0;l<P->nFaces;l++)if(P->face[l]==index)jj=l;
      if(jj>-1)
       {
#ifdef MFALLOWVERBOSE
        if(verbose){printf("   face %d has thisP index\n",jj);fflush(stdout);}
#endif

        P->nFaceV[jj]++;
       }
     }
   }

#define MFREMOVEREDUNDANTFACES
#ifdef MFREMOVEREDUNDANTFACES

#ifdef MFALLOWVERBOSE
  if(verbose){printf("Now remove faces with too few points\n");fflush(stdout);}
#endif

  for(i=0;i<P->nFaces;i++)
   {

#ifdef MFALLOWVERBOSE
    if(verbose){printf(" Face %d contains %d vertices\n",i,P->nFaceV[i]);fflush(stdout);
               printf("This Polytope now has %d faces\n",P->nFaces);}
#endif

    if(P->nFaceV[i]<P->k)
     {

#ifdef MFALLOWVERBOSE
      if(verbose){printf(" Removing face %d\n",P->face[i]);fflush(stdout);}
#endif

      for(j=0;j<P->n;j++)
       {

#ifdef MFALLOWVERBOSE
        if(verbose)
         {
          printf(" Vertex %d has %d indices [",j,P->nIndices[j]);
          for(jj=0;jj<P->nIndices[j];jj++)
           {
            if(jj>0)printf(",");
            printf("%d",(P->indices[j])[jj]);
           }
          printf("]\n");
          fflush(stdout);
         }
#endif

        for(jj=0;jj<P->nIndices[j];jj++)
         {
          if((P->indices[j])[jj]==P->face[i])
           {

#ifdef MFALLOWVERBOSE
            if(verbose){printf("   index %d is the face being removed, moving index set down\n",jj);fflush(stdout);}
#endif

            for(l=jj;l<P->nIndices[j]-1;l++)
             {
              (P->indices[j])[l]=(P->indices[j])[l+1];
             }
            P->nIndices[j]--;
            jj--;
           }
         }
#ifdef MFALLOWVERBOSE
        if(verbose)
         {
          printf(" >Vertex %d now has %d indices [",j,P->nIndices[j]);
          for(jj=0;jj<P->nIndices[j];jj++)
           {
            if(jj>0)printf(",");
            printf("%d",(P->indices[j])[jj]);
           }
          printf("]\n");
          fflush(stdout);
         }
#endif

       }
      P->nFaces--;
      P->face[i]=P->face[P->nFaces];
      P->nFaceV[i]=P->nFaceV[P->nFaces];
      if(P->faceN[i]!=NULL)
       {
        MFFreeKVector(P->faceN[i],e);
        P->faceN[i]=NULL;
       }
      P->faceN[i]=P->faceN[P->nFaces];
      P->faceN[P->nFaces]=NULL;
      P->faceO[i]=P->faceO[P->nFaces];
      i--;
     }
   }
#endif

#ifdef MFALLOWVERBOSE
  if(verbose)printf("done MFPolytopeUpdateFaceList\n");
#endif

  return;
 }

int MFPolytopeVerticesOnSameEdge(MFPolytope P,int i,int j, MFErrorHandler e)
 {
  static char RoutineName[]={"MFPolytopeVerticesOnSameEdge"};
  int n1,n2;
/*int *inter1,*inter2;*/
  int *inter;
  int n,m,result;

#ifdef MFNOCONFIDENCE
  if(P==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"First argument is NULL");
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  n1=P->nIndices[i];
  n2=P->nIndices[j];
/*inter1=P->indices[i];
  inter2=P->indices[j];*/
  m=n1;if(n2>m)m=n2;
  inter=(int*)malloc(m*sizeof(int));

#ifndef MFNOSAFETYNET
  if(inter==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",m*sizeof(int));
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return -1;
   }
#endif

  n=MFPolytopeIntersectIndexSets(P,i,j,inter,e);

  if(inter!=NULL)free(inter);
  free(inter);
  result=n>=P->k-1;
  return result;
 }

int MFPolytopeVerticesOnEdgeWithIndexLessThan(MFPolytope P,int i,int j,int edge, MFErrorHandler e)
 {
  static char RoutineName[]={"MFPolytopeVerticesOnEdgeWithIndexLessThan"};
  int n1,n2;
  int *inter;
  int n,m,result;

  n1=P->nIndices[i];
  n2=P->nIndices[j];
  m=n1;if(n2>m)m=n2;
  inter=(int*)malloc(m*sizeof(int));

#ifndef MFNOSAFETYNET
  if(inter==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",m*sizeof(int));
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return -1;
   }
#endif

  n=MFPolytopeIntersectIndexSets(P,i,j,inter,e);

  result=1;
  for(m=0;m<n;m++)if(inter[m]>=edge)result=0;

  if(inter!=NULL)free(inter);
  return result;
 }

int MFPolytopeVertexLargestIndex(MFPolytope P,int i, MFErrorHandler e)
 {
  static char RoutineName[]={"MFPolytopeVertexLargestIndex"};
  int j;

#ifdef MFNOCONFIDENCE
  if(P==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"First argument is NULL");
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(i<0 || i>=P->n)
   {
    sprintf(MFPolytopeErrorMsg,"Invalid vertex number %d (must be in [0,%d) )",i,P->n);
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(P->indices==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"Polytope has NULL vertex index list");
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  return (P->indices[i])[P->nIndices[i]-1];
 }

int MFPolytopeVertexSmallestIndex(MFPolytope P,int i, MFErrorHandler e)
 {
  static char RoutineName[]={"MFPolytopeVertexSmallestIndex"};
  int j;

#ifdef MFNOCONFIDENCE
  if(P==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"First argument is NULL");
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(i<0 || i>=P->n)
   {
    sprintf(MFPolytopeErrorMsg,"Invalid vertex number %d (must be in [0,%d) )",i,P->n);
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(P->indices==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"Polytope has NULL vertex index list");
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  return (P->indices[i])[0];
 }

int MFPolytopeLargestVertexIndex(MFPolytope P, MFErrorHandler e)
 {
  static char RoutineName[]={"MFPolytopeLargestVertexIndex"};
  int i,n;
  int t;
  int result;

#ifdef MFNOCONFIDENCE
  if(P==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"First argument is NULL");
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  n=MFPolytopeNumberOfVertices(P,e);
  if(n>0)result=MFPolytopeVertexLargestIndex(P,0,e);
   else  result=-1;

  for(i=1;i<n;i++)
   {
    t=MFPolytopeVertexLargestIndex(P,i,e);
    if(t>result)result=t;
   }

  return result;
 }

int MFPolytopeSmallestVertexIndex(MFPolytope P, MFErrorHandler e)
 {
  static char RoutineName[]={"MFPolytopeSmallestVertexIndex"};
  int i,n;
  int t;
  int result;

#ifdef MFNOCONFIDENCE
  if(P==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"First argument is NULL");
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  n=MFPolytopeNumberOfVertices(P,e);
  if(n>0)result=MFPolytopeVertexSmallestIndex(P,0,e);
   else  result=-1;

  for(i=1;i<n;i++)
   {
    t=MFPolytopeVertexSmallestIndex(P,i,e);
    if(t<result)result=t;
   }

  return result;
 }

int MFPolytopeInterior(MFPolytope P, MFKVector s, MFErrorHandler e)
 {
  static char RoutineName[]={"MFPolytopeInterior"};
  int face;

#ifdef MFNOCONFIDENCE
  if(P==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"First argument is NULL");
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    return -1;
   }

  if(s==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"Second argument is NULL");
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    return -1;
   }
#endif

  if(P->n==0)return 0;
  for(face=0;face<P->nFaces;face++)
   {
    if(MFKVDot(s,P->faceN[face],e)-P->faceO[face]>0.)
     {
      return 0;
     }
   }
  return 1;
 }

int MFPolytopeTotallyInterior(MFPolytope P, MFKVector s, double eps, MFErrorHandler e)
 {
  static char RoutineName[]={"MFPolytopeTotallyInterior"};
  int face;

#ifdef MFNOCONFIDENCE
  if(P==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"First argument is NULL");
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(s==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"Second argument is NULL");
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  if(P->n==0)
   {
    return 0;
   }
  for(face=0;face<P->nFaces;face++)
   {
     if(MFKVDot(s,P->faceN[face],e)-P->faceO[face]>-eps)
      {
       return 0;
      }
    }
  return 1;
 }

double MFPolytopeLargestRadiusOfVertex(MFPolytope P, MFErrorHandler e)
 {
  static char RoutineName[]={"MFPolytopeLargestRadiusOfVertex"};
  double r,result;
  int i;

#ifdef MFNOCONFIDENCE
  if(P==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"First argument is NULL");
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  result=0;
  for(i=0;i<P->n;i++)
   if((r=MFPolytopeRadiusOfVertex(P,i,e))>result)result=r;

  return result;
 }

void MFWritePolytope(FILE *fid,MFPolytope P, MFErrorHandler e)
 {
  static char RoutineName[]={"MFWritePolytope"};
  int i,j;

  fprintf(fid,"%s\n","Polytope");
  if(P==NULL)
   {
    fprintf(fid,"-1\n");
    return;
   }
  fprintf(fid,"%d\n",P->k);

  fprintf(fid,"%d %d\n",P->n,P->m);
  for(i=0;i<P->n*P->k;i++)
   {
    if(i>0)fprintf(fid," ");
    fprintf(fid,"%le",P->v[i]);
   }
  fprintf(fid,"\n");
  for(i=0;i<P->n;i++)
   {
    if(i>0)fprintf(fid," ");
    fprintf(fid,"%d",P->nIndices[i]);
   }
  fprintf(fid,"\n");
  for(i=0;i<P->n;i++)
   {
    if(i>0)fprintf(fid," ");
    fprintf(fid,"%d",P->mIndices[i]);
   }
  fprintf(fid,"\n");
  for(i=0;i<P->n;i++)
   {
    for(j=0;j<P->nIndices[i];j++)
     {
      if(j>0)fprintf(fid," ");
      fprintf(fid,"%d",(P->indices[i])[j]);
     }
    fprintf(fid,"\n");
   }

  fprintf(fid,"%d %d\n",P->nFaces,P->mFaces);
  for(i=0;i<P->nFaces;i++)
   {
    if(i>0)fprintf(fid," ");
    fprintf(fid,"%d",P->face[i]);
   }
  fprintf(fid,"\n");
  for(i=0;i<P->nFaces;i++)
   {
    if(i>0)fprintf(fid," ");
    fprintf(fid,"%d",P->nFaceV[i]);
   }
  fprintf(fid,"\n");
  for(i=0;i<P->nFaces;i++)
   {
    MFWriteKVector(fid,P->faceN[i],e);
   }
  fprintf(fid,"\n");
  for(i=0;i<P->nFaces;i++)
   {
    if(i>0)fprintf(fid," ");
    fprintf(fid,"%le",P->faceO[i]);
   }
  fprintf(fid,"\n");

  return;
 }

MFPolytope MFReadPolytope(FILE *fid, MFErrorHandler e)
 {
  static char RoutineName[]={"MFReadPolytope"};
  int i,j;
  MFPolytope P;
  char tag[100]="";
  int k=0;

  fscanf(fid,"%s\n",tag);

#ifdef MFNOCONFIDENCE
  if(strcmp(tag,"Polytope"))
   {
    sprintf(MFPolytopeErrorMsg,"Next Object is not a Polytope! (%s)\n",RoutineName,tag);
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  fscanf(fid,"%d\n",&k);
  if(k<0)
   {
    return NULL;
   }

  P=(MFPolytope)malloc(sizeof(struct MFPolytopeSt));

#ifndef MFNOSAFETYNET
  if(P==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct MFPolytopeSt));
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  P->k=k;

  fscanf(fid,"%d %d\n",&(P->n),&(P->m));

  if(P->n>0)
   {
    P->v=(double*)malloc(P->n*P->k*sizeof(double));

#ifndef MFNOSAFETYNET
    if(P->v==NULL)
     {
      sprintf(MFPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",P->n*P->k*sizeof(double));
      MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return NULL;
     }
#endif
    for(i=0;i<P->n*P->k;i++)
     {
      if(i>0)fscanf(fid," ");
      fscanf(fid,"%le",&(P->v[i]));
     }
    fscanf(fid,"\n");

    P->nIndices=(int*)malloc(P->n*sizeof(int));

#ifndef MFNOSAFETYNET
    if(P->nIndices==NULL)
     {
      sprintf(MFPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",P->n*sizeof(int));
      MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return NULL;
     }
#endif
    for(i=0;i<P->n;i++)
     {
      if(i>0)fscanf(fid," ");
      fscanf(fid,"%d",&(P->nIndices[i]));
     }
    fscanf(fid,"\n");

    P->mIndices=(int*)malloc(P->n*sizeof(int));

#ifndef MFNOSAFETYNET
    if(P->mIndices==NULL)
     {
      sprintf(MFPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",P->n*sizeof(int));
      MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return NULL;
     }
#endif
    for(i=0;i<P->n;i++)
     {
      if(i>0)fscanf(fid," ");
      fscanf(fid,"%d",&(P->mIndices[i]));
     }
    fscanf(fid,"\n");

    P->indices=(int**)malloc(P->n*sizeof(int*));

#ifndef MFNOSAFETYNET
    if(P->indices==NULL)
     {
      sprintf(MFPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",P->n*sizeof(int*));
      MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return NULL;
     }
#endif
  
    for(i=0;i<P->n;i++)
     {
      P->indices[i]=(int*)malloc(P->mIndices[i]*sizeof(int));

#ifndef MFNOSAFETYNET
      if(P->indices[i]==NULL)
       {
        sprintf(MFPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",P->mIndices[i]*sizeof(int));
        MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
        MFErrorHandlerOutOfMemory(e);
        return NULL;
       }
#endif
      for(j=0;j<P->nIndices[i];j++)
       {
        if(j>0)fscanf(fid," ");
        fscanf(fid,"%d",&((P->indices[i])[j]));
       }
      fscanf(fid,"\n");
     }
   }else{
    P->v=NULL;
    P->nIndices=NULL;
    P->mIndices=NULL;
    P->indices=NULL;
   }

  fscanf(fid,"%d %d\n",&(P->nFaces),&(P->mFaces));

  if(P->nFaces>0)
   {
    P->face=(int*)malloc(P->mFaces*sizeof(int));

#ifndef MFNOSAFETYNET
    if(P->face==NULL)
     {
      sprintf(MFPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",P->mFaces*sizeof(int));
      MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return NULL;
     }
#endif
    for(i=0;i<P->nFaces;i++)
     {
      if(i>0)fscanf(fid," ");
      fscanf(fid,"%d",&(P->face[i]));
     }
    fscanf(fid,"\n");

    P->nFaceV=(int*)malloc(P->mFaces*sizeof(int));

#ifndef MFNOSAFETYNET
    if(P->nFaceV==NULL)
     {
      sprintf(MFPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",P->mFaces*sizeof(int));
      MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return NULL;
     }
#endif
    for(i=0;i<P->nFaces;i++)
     {
      if(i>0)fscanf(fid," ");
      fscanf(fid,"%d",&(P->nFaceV[i]));
     }
    fscanf(fid,"\n");

    P->faceN=(MFKVector*)malloc(P->mFaces*sizeof(MFKVector));

#ifndef MFNOSAFETYNET
    if(P->faceN==NULL)
     {
      sprintf(MFPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",P->mFaces*sizeof(MFKVector));
      MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return NULL;
     }
#endif
    for(i=0;i<P->mFaces;i++)P->faceN[i]=NULL;
    for(i=0;i<P->nFaces;i++)
     {
      P->faceN[i]=MFReadKVector(fid,e);
     }
    fscanf(fid,"\n");

    P->faceO=(double*)malloc(P->mFaces*sizeof(double));

#ifndef MFNOSAFETYNET
    if(P->faceO==NULL)
     {
      sprintf(MFPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",P->mFaces*sizeof(double));
      MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return NULL;
     }
#endif
    for(i=0;i<P->nFaces;i++)
     {
      if(i>0)fscanf(fid," ");
      fscanf(fid,"%le",&(P->faceO[i]));
     }
    fscanf(fid,"\n");
   }else{
    P->face=NULL;
    P->nFaceV=NULL;
    P->faceN=NULL;
    P->faceO=NULL;
   }

  return P;
 }

int MFTestPolytope(MFPolytope P, MFErrorHandler e)
 {
  static char RoutineName[]={"MFTestPolytope"};
  int coordinate;
  int vertex;
  int face;
  MFKVector s;
  int good;

#ifdef MFNOCONFIDENCE
  if(P==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"First argument is NULL");
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  s=MFCreateKVector(P->k,e);

  good=1;
  for(vertex=0;vertex<P->n;vertex++)
   {
    for(coordinate=0;coordinate<P->k;coordinate++)MFKVSetC(s,coordinate,P->v[coordinate+vertex*P->k],e);
    for(face=0;face<P->nFaces;face++)
     {
      if(MFKVDot(s,P->faceN[face],e)>P->faceO[face]+1.e-5)
       {
        printf("ERROR, Here's a polytope vertex (%d) that lies completely outside the planes (%d)!\n",vertex,face);
        printf("       ");
        MFPrintKVector(stdout,s,e);
        printf(".");
        MFPrintKVector(stdout,P->faceN[face],e);
        printf("= %le > %le\n",MFKVDot(s,P->faceN[face],e),P->faceO[face]);
        good=0;
       }
     }
   }

  if(!good)
   {
    MFPrintPolytope(stdout,P,e);fflush(stdout);
    abort();
   }
  MFFreeKVector(s,e);

  return 1;
 }

void MFPolytopeMergeCloseVertices(MFPolytope P,double eps, MFErrorHandler e)
 {
  static char RoutineName[]={"MFPolytopeMergeCloseVertices"};
  int i,j,l;
  double dist;
  int nElim;
  int verbose=0;

#ifdef MFNOCONFIDENCE
  if(P==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"First argument is NULL");
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

#ifdef MFALLOWVERBOSE
  if(verbose)printf("MFPolytopeMergeCloseVertices\n");
#endif

  nElim=0;
  loop:
  for(i=0;i<P->n;i++)
   {
    for(j=0;j<P->n;j++)
     {
      if(i!=j)
       {
        dist=0.;
        for(l=0;l<P->k;l++)
         dist+=fabs(P->v[l+i*P->k]-P->v[l+j*P->k]);
        if(dist/P->k<eps)
         {
          nElim++;

#ifdef MFALLOWVERBOSE
          if(verbose)
           {
            printf("MFPolytopeMergeCloseVertices: Merging vertices %d and %d\n",i,j);
            printf("  before:");
            MFPrintPolytope(stdout,P,e);
           }
#endif

          MFPolytopeMergeVertices(P,i,j,e);

#ifdef MFALLOWVERBOSE
          if(verbose)
           {
            printf("  after:");
            MFPrintPolytope(stdout,P,e);
           }
#endif

          goto loop;
         }
       }
     }
   }

#ifdef MFALLOWVERBOSE
  if(verbose&&nElim>0)printf("MFPolytopeMergeCloseVertices: Merged %d vertices\n",nElim);
#endif

  return;
 }

void MFPolytopeMergeVertices(MFPolytope P,int i,int j, MFErrorHandler e)
 {
  static char RoutineName[]={"MFPolytopeMergeVertices"};
  int n1,n2,n;
  int l;
  int *Union;

#ifdef MFNOCONFIDENCE
  if(P==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"First argument is NULL");
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  n1=P->nIndices[i];
  n2=P->nIndices[j];
  n=n1+n2;
  Union=(int*)malloc(n*sizeof(int));

#ifndef MFNOSAFETYNET
  if(Union==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(int));
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  n=MFPolytopeUnionOfIndexSets(P,i,j,Union,e);
  P->nIndices[i]=n;
  printf("   free P(0x%8.8x)->indices[%d]=0x%8.8x\n",P,i,P->indices[i]);fflush(stdout);
  free(P->indices[i]);
  P->indices[i]=Union;
  if(j<P->n-1)
   {
    printf("%s replace P(0x%8.8x)->indices[%d]=0x%8.8x with P(0x%8.8x)->indices[%d]=0x%8.8x\n",RoutineName,P,j,P->indices[j],P,P->n-1,P->indices[P->n-1]);fflush(stdout);
    printf("   free P(0x%8.8x)->indices[%d]=0x%8.8x\n",P,j,P->indices[j]);fflush(stdout);
    free(P->indices[j]);
    P->nIndices[j]=P->nIndices[P->n-1];
    P->indices[j]=P->indices[P->n-1];
    for(l=0;l<P->k;l++)
      P->v[l+P->k*j]=P->v[l+P->k*(P->n-1)];
   }
  (P->n)--;

/*if(!MFTestVertexIndexOrder(P->nIndices[i],P->indices[i]))
   {
    printf("Error: MFPolytopeMergeVertices at exit vertex indices not sorted\n");
    fflush(stdout);
   }*/

  return;
 }

int MFPolytopeUnionOfIndexSets(MFPolytope P,int v1,int v2,int *inter, MFErrorHandler e)
 {
  static char RoutineName[]={"MFPolytopeUnionOfIndexSets"};

/* Returns the number of entries two index sets have in common  */
/*   Assumes that the index sets are sorted in increasing order */

  int n;
  int n1,n2;
  int *inter1;
  int *inter2;
  int i1,i2;
  int i;
  int verbose=0;

#ifdef MFNOCONFIDENCE
  if(P==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"First argument is NULL");
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  n1=P->nIndices[v1];
  inter1=P->indices[v1];
  n2=P->nIndices[v2];
  inter2=P->indices[v2];

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("%s\n",RoutineName);
    printf("set 1: [");
    for(i=0;i<n1;i++){if(i>0)printf(",");printf("%d",inter1[i]);}
    printf("]\n");
    printf("set 2: [");
    for(i=0;i<n2;i++){if(i>0)printf(",");printf("%d",inter2[i]);}
    printf("]\n");
   }
#endif

  n=0;
  i1=0;i2=0;
  while(i1<n1||i2<n2)
   {
    if(i1>=n1)
     {

#ifdef MFALLOWVERBOSE
      if(verbose)printf("take 2 %d (out of 1)\n",inter2[i2]);
#endif

      inter[n]=inter2[i2];
      n++;
      i2++;
    }else if(i2>=n2)
     {

#ifdef MFALLOWVERBOSE
      if(verbose)printf("take 1 %d (out of 2)\n",inter1[i2]);
#endif

      inter[n]=inter1[i1];
      n++;
      i1++;
     }else if(inter1[i1]<inter2[i2])
     {

#ifdef MFALLOWVERBOSE
      if(verbose)printf("take 1 %d (smaller than 2 %d)\n",inter1[i1],inter2[i2]);
#endif

      inter[n]=inter1[i1];
      n++;
      i1++;
     }else if(inter2[i2]<inter1[i1])
     {

#ifdef MFALLOWVERBOSE
      if(verbose)printf("take 2 %d (smaller than 1 %d)\n",inter2[i1],inter1[i2]);
#endif

      inter[n]=inter2[i2];
      n++;
      i2++;
     }else{

#ifdef MFALLOWVERBOSE
      if(verbose)printf("take 2 %d (same as 1 %d)\n",inter2[i1],inter1[i2]);
#endif

      inter[n]=inter2[i2];
      n++;
      i1++;
      i2++;
     }
   }

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("union: [");
    for(i=0;i<n;i++){if(i>0)printf(",");printf("%d",inter[i]);}
    printf("]\n");fflush(stdout);
   }
#endif

  return n;
 }

int MFTestVertexIndexOrder(int n,int *s, MFErrorHandler e)
 {
  static char RoutineName[]={"MFTestVertexIndexOrder"};
  int i;

  for(i=1;i<n;i++)
   if(s[i-1]>=s[i])
    {
     return 0;
    }
  return 1;
 }

int MFPolytopeTestVertexList(MFPolytope P, MFErrorHandler e)
 {
  static char RoutineName[]={"MFPolytopeTestVertexList"};
  int i,j,k,n,m;
  int same;
  int good;

#ifdef MFNOCONFIDENCE
  if(P==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"First argument is NULL");
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif
  
  n=P->n;
  good=1;
  for(i=0;i<n-1;i++)
   {
    for(j=i+1;j<n;j++)
     {
      if(P->nIndices[i]==P->nIndices[j])
       {
        m=P->nIndices[i];
        same=1;
        for(k=0;k<m;k++)
         {
          if((P->indices[i])[k]!=(P->indices[j])[k])same=0;
         }
        if(same)
         {
          printf("Problem with Polytope 0x%8.8x, two vertices %d and %d with the same index set (",P,i,j);
          for(k=0;k<m;k++)
           {
            if(k>0)printf(",");
            printf("%d",(P->indices[i])[k]);
           }
          printf("),(");
          for(k=0;k<m;k++)
           {
            if(k>0)printf(",");
            printf("%d",(P->indices[j])[k]);
           }
          printf(")\n");
          fflush(stdout);
          good=0;
         }
       }
     }
   }
  
  for(i=0;i<n;i++)
   {
    m=P->nIndices[i];
    for(j=0;j<m-1;j++)
     {
      for(k=j+1;k<m;k++)
       {
        if((P->indices[i])[j]==(P->indices[i])[k])
         {
          printf("Problem with Polytope 0x%8.8x, vertex %d has a repeated index (",P,i);
          for(k=0;k<m;k++)
           {
            if(k>0)printf(",");
            printf("%d",(P->indices[i])[k]);
           }
          printf(")\n");
         }
       }
     }
   }
  if(!good)
   {
    MFPrintPolytope(stdout,P,e);fflush(stdout);
   }

  return good;
 }

#define READY
#ifdef READY
MFPolytope MFCreateSectionOfPolytope(MFPolytope P, int index, MFKVector nrm, double on, MFErrorHandler e)
 {
  static char RoutineName[]={"MFCreateSectionOfPolytope"};
  MFPolytope result;
  int i,j;
  int k;
  int n,m;
  double *d;
  double t;
  int *inter;

  result=(MFPolytope)malloc(sizeof(struct MFPolytopeSt));
  inter=(int*)malloc(P->k*sizeof(int));

  k=P->k;

  result->k=P->k;
  result->R=P->R;

  d=(double*)malloc(P->n*sizeof(double));

  t=fabs(on);
  for(j=0;j<P->k;j++)if(fabs(MFKV_C(nrm,j,e))>t)t=fabs(MFKV_C(nrm,j,e));

/* Mark the vertices of P. */

  for(i=0;i<P->n;i++)
   {
    d[i]=-on;
    for(j=0;j<P->k;j++)d[i]+=P->v[j+P->k*i]*MFKV_C(nrm,j,e);
    d[i]=d[i]/t;
   }

/* Count the number of edges that cross (the number of vertices in the result). */

  result->n=0;
  for(i=0;i<P->n;i++)
   {
    for(j=0;j<i;j++)
     {
      if(MFPolytopeIntersectIndexSets(P,i,j,inter,e)==P->k-1 && d[i]*d[j]<=0.)result->n++;
     }
   }

/* Allocate space for the new vertices. */

  result->m=result->n;
  result->v=(double*)malloc(result->n*result->k*sizeof(double));
  result->nIndices=(int*)malloc(result->n*sizeof(int));
  result->mIndices=(int*)malloc(result->n*sizeof(int));
  result->indices=(int**)malloc(result->n*sizeof(int*));

  for(i=0;i<result->n;i++)
   {
    result->nIndices[i]=result->k;
    result->mIndices[i]=result->k;
    result->indices[i]=(int*)malloc(result->k*sizeof(int));
   }

/* Fill in the values. */
 
  n=0;
  for(i=0;i<P->n;i++)
   {
    for(j=0;j<i;j++)
     {
      if(MFPolytopeIntersectIndexSets(P,i,j,inter,e)==P->k-1 && d[i]*d[j]<=0.)
       {
        t=d[i]/(d[i]-d[j]);
        m=MFPolytopeIntersectIndexSets(P,i,j,result->indices[n],e);

        result->indices[n][result->k-1]=index;

        for(m=0;m<result->n;m++)result->v[m+n*k]=(1-t)*P->v[m+i*k]+t*P->v[m+j*k];
        n++;
       }
     }
   }

  if(result->n>0)
   {
    result->nFaces=P->nFaces+1;
    result->mFaces=result->nFaces;
    result->face=(int*)malloc(result->mFaces*sizeof(int));
    for(i=0;i<P->nFaces;i++)result->face[i]=P->face[i];
    result->face[P->nFaces]=index;
   }else{
    result->nFaces=0;
    result->mFaces=0;
    result->face=(int*)NULL;
   }

  free(d);
  free(inter);

  return result;
 }
#endif

int MFPolytopeClosestFace(MFPolytope P, MFKVector s, double *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFPolytopeClosestFace"};
  int face;
  int result;
  double t;

#ifdef MFNOCONFIDENCE
  if(P==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"First argument is NULL");
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    return -1;
   }
  
  if(s==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"Second argument is NULL");
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    return -1;
   }
#endif

  *d=0.;
  if(P->n==0)return -1;
  result=-1;
  for(face=0;face<P->nFaces;face++)
   {
    t=MFKVDot(s,P->faceN[face],e)-P->faceO[face];
    if(face==0||t>*d)
     {
      *d=t;
      result=P->face[face];
     }
   }
  return result;
 }

int MFPolytopeGetVertexMark(MFPolytope P,int i, MFErrorHandler e)
 {
  static char RoutineName[]={"MFPolytopeGetVertexMark"};

#ifdef MFNOCONFIDENCE
  if(P==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"First argument is NULL");
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    return -1;
   }
  if(i<0 || i>=P->n)
   {
    sprintf(MFPolytopeErrorMsg,"Invalid vertex number %d (must be in [0,%d) )",i,P->n);
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    return -1;
   }
#endif

  return P->mark[i];
 }

void MFPolytopeSetVertexMark(MFPolytope P,int i, int mark, MFErrorHandler e)
 {
  static char RoutineName[]={"MFPolytopeSetVertexMark"};

#ifdef MFNOCONFIDENCE
  if(P==NULL)
   {
    sprintf(MFPolytopeErrorMsg,"First argument is NULL");
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    return;
   }
  if(i<0 || i>=P->n)
   {
    sprintf(MFPolytopeErrorMsg,"Invalid vertex number %d (must be in [0,%d) )",i,P->n);
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  P->mark[i]=mark;
 }

void MFRecomputePolytopeFromFaces(MFPolytope P, MFErrorHandler e)
 {
  static char RoutineName[]={"MFRecomputePolytopeFromFaces"};
  int i,j,k,l,m;
  int doit,newV;
  double det,vx,vy,vz;
  double a[9],b[3];
  int verbose=0;

  if(P->k!=3)abort();

  P->n=0;
  P->R=1.e100;
  if(P->k==3)
   {
    for(i=0;i<P->nFaces;i++)P->nFaceV[i]=0;

    for(i=0;i<P->nFaces-2;i++)
     {
      for(j=i+1;j<P->nFaces-1;j++)
       {
        for(k=j+1;k<P->nFaces;k++)
         {
          a[0]=MFKV_C(P->faceN[i],0,e); a[3]=MFKV_C(P->faceN[i],1,e); a[6]=MFKV_C(P->faceN[i],2,e);b[0]=P->faceO[i];
          a[1]=MFKV_C(P->faceN[j],0,e); a[4]=MFKV_C(P->faceN[j],1,e); a[7]=MFKV_C(P->faceN[j],2,e);b[1]=P->faceO[j];
          a[2]=MFKV_C(P->faceN[k],0,e); a[5]=MFKV_C(P->faceN[k],1,e); a[8]=MFKV_C(P->faceN[k],2,e);b[2]=P->faceO[k];
          MFSolveFull(3,a,b,e);
          vx=b[0]; vy=b[1]; vz=b[2];
          doit=1;

#ifdef MFALLOWVERBOSE
          if(verbose){printf("Vertex defined by faces %d, %d and %d has coordinates (%lf,%lf,%lf). Det=%le\n",i,j,k,vx,vy,vz,det);}
#endif

          m=0;
          for(l=0;l<P->nFaces;l++)
           {
            if(MFKV_C(P->faceN[l],0,e)*vx+MFKV_C(P->faceN[l],1,e)*vy+MFKV_C(P->faceN[l],2,e)*vz-P->faceO[l]>1.e-7)
             {
              doit=0;

#ifdef MFALLOWVERBOSE
              if(verbose){printf("lies right of face %d,  %le (skipping)\n",l,MFKV_C(P->faceN[l],0,e)*vx+MFKV_C(P->faceN[l],1,e)*vy+MFKV_C(P->faceN[l],2,e)*vz-P->faceO[l]);}
#endif

             }else if(fabs(MFKV_C(P->faceN[l],0,e)*vx+MFKV_C(P->faceN[l],1,e)*vy+MFKV_C(P->faceN[l],2,e)*vz-P->faceO[l])<1.e-7)
             {

#ifdef MFALLOWVERBOSE
              if(verbose){printf("lies on       face %d,  %le\n",l,MFKV_C(P->faceN[l],0,e)*vx+MFKV_C(P->faceN[l],1,e)*vy+MFKV_C(P->faceN[l],2,e)*vz-P->faceO[l]);}
#endif

              if(l!=i && l!=j && l!=k && ( l<i || l<j || l<k))
               {
                doit=0;

#ifdef MFALLOWVERBOSE
                if(verbose){printf("  already been done! (skipping)\n");}
#endif

               }
              m++;
             }else{

#ifdef MFALLOWVERBOSE
              if(verbose){printf("lies left of  face %d, %le\n",l,MFKV_C(P->faceN[l],0,e)*vx+MFKV_C(P->faceN[l],1,e)*vy+MFKV_C(P->faceN[l],2,e)*vz-P->faceO[l]);}
#endif

             }
           }
          if(doit)
           {
            newV=MFPolytopeAddVertex(P,e);
            P->v[0+3*newV]=vx;
            P->v[1+3*newV]=vy;
            P->v[2+3*newV]=vz;
            P->nIndices[newV]=m;
            P->mIndices[newV]=m;
            P->mark[newV]=0;
            P->indices[newV]=(int*)malloc(m*sizeof(int));

            m=0;
            for(l=0;l<P->nFaces;l++)
             {
              if(fabs(MFKV_C(P->faceN[l],0,e)*vx+MFKV_C(P->faceN[l],1,e)*vy+MFKV_C(P->faceN[l],2,e)*vz-P->faceO[l])<1.e-7)
               {
                (P->indices[newV])[m]=P->face[l];m++;
                P->nFaceV[l]++;
               }
             }
           }
         }
       }
     }
   }

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("Recomputed:\n");
    MFPrintPolytope(stdout,P,e);
   }
#endif

#ifdef MFNOCONFIDENCE
  if(!MFPolytopeTestVertexList(P,e))
   {
    sprintf(MFPolytopeErrorMsg,"Polytope vertex list is bad.");
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  return;
 }

void MFClipPolytope(MFPolytope P,int index, double *din, int mark, MFErrorHandler e)
 {
  static char RoutineName[]={"MFClipPolytope"};
  int i,j,k,l;
  int n;
  int flag;
  int verbose=0;
  int *inter;
  int n1,n2;
  double *vi,*vj,*v;
  double t;
  int nold,m,newV;
  int ln,li,lj;
  double *d,dmax;
  int allIn;

  allIn=1;
  for(i=0;i<P->n;i++)if(din[i]>0.)allIn=0;

  if(allIn)return;

  if(0){printf("before:\n");MFPrintPolytope(stdout,P,e);fflush(stdout);}

  d=(double*)malloc((P->m)*sizeof(double));

#ifndef MFNOSAFETYNET
  if(d==NULL&&P->m>0)
   {
    sprintf(MFPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",P->m*sizeof(double));
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  dmax=fabs(din[0]);for(i=1;i<P->n;i++)if(fabs(din[i])>dmax)dmax=fabs(din[i]);

  for(i=0;i<P->n;i++)d[i]=din[i]/dmax;
  if(0){for(i=0;i<P->n;i++)printf("  vertex %d is side %lf\n",i,d[i]);fflush(stdout);}

/* Add vertices on edges that cross the surface */

  nold=P->n;
  for(i=0;i<nold;i++)
   {
    n1=P->nIndices[i];
    for(j=i+1;j<nold;j++)
     {

/* This number should be relative to the longest edge of the polytope */

      if(fabs(d[i])>1.e-14 && fabs(d[i])>1.e-14 &&
            (d[i]<0&&d[j]>0 || d[i]>0 && d[j]<0))
       {
        n2=P->nIndices[j];
        m=n1;if(n2>m)m=n2;

        if(m>0)
         {
          inter=(int*)malloc(m*sizeof(int));

#ifndef MFNOSAFETYNET
          if(inter==NULL)
           {
            sprintf(MFPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",m*sizeof(int));
            MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
            MFErrorHandlerOutOfMemory(e);
            return;
           }
#endif

          n=MFPolytopeIntersectIndexSets(P,i,j,inter,e);
         }else{
          inter=NULL;
          n=0;
         }

        if(0){printf("Vertices %d (%lf) and %d (%lf) have %d indices in common\n",i,d[i],j,d[j],n);fflush(stdout);}
        if(n<(P->k)-1)
         {
          if(inter!=NULL)free(inter);
          inter=NULL;
         }else{
          newV=MFPolytopeAddVertex(P,e);
          t=-d[i]/(d[j]-d[i]);
          if(0){printf("  add vertex, t=%lf\n",t);fflush(stdout);}
          li=P->k*i;
          lj=P->k*j;
          ln=P->k*newV;
          for(l=0;l<P->k;l++)
            P->v[l+ln]=P->v[l+li]+t*(P->v[l+lj]-P->v[l+li]);

          P->nIndices[newV]=n;
          P->mIndices[newV]=m;
          P->mark[newV]=mark;
          P->indices[newV]=inter;
          d=(double*)realloc((void*)d,(P->m)*sizeof(double));

#ifndef MFNOSAFETYNET
          if(d==NULL)
           {
            sprintf(MFPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",P->m*sizeof(double));
            MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
            MFErrorHandlerOutOfMemory(e);
            return;
           }
#endif

          d[newV]=-100.;
          MFPolytopeAddVertexIndexIntoSet(P,newV,index,e);
         }
       }
     }
   }

  if(0){
  printf("after adding vertices:\n");MFPrintPolytope(stdout,P,e);fflush(stdout);
  for(i=0;i<P->n;i++)printf("  vertex %d is side %lf\n",i,d[i]);fflush(stdout);}

/* remove vertices that are on the wrong side */

  for(i=0;i<P->n;i++)
   {
    if(d[i]>0)
     {
      if(i<P->n-1)
       {
        for(j=0;j<P->k;j++)P->v[j+P->k*i]=P->v[j+P->k*(P->n-1)];
        if(P->indices[i]!=NULL)free(P->indices[i]);
        P->indices[i]=P->indices[P->n-1];
        P->nIndices[i]=P->nIndices[P->n-1];
        P->mIndices[i]=P->mIndices[P->n-1];
        P->mark[i]=P->mark[P->n-1];
        d[i]=d[P->n-1];
       }else{
        if(P->indices[i]!=NULL)free(P->indices[i]);
       }
      P->indices[P->n-1]=NULL;
      P->nIndices[P->n-1]=0;
      P->mIndices[P->n-1]=0;
      P->mark[P->n-1]=0;
      d[P->n-1]=0.;
      (P->n)--;
      i--;
     }else if(d[i]==0) MFPolytopeAddVertexIndexIntoSet(P,i,index,e);
   }

  free(d);

  if(0){printf("after adding and removing vertices:\n");MFPrintPolytope(stdout,P,e);fflush(stdout);}

/* Now remove redundant half spaces */
/*   For each appearing index, count occurances. */
/*    If it appears fewer than k-1 times, remove it */
/*    else If dim(span(vertices))<k-1 times, remove it */

  if(P->nFaces>=P->mFaces)
   {
    P->mFaces+=MFPolytopeIncrementToAllocateToFaceList;
    P->face=(int*)realloc((void*)P->face,P->mFaces*sizeof(int));

#ifndef MFNOSAFETYNET
    if(P->face==NULL)
     {
      sprintf(MFPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",P->mFaces*sizeof(int));
      MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

    P->nFaceV=(int*)realloc((void*)P->nFaceV,P->mFaces*sizeof(int));

#ifndef MFNOSAFETYNET
    if(P->nFaceV==NULL)
     {
      sprintf(MFPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",P->mFaces*sizeof(int));
      MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

    P->faceN=(MFKVector*)realloc((void*)(P->faceN),P->mFaces*sizeof(MFKVector));

#ifndef MFNOSAFETYNET
    if(P->faceN==NULL)
     {
      sprintf(MFPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",P->mFaces*sizeof(MFKVector));
      MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

    for(i=P->nFaces;i<P->mFaces;i++)P->faceN[i]=NULL;
  
    P->faceO=(double*)realloc((void*)P->faceO,P->mFaces*sizeof(double));

#ifndef MFNOSAFETYNET
    if(P->faceO==NULL)
     {
      sprintf(MFPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",P->mFaces*sizeof(double));
      MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

   }

/* This must be to guard against something */
  P->face[P->nFaces]=index;
  if(P->faceN[P->nFaces]!=NULL)
   {
    MFFreeKVector(P->faceN[P->nFaces],e);
    P->faceN[P->nFaces]=NULL;
   }
  P->faceO[P->nFaces]=0.;
  P->nFaceV[P->nFaces]=0;
  P->nFaces++;

#ifdef MFALLOWVERBOSE
  if(verbose){MFPrintPolytope(stdout,P,e);
              printf("\n%s   -- call UpdateFaceList, P has %d vertices \n",RoutineName,P->n);}
#endif

#ifdef MFNOCONFIDENCE
  if(!MFPolytopeTestVertexList(P,e))
   {
    sprintf(MFPolytopeErrorMsg,"Vertex list is bad.");
    MFSetError(e,12,RoutineName,MFPolytopeErrorMsg,__LINE__,__FILE__);
   }
#endif

  MFPolytopeUpdateFaceList(P,e);

  if(0){printf("after:\n");MFPrintPolytope(stdout,P,e);fflush(stdout);}

  return;
 }

#ifdef __cplusplus
}
#endif
