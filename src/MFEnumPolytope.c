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

static char *id="@(#) $Id: MFEnumPolytope.c,v 1.3 2007/02/13 01:22:33 mhender Exp $";

static char MFEnumPolytopeErrorMsg[256]="";

#include <MFEnumPolytope.h>
#include <MFEnumDualPolytope.h>
#include <MFPolytope.h>
#include <MFAtlas.h>
#include <MFAtlasFriends.h>
#include <MFChart.h>
#include <MFKVector.h>
#include <MFPrint.h>
#include <MFPrintEnum.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <MFErrorHandler.h>

#ifdef __cplusplus
 extern "C" {
#endif

struct MFCellSt {
             int nIndices;
             int *indices;
             int mFaceCells;
             int nFaceCells;
             int *faceCells;
             int nCellFaces;
             int mCellFaces;
             int *cellFaces;
            };

typedef struct MFCellSt *MFCell;

struct MFEnumPolytopeSt {
                  MFKVector *v;
                  int k;
                  int n;
                  int *nCells;
                  int *mCells;
                  MFCell **cells;
                 };

#define EPNumberOfCellsToAllocate 1

/* Routines special for access to other data structs that aren't public */

MFChart MFAtlasChart(MFAtlas,int,MFErrorHandler);
int MFPolytopeIntersectIndexSets(MFPolytope,int,int,int*,MFErrorHandler);

static int MFIntersectIndices(int,int*,int,int*,int*,MFErrorHandler);
static void MFCreateZeroCellV(MFCell*,MFErrorHandler,int,...);
static void MFCreateZeroCell(MFCell*,int,int*,MFErrorHandler);
static void MFCreateZeroCellWithEdgeList(MFCell*,int,int*,MFErrorHandler);
static void MFFreeCell(MFCell,MFErrorHandler);
static void MFCreateEmptyKComplex(int,int,MFEnumPolytope*,int,MFErrorHandler);
static int MFIndicesEqual(int,int*,int,int*,MFErrorHandler);
static void EPSortPolytopeIndices(int,int*,MFErrorHandler);
static void EPCollapseIndices(int*,int*,MFErrorHandler);
static int MFEnumPolytopeVerticesOnCell(MFEnumPolytope,int,int,int,int**,MFErrorHandler);
static int MFEnumPolytopeIsCellExterior(MFAtlas,int,MFEnumPolytope,int,int,MFErrorHandler);

static int MFTestEnumPolytope(MFEnumPolytope,int,MFErrorHandler);

MFEnumPolytope MFEnumeratePolytope(MFPolytope P, MFErrorHandler e)
 {
  static char RoutineName[]={"MFEnumeratePolytope"};
  int    i,j,n;
  MFEnumPolytope EP=NULL;
  int *indices;
  MFKVector s;
  int III;
  int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s\n",RoutineName);fflush(stdout);}
#endif

/* In 2d, boundary complex is a polygon, put edges head to tail, so interior is
                   to left */
/* In 3d, boundary complex is a polyhedron, put edges of each 2-cell head to
                   tail with outward normal */

  MFCreateEmptyKComplex(MFPolytopeDimension(P,e),MFPolytopeDimension(P,e),&EP,MFPolytopeNumberOfVertices(P,e),e);
  if(MFPolytopeNumberOfFaces(P,e)>0)
    indices=(int*)malloc(MFPolytopeNumberOfFaces(P,e)*sizeof(int));
   else
    indices=NULL;

#ifndef MFNOSAFETYNET
  if(MFPolytopeNumberOfFaces(P,e)> 0 && indices==NULL)
   {
    sprintf(MFEnumPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",MFPolytopeNumberOfFaces(P,e)*sizeof(int));
    MFSetError(e,12,RoutineName,MFEnumPolytopeErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  for(III=0;III<MFPolytopeNumberOfFaces(P,e);III++)indices[III]=-1;

  for(i=0;i<MFPolytopeNumberOfVertices(P,e);i++)
   {
    n=MFPolytopeNumberOfVertexIndices(P,i,e);
    for(j=0;j<n;j++)indices[j]=MFPolytopeVertexIndex(P,i,j,e);
    MFCreateZeroCell(&((EP->cells[0])[i]),n,indices,e);

    s=MFCreateKVector(MFPolytopeDimension(P,e),e);
    MFPolytopeVertex(P,i,s,e);
    EP->v[i]=s;
   }

#ifdef MFALLOWVERBOSE
  if(verbose){printf("  call MFEnumerate\n");fflush(stdout);}
#endif

  MFEnumerate(EP,e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("  back from MFEnumerate\n");fflush(stdout);}
#endif

  if(indices!=NULL)free(indices);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("  done %s\n",RoutineName);fflush(stdout);}
#endif

  return EP;
 }

void MFEnumerate(MFEnumPolytope EP, MFErrorHandler e)
 {
  static char RoutineName[]={"MFEnumerate"};
  int    intersection[1000]={0}; /* ??? */

  int    nv,ne,ni;
  int    iv,jv;
  int     n1, n2;
  int    *i1,*i2;
  int    d,k;
  int    mf;
  int    i,n;
  int    alreadyThere;
  int    j,f;
  MFCell tmp;
  int III;
  int    verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("\n\n\n\nMFEnumerate\n");fflush(stdout);}
#endif

/* Pairwise intersection of vertices */

  k=EP->k;

  for(d=1;d<k+1;d++)
   {

#ifdef MFALLOWVERBOSE
    if(verbose){printf("  intersecting %d cells to build %d cells\n",d-1,d);fflush(stdout);}
#endif

    nv=EP->nCells[d-1];

#ifdef MFALLOWVERBOSE
    if(verbose){printf("  There are %d %d-cells\n",nv,d-1);fflush(stdout);}
#endif

    if(EP->mCells[d]==0)
     {
      EP->mCells[d]=EPNumberOfCellsToAllocate;
      EP->cells[d]=(MFCell*)realloc(EP->cells[d],EP->mCells[d]*sizeof(MFCell));

#ifndef MFNOSAFETYNET
      if(EP->cells[d]==NULL)
       {
        sprintf(MFEnumPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",EP->mCells[d]*sizeof(MFCell));
        MFSetError(e,12,RoutineName,MFEnumPolytopeErrorMsg,__LINE__,__FILE__);
        MFErrorHandlerOutOfMemory(e);
        return;
       }
#endif

     }

#ifdef MFALLOWVERBOSE
    if(verbose)
     {
      printf("  There is space for %d %d-cells\n",EP->mCells[d],d);fflush(stdout);
      printf("  There are %d %d-cells\n",nv,d);fflush(stdout);
     }
#endif

    ne=0;
    for(iv=0;iv<nv-1;iv++)
     {
      for(jv=iv+1;jv<nv;jv++)
       {
        n1=((EP->cells[d-1])[iv])->nIndices;
        i1=((EP->cells[d-1])[iv])->indices;
        n2=((EP->cells[d-1])[jv])->nIndices;
        i2=((EP->cells[d-1])[jv])->indices;
        ni=MFIntersectIndices(n1,i1,n2,i2,intersection,e);

#ifdef MFALLOWVERBOSE
        if(verbose)
         {
          fprintf(stdout,"pair %2d,%2d ",iv,jv);
          fprintf(stdout,"      indices [");
          for(i=0;i<n1;i++)
           {
            if(i>0)fprintf(stdout," ");
            fprintf(stdout,"%2d",i1[i]);
           }
          fprintf(stdout,"], and [");
          for(i=0;i<n2;i++)
           {
            if(i>0)fprintf(stdout," ");
            fprintf(stdout,"%2d",i2[i]);
           }
          fprintf(stdout,"]      intersection [");
          for(i=0;i<ni;i++)
           {
            if(i>0)fprintf(stdout," ");
            fprintf(stdout,"%2d",intersection[i]);
           }
          fprintf(stdout,"]\n");fflush(stdout);
         }
#endif

        if(ni>=k-d)
         {

#ifdef MFALLOWVERBOSE
          if(verbose){printf("  Intersection large enough\n");fflush(stdout);}
#endif

          alreadyThere=0;
          n=-1;
          for(i=0;i<ne;i++)
           {
            n2=((EP->cells[d])[i])->nIndices;
            i2=((EP->cells[d])[i])->indices;
            if(MFIndicesEqual(ni,intersection,n2,i2,e)){alreadyThere=1;n=i;}
           }
          if(!alreadyThere)
           {

#ifdef MFALLOWVERBOSE
            if(verbose){printf("  Found a new %d-cell\n",d);fflush(stdout);}
#endif

            if(EP->mCells[d]==0)
             {
              EP->mCells[d]=EPNumberOfCellsToAllocate;

#ifdef MFALLOWVERBOSE
              if(verbose){printf("  increase size of EP->cells[%d] to %d\n",d,EP->mCells[d]);fflush(stdout);}
#endif

              EP->cells[d]=(MFCell*)realloc(EP->cells[d],EP->mCells[d]*sizeof(MFCell));

#ifndef MFNOSAFETYNET
              if(EP->cells[d]==NULL)
               {
                sprintf(MFEnumPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",EP->mCells[d]*sizeof(MFCell));
                MFSetError(e,12,RoutineName,MFEnumPolytopeErrorMsg,__LINE__,__FILE__);
                MFErrorHandlerOutOfMemory(e);
                return;
               }
#endif

             }else if(ne>=EP->mCells[d])
             {
              EP->mCells[d]+=EPNumberOfCellsToAllocate;

#ifdef MFALLOWVERBOSE
              if(verbose){printf("  allocate %d EP->cells[%d]\n",EP->mCells[d],d);fflush(stdout);}
#endif

              EP->cells[d]=(MFCell*)realloc(EP->cells[d],EP->mCells[d]*sizeof(MFCell));

#ifndef MFNOSAFETYNET
              if(EP->cells[d]==NULL)
               {
                sprintf(MFEnumPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",EP->mCells[d]*sizeof(MFCell));
                MFSetError(e,12,RoutineName,MFEnumPolytopeErrorMsg,__LINE__,__FILE__);
                MFErrorHandlerOutOfMemory(e);
                return;
               }
#endif

             }

#ifdef MFALLOWVERBOSE
            if(verbose)
             {
              printf("  This will be %d-cell number %d, address=0x%8.8x, location=0x0%8.8x\n",d,ne,(EP->cells[d])+ne,(EP->cells[d])[ne]);fflush(stdout);
              printf("  sizeof(struct MFCellSt)=%d\n",sizeof(struct MFCellSt));fflush(stdout);
             }
#endif

            tmp=(MFCell)malloc(sizeof(struct MFCellSt));

#ifdef MFALLOWVERBOSE
            if(verbose){printf("  done malloc\n");fflush(stdout);}
#endif

            (EP->cells[d])[ne]=tmp;

#ifdef MFALLOWVERBOSE
            if(verbose){printf("  done assignment\n");fflush(stdout);}
#endif


#ifndef MFNOSAFETYNET
            if(EP->cells[d][ne]==NULL)
             {
              sprintf(MFEnumPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct MFCellSt));
              MFSetError(e,12,RoutineName,MFEnumPolytopeErrorMsg,__LINE__,__FILE__);
              MFErrorHandlerOutOfMemory(e);
              return;
             }
#endif

#ifdef MFALLOWVERBOSE
            if(verbose){printf("  done allocating cell\n");fflush(stdout);}
#endif

            ((EP->cells[d])[ne])->nIndices=ni;

#ifdef MFALLOWVERBOSE
            if(verbose){printf("  done setting number of indices\n");fflush(stdout);}
#endif

            if(ni>0)
             {

#ifdef MFALLOWVERBOSE
              if(verbose){printf("  allocating space for %d indices\n",ni);fflush(stdout);}
#endif

              ((EP->cells[d])[ne])->indices=(int*)malloc(ni*sizeof(int));

#ifndef MFNOSAFETYNET
              if((EP->cells[d][ne])->indices==NULL)
               {
                sprintf(MFEnumPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",ni*sizeof(int));
                MFErrorHandlerOutOfMemory(e);
                return;
               }
#endif

#ifdef MFALLOWVERBOSE
              if(verbose){printf("  copying indices\n");fflush(stdout);}
#endif

              for(i=0;i<ni;i++)
                ((EP->cells[d])[ne])->indices[i]=intersection[i];
             }else
              ((EP->cells[d])[ne])->indices=NULL;

#ifdef MFALLOWVERBOSE
            if(verbose){printf("  done indices\n");fflush(stdout);}
#endif

            ((EP->cells[d])[ne])->mCellFaces=2;
            ((EP->cells[d])[ne])->nCellFaces=2;

#ifdef MFALLOWVERBOSE
            if(verbose){printf("  allocating space for cell faces\n");fflush(stdout);}
#endif

            ((EP->cells[d])[ne])->cellFaces=(int*)malloc(((EP->cells[d])[ne])->mCellFaces*sizeof(int));

#ifndef MFNOSAFETYNET
            if((EP->cells[d][ne])->cellFaces==NULL)
             {
              sprintf(MFEnumPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",((EP->cells[d])[ne])->mCellFaces*sizeof(int));
              MFSetError(e,12,RoutineName,MFEnumPolytopeErrorMsg,__LINE__,__FILE__);
              MFErrorHandlerOutOfMemory(e);
              return;
             }
#endif

#ifdef MFALLOWVERBOSE
            if(verbose){printf("  set cell faces\n");fflush(stdout);}
#endif

            ((EP->cells[d])[ne])->cellFaces[0]=iv;
            ((EP->cells[d])[ne])->cellFaces[1]=jv;

            ((EP->cells[d])[ne])->mFaceCells=0;
            ((EP->cells[d])[ne])->nFaceCells=0;
            ((EP->cells[d])[ne])->faceCells=NULL;

#ifdef MFALLOWVERBOSE
            if(verbose){printf("  done creating %d-cell number %d\n",d,ne);fflush(stdout);}
#endif
            n=ne;
            ne++;
           }else{

#ifdef MFALLOWVERBOSE
            if(verbose){printf("  Found an old %d-cell\n",d);fflush(stdout);}
#endif
            alreadyThere=0;
            mf=((EP->cells[d])[n])->nCellFaces;
            for(i=0;i<mf;i++)
              if(((EP->cells[d])[n])->cellFaces[i]==iv)alreadyThere=1;
            if(!alreadyThere)
             {
              if(((EP->cells[d])[n])->nCellFaces>=((EP->cells[d])[n])->mCellFaces)
               {
                ((EP->cells[d])[n])->mCellFaces++;
                ((EP->cells[d])[n])->cellFaces=(int*)realloc(((EP->cells[d])[n])->cellFaces,((EP->cells[d])[n])->mCellFaces*sizeof(int));

#ifndef MFNOSAFETYNET
                if((EP->cells[d][n])->cellFaces==NULL)
                 {
                  sprintf(MFEnumPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",((EP->cells[d])[n])->mCellFaces*sizeof(int));
                  MFSetError(e,12,RoutineName,MFEnumPolytopeErrorMsg,__LINE__,__FILE__);
                  MFErrorHandlerOutOfMemory(e);
                  return;
                 }
#endif
               }

#ifdef MFALLOWVERBOSE
              if(verbose){printf("  Add a cell face\n",d);fflush(stdout);}
#endif
              ((EP->cells[d])[n])->cellFaces[((EP->cells[d])[n])->nCellFaces]=iv;
              (((EP->cells[d])[n])->nCellFaces)++;
             }
            alreadyThere=0;
            mf=((EP->cells[d])[n])->nCellFaces;
            for(i=0;i<mf;i++)
              if(((EP->cells[d])[n])->cellFaces[i]==jv)alreadyThere=1;
            if(!alreadyThere)
             {
              if(((EP->cells[d])[n])->nCellFaces>=((EP->cells[d])[n])->mCellFaces)
               {
                ((EP->cells[d])[n])->mCellFaces++;
                ((EP->cells[d])[n])->cellFaces=(int*)realloc(((EP->cells[d])[n])->cellFaces,((EP->cells[d])[n])->mCellFaces*sizeof(int));

#ifndef MFNOSAFETYNET
                if((EP->cells[d][n])->cellFaces==NULL)
                 {
                  sprintf(MFEnumPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",((EP->cells[d])[n])->mCellFaces*sizeof(int));
                  MFSetError(e,12,RoutineName,MFEnumPolytopeErrorMsg,__LINE__,__FILE__);
                  MFErrorHandlerOutOfMemory(e);
                  return;
                 }
#endif

               }
              ((EP->cells[d])[n])->cellFaces[((EP->cells[d])[n])->nCellFaces]=jv;
              (((EP->cells[d])[n])->nCellFaces)++;
             }
           }

/* Face Cells of one lower dimension */

          mf=((EP->cells[d])[n])->nCellFaces;
          for(i=0;i<mf;i++)
           {
            f=((EP->cells[d])[n])->cellFaces[i];
            if(((EP->cells[d-1])[f])->nFaceCells>0)
             {
              alreadyThere=0;
              for(j=0;j<((EP->cells[d-1])[f])->nFaceCells;j++)
                if(((EP->cells[d-1])[f])->faceCells[j]==n)alreadyThere=1;
              if(!alreadyThere)
               {
                j=((EP->cells[d-1])[f])->nFaceCells+1;
                if(j>=((EP->cells[d-1])[f])->mFaceCells)
                 {
                  ((EP->cells[d-1])[f])->mFaceCells=j;
                  ((EP->cells[d-1])[f])->faceCells=(int*)realloc(((EP->cells[d-1])[f])->faceCells,((EP->cells[d-1])[f])->mFaceCells*sizeof(int));

#ifndef MFNOSAFETYNET
                  if((EP->cells[d-1][f])->faceCells==NULL)
                   {
                    sprintf(MFEnumPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",((EP->cells[d-1])[f])->mFaceCells*sizeof(int));
                    MFSetError(e,12,RoutineName,MFEnumPolytopeErrorMsg,__LINE__,__FILE__);
                    MFErrorHandlerOutOfMemory(e);
                    return;
                   }
#endif
                 }
                ((EP->cells[d-1])[f])->nFaceCells=j;
                ((EP->cells[d-1])[f])->faceCells[j-1]=n;
               }
             }else{
              if(0>=((EP->cells[d-1])[f])->mFaceCells)
               {
                ((EP->cells[d-1])[f])->mFaceCells=1;
                ((EP->cells[d-1])[f])->faceCells=(int*)malloc(((EP->cells[d-1])[f])->mFaceCells*sizeof(int));

#ifndef MFNOSAFETYNET
                if((EP->cells[d-1][f])->faceCells==NULL)
                 {
                  sprintf(MFEnumPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",((EP->cells[d-1])[f])->mFaceCells*sizeof(int));
                  MFSetError(e,12,RoutineName,MFEnumPolytopeErrorMsg,__LINE__,__FILE__);
                  MFErrorHandlerOutOfMemory(e);
                  return;
                 }
#endif
               }
              ((EP->cells[d-1])[f])->nFaceCells=1;
              ((EP->cells[d-1])[f])->faceCells[0]=n;
             }
           }
         }
       }
      EP->nCells[d]=ne;
     }

#ifdef MFALLOWVERBOSE
    if(verbose){printf("  done intersecting %d cells to build %d cells\n",d-1,d);fflush(stdout);}
#endif

   }

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    MFPrintEnumPolytope(stdout,EP,e);
    printf("done %s\n",RoutineName);fflush(stdout);
   }
#endif

  return;
 }

int MFIntersectIndices(int n1,int *inter1,int n2,int *inter2,int *inter, MFErrorHandler e)
 {
  static char RoutineName[]={"MFIntersectIndices"};
  int n;
  int i1,i2;

/*if(!MFTestVertexIndexOrder(n1,inter1,e))
   {
    printf("Error: MFIntersectIndices at entry, indices not sorted\n");
    fflush(stdout);
   }
  if(!MFTestVertexIndexOrder(n2,inter2,e))
   {
    printf("Error: MFIntersectIndices at entry, indices not sorted\n");
    fflush(stdout);
   }*/

  n=0;
  i1=0;i2=0;
  while(i1<n1&&i2<n2)
   {
    if(inter1[i1]==inter2[i2])
     {
      if(inter!=NULL)inter[n]=inter1[i1];
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

/*if(inter!=NULL && !MFTestVertexIndexOrder(n,inter,e))
   {
    printf("Error: MFIntersectIndices at exit, indices not sorted\n");
    fflush(stdout);
   }*/

  return n;
 }

void MFCreateZeroCellV(MFCell *cell,MFErrorHandler e, int n,...)
 {
  static char RoutineName[]={"MFCreateZeroCellV"};
  va_list indices;
  int i,j;

  *cell=(MFCell)malloc(sizeof(struct MFCellSt));

#ifndef MFNOSAFETYNET
  if(*cell==NULL)
   {
    sprintf(MFEnumPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct MFCellSt));
    MFSetError(e,12,RoutineName,MFEnumPolytopeErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  (*cell)->nIndices=n;
  (*cell)->indices=(int*)malloc(n*sizeof(int));

#ifndef MFNOSAFETYNET
  if((*cell)->indices==NULL)
   {
    sprintf(MFEnumPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(int));
    MFSetError(e,12,RoutineName,MFEnumPolytopeErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  (*cell)->mCellFaces=0;
  (*cell)->nCellFaces=0;
  (*cell)->cellFaces=NULL;
  (*cell)->mFaceCells=0;
  (*cell)->nFaceCells=0;
  (*cell)->faceCells=NULL;

  va_start(indices,n);
  for(i=0;i<n;i++)
    (*cell)->indices[i]=va_arg(indices,int);
  va_end(indices);

  return;
 }

void MFCreateZeroCell(MFCell *cell,int n,int* indices, MFErrorHandler e)
 {
  static char RoutineName[]={"MFCreateZeroCell"};
  int i,j;

  *cell=(MFCell)malloc(sizeof(struct MFCellSt));

#ifndef MFNOSAFETYNET
  if(*cell==NULL)
   {
    sprintf(MFEnumPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct MFCellSt));
    MFSetError(e,12,RoutineName,MFEnumPolytopeErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  (*cell)->nIndices=n;
  if(n>0)
   {
    (*cell)->indices=(int*)malloc(n*sizeof(int));

#ifndef MFNOSAFETYNET
    if((*cell)->indices==NULL)
     {
      sprintf(MFEnumPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(int));
      MFSetError(e,12,RoutineName,MFEnumPolytopeErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

   }else{
    (*cell)->indices=NULL;
   }
  (*cell)->mCellFaces=0;
  (*cell)->nCellFaces=0;
  (*cell)->cellFaces=NULL;
  (*cell)->mFaceCells=0;
  (*cell)->nFaceCells=0;
  (*cell)->faceCells=NULL;

  for(i=0;i<n;i++)(*cell)->indices[i]=indices[i];

  return;
 }

void  MFCreateEmptyKComplex(int k, int n, MFEnumPolytope *EP,int m, MFErrorHandler e)
 {
  static char RoutineName[]={"MFCreateEmptyKComplex"};
  int i,j;

  *EP=(MFEnumPolytope)malloc(sizeof(struct MFEnumPolytopeSt));

#ifndef MFNOSAFETYNET
  if(*EP==NULL)
   {
    sprintf(MFEnumPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct MFEnumPolytopeSt));
    MFSetError(e,12,RoutineName,MFEnumPolytopeErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  (*EP)->k=k;
  (*EP)->n=n;
  (*EP)->nCells=(int*)malloc((k+1)*sizeof(int));

#ifndef MFNOSAFETYNET
  if((*EP)->nCells==NULL)
   {
    sprintf(MFEnumPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",(k+1)*sizeof(int));
    MFSetError(e,12,RoutineName,MFEnumPolytopeErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  (*EP)->mCells=(int*)malloc((k+1)*sizeof(int));

#ifndef MFNOSAFETYNET
  if((*EP)->mCells==NULL)
   {
    sprintf(MFEnumPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",(k+1)*sizeof(int));
    MFSetError(e,12,RoutineName,MFEnumPolytopeErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  (*EP)->cells=(MFCell**)malloc((k+1)*sizeof(MFCell*));

#ifndef MFNOSAFETYNET
  if((*EP)->cells==NULL)
   {
    sprintf(MFEnumPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",(k+1)*sizeof(MFCell*));
    MFSetError(e,12,RoutineName,MFEnumPolytopeErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  (*EP)->mCells[0]=m;
  (*EP)->nCells[0]=m;
  if(m>0)
   {
    (*EP)->cells[0]=(MFCell*)malloc((*EP)->mCells[0]*sizeof(MFCell));

#ifndef MFNOSAFETYNET
    if((*EP)->cells[0]==NULL)
     {
      sprintf(MFEnumPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",((*EP)->mCells[0])*sizeof(MFCell));
      MFSetError(e,12,RoutineName,MFEnumPolytopeErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

    (*EP)->v=(MFKVector*)malloc((*EP)->mCells[0]*sizeof(MFKVector));

#ifndef MFNOSAFETYNET
    if((*EP)->v==NULL)
     {
      sprintf(MFEnumPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",((*EP)->mCells[0])*sizeof(MFKVector));
      MFSetError(e,12,RoutineName,MFEnumPolytopeErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

    for(j=0;j<(*EP)->mCells[0];j++)
     {
      ((*EP)->cells[0])[j]=NULL;
      (*EP)->v[j]=NULL;
     }
   }else{
    (*EP)->mCells[0]=EPNumberOfCellsToAllocate;
    (*EP)->nCells[0]=0;
    (*EP)->cells[0]=(MFCell*)malloc((*EP)->mCells[0]*sizeof(MFCell));

#ifndef MFNOSAFETYNET
    if((*EP)->cells[0]==NULL)
     {
      sprintf(MFEnumPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",((*EP)->mCells[0])*sizeof(MFCell));
      MFSetError(e,12,RoutineName,MFEnumPolytopeErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

    (*EP)->v=(MFKVector*)malloc((*EP)->mCells[0]*sizeof(MFKVector));

#ifndef MFNOSAFETYNET
    if((*EP)->v==NULL)
     {
      sprintf(MFEnumPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",((*EP)->mCells[0])*sizeof(MFKVector));
      MFSetError(e,12,RoutineName,MFEnumPolytopeErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

    for(j=0;j<(*EP)->mCells[0];j++)
     {
      ((*EP)->cells[0])[j]=NULL;
      (*EP)->v[j]=NULL;
     }
   }

  for(i=1;i<=k;i++)
   {
    (*EP)->mCells[i]=EPNumberOfCellsToAllocate;
    (*EP)->nCells[i]=0;
    (*EP)->cells[i]=(MFCell*)malloc((*EP)->mCells[i]*sizeof(MFCell));

#ifndef MFNOSAFETYNET
    if((*EP)->cells[i]==NULL)
     {
      sprintf(MFEnumPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",((*EP)->mCells[i])*sizeof(MFCell));
      MFSetError(e,12,RoutineName,MFEnumPolytopeErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

    for(j=0;j<(*EP)->mCells[i];j++)((*EP)->cells[i])[j]=NULL;
   }

  return;
 }

int MFIndicesEqual(int n1,int *inter1,int n2,int *inter2, MFErrorHandler e)
 {
  static char RoutineName[]={"MFIndicesEqual"};
  int result;
  int i;

  result=0;
  if(n1!=n2)
   {
    result=0;
   }else{
    result=1;
    for(i=0;i<n1;i++)
      if(inter1[i]!=inter2[i])result=0;
   } 
  return result;
 }

void MFFreeEnumPolytope(MFEnumPolytope EP, MFErrorHandler e)
 {
  static char RoutineName[]={"MFFreeEnumPolytope"};
  int i,j;
  int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("MFFreeEnumPolytope\n");fflush(stdout);}
#endif

  if(EP==NULL)return;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("  delete vertices\n");fflush(stdout);}
#endif

  for(i=0;i<EP->nCells[0];i++)
   {

#ifdef MFALLOWVERBOSE
    if(verbose){printf("  %d \n",i);fflush(stdout);}
#endif

    if(EP->v[i]!=NULL)MFFreeKVector(EP->v[i],e);
   }

#ifdef MFALLOWVERBOSE
  if(verbose){printf("  delete vertex list\n");fflush(stdout);}
#endif

  if(EP->v!=NULL)
   {
    free(EP->v);
   }

#ifdef MFALLOWVERBOSE
  if(verbose){printf("  delete cells\n");fflush(stdout);}
#endif

  if(EP->nCells!=NULL)
   {
    for(j=0;j<=EP->k;j++)
     {

#ifdef MFALLOWVERBOSE
      if(verbose){printf("    dimension %d\n",j);fflush(stdout);}
#endif

      for(i=0;i<EP->nCells[j];i++)
       {

#ifdef MFALLOWVERBOSE
        if(verbose){printf("      %d\n",i);fflush(stdout);}
#endif

        if((EP->cells[j])[i]!=NULL)MFFreeCell((EP->cells[j])[i],e);
       }

#ifdef MFALLOWVERBOSE
      if(verbose){printf("   delete list of cells of dimension %d\n",j);fflush(stdout);}
#endif

      if(EP->cells[j]!=NULL)
       {
        free(EP->cells[j]);
       }
     }
   }

#ifdef MFALLOWVERBOSE
  if(verbose){printf("   delete cell counts\n");fflush(stdout);}
#endif

  if(EP->nCells!=NULL)
   {
    free(EP->nCells);
   }
  if(EP->mCells!=NULL)
   {
    free(EP->mCells);
   }

#ifdef MFALLOWVERBOSE
  if(verbose){printf("   delete list of list of cells\n");fflush(stdout);}
#endif

  if(EP->cells!=NULL)
   {
    free(EP->cells);
   }
     
#ifdef MFALLOWVERBOSE
  if(verbose){printf("   free data structure\n");fflush(stdout);}
#endif

  free(EP);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("done MFFreeEnumPolytope\n");fflush(stdout);}
#endif

  return;
 }

void MFFreeCell(MFCell C, MFErrorHandler e)
 {
  static char RoutineName[]={"MFFreeCell"};
  int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("MFFreeCell\n");fflush(stdout);}
#endif

  if(C==NULL)return;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("   delete indices\n");fflush(stdout);}
#endif

  if(C->indices!=NULL)
   {
    free(C->indices);
   }

#ifdef MFALLOWVERBOSE
  if(verbose){printf("   delete faceCells\n");fflush(stdout);}
#endif

  if(C->faceCells!=NULL)
   {
    free(C->faceCells);
   }

#ifdef MFALLOWVERBOSE
  if(verbose){printf("   delete cellFaces\n");fflush(stdout);}
#endif

  if(C->cellFaces!=NULL)
   {
    free(C->cellFaces);
   }

#ifdef MFALLOWVERBOSE
  if(verbose){printf("   free data structure\n");fflush(stdout);}
#endif

  free(C);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("done MFFreeCell\n");fflush(stdout);}
#endif

  return;
 }

long MFEnumPolytopeNumberOfVertices(MFEnumPolytope EP, MFErrorHandler e)
 {
  static char RoutineName[]={"MFEnumPolytopeNumberOfVertices"};
  return EP->nCells[0];
 }

int MFEnumPolytopeNumberOfCells(MFEnumPolytope EP,int dim, MFErrorHandler e)
 {
  static char RoutineName[]={"MFEnumPolytopeNumberOfCells"};
  if(dim>EP->n)return 0;
  return EP->nCells[dim];
 }

int MFEnumPolytopeNumberOfCellFaces(MFEnumPolytope EP,int dim,long cell, MFErrorHandler e)
 {
  static char RoutineName[]={"MFEnumPolytopeNumberOfCellFaces"};
  return ((EP->cells[dim])[cell])->nCellFaces;
 }

long MFEnumPolytopeCellFace(MFEnumPolytope EP,int dim,long cell,int face, MFErrorHandler e)
 {
  static char RoutineName[]={"MFEnumPolytopeCellFace"};
  int result;

#ifdef MFNOCONFIDENCE
  if(dim<0||dim>EP->k)fprintf(stderr,"Error at line %d in file %s, you asked for a cell face of a %d-cell, but k=%d\n",__LINE__,__FILE__,dim,EP->k);
  if(cell<0||cell>=EP->nCells[dim])fprintf(stderr,"Error at line %d in file %s, you asked for a cell face of %d-cell %d, but there are %d cells of this dimension \n",__LINE__,__FILE__,dim,cell,EP->nCells[dim]);
  if(face<0||face>=((EP->cells[dim])[cell])->nCellFaces)fprintf(stderr,"Error at line %d in file %s, you asked for cell face %d of %d-cell %d, but there are %d cell faces\n",__LINE__,__FILE__,face,dim,cell,((EP->cells[dim])[cell])->nCellFaces);
#endif

/* Signs: if cell<0 return -( (...[-cell])->cellFaces[face]) ? */

  result=((EP->cells[dim])[cell])->cellFaces[face];

  return result;
 }

int MFEnumPolytopeNumberOfFaceCells(MFEnumPolytope EP,int dim,long face, MFErrorHandler e)
 {
  static char RoutineName[]={"MFEnumPolytopeNumberOfFaceCells"};
  return ((EP->cells[dim])[face])->nFaceCells;
 }

int MFEnumPolytopeNumberOfCellIndices(MFEnumPolytope EP,int dim,long cell, MFErrorHandler e)
 {
  static char RoutineName[]={"MFEnumPolytopeNumberOfCellIndices"};

#ifdef MFNOCONFIDENCE
  if(EP==NULL)
   {
    printf("EP NULL in MFEnumPolytopeNumberOfCellIndices\n");
    fflush(stdout);abort();
   }
  if(dim<0 || dim>EP->k)
   {
    printf("In MFEnumPolytopeNumberOfCellIndices, dim (%d) is invalid. Must be positice and less than or equal to %d.\n",dim,EP->k);
    fflush(stdout);abort();
   }
  if(cell<0 || !(cell<EP->nCells[dim]))
   {
    printf("In MFEnumPolytopeNumberOfCellIndices, cell (%d) is invalid. Must be positice and less than %d.\n",cell,EP->nCells[dim]);
    fflush(stdout);abort();
   }
#endif

  return ((EP->cells[dim])[cell])->nIndices;
 }

int MFEnumPolytopeCellIndex(MFEnumPolytope EP,int dim,long cell,int i, MFErrorHandler e)
 {
  static char RoutineName[]={"MFEnumPolytopeCellIndex"};
  return ((EP->cells[dim])[cell])->indices[i];
 }

long MFEnumPolytopeFaceCell(MFEnumPolytope EP,int dim,long face,int cell, MFErrorHandler e)
 {
  static char RoutineName[]={"MFEnumPolytopeFaceCell"};
  return ((EP->cells[dim])[face])->faceCells[cell];
 }

void EPSortPolytopeIndices(int n,int *list, MFErrorHandler e)
 {
  static char RoutineName[]={"EPSortPolytopeIndices"};
  int i,j,t;

  for(i=0;i<n-1;i++)
   for(j=i;j<n;j++)
    {
     if(list[j]<list[i])
      {
       t=list[j];
       list[j]=list[i];
       list[i]=t;
      }
    }
 }

void EPCollapseIndices(int *n,int *list, MFErrorHandler e)
 {
  static char RoutineName[]={"EPCollapseIndices"};
  int i,j;

  EPSortPolytopeIndices(*n,list,e);
  for(i=0;i<*n-1;i++)
   {
    if(list[i]==list[i+1])
     {
      (*n)--;
      for(j=i+1;j<*n;j++)list[j]=list[j+1];
      i--;
     }
   }

  return;
 }

void MFCreateZeroCellWithEdgeList(MFCell *cell,int n,int *edges, MFErrorHandler e)
 {
  static char RoutineName[]={"MFCreateZeroCellWithEdgeList"};
  int i,j;

  *cell=(MFCell)malloc(sizeof(struct MFCellSt));

#ifndef MFNOSAFETYNET
  if(*cell==NULL)
   {
    sprintf(MFEnumPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct MFCellSt));
    MFSetError(e,12,RoutineName,MFEnumPolytopeErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  (*cell)->nIndices=0;
  (*cell)->indices=NULL;
  (*cell)->mCellFaces=0;
  (*cell)->nCellFaces=0;
  (*cell)->cellFaces=NULL;
  (*cell)->mFaceCells=n;
  (*cell)->nFaceCells=n;
  (*cell)->faceCells=(int*)malloc(n*sizeof(int));

#ifndef MFNOSAFETYNET
  if((*cell)->faceCells==NULL)
   {
    sprintf(MFEnumPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(int));
    MFSetError(e,12,RoutineName,MFEnumPolytopeErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  for(i=0;i<n;i++)(*cell)->faceCells[i]=edges[i];

  return;
 }

int MFEnumPolytopeDimension(MFEnumPolytope EP, MFErrorHandler e)
 {
  static char RoutineName[]={"MFEnumPolytopeDimension"};
  return EP->k;
 }

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

int MFTestEnumPolytope(MFEnumPolytope EP,int closed, MFErrorHandler e)
 {
  static char RoutineName[]={"MFTestEnumPolytope"};
  int k;
  int d;
  int nCells,nFaces,nAdj;
  long cell,face,adj;
  long cellface,adjcell;
  int nFaceAdjs,nAdjFaces;
  long faceadj,adjface;
  int i;
  int present,presentTwice,passed;

  k=MFEnumPolytopeDimension(EP,e);

  passed=1;
  for(d=0;d<k+1;d++)
   {
    nCells=MFEnumPolytopeNumberOfCells(EP,d,e);
    for(cell=0;cell<nCells;cell++)
     {
      if(d-1>-1)
       {
        nFaces=MFEnumPolytopeNumberOfCellFaces(EP,d,cell,e);
        if(closed && nFaces<d+1)
         {
          passed=0;
          printf("MFTestEnumPolytope:\n %d-cell %ld has %d faces, fewer than %d\n",d,cell,nFaces,d+1);
         }
        for(face=0;face<nFaces;face++)
         {
          cellface=MFEnumPolytopeCellFace(EP,d,cell,face,e);
          presentTwice=0;
          for(i=0;i<nFaces;i++)
           {
            if(i!=face && cellface==MFEnumPolytopeCellFace(EP,d,cell,i,e))presentTwice=1;
           }
          if(presentTwice)
           {
            passed=0;
            printf("MFTestEnumPolytope:\n A face of %d-cell %ld (%d-cell %ld) appears twice in the face list.\n",d,cell,d-1,cellface);
            printf("The face list of %d-cell %ld is [",d,cell);
            for(i=0;i<nFaces;i++)
             {
              if(i>0)printf(",");
              printf("%d",MFEnumPolytopeCellFace(EP,d,cell,i,e));
             }
            printf("]\n");
           }
          nFaceAdjs=MFEnumPolytopeNumberOfFaceCells(EP,d-1,cellface,e);
          present=0;
          for(faceadj=0;faceadj<nFaceAdjs;faceadj++)
           {
            if(MFEnumPolytopeFaceCell(EP,d-1,cellface,faceadj,e)==cell)present=1;
           }
          if(!present)
           {
            passed=0;
            printf("MFTestEnumPolytope:\n Face of %d-cell %ld (%d-cell %ld) does not list %d-cell %ld as a cell it is contained in\n",d,cell,d-1,cellface,d,cell);
            printf("%d-cell %ld is contained in cells (%d-cells) [",d-1,cellface,d);
            for(faceadj=0;faceadj<nFaceAdjs;faceadj++)
             {
              if(faceadj>0)printf(",");
              printf("%d",MFEnumPolytopeFaceCell(EP,d-1,cellface,faceadj,e));
             }
            printf("]\n");
           }
         }
       }
      if(d+1<k)
       {
        nAdj=MFEnumPolytopeNumberOfFaceCells(EP,d,cell,e);
        for(adj=0;adj<nAdj;adj++)
         {
          adjcell=MFEnumPolytopeFaceCell(EP,d,cell,adj,e);
          presentTwice=0;
          for(i=0;i<nAdj;i++)
           {
            if(i!=adj && adjcell==MFEnumPolytopeFaceCell(EP,d,cell,i,e))presentTwice=1;
           }
          if(presentTwice)
           {
            passed=0;
            printf("MFTestEnumPolytope:\n A containing cell of %d-cell %ld (%d-cell %ld) appears twice in the list of containing cells.\n",d,cell,d+1,adjcell);
            printf("The containing cell list of %d-cell %ld is [",d,cell);
            for(i=0;i<nAdj;i++)
             {
              if(i>0)printf(",");
              printf("%d",MFEnumPolytopeFaceCell(EP,d,cell,i,e));
             }
            printf("]\n");
           }
          nAdjFaces=MFEnumPolytopeNumberOfCellFaces(EP,d+1,adjcell,e);
          present=0;
          for(adjface=0;adjface<nAdjFaces;adjface++)
           {
            if(MFEnumPolytopeCellFace(EP,d+1,adjcell,adjface,e)==cell)present=1;
           }
          if(!present)
           {
            passed=0;
            printf("MFTestEnumPolytope:\n Containing cell of %d-cell %ld (%d-cell %ld) does not list %d-cell %ld as a face\n",d,cell,d+1,adjcell,d,cell);
            printf("%d-cell %ld has faces (%d-cells) [",d+1,adjcell,d);
            for(adjface=0;adjface<nAdjFaces;adjface++)
             {
              if(adjface>0)printf(",");
              printf("%d",MFEnumPolytopeCellFace(EP,d+1,adjcell,adjface,e));
             }
            printf("]\n");
           }
         }
       }
     }
   }
  return passed;
 }

int MFEnumPolytopeIsCellExterior(MFAtlas A,int chart,MFEnumPolytope P,int dim,int cell, MFErrorHandler e)
 {
  static char RoutineName[]={"MFEnumPolytopeIsCellExterior"};
  int hasHyperCube,hasNotInBoth;
  int indicesOneIn;
  int result;
  int nV;
  int *Vs;
  MFCell c;
  int i;
  int *indices;
  int imin;
  double rmin,r,vij;
  MFPolytope PP;
  MFKVector vi;
  MFKVector vj;
  MFKVector vk;
  MFKVector v,diff;
  MFKVector diffi,diffj;
  double s,t,R;
  double a11,a12,a22,b1,b2;
  int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("MFEnumPolytopeIsCellExterior P=0x%8.8x--\n",P,e);fflush(stdout);}
#endif

  if(dim==MFEnumPolytopeDimension(P,e))return 0;

  PP=MFChartPolytope(MFAtlasChart(A,chart,e),e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("  %d-cell %d\n",dim,cell);fflush(stdout);}
#endif

#ifdef MFALLOWVERBOSE
  if(verbose&&dim==1)
   {
    indices=NULL;
    nV=0;
    nV=MFEnumPolytopeVerticesOnCell(P,dim,cell,nV,&indices,e);
    printf("There are %d vertices lying on edge %d [",nV,cell);
    for(i=0;i<nV;i++)
     {
      if(i>0)printf(",");
      printf("%d",indices[i]);
     }
    printf("]\n");fflush(stdout);
    free(indices);
   }
#endif

  c=(P->cells[dim])[cell];

  result=0;
  for(i=0;i<c->nIndices;i++)
   {
    if(MFAtlasIsHalfSpaceHyperCube(A,c->indices[i],e))
     {
      result=1;

#ifdef MFALLOWVERBOSE
      if(verbose){printf(" Cell lies on a hypercube face\n");fflush(stdout);}
#endif

      return result;
     }
    if(!MFAtlasIsHalfSpaceInBothPolytopes(A,c->indices[i],e))
     {
      result=1;

#ifdef MFALLOWVERBOSE
      if(verbose){printf(" Cell does not lie on a face shared by an adjacent polytope.\n");fflush(stdout);}
#endif

      return result;
     }
   }

  result=0;

  nV=0;
  indices=NULL;
  nV=MFEnumPolytopeVerticesOnCell(P,dim,cell,nV,&indices,e);

  indicesOneIn=0;

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("There are %d vertices lying on the cell [",nV);
    for(i=0;i<nV;i++)
     {
      if(i>0)printf(",");
      printf("%d",indices[i]);
     }
    printf("]\n");
   }
#endif

  for(i=0;i<nV;i++)
   {

#ifdef MFALLOWVERBOSE
    if(verbose)
     {
      printf("  Vertex %d is ",indices[i]);
      vi=MFCreateKVector(MFAtlasK(A,e),e);
      MFPolytopeVertex(PP,indices[i],vi,e);
      MFPrintKVector(stdout,vi,e);
      MFFreeKVector(vi,e);
      printf(" radius %lf",MFPolytopeRadiusOfVertex(PP,indices[i],e));
     }
#endif

    if(MFPolytopeRadiusOfVertex(PP,indices[i],e)<MFAtlasChartRadius(A,chart,e),e)
     {

#ifdef MFALLOWVERBOSE
      if(verbose)printf(" and is in chart %d (radius %lf)\n",chart,MFAtlasChartRadius(A,chart,e));
#endif

      indicesOneIn=1;
     }
   }

  if(indicesOneIn)
   {

#ifdef MFALLOWVERBOSE
    if(verbose)printf(" One vertex is interior to the chart\n");
#endif

    result=0;
    if(indices!=NULL)free(indices);
    return result;
   }else{

#ifdef MFALLOWVERBOSE
    if(verbose)printf(" All vertices are exterior to the chart\n");
#endif

   }

  if(dim==0)
   {
    result=0;
    if(MFPolytopeRadiusOfVertex(MFChartPolytope(MFAtlasChart(A,chart,e),e),cell,e)>MFAtlasChartRadius(A,chart,e))result=1;
/* }else if(dim==1)
   {
    vi=MFCreateKVector(MFAtlasK(A,e),e);
    vj=MFCreateKVector(MFAtlasK(A,e),e);
    diff=MFCreateKVector(MFAtlasK(A,e),e);
    v=MFCreateKVector(MFAtlasK(A,e),e);
    MFPolytopeVertex(PP,indices[0],vi,e);
    MFPolytopeVertex(PP,indices[1],vj,e);
    MFKVDiff(vi,vj,diff,e);
    t=MFKVDot(vi,diff,e)/MFKVNorm(diff,e);
    if(t<0.)t=0.;
    if(t>1.)t=1.;
    MFKVScaleMul(t,diff,vj,e);
    MFKVAdd(vi,vj,v,e);
    R=MFKVNorm(v,e);
    
    MFFreeKVector(vi,e);
    MFFreeKVector(vj,e);
    MFFreeKVector(diff,e);
    MFFreeKVector(v,e);
    result=0;
    if(R>MFAtlasChartRadius(A,chart,e))result=1;
   }else if(dim==2)
   {
    vi=MFCreateKVector(MFAtlasK(A,e),e);
    vj=MFCreateKVector(MFAtlasK(A,e),e);
    vk=MFCreateKVector(MFAtlasK(A,e),e);
    diffi=MFCreateKVector(MFAtlasK(A,e),e);
    diffj=MFCreateKVector(MFAtlasK(A,e),e);
    v=MFCreateKVector(MFAtlasK(A,e),e);
    MFPolytopeVertex(PP,indices[0],vi,e);
    MFPolytopeVertex(PP,indices[1],vj,e);
    MFPolytopeVertex(PP,indices[2],vk,e);
    MFKVDiff(vi,vj,diffi,e);
    MFKVDiff(vi,vk,diffj,e);
    a11=MFKVNorm(diffi,e);
    a12=MFKVDot(diffi,diffj,e);
    a22=MFKVNorm(diffj,e);
    b1=-MFKVDot(vi,diffi,e);
    b2=-MFKVDot(vi,diffj,e);
    s=( b1*a22-b2*a12)/(a11*a22-a12*a12);
    t=(-b1*a12+b2*a22)/(a11*a22-a12*a12);
    if(s<0.)s=0.;
    if(t<0.)t=0.;
    MFKVScaleMul(t,diffi,vj,e);
    MFKVScaleMul(s,diffj,vk,e);
    MFKVAdd(vi,vj,diffi,e);
    MFKVAdd(diffi,diffj,v,e);
    R=MFKVNorm(v,e);

    MFFreeKVector(vi,e);
    MFFreeKVector(vj,e);
    MFFreeKVector(vk,e);
    MFFreeKVector(diffi,e);
    MFFreeKVector(diffj,e);
    MFFreeKVector(v,e);
    result=0;
    if(R>MFAtlasChartRadius(A,chart,e))result=1;
*/
    if(indices!=NULL)free(indices);
    return result;
   }

  result=0;
  if(0 && nV>1)
   {

#ifdef MFALLOWVERBOSE
    if(verbose)printf("Punt -- Trying to eliminate a cell of dimension %d\n",dim);
#endif

    imin=0;
    rmin=MFPolytopeRadiusOfVertex(PP,indices[0],e);

#ifdef MFALLOWVERBOSE
    if(verbose)printf("   Vertex %d is radius %lf\n",indices[0],rmin);
#endif

    for(i=1;i<nV;i++)
     {
      r=MFPolytopeRadiusOfVertex(PP,indices[i],e);
      if(r<rmin)
       {
        rmin=r;
        imin=i;
       }
     }

#ifdef MFALLOWVERBOSE
    if(verbose)printf(" The closest vertex to the center is Vertex %d, radius %lf, ball is radius %lf\n",indices[imin],rmin,MFAtlasChartRadius(A,chart,e));
#endif

    vi=MFCreateKVector(MFAtlasK(A,e),e);
    vj=MFCreateKVector(MFAtlasK(A,e),e);
    MFPolytopeVertex(PP,indices[imin],vi,e);
    r=rmin*rmin;
    result=0;
    for(i=0;i<nV;i++)
     {
      if(i!=imin)
       {
        MFPolytopeVertex(PP,indices[i],vj,e);
        vij=MFKVDot(vi,vj,e);
        if(vij>=r)result=1;

#ifdef MFALLOWVERBOSE
        if(verbose)
         {
          printf("  v[%d]=",imin);
          MFPrintKVector(stdout,vi,e);
          printf("\n");
          printf("  v[%d]=",i);
          MFPrintKVector(stdout,vj,e);
          printf("\n");
          printf("    v[%d].v[%d]=%lf, v[%d].v[%d]=%lf\n",indices[i],indices[imin],vij,indices[imin],indices[imin],r);
          printf("    v[%d].(v[%d]-v[%d])=%lf, v[%d].(v[%d]-v[%d])=%lf\n",indices[i],indices[i],indices[imin],r-vij,indices[imin],indices[imin],indices[i],MFKVDot(vj,vj,e)-vij);
         }
#endif

       }
     }
    if(result&&rmin<MFAtlasChartRadius(A,chart,e))result=0;

#ifdef MFALLOWVERBOSE
    if(verbose)printf("   result is %d\n",result);
#endif

    MFFreeKVector(vi,e);
    MFFreeKVector(vj,e);
   }
   
  if(indices!=NULL)free(indices);

  return result;
 }

int MFEnumPolytopeVerticesOnCell(MFEnumPolytope P,int dim,int cell,int nV,int **indices, MFErrorHandler e)
 {
  static char RoutineName[]={"MFEnumPolytopeVerticesOnCell"};
  int nF,f;
  int n;
  MFCell c;
  int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("Counting the vertices contained on %d-cell %d, nV=%d\n",dim,cell,nV);fflush(stdout);}
#endif

  if(dim==0)
   {

#ifdef MFALLOWVERBOSE
    if(verbose){printf("   Single vertex %d, checking to see if it is in the list nV=%d\n",cell,nV);fflush(stdout);}
#endif

    if(*indices!=NULL)
     {
      n=0;
      for(f=0;f<nV;f++)if(((*indices)[f])==cell)n=1;
      if(n==1)
       {

#ifdef MFALLOWVERBOSE
        if(verbose){printf("     Already there, returning\n");fflush(stdout);}
#endif

        return nV;
       }
     }

#ifdef MFALLOWVERBOSE
    if(verbose){printf("     Not there, adding it\n");fflush(stdout);}
#endif

    if(*indices==NULL)
     {

#ifdef MFALLOWVERBOSE
      if(verbose){printf("     *indices is NULL, mallocing it\n");fflush(stdout);}
#endif

      *indices=(int*)malloc(10*sizeof(int));

#ifndef MFNOSAFETYNET
      if(*indices==NULL)
       {
        sprintf(MFEnumPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",10*sizeof(int));
        MFSetError(e,12,RoutineName,MFEnumPolytopeErrorMsg,__LINE__,__FILE__);
        MFErrorHandlerOutOfMemory(e);
        return -1;
       }
#endif

     }else{

#ifdef MFALLOWVERBOSE
      if(verbose){printf("     *indices is not NULL, reallocing it\n");fflush(stdout);}
#endif

      *indices=(int*)realloc((void*)(*indices),(nV+10)*sizeof(int));

#ifndef MFNOSAFETYNET
      if(*indices==NULL)
       {
        sprintf(MFEnumPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",(nV+10)*sizeof(int));
        MFSetError(e,12,RoutineName,MFEnumPolytopeErrorMsg,__LINE__,__FILE__);
        MFErrorHandlerOutOfMemory(e);
        return -1;
       }
#endif

     }

    (*indices)[nV]=cell;
    return nV+1;
   }


  c=(P->cells[dim])[cell];

  nF=c->nCellFaces;

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("   Pushing down to check the faces of this %d-cell [",dim);
    for(f=0;f<nF;f++)
     {
      if(f>0)printf(",");
      printf("%d",c->cellFaces[f]);
     }
    printf("]\n");
    fflush(stdout);
    printf("     There are %d faces\n",nF);fflush(stdout);
   }
#endif

  n=nV;
  for(f=0;f<nF;f++)
   {

#ifdef MFALLOWVERBOSE
    if(verbose){printf("     Check face %d, %d-cell %d\n",f,dim-1,c->faceCells[f]);fflush(stdout);}
#endif

    n=MFEnumPolytopeVerticesOnCell(P,dim-1,c->cellFaces[f],n,indices,e);
   }

  return n;
 }

MFEnumPolytope MFEnumerateAtlas(MFAtlas A, MFErrorHandler e)
 {
  static char RoutineName[]={"MFEnumerateAtlas"};
  int chart,nCharts;
  int dim;
  int k;
  int i,l;
  MFCell mfcell;
  MFCell cellFace;
  MFCell faceCell;
  int chart2;
  long nCells2,cell2;
  MFCell mfcell2;
  long cell,nCells;
  MFPolytope P;
  MFEnumPolytope *EPS;
  MFEnumPolytope Complex=NULL;
  int m;
  int indices[10000]={0}; /*!!!*/
  int hasHyperCube=0;
  int hasNotInBoth=0;
  int found;
  int failed;
  int verbose=0;
  int verboseResult=0;
  int verboseTerse=0;

  int nSame;
  int mSame;
  int maxSame;
  MFCell mfcellClosest;
  int cellClosest;

/* In 2d, edge dual takes you to neighboring polygon, and then the next edge CCW
                   which gives the next neighbor (sign?)*/
/* In 3d, face dual takes you to neighboring polyhedron. Edge on face in primal
                   walks to next face whcih gives next neighbor */

#ifdef MFALLOWVERBOSE
  if(verbose){printf("In MFEnumeratePolytopes\n");fflush(stdout);}
#endif

  nCharts=MFAtlasNumberOfCharts(A,e);
  k=MFAtlasK(A,e);
  EPS=(MFEnumPolytope*)malloc(nCharts*sizeof(MFEnumPolytope));

#ifndef MFNOSAFETYNET
  if(EPS==NULL)
   {
    sprintf(MFEnumPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",nCharts*sizeof(MFEnumPolytope));
    MFSetError(e,12,RoutineName,MFEnumPolytopeErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

#ifdef MFALLOWVERBOSE
  if(verbose){printf("Atlas is dimension %d, with %d charts.\n",k,nCharts);fflush(stdout);}
#endif

  MFCreateEmptyKComplex(k,MFAtlasN(A,e),&Complex,0,e);

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("In MFEnumeratePolytopes\n");fflush(stdout);
    for(chart=0;chart<nCharts;chart++)
     {
      P=MFChartPolytope(MFAtlasChart(A,chart,e),e);
      printf("\nPolytope %d\n\n",chart);
      MFPrintPolytopeTerse(stdout,P,e);
     }
    printf("\nFaces\n\n");
    MFPrintAtlasFaceList(stdout,A,e);
    fflush(stdout);
   }
#endif

  for(chart=0;chart<nCharts;chart++)
   {

#ifdef MFALLOWVERBOSE
    if(verbose){printf("chart %d/%d\n",chart,nCharts);fflush(stdout);}
#endif

    if(MFAtlasChart(A,chart,e)==NULL)
     {
      {printf("chart %d/%d is NULL! ,atlas 0x%8.8x\n",chart,nCharts,A);fflush(stdout);}
       continue;
     }
    P=MFChartPolytope(MFAtlasChart(A,chart,e),e);

/*  for(face=0;face<MFPolytopeNumberOfFaces(P,e);face++)
     {
      face1=MFPolytopeFaceIndex(P,face);
     Merge all vertices on a face if the face is not also a face of
     the polytope across the face 
     } */

/*  MFPolytopeMergeCloseVertices(P,.05*MFAtlasChartRadius(A,chart,e),e);
    MFPolytopeUpdateFaceList(P,e);*/

#ifdef MFALLOWVERBOSE
    if(verbose){printf("call MFEnumeratePolytope, chart %d\n",chart);fflush(stdout);}
#endif

    if(MFPolytopeNumberOfVertices(P,e)>0)
      EPS[chart]=MFEnumeratePolytope(P,e);
     else
      EPS[chart]=NULL;

#ifdef MFALLOWVERBOSE
    if(verbose){printf("call MFTestEnumPolytope, chart %d\n",chart);fflush(stdout);}
#endif

    if(0 && EPS[chart]!=NULL && (failed=!MFTestEnumPolytope(EPS[chart],1,e)))
     {
      printf("----------------------------------------------------\n");
      printf("----------------------------------------------------\n");
      printf("\n%s EPS[%d] failed tests:\n\n",RoutineName,chart);
      printf("P[%d]=0x%8.8x\n",chart,P,e);
      MFPrintPolytope(stdout,P,e);
      printf("\n");
      printf("EPS[%d]=0x%8.8x\n",chart,EPS[chart]);
      MFPrintEnumPolytope(stdout,EPS[chart],e);
      printf("\n");
      printf("----------------------------------------------------\n");
      printf("----------------------------------------------------\n");
      fflush(stdout);
      abort();
     }

#ifdef MFALLOWVERBOSE
    if(verbose)
     {
      printf("  Chart %d:\n",chart);
      printf("EPS[%d]=0x%8.8x\n",chart,EPS[chart]);
      printf("\nCenter:");fflush(stdout);
      if(MFChartCenter(MFAtlasChart(A,chart,e),e)==NULL)
       {
        printf("NULL\n");fflush(stdout);abort();
       }else
        MFPrintNVector(stdout,MFChartCenter(MFAtlasChart(A,chart,e),e),e);
      fflush(stdout);
      printf("\n");fflush(stdout);
      printf("\nPolytope:\n\n");fflush(stdout);
      MFPrintPolytope(stdout,MFChartPolytope(MFAtlasChart(A,chart,e),e),e);fflush(stdout);
      printf("\n");fflush(stdout);
      printf("\nEnumPolytope:\n\n");fflush(stdout);
      if(EPS[chart]!=NULL){MFPrintEnumPolytope(stdout,EPS[chart],e);fflush(stdout);}
      printf("\n");
      fflush(stdout);
     }
#endif

/*  if(failed)abort();*/

#ifdef MFALLOWVERBOSE
    if(verbose){printf("In middle of MFEnumeratePolytopes\n");fflush(stdout);}
#endif

    if(EPS[chart]!=NULL)
     {
    for(dim=0;dim<k+1;dim++)
     {
      nCells=MFEnumPolytopeNumberOfCells(EPS[chart],dim,e);

#ifdef MFALLOWVERBOSE
      if(verbose){printf("    There are %d cells of dimension %d\n",nCells,dim);fflush(stdout);}
#endif

      for(cell=0;cell<nCells;cell++)
       {
        mfcell=((EPS[chart])->cells[dim])[cell];

#ifdef MFALLOWVERBOSE
        if(verbose)
         {
          printf("      %d [",cell);
          for(i=0;i<mfcell->nIndices;i++)
           {
            if(i>0)printf(",");
            printf("%d",mfcell->indices[i]);
           }
          printf("]");
          fflush(stdout);
         }
#endif

        if(!MFEnumPolytopeIsCellExterior(A,chart,EPS[chart],dim,cell,e))
         {
          indices[0]=chart;
          m=1;
          for(i=0;i<mfcell->nIndices;i++)
           {
            if((indices[m]=MFAtlasLeftPolytope(A,mfcell->indices[i],e))==chart)
                indices[m]=MFAtlasRightPolytope(A,mfcell->indices[i],e);
            m++;
           }
          EPSortPolytopeIndices(m,indices,e);
          mfcell->indices=(int*)realloc((void*)mfcell->indices,(m+1)*sizeof(int));

#ifndef MFNOSAFETYNET
          if(mfcell->indices==NULL)
           {
            sprintf(MFEnumPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",(m+1)*sizeof(int));
            MFSetError(e,12,RoutineName,MFEnumPolytopeErrorMsg,__LINE__,__FILE__);
            MFErrorHandlerOutOfMemory(e);
            return NULL;
           }
#endif

          mfcell->nIndices=m;
          for(i=0;i<m;i++)mfcell->indices[i]=indices[i];

#ifdef MFALLOWVERBOSE
          if(verbose)
           {
            printf("        Polytopes [");
            for(i=0;i<m;i++)
             {
              if(i>0)printf(",");
              printf("%d",mfcell->indices[i]);
             }
            printf("]\n");
            fflush(stdout);
           }
#endif
          mfcell->indices[m]=-1;
         }else{
          mfcell->nIndices=0;
          if(mfcell->indices!=NULL)
           {
            free(mfcell->indices);
           }
          mfcell->indices=NULL;

#ifdef MFALLOWVERBOSE
          if(verbose){printf("        Not a real %d-cell\n",dim);fflush(stdout);}
#endif

         }
       }

#ifdef MFALLOWVERBOSE
      if(verbose){printf("done print dim %d in middle of MFEnumeratePolytopes\n",dim);fflush(stdout);}
#endif

     }
#ifdef MFALLOWVERBOSE
      if(verbose){printf("done print in middle of MFEnumeratePolytopes\n");fflush(stdout);}
#endif
     }
   }

/* Now have for each cell of each polytope a list of the k-cells to which it belongs */

#ifdef MFALLOWVERBOSE
  if(verbose){printf("\nSecond step: assign global numbers to the cells in the polytope of each chart.\n\n");fflush(stdout);}
#endif

  for(i=0;i<k+1;i++)Complex->nCells[i]=0;

  for(chart=0;chart<nCharts;chart++)
   {

#ifdef MFALLOWVERBOSE
    if(verbose){printf("  Chart %d - \n",chart);fflush(stdout);}
#endif

    if(EPS[chart]!=NULL)
     {
    for(dim=0;dim<k+1;dim++)
     {
      nCells=MFEnumPolytopeNumberOfCells(EPS[chart],dim,e);

#ifdef MFALLOWVERBOSE
      if(verbose){printf("    There are %d cells of dimension %d \n",nCells,dim);fflush(stdout);}
#endif

      found=0;
      for(cell=0;cell<nCells;cell++)
       {
        mfcell=((EPS[chart])->cells[dim])[cell];
        if(mfcell->nIndices>0 && mfcell->indices[mfcell->nIndices]==-1)
         {
          found=1;

#ifdef MFALLOWVERBOSE
          if(verbose)
           {
            printf("     %d-cell %d is unassigned, %d, signature [",dim,cell,mfcell->indices[mfcell->nIndices]);
            for(i=0;i<mfcell->nIndices;i++)
             {
              if(i>0)printf(",");
              printf("%d",mfcell->indices[i]);
             }
            printf("]\n");
            fflush(stdout);
           }
#endif

          for(i=0;i<mfcell->nIndices;i++)
           {
            chart2=mfcell->indices[i];

#ifdef MFALLOWVERBOSE
            if(verbose){printf("       Polytope %d shares this cell\n",chart2);fflush(stdout);}
#endif

            if(chart2!=chart && EPS[chart2]!=NULL)
             {
              nCells2=MFEnumPolytopeNumberOfCells(EPS[chart2],dim,e);
              maxSame=0;
              mSame=0;
              mfcellClosest=NULL;
              cellClosest=-1;
              for(cell2=0;cell2<nCells2;cell2++)
               {
                mfcell2=((EPS[chart2])->cells[dim])[cell2];
                if(mfcell2->nIndices==0)
                 {

#ifdef MFALLOWVERBOSE
                  if(verbose){printf("         cell %d is not a real cell\n",cell2);fflush(stdout);}
#endif
                 }else{
                  nSame=MFIntersectIndices(mfcell->nIndices,mfcell->indices,mfcell2->nIndices,mfcell2->indices,NULL,e);

#ifdef MFALLOWVERBOSE
                  if(verbose)
                   {
                    printf("         cell %d, %d, signature [",cell2,mfcell2->indices[mfcell2->nIndices]);
                    for(l=0;l<mfcell->nIndices;l++)
                     {
                      if(l>0)printf(",");
                      printf("%d",mfcell2->indices[l]);
                     }
                    printf("] matches %d, max so far is %d\n",nSame,maxSame);
                    fflush(stdout);
                   }
#endif

/* Here don't want to require match. Know that one of these vertices is the one */
                  if(nSame>maxSame)
                   {
                    maxSame=nSame;
                    mSame=1;
                    mfcellClosest=mfcell2;
                    cellClosest=cell2;
                   }else if(nSame==maxSame) mSame++;
                 }
               }

#ifdef MFALLOWVERBOSE
              if(verbose){printf("       blah %d \n",chart2);fflush(stdout);}
#endif

/*            if(mfcellClosest!=NULL && mSame==1 )*/
              if(mfcellClosest!=NULL && maxSame==mfcell->nIndices )
               {
                mfcell2=mfcellClosest;

#ifdef MFALLOWVERBOSE
                if(verbose)
                 {
                  printf("         match! with cell %d, chart %d signature [",cellClosest,chart2);
                  for(l=0;l<mfcell2->nIndices;l++)
                   {
                    if(l>0)printf(",");
                    printf("%d",mfcell2->indices[l]);
                   }
                  printf("]");
                  printf(" %d matches with chart %d cell %d [",maxSame,chart,cell);
                  for(l=0;l<mfcell->nIndices;l++)
                   {
                    if(l>0)printf(",");
                    printf("%d",mfcell->indices[l]);
                   }
                  printf("]\n");
                  fflush(stdout);
                 }
#endif
                if(mfcell->indices[mfcell->nIndices]==-1
                && mfcell2->indices[mfcell2->nIndices]==-1)
                 {
                  mfcell->indices[mfcell->nIndices]=Complex->nCells[dim];

#ifdef MFALLOWVERBOSE
                  if(verbose)printf("     unassigned %d-cell %d on chart %d will be assigned index %d\n",dim,cell,chart,mfcell->indices[mfcell->nIndices]);
#endif

                  Complex->nCells[dim]++;
                 }
                if(mfcell2->indices[mfcell2->nIndices]==-1)
                 {
                  mfcell2->indices[mfcell2->nIndices]=mfcell->indices[mfcell->nIndices];

#ifdef MFALLOWVERBOSE
                  if(verbose)printf("     matching %d-cell %d on chart %d will also be assigned index %d\n",dim,cellClosest,chart2,mfcell2->indices[mfcell2->nIndices]);
#endif
                 }else{
                  mfcell->indices[mfcell->nIndices]=mfcell2->indices[mfcell2->nIndices];

#ifdef MFALLOWVERBOSE
                  if(verbose)printf("     unassigned %d-cell %d on chart %d will be assigned index %d from matching cell\n",dim,cell,chart,mfcell2->indices[mfcell2->nIndices]);
#endif
                 }
               }
             } /* if chart!=chart2 */

#ifdef MFALLOWVERBOSE
            if(verbose){printf("       blah2 %d \n",chart2);fflush(stdout);}
#endif
           } /* for Indices */
          if(mfcell->indices[mfcell->nIndices]==-1 && dim==k)
           {
            mfcell->indices[mfcell->nIndices]=Complex->nCells[dim];

#ifdef MFALLOWVERBOSE
            if(verbose)printf("     unassigned %d-cell %d on chart %d matches nothing. Will be assigned index %d\n",dim,cell,chart,mfcell->indices[mfcell->nIndices]);
#endif
            Complex->nCells[dim]++;
           }
         }   /* if unassigned */
       }     /* for cell */

#ifdef MFALLOWVERBOSE
      if(verbose&&!found){printf("        None become new ones\n");fflush(stdout);}
#endif
     } /* for dim */
     }
   }   /* for chart */

/* At this point we have found all the vertices which are simple. */

/* Now have for each cell of each polytope a global index, and a count of the cells of each dimension */
/* May have unmatched cells (index -1) */

/* If k==2 can walk around 0-cells. 
                           1-cells are already identified, 
                           2-cells are unique */

/* If k>2 can walk around (k-2)-cells. 
                          (k-1)-cells are already identified, 
                           k-cells are unique */

#ifdef MFALLOWVERBOSE
  if(verbose){printf("\nThird step: Allocate space for cells in the Complex.\n\n");fflush(stdout);}
#endif

  for(dim=0;dim<k+1;dim++)
   {
    if(dim==0)
     {
      Complex->mCells[dim]=Complex->nCells[dim];
      Complex->cells[dim]=(MFCell*)realloc((void*)(Complex->cells[dim]),Complex->mCells[dim]*sizeof(MFCell));

#ifndef MFNOSAFETYNET
      if(Complex->cells[dim]==NULL && Complex->mCells[dim]!=0)
       {
        sprintf(MFEnumPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",Complex->mCells[dim]*sizeof(MFCell));
        MFSetError(e,12,RoutineName,MFEnumPolytopeErrorMsg,__LINE__,__FILE__);
        MFErrorHandlerOutOfMemory(e);
        return NULL;
       }
#endif

      Complex->v=(MFKVector*)realloc(Complex->v,Complex->mCells[dim]*sizeof(MFKVector));

#ifndef MFNOSAFETYNET
      if(Complex->v==NULL && Complex->mCells[dim]!=0)
       {
        sprintf(MFEnumPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",Complex->mCells[dim]*sizeof(MFKVector));
        MFSetError(e,12,RoutineName,MFEnumPolytopeErrorMsg,__LINE__,__FILE__);
        MFErrorHandlerOutOfMemory(e);
        return NULL;
       }
#endif

      for(i=0;i<Complex->mCells[dim];i++)Complex->v[i]=NULL;
     }else{
      Complex->mCells[dim]=Complex->nCells[dim];
      Complex->cells[dim]=(MFCell*)realloc((void*)Complex->cells[dim],Complex->mCells[dim]*sizeof(MFCell));

#ifndef MFNOSAFETYNET
      if(Complex->cells[dim]==NULL && Complex->mCells[dim]!=0)
       {
        sprintf(MFEnumPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",Complex->mCells[dim]*sizeof(MFCell));
        MFSetError(e,12,RoutineName,MFEnumPolytopeErrorMsg,__LINE__,__FILE__);
        MFErrorHandlerOutOfMemory(e);
        return NULL;
       }
#endif

     }
#ifdef MFALLOWVERBOSE
    if(verbose){printf("    There are %d cells of dimension %d\n",Complex->nCells[dim],dim);fflush(stdout);}
#endif

    for(i=0;i<Complex->nCells[dim];i++)
     {
      mfcell=(MFCell)malloc(sizeof(struct MFCellSt));

#ifndef MFNOSAFETYNET
      if(mfcell==NULL)
       {
        sprintf(MFEnumPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct MFCellSt));
        MFSetError(e,12,RoutineName,MFEnumPolytopeErrorMsg,__LINE__,__FILE__);
        MFErrorHandlerOutOfMemory(e);
        return NULL;
       }
#endif

      (Complex->cells[dim])[i]=mfcell;
      mfcell->nIndices=0;
      mfcell->indices=NULL;
      mfcell->nFaceCells=0;
      mfcell->mFaceCells=0;
      mfcell->faceCells=NULL;
      mfcell->nCellFaces=0;
      mfcell->mCellFaces=0;
      mfcell->cellFaces=NULL;
     }
   }

#ifdef MFALLOWVERBOSE
  if(verbose){printf("\nFourth step: Place the k-cells from the polytopes into the Complex.\n\n");fflush(stdout);}
#endif

  for(chart=0;chart<nCharts;chart++)
   {
    if(EPS[chart]!=NULL)
     {

#ifdef MFALLOWVERBOSE
    if(verboseTerse)
     {
      printf("\n  Chart %d:\n\n",chart);
      MFPrintEnumPolytope(stdout,EPS[chart],e);
     }
    if(verbose){printf("    %d-cell %d\n",k,chart);fflush(stdout);}
#endif

    mfcell2=((EPS[chart])->cells[k])[0];
    if(Complex->cells[k]!=NULL)
      mfcell=(Complex->cells[k])[chart];
     else
      mfcell=NULL;

    if(chart<Complex->nCells[k] && mfcell!=NULL)
     {

#ifdef MFALLOWVERBOSE
      if(verbose){printf("        mfcell is (Complex->cells[%d])[%d] = 0x%8.8x\n",k,chart,mfcell);fflush(stdout);}
#endif

      mfcell->nFaceCells=0;
      mfcell->mFaceCells=0;
      mfcell->faceCells=NULL;
  
      m=0;
      for(i=0;i<mfcell2->nCellFaces;i++)
       {
        cellFace=(EPS[chart]->cells[k-1])[mfcell2->cellFaces[i]];
        if(cellFace->nIndices>0 && cellFace->indices[cellFace->nIndices]>-1 )m++;
       }
      if(m==0)
       {
        mfcell->nCellFaces=0;
        mfcell->mCellFaces=0;
        mfcell->cellFaces=NULL;

#ifdef MFALLOWVERBOSE
        if(verbose){printf("           has no cell faces\n");fflush(stdout);}
#endif

       }else{
        mfcell->nCellFaces=0;
        mfcell->mCellFaces=m;
        mfcell->cellFaces=(int*)malloc(mfcell->mCellFaces*sizeof(int));

#ifndef MFNOSAFETYNET
        if(mfcell->cellFaces==NULL)
         {
          sprintf(MFEnumPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",mfcell->mCellFaces*sizeof(int));
          MFSetError(e,12,RoutineName,MFEnumPolytopeErrorMsg,__LINE__,__FILE__);
          MFErrorHandlerOutOfMemory(e);
          return NULL;
         }
#endif
  
        for(i=0;i<mfcell2->nCellFaces;i++)
         {
          cellFace=(EPS[chart]->cells[k-1])[mfcell2->cellFaces[i]];
          if(cellFace->nIndices>0 && cellFace->indices[cellFace->nIndices]>-1)
           {
            mfcell->cellFaces[mfcell->nCellFaces]=cellFace->indices[cellFace->nIndices];
            (mfcell->nCellFaces)++;
           }
         }

#ifdef MFALLOWVERBOSE
        if(verbose)
         {
          printf("        will have %d cell faces [",mfcell->nCellFaces);
          for(i=0;i<mfcell->nCellFaces;i++)
           {
            if(i>0)printf(",");
            printf("%d",mfcell->cellFaces[i]);
           }
          printf("]\n");
         }
#endif

       }
     }
     }
   }

#ifdef MFALLOWVERBOSE
  if(verbose){printf("\nFifth step: Merge the cell lists from the polytopes into the Complex.\n\n");fflush(stdout);}
#endif

  for(chart=0;chart<nCharts;chart++)
   {

#ifdef MFALLOWVERBOSE
    if(verbose){printf("    chart %d\n",chart);fflush(stdout);}
#endif

    if(EPS[chart]!=NULL)
     {
    for(dim=0;dim<k;dim++)
     {
      nCells=MFEnumPolytopeNumberOfCells(EPS[chart],dim,e);

#ifdef MFALLOWVERBOSE
      if(verbose){printf("      There are %d cells of dimension %d\n",nCells,dim);fflush(stdout);}
#endif

      for(cell=0;cell<nCells;cell++)
       {
        mfcell=(EPS[chart]->cells[dim])[cell];
        if(mfcell->nIndices>0 && mfcell->indices[mfcell->nIndices]!=-1)
         {
          i=mfcell->indices[mfcell->nIndices];
          mfcell2=Complex->cells[dim][i];

          m=0;
          for(i=0;i<mfcell->nCellFaces;i++)
           {
            cellFace=(EPS[chart]->cells[dim-1])[mfcell->cellFaces[i]];
            if(cellFace->nIndices>0 && cellFace->indices[cellFace->nIndices]>-1 )m++;
           }

#ifdef MFALLOWVERBOSE
          if(verbose){printf("        cell %d is global %d-cell %d\n",cell,dim,mfcell->indices[mfcell->nIndices]);fflush(stdout);}
#endif

          if(m==0)
           {

#ifdef MFALLOWVERBOSE
            if(verbose){printf("          has no new cellFaces to contribute\n");fflush(stdout);}
#endif

           }else{

#ifdef MFALLOWVERBOSE
            if(verbose)
             {
              printf("          has %d new cellFaces to contribute [",m);
              for(i=0;i<mfcell->nCellFaces;i++)
               {
                cellFace=((EPS[chart])->cells[dim-1])[mfcell->cellFaces[i]];
                if(i>0)printf(",");
                if(cellFace->nIndices>0 && cellFace->indices[cellFace->nIndices]>-1)
                  printf("%d",cellFace->indices[cellFace->nIndices]);
                 else
                  printf("--");
               }
              printf("]\n");
              fflush(stdout);
             }
#endif

            if(mfcell2->mCellFaces==0)
             {
              mfcell2->mCellFaces=m;
              if(mfcell2->mCellFaces>-1)
               {
                mfcell2->cellFaces=(int*)malloc(mfcell2->mCellFaces*sizeof(int));

#ifndef MFNOSAFETYNET
                if(mfcell2->cellFaces==NULL)
                 {
                  sprintf(MFEnumPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",mfcell->mCellFaces*sizeof(int));
                  MFSetError(e,12,RoutineName,MFEnumPolytopeErrorMsg,__LINE__,__FILE__);
                  MFErrorHandlerOutOfMemory(e);
                  return NULL;
                 }
#endif

               }
             }else{
              mfcell2->mCellFaces+=m;
              if(mfcell2->mCellFaces>-1)
               {
                mfcell2->cellFaces=(int*)realloc((void*)(mfcell2->cellFaces),mfcell2->mCellFaces*sizeof(int));

#ifndef MFNOSAFETYNET
                if(mfcell2->cellFaces==NULL)
                 {
                  sprintf(MFEnumPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",mfcell2->mCellFaces*sizeof(int));
                  MFSetError(e,12,RoutineName,MFEnumPolytopeErrorMsg,__LINE__,__FILE__);
                  MFErrorHandlerOutOfMemory(e);
                  return NULL;
                 }
#endif

               }
             }
    
            for(i=0;i<mfcell->nCellFaces;i++)
             {
              cellFace=((EPS[chart])->cells[dim-1])[mfcell->cellFaces[i]];
              if(cellFace->nIndices>0 && cellFace->indices[cellFace->nIndices]>-1)
               {
                mfcell2->cellFaces[mfcell2->nCellFaces]=cellFace->indices[cellFace->nIndices];

#ifdef MFALLOWVERBOSE
                if(verbose){printf("             add cellFace %d\n",mfcell2->cellFaces[mfcell2->nCellFaces]);fflush(stdout);}
#endif

                mfcell2->nCellFaces++;
               }
             }

#ifdef MFALLOWVERBOSE
            if(verbose)
             {
              printf("        before collapse: cell faces [");
              for(i=0;i<mfcell2->nCellFaces;i++)
               {
                if(i>0)printf(",");
                printf("%d",mfcell2->cellFaces[i]);
               }
              printf("]\n");
              fflush(stdout);
             }
#endif

            EPCollapseIndices(&(mfcell2->nCellFaces),mfcell2->cellFaces,e);

#ifdef MFALLOWVERBOSE
            if(verbose)
             {
              printf("        after collapse: cell faces [");
              for(i=0;i<mfcell2->nCellFaces;i++)
               {
                if(i>0)printf(",");
                printf("%d",mfcell2->cellFaces[i]);
               }
              printf("]\n");
              fflush(stdout);
             }
#endif
           }

          m=0;
          for(i=0;i<mfcell->nFaceCells;i++)
           {
            faceCell=(EPS[chart]->cells[dim+1])[mfcell->faceCells[i]];
            if(faceCell->nIndices>0 && faceCell->indices[faceCell->nIndices]>-1)m++;
           }

          if(m==0)
           {

#ifdef MFALLOWVERBOSE
            if(verbose){printf("        has no new faceCells to contribute\n");fflush(stdout);}
#endif

           }else{

#ifdef MFALLOWVERBOSE
            if(verbose)
             {
              printf("          has %d new faceCells to contribute [",m);
              for(i=0;i<mfcell->nFaceCells;i++)
               {
                faceCell=(EPS[chart]->cells[dim+1])[mfcell->faceCells[i]];
                if(i>0)printf(",");
                if(faceCell->nIndices>0 && faceCell->indices[faceCell->nIndices]>-1)
                  printf("%d",faceCell->indices[faceCell->nIndices]);
                 else
                  printf("--");
               }
              printf("]\n");
              fflush(stdout);
             }
#endif

            if(mfcell2->mFaceCells==0)
             {
              mfcell2->mFaceCells=m+1;
              if(mfcell2->mFaceCells>-1)
               {
                mfcell2->faceCells=(int*)malloc(mfcell2->mFaceCells*sizeof(int));

#ifndef MFNOSAFETYNET
                if(mfcell2->faceCells==NULL)
                 {
                  sprintf(MFEnumPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",mfcell2->mFaceCells*sizeof(int));
                  MFSetError(e,12,RoutineName,MFEnumPolytopeErrorMsg,__LINE__,__FILE__);
                  MFErrorHandlerOutOfMemory(e);
                  return NULL;
                 }
#endif

               }
             }else{
              mfcell2->mFaceCells+=m+1;
              if(mfcell2->mFaceCells>-1)
               {
                mfcell2->faceCells=(int*)realloc((void*)(mfcell2->faceCells),mfcell2->mFaceCells*sizeof(int));

#ifndef MFNOSAFETYNET
                if(mfcell2->faceCells==NULL)
                 {
                  sprintf(MFEnumPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",mfcell2->mFaceCells*sizeof(int));
                  MFSetError(e,12,RoutineName,MFEnumPolytopeErrorMsg,__LINE__,__FILE__);
                  MFErrorHandlerOutOfMemory(e);
                  return NULL;
                 }
#endif

               }
             }
  
            for(i=0;i<mfcell->nFaceCells;i++)
             {
              faceCell=((EPS[chart])->cells[dim+1])[mfcell->faceCells[i]];
              if(faceCell->nIndices>0 && faceCell->indices[faceCell->nIndices]>-1)
               {
                mfcell2->faceCells[mfcell2->nFaceCells]=faceCell->indices[faceCell->nIndices];

#ifdef MFALLOWVERBOSE
                if(verbose){printf("             add faceCell %d\n",mfcell2->faceCells[mfcell2->nFaceCells]);fflush(stdout);}
#endif

                mfcell2->nFaceCells++;
               }
             }

#ifdef MFALLOWVERBOSE
            if(verbose)
             {
              printf("        before collapse: face cells [");
              for(i=0;i<mfcell2->nFaceCells;i++)
               {
                if(i>0)printf(",");
                printf("%d",mfcell2->faceCells[i]);
               }
              printf("]\n");
              fflush(stdout);
             }
#endif

            EPCollapseIndices(&(mfcell2->nFaceCells),mfcell2->faceCells,e);

#ifdef MFALLOWVERBOSE
            if(verbose)
             {
              printf("        after collapse: face cells [");
              for(i=0;i<mfcell2->nFaceCells;i++)
               {
                if(i>0)printf(",");
                printf("%d",mfcell2->faceCells[i]);
               }
              printf("]\n");
              fflush(stdout);
             }
#endif
           }
         }
       }
     }
     }  /*if(EPS[chart]!=NULL)*/
   }

  for(chart=0;chart<nCharts;chart++)
   if(EPS[chart]!=NULL)MFFreeEnumPolytope(EPS[chart],e);
  free(EPS);

#ifdef MFALLOWVERBOSE
  if((failed=!MFTestEnumPolytope(Complex,0,e)) || verboseResult)
   {
    printf("----------------------------------------------------\n");
    printf("----------------------------------------------------\n");
    printf("\nComplex:\n\n");
    MFPrintEnumPolytope(stdout,Complex,e);
    printf("\n");
    printf("----------------------------------------------------\n");
    printf("----------------------------------------------------\n");
    fflush(stdout);
   }

  if(verboseTerse)
   {
    int halfspace;

    printf("\nHalf Spaces:\n\n");
    for(halfspace=0;halfspace<MFAtlasNumberOfHalfSpaces(A,e);halfspace++)
     {
      printf("      %d       [%d,%d]\n",halfspace,MFAtlasLeftPolytope(A,halfspace,e),MFAtlasRightPolytope(A,halfspace,e));
     }
    fflush(stdout);
   }
#endif

  return Complex;
 }

MFKVector MFEnumPolytopeVertex(MFEnumPolytope P,long i, MFErrorHandler e)
 {
  static char RoutineName[]={"MFEnumPolytopeVertex"};
  return P->v[i];

#ifdef __cplusplus
}
#endif
 }
