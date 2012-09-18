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

static char *id="@(#) $Id: MFEnumDualPolytope.c,v 1.3 2007/02/13 01:22:33 mhender Exp $";

static char MFEnumDualPolytopeErrorMsg[256]="";

#include <MFAtlas.h>
#include <MFEnumDualPolytope.h>
#include <MFEnumPolytope.h>
#include <MFPolytope.h>
#include <MFAtlas.h>
#include <MFChart.h>
#include <MFNVector.h>
#include <MFPrint.h>
#include <MFPrintEnum.h>
#include <stdio.h>
#include <stdlib.h>
#include <MFErrorHandler.h>

#ifdef __cplusplus
 extern "C" {
#endif

struct MFDualCellSt {
             int nIndices;
             int *indices;
             int mFaceCells;
             int nFaceCells;
             int *faceCells;
             int nCellFaces;
             int mCellFaces;
             int *cellFaces;
            };

typedef struct MFDualCellSt *MFDualCell;

struct MFEnumDualPolytopeSt {
                  MFNVector *v;
                  int n;
                  int *nCells;
                  int *mCells;
                  MFDualCell **cells;
                 };
#define EPNumberOfCellsToAllocate 1

MFChart MFAtlasChart(MFAtlas,int,MFErrorHandler);

int MFDualIntersectIndices(int,int*,int,int*,int*,MFErrorHandler);
void MFCreateDualZeroCellV(MFDualCell*,MFErrorHandler,int,...);
void MFCreateDualZeroCell(MFDualCell*,int,int*,MFErrorHandler);
void MFCreateDualZeroCellWithEdgeList(MFDualCell*,int,int*,MFErrorHandler);
void MFFreeDualCell(MFDualCell,MFErrorHandler);
void MFCreateEmptyDualNComplex(int,MFEnumDualPolytope*,int,MFErrorHandler);
int MFDualIndicesEqual(int,int*,int,int*,MFErrorHandler);
void EPSortDualPolytopeIndices(int,int*,MFErrorHandler);
void EPDualCollapseIndices(int*,int*,MFErrorHandler);
int MFEnumDualPolytopeVerticesOnDualCell(MFEnumDualPolytope,int,int,int,int**,MFErrorHandler);
int MFEnumDualPolytopeIsCellExterior(MFAtlas,int,MFEnumDualPolytope,int,int,MFErrorHandler);

MFEnumDualPolytope MFCreateDualOfEnumPolytope(MFEnumPolytope,MFErrorHandler);
int MFTestEnumDualPolytope(MFEnumDualPolytope,int,MFErrorHandler);

struct MFCellSt;

void MFEnumerateDual(MFEnumDualPolytope EP, MFErrorHandler e)
 {
  static char RoutineName[]={"MFEnumerateDual"};
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
  int    verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose)fprintf(stdout,"MFEnumerateDual\n");
#endif

/* Pairwise intersection of vertices */

  k=EP->n;

  for(d=1;d<k+1;d++)
   {

#ifdef MFALLOWVERBOSE
    if(verbose)fprintf(stdout,"  intersecting %d cells to build %d cells\n",d-1,d);
#endif

    nv=EP->nCells[d-1];

#ifdef MFALLOWVERBOSE
    if(verbose)fprintf(stdout,"  There are %d %d-cells\n",nv,d-1);
#endif

    if(EP->mCells[d]==0)
     {
      EP->mCells[d]=EPNumberOfCellsToAllocate;
      EP->cells[d]=(MFDualCell*)realloc((void*)(EP->cells[d]),EP->mCells[d]*sizeof(MFDualCell));

#ifndef MFNOSAFETYNET
      if(EP->cells[d]==NULL)
       {
        sprintf(MFEnumDualPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",EP->mCells[d]*sizeof(MFDualCell));
        MFSetError(e,12,RoutineName,MFEnumDualPolytopeErrorMsg,__LINE__,__FILE__);
        MFErrorHandlerOutOfMemory(e);
        return;
       }
#endif

     }

#ifdef MFALLOWVERBOSE
    if(verbose)fprintf(stdout,"  There is space for %d %d-cells\n",EP->mCells[d],d);
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
        ni=MFDualIntersectIndices(n1,i1,n2,i2,intersection,e);

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
          alreadyThere=0;
          n=-1;
          for(i=0;i<ne;i++)
           {
            n2=((EP->cells[d])[i])->nIndices;
            i2=((EP->cells[d])[i])->indices;
            if(MFDualIndicesEqual(ni,intersection,n2,i2,e)){alreadyThere=1;n=i;}
           }
          if(!alreadyThere)
           {

#ifdef MFALLOWVERBOSE
            if(verbose)fprintf(stdout,"  Found a new %d-cell\n",d);
#endif
            if(EP->mCells[d]==0)
             {
              EP->mCells[d]=EPNumberOfCellsToAllocate;
              EP->cells[d]=(MFDualCell*)realloc((void*)(EP->cells[d]),EP->mCells[d]*sizeof(MFDualCell));

#ifndef MFNOSAFETYNET
              if(EP->cells[d]==NULL)
               {
                sprintf(MFEnumDualPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",EP->mCells[d]*sizeof(MFDualCell));
                MFSetError(e,12,RoutineName,MFEnumDualPolytopeErrorMsg,__LINE__,__FILE__);
                MFErrorHandlerOutOfMemory(e);
                return;
               }
#endif
             }else if(ne>=EP->mCells[d])
             {
              EP->mCells[d]+=EPNumberOfCellsToAllocate;
              EP->cells[d]=(MFDualCell*)realloc((void*)(EP->cells[d]),EP->mCells[d]*sizeof(MFDualCell));

#ifndef MFNOSAFETYNET
              if(EP->cells[d]==NULL)
               {
                sprintf(MFEnumDualPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",EP->mCells[d]*sizeof(MFDualCell));
                MFSetError(e,12,RoutineName,MFEnumDualPolytopeErrorMsg,__LINE__,__FILE__);
                MFErrorHandlerOutOfMemory(e);
                return;
               }
#endif
             }

#ifdef MFALLOWVERBOSE
            if(verbose)fprintf(stdout,"  This will be %d-cell number %d\n",d,ne);
#endif
            (EP->cells[d])[ne]=(MFDualCell)malloc(sizeof(struct MFDualCellSt));

#ifndef MFNOSAFETYNET
            if(EP->cells[d][ne]==NULL)
             {
              sprintf(MFEnumDualPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct MFDualCellSt));
              MFSetError(e,12,RoutineName,MFEnumDualPolytopeErrorMsg,__LINE__,__FILE__);
              MFErrorHandlerOutOfMemory(e);
              return;
             }
#endif

            ((EP->cells[d])[ne])->nIndices=ni;
            if(ni>0)
             {
              ((EP->cells[d])[ne])->indices=(int*)malloc(ni*sizeof(int));

#ifndef MFNOSAFETYNET
              if((EP->cells[d][ne])->indices==NULL)
               {
                sprintf(MFEnumDualPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",ni*sizeof(int));
                MFSetError(e,12,RoutineName,MFEnumDualPolytopeErrorMsg,__LINE__,__FILE__);
                MFErrorHandlerOutOfMemory(e);
                return;
               }
#endif
              for(i=0;i<ni;i++)
                ((EP->cells[d])[ne])->indices[i]=intersection[i];
             }else
              ((EP->cells[d])[ne])->indices=NULL;

            ((EP->cells[d])[ne])->mCellFaces=2;
            ((EP->cells[d])[ne])->nCellFaces=2;
            ((EP->cells[d])[ne])->cellFaces=(int*)malloc(((EP->cells[d])[ne])->mCellFaces*sizeof(int));

#ifndef MFNOSAFETYNET
            if((EP->cells[d][ne])->cellFaces==NULL)
             {
              sprintf(MFEnumDualPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",((EP->cells[d])[ne])->mCellFaces*sizeof(int));
              MFSetError(e,12,RoutineName,MFEnumDualPolytopeErrorMsg,__LINE__,__FILE__);
              MFErrorHandlerOutOfMemory(e);
              return;
             }
#endif
            ((EP->cells[d])[ne])->cellFaces[0]=iv;
            ((EP->cells[d])[ne])->cellFaces[1]=jv;

            ((EP->cells[d])[ne])->mFaceCells=0;
            ((EP->cells[d])[ne])->nFaceCells=0;
            ((EP->cells[d])[ne])->faceCells=NULL;

            n=ne;
            ne++;
           }else{
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
                  sprintf(MFEnumDualPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",((EP->cells[d])[n])->mCellFaces*sizeof(int));
                  MFSetError(e,12,RoutineName,MFEnumDualPolytopeErrorMsg,__LINE__,__FILE__);
                  MFErrorHandlerOutOfMemory(e);
                  return;
                 }
#endif
               }
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
                  sprintf(MFEnumDualPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",((EP->cells[d])[n])->mCellFaces*sizeof(int));
                  MFSetError(e,12,RoutineName,MFEnumDualPolytopeErrorMsg,__LINE__,__FILE__);
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
                    sprintf(MFEnumDualPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",((EP->cells[d-1])[f])->mFaceCells*sizeof(int));
                    MFSetError(e,12,RoutineName,MFEnumDualPolytopeErrorMsg,__LINE__,__FILE__);
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
                  sprintf(MFEnumDualPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",((EP->cells[d-1])[f])->mFaceCells*sizeof(int));
                  MFSetError(e,12,RoutineName,MFEnumDualPolytopeErrorMsg,__LINE__,__FILE__);
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
   }

#ifdef MFALLOWVERBOSE
  if(verbose)MFPrintEnumDualPolytope(stdout,EP,e);
#endif

  return;
 }

int MFDualIntersectIndices(int n1,int *inter1,int n2,int *inter2,int *inter, MFErrorHandler e)
 {
  static char RoutineName[]={"MFDualIntersectIndices"};
  int n;
  int i1,i2;

/*if(!MFTestVertexIndexOrder(n1,inter1))
   {
    printf("Error: MFDualIntersectIndices at entry, indices not sorted\n");
    fflush(stdout);
   }
  if(!MFTestVertexIndexOrder(n2,inter2))
   {
    printf("Error: MFDualIntersectIndices at entry, indices not sorted\n");
    fflush(stdout);
   }*/

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

/*if(!MFTestVertexIndexOrder(n,inter))
   {
    printf("Error: MFDualIntersectIndices at exit, indices not sorted\n");
    fflush(stdout);
   }*/

  return n;
 }

void MFCreateDualZeroCellV(MFDualCell *cell,MFErrorHandler e, int n,...)
 {
  static char RoutineName[]={"MFCreateDualZeroCellV"};
  va_list indices;
  int i,j;

  *cell=(MFDualCell)malloc(sizeof(struct MFDualCellSt));

#ifndef MFNOSAFETYNET
  if(*cell==NULL)
   {
    sprintf(MFEnumDualPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct MFDualCellSt));
    MFSetError(e,12,RoutineName,MFEnumDualPolytopeErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  (*cell)->nIndices=n;
  (*cell)->indices=(int*)malloc(n*sizeof(int));

#ifndef MFNOSAFETYNET
  if((*cell)->indices==NULL)
   {
    sprintf(MFEnumDualPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(int));
    MFSetError(e,12,RoutineName,MFEnumDualPolytopeErrorMsg,__LINE__,__FILE__);
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

void MFCreateDualZeroCell(MFDualCell *cell,int n,int* indices, MFErrorHandler e)
 {
  static char RoutineName[]={"MFCreateDualZeroCell"};
  int i,j;

  *cell=(MFDualCell)malloc(sizeof(struct MFDualCellSt));

#ifndef MFNOSAFETYNET
  if(*cell==NULL)
   {
    sprintf(MFEnumDualPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct MFDualCellSt));
    MFSetError(e,12,RoutineName,MFEnumDualPolytopeErrorMsg,__LINE__,__FILE__);
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
      sprintf(MFEnumDualPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(int));
      MFSetError(e,12,RoutineName,MFEnumDualPolytopeErrorMsg,__LINE__,__FILE__);
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

void  MFCreateEmptyNComplex(int n, MFEnumDualPolytope *EP,int m, MFErrorHandler e)
 {
  static char RoutineName[]={"MFCreateEmptyNComplex"};
  int i,j;

  *EP=(MFEnumDualPolytope)malloc(sizeof(struct MFEnumDualPolytopeSt));

#ifndef MFNOSAFETYNET
  if(*EP==NULL)
   {
    sprintf(MFEnumDualPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct MFEnumDualPolytopeSt));
    MFSetError(e,12,RoutineName,MFEnumDualPolytopeErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  (*EP)->n=n;
  (*EP)->nCells=(int*)malloc((n+1)*sizeof(int));

#ifndef MFNOSAFETYNET
  if((*EP)->nCells==NULL)
   {
    sprintf(MFEnumDualPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",(n+1)*sizeof(int));
    MFSetError(e,12,RoutineName,MFEnumDualPolytopeErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  (*EP)->mCells=(int*)malloc((n+1)*sizeof(int));

#ifndef MFNOSAFETYNET
  if((*EP)->mCells==NULL)
   {
    sprintf(MFEnumDualPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",(n+1)*sizeof(int));
    MFSetError(e,12,RoutineName,MFEnumDualPolytopeErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  (*EP)->cells=(MFDualCell**)malloc((n+1)*sizeof(MFDualCell*));

#ifndef MFNOSAFETYNET
  if((*EP)->cells==NULL)
   {
    sprintf(MFEnumDualPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",(n+1)*sizeof(MFDualCell*));
    MFSetError(e,12,RoutineName,MFEnumDualPolytopeErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  (*EP)->mCells[0]=m;
  (*EP)->nCells[0]=m;
  if(m>0)
   {
    (*EP)->cells[0]=(MFDualCell*)malloc((*EP)->mCells[0]*sizeof(MFDualCell));

#ifndef MFNOSAFETYNET
    if((*EP)->cells[0]==NULL)
     {
      sprintf(MFEnumDualPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",((*EP)->mCells[0])*sizeof(MFDualCell));
      MFSetError(e,12,RoutineName,MFEnumDualPolytopeErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

    (*EP)->v=(MFNVector*)malloc((*EP)->mCells[0]*sizeof(MFNVector));

#ifndef MFNOSAFETYNET
    if((*EP)->v==NULL)
     {
      sprintf(MFEnumDualPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",((*EP)->mCells[0])*sizeof(MFNVector));
      MFSetError(e,12,RoutineName,MFEnumDualPolytopeErrorMsg,__LINE__,__FILE__);
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
    (*EP)->cells[0]=(MFDualCell*)malloc((*EP)->mCells[0]*sizeof(MFDualCell));

#ifndef MFNOSAFETYNET
    if((*EP)->cells[0]==NULL)
     {
      sprintf(MFEnumDualPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",((*EP)->mCells[0])*sizeof(MFDualCell));
      MFSetError(e,12,RoutineName,MFEnumDualPolytopeErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

    if((*EP)->cells[0]==NULL)fprintf(stderr,"Line %d in file %s - out of space allocating a block of %d bytes\n",__LINE__,__FILE__,((*EP)->mCells[0])*sizeof(MFDualCell));
    (*EP)->v=(MFNVector*)malloc((*EP)->mCells[0]*sizeof(MFNVector));

#ifndef MFNOSAFETYNET
    if((*EP)->v==NULL)
     {
      sprintf(MFEnumDualPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",((*EP)->mCells[0])*sizeof(MFNVector));
      MFSetError(e,12,RoutineName,MFEnumDualPolytopeErrorMsg,__LINE__,__FILE__);
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

  for(i=1;i<=n;i++)
   {
    (*EP)->mCells[i]=EPNumberOfCellsToAllocate;
    (*EP)->nCells[i]=0;
    (*EP)->cells[i]=(MFDualCell*)malloc((*EP)->mCells[i]*sizeof(MFDualCell));

#ifndef MFNOSAFETYNET
    if((*EP)->cells[i]==NULL)
     {
      sprintf(MFEnumDualPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",((*EP)->mCells[i])*sizeof(MFDualCell));
      MFSetError(e,12,RoutineName,MFEnumDualPolytopeErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

    for(j=0;j<(*EP)->mCells[i];j++)((*EP)->cells[i])[j]=NULL;
   }

  return;
 }

int MFDualIndicesEqual(int n1,int *inter1,int n2,int *inter2, MFErrorHandler e)
 {
  static char RoutineName[]={"MFDualIndicesEqual"};
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

void MFFreeEnumDualPolytope(MFEnumDualPolytope EP, MFErrorHandler e)
 {
  static char RoutineName[]={"MFFreeEnumDualPolytope"};
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

    if(EP->v[i]!=NULL)MFFreeNVector(EP->v[i],e);
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
    for(j=0;j<EP->n+1;j++)
     {

#ifdef MFALLOWVERBOSE
      if(verbose){printf("    dimension %d\n",j);fflush(stdout);}
#endif

      for(i=0;i<EP->nCells[j];i++)
       {

#ifdef MFALLOWVERBOSE
        if(verbose){printf("      %d\n",i);fflush(stdout);}
#endif

        if((EP->cells[j])[i]!=NULL)MFFreeDualCell((EP->cells[j])[i],e);
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

void MFFreeDualCell(MFDualCell C, MFErrorHandler e)
 {
  static char RoutineName[]={"MFFreeDualCell"};
  int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("MFFreeDualCell\n");fflush(stdout);}
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
  if(verbose){printf("done MFFreeDualCell\n");fflush(stdout);}
#endif

  return;
 }

long MFEnumDualPolytopeNumberOfVertices(MFEnumDualPolytope EP, MFErrorHandler e)
 {
  static char RoutineName[]={"MFEnumDualPolytopeNumberOfVertices"};

  return EP->nCells[0];
 }

int MFEnumDualPolytopeNumberOfCells(MFEnumDualPolytope EP,int dim, MFErrorHandler e)
 {
  static char RoutineName[]={"MFEnumDualPolytopeNumberOfCells"};

  if(dim>EP->n)return 0;

  return EP->nCells[dim];
 }

int MFEnumDualPolytopeNumberOfCellFaces(MFEnumDualPolytope EP,int dim,long cell, MFErrorHandler e)
 {
  static char RoutineName[]={"MFEnumDualPolytopeNumberOfCellFaces"};

  return ((EP->cells[dim])[cell])->nCellFaces;
 }

long MFEnumDualPolytopeCellFace(MFEnumDualPolytope EP,int dim,long cell,int face, MFErrorHandler e)
 {
  static char RoutineName[]={"MFEnumDualPolytopeCellFace"};

#ifdef MFNOCONFIDENCE
  if(dim<0||dim>EP->n)fprintf(stderr,"Error at line %d in file %s, you asked for a cell face of a %d-cell, but k=%d\n",__LINE__,__FILE__,dim,EP->n);
  if(cell<0||cell>=EP->nCells[dim])fprintf(stderr,"Error at line %d in file %s, you asked for a cell face of %d-cell %d, but there are %d cells of this dimension \n",__LINE__,__FILE__,dim,cell,EP->nCells[dim]);
  if(face<0||face>=((EP->cells[dim])[cell])->nCellFaces){fprintf(stderr,"Error at line %d in file %s, you asked for cell face %d of %d-cell %d, but there are %d cell faces\n",__LINE__,__FILE__,face,dim,cell,((EP->cells[dim])[cell])->nCellFaces);abort();}
#endif

  return ((EP->cells[dim])[cell])->cellFaces[face];
 }

int MFEnumDualPolytopeNumberOfFaceCells(MFEnumDualPolytope EP,int dim,long face, MFErrorHandler e)
 {
  static char RoutineName[]={"MFEnumDualPolytopeNumberOfFaceCells"};

  return ((EP->cells[dim])[face])->nFaceCells;
 }

int MFEnumDualPolytopeNumberOfCellIndices(MFEnumDualPolytope EP,int dim,long cell, MFErrorHandler e)
 {
  static char RoutineName[]={"MFEnumDualPolytopeNumberOfCellIndices"};

  return ((EP->cells[dim])[cell])->nIndices;
 }

int MFEnumDualPolytopeCellIndex(MFEnumDualPolytope EP,int dim,long cell,int i, MFErrorHandler e)
 {
  static char RoutineName[]={"MFEnumDualPolytopeCellIndex"};

  return ((EP->cells[dim])[cell])->indices[i];
 }

long MFEnumDualPolytopeFaceCell(MFEnumDualPolytope EP,int dim,long face,int cell, MFErrorHandler e)
 {
  static char RoutineName[]={"MFEnumDualPolytopeFaceCell"};

  return ((EP->cells[dim])[face])->faceCells[cell];
 }

void EPSortDualPolytopeIndices(int n,int *list, MFErrorHandler e)
 {
  static char RoutineName[]={"EPSortDualPolytopeIndices"};
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

void EPDualCollapseIndices(int *n,int *list, MFErrorHandler e)
 {
  static char RoutineName[]={"EPDualCollapseIndices"};
  int i,j;
  EPSortDualPolytopeIndices(*n,list,e);
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

void MFCreateDualZeroCellWithEdgeList(MFDualCell *cell,int n,int *edges, MFErrorHandler e)
 {
  static char RoutineName[]={"MFCreateDualZeroCellWithEdgeList"};
  int i,j;

  *cell=(MFDualCell)malloc(sizeof(struct MFDualCellSt));

#ifndef MFNOSAFETYNET
  if(*cell==NULL)
   {
    sprintf(MFEnumDualPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct MFDualCellSt));
    MFSetError(e,12,RoutineName,MFEnumDualPolytopeErrorMsg,__LINE__,__FILE__);
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
    sprintf(MFEnumDualPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(int));
    MFSetError(e,12,RoutineName,MFEnumDualPolytopeErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  for(i=0;i<n;i++)(*cell)->faceCells[i]=edges[i];

  return;
 }

int MFEnumDualPolytopeDimension(MFEnumDualPolytope EP, MFErrorHandler e)
 {
  static char RoutineName[]={"MFEnumDualPolytopeDimension"};
  return EP->n;
 }

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

int MFTestEnumDualPolytope(MFEnumDualPolytope EP,int closed, MFErrorHandler e)
 {
  static char RoutineName[]={"MFTestDualEnumPolytope"};
  int k;
  int d;
  int nCells,nFaces,nAdj;
  long cell,face,adj;
  long cellface,adjcell;
  int nFaceAdjs,nAdjFaces;
  long faceadj,adjface;
  int i;
  int present,presentTwice,passed;
  int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose)printf("MFTestEnumDualPolytope:\n");
#endif

  k=MFEnumDualPolytopeDimension(EP,e);

  passed=1;
  for(d=0;d<k+1;d++)
   {
    nCells=MFEnumDualPolytopeNumberOfCells(EP,d,e);
    for(cell=0;cell<nCells;cell++)
     {
      if(d-1>-1)
       {
        nFaces=MFEnumDualPolytopeNumberOfCellFaces(EP,d,cell,e);
        if(closed && nFaces<d+1)
         {
          passed=0;
          printf("  %d-cell %ld has %d faces, fewer than %d\n",d,cell,nFaces,d+1);
         }
        for(face=0;face<nFaces;face++)
         {
          cellface=MFEnumDualPolytopeCellFace(EP,d,cell,face,e);
          presentTwice=0;
          for(i=0;i<nFaces;i++)
           {
            if(i!=face && cellface==MFEnumDualPolytopeCellFace(EP,d,cell,i,e))presentTwice=1;
           }
          if(presentTwice)
           {
            passed=0;
            printf("  A face of %d-cell %ld (%d-cell %ld) appears twice in the face list.\n",d,cell,d-1,cellface);
            printf("The face list of %d-cell %ld is [",d,cell);
            for(i=0;i<nFaces;i++)
             {
              if(i>0)printf(",");
              printf("%d",MFEnumDualPolytopeCellFace(EP,d,cell,i,e));
             }
            printf("]\n");
           }
          nFaceAdjs=MFEnumDualPolytopeNumberOfFaceCells(EP,d-1,cellface,e);
          present=0;
          for(faceadj=0;faceadj<nFaceAdjs;faceadj++)
           {
            if(MFEnumDualPolytopeFaceCell(EP,d-1,cellface,faceadj,e)==cell)present=1;
           }
          if(!present)
           {
            passed=0;
            printf("  Face of %d-cell %ld (%d-cell %ld) does not list %d-cell %ld as a cell it is contained in\n",d,cell,d-1,cellface,d,cell);
            printf("%d-cell %ld is contained in cells (%d-cells) [",d-1,cellface,d);
            for(faceadj=0;faceadj<nFaceAdjs;faceadj++)
             {
              if(adjface>0)printf(",");
              printf("%d",MFEnumDualPolytopeCellFace(EP,d+1,cellface,faceadj,e));
             }
            printf("]\n");
           }
         }
       }
      if(d+1<k)
       {
        nAdj=MFEnumDualPolytopeNumberOfFaceCells(EP,d,cell,e);
        for(adj=0;adj<nAdj;adj++)
         {
          adjcell=MFEnumDualPolytopeFaceCell(EP,d,cell,adj,e);
          presentTwice=0;
          for(i=0;i<nAdj;i++)
           {
            if(i!=adj && adjcell==MFEnumDualPolytopeFaceCell(EP,d,cell,i,e))presentTwice=1;
           }
          if(presentTwice)
           {
            passed=0;
            printf("  A containing cell of %d-cell %ld (%d-cell %ld) appears twice in the list of containing cells.\n",d,cell,d+1,adjcell);
            printf("The containing cell list of %d-cell %ld is [",d,cell);
            for(i=0;i<nAdj;i++)
             {
              if(i>0)printf(",");
              printf("%d",MFEnumDualPolytopeFaceCell(EP,d,cell,i,e));
             }
            printf("]\n");
           }
          nAdjFaces=MFEnumDualPolytopeNumberOfCellFaces(EP,d+1,adjcell,e);
          present=0;
          for(adjface=0;adjface<nAdjFaces;adjface++)
           {
            if(MFEnumDualPolytopeCellFace(EP,d+1,adjcell,adjface,e)==cell)present=1;
           }
          if(!present)
           {
            passed=0;
            printf("  Containing cell of %d-cell %ld (%d-cell %ld) does not list %d-cell %ld as a face\n",d,cell,d+1,adjcell,d,cell);
            printf("%d-cell %ld has faces (%d-cells) [",d+1,adjcell,d);
            for(adjface=0;adjface<nAdjFaces;adjface++)
             {
              if(adjface>0)printf(",");
              printf("%d",MFEnumDualPolytopeCellFace(EP,d+1,adjcell,adjface,e));
             }
            printf("]\n");
           }
         }
       }
     }
   }

  if(passed)printf("  passed without problems\n");
  return passed;
 }

int MFEnumDualPolytopeVerticesOnDualCell(MFEnumDualPolytope P,int dim,int cell,int nV,int **indices, MFErrorHandler e)
 {
  static char RoutineName[]={"MFEnumDualPolytopeVerticesOnDualCell"};
  int nF,f;
  int n;
  MFDualCell c;
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
        sprintf(MFEnumDualPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",10*sizeof(int));
        MFSetError(e,12,RoutineName,MFEnumDualPolytopeErrorMsg,__LINE__,__FILE__);
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
        sprintf(MFEnumDualPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",(nV+10)*sizeof(int));
        MFSetError(e,12,RoutineName,MFEnumDualPolytopeErrorMsg,__LINE__,__FILE__);
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
   }
  if(verbose){printf("     There are %d faces\n",nF);fflush(stdout);}
#endif

  n=nV;
  for(f=0;f<nF;f++)
   {

#ifdef MFALLOWVERBOSE
    if(verbose){printf("     Check face %d, %d-cell %d\n",f,dim-1,c->faceCells[f]);fflush(stdout);}
#endif

    n=MFEnumDualPolytopeVerticesOnDualCell(P,dim-1,c->cellFaces[f],n,indices,e);
   }

  return n;
 }

MFEnumDualPolytope MFEnumDualOfAtlas(MFAtlas A, MFErrorHandler e)
 {
  static char RoutineName[]={"MFEnumDualOfAtlas"};
  MFEnumPolytope Complex;
  MFEnumDualPolytope DualComplex;
  int nCharts;
  int chart;
  int failed;
  int verboseResult=0;
  int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("In %s\n",RoutineName);fflush(stdout);}
  if(verbose){printf("  Call MFEnumerateAtlas\n");fflush(stdout);}
#endif

  Complex=MFEnumerateAtlas(A,e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("  Call MFCreateDualOfEnumPolytope\n");fflush(stdout);}
#endif

  DualComplex=MFCreateDualOfEnumPolytope(Complex,e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("  Call MFCreateDualOfEnumPolytope\n");fflush(stdout);}
#endif

  nCharts=MFAtlasNumberOfCharts(A,e);
  for(chart=0;chart<nCharts;chart++)
   {
    DualComplex->v[chart]=MFChartCenter(MFAtlasChart(A,chart,e),e);
    MFRefNVector(DualComplex->v[chart],e);
   }

#ifdef MFALLOWVERBOSE
  if(verbose){printf("  Test\n");fflush(stdout);}
  if(0 && (failed=!MFTestEnumDualPolytope(DualComplex,1,e)) || verboseResult )
   {
    if(!verboseResult)
     {
      printf("----------------------------------------------------\n");
      printf("----------------------------------------------------\n");
     }
    printf("\nDual Complex failed tests:\n\n");
    MFPrintEnumDualPolytope(stdout,DualComplex,e);
    printf("\n");
    printf("----------------------------------------------------\n");
    printf("----------------------------------------------------\n");
    fflush(stdout);
   }
#endif

#ifdef MFALLOWVERBOSE
  if(verbose){printf("  Call MFFreeEnumPolytope\n");fflush(stdout);}
#endif

  MFFreeEnumPolytope(Complex,e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("Done %s\n",RoutineName);fflush(stdout);}
#endif

  return DualComplex;
 }

MFNVector MFEnumDualPolytopeVertex(MFEnumDualPolytope P,long i, MFErrorHandler e)
 {
  static char RoutineName[]={"MFEnumDualPolytopeVertex"};

  return P->v[i];
 }

#ifdef HELLIFIKNOWWHYTHEREARETWOOFTHESE
int MFTestEnumDualPolytope(MFEnumDualPolytope EP,int closed, MFErrorHandler e)
 {
  static char RoutineName[]={"MFTestEnumDualPolytope"};
  int k;
  int d;
  int nCells,nFaces,nAdj;
  long cell,face,adj;
  long cellface,adjcell;
  int nFaceAdjs,nAdjFaces;
  long faceadj,adjface;
  int i;
  int present,presentTwice,passed;
  int verbose=0;

  k=MFEnumDualPolytopeDimension(EP,e);

#ifdef MFALLOWVERBOSE
  if(verbose)printf("MFTestEnumDualPolytope:\n");
#endif

  passed=1;
  for(d=0;d<k+1;d++)
   {
    nCells=MFEnumDualPolytopeNumberOfCells(EP,d,e);
    for(cell=0;cell<nCells;cell++)
     {
      if(d-1>-1)
       {
        nFaces=MFEnumDualPolytopeNumberOfCellFaces(EP,d,cell,e);
        if(closed && nFaces<d+1)
         {
          passed=0;
          printf("  %d-cell %ld has %d faces, fewer than %d\n",d,cell,nFaces,d+1);
          for(face=0;face<nFaces;face++)
           {
            cellface=MFEnumDualPolytopeCellFace(EP,d,cell,face);
            printf("     face %d is %d-cell %d\n",face,d-1,cellface);
           }
         }
        for(face=0;face<nFaces;face++)
         {
          cellface=MFEnumDualPolytopeCellFace(EP,d,cell,face,e);
          presentTwice=0;
          for(i=0;i<nFaces;i++)
           {
            if(i!=face && cellface==MFEnumDualPolytopeCellFace(EP,d,cell,i,e))presentTwice=1;
           }
          if(presentTwice)
           {
            passed=0;
            printf("  A face of %d-cell %ld (%d-cell %ld) appears twice in the face list.\n",d,cell,d-1,cellface);
            printf("The face list of %d-cell %ld is [",d,cell);
            for(i=0;i<nFaces;i++)
             {
              if(i>0)printf(",");
              printf("%d",MFEnumDualPolytopeCellFace(EP,d,cell,i,e));
             }
            printf("]\n");
           }
          nFaceAdjs=MFEnumDualPolytopeNumberOfFaceCells(EP,d-1,cellface,e);
          present=0;
          for(faceadj=0;faceadj<nFaceAdjs;faceadj++)
           {
            if(MFEnumDualPolytopeFaceCell(EP,d-1,cellface,faceadj,e)==cell)present=1;
           }
          if(!present)
           {
            passed=0;
            printf("  Face of %d-cell %ld (%d-cell %ld) does not list %d-cell %ld as a cell it is contained in\n",d,cell,d-1,cellface,d,cell);
            printf("%d-cell %ld is contained in cells (%d-cells) [",d-1,cellface,d);
            for(faceadj=0;faceadj<nFaceAdjs;faceadj++)
             {
              if(faceadj>0)printf(",");
              printf("%d",MFEnumDualPolytopeFaceCell(EP,d-1,cellface,faceadj,e));
             }
            printf("]\n");
           }
         }
       }
      if(d+1<k)
       {
        nAdj=MFEnumDualPolytopeNumberOfFaceCells(EP,d,cell,e);
        for(adj=0;adj<nAdj;adj++)
         {
          adjcell=MFEnumDualPolytopeFaceCell(EP,d,cell,adj,e);
          presentTwice=0;
          for(i=0;i<nAdj;i++)
           {
            if(i!=adj && adjcell==MFEnumDualPolytopeFaceCell(EP,d,cell,i,e))presentTwice=1;
           }
          if(presentTwice)
           {
            passed=0;
            printf("  A containing cell of %d-cell %ld (%d-cell %ld) appears twice in the list of containing cells.\n",d,cell,d+1,adjcell);
            printf("The containing cell list of %d-cell %ld is [",d,cell);
            for(i=0;i<nAdj;i++)
             {
              if(i>0)printf(",");
              printf("%d",MFEnumDualPolytopeFaceCell(EP,d,cell,i,e));
             }
            printf("]\n");
           }
          nAdjFaces=MFEnumDualPolytopeNumberOfCellFaces(EP,d+1,adjcell,e);
          present=0;
          for(adjface=0;adjface<nAdjFaces;adjface++)
           {
            if(MFEnumDualPolytopeCellFace(EP,d+1,adjcell,adjface,e)==cell)present=1;
           }
          if(!present)
           {
            passed=0;
            printf("  Containing cell of %d-cell %ld (%d-cell %ld) does not list %d-cell %ld as a face\n",d,cell,d+1,adjcell,d,cell);
            printf("%d-cell %ld has faces (%d-cells) [",d+1,adjcell,d);
            for(adjface=0;adjface<nAdjFaces;adjface++)
             {
              if(adjface>0)printf(",");
              printf("%d",MFEnumDualPolytopeCellFace(EP,d+1,adjcell,adjface,e));
             }
            printf("]\n");
           }
         }
       }
     }
   }
  return passed;
 }
#endif

MFEnumDualPolytope MFCreateDualOfEnumPolytope(MFEnumPolytope EP, MFErrorHandler e)
 {
  static char RoutineName[]={"MFCreateDualOfEnumPolytope"};
  MFEnumDualPolytope dual=NULL;
  int dim;
  int n;
  int i,j;
  int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("MFCreateDualOfEnumPolytope - initial enumerated complex:\n");fflush(stdout);}
  if(verbose)MFPrintEnumPolytope(stdout,EP,e);
#endif
 
  n=MFEnumPolytopeDimension(EP,e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("   The complex has %d %d-cells\n",MFEnumPolytopeNumberOfCells(EP,n,e),n);fflush(stdout);}
#endif

  MFCreateEmptyNComplex(n,&dual,MFEnumPolytopeNumberOfCells(EP,n,e),e);

  for(dim=0;dim<n+1;dim++)
   {

#ifdef MFALLOWVERBOSE
    if(verbose){printf("   dimension %d in complex is dimension %d in dual\n",dim,n-dim);fflush(stdout);}
#endif

    dual->nCells[n-dim]=MFEnumPolytopeNumberOfCells(EP,dim,e);
    dual->mCells[n-dim]=dual->nCells[n-dim];

#ifdef MFALLOWVERBOSE
    if(verbose){printf("      There will be %d %d-cells in the dual\n",dual->nCells[n-dim],n-dim);fflush(stdout);}
#endif

    if(dual->nCells[n-dim]>0)
     {
      dual->cells[n-dim]=(MFDualCell*)realloc((void*)(dual->cells[n-dim]),dual->mCells[n-dim]*sizeof(MFDualCell));

#ifndef MFNOSAFETYNET
      if(dual->cells[n-dim]==NULL)
       {
        sprintf(MFEnumDualPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",dual->mCells[n-dim]*sizeof(MFDualCell));
        MFSetError(e,12,RoutineName,MFEnumDualPolytopeErrorMsg,__LINE__,__FILE__);
        MFErrorHandlerOutOfMemory(e);
        return NULL;
       }
#endif

      for(i=0;i<dual->nCells[n-dim];i++)
       {
        (dual->cells[n-dim])[i]=(MFDualCell)malloc(sizeof(struct MFDualCellSt));

#ifndef MFNOSAFETYNET
        if((dual->cells[n-dim])[i]==NULL)
         {
          sprintf(MFEnumDualPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct MFDualCellSt));
          MFSetError(e,12,RoutineName,MFEnumDualPolytopeErrorMsg,__LINE__,__FILE__);
          MFErrorHandlerOutOfMemory(e);
          return NULL;
         }
#endif
        
        ((dual->cells[n-dim])[i])->nIndices=MFEnumPolytopeNumberOfCellIndices(EP,dim,i,e);

#ifdef MFALLOWVERBOSE
        if(verbose){printf("   Setting the indices of dual %d-cell %d (there are %d of them)\n",n-dim,i,((dual->cells[n-dim])[i])->nIndices);fflush(stdout);}
#endif

        if(((dual->cells[n-dim])[i])->nIndices>0)
         {
          ((dual->cells[n-dim])[i])->indices=(int*)malloc((((dual->cells[n-dim])[i])->nIndices)*sizeof(int));

#ifndef MFNOSAFETYNET
          if(((dual->cells[n-dim])[i])->indices==NULL)
           {
            sprintf(MFEnumDualPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",(((dual->cells[n-dim])[i])->nIndices)*sizeof(int));
            MFSetError(e,12,RoutineName,MFEnumDualPolytopeErrorMsg,__LINE__,__FILE__);
            MFErrorHandlerOutOfMemory(e);
            return NULL;
           }
#endif

          for(j=0;j<((dual->cells[n-dim])[i])->nIndices;j++)
           ((dual->cells[n-dim])[i])->indices[j]=MFEnumPolytopeCellIndex(EP,dim,i,j,e);
         }else{
          ((dual->cells[n-dim])[i])->indices=NULL;
         }

        ((dual->cells[n-dim])[i])->nFaceCells=MFEnumPolytopeNumberOfCellFaces(EP,dim,i,e);
        ((dual->cells[n-dim])[i])->mFaceCells=((dual->cells[n-dim])[i])->nFaceCells;

#ifdef MFALLOWVERBOSE
        if(verbose){printf("   Setting the facecells of dual %d-cell %d (there are %d of them)\n",n-dim,i,((dual->cells[n-dim])[i])->nFaceCells);fflush(stdout);}
#endif

        if(((dual->cells[n-dim])[i])->nFaceCells>0)
         {
          ((dual->cells[n-dim])[i])->faceCells=(int*)malloc((((dual->cells[n-dim])[i])->nFaceCells)*sizeof(int));

#ifndef MFNOSAFETYNET
          if(((dual->cells[n-dim])[i])->faceCells==NULL)
           {
            sprintf(MFEnumDualPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",(((dual->cells[n-dim])[i])->nFaceCells)*sizeof(int));
            MFSetError(e,12,RoutineName,MFEnumDualPolytopeErrorMsg,__LINE__,__FILE__);
            MFErrorHandlerOutOfMemory(e);
            return NULL;
           }
#endif

          for(j=0;j<((dual->cells[n-dim])[i])->nFaceCells;j++)
           ((dual->cells[n-dim])[i])->faceCells[j]=MFEnumPolytopeCellFace(EP,dim,i,j,e);
         }else{
          ((dual->cells[n-dim])[i])->faceCells=NULL;
         }

        ((dual->cells[n-dim])[i])->nCellFaces=MFEnumPolytopeNumberOfFaceCells(EP,dim,i,e);
        ((dual->cells[n-dim])[i])->mCellFaces=((dual->cells[n-dim])[i])->nCellFaces;

#ifdef MFALLOWVERBOSE
        if(verbose){printf("   Setting the cellFaces of dual %d-cell %d (there are %d of them)\n",n-dim,i,((dual->cells[n-dim])[i])->nCellFaces);fflush(stdout);}
#endif

        if(((dual->cells[n-dim])[i])->nCellFaces>0)
         {
          ((dual->cells[n-dim])[i])->cellFaces=(int*)malloc((((dual->cells[n-dim])[i])->nCellFaces)*sizeof(int));

#ifndef MFNOSAFETYNET
          if(((dual->cells[n-dim])[i])->cellFaces==NULL)
           {
            sprintf(MFEnumDualPolytopeErrorMsg,"Out of memory, trying to allocate %d bytes",(((dual->cells[n-dim])[i])->nCellFaces)*sizeof(int));
            MFSetError(e,12,RoutineName,MFEnumDualPolytopeErrorMsg,__LINE__,__FILE__);
            MFErrorHandlerOutOfMemory(e);
            return NULL;
           }
#endif

          for(j=0;j<((dual->cells[n-dim])[i])->nCellFaces;j++)
           ((dual->cells[n-dim])[i])->cellFaces[j]=MFEnumPolytopeFaceCell(EP,dim,i,j,e);
         }else{
          ((dual->cells[n-dim])[i])->cellFaces=NULL;
         }
       }
     }
   }

  return dual;
 }

#ifdef __cplusplus
}
#endif
