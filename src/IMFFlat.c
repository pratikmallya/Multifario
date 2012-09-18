/*
 *  PROGRAM NAME:  multifario
 *
 *  (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
 *  CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
 *
 *  Please refer to the LICENSE file in the top directory
 *
 *      author: Mike Henderson mhender@watson.ibm.com
 *      date:   September 15, 2003
 */

static char *id="@(#) $Id: IMFFlat.c,v 1.6 2011/07/21 17:42:46 mhender Exp $";

static char MFFlatMFErrorHandlerMsg[256]="";

#include <MFImplicitMF.h>
#include <IMFExpansionSpace.h>
#include <IMFExpansion.h>
#include <IMFExpansionPt.h>
#include <MFErrorHandler.h>
#include <stdlib.h>
#include <string.h>

#ifdef __cplusplus
 extern "C" {
#endif

static void MFFreeNSpaceData(void*,MFErrorHandler);
static int MFProjectNSpace(int,int,MFNVector,MFNKMatrix,MFNVector,void*,int*,MFErrorHandler);
static int MFTangentNSpace(int,int,MFNVector,MFNKMatrix,void*,MFErrorHandler);
static double MFScaleNSpace(int,int,MFNVector,MFNKMatrix,void*,MFErrorHandler);
static void MFWriteNSpaceData(FILE*,void*,MFErrorHandler);
static MFImplicitMF MFReadNSpaceMF(FILE*,MFErrorHandler);

static int MFFlatMFProjectToSave(MFNVector,double*,void*,MFErrorHandler);
static int MFFlatMFProjectToDraw(MFNVector,double*,void*,MFErrorHandler);
static int MFFlatMFProjectForBB(MFNVector,double*,void*,MFErrorHandler);

MFNVector MFNVectorFactory(MFImplicitMF,MFErrorHandler);
MFNKMatrix MFNKMatrixFactory(MFImplicitMF,MFErrorHandler);

struct IMFFlatData
 {
  int n;
  int k;
 };

/*! \fn MFImplicitMF IMFCreateFlat(int n, int k, MFErrorHandler e)
 *  \brief A flat k-manifold embedded in Euclidean n-space
 *  
 *  \param n The dimension of the embedding space.
 *  \param k The dimension of the manifold.
 *  \param e An error handler.
 *  \returns A new MFImplicitMF.
 */
MFImplicitMF IMFCreateFlat(int n, int k, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFCreateFlat"};
  MFImplicitMF nSpace;
  struct IMFFlatData *data;
  MFNSpace space;
  int verbose=0;

  if(verbose)printf("%s, n=%d, k=%d\n",RoutineName,n,k);
  nSpace=MFIMFCreateBaseClass(n,k,"IMFFlat",e);

  space=IMFCreateExpansionSpace(n,e);
  MFIMFSetSpace(nSpace,space,e);
  MFFreeNSpace(space,e);

  MFIMFSetFreeData(nSpace,MFFreeNSpaceData,e);
  MFIMFSetProjectForSave(nSpace,MFFlatMFProjectToSave,e);
  MFIMFSetProjectForDraw(nSpace,MFFlatMFProjectToDraw,e);
  MFIMFSetProjectForBB(nSpace,MFFlatMFProjectForBB,e);

  MFIMFSetVectorFactory(nSpace,MFNVectorFactory,e);
  MFIMFSetMatrixFactory(nSpace,MFNKMatrixFactory,e);

  data=(struct IMFFlatData*)malloc(sizeof(struct IMFFlatData));

#ifndef MFNOSAFETYNET
  if(data==NULL)
   {
    sprintf(MFFlatMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct IMFFlatData));
    MFSetError(e,12,RoutineName,MFFlatMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  data->n=n;
  data->k=k;
  MFIMFSetData(nSpace,(void*)data,e);
  if(verbose)printf("done %s\n",RoutineName);

  return nSpace;
 }

void MFFreeNSpaceData(void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFFreeNSpaceData"};
  struct IMFFlatData *data;

  data=(struct IMFFlatData*)d;
  if(data!=NULL)free(data);
  return;
 }

int MFFlatMFProjectToSave(MFNVector u, double *x, void *d, MFErrorHandler e)
 {
  int i;
  struct IMFFlatData *data;

  data=(struct IMFFlatData*)d;

  if(x==NULL)return data->n;

  for(i=0;i<data->n;i++)x[i]=(IMFExpansionU((IMFExpansionNVGetE(u,e)),e))[i];

  return 0;
 }

int MFFlatMFProjectForBB(MFNVector u, double *x, void *d, MFErrorHandler e)
 {
  int i;
  struct IMFFlatData *data;

  data=(struct IMFFlatData*)d;

  if(x==NULL)return data->n;
/*if(x==NULL)return data->n+1;*/

  for(i=0;i<data->n;i++)x[i]=(IMFExpansionU(IMFExpansionNVGetE(u,e),e))[i];
/*x[data->n]=.005*IMFExpansionNVGetT(u,e);*/

  return 0;
 }

int MFFlatMFProjectToDraw(MFNVector u, double *x, void *d, MFErrorHandler e)
 {
  int i,j;
  struct IMFFlatData *data;

  data=(struct IMFFlatData*)d;

  if(x==NULL)return data->n;

  for(i=0;i<data->n;i++)
   {
    x[i]=MFNV_C(u,i,e);
    if(x[i]!=x[i]){printf("NaN in FlatProjectToDraw\n");fflush(stdout);}
   }

  return 0;
 }

#ifdef __cplusplus
}
#endif
