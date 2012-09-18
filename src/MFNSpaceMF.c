/*
 *  PROGRAM NAME:  multifario
 *
 *  (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
 *  CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
 *
 *  Please refer to the LICENSE file in the top directory
 *
 *      author: Mike Henderson mhender@watson.ibm.com
 *      date:   March 1, 1999
 */

static char *id="@(#) $Id: MFNSpaceMF.c,v 1.3 2007/02/13 01:22:34 mhender Exp $";

static char MFNSpaceMFErrorHandlerMsg[256]="";

#include <MFImplicitMF.h>
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

static int MFNSpaceMFProjectToSave(MFNVector,double*,void*,MFErrorHandler);
static int MFNSpaceMFProjectToDraw(MFNVector,double*,void*,MFErrorHandler);
static int MFNSpaceMFProjectForBB(MFNVector,double*,void*,MFErrorHandler);

MFNVector MFNVectorFactory(MFImplicitMF,MFErrorHandler);
MFNKMatrix MFNKMatrixFactory(MFImplicitMF,MFErrorHandler);

struct MFIMFNSpaceData
 {
  int n;
  double R;
 };

/*! \fn MFImplicitMF MFIMFCreateNSpaceWithRadius(int k,double R);
 *  \brief Creates Euclidean k-space embedded in k-space. This is the same as MFIMFCreateNSpace(k) and MFSetR(.,R)
 *
 *  \param   k The dimension of the manifold and the embedding space.
 *  \param   R The radius.
 *  \returns An implicitly defined manifold.
 */
MFImplicitMF MFIMFCreateNSpaceWithRadius(int n, double R, MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFCreateNSpaceWithRadius"};
  MFImplicitMF nSpace;
  struct MFIMFNSpaceData *data;
  MFNSpace space;

  nSpace=MFIMFCreateBaseClass(n,n,"NSpaceMF",e);

  space=MFCreateNSpace(n,e);
  MFIMFSetSpace(nSpace,space,e);
  MFFreeNSpace(space,e);

  MFIMFSetFreeData(nSpace,MFFreeNSpaceData,e);
  MFIMFSetProject(nSpace,MFProjectNSpace,e);
  MFIMFSetTangent(nSpace,MFTangentNSpace,e);
  MFIMFSetScale(nSpace,MFScaleNSpace,e);
  MFIMFSetR(nSpace,R,e);
  MFIMFSetWriteData(nSpace,MFWriteNSpaceData,e);
  MFIMFSetProjectForSave(nSpace,MFNSpaceMFProjectToSave,e);
  MFIMFSetProjectForDraw(nSpace,MFNSpaceMFProjectToDraw,e);
  MFIMFSetProjectForBB(nSpace,MFNSpaceMFProjectForBB,e);

  MFIMFSetVectorFactory(nSpace,MFNVectorFactory,e);
  MFIMFSetMatrixFactory(nSpace,MFNKMatrixFactory,e);

  data=(struct MFIMFNSpaceData*)malloc(sizeof(struct MFIMFNSpaceData));

#ifndef MFNOSAFETYNET
  if(data==NULL)
   {
    sprintf(MFNSpaceMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct MFIMFNSpaceData));
    MFSetError(e,12,RoutineName,MFNSpaceMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  data->n=n;
  data->R=R;
  MFIMFSetData(nSpace,(void*)data,e);

  return nSpace;
 }

/*! \fn MFImplicitMF MFIMFCreateNSpace(int k);
 *  \brief Creates Euclidean k-space embedded in k-space.
 *
 *  \returns An implicitly defined manifold.
 */
MFImplicitMF MFIMFCreateNSpace(int n, MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFCreateNSpace"};
  MFImplicitMF nSpace;

  nSpace=MFIMFCreateNSpaceWithRadius(n,-1.,e);

  return nSpace;
 }

void MFFreeNSpaceData(void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFFreeNSpaceData"};
  struct MFIMFNSpaceData *data;

  data=(struct MFIMFNSpaceData*)d;
  if(data!=NULL)free(data);
  return;
 }

int MFProjectNSpace(int n,int k,MFNVector vu0,MFNKMatrix mPhi,MFNVector vu,void *d,int *index, MFErrorHandler e)
 {
  static char RoutineName[]={"MFProjectNSpace"};
  static int i;
  double *u,*u0,*Phi;
  struct MFIMFNSpaceData *data;

  data=(struct MFIMFNSpaceData*)d;

  u0=MFNV_CStar(vu0,e);
  Phi=MFNKM_CStar(mPhi,e);
  u=MFNV_CStar(vu,e);

  for(i=0;i<n;i++)u[i]=u0[i];
  *index=0;
  return 1;
 }

int MFTangentNSpace(int n,int k,MFNVector vu,MFNKMatrix mPhi,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFTangentNSpace"};
  static int i;
  double *u,*Phi;
  struct MFIMFNSpaceData *data;

  data=(struct MFIMFNSpaceData*)d;

  Phi=MFNKM_CStar(mPhi,e);
  u=MFNV_CStar(vu,e);

  for(i=0;i<n*n;i++)Phi[i]=0.;
  for(i=0;i<n;i++)Phi[i*(n+1)]=1.;

  return 1;
 }

double MFScaleNSpace(int n,int k,MFNVector u,MFNKMatrix Phi,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFScaleNSpace"};
  struct MFIMFNSpaceData *data;

  data=(struct MFIMFNSpaceData*)d;
  return .1;
 }

void MFWriteNSpaceData(FILE *fid,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFWriteNSpaceData"};
  struct MFIMFNSpaceData *data;

  data=(struct MFIMFNSpaceData*)d;

  fprintf(fid,"%d %lf\n",data->n,data->R);
  return;
 }

MFImplicitMF MFReadNSpaceMF(FILE *fid, MFErrorHandler e)
 {
  static char RoutineName[]={"MFReadNSpaceMF"};
  int n=0;
  double R=0.;
  MFImplicitMF nSpace;

  fscanf(fid,"%d %lf\n",&n,&R);
  nSpace=MFIMFCreateNSpaceWithRadius(n,R,e);

  return nSpace;
 }

int MFNSpaceMFProjectToSave(MFNVector u, double *x, void *d, MFErrorHandler e)
 {
  int i;
  struct MFIMFNSpaceData *data;

  data=(struct MFIMFNSpaceData*)d;

  if(x==NULL)return data->n;

  for(i=0;i<data->n;i++)x[i]=MFNV_C(u,i,e);

  return 0;
 }

int MFNSpaceMFProjectToDraw(MFNVector u, double *x, void *d, MFErrorHandler e)
 {
  int i;
  struct MFIMFNSpaceData *data;

  data=(struct MFIMFNSpaceData*)d;

  if(x==NULL)return data->n;

  for(i=0;i<data->n;i++)x[i]=MFNV_C(u,i,e);

  return 0;
 }

int MFNSpaceMFProjectForBB(MFNVector u, double *x, void *d, MFErrorHandler e)
 {
  int i;
  struct MFIMFNSpaceData *data;

  data=(struct MFIMFNSpaceData*)d;

  if(x==NULL)return data->n;

  for(i=0;i<data->n;i++)x[i]=MFNV_C(u,i,e);

  return 0;
 }

/*! @} */

/*! @} */

#ifdef __cplusplus
}
#endif
