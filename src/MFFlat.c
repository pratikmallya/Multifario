/*
 *  PROGRAM NAME:  multifario
 *
 *  (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
 *  CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
 *
 *  Please refer to the LICENSE file in the top directory
 *
 *      author: Mike Henderson mhender@watson.ibm.com
 *      date:   April 27, 1999
 */

static char *id="@(#) $Id: MFFlat.c,v 1.4 2011/07/21 17:42:46 mhender Exp $";

static char MFFlatMFErrorHandlerMsg[256]="";

#include <MFImplicitMF.h>
#include <MFErrorHandler.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef __cplusplus
 extern "C" {
#endif

static void MFFreeFlatData(void*,MFErrorHandler);
static int MFProjectFlat(int,int,MFNVector,MFNKMatrix,MFNVector,void*,int*,MFErrorHandler);
static int MFTangentFlat(int,int,MFNVector,MFNKMatrix,void*,MFErrorHandler);
static double MFScaleFlat(int,int,MFNVector,MFNKMatrix,void*,MFErrorHandler);
static void MFWriteFlatData(FILE*,void*,MFErrorHandler);
static MFImplicitMF MFReadFlat(FILE*,MFErrorHandler);

static int MFFlatProjectToSave(MFNVector,double*,void *d,MFErrorHandler);
static int MFFlatProjectToDraw(MFNVector,double*,void *d,MFErrorHandler);
static int MFFlatProjectForBB(MFNVector,double*,void *d,MFErrorHandler);

MFNVector MFNVectorFactory(MFImplicitMF,MFErrorHandler);
MFNKMatrix MFNKMatrixFactory(MFImplicitMF,MFErrorHandler);

struct MFFlatData
 {
  int n;
  int k;
  double *o;
  double *v;
 };

/*! \fn MFImplicitMF MFIMFCreateFlat(int n, int k, double *o,double *v, MFErrorHandler e);
 *  \brief Creates a flat k-space embedded in n-space. o is an array of length n, and v is an array of length n*k giving 
 *         the tangent plane. The jth basis vector is v[i+n*j]. The manifold is o[i]+sum_j v[i+n*j].s[j].
 *
 *  \param n The dimension of the embedding space.
 *  \param k The dimension of the manifold.
 *  \param o A point on the manifold.
 *  \param v The (constant) tangent space.  The jth basis vector is v[i+n*j].
 *  \param e An error handler.
 *  \returns An implicitly defined manifold.
 */
MFImplicitMF MFIMFCreateFlat(int n, int k, double *o,double *v, MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFCreateFlat"};
  MFImplicitMF flat;
  MFNSpace space;
  struct MFFlatData *data;
  int i;

  flat=MFIMFCreateBaseClass(n,k,"Flat",e);

  space=MFCreateNSpace(n,e);
  MFIMFSetSpace(flat,space,e);
  MFFreeNSpace(space,e);

  data=(struct MFFlatData*)malloc(sizeof(struct MFFlatData));

#ifndef MFNOSAFETYNET
  if(data==NULL)
   {
    sprintf(MFFlatMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct MFFlatData));
    MFSetError(e,12,RoutineName,MFFlatMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(flat);
    return NULL;
   }
#endif

  data->n=n;
  data->k=k;
  data->o=(double*)malloc(n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(data->o==NULL)
   {
    sprintf(MFFlatMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(double));
    MFSetError(e,12,RoutineName,MFFlatMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(flat);
    free(data);
    return NULL;
   }
#endif

  for(i=0;i<n;i++)data->o[i]=o[i];
  data->v=(double*)malloc(k*n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(data->v==NULL)
   {
    sprintf(MFFlatMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*k*sizeof(double));
    MFSetError(e,12,RoutineName,MFFlatMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(flat);
    free(data->o);
    free(data);
    return NULL;
   }
#endif

  for(i=0;i<k*n;i++)data->v[i]=v[i];

  MFIMFSetData(flat,(void*)data,e);
  MFIMFSetFreeData(flat,MFFreeFlatData,e);
  MFIMFSetProject(flat,MFProjectFlat,e);
  MFIMFSetTangent(flat,MFTangentFlat,e);
  MFIMFSetScale(flat,MFScaleFlat,e);
  MFIMFSetWriteData(flat,MFWriteFlatData,e);
  MFIMFSetProjectForSave(flat,MFFlatProjectToSave,e);
  MFIMFSetProjectForDraw(flat,MFFlatProjectToDraw,e);
  MFIMFSetProjectForBB(flat,MFFlatProjectForBB,e);

  MFIMFSetVectorFactory(flat,MFNVectorFactory,e);
  MFIMFSetMatrixFactory(flat,MFNKMatrixFactory,e);

  return flat;
 }

void MFFreeFlatData(void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFFreeFlatData"};
  struct MFFlatData *data;

  if(d==NULL)return;
  data=(struct MFFlatData*)d;
  if(data->o!=NULL)free(data->o);
  if(data->v!=NULL)free(data->v);
  free(data);
  return;
 }

int MFProjectFlat(int n,int k,MFNVector vu0,MFNKMatrix mPhi,MFNVector vu,void *d,int *index, MFErrorHandler e)
 {
  static char RoutineName[]={"MFProjectFlat"};
  int i,j,l;
  struct MFFlatData *data;
  double *u0, *u;
  double s;

  u0=MFNV_CStar(vu0,e);
  u=MFNV_CStar(vu,e);

  data=(struct MFFlatData*)d;

  for(i=0;i<data->n;i++)
   {
    u[i]=data->o[i];
    for(j=0;j<data->k;j++)
     {
      s=0.;for(l=0;l<data->n;l++)s+=(u0[l]-data->o[l])*(data->v)[l+(data->n)*j];
      u[i]+=s*(data->v)[i+(data->n)*j];
     }
   }
  *index=0;
  return 1;
 }

int MFTangentFlat(int n,int k,MFNVector vu,MFNKMatrix mPhi,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFTangentFlat"};
  int i;
  struct MFFlatData *data;
  double *Phi;

  data=(struct MFFlatData*)d;

  Phi=MFNKM_CStar(mPhi,e);

  for(i=0;i<n*k;i++)Phi[i]=data->v[i];

  return 1;
 }

double MFScaleFlat(int n,int k,MFNVector vu,MFNKMatrix mPhi,void *d, MFErrorHandler e)
 {
  return .3;
 }

void MFWriteFlatData(FILE *fid,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFWriteFlatData"};
  struct MFFlatData *data;
  int i;

  data=(struct MFFlatData*)d;

  fprintf(fid,"%d %d\n",data->n,data->k);
  for(i=0;i<data->n;i++)fprintf(fid,"%lf",data->o[i]);
  for(i=0;i<data->k*data->n;i++)fprintf(fid,"%lf",data->v[i]);

  fflush(fid);
  return;
 }

MFImplicitMF MFReadFlat(FILE *fid, MFErrorHandler e)
 {
  static char RoutineName[]={"MFReadFlat"};
  int n=0;
  int k=0;
  double *o;
  double *v;
  MFImplicitMF flat;
  struct MFFlatData *data;
  int i;

  fscanf(fid,"%d %d\n",&n,&k);

  o=(double*)malloc(n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(o==NULL)
   {
    sprintf(MFFlatMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(double));
    MFSetError(e,12,RoutineName,MFFlatMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(flat);
    free(data);
    return NULL;
   }
#endif

  for(i=0;i<n;i++)fscanf(fid,"%lf",o+i);
  v=(double*)malloc(k*n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(v==NULL)
   {
    sprintf(MFFlatMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*k*sizeof(double));
    MFSetError(e,12,RoutineName,MFFlatMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(o);
    return NULL;
   }
#endif

  for(i=0;i<k*n;i++)fscanf(fid,"%lf",v+i);

  flat=MFIMFCreateFlat(n,k,o,v,e);
  free(o);
  free(v);

  return flat;
 }

int MFFlatProjectToSave(MFNVector u, double *x, void *d, MFErrorHandler e)
 {
  struct MFFlatData *data;
  int i;

  data=(struct MFFlatData*)d;
  if(x==NULL)return data->n;

  for(i=0;i<data->n;i++)x[i]=MFNV_C(u,i,e);

  return 0;
 }

int MFFlatProjectToDraw(MFNVector u, double *x, void *d, MFErrorHandler e)
 {
  struct MFFlatData *data;
  int i;

  data=(struct MFFlatData*)d;
  if(x==NULL)return data->n;

  for(i=0;i<data->n;i++)x[i]=MFNV_C(u,i,e);

  return 0;
 }

int MFFlatProjectForBB(MFNVector u, double *x, void *d, MFErrorHandler e)
 {
  struct MFFlatData *data;
  int i;

  data=(struct MFFlatData*)d;
  if(x==NULL)return data->n;

  for(i=0;i<data->n;i++)x[i]=MFNV_C(u,i,e);

  return 0;
 }

/*! @} */

/*! @} */

#ifdef __cplusplus
}
#endif
