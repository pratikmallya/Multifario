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

static char *id="@(#) $Id: IMFExpansionPt.c,v 1.7 2011/07/21 17:42:46 mhender Exp $";

static char MFNVectorErrorMsg[256]="";

#include <MFErrorHandler.h>
#include <MFNVector.h>
#include <IMFExpansion.h>
#include <IMFExpansionSpace.h>
#include <IMFExpansionPt.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#ifndef DBL_QNAN
#define DBL_QNAN 1.e200
#endif

#ifdef __cplusplus
 extern "C" {
#endif

struct IMFExpansionNVectorData
 {
  int nC;
  int k;
  IMFExpansion E;
  double t;
  MFNVector sigma;
  MFKVector s0;
  int chart0;
  int prevChart;
  int type;
 };

static void IMFFreeExpansionNVectorData(void*,MFErrorHandler);
static int IMFExpansionNVGetNC(void*,MFErrorHandler);
static double IMFExpansionNVGetC(int,void*,MFErrorHandler);
static void IMFExpansionNVSetC(int,double,void*,MFErrorHandler);
static void IMFExpansionNVDiff(void*,void*,void*,MFErrorHandler);
static void IMFExpansionNVAdd(void*,void*,void*,MFErrorHandler);
static void IMFPrintExpansionNVector(FILE*,void*,MFErrorHandler);
static MFNVector MFCloneExpansionNVector(void*,MFErrorHandler);

MFNVector IMFCreateExpansionNVector(IMFExpansion E,double t, MFNVector sigma, int prevChart, int type, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFCreateExpansionNVector"};
  struct IMFExpansionNVectorData *data;
  int i;
  MFNVector thisExpansion;

  data=(struct IMFExpansionNVectorData*)malloc(sizeof(struct IMFExpansionNVectorData));

#ifndef MFNOSAFETYNET
  if(data==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct IMFExpansionNVectorData));
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return thisExpansion;
   }
#endif

  data->nC=IMFExpansionN(E,e);
  data->k=IMFExpansionK(E,e);
  data->E=E;
  IMFRefExpansion(E,e);

  data->t=t;
  data->sigma=MFCloneNVector(sigma,e);
  data->s0=NULL;
  data->chart0=-1;
  data->prevChart=prevChart;
  data->type=type;

  thisExpansion=MFCreateNVectorBaseClass("IMFExpansionVector",e);

  MFNVectorSetData(thisExpansion,data,e);
  MFNVectorSetFreeData(thisExpansion,IMFFreeExpansionNVectorData,e);
  MFNVectorSetGetNC(thisExpansion,IMFExpansionNVGetNC,e);
  MFNVectorSetGetC(thisExpansion,IMFExpansionNVGetC,e);
  MFNVectorSetSetC(thisExpansion,IMFExpansionNVSetC,e);
  MFNVectorSetDiff(thisExpansion,IMFExpansionNVDiff,e);
  MFNVectorSetAdd(thisExpansion,IMFExpansionNVAdd,e);
  MFNVectorSetClone(thisExpansion,MFCloneExpansionNVector,e);
  MFNVectorSetPrint(thisExpansion,IMFPrintExpansionNVector,e);

  return thisExpansion;
 }

void IMFFreeExpansionNVectorData(void *data, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFFreeExpansionNVectorData"};
  struct IMFExpansionNVectorData *thisExpansion;

  if(data==(void*)0x0087dc00){printf("In %s\n",RoutineName);fflush(stdout);}

  thisExpansion=(struct IMFExpansionNVectorData*)data;

  if(thisExpansion==NULL)return;
  IMFFreeExpansion(thisExpansion->E,e);
  if(thisExpansion->sigma!=NULL)MFFreeNVector(thisExpansion->sigma,e);
  if(thisExpansion->s0!=NULL)MFFreeKVector(thisExpansion->s0,e);
  free(thisExpansion);

  return;
 }

int IMFExpansionNVGetNC(void *data, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFNV_NC"};
  struct IMFExpansionNVectorData *thisExpansion;

  thisExpansion=(struct IMFExpansionNVectorData*)data;

#ifdef MFNOCONFIDENCE
  if(thisExpansion==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector Data (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return -1;
   }
#endif

  return thisExpansion->nC;
 }

double IMFExpansionNVGetC(int i,void *data, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFExpansionNVGetC"};
  struct IMFExpansionNVectorData *thisExpansion;
  double result;

  thisExpansion=(struct IMFExpansionNVectorData*)data;

#ifdef MFNOCONFIDENCE
  if(thisExpansion==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector Data (argument 2) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return DBL_QNAN;
   }

  if(i<0|| !(i<thisExpansion->nC))
   {
    sprintf(MFNVectorErrorMsg,"Coordinate %d (argument 1) is illegal. Must be in 0 to %d",i,thisExpansion->nC-1);
    MFSetError(e,8,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return DBL_QNAN;
   }
#endif

  result=IMFExpansionU(thisExpansion->E,e)[i];

  return result;
 }

void IMFExpansionNVSetC(int i,double vl,void *data, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFExpansionNVSetC"};
  struct IMFExpansionNVectorData *thisExpansion;

  thisExpansion=(struct IMFExpansionNVectorData*)data;

#ifdef MFNOCONFIDENCE
  if(thisExpansion==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector Data (argument 3) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(i<0|| !(i<thisExpansion->nC))
   {
    sprintf(MFNVectorErrorMsg,"Coordinate %d (argument 1) is illegal. Must be in 0 to %d",i,thisExpansion->nC-1);
    MFSetError(e,8,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  IMFExpansionU(thisExpansion->E,e)[i]=vl;

  return;
 }

void IMFExpansionNVDiff(void *adata,void *bdata, void *cdata, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFExpansionNVDiff"};
  int i,n;
  struct IMFExpansionNVectorData *a;
  struct IMFExpansionNVectorData *b;
  struct IMFExpansionNVectorData *c;

  a=(struct IMFExpansionNVectorData*)adata;
  b=(struct IMFExpansionNVectorData*)bdata;
  c=(struct IMFExpansionNVectorData*)cdata;

#ifdef MFNOCONFIDENCE
  if(a==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector Data for a (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(b==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector Data for b (argument 2) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(c==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector Data for c (argument 3) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(a->nC!=b->nC || a->nC!=c->nC || b->nC!=c->nC)
   {
    sprintf(MFNVectorErrorMsg,"Vectors must all be the same length a=%d, b=%d, c=%d",a->nC,b->nC,c->nC);
    MFSetError(e,4,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  n=a->nC;

  for(i=0;i<n;i++)
    IMFExpansionU(c->E,e)[i]=IMFExpansionU(a->E,e)[i]-IMFExpansionU(b->E,e)[i];

  return;
 }

void IMFExpansionNVAdd(void *adata,void *bdata,void *cdata, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFExpansionNVAdd"};
  int i,n;
  struct IMFExpansionNVectorData *a;
  struct IMFExpansionNVectorData *b;
  struct IMFExpansionNVectorData *c;

  a=(struct IMFExpansionNVectorData*)adata;
  b=(struct IMFExpansionNVectorData*)bdata;
  c=(struct IMFExpansionNVectorData*)cdata;

#ifdef MFNOCONFIDENCE
  if(a==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector Data for a (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(b==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector Data for b (argument 2) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(c==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector Data for c (argument 3) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(a->nC!=b->nC || a->nC!=c->nC || b->nC!=c->nC)
   {
    sprintf(MFNVectorErrorMsg,"Vectors must all be the same length a=%d, b=%d, c=%d",a->nC,b->nC,c->nC);
    MFSetError(e,4,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  n=a->nC;

  for(i=0;i<n;i++)
    IMFExpansionU(c->E,e)[i]=IMFExpansionU(a->E,e)[i]-IMFExpansionU(b->E,e)[i];

  return;
 }

void IMFPrintExpansionNVector(FILE *fid,void *data, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFPrintExpansionNVector"};
  int i;
  struct IMFExpansionNVectorData *u;

  u=(struct IMFExpansionNVectorData*)data;

#ifdef MFNOCONFIDENCE
  if(fid==NULL)
   {
    sprintf(MFNVectorErrorMsg,"fid (argument 1) is NULL.");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(data==NULL)
   {
    sprintf(MFNVectorErrorMsg,"data (argument 2) is NULL.");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  if(u->nC<1)return;

  if(u->nC<10)
   {
    fprintf(fid,"(%lf",IMFExpansionU(u->E,e)[0]);
    for(i=1;i<u->nC;i++)
      fprintf(fid,",%lf",IMFExpansionU(u->E,e)[i]);
    fprintf(fid,")");
   }else{
    fprintf(fid,"(%lf",IMFExpansionU(u->E,e)[0]);
    for(i=1;i<5;i++)
      fprintf(fid,",%lf",IMFExpansionU(u->E,e)[i]);
    fprintf(fid,",...");
    for(i=u->nC-5;i<u->nC;i++)
      fprintf(fid,",%lf",IMFExpansionU(u->E,e)[i]);
    fprintf(fid,")");
   }
  fprintf(fid,"\n");

  return;
 }

IMFExpansion IMFExpansionNVGetE(MFNVector thisExpansion, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFExpansionNVGetE"};
  struct IMFExpansionNVectorData *data;

  data=(struct IMFExpansionNVectorData*)MFNVectorGetData(thisExpansion,e);

  return data->E;
 }

double IMFExpansionNVGetT(MFNVector thisExpansion, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFExpansionNVGetE"};
  struct IMFExpansionNVectorData *data;

  data=(struct IMFExpansionNVectorData*)MFNVectorGetData(thisExpansion,e);

  return data->t;
 }

MFNVector IMFExpansionNVGetSigma(MFNVector thisExpansion, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFExpansionNVGetSigma"};
  struct IMFExpansionNVectorData *data;

  data=(struct IMFExpansionNVectorData*)MFNVectorGetData(thisExpansion,e);

  return data->sigma;
 }

MFNVector MFCloneExpansionNVector(void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFCloneExpansionNVector"};
  struct IMFExpansionNVectorData *data;
  IMFExpansion Eclone;
  MFNVector result;
  MFKVector tmp;
  int i;

  data=(struct IMFExpansionNVectorData*)d;

  Eclone=IMFCloneExpansion(data->E,e);
  result=IMFCreateExpansionNVector(Eclone,data->t,data->sigma,data->prevChart,data->type,e);
  IMFFreeExpansion(Eclone,e);
  IMFExpansionNVSetChart0(result,data->chart0,e);
  if(data->s0!=NULL)
    IMFExpansionNVSetS0(result,data->s0,e);

  return result;
 }

int IMFExpansionNVGetType(MFNVector thisExpansion, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFExpansionNVGetE"};
  struct IMFExpansionNVectorData *data;

  data=(struct IMFExpansionNVectorData*)MFNVectorGetData(thisExpansion,e);

  return data->type;
 }

void IMFExpansionNVSetT(MFNVector thisExpansion, double t, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFExpansionNVGetE"};
  struct IMFExpansionNVectorData *data;

  data=(struct IMFExpansionNVectorData*)MFNVectorGetData(thisExpansion,e);

  data->t=t;

  return;
 }

void IMFExpansionNVSetType(MFNVector thisExpansion, int type, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFExpansionNVGetE"};
  struct IMFExpansionNVectorData *data;

  data=(struct IMFExpansionNVectorData*)MFNVectorGetData(thisExpansion,e);

  data->type=type;

  return;
 }

MFKVector IMFExpansionNVGetS0(MFNVector thisExpansion, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFExpansionNVGetS0"};
  struct IMFExpansionNVectorData *data;

  data=(struct IMFExpansionNVectorData*)MFNVectorGetData(thisExpansion,e);

  return data->s0;
 }

void IMFExpansionNVSetS0(MFNVector thisExpansion, MFKVector S, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFExpansionNVSetS0"};
  struct IMFExpansionNVectorData *data;

  data=(struct IMFExpansionNVectorData*)MFNVectorGetData(thisExpansion,e);

  MFRefKVector(S,e);
  if(data->s0!=NULL)MFFreeKVector(data->s0,e);
  data->s0=S;

  return;
 }

int IMFExpansionNVGetChart0(MFNVector thisExpansion, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFExpansionNVGetChart0"};
  struct IMFExpansionNVectorData *data;

  data=(struct IMFExpansionNVectorData*)MFNVectorGetData(thisExpansion,e);

  return data->chart0;
 }

void IMFExpansionNVSetChart0(MFNVector thisExpansion, int chart0, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFExpansionNVSetChart0"};
  struct IMFExpansionNVectorData *data;

  data=(struct IMFExpansionNVectorData*)MFNVectorGetData(thisExpansion,e);

  data->chart0=chart0;
 }

int IMFExpansionNVGetPrev(MFNVector thisExpansion, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFExpansionNVGetPrev"};
  struct IMFExpansionNVectorData *data;

  data=(struct IMFExpansionNVectorData*)MFNVectorGetData(thisExpansion,e);

  return data->prevChart;
 }

void IMFExpansionNVSetPrev(MFNVector thisExpansion, int prevChart, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFExpansionNVSetPrev"};
  struct IMFExpansionNVectorData *data;

  data=(struct IMFExpansionNVectorData*)MFNVectorGetData(thisExpansion,e);

  data->prevChart=prevChart;
 }

#ifdef __cplusplus
}
#endif
