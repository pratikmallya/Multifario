/* 
 *  PROGRAM NAME:  multifario
 *
 *  (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
 *  CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
 *
 *  Please refer to the LICENSE file in the top directory
 *
 *
 *      author: Mike Henderson mhender@watson.ibm.com
 *      date:   November 11, 1997
 *              February 2, 1999   Ported to C
 */

static char *id="@(#) $Id: MFDenseNVector.c,v 1.6 2011/07/21 17:42:46 mhender Exp $";

#include <MFErrorHandler.h>
#include <MFNVector.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef DBL_QNAN
#define DBL_QNAN 1.e300
#endif

#ifdef __cplusplus
 extern "C" {
#endif

struct MFDenseNVectorData
 {
  int nC;
  double *data;
  int wrapped;
 };

void MFFreeDenseNVectorData(void*,MFErrorHandler);
int MFDenseNVGetNC(void*,MFErrorHandler);
double MFDenseNVGetC(int,void*,MFErrorHandler);
void MFDenseNVSetC(int,double,void*,MFErrorHandler);
void MFDenseNVDiff(void*,void*,void*,MFErrorHandler);
void MFDenseNVAdd(void*,void*,void*,MFErrorHandler);
void MFWriteDenseNVectorData(FILE*,void*,MFErrorHandler);
MFNVector MFReadDenseNVector(FILE*,MFErrorHandler);
void MFPrintDenseNVector(FILE*,void*,MFErrorHandler);
MFNVector MFCloneDenseNVector(void*,MFErrorHandler);

static char MFNVectorErrorMsg[256]="";

MFNVector MFCreateNVectorWithData(int n,double *vl, MFErrorHandler e)
 {
  static char RoutineName[]={"MFCreateNVectorWithData"};
  struct MFDenseNVectorData *data;
  int i;
  MFNVector thisVector;

#ifdef MFNOCONFIDENCE
  if(n<1)
   {
    sprintf(MFNVectorErrorMsg,"Length of Vector %d (argument 1) is Illegal. Must be positive.",n);
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return NULL;
   }
#endif

  data=(struct MFDenseNVectorData*)malloc(sizeof(struct MFDenseNVectorData)); /*done*/

#ifndef MFNOSAFETYNET
  if(data==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct MFDenseNVectorData));
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return thisVector;
   }
#endif

  data->nC=n;
  data->data=(double*)malloc(n*sizeof(double));


#ifndef MFNOSAFETYNET
  if(n>0&&data->data==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(double));
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return thisVector;
   }
#endif

  if(n>0&&vl!=NULL)
   {
    for(i=0;i<n;i++)(data->data)[i]=vl[i];
   }else if(vl==NULL)
   {
    for(i=0;i<n;i++)(data->data)[i]=0.;
   }

  data->wrapped=0;

  thisVector=MFCreateNVectorBaseClass("DENSE",e);

  MFNVectorSetData(thisVector,data,e);
  MFNVectorSetFreeData(thisVector,MFFreeDenseNVectorData,e);
  MFNVectorSetWriteData(thisVector,MFWriteDenseNVectorData,e);
  MFNVectorSetGetNC(thisVector,MFDenseNVGetNC,e);
  MFNVectorSetGetC(thisVector,MFDenseNVGetC,e);
  MFNVectorSetSetC(thisVector,MFDenseNVSetC,e);
  MFNVectorSetDiff(thisVector,MFDenseNVDiff,e);
  MFNVectorSetAdd(thisVector,MFDenseNVAdd,e);
  MFNVectorSetClone(thisVector,MFCloneDenseNVector,e);
  MFNVectorSetPrint(thisVector,MFPrintDenseNVector,e);

  return thisVector;
 }

MFNVector MFCreateWrappedNVector(int n,double *vl, MFErrorHandler e)
 {
  static char RoutineName[]={"MFCreateWrappedNVector"};
  struct MFDenseNVectorData *data;
  int i;
  MFNVector thisVector;

  thisVector=NULL;

#ifdef MFNOCONFIDENCE
  if(n<1)
   {
    sprintf(MFNVectorErrorMsg,"Length of Vector %d (argument 1) is Illegal. Must be positive.",n);
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return NULL;
   }

  if(vl==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to data (argument 2) is Illegal. Must be non-NULL.");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return thisVector;
   }
#endif

  data=(struct MFDenseNVectorData*)malloc(sizeof(struct MFDenseNVectorData));

#ifndef MFNOSAFETYNET
  if(data==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct MFDenseNVectorData));
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return thisVector;
   }
#endif

  data->nC=n;
  data->data=vl;
  data->wrapped=1;

  thisVector=MFCreateNVectorBaseClass("DENSE",e);

  MFNVectorSetData(thisVector,data,e);
  MFNVectorSetFreeData(thisVector,MFFreeDenseNVectorData,e);
  MFNVectorSetWriteData(thisVector,MFWriteDenseNVectorData,e);
  MFNVectorSetGetNC(thisVector,MFDenseNVGetNC,e);
  MFNVectorSetGetC(thisVector,MFDenseNVGetC,e);
  MFNVectorSetSetC(thisVector,MFDenseNVSetC,e);
  MFNVectorSetDiff(thisVector,MFDenseNVDiff,e);
  MFNVectorSetAdd(thisVector,MFDenseNVAdd,e);
  MFNVectorSetClone(thisVector,MFCloneDenseNVector,e);
  MFNVectorSetPrint(thisVector,MFPrintDenseNVector,e);

  return thisVector;
 }

MFNVector MFCreateNVector(int n, MFErrorHandler e)
 {
  static char RoutineName[]={"MFCreateNVector"};
  int i;
  MFNVector thisVector;

  thisVector=MFCreateNVectorWithData(n,NULL,e);

  return thisVector;
 }

void MFFreeDenseNVectorData(void *data, MFErrorHandler e)
 {
  static char RoutineName[]={"MFFreeDenseNVectorData"};
  struct MFDenseNVectorData *thisVector;
  int verbose=0;

  thisVector=(struct MFDenseNVectorData*)data;

  if(verbose){printf("%s line %d, data=0x%8.8x\n",RoutineName,__LINE__,data);fflush(stdout);}

  if(thisVector==NULL)return;
  if(!thisVector->wrapped&&thisVector->data!=NULL)
   {
    if(verbose){printf("%s line %d, data->data=0x%8.8x\n",RoutineName,__LINE__,thisVector->data);fflush(stdout);}
    free(thisVector->data);
   }
  if(verbose){printf("%s line %d\n",RoutineName,__LINE__);fflush(stdout);}
  free(thisVector);
  if(verbose){printf("%s line %d\n",RoutineName,__LINE__);fflush(stdout);}

  return;
 }

int MFDenseNVGetNC(void *data, MFErrorHandler e)
 {
  static char RoutineName[]={"MFNV_NC"};
  struct MFDenseNVectorData *thisVector;

  thisVector=(struct MFDenseNVectorData*)data;

#ifdef MFNOCONFIDENCE
  if(thisVector==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector Data (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return -1;
   }
#endif

  return thisVector->nC;
 }

double MFDenseNVGetC(int i,void *data, MFErrorHandler e)
 {
  static char RoutineName[]={"MFDenseNVGetC"};
  struct MFDenseNVectorData *thisVector;
  double result;

  thisVector=(struct MFDenseNVectorData*)data;

#ifdef MFNOCONFIDENCE
  if(thisVector==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector Data (argument 2) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return DBL_QNAN;
   }

  if(i<0|| !(i<thisVector->nC))
   {
    sprintf(MFNVectorErrorMsg,"Coordinate %d (argument 1) is illegal. Must be in 0 to %d",i,thisVector->nC-1);
    MFSetError(e,8,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return DBL_QNAN;
   }
#endif

  result=thisVector->data[i];

  return result;
 }

void MFDenseNVSetC(int i,double vl,void *data, MFErrorHandler e)
 {
  static char RoutineName[]={"MFDenseNVSetC"};
  struct MFDenseNVectorData *thisVector;

  thisVector=(struct MFDenseNVectorData*)data;

#ifdef MFNOCONFIDENCE
  if(thisVector==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector Data (argument 3) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(i<0|| !(i<thisVector->nC))
   {
    sprintf(MFNVectorErrorMsg,"Coordinate %d (argument 1) is illegal. Must be in 0 to %d",i,thisVector->nC-1);
    MFSetError(e,8,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  thisVector->data[i]=vl;

  return;
 }

void MFDenseNVDiff(void *adata,void *bdata, void *cdata, MFErrorHandler e)
 {
  static char RoutineName[]={"MFDenseNVDiff"};
  int i,n;
  struct MFDenseNVectorData *a;
  struct MFDenseNVectorData *b;
  struct MFDenseNVectorData *c;

  a=(struct MFDenseNVectorData*)adata;
  b=(struct MFDenseNVectorData*)bdata;
  c=(struct MFDenseNVectorData*)cdata;

#ifdef MFNOCONFIDENCE
  if(a==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector Data for a (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  b=(struct MFDenseNVectorData*)bdata;
  if(b==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector Data for b (argument 2) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  c=(struct MFDenseNVectorData*)cdata;
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
    c->data[i]=a->data[i]-b->data[i];

  return;
 }

void MFDenseNVAdd(void *adata,void *bdata,void *cdata, MFErrorHandler e)
 {
  static char RoutineName[]={"MFDenseNVAdd"};
  int i,n;
  struct MFDenseNVectorData *a;
  struct MFDenseNVectorData *b;
  struct MFDenseNVectorData *c;

  a=(struct MFDenseNVectorData*)adata;
  b=(struct MFDenseNVectorData*)bdata;
  c=(struct MFDenseNVectorData*)cdata;

#ifdef MFNOCONFIDENCE
  if(a==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector Data for a (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  b=(struct MFDenseNVectorData*)bdata;
  if(b==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector Data for b (argument 2) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  c=(struct MFDenseNVectorData*)cdata;
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

  for(i=0;i<n;i++)c->data[i]=a->data[i]+b->data[i];

  return;
 }

double *MFNV_CStar(MFNVector thisVector, MFErrorHandler e)
 {
  static char RoutineName[]={"MFNV_CStar"};
  struct MFDenseNVectorData *u;

#ifdef MFNOCONFIDENCE
  if(thisVector==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return NULL;
   }

  if(strcmp(MFNVGetId(thisVector,e),"DENSE"))
   {
    sprintf(MFNVectorErrorMsg,"Trying to get dense arrays from non-dense Vector type \"%s\"",MFNVGetId(thisVector,e));
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return NULL;
   }
#endif

  u=(struct MFDenseNVectorData*)MFNVectorGetData(thisVector,e);

  return u->data;
 }

void MFWriteDenseNVectorData(FILE *fid,void *data, MFErrorHandler e)
 {
  static char RoutineName[]={"MFWriteDenseNVectorData"};
  int i;
  struct MFDenseNVectorData *u;

  u=(struct MFDenseNVectorData*)data;

#ifdef MFNOCONFIDENCE
  if(u==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector Data (argument 2) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  fprintf(fid,"%d\n",u->nC);
  for(i=0;i<u->nC;i++)
   {
    if(i>0)fprintf(fid," ");
    fprintf(fid,"%lf",u->data[i]);
   }
  fprintf(fid,"\n");

  return;
 }

MFNVector MFReadDenseNVector(FILE *fid, MFErrorHandler e)
 {
  static char RoutineName[]={"MFReadDenseNVector"};
  int n;
  double *vl;
  MFNVector thisVector;
  int i;

  fscanf(fid,"%d\n",&n);

  vl=(double*)malloc(n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(n>0&&vl==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(double));
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  for(i=0;i<n;i++)
   {
    if(i>0)fscanf(fid," ");
    fscanf(fid,"%lf",&(vl[i]));
   }
  fscanf(fid,"\n");

  thisVector=MFCreateNVectorWithData(n,vl,e);

  free(vl);

  return thisVector;
 }

MFNVector MFCloneDenseNVector(void *data, MFErrorHandler e)
 {
  static char RoutineName[]={"MFCloneDenseNVector"};
  struct MFDenseNVectorData *u;

  u=(struct MFDenseNVectorData*)data;

  return MFCreateNVectorWithData(u->nC,u->data,e);
 }

void MFPrintDenseNVector(FILE *fid,void *data, MFErrorHandler e)
 {
  static char RoutineName[]={"MFPrintDenseNVector"};
  int i;
  struct MFDenseNVectorData *u;

  u=(struct MFDenseNVectorData*)data;

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
    fprintf(fid,"(%lf",(u->data)[0]);
    for(i=1;i<u->nC;i++)
      fprintf(fid,",%lf",(u->data)[i]);
    fprintf(fid,")");
   }else{
    fprintf(fid,"(%lf",(u->data)[0]);
    for(i=1;i<5;i++)
      fprintf(fid,",%lf",(u->data)[i]);
    fprintf(fid,",...");
    for(i=u->nC-5;i<u->nC;i++)
      fprintf(fid,",%lf",(u->data)[i]);
    fprintf(fid,")");
   }

  return;
 }

#ifdef __cplusplus
}
#endif
