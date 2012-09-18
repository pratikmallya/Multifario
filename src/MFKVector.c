/* 
 *  PROGRAM NAME:  multifario
 *
 *  (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
 *  CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
 *
 *  Please refer to the LICENSE file in the top directory
 *
 *      author: Mike Henderson mhender@watson.ibm.com
 *      date:   November 11, 1997
 *              February 2, 1999   Ported to C
 */

static char *id="@(#) $Id: MFKVector.c,v 1.4 2008/05/02 12:44:18 mhender Exp $";

#include <MFErrorHandler.h>
#include <MFKVector.h>
#include <float.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifndef DBL_QNAN
#define DBL_QNAN 1.e200
#endif

#ifdef __cplusplus
 extern "C" {
#endif

void MFSetError(MFErrorHandler,int,char*,char*,int,char*);
static char MFKVectorErrorMsg[256]="";

struct MFKVectorSt
 {
  int nC;
  double *data;
  int nRefs;
 };

MFKVector MFCreateKVector(int n, MFErrorHandler e)
 {
  static char RoutineName[]={"MFCreateKVector"};
  MFKVector thisVector;
  int i;

#ifdef MFNOCONFIDENCE
  if(n<1)
   {
    sprintf(MFKVectorErrorMsg,"Length of Vector %d (argument 1) is Illegal. Must be positive.",n);
    MFSetError(e,12,RoutineName,MFKVectorErrorMsg,__LINE__,__FILE__);
    return NULL;
   }
#endif

  thisVector=(MFKVector)malloc(sizeof(struct MFKVectorSt));

#ifdef MFNOSAFETYNET
  if(thisVector==NULL)
   {
    sprintf(MFKVectorErrorMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct MFKVectorSt));
    MFSetError(e,12,RoutineName,MFKVectorErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  thisVector->nC=n;

  thisVector->data=(double*)malloc(n*sizeof(double));

#ifdef MFNOSAFETYNET
  if(thisVector->data==NULL)
   {
    sprintf(MFKVectorErrorMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(double));
    MFSetError(e,12,RoutineName,MFKVectorErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  for(i=0;i<n;i++)thisVector->data[i]=0.;
  thisVector->nRefs=1;

  return thisVector;
 }

MFKVector MFCreateKVectorWithData(int n,double *vl, MFErrorHandler e)
 {
  static char RoutineName[]={"MFCreateKVectorWithFullData"};
  int i;
  MFKVector thisVector;

  if(n<1)
   {
    sprintf(MFKVectorErrorMsg,"Length of Vector %d (argument 1) is Illegal. Must be positive.",n);
    MFSetError(e,4,RoutineName,MFKVectorErrorMsg,__LINE__,__FILE__);
    return NULL;
   }

  thisVector=(MFKVector)malloc(sizeof(struct MFKVectorSt));

#ifdef MFNOSAFETYNET
  if(thisVector==NULL)
   {
    sprintf(MFKVectorErrorMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct MFKVectorSt));
    MFSetError(e,12,RoutineName,MFKVectorErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return thisVector;
   }
#endif

  thisVector->data=(double*)malloc(n*sizeof(double));

#ifdef MFNOSAFETYNET
  if(thisVector->data==NULL)
   {
    sprintf(MFKVectorErrorMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(double));
    MFSetError(e,12,RoutineName,MFKVectorErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return thisVector;
   }
#endif

#ifdef MFNOCONFIDENCE
  if(vl==NULL)
   {
    sprintf(MFKVectorErrorMsg,"The pointer to the array of coordinates (argument 2) is NULL");
    MFSetError(e,12,RoutineName,MFKVectorErrorMsg,__LINE__,__FILE__);
    thisVector->data=NULL;
    thisVector->nRefs=1;

    return thisVector;
   }
#endif

  for(i=0;i<n;i++)thisVector->data[i]=vl[i];
  thisVector->nC=n;

  thisVector->nRefs=1;

  return thisVector;
 }

void MFFreeKVector(MFKVector thisVector, MFErrorHandler e)
 {
  static char RoutineName[]={"MFFreeKVector"};

#ifdef MFNOCONFIDENCE
  if(thisVector==NULL)
   {
    sprintf(MFKVectorErrorMsg,"Pointer to Vector (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFKVectorErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  thisVector->nRefs--;

  if(thisVector->nRefs<1)
   {
    if(thisVector->data!=NULL){free(thisVector->data);thisVector->data=NULL;}
    free(thisVector);
   }
  return;
 }

int MFKV_NC(MFKVector thisVector, MFErrorHandler e)
 {
  static char RoutineName[]={"MFKV_NC"};

#ifdef MFNOCONFIDENCE
  if(thisVector==NULL)
   {
    sprintf(MFKVectorErrorMsg,"Pointer to Vector (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFKVectorErrorMsg,__LINE__,__FILE__);
    return -1;
   }
#endif

  return thisVector->nC;
 }

double MFKV_C(MFKVector thisVector,int i, MFErrorHandler e)
 {
  static char RoutineName[]={"MFKV_C"};
  int j;

#ifdef MFNOCONFIDENCE
  if(thisVector==NULL)
   {
    sprintf(MFKVectorErrorMsg,"Pointer to Vector (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFKVectorErrorMsg,__LINE__,__FILE__);
    return DBL_QNAN;
   }
  if(i<0|| !(i<thisVector->nC))
   {
    sprintf(MFKVectorErrorMsg,"Coordinate %d (argument 2) is illegal. Must be in 0 to %d",i,thisVector->nC-1);
    MFSetError(e,8,RoutineName,MFKVectorErrorMsg,__LINE__,__FILE__);
    return DBL_QNAN;
   }
#endif

  return thisVector->data[i];
 }

void MFKVSetC(MFKVector thisVector,int i,double vl, MFErrorHandler e)
 {
  static char RoutineName[]={"MFKVSetC"};
  int j;

#ifdef MFNOCONFIDENCE
  if(thisVector==NULL)
   {
    sprintf(MFKVectorErrorMsg,"Pointer to Vector (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFKVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(i<0|| !(i<thisVector->nC))
   {
    sprintf(MFKVectorErrorMsg,"Coordinate %d (argument 2) is illegal. Must be in 0 to %d",i,thisVector->nC-1);
    MFSetError(e,8,RoutineName,MFKVectorErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  thisVector->data[i]=vl;
  return;
 }

void MFRefKVector(MFKVector thisVector, MFErrorHandler e)
 {
  static char RoutineName[]={"MFRefKVector"};

#ifdef MFNOCONFIDENCE
  if(thisVector==NULL)
   {
    sprintf(MFKVectorErrorMsg,"Pointer to Vector (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFKVectorErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  thisVector->nRefs++;
  return;
 }

double MFKVDot(MFKVector a,MFKVector b, MFErrorHandler e)
 {
  static char RoutineName[]={"MFKVDot"};
  double d;
  int i;

#ifdef MFNOCONFIDENCE
  if(a==NULL)
   {
    sprintf(MFKVectorErrorMsg,"Vector a (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFKVectorErrorMsg,__LINE__,__FILE__);
    return DBL_QNAN;
   }

  if(b==NULL)
   {
    sprintf(MFKVectorErrorMsg,"Vector b (argument 2) is NULL");
    MFSetError(e,12,RoutineName,MFKVectorErrorMsg,__LINE__,__FILE__);
    return DBL_QNAN;
   }

  if(a->nC!=b->nC)
   {
    sprintf(MFKVectorErrorMsg,"Vectors must be the same length a=%d, b=%d",a->nC,b->nC);
    MFSetError(e,12,RoutineName,MFKVectorErrorMsg,__LINE__,__FILE__);
    return DBL_QNAN;
   }
#endif

  d=0;
  for(i=0;i<a->nC;i++)d+=a->data[i]*b->data[i];

  return d;
 }

double MFKVNorm(MFKVector a, MFErrorHandler e)
 {
  static char RoutineName[]={"MFKVNorm"};
  double d;
  int i;

#ifdef MFNOCONFIDENCE
  if(a==NULL)
   {
    sprintf(MFKVectorErrorMsg,"Vector a (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFKVectorErrorMsg,__LINE__,__FILE__);
    return -1.;
   }
#endif

  d=0;
  for(i=0;i<a->nC;i++)d+=a->data[i]*a->data[i];

  return sqrt(d);
 }

void MFKVScale(double s,MFKVector a, MFErrorHandler e)
 {
  static char RoutineName[]={"MFKVScale"};
  int i;

  for(i=0;i<a->nC;i++)a->data[i]*=s;

  return;
 }

void MFKVScaleMul(double s,MFKVector a,MFKVector b, MFErrorHandler e)
 {
  static char RoutineName[]={"MFKVScaleMul"};
  int i;


#ifdef MFNOCONFIDENCE
  if(a==NULL)
   {
    sprintf(MFKVectorErrorMsg,"Vector a (argument 2) is NULL");
    MFSetError(e,12,RoutineName,MFKVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(b==NULL)
   {
    sprintf(MFKVectorErrorMsg,"Vector a (argument 3) is NULL");
    MFSetError(e,12,RoutineName,MFKVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(a->nC!=b->nC)
   {
    sprintf(MFKVectorErrorMsg,"Vectors must be the same length a=%d, b=%d",a->nC,b->nC);
    MFSetError(e,12,RoutineName,MFKVectorErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  for(i=0;i<a->nC;i++)b->data[i]=s*a->data[i];

  return;
 }

double *MFKV_CStar(MFKVector thisVector, MFErrorHandler e)
 {
  static char RoutineName[]={"MFKV_CStar"};


#ifdef MFNOCONFIDENCE
  if(thisVector==NULL)
   {
    sprintf(MFKVectorErrorMsg,"Pointer to Vector (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFKVectorErrorMsg,__LINE__,__FILE__);
    return NULL;
   }
#endif

  return thisVector->data;
 }

void MFKVDiff(MFKVector a,MFKVector b, MFKVector c, MFErrorHandler e)
 {
  static char RoutineName[]={"MFKVDiff"};
  MFKVector result;
  int i,n;

#ifdef MFNOCONFIDENCE
  if(a==NULL)
   {
    sprintf(MFKVectorErrorMsg,"Vector a (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFKVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(b==NULL)
   {
    sprintf(MFKVectorErrorMsg,"Vector b (argument 2) is NULL");
    MFSetError(e,12,RoutineName,MFKVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(c==NULL)
   {
    sprintf(MFKVectorErrorMsg,"Vector c (argument 3) is NULL");
    MFSetError(e,12,RoutineName,MFKVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(a->nC!=b->nC || a->nC!=c->nC || b->nC!=c->nC)
   {
    sprintf(MFKVectorErrorMsg,"Vectors must all be the same length a=%d, b=%d, c=%d",a->nC,b->nC,c->nC);
    MFSetError(e,12,RoutineName,MFKVectorErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  n=a->nC;

  for(i=0;i<n;i++)
    c->data[i]=a->data[i]-b->data[i];

  return;
 }

void MFKVAdd(MFKVector a,MFKVector b, MFKVector c, MFErrorHandler e)
 {
  static char RoutineName[]={"MFKVAdd"};
  MFKVector result;
  int i,n;

#ifdef MFNOCONFIDENCE
  if(a==NULL)
   {
    sprintf(MFKVectorErrorMsg,"Vector a (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFKVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(b==NULL)
   {
    sprintf(MFKVectorErrorMsg,"Vector b (argument 2) is NULL");
    MFSetError(e,12,RoutineName,MFKVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(c==NULL)
   {
    sprintf(MFKVectorErrorMsg,"Vector c (argument 3) is NULL");
    MFSetError(e,12,RoutineName,MFKVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(a->nC!=b->nC || a->nC!=c->nC || b->nC!=c->nC)
   {
    sprintf(MFKVectorErrorMsg,"Vectors must all be the same length a=%d, b=%d, c=%d",a->nC,b->nC,c->nC);
    MFSetError(e,12,RoutineName,MFKVectorErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  n=a->nC;

  for(i=0;i<n;i++)
    c->data[i]=a->data[i]+b->data[i];

  return;
 }

void MFWriteKVector(FILE *fid,MFKVector s, MFErrorHandler e)
 {
  static char RoutineName[]={"MFWriteKVector"};
  int i;

  fprintf(fid,"%s\n","KVector");
  fprintf(fid,"%d %d\n",s->nC,s->nRefs);
  for(i=0;i<s->nC;i++)
   {
    if(i>0)fprintf(fid," ");
    fprintf(fid,"%lf",s->data[i]);
   }
  fprintf(fid,"\n");

  return;
 }

MFKVector MFReadKVector(FILE *fid, MFErrorHandler e)
 {
  static char RoutineName[]={"MFReadKVector"};
  int i;
  MFKVector s;
  char tag[100]="";

  fscanf(fid,"%s\n",tag);

#ifdef MFNOCONFIDENCE
  if(strcmp(tag,"KVector"))
   {
    sprintf(MFKVectorErrorMsg,"Next Object is not a KVector! (%s)\n",RoutineName,tag);
    MFSetError(e,12,RoutineName,MFKVectorErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  s=(MFKVector)malloc(sizeof(struct MFKVectorSt));
#ifdef MFNOSAFETYNET
  if(s==NULL)
   {
    sprintf(MFKVectorErrorMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct MFKVectorSt));
    MFSetError(e,12,RoutineName,MFKVectorErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  fscanf(fid,"%d %d\n",&(s->nC),&(s->nRefs));

  s->data=(double*)malloc(s->nC*sizeof(double)); /*done*/

#ifdef MFNOSAFETYNET
  if(s->data==NULL)
   {
    sprintf(MFKVectorErrorMsg,"Out of memory, trying to allocate %d bytes",s->nC*sizeof(double));
    MFSetError(e,12,RoutineName,MFKVectorErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif
  
  for(i=0;i<s->nC;i++)
   {
    if(i>0)fscanf(fid," ");
    fscanf(fid,"%lf",&(s->data[i]));
   }
  fscanf(fid,"\n");

  return s;
 }

int MFKVectorGetNRefs(MFKVector u, MFErrorHandler e)
 {
  static char RoutineName[]={"MFKVectorGetNRefs"};

  return u->nRefs;
 }

#ifdef __cplusplus
}
#endif
