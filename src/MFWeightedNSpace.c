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

static char *id="@(#) $Id: MFWeightedNSpace.c,v 1.3 2007/02/13 01:22:34 mhender Exp $";

static char MFWeightedNSpaceErrorMsg[256]="";

#include <MFNSpace.h>
#include <MFNVector.h>
#include <MFErrorHandler.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef __cplusplus
 extern "C" {
#endif

struct MFWeightedNSpaceData
 {
  int n;
  int l;
  double s;
  int *e;
 };

static void MFFreeWeightedNSpaceData(void*,MFErrorHandler);
static double MFWeightedNSpaceDistance(MFNSpace,MFNVector,MFNVector,void*,MFErrorHandler);
static void MFWeightedNSpaceDirection(MFNSpace,MFNVector,MFNVector,MFNVector,void*,MFErrorHandler);
static void MFWeightedNSpaceAdd(MFNSpace,MFNVector,MFNVector,MFNVector,void*,MFErrorHandler);
static void MFWeightedNSpaceScale(MFNSpace,double,MFNVector,MFNVector,void*,MFErrorHandler);
static double MFWeightedNSpaceInner(MFNSpace,MFNVector,MFNVector,void*,MFErrorHandler);
static void MFWriteWeightedNSpaceData(FILE*,MFNSpace,void*,MFErrorHandler);
void MFReadWeightedNSpaceData(FILE*,MFNSpace,MFErrorHandler);

MFNSpace MFCreateNSpace(int n, MFErrorHandler e)
 {
  static char RoutineName[]={"MFCreateNSpace"};
  MFNSpace thisSpace;
  struct MFWeightedNSpaceData *data;

  thisSpace=MFCreateNSpaceBaseClass("WeightedNSpace",e);
  MFNSpaceSetDistance(thisSpace,MFWeightedNSpaceDistance,e);
  MFNSpaceSetInnerProduct(thisSpace,MFWeightedNSpaceInner,e);
  MFNSpaceSetDirection(thisSpace,MFWeightedNSpaceDirection,e);
  MFNSpaceSetAdd(thisSpace,MFWeightedNSpaceAdd,e);
  MFNSpaceSetScale(thisSpace,MFWeightedNSpaceScale,e);
  MFNSpaceSetFreeData(thisSpace,MFFreeWeightedNSpaceData,e);
  MFNSpaceSetWriteData(thisSpace,MFWriteWeightedNSpaceData,e);

  data=(struct MFWeightedNSpaceData*)malloc(sizeof(struct MFWeightedNSpaceData));

#ifndef MFNOSAFETYNET
  if(data==NULL)
   {
    sprintf(MFWeightedNSpaceErrorMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct MFWeightedNSpaceData));
    MFSetError(e,12,RoutineName,MFWeightedNSpaceErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif
  
  data->n=n;
  data->l=0;
  data->s=1.;
  data->e=(int*)NULL;

  MFNSpaceSetData(thisSpace,(void*)data,e);

  return thisSpace;
 }

MFNSpace MFCreateWeightedNSpace(int n,int l,double s, MFErrorHandler e)
 {
  static char RoutineName[]={"MFCreateWeightedNSpace"};
  MFNSpace thisSpace;
  struct MFWeightedNSpaceData *data;

  thisSpace=MFCreateNSpaceBaseClass("WeightedNSpace",e);
  MFNSpaceSetDistance(thisSpace,MFWeightedNSpaceDistance,e);
  MFNSpaceSetInnerProduct(thisSpace,MFWeightedNSpaceInner,e);
  MFNSpaceSetDirection(thisSpace,MFWeightedNSpaceDirection,e);
  MFNSpaceSetAdd(thisSpace,MFWeightedNSpaceAdd,e);
  MFNSpaceSetScale(thisSpace,MFWeightedNSpaceScale,e);
  MFNSpaceSetFreeData(thisSpace,MFFreeWeightedNSpaceData,e);
  MFNSpaceSetWriteData(thisSpace,MFWriteWeightedNSpaceData,e);

  data=(struct MFWeightedNSpaceData*)malloc(sizeof(struct MFWeightedNSpaceData));

#ifndef MFNOSAFETYNET
  if(data==NULL)
   {
    sprintf(MFWeightedNSpaceErrorMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct MFWeightedNSpaceData));
    MFSetError(e,12,RoutineName,MFWeightedNSpaceErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  data->n=n;
  data->l=l;
  data->s=s;
  data->e=(int*)NULL;

  MFNSpaceSetData(thisSpace,(void*)data,e);

  return thisSpace;
 }

MFNSpace MFCreateNSpaceWithExponents(int n, int *exp, MFErrorHandler e)
 {
  static char RoutineName[]={"MFCreateNSpaceWithExponents"};
  MFNSpace thisSpace;
  struct MFWeightedNSpaceData *data;
  int i;

  thisSpace=MFCreateNSpaceBaseClass("WeightedNSpace",e);
  MFNSpaceSetDistance(thisSpace,MFWeightedNSpaceDistance,e);
  MFNSpaceSetInnerProduct(thisSpace,MFWeightedNSpaceInner,e);
  MFNSpaceSetDirection(thisSpace,MFWeightedNSpaceDirection,e);
  MFNSpaceSetAdd(thisSpace,MFWeightedNSpaceAdd,e);
  MFNSpaceSetScale(thisSpace,MFWeightedNSpaceScale,e);
  MFNSpaceSetFreeData(thisSpace,MFFreeWeightedNSpaceData,e);
  MFNSpaceSetWriteData(thisSpace,MFWriteWeightedNSpaceData,e);

  data=(struct MFWeightedNSpaceData*)malloc(sizeof(struct MFWeightedNSpaceData));

#ifndef MFNOSAFETYNET
  if(data==NULL)
   {
    sprintf(MFWeightedNSpaceErrorMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct MFWeightedNSpaceData));
    MFSetError(e,12,RoutineName,MFWeightedNSpaceErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  data->n=n;
  data->l=0;
  data->s=1.;
  for(i=0;i<n;i++)data->e[i]=exp[i];

  MFNSpaceSetData(thisSpace,(void*)data,e);

  return thisSpace;
 }

MFNSpace MFCreateWeightedNSpaceWithExponents(int n,int l,double s, int *exp, MFErrorHandler e)
 {
  static char RoutineName[]={"MFCreateWeightedNSpaceWithExponents"};
  MFNSpace thisSpace;
  struct MFWeightedNSpaceData *data;
  int i;

  thisSpace=MFCreateNSpaceBaseClass("WeightedNSpace",e);
  MFNSpaceSetDistance(thisSpace,MFWeightedNSpaceDistance,e);
  MFNSpaceSetInnerProduct(thisSpace,MFWeightedNSpaceInner,e);
  MFNSpaceSetDirection(thisSpace,MFWeightedNSpaceDirection,e);
  MFNSpaceSetAdd(thisSpace,MFWeightedNSpaceAdd,e);
  MFNSpaceSetScale(thisSpace,MFWeightedNSpaceScale,e);
  MFNSpaceSetFreeData(thisSpace,MFFreeWeightedNSpaceData,e);
  MFNSpaceSetWriteData(thisSpace,MFWriteWeightedNSpaceData,e);

  data=(struct MFWeightedNSpaceData*)malloc(sizeof(struct MFWeightedNSpaceData));

#ifndef MFNOSAFETYNET
  if(data==NULL)
   {
    sprintf(MFWeightedNSpaceErrorMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct MFWeightedNSpaceData));
    MFSetError(e,12,RoutineName,MFWeightedNSpaceErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  data->n=n;
  data->n=n;
  data->l=l;
  data->s=s;
  for(i=0;i<n;i++)data->e[i]=exp[i];

  MFNSpaceSetData(thisSpace,(void*)data,e);

  return thisSpace;
 }

void MFFreeWeightedNSpaceData(void *data, MFErrorHandler e)
 {
  static char RoutineName[]={"MFFreeWeightedNSpaceData"};
  struct MFWeightedNSpaceData *d;

  d=(struct MFWeightedNSpaceData*)data;
  if(d!=(struct MFWeightedNSpaceData*)NULL)free(d);

  return;
 }

double MFWeightedNSpaceDistance(MFNSpace thisSpace,MFNVector v0,MFNVector v1,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFWeightedNSpaceDistance"};
  double result;
  int i;
  struct MFWeightedNSpaceData *data;

  data=(struct MFWeightedNSpaceData*)d;

  if(data->e==(int*)NULL)
   {
    result=0.;
    for(i=0;i<data->l;i++)result+=pow(MFNV_C(v0,i,e)-MFNV_C(v1,i,e),2)*data->s;
    for(i=data->l;i<data->n;i++)result+=pow(MFNV_C(v0,i,e)-MFNV_C(v1,i,e),2);
    result=sqrt(result);
   }else{
    result=0.;
    for(i=0;i<data->l;i++)
      result+=pow(MFNV_C(v0,i,e)-MFNV_C(v1,i,e),2*data->e[i])*data->s;
    for(i=data->l;i<data->n;i++)
      result+=pow(MFNV_C(v0,i,e)-MFNV_C(v1,i,e),2*data->e[i]);
    result=sqrt(result);
   }

  return result;
 }

void MFWeightedNSpaceDirection(MFNSpace thisSpace,MFNVector v0,MFNVector v1,MFNVector diff,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFWeightedNSpaceDirection"};
  int i;
  struct MFWeightedNSpaceData *data;
  data=(struct MFWeightedNSpaceData*)d;

  for(i=0;i<data->n;i++)MFNVSetC(diff,i,MFNV_C(v1,i,e)-MFNV_C(v0,i,e),e);

  return;
 }

double MFWeightedNSpaceInner(MFNSpace thisSpace,MFNVector v0,MFNVector v1,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFWeightedNSpaceInner"};
  double result;
  int i;
  int verbose=0;
  struct MFWeightedNSpaceData *data;

  data=(struct MFWeightedNSpaceData*)d;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s 0x%8.8x\n",RoutineName,thisSpace);fflush(stdout);}
#endif

  if(data->e==(int*)NULL)
   {

#ifdef MFALLOWVERBOSE
    if(verbose)
     {
      printf("%s - e==NULL\n",RoutineName);fflush(stdout);
      printf("s==%lf\n",data->s);fflush(stdout);
      printf("l==%d\n",data->l);fflush(stdout);
     }
#endif

    result=0.;
    for(i=0;i<data->l;i++)
     {
      result+=MFNV_C(v0,i,e)*MFNV_C(v1,i,e)*data->s;

#ifdef MFALLOWVERBOSE
      if(verbose){printf("   %lf * %lf * %lf += %lf\n",MFNV_C(v0,i,e),MFNV_C(v1,i,e),data->s,result);fflush(stdout);}
#endif

     }
    for(i=data->l;i<data->n;i++)
     {
      result+=MFNV_C(v0,i,e)*MFNV_C(v1,i,e);

#ifdef MFALLOWVERBOSE
      if(verbose){printf("   %lf * %lf += %lf\n",MFNV_C(v0,i,e),MFNV_C(v1,i,e),result);fflush(stdout);}
#endif

     }
   }else{

#ifdef MFALLOWVERBOSE
    if(verbose)
     {
      printf("%s - e!=NULL\n",RoutineName);fflush(stdout);
      printf("s==%lf\n",data->s);fflush(stdout);
      printf("l==%d\n",data->l);fflush(stdout);
     }
#endif
    result=0.;
    for(i=0;i<data->l;i++)
      result+=pow(MFNV_C(v0,i,e),data->e[i])*pow(MFNV_C(v1,i,e),data->e[i])*data->s;
    for(i=data->l;i<data->n;i++)
      result+=pow(MFNV_C(v0,i,e)*MFNV_C(v1,i,e),data->e[i]);
   }

  return result;
 }

void MFWriteWeightedNSpaceData(FILE *fid,MFNSpace thisSpace,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFWriteWeightedNSpaceData"};
  int i;
  struct MFWeightedNSpaceData *data;

  data=(struct MFWeightedNSpaceData*)d;

  fprintf(fid,"%d %d %lf\n",data->n,data->l,data->s);
  if(data->e==(int*)NULL)
   {
    fprintf(fid,"0\n");
   }else{
    fprintf(fid,"1\n");
    for(i=0;i<data->n;i++)
      fprintf(fid,"%d\n",data->e[i]);
   }

  return;
 }

void MFReadWeightedNSpaceData(FILE *fid,MFNSpace thisSpace, MFErrorHandler e)
 {
  static char RoutineName[]={"MFReadWeightedNSpaceData"};
  int i,n;
  char tag[100]="";
  struct MFWeightedNSpaceData *data;

  MFNSpaceSetDistance(thisSpace,MFWeightedNSpaceDistance,e);
  MFNSpaceSetInnerProduct(thisSpace,MFWeightedNSpaceInner,e);
  MFNSpaceSetDirection(thisSpace,MFWeightedNSpaceDirection,e);
  MFNSpaceSetAdd(thisSpace,MFWeightedNSpaceAdd,e);
  MFNSpaceSetScale(thisSpace,MFWeightedNSpaceScale,e);
  MFNSpaceSetFreeData(thisSpace,MFFreeWeightedNSpaceData,e);
  MFNSpaceSetWriteData(thisSpace,MFWriteWeightedNSpaceData,e);

  data=(struct MFWeightedNSpaceData*)malloc(sizeof(struct MFWeightedNSpaceData));

#ifndef MFNOSAFETYNET
  if(data==NULL)
   {
    sprintf(MFWeightedNSpaceErrorMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct MFWeightedNSpaceData));
    MFSetError(e,12,RoutineName,MFWeightedNSpaceErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  fscanf(fid,"%d %d %lf\n",&(data->n),&(data->l),&(data->s));
  fscanf(fid,"%d\n",&n);
  if(n==0)
    data->e=(int*)NULL;
   else{
    data->e=(int*)malloc(data->n*sizeof(int));

#ifndef MFNOSAFETYNET
  if(data->e==NULL)
   {
    sprintf(MFWeightedNSpaceErrorMsg,"Out of memory, trying to allocate %d bytes",data->n*sizeof(int));
    MFSetError(e,12,RoutineName,MFWeightedNSpaceErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

    for(i=0;i<data->n;i++)
              fscanf(fid,"%d\n",&(data->e[i]));
   }
  MFNSpaceSetData(thisSpace,(void*)data,e);

  return;
 }

void MFWeightedNSpaceAdd(MFNSpace thisSpace,MFNVector v0,MFNVector v1,MFNVector sum,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFWeightedNSpaceAdd"};
  int i;
  struct MFWeightedNSpaceData *data;

  data=(struct MFWeightedNSpaceData*)d;

  for(i=0;i<data->n;i++)MFNVSetC(sum,i,MFNV_C(v1,i,e)+MFNV_C(v0,i,e),e);

  return;
 }

void MFWeightedNSpaceScale(MFNSpace thisSpace,double s, MFNVector v,MFNVector prod,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFWeightedNSpaceScale"};
  int i;
  struct MFWeightedNSpaceData *data;

  data=(struct MFWeightedNSpaceData*)d;

  for(i=0;i<data->n;i++)MFNVSetC(prod,i,s*MFNV_C(v,i,e),e);

  return;
 }


#ifdef __cplusplus
}
#endif
