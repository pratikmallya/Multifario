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
 *              April 29, 2002     Made into a base class
 */

static char *id="@(#) $Id: MFNVector.c,v 1.7 2011/07/21 17:42:46 mhender Exp $";

#include <multifarioConfig.h>
#include <MFErrorHandler.h>
#include <MFNVector.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef DBL_QNAN
#define DBL_QNAN 1.e200
#endif

static char MFNVectorErrorMsg[256]="";

#ifdef __cplusplus
 extern "C" {
#endif

MFNVector MFReadDenseNVector(FILE*);
MFNVector MFReadAUTOBVNVector(FILE*);

struct MFNVectorSt
 {
  char *type;

  void *data;
  void (*freedata)(void*,MFErrorHandler);
  void (*writedata)(FILE*,void*,MFErrorHandler);

  MFNVector (*clone)(void*,MFErrorHandler);
  int (*getNC)(void*,MFErrorHandler);
  double (*getC)(int,void*,MFErrorHandler);
  void (*setC)(int,double,void*,MFErrorHandler);
  void (*add)(void*,void*,void*,MFErrorHandler);
  void (*diff)(void*,void*,void*,MFErrorHandler);
  void (*print)(FILE*,void*,MFErrorHandler);

  int index;
  int index2;
  int nRefs;
 };

MFNVector MFCreateNVectorBaseClass(char *type, MFErrorHandler e)
 {
  static char RoutineName[]={"MFCreateNVectorBaseClass"};
  MFNVector thisVector;
  int i;

  thisVector=(MFNVector)malloc(sizeof(struct MFNVectorSt)); /*done*/
  if(thisVector==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct MFNVectorSt));
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return thisVector;
   }

  thisVector->data=NULL;
  thisVector->freedata=NULL;
  thisVector->writedata=NULL;
  thisVector->clone=NULL;
  thisVector->getNC=NULL;
  thisVector->getC=NULL;
  thisVector->setC=NULL;
  thisVector->add=NULL;
  thisVector->diff=NULL;
  thisVector->print=NULL;

  thisVector->nRefs=1;

  thisVector->index=0;
  thisVector->index2=0;

  if(type!=NULL)
   {
    thisVector->type=(char*)malloc((strlen(type)+1)*sizeof(char));

#ifdef MFNOCONFIDENCE
    if(thisVector->type==NULL)
     {
      sprintf(MFNVectorErrorMsg,"Out of memory, trying to allocate %d bytes",(strlen(type)+1)*sizeof(char));
      MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

    strcpy(thisVector->type,type);
   }else
   thisVector->type=NULL;

  return thisVector;
 }

void MFFreeNVector(MFNVector thisVector, MFErrorHandler e)
 {
  static char RoutineName[]={"MFFreeNVector"};
  int i;

#ifdef MFNOCONFIDENCE
  if(thisVector==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  thisVector->nRefs--;

  if(thisVector->nRefs<1)
   {
    if(thisVector->type!=NULL){free(thisVector->type);thisVector->type=NULL;}

    if(thisVector->freedata!=NULL&&thisVector->data!=NULL){(thisVector->freedata)(thisVector->data,e);thisVector->data=NULL;}

    free(thisVector);
   }
  return;
 }

int MFNV_NC(MFNVector thisVector, MFErrorHandler e)
 {
  static char RoutineName[]={"MFNV_NC"};
  int result;

#ifdef MFNOCONFIDENCE
  if(thisVector==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return -1;
   }

  if(thisVector->getNC==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Routine to get the number of coordinates has not been provided");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return -1;
   }
#endif

  result=(thisVector->getNC)(thisVector->data,e);

  return result;
 }

double MFNV_C(MFNVector thisVector,int i, MFErrorHandler e)
 {
  static char RoutineName[]={"MFNV_C"};
  int j;
  double result;

#ifdef MFNOCONFIDENCE
  if(thisVector==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return DBL_QNAN;
   }

  if(thisVector->getC==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Routine to get a coordinate (vector type %s) has not been provided",thisVector->type);
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return DBL_QNAN;
   }
#endif

  result=(thisVector->getC)(i,thisVector->data,e);

  return result;
 }

void MFNVSetC(MFNVector thisVector,int i,double vl, MFErrorHandler e)
 {
  static char RoutineName[]={"MFNVSetC"};
  int rc;
  int j;

#ifdef MFNOCONFIDENCE
  if(thisVector==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(thisVector->setC==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Routine to set a coordinate (vector type %s) has not been provided",thisVector->type);
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  (thisVector->setC)(i,vl,thisVector->data,e);

  return;
 }

void MFPrintNVector(FILE *fid, MFNVector thisVector, MFErrorHandler e)
 {
  static char RoutineName[]={"MFPrintNVector"};
  int j;
  double result;

#ifdef MFNOCONFIDENCE
  if(thisVector==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(thisVector->print==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Routine to print an NVector of type %s has not been provided",thisVector->print);
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  (thisVector->print)(fid,thisVector->data,e);

  return;
 }

void MFRefNVector(MFNVector thisVector, MFErrorHandler e)
 {
  static char RoutineName[]={"MFRefNVector"};

#ifdef MFNOCONFIDENCE
  if(thisVector==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  thisVector->nRefs++;
  return;
 }

void MFNVDiff(MFNVector a,MFNVector b, MFNVector c, MFErrorHandler e)
 {
  static char RoutineName[]={"MFNVDiff"};

#ifdef MFNOCONFIDENCE
  if(a==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Vector a (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(b==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Vector b (argument 2) is NULL");
    MFSetError(e,4,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(c==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Vector c (argument 3) is NULL");
    MFSetError(e,4,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(a->diff==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Routine to take the difference of two vectors not been provided");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(strcmp(a->type,b->type))
   {
    sprintf(MFNVectorErrorMsg,"Can't take the difference of vectors of different type!");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(strcmp(b->type,c->type))
   {
    sprintf(MFNVectorErrorMsg,"The difference of two vectors must be the same type as the vectors!");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  (a->diff)(a->data,b->data,c->data,e);

  return;
 }

void MFNVAdd(MFNVector a,MFNVector b, MFNVector c, MFErrorHandler e)
 {
  static char RoutineName[]={"MFNVAdd"};

#ifdef MFNOCONFIDENCE
  if(a==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Vector a (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(b==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Vector b (argument 2) is NULL");
    MFSetError(e,4,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(c==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Vector c (argument 3) is NULL");
    MFSetError(e,4,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(a->add==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Routine to add two vectors not been provided");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(strcmp(a->type,b->type))
   {
    sprintf(MFNVectorErrorMsg,"Can't add vectors of different type!");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(strcmp(b->type,c->type))
   {
    sprintf(MFNVectorErrorMsg,"The sum of two vectors must be the same type as the vectors!");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  (a->add)(a->data,b->data,c->data,e);

  return;
 }

void MFWriteNVector(FILE *fid,MFNVector u, MFErrorHandler e)
 {
  static char RoutineName[]={"MFWriteNVector"};
  int i;

  fprintf(fid,"%s\n","NVector");
  if(u->type!=NULL)
    fprintf(fid,"%d %s\n",strlen(u->type),u->type);
   else
    fprintf(fid,"%d \n",0);
  if(u->writedata!=NULL) (u->writedata)(fid,u->data,e);
  fprintf(fid,"%d %d %d\n",u->index,u->index2,u->nRefs);

  return;
 }

MFNVector MFReadNVector(FILE *fid, MFErrorHandler e)
 {
  static char RoutineName[]={"MFReadNVector"};
  int n;
  MFNVector u;
  char tag[100]="";

  fscanf(fid,"%s\n",tag);

#ifdef MFNOCONFIDENCE
  if(strcmp(tag,"NVector"))
   {
    sprintf(MFNVectorErrorMsg,"Next Object is not a NVector! (%s)\n",RoutineName,tag);
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  u=MFCreateNVectorBaseClass(NULL,e);

  fscanf(fid,"%d",&n);
  if(n>0)
   {
    fscanf(fid,"%s\n",tag);
#ifdef MFNOCONFIDENCE
   }else{
    sprintf(MFNVectorErrorMsg,"Unknown input Vector type, type length 0");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return NULL;
#endif
   }

  if(!strcmp(tag,"DENSE"))u=MFReadDenseNVector(fid);
   else if(!strcmp(u->type,"LOCA"))n=n;
   else if(!strcmp(u->type,"AUTO"))u=MFReadAUTOBVNVector(fid);
   else{
    sprintf(MFNVectorErrorMsg,"Unknown input Vector type \"%s\"",u->type);
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return NULL;
   }

  fscanf(fid,"%d %d %d\n",&(u->index),&(u->index2),&(u->nRefs));

  return u;
 }

char *MFNVGetType(MFNVector u, MFErrorHandler e)
 {
  static char RoutineName[]={"MFNVGetType"};

#ifdef MFNOCONFIDENCE
  if(u==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  return u->type;
 }

int MFNVGetIndex(MFNVector u, MFErrorHandler e)
 {
  static char RoutineName[]={"MFNVGetIndex"};

#ifdef MFNOCONFIDENCE
  if(u==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  return u->index;
 }

void MFNVSetIndex(MFNVector u, int index, MFErrorHandler e)
 {
  static char RoutineName[]={"MFNVSetIndex"};

  u->index=index;
 }

int MFNVGetIndex2(MFNVector u, MFErrorHandler e)
 {
  static char RoutineName[]={"MFNVGetIndex2"};

#ifdef MFNOCONFIDENCE
  if(u==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  return u->index2;
 }

void MFNVSetIndex2(MFNVector u, int index, MFErrorHandler e)
 {
  static char RoutineName[]={"MFNVSetIndex2"};

#ifdef MFNOCONFIDENCE
  if(u==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  u->index2=index;
 }

MFNVector MFCloneNVector(MFNVector thisVector, MFErrorHandler e)
 {
  static char RoutineName[]={"MFCloneNVector"};

#ifdef MFNOCONFIDENCE
  if(thisVector==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return NULL;
   }

  if(thisVector->clone==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Clone routine not provided");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return NULL;
   }
#endif

  return (thisVector->clone)(thisVector->data,e);
 }

char *MFNVGetId(MFNVector thisVector, MFErrorHandler e)
 {
  static char RoutineName[]={"MFNVGetId"};

#ifdef MFNOCONFIDENCE
  if(thisVector==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return NULL;
   }
#endif

  return thisVector->type;
 }

void MFNVectorSetPrint(MFNVector thisVector, void (*print)(FILE*,void*,MFErrorHandler), MFErrorHandler e)
 {
  static char RoutineName[]={"MFNVectorSetPrint"};

#ifdef MFNOCONFIDENCE
  if(thisVector==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  thisVector->print=print;

  return;
 }

void MFNVectorSetData(MFNVector thisVector, void *data, MFErrorHandler e)
 {
  static char RoutineName[]={"MFNVectorSetData"};

#ifdef MFNOCONFIDENCE
  if(thisVector==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  thisVector->data=data;

  return;
 }

void *MFNVectorGetData(MFNVector thisVector, MFErrorHandler e)
 {
  static char RoutineName[]={"MFNVectorGetData"};

#ifdef MFNOCONFIDENCE
  if(thisVector==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  return thisVector->data;
 }

void MFNVectorSetFreeData(MFNVector thisVector, void (*freedata)(void*,MFErrorHandler), MFErrorHandler e)
 {
  static char RoutineName[]={"MFNVectorSetFreeData"};

#ifdef MFNOCONFIDENCE
  if(thisVector==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  thisVector->freedata=freedata;

  return;
 }

void MFNVectorSetWriteData(MFNVector thisVector, void (*writedata)(FILE*,void*,MFErrorHandler), MFErrorHandler e)
 {
  static char RoutineName[]={"MFNVectorSetWriteData"};

#ifdef MFNOCONFIDENCE
  if(thisVector==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  thisVector->writedata=writedata;

  return;
 }

void MFNVectorSetClone(MFNVector thisVector, MFNVector (*clone)(void*,MFErrorHandler), MFErrorHandler e)
 {
  static char RoutineName[]={"MFNVectorSetClone"};

#ifdef MFNOCONFIDENCE
  if(thisVector==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  thisVector->clone=clone;

  return;
 }

void MFNVectorSetGetNC(MFNVector thisVector, int (*getNC)(void*,MFErrorHandler), MFErrorHandler e)
 {
  static char RoutineName[]={"MFNVectorSetGetNC"};

#ifdef MFNOCONFIDENCE
  if(thisVector==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  thisVector->getNC=getNC;

  return;
 }

void MFNVectorSetGetC(MFNVector thisVector, double (*getC)(int,void*,MFErrorHandler), MFErrorHandler e)
 {
  static char RoutineName[]={"MFNVectorSetGetC"};

#ifdef MFNOCONFIDENCE
  if(thisVector==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  thisVector->getC=getC;

  return;
 }

void MFNVectorSetSetC(MFNVector thisVector, void (*setC)(int,double,void*,MFErrorHandler), MFErrorHandler e)
 {
  static char RoutineName[]={"MFNVectorSetSetC"};

#ifdef MFNOCONFIDENCE
  if(thisVector==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  thisVector->setC=setC;

  return;
 }

void MFNVectorSetAdd(MFNVector thisVector, void (*add)(void*,void*,void*,MFErrorHandler), MFErrorHandler e)
 {
  static char RoutineName[]={"MFNVectorSetAdd"};

#ifdef MFNOCONFIDENCE
  if(thisVector==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  thisVector->add=add;

  return;
 }

void MFNVectorSetDiff(MFNVector thisVector, void (*diff)(void*,void*,void*,MFErrorHandler), MFErrorHandler e)
 {
  static char RoutineName[]={"MFNVectorSetDiff"};

#ifdef MFNOCONFIDENCE
  if(thisVector==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  thisVector->diff=diff;

  return;
 }

#ifdef __cplusplus
}
#endif
