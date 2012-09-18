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

static char *id="@(#) $Id: MFContinuationMethod.c,v 1.3 2007/02/13 01:22:33 mhender Exp $";

char MFContinuationMethodErrorMsg[256];

#include <MFAtlas.h>
#include <MFContinuationMethod.h>
#include <MFErrorHandler.h>
#include <string.h>

#ifdef __cplusplus
 extern "C" {
#endif

struct MFContinuationMethodST{
                              char *type;
                              void *parmBlock;
                              void (*ExtendAtlasMultipleWithTangents)(struct MFContinuationMethodST*,MFAtlas,MFImplicitMF,MFNRegion,int,MFNVector*,MFNKMatrix*,MFErrorHandler);
                              void (*CloseAtlas)(struct MFContinuationMethodST*,MFAtlas,MFErrorHandler);
                              void (*FlushAtlas)(struct MFContinuationMethodST*,MFAtlas,MFErrorHandler);
                              void (*FreeParmBlock)(void*,MFErrorHandler);
                              int nRefs;
                             };

MFContinuationMethod MFCreateContinuationMethodBaseClass(char *type, MFErrorHandler e)
 {
  static char RoutineName[]={"MFCreateContinuationMethodBaseClass"};
  MFContinuationMethod result;

  result=(struct MFContinuationMethodST*)malloc(sizeof(struct MFContinuationMethodST));

#ifndef MFNOSAFETYNET
  if(result==NULL)
   {
    sprintf(MFContinuationMethodErrorMsg,"Out of space trying to allocate %d bytes.",sizeof(struct MFContinuationMethodST));
    MFSetError(e,12,RoutineName,MFContinuationMethodErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  result->type=(char*)malloc((strlen(type)+1)*sizeof(char));

#ifndef MFNOSAFETYNET
  if(result->type==NULL)
   {
    sprintf(MFContinuationMethodErrorMsg,"Out of space trying to allocate %d bytes.",10*sizeof(char));
    MFSetError(e,12,RoutineName,MFContinuationMethodErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(result);
    return NULL;
   }
#endif

  strcpy(result->type,type);

  result->ExtendAtlasMultipleWithTangents=NULL;
  result->CloseAtlas=NULL;
  result->FlushAtlas=NULL;
  result->nRefs=1;

  return result;
 }

void MFFreeContinuationMethod(MFContinuationMethod thisMethod, MFErrorHandler e)
 {
  static char RoutineName[]={"MFFreeContinuationMethod"};

#ifdef MFNOCONFIDENCE
  if(thisMethod==NULL)
   {
    sprintf(MFContinuationMethodErrorMsg,"Pointer to ContinuationMethod (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFContinuationMethodErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  thisMethod->nRefs--;

  if(thisMethod->nRefs<1)
   {
    if(thisMethod->type!=NULL)free(thisMethod->type);
    if(thisMethod->parmBlock!=NULL && thisMethod->FreeParmBlock!=NULL)thisMethod->FreeParmBlock(thisMethod->parmBlock,e);
    free(thisMethod);
   }

  return;
 }

void MFRefContinuationMethod(MFContinuationMethod thisMethod, MFErrorHandler e)
 {
  static char RoutineName[]={"MFRefContinuationMethod"};

#ifdef MFNOCONFIDENCE
  if(thisMethod==NULL)
   {
    sprintf(MFContinuationMethodErrorMsg,"Pointer to ContinuationMethod (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFContinuationMethodErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  thisMethod->nRefs++;

  return;
 }

void MFContinuationMethodExtendAtlasMultipleWithTangents(MFContinuationMethod thisMethod,MFAtlas A,MFImplicitMF M,MFNRegion Omega,int n,MFNVector *u0,MFNKMatrix *Phi0, MFErrorHandler e)
 {
  static char RoutineName[]={"MFContinuationMethodExtendAtlasMultipleWithTangents"};

#ifdef MFNOCONFIDENCE
  if(thisMethod==NULL)
   {
    sprintf(MFContinuationMethodErrorMsg,"Pointer to ContinuationMethod (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFContinuationMethodErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(thisMethod->ExtendAtlasMultipleWithTangents==NULL)
   {
    sprintf(MFContinuationMethodErrorMsg,"This ContinuationMethod has no ComputeAtlasMultipleWithTangents");
    MFSetError(e,12,RoutineName,MFContinuationMethodErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  thisMethod->ExtendAtlasMultipleWithTangents(thisMethod,A,M,Omega,n,u0,Phi0,e);

  return;
 }

void MFFlushAtlas(MFContinuationMethod thisMethod,MFAtlas A, MFErrorHandler e)
 {
  static char RoutineName[]={"MFContinuationMethodFlushAtlas"};

#ifdef MFNOCONFIDENCE
  if(thisMethod==NULL)
   {
    sprintf(MFContinuationMethodErrorMsg,"Pointer to ContinuationMethod (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFContinuationMethodErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  if(thisMethod->FlushAtlas==NULL)return;

  thisMethod->FlushAtlas(thisMethod,A,e);

  return;
 }

void MFCloseAtlas(MFContinuationMethod thisMethod,MFAtlas A, MFErrorHandler e)
 {
  static char RoutineName[]={"MFContinuationMethodCloseAtlas"};

#ifdef MFNOCONFIDENCE
  if(thisMethod==NULL)
   {
    sprintf(MFContinuationMethodErrorMsg,"Pointer to ContinuationMethod (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFContinuationMethodErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  if(thisMethod->CloseAtlas==NULL)return;

  thisMethod->CloseAtlas(thisMethod,A,e);

  return;
 }

MFAtlas MFComputeAtlas(MFContinuationMethod thisMethod, MFImplicitMF M, MFNRegion Omega, MFNVector u0, MFErrorHandler e)
 {
  static char RoutineName[]={"MFComputeAtlas"};
  MFAtlas S;

#ifdef MFNOCONFIDENCE
  if(thisMethod==NULL)
   {
    sprintf(MFContinuationMethodErrorMsg,"Pointer to ContinuationMethod (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFContinuationMethodErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(thisMethod->ExtendAtlasMultipleWithTangents==NULL)
   {
    sprintf(MFContinuationMethodErrorMsg,"This ContinuationMethod has no ComputeAtlas");
    MFSetError(e,12,RoutineName,MFContinuationMethodErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  S=MFCreateAtlas(M,e);

  MFExtendAtlasMultipleWithTangents(S,thisMethod,M,Omega,1,&u0,NULL,e);

  return S;
 }

MFAtlas MFComputeAtlasWithTangent(MFContinuationMethod thisMethod, MFImplicitMF M, MFNRegion Omega, MFNVector u0, MFNKMatrix Phi0, MFErrorHandler e)
 {
  static char RoutineName[]={"MFComputeAtlasWithTangent"};
  MFAtlas S;

#ifdef MFNOCONFIDENCE
  if(thisMethod==NULL)
   {
    sprintf(MFContinuationMethodErrorMsg,"Pointer to ContinuationMethod (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFContinuationMethodErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(thisMethod->ExtendAtlasMultipleWithTangents==NULL)
   {
    sprintf(MFContinuationMethodErrorMsg,"This ContinuationMethod has no ComputeAtlasWithTangent");
    MFSetError(e,12,RoutineName,MFContinuationMethodErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  S=MFCreateAtlas(M,e);

  MFExtendAtlasMultipleWithTangents(S,thisMethod,M,Omega,1,&u0,&Phi0,e);
  
  return S;
 }


MFAtlas MFComputeAtlasMultiple(MFContinuationMethod thisMethod, MFImplicitMF M, MFNRegion Omega, int m, MFNVector *u0, MFErrorHandler e)
 {
  static char RoutineName[]={"MFComputeAtlasMultiple"};
  MFAtlas S;

#ifdef MFNOCONFIDENCE
  if(thisMethod==NULL)
   {
    sprintf(MFContinuationMethodErrorMsg,"Pointer to ContinuationMethod (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFContinuationMethodErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(thisMethod->ExtendAtlasMultipleWithTangents==NULL)
   {
    sprintf(MFContinuationMethodErrorMsg,"This ContinuationMethod has no ComputeAtlasMultiple");
    MFSetError(e,12,RoutineName,MFContinuationMethodErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  S=MFCreateAtlas(M,e);

  MFExtendAtlasMultipleWithTangents(S,thisMethod,M,Omega,m,u0,NULL,e);

  return S;
 }

MFAtlas MFComputeAtlasMultipleWithTangents(MFContinuationMethod thisMethod, MFImplicitMF M, MFNRegion Omega, int m, MFNVector *u0, MFNKMatrix *Phi0, MFErrorHandler e)
 {
  static char RoutineName[]={"MFComputeAtlasMultipleWithTangents"};
  MFAtlas S;

#ifdef MFNOCONFIDENCE
  if(thisMethod==NULL)
   {
    sprintf(MFContinuationMethodErrorMsg,"Pointer to ContinuationMethod (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFContinuationMethodErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(thisMethod->ExtendAtlasMultipleWithTangents==NULL)
   {
    sprintf(MFContinuationMethodErrorMsg,"This ContinuationMethod has no ComputeAtlasMultipleWithTangents");
    MFSetError(e,12,RoutineName,MFContinuationMethodErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  S=MFCreateAtlas(M,e);

  MFExtendAtlasMultipleWithTangents(S,thisMethod,M,Omega,m,u0,Phi0,e);

  return S;
 }

void MFExtendAtlas(MFAtlas S, MFContinuationMethod thisMethod, MFImplicitMF M, MFNRegion Omega, MFNVector u0, MFErrorHandler e)
 {
  static char RoutineName[]={"MFExtendAtlas"};

#ifdef MFNOCONFIDENCE
  if(thisMethod==NULL)
   {
    sprintf(MFContinuationMethodErrorMsg,"Pointer to ContinuationMethod (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFContinuationMethodErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(thisMethod->ExtendAtlasMultipleWithTangents==NULL)
   {
    sprintf(MFContinuationMethodErrorMsg,"This ContinuationMethod has no ExtendAtlas");
    MFSetError(e,12,RoutineName,MFContinuationMethodErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  MFExtendAtlasMultipleWithTangents(S,thisMethod,M,Omega,1,&u0,NULL,e);

  return;
 }

void MFExtendAtlasMultiple(MFAtlas S, MFContinuationMethod thisMethod, MFImplicitMF M, MFNRegion Omega, int m, MFNVector *u0, MFErrorHandler e)
 {
  static char RoutineName[]={"MFExtendAtlasMultiple"};

#ifdef MFNOCONFIDENCE
  if(thisMethod==NULL)
   {
    sprintf(MFContinuationMethodErrorMsg,"Pointer to ContinuationMethod (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFContinuationMethodErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(thisMethod->ExtendAtlasMultipleWithTangents==NULL)
   {
    sprintf(MFContinuationMethodErrorMsg,"This ContinuationMethod has no ExtendAtlasMultiple");
    MFSetError(e,12,RoutineName,MFContinuationMethodErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  MFExtendAtlasMultipleWithTangents(S,thisMethod,M,Omega,m,u0,NULL,e);

  return;
 }

void MFExtendAtlasWithTangent(MFAtlas S, MFContinuationMethod thisMethod, MFImplicitMF M, MFNRegion Omega, MFNVector u0,MFNKMatrix Phi0, MFErrorHandler e)
 {
  static char RoutineName[]={"MFExtendAtlasWithTangent"};

#ifdef MFNOCONFIDENCE
  if(thisMethod==NULL)
   {
    sprintf(MFContinuationMethodErrorMsg,"Pointer to ContinuationMethod (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFContinuationMethodErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(thisMethod->ExtendAtlasMultipleWithTangents==NULL)
   {
    sprintf(MFContinuationMethodErrorMsg,"This ContinuationMethod has no ExtendAtlasWithTangent");
    MFSetError(e,12,RoutineName,MFContinuationMethodErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  MFExtendAtlasMultipleWithTangents(S,thisMethod,M,Omega,1,&u0,&Phi0,e);

  return;
 }

void MFExtendAtlasMultipleWithTangents(MFAtlas S, MFContinuationMethod thisMethod, MFImplicitMF M, MFNRegion Omega, int m, MFNVector *u0, MFNKMatrix *Phi0, MFErrorHandler e)
 {
  static char RoutineName[]={"MFExtendAtlasMultipleWithTangents"};

#ifdef MFNOCONFIDENCE
  if(thisMethod==NULL)
   {
    sprintf(MFContinuationMethodErrorMsg,"Pointer to ContinuationMethod (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFContinuationMethodErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(thisMethod->ExtendAtlasMultipleWithTangents==NULL)
   {
    sprintf(MFContinuationMethodErrorMsg,"This ContinuationMethod has no ExtendAtlasWithTangent");
    MFSetError(e,12,RoutineName,MFContinuationMethodErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  thisMethod->ExtendAtlasMultipleWithTangents(thisMethod,S,M,Omega,m,u0,Phi0,e);

  return;
 }

char *MFContinuationMethodGetType(MFContinuationMethod thisMethod, MFErrorHandler e)
 {
  static char RoutineName[]={"MFContinuationMethodGetType"};

#ifdef MFNOCONFIDENCE
  if(thisMethod==NULL)
   {
    sprintf(MFContinuationMethodErrorMsg,"Pointer to ContinuationMethod (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFContinuationMethodErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  return thisMethod->type;
 }

void *MFContinuationMethodGetParmBlock(MFContinuationMethod thisMethod, MFErrorHandler e)
 {
  static char RoutineName[]={"MFContinuationMethodGetParmBlock"};

#ifdef MFNOCONFIDENCE
  if(thisMethod==NULL)
   {
    sprintf(MFContinuationMethodErrorMsg,"Pointer to ContinuationMethod (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFContinuationMethodErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  return thisMethod->parmBlock;
 }

void MFContinuationMethodSetParmBlock(MFContinuationMethod thisMethod, void *parmBlock, MFErrorHandler e)
 {
  static char RoutineName[]={"MFContinuationMethodSetParmBlock"};

#ifdef MFNOCONFIDENCE
  if(thisMethod==NULL)
   {
    sprintf(MFContinuationMethodErrorMsg,"Pointer to ContinuationMethod (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFContinuationMethodErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  thisMethod->parmBlock=parmBlock;
  return;
 }

void MFContinuationMethodSetFreeParmBlock(MFContinuationMethod thisMethod, void (*FreeParmBlock)(void*,MFErrorHandler), MFErrorHandler e)
 {
  static char RoutineName[]={"MFContinuationMethodSetFreeParmBlock"};

#ifdef MFNOCONFIDENCE
  if(thisMethod==NULL)
   {
    sprintf(MFContinuationMethodErrorMsg,"Pointer to ContinuationMethod (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFContinuationMethodErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  thisMethod->FreeParmBlock=FreeParmBlock;
  return;
 }

void MFContinuationMethodSetExtendAtlasMultipleWithTangents(MFContinuationMethod thisMethod, void (*ExtendAtlasMultipleWithTangents)(struct MFContinuationMethodST*,MFAtlas,MFImplicitMF,MFNRegion,int,MFNVector*,MFNKMatrix*,MFErrorHandler), MFErrorHandler e)
 {
  static char RoutineName[]={"MFContinuationMethodSetExtendAtlasMultipleWithTangents"};

#ifdef MFNOCONFIDENCE
  if(thisMethod==NULL)
   {
    sprintf(MFContinuationMethodErrorMsg,"Pointer to ContinuationMethod (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFContinuationMethodErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  thisMethod->ExtendAtlasMultipleWithTangents=ExtendAtlasMultipleWithTangents;
  return;
 }

void MFContinuationMethodSetCloseAtlas(MFContinuationMethod thisMethod, void (*CloseAtlas)(struct MFContinuationMethodST*,MFAtlas,MFErrorHandler), MFErrorHandler e)
 {
  static char RoutineName[]={"MFContinuationMethodSetCloseAtlas"};

#ifdef MFNOCONFIDENCE
  if(thisMethod==NULL)
   {
    sprintf(MFContinuationMethodErrorMsg,"Pointer to ContinuationMethod (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFContinuationMethodErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  thisMethod->CloseAtlas=CloseAtlas;
  return;
 }

void MFContinuationMethodSetFlushAtlas(MFContinuationMethod thisMethod, void (*FlushAtlas)(struct MFContinuationMethodST*,MFAtlas,MFErrorHandler), MFErrorHandler e)
 {
  static char RoutineName[]={"MFContinuationMethodSetFlushAtlas"};

#ifdef MFNOCONFIDENCE
  if(thisMethod==NULL)
   {
    sprintf(MFContinuationMethodErrorMsg,"Pointer to ContinuationMethod (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFContinuationMethodErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  thisMethod->FlushAtlas=FlushAtlas;
  return;
 }


#ifdef __cplusplus
}
#endif
