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

static char *id="@(#) $Id: MFNSpace.c,v 1.5 2007/07/06 12:38:36 mhender Exp $";

#include <multifarioConfig.h>

static char MFNSpaceErrorMsg[256]="";

#include <MFNSpace.h>
#include <MFNVector.h>
#include <MFNKMatrix.h>
#include <MFErrorHandler.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#ifdef __cplusplus
 extern "C" {
#endif

struct MFNSpaceSt
 {
  int nRefs;
  char *id;
  void *d;

  double (*distance)(MFNSpace,MFNVector,MFNVector,void*,MFErrorHandler);
  double (*inner)(MFNSpace,MFNVector,MFNVector,void*,MFErrorHandler);
  void (*direction)(MFNSpace,MFNVector,MFNVector,MFNVector,void*,MFErrorHandler);
  void (*add)(MFNSpace,MFNVector,MFNVector,MFNVector,void*,MFErrorHandler);
  void (*scale)(MFNSpace,double,MFNVector,MFNVector,void*,MFErrorHandler);
  void (*writedata)(FILE*,MFNSpace,void*,MFErrorHandler);
  void (*freedata)(void*,MFErrorHandler);
 };

void MFRefNSpace(MFNSpace thisSpace, MFErrorHandler e)
 {
  static char RoutineName[]={"MFRefNSpace"};

  thisSpace->nRefs++;

  return;
 }

void MFFreeNSpace(MFNSpace thisSpace, MFErrorHandler e)
 {
  static char RoutineName[]={"MFFreeNSpace"};

  thisSpace->nRefs--;
  if(thisSpace->nRefs<1)
   {
    if(thisSpace->id!=NULL)free(thisSpace->id);
    if(thisSpace->freedata!=NULL && thisSpace->d!=NULL)thisSpace->freedata(thisSpace->d,e);
    free(thisSpace);
   }
  return;
 }

double MFNSpaceDistance(MFNSpace thisSpace,MFNVector v0,MFNVector v1, MFErrorHandler e)
 {
  static char RoutineName[]={"MFNSpaceDistance"};
  double result;
  int i;

  result=thisSpace->distance(thisSpace,v0,v1,thisSpace->d,e);

  return result;
 }

void MFNSpaceDirection(MFNSpace thisSpace,MFNVector v0,MFNVector v1,MFNVector diff, MFErrorHandler e)
 {
  static char RoutineName[]={"MFNSpaceDirection"};
  double d;
  int i;

  thisSpace->direction(thisSpace,v0,v1,diff,thisSpace->d,e);

  return;
 }

double MFNSpaceInner(MFNSpace thisSpace,MFNVector v0,MFNVector v1, MFErrorHandler e)
 {
  static char RoutineName[]={"MFNSpaceInner"};
  double result;
  int i;

  result=thisSpace->inner(thisSpace,v0,v1,thisSpace->d,e);

  return result;
 }

void MFWriteNSpace(FILE *fid,MFNSpace u, MFErrorHandler e)
 {
  static char RoutineName[]={"MFWriteNSpace"};
  int i;

  fprintf(fid,"%s\n","NSpace");
  fprintf(fid,"%d\n",strlen(u->id));
  fprintf(fid,"%s\n",u->id);
  fprintf(fid,"%d\n",u->nRefs);
  u->writedata(fid,u,u->d,e);

  return;
 }

void MFReadWeightedNSpaceData(FILE *fid,MFNSpace,MFErrorHandler);
#ifdef MFTPBVP
void MFReadTPBVPNSpaceData(FILE *fid,MFNSpace,MFErrorHandler);
#endif
/*void MFReadAUTONSpaceData(FILE *fid,MFNSpace,MFErrorHandler);*/

MFNSpace MFReadNSpace(FILE *fid, MFErrorHandler e)
 {
  static char RoutineName[]={"MFReadNSpace"};
  int i;
  MFNSpace Omega;
  char tag[100]="";
  int n=0;
  char *id;
  int verbose=0;

  fscanf(fid,"%s\n",tag);
  if(strcmp(tag,"NSpace"))
   {
    sprintf(MFNSpaceErrorMsg,"Next Object is not a NSpace! (%s) -- %s\n",RoutineName,tag);
    MFSetError(e,12,RoutineName,MFNSpaceErrorMsg,__LINE__,__FILE__);
    return NULL;
   }
  if(verbose)
   {
    printf("%s\n",RoutineName);
    printf("  tag %s\n",tag);fflush(stdout);
   }

  fscanf(fid,"%d\n",&n);
  if(verbose){printf("  len(id) %d\n",n);fflush(stdout);}

  id=(char*)malloc((n+1)*sizeof(char));

#ifndef MFNOSAFETYNET
  if(id==NULL)
   {
    sprintf(MFNSpaceErrorMsg,"Out of memory, trying to allocate %d bytes",(n+1)*sizeof(char));
    MFSetError(e,12,RoutineName,MFNSpaceErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  fscanf(fid,"%s\n",id);
  if(verbose){printf("  id %s\n",id);fflush(stdout);}

  Omega= MFCreateNSpaceBaseClass(id,e);
  fscanf(fid,"%d\n",&(Omega->nRefs));

  if(!strcmp(id,"WeightedNSpace"))
   {
    MFReadWeightedNSpaceData(fid,Omega,e);

#ifdef MFTPBVP
  }else if(!strcmp(id,"TPBVPNSpace"))
   {
    MFReadTPBVPNSpaceData(fid,Omega,e);
#endif

  }else if(!strcmp(id,"AUTONSpace"))
   {
/*  MFReadAUTONSpaceData(fid,Omega,e);*/

   }else{
    sprintf(MFNSpaceErrorMsg,"Object is not a known NSpace! (%s)\n",RoutineName,tag);
    MFSetError(e,12,RoutineName,MFNSpaceErrorMsg,__LINE__,__FILE__);
    free(id);
    return NULL;
   }
  free(id);

  return Omega;
 }

void MFNSpaceAdd(MFNSpace thisSpace,MFNVector v0,MFNVector v1,MFNVector sum, MFErrorHandler e)
 {
  static char RoutineName[]={"MFNSpaceAdd"};
  double d;
  int i;

#ifdef MFNOCONFIDENCE
  if(thisSpace==NULL)
   {
    sprintf(MFNSpaceErrorMsg,"NSpace, (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNSpaceErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(thisSpace->add==NULL)
   {
    sprintf(MFNSpaceErrorMsg,"NSpace, type %s. add called, but not set.",thisSpace->id);
    MFSetError(e,12,RoutineName,MFNSpaceErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  thisSpace->add(thisSpace,v0,v1,sum,thisSpace->d,e);

  return;
 }

double MFNSpaceTangentDistance(MFNSpace thisSpace, MFNKMatrix Phi0, MFNKMatrix Phi1, MFErrorHandler e)
 {
  static char RoutineName[]={"MFNSpaceTangentDistance"};
  int i,j,n,k;
  double id[9];
  MFNVector phi0,phi1;

  k=MFNKMatrixK(Phi0,e);
  n=MFNKMatrixN(Phi0,e);

  if(k>2)
   {

#ifdef MFNOCONFIDENCE
    sprintf(MFNSpaceErrorMsg,"Routine only implemented for k<=2, thisSpace test will not be used.");
    MFSetError(e,4,RoutineName,MFNSpaceErrorMsg,__LINE__,__FILE__);
#endif

    return 0.;
   }


  for(i=0;i<k;i++)
   {
    phi0=MFMColumn(Phi0,i,e);
    for(j=0;j<k;j++)
     {
      phi1=MFMColumn(Phi1,j,e);
      id[i+k*j]=MFNSpaceInner(thisSpace,phi0,phi1,e);
      MFFreeNVector(phi1,e);
     }
    MFFreeNVector(phi0,e);
   }

  if(k==1)return fabs(id[0]-1.);
   else if(k==2)return sqrt(fabs(id[0]*id[3]-id[1]*id[2]-1.));

  return 0.;
 }

void MFNSpaceScale(MFNSpace thisSpace,double s,MFNVector v,MFNVector prod, MFErrorHandler e)
 {
  static char RoutineName[]={"MFNSpaceScale"};
  double d;
  int i;

#ifdef MFNOCONFIDENCE
    if(thisSpace==NULL)
     {
      sprintf(MFNSpaceErrorMsg,"NSpace, (argument 1) is NULL");
      MFSetError(e,12,RoutineName,MFNSpaceErrorMsg,__LINE__,__FILE__);
      return;
     }

    if(thisSpace->add==NULL)
     {
      sprintf(MFNSpaceErrorMsg,"NSpace, type %s. scale called, but not set.",thisSpace->id);
      MFSetError(e,12,RoutineName,MFNSpaceErrorMsg,__LINE__,__FILE__);
      return;
     }
#endif

  thisSpace->scale(thisSpace,s,v,prod,thisSpace->d,e);

  return;
 }

MFNSpace MFCreateNSpaceBaseClass(char *id, MFErrorHandler e)
 {
  static char RoutineName[]={"MFCreateNSpaceBaseClass"};
  MFNSpace thisSpace;

  thisSpace=(MFNSpace)malloc(sizeof(struct MFNSpaceSt));

#ifndef MFNOSAFETYNET
  if(thisSpace==NULL)
   {
    sprintf(MFNSpaceErrorMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct MFNSpaceSt));
    MFSetError(e,12,RoutineName,MFNSpaceErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  thisSpace->id=(char*)malloc((strlen(id)+1)*sizeof(char));

#ifndef MFNOSAFETYNET
  if(thisSpace->id==NULL)
   {
    sprintf(MFNSpaceErrorMsg,"Out of memory, trying to allocate %d bytes",(strlen(id)+1)*sizeof(char));
    MFSetError(e,12,RoutineName,MFNSpaceErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  strcpy(thisSpace->id,id);
  thisSpace->nRefs=1;
  thisSpace->d=NULL;
  thisSpace->distance=NULL;
  thisSpace->inner=NULL;
  thisSpace->direction=NULL;
  thisSpace->add=NULL;
  thisSpace->scale=NULL;
  thisSpace->writedata=NULL;
  thisSpace->freedata=NULL;

  return thisSpace;
 }

void MFNSpaceSetDistance(MFNSpace thisSpace,double (*distance)(MFNSpace,MFNVector,MFNVector,void*,MFErrorHandler), MFErrorHandler e)
 {
  static char RoutineName[]={"MFNSpaceSetDistance"};

#ifdef MFNOCONFIDENCE
  if(thisSpace==NULL)
   {
    sprintf(MFNSpaceErrorMsg,"NSpace, (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNSpaceErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  thisSpace->distance=distance;
  return;
 }

void MFNSpaceSetInnerProduct(MFNSpace thisSpace,double (*inner)(MFNSpace,MFNVector,MFNVector,void*,MFErrorHandler), MFErrorHandler e)
 {
  static char RoutineName[]={"MFNSpaceSetDistance"};

#ifdef MFNOCONFIDENCE
  if(thisSpace==NULL)
   {
    sprintf(MFNSpaceErrorMsg,"NSpace, (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNSpaceErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  thisSpace->inner=inner;
  return;
 }

void MFNSpaceSetDirection(MFNSpace thisSpace,void (*direction)(MFNSpace,MFNVector,MFNVector,MFNVector,void*,MFErrorHandler), MFErrorHandler e)
 {
  static char RoutineName[]={"MFNSpaceSetDistance"};

#ifdef MFNOCONFIDENCE
  if(thisSpace==NULL)
   {
    sprintf(MFNSpaceErrorMsg,"NSpace, (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNSpaceErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  thisSpace->direction=direction;
  return;
 }

void MFNSpaceSetAdd(MFNSpace thisSpace,void (*add)(MFNSpace,MFNVector,MFNVector,MFNVector,void*,MFErrorHandler), MFErrorHandler e)
 {
  static char RoutineName[]={"MFNSpaceSetDistance"};

#ifdef MFNOCONFIDENCE
  if(thisSpace==NULL)
   {
    sprintf(MFNSpaceErrorMsg,"NSpace, (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNSpaceErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  thisSpace->add=add;
  return;
 }

void MFNSpaceSetScale(MFNSpace thisSpace,void (*scale)(MFNSpace,double,MFNVector,MFNVector,void*,MFErrorHandler), MFErrorHandler e)
 {
  static char RoutineName[]={"MFNSpaceSetDistance"};

#ifdef MFNOCONFIDENCE
  if(thisSpace==NULL)
   {
    sprintf(MFNSpaceErrorMsg,"NSpace, (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNSpaceErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  thisSpace->scale=scale;
  return;
 }

void MFNSpaceSetWriteData(MFNSpace thisSpace,void (*writedata)(FILE*,MFNSpace,void*,MFErrorHandler), MFErrorHandler e)
 {
  static char RoutineName[]={"MFNSpaceSetDistance"};

#ifdef MFNOCONFIDENCE
  if(thisSpace==NULL)
   {
    sprintf(MFNSpaceErrorMsg,"NSpace, (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNSpaceErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  thisSpace->writedata=writedata;
  return;
 }

void MFNSpaceSetFreeData(MFNSpace thisSpace,void (*freedata)(void *,MFErrorHandler), MFErrorHandler e)
 {
  static char RoutineName[]={"MFNSpaceSetDistance"};

#ifdef MFNOCONFIDENCE
  if(thisSpace==NULL)
   {
    sprintf(MFNSpaceErrorMsg,"NSpace, (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNSpaceErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  thisSpace->freedata=freedata;
  return;
 }

void MFNSpaceSetData(MFNSpace thisSpace,void *data, MFErrorHandler e)
 {
  static char RoutineName[]={"MFNSpaceSetDistance"};

#ifdef MFNOCONFIDENCE
  if(thisSpace==NULL)
   {
    sprintf(MFNSpaceErrorMsg,"NSpace, (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNSpaceErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  thisSpace->d=data;
  return;
 }

void *MFNSpaceGetData(MFNSpace thisSpace, MFErrorHandler e)
 {
  static char RoutineName[]={"MFNSpaceGetData"};

#ifdef MFNOCONFIDENCE
  if(thisSpace==NULL)
   {
    sprintf(MFNSpaceErrorMsg,"NSpace, (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNSpaceErrorMsg,__LINE__,__FILE__);
    return NULL;
   }
#endif

  return thisSpace->d;
 }

char *MFNSpaceGetId(MFNSpace thisSpace, MFErrorHandler e)
 {
  static char RoutineName[]={"MFNSpaceGetId"};

#ifdef MFNOCONFIDENCE
  if(thisSpace==NULL)
   {
    sprintf(MFNSpaceErrorMsg,"NSpace, (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNSpaceErrorMsg,__LINE__,__FILE__);
    return NULL;
   }
#endif

  return thisSpace->id;
 }


#ifdef __cplusplus
}
#endif
