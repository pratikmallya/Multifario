/*
 *  PROGRAM NAME:  multifario
 *
 *  (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
 *  CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
 *
 *  Please refer to the LICENSE file in the top directory
 *
 *      author: Mike Henderson mhender@watson.ibm.com
 *      date:   July 15, 2002 Modified MFAUTO
 */

static char *id="@(#) $Id: MFAUTONSpace.c,v 1.10 2009/04/09 15:01:33 mhender Exp $";

#include <multifarioConfig.h>

char MFAUTONSpaceErrorMsg[256]="";

#include <MFNSpace.h>
#include <MFNVector.h>
#include <MFErrorHandler.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <MFAUTO.h>

#ifdef __cplusplus
 extern "C" {
#endif

struct MFAUTONSpaceData
 {
   integer nicp;
   integer *icp;
   doublereal *thu;
   doublereal *thl;
   iap_type *iap;
 };

static void MFFreeAUTONSpaceData(void*,MFErrorHandler);
static double MFAUTONSpaceDistance(MFNSpace,MFNVector,MFNVector,void*,MFErrorHandler);
static void MFAUTONSpaceDirection(MFNSpace,MFNVector,MFNVector,MFNVector,void*,MFErrorHandler);
static void MFAUTONSpaceAdd(MFNSpace,MFNVector,MFNVector,MFNVector,void*,MFErrorHandler);
static void MFAUTONSpaceScale(MFNSpace,double,MFNVector,MFNVector,void*,MFErrorHandler);
static double MFAUTONSpaceInner(MFNSpace,MFNVector,MFNVector,void*,MFErrorHandler);
void MFPrintAUTOBVNVectorFull(FILE*,MFNVector,MFErrorHandler);

doublereal *MFAUTONSpaceGetThl(MFNSpace,MFErrorHandler);
doublereal *MFAUTONSpaceGetThu(MFNSpace,MFErrorHandler);

MFNSpace MFCreateAUTONSpace(MFAUTOTPBVP tpbvp, doublereal *thu, doublereal *thl, MFErrorHandler e)
 {
  static char RoutineName[]={"MFCreateAUTONSpace"};
  MFNSpace this;
  struct MFAUTONSpaceData *data;
  int i;

  this=MFCreateNSpaceBaseClass("AUTONSpace",e);
  MFNSpaceSetDistance(this,MFAUTONSpaceDistance,e);
  MFNSpaceSetInnerProduct(this,MFAUTONSpaceInner,e);
  MFNSpaceSetDirection(this,MFAUTONSpaceDirection,e);
  MFNSpaceSetAdd(this,MFAUTONSpaceAdd,e);
  MFNSpaceSetScale(this,MFAUTONSpaceScale,e);
  MFNSpaceSetFreeData(this,MFFreeAUTONSpaceData,e);

  data=malloc(sizeof(struct MFAUTONSpaceData));

#ifndef MFNOSAFETYNET
  if(data==NULL)
   {
    sprintf(MFAUTONSpaceErrorMsg,"Out of memory trying to allocate %d bytes.",sizeof(struct MFAUTONSpaceData));
    MFSetError(e,12,RoutineName,MFAUTONSpaceErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  data->nicp=MFAUTOTPBVPGetNPAR(tpbvp,e);
  data->icp=malloc(2*NPARX*sizeof(integer));

#ifndef MFNOSAFETYNET
  if(data->icp==NULL)
   {
    sprintf(MFAUTONSpaceErrorMsg,"Out of memory trying to allocate %d bytes.",2*NPARX*sizeof(integer));
    MFSetError(e,12,RoutineName,MFAUTONSpaceErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  for(i=0;i<2*NPARX;i++)data->icp[i]=i;
  for(i=0;i<data->nicp;i++)data->icp[i]=(MFAUTOTPBVPGetICP(tpbvp,e))[i];

  data->thu=thu;
  data->thl=malloc(NPARX*sizeof(doublereal));

#ifndef MFNOSAFETYNET
  if(data->thl==NULL)
   {
    sprintf(MFAUTONSpaceErrorMsg,"Out of memory trying to allocate %d bytes.",NPARX*sizeof(doublereal));
    MFSetError(e,12,RoutineName,MFAUTONSpaceErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  for(i=0;i<NPARX;i++)data->thl[i]=0.;
  for(i=0;i<data->nicp;i++)data->thl[i]=thl[i];

  data->iap=malloc(sizeof(iap_type));

#ifndef MFNOSAFETYNET
  if(data->iap==NULL)
   {
    sprintf(MFAUTONSpaceErrorMsg,"Out of memory trying to allocate %d bytes.",sizeof(iap_type));
    MFSetError(e,12,RoutineName,MFAUTONSpaceErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  data->iap->ndim=MFAUTOTPBVPGetNDIM(tpbvp,e);
  data->iap->ntst=MFAUTOTPBVPGetNTST(tpbvp,e);
  data->iap->ncol=MFAUTOTPBVPGetNCOL(tpbvp,e);
  data->iap->nfpr=MFAUTOTPBVPGetNFPR(tpbvp,e);
  MFNSpaceSetData(this,(void*)data,e);

  return this;
 }

void MFFreeAUTONSpaceData(void *data, MFErrorHandler e)
 {
  static char RoutineName[]={"MFFreeNSpaceData"};
  struct MFAUTONSpaceData *d;

  d=(struct MFAUTONSpaceData*)data;
  if(d!=NULL)
   {
    if(d->thl!=NULL)free(d->thl);
    if(d->icp!=NULL)free(d->icp);
    if(d->iap!=NULL)free(d->iap);
    free(d);
   }

  return;
 }

double MFAUTONSpaceDistance(MFNSpace this,MFNVector v0,MFNVector v1,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFNSpaceDistance"};
  double ss;
  MFNVector diff;

  int i;
  diff=MFCloneNVector(v0,e);
  MFAUTONSpaceDirection(this,v0,v1,diff,d,e);
  ss=sqrt(MFAUTONSpaceInner(this,diff,diff,d,e));
  MFFreeNVector(diff,e);
  return ss;
 }

void MFAUTONSpaceDirection(MFNSpace this,MFNVector v0,MFNVector v1,MFNVector diff,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFNSpaceDirection"};
  int i,j;
  struct MFAUTONSpaceData *data;
  long ntst,ncol,ndim,npar,nfpr;
  long *icp;
  double **u0;
  double *par0;
  double **u1;
  double *par1;
  double **udir;
  double *pardir;

  data=(struct MFAUTONSpaceData*)d;

  ntst=MFAUTOBVNVGetNtst(v0,e);
  ncol=MFAUTOBVNVGetNcol(v0,e);
  ndim=MFAUTOBVNVGetNdim(v0,e);
  npar=MFAUTOBVNVGetNpar(v0,e);
  nfpr=data->iap->nfpr;
  icp=data->icp;

  u0=MFAUTOBVNVGetU(v0,e);
  par0=MFAUTOBVNVGetPar(v0,e);
  u1=MFAUTOBVNVGetU(v1,e);
  par1=MFAUTOBVNVGetPar(v1,e);
  udir=MFAUTOBVNVGetU(diff,e);
  pardir=MFAUTOBVNVGetPar(diff,e);

  for(i=0;i<ntst;i++)
    for(j=0;j<ndim*ncol;j++)
      udir[i][j]=u1[i][j]-u0[i][j];
  for(j=0;j<ndim;j++)
     udir[ntst][j]=u1[ntst][j]-u0[ntst][j];

  for(i=0;i<npar;i++)
    pardir[i]=par1[i]-par0[i];

  for(i=0;i<nfpr;i++)
    pardir[icp[i]]=par1[icp[i]]-par0[icp[i]];

  return;
 }

void MFAUTONSpaceAdd(MFNSpace this,MFNVector v0,MFNVector v1,MFNVector sum,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFNAUTOSpaceAdd"};
  struct MFAUTONSpaceData *data;
  int i,j;
  long ntst,ncol,ndim,nfpr,*icp;
  double **u0;
  double *par0;
  double **u1;
  double *par1;
  double **usum;
  double *parsum;

  /* This function now calls the AUTO routines to compute the
     distance between two vectors.  I got this code from
     scaleb.*/

  data=(struct MFAUTONSpaceData*)d;

  ntst=MFAUTOBVNVGetNtst(v0,e);
  ncol=MFAUTOBVNVGetNcol(v0,e);
  ndim=MFAUTOBVNVGetNdim(v0,e);
  nfpr=data->iap->nfpr;
  icp=MFAUTOBVNVGetICP(v0,e);

  u0=MFAUTOBVNVGetU(v0,e);
  par0=MFAUTOBVNVGetPar(v0,e);
  u1=MFAUTOBVNVGetU(v1,e);
  par1=MFAUTOBVNVGetPar(v1,e);
  usum=MFAUTOBVNVGetU(sum,e);
  parsum=MFAUTOBVNVGetPar(sum,e);

  for(i=0;i<ntst+1;i++)
   {
    for(j=0;j<ndim*ncol;j++)
     {
      usum[i][j]=u0[i][j]+u1[i][j];
     }
   }
  for(i=0;i<nfpr;i++)
    parsum[icp[i]]=par0[icp[i]]+par1[icp[i]];

  return;
 }

double MFAUTONSpaceInner(MFNSpace this,MFNVector v0,MFNVector v1,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTONSpaceInner"};
  int i;
  int it,icol;
  long ntst,ncol,ndim,npar,nfpr;
  struct MFAUTONSpaceData *data;
  double **u0;
  double *par0;
  double **u1;
  double *par1;
  double *dtm;
  double ss;
  long *icp;
  double t;
  int idim;

  /* This function now calls the AUTO routines to compute the
     product of two vectors.  I got this code from scaleb.*/

  data=(struct MFAUTONSpaceData*)d;

  ntst=MFAUTOBVNVGetNtst(v0,e);
  ncol=MFAUTOBVNVGetNcol(v0,e);
  ndim=MFAUTOBVNVGetNdim(v0,e);
  npar=MFAUTOBVNVGetNpar(v0,e);
  nfpr=MFAUTOBVNVGetNfpr(v0,e);
  icp=MFAUTOBVNVGetICP(v0,e);

  u0=MFAUTOBVNVGetU(v0,e);
  par0=MFAUTOBVNVGetPar(v0,e);
  u1=MFAUTOBVNVGetU(v1,e);
  par1=MFAUTOBVNVGetPar(v1,e);

  dtm=MFAUTOBVNVGetDt(v0,e);
  
  ss=0.;
  for(i=0;i<nfpr;i++)
    ss+=data->thl[i]*par0[icp[i]]*par1[icp[i]];

  ss += rinpr(data->iap, ndim, u0, u1, dtm, data->thu);

  return ss;
 }

void MFAUTONSpaceScale(MFNSpace this,double s,MFNVector v,MFNVector prod,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFNAUTOSpaceScale"};
  int i,j;
  struct MFAUTONSpaceData *data;
  double **vu,**pu;
  double  *vr, *pr;
  int ntst,ncol,ndim,nfpr;
  long *icp;
  int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("%s\n",RoutineName);
    printf("s: %lf\n",s);
    printf("v:\n");
    MFPrintNVector(stdout,v);fflush(stdout);
   }
#endif

  ntst=MFAUTOBVNVGetNtst(v,e);
  ncol=MFAUTOBVNVGetNcol(v,e);
  ndim=MFAUTOBVNVGetNdim(v,e);
  nfpr=MFAUTOBVNVGetNfpr(v,e);
  icp=MFAUTOBVNVGetICP(v,e);

  vu=MFAUTOBVNVGetU(v,e);
  vr=MFAUTOBVNVGetPar(v,e);
  pu=MFAUTOBVNVGetU(prod,e);
  pr=MFAUTOBVNVGetPar(prod,e);

  for(i=0;i<ntst+1;i++)for(j=0;j<ncol*ndim;j++)(pu[i])[j]=s*(vu[i])[j];

  for(i=0;i<NPARX;i++)pr[i]=s*vr[i];

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("\ns*v:\n");
    MFPrintNVector(stdout,prod);
    printf("\n");
   }
#endif

  return;
 }

doublereal *MFAUTONSpaceGetThl(MFNSpace this, MFErrorHandler e)
 {
  static char RoutineName[]={"MFNAUTONSpaceGetThl"};

  struct MFAUTONSpaceData *data;
  if(strcmp(MFNSpaceGetId(this,e),"AUTONSpace"))
   {
    sprintf(MFAUTONSpaceErrorMsg,"NSpace, (argument 1) is not an AUTO NSpace");
    MFSetError(e,12,RoutineName,MFAUTONSpaceErrorMsg,__LINE__,__FILE__);
    return NULL;
   }

  data=(struct MFAUTONSpaceData*)MFNSpaceGetData(this,e);

  return data->thl;
 }

doublereal *MFAUTONSpaceGetThu(MFNSpace this, MFErrorHandler e)
 {
  static char RoutineName[]={"MFNAUTONSpaceGetThu"};

  struct MFAUTONSpaceData *data;
  if(strcmp(MFNSpaceGetId(this,e),"AUTONSpace"))
   {
    sprintf(MFAUTONSpaceErrorMsg,"NSpace, (argument 1) is not an AUTO NSpace");
    MFSetError(e,12,RoutineName,MFAUTONSpaceErrorMsg,__LINE__,__FILE__);
    return NULL;
   }

  data=(struct MFAUTONSpaceData*)MFNSpaceGetData(this,e);

  return data->thu;
 }

#ifdef __cplusplus
}
#endif
