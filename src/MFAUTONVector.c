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

static char *id="@(#) $Id: MFAUTONVector.c,v 1.10 2009/04/09 15:01:33 mhender Exp $";

#include <multifarioConfig.h>

static char MFNVectorErrorMsg[256]="";

#include <MFErrorHandler.h>
#include <MFNVector.h>
#include <MFAUTO.h>
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

struct MFAUTOBVNVectorData
 {
  int nC;

  int ntst;
  int ndim;
  int ncol;
  int npar;
  int nfpr;
  int nicp;
  long *icp;
  doublereal **u;
  doublereal *par;
  doublereal *t;
  doublereal *dt;
  doublereal *thu;
  char type[4];

  double R;
  int nit;

/* Bifurcation detection values */

  int nuz;
  doublereal *uzbv; /* Computed in fnuzbv */
  doublereal lpbv;  /* Computed in fnlpbv */
  doublereal bpbv;  /* Computed in fnbpbv */
  doublereal spbv;  /* Computed in fnspbv */

  doublecomplex *ev;  /* Computed in fnspbv (stability of periodic orbits) */
  doublereal **p0;  /* Computed in solvbv */
  doublereal **p1;  /* Computed in solvbv */
 };

static void MFFreeAUTOBVNVectorData(void*,MFErrorHandler);
static int MFAUTOBVNVGetNC(void*,MFErrorHandler);
static doublereal MFAUTOBVNVGetC(int,void*,MFErrorHandler);
static void MFAUTOBVNVSetC(int,double,void*,MFErrorHandler);
static void MFAUTOBVNVDiff(void*,void*,void*,MFErrorHandler);
static void MFAUTOBVNVAdd(void*,void*,void*,MFErrorHandler);
static void MFWriteAUTOBVNVectorData(FILE*,void*,MFErrorHandler);
MFNVector MFReadAUTOBVNVector(FILE*,MFErrorHandler);
static MFNVector MFCloneAUTOBVNVector(void*,MFErrorHandler);

MFNVector MFCreateAUTOBVNVector(int,int,int,int,int,int,int,long*,MFErrorHandler);
void MFAUTOBVNVCopyDataValues(MFNVector,MFNVector,MFErrorHandler);

doublereal **MFAUTOBVNVGetU(MFNVector,MFErrorHandler);
doublereal *MFAUTOBVNVGetPar(MFNVector,MFErrorHandler);
int MFAUTOBVNVGetNtst(MFNVector,MFErrorHandler);
int MFAUTOBVNVGetNcol(MFNVector,MFErrorHandler);
int MFAUTOBVNVGetNdim(MFNVector,MFErrorHandler);
int MFAUTOBVNVGetNpar(MFNVector,MFErrorHandler);
int MFAUTOBVNVGetNfpr(MFNVector,MFErrorHandler);
long *MFAUTOBVNVGetICP(MFNVector,MFErrorHandler);

void MFAUTOBVNVSetT(MFNVector,doublereal*,MFErrorHandler);

void MFAUTOBVNVSetNUZ(MFNVector,int,MFErrorHandler);
void MFAUTOBVNVSetUZBV(MFNVector,int,doublereal,MFErrorHandler);
void MFAUTOBVNVSetLPBV(MFNVector,doublereal,MFErrorHandler);
void MFAUTOBVNVSetBPBV(MFNVector,doublereal,MFErrorHandler);
void MFAUTOBVNVSetSPBV(MFNVector,doublereal,MFErrorHandler);
void MFAUTOBVNVSetEV(MFNVector,doublecomplex*,MFErrorHandler);
void MFAUTOBVNVSetP0(MFNVector,doublereal**,MFErrorHandler);
void MFAUTOBVNVSetP1(MFNVector,doublereal**,MFErrorHandler);

int MFAUTOBVNVGetNUZ(MFNVector,MFErrorHandler);
doublereal MFAUTOBVNVGetUZBV(MFNVector,int,MFErrorHandler);
doublereal MFAUTOBVNVGetLPBV(MFNVector,MFErrorHandler);
doublereal MFAUTOBVNVGetBPBV(MFNVector,MFErrorHandler);
doublereal MFAUTOBVNVGetSPBV(MFNVector,MFErrorHandler);
doublecomplex *MFAUTOBVNVGetEV(MFNVector,MFErrorHandler);
doublereal** MFAUTOBVNVGetP0(MFNVector,MFErrorHandler);
doublereal** MFAUTOBVNVGetP1(MFNVector,MFErrorHandler);

void MFPrintAUTOBVNVectorFull(FILE*,MFNVector,MFErrorHandler);
void MFPrintAUTOBVNVectorDataFull(FILE*,void*,MFErrorHandler);
void MFPrintAUTOBVNVector(FILE*,void*,MFErrorHandler);

void MFFreeAUTOBVNVectorData(void *data, MFErrorHandler e)
 {
  static char RoutineName[]={"MFFreeAUTOBVNVectorData"};
  int i;
  struct MFAUTOBVNVectorData *this;
  int verbose=0;

  this=(struct MFAUTOBVNVectorData*)data;

  if(this->u!=NULL)
   {
    free_dmatrix(this->u);
    this->u=NULL;
   }
  if(this->par!=NULL)
   {
    free(this->par);
    this->par=NULL;
   }
  if(this->t!=NULL)
   {
   free(this->t);
    this->t=NULL;
   }
  if(this->dt!=NULL)
   {
    free(this->dt);
    this->dt=NULL;
   }
  if(this->ev!=NULL)
   {
    free(this->ev);
    this->ev=NULL;
   }
  if(this->p0!=NULL)
   {
    free_dmatrix(this->p0);
    this->p0=NULL;
   }
  if(this->p1!=NULL)
   {
    free_dmatrix(this->p1);
    this->p1=NULL;
   }

  if(this->uzbv!=NULL)free(this->uzbv);

  free(this);

  return;
 }

int MFAUTOBVNVGetNC(void *data, MFErrorHandler e)
 {
  static char RoutineName[]={"MFNV_NC"};
  struct MFAUTOBVNVectorData *this;

  this=(struct MFAUTOBVNVectorData*)data;

#ifdef MFNOCONFIDENCE
  if(this==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector Data (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return -1;
   }
#endif

  return this->nC;
 }

double MFAUTOBVNVGetC(int i,void *data, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOBVNVGetC"};
  struct MFAUTOBVNVectorData *this;
  int j;
  double result;

  this=(struct MFAUTOBVNVectorData*)data;

#ifdef MFNOCONFIDENCE
  if(this==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector Data (argument 2) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return DBL_QNAN;
   }

  if(i<0|| !(i<this->nC))
   {
    sprintf(MFNVectorErrorMsg,"Coordinate %d (argument 1) is illegal. Must be in 0 to %d",i,this->nC-1);
    MFSetError(e,8,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return DBL_QNAN;
   }
#endif

  if(i<(this->ntst*this->ncol+1)*this->ndim)
   {
    j=i%(this->ntst+1);
    i=i/(this->ntst+1);
    result=(this->u[i])[j];
   }else{
    i=i-(this->ntst*this->ncol+1)*this->ndim;
    result=this->par[this->icp[i]];
   }

  return result;
 }

void MFAUTOBVNVSetC(int i,doublereal vl,void *data, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOBVNVSetC"};
  int j;
  struct MFAUTOBVNVectorData *this;

  this=(struct MFAUTOBVNVectorData*)data;

#ifdef MFNOCONFIDENCE
  if(this==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector Data (argument 3) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(i<0|| !(i<this->nC))
   {
    sprintf(MFNVectorErrorMsg,"Coordinate %d (argument 1) is illegal. Must be in 0 to %d",i,this->nC-1);
    MFSetError(e,8,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  if(i<(this->ntst*this->ncol+1)*this->ndim)
   {
    j=i%(this->ntst+1);
    i=i/(this->ntst+1);
    (this->u[i])[j]=vl;
   }else{
    i=i-(this->ntst*this->ncol+1)*this->ndim;
    this->par[this->icp[i]]=vl;
   }

  return;
 }

void MFAUTOBVNVDiff(void *adata,void *bdata, void *cdata, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOBVNVDiff"};
  int i,ni;
  int j,nj;
  struct MFAUTOBVNVectorData *a;
  struct MFAUTOBVNVectorData *b;
  struct MFAUTOBVNVectorData *c;

  a=(struct MFAUTOBVNVectorData*)adata;

#ifndef MFNOSAFETYNET
  if(a==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector Data for a (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  b=(struct MFAUTOBVNVectorData*)bdata;

#ifndef MFNOSAFETYNET
  if(b==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector Data for b (argument 2) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  c=(struct MFAUTOBVNVectorData*)cdata;

#ifndef MFNOSAFETYNET
  if(c==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector Data for c (argument 3) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

#ifdef MFNOCONFIDENCE
  if(a->nC!=b->nC || a->nC!=c->nC || b->nC!=c->nC)
   {
    sprintf(MFNVectorErrorMsg,"Vectors must all be the same length a=%d, b=%d, c=%d",a->nC,b->nC,c->nC);
    MFSetError(e,4,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  ni=a->ntst+1;
  nj=a->ncol*a->ndim;

  for(i=0;i<ni;i++)
   for(j=0;j<nj;j++)
    {
     (c->u[i])[j]=(a->u[i])[j]-(b->u[i])[j];
    }
  for(i=0;i<a->npar;i++)
   {
    c->par[i]=a->par[i]-b->par[i];
   }
  if(a->t!=NULL)
   {
    for(i=0;i<a->ntst+1;i++)
       c->t[i]=a->t[i];
    for(i=0;i<a->ntst;i++)
       c->dt[i]=a->dt[i];
   }else if(b->t!=NULL)
   {
    for(i=0;i<b->ntst+1;i++)
       c->t[i]=b->t[i];
    for(i=0;i<b->ntst;i++)
       c->dt[i]=b->dt[i];
   }

  return;
 }

void MFAUTOBVNVAdd(void *adata,void *bdata,void *cdata, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOBVNVAdd"};
  int i,ni;
  int j,nj;
  struct MFAUTOBVNVectorData *a;
  struct MFAUTOBVNVectorData *b;
  struct MFAUTOBVNVectorData *c;

  a=(struct MFAUTOBVNVectorData*)adata;

#ifndef MFSAFETYNET
  if(a==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector Data for a (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  b=(struct MFAUTOBVNVectorData*)bdata;

#ifdef MFSAFETYNET
  if(b==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector Data for b (argument 2) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  c=(struct MFAUTOBVNVectorData*)cdata;

#ifdef MFSAFETYNET
  if(c==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector Data for c (argument 3) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

#ifdef MFNOCONFIDENCE
  if(a->nC!=b->nC || a->nC!=c->nC || b->nC!=c->nC)
   {
    sprintf(MFNVectorErrorMsg,"Vectors must all be the same length a=%d, b=%d, c=%d",a->nC,b->nC,c->nC);
    MFSetError(e,4,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  ni=a->ntst+1;
  nj=a->ncol*a->ndim;

  for(i=0;i<ni;i++)
   for(j=0;j<nj;j++)
     (c->u[i])[j]=(a->u[i])[j]+(b->u[i])[j];
  for(i=0;i<a->npar;i++)
     c->par[i]=a->par[i]+b->par[i];
  if(a->t!=NULL)
   {
    for(i=0;i<a->ntst+1;i++)
       c->t[i]=a->t[i];
    for(i=0;i<a->ntst;i++)
       c->dt[i]=a->dt[i];
   }

  return;
 }

void MFWriteAUTOBVNVectorData(FILE *fid,void *data, MFErrorHandler e)
 {
  static char RoutineName[]={"MFWriteAUTOBVNVectorData"};
  int i,j;
  struct MFAUTOBVNVectorData *u;

#ifdef MFNOCONFIDENCE
  if(fid==NULL)
   {
    sprintf(MFNVectorErrorMsg,"fid (argument 1) is NULL.");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(data==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector Data (argument 2) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  u=(struct MFAUTOBVNVectorData*)data;

  fprintf(fid,"%d %d %d %d %d %d\n",u->ntst,u->ncol,u->ndim,u->npar,u->nuz,u->nfpr);
  for(i=0;i<u->nfpr-1;i++)fprintf(fid,"%d ",u->icp[i]);
  fprintf(fid,"%d\n",u->icp[u->nfpr-1]);
  for(i=0;i<u->ntst+1;i++)
   {
    for(j=0;j<u->ncol*u->ndim;j++)
     {
      if(i>0||j>0)fprintf(fid," ");
      fprintf(fid,"%lf",(u->u[i])[j]);
     }
   }
  fprintf(fid,"\n");

  for(i=0;i<u->npar;i++)
   {
    if(i>0)fprintf(fid," ");
    fprintf(fid,"%lf",u->par[i]);
   }
  fprintf(fid,"\n");

  if(u->t!=NULL)
   {
    fprintf(fid,"%d\n",1);
    for(i=0;i<u->ntst+1;i++)
     {
      if(i>0)fprintf(fid," ");
      fprintf(fid,"%lf",u->t[i]);
     }
    fprintf(fid,"\n");
   }else{
    fprintf(fid,"%d\n",0);
   }

  if(u->ev!=NULL)
   {
    fprintf(fid,"%d\n",1);
    for(i=0;i<u->ndim;i++)
     {
      if(i>0)fprintf(fid," ");
      fprintf(fid,"%lf",u->ev[i]);
     }
    fprintf(fid,"\n");
   }else{
    fprintf(fid,"%d %d\n",0,0);
   }

  return;
 }

MFNVector MFReadAUTOBVNVector(FILE *fid, MFErrorHandler e)
 {
  static char RoutineName[]={"MFReadAUTOBVNVector"};
  MFNVector this;
  int ntst,ncol,ndim,npar,nuz,nfpr,nicp;
  int i,j,there;
  struct MFAUTOBVNVectorData *data;
  long *icp;

  fscanf(fid,"%d %d %d %d %d %d\n",&ntst,&ncol,&ndim,&npar,&nuz,&nfpr,&nicp);
  icp=malloc(nfpr*sizeof(long));

#ifndef MFNOSAFETYNET
  if(icp==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Out of memory trying to allocate %d bytes.",nfpr*sizeof(long));
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  for(i=0;i<nicp-1;i++)fscanf(fid,"%d",icp+i);
  fscanf(fid,"%d\n",icp+nicp-1);

  this=MFCreateAUTOBVNVector(ntst,ncol,ndim,npar,nuz,nfpr,nicp,icp,e);
  data=(struct MFAUTOBVNVectorData*)MFNVectorGetData(this,e);

  for(i=0;i<ntst+1;i++)
   {
    for(j=0;j<ncol*ndim;j++)
     {
      if(i>0)fscanf(fid," ");
      fscanf(fid,"%lf",&((data->u[i])[j]));
     }
   }
  fscanf(fid,"\n");

  for(i=0;i<data->npar;i++)
   {
    if(i>0)fscanf(fid," ");
    fscanf(fid,"%lf",&(data->par[i]));
   }
  fscanf(fid,"\n");

  fscanf(fid,"%d\n",&there);
  if(there)
   {
    data->t=malloc((data->ntst+1)*sizeof(double));

#ifndef MFNOSAFETYNET
    if(data->t==NULL)
     {
      sprintf(MFNVectorErrorMsg,"Out of memory trying to allocate %d bytes.",(data->ntst+1)*sizeof(double));
      MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

    for(i=0;i<data->ntst+1;i++)
     {
      if(i>0)fscanf(fid," ");
      fscanf(fid,"%lf",&(data->t[i]));
     }
    fscanf(fid,"\n");
    data->dt=malloc((data->ntst+1)*sizeof(double));

#ifndef MFNOSAFETYNET
    if(data->dt==NULL)
     {
      sprintf(MFNVectorErrorMsg,"Out of memory trying to allocate %d bytes.",(data->ntst+1)*sizeof(double));
      MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

    for(i=0;i<data->ntst;i++)data->dt[i]=data->t[i+1]-data->t[i];
   }

  fscanf(fid,"%d\n",&there);
  if(there)
   {
    data->ev=malloc(data->ndim*sizeof(doublecomplex));

#ifndef MFNOSAFETYNET
    if(data->ev==NULL)
     {
      sprintf(MFNVectorErrorMsg,"Out of memory trying to allocate %d bytes.",data->ndim*sizeof(doublecomplex));
      MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

    for(i=0;i<data->ndim;i++)
     {
      if(i>0)fscanf(fid," ");
      fscanf(fid,"%lf",&(data->ev[i]));
     }
    fscanf(fid,"\n");
   }

  MFNVectorSetData(this,data,e);

  return this;
 }

MFNVector MFCloneAUTOBVNVector(void *data, MFErrorHandler e)
 {
  static char RoutineName[]={"MFCloneAUTOBVNVector"};
  struct MFAUTOBVNVectorData *u;
  struct MFAUTOBVNVectorData *v;
  MFNVector clone;
  int i,j;

  u=(struct MFAUTOBVNVectorData*)data;

  clone=MFCreateAUTOBVNVector(u->ntst,u->ncol,u->ndim,u->npar,u->nuz,u->nfpr,u->nicp,u->icp,e);

  v=(struct MFAUTOBVNVectorData*)MFNVectorGetData(clone,e);

  for(i=0;i<u->ntst+1;i++)
   for(j=0;j<u->ncol*u->ndim;j++)
     (v->u[i])[j]=(u->u[i])[j];
  for(i=0;i<u->nfpr;i++)v->par[i]=u->par[i];
  for(i=0;i<u->ntst+1;i++)v->t[i]=u->t[i];
  for(i=0;i<u->ntst;i++)v->dt[i]=u->dt[i];

  v->thu=u->thu;

  v->R=u->R;
  v->nit=u->nit;

  return clone;
 }

void MFPrintAUTOBVNVector(FILE *fid,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFPrintAUTOBVNVector"};
  double unorm;
  int i;
  int it,icol,idim;
  struct MFAUTOBVNVectorData *data;
  iap_type iap;
  int j;
  int branch,point,label;
  char type[3];
  double umaxnorm;
  double ul2norm;
  int n;
  int ic;

  if(0)
   {
    MFPrintAUTOBVNVectorDataFull(fid,d,e);
    return;
   }

  data=(struct MFAUTOBVNVectorData*)d;

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

  iap.ndim = data->ndim;
  iap.ntst = data->ntst;
  iap.ncol = data->ncol;

  branch=1;
  point=1;
  label=1;
  fprintf(fid," %d %d %s %d",branch,point,data->type,label);fflush(stdout);
  for(i=0;i<data->npar;i++)fprintf(fid," %le",data->par[data->icp[i]]);

  ul2norm=sqrt(rinpr(&iap, data->ndim, data->u, data->u, data->dt, data->thu));
  fprintf(fid," %le",ul2norm);

  for(i=0;i<data->ndim;i++)
   {
    umaxnorm=0.;
    for(j=0;j<data->ntst+1;j++)if(fabs((data->u)[j][i])>umaxnorm)umaxnorm=fabs((data->u)[j][i]);
    fprintf(fid," %le",umaxnorm);
   }

  return;
 }

MFNVector MFCreateAUTOBVNVector(int ntst,int ncol, int ndim, int npar, int nuz, int nfpr, int nicp, long *icp, MFErrorHandler e)
 {
  static char RoutineName[]={"MFCreateAUTOBVNVector"};
  MFNVector this;
  int i,j;
  struct MFAUTOBVNVectorData *data;

  int verbose=0;

  data=malloc(sizeof(struct MFAUTOBVNVectorData));

#ifndef MFNOSAFETYNET
  if(data==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct MFAUTOBVNVectorData));
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  data->nC=(ntst*ncol+1)*ndim+npar;
  data->ntst=ntst;
  data->ndim=ndim;
  data->ncol=ncol;
  data->npar=npar;
  data->nfpr=nfpr;
  data->nicp=nicp;
  data->icp=icp;
  data->nuz=nuz;
  strcpy(data->type,"  ");

  data->u=dmatrix(data->ntst+1,data->ncol*data->ndim);

#ifndef MFNOSAFETYNET
  if(data->u==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Out of memory, trying to allocate %d bytes",(data->ntst+1)*data->ncol*data->ndim*sizeof(double));
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  data->par=malloc(NPARX*sizeof(double));

#ifndef MFNOSAFETYNET
    if(data->par==NULL)
     {
      sprintf(MFNVectorErrorMsg,"Out of memory, trying to allocate %d bytes",NPARX*sizeof(double));
      MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

  data->R=1.;
  data->nit=-1;

  data->t=malloc((ntst+1)*sizeof(double));

#ifndef MFNOSAFETYNET
  if(data->t==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Out of memory, trying to allocate %d bytes",(ntst+1)*sizeof(double));
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  data->dt=malloc((ntst+1)*sizeof(double));

#ifndef MFNOSAFETYNET
  if(data->dt==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Out of memory, trying to allocate %d bytes",(ntst+1)*sizeof(double));
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  data->thu=NULL;

  for(i=0;i<data->ntst+1;i++)
   for(j=0;j<data->ncol*data->ndim;j++)
     (data->u[i])[j]=0.;
  for(i=0;i<data->npar;i++)data->par[i]=0.;
  for(i=0;i<data->ntst+1;i++)data->t[i]=1.*i/data->ntst;
  for(i=0;i<data->ntst  ;i++)data->dt[i]=data->t[i+1]-data->t[i];

  data->ev=malloc(ndim*sizeof(doublecomplex));

#ifndef MFNOSAFETYNET
  if(data->ev==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Out of memory, trying to allocate %d bytes",ndim*sizeof(doublecomplex));
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  for(i=0;i<data->ndim;i++){data->ev[i].r=0.;data->ev[i].i=0.;}
  data->p0=dmatrix(ndim,ndim);

#ifndef MFNOSAFETYNET
  if(data->p0==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Out of memory, trying to allocate %d bytes",ndim*ndim*sizeof(double));
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  for(i=0;i<ndim;i++)for(j=0;j<ndim;j++)data->p0[i][j]=0.;

  data->p1=dmatrix(ndim,ndim);

#ifndef MFNOSAFETYNET
  if(data->p1==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Out of memory, trying to allocate %d bytes",ndim*ndim*sizeof(double));
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  for(i=0;i<ndim;i++)for(j=0;j<ndim;j++)data->p1[i][j]=0.;

  data->uzbv=malloc(data->nuz*sizeof(double));

#ifndef MFNOSAFETYNET
  if(data->uzbv==NULL && data->nuz>0)
   {
    sprintf(MFNVectorErrorMsg,"Out of memory, trying to allocate %d bytes",data->nuz*sizeof(double));
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  data->lpbv=1.;
  data->bpbv=1.;
  data->spbv=1.;
  for(i=0;i<nuz;i++)data->uzbv[i]=1.;

  this=MFCreateNVectorBaseClass("AUTO",e);
  MFNVectorSetData(this,data,e);
  MFNVectorSetFreeData(this,MFFreeAUTOBVNVectorData,e);

  MFNVectorSetWriteData(this,MFWriteAUTOBVNVectorData,e);
  MFNVectorSetGetNC(this,MFAUTOBVNVGetNC,e);
  MFNVectorSetGetC(this,MFAUTOBVNVGetC,e);
  MFNVectorSetSetC(this,MFAUTOBVNVSetC,e);
  MFNVectorSetDiff(this,MFAUTOBVNVDiff,e);
  MFNVectorSetAdd(this,MFAUTOBVNVAdd,e);
  MFNVectorSetClone(this,MFCloneAUTOBVNVector,e);
  MFNVectorSetPrint(this,MFPrintAUTOBVNVector,e);

  return this;
 }

doublereal **MFAUTOBVNVGetU(MFNVector this, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOBVNVGetU"};
  struct MFAUTOBVNVectorData *data;

#ifdef MFNOCONFIDENCE
  if(this==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return NULL;
   }

  if(strcmp(MFNVGetId(this,e),"AUTO"))
   {
    sprintf(MFNVectorErrorMsg,"Trying to get arrays from non-AUTO Vector type \"%s\"",MFNVGetId(this,e));
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return NULL;
   }
#endif

  data=(struct MFAUTOBVNVectorData*)MFNVectorGetData(this,e);

  return(data->u);
 }

doublereal *MFAUTOBVNVGetPar(MFNVector this, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOBVNVGetPar"};
  struct MFAUTOBVNVectorData *data;

#ifdef MFNOCONFIDENCE
  if(this==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return NULL;
   }

  if(strcmp(MFNVGetId(this,e),"AUTO"))
   {
    sprintf(MFNVectorErrorMsg,"Trying to get arrays from non-AUTO Vector type \"%s\"",MFNVGetId(this,e));
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return NULL;
   }
#endif

  data=(struct MFAUTOBVNVectorData*)MFNVectorGetData(this,e);

  return(data->par);
 }

doublereal *MFAUTOBVNVGetT(MFNVector this, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOBVNVGetT"};
  struct MFAUTOBVNVectorData *data;

#ifdef MFNOCONFIDENCE
  if(this==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return NULL;
   }

  if(strcmp(MFNVGetId(this,e),"AUTO"))
   {
    sprintf(MFNVectorErrorMsg,"Trying to get arrays from non-AUTO Vector type \"%s\"",MFNVGetId(this,e));
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return NULL;
   }
#endif

  data=(struct MFAUTOBVNVectorData*)MFNVectorGetData(this,e);

  return(data->t);
 }

void MFAUTOBVNVSetT(MFNVector this,doublereal *t, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOBVNVSetT"};
  struct MFAUTOBVNVectorData *data;
  int i;

#ifdef MFNOCONFIDENCE
  if(this==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(strcmp(MFNVGetId(this,e),"AUTO"))
   {
    sprintf(MFNVectorErrorMsg,"Trying to get arrays from non-AUTO Vector type \"%s\"",MFNVGetId(this,e));
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  data=(struct MFAUTOBVNVectorData*)MFNVectorGetData(this,e);

  for(i=0;i<data->ntst+1;i++)data->t[i]=t[i];
  for(i=0;i<data->ntst;i++)data->dt[i]=t[i+1]-t[i];
  return;
 }


doublereal *MFAUTOBVNVGetDt(MFNVector this, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOBVNVGetDt"};
  struct MFAUTOBVNVectorData *data;

#ifdef MFNOCONFIDENCE
  if(this==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return NULL;
   }

  if(strcmp(MFNVGetId(this,e),"AUTO"))
   {
    sprintf(MFNVectorErrorMsg,"Trying to get arrays from non-AUTO Vector type \"%s\"",MFNVGetId(this,e));
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return NULL;
   }
#endif

  data=(struct MFAUTOBVNVectorData*)MFNVectorGetData(this,e);

  return(data->dt);
 }

doublecomplex *MFAUTOBVNVGetEV(MFNVector this, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOBVNVGetEV"};
  struct MFAUTOBVNVectorData *data;

#ifdef MFNOCONFIDENCE
  if(this==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return NULL;
   }

  if(strcmp(MFNVGetId(this,e),"AUTO"))
   {
    sprintf(MFNVectorErrorMsg,"Trying to get arrays from non-AUTO Vector type \"%s\"",MFNVGetId(this,e));
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return NULL;
   }
#endif

  data=(struct MFAUTOBVNVectorData*)MFNVectorGetData(this,e);

  return(data->ev);
 }

void MFAUTOBVNVSetEV(MFNVector this, doublecomplex *ev, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOBVNVSetEV"};
  struct MFAUTOBVNVectorData *data;
  int i;


#ifdef MFNOCONFIDENCE
  if(this==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(strcmp(MFNVGetId(this,e),"AUTO"))
   {
    sprintf(MFNVectorErrorMsg,"Trying to get arrays from non-AUTO Vector type \"%s\"",MFNVGetId(this,e));
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  data=(struct MFAUTOBVNVectorData*)MFNVectorGetData(this,e);

  if(data->ev!=NULL)free(data->ev);
  data->ev=ev;

  return;
 }

void MFAUTOBVNVSetNUZ(MFNVector this,int nuz, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOBVNVSetNUZ"};
  struct MFAUTOBVNVectorData *data;
  int i;

#ifdef MFNOCONFIDENCE
  if(this==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(strcmp(MFNVGetId(this,e),"AUTO"))
   {
    sprintf(MFNVectorErrorMsg,"Trying to set number of user zeroes in a non-AUTO Vector type \"%s\"",MFNVGetId(this,e));
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  data=(struct MFAUTOBVNVectorData*)MFNVectorGetData(this,e);

  data->nuz=nuz;

  return;
 }

void MFAUTOBVNVSetUZBV(MFNVector this,int i,doublereal uzbv, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOBVNVSetUZBV"};
  struct MFAUTOBVNVectorData *data;

#ifdef MFNOCONFIDENCE
  if(this==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(strcmp(MFNVGetId(this,e),"AUTO"))
   {
    sprintf(MFNVectorErrorMsg,"Trying to get a user zero value in a non-AUTO Vector type \"%s\"",MFNVGetId(this,e));
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  data=(struct MFAUTOBVNVectorData*)MFNVectorGetData(this,e);

  if(i<0)
   {
    sprintf(MFNVectorErrorMsg,"Trying to set an invalid user zero value %d must be in the range [0,%d)",i,data->nuz);
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }
  if(data->uzbv==NULL||i>=data->nuz)
   {
    data->nuz=i+1;
    data->uzbv=realloc(data->uzbv,data->nuz*sizeof(double));

#ifndef MFNOSAFETYNET
    if(data->uzbv==NULL && data->nuz>0)
     {
      sprintf(MFNVectorErrorMsg,"Out of memory, trying to allocate %d bytes",data->nuz*sizeof(double));
      MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

   }

  data->uzbv[i]=uzbv;

  return;
 }

void MFAUTOBVNVSetLPBV(MFNVector this,doublereal lpbv, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOBVNVSetLPBV"};
  struct MFAUTOBVNVectorData *data;
  int i;

#ifdef MFNOCONFIDENCE
  if(this==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(strcmp(MFNVGetId(this,e),"AUTO"))
   {
    sprintf(MFNVectorErrorMsg,"Trying to set the Limit Point value in a non-AUTO Vector type \"%s\"",MFNVGetId(this,e));
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  data=(struct MFAUTOBVNVectorData*)MFNVectorGetData(this,e);

  data->lpbv=lpbv;

  return;
 }

void MFAUTOBVNVSetBPBV(MFNVector this,doublereal bpbv, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOBVNVSetBPBV"};
  struct MFAUTOBVNVectorData *data;
  int i;

#ifdef MFNOCONFIDENCE
  if(this==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(strcmp(MFNVGetId(this,e),"AUTO"))
   {
    sprintf(MFNVectorErrorMsg,"Trying to set the Bifurcation Point value in a non-AUTO Vector type \"%s\"",MFNVGetId(this,e));
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  data=(struct MFAUTOBVNVectorData*)MFNVectorGetData(this,e);

  data->bpbv=bpbv;

  return;
 }

void MFAUTOBVNVSetSPBV(MFNVector this,doublereal spbv, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOBVNVSetSPBV"};
  struct MFAUTOBVNVectorData *data;
  int i;

#ifdef MFNOCONFIDENCE
  if(this==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(strcmp(MFNVGetId(this,e),"AUTO"))
   {
    sprintf(MFNVectorErrorMsg,"Trying to set the Periodic Bifurcation Point value in a non-AUTO Vector type \"%s\"",MFNVGetId(this,e));
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  data=(struct MFAUTOBVNVectorData*)MFNVectorGetData(this,e);

  data->spbv=spbv;

  return;
 }

void MFAUTOBVNVSetP0(MFNVector this,doublereal **p0, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOBVNVSetP0"};
  struct MFAUTOBVNVectorData *data;
  int i;

#ifdef MFNOCONFIDENCE
  if(this==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(strcmp(MFNVGetId(this,e),"AUTO"))
   {
    sprintf(MFNVectorErrorMsg,"Trying to set P0 in a non-AUTO Vector type \"%s\"",MFNVGetId(this,e));
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  data=(struct MFAUTOBVNVectorData*)MFNVectorGetData(this,e);

  data->p0=p0;

  return;
 }

void MFAUTOBVNVSetP1(MFNVector this,doublereal **p1, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOBVNVSetP1"};
  struct MFAUTOBVNVectorData *data;
  int i;

#ifdef MFNOCONFIDENCE
  if(this==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(strcmp(MFNVGetId(this,e),"AUTO"))
   {
    sprintf(MFNVectorErrorMsg,"Trying to set P1 in a non-AUTO Vector type \"%s\"",MFNVGetId(this,e));
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  data=(struct MFAUTOBVNVectorData*)MFNVectorGetData(this,e);

  data->p1=p1;

  return;
 }

int MFAUTOBVNVGetNUZ(MFNVector this, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOBVNVGetNUZ"};
  struct MFAUTOBVNVectorData *data;
  int i;

#ifdef MFNOCONFIDENCE
  if(this==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return -1;
   }

  if(strcmp(MFNVGetId(this,e),"AUTO"))
   {
    sprintf(MFNVectorErrorMsg,"Trying to get the number of user zeroes from non-AUTO Vector type \"%s\"",MFNVGetId(this,e));
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return -1;
   }
#endif

  data=(struct MFAUTOBVNVectorData*)MFNVectorGetData(this,e);

  return data->nuz;
 }

doublereal MFAUTOBVNVGetUZBV(MFNVector this,int i, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOBVNVGetUZBV"};
  struct MFAUTOBVNVectorData *data;

#ifdef MFNOCONFIDENCE
  if(this==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return 0.;
   }

  if(strcmp(MFNVGetId(this,e),"AUTO"))
   {
    sprintf(MFNVectorErrorMsg,"Trying to get user zero from non-AUTO Vector type \"%s\"",MFNVGetId(this,e));
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return 0.;
   }
#endif

  data=(struct MFAUTOBVNVectorData*)MFNVectorGetData(this,e);

  if(i<0 || i>=data->nuz)
   {
    sprintf(MFNVectorErrorMsg,"Trying to get an invalid user zero value %d must be in the range [0,%d)",i,data->nuz);
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return 0.;
   }

  return data->uzbv[i];
 }

doublereal MFAUTOBVNVGetLPBV(MFNVector this, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOBVNVGetLPBV"};
  struct MFAUTOBVNVectorData *data;
  int i;

  if(this==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return 0.;
   }

  if(strcmp(MFNVGetId(this,e),"AUTO"))
   {
    sprintf(MFNVectorErrorMsg,"Trying to get limit point bifurcation value from a non-AUTO Vector type \"%s\"",MFNVGetId(this,e));
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return 0.;
   }
  data=(struct MFAUTOBVNVectorData*)MFNVectorGetData(this,e);

  return data->lpbv;
 }

doublereal MFAUTOBVNVGetBPBV(MFNVector this, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOBVNVGetBPBV"};
  struct MFAUTOBVNVectorData *data;
  int i;

#ifdef MFNOCONFIDENCE
  if(this==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return 0.;
   }

  if(strcmp(MFNVGetId(this,e),"AUTO"))
   {
    sprintf(MFNVectorErrorMsg,"Trying to get bifurcation point bifurcation value from a non-AUTO Vector type \"%s\"",MFNVGetId(this,e));
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return 0.;
   }
#endif

  data=(struct MFAUTOBVNVectorData*)MFNVectorGetData(this,e);

  return data->bpbv;
 }

doublereal MFAUTOBVNVGetSPBV(MFNVector this, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOBVNVGetSPBV"};
  struct MFAUTOBVNVectorData *data;
  int i;

#ifdef MFNOCONFIDENCE
  if(this==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return 0.;
   }

  if(strcmp(MFNVGetId(this,e),"AUTO"))
   {
    sprintf(MFNVectorErrorMsg,"Trying to get periodic bifurcation value from a non-AUTO Vector type \"%s\"",MFNVGetId(this,e));
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return 0.;
   }
#endif

  data=(struct MFAUTOBVNVectorData*)MFNVectorGetData(this,e);

  return data->spbv;
 }

doublereal** MFAUTOBVNVGetP0(MFNVector this, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOBVNVGetP0"};
  struct MFAUTOBVNVectorData *data;
  int i;

#ifdef MFNOCONFIDENCE
  if(this==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return NULL;
   }

  if(strcmp(MFNVGetId(this,e),"AUTO"))
   {
    sprintf(MFNVectorErrorMsg,"Trying to get P0 from a non-AUTO Vector type \"%s\"",MFNVGetId(this,e));
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return NULL;
   }
#endif

  data=(struct MFAUTOBVNVectorData*)MFNVectorGetData(this,e);

  return data->p0;
 }

doublereal** MFAUTOBVNVGetP1(MFNVector this, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOBVNVGetP1"};
  struct MFAUTOBVNVectorData *data;
  int i;

#ifdef MFNOCONFIDENCE
  if(this==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return NULL;
   }

  if(strcmp(MFNVGetId(this,e),"AUTO"))
   {
    sprintf(MFNVectorErrorMsg,"Trying to get P1 from a non-AUTO Vector type \"%s\"",MFNVGetId(this,e));
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return NULL;
   }
#endif

  data=(struct MFAUTOBVNVectorData*)MFNVectorGetData(this,e);

  return data->p1;
 }

void MFAUTOBVNVCopyDataValues(MFNVector this,MFNVector that, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOBVNVCopyDataValues"};
  struct MFAUTOBVNVectorData *thisdata;
  struct MFAUTOBVNVectorData *thatdata;
  int i,j;

#ifdef MFNOCONFIDENCE
  if(this==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(strcmp(MFNVGetId(this,e),"AUTO"))
   {
    sprintf(MFNVectorErrorMsg,"Trying to copy data values from a non-AUTO Vector type \"%s\"",MFNVGetId(this,e));
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  thisdata=(struct MFAUTOBVNVectorData*)MFNVectorGetData(this,e);

#ifdef MFNOCONFIDENCE
  if(strcmp(MFNVGetId(that,e),"AUTO"))
   {
    sprintf(MFNVectorErrorMsg,"Trying to copy data values from a non-AUTO Vector type \"%s\"",MFNVGetId(that,e));
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  thatdata=(struct MFAUTOBVNVectorData*)MFNVectorGetData(that,e);

#ifdef MFNOCONFIDENCE
  if(thisdata->nC!=thatdata->nC || thisdata->ntst!=thatdata->ntst || thisdata->ndim!=thatdata->ndim
  || thisdata->ncol!=thatdata->ncol || thisdata->npar!=thatdata->npar)
   {
    sprintf(MFNVectorErrorMsg,"Trying to copy data values from an incompatible AUTO Vector.\n   nC %d and %d\n  ndim %d and %d\n   ncol %d and %d\n   ntst %d and %d\n   npar %d and %d.",thisdata->nC,thatdata->nC,thisdata->ndim,thatdata->ndim,thisdata->ncol,thatdata->ncol,thisdata->ntst,thatdata->ntst,thisdata->npar,thatdata->npar);
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  for(i=0;i<thisdata->ntst+1;i++)
    for(j=0;j<thisdata->ndim*thisdata->ncol;j++)
      (thatdata->u[i])[j]=(thisdata->u[i])[j];

  for(i=0;i<thisdata->npar;i++)thatdata->par[i]=thisdata->par[i];
  for(i=0;i<thisdata->ntst+1;i++)thatdata->t[i]=thisdata->t[i];
  for(i=0;i<thisdata->ntst  ;i++)thatdata->dt[i]=thisdata->dt[i];

/* Bifurcation detection values */
  if(thatdata->nuz!=thisdata->nuz)
   {
    thatdata->nuz=thisdata->nuz;
    if(thatdata->uzbv!=NULL)free(thatdata->uzbv);
    thatdata->uzbv=malloc(thatdata->nuz*sizeof(double));

#ifndef MFNOSAFETYNET
    if(thatdata->uzbv==NULL && thatdata->nuz>0)
     {
      sprintf(MFNVectorErrorMsg,"Out of memory, trying to allocate %d bytes",thatdata->nuz*sizeof(double));
      MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

   }

  for(i=0;i<thisdata->nuz;i++)thatdata->uzbv[i]=thisdata->uzbv[i];
  thatdata->lpbv=thisdata->lpbv;
  thatdata->bpbv=thisdata->bpbv;
  thatdata->spbv=thisdata->spbv;

  for(i=0;i<thisdata->ndim;i++)thatdata->ev[i]=thisdata->ev[i];

  for(i=0;i<thisdata->ndim;i++)
    for(j=0;j<thisdata->ndim;j++)
      (thatdata->p0[i])[j]=(thisdata->p0[i])[j];

  for(i=0;i<thisdata->ndim;i++)
    for(j=0;j<thisdata->ndim;j++)
      (thatdata->p1[i])[j]=(thisdata->p1[i])[j];

  return;
 }

int MFAUTOBVNVGetNtst(MFNVector this, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOBVNVGetNtst"};
  struct MFAUTOBVNVectorData *data;
  int i,j;

#ifdef MFNOCONFIDENCE
  if(this==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(strcmp(MFNVGetId(this,e),"AUTO"))
   {
    sprintf(MFNVectorErrorMsg,"Trying to copy data values from a non-AUTO Vector type \"%s\"",MFNVGetId(this,e));
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
  }
#endif

  data=(struct MFAUTOBVNVectorData*)MFNVectorGetData(this,e);

  return data->ntst;
 }

int MFAUTOBVNVGetNcol(MFNVector this, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOBVNVGetNcol"};
  struct MFAUTOBVNVectorData *data;
  int i,j;

#ifdef MFNOCONFIDENCE
  if(this==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(strcmp(MFNVGetId(this,e),"AUTO"))
   {
    sprintf(MFNVectorErrorMsg,"Trying to copy data values from a non-AUTO Vector type \"%s\"",MFNVGetId(this,e));
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
  }
#endif

  data=(struct MFAUTOBVNVectorData*)MFNVectorGetData(this,e);

  return data->ncol;
 }

int MFAUTOBVNVGetNdim(MFNVector this, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOBVNVGetNdim"};
  struct MFAUTOBVNVectorData *data;
  int i,j;

#ifdef MFNOCONFIDENCE
  if(this==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(strcmp(MFNVGetId(this,e),"AUTO"))
   {
    sprintf(MFNVectorErrorMsg,"Trying to copy data values from a non-AUTO Vector type \"%s\"",MFNVGetId(this,e));
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
  }
#endif

  data=(struct MFAUTOBVNVectorData*)MFNVectorGetData(this,e);

  return data->ndim;
 }

int MFAUTOBVNVGetNpar(MFNVector this, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOBVNVGetNpar"};
  struct MFAUTOBVNVectorData *data;
  int i,j;

#ifdef MFNOCONFIDENCE
  if(this==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(strcmp(MFNVGetId(this,e),"AUTO"))
   {
    sprintf(MFNVectorErrorMsg,"Trying to copy data values from a non-AUTO Vector type \"%s\"",MFNVGetId(this,e));
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
  }
#endif

  data=(struct MFAUTOBVNVectorData*)MFNVectorGetData(this,e);

  return data->npar;
 }

double MFAUTOBVNVParmDot(MFNVector u,MFNVector v, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOBVNVParmDot"};
  struct MFAUTOBVNVectorData *udata;
  struct MFAUTOBVNVectorData *vdata;
  int i;
  double result;

  udata=(struct MFAUTOBVNVectorData*)MFNVectorGetData(u,e);
  vdata=(struct MFAUTOBVNVectorData*)MFNVectorGetData(v,e);

  result=0;
  for(i=0;i<udata->npar;i++)
    result+=udata->par[i]*vdata->par[i];

  return result;
 }

void MFPrintAUTOBVNVectorFull(FILE *fid,MFNVector u0, MFErrorHandler e)
 {
  static char RoutineName[]={"MFPrintAUTOBVNVectorFull"};

  MFPrintAUTOBVNVectorDataFull(fid,MFNVectorGetData(u0,e),e);

  return;
 }

void MFPrintAUTOBVNVectorDataFull(FILE *fid,void *u0, MFErrorHandler e)
 {
  static char RoutineName[]={"MFPrintAUTOBVNVectorDataFull"};
  double unorm;
  int i;
  int it,icol,idim;
  struct MFAUTOBVNVectorData *u;

#ifdef MFNOCONFIDENCE
  if(fid==NULL)
   {
    sprintf(MFNVectorErrorMsg,"fid (argument 1) is NULL.");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  u=(struct MFAUTOBVNVectorData*)u0;
  fprintf(fid,"\n");
  fprintf(fid,"          t          dt");
  for(idim=0;idim<u->ndim;idim++)fprintf(fid,"       u[%d]",idim);
  fprintf(fid,"\n");
  for(it=0;it<u->ntst;it++)
   {
    for(icol=0;icol<u->ncol;icol++)
     {
      fprintf(fid,"%2d %2d ",it,icol);
      fprintf(fid,"%10.7lf ",(u->t)[it]+icol*(u->dt)[it]/u->ncol);
      fprintf(fid,"%lf ",(u->dt)[it]);
      for(idim=0;idim<u->ndim;idim++)fprintf(fid,"%10.7lf ",(u->u)[it][idim+u->ndim*icol]);
      fprintf(fid,"\n");fflush(fid);
     }
   }
  it=u->ntst;
  icol=0;
  fprintf(fid,"%2d %2d ",it,icol);
  fprintf(fid,"%10.7lf ",(u->t)[it]+icol*(u->dt)[it]/u->ncol);
  fprintf(fid,"         ");
  for(idim=0;idim<u->ndim;idim++)fprintf(fid,"%10.7lf ",(u->u)[it][idim+u->ndim*icol]);
  fprintf(fid,"\n");fflush(fid);

  fprintf(fid,"                          ");
  for(i=0;i<u->nfpr;i++)
   {
    if(i>0)fprintf(fid,",");
    fprintf(fid,"%10.7lf",u->par[i]);
   }
  fprintf(fid,"\n");fflush(stdout);

  return;
 }

void MFAUTOBVNVSetType(MFNVector this, char *type, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOBVNVSetType"};
  struct MFAUTOBVNVectorData *data;
  int i;

#ifdef MFNOCONFIDENCE
  if(this==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(strcmp(MFNVGetId(this,e),"AUTO"))
   {
    sprintf(MFNVectorErrorMsg,"Trying to get arrays from non-AUTO Vector type \"%s\"",MFNVGetId(this,e));
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  data=(struct MFAUTOBVNVectorData*)MFNVectorGetData(this,e);

  strncpy(data->type,type,2);
  data->type[3]=0x0;

  return;
 }

int MFAUTOBVNVGetNit(MFNVector this, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOBVNVGetNit"};
  struct MFAUTOBVNVectorData *data;
  int i,j;

#ifdef MFNOCONFIDENCE
  if(this==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(strcmp(MFNVGetId(this,e),"AUTO"))
   {
    sprintf(MFNVectorErrorMsg,"Trying to copy data values from a non-AUTO Vector type \"%s\"",MFNVGetId(this,e));
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  data=(struct MFAUTOBVNVectorData*)MFNVectorGetData(this,e);

  return data->nit;
 }

double MFAUTOBVNVGetR(MFNVector this, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOBVNVGetR"};
  struct MFAUTOBVNVectorData *data;
  int i,j;

#ifdef MFNOCONFIDENCE
  if(this==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return 0.;
   }

  if(strcmp(MFNVGetId(this,e),"AUTO"))
   {
    sprintf(MFNVectorErrorMsg,"Trying to copy data values from a non-AUTO Vector type \"%s\"",MFNVGetId(this,e));
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return 0.;
   }
#endif

  data=(struct MFAUTOBVNVectorData*)MFNVectorGetData(this,e);

  return data->R;
 }

void MFAUTOBVNVSetR(MFNVector this, double R, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOBVNVSetR"};
  struct MFAUTOBVNVectorData *data;
  int i;

#ifdef MFNOCONFIDENCE
  if(this==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(strcmp(MFNVGetId(this,e),"AUTO"))
   {
    sprintf(MFNVectorErrorMsg,"Trying to get arrays from non-AUTO Vector type \"%s\"",MFNVGetId(this,e));
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  data=(struct MFAUTOBVNVectorData*)MFNVectorGetData(this,e);

  data->R=R;

  return;
 }

void MFAUTOBVNVSetNit(MFNVector this, int nit, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOBVNVSetNit"};
  struct MFAUTOBVNVectorData *data;
  int i;

#ifdef MFNOCONFIDENCE
  if(this==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(strcmp(MFNVGetId(this,e),"AUTO"))
   {
    sprintf(MFNVectorErrorMsg,"Trying to get arrays from non-AUTO Vector type \"%s\"",MFNVGetId(this,e));
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  data=(struct MFAUTOBVNVectorData*)MFNVectorGetData(this,e);

  data->nit=nit;

  return;
 }

long *MFAUTOBVNVGetICP(MFNVector this, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOBVNVGetNpar"};
  struct MFAUTOBVNVectorData *data;
  int i,j;

#ifdef MFNOCONFIDENCE
  if(this==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(strcmp(MFNVGetId(this,e),"AUTO"))
   {
    sprintf(MFNVectorErrorMsg,"Trying to copy data values from a non-AUTO Vector type \"%s\"",MFNVGetId(this,e));
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
  }
#endif

  data=(struct MFAUTOBVNVectorData*)MFNVectorGetData(this,e);

  return data->icp;
 }

int MFAUTOBVNVGetNfpr(MFNVector this, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOBVNVGetNfpr"};
  struct MFAUTOBVNVectorData *data;
  int i,j;

#ifdef MFNOCONFIDENCE
  if(this==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(strcmp(MFNVGetId(this,e),"AUTO"))
   {
    sprintf(MFNVectorErrorMsg,"Trying to copy data values from a non-AUTO Vector type \"%s\"",MFNVGetId(this,e));
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
  }
#endif

  data=(struct MFAUTOBVNVectorData*)MFNVectorGetData(this,e);

  return data->nfpr;
 }

void MFAUTOBVNVSetThu(MFNVector this, doublereal *thu, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOBVNVSetNit"};
  struct MFAUTOBVNVectorData *data;
  int i;

#ifdef MFNOCONFIDENCE
  if(this==NULL)
   {
    sprintf(MFNVectorErrorMsg,"Pointer to Vector (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(strcmp(MFNVGetId(this,e),"AUTO"))
   {
    sprintf(MFNVectorErrorMsg,"Trying to get arrays from non-AUTO Vector type \"%s\"",MFNVGetId(this,e));
    MFSetError(e,12,RoutineName,MFNVectorErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  data=(struct MFAUTOBVNVectorData*)MFNVectorGetData(this,e);

  data->thu=thu;

  return;
 }

#ifdef __cplusplus
}
#endif
