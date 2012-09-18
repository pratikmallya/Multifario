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

static char *id="@(#) $Id: MFTPBVPNSpace.c,v 1.3 2007/02/13 01:22:34 mhender Exp $";

static char MFTPBVPNSpaceErrorMsg[256]="";

#include <MFNSpace.h>
#include <MFNVector.h>
#include <MFErrorHandler.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef __cplusplus
 extern "C" {
#endif

struct MFTPBVPNSpaceData
 {
  int n;
  int nx;
  int nu;
  int np;
  int *parmPeriodic;
  double *parmPeriod;
 };

static void MFFreeTPBVPNSpaceData(void*,MFErrorHandler);
static double MFTPBVPNSpaceDistance(MFNSpace,MFNVector,MFNVector,void*,MFErrorHandler);
static void MFTPBVPNSpaceDirection(MFNSpace,MFNVector,MFNVector,MFNVector,void*,MFErrorHandler);
static void MFTPBVPNSpaceAdd(MFNSpace,MFNVector,MFNVector,MFNVector,void*,MFErrorHandler);
static void MFTPBVPNSpaceScale(MFNSpace,double,MFNVector,MFNVector,void*,MFErrorHandler);
static double MFTPBVPNSpaceInner(MFNSpace,MFNVector,MFNVector,void*,MFErrorHandler);
static void MFWriteTPBVPNSpaceData(FILE*,MFNSpace,void*,MFErrorHandler);
static void MFReadTPBVPNSpaceData(FILE*,MFNSpace,MFErrorHandler);

MFNSpace MFCreateTPBVPNSpace(int nx, int nu, int np, MFErrorHandler e)
 {
  static char RoutineName[]={"MFCreateTPBVPNSpace"};
  MFNSpace thisSpace;
  struct MFTPBVPNSpaceData *data;
  int i;

  thisSpace=MFCreateNSpaceBaseClass("TPBVPNSpace",e);
  MFNSpaceSetDistance(thisSpace,MFTPBVPNSpaceDistance,e);
  MFNSpaceSetInnerProduct(thisSpace,MFTPBVPNSpaceInner,e);
  MFNSpaceSetDirection(thisSpace,MFTPBVPNSpaceDirection,e);
  MFNSpaceSetAdd(thisSpace,MFTPBVPNSpaceAdd,e);
  MFNSpaceSetScale(thisSpace,MFTPBVPNSpaceScale,e);
  MFNSpaceSetFreeData(thisSpace,MFFreeTPBVPNSpaceData,e);
  MFNSpaceSetWriteData(thisSpace,MFWriteTPBVPNSpaceData,e);

  data=(struct MFTPBVPNSpaceData*)malloc(sizeof(struct MFTPBVPNSpaceData));

#ifndef MFNOSAFETYNET
  if(data==NULL)
   {
    sprintf(MFTPBVPNSpaceErrorMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct MFTPBVPNSpaceData));
    MFSetError(e,12,RoutineName,MFTPBVPNSpaceErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  data->n=nx*nu+np+nx;
  data->nx=nx;
  data->nu=nu;
  data->np=np;
  data->parmPeriodic=(int*)malloc(np*sizeof(int));

#ifndef MFNOSAFETYNET
  if(data->parmPeriodic==NULL)
   {
    sprintf(MFTPBVPNSpaceErrorMsg,"Out of memory, trying to allocate %d bytes",np*sizeof(int));
    MFSetError(e,12,RoutineName,MFTPBVPNSpaceErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  data->parmPeriod=(double*)malloc(np*sizeof(double));

#ifndef MFNOSAFETYNET
  if(data->parmPeriod==NULL)
   {
    sprintf(MFTPBVPNSpaceErrorMsg,"Out of memory, trying to allocate %d bytes",np*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPNSpaceErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  for(i=0;i<np;i++)
   {
    data->parmPeriodic[i]=0;
    data->parmPeriod[i]=0.;
   }
  
  MFNSpaceSetData(thisSpace,(void*)data,e);

  return thisSpace;
 }

void MFFreeTPBVPNSpaceData(void *data, MFErrorHandler e)
 {
  static char RoutineName[]={"MFFreeTPBVPNSpace"};
  struct MFTPBVPNSpaceData *d;

  d=(struct MFTPBVPNSpaceData*)data;
  if(d!=NULL)
   {
    if(d->parmPeriodic!=NULL)free(d->parmPeriodic);
    if(d->parmPeriod!=NULL)free(d->parmPeriod);
    free(d);
   }

  return;
 }

double MFTPBVPNSpaceDistance(MFNSpace thisSpace,MFNVector v0,MFNVector v1,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFNSpaceDistance"};
  double result;
  int i,j,o;
  int nx,nu,np;
  double t;
  double *u0;
  double *u1;
  struct MFTPBVPNSpaceData *data;

  data=(struct MFTPBVPNSpaceData*)d;
  nx=data->nx;
  nu=data->nu;
  np=data->np;

  result=0.;
  o=nx*nu+np;

  u0=MFNV_CStar(v0,e);
  u1=MFNV_CStar(v1,e);

  for(j=0;j<nu;j++)
    result+=(u1[j]-u0[j])*(u1[j]-u0[j])*u0[o]/2.;
  for(i=1;i<nx-1;i++)
   {
    for(j=0;j<nu;j++)
      result+=(u1[j+nu*i]-u0[j+nu*i])*(u1[j+nu*i]-u0[j+nu*i])*u0[o+i];
   }
  for(j=0;j<nu;j++)
    result+=(u1[j+(nx-1)*nu]-u0[j+(nx-1)*nu])*(u1[j+(nx-1)*nu]-u0[j+(nx-1)*nu])*u0[o+nx-1]/2.;

  for(i=0;i<np;i++)
   {
    if(!data->parmPeriodic[i])
     {
      result+=(u1[nx*nu+i]-u0[nx*nu+i])*(u1[nx*nu+i]-u0[nx*nu+i]);
     }else{
      t=fabs(u1[nx*nu+i]-u0[nx*nu+i]);
      if(t>data->parmPeriod[i]/2.)t=t-data->parmPeriod[i]/2.;
      result+=t*t;
     }
   }
  result=sqrt(result);

  return result;
 }

void MFTPBVPNSpaceDirection(MFNSpace thisSpace,MFNVector v0,MFNVector v1,MFNVector diff,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFNSpaceDirection"};
  int i;
  int nx,nu,np;
  struct MFTPBVPNSpaceData *data;
  double t;
  double *u0;
  double *u1;
  double *di;

  data=(struct MFTPBVPNSpaceData*)d;
  nx=data->nx;
  nu=data->nu;
  np=data->np;

  u0=MFNV_CStar(v0,e);
  u1=MFNV_CStar(v1,e);
  di=MFNV_CStar(diff,e);

  for(i=0;i<nx*nu;i++)di[i]=u1[i]-u0[i];

  for(i=0;i<np;i++)
   {
    if(!data->parmPeriodic[i])
     {
      di[nx*nu+i]=u1[nx*nu+i]-u0[nx*nu+i];
     }else{
      t=fabs(u1[nx*nu+i]-u0[nx*nu+i]);
      if(t<data->parmPeriod[i]/2.)
       {
        di[nx*nu+i]=u1[nx*nu+i]-u0[nx*nu+i];
       }else{
        di[nx*nu+i]=u0[nx*nu+i]-u1[nx*nu+i];
       }
     }
   }

  for(i=nx*nu+np;i<nx*nu+np+nx;i++)di[i]=u1[i];

  return;
 }

double MFTPBVPNSpaceInner(MFNSpace thisSpace,MFNVector v0,MFNVector v1,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFNSpaceInner"};
  double result;
  int i,j,o;
  int nx,nu,np,n;
  int verbose;
  struct MFTPBVPNSpaceData *data;
  double *u0;
  double *u1;

  data=(struct MFTPBVPNSpaceData*)d;
  n=data->n;
  nx=data->nx;
  nu=data->nu;
  np=data->np;

  verbose=0;

  if(verbose){printf("%s, n=%d, nx=%d, nu=%d, np=%d, n should be %d\n",RoutineName,n,nx,nu,np,nx*nu+np+nx);fflush(stdout);}

  data=(struct MFTPBVPNSpaceData*)d;
  u0=MFNV_CStar(v0,e);
  u1=MFNV_CStar(v1,e);

  o=nx*nu+np;
  result=0.;
  for(j=0;j<nu;j++)
    result+=u0[j]*u1[j]*u0[o]/2;
  for(i=1;i<nx-1;i++)
   {
    for(j=0;j<nu;j++)
      result+=u0[j+nu*i]*u1[j+nu*i]*u0[o+i];
   }
  for(j=0;j<nu;j++)
    result+=u0[j+(nx-1)*nu]*u1[j+(nx-1)*nu]*u0[o+nx-1]/2;

  o=nx*nu;
  for(i=0;i<np;i++)
    result+=u0[o+i]*u1[o+i];

  return result;
 }

void MFWriteTPBVPNSpaceData(FILE *fid,MFNSpace thisSpace,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFWriteTPBVPNSpaceData"};
  int i;
  struct MFTPBVPNSpaceData *data;

  data=(struct MFTPBVPNSpaceData*)d;

  fprintf(fid,"%d %d %d %d\n",data->n,data->nx,data->nu,data->np);

  return;
 }

void MFReadTPBVPNSpaceData(FILE *fid, MFNSpace thisSpace, MFErrorHandler e)
 {
  static char RoutineName[]={"MFTPBVPReadNSpaceData"};
  int i;
  struct MFTPBVPNSpaceData *data;

  data=(struct MFTPBVPNSpaceData*)malloc(sizeof(struct MFTPBVPNSpaceData));

#ifndef MFNOSAFETYNET
  if(data==NULL)
   {
    sprintf(MFTPBVPNSpaceErrorMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct MFTPBVPNSpaceData));
    MFSetError(e,12,RoutineName,MFTPBVPNSpaceErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif


  fscanf(fid,"%d %d %d %d\n",&(data->n),&(data->nx),&(data->nu),&(data->np));
  MFNSpaceSetData(thisSpace,data,e);

  return;
 }

void MFTPBVPNSpaceAdd(MFNSpace thisSpace,MFNVector v0,MFNVector v1,MFNVector sum,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFNSpaceAdd"};
  int i;
  int nx,nu,np;
  struct MFTPBVPNSpaceData *data;
  double *u0,*u1,*v;

  data=(struct MFTPBVPNSpaceData*)d;
  nx=data->nx;
  nu=data->nu;
  np=data->np;
  u0=MFNV_CStar(v0,e);
  u1=MFNV_CStar(v1,e);
  v =MFNV_CStar(sum,e);

  for(i=0;i<nx*nu+np;i++)
    v[i]=u1[i]+u0[i];
  for(i=nx*nu+np;i<nx*nu+np+nx;i++)
    v[i]=u1[i];

  return;
 }

void MFTPBVPNSpaceScale(MFNSpace thisSpace,double s, MFNVector v,MFNVector prod,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFNSpaceScale"};
  int i;
  int nx,nu,np;
  struct MFTPBVPNSpaceData *data;
  double *u,*p;

#ifdef MFNOCONFIDENCE
  if(thisSpace==NULL)
   {
    sprintf(MFTPBVPNSpaceErrorMsg,"Space (arg 1) is NULL!");
    MFSetError(e,12,RoutineName,MFTPBVPNSpaceErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(v==NULL)
   {
    sprintf(MFTPBVPNSpaceErrorMsg,"v (arg 4) is NULL!");
    MFSetError(e,12,RoutineName,MFTPBVPNSpaceErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(prod==NULL)
   {
    sprintf(MFTPBVPNSpaceErrorMsg,"prod (arg 5) is NULL!");
    MFSetError(e,12,RoutineName,MFTPBVPNSpaceErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  data=(struct MFTPBVPNSpaceData*)d;
  nx=data->nx;
  nu=data->nu;
  np=data->np;
  u=MFNV_CStar(v,e);
  p=MFNV_CStar(prod,e);

  for(i=0;i<nx*nu+np;i++)p[i]=s*u[i];
  for(i=nx*nu+np;i<nx*nu+np+nx;i++)p[i]=u[i];

  return;
 }

void MFTPBVPNSpaceSetPeriodicParameter(MFNSpace thisSpace,int p,double T, MFErrorHandler e)
 {
  struct MFTPBVPNSpaceData *data;
  static char RoutineName[]={"MFTPBVPNSpaceSetPeriodicParameter"};

#ifdef MFNOCONFIDENCE
  if(strcmp("TPBVPNSpace",MFNSpaceGetId(thisSpace,e)))
   {
    sprintf(MFTPBVPNSpaceErrorMsg,"Space not a TPBVPNSpace!");
    MFSetError(e,12,RoutineName,MFTPBVPNSpaceErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  data=(struct MFTPBVPNSpaceData*)MFNSpaceGetData(thisSpace,e);

#ifdef MFNOCONFIDENCE
  if(p<0||p>=data->np)
   {
    sprintf(MFTPBVPNSpaceErrorMsg,"Parameter %d (argument 2) is invalid must be in range [0,%d]!",p,data->np);
    MFSetError(e,12,RoutineName,MFTPBVPNSpaceErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  data->parmPeriodic[p]=1;
  data->parmPeriod[p]=T;

  return;
 }

#ifdef __cplusplus
}
#endif
