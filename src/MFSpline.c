/*
 *  PROGRAM NAME:  multifario
 *
 *  (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
 *  CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
 *
 *  Please refer to the LICENSE file in the top directory
 *
 *      author: Mike Henderson mhender@watson.ibm.com
 *      date:   August 26, 1999 (MFExpression)
 *              February 25, 2003 changed name and added subroutine support
 */

static char *id="@(#) $Id: MFSpline.c,v 1.2 2011/07/21 17:42:46 mhender Exp $";

static char MFSplineErrorMsg[256]="";

extern double MFEpsilon;

#include <multifarioConfig.h>
#include <MFImplicitMF.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <ExpCmp.h>
#include <MFFortran.h>
#include <MFErrorHandler.h>
#define round(x) ((x)>=0?(long)((x)+0.5):(long)((x)-0.5))

#ifndef MAX
#define MAX(X, Y) ((X) < (Y) ? (Y) : (X))
#endif

#ifndef DBL_QNAN
#define DBL_QNAN 1.e200
#endif

#ifdef __cplusplus
 extern "C" {
#endif

static void MFFreeSplineData(void*,MFErrorHandler);
static int MFProjectSpline(int,int,MFNVector,MFNKMatrix,MFNVector,void*,int*,MFErrorHandler);
static int MFTangentSpline(int,int,MFNVector,MFNKMatrix,void*,MFErrorHandler);
static double MFScaleSpline(int,int,MFNVector,MFNKMatrix,void*,MFErrorHandler);
static void MFWriteSplineData(FILE*,void*,MFErrorHandler);
static MFImplicitMF MFReadSpline(FILE*,MFErrorHandler);
static double *MFSplineNV_CStar(MFNVector,MFErrorHandler);

static int MFSplineProjectToSave(MFNVector,double*,void*,MFErrorHandler);
static void MFSplineSetStability(MFImplicitMF,MFNVector,MFNKMatrix,void*,MFErrorHandler);
int MFStopSpline(MFImplicitMF,MFNVector,MFNKMatrix,MFNVector,MFNKMatrix,void*,MFErrorHandler);
static double MFFindS(double*,void*,MFErrorHandler);
MFNVector MFCreateSplineNVector(int,MFErrorHandler);
MFNSpace MFCreateSplineNSpace(MFErrorHandler);

struct MFFitpackSplineData
 {
  int    npts;
  int    n;
  double *temp;
  double **yp;
  double  sigma;
  double **x;
  double *s;
  double  p;
  int    periodic;
 };

void CALLFCURV2(double *t,int *n,double *x,double *y,double *yp,double *sigma, double *ans);
void CALLFCURVP2(double *t,int *n,double *x,double *y,double *p,double *yp,double *sigma, double *ans);

static void spline(double s,double *u, struct MFFitpackSplineData *data);
static void dspline(double s, double *du, struct MFFitpackSplineData *data);
static void ddspline(double s, double *ddu, struct MFFitpackSplineData *data);

MFNVector MFSplineNVectorFactory(MFImplicitMF,MFErrorHandler);
MFNKMatrix MFNKMatrixFactory(MFImplicitMF,MFErrorHandler);
void MFSplineNVectorSetS(MFNVector,double,MFErrorHandler);
double MFSplineNVectorGetS(MFNVector,MFErrorHandler);

MFImplicitMF MFIMFCreatePeriodicSpline(int npts, int n, double **x, MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFCreatePeriodicSpline"};
  MFImplicitMF thisMF;
  MFNSpace space;
  double ds;
  int i,j,k;
  int ierr;
  struct MFFitpackSplineData *data;
  double t;
  double z;
  double u[4];

  k=1;

  thisMF=MFIMFCreateBaseClass(n,k,"PeriodicSpline",e);

  space=MFCreateSplineNSpace(e);
  MFIMFSetSpace(thisMF,space,e);
  MFFreeNSpace(space,e);

  data=(struct MFFitpackSplineData*)malloc(sizeof(struct MFFitpackSplineData));

#ifndef MFNOSAFETYNET
  if(data==NULL)
   {
    sprintf(MFSplineErrorMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct MFFitpackSplineData));
    MFSetError(e,12,RoutineName,MFSplineErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(thisMF);
    return NULL;
   }
#endif

  data->npts=npts;
  data->n   =n;
  data->periodic=1;

  data->temp=(double*)malloc(2*((data->npts+1)+1)*sizeof(double));

#ifndef MFNOSAFETYNET
  if(data->temp==NULL)
   {
    sprintf(MFSplineErrorMsg,"Out of memory, trying to allocate %d bytes",(2*(data->npts+1)+1)*sizeof(double));
    MFSetError(e,12,RoutineName,MFSplineErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(thisMF);
    return NULL;
   }
#endif

  data->x=(double**)malloc(n*sizeof(double*));

#ifndef MFNOSAFETYNET
  if(data->x==NULL)
   {
    sprintf(MFSplineErrorMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(double));
    MFSetError(e,12,RoutineName,MFSplineErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(thisMF);
    return NULL;
   }
#endif

  data->yp=(double**)malloc(n*sizeof(double*));

#ifndef MFNOSAFETYNET
  if(data->yp==NULL)
   {
    sprintf(MFSplineErrorMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(double));
    MFSetError(e,12,RoutineName,MFSplineErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(thisMF);
    return NULL;
   }
#endif

  for(i=0;i<n;i++)
   {
    data->x[i]=(double*)malloc((data->npts+2)*sizeof(double));

#ifndef MFNOSAFETYNET
  if(data->x[i]==NULL)
   {
    sprintf(MFSplineErrorMsg,"Out of memory, trying to allocate %d bytes",(data->npts+2)*sizeof(double));
    MFSetError(e,12,RoutineName,MFSplineErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(thisMF);
    return NULL;
   }
#endif

    data->yp[i]=(double*)malloc((data->npts+2)*sizeof(double));

#ifndef MFNOSAFETYNET
  if(data->yp[i]==NULL)
   {
    sprintf(MFSplineErrorMsg,"Out of memory, trying to allocate %d bytes",(data->npts)*sizeof(double));
    MFSetError(e,12,RoutineName,MFSplineErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(thisMF);
    return NULL;
   }
#endif

    for(j=0;j<npts;j++)data->x[i][j]=x[i][j];
   }

  data->s=(double*)malloc((data->npts+1)*sizeof(double));

#ifndef MFNOSAFETYNET
  if(data->s==NULL)
   {
    sprintf(MFSplineErrorMsg,"Out of memory, trying to allocate %d bytes",npts*sizeof(double));
    MFSetError(e,12,RoutineName,MFSplineErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(thisMF);
    return NULL;
   }
#endif

  i=0;
  data->s[i]=0.;
  for(i=1;i<data->npts;i++)
   {
    ds=0.;
    for(j=0;j<n;j++)ds+=(data->x[j][i]-data->x[j][i-1])*(data->x[j][i]-data->x[j][i-1]);
    data->s[i]=data->s[i-1]+sqrt(ds);
   }

  ds=0.;
  for(j=0;j<n;j++)ds+=(data->x[j][0]-data->x[j][data->npts-1])*(data->x[j][0]-data->x[j][data->npts-1]);
  data->p=data->s[data->npts-1]+sqrt(ds);
  for(i=0;i<data->npts;i++)data->s[i]=data->s[i]/data->p;
  data->p=1.;
  data->sigma=0.;


/*
   ierr contains an error flag,
        = 0 for normal return,
        = 1 if n is less than 2,
        = 2 if p is less than or equal to x(n)-x(1),
        = 3 if x-values are not strictly increasing.
 */

  for(j=0;j<n;j++)
   {
    CALLCURVP1(&(data->npts),data->s,data->x[j],&(data->p),data->yp[j],data->temp,&(data->sigma),&ierr);
   }

/*
  for(i=0;i<data->npts;i++)
   {
    printf("%2d %14.7le (",i,data->s[i]);
    for(j=0;j<data->n;j++)
     {
      if(j>0)printf(",");
      printf("%14.7e",data->x[j][i]);
     }
    printf(") (");fflush(stdout);
    for(j=0;j<data->n;j++)
     {
      if(j>0)printf(",");
      printf("%14.7le",data->yp[j][i]);
     }
    printf(")\n");fflush(stdout);
   }
  printf("\n");fflush(stdout);
 
  for(t=0.;t<1.05;t+=.1)
   {
    printf("  %14.7e (",t);
    for(j=0;j<data->n;j++)
     {
      if(j>0)printf(",");
      CALLFCURVP2(&t,&(data->npts),data->s,data->x[j],&(data->p),data->yp[j],&(data->sigma),&z);
      printf("%14.7e",z);
     }
    printf(")\n");fflush(stdout);
   }
  printf("\n");fflush(stdout);
 */

  MFIMFSetData(thisMF,(void*)data,e);
  MFIMFSetFreeData(thisMF,MFFreeSplineData,e);
  MFIMFSetProject(thisMF,MFProjectSpline,e);
  MFIMFSetTangent(thisMF,MFTangentSpline,e);
/*MFIMFSetScale(thisMF,MFScaleSpline,e);*/

  MFIMFSetProjectForSave(thisMF,MFSplineProjectToSave,e);
  MFIMFSetProjectForDraw(thisMF,MFSplineProjectToSave,e);
  MFIMFSetProjectForBB(thisMF,MFSplineProjectToSave,e);

  MFIMFSetVectorFactory(thisMF,MFSplineNVectorFactory,e);
  MFIMFSetMatrixFactory(thisMF,MFNKMatrixFactory,e);

  printf("done %s\n",RoutineName);fflush(stdout);
  return thisMF;
 }

MFImplicitMF MFIMFCreateSpline(int npts, int n, double **x, MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFCreateSpline"};
  MFImplicitMF thisMF;
  MFNSpace space;
  double ds;
  int ierr;
  int i,j,k;
  struct MFFitpackSplineData *data;

  k=1;

  thisMF=MFIMFCreateBaseClass(n,k,"PeriodicSpline",e);

  space=MFCreateSplineNSpace(e);
  MFIMFSetSpace(thisMF,space,e);
  MFFreeNSpace(space,e);

  data=(struct MFFitpackSplineData*)malloc(sizeof(struct MFFitpackSplineData));

#ifndef MFNOSAFETYNET
  if(data==NULL)
   {
    sprintf(MFSplineErrorMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct MFFitpackSplineData));
    MFSetError(e,12,RoutineName,MFSplineErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(thisMF);
    return NULL;
   }
#endif

  data->npts=npts;
  data->n   =n;
  data->periodic=0;

  data->temp=(double*)malloc(2*npts*sizeof(double));

#ifndef MFNOSAFETYNET
  if(data->temp==NULL)
   {
    sprintf(MFSplineErrorMsg,"Out of memory, trying to allocate %d bytes",2*npts*sizeof(double));
    MFSetError(e,12,RoutineName,MFSplineErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(thisMF);
    return NULL;
   }
#endif

  data->x=(double**)malloc(n*sizeof(double*));

#ifndef MFNOSAFETYNET
  if(data->x==NULL)
   {
    sprintf(MFSplineErrorMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(double));
    MFSetError(e,12,RoutineName,MFSplineErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(thisMF);
    return NULL;
   }
#endif

  data->yp=(double**)malloc(n*sizeof(double*));

#ifndef MFNOSAFETYNET
  if(data->yp==NULL)
   {
    sprintf(MFSplineErrorMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(double));
    MFSetError(e,12,RoutineName,MFSplineErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(thisMF);
    return NULL;
   }
#endif

  for(i=0;i<n;i++)
   {
    data->x[i]=(double*)malloc(npts*sizeof(double));

#ifndef MFNOSAFETYNET
  if(data->x[i]==NULL)
   {
    sprintf(MFSplineErrorMsg,"Out of memory, trying to allocate %d bytes",npts*sizeof(double));
    MFSetError(e,12,RoutineName,MFSplineErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(thisMF);
    return NULL;
   }
#endif

    data->yp[i]=(double*)malloc(npts*sizeof(double));

#ifndef MFNOSAFETYNET
  if(data->yp[i]==NULL)
   {
    sprintf(MFSplineErrorMsg,"Out of memory, trying to allocate %d bytes",npts*sizeof(double));
    MFSetError(e,12,RoutineName,MFSplineErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(thisMF);
    return NULL;
   }
#endif

    for(j=0;j<npts;j++)data->x[i][j]=x[i][j];
   }

  data->s=(double*)malloc(npts*sizeof(double));

#ifndef MFNOSAFETYNET
  if(data->s==NULL)
   {
    sprintf(MFSplineErrorMsg,"Out of memory, trying to allocate %d bytes",npts*sizeof(double));
    MFSetError(e,12,RoutineName,MFSplineErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(thisMF);
    return NULL;
   }
#endif

  data->s[0]=0.;
  for(i=1;i<npts;i++)
   {
    ds=0.;
    for(j=0;j<n;j++)ds+=(x[j][i]-x[j][i-1])*(x[j][i]-x[j][i-1]);
    data->s[i]=data->s[i-1]+sqrt(ds);
   }
  data->p=data->s[npts-1];

  for(i=0;i<npts;i++)data->s[i]=data->s[i]/data->p;
  data->p=1.;

  data->sigma=0.001;

  for(j=0;j<n;j++)CALLCURV1(&(data->npts),data->s,data->x[j],&(data->p),data->yp[j],data->temp,&(data->sigma),&ierr);

  MFIMFSetData(thisMF,(void*)data,e);
  MFIMFSetFreeData(thisMF,MFFreeSplineData,e);
  MFIMFSetProject(thisMF,MFProjectSpline,e);
  MFIMFSetTangent(thisMF,MFTangentSpline,e);
  MFIMFSetScale(thisMF,MFScaleSpline,e);
  MFIMFSetStop(thisMF,MFStopSpline,e);

  MFIMFSetProjectForSave(thisMF,MFSplineProjectToSave,e);
  MFIMFSetProjectForDraw(thisMF,MFSplineProjectToSave,e);
  MFIMFSetProjectForBB(thisMF,MFSplineProjectToSave,e);

  MFIMFSetVectorFactory(thisMF,MFSplineNVectorFactory,e);
  MFIMFSetMatrixFactory(thisMF,MFNKMatrixFactory,e);

  return thisMF;
 }

void MFFreeSplineData(void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFFreeSplineData"};
  int i;
  struct MFFitpackSplineData *data;

  data=(struct MFFitpackSplineData*)d;

  if(data->temp!=NULL)free(data->temp);
  if(data->s   !=NULL)free(data->s);

  if(data->x!=NULL)
   {
    for(i=0;i<data->n;i++)
         if(data->x[i]!=NULL)free(data->x[i]);
    free(data->x);
   }

  if(data->yp!=NULL)
   {
    for(i=0;i<data->n;i++)
         if(data->yp[i]!=NULL)free(data->yp[i]);
    free(data->yp);
   }

  free(data);
  return;
 }

int MFProjectSpline(int n,int k,MFNVector vu0,MFNKMatrix mPhi,MFNVector vu,void *d,int *index, MFErrorHandler e)
 {
  static char RoutineName[]={"MFProjectSpline"};
  double s;
  int    rc;
  struct MFFitpackSplineData *data;
  double *u0;
  double *u ;
  double *phi;
  double f,df;
  int i;
  int itimes;
  int npts;
  static double *du=NULL;
  MFNVector PhiC;
  int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("In %s\n",RoutineName);fflush(stdout);}
#endif

  data=(struct MFFitpackSplineData*)d;

  u0=MFSplineNV_CStar(vu0,e);
  u =MFSplineNV_CStar(vu ,e);
  PhiC=MFMColumn(mPhi,0,e);
  phi=MFSplineNV_CStar(PhiC,e);
  du=(double*)realloc((void*)du,(data->n)*sizeof(double));

/*  solve f(s)=(u(s)-u0).phi=0; */

  s=MFSplineNVectorGetS(vu0,e);

  spline(s,u,data);
  f=0.;for(i=0;i<n;i++) f+=(u[i]-u0[i])*phi[i];

  itimes=0;
  while(fabs(f)>1.e-6&&itimes<10)
   {
#ifdef MFALLOWVERBOSE
     if(verbose){printf("%d, s=%f, f=%lf, df=%lf\n",itimes,s,f,df);fflush(stdout);}
#endif

      spline(s,u,data);
     dspline(s,du,data);
     f=0.;for(i=0;i<n;i++) f+=(u[i]-u0[i])*phi[i];
    df=0.;for(i=0;i<n;i++)df+=du[i]*phi[i];
    s=s-f/df;

     spline(s,u,data);
    dspline(s,du,data);
    itimes++;
   }
#ifdef MFALLOWVERBOSE
  if(verbose){printf("%d, s=%f, f=%lf, df=%lf\n",itimes,s,f,df);fflush(stdout);}
#endif
  rc=1;
  if(itimes>9)rc=0;

  *index=0;

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("u0=");MFPrintNVector(stdout,vu0,e);printf("\n");fflush(stdout);
    printf("u =");MFPrintNVector(stdout,vu,e);printf(" s=%10.7f\n",s);fflush(stdout);
    printf("done %s, rc=%d\n",RoutineName,rc);fflush(stdout);
   }
#endif

  MFSplineNVectorSetS(vu,s,e);
  MFFreeNVector(PhiC,e);
  return rc;
 }

int MFTangentSpline(int n,int k,MFNVector vu,MFNKMatrix mPhi,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFTangentSpline"};
  double s;
  double t;
  double *u;
  double *Phi;
  double *u0;
  int    index;
  MFNVector vu0;
  int i;
  struct MFFitpackSplineData *data;
  int verbose=0;
  MFNVector PhiC;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("In %s\n",RoutineName);fflush(stdout);}
#endif

  data=(struct MFFitpackSplineData*)d;


  PhiC=MFMColumn(mPhi,0,e);
  Phi=MFSplineNV_CStar(PhiC,e);
  u =MFSplineNV_CStar(vu ,e);

  vu0=MFCloneNVector(vu,e);
  u0=MFSplineNV_CStar(vu0,e);
  t =MFSplineNVectorGetS(vu0,e);

  dspline(t,Phi,data);

  t=0.;
  for(i=0;i<n;i++)t+=Phi[i]*Phi[i];
  if(fabs(t)>1.e-15)t=1./sqrt(t);

  for(i=0;i<n;i++)Phi[i]=Phi[i]*t;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("done %s\n",RoutineName);fflush(stdout);}
#endif

  MFFreeNVector(vu0,e);
  MFFreeNVector(PhiC,e);
  return 1;
 }

double MFScaleSpline(int n,int k,MFNVector vu,MFNKMatrix mPhi,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFScaleSpline"};

  double s;
  double t;
  int i;
   struct MFFitpackSplineData *data;
  int npts;
  int verbose=0;
  double *u;
  double *u0;
  int    index;
  MFNVector vu0;
  double r;
  static double *ddF=NULL;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s\n",RoutineName);fflush(stdout);}
#endif

  data=(struct MFFitpackSplineData*)d;

  ddF   =(double*)realloc(ddF ,(data->n)*sizeof(double));

  u =MFSplineNV_CStar(vu,e);

  vu0=MFCloneNVector(vu,e);
   u0=MFSplineNV_CStar(vu0,e);
   t =MFSplineNVectorGetS(vu0,e);
  ddspline(t,ddF,data);

  t=0.;
  for(i=0;i<n;i++)t+=ddF[i]*ddF[i];
  t=sqrt(t);
 
  r=sqrt(2*MFEpsilon/t);

#ifdef MFALLOWVERBOSE
 if(verbose)
  {
   printf("    |u''|=%lf\n",t);fflush(stdout);
   printf("     r   =%lf\n",r);fflush(stdout);
   printf("done %s\n",RoutineName);fflush(stdout);
  }
#endif

  MFFreeNVector(vu0,e);

  return r;
 }

int MFStopSpline(MFImplicitMF thisMF,MFNVector u0,MFNKMatrix Phi0,MFNVector u1,MFNKMatrix Phi1,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFStopSpline"};
  struct MFFitpackSplineData *data;
  int npts;
  double s0,s1;
  int verbose=0;

  data=(struct MFFitpackSplineData*)d;

  if(data->periodic)return 0;

  s0=MFFindS(MFSplineNV_CStar(u0,e),d,e);
  s1=MFFindS(MFSplineNV_CStar(u1,e),d,e);

  if(verbose)printf("\n\n***%s\n",RoutineName);fflush(stdout);

  return 0;
 }

int MFSplineProjectToSave(MFNVector u, double *x, void *d, MFErrorHandler e)
 {
  int i;
  struct MFFitpackSplineData *data;

  data=(struct MFFitpackSplineData*)d;

  if(x==NULL)return data->n;
  for(i=0;i<data->n;i++)x[i]=MFNV_C(u,i,e);

  return 0;
 }

/*! @} */

/*! @} */

static double MFFindS(double *u, void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFFindS"};
  int    i,j;
  double  s;
  double dist;
  double dmin;
  int    jmin;
  double  tmp;
  struct MFFitpackSplineData *data;
  int N,M;
  double sL,sR;
  double dL,dR;
  double ss;
  int verbose=0;
  static double *us;

  if(verbose){printf("In %s\n",RoutineName);fflush(stdout);}
  if(verbose){printf("     find s for u=(%le,%le,%le,%le)\n",u[0],u[1],u[2],u[3]);fflush(stdout);}

  us=(double*)realloc((void*)us,N*sizeof(double));

  data=(struct MFFitpackSplineData*)d;

  N=data->n;
  M=data->npts;

  jmin=0;
  j=0;dist=0;for(i=0;i<N;i++)dist+=(u[i]-data->x[i][j])*(u[i]-data->x[i][j]);dist=sqrt(dist)/N;
  if(verbose){printf("  %2d   s=%10.7f,  d=%10.7f u=(%le,%le,%le,%le)\n",j,data->s[j],dist,data->x[0][j],data->x[1][j],data->x[2][j],data->x[3][j]);}
  dmin=dist;
  for(j=1;j<M;j++)
   {
    dist=0;for(i=0;i<N;i++)dist+=(u[i]-data->x[i][j])*(u[i]-data->x[i][j]);dist=sqrt(dist)/N;
    if(0&&verbose){printf("  %2d   s=%10.7f,  d=%10.7f u=(%le,%le,%le,%le)\n",j,data->s[j],dist,data->x[0][j],data->x[1][j],data->x[2][j],data->x[3][j]);}
    if(dist<dmin){dmin=dist;jmin=j;}
   }
  if(verbose){printf("     dmin=%10.7f jmin=%d\n",dmin,jmin);fflush(stdout);}

  if(jmin==0)
   {
    sL=data->s[M-1]-1.;
    j=M-1;dist=0;for(i=0;i<N;i++)dist+=(u[i]-data->x[i][j])*(u[i]-data->x[i][j]);dist=sqrt(dist)/N;
    dL=dist;
    sR=data->s[1];
    j=1;dist=0;for(i=0;i<N;i++)dist+=(u[i]-data->x[i][j])*(u[i]-data->x[i][j]);dist=sqrt(dist)/N;
    dR=dist;
    if(verbose){printf("     use interval between knots %d and %d\n",M-1,jmin+1);fflush(stdout);}
   }else if(jmin==M-1)
   {
    sL=data->s[M-2];
    j=jmin-1;dist=0;for(i=0;i<N;i++)dist+=(u[i]-data->x[i][j])*(u[i]-data->x[i][j]);dist=sqrt(dist)/N;
    dL=dist;
    sR=data->s[0]+1.;
    j=0;dist=0;for(i=0;i<N;i++)dist+=(u[i]-data->x[i][j])*(u[i]-data->x[i][j]);dist=sqrt(dist)/N;
    dR=dist;
    if(verbose){printf("     use interval between knots %d and %d\n",jmin-1,0);fflush(stdout);}
   }else{
    sL=data->s[jmin-1];
    j=jmin-1;dist=0;for(i=0;i<N;i++)dist+=(u[i]-data->x[i][j])*(u[i]-data->x[i][j]);dist=sqrt(dist)/N;
    dL=dist;
    sR=data->s[jmin+1];
    j=jmin+1;dist=0;for(i=0;i<N;i++)dist+=(u[i]-data->x[i][j])*(u[i]-data->x[i][j]);dist=sqrt(dist)/N;
    dR=dist;
    if(verbose){printf("     use interval between knots %d and %d\n",jmin-1,jmin+1);fflush(stdout);}
   }

  while(fabs(sR-sL)>1.e-7)
   {
    s=(sL+sR)/2;
    ss=s;if(ss<0.)ss+=1.;if(ss>1.)ss-=1.;
    if(verbose){printf("\nss=%f\n",ss);fflush(stdout);}
    spline(ss,us,data);
    dist=0;for(i=0;i<N;i++)dist+=(u[i]-us[i])*(u[i]-us[i]);dist=sqrt(dist)/N;

    if(verbose){printf("     sL=%10.7lf, s=%10.7lf, sR=%10.7lf\n",sL,s,sR);fflush(stdout);}
    if(verbose){printf("     dL=%10.7lf, d=%10.7lf, dR=%10.7lf\n\n",dL,dist,dR);fflush(stdout);}

    if(dist<dL)
     {
      sL=s;
      dL=dist;
     }else{
      sR=s;
      dR=dist;
     }
   }
  ss=s;if(ss<0.)ss+=1.;if(ss>1.)ss-=1.;
  if(1||verbose){printf("     FindS s=%lf, error=%14.7le\n",s,dist);fflush(stdout);}

  return ss;
 }

static void spline(double s,double *u, struct MFFitpackSplineData *data)
 {
  static char RoutineName[]={"spline"};
  int i,j;
  int verbose=0;
  double z=0.;

  int ierr;
  int m;

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("In %s\n",RoutineName);
    printf("  n   =%d\n",data->n);
    printf("  npts=%d\n",data->npts);
    fflush(stdout);
   }
#endif

  if(data->periodic)
   {
    for(i=0;i<data->n;i++){CALLFCURVP2(&s,&(data->npts),data->s,data->x[i],&(data->p),data->yp[i],&(data->sigma),&z);u[i]=z;}
   }else{
    for(i=0;i<data->n;i++){CALLFCURV2(&s,&(data->npts),data->s,data->x[i],data->yp[i],&(data->sigma),&z);u[i]=z;}
   }

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("  u(%14.7le):  (",s);
    for(i=0;i<data->n;i++)
     {
      if(i>0)printf(",");
      printf("%14.7le",u[i]);
     }
    printf(")\n");fflush(stdout);
  
    printf("done %s\n",RoutineName);fflush(stdout);
   }
#endif

  return;
 }

static void dspline(double s, double *du, struct MFFitpackSplineData *data)
 {
  static char RoutineName[]={"dspline"};
  int i;
  int verbose=0;
  double z;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("In %s\n",RoutineName);fflush(stdout);}
#endif

  for(i=0;i<data->n;i++)
   {
    CALLFDCURVP2(&s,&(data->npts),data->s,data->x[i],&(data->p),data->yp[i],&(data->sigma),&z);
    du[i]=z;
   }

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf(" du(%14.7le):  (",s);
    for(i=0;i<data->n;i++)
     {
      if(i>0)printf(",");
      printf("%14.7le",du[i]);
     }
    printf(")\n");fflush(stdout);
  
    printf("done %s\n",RoutineName);fflush(stdout);
  }
#endif

  return;
 }

static void ddspline(double s, double *ddu, struct MFFitpackSplineData *data)
 {
  static char RoutineName[]={"ddspline"};
  double t;
  int i;
  static double *Fplus =NULL;
  static double *F     =NULL;
  static double *Fminus=NULL;
  double eps=1.e-2;
  int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("In %s\n",RoutineName);fflush(stdout);}
#endif

  Fplus =(double*)realloc(Fplus ,(data->n)*sizeof(double));
  F     =(double*)realloc(Fplus ,(data->n)*sizeof(double));
  Fminus=(double*)realloc(Fminus,(data->n)*sizeof(double));

  t=s+eps;spline(t,Fplus ,data);
  t=s    ;spline(t,F     ,data);
  t=s-eps;spline(t,Fminus,data);

  if(verbose){printf("  difference du\n");fflush(stdout);}
  for(i=0;i<data->n;i++)
   {
    ddu[i]=.25*(Fplus[i]-2*F[i]+Fminus[i])/eps/eps;
    if(verbose){printf("  %d ddu=%14.7le f+=%14.7le, f0=%14.7le f-=%14.7le\n",i,ddu[i],Fplus[i],F[i],Fminus[i]);fflush(stdout);}
   }

#ifdef MFALLOWVERBOSE
  if(verbose){printf("done %s\n",RoutineName);fflush(stdout);}
#endif

  return;
 }

struct MFSplineNVectorData
 {
  int n;
  double *u;
  double  s;
 };

static void MFFreeSplineNVectorData(void*,MFErrorHandler);
static int MFSplineNVGetNC(void*,MFErrorHandler);
static double MFSplineNVGetC(int,void*,MFErrorHandler);
static void MFSplineNVSetC(int,double,void*,MFErrorHandler);
static void MFSplineNVDiff(void*,void*,void*,MFErrorHandler);
static void MFSplineNVAdd(void*,void*,void*,MFErrorHandler);
static MFNVector MFCloneSplineNVector(void*,MFErrorHandler);


void MFPrintSplineNVector(FILE*,void*,MFErrorHandler);

void MFFreeSplineNVectorData(void *data, MFErrorHandler e)
 {
  static char RoutineName[]={"MFFreeSplineNVectorData"};
  int i;
  struct MFSplineNVectorData *this;
  int verbose=0;

  this=(struct MFSplineNVectorData*)data;

  if(this->u!=NULL)
   {
    free(this->u);
    this->u=NULL;
   }

  free(this);

  return;
 }

int MFSplineNVGetNC(void *data, MFErrorHandler e)
 {
  static char RoutineName[]={"MFSplineNVGetNC"};
  struct MFSplineNVectorData *this;

  this=(struct MFSplineNVectorData*)data;

#ifdef MFNOCONFIDENCE
  if(this==NULL)
   {
    sprintf(MFSplineErrorMsg,"Pointer to Vector Data (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFSplineErrorMsg,__LINE__,__FILE__);
    return -1;
   }
#endif

  return this->n;
 }

double MFSplineNVGetC(int i,void *data, MFErrorHandler e)
 {
  static char RoutineName[]={"MFSplineNVGetC"};
  struct MFSplineNVectorData *this;
  int j;
  double result;

  this=(struct MFSplineNVectorData*)data;

#ifdef MFNOCONFIDENCE
  if(this==NULL)
   {
    sprintf(MFSplineErrorMsg,"Pointer to Vector Data (argument 2) is NULL");
    MFSetError(e,12,RoutineName,MFSplineErrorMsg,__LINE__,__FILE__);
    return DBL_QNAN;
   }

  if(i<0|| !(i<this->n))
   {
    sprintf(MFSplineErrorMsg,"Coordinate %d (argument 1) is illegal. Must be in 0 to %d",i,this->n-1);
    MFSetError(e,8,RoutineName,MFSplineErrorMsg,__LINE__,__FILE__);
    return DBL_QNAN;
   }
#endif

  return (this->u)[i];
 }

void MFSplineNVSetC(int i,double vl,void *data, MFErrorHandler e)
 {
  static char RoutineName[]={"MFSplineNVSetC"};
  int j;
  struct MFSplineNVectorData *this;

  this=(struct MFSplineNVectorData*)data;

#ifdef MFNOCONFIDENCE
  if(this==NULL)
   {
    sprintf(MFSplineErrorMsg,"Pointer to Vector Data (argument 3) is NULL");
    MFSetError(e,12,RoutineName,MFSplineErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(i<0|| !(i<this->n))
   {
    sprintf(MFSplineErrorMsg,"Coordinate %d (argument 1) is illegal. Must be in 0 to %d",i,this->n-1);
    MFSetError(e,8,RoutineName,MFSplineErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

   this->u[i]=vl;

  return;
 }

void MFSplineNVDiff(void *adata,void *bdata, void *cdata, MFErrorHandler e)
 {
  static char RoutineName[]={"MFSplineNVDiff"};
  int i,ni;
  int j,nj;
  struct MFSplineNVectorData *a;
  struct MFSplineNVectorData *b;
  struct MFSplineNVectorData *c;

  a=(struct MFSplineNVectorData*)adata;

#ifndef MFNOSAFETYNET
  if(a==NULL)
   {
    sprintf(MFSplineErrorMsg,"Pointer to Vector Data for a (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFSplineErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  b=(struct MFSplineNVectorData*)bdata;

#ifndef MFNOSAFETYNET
  if(b==NULL)
   {
    sprintf(MFSplineErrorMsg,"Pointer to Vector Data for b (argument 2) is NULL");
    MFSetError(e,12,RoutineName,MFSplineErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  c=(struct MFSplineNVectorData*)cdata;

#ifndef MFNOSAFETYNET
  if(c==NULL)
   {
    sprintf(MFSplineErrorMsg,"Pointer to Vector Data for c (argument 3) is NULL");
    MFSetError(e,12,RoutineName,MFSplineErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

#ifdef MFNOCONFIDENCE
  if(a->n!=b->n || a->n!=c->n || b->n!=c->n)
   {
    sprintf(MFSplineErrorMsg,"Vectors must all be the same length a=%d, b=%d, c=%d",a->n,b->n,c->n);
    MFSetError(e,4,RoutineName,MFSplineErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  for(i=0;i<a->n;i++)
    {
     c->u[i]=a->u[i]-b->u[i];
    }
  c->s=0.;

  return;
 }

void MFSplineNVAdd(void *adata,void *bdata,void *cdata, MFErrorHandler e)
 {
  static char RoutineName[]={"MFSplineNVAdd"};
  int i,ni;
  int j,nj;
  struct MFSplineNVectorData *a;
  struct MFSplineNVectorData *b;
  struct MFSplineNVectorData *c;

  a=(struct MFSplineNVectorData*)adata;

#ifndef MFSAFETYNET
  if(a==NULL)
   {
    sprintf(MFSplineErrorMsg,"Pointer to Vector Data for a (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFSplineErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  b=(struct MFSplineNVectorData*)bdata;

#ifdef MFSAFETYNET
  if(b==NULL)
   {
    sprintf(MFSplineErrorMsg,"Pointer to Vector Data for b (argument 2) is NULL");
    MFSetError(e,12,RoutineName,MFSplineErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  c=(struct MFSplineNVectorData*)cdata;

#ifdef MFSAFETYNET
  if(c==NULL)
   {
    sprintf(MFSplineErrorMsg,"Pointer to Vector Data for c (argument 3) is NULL");
    MFSetError(e,12,RoutineName,MFSplineErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

#ifdef MFNOCONFIDENCE
  if(a->n!=b->n || a->n!=c->n || b->n!=c->n)
   {
    sprintf(MFSplineErrorMsg,"Vectors must all be the same length a=%d, b=%d, c=%d",a->n,b->n,c->n);
    MFSetError(e,4,RoutineName,MFSplineErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  for(i=0;i<ni;i++)
     c->u[i]=a->u[i]+b->u[i];
  c->s=0.;

  return;
 }

MFNVector MFCloneSplineNVector(void *data, MFErrorHandler e)
 {
  static char RoutineName[]={"MFCloneSplineNVector"};
  struct MFSplineNVectorData *u;
  struct MFSplineNVectorData *v;
  MFNVector clone;
  int i,j;

  u=(struct MFSplineNVectorData*)data;

  clone=MFCreateSplineNVector(u->n,e);

  v=(struct MFSplineNVectorData*)MFNVectorGetData(clone,e);

  for(i=0;i<u->n;i++)v->u[i]=u->u[i];
  v->s=u->s;

  return clone;
 }

void MFPrintSplineNVector(FILE *fid,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFPrintSplineNVector"};
  int i;
  struct MFSplineNVectorData *data;

  data=(struct MFSplineNVectorData*)d;

#ifdef MFNOCONFIDENCE
  if(fid==NULL)
   {
    sprintf(MFSplineErrorMsg,"fid (argument 1) is NULL.");
    MFSetError(e,12,RoutineName,MFSplineErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(data==NULL)
   {
    sprintf(MFSplineErrorMsg,"data (argument 2) is NULL.");
    MFSetError(e,12,RoutineName,MFSplineErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  fprintf(fid,"(");
  for(i=0;i<data->n;i++)
   {
    if(i>0)fprintf(fid,",");
    fprintf(fid,"%le",data->u[i]);
   }
  fprintf(fid,")");

  return;
 }

MFNVector MFCreateSplineNVector(int n,MFErrorHandler e)
 {
  static char RoutineName[]={"MFCreateSplineNVector"};
  MFNVector this;
  int i,j;
  struct MFSplineNVectorData *data;

  int verbose=0;

  data=malloc(sizeof(struct MFSplineNVectorData));

#ifndef MFNOSAFETYNET
  if(data==NULL)
   {
    sprintf(MFSplineErrorMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct MFSplineNVectorData));
    MFSetError(e,12,RoutineName,MFSplineErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  data->n=n;
  data->u=(double*)malloc(n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(data->u==NULL)
   {
    sprintf(MFSplineErrorMsg,"Out of memory, trying to allocate %d bytes",(data->n)*sizeof(double));
    MFSetError(e,12,RoutineName,MFSplineErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  for(i=0;i<n;i++)data->u[i]=0.;
  data->s=0.;

  this=MFCreateNVectorBaseClass("Spline",e);
  MFNVectorSetData(this,data,e);
  MFNVectorSetFreeData(this,MFFreeSplineNVectorData,e);

  MFNVectorSetGetNC(this,MFSplineNVGetNC,e);
  MFNVectorSetGetC(this,MFSplineNVGetC,e);
  MFNVectorSetSetC(this,MFSplineNVSetC,e);
  MFNVectorSetDiff(this,MFSplineNVDiff,e);
  MFNVectorSetAdd(this,MFSplineNVAdd,e);
  MFNVectorSetClone(this,MFCloneSplineNVector,e);
  MFNVectorSetPrint(this,MFPrintSplineNVector,e);

  return this;
 }
void MFSplineNVectorSetS(MFNVector this,double s,MFErrorHandler e)
 {
  struct MFSplineNVectorData *data;
  data=(struct MFSplineNVectorData*)MFNVectorGetData(this,e);
  data->s=s;
  return;
 }

double MFSplineNVectorGetS(MFNVector this,MFErrorHandler e)
 {
  struct MFSplineNVectorData *data;
  data=(struct MFSplineNVectorData*)MFNVectorGetData(this,e);

  return data->s;
 }


MFNVector MFSplineNVectorFactory(MFImplicitMF M,MFErrorHandler e)
 {
  return MFCreateSplineNVector(MFIMF_N(M,e),e);
 }

static double *MFSplineNV_CStar(MFNVector u,MFErrorHandler e)
 {
  struct MFSplineNVectorData *data;

  data=(struct MFSplineNVectorData*)MFNVectorGetData(u,e);

  return data->u;
 }

static double MFSplineNSpaceDistance(MFNSpace,MFNVector,MFNVector,void*,MFErrorHandler);
static void MFSplineNSpaceDirection(MFNSpace,MFNVector,MFNVector,MFNVector,void*,MFErrorHandler);
static void MFSplineNSpaceAdd(MFNSpace,MFNVector,MFNVector,MFNVector,void*,MFErrorHandler);
static void MFSplineNSpaceScale(MFNSpace,double,MFNVector,MFNVector,void*,MFErrorHandler);
static double MFSplineNSpaceInner(MFNSpace,MFNVector,MFNVector,void*,MFErrorHandler);

MFNSpace MFCreateSplineNSpace(MFErrorHandler e)
 {
  static char RoutineName[]={"MFCreateSplineNSpace"};
  MFNSpace this;
  struct MFSplineNSpaceData *data;

  this=MFCreateNSpaceBaseClass("SplineNSpace",e);
  MFNSpaceSetDistance(this,MFSplineNSpaceDistance,e);
  MFNSpaceSetInnerProduct(this,MFSplineNSpaceInner,e);
  MFNSpaceSetDirection(this,MFSplineNSpaceDirection,e);
  MFNSpaceSetAdd(this,MFSplineNSpaceAdd,e);
  MFNSpaceSetScale(this,MFSplineNSpaceScale,e);

  return this;
 }

double MFSplineNSpaceDistance(MFNSpace this,MFNVector v0,MFNVector v1,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFNSpaceDistance"};
  double ss;
  MFNVector diff;

  int i;
  diff=MFCloneNVector(v0,e);
  MFSplineNSpaceDirection(this,v0,v1,diff,d,e);
  ss=sqrt(MFSplineNSpaceInner(this,diff,diff,d,e));
  MFFreeNVector(diff,e);
  return ss;
 }

void MFSplineNSpaceDirection(MFNSpace this,MFNVector v0,MFNVector v1,MFNVector diff,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFNSpaceDirection"};
  int i,j;
  double *u0,*u1,*dl;

  u0=MFSplineNV_CStar(v0,e);
  u1=MFSplineNV_CStar(v1,e);
  dl=MFSplineNV_CStar(diff,e);

  for(i=0;i<MFNV_NC(v0,e);i++)dl[i]=u1[i]-u0[i];

  return;
 }

void MFSplineNSpaceAdd(MFNSpace this,MFNVector v0,MFNVector v1,MFNVector sum,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFNSplineSpaceAdd"};
  int i;
  double *u0;
  double *u1;
  double *usum;

  u0=MFSplineNV_CStar(v0,e);
  u1=MFSplineNV_CStar(v1,e);
  usum=MFSplineNV_CStar(sum,e);

  for(i=0;i<MFNV_NC(v0,e);i++)usum[i]=u0[i]+u1[i];

  return;
 }

double MFSplineNSpaceInner(MFNSpace this,MFNVector v0,MFNVector v1,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFSplineNSpaceInner"};
  int i;
  double *u0;
  double *u1;
  double prod;

  u0=MFSplineNV_CStar(v0,e);
  u1=MFSplineNV_CStar(v1,e);

  prod=0.;for(i=0;i<MFNV_NC(v0,e);i++)prod+=u0[i]*u1[i];

  return prod;
 }

void MFSplineNSpaceScale(MFNSpace this,double s,MFNVector v,MFNVector prod,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFNSplineSpaceScale"};
  int i;
  double *u;
  double *p;

  u=MFSplineNV_CStar(v,e);
  p=MFSplineNV_CStar(prod,e);

  for(i=0;i<MFNV_NC(v,e);i++)p[i]=s*u[i];

  return;
 }

#ifdef __cplusplus
}
#endif
