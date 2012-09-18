/*
 *  PROGRAM NAME:  multifario
 *
 *  (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
 *  CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
 *
 *  Please refer to the LICENSE file in the top directory
 *
 *      author: Mike Henderson mhender@watson.ibm.com
 *      date:   December 9, 2002   modified InvMF
 */

static char *id="@(#) $Id: IMFSphereOnExpansion.c,v 1.8 2011/07/21 17:42:46 mhender Exp $";

static char IMFSphereOnExpansionErrorMsg[256];

#include <MFAtlas.h>
#include <IMFExpansion.h>
#include <MFErrorHandler.h>
#include <IMFFlow.h>
#include <IMFNSpace.h>
#include <MFFortran.h>
#include <MFPrint.h>
#include <IMF.h>
#include <math.h>

#ifdef __cplusplus
 extern "C" {
#endif

#ifndef MAX
#define MAX(X, Y) ((X) < (Y) ? (Y) : (X))
#endif
#define round(x) ((x)>=0?(long)((x)+0.5):(long)((x)-0.5))

struct IMFSphereOnEData
    {
     int n;
     int k;
     IMFExpansion E;
     IMFFlow F;
     MFKVector p0;
     double eps;
     double R;
    };

void MFFreeIMFSphereOnEData(void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFFreeIMFSphereOnEData"};

  struct IMFSphereOnEData *data;
  data=(struct IMFSphereOnEData*)d;
  IMFFreeExpansion(data->E,e);
  IMFFreeFlow(data->F,e);
  free(data);
 }

static int MFProjectIMFSphereOnE(int,int,MFNVector,MFNKMatrix,MFNVector,void*,int*,MFErrorHandler);
static int MFTangentIMFSphereOnE(int,int,MFNVector,MFNKMatrix,void*,MFErrorHandler);
static int MFIMFSphereOnEProjectForBB(MFNVector,double*,void*,MFErrorHandler);
static double MFScaleIMFSphereOnE(int,int,MFNVector,MFNKMatrix,void*,MFErrorHandler);
static MFNVector MFSOENVectorFactory(MFImplicitMF,MFErrorHandler);
static MFNKMatrix MFSOENKMatrixFactory(MFImplicitMF,MFErrorHandler);

/*! \fn MFImplicitMF IMFCreateSphereOnExpansion(IMFExpansion E, IMFFlow F, MFKVector p0,double eps, double R, double r, MFErrorHandler e);
 *  \brief Creates a manifold which is the intersection of a sphere and the surface defined by an expansion.
 *
 *  \param E The expansion.
 *  \param F The flow.
 *  \param p0 The parameters for the flow.
 *  \param eps The tolerance allowed for the error on the conditions that the manifold lie on the expansion.
 *  \param R The raidus of the sphere.
 *  \param r The default scale associated with the new manifold (MFIMFSetR).
 *  \returns A new MFImplicitMF.
 */
MFImplicitMF IMFCreateSphereOnExpansion(IMFExpansion E, IMFFlow F, MFKVector p0,double eps, double R, double r, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFCreateSphereOnExpansion"};

/*   Manifold c with
          c-E(s)   =0
         |c-E(0)|^2=R^2               */

  MFImplicitMF thisSphere;
  int n,k;
  struct IMFSphereOnEData *data;
  MFNSpace space;

  n=IMFExpansionN(E,e);
  k=IMFExpansionK(E,e);
  thisSphere=MFIMFCreateBaseClass(n+k,k-1,"InvariantManifoldInitialCurve",e);
  space=InvMFCreateNSpace(n,e);
  MFIMFSetSpace(thisSphere,space,e);
  MFFreeNSpace(space,e);

  data=(struct IMFSphereOnEData*)malloc(sizeof(struct IMFSphereOnEData));

#ifndef MFNOSAFETYNET
  if(data==NULL)
   {
    sprintf(IMFSphereOnExpansionErrorMsg,"Out of memory trying to allocate %d bytes\n",sizeof(struct IMFSphereOnEData));
    MFSetError(e,12,RoutineName,IMFSphereOnExpansionErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  data->E=E;IMFRefExpansion(E,e);
  data->F=F;IMFRefFlow(F,e);
  data->p0=p0;if(p0!=NULL)MFRefKVector(p0,e);
  data->eps=eps;
  data->R=R;
  data->n=n;
  data->k=k;
  MFIMFSetData(thisSphere,(void*)data,e);

  MFIMFSetFreeData(thisSphere,MFFreeIMFSphereOnEData,e);
  MFIMFSetProject(thisSphere,MFProjectIMFSphereOnE,e);
  MFIMFSetTangent(thisSphere,MFTangentIMFSphereOnE,e);
  MFIMFSetScale(thisSphere,MFScaleIMFSphereOnE,e);
  MFIMFSetR(thisSphere,r,e);

  MFIMFSetProjectForBB(thisSphere,MFIMFSphereOnEProjectForBB,e);
  MFIMFSetProjectForDraw(thisSphere,MFIMFSphereOnEProjectForBB,e);

  MFIMFSetVectorFactory(thisSphere,MFSOENVectorFactory,e);
  MFIMFSetMatrixFactory(thisSphere,MFSOENKMatrixFactory,e);

  return thisSphere;
 }

int MFProjectIMFSphereOnE(int n,int k,MFNVector vu0,MFNKMatrix mPhi,MFNVector vu,void *d,int *index, MFErrorHandler e)
 {
  static char RoutineName[]={"MFProjectIMFSphereOnE"};
  double err;
  double *A=NULL; /* [n*n];*/
  double *b=NULL; /* [n];  */
  int *ipvt=NULL;    /* [n];  */
  int one=1;
  int info=0;
  int itimes;
  int i,j,l,m;
  char trans='N';
  struct IMFSphereOnEData *data;
  double f;
  double *df=NULL;
  double t,dnorm,unorm;
  double *u0,*Phi,*u,*EPhi;

  double R;
  MFNVector vM0;
  MFNVector vM;
  MFNVector vdM;
  double *M0;
  double *M;
  double *dM;
  double *s;
  MFKVector vs;

  int verbose=0;

  data=(struct IMFSphereOnEData*)d;

  m=IMFExpansionN(data->E,e);
  u0=MFNV_CStar(vu0,e);
  Phi=MFNKM_CStar(mPhi,e);
  u=MFNV_CStar(vu,e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("In %s. n=%d, m=%d, k=%d\n",RoutineName,n,m,k);fflush(stdout);}
#endif

  A=(double*)malloc(n*n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(A==NULL)
   {
    sprintf(IMFSphereOnExpansionErrorMsg,"Out of memory trying to allocate %d bytes\n",n*n*sizeof(double));
    MFSetError(e,12,RoutineName,IMFSphereOnExpansionErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0;
   }
#endif

  b=(double*)malloc(n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(b==NULL)
   {
    sprintf(IMFSphereOnExpansionErrorMsg,"Out of memory trying to allocate %d bytes\n",n*sizeof(double));
    MFSetError(e,12,RoutineName,IMFSphereOnExpansionErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0;
   }
#endif

  ipvt=(int*)malloc(n*sizeof(int));

#ifndef MFNOSAFETYNET
  if(ipvt==NULL)
   {
    sprintf(IMFSphereOnExpansionErrorMsg,"Out of memory trying to allocate %d bytes\n",n*sizeof(int));
    MFSetError(e,12,RoutineName,IMFSphereOnExpansionErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0;
   }
#endif

  df=(double*)malloc(n*(n-k)*sizeof(double));

#ifndef MFNOSAFETYNET
  if(df==NULL)
   {
    sprintf(IMFSphereOnExpansionErrorMsg,"Out of memory trying to allocate %d bytes\n",n*(n-k)*sizeof(double));
    MFSetError(e,12,RoutineName,IMFSphereOnExpansionErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0;
   }
#endif

  for(i=0;i<n;i++)u[i]=u0[i];

  vs=MFCreateKVector(k+1,e);
  s=MFKV_CStar(vs,e);

  vM0=MFCreateNVector(m,e);
  M0=MFNV_CStar(vM0,e);

  vdM=MFCreateNVector(m,e);
  dM=MFNV_CStar(vdM,e);

  R=data->R;
  for(i=0;i<k+1;i++)s[i]=0.;
  IMFEvaluateExpansion(data->E,vs,vM0,e);

  for(i=0;i<k+1;i++)s[i]=u[m+i];

  err=1.;
  itimes=0;
  while(err>1.e-10)
   {
    for(i=0;i<k+1;i++)u[m+i]=s[i];
    IMFEvaluateExpansion(data->E,vs,vu,e);

    f=-R*R;for(i=0;i<m;i++)f+=(u[i]-M0[i])*(u[i]-M0[i]);

    b[0]=-f;
    for(j=1;j<k+1;j++)
      {b[j]=0.;for(i=0;i<k+1;i++)b[j]-=Phi[i+m+n*(j-1)]*(u[i+m]-u0[i+m]);}

    for(j=0;j<k+1;j++)
     {
      IMFEvaluateExpansionDerivative(data->E,vs,j,vdM,e);
      A[0+(k+1)*j]=0.;for(i=0;i<m;i++)A[0+(k+1)*j]+=2*(u[i]-M0[i])*dM[i];
  
      for(l=1;l<k+1;l++)
        A[l+(k+1)*j]=Phi[j+m+n*(l-1)];
     }

#ifdef MFALLOWVERBOSE
    if(verbose)
     {
      printf("System for Newton Correction\n");
      IMFPrintFull(k+1,A,b,e);
     }
#endif

    err=0.;
    for(i=0;i<k+1;i++)err+=b[i]*b[i];
    err=sqrt(err);

#ifdef MFALLOWVERBOSE
    if(verbose)printf("  %d  %le\n",itimes,err);fflush(stdout);
#endif

    MFSolveFull(k+1,A,b,e);

    t=1.;
    dnorm=fabs(b[0]);
    for(i=1;i<k+1;i++)
      if(fabs(b[i])>dnorm)dnorm=fabs(b[i]);
    unorm=fabs(s[0]);
    for(i=1;i<k+1;i++)
      if(fabs(s[i])>unorm)unorm=fabs(s[i]);
    if(dnorm>.1*(unorm+1))t=.1*(unorm+1)/dnorm; /* Damping */

    for(i=0;i<k+1;i++){s[i]=s[i]+t*b[i];u[i+m]=s[i];}

    itimes++;
    if(itimes>10)
     {
      free(A);
      free(b);
      free(ipvt);
      free(df);
      MFFreeKVector(vs,e);
      MFFreeNVector(vM0,e);
      MFFreeNVector(vdM,e);

#ifdef MFALLOWVERBOSE
      if(verbose){printf("Done %s -- too many iterations\n",RoutineName );fflush(stdout);}
#endif

      return 0;
     }
   }
  *index=0;
  for(i=0;i<k+1;i++)u[i+m]=s[i];

#ifdef MFALLOWVERBOSE
  if(verbose){printf("Result is ");MFPrintNVector(stdout,vu,e);printf("\n");fflush(stdout);}
#endif

  free(A);
  free(b);
  free(ipvt);
  free(df);
  MFFreeKVector(vs,e);
  MFFreeNVector(vM0,e);
  MFFreeNVector(vdM,e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("Done %s\n",RoutineName);fflush(stdout);}
#endif

  return 1;
 }

int MFTangentIMFSphereOnE(int n,int k,MFNVector vu,MFNKMatrix mPhi,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFTangentIMFSphereOnE"};
  double *A=NULL; /* [n*n];*/
  double *b=NULL; /* [n];  */
  int *ipvt=NULL;    /* [n];  */
  int one=1;
  int info=0;
  int itimes;
  int i,j,l,m;
  char trans='N';
  struct IMFSphereOnEData *data;
  double t;
  double *u0,*Phi,*u;

  char jobvl='N';                   /* No left  eigenvectors */
  char jobvr='A';                   /*    right eigenvectors */
  int lda;                        /* Leading dimension of A */
  double *sv;    /* n */
  double *U;    /* n*n */
  int ldU;
  double *vr;    /* n*k */
  double *V;    /* n*n */
  int ldV;
  double *work;                     /* Work array */
  int lwork;                     /* Length of work array */
  int msv;

  double R;
  MFNVector vM0;
  MFNVector vM;
  MFNVector vdM;
  MFKVector vs;
  double *M0;
  double *M;
  double *dM;
  double *s;

  int verbose=0;

  data=(struct IMFSphereOnEData*)d;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("In %s\n",RoutineName);fflush(stdout);}
#endif

  m=IMFExpansionN(data->E,e);

  Phi=MFNKM_CStar(mPhi,e);
  u=MFNV_CStar(vu,e);

  A=(double*)malloc(n*n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(A==NULL)
   {
    sprintf(IMFSphereOnExpansionErrorMsg,"Out of memory trying to allocate %d bytes\n",n*n*sizeof(double));
    MFSetError(e,12,RoutineName,IMFSphereOnExpansionErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0;
   }
#endif

  b=(double*)malloc(n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(b==NULL)
   {
    sprintf(IMFSphereOnExpansionErrorMsg,"Out of memory trying to allocate %d bytes\n",n*sizeof(double));
    MFSetError(e,12,RoutineName,IMFSphereOnExpansionErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0;
   }
#endif

  ipvt=(int*)malloc(n*sizeof(int));

#ifndef MFNOSAFETYNET
  if(ipvt==NULL)
   {
    sprintf(IMFSphereOnExpansionErrorMsg,"Out of memory trying to allocate %d bytes\n",n*sizeof(int));
    MFSetError(e,12,RoutineName,IMFSphereOnExpansionErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0;
   }
#endif

  vs=MFCreateKVector(k+1,e);
  s=MFKV_CStar(vs,e);
  vM0=MFCreateNVector(m,e);
  M0=MFNV_CStar(vM0,e);
  vdM=MFCreateNVector(m,e);
  dM=MFNV_CStar(vdM,e);

  R=data->R;
  for(i=0;i<k+1;i++)s[i]=0.;
  IMFEvaluateExpansion(data->E,vs,vM0,e);

  for(i=0;i<k+1;i++)s[i]=u[m+i];

  for(i=0;i<n;i++)b[i]=0.;

  for(j=0;j<k+1;j++)
   {
    IMFEvaluateExpansionDerivative(data->E,vs,j,vdM,e);
    A[0+(k+1)*j]=0.;for(i=0;i<m;i++)A[0+(k+1)*j]+=2*dM[i]*(u[i]-M0[i]);
 
    for(l=1;l<k+1;l++)A[l+(k+1)*j]=0.;
   }

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("  u: (%lf",u[0]);for(i=1;i<m+k+1;i++)printf(",%lf",u[i]);printf(")\n");fflush(stdout);
    printf("System for Tangent\n");
    IMFPrintFull(k+1,A,b,e);
   }
#endif

  lda=k+1;
  ldU=k+1;

  U=(double*)malloc((k+1)*(k+1)*sizeof(double));

#ifndef MFNOSAFETYNET
  if(U==NULL)
   {
    sprintf(IMFSphereOnExpansionErrorMsg,"Out of memory trying to allocate %d bytes\n",(k+1)*(k+1)*sizeof(double));
    MFSetError(e,12,RoutineName,IMFSphereOnExpansionErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0;
   }
#endif

  vr=(double*)malloc((k+1)*k*sizeof(double));

#ifndef MFNOSAFETYNET
  if(vr==NULL)
   {
    sprintf(IMFSphereOnExpansionErrorMsg,"Out of memory trying to allocate %d bytes\n",(k+1)*k*sizeof(double));
    MFSetError(e,12,RoutineName,IMFSphereOnExpansionErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0;
   }
#endif

  ldV=(k+1);
  V=(double*)malloc((k+1)*(k+1)*sizeof(double));

#ifndef MFNOSAFETYNET
  if(V==NULL)
   {
    sprintf(IMFSphereOnExpansionErrorMsg,"Out of memory trying to allocate %d bytes\n",(k+1)*(k+1)*sizeof(double));
    MFSetError(e,12,RoutineName,IMFSphereOnExpansionErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0;
   }
#endif

  sv=(double*)malloc((k+1)*sizeof(double));

#ifndef MFNOSAFETYNET
  if(sv==NULL)
   {
    sprintf(IMFSphereOnExpansionErrorMsg,"Out of memory trying to allocate %d bytes\n",(k+1)*sizeof(double));
    MFSetError(e,12,RoutineName,IMFSphereOnExpansionErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0;
   }
#endif

  i=-1;
  info=0;
  jobvl='N';
  jobvr='A';
  CALLDGESVD(&jobvl,&jobvr,&lda,&lda,A,&lda,sv,U,&ldU,V,&ldV,&t,&i,&info);
  lwork=MAX(4*(k-1),5*(k-1)-4);
  lwork=round(t);

  work=(double*)malloc(lwork*sizeof(double));

#ifndef MFNOSAFETYNET
  if(work==NULL)
   {
    sprintf(IMFSphereOnExpansionErrorMsg,"Out of memory trying to allocate %d bytes\n",lwork*sizeof(double));
    MFSetError(e,12,RoutineName,IMFSphereOnExpansionErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0;
   }
#endif

  info=0;
  jobvl='N';
  jobvr='A';
  CALLDGESVD(&jobvl,&jobvr,&lda,&lda,A,&lda,sv,U,&ldU,V,&ldV,work,&lwork,&info);

  msv=0;
  for(i=0;i<k+1;i++)
   {

#ifdef MFALLOWVERBOSE
    if(verbose){printf("Singular value %d is %lf\n",i,sv[i]);fflush(stdout);}
#endif

    if(fabs(sv[i])<1.e-7)
     {
      for(j=0;j<k+1;j++)vr[j+(k+1)*msv]=V[i+(k+1)*j];

#ifdef MFALLOWVERBOSE
      if(verbose){printf(" null vector (%lf",vr[0+(k+1)*msv]);for(j=1;j<k+1;j++)printf(",%lf",vr[j+(k+1)*msv]);printf(")\n");fflush(stdout);}
#endif

      msv++;
     }
   }

  t=0.;
  for(j=0;j<k+1;j++)t+=vr[j]*vr[j];
  if(fabs(t)>1.e-15)
   {
    t=1./sqrt(t);
    for(j=0;j<k+1;j++){vr[j]=vr[j]*t;}
   }

  for(i=1;i<msv;i++)
   {
    for(l=0;l<i;l++)
     {
      t=0.;
      for(j=0;j<k+1;j++)t+=vr[j+(k+1)*l]*vr[j+(k+1)*i];
      for(j=0;j<k+1;j++)vr[j+(k+1)*i]=vr[j+(k+1)*i]-t*vr[j+(k+1)*l];
     }

    t=0.;
    for(j=0;j<k+1;j++)t+=vr[j+(k+1)*i]*vr[j+(k+1)*i];
    if(fabs(t)>1.e-15)
     {
      t=1./sqrt(t);
      for(j=0;j<k+1;j++){vr[j+(k+1)*i]=vr[j+(k+1)*i]*t;}
     }
   }

  for(j=0;j<n*k;j++)Phi[j]=0.;
  for(j=0;j<k;j++)
   {
    for(l=0;l<k+1;l++)
     {
      IMFEvaluateExpansionDerivative(data->E,vs,l,vdM,e);
      for(i=0;i<m;i++)
        Phi[i+n*j]+=dM[i]*vr[l+(k+1)*j];
     }
    for(i=0;i<k+1;i++)Phi[i+m+n*j]=vr[i+(k+1)*j];
   }

#ifdef MFALLOWVERBOSE
  if(verbose){printf("Tangent is\n");MFPrintNKMatrix(stdout,mPhi,e);fflush(stdout);}
#endif

  free(A);
  free(b);
  free(ipvt);
  MFFreeKVector(vs,e);
  MFFreeNVector(vM0,e);
  MFFreeNVector(vdM,e);
  free(U);
  free(V);
  free(sv);
  free(vr);
  free(work);

  return 1;
 }

int MFIMFSphereOnEProjectForBB(MFNVector u, double *x, void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFSphereOnEProjectForBB"};
  int i;
  struct IMFSphereOnEData *data;

  data=(struct IMFSphereOnEData*)d;

 if(x==NULL)return IMFExpansionN(data->E,e);
  for(i=0;i<IMFExpansionN(data->E,e);i++)x[i]=MFNV_C(u,i,e);

  return 0;
 }

IMFExpansion IMFSphereOnExpansionGetLocal(MFAtlas c, MFNVector vu, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFSphereOnExpansionGetLocal"};
  double *u;
  double *du;
  double *ddu;
  double *Du;
  double *DDu;
  MFNKMatrix mPhi;
  IMFExpansion result;
  int i,j,l,n,k,m;
  double *A=NULL;
  double *b=NULL;
  MFNVector vM0,vdM,vdMi,vdMj;
  MFKVector vs;
  double *s,*M0,*dM,*dMi,*dMj;
  int ii,jj,ll,kk;

  struct IMFSphereOnEData *data;

  data=(struct IMFSphereOnEData*)MFIMFGetData(MFAtlasMF(c,e),e);

/*
  second derivative is?
  third derivative is?
*/

  n=MFAtlasN(c,e);
  k=MFAtlasK(c,e);
  m=IMFExpansionN(data->E,e);

  result=IMFCreateExpansion(m,k,e);
  u=MFNV_CStar(vu,e);
  mPhi=MFIMFTangentSpace(MFAtlasMF(c,e),vu,e);

  du=MFNKM_CStar(mPhi,e);
  ddu=(double*)malloc(n*n*k*sizeof(double));

#ifndef MFNOSAFETYNET
  if(ddu==NULL)
   {
    sprintf(IMFSphereOnExpansionErrorMsg,"Out of memory trying to allocate %d bytes\n",n*n*k*sizeof(double));
    MFSetError(e,12,RoutineName,IMFSphereOnExpansionErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  Du=(double*)malloc(m*k*sizeof(double));

#ifndef MFNOSAFETYNET
  if(Du==NULL)
   {
    sprintf(IMFSphereOnExpansionErrorMsg,"Out of memory trying to allocate %d bytes\n",m*k*sizeof(double));
    MFSetError(e,12,RoutineName,IMFSphereOnExpansionErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  DDu=(double*)malloc(m*m*k*sizeof(double));

#ifndef MFNOSAFETYNET
  if(DDu==NULL)
   {
    sprintf(IMFSphereOnExpansionErrorMsg,"Out of memory trying to allocate %d bytes\n",m*m*k*sizeof(double));
    MFSetError(e,12,RoutineName,IMFSphereOnExpansionErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  A=(double*)malloc(n*n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(A==NULL)
   {
    sprintf(IMFSphereOnExpansionErrorMsg,"Out of memory trying to allocate %d bytes\n",n*n*sizeof(double));
    MFSetError(e,12,RoutineName,IMFSphereOnExpansionErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  b=(double*)malloc(n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(b==NULL)
   {
    sprintf(IMFSphereOnExpansionErrorMsg,"Out of memory trying to allocate %d bytes\n",n*sizeof(double));
    MFSetError(e,12,RoutineName,IMFSphereOnExpansionErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  vs=MFCreateKVector(k+1,e);
  s=MFKV_CStar(vs,e);
  vM0=MFCreateNVector(m,e);
  M0=MFNV_CStar(vM0,e);
  vdM=MFCreateNVector(m,e);
  dM=MFNV_CStar(vdM,e);
  vdMi=MFCreateNVector(m,e);
  dMi=MFNV_CStar(vdMi,e);
  vdMj=MFCreateNVector(m,e);
  dMj=MFNV_CStar(vdMj,e);

  for(i=0;i<k+1;i++)s[i]=0.;
  IMFEvaluateExpansion(data->E,vs,vM0,e);

  for(i=0;i<k+1;i++)s[i]=u[m+i];

  for(ii=0;ii<k;ii++)
   {
    for(jj=0;jj<k;jj++)
     {
      b[0]=0.;
      for(ll=0;ll<k+1;ll++)
       {
        for(kk=0;kk<k+1;kk++)
         {
          IMFEvaluateExpansionDerivative(data->E,vs,ll,vdMi,e);
          IMFEvaluateExpansionDerivative(data->E,vs,kk,vdMj,e);
          for(i=0;i<m;i++)b[0]-=2*dMi[i]*du[ll+m+n*ii]*dMj[i]*du[kk+m+n*jj];
         }
       }
      for(i=1;i<k+1;i++)b[i]=0.;

      for(j=0;j<k+1;j++)
       {
        A[0+(k+1)*j]=0.;
        IMFEvaluateExpansionDerivative(data->E,vs,j,vdM,e);
        for(i=0;i<m;i++)A[0+(k+1)*j]+=2*dM[i]*(u[i]-M0[i]);
    
        for(l=1;l<k+1;l++)A[l+(k+1)*j]=du[j+m+n*(l-1)];
       }

      MFSolveFull(k+1,A,b,e);
      for(i=0;i<k+1;i++)ddu[i+m+n*(ii+k*jj)]=b[i];
     }
   }

  for(ii=0;ii<k;ii++)
   {
    for(jj=0;jj<k;jj++)
     {
      for(i=0;i<m;i++)
       {
        ddu[i+n*(ii+k*jj)]=0.;
        for(j=0;j<k+1;j++)
         {
          for(ll=0;ll<k+1;ll++)
           {
            IMFEvaluateExpansionDerivative(data->E,vs,ll,vdM,e);
            ddu[i+n*(ii+k*jj)]+=dM[i]*ddu[ll+m+n*(ii+k*jj)];
            for(kk=0;kk<k+1;kk++)
             {
              IMFEvaluateExpansionSecondDerivative(data->E,vs,ll,kk,vdM,e);
              ddu[i+n*(ii+k*jj)]+=dM[i]*du[ll+m+n*ii]*du[kk+m+n*jj];
             }
           }
         }
       }
     }
   }

  for(j=0;j<k;j++)
   {
    for(i=0;i<m;i++)
      Du[i+m*j]=du[i+n*j];
    for(l=0;l<k;l++)
     {
      for(i=0;i<m;i++)
        DDu[i+m*(j+k*l)]=du[i+n*(j+k*l)];
     }
   }

  IMFExpansionSetDerivatives(result,u,Du,DDu,NULL,e);

  MFFreeNKMatrix(mPhi,e);
  MFFreeNVector(vM0,e);
  MFFreeNVector(vdM,e);
  MFFreeNVector(vdMi,e);
  MFFreeNVector(vdMj,e);
  MFFreeKVector(vs,e);
  free(Du);
  free(DDu);
  free(ddu);
  free(A);
  free(b);
  return result;
 }

double MFScaleIMFSphereOnE(int n,int k,MFNVector vu,MFNKMatrix mPhi,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFScaleIMFSphereOnE"};
  struct IMFSphereOnEData *data;
  double R,Rf;

  data=(struct IMFSphereOnEData*)d;

  R=IMFExpansionR(data->E,data->eps,e);
  Rf=IMFFlowR(data->F,data->eps,vu,data->p0,mPhi,e);
  if(Rf<R)R=Rf;

  if(R>data->R||R!=R)R=data->R;

  return R;
 }

MFNVector MFSOENVectorFactory(MFImplicitMF SOE,MFErrorHandler e)
 {
  static char RoutineName[]={"MFSOENVectorFactory"};
  struct IMFSphereOnEData *data;

  data=(struct IMFSphereOnEData *)MFIMFGetData(SOE,e);

  return MFCreateNVector(data->n,e);
 }

MFNKMatrix MFSOENKMatrixFactory(MFImplicitMF SOE,MFErrorHandler e)
 {
  static char RoutineName[]={"MFSOENKMatrixFactory"};
  struct IMFSphereOnEData *data;

  data=(struct IMFSphereOnEData *)MFIMFGetData(SOE,e);

  return MFCreateNKMatrixWithData(data->n,data->k,NULL,e);
 }

#ifdef __cplusplus
}
#endif
