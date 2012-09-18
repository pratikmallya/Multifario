/* 
 *  PROGRAM NAME:  multifario
 *
 *  (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
 *  CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
 *
 *  Please refer to the LICENSE file in the top directory
 *
 *      author: Mike Henderson mhender@watson.ibm.com
 *      date:   February 22, 1999
 */

static char *id="@(#) $Id: MFPendula.c,v 1.5 2011/07/21 17:42:46 mhender Exp $";

static char MFPendulaMFErrorHandlerMsg[256]="";

extern double MFEpsilon;

#include <MFImplicitMF.h>
#include <MFNKMatrix.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <MFFortran.h>

#define GU GuPendula
#define GUU GuuPendula

#ifdef __cplusplus
 extern "C" {
#endif

/* Space for LU decompositions */

static double *A,*B,*C,*D,*R,*S,*X,*Y;
static int   *Pivots;

static int nBorders,nBandsU,nBandsL;
static int ldA,ldB,ldC,ldD;

static void MFFreePendulaData(void*,MFErrorHandler);
static int MFProjectPendula(int,int,MFNVector,MFNKMatrix,MFNVector,void*,int*,MFErrorHandler);
static int MFTangentPendula(int,int,MFNVector,MFNKMatrix,void*,MFErrorHandler);
static int MFTangentPendulaWithGuess(int,int,MFNVector,MFNKMatrix,MFNKMatrix,void*,MFErrorHandler);
static double MFScalePendula(int,int,MFNVector,MFNKMatrix,void*,MFErrorHandler);
static void MFPendulaInitializeSolver(int,MFErrorHandler);
static double GPendula(int,double*,int,double,double,double,MFErrorHandler);
static double GuPendula(int,int,double*,int,double,double,double,MFErrorHandler);
static double GuuPendula(int,int,int,double*,int,double,double,double,MFErrorHandler);
static double TGuPendula(int,int,double*,int,double,double,double,MFErrorHandler);
static double TGuuPendula(int,int,int,double*,int,double,double,double,MFErrorHandler);
static int MFTangentPendulaSingle(int,int,int,double*,double*,double*,double*,void*,MFErrorHandler);
static void MFCurvaturePendula(int,int,double*,double*,double*,double*,double*,double*,void*,MFErrorHandler);
static void MFWritePendulaData(FILE*,void*,MFErrorHandler);
static MFImplicitMF MFReadPendula(FILE*,MFErrorHandler);

static int MFPendulaProjectToSave(MFNVector,double*,void*,MFErrorHandler);
static int MFPendulaProjectToDraw(MFNVector,double*,void*,MFErrorHandler);
static int MFPendulaProjectForBB(MFNVector,double*,void*,MFErrorHandler);

MFNVector MFNVectorFactory(MFImplicitMF,MFErrorHandler);
MFNKMatrix MFNKMatrixFactory(MFImplicitMF,MFErrorHandler);

struct MFPendulaData
 {
  int nt;
  double kappa;
  double gamma;
  double R;
 };

MFImplicitMF MFIMFCreatePendula(int nt,double kappa,double gamma,int windingno, MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFCreatePendula"};
  MFImplicitMF pendula;
  int i;
  MFNSpace space;
  struct MFPendulaData *data;

  pendula=MFIMFCreateBaseClass(2*nt+3,2,"Pendula",e);

  space=MFCreateWeightedNSpace(2*nt+3,2*nt,1./nt,e);
  MFIMFSetSpace(pendula,space,e);
  MFFreeNSpace(space,e);

  data=(struct MFPendulaData*)malloc(sizeof(struct MFPendulaData)); /*done*/

#ifndef MFNOSAFETYNET
  if(data==NULL)
   {
    sprintf(MFPendulaMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct MFPendulaData));
    MFSetError(e,12,RoutineName,MFPendulaMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  data->nt=nt;
  data->kappa=kappa;
  data->gamma=gamma;
  data->R=(double)windingno;

  MFIMFSetData(pendula,(void*)data,e);
  MFIMFSetFreeData(pendula,MFFreePendulaData,e);
  MFIMFSetProject(pendula,MFProjectPendula,e);
  MFIMFSetTangent(pendula,MFTangentPendula,e);
  MFIMFSetTangentWithGuess(pendula,MFTangentPendulaWithGuess,e);
  MFIMFSetScale(pendula,MFScalePendula,e);
  MFIMFSetWriteData(pendula,MFWritePendulaData,e);
  MFIMFSetProjectForSave(pendula,MFPendulaProjectToSave,e);
  MFIMFSetProjectForDraw(pendula,MFPendulaProjectToDraw,e);
  MFIMFSetProjectForBB(pendula,MFPendulaProjectForBB,e);

  MFIMFSetVectorFactory(pendula,MFNVectorFactory,e);
  MFIMFSetMatrixFactory(pendula,MFNKMatrixFactory,e);

  MFPendulaInitializeSolver(nt,e);

  return pendula;
 }

void MFFreePendulaData(void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFFreePendulaData"};
  free((struct MFPendulaData*)d);
  return;
 }

double MFScalePendula(int n,int k,MFNVector vu,MFNKMatrix mPhi, void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFScalePendula"};
  double *a00,*a01,*a11;
  double A,B,C,l;
  double r;
  int i,nt;
  double *u,*Phi;
  int verbose=0;

  Phi=MFNKM_CStar(mPhi,e);
  u=MFNV_CStar(vu,e);

  nt=(n-3)/2;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("  MFScalePendula \n");fflush(stdout);}
#endif

  a00=(double*)malloc(n*sizeof(double)); /*done*/

#ifndef MFNOSAFETYNET
  if(a00==NULL)
   {
    sprintf(MFPendulaMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(double));
    MFSetError(e,12,RoutineName,MFPendulaMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return -1.;
   }
#endif

  a01=(double*)malloc(n*sizeof(double)); /*done*/

#ifndef MFNOSAFETYNET
  if(a01==NULL)
   {
    sprintf(MFPendulaMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(double));
    MFSetError(e,12,RoutineName,MFPendulaMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return -1.;
   }
#endif

  a11=(double*)malloc(n*sizeof(double)); /*done*/

#ifndef MFNOSAFETYNET
  if(a11==NULL)
   {
    sprintf(MFPendulaMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(double));
    MFSetError(e,12,RoutineName,MFPendulaMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return -1.;
   }
#endif

  MFCurvaturePendula(n,k,u,Phi,Phi+n,a00,a01,a11,d,e);

  A=0.;
  for(i=0;i<n-3;i++)A+=a00[i]*a00[i]/nt;
  for(i=n-3;i<n;i++)A+=a00[i]*a00[i];
  A=sqrt(A);

#ifdef MFALLOWVERBOSE
   if(verbose){printf("  |a_00|=%lf\n",A);fflush(stdout);}
#endif

  B=0.;
  for(i=0;i<n-3;i++)B+=a01[i]*a01[i]/nt;
  for(i=n-3;i<n;i++)B+=a01[i]*a01[i];
  B=sqrt(B);

#ifdef MFALLOWVERBOSE
   if(verbose){printf("  |a_01|=|a_10|=%lf\n",B);fflush(stdout);}
#endif

  C=0.;
  for(i=0;i<n-3;i++)C+=a11[i]*a11[i]/nt;
  for(i=n-3;i<n;i++)C+=a11[i]*a11[i];
  C=sqrt(C);

#ifdef MFALLOWVERBOSE
   if(verbose){printf("  |a_11|=%lf\n",C);fflush(stdout);}
#endif

  free(a00);
  free(a01);
  free(a11);

  if(A+C>0)
    l=.5*(A+C)+.5*sqrt((A-C)*(A-C)+4*B*B);
   else
    l=-.5*(A+C)+.5*sqrt((A-C)*(A-C)+4*B*B);
  r=.6*sqrt(2*MFEpsilon/l);

#ifdef MFALLOWVERBOSE
   if(verbose){printf("%s r=%lf\n",RoutineName,r);fflush(stdout);}
#endif

  if(r>.4)r=.4;

  return r;
 }

double GPendula(int i,double *point,int nt, double kappa, double gamma, double r, MFErrorHandler e)
 {
  static char RoutineName[]={"GPendula"};
  double w,Ic,In;
  double x,h;
  double W;
  double result;
  int n;
  int j,m;
  int jm,j0,jp;
  double f,fm,f0,fp;

#ifdef MFNOCONFIDENCE
  if(i<0 || i>2*nt+1)
   {
    sprintf(MFPendulaMFErrorHandlerMsg,"i<0 || i>n-k, in G i=%d, n-k=%d\n",i,2*nt+1);
    MFSetError(e,12,RoutineName,MFPendulaMFErrorHandlerMsg,__LINE__,__FILE__);
    return 0.;
   }
#endif

  w=2*3.1415926/(10*point[2*nt]);
  In=point[2*nt+1];
  Ic=3.1415926*kappa*point[2*nt+2];
  Ic=point[2*nt+2];

/*
 
   EVALUATES G (U,L)
                                 I
      U IS:
             1     2     3        N-4  N-3   N-2  N-1 N
 
           P1(1),P2(1),P1(2),...,P1(M),P2(M), T , In, Ic
 
*/

  h=r*2.*3.1415926/nt;

  j=i/2;
  m=i%2;

  x=j*h;

  jm=j-1;
  j0=j;
  jp=j+1;

  if(j0== 0)jm=nt-1;
  if(j0==nt-1)jp=0;

  result=0.;
  if(i<2*nt)
   {
    switch(m)
     {
      case 0:
       fp=w+sin(w*(x+h)+point[2*jp])+kappa*(point[2*jp]-point[2*jp+1])-In-Ic;
       f0=w+sin(w*x  +point[2*j0])+kappa*(point[2*j0]-point[2*j0+1])-In-Ic;
       fm=w+sin(w*(x-h)+point[2*jm])+kappa*(point[2*jm]-point[2*jm+1])-In-Ic;
       f=.25*(fp+2*f0+fm);
       result=w*w*(point[2*jp]-2.*point[2*j0]+point[2*jm])+gamma*w*(point[2*jp]-point[2*jm])*h+f*h*h;
       break;
      case 1:
       fp=w+sin(w*(x+h)+point[2*jp+1])+kappa*(point[2*jp+1]-point[2*jp])-In+Ic;
       f0=w+sin(w*(x  )+point[2*j0+1])+kappa*(point[2*j0+1]-point[2*j0])-In+Ic;
       fm=w+sin(w*(x-h)+point[2*jm+1])+kappa*(point[2*jm+1]-point[2*jm])-In+Ic;
       f=.25*(fp+2*f0+fm);
       result=w*w*(point[2*jp+1]-2.*point[2*j0+1]+point[2*jm+1])+gamma*w*(point[2*jp+1]-point[2*jm+1])*h+f*h*h;
       break;
     }
   }else{
    result=.5*(point[2*(nt-1)]+point[2*(nt-1)+1]);    /* Phase Constraint */
   }

  return(result);
 }

void MFPendulaInitializeSolver(int nt, MFErrorHandler e)
 {
  static char RoutineName[]={"MFPendulaInitializeSolver"};
  int nEqs,n;
  int k;

  nEqs=2;
  k=2;
  n=nEqs*nt+3;
  nBorders=nEqs+1;
  nBandsU=2*nEqs-1;
  nBandsL=2*nEqs-1;

  ldA=2*nBandsL+nBandsU+1;
  ldB=n-k-nBorders;
  ldC=nBorders+k;
  ldD=nBorders+k;

  A=(double*)malloc(ldA*(n-k-nBorders)*sizeof(double)); /*done*/

#ifndef MFNOSAFETYNET
  if(A==NULL)
   {
    sprintf(MFPendulaMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",ldA*(n-k-nBorders)*sizeof(double));
    MFSetError(e,12,RoutineName,MFPendulaMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  B=(double*)malloc(ldB*(nBorders+k)*sizeof(double)); /*done*/

#ifndef MFNOSAFETYNET
  if(B==NULL)
   {
    sprintf(MFPendulaMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",ldB*(nBorders+k)*sizeof(double));
    MFSetError(e,12,RoutineName,MFPendulaMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  C=(double*)malloc(ldC*(n-k-nBorders)*sizeof(double)); /*done*/

#ifndef MFNOSAFETYNET
  if(C==NULL)
   {
    sprintf(MFPendulaMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",ldC*(n-k-nBorders)*sizeof(double));
    MFSetError(e,12,RoutineName,MFPendulaMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  D=(double*)malloc(ldD*(nBorders+k)*sizeof(double)); /*done*/

#ifndef MFNOSAFETYNET
  if(D==NULL)
   {
    sprintf(MFPendulaMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",ldD*(nBorders+k)*sizeof(double));
    MFSetError(e,12,RoutineName,MFPendulaMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  R=(double*)malloc((n-k-nBorders)*sizeof(double)); /*done*/

#ifndef MFNOSAFETYNET
  if(R==NULL)
   {
    sprintf(MFPendulaMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",(n-k-nBorders)*sizeof(double));
    MFSetError(e,12,RoutineName,MFPendulaMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  S=(double*)malloc((nBorders+k)*sizeof(double)); /*done*/

#ifndef MFNOSAFETYNET
  if(S==NULL)
   {
    sprintf(MFPendulaMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",(nBorders+k)*sizeof(double));
    MFSetError(e,12,RoutineName,MFPendulaMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  X=(double*)malloc((nBorders+k+nBandsL+1)*(nBorders+k+nBandsL+1)*sizeof(double)); /*done*/

#ifndef MFNOSAFETYNET
  if(X==NULL)
   {
    sprintf(MFPendulaMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",(nBorders+k+nBandsL+1)*(nBorders+k+nBandsL+1)*sizeof(double));
    MFSetError(e,12,RoutineName,MFPendulaMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  Y=(double*)malloc((nBorders+k+nBandsL+1)*sizeof(double)); /*done*/

#ifndef MFNOSAFETYNET
  if(Y==NULL)
   {
    sprintf(MFPendulaMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",(nBorders+k+nBandsL+1)*sizeof(double));
    MFSetError(e,12,RoutineName,MFPendulaMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  Pivots=(int*)malloc(n*sizeof(int));

#ifndef MFNOSAFETYNET
  if(Pivots==NULL)
   {
    sprintf(MFPendulaMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(int));
    MFSetError(e,12,RoutineName,MFPendulaMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  return;
 }

int MFProjectPendula(int n,int k,MFNVector vu0,MFNKMatrix mPhi,MFNVector vu,void *d,int *index, MFErrorHandler e)
 {
  static char RoutineName[]={"MFProjectPendula"};
  static double eps=.001;
  static double emax;
  static double error;
  static double delta;
  static double factor;
  static int i,j;
  static int i0,i1;
  static int itimes;
  static int nBands,bandW;
  static int ierr;
  static double *tangentA;
  static double *tangentB;
  double *u,*u0,*Phi;

  static int nt;
  static int na,nb;
  static double kappa,gamma,r;

  static int verbose=0;
  static int test=0;

  u0=MFNV_CStar(vu0,e);
  Phi=MFNKM_CStar(mPhi,e);
  u=MFNV_CStar(vu,e);

  nt=((struct MFPendulaData*)d)->nt;
  kappa=((struct MFPendulaData*)d)->kappa;
  gamma=((struct MFPendulaData*)d)->gamma;
  r=((struct MFPendulaData*)d)->R;

  verbose=0;
  test=0;

  tangentA=Phi;
  tangentB=Phi+n;

  nBands=nBandsL+nBandsU+1;

  for(i=0;i<n;i++)u[i]=u0[i];

/* Test basis for the tangent Space */

#ifdef MFALLOWVERBOSE
  if(test)
   {
    error=0.;
    for(i=0;i<2*nt;i++)
      error+=tangentA[i]*tangentA[i]/nt;
    for(i=2*nt;i<n;i++)
      error+=tangentA[i]*tangentA[i];
    printf("direction A.direction A: %e\n",error);

    error=0.;
    for(i=0;i<2*nt;i++)
      error+=tangentA[i]*tangentB[i]/nt;
    for(i=2*nt;i<n;i++)
      error+=tangentA[i]*tangentB[i];
    printf("direction A.direction B: %e\n",error);

    error=0.;
    for(i=0;i<2*nt;i++)
      error+=tangentB[i]*tangentB[i]/nt;
    for(i=2*nt;i<n;i++)
      error+=tangentB[i]*tangentB[i];
    printf("direction B.direction B: %e\n",error);
   }
#endif

  error=1.;
  itimes=0;
  while(error>1.e-10)
   {
    for(j=0;j<n-k-nBorders;j++)
     {
      i0=j-nBandsU;
      if(i0<0)i0=0;
      i1=j+nBandsL;
      if(i1>n-k-nBorders-1)i1=n-k-nBorders-1;
      for(i=i0;i<i1+1;i++)
        A[i-j+nBands-1+ldA*j]=GU(i,j,u,nt,kappa,gamma,r,e);
     }

    for(i=0;i<n-k-nBorders;i++)
     {
      for(j=0;j<nBorders+k;j++)
        B[i+ldB*j]=GU(i,n-k-nBorders+j,u,nt,kappa,gamma,r,e);
     }

    for(i=0;i<nBorders;i++)
     {
      for(j=0;j<n-k-nBorders;j++)
        C[i+ldC*j]=GU(n-k-nBorders+i,j,u,nt,kappa,gamma,r,e);
     }
    for(j=0;j<n-k-nBorders;j++)
     {
      if(j<2*nt)
       {
        C[nBorders+0+ldC*j]=tangentA[j]/nt;
        C[nBorders+1+ldC*j]=tangentB[j]/nt;
       }else{
        C[nBorders+0+ldC*j]=tangentA[j];
        C[nBorders+1+ldC*j]=tangentB[j];
       }
     }

    for(i=0;i<nBorders;i++)
     {
      for(j=0;j<nBorders+k;j++)
        D[i+ldD*j]=GU(n-k-nBorders+i,n-k-nBorders+j,u,nt,kappa,gamma,r,e);
     }
    for(j=0;j<nBorders+k;j++)
     {
      if(n-k-nBorders+j<2*nt)
       {
        D[nBorders+0+ldD*j]=tangentA[n-k-nBorders+j]/nt;
        D[nBorders+1+ldD*j]=tangentB[n-k-nBorders+j]/nt;
       }else{
        D[nBorders+0+ldD*j]=tangentA[n-k-nBorders+j];
        D[nBorders+1+ldD*j]=tangentB[n-k-nBorders+j];
       }
     }

    for(i=0;i<n-k-nBorders;i++)
      R[i]=-GPendula(i,u,nt,kappa,gamma,r,e);

    for(i=0;i<nBorders;i++)
      S[i]=-GPendula(n-k-nBorders+i,u,nt,kappa,gamma,r,e);

    S[nBorders+0]=0.;
    S[nBorders+1]=0.;
    for(i=0;i<n;i++)
     {
      if(i<2*nt)
       {
        S[nBorders+0]+=tangentA[i]*(u0[i]-u[i])/nt;
        S[nBorders+1]+=tangentB[i]*(u0[i]-u[i])/nt;
       }else{
        S[nBorders+0]+=tangentA[i]*(u0[i]-u[i]);
        S[nBorders+1]+=tangentB[i]*(u0[i]-u[i]);
       }
     }

/* Calculate the size of the right hand side. */

    error=0.;
    for(i=0;i<n-k-nBorders;i++)
     {
      if(i<2*nt)
       {
        error+=R[i]*R[i]/nt;
       }else{
        error+=R[i]*R[i];
       }
     }
    for(i=0;i<nBorders+k;i++)
     {
      if(i+n-k-nBorders<2*nt)
       {
        error+=S[i]*S[i]/nt;
       }else{
        error+=S[i]*S[i];
       }
     }
    error=sqrt(error);

/* Compute the correction. */

#ifdef MFALLOWVERBOSE
    if(verbose)printf("\n                   Size of rhs: %e\n",error);
#endif

    na=n-k-nBorders;
    nb=nBorders+k;
    CALLDBOFA(A,&ldA,&na,&nBandsL,&nBandsU,B,&ldB,&nb,C,&ldC,D,&ldD,X,NULL,Pivots,&ierr);

#ifndef MFNOSAFETYNET
    if(ierr!=0)
     {
      sprintf(MFPendulaMFErrorHandlerMsg," Problem with factor, zero on diagonal %d",ierr);
      MFSetError(e,12,RoutineName,MFPendulaMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return 0;
     }
#endif

    ierr=0;
    CALLDBOSL(A,&ldA,&na,&nBandsL,&nBandsU,B,&ldB,&nb,C,&ldC,D,&ldD,X,Y,R,S,NULL,Pivots,&ierr);

/* Is GPendula(u+eps*delta)-GPendula(u) ~ -eps*GPendula(u)  */
/*    P(u+eps*delta)-P(u) ~ -eps*P(u)  */
/*    N(u+eps*delta)-N(u) ~ -eps*N(u)? */

/* Calculate the size of the Correction. */

    delta=0.;
    for(i=0;i<n-k-nBorders;i++)
     {
      if(i<2*nt)
       {
        delta+=R[i]*R[i]/nt;
       }else{
        delta+=R[i]*R[i];
       }
     }
    for(i=0;i<nBorders+k;i++)
     {
      if(i+n-k-nBorders<2*nt)
       {
        delta+=S[i]*S[i]/nt;
       }else{
        delta+=S[i]*S[i];
       }
     }
    delta=sqrt(delta);

#ifdef MFALLOWVERBOSE
    if(verbose)printf("\n                   Size of correction: %e\n",delta);
#endif

/* Apply the Correction */

    factor=1.;

    for(i=0;i<n-k-nBorders;i++)
      u[i]+=factor*R[i];
    for(i=0;i<nBorders+k;i++)
      u[n-k-nBorders+i]+=factor*S[i];

    itimes++;

#ifndef MFNOSAFETYNET
    if(itimes>20)
     {
      sprintf(MFPendulaMFErrorHandlerMsg," Too many refine iterations %d",itimes);
      MFSetError(e,12,RoutineName,MFPendulaMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return 0;
     }
#endif

   }

#ifdef MFALLOWVERBOSE
  if(verbose)printf("\n");
#endif

  *index=0;
  return 1;
 }

int MFTangentPendula(int n,int k,MFNVector vu,MFNKMatrix mPhi,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFTangentPendula"};
  int i;
  double *Phi0;
  double *u,*Phi;

  Phi=MFNKM_CStar(mPhi,e);
  u=MFNV_CStar(vu,e);

  Phi0=(double*)malloc(n*k*sizeof(double)); /*done*/

#ifndef MFNOSAFETYNET
  if(Phi0==NULL)
   {
    sprintf(MFPendulaMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*k*sizeof(double));
    MFSetError(e,12,RoutineName,MFPendulaMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0;
   }
#endif

  for(i=0;i<n*k;i++)Phi0[i]=0.;
  Phi0[n-2]=1.;
  Phi0[n-1+n]=1.;
   
  MFTangentPendulaSingle(n,k,0,u,Phi0,Phi0+n,Phi,d,e);
  MFTangentPendulaSingle(n,k,1,u,Phi ,Phi0+n,Phi,d,e);
  free(Phi0);
  return 1;
 }

int MFTangentPendulaWithGuess(int n,int k,MFNVector vu,MFNKMatrix mPhi0,MFNKMatrix mPhi,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFTangentPendulaWithGuess"};
  double *u,*Phi,*Phi0;

  Phi=MFNKM_CStar(mPhi,e);
  Phi0=MFNKM_CStar(mPhi0,e);
  u=MFNV_CStar(vu,e);

  MFTangentPendulaSingle(n,k,0,u,Phi0,Phi0+n,Phi,d,e);
  MFTangentPendulaSingle(n,k,1,u,Phi ,Phi0+n,Phi,d,e);
  return 1;
 }

int MFTangentPendulaSingle(int n,int k,int t, double *u,double *tangentA,double *tangentB,double *Phi, void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFTangentPendulaSingle"};
  static double *result;
  static int i,j;
  static int i0,i1;
  static int verbose=0;
  static int test=0;
  static int nBands,bandW;
  static int na,nb;
  static int ierr;
  static double s;
  static int nt;
  static double kappa,gamma,r;

  nt=((struct MFPendulaData*)d)->nt;
  kappa=((struct MFPendulaData*)d)->kappa;
  gamma=((struct MFPendulaData*)d)->gamma;
  r=((struct MFPendulaData*)d)->R;

#ifdef MFALLOWVERBOSE
  if(verbose||test)printf("surfaceTangentBordered direction %d\n",t);

  if(verbose)
   {
    printf("TangentA:\n");
    printf("[");
    if(n<10)
     {
      for(i=0;i<n;i++)
       {
        if(i>0)printf(",");
        printf("%7.3lf",tangentA[i]);
       }
     }else{
      for(i=0;i<4;i++)
       {
        if(i>0)printf(",");
        printf("%7.3lf",tangentA[i]);
       }
      printf(",...");
      for(i=n-3;i<n;i++)
       {
        printf(",%7.3lf",tangentA[i]);
       }
     }
    printf("]\n");fflush(stdout);
   
    printf("TangentB:\n");
    printf("[");
    if(n<10)
     {
      for(i=0;i<n;i++)
       {
        if(i>0)printf(",");
        printf("%7.3lf",tangentB[i]);
       }
     }else{
      for(i=0;i<4;i++)
       {
        if(i>0)printf(",");
        printf("%7.3lf",tangentB[i]);
       }
      printf(",...");
      for(i=n-3;i<n;i++)
       {
        printf(",%7.3lf",tangentB[i]);
       }
     }
    printf("]\n");fflush(stdout);
   }

  if(test)
   {
    double e;

    e=0.;
    for(i=0;i<2*nt;i++)
      e+=tangentA[i]*tangentA[i]/nt;
    for(i=2*nt;i<n;i++)
      e+=tangentA[i]*tangentA[i];
    printf("direction A.direction A: %e\n",e);

    e=0.;
    for(i=0;i<2*nt;i++)
      e+=tangentA[i]*tangentB[i]/nt;
    for(i=2*nt;i<n;i++)
      e+=tangentA[i]*tangentB[i];
    printf("direction A.direction B: %e\n",e);

    e=0.;
    for(i=0;i<2*nt;i++)
      e+=tangentB[i]*tangentB[i]/nt;
    for(i=2*nt;i<n;i++)
      e+=tangentB[i]*tangentB[i];
    printf("direction B.direction B: %e\n",e);
   }
#endif

  nBands=nBandsL+nBandsU+1;

  for(j=0;j<n-k-nBorders;j++)
   {
    i0=j-nBandsU;
    if(i0<0)i0=0;
    i1=j+nBandsL;
    if(i1>n-k-nBorders-1)i1=n-k-nBorders-1;
    for(i=i0;i<i1+1;i++)
      A[i-j+nBands-1+ldA*j]=GU(i,j,u,nt,kappa,gamma,r,e);
   }

  for(i=0;i<n-k-nBorders;i++)
   {
    for(j=0;j<nBorders+k;j++)
      B[i+ldB*j]=GU(i,n-k-nBorders+j,u,nt,kappa,gamma,r,e);
   }

  for(i=0;i<nBorders;i++)
   {
    for(j=0;j<n-k-nBorders;j++)
      C[i+ldC*j]=GU(n-k-nBorders+i,j,u,nt,kappa,gamma,r,e);
   }
  for(j=0;j<n-k-nBorders;j++)
   {
    if(j<2*nt)
     {
      C[nBorders+0+ldC*j]=tangentA[j]/nt;
      C[nBorders+1+ldC*j]=tangentB[j]/nt;
     }else{
      C[nBorders+0+ldC*j]=tangentA[j];
      C[nBorders+1+ldC*j]=tangentB[j];
     }
   }

  for(i=0;i<nBorders;i++)
   {
    for(j=0;j<nBorders+k;j++)
      D[i+ldD*j]=GU(n-k-nBorders+i,n-k-nBorders+j,u,nt,kappa,gamma,r,e);
   }
  for(j=0;j<nBorders+k;j++)
   {
    if(n-k-nBorders+j<2*nt)
     {
      D[nBorders+0+ldD*j]=tangentA[n-k-nBorders+j]/nt;
      D[nBorders+1+ldD*j]=tangentB[n-k-nBorders+j]/nt;
     }else{
      D[nBorders+0+ldD*j]=tangentA[n-k-nBorders+j];
      D[nBorders+1+ldD*j]=tangentB[n-k-nBorders+j];
     }
   }
  for(i=0;i<n-k-nBorders;i++)R[i]=0.;

  for(i=0;i<nBorders+k;i++)S[i]=0.;
  S[nBorders+t]=1.;


/* Calculate the size of the right hand side. */

/* Compute the tangent. */

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf(" System:\n\n");
    for(i=0;i<n-k-nBorders;i++)
     {
      printf("[");
      for(j=0;j<n-k-nBorders;j++)
       {
        if(j>i-nBandsL-1 && j<i+nBandsU+1)
         {
          printf("%7.3lf",A[i-j+nBands-1+ldA*j]);
         }else{
          printf("       ");
         }
       }
      printf(" | ");
      for(j=0;j<nBorders+k;j++)printf("%7.3lf",B[i+ldB*j]);
      printf("]");
      printf(" [%7.3lf]\n",R[i]);
     }
  
    printf("[");
    for(j=0;j<n-k-nBorders;j++)printf("--------");
    printf(" + ");
    for(j=0;j<nBorders+k;j++)printf("--------");
    printf("]");
    printf(" [-------]\n");

    for(i=0;i<nBorders+k;i++)
     {
      printf("[");
      for(j=0;j<n-k-nBorders;j++)printf("%7.3lf",C[i+ldC*j]);
      printf(" | ");
      for(j=0;j<nBorders+k;j++)printf("%7.3lf",D[i+ldC*j]);
      printf("]");
      printf(" [%7.3lf]\n",S[i]);
     }
   }
#endif

  na=n-k-nBorders;
  nb=nBorders+k;
  CALLDBOFA(A,&ldA,&na,&nBandsL,&nBandsU,B,&ldB,&nb,C,&ldC,D,&ldD,X,NULL,Pivots,&ierr);

#ifdef MFNOCONFIDENCE
  if(ierr!=0)
   {
    sprintf(MFPendulaMFErrorHandlerMsg," Problem with factor, zero on diagonal %d",ierr);
    MFSetError(e,12,RoutineName,MFPendulaMFErrorHandlerMsg,__LINE__,__FILE__);
    return;
   }
#endif

  ierr=0;
  CALLDBOSL(A,&ldA,&na,&nBandsL,&nBandsU,B,&ldB,&nb,C,&ldC,D,&ldD,X,Y,R,S,NULL,Pivots,&ierr);

/* Extract the result */

  for(i=0;i<n-k-nBorders;i++)
    Phi[i+n*t]=R[i];
  for(i=0;i<nBorders+k;i++)
    Phi[n-k-nBorders+i+n*t]=S[i];

/* Normalize */

  s=0.;
  for(i=0;i<2*nt;i++)s+=Phi[i+n*t]*Phi[i+n*t]/nt;
  for(i=2*nt;i<n;i++)s+=Phi[i+n*t]*Phi[i+n*t];

#ifdef MFNOCONFIDENCE
  if(fabs(s)<.00001)
   {
    sprintf(MFPendulaMFErrorHandlerMsg,"Warning, result small in norm -- MFPendula Single Tangent\n");
    MFSetError(e,12,RoutineName,MFPendulaMFErrorHandlerMsg,__LINE__,__FILE__);
    return;
   }
#endif
  s=1./sqrt(s);

  for(i=0;i<n;i++)Phi[i+n*t]=s*Phi[i+n*t];

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("Tangent (t=%d):\n",t);
    printf("[");
    if(n<10)
     {
      for(i=0;i<n;i++)
       {
        if(i>0)printf(",");
        printf("%7.3lf",Phi[i+n*t]);
       }
     }else{
      for(i=0;i<4;i++)
       {
        if(i>0)printf(",");
        printf("%7.3lf",Phi[i+n*t]);
       }
      printf(",...");
      for(i=n-3;i<n;i++)
       {
        printf(",%7.3lf",Phi[i+n*t]);
       }
     }
    printf("]\n");fflush(stdout);
    printf("done surfaceTangentBordered\n");
   }
#endif

  return 1;
 }

double GuPendula(int i,int j, double *point,int nt, double kappa, double gamma, double r, MFErrorHandler e)
 {
  static char RoutineName[]={"GuPendula"};
  double w,Ic,In;
  double dwdT;
  double x,h;
  double W;
  double result;
  int n;
  int I,m;
  int Im,I0,Ip;
  double f,fm,f0,fp;

#ifdef MFNOCONFIDENCE
  if(i<0 || i>2*nt+1)
   {
    sprintf(MFPendulaMFErrorHandlerMsg,"i<0 || i>n-k, in G i=%d, n-k=%d\n",i,2*nt+1);
    MFSetError(e,12,RoutineName,MFPendulaMFErrorHandlerMsg,__LINE__,__FILE__);
    return 0.;
   }

  if(j<0 || j>2*nt+3)
   {
    sprintf(MFPendulaMFErrorHandlerMsg,"j<0 || j>n, in Guu i=%d, j=%d, n=%d\n",i,j,2*nt+3);
    MFSetError(e,12,RoutineName,MFPendulaMFErrorHandlerMsg,__LINE__,__FILE__);
    return 0.;
   }
#endif

  w=2*3.1415926/(10*point[2*nt]);
  dwdT=-10*w*w/(2*3.1415926);
  In=point[2*nt+1];
  Ic=3.1415926*kappa*point[2*nt+2];
  Ic=point[2*nt+2];

/*
 
   EVALUATES GU(U,L)
                                 I
      U IS:
             1     2     3        N-4  N-3   N-2  N-1 N
 
           P1(1),P2(1),P1(2),...,P1(M),P2(M), T , In, Ic
 
*/

  h=r*2.*3.1415926/nt;

  I=i/2;
  m=i%2;

  x=I*h;

  Im=I-1;
  I0=I;
  Ip=I+1;

  if(I0== 0)Im=nt-1;
  if(I0==nt-1)Ip=0;

  result=0.;
  if(i<2*nt)
   {
    result=0.;
    switch(m)
     {
      case 0:
       if(j==2*Ip)
        {
         fp=cos(w*(x+h)+point[2*Ip])+kappa;
         f0=0.;
         fm=0.;
         f=.25*(fp+2*f0+fm);
         result=w*w+gamma*w*h+f*h*h;
        }else if(j== 2*Ip+1)
        {
         fp=-kappa;
         f0=0.;
         fm=0.;
         f=.25*(fp+2*f0+fm);
         result=f*h*h;
        }else if(j== 2*I0)
        {
         fp=0.;
         f0=cos(w*x  +point[2*I0])+kappa;
         fm=0.;
         f=.25*(fp+2*f0+fm);
         result=-2.*w*w+f*h*h;
        }else if(j== 2*I0+1)
        {
         fp=0.;
         f0=-kappa;
         fm=0.;
         f=.25*(fp+2*f0+fm);
         result=f*h*h;
        }else if(j== 2*Im)
        {
         fp=0.;
         f0=0.;
         fm=cos(w*(x-h)+point[2*Im])+kappa;
         f=.25*(fp+2*f0+fm);
         result=w*w-gamma*w*h+f*h*h;
        }else if(j== 2*Im+1)
        {
         fp=0.;
         f0=0.;
         fm=-kappa;
         f=.25*(fp+2*f0+fm);
         result=f*h*h;
        }else if(j== 2*nt)
        {   /* w */
         fp=1.+(x+h)*cos(w*(x+h)+point[2*Ip]);
         f0=1.+ x   *cos(w*x  +point[2*I0]);
         fm=1.+(x-h)*cos(w*(x-h)+point[2*Im]);
         f=.25*(fp+2*f0+fm);
         result=2.*w*(point[2*Ip]-2.*point[2*I0]+point[2*Im])+gamma*(point[2*Ip]-point[2*Im])*h+f*h*h;
         result=result*dwdT;
        }else if(j== 2*nt+1)
        {   /* In */
         fp=-1.;
         f0=-1.;
         fm=-1.;
         f=.25*(fp+2*f0+fm);
         result=f*h*h;
        }else if(j== 2*nt+2)
        {   /* Ic */
         fp=-3.1415926*kappa;
         f0=-3.1415926*kappa;
         fm=-3.1415926*kappa;
         fp=-1.;
         f0=-1.;
         fm=-1.;
         f=.25*(fp+2*f0+fm);
         result=f*h*h;
        }else
         result=0.;
       break;
      case 1:
       if(j== 2*Ip)
        {
         fp=-kappa;
         f0=0.;
         fm=0.;
         f=.25*(fp+2*f0+fm);
         result=f*h*h;
        }else if(j== 2*Ip+1)
        {
         fp=cos(w*(x+h)+point[2*Ip+1])+kappa;
         f0=0.;
         fm=0.;
         f=.25*(fp+2*f0+fm);
         result=w*w+gamma*w*h+f*h*h;
        }else if(j== 2*I0)
        {
         fp=0.;
         f0=-kappa;
         fm=0.;
         f=.25*(fp+2*f0+fm);
         result=f*h*h;
        }else if(j== 2*I0+1)
        {
         fp=0.;
         f0=cos(w*(x  )+point[2*I0+1])+kappa;
         fm=0.;
         f=.25*(fp+2*f0+fm);
         result=-2.*w*w+f*h*h;
        }else if(j== 2*Im)
        {
         fp=0.;
         f0=0.;
         fm=-kappa;
         f=.25*(fp+2*f0+fm);
         result=f*h*h;
        }else if(j== 2*Im+1)
        {
         fp=0.;
         f0=0.;
         fm=cos(w*(x-h)+point[2*Im+1])+kappa;
         f=.25*(fp+2*f0+fm);
         result=w*w-gamma*w*h+f*h*h;
        }else if(j== 2*nt)
        {   /* w  */
         fp=1.+(x+h)*cos(w*(x+h)+point[2*Ip+1]);
         f0=1.+ x   *cos(w*(x  )+point[2*I0+1]);
         fm=1.+(x-h)*cos(w*(x-h)+point[2*Im+1]);
         f=.25*(fp+2*f0+fm);
         result=2*w*(point[2*Ip+1]-2.*point[2*I0+1]+point[2*Im+1])+gamma*(point[2*Ip+1]-point[2*Im+1])*h+f*h*h;
         result=result*dwdT;
        }else if(j== 2*nt+1)
        {   /* In */
         fp=-1.;
         f0=-1.;
         fm=-1.;
         f=.25*(fp+2*f0+fm);
         result=f*h*h;
        }else if(j== 2*nt+2)
        {   /* Ic */
         fp=1.;
         f0=1.;
         fm=1.;
         f=.25*(fp+2*f0+fm);
         result=f*h*h;
        }else
         result=0.;
       break;
     }
   }else{
    /* Phase Constraint */
    if(j==2*(nt-1))result=.5;
     else if(j==2*(nt-1)+1)result=.5;
     else result=0.;
   }

  return(result);
 }

double TGuuPendula(int i,int j,int k, double *point,int nt, double kappa, double gamma, double r, MFErrorHandler e)
 {
  static char RoutineName[]={"TGuuPendula"};
  double eps=.001;
  double t;
  double gm,gp;
  double test;
  double result;

  t=point[k];

  point[k]=t+eps/2;
  gp=TGuPendula(i,j,point,nt,kappa,gamma,r,e);

  point[k]=t-eps/2;
  gm=TGuPendula(i,j,point,nt,kappa,gamma,r,e);

  point[k]=t;

  result=(gp-gm)/eps;
  test=GuuPendula(i,j,k,point,nt,kappa,gamma,r,e);
  
  return result;
 }

double TGuPendula(int i,int j, double *point,int nt, double kappa, double gamma, double r, MFErrorHandler e)
 {
  static char RoutineName[]={"TGuPendula"};
  double eps=.001;
  double t;
  double gm,gp;
  double test;
  double result;

  t=point[j];

  point[j]=t+eps/2;
  gp=GPendula(i,point,nt,kappa,gamma,r,e);

  point[j]=t-eps/2;
  gm=GPendula(i,point,nt,kappa,gamma,r,e);

  point[j]=t;

  result=(gp-gm)/eps;
  test=GuPendula(i,j,point,nt,kappa,gamma,r,e);

#ifdef MFNOCONFIDENCE
  if(fabs(test-result)>eps)
   {printf("Large error in Gu, %d, %d, diff=%lf, exact=%lf\n",i,j,result,test);fflush(stdout);}
#endif
  
  return result;
 }

double GuuPendula(int i,int j,int k, double *point,int nt, double kappa, double gamma, double r, MFErrorHandler e)
 {
  static char RoutineName[]={"GuuPendula"};
  double w,Ic,In;
  double dwdT,ddwdTdw;
  double x,h;
  double W;
  double result;
  int n;
  int I,m;
  int Im,I0,Ip;
  double f,fm,f0,fp;


#ifdef MFNOCONFIDENCE
  if(i<0 || i>2*nt+1)
   {
    sprintf(MFPendulaMFErrorHandlerMsg,"i<0 || i>n-k, in G i=%d, n-k=%d\n",i,2*nt+1);
    MFSetError(e,12,RoutineName,MFPendulaMFErrorHandlerMsg,__LINE__,__FILE__);
    return 0.;
   }

  if(j<0 || j>2*nt+3)
   {
    sprintf(MFPendulaMFErrorHandlerMsg,"j<0 || j>n, in Guu i=%d, j=%d, k=%d, n=%d\n",i,j,k,2*nt+3);
    MFSetError(e,12,RoutineName,MFPendulaMFErrorHandlerMsg,__LINE__,__FILE__);
    return 0.;
   }

  if(k<0 || k>2*nt+3)
   {
    sprintf(MFPendulaMFErrorHandlerMsg,"k<0 || k>n, in Guu i=%d, j=%d, k=%d, n=%d\n",i,j,k,2*nt+3);
    MFSetError(e,12,RoutineName,MFPendulaMFErrorHandlerMsg,__LINE__,__FILE__);
    return 0.;
   }
#endif

  if(j>i+nBandsU && j!=2*nt && k!=2*nt)return 0.;
  if(j<i-nBandsL && j!=2*nt && k!=2*nt)return 0.;
  if(j!=k && j!=2*nt && k!=2*nt)return 0.;

  w=2*3.1415926/(10*point[2*nt]);
  dwdT=-10*w*w/(2*3.1415926);
  ddwdTdw=-2*10*w/(2*3.1415926);
  In=point[2*nt+1];
  Ic=3.1415926*kappa*point[2*nt+2];
  Ic=point[2*nt+2];

/*
 
   EVALUATES GUU(U,L)
                                 I
      U IS:
             1     2     3        N-4  N-3   N-2  N-1 N
 
           P1(1),P2(1),P1(2),...,P1(M),P2(M), T , In, Ic
 
*/

  h=r*2.*3.1415926/nt;

  I=i/2;
  m=i%2;

  x=I*h;

  Im=I-1;
  I0=I;
  Ip=I+1;

  if(I0== 0)Im=nt-1;
  if(I0==nt-1)Ip=0;

  result=0.;
  if(i<2*nt)
   {
    result=0.;
    switch(m)
     {
      case 0:
       if(j==2*Ip)
        {
         if(k==2*Ip)
          {
           fp=-sin(w*(x+h)+point[2*Ip]);
           f0=0.;
           fm=0.;
           f=.25*(fp+2*f0+fm);
           result=f*h*h;
          }else if(k==2*nt)
          {
           fp=-(x+h)*sin(w*(x+h)+point[2*Ip]);
           f0=0.;
           fm=0.;
           f=.25*(fp+2*f0+fm);
           result=2*w+gamma*h+f*h*h;
           result=result*dwdT;
          }else result=0.;
        }else if(j== 2*I0)
        {
         if(k==2*I0)
          {
           fp=0.;
           f0=-sin(w*x  +point[2*I0]);
           fm=0.;
           f=.25*(fp+2*f0+fm);
           result=f*h*h;
          }else if(k==2*nt)
          {
           fp=0.;
           f0=-x*sin(w*x  +point[2*I0]);
           fm=0.;
           f=.25*(fp+2*f0+fm);
           result=-4*w+f*h*h;
           result=result*dwdT;
          }else result=0.;
        }else if(j== 2*Im)
        {
         if(k==2*Im)
          {
           fp=0.;
           f0=0.;
           fm=-sin(w*(x-h)+point[2*Im]);
           f=.25*(fp+2*f0+fm);
           result=f*h*h;
          }else if(k==2*nt)
          {
           fp=0.;
           f0=0.;
           fm=-(x-h)*sin(w*(x-h)+point[2*Im]);
           f=.25*(fp+2*f0+fm);
           result=2*w-gamma*h+f*h*h;
           result=result*dwdT;
          }else result=0.;
        }else if(j== 2*nt)
        {   /* w */
         if(k==2*Ip)
          {
           fp=-(x+h)*sin(w*(x+h)+point[2*Ip]);
           f0=0.;
           fm=0.;
           f=.25*(fp+2*f0+fm);
           result=2.*w+gamma*h+f*h*h;
           result=result*dwdT;
          }else if(k==2*I0)
          {
           fp=0.;
           f0=-x*sin(w*x  +point[2*I0]);
           fm=0.;
           f=.25*(fp+2*f0+fm);
           result=-4*w+f*h*h;
           result=result*dwdT;
          }else if(k==2*Im)
          {
           fp=0.;
           f0=0.;
           fm=-(x-h)*sin(w*(x-h)+point[2*Im]);
           f=.25*(fp+2*f0+fm);
           result=2.*w-gamma*h+f*h*h;
           result=result*dwdT;
          }else if(k==2*nt)
          {
           fp=-(x+h)*(x+h)*sin(w*(x+h)+point[2*Ip]);
           f0=- x   * x   *sin(w*x  +point[2*I0]);
           fm=-(x-h)*(x-h)*sin(w*(x-h)+point[2*Im]);
           f=.25*(fp+2*f0+fm);
           result=2*(point[2*Ip]-2.*point[2*I0]+point[2*Im])+f*h*h;
           result=result*dwdT*dwdT;
           fp=1.+(x+h)*cos(w*(x+h)+point[2*Ip]);
           f0=1.+ x   *cos(w*x  +point[2*I0]);
           fm=1.+(x-h)*cos(w*(x-h)+point[2*Im]);
           f=.25*(fp+2*f0+fm);
           result+=dwdT*ddwdTdw*(2.*w*(point[2*Ip]-2.*point[2*I0]+point[2*Im])+gamma*(point[2*Ip]-point[2*Im])*h+f*h*h);
          }else result=0.;
        }else result=0.;
       break;
      case 1:
       if(j== 2*Ip+1)
        {
         if(k== 2*Ip+1)
          {
           fp=-sin(w*(x+h)+point[2*Ip+1]);
           f0=0.;
           fm=0.;
           f=.25*(fp+2*f0+fm);
           result=f*h*h;
          }else if(k== 2*nt)
          {
           fp=-(x+h)*sin(w*(x+h)+point[2*Ip+1]);
           f0=0.;
           fm=0.;
           f=.25*(fp+2*f0+fm);
           result=2*w+gamma*h+f*h*h;
           result=result*dwdT;
          }else result=0.;
        }else if(j== 2*I0+1)
        {
         if(k== 2*I0+1)
          {
           fp=0.;
           f0=-sin(w*(x  )+point[2*I0+1]);
           fm=0.;
           f=.25*(fp+2*f0+fm);
           result=f*h*h;
          }else if(k== 2*nt)
          {
           fp=0.;
           f0=-x*sin(w*(x  )+point[2*I0+1]);
           fm=0.;
           f=.25*(fp+2*f0+fm);
           result=-4*w+f*h*h;
           result=result*dwdT;
          }else result=0.;
        }else if(j== 2*Im+1)
        {
         if(k== 2*Im+1)
          {
           fp=0.;
           f0=0.;
           fm=-sin(w*(x-h)+point[2*Im+1]);
           f=.25*(fp+2*f0+fm);
           result=f*h*h;
          }else if(k== 2*nt)
          {
           fp=0.;
           f0=0.;
           fm=-(x-h)*sin(w*(x-h)+point[2*Im+1]);
           f=.25*(fp+2*f0+fm);
           result=2*w-gamma*h+f*h*h;
           result=result*dwdT;
          }else result=0.;
        }else if(j== 2*nt)
        {   /* w  */
         if(k== 2*Ip+1)
          {
           fp=-(x+h)*sin(w*(x+h)+point[2*Ip+1]);
           f0=0.;
           fm=0.;
           f=.25*(fp+2*f0+fm);
           result=2*w+gamma*h+f*h*h;
           result=result*dwdT;
          }else if(k== 2*I0+1)
          {
           fp=0.;
           f0=-x*sin(w*(x  )+point[2*I0+1]);
           fm=0.;
           f=.25*(fp+2*f0+fm);
           result=-4*w+f*h*h;
           result=result*dwdT;
          }else if(k== 2*Im+1)
          {
           fp=0.;
           f0=0.;
           fm=-(x-h)*sin(w*(x-h)+point[2*Im+1]);
           f=.25*(fp+2*f0+fm);
           result=2*w-gamma*h+f*h*h;
           result=result*dwdT;
          }else if(k== 2*nt)
          {
           fp=-(x+h)*(x+h)*sin(w*(x+h)+point[2*Ip+1]);
           f0=- x   * x   *sin(w*(x  )+point[2*I0+1]);
           fm=-(x-h)*(x-h)*sin(w*(x-h)+point[2*Im+1]);
           f=.25*(fp+2*f0+fm);
           result=2*(point[2*Ip+1]-2.*point[2*I0+1]+point[2*Im+1])+f*h*h;
           result=result*dwdT*dwdT;
           fp=1.+(x+h)*cos(w*(x+h)+point[2*Ip+1]);
           f0=1.+ x   *cos(w*x  +point[2*I0+1]);
           fm=1.+(x-h)*cos(w*(x-h)+point[2*Im+1]);
           f=.25*(fp+2*f0+fm);
           result+=dwdT*ddwdTdw*(2.*w*(point[2*Ip+1]-2.*point[2*I0+1]+point[2*Im+1])+gamma*(point[2*Ip+1]-point[2*Im+1])*h+f*h*h);
          }else result=0.;
        }else
         result=0.;
       break;
     }
   }else{
    /* Phase Constraint */
    result=0.;
   }

  return(result);
 }

void MFCurvaturePendula(int n,int k,double *u,double *tangentA,double *tangentB,double *a00, double *a01, double *a11, void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFCurvaturePendula"};
  static double *result;
  static double err;
  static int i,j,l;
  static int i0,i1;
  static int verbose=0;
  static int test=0;
  static int nBands,bandW;
  static int na,nb;
  static int ierr;
  static double scl;
  static int nt;
  static double kappa,gamma,r;

  nt=((struct MFPendulaData*)d)->nt;
  kappa=((struct MFPendulaData*)d)->kappa;
  gamma=((struct MFPendulaData*)d)->gamma;
  r=((struct MFPendulaData*)d)->R;

  nBands=nBandsL+nBandsU+1;

  for(j=0;j<n-k-nBorders;j++)
   {
    i0=j-nBandsU;
    if(i0<0)i0=0;
    i1=j+nBandsL;
    if(i1>n-k-nBorders-1)i1=n-k-nBorders-1;
    for(i=i0;i<i1+1;i++)
      A[i-j+nBands-1+ldA*j]=GU(i,j,u,nt,kappa,gamma,r,e);
   }

  for(i=0;i<n-k-nBorders;i++)
   {
    for(j=0;j<nBorders+k;j++)
      B[i+ldB*j]=GU(i,n-k-nBorders+j,u,nt,kappa,gamma,r,e);
   }

  for(i=0;i<nBorders;i++)
   {
    for(j=0;j<n-k-nBorders;j++)
      C[i+ldC*j]=GU(n-k-nBorders+i,j,u,nt,kappa,gamma,r,e);
   }
  for(j=0;j<n-k-nBorders;j++)
   {
    if(j<2*nt)
     {
      C[nBorders+0+ldC*j]=tangentA[j]/nt;
      C[nBorders+1+ldC*j]=tangentB[j]/nt;
     }else{
      C[nBorders+0+ldC*j]=tangentA[j];
      C[nBorders+1+ldC*j]=tangentB[j];
     }
   }

  for(i=0;i<nBorders;i++)
   {
    for(j=0;j<nBorders+k;j++)
      D[i+ldD*j]=GU(n-k-nBorders+i,n-k-nBorders+j,u,nt,kappa,gamma,r,e);
   }
  for(j=0;j<nBorders+k;j++)
   {
    if(n-k-nBorders+j<2*nt)
     {
      D[nBorders+0+ldD*j]=tangentA[n-k-nBorders+j]/nt;
      D[nBorders+1+ldD*j]=tangentB[n-k-nBorders+j]/nt;
     }else{
      D[nBorders+0+ldD*j]=tangentA[n-k-nBorders+j];
      D[nBorders+1+ldD*j]=tangentB[n-k-nBorders+j];
     }
   }

/* Calculate the size of the right hand side. */

/* Compute the tangent. */

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf(" System:\n\n");
    for(i=0;i<n-k-nBorders;i++)
     {
      printf("[");
      for(j=0;j<n-k-nBorders;j++)
       {
        if(j>i-nBandsL-1 && j<i+nBandsU+1)
         {
          printf("%7.3lf",A[i-j+nBands-1+ldA*j]);
         }else{
          printf("       ");
         }
       }
      printf(" | ");
      for(j=0;j<nBorders+k;j++)printf("%7.3lf",B[i+ldB*j]);
      printf("]");
      printf(" [%7.3lf]\n",R[i]);
     }
  
    printf("[");
    for(j=0;j<n-k-nBorders;j++)printf("--------");
    printf(" + ");
    for(j=0;j<nBorders+k;j++)printf("--------");
    printf("]");
    printf(" [-------]\n");

    for(i=0;i<nBorders+k;i++)
     {
      printf("[");
      for(j=0;j<n-k-nBorders;j++)printf("%7.3lf",C[i+ldC*j]);
      printf(" | ");
      for(j=0;j<nBorders+k;j++)printf("%7.3lf",D[i+ldC*j]);
      printf("]");
      printf(" [%7.3lf]\n",S[i]);
     }
   }
#endif

  na=n-k-nBorders;
  nb=nBorders+k;
  CALLDBOFA(A,&ldA,&na,&nBandsL,&nBandsU,B,&ldB,&nb,C,&ldC,D,&ldD,X,NULL,Pivots,&ierr);

#ifndef MFNOCONFIDENCE
  if(ierr!=0)
   {
    sprintf(MFPendulaMFErrorHandlerMsg," Problem with factor, zero on diagonal %d",ierr);
    MFSetError(e,12,RoutineName,MFPendulaMFErrorHandlerMsg,__LINE__,__FILE__);
    return;
   }
#endif

/* -GuuPhiSPhiT */

  for(i=0;i<n-k-nBorders;i++)
   {
    R[i]=0.;
    for(j=0;j<n;j++)
     {
      R[i]+=GUU(i,j,j,u,nt,kappa,gamma,r,e)*tangentA[j]*tangentA[j];
      for(l=j+1;l<n;l++)
        R[i]+=2*GUU(i,j,l,u,nt,kappa,gamma,r,e)*tangentA[j]*tangentA[l];
     }
   }

  for(i=0;i<nBorders;i++)
   {
    S[i]=0.;
    for(j=0;j<n;j++)
     {
      S[i]+=GUU(n-k-nBorders+i,j,j,u,nt,kappa,gamma,r,e)*tangentA[j]*tangentA[j];
      for(l=j+1;l<n;l++)
        S[i]+=2*GUU(n-k-nBorders+i,j,l,u,nt,kappa,gamma,r,e)*tangentA[j]*tangentA[l];
     }
   }
  for(i=nBorders;i<nBorders+k;i++)S[i]=0.;

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    err=0.;
    for(i=0;i<n-k-nBorders;i++)err+=R[i]*R[i]/nt;
    for(i=0;i<nBorders;i++)err+=S[i]*S[i]/nt;
    for(i=nBorders;i<nBorders+k;i++){if(i+n-k-nBorders==2*nt)err+=1.*(S[i]*S[i]);else err+=S[i]*S[i];}
    err=sqrt(err);
    printf("Norm of rhs for a00 is %lf\n",err);fflush(stdout);
   }
#endif

  ierr=0;
  CALLDBOSL(A,&ldA,&na,&nBandsL,&nBandsU,B,&ldB,&nb,C,&ldC,D,&ldD,X,Y,R,S,NULL,Pivots,&ierr);

/* Extract the result */

  for(i=0;i<n-k-nBorders;i++)
    a00[i]=R[i];
  for(i=0;i<nBorders+k;i++)
    a00[n-k-nBorders+i]=S[i];

  for(i=0;i<n-k-nBorders;i++)
   {
    R[i]=0.;
    for(j=0;j<n;j++)
     {
      R[i]+=GUU(i,j,j,u,nt,kappa,gamma,r,e)*tangentA[j]*tangentB[j];
      for(l=j+1;l<n;l++)
        R[i]+=GUU(i,j,l,u,nt,kappa,gamma,r,e)*(tangentA[j]*tangentB[l]+tangentA[l]*tangentB[j]);
     }
   }

  for(i=0;i<nBorders;i++)
   {
    S[i]=0.;
    for(j=0;j<n;j++)
     {
      S[i]+=GUU(n-k-nBorders+i,j,j,u,nt,kappa,gamma,r,e)*tangentA[j]*tangentB[j];
      for(l=j+1;l<n;l++)
        S[i]+=GUU(n-k-nBorders+i,j,l,u,nt,kappa,gamma,r,e)*(tangentA[j]*tangentB[l]+tangentA[l]*tangentB[j]);
     }
   }
  for(i=nBorders;i<nBorders+k;i++)S[i]=0.;

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    err=0.;
    for(i=0;i<n-k-nBorders;i++)err+=R[i]*R[i]/nt;
    for(i=0;i<nBorders;i++)err+=S[i]*S[i]/nt;
    for(i=nBorders;i<nBorders+k;i++){if(i+n-k-nBorders==2*nt)err+=1.*(S[i]*S[i]);else err+=S[i]*S[i];}
    err=sqrt(err);

    printf("Norm of rhs for a01 is %lf\n",e);fflush(stdout);
   }
#endif

  ierr=0;
  CALLDBOSL(A,&ldA,&na,&nBandsL,&nBandsU,B,&ldB,&nb,C,&ldC,D,&ldD,X,Y,R,S,NULL,Pivots,&ierr);

/* Extract the result */

  for(i=0;i<n-k-nBorders;i++)
    a01[i]=R[i];
  for(i=0;i<nBorders+k;i++)
    a01[n-k-nBorders+i]=S[i];

/* -GuuPhiSPhiT */

  for(i=0;i<n-k-nBorders;i++)
   {
    R[i]=0.;
    for(j=0;j<n;j++)
     {
      R[i]+=GUU(i,j,j,u,nt,kappa,gamma,r,e)*tangentB[j]*tangentB[j];
      for(l=j+1;l<n;l++)
       R[i]+=2*GUU(i,j,l,u,nt,kappa,gamma,r,e)*tangentB[j]*tangentB[l];
     }
   }

  for(i=0;i<nBorders;i++)
   {
    S[i]=0.;
    for(j=0;j<n;j++)
     {
      S[i]+=GUU(n-k-nBorders+i,j,j,u,nt,kappa,gamma,r,e)*tangentB[j]*tangentB[j];
      for(l=j+1;l<n;l++)
       S[i]+=2*GUU(n-k-nBorders+i,j,l,u,nt,kappa,gamma,r,e)*tangentB[j]*tangentB[l];
     }
   }
  for(i=nBorders;i<nBorders+k;i++)S[i]=0.;

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    err=0.;
    for(i=0;i<n-k-nBorders;i++)err+=R[i]*R[i]/nt;
    for(i=0;i<nBorders;i++)err+=S[i]*S[i]/nt;
    for(i=nBorders;i<nBorders+k;i++){if(i+n-k-nBorders==2*nt)err+=1.*(S[i]*S[i]);else err+=S[i]*S[i];}
    err=sqrt(err);

    printf("Norm of rhs for a11 is %lf\n",err);fflush(stdout);
   }
#endif

  ierr=0;
  CALLDBOSL(A,&ldA,&na,&nBandsL,&nBandsU,B,&ldB,&nb,C,&ldC,D,&ldD,X,Y,R,S,NULL,Pivots,&ierr);

/* Extract the result */

  for(i=0;i<n-k-nBorders;i++)
    a11[i]=R[i];
  for(i=0;i<nBorders+k;i++)
    a11[n-k-nBorders+i]=S[i];

  return;
 }

void MFWritePendulaData(FILE *fid,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFWritePendulaData"};
  struct MFPendulaData *data;
  int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s\n",RoutineName);fflush(stdout);}
#endif

  data=(struct MFPendulaData*)d;
  fprintf(fid,"%d %lf %lf %lf\n",data->nt,data->kappa,data->gamma,data->R);
  fflush(fid);

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("  nt=%d\n",data->nt);
    printf("  kappa=%lf\n",data->kappa);
    printf("  gamma=%lf\n",data->gamma);
    printf("  R=%lf\n",data->R);fflush(stdout);
   }
#endif

  return;
 }

MFImplicitMF MFReadPendula(FILE *fid, MFErrorHandler e)
 {
  static char RoutineName[]={"MFReadPendula"};
  MFImplicitMF pendula;
  int nt=0;
  double kappa=0.;
  double gamma=0.;
  double R=0.;
  int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s\n",RoutineName);fflush(stdout);}
#endif

  fscanf(fid,"%d %lf %lf %lf\n",&nt,&kappa,&gamma,&R);
  pendula=MFIMFCreatePendula(nt,kappa,gamma,1,e);

  return pendula;
 }

int MFPendulaProjectToSave(MFNVector u, double *x, void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFPendulaProjectToSave"};
  int i,n;
  struct MFPendulaData *data;

  data=(struct MFPendulaData*)d;

  n=2*(data->nt)+3;

  if(x==NULL)return n;
  for(i=0;i<n;i++)x[i]=MFNV_C(u,i,e);

  return 0;
 }

int MFPendulaProjectToDraw(MFNVector u, double *x, void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFPendulaProjectToDraw"};
  int i,n;
  struct MFPendulaData *data;

  data=(struct MFPendulaData*)d;

  n=2*(data->nt)+3;

  if(x==NULL)return 3;
  x[0]=MFNV_C(u,n-3,e);
  x[1]=MFNV_C(u,n-2,e);
  x[2]=MFNV_C(u,n-1,e);

  return 0;
 }

int MFPendulaProjectForBB(MFNVector u, double *x, void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFPendulaProjectForBB"};
  int i,n;
  struct MFPendulaData *data;

  data=(struct MFPendulaData*)d;

  n=2*(data->nt)+3;

  if(x==NULL)return 3;
  x[0]=MFNV_C(u,n-3,e);
  x[1]=MFNV_C(u,n-2,e);
  x[2]=MFNV_C(u,n-1,e);

  return 0;
 }

#ifdef __cplusplus
}
#endif
