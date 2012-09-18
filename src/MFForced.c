/* 
 *  PROGRAM NAME:  multifario
 *
 *  (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
 *  CORPORATION 12/1/2005.  ALL RIGHTS RESERVED.
 *
 *  Please refer to the LICENSE file in the top directory
 *
 *      author: Mike Henderson mhender@watson.ibm.com
 *      date:   May 5, 2005
 */

static char *id="@(#) $Id: MFForced.c,v 1.5 2011/07/21 17:42:46 mhender Exp $";

#include <MFImplicitMF.h>
#include <MFNKMatrix.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <MFFortran.h>
#include <MFNSpace.h>
#include <MFNVector.h>
#include <MFErrorHandler.h>

#define TWOPI 6.2831853071795862320

#define GU GuForcedOscillator
#define GUU GuuForcedOscillator

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))

#ifdef __cplusplus
 extern "C" {
#endif

static char MFForcedOscillatorMFErrorHandlerMsg[256]="";

extern double MFEpsilon;

/* Space for LU decompositions */

double *A,*B,*C,*D,*R,*S,*X,*Y;
int   *Pivots;

int nBorders,nBandsU,nBandsL;
int ldA,ldB,ldC,ldD;

MFNSpace MFCreateForcedOscillatorNSpace(MFErrorHandler);

double *MFKV_CStar(MFKVector,MFErrorHandler);
double *MFNV_CStar(MFNVector,MFErrorHandler);
double *MFNKM_CStar(MFNKMatrix,MFErrorHandler);

double TGuForcedOscillator(int,int,double*,int,double,MFErrorHandler);
double TGuuForcedOscillator(int,int,int,double*,int,double,MFErrorHandler);
double FGuuForcedOscillator(int,int,int,double*,int,double,MFErrorHandler);

void MFFreeForcedOscillatorData(void*,MFErrorHandler);
int MFProjectForcedOscillator(int,int,MFNVector,MFNKMatrix,MFNVector,void*,int*,MFErrorHandler);
int MFTangentForcedOscillator(int,int,MFNVector,MFNKMatrix,void*,MFErrorHandler);
int MFTangentForcedOscillatorWithGuess(int,int,MFNVector,MFNKMatrix,MFNKMatrix,void*,MFErrorHandler);
double MFScaleForcedOscillator(int,int,MFNVector,MFNKMatrix,void*,MFErrorHandler);
void MFForcedOscillatorInitializeSolver(int,MFErrorHandler);
double GForcedOscillator(int,double*,int,double,MFErrorHandler);
double GuForcedOscillator(int,int,double*,int,double,MFErrorHandler);
double GuuForcedOscillator(int,int,int,double*,int,double,MFErrorHandler);
int MFTangentForcedOscillatorSingle(int,int,int,double*,double*,double*,double*,void*,MFErrorHandler);
void MFCurvatureForcedOscillator(int,int,double*,double*,double*,double*,double*,double*,void*,MFErrorHandler);
void MFWriteForcedOscillatorData(FILE*,void*,MFErrorHandler);
MFImplicitMF MFReadForcedOscillator(FILE*,MFErrorHandler);

int MFForcedOscillatorProjectToSave(MFNVector,double*,void*,MFErrorHandler);
int MFForcedOscillatorProjectToDraw(MFNVector,double*,void*,MFErrorHandler);
int MFForcedOscillatorProjectForBB(MFNVector,double*,void*,MFErrorHandler);

MFNVector MFNVectorFactory(MFImplicitMF,MFErrorHandler);
MFNKMatrix MFNKMatrixFactory(MFImplicitMF,MFErrorHandler);

struct MFForcedOscillatorData
 {
  int nt;
  double R;
  MFNSpace space;
 };

MFImplicitMF MFIMFCreateForcedOscillator(int nt,int windingno,MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFCreateForcedOscillator"};
  MFImplicitMF forced;
  int i;
  MFNSpace space;
  struct MFForcedOscillatorData *data;

  forced=MFIMFCreateBaseClass(nt+2,2,"ForcedOscillator",e);

  space=MFCreateForcedOscillatorNSpace(e);
  MFIMFSetSpace(forced,space,e);
  MFFreeNSpace(space,e);

  data=(struct MFForcedOscillatorData*)malloc(sizeof(struct MFForcedOscillatorData));
#ifndef MFNPOSAFETYNET
  if(data==NULL)
   {
    sprintf(MFForcedOscillatorMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct MFForcedOscillatorData));
    MFSetError(e,12,RoutineName,MFForcedOscillatorMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  data->nt=nt;
  data->R=(double)windingno;
  data->space=space;MFRefNSpace(space,e);

  MFIMFSetData(forced,(void*)data,e);
  MFIMFSetFreeData(forced,MFFreeForcedOscillatorData,e);
  MFIMFSetSpace(forced,space,e);
  MFIMFSetProject(forced,MFProjectForcedOscillator,e);
  MFIMFSetTangent(forced,MFTangentForcedOscillator,e);
  MFIMFSetTangentWithGuess(forced,MFTangentForcedOscillatorWithGuess,e);
  MFIMFSetScale(forced,MFScaleForcedOscillator,e);
  MFIMFSetWriteData(forced,MFWriteForcedOscillatorData,e);
  MFIMFSetProjectForSave(forced,MFForcedOscillatorProjectToSave,e);
  MFIMFSetProjectForDraw(forced,MFForcedOscillatorProjectToDraw,e);
  MFIMFSetProjectForBB(forced,MFForcedOscillatorProjectForBB,e);

  MFIMFSetVectorFactory(forced,MFNVectorFactory,e);
  MFIMFSetMatrixFactory(forced,MFNKMatrixFactory,e);

  MFForcedOscillatorInitializeSolver(nt,e);

  return forced;
 }

void MFFreeForcedOscillatorData(void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFFreeForcedOscillatorData"};
  MFFreeNSpace(((struct MFForcedOscillatorData*)d)->space,e);
  free((struct MFForcedOscillatorData*)d);
  return;
 }

extern double MFEpsilon;

double MFScaleForcedOscillator(int n,int k,MFNVector vu,MFNKMatrix mPhi, void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFScaleForcedOscillator"};
  double *a00,*a01,*a11;
  double A,B,C,l;
  double r;
  int i,nt;
  double *u,*Phi;
  int verbose;

  verbose=0;

  Phi=MFNKM_CStar(mPhi,e);
  u=MFNV_CStar(vu,e);

  nt=n-2;

#ifdef MFALLOWVERBOSE
    if(verbose){printf("  MFScaleForcedOscillator \n");fflush(stdout);}
#endif

    a00=(double*)malloc(n*sizeof(double));
#ifndef MFNPOSAFETYNET
    if(a00==NULL)
     {
      sprintf(MFForcedOscillatorMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(double));
      MFSetError(e,12,RoutineName,MFForcedOscillatorMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return -1.;
     }
#endif

    a01=(double*)malloc(n*sizeof(double));
#ifndef MFNPOSAFETYNET
    if(a01==NULL)
     {
      sprintf(MFForcedOscillatorMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(double));
      MFSetError(e,12,RoutineName,MFForcedOscillatorMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return -1.;
     }
#endif

    a11=(double*)malloc(n*sizeof(double));
#ifndef MFNPOSAFETYNET
    if(a11==NULL)
     {
      sprintf(MFForcedOscillatorMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(double));
      MFSetError(e,12,RoutineName,MFForcedOscillatorMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return -1.;
     }
#endif

  MFCurvatureForcedOscillator(n,k,u,Phi,Phi+n,a00,a01,a11,d,e);

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
  r=2*sqrt(2*MFEpsilon/l);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s r=%lf\n",RoutineName,r);fflush(stdout);}
#endif

  if(r>.5)r=.5;

  printf("%s, r=%lf\n",RoutineName,r);

  return r;
 }

void MFForcedOscillatorInitializeSolver(int nt, MFErrorHandler e)
 {
  static char RoutineName[]={"MFForcedOscillatorInitializeSolver"};
  int nEqs,n;
  int k;
  int verbose=0;

  nEqs=1;
  k=2;
  n=nEqs*nt+2;
  nBorders=nEqs;
  nBandsU=2*nEqs-1;
  nBandsL=2*nEqs-1;

  ldA=2*nBandsL+nBandsU+1;
  ldB=n-k-nBorders;
  ldC=nBorders+k;
  ldD=nBorders+k;

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("%s\n",RoutineName);fflush(stdout);
    printf("   There are %d mesh point\n",nt);fflush(stdout);
    printf("             %d equations\n",nEqs);fflush(stdout);
    printf("             %d parameters\n",k);fflush(stdout);
    printf("         n=  %d \n",n);fflush(stdout);
    printf("    nBorders=%d \n",nBorders);fflush(stdout);
    printf("    nBandsU =%d \n",nBandsU);fflush(stdout);
    printf("    nBandsL =%d \n",nBandsL);fflush(stdout);
    printf("    ldA     =%d \n",ldA);fflush(stdout);
    printf("    ldB     =%d \n",ldB);fflush(stdout);
    printf("    ldC     =%d \n",ldC);fflush(stdout);
    printf("    ldD     =%d \n",ldD);fflush(stdout);
   }
#endif

  A=(double*)malloc(ldA*(n-k-nBorders)*sizeof(double));
#ifndef MFNPOSAFETYNET
  if(A==NULL)
   {
    sprintf(MFForcedOscillatorMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",ldA*(n-k-nBorders)*sizeof(double));
    MFSetError(e,12,RoutineName,MFForcedOscillatorMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  B=(double*)malloc(ldB*(nBorders+k)*sizeof(double));
#ifndef MFNPOSAFETYNET
  if(B==NULL)
   {
    sprintf(MFForcedOscillatorMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",ldB*(nBorders+k)*sizeof(double));
    MFSetError(e,12,RoutineName,MFForcedOscillatorMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  C=(double*)malloc(ldC*(n-k-nBorders)*sizeof(double));
#ifndef MFNPOSAFETYNET
  if(C==NULL)
   {
    sprintf(MFForcedOscillatorMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",ldC*(n-k-nBorders)*sizeof(double));
    MFSetError(e,12,RoutineName,MFForcedOscillatorMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  D=(double*)malloc(ldD*(nBorders+k)*sizeof(double));
#ifndef MFNPOSAFETYNET
  if(D==NULL)
   {
    sprintf(MFForcedOscillatorMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",ldD*(nBorders+k)*sizeof(double));
    MFSetError(e,12,RoutineName,MFForcedOscillatorMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  R=(double*)malloc((n-k-nBorders)*sizeof(double));
#ifndef MFNPOSAFETYNET
  if(R==NULL)
   {
    sprintf(MFForcedOscillatorMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",(n-k-nBorders)*sizeof(double));
    MFSetError(e,12,RoutineName,MFForcedOscillatorMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  S=(double*)malloc((nBorders+k)*sizeof(double));
#ifndef MFNPOSAFETYNET
  if(S==NULL)
   {
    sprintf(MFForcedOscillatorMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",(nBorders+k)*sizeof(double));
    MFSetError(e,12,RoutineName,MFForcedOscillatorMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  X=(double*)malloc((nBorders+k+nBandsL+1)*(nBorders+k+nBandsL+1)*sizeof(double));
#ifndef MFNPOSAFETYNET
  if(X==NULL)
   {
    sprintf(MFForcedOscillatorMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",(nBorders+k+nBandsL+1)*(nBorders+k+nBandsL+1)*sizeof(double));
    MFSetError(e,12,RoutineName,MFForcedOscillatorMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  Y=(double*)malloc((nBorders+k+nBandsL+1)*sizeof(double));
#ifndef MFNPOSAFETYNET
  if(X==NULL)
   {
    sprintf(MFForcedOscillatorMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",(nBorders+k+nBandsL+1)*sizeof(double));
    MFSetError(e,12,RoutineName,MFForcedOscillatorMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  Pivots=(int*)malloc(n*sizeof(int));
#ifndef MFNPOSAFETYNET
  if(Pivots==NULL)
   {
    sprintf(MFForcedOscillatorMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(int));
    MFSetError(e,12,RoutineName,MFForcedOscillatorMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  return;
 }

static int firsttime=1;
static FILE *ofid=NULL;
static int ocnt=0;

int MFProjectForcedOscillator(int n,int k,MFNVector vu0,MFNKMatrix mPhi,MFNVector vu,void *d,int *index, MFErrorHandler eh)
 {
  static char RoutineName[]={"MFProjectForcedOscillator"};
  static double eps=.001;
  static double e;
  static double emax;
  static double error;
  static double delta;
  static double factor;
  static int i,j;
  static int i0,i1;
  static int itimes;
  static int verbose,test;
  static int nBands,bandW;
  static int ierr;
  static double *tangentA;
  static double *tangentB;
  double *u,*u0,*Phi;

  static int nt;
  static int na,nb;
  static double r;
  int ip,im;

  u0=MFNV_CStar(vu0,eh);
  Phi=MFNKM_CStar(mPhi,eh);
  u=MFNV_CStar(vu,eh);

  nt=((struct MFForcedOscillatorData*)d)->nt;
  r=((struct MFForcedOscillatorData*)d)->R;

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
    e=0.;
    for(i=0;i<nt;i++)
     {
      ip=i+1;if(ip>=nt)ip=0;
      e+=(tangentA[ip]+tangentA[i])*(tangentA[ip]+tangentA[i])/4/nt;
      e+=(tangentA[ip]-tangentA[i])*(tangentA[ip]-tangentA[i])*nt;
     }
    for(i=nt;i<n;i++)
      e+=tangentA[i]*tangentA[i];
    printf("direction A.direction A: %e\n",e);

    e=0.;
    for(i=0;i<nt;i++)
     {
      ip=i+1;if(ip>=nt)ip=0;
      e+=(tangentA[ip]+tangentA[i])*(tangentB[ip]+tangentB[i])/4/nt;
      e+=(tangentA[ip]-tangentA[i])*(tangentB[ip]-tangentB[i])*nt;
     }
    for(i=nt;i<n;i++)
      e+=tangentA[i]*tangentB[i];
    printf("direction A.direction B: %e\n",e);

    e=0.;
    for(i=0;i<nt;i++)
     {
      ip=i+1;if(ip>=nt)ip=0;
      e+=(tangentB[ip]+tangentB[i])*(tangentB[ip]+tangentB[i])/4/nt;
      e+=(tangentB[ip]-tangentB[i])*(tangentB[ip]-tangentB[i])*nt;
     }
    for(i=nt;i<n;i++)
      e+=tangentB[i]*tangentB[i];
    printf("direction B.direction B: %e\n",e);
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
        A[i-j+nBands-1+ldA*j]=GU(i,j,u,nt,r,eh);
     }

    for(i=0;i<n-k-nBorders;i++)
     {
      for(j=0;j<nBorders+k;j++)
        B[i+ldB*j]=GU(i,n-k-nBorders+j,u,nt,r,eh);
     }

    for(i=0;i<nBorders;i++)
     {
      for(j=0;j<n-k-nBorders;j++)
        C[i+ldC*j]=GU(n-k-nBorders+i,j,u,nt,r,eh);
     }
    for(j=0;j<n-k-nBorders;j++)
     {
      i=j;
      im=i-1;if(im<0)im=nt-1;
      ip=i+1;if(ip>=nt)ip=0;
      if(j<nt)
       {
        C[nBorders+0+ldC*j] =(tangentA[ip]+2*tangentA[i]+tangentA[im])/4/nt;
        C[nBorders+0+ldC*j]-=(tangentA[ip]-2*tangentA[i]+tangentA[im])*nt;
        C[nBorders+1+ldC*j] =(tangentB[ip]+2*tangentB[i]+tangentB[im])/4/nt;
        C[nBorders+1+ldC*j]-=(tangentB[ip]-2*tangentB[i]+tangentB[im])*nt;
       }else{
        C[nBorders+0+ldC*j]=tangentA[j];
        C[nBorders+1+ldC*j]=tangentB[j];
       }
     }

    for(i=0;i<nBorders;i++)
     {
      for(j=0;j<nBorders+k;j++)
        D[i+ldD*j]=GU(n-k-nBorders+i,n-k-nBorders+j,u,nt,r,eh);
     }
    for(j=0;j<nBorders+k;j++)
     {
      i=n-k-nBorders+j;
      im=i-1;if(im<0)im=nt-1;
      ip=i+1;if(ip>=nt)ip=0;
      if(n-k-nBorders+j<nt)
       {
        D[nBorders+0+ldC*j] =(tangentA[ip]+2*tangentA[i]+tangentA[im])/4/nt;
        D[nBorders+0+ldC*j]-=(tangentA[ip]-2*tangentA[i]+tangentA[im])*nt;
        D[nBorders+1+ldC*j] =(tangentB[ip]+2*tangentB[i]+tangentB[im])/4/nt;
        D[nBorders+1+ldC*j]-=(tangentB[ip]-2*tangentB[i]+tangentB[im])*nt;
       }else{
        D[nBorders+0+ldD*j]=tangentA[n-k-nBorders+j];
        D[nBorders+1+ldD*j]=tangentB[n-k-nBorders+j];
       }
     }

    for(i=0;i<n-k-nBorders;i++)
      R[i]=-GForcedOscillator(i,u,nt,r,eh);

    for(i=0;i<nBorders;i++)
      S[i]=-GForcedOscillator(n-k-nBorders+i,u,nt,r,eh);

    S[nBorders+0]=0.;
    S[nBorders+1]=0.;
    for(i=0;i<n;i++)
     {
      ip=i+1;if(ip>=nt)ip=0;
      if(i<nt)
       {
        S[nBorders+0]+=(tangentA[ip]+tangentA[i])*(u0[ip]-u[ip]+u0[i]-u[i])/nt/4;
        S[nBorders+0]+=(tangentA[ip]-tangentA[i])*(u0[ip]-u[ip]-u0[i]+u[i])*nt;
        S[nBorders+1]+=(tangentB[ip]+tangentB[i])*(u0[ip]-u[ip]+u0[i]-u[i])/nt/4;
        S[nBorders+1]+=(tangentB[ip]-tangentB[i])*(u0[ip]-u[ip]-u0[i]+u[i])*nt;
       }else{
        S[nBorders+0]+=tangentA[i]*(u0[i]-u[i]);
        S[nBorders+1]+=tangentB[i]*(u0[i]-u[i]);
       }
     }

/* Calculate the size of the right hand side. */

    error=0.;
    for(i=0;i<n-k-nBorders;i++)
     {
      if(i<nt-1)
       {
        error+=(R[i+1]+R[i])*(R[i+1]+R[i])/4/nt;
        error+=(R[i+1]-R[i])*(R[i+1]-R[i])*nt;
       }else if(i<nt)
       {
        error+=(S[0]+R[i])*(S[0]+R[i])/4/nt;
        error+=(S[0]-R[i])*(S[0]-R[i])*nt;
       }else{
        error+=R[i]*R[i];
       }
     }
    for(i=0;i<nBorders+k;i++)
     {
      if(i+n-k-nBorders<nt-1)
       {
        error+=(S[i+1]+S[i])*(S[i+1]+S[i])/4/nt;
        error+=(S[i+1]-S[i])*(S[i+1]-S[i])*nt;
       }else if(i+n-k-nBorders<nt)
       {
        error+=(R[0]+S[i])*(R[0]+S[i])/4/nt;
        error+=(R[0]-S[i])*(R[0]-S[i])*nt;
       }else{
        error+=S[i]*S[i];
       }
     }
    error=sqrt(error);

#ifdef MFALLOWVERBOSE
    if(verbose)printf("\n                   Size of rhs: %e\n",error);
#endif

    na=n-k-nBorders;
    nb=nBorders+k;
    CALLDBOFA(A,&ldA,&na,&nBandsL,&nBandsU,B,&ldB,&nb,C,&ldC,D,&ldD,X,NULL,Pivots,&ierr);
    if(ierr!=0)
     {
      printf(" Problem with factor, zero on diagonal %d\n",ierr);
      return 0;
     }
    ierr=0;
    CALLDBOSL(A,&ldA,&na,&nBandsL,&nBandsU,B,&ldB,&nb,C,&ldC,D,&ldD,X,Y,R,S,NULL,Pivots,&ierr);

/* Is GForcedOscillator(u+eps*delta)-GForcedOscillator(u) ~ -eps*GForcedOscillator(u)  */
/*    P(u+eps*delta)-P(u) ~ -eps*P(u)  */
/*    N(u+eps*delta)-N(u) ~ -eps*N(u)? */

/* Calculate the size of the Correction. */

    delta=0.;
    for(i=0;i<n-k-nBorders;i++)
     {
      if(i<nt)
       {
        delta+=R[i]*R[i]/nt;
       }else{
        delta+=R[i]*R[i];
       }
     }
    for(i=0;i<nBorders+k;i++)
     {
      if(i+n-k-nBorders<nt)
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
/*  if(delta>MFEpsilon)factor=MFEpsilon/delta;*/
/*  if(delta>MFEpsilon)return 0;*/

    for(i=0;i<n-k-nBorders;i++)
      u[i]+=factor*R[i];
    for(i=0;i<nBorders+k;i++)
      u[n-k-nBorders+i]+=factor*S[i];

#ifdef MFALLOWVERBOSE
    if(verbose)printf(" %d error %e size of correction %e (%le)\n",itimes,error,delta*factor,delta);
#endif

    itimes++;
    if(itimes>6&&!firsttime || itimes>100)
     {
      printf(" Too many refine iterations %d\n",itimes);fflush(stdout);
      *index=0;
      return 0;
     }
   }

#ifdef MFALLOWVERBOSE
  if(verbose)printf("\n");
#endif

  *index=0;

#ifdef WRITESOLUTIONS
  if(ofid==NULL)
   {
    ofid=fopen("junk.dx","w");
    fprintf(ofid,"object \"lines\" class gridconnections counts %d\n\n",nt);
   }
  fprintf(ofid,"#  e=%lf B=%lf\n",u[nt],u[nt+1]);
  fprintf(ofid,"object \"p%d\" class array type float rank 1 shape 2 items %d data follows\n",ocnt,nt);
  for(i=0;i<nt;i++)
    fprintf(ofid,"    %lf   %lf\n",r*3.1415926*i/nt,u[i]);
  fprintf(ofid,"\n");
 
  fprintf(ofid,"object \"obj%d\" class field\n",ocnt);
  fprintf(ofid,"component \"positions\" value \"p%d\"\n",ocnt);
  fprintf(ofid,"component \"connections\" value \"lines\"\n\n");fflush(ofid);
#endif

  ocnt++;
  firsttime=0;
  return 1;
 }

void ForcedWrite(MFNVector U, MFErrorHandler e)
 {
#ifdef WRITESOLUTIONS
  double *u;
  int i,nt;

  nt=MFNV_NC(U,e)-2;

  u=MFNV_CStar(U,e);
  if(ofid==NULL)
   {
    ofid=fopen("junk.dx","w");
    fprintf(ofid,"object \"lines\" class gridconnections counts %d\n\n",nt);
   }
  fprintf(ofid,"#  e=%lf B=%lf\n",u[nt],u[nt+1]);
  fprintf(ofid,"object \"p%d\" class array type float rank 1 shape 2 items %d data follows\n",ocnt,nt);
  for(i=0;i<nt;i++)
    fprintf(ofid,"    %lf   %lf\n",4*3.1415926*i/nt,u[i]);
  fprintf(ofid,"\n");

  fprintf(ofid,"object \"obj%d\" class field\n",ocnt);
  fprintf(ofid,"component \"positions\" value \"p%d\"\n",ocnt);
  fprintf(ofid,"component \"connections\" value \"lines\"\n\n");fflush(ofid);

  ocnt++;
#endif
  return;
 }

void ForcedCloser(MFErrorHandler e)
 {
#ifdef WRITESOLUTIONS
  int i;

  fprintf(ofid,"object \"default\" class group\n");
  for(i=0;i<ocnt;i++)
    fprintf(ofid," member %d \"obj%d\"\n",i,i);
  fprintf(ofid,"end\n");

  fclose(ofid);
  ofid=NULL;
#endif
  return;
 }

int MFTangentForcedOscillator(int n,int k,MFNVector vu,MFNKMatrix mPhi,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFTangentForcedOscillator"};
  int i;
  double *Phi0;
  double *u,*Phi;

  Phi=MFNKM_CStar(mPhi,e);
  u=MFNV_CStar(vu,e);

  Phi0=(double*)malloc(n*k*sizeof(double));
#ifndef MFNPOSAFETYNET
  if(Phi0==NULL)
   {
    sprintf(MFForcedOscillatorMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*k*sizeof(double));
    MFSetError(e,12,RoutineName,MFForcedOscillatorMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0;
   }
#endif

  for(i=0;i<n*k;i++)Phi0[i]=0.;
  Phi0[n-2]=1.;
  Phi0[n-1+n]=1.;
   
  MFTangentForcedOscillatorSingle(n,k,0,u,Phi0,Phi0+n,Phi,d,e);
  MFTangentForcedOscillatorSingle(n,k,1,u,Phi ,Phi0+n,Phi,d,e);
  free(Phi0);
  MFGramSchmidt(((struct MFForcedOscillatorData*)d)->space,mPhi,e);

  if(0)
   {
    double e;
    double *tangentA,*tangentB;
    int i,ip;
    int nt;

    printf("%s, n=%d\n",RoutineName,n);fflush(stdout);
    nt=n-2;
    tangentA=Phi;
    tangentB=Phi+n;
    e=0.;
    for(i=0;i<nt;i++)
     {
      ip=i+1;if(ip>=nt)ip=0;
      e+=(tangentA[ip]+tangentA[i])*(tangentA[ip]+tangentA[i])/4/nt;
      e+=(tangentA[ip]-tangentA[i])*(tangentA[ip]-tangentA[i])*nt;
     }
    for(i=nt;i<n;i++)
      e+=tangentA[i]*tangentA[i];
    printf("direction A.direction A: %e\n",e);

    e=0.;
    for(i=0;i<nt;i++)
     {
      ip=i+1;if(ip>=nt)ip=0;
      e+=(tangentA[ip]+tangentA[i])*(tangentB[ip]+tangentB[i])/4/nt;
      e+=(tangentA[ip]-tangentA[i])*(tangentB[ip]-tangentB[i])*nt;
     }
    for(i=nt;i<n;i++)
      e+=tangentA[i]*tangentB[i];
    printf("direction A.direction B: %e\n",e);

    e=0.;
    for(i=0;i<nt;i++)
     {
      ip=i+1;if(ip>=nt)ip=0;
      e+=(tangentB[ip]+tangentB[i])*(tangentB[ip]+tangentB[i])/4/nt;
      e+=(tangentB[ip]-tangentB[i])*(tangentB[ip]-tangentB[i])*nt;
     }
    for(i=nt;i<n;i++)
      e+=tangentB[i]*tangentB[i];
    printf("direction B.direction B: %e\n",e);
   }

  return 1;
 }

int MFTangentForcedOscillatorWithGuess(int n,int k,MFNVector vu,MFNKMatrix mPhi0,MFNKMatrix mPhi,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFTangentForcedOscillatorWithGuess"};
  double *u,*Phi,*Phi0;

  Phi=MFNKM_CStar(mPhi,e);
  Phi0=MFNKM_CStar(mPhi0,e);
  u=MFNV_CStar(vu,e);

  MFTangentForcedOscillatorSingle(n,k,0,u,Phi0,Phi0+n,Phi,d,e);
  MFTangentForcedOscillatorSingle(n,k,1,u,Phi ,Phi0+n,Phi,d,e);
  MFGramSchmidt(((struct MFForcedOscillatorData*)d)->space,mPhi,e);

  return 1;
 }

int MFTangentForcedOscillatorSingle(int n,int k,int t, double *u,double *tangentA,double *tangentB,double *Phi, void *d, MFErrorHandler eh)
 {
  static double *result;
  static int i,j;
  static int i0,i1;
  static int verbose,test;
  static int nBands,bandW;
  static int na,nb;
  static int ierr;
  static double s;
  static double e;
  static int nt;
  static double r;
  int ip,im;

  nt=((struct MFForcedOscillatorData*)d)->nt;
  r=((struct MFForcedOscillatorData*)d)->R;

  verbose=0;
  test=0;

#ifdef MFALLOWVERBOSE
  if(verbose||test)printf("surfaceTangentBordered direction %d\n",t);
#endif

#ifdef MFALLOWVERBOSE
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
    e=0.;
    for(i=0;i<nt;i++)
     {
      ip=i+1;if(ip>=nt)ip=0;
      e+=(tangentA[ip]+tangentA[i])*(tangentA[ip]+tangentA[i])/4/nt;
      e+=(tangentA[ip]-tangentA[i])*(tangentA[ip]-tangentA[i])*nt;
     }
    for(i=nt;i<n;i++)
      e+=tangentA[i]*tangentA[i];
    printf("direction A.direction A: %e\n",e);

    e=0.;
    for(i=0;i<nt;i++)
     {
      e+=(tangentA[ip]+tangentA[i])*(tangentB[ip]+tangentB[i])/4/nt;
      e+=(tangentA[ip]-tangentA[i])*(tangentB[ip]-tangentB[i])*nt;
     }
    for(i=nt;i<n;i++)
      e+=tangentA[i]*tangentB[i];
    printf("direction A.direction B: %e\n",e);

    e=0.;
    for(i=0;i<nt;i++)
     {
      e+=(tangentB[ip]+tangentB[i])*(tangentB[ip]+tangentB[i])/4/nt;
      e+=(tangentB[ip]-tangentB[i])*(tangentB[ip]-tangentB[i])*nt;
     }
    for(i=nt;i<n;i++)
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
      A[i-j+nBands-1+ldA*j]=GU(i,j,u,nt,r,eh);
   }

  for(i=0;i<n-k-nBorders;i++)
   {
    for(j=0;j<nBorders+k;j++)
      B[i+ldB*j]=GU(i,n-k-nBorders+j,u,nt,r,eh);
   }

  for(i=0;i<nBorders;i++)
   {
    for(j=0;j<n-k-nBorders;j++)
      C[i+ldC*j]=GU(n-k-nBorders+i,j,u,nt,r,eh);
   }
  for(j=0;j<n-k-nBorders;j++)
   {
    i=j;
    im=i-1;if(im<0)im=nt-1;
    ip=i+1;if(ip>=nt)ip=0;
    if(j<nt)
     {
      C[nBorders+0+ldC*j] =(tangentA[ip]+2*tangentA[i]+tangentA[im])/4/nt;
      C[nBorders+0+ldC*j]-=(tangentA[ip]-2*tangentA[i]+tangentA[im])*nt;
      C[nBorders+1+ldC*j] =(tangentB[ip]+2*tangentB[i]+tangentB[im])/4/nt;
      C[nBorders+1+ldC*j]-=(tangentB[ip]-2*tangentB[i]+tangentB[im])*nt;
     }else{
      C[nBorders+0+ldC*j]=tangentA[j];
      C[nBorders+1+ldC*j]=tangentB[j];
     }
   }

  for(i=0;i<nBorders;i++)
   {
    for(j=0;j<nBorders+k;j++)
      D[i+ldD*j]=GU(n-k-nBorders+i,n-k-nBorders+j,u,nt,r,eh);
   }
  for(j=0;j<nBorders+k;j++)
   {
    i=n-k-nBorders+j;
    im=i-1;if(im<0)im=nt-1;
    ip=i+1;if(ip>=nt)ip=0;
    if(n-k-nBorders+j<nt)
     {
      D[nBorders+0+ldC*j] =(tangentA[ip]+2*tangentA[i]+tangentA[im])/4/nt;
      D[nBorders+0+ldC*j]-=(tangentA[ip]-2*tangentA[i]+tangentA[im])*nt;
      D[nBorders+1+ldC*j] =(tangentB[ip]+2*tangentB[i]+tangentB[im])/4/nt;
      D[nBorders+1+ldC*j]-=(tangentB[ip]-2*tangentB[i]+tangentB[im])*nt;
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
    for(j=0;j<n-k-nBorders;j++)printf("-------");
    printf(" + ");
    for(j=0;j<nBorders+k;j++)printf("-------");
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
  if(ierr!=0)
   {
    printf(" Problem with factor, zero on diagonal %d\n",ierr);
    exit(8);
   }
  ierr=0;
  CALLDBOSL(A,&ldA,&na,&nBandsL,&nBandsU,B,&ldB,&nb,C,&ldC,D,&ldD,X,Y,R,S,NULL,Pivots,&ierr);

/* Extract the result */

  for(i=0;i<n-k-nBorders;i++)
    Phi[i+n*t]=R[i];
  for(i=0;i<nBorders+k;i++)
    Phi[n-k-nBorders+i+n*t]=S[i];

/* Normalize */

  s=0.;
  for(i=0;i<nt;i++)
   {
    ip=i+1;if(ip>=nt)ip=0;
    s+=(Phi[ip+n*t]+Phi[i+n*t])*(Phi[ip+n*t]+Phi[i+n*t])/4/nt;
    s+=(Phi[ip+n*t]-Phi[i+n*t])*(Phi[ip+n*t]-Phi[i+n*t])*nt;
   }
  if(fabs(s)<.00001){printf("Warning, result small in norm -- MFForcedOscillator Single Tangent\n");fflush(stdout);}
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

double TGuuForcedOscillator(int i,int j,int k, double *point,int nt, double r, MFErrorHandler e)
 {
  static char RoutineName[]={"TGuuForcedOscillator"};
  double eps=.001;
  double t;
  double gm,gp;
  double test;
  double result;
  printf("%s\n",RoutineName);fflush(stdout);

  t=point[k];

  point[k]=t+eps/2;
  gp=TGuForcedOscillator(i,j,point,nt,r,e);

  point[k]=t-eps/2;
  gm=TGuForcedOscillator(i,j,point,nt,r,e);

  point[k]=t;

  result=(gp-gm)/eps;
  test=GuuForcedOscillator(i,j,k,point,nt,r,e);
  if(fabs(test-result)>eps)
   {printf("Large error in Guu, %d, %d, %d, diff=%lf, exact=%lf\n",i,j,k,result,test);fflush(stdout);}
  
  return result;
 }

double FGuuForcedOscillator(int i,int j,int k, double *point,int nt, double r, MFErrorHandler e)
 {
  static char RoutineName[]={"FGuuForcedOscillator"};
  double eps=.001;
  double t;
  double gm,gp;
  double test;
  double result;

  t=point[k];

  point[k]=t+eps/2;
  gp=GuForcedOscillator(i,j,point,nt,r,e);

  point[k]=t-eps/2;
  gm=GuForcedOscillator(i,j,point,nt,r,e);

  point[k]=t;

  result=(gp-gm)/eps;
  
  return result;
 }

double TGuForcedOscillator(int i,int j, double *point,int nt, double r, MFErrorHandler e)
 {
  static char RoutineName[]={"TGuForcedOscillator"};
  double eps=.001;
  double t;
  double gm,gp;
  double test;
  double result;
  printf("%s\n",RoutineName);fflush(stdout);

  t=point[j];

  point[j]=t+eps/2;
  gp=GForcedOscillator(i,point,nt,r,e);

  point[j]=t-eps/2;
  gm=GForcedOscillator(i,point,nt,r,e);

  point[j]=t;

  result=(gp-gm)/eps;
  test=GuForcedOscillator(i,j,point,nt,r,e);
  if(fabs(test-result)>eps)
   {printf("Large error in Gu, %d, %d, diff=%lf, exact=%lf\n",i,j,result,test);fflush(stdout);}
  
  return result;
 }

void MFCurvatureForcedOscillator(int n,int k,double *u,double *tangentA,double *tangentB,double *a00, double *a01, double *a11, void *d, MFErrorHandler eh)
 {
  static double *result;
  static int i,j,l;
  static int i0,i1;
  static int verbose,test;
  static int nBands,bandW;
  static int na,nb;
  static int ierr;
  static double e;
  static double scl;
  static int nt;
  static double r;
  int ip,im;

  nt=((struct MFForcedOscillatorData*)d)->nt;
  r=((struct MFForcedOscillatorData*)d)->R;

  verbose=0;
  nBands=nBandsL+nBandsU+1;

  for(j=0;j<n-k-nBorders;j++)
   {
    i0=j-nBandsU;
    if(i0<0)i0=0;
    i1=j+nBandsL;
    if(i1>n-k-nBorders-1)i1=n-k-nBorders-1;
    for(i=i0;i<i1+1;i++)
      A[i-j+nBands-1+ldA*j]=GU(i,j,u,nt,r,eh);
   }

  for(i=0;i<n-k-nBorders;i++)
   {
    for(j=0;j<nBorders+k;j++)
      B[i+ldB*j]=GU(i,n-k-nBorders+j,u,nt,r,eh);
   }

  for(i=0;i<nBorders;i++)
   {
    for(j=0;j<n-k-nBorders;j++)
      C[i+ldC*j]=GU(n-k-nBorders+i,j,u,nt,r,eh);
   }
  for(j=0;j<n-k-nBorders;j++)
   {
    i=j;
    im=i-1;if(im<0)im=nt-1;
    ip=i+1;if(ip>=nt)ip=0;
    if(j<nt)
     {
      C[nBorders+0+ldC*j] =(tangentA[ip]+2*tangentA[i]+tangentA[im])/4/nt;
      C[nBorders+0+ldC*j]-=(tangentA[ip]-2*tangentA[i]+tangentA[im])*nt;
      C[nBorders+1+ldC*j] =(tangentB[ip]+2*tangentB[i]+tangentB[im])/4/nt;
      C[nBorders+1+ldC*j]-=(tangentB[ip]-2*tangentB[i]+tangentB[im])*nt;
     }else{
      C[nBorders+0+ldC*j]=tangentA[j];
      C[nBorders+1+ldC*j]=tangentB[j];
     }
   }

  for(i=0;i<nBorders;i++)
   {
    for(j=0;j<nBorders+k;j++)
      D[i+ldD*j]=GU(n-k-nBorders+i,n-k-nBorders+j,u,nt,r,eh);
   }
  for(j=0;j<nBorders+k;j++)
   {
    i=n-k-nBorders+j;
    im=i-1;if(im<0)im=nt-1;
    ip=i+1;if(ip>=nt)ip=0;
    if(n-k-nBorders+j<nt)
     {
      D[nBorders+0+ldC*j] =(tangentA[ip]+2*tangentA[i]+tangentA[im])/4/nt;
      D[nBorders+0+ldC*j]-=(tangentA[ip]-2*tangentA[i]+tangentA[im])*nt;
      D[nBorders+1+ldC*j] =(tangentB[ip]+2*tangentB[i]+tangentB[im])/4/nt;
      D[nBorders+1+ldC*j]-=(tangentB[ip]-2*tangentB[i]+tangentB[im])*nt;
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
    for(j=0;j<n-k-nBorders;j++)printf("-------");
    printf(" + ");
    for(j=0;j<nBorders+k;j++)printf("-------");
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
  if(ierr!=0)
   {
    printf(" Problem with factor, zero on diagonal %d\n",ierr);
    exit(8);
   }

/* -GuuPhiSPhiT */

  for(i=0;i<n-k-nBorders;i++)
   {
    R[i]=0.;
    for(j=MAX(0,i-1);j<MIN(n,i+1);j++)
     {
      R[i]+=GUU(i,j,j,u,nt,r,eh)*tangentA[j]*tangentA[j];
      R[i]+=2*GUU(i,j,j+1,u,nt,r,eh)*tangentA[j]*tangentA[j+1];
     }
   }

  for(i=0;i<nBorders;i++)
   {
    S[i]=0.;
    for(j=MAX(0,i-1);j<MIN(n,i+1);j++)
     {
      S[i]+=GUU(n-k-nBorders+i,j,j,u,nt,r,eh)*tangentA[j]*tangentA[j];
      S[i]+=2*GUU(n-k-nBorders+i,j,j+1,u,nt,r,eh)*tangentA[j]*tangentA[j+1];
     }
   }
  for(i=nBorders;i<nBorders+k;i++)S[i]=0.;

  e=0.;
  for(i=0;i<n-k-nBorders;i++)e+=R[i]*R[i]/nt;
  for(i=0;i<nBorders;i++)e+=S[i]*S[i]/nt;
  for(i=nBorders;i<nBorders+k;i++){if(i+n-k-nBorders==nt)e+=1.*(S[i]*S[i]);else e+=S[i]*S[i];}
  e=sqrt(e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("Norm of rhs for a00 is %lf\n",eh);fflush(stdout);}
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
    for(j=MAX(0,i-1);j<MIN(n,i+1);j++)
     {
      R[i]+=GUU(i,j,j,u,nt,r,eh)*tangentA[j]*tangentB[j];
      R[i]+=GUU(i,j,j+1,u,nt,r,eh)*(tangentA[j]*tangentB[j+1]+tangentA[j+1]*tangentB[j]);
     }
   }

  for(i=0;i<nBorders;i++)
   {
    S[i]=0.;
    for(j=MAX(0,i-1);j<MIN(n,i+1);j++)
     {
      S[i]+=GUU(n-k-nBorders+i,j,j,u,nt,r,eh)*tangentA[j]*tangentB[j];
      S[i]+=GUU(n-k-nBorders+i,j,j+1,u,nt,r,eh)*(tangentA[j]*tangentB[j+1]+tangentA[j+1]*tangentB[j]);
     }
   }
  for(i=nBorders;i<nBorders+k;i++)S[i]=0.;

  e=0.;
  for(i=0;i<n-k-nBorders;i++)e+=R[i]*R[i]/nt;
  for(i=0;i<nBorders;i++)e+=S[i]*S[i]/nt;
  for(i=nBorders;i<nBorders+k;i++){if(i+n-k-nBorders==nt)e+=1.*(S[i]*S[i]);else e+=S[i]*S[i];}
  e=sqrt(e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("Norm of rhs for a01 is %lf\n",eh);fflush(stdout);}
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
    for(j=MAX(0,i-1);j<MIN(n,i+1);j++)
     {
      R[i]+=GUU(i,j,j,u,nt,r,eh)*tangentB[j]*tangentB[j];
      R[i]+=2*GUU(i,j,j+1,u,nt,r,eh)*tangentB[j]*tangentB[j+1];
     }
   }

  for(i=0;i<nBorders;i++)
   {
    S[i]=0.;
    for(j=MAX(0,i-1);j<MIN(n,i+1);j++)
     {
      S[i]+=GUU(n-k-nBorders+i,j,j,u,nt,r,eh)*tangentB[j]*tangentB[j];
      S[i]+=2*GUU(n-k-nBorders+i,j,j+1,u,nt,r,eh)*tangentB[j]*tangentB[j+1];
     }
   }
  for(i=nBorders;i<nBorders+k;i++)S[i]=0.;

  e=0.;
  for(i=0;i<n-k-nBorders;i++)e+=R[i]*R[i]/nt;
  for(i=0;i<nBorders;i++)e+=S[i]*S[i]/nt;
  for(i=nBorders;i<nBorders+k;i++){if(i+n-k-nBorders==nt)e+=1.*(S[i]*S[i]);else e+=S[i]*S[i];}
  e=sqrt(e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("Norm of rhs for a11 is %lf\n",e);fflush(stdout);}
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

void MFWriteForcedOscillatorData(FILE *fid,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFWriteForcedOscillatorData"};
  struct MFForcedOscillatorData *data;
  int verbose;

  verbose=0;

  if(verbose){printf("%s\n",RoutineName);fflush(stdout);}

  data=(struct MFForcedOscillatorData*)d;
  fprintf(fid,"%d %lf\n",data->nt,data->R);
  fflush(fid);

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("  nt=%d\n",data->nt);
    printf("  R=%lf\n",data->R);fflush(stdout);
   }
#endif

  return;
 }

MFImplicitMF MFReadForcedOscillator(FILE *fid, MFErrorHandler e)
 {
  static char RoutineName[]={"MFReadForcedOscillator"};
  MFImplicitMF forced;
  int nt=0;
  double R=0.;
  int verbose;

  verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s\n",RoutineName);fflush(stdout);}
#endif

  fscanf(fid,"%d %lf\n",&nt,&R);
  forced=MFIMFCreateForcedOscillator(nt,1,e);

  return forced;
 }

int MFForcedOscillatorProjectToSave(MFNVector u, double *x, void *d, MFErrorHandler e)
 {
  int i,n;
  struct MFForcedOscillatorData *data;

  data=(struct MFForcedOscillatorData*)d;

  n=data->nt+2;

  if(x==NULL)return n;
  for(i=0;i<n;i++)x[i]=MFNV_C(u,i,e);

  return 0;
 }

int MFForcedOscillatorProjectForBB(MFNVector u, double *x, void *d, MFErrorHandler e)
 {
  return MFForcedOscillatorProjectToDraw(u,x,d,e);

  return 0;
 }

int MFForcedOscillatorProjectToDraw(MFNVector u, double *x, void *d, MFErrorHandler e)
 {
  int i,n;
  double h;
  struct MFForcedOscillatorData *data;

  data=(struct MFForcedOscillatorData*)d;

  n=data->nt+2;

  if(x==NULL)return 4;
  x[0]=MFNV_C(u,n-2,e);
  x[1]=MFNV_C(u,n-1,e);
  x[2]=MFNV_C(u,0,e);
  h=data->R*3.1415926/data->nt;
  x[3]=(MFNV_C(u,1,e)-MFNV_C(u,0,e))/h;

  return 0;
 }

double GForcedOscillator(int i,double *point,int nt, double r, MFErrorHandler eh)
 {
  static char RoutineName[]={"GForcedOscillator"};
  double e,B;
  double t,h;
  double up,um,vp,vm;
  double result;
  int im,i0,ip;
  int n,k;

/*
 
   EVALUATES G (U,L)
                                 I
      U IS:
             1     2     3         N-4    N-3   N-2  N-1 N
 
           P1(1),P1(2),P1(3),...,P1(M-1),P1(M),  w ,  e, B

      N=nt+2
      k=2
 
*/
  n=nt+2;
  k=2;

  if(i<0 || i>=n-k)
   {
    printf("i<0 || i>=n-k in GForcedOscillator i=%d, n=%d\n",i,n-k);fflush(stdout);
    return(0.);
   }

  e=point[nt+0];
  B=.05*point[nt+1];

  h=r*3.1415926/nt;
  t=i*h;

  im=i-1;
  i0=i;
  ip=i+1;

  if(i0== 0)im=nt-1;
  if(i0==nt-1)ip=0;

  result=0.;

  up=(point[ip]-point[i0])/h;
  um=(point[i0]-point[im])/h;
  vp=(point[ip]+point[i0])/2;
  vm=(point[i0]+point[im])/2;

  result=(up-um)/h
         +(e-B)*(up+um)/2*(up+um)/2*(up+um)/2
         -(e/2-B)*(up+um)/2
         +((1+B*sin(2*(t+h/2)))*vp+(1+B*sin(2*(t-h/2)))*vm)/2;
  result=result*h*h;
  return(result);
 }

double GuForcedOscillator(int i,int j, double *point,int nt, double r, MFErrorHandler eh)
 {
  static char RoutineName[]={"GuForcedOscillator"};
  double e,B;
  double t,h;
  double up,um,vp,vm;
  double dupdip,dupdi0,dupdim;
  double dumdip,dumdi0,dumdim;
  double dvpdip,dvpdi0,dvpdim;
  double dvmdip,dvmdi0,dvmdim;
  double result;
  int im,i0,ip;

  if(i<0 || i>nt)
   {
    printf("i<0 || i>n-k, in Gu i=%d, n-k=%d\n",i,nt);
    return(0.);
   }

  if(j<0 || j>nt+2)
   {
    printf("j<0 || j>n, in Gu i=%d, j=%d, n=%d\n",i,j,nt+2);
    return(0.);
   }

  e=point[nt+0];
  B=0.05*point[nt+1];

/*
 
   EVALUATES GU(U,L)
                                 I
      U IS:
             1     2     3        N-4  N-3   N-2  N-1 N
 
           P1(1),P2(1),P1(2),...,P1(M),P2(M), T , In, Ic
 
*/

  h=r*3.1415926/nt;
  t=i*h;

  im=i-1;
  i0=i;
  ip=i+1;

  if(i0== 0)im=nt-1;
  if(i0==nt-1)ip=0;

  up=(point[ip]-point[i0])/h; dupdip= 1/h; dupdi0=-1/h; dupdim= 0.;
  um=(point[i0]-point[im])/h; dumdip= 0.; dumdi0=1/h; dumdim= -1/h;
  vp=(point[ip]+point[i0])/2; dvpdip= .5; dvpdi0=.5; dvpdim= 0.;
  vm=(point[i0]+point[im])/2; dvmdip= 0.; dvmdi0=.5; dvmdim= .5;

  if(j==ip)
   {
    result=(dupdip-dumdip)/h
       +3*(e-B)*(up+um)/2*(up+um)/2*(dupdip+dumdip)/2
         -(e/2-B)*(dupdip+dumdip)/2
       +((1+B*sin(2*(t+h/2)))*dvpdip+(1+B*sin(2*(t-h/2)))*dvmdip)/2;
   }else if(j==i0)
   {
    result=(dupdi0-dumdi0)/h
       +3*(e-B)*(up+um)/2*(up+um)/2*(dupdi0+dumdi0)/2
         -(e/2-B)*(dupdi0+dumdi0)/2
       +((1+B*sin(2*(t+h/2)))*dvpdi0+(1+B*sin(2*(t-h/2)))*dvmdi0)/2;
   }else if(j==im)
   {
    result=(dupdim-dumdim)/h
       +3*(e-B)*(up+um)/2*(up+um)/2*(dupdim+dumdim)/2
         -(e/2-B)*(dupdim+dumdim)/2
       +((1+B*sin(2*(t+h/2)))*dvpdim+(1+B*sin(2*(t-h/2)))*dvmdim)/2;
   }else if(j==nt  )
   {   /* e */
    result=(up+um)/2*(up+um)/2*(up+um)/2-(up+um)/4;
   }else if(j==nt+1)
   {   /* B */
    result=-(up+um)/2*(up+um)/2*(up+um)/2+(up+um)/2
         +(sin(2*(t+h/2))*vp+sin(2*(t-h/2))*vm)/2;
    result=.05*result;
   }else
    result=0.;

  result=result*h*h;
  
  return(result);
 }

double GuuForcedOscillator(int i,int j,int k, double *point,int nt, double r, MFErrorHandler eh)
 {
  static char RoutineName[]={"GuuForcedOscillator"};
  double e,B;
  double up,um,vp,vm;
  double dupdip,dupdi0,dupdim,dupdw;
  double dumdip,dumdi0,dumdim,dumdw;
  double dvpdip,dvpdi0,dvpdim,dvpdw;
  double dvmdip,dvmdi0,dvmdim,dvmdw;
  double t,h;
  double result;
  int im,i0,ip;

  if(i<0 || i>nt)
   {
    printf("i<0 || i>n-k, in Guu i=%d, n-k=%d\n",i,nt);
    abort();
    return(0.);
   }

  if(j<0 || j>nt+2)
   {
    printf("j<0 || j>n, in Guu i=%d, j=%d, k=%d, n=%d\n",i,j,k,nt+2);
    return(0.);
   }

  if(k<0 || k>nt+2)
   {
    printf("k<0 || k>n, in Guu i=%d, j=%d, k=%d, n=%d\n",i,j,k,nt+2);
    return(0.);
   }

  e=point[nt  ];
  B=0.05*point[nt+1];

/*
 
   EVALUATES GUU(U,L)
                                 I
      U IS:
             1     2     3        N-4  N-3   N-2  N-1 N
 
           P1(1),P2(1),P1(2),...,P1(M),P2(M), T , In, Ic
 
*/

  h=r*3.1415926/nt;
  t=i*h;

  im=i-1;
  i0=i;
  ip=i+1;

  if(i0== 0)im=nt-1;
  if(i0==nt-1)ip=0;

  up=(point[ip]-point[i0])/h; dupdip= 1/h; dupdi0=-1/h; dupdim= 0.;
  um=(point[i0]-point[im])/h; dumdip= 0.; dumdi0=1/h; dumdim= -1/h;
  vp=(point[ip]+point[i0])/2; dvpdip= .5; dvpdi0=.5; dvpdim= 0.;
  vm=(point[i0]+point[im])/2; dvmdip= 0.; dvmdi0=.5; dvmdim= .5;
  if(j==ip)
   {
    if(k==ip)
     {
/* ip,ip */
      result=6*(e-B)*(up+um)/2*(dupdip+dumdip)/2*(dupdip+dumdip)/2;
     }else if(k==i0)
     {
/* i0,ip */
      result=6*(e-B)*(up+um)/2*(dupdip+dumdip)/2*(dupdi0+dumdi0)/2;
     }else if(k==im)
     {
/* im,ip */
      result=6*(e-B)*(up+um)/2*(dupdip+dumdip)/2*(dupdim+dumdim)/2;
     }else if(k==nt)
     {
/* ip,e */
      result=3*(up+um)/2*(up+um)/2*(dupdip+dumdip)/2+(dupdip+dumdip)/4;
     }else if(k==nt+1)
     {
/* ip,B */
      result=-3*(up+um)/2*(up+um)/2*(dupdip+dumdip)/2+(dupdip+dumdip)/2
       +(sin(2*(t+h/2))*dvpdip+sin(2*(t-h/2))*dvmdip)/2;
      result=.05*result;
     }else result=0.;
     }else if(j==i0)
     {
      if(k==ip)
       {
/* ip,i0 */
        result=6*(e-B)*(up+um)/2*(dupdi0+dumdi0)/2*(dupdip+dumdip)/2;
       }else if(k==i0)
       {
/* i0,i0 */
        result=6*(e-B)*(up+um)/2*(dupdi0+dumdi0)/2*(dupdi0+dumdi0)/2;
       }else if(k==im)
       {
/* im,i0 */
        result=6*(e-B)*(up+um)/2*(dupdi0+dumdi0)/2*(dupdim+dumdim)/2;
       }else if(k==nt)
       {
/* i0,e */
        result=3*(up+um)/2*(up+um)/2*(dupdi0+dumdi0)/2+(dupdi0+dumdi0)/4;
       }else if(k==nt+1)
       {
/* i0,B */
        result=-3*(up+um)/2*(up+um)/2*(dupdi0+dumdi0)/2+(dupdi0+dumdi0)/2
         +(sin(2*(t+h/2))*dvpdi0+sin(2*(t-h/2))*dvmdi0)/2;
        result=.05*result;
       }else result=0.;
     }else if(j== im)
     {
      if(k==ip)
       {
/* ip,im */
        result=6*(e-B)*(up+um)/2*(dupdip+dumdip)/2*(dupdim+dumdim)/2;
       }else if(k==i0)
       {
/* i0,im */
        result=6*(e-B)*(up+um)/2*(dupdi0+dumdi0)/2*(dupdim+dumdim)/2;
       }else if(k==im)
       {
/* im,im */
        result=6*(e-B)*(up+um)/2*(dupdim+dumdim)/2*(dupdim+dumdim)/2;
       }else if(k==nt)
       {
/* im,e */
        result=3*(up+um)/2*(up+um)/2*(dupdim+dumdim)/2+(dupdim+dumdim)/4;
       }else if(k==nt+1)
       {
/* im,B */
        result=-3*(up+um)/2*(up+um)/2*(dupdim+dumdim)/2+(dupdim+dumdim)/2
         +(sin(2*(t+h/2))*dvpdim+sin(2*(t-h/2))*dvmdim)/2;
        result=.05*result;
       }else result=0.;
     }else if(j==nt)
     {   
      if(k==ip)
       {
/* ip,e */
        result=3*(up+um)/2*(up+um)/2*(dupdip+dumdip)/2+(dupdip+dumdip)/4;
       }else if(k==i0)
       {
/* i0,e */
        result=3*(up+um)/2*(up+um)/2*(dupdi0+dumdi0)/2+(dupdi0+dumdi0)/4;
       }else if(k==im)
       {
/* im,e */
        result=3*(up+um)/2*(up+um)/2*(dupdim+dumdim)/2+(dupdim+dumdim)/4;
       }else if(k==nt)
       {
/* e,e */
        result=0.;
       }else if(k==nt+1)
       {
/* B,e */
        result=0.;
        result=.05*result;
       }else result=0.;
     }else if(j==nt+1)
     {
      if(k==ip)
       {
/* ip,B */
        result=-3*(up+um)/2*(up+um)/2*(dupdip+dumdip)/2+(dupdip+dumdip)/2
         +(sin(2*(t+h/2))*dvpdip+sin(2*(t-h/2))*dvmdip)/2;
        result=.05*result;
       }else if(k==i0)
       {
/* i0,B */
        result=-3*(up+um)/2*(up+um)/2*(dupdi0+dumdi0)/2+(dupdi0+dumdi0)/2
         +(sin(2*(t+h/2))*dvpdi0+sin(2*(t-h/2))*dvmdi0)/2;
        result=.05*result;
       }else if(k==im)
       {
/* im,B */
        result=-3*(up+um)/2*(up+um)/2*(dupdim+dumdim)/2+(dupdim+dumdim)/2
         +(sin(2*(t+h/2))*dvpdim+sin(2*(t-h/2))*dvmdim)/2;
        result=.05*result;
       }else if(k==nt)
       {
/* e,B */
        result=0.;
        result=.05*result;
       }else if(k==nt+1)
       {
/* B,B */
        result=0.;
        result=.05*.05*result;
       }else result=0.;
     }else result=0.;
  result=result*h*h;

  return(result);
 }

void ForcedIntegrate(double x, int p, int q, int nt, MFErrorHandler e)
 {
  double xp,x0,xm;
  double f,fp;
  double h;

/*h=q*2*3.1415926/nt;

  up=(xp-x0)/h;
  um=(x0-xm)/h;
  vp=(xp+x0)/2;
  vm=(x0+xm)/2;
  f=(up-um)/h
         +(B-e)*(up+um)/2*(up+um)/2*(up+um)/2
         -(B-e/2)*(up+um)/2
         +((1+B*sin(2*(t+h/2)))*vp+(1+B*sin(2*(t-h/2)))*vm)/2;*/
 } 

double MFForcedOscillatorNSpaceDistance(MFNSpace,MFNVector,MFNVector,void*,MFErrorHandler);
double MFForcedOscillatorNSpaceInner(MFNSpace,MFNVector,MFNVector,void*,MFErrorHandler);
void MFForcedOscillatorNSpaceDirection(MFNSpace,MFNVector,MFNVector,MFNVector,void*,MFErrorHandler);
void MFForcedOscillatorNSpaceAdd(MFNSpace,MFNVector,MFNVector,MFNVector,void*,MFErrorHandler);
void MFForcedOscillatorNSpaceScale(MFNSpace,double,MFNVector,MFNVector,void*,MFErrorHandler);

MFNSpace MFCreateForcedOscillatorNSpace(MFErrorHandler e)
 {
  static char RoutineName[]={"MFCreateForcedOscillatorNSpace"};
  MFNSpace thisOscillator;

  thisOscillator=MFCreateNSpaceBaseClass("ForcedOscillatorNSpace",e);
  MFNSpaceSetDistance(thisOscillator,MFForcedOscillatorNSpaceDistance,e);
  MFNSpaceSetInnerProduct(thisOscillator,MFForcedOscillatorNSpaceInner,e);
  MFNSpaceSetDirection(thisOscillator,MFForcedOscillatorNSpaceDirection,e);
  MFNSpaceSetAdd(thisOscillator,MFForcedOscillatorNSpaceAdd,e);
  MFNSpaceSetScale(thisOscillator,MFForcedOscillatorNSpaceScale,e);

  return thisOscillator;
 }

double MFForcedOscillatorNSpaceDistance(MFNSpace thisOscillator,MFNVector v0,MFNVector v1,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFForcedOscillatorNSpaceDistance"};
  double result;
  int i,ip,n;
  double dt;
  double xp,x;

  n=MFNV_NC(v0,e);
  dt=1./(n-2);
  result=0.;
  for(i=0;i<n-2;i++)
   {
    ip=i+1;if(ip>=n-2)ip=0;
    xp=MFNV_CStar(v1,e)[ip]-MFNV_CStar(v0,e)[ip];
    x =MFNV_CStar(v1,e)[i ]-MFNV_CStar(v0,e)[i ];
    result+=(xp+x)*(xp+x)/4.*dt+(xp-x)*(xp-x)/dt;
   }
  result+=(MFNV_CStar(v1,e)[n-2]-MFNV_CStar(v0,e)[n-2])*(MFNV_CStar(v1,e)[n-2]-MFNV_CStar(v0,e)[n-2]);
  result+=(MFNV_CStar(v1,e)[n-1]-MFNV_CStar(v0,e)[n-1])*(MFNV_CStar(v1,e)[n-1]-MFNV_CStar(v0,e)[n-1]);
  result=sqrt(result);

  return result;
 }

double MFForcedOscillatorNSpaceInner(MFNSpace thisOscillator,MFNVector v0,MFNVector v1,void *d,MFErrorHandler e)
 {
  static char RoutineName[]={"MFForcedOscillatorNSpaceInner"};
  double result;
  int i,ip,n;
  double dt;
  double xp,x,yp,y;

  n=MFNV_NC(v0,e);
  dt=1./(n-2);
  result=0.;
  for(i=0;i<n-2;i++)
   {
    ip=i+1;if(ip>=n-2)ip=0;
    xp=MFNV_CStar(v1,e)[ip];
    x =MFNV_CStar(v1,e)[i ];
    yp=MFNV_CStar(v0,e)[ip];
    y =MFNV_CStar(v0,e)[i ];
    result+=(xp+x)*(yp+y)/4.*dt+(xp-x)*(yp-y)/dt;
   }
  result+=MFNV_CStar(v1,e)[n-2]*MFNV_CStar(v0,e)[n-2];
  result+=MFNV_CStar(v1,e)[n-1]*MFNV_CStar(v0,e)[n-1];

  return result;
 }

void MFForcedOscillatorNSpaceDirection(MFNSpace thisOscillator,MFNVector v0,MFNVector v1,MFNVector diff,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFForcedOscillatorNSpaceDirection"};
  int i,n;

  n=MFNV_NC(v0,e);
  for(i=0;i<n;i++)MFNV_CStar(diff,e)[i]=MFNV_CStar(v1,e)[i]-MFNV_CStar(v0,e)[i];
 }

void MFForcedOscillatorNSpaceAdd(MFNSpace thisOscillator,MFNVector v0,MFNVector v1,MFNVector sum,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFForcedOscillatorNSpaceAdd"};
  int i,n;

  n=MFNV_NC(v0,e);
  for(i=0;i<n;i++)MFNV_CStar(sum,e)[i]=MFNV_CStar(v1,e)[i]+MFNV_CStar(v0,e)[i];

  return;
 }

void MFForcedOscillatorNSpaceScale(MFNSpace thisOscillator,double s, MFNVector v,MFNVector prod,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFForcedOscillatorNSpaceScale"};
  int i,n;

  n=MFNV_NC(v,e);
  for(i=0;i<n;i++)MFNV_CStar(prod,e)[i]=s*MFNV_CStar(v,e)[i];

  return;
 }

#ifdef __cplusplus
}
#endif
