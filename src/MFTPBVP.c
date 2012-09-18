/* 
 *  PROGRAM NAME:  multifario
 *
 *  (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
 *  CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
 *
 *  Please refer to the LICENSE file in the top directory
 *
 *      author: Mike Henderson mhender@watson.ibm.com
 *      date:   July 27, 2000
 */

static char *id="@(#) $Id: MFTPBVP.c,v 1.12 2011/07/21 17:42:46 mhender Exp $";

#include <multifarioConfig.h>
#include <MFImplicitMF.h>
#include <MFNSpace.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <MFFortran.h>
#include <MFErrorHandler.h>
#include <MFPrint.h>
#include <MFTPBVP.h>

#define round(x) ((x)>=0?(long)((x)+0.5):(long)((x)-0.5))

#ifdef __cplusplus
 extern "C" {
#endif

static double epsilon;

static char MFTPBVPMFErrorHandlerMsg[256]="";

/*!
 * \brief A function defining the right hand side of a two point boundary value problem.
 */
/*typedef void (*MFTPBVPFFUNCTION)(double,int,double*,int,double*,double*,double*,double*);*/

/*!
 * \brief A function defining the boundary value equations of a two point boundary value problem.
 */
/*typedef void (*MFTPBVPAFUNCTION)(int,int,double*,double*,int,double*,double*,double*,double*,double*);*/

/*!
 * \brief A function defining the integral piece of an integral constraint of a two point boundary value problem.
 */
/*typedef void (*MFTPBVPLFUNCTION)(int,double,int,double*,int,double*,double*,double*,double*);*/

/*!
 * \brief A function defining the non-vector piece of an integral constraint of a two point boundary value problem.
 */
/*typedef void (*MFTPBVPMFUNCTION)(int,int,double*,double*,double*);*/


struct MFTPBVPData
 {
  int nx;
  int nu;
  int np;
  int nbc;
  int nic;
  int k;
  MFTPBVPFFUNCTION f;
  MFTPBVPFFUNCTION fu;
  MFTPBVPFFUNCTION fl;
  MFTPBVPAFUNCTION a;
  MFTPBVPAFUNCTION au;
  MFTPBVPAFUNCTION al;
  MFTPBVPLFUNCTION l;
  MFTPBVPLFUNCTION lu;
  MFTPBVPLFUNCTION ll;
  MFTPBVPMFUNCTION m;
  MFTPBVPMFUNCTION ml;
  MFNSpace space;
 };

/* u'(x) = f(u(x),x,l) */

/* a(u(0),u(1),l) = 0 */

/* int l(u(x),x,l) dx + m(l) = 0 */

extern double MFEpsilon;

/*! \fn MFImplicitMF MFIMFCreateTPBVP(int k, int nx,int nu,int np, MFTPBVPFFUNCTION f, MFTPBVPFFUNCTION fu, MFTPBVPFFUNCTION fl, int nbc, MFTPBVPAFUNCTION a, MFTPBVPAFUNCTION au, MFTPBVPAFUNCTION al, int nic, MFTPBVPLFUNCTION l, MFTPBVPLFUNCTION lu, MFTPBVPLFUNCTION ll, MFTPBVPMFUNCTION m, MFTPBVPMFUNCTION ml);
 *  \brief Creates a manifold which is the solution manifold of a two point boundary value problem with integral constraints.
 *         Keller's second order box scheme is used.
 *
 *  \param k The number of degrees of freedom (the dimension of the solution manifold).
 *  \param nx The number of mesh intervals to use in the discretization.
 *  \param nu The number of functions defined on the mesh.
 *  \param np The number of scalar parameters.
 *  \param f  The right hand siade of the ODE's u'=f(u,p).
 *  \param fu The first derivatives of the right hand side with respect to u.
 *  \param fl The first derivatives of the right hand side with respect to the parameters l.
 *  \param nbc The number of boundary conditions.
 *  \param a  The function which defines the boundary conditions a(u(0),u(1),p)=0
 *  \param au The derivative of the boundary conditions with respect to u.
 *  \param al The derivative of the boundary conditions with respect to the parameters l.
 *  \param nic The number of integral conditions
 *  \param l The function which defines the integral part of the integral conditions int_0^1 l(u(t),p) dt + m(l)=0
 *  \param lu The derivative of the integral part of the integral conditions with respect to u.
 *  \param ll The derivative of the integral part of the integral conditions with respect to the parameters l.
 *  \param m The function which defines the non-integral part of the integral conditions int_0^1 l(u(t),p) dt + m(l)=0
 *  \param ml The derivative of the non-integral part of the integral conditions with respect to the parameters l.
 *  \returns An implicitly defined manifold.
 */
MFImplicitMF MFIMFCreateTPBVP(int k, int nx,int nu,int np,
     MFTPBVPFFUNCTION f,
     MFTPBVPFFUNCTION fu,
     MFTPBVPFFUNCTION fl,
     int nbc,
     MFTPBVPAFUNCTION a,
     MFTPBVPAFUNCTION au,
     MFTPBVPAFUNCTION al,
     int nic,
     MFTPBVPLFUNCTION l,
     MFTPBVPLFUNCTION lu,
     MFTPBVPLFUNCTION ll,
     MFTPBVPMFUNCTION m,
     MFTPBVPMFUNCTION ml, MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFCreateTPBVP"};
  MFImplicitMF thisBVP;
  struct MFTPBVPData *data;
  MFNSpace space;

  thisBVP=MFIMFCreateBaseClass(nx*nu+np+nx,k,"TPBVP",e);

  space=MFCreateTPBVPNSpace(nx,nu,np,e);
  MFIMFSetSpace(thisBVP,space,e);
  MFFreeNSpace(space,e);

  data=(struct MFTPBVPData *)malloc(sizeof(struct MFTPBVPData)); /*done*/

#ifndef MFNOSAFETYNET
  if(data==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct MFTPBVPData));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  data->np=np;
  data->nx=nx;
  data->nu=nu;
  data->nbc=nbc;
  data->nic=nic;
  data->k=k;
  data->f=f;
  data->fu=fu;
  data->fl=fl;
  data->a=a;
  data->au=au;
  data->al=al;
  data->l=l;
  data->lu=lu;
  data->ll=ll;
  data->m=m;
  data->ml=ml;
  data->space=space;
  epsilon=1.e-7;

  MFIMFSetData(thisBVP,(void*)data,e);
  MFIMFSetFreeData(thisBVP,MFFreeTPBVPData,e);
  MFIMFSetProject(thisBVP,MFProjectTPBVP,e);
  MFIMFSetTangent(thisBVP,MFTangentTPBVP,e);
  MFIMFSetTangentWithGuess(thisBVP,MFTangentTPBVPWithGuess,e);
  MFIMFSetScale(thisBVP,MFScaleTPBVP,e);
  MFIMFSetWriteData(thisBVP,MFWriteTPBVPData,e);
  MFIMFSetSingular(thisBVP,MFSingularTPBVP,e);
  MFIMFSetProjectForSave(thisBVP,MFTPBVPProjectToSave,e);
  MFIMFSetProjectForDraw(thisBVP,MFTPBVPProjectToDraw,e);
  MFIMFSetProjectForBB(thisBVP,MFTPBVPProjectForBB,e);

  MFIMFSetVectorFactory(thisBVP,MFNVectorFactory,e);
  MFIMFSetMatrixFactory(thisBVP,MFNKMatrixFactory,e);

  return thisBVP;
 }

void MFFreeTPBVPData(void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFFreeTPBVPData"};

  free(d);
 }

void MFWriteTPBVPData(FILE *fid,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFWriteTPBVPData"};
  struct MFTPBVPData *data;

  data=(struct MFTPBVPData*)d;

  fprintf(fid,"%d %d %d %d %d %d\n",data->k,data->nx,data->nu,data->np,data->nbc,data->nic);
  return;
 }

MFImplicitMF MFReadTPBVP(FILE *fid, MFErrorHandler e)
 {
  static char RoutineName[]={"MFReadTPBVP"};
  MFImplicitMF result;
  int k,nx,nu,np,nic,nbc;

  fscanf(fid,"%d %d %d %d %d %d\n",&k,&nx,&nu,&np,&nbc,&nic);

  result=MFIMFCreateTPBVP(k,nx,nu,np,
     NULL,
     NULL,
     NULL,
     nbc,
     NULL,
     NULL,
     NULL,
     nic,
     NULL,
     NULL,
     NULL,
     NULL,
     NULL,e);

  return result;
 }

int MFProjectTPBVP(int n,int k,MFNVector vu0,MFNKMatrix mPhi,MFNVector vu,void *d,int *index, MFErrorHandler e)
 {
  static char RoutineName[]={"MFProjectTPBVP"};
  static double *A=NULL;
  static double *B=NULL;
  static double *C=NULL;
  static double *D=NULL;
  static double *r=NULL;
  static double *s=NULL;
  static double *r0=NULL;
  static double *s0=NULL;
  static double *X0=NULL;
  int itimes;
  int nx,nu,np,nbc,nic;
  int i,j,J;
  double error,delta;
  double *u, *u0, *Phi;
  double tol;

  int verbose=0;
  int verbosest=0;

  u0=MFNV_CStar(vu0,e);
  Phi=MFNKM_CStar(mPhi,e);
  u=MFNV_CStar(vu,e);

  nx=((struct MFTPBVPData*)(d))->nx;
  nu=((struct MFTPBVPData*)(d))->nu;
  np=((struct MFTPBVPData*)(d))->np;
  nbc=((struct MFTPBVPData*)(d))->nbc;
  nic=((struct MFTPBVPData*)(d))->nic;
  X0=(double*)realloc((void*)X0,(2*nu+np)*(2*nu+np)*sizeof(double));

#ifndef MFNOSAFETYNET
  if(X0==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",(2*nu+np)*(2*nu+np)*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0;
   }
#endif

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s nx=%d, nu=%d, np=%d, nic=%d\n",RoutineName,nx,nu,np,nic);fflush(stdout);}
  if(verbosest)
   {
    printf("    u0:[");fflush(stdout);
    for(i=0;i<nx;i++)
     {
      if(i>0)printf(",\n        ");
      printf("(");
      for(j=0;j<nu;j++)
       {
        if(j>0)printf(",");
        printf("%lf",u0[j+nu*i]);fflush(stdout);
       }
      printf(")");fflush(stdout);
     }
    printf("]\n");fflush(stdout);
    printf("l:[(");fflush(stdout);
    for(i=0;i<np;i++)
     {
      if(i>0)printf(",");
      printf("%lf",u0[nx*nu+i]);fflush(stdout);
     }
    printf(")]\n");fflush(stdout);
    printf("x:[");fflush(stdout);
    for(i=0;i<nx;i++)
     {
      if(i>0)printf(",");
      if(i%6==5)printf("\n");
      printf("%lf",u0[nx*nu+np+i]);fflush(stdout);
     }
    printf("]\n");fflush(stdout);
   }
#endif

  for(i=0;i<n;i++)u[i]=u0[i];

/*MFTPBVPTestJacobian(nx,nu,np,nbc,nic,n,k,u0,Phi,d,e);*/

  itimes=0;
  error=1.;
  tol=1.e-7;
  while(error>tol && itimes<10)
   {
    error=MFTPBVPGetRes(n,k,d,u0,u,Phi,&r,&s,e);
#ifdef MFALLOWVERBOSE
    if(verbose)printf(" %d error %le\n",itimes,error);
    if(verbosest)
     {
      printf("residual:\n");
      for(i=0;i<(nx-1)*nu;i++)
        printf("  %4d %le\n",i,r[i]);
      printf("  ----------------------\n");
      for(i=0;i<nbc+nic+k;i++)
        printf("  %4d %le\n",i+(nx-1)*nu,s[i]);
     }
#endif
    if(itimes==0||error>1.e-10)
     {
      MFTPBVPGetJac(n,k,d,u,u0,Phi,&A,&B,&C,&D,e);
      error=MFTPBVPGetRes(n,k,d,u0,u,Phi,&r,&s,e);
      MFSolveBordered((nx-1)*nu,nbc+nic+k,nu-1,2*nu-1,A,B,C,D,r,s,NULL,X0,e);
      for(i=0;i<(nx-1)*nu;i++)delta+=r[i]*r[i];
      for(i=0;i<nbc+nic+k;i++)delta+=s[i]*s[i];

      for(i=0;i<(nx-1)*nu;i++)u[i]=u[i]-r[i];
      for(i=0;i<nu+np;i++)u[(nx-1)*nu+i]=u[(nx-1)*nu+i]-s[i];
     }
    itimes++;
   }

  *index=0;
  if(0)
   {
    static double *A=NULL;
    static double *B=NULL;
    static double *C=NULL;
    static double *D=NULL;
    static double *Y=NULL;
    static int *Pivots=NULL;
    int nrhs;
    int ierr;
    int i,j,l;
    char trans='N';

    A=(double*)realloc((void*)A,nu* nu*sizeof(double));

#ifndef MFNOSAFETYNET
    if(A==NULL)
     {
      sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",nu*nu*sizeof(double));
      MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return 0;
     }
#endif

    B=(double*)realloc((void*)B,nu*(nu+np)*sizeof(double));

#ifndef MFNOSAFETYNET
    if(B==NULL)
     {
      sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",nu*(2*nu+np)*sizeof(double));
      MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      free(A);
      return 0;
     }
#endif

    C=(double*)realloc((void*)C,(nu+np)*nu*sizeof(double));

#ifndef MFNOSAFETYNET
    if(C==NULL)
     {
      sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",(2*nu+np)*nu*sizeof(double));
      MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      free(A);
      free(B);
      return 0;
     }
#endif

    D=(double*)realloc((void*)D,(nu+np)*(nu+np)*sizeof(double));

#ifndef MFNOSAFETYNET
    if(D==NULL)
     {
      sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",(nu+np)*(nu+np)*sizeof(double));
      MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      free(A);
      free(B);
      free(C);
      return 0;
     }
#endif

    Y =(double*)realloc((void*)Y,(nu+np)*(nu+np)*sizeof(double));

#ifndef MFNOSAFETYNET
    if(Y==NULL)
     {
      sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",(nu+np)*(nu+np)*sizeof(double));
      MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      free(A);
      free(B);
      free(C);
      free(D);
      return 0;
     }
#endif

    Pivots=(int*)realloc((void*)Pivots,nu*sizeof(int));

#ifndef MFNOSAFETYNET
    if(Pivots==NULL)
     {
      sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",nu*sizeof(int));
      MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      free(A);
      free(B);
      free(C);
      free(D);
      free(Y);
      return 0;
     }
#endif

    for(i=0;i<nu;i++)
     {
      for(j=0;j<nu;j++)A[i+nu*j]=X0[i+(nbc+nic+k+nu)*j];
      for(j=0;j<nu+np;j++)B[i+nu*j]=X0[i+(nbc+nic+k+nu)*(j+nu)];
     }

    for(i=0;i<nu+np;i++)
     {
      for(j=0;j<nu;j++)C[i+(nu+np)*j]=X0[i+nu+(nbc+nic+k+nu)*j];
      for(j=0;j<nu+np;j++)D[i+(nu+np)*j]=X0[i+nu+(nbc+nic+k+nu)*(j+nu)];
     }

    CALLDGETRF(&nu,&nu,A,&nu,Pivots,&ierr);

    nrhs=nu+np;
    CALLDGETRS(&trans,&nu,&nrhs,A,&nu,Pivots,B,&nu,&ierr);

    for(i=0;i<nu+np;i++)
     {
      for(j=0;j<nu+np;j++)
       {
        Y[i+(nu+np)*j]=D[i+(nu+np)*j];
        for(l=0;l<nu;l++)
          Y[i+(nu+np)*j]-=C[i+(nu+np)*l]*B[l+nu*j];
       }
     }
    *index=MFAtlasDetSign(nbc+nic+k,Y,e);
   }else{
    *index=MFAtlasDetSign(2*nu+np,X0,e);
   }

#ifdef MFALLOWVERBOSE
  if(verbose){printf("done %s, index=%d (%le,%le)\n",RoutineName,*index,u[nx*nu],u[nx*nu+1]);fflush(stdout);}
#endif

  return error<tol;
 }

int MFSolveBordered(int n1,int n2,int nl,int nu,double *A,double *B,double *C,double *D,double *R,double *S,double *Z,double *X0, MFErrorHandler e)
 {
  static char RoutineName[]={"MFSolveBordered"};
  int i,j;
  int ierr;
  double *X=NULL;
  static double *W=NULL;
  static double *Y=NULL;
  static int *Pivots=NULL;
  static int ldA;
  static int ldX,result;
  int rc;

#ifdef HAVE_LAPACK

  ldA=2*nl+nu+1;
  ldX=n2+nl+1;

  if(Z==NULL)
   {
    W=(double*)realloc((void*)W,(ldA+n2)*(ldA+n2)*sizeof(double));

#ifndef MFNOSAFETYNET
    if(W==NULL)
     {
      sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",(ldA+n2)*(ldA+n2)*sizeof(double));
      MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return 0;
     }
#endif

    X=W;
   }else{
    X=Z;
   }
  Y=(double*)realloc((void*)Y,(ldA+n2)*sizeof(double));

#ifndef MFNOSAFETYNET
  if(Y==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",(ldA+n2)*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0;
   }
#endif

  Pivots=(int*)realloc((void*)Pivots,(n1+n2)*sizeof(int));

#ifndef MFNOSAFETYNET
    if(Pivots==NULL)
     {
      sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",(ldA+n2)*sizeof(double));
      MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return 0;
     }
#endif


  CALLDBOFA(A,&ldA,&n1,&nl,&nu,B,&n1,&n2,C,&n2,D,&n2,X,X0,Pivots,&ierr);
  rc=ierr;
  ierr=0;
  CALLDBOSL(A,&ldA,&n1,&nl,&nu,B,&n1,&n2,C,&n2,D,&n2,X,Y,R,S,X0,Pivots,&ierr);
  if(ierr!=0)printf("dbosl returned ierr=%d\n",ierr);
  if(rc==0)rc=ierr;

  return rc;
#else
  sprintf(MFTPBVPMFErrorHandlerMsg,"The TPBVP manifold requires dgesvd from LAPACK.");
  MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
  return 0;
#endif
 }

double MFGetBandedMatrixElement(int i,int j,double *A,int n,int nBandsL,int nBandsU, MFErrorHandler e)
 {
  static char RoutineName[]={"MFGetBandedMatrixElement"};
  int ldA=0;
  int nBands=0;

  ldA=2*nBandsL+nBandsU+1;
  nBands=nBandsL+nBandsU+1;

  if(i-j+nBands-1>-1 && i-j+nBands-1<ldA)
   return A[i-j+nBands-1+ldA*j];
  else
   return 0.;
 }

void MFSetBandedMatrixElement(int i,int j,double *A,double c,int n,int nBandsL,int nBandsU, MFErrorHandler e)
 {
  static char RoutineName[]={"MFSetBandedMatrixElement"};
  int ldA=0;
  int nBands=0;

  ldA=2*nBandsL+nBandsU+1;
  nBands=nBandsL+nBandsU+1;

  A[i-j+nBands-1+ldA*j]=c;
  return;
 }

void MFIncrementBandedMatrixElement(int i,int j,double *A,double c,int n,int nBandsL,int nBandsU, MFErrorHandler e)
 {
  static char RoutineName[]={"MFIncrementBandedMatrixElement"};
  int ldA=0;
  int nBands=0;

  ldA=2*nBandsL+nBandsU+1;
  nBands=nBandsL+nBandsU+1;
  A[i-j+nBands-1+ldA*j]+=c;
  return;
 }

void MFPrintBorderedBandedMatrix(FILE *fid, int n1,int n2,int nBandsL,int nBandsU,double *A,double *B,double *C,double *D, MFErrorHandler e)
 {
  static char RoutineName[]={"MFPrintBorderedBandedMatrix"};
  int i,j;
  int ldA=0;
  int nBands=0;

  ldA=2*nBandsL+nBandsU+1;
  nBands=nBandsL+nBandsU+1;

  fprintf(fid,"    ");
  for(i=0;i<n1;i++)fprintf(fid," %4d ",i);
  for(i=0;i<n2;i++)fprintf(fid," %4d ",i+n1);
  fprintf(fid,"\n");fflush(stdout);
  for(i=0;i<n1;i++)
   {
    fprintf(fid,"%2d [",i);
    for(j=0;j<n1;j++)
     {
      if(j>0)fprintf(fid," ");
      if(j<i-nBandsL||j>i+nBandsU)fprintf(fid,"     ");
       else fprintf(fid,"%5.2lf",A[i-j+nBands-1+ldA*j]);
     }
    fprintf(fid,"|");
    for(j=0;j<n2;j++)
     {
      if(j>0)fprintf(fid," ");
      fprintf(fid,"%5.2lf",B[i+n1*j]);
     }
    fprintf(fid,"]\n");
   }

  fprintf(fid,"   [");
  for(j=0;j<n1;j++)
   {
    if(j>0)fprintf(fid,"-");
    fprintf(fid,"-----");
   }
  fprintf(fid,"+");
  for(j=0;j<n2;j++)
   {
    if(j>0)fprintf(fid,"-");
    fprintf(fid,"-----");
   }
  fprintf(fid,"]\n");

  for(i=0;i<n2;i++)
   {
    fprintf(fid,"%2d [",i+n1);
    for(j=0;j<n1;j++)
     {
      if(j>0)fprintf(fid," ");
      fprintf(fid,"%5.2lf",C[i+n2*j]);
     }
    fprintf(fid,"|");
    for(j=0;j<n2;j++)
     {
      if(j>0)fprintf(fid," ");
      fprintf(fid,"%5.2lf",D[i+n2*j]);
     }
    fprintf(fid,"]\n");
   }
  fflush(stdout);

  return;
 }

void MFTPBVPGetJac(int n, int k, void *d,double *u,double *u0,double *Phi,double **pA,double **pB,double **pC,double **pD, MFErrorHandler e)
 {
  static char RoutineName[]={"MFTPBVPGetJac"};
  double *A;
  double *B;
  double *C;
  double *D;
  int nx,nu,np,nic,nbc;
  MFTPBVPFFUNCTION f,fu,fl;
  MFTPBVPAFUNCTION a,au,al;
  MFTPBVPLFUNCTION l,lu,ll;
  MFTPBVPMFUNCTION m,ml;
  int i,j,J,jj,o;
  double t;
  static double *U=NULL;
  static double *UL=NULL;
  static double *UR=NULL;
  static double *P=NULL;
  static double *Uref=NULL;
  static double *ULref=NULL;
  static double *URref=NULL;
  static double *Pref=NULL;
  static double *fv=NULL;
  static double *flv=NULL;
  static double *av=NULL;
  static double *lv=NULL;
  static double *llv=NULL;
  static double *mv=NULL;
  static int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s\n",RoutineName);fflush(stdout);}
#endif

  nx=((struct MFTPBVPData*)(d))->nx;
  nu=((struct MFTPBVPData*)(d))->nu;
  nbc=((struct MFTPBVPData*)(d))->nbc;
  np=((struct MFTPBVPData*)(d))->np;
  nic=((struct MFTPBVPData*)(d))->nic;
  k=((struct MFTPBVPData*)(d))->k;
  f=((struct MFTPBVPData*)(d))->f;
  a=((struct MFTPBVPData*)(d))->a;
  l=((struct MFTPBVPData*)(d))->l;
  m=((struct MFTPBVPData*)(d))->m;
  fu=((struct MFTPBVPData*)(d))->fu;
  au=((struct MFTPBVPData*)(d))->au;
  lu=((struct MFTPBVPData*)(d))->lu;
  fl=((struct MFTPBVPData*)(d))->fl;
  al=((struct MFTPBVPData*)(d))->al;
  ll=((struct MFTPBVPData*)(d))->ll;
  ml=((struct MFTPBVPData*)(d))->ml;

  *pA=(double*)realloc((void*)(*pA),(nx-1)*nu*(4*nu+nu+1)*sizeof(double));

#ifndef MFNOSAFETYNET
    if(*pA==NULL)
     {
      sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",(nx-1)*nu*sizeof(double));
      MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

  *pB=(double*)realloc((void*)(*pB),(nx-1)*nu*(nbc+nic+k)*sizeof(double));

#ifndef MFNOSAFETYNET
    if(*pB==NULL)
     {
      sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",(nx-1)*nu*(nbc+nic+k)*sizeof(double));
      MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

  *pC=(double*)realloc((void*)(*pC),(nx-1)*nu*(nbc+nic+k)*sizeof(double));

#ifndef MFNOSAFETYNET
    if(*pC==NULL)
     {
      sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",(nx-1)*nu*(nbc+nic+k)*sizeof(double));
      MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

  *pD=(double*)realloc((void*)(*pD),(nbc+nic+k)*(nbc+nic+k)*sizeof(double));

#ifndef MFNOSAFETYNET
    if(*pD==NULL)
     {
      sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",(nbc+nic+k-1)*(nbc+nic+k)*sizeof(double));
      MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

  fv=(double*)realloc((void*)fv,nu*(nu+np)*sizeof(double));

#ifndef MFNOSAFETYNET
    if(fv==NULL)
     {
      sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",nu*(nu+np)*sizeof(double));
      MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

  flv=(double*)realloc((void*)flv,nu*(nu+np)*sizeof(double));

#ifndef MFNOSAFETYNET
    if(flv==NULL)
     {
      sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",nu*(nu+np)*sizeof(double));
      MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

  av=(double*)realloc((void*)av,nbc*(2*nu+np)*sizeof(double));

#ifndef MFNOSAFETYNET
    if(av==NULL)
     {
      sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",nbc*(2*nu+np)*sizeof(double));
      MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

  if(nic>0)
   {
    lv=(double*)realloc((void*)lv,nic*(nu+np)*sizeof(double));

#ifndef MFNOSAFETYNET
    if(lv==NULL)
     {
      sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",nic*(nu+np)*sizeof(double));
      MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

    llv=(double*)realloc((void*)llv,nic*(nu+np)*sizeof(double));

#ifndef MFNOSAFETYNET
    if(llv==NULL)
     {
      sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",nic*(nu+np)*sizeof(double));
      MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

    mv=(double*)realloc((void*)mv,nic*(nu+np)*sizeof(double));

#ifndef MFNOSAFETYNET
    if(mv==NULL)
     {
      sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",nic*(nu+np)*sizeof(double));
      MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

   }

  U=(double*)realloc((void*)U,(nu+np+1)*sizeof(double));

#ifndef MFNOSAFETYNET
    if(U==NULL)
     {
      sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",(nu+np+1)*sizeof(double));
      MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

  UL=(double*)realloc((void*)UL,(nu+np+1)*sizeof(double));

#ifndef MFNOSAFETYNET
    if(UL==NULL)
     {
      sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",(nu+np+1)*sizeof(double));
      MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

  UR=(double*)realloc((void*)UR,(nu+np+1)*sizeof(double));

#ifndef MFNOSAFETYNET
    if(UR==NULL)
     {
      sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",(nu+np+1)*sizeof(double));
      MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

  P=(double*)realloc((void*)P,(nu+np+1)*sizeof(double));

#ifndef MFNOSAFETYNET
    if(P==NULL)
     {
      sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",(nu+np+1)*sizeof(double));
      MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

  Uref=(double*)realloc((void*)Uref,(nu+np+1)*sizeof(double));

#ifndef MFNOSAFETYNET
    if(Uref==NULL)
     {
      sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",(nu+np+1)*sizeof(double));
      MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

  URref=(double*)realloc((void*)URref,(nu+np+1)*sizeof(double));

#ifndef MFNOSAFETYNET
    if(URref==NULL)
     {
      sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",(nu+np+1)*sizeof(double));
      MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

  ULref=(double*)realloc((void*)ULref,(nu+np+1)*sizeof(double));

#ifndef MFNOSAFETYNET
    if(ULref==NULL)
     {
      sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",(nu+np+1)*sizeof(double));
      MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

  Pref=(double*)realloc((void*)Pref,(nu+np+1)*sizeof(double));

#ifndef MFNOSAFETYNET
    if(Pref==NULL)
     {
      sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",(nu+np+1)*sizeof(double));
      MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

  A=*pA;
  B=*pB;
  C=*pC;
  D=*pD;

  for(i=0;i<(nx-1)*nu*(4*nu+nu+1);i++)A[i]=0.;
  for(i=0;i<(nx-1)*nu*(nbc+nic+k);i++)B[i]=0.;
  for(i=0;i<(nx-1)*nu*(nbc+nic+k);i++)C[i]=0.;
  for(i=0;i<(nbc+nic+k)*(nbc+nic+k);i++)D[i]=0.;

  for(i=0;i<nx-1;i++)
   {
    for(j=0;j<nu;j++)U[j]=.5*(u[j+nu*(i+1)]+u[j+nu*i]);
    for(j=0;j<np;j++)P[j]=u[nx*nu+j];
    t=.5*(u[nx*nu+np+i+1]+u[nx*nu+np+i]);
    for(j=0;j<nu;j++)Uref[j]=.5*(u0[j+nu*(i+1)]+u0[j+nu*i]);
    for(j=0;j<np;j++)Pref[j]=u0[nx*nu+j];

    fu(t,nu,U,np,P,Uref,Pref,fv,e);
    fl(t,nu,U,np,P,Uref,Pref,flv,e);

    for(j=0;j<nu;j++)
     {
      for(J=0;J<nu;J++)
       {
        MFSetBandedMatrixElement(j+nu*i,J+nu*i,A,0.,(nx-1)*nu,nu-1,2*nu-1,e);
        if(j==J)MFIncrementBandedMatrixElement(j+nu*i,J+nu*i,A,-1.,(nx-1)*nu,nu-1,2*nu-1,e);
        MFIncrementBandedMatrixElement(j+nu*i,J+nu*i,A,-.5*(u[nx*nu+np+i+1]-u[nx*nu+np+i])*fv[j+nu*J],(nx-1)*nu,nu-1,2*nu-1,e);
        if(i!=nx-2)
         {
          MFSetBandedMatrixElement(j+nu*i,J+nu*(i+1),A,0.,(nx-1)*nu,nu-1,2*nu-1,e);
          if(j==J)MFIncrementBandedMatrixElement(j+nu*i,J+nu*(i+1),A,1.,(nx-1)*nu,nu-1,2*nu-1,e);
          MFIncrementBandedMatrixElement(j+nu*i,J+nu*(i+1),A,-.5*(u[nx*nu+np+i+1]-u[nx*nu+np+i])*fv[j+nu*J],(nx-1)*nu,nu-1,2*nu-1,e);
          B[j+nu*i+(nx-1)*nu*J]=0.;
         }else{
          if(j==J)B[j+nu*i+(nx-1)*nu*J]=1.;
          B[j+nu*i+(nx-1)*nu*J]+=-.5*(u[nx*nu+np+i+1]-u[nx*nu+np+i])*fv[j+nu*J];
         }
       }
      for(J=0;J<np;J++)
       {
        B[j+nu*i+(nx-1)*nu*(nu+J)]=-(u[nx*nu+np+i+1]-u[nx*nu+np+i])*flv[j+nu*J];
       }
     }
   }

  for(j=0;j<nu;j++)UL[j]=.5*(u[j+nu]+u[j]);
  for(j=0;j<nu;j++)UR[j]=.5*(u[j+nu*(nx-1)]+u[j+nu*(nx-2)]);
  for(j=0;j<np;j++)P[j]=u[nx*nu+j];
  for(j=0;j<nu;j++)ULref[j]=.5*(u0[j+nu]+u0[j]);
  for(j=0;j<nu;j++)URref[j]=.5*(u0[j+nu*(nx-1)]+u0[j+nu*(nx-2)]);
  for(j=0;j<np;j++)Pref[j]=u0[nx*nu+j];
  au(nbc,nu,UL,UR,np,P,ULref,URref,Pref,av,e);
  for(j=0;j<nbc;j++)
   {
    for(i=2*nu;i<nu*(nx-2);i++)C[j+(nbc+nic+k)*i]=0.;
    for(J=0;J<nu;J++)
     {
      C[j+(nbc+nic+k)*J]=.5*av[j+nbc*J];
      C[j+(nbc+nic+k)*(J+nu)]=.5*av[j+nbc*J];
      C[j+(nbc+nic+k)*(J+nu*(nx-2))]=.5*av[j+nbc*(J+nu)];
      D[j+(nbc+nic+k)*J]=.5*av[j+nbc*(J+nu)];
     }
   }
  for(j=0;j<nu;j++)UL[j]=.5*(u[j+nu]+u[j]);
  for(j=0;j<nu;j++)UR[j]=.5*(u[j+nu*(nx-1)]+u[j+nu*(nx-2)]);
  for(j=0;j<np;j++)P[j]=u[nx*nu+j];
  for(j=0;j<nu;j++)ULref[j]=.5*(u0[j+nu]+u0[j]);
  for(j=0;j<nu;j++)URref[j]=.5*(u0[j+nu*(nx-1)]+u0[j+nu*(nx-2)]);
  for(j=0;j<np;j++)Pref[j]=u0[nx*nu+j];
  al(nbc,nu,UL,UR,np,P,ULref,URref,Pref,av,e);
  for(j=0;j<nbc;j++)
   {
    for(J=0;J<np;J++)
     {
      D[j+(nbc+nic+k)*(nu+J)]=av[j+nbc*J];
     }
   }

  if(nic>0)
   {
    for(j=0;j<np;j++)P[j]=u[nx*nu+j];
    for(j=0;j<np;j++)Pref[j]=u0[nx*nu+j];
    ml(nic,np,P,Pref,mv,e);
    for(i=0;i<nx-1;i++)
     {
      for(j=0;j<nu;j++)U[j]=u[j+nu*i];
      for(j=0;j<np;j++)P[nu+j]=u[nx*nu+j];
      for(j=0;j<nu;j++)Uref[j]=u0[j+nu*i];
      for(j=0;j<np;j++)Pref[nu+j]=u0[nx*nu+j];
      t=u[nx*nu+np+i];
      lu(nic,t,nu,U,np,P,Uref,Pref,lv,e);
      ll(nic,t,nu,U,np,P,Uref,Pref,llv,e);
      if(i==0)
       {
        for(j=0;j<nic;j++)
         {
          for(J=0;J<nu;J++)
            C[nbc+j+(nbc+nic+k)*J]=(u[nx*nu+np+1]-u[nx*nu+np+0])*lv[j+nic*J]/2;
          for(J=0;J<np;J++)
            D[nbc+j+(nbc+nic+k)*(nu+J)]=mv[j+nic*J]
                             +(u[nx*nu+np+1]-u[nx*nu+np+0])*llv[j+nic*J]/2;
         }
  
        for(jj=0;jj<nu;jj++)U[jj]=u[jj+nu*(nx-1)];
        t=u[nx*nu+np+nx-1];
        for(jj=0;jj<nu;jj++)Uref[jj]=u[jj+nu*(nx-1)];
        lu(nic,t,nu,U,np,P,Uref,Pref,lv,e);
        ll(nic,t,nu,U,np,P,Uref,Pref,llv,e);
  
        for(j=0;j<nic;j++)
         {
          for(J=0;J<nu;J++)
            D[nbc+j+(nbc+nic+k)*J]=(u[nx*nu+np+nx]-u[nx*nu+np+nx-1])*lv[j+nic*J]/2;
          for(J=0;J<np;J++)
            D[nbc+j+(nbc+nic+k)*(nu+J)]+=
                              (u[nx*nu+np+nx]-u[nx*nu+np+nx-1])*llv[j+nic*J]/2;
         }
       }else{
        for(j=0;j<nic;j++)
         {
          for(J=0;J<nu;J++)
            C[nbc+j+(nbc+nic+k)*(J+nu*i)]=(u[nx*nu+np+i+1]-u[nx*nu+np+i-1])*lv[j+nic*J];
          for(J=0;J<np;J++)
            D[nbc+j+(nbc+nic+k)*(nu+J)]+=(u[nx*nu+np+i+1]-u[nx*nu+np+i-1])*llv[j+nic*J];
         }
       }
     }
   }

/* Arclength */

  for(j=0;j<k;j++)
   {
    o=nx*nu+np+n*j;
    for(jj=0;jj<nu;jj++)
     {
      if(Phi!=NULL)
        C[nbc+nic+j+(nbc+nic+k)*jj]=Phi[jj+n*j]*Phi[o]/2;
       else
        C[nbc+nic+j+(nbc+nic+k)*jj]=0.;
     }
    for(i=1;i<nx-1;i++)
     {
      for(jj=0;jj<nu;jj++)
       {
        if(Phi!=NULL)
          C[nbc+nic+j+(nbc+nic+k)*(jj+nu*i)]=Phi[jj+nu*i+n*j]*Phi[o+i];
         else
          C[nbc+nic+j+(nbc+nic+k)*(jj+nu*i)]=0.;
       }
     }
    for(jj=0;jj<nu;jj++)
     {
      if(Phi!=NULL)
        D[nbc+nic+j+(nbc+nic+k)*jj]=Phi[jj+(nx-1)*nu+n*j]*Phi[o+nx-1]/2;
       else
        D[nbc+nic+j+(nbc+nic+k)*jj]=0.;
     }
    for(i=0;i<np;i++)
     {
      if(Phi!=NULL)
        D[nbc+nic+j+(nbc+nic+k)*(i+nu)]=Phi[i+nx*nu+n*j];
       else
        D[nbc+nic+j+(nbc+nic+k)*(i+nu)]=0.;
     }
   }

/*MFPrintBorderedBandedMatrix(stdout,(nx-1)*nu,nbc+nic+k,nu-1,2*nu-1,A,B,C,D);*/

#ifdef MFALLOWVERBOSE
  if(verbose){printf("done %s\n",RoutineName);fflush(stdout);}
#endif

  return;
 }

double MFTPBVPGetRes(int n, int k,void *d,double *u0,double *u,double *Phi,double **pr,double **ps, MFErrorHandler e)
 {
  static char RoutineName[]={"MFTPBVPGetRes"};
  int nx,nu,np,nic,nbc;
  MFTPBVPFFUNCTION f,fu,fl;
  MFTPBVPAFUNCTION a,au,al;
  MFTPBVPLFUNCTION l,lu,ll;
  MFTPBVPMFUNCTION m,ml;
  int i,j,J,jj,o;
  double error;
  double emax;
  int iemax;
  static double *U=NULL;
  static double *UL=NULL;
  static double *UR=NULL;
  static double *P=NULL;
  static double *Uref=NULL;
  static double *ULref=NULL;
  static double *URref=NULL;
  static double *Pref=NULL;
  static double *fv=NULL;
  static double *av=NULL;
  static double *bv=NULL;
  static double *cv=NULL;
  static double *lv=NULL;
  static double *llv=NULL;
  static double *mv=NULL;
  double *r;
  double *s;
  double t;
  int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s\n",RoutineName);fflush(stdout);}
#endif

  nx=((struct MFTPBVPData*)(d))->nx;
  nu=((struct MFTPBVPData*)(d))->nu;
  np=((struct MFTPBVPData*)(d))->np;
  nbc=((struct MFTPBVPData*)(d))->nbc;
  nic=((struct MFTPBVPData*)(d))->nic;
  k=((struct MFTPBVPData*)(d))->k;
  f=((struct MFTPBVPData*)(d))->f;
  a=((struct MFTPBVPData*)(d))->a;
  l=((struct MFTPBVPData*)(d))->l;
  m=((struct MFTPBVPData*)(d))->m;
  fu=((struct MFTPBVPData*)(d))->fu;
  au=((struct MFTPBVPData*)(d))->au;
  lu=((struct MFTPBVPData*)(d))->lu;
  fl=((struct MFTPBVPData*)(d))->fl;
  al=((struct MFTPBVPData*)(d))->al;
  ll=((struct MFTPBVPData*)(d))->ll;
  ml=((struct MFTPBVPData*)(d))->ml;

  *pr=(double*)realloc((void*)(*pr),nx*nu*sizeof(double));

#ifndef MFNOSAFETYNET
    if(*pr==NULL)
     {
      sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",nx*nu*sizeof(double));
      MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return 1.;
     }
#endif

  *ps=(double*)realloc((void*)(*ps),(nbc+nic+k)*sizeof(double));

#ifndef MFNOSAFETYNET
    if(*ps==NULL)
     {
      sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",(nbc+nic+k)*sizeof(double));
      MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return -1.;
     }
#endif

  fv=(double*)realloc((void*)fv,nu*(nu+np)*sizeof(double));

#ifndef MFNOSAFETYNET
    if(fv==NULL)
     {
      sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",nu*(nu+np)*sizeof(double));
      MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return -1.;
     }
#endif

  av=(double*)realloc((void*)av,nbc*(2*nu+np)*sizeof(double));

#ifndef MFNOSAFETYNET
    if(av==NULL)
     {
      sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",(2*nu+np)*sizeof(double));
      MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return -1.;
     }
#endif

  bv=(double*)realloc((void*)bv,nu*(nu+np)*sizeof(double));

#ifndef MFNOSAFETYNET
    if(bv==NULL)
     {
      sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",nu*(nu+np)*sizeof(double));
      MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return -1.;
     }
#endif

  cv=(double*)realloc((void*)cv,nu*(nu+np)*sizeof(double));

#ifndef MFNOSAFETYNET
    if(cv==NULL)
     {
      sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",nu*(nu+np)*sizeof(double));
      MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return -1.;
     }
#endif

  if(nic>0)
   {
    lv=(double*)realloc((void*)lv,nic*(nu+np)*sizeof(double));

#ifndef MFNOSAFETYNET
    if(lv==NULL)
     {
      sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",nic*(nu+np)*sizeof(double));
      MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return -1.;
     }
#endif

    llv=(double*)realloc((void*)llv,nic*(nu+np)*sizeof(double));

#ifndef MFNOSAFETYNET
    if(llv==NULL)
     {
      sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",nic*(nu+np)*sizeof(double));
      MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return -1.;
     }
#endif

    mv=(double*)realloc((void*)mv,nic*(nu+np)*sizeof(double));

#ifndef MFNOSAFETYNET
    if(mv==NULL)
     {
      sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",nic*(nu+np)*sizeof(double));
      MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return -1.;
     }
#endif

   }
  U=(double*)realloc((void*)U,nu*sizeof(double));

#ifndef MFNOSAFETYNET
    if(U==NULL)
     {
      sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",nu*sizeof(double));
      MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return -1.;
     }
#endif

  UL=(double*)realloc((void*)UL,nu*sizeof(double));

#ifndef MFNOSAFETYNET
    if(UL==NULL)
     {
      sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",nu*sizeof(double));
      MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return -1.;
     }
#endif

  UR=(double*)realloc((void*)UR,nu*sizeof(double));

#ifndef MFNOSAFETYNET
    if(UR==NULL)
     {
      sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",nu*sizeof(double));
      MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return -1.;
     }
#endif

  P=(double*)realloc((void*)P,np*sizeof(double));

#ifndef MFNOSAFETYNET
    if(P==NULL)
     {
      sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",np*sizeof(double));
      MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return -1.;
     }
#endif

  Uref=(double*)realloc((void*)Uref,nu*sizeof(double));

#ifndef MFNOSAFETYNET
    if(Uref==NULL)
     {
      sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",nu*sizeof(double));
      MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return -1.;
     }
#endif

  ULref=(double*)realloc((void*)ULref,nu*sizeof(double));

#ifndef MFNOSAFETYNET
    if(ULref==NULL)
     {
      sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",nu*sizeof(double));
      MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return -1.;
     }
#endif

  URref=(double*)realloc((void*)URref,nu*sizeof(double));

#ifndef MFNOSAFETYNET
    if(URref==NULL)
     {
      sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",nu*sizeof(double));
      MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return -1.;
     }
#endif

  Pref=(double*)realloc((void*)Pref,np*sizeof(double));

#ifndef MFNOSAFETYNET
    if(Pref==NULL)
     {
      sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",np*sizeof(double));
      MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return -1.;
     }
#endif


  r=*pr;
  s=*ps;

  error=0.;
  emax=0.;
  iemax=-1;

/* EQ */

  for(i=0;i<nx-1;i++)
   {

#ifdef MFALLOWVERBOSE
    if(verbose)
     {
      printf("  %d,     u=(%lf,",i,u[0+nu*i]);
      for(j=1;j<nu;j++)printf(",%lf",u[j+nu*i]);
      printf(")   x=%lf\n",u[nx*nu+np+i]);fflush(stdout);
     }
#endif

    for(j=0;j<nu;j++)U[j]=.5*(u[j+nu*(i+1)]+u[j+nu*i]);
    for(j=0;j<np;j++)P[j]=u[nx*nu+j];
    t=.5*(u[nx*nu+np+i+1]+u[nx*nu+np+i]);
    for(j=0;j<nu;j++)Uref[j]=.5*(u0[j+nu*(i+1)]+u0[j+nu*i]);
    for(j=0;j<np;j++)Pref[j]=u0[nx*nu+j];

#ifdef MFALLOWVERBOSE
    if(verbose)
     {
      printf("                     U.%d      u=(%lf,",i,U[0]);
      for(j=1;j<nu;j++)printf(",%lf",U[j]);
      printf(")   x=%lf\n",t);fflush(stdout);
     }
#endif

    f(t,nu,U,np,P,Uref,Pref,fv,e);

    for(j=0;j<nu;j++)
     {
      r[j+nu*i]=u[j+nu*(i+1)]-u[j+nu*i]-(u[nx*nu+np+i+1]-u[nx*nu+np+i])*fv[j];
      error+=r[j+nu*i]*r[j+nu*i];
      if(fabs(r[j+nu*i])>emax)
       {
        emax=fabs(r[j+nu*i]);
        iemax=j+nu*i;
       }

#ifdef MFALLOWVERBOSE
      if(verbose){printf("                     EQ%d.%d, %lf = %lf - %lf - (%lf-%lf)*%lf \n",i,j,r[j+nu*i],u[j+nu*(i+1)],u[j+nu*i],u[nx*nu+np+i+1],u[nx*nu+np+i],fv[j]);fflush(stdout);}
#endif

     }
   }

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("  %d,     u=(%lf,",nx-1,u[0+nu*(nx-1)]);
    for(j=1;j<nu;j++)printf(",%lf",u[j+nu*(nx-1)]);
    printf(")   x=%lf\n",u[nx*nu+np+nx-1]);fflush(stdout);
   }
#endif

  for(j=0;j<nu;j++)UL[j]=.5*(u[j+nu]+u[j]);
  for(j=0;j<nu;j++)UR[j]=.5*(u[j+nu*(nx-1)]+u[j+nu*(nx-2)]);
  for(j=0;j<np;j++)P[j]=u[nx*nu+j];
  for(j=0;j<nu;j++)ULref[j]=.5*(u0[j+nu]+u0[j]);
  for(j=0;j<nu;j++)URref[j]=.5*(u0[j+nu*(nx-1)]+u0[j+nu*(nx-2)]);
  for(j=0;j<np;j++)Pref[j]=u0[nx*nu+j];
  a(nbc,nu,UL,UR,np,P,ULref,URref,Pref,av,e);

/* BC */

  for(j=0;j<nbc;j++)
   {
    s[j]=av[j];

#ifdef MFALLOWVERBOSE
    if(verbose){printf("                     BC%d, %lf\n",j,s[j]);fflush(stdout);}
#endif

    error+=s[j]*s[j];
    if(fabs(s[j])>emax)
     {
      emax=fabs(s[j]);
      iemax=(nx-1)*nu+j;
     }
   }

/* IC */

  if(nic>0)
   {
    m(nic,np,P,Pref,mv,e);
    for(i=0;i<nx-1;i++)
     {
      for(j=0;j<nu;j++)U[j]=u[j+nu*i];
      for(j=0;j<np;j++)P[j]=u[nx*nu+j];
      t=u[nx*nu+np+i];
      for(j=0;j<nu;j++)Uref[j]=u0[j+nu*i];
      for(j=0;j<np;j++)Pref[j]=u0[nx*nu+j];
      l(nic,t,nu,U,np,P,Uref,Pref,lv,e);
      if(i==0)
       {
        m(nic,np,P,Pref,mv,e);
        for(j=0;j<nic;j++)
         s[nbc+j]=mv[j]+lv[j]*(u[nx*nu+np+1]-u[nx*nu+np+0])/2;
  
        for(j=0;j<nu;j++)U[j]=u[j+nu*(nx-1)];
        for(j=0;j<nu;j++)Uref[j]=u0[j+nu*(nx-1)];
        t=u[nx*nu+np+nx-1];
        l(nic,t,nu,U,np,P,Uref,Pref,lv,e);
  
        for(j=0;j<nic;j++)
         s[nbc+j]+=lv[j]*(u[nx*nu+np+nx]-u[nx*nu+np+nx-1])/2;
       }else{
        for(j=0;j<nic;j++)
         s[nbc+j]+=lv[j]*(u[nx*nu+np+i+1]-u[nx*nu+np+i-1]);
       }
     }
    for(j=0;j<nic;j++)
     {
      error+=s[nbc+j]*s[nbc+j];
      if(fabs(s[nbc+j])>emax)
       {
        emax=fabs(s[nbc+j]);
        iemax=(nx-1)*nu+nu+j;
       }

#ifdef MFALLOWVERBOSE
      if(verbose){printf("                     IC%d, %le\n",j,s[nbc+j]);fflush(stdout);}
#endif

     }
   }

/* Arclength */

  for(j=0;j<k;j++)
   {
    o=nx*nu+np+n*j;
    s[nbc+nic+j]=0.;
    if(Phi!=NULL)
     {
      for(jj=0;jj<nu;jj++)
        s[nbc+nic+j]+=Phi[jj+n*j]*(u[jj]-u0[jj])*Phi[o]/2;
      for(i=1;i<nx-1;i++)
       {
        for(jj=0;jj<nu;jj++)
         {
          s[nbc+nic+j]+=Phi[jj+nu*i+n*j]*(u[jj+nu*i]-u0[jj+nu*i])*Phi[o+i];
         }
       }
      for(jj=0;jj<nu;jj++)
        s[nbc+nic+j]+=Phi[jj+nu*(nx-1)+n*j]*(u[jj+nu*(nx-1)]-u0[jj+nu*(nx-1)])*Phi[o+nx-1]/2;
      for(i=0;i<np;i++)
       {
        s[nbc+nic+j]+=Phi[nx*nu+i+n*j]*(u[nx*nu+i]-u0[nx*nu+i]);
       }
      error+=s[nbc+nic+j]*s[nbc+nic+j];
      if(fabs(s[nbc+nic+j])>emax)
       {
        emax=fabs(s[nbc+nic+j]);
        iemax=(nx-1)*nu+nbc+nic+j;
       }
     }

#ifdef MFALLOWVERBOSE
    if(verbose){printf("                     TC%d, %lf\n",j,s[nbc+nic+j]);fflush(stdout);}
#endif

   }

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("Maximum error %lf at component %d\n",emax,iemax);
    printf("done %s\n",RoutineName);fflush(stdout);
   }
#endif

  return sqrt(error);
 }

int MFTangentTPBVPWithGuess(int n,int k,MFNVector u,MFNKMatrix Phi1,MFNKMatrix Phi,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFTangentTPBVPWithGuess"};
  int rc;

  rc=MFTangentTPBVP(n,k,u,Phi,d,e);

  return rc;
 }

double MFScaleTPBVP(int n,int k,MFNVector vu,MFNKMatrix mPhi,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFScaleTPBVP"};
  int i,j;
  double R;
  double ev,evmax;
  int info;
  char jobvl;
  char jobvr;
  static double *a=NULL;
  int lda;
  static double *wr=NULL;
  static double *wi=NULL;
  int ldvl;
  static double *vl=NULL;
  int ldvr;
  static double *vr=NULL;
  static double *work=NULL;
  int lwork;
  double *u, *Phi;
  int verbose=0;

#ifdef HAVE_LAPACK

  Phi=MFNKM_CStar(mPhi,e);
  u=MFNV_CStar(vu,e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s\n",RoutineName);fflush(stdout);}
#endif

  lda=k;
  a=(double*)realloc((void*)a,k*k*sizeof(double));

#ifndef MFNOSAFETYNET
  if(a==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",k*k*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return -1.;
   }
#endif

  wr=(double*)realloc((void*)wr,k*sizeof(double));

#ifndef MFNOSAFETYNET
  if(wr==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",k*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return -1.;
   }
#endif

  wi=(double*)realloc((void*)wi,k*sizeof(double));

#ifndef MFNOSAFETYNET
  if(wi==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",k*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return -1.;
   }
#endif

  ldvl=k;
  vl=(double*)realloc((void*)vl,k*k*sizeof(double));

#ifndef MFNOSAFETYNET
  if(vl==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",k*k*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return -1.;
   }
#endif

  ldvr=k;
  vr=(double*)realloc((void*)vr,k*k*sizeof(double));

#ifndef MFNOSAFETYNET
  if(vr==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",k*k*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return -1.;
   }
#endif

  lwork=3*k;
  work=(double*)realloc((void*)work,lwork*sizeof(double));

#ifndef MFNOSAFETYNET
  if(work==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",lwork*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return -1.;
   }
#endif

  for(i=0;i<k;i++)
   for(j=0;j<k;j++)
    {
     a[i+k*j]=MFTPBVPCurvatureSingle(n,i,j,u,u,Phi,d,e);

#ifdef MFALLOWVERBOSE
     if(verbose){printf("a[%d,%d]=%lf\n",i,j,a[i+k*j]);fflush(stdout);}
#endif

    }

  info=0;
  jobvl='N';
  jobvr='N';

  CALLDGEEV(&jobvl,&jobvr,&k,a,&lda,wr,wi,vl,&ldvl,vr,&ldvr,work,&lwork,&info);

  evmax=sqrt(wr[0]*wr[0]+wi[0]*wi[0]);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("|ev[0]|=%lf\n",evmax);fflush(stdout);}
#endif

  for(i=1;i<k;i++)
   {
    ev=sqrt(wr[i]*wr[i]+wi[i]*wi[i]);

#ifdef MFALLOWVERBOSE
    if(verbose){printf("|ev[%d]|=%lf\n",i,ev);fflush(stdout);}
#endif

    if(ev>evmax)evmax=ev;
   }

  R=sqrt(2*MFEpsilon/evmax);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("R=%lf\n",R);fflush(stdout);}
#endif

  if(R!=R)R=1.;

/* R=.1;
  if(R<.05)R=.05;
  R=3.*R;*/
  return R;
#else
  sprintf(MFTPBVPMFErrorHandlerMsg,"The TPBVP manifold requires dgesvd from LAPACK.");
  MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
  return -1.;
#endif
 }

int MFSolveFull(int n,double *A,double *b, MFErrorHandler e)
 {
  static char RoutineName[]={"MFSolveFull"};
  int ierr;
  static int *Pivots=NULL;
  static int nPivots=0;
  char trans='N';
  int one=1;

#ifdef HAVE_LAPACK

  if(n>nPivots)
   {
    Pivots=(int*)realloc((void*)Pivots,n*sizeof(int));

#ifndef MFNOSAFETYNET
    if(Pivots==NULL)
     {
      sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(double));
      MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return -1;
     }
#endif

    nPivots=n;
   }

  CALLDGETRF(&n,&n,A,&n,Pivots,&ierr);

  if(ierr!=0)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Problem with factor, zero on diagonal %d\n",ierr);
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    return 0;
   }
  ierr=0;

  CALLDGETRS(&trans,&n,&one,A,&n,Pivots,b,&n,&ierr);

  return 1;
#else
  sprintf(MFTPBVPMFErrorHandlerMsg,"MFTPBVP requires LAPACK");
  MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
  return 0;
#endif
 }

void MFPrintBorderedBandedMatrixMinusFull(FILE *fid, int n1,int n2,int nBandsL,int nBandsU,double *A,double *B,double *C,double *D,double *L, MFErrorHandler e)
 {
  static char RoutineName[]={"MFPrintBorderedBandedMatrixMinusFull"};
  int i,j;
  int ldA=0;
  int nBands=0;

  ldA=2*nBandsL+nBandsU+1;
  nBands=nBandsL+nBandsU+1;

  fprintf(fid,"    ");
  for(i=0;i<n1;i++)fprintf(fid," %4d ",i);
  for(i=0;i<n2;i++)fprintf(fid," %4d ",i+n1);
  fprintf(fid,"\n");fflush(fid);
  for(i=0;i<n1;i++)
   {
    fprintf(fid,"%2d [",i);
    for(j=0;j<n1;j++)
     {
      if(j>0)fprintf(fid," ");
      if(j<i-nBandsL||j>i+nBandsU)fprintf(fid,"%6.0le",fabs(-L[i+(n1+n2)*j]));
       else fprintf(fid,"%6.0le",fabs(A[i-j+nBands-1+ldA*j]-L[i+(n1+n2)*j]));
     }
    fprintf(fid,"|");
    for(j=0;j<n2;j++)
     {
      if(j>0)fprintf(fid," ");
      fprintf(fid,"%6.0le",fabs(B[i+n1*j]-L[i+(n1+n2)*(n1+j)]));
     }
    fprintf(fid,"]\n");
   }

  fprintf(fid,"   [");
  for(j=0;j<n1;j++)
   {
    if(j>0)fprintf(fid,"-");
    fprintf(fid,"-----");
   }
  fprintf(fid,"+");
  for(j=0;j<n2;j++)
   {
    if(j>0)fprintf(fid,"-");
    fprintf(fid,"-----");
   }
  fprintf(fid,"]\n");

  for(i=0;i<n2;i++)
   {
    fprintf(fid,"%2d [",i+n1);
    for(j=0;j<n1;j++)
     {
      if(j>0)fprintf(fid," ");
      fprintf(fid,"%6.0le",fabs(C[i+n2*j]-L[i+n1+(n1+n2)*j]));
     }
    fprintf(fid,"|");
    for(j=0;j<n2;j++)
     {
      if(j>0)fprintf(fid," ");
      fprintf(fid,"%6.0le",fabs(D[i+n2*j]-L[i+n1+(n1+n2)*(j+n1)]));
     }
    fprintf(fid,"]\n");
   }
  fflush(fid);

  return;
 }

int MFTPBVPGetNX(MFImplicitMF M, MFErrorHandler e)
 {
  static char RoutineName[]={"MFTPBVPGetNX"};

  return ((struct MFTPBVPData*)(MFIMFGetData(M,e)))->nx;
 }

int MFTPBVPGetNU(MFImplicitMF M, MFErrorHandler e)
 {
  static char RoutineName[]={"MFTPBVPGetNU"};

  return ((struct MFTPBVPData*)(MFIMFGetData(M,e)))->nu;
 }

int MFTPBVPGetNP(MFImplicitMF M, MFErrorHandler e)
 {
  static char RoutineName[]={"MFTPBVPGetNP"};

  return ((struct MFTPBVPData*)(MFIMFGetData(M,e)))->np;
 }

int MFTPBVPGetNIC(MFImplicitMF M, MFErrorHandler e)
 {
  static char RoutineName[]={"MFTPBVPGetNIC"};

  return ((struct MFTPBVPData*)(MFIMFGetData(M,e)))->nic;
 }

int MFTPBVPGetNBC(MFImplicitMF M, MFErrorHandler e)
 {
  static char RoutineName[]={"MFTPBVPGetNBC"};

  return ((struct MFTPBVPData*)(MFIMFGetData(M,e)))->nbc;
 }

/*! \fn MFNVector MFTPBVPIntegrateForInitialSolution(MFImplicitMF M,double *u0,double *p,double *x);
 *  \brief Solves an initial value problem in place of a MFTPBVPMF, and returns the solution. This may be useful in
 *         constructing initial guesses.
 *
 *  \param M An MFTPBVPMF
 *  \param u0 An array of length nu with the initial condition.
 *  \param p  An array of length nl with the parameter values.
 *  \param x  The nx+1 mesh points on [0,1].
 *  \returns A solution (u(t),p).
 */
MFNVector MFTPBVPIntegrateForInitialSolution(MFImplicitMF M, double *u0,double *p0, double *x, MFErrorHandler e)
 {
  static char RoutineName[]={"MFTPBVPIntegrateForInitialSolution"};
  int i,j,l;
  int nx,nu,np;
  MFNVector result;
  double error,err,z;
  int itimes;
  double t;
  double *res;
  static double *u=NULL;
  static double *up=NULL;
  static double *U=NULL;
  static double *f=NULL;
  static double *fp=NULL;
  static double *fu=NULL;
  static double *b=NULL;
  static double *A=NULL;
  struct MFTPBVPData *data;

  nx=MFTPBVPGetNX(M,e);
  nu=MFTPBVPGetNU(M,e);
  np=MFTPBVPGetNP(M,e);

  u=(double*)realloc((void*)u,nu*sizeof(double));

#ifndef MFNOSAFETYNET
    if(u==NULL)
     {
      sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",nu*sizeof(double));
      MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return NULL;
     }
#endif

  up=(double*)realloc((void*)up,nu*sizeof(double));

#ifndef MFNOSAFETYNET
    if(up==NULL)
     {
      sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",nu*sizeof(double));
      MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return NULL;
     }
#endif

  U=(double*)realloc((void*)U,nu*sizeof(double));

#ifndef MFNOSAFETYNET
    if(U==NULL)
     {
      sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",nu*sizeof(double));
      MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return NULL;
     }
#endif

  f=(double*)realloc((void*)f,nu*sizeof(double));

#ifndef MFNOSAFETYNET
    if(f==NULL)
     {
      sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",nu*sizeof(double));
      MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
     }
#endif

  fp=(double*)realloc((void*)fp,nu*sizeof(double));

#ifndef MFNOSAFETYNET
    if(fp==NULL)
     {
      sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",nu*sizeof(double));
      MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return NULL;
     }
#endif

  fu=(double*)realloc((void*)fu,nu*nu*sizeof(double));

#ifndef MFNOSAFETYNET
    if(fu==NULL)
     {
      sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",nu*sizeof(double));
      MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return NULL;
     }
#endif

  b=(double*)realloc((void*)b,nu*sizeof(double));

#ifndef MFNOSAFETYNET
    if(b==NULL)
     {
      sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",nu*sizeof(double));
      MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return NULL;
     }
#endif

  A=(double*)realloc((void*)A,nu*nu*sizeof(double));

#ifndef MFNOSAFETYNET
    if(A==NULL)
     {
      sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",nu*nu*sizeof(double));
      MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return NULL;
     }
#endif


  data=(struct MFTPBVPData*)MFIMFGetData(M,e);
  t=0.;
  (data->f)(t,nu,u0,np,p0,u0,p0,f,e);

  result=MFIMFVectorFactory(M,e);
  res=MFNV_CStar(result,e);

  for(j=0;j<nu;j++)up[j]=u0[j]-(x[1]-x[0])*f[j]/2;
  for(j=0;j<nu;j++)res[j]=up[j];
  for(j=0;j<nu;j++)up[j]=u0[j]+(x[1]-x[0])*f[j]/2;
  for(j=0;j<nu;j++)res[nu+j]=up[j];
  for(j=0;j<nu;j++)u[j]=up[j];

  for(i=0;i<nx-1;i++)
   {
/*  (u[i+1]-u[i])-(x[i+1]-x[i])*f((u[i+1]+u[i])/2)=0
i=1 (u[1]-u[0])-(x[1]-x[0])*f((u[1]+u[0])/2)=0        */

    t=.5*(x[i+1]+x[i]);
    itimes=0;
    error=1.;
    while(error>1.e-10 && itimes<50)
     {
      for(j=0;j<nu;j++)U[j]=.5*(u[j]+up[j]);

      (data->f)(t,nu,U,np,p0,U,p0,f,e);
      (data->fu)(t,nu,U,np,p0,U,p0,fu,e);

      error=0.;
      for(j=0;j<nu;j++)
       {
        b[j]=-(u[j]-up[j])+(x[i+1]-x[i])*f[j];
        error+=b[j]*b[j];
        for(l=0;l<nu;l++)
         {
          if(j==l)A[j+nu*l]=1.-(x[i+1]-x[i])*fu[j+nu*l];
           else A[j+nu*l]=-(x[i+1]-x[i])*fu[j+nu*l];
         }
       }
      error=sqrt(error);

      MFSolveFull(nu,A,b,e);

      for(j=0;j<nu;j++)u[j]=u[j]+b[j];
      itimes++;
     }
    for(j=0;j<nu;j++)
     {
      MFNVSetC(result,j+nu*(i+1),up[j],e);
      up[j]=u[j];
     }
   }
  for(j=0;j<np;j++)MFNVSetC(result,j+nx*nu,p0[j],e);
  for(j=0;j<nx;j++)MFNVSetC(result,j+nx*nu+np,x[j],e);

  MFNVSetIndex(result,0,e);

  return result;
 }

double MFTPBVPCurvatureSingle(int n,int it,int jt, double *u,double *u0,double *Phi,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFTPBVPCurvatureSingle"};
  static double *A=NULL;
  static double *B=NULL;
  static double *C=NULL;
  static double *D=NULL;
  static double *r=NULL;
  static double *s=NULL;
  static double *r1=NULL;
  static double *s1=NULL;
  static double *du=NULL;
  double result;
  int k,nx,nu,np,nbc,nic;
  int i,j;
  double err;
  static int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s\n",RoutineName);fflush(stdout);}
#endif

  nx=((struct MFTPBVPData*)(d))->nx;
  nu=((struct MFTPBVPData*)(d))->nu;
  np=((struct MFTPBVPData*)(d))->np;
  nbc=((struct MFTPBVPData*)(d))->nbc;
  nic=((struct MFTPBVPData*)(d))->nic;
  k=((struct MFTPBVPData*)(d))->k;

  r=(double*)realloc((void*)r,((nx-1)*nu+nu+np)*sizeof(double));

#ifndef MFNOSAFETYNET
  if(r==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",((nx-1)*nu+nu+np)*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return -1.;
   }
#endif

  s=(double*)realloc((void*)s,(nbc+nic+k)*sizeof(double));

#ifndef MFNOSAFETYNET
  if(s==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",(nbc+nic+k)*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return -1.;
   }
#endif

  r1=(double*)realloc((void*)r1,((nx-1)*nu+nu+np)*sizeof(double));

#ifndef MFNOSAFETYNET
  if(r1==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",((nx-1)*nu+nu+np)*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return -1.;
   }
#endif

  s1=(double*)realloc((void*)s1,(nbc+nic+k)*sizeof(double));

#ifndef MFNOSAFETYNET
  if(s1==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",(nbc+nic+k)*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return -1.;
   }
#endif

  du=(double*)realloc((void*)du,n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(du==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return -1.;
   }
#endif


  MFTPBVPGetJac(n,k,d,u,u0,Phi,&A,&B,&C,&D,e);

#ifdef MFALLOWVERBOSE
  if(0&&verbose)
   {
    int l;

    MFPrintBorderedBandedMatrixByBlock(stdout,nx,nu,np,nbc,nic,k,A,B,C,D,e);
    for(i=0;i<nx;i++)
     {
      for(l=0;l<k;l++)
       {
        if(i==0&&l==0)printf("      Phi =[(%le",Phi[0+nu*i+n*l]);
         else printf("           (%le",Phi[0+nu*i+n*l]);
        for(j=1;j<nu;j++)printf(",%le",Phi[j+nu*i+n*l]);printf(")");
       }
       printf("]\n");fflush(stdout);
     }
    for(l=0;l<k;l++)
     {
      if(l==0)printf("           [(%le",Phi[0+nu*nx+n*l]);
       else printf("           (%le",Phi[0+nu*nx+n*l]);
      for(j=1;j<np;j++)printf(",%le",Phi[j+nu*nx+n*l]);printf(")");
     }
     printf("]\n");fflush(stdout);
   }
#endif


  err=1.e-5;
  for(i=nx*nu+np;i<n;i++)du[i]=u[i];
  for(i=0;i<nx*nu+np;i++)du[i]=u[i]+Phi[i+k*it]*err/2+Phi[i+k*jt]*err/2;
  MFTPBVPGetRes(n,k,d,u,du,Phi,&r1,&s1,e);
  for(i=0;i<(nx-1)*nu;i++)r[i]=r1[i];
  for(i=0;i<nbc+nic+k;i++)s[i]=s1[i];

  for(i=0;i<nx*nu+np;i++)du[i]=u[i]+Phi[i+k*it]*err/2-Phi[i+k*jt]*err/2;
  MFTPBVPGetRes(n,k,d,u,du,Phi,&r1,&s1,e);
  for(i=0;i<(nx-1)*nu;i++)r[i]-=r1[i];
  for(i=0;i<nbc+nic+k;i++)s[i]-=s1[i];

  for(i=0;i<nx*nu+np;i++)du[i]=u[i]-Phi[i+k*it]*err/2+Phi[i+k*jt]*err/2;
  MFTPBVPGetRes(n,k,d,u,du,Phi,&r1,&s1,e);
  for(i=0;i<(nx-1)*nu;i++)r[i]-=r1[i];
  for(i=0;i<nbc+nic+k;i++)s[i]-=s1[i];

  for(i=0;i<nx*nu+np;i++)du[i]=u[i]-Phi[i+k*it]*err/2-Phi[i+k*jt]*err/2;
  MFTPBVPGetRes(n,k,d,u,du,Phi,&r1,&s1,e);
  for(i=0;i<(nx-1)*nu;i++)r[i]=-(r[i]+r1[i])/err/err;
  for(i=0;i<nbc+nic+k;i++)s[i]=-(s[i]+s1[i])/err/err;

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("  Rhs= - GuuPhiPhi\n");
    for(i=0;i<nx-1;i++)
     {
      if(i==0)printf("  Rhs    r=(%le",r[0+nu*i]);
       else printf("           (%le",r[0+nu*i]);
      for(j=1;j<nu;j++)printf(",%le",r[j+nu*i]);printf(")\n");fflush(stdout);
     }
    printf("         s=(%le",s[0]);for(j=1;j<nu;j++)printf(",%le",s[j]);printf(")\n");fflush(stdout);
    printf("         p=(%le",s[nu]);for(j=nu+1;j<nbc+nic+k;j++)printf(",%le",s[j]);printf(")\n");fflush(stdout);
   }
#endif

  MFSolveBordered((nx-1)*nu,nbc+nic+k,nu-1,2*nu-1,A,B,C,D,r,s,NULL,NULL,e);

/* Extract the result */

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    for(i=0;i<nx-1;i++)
     {
      if(i==0)printf("         r=(%le",r[0+nu*i]);
       else printf("           (%le",r[0+nu*i]);
      for(j=1;j<nu;j++)printf(",%le",r[j+nu*i]);printf(")\n");fflush(stdout);
     }
    printf("         s=(%le",s[0]);for(j=1;j<nu;j++)printf(",%le",s[j]);printf(")\n");fflush(stdout);
    printf("         p=(%le",s[nu]);for(j=nu+1;j<nbc+nic+k;j++)printf(",%le",s[j]);printf(")\n");fflush(stdout);
   }
#endif

  result=0;
  for(j=0;j<nu;j++)result+=r[j+nu*0]*r[j+nu*0]*u[nx*nu+np+0]/2.;
  for(i=1;i<nx-1;i++)for(j=0;j<nu;j++)result+=r[j+nu*i]*r[j+nu*i]*u[nx*nu+np+i];
  for(i=0;i<nu;i++)result+=s[i]*s[i]*u[nx*nu+np+nx-1]/2.;
  for(i=nu;i<nbc+nic+k;i++)result+=s[i]*s[i];
  result=sqrt(result);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("done %s, result=%lf\n",RoutineName,result);fflush(stdout);}
#endif

  return result;
 }

/*! \fn MFNVector MFTPBVPIntegrateForTangent(MFImplicitMF M,MFNVector u,double *du0,double *dp);
 *  \brief Solves an initial value problem of the linearization of a MFTPBVPMF, and returns the solution. This may be useful in
 *         constructing initial approximations of the columns of the basis for the tangent space.
 *
 *  \param M An MFTPBVPMF
 *  \param u A solution (u(t),p) which is the "point" at which the variational equations are written.
 *  \param du0 The initial perturbation (at x=0).
 *  \param dp The perturbation of the parameters.
 *  \returns A solution (du(t),dp) that might be used as a basis vector for the tangent space.
 */
MFNVector MFTPBVPIntegrateForTangent(MFImplicitMF M, MFNVector u, double *du0,double *dp, MFErrorHandler e)
 {
  static char RoutineName[]={"MFTPBVPIntegrateForTangent"};
  int i,j,l;
  int n,nx,nu,np;
  MFNVector result;
  double t;
  static double *a=NULL;
  static double *b=NULL;
  static double *U=NULL;
  static double *P=NULL;
  static double *fu=NULL;
  static double *fl=NULL;
  double dot;
  double x,xm;
  struct MFTPBVPData *data;

  data=(struct MFTPBVPData*)MFIMFGetData(M,e);

  nx=MFTPBVPGetNX(M,e);
  nu=MFTPBVPGetNU(M,e);
  np=MFTPBVPGetNP(M,e);

/*  (du[i+1]-du[i])-(x[i+1]-x[i])*fu((u[i+1]+u[i])/2)(du[i+1]+du[i])/2
               -sum (x[i+1]-x[i])*fl*dp=0 */

  result=MFCreateNVector(MFIMF_N(M,e),e);
  for(j=0;j<nu;j++)MFNVSetC(result,j,du0[j],e);
  for(j=0;j<np;j++)MFNVSetC(result,nx*nu+j,dp[j],e);
  for(j=0;j<nx;j++)MFNVSetC(result,nx*nu+np+j,MFNV_C(u,nx*nu+np+j,e),e);

  U=(double*)realloc((void*)U,nu*sizeof(double));

#ifndef MFNOSAFETYNET
  if(U==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",nu*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  P=(double*)realloc((void*)P,np*sizeof(double));

#ifndef MFNOSAFETYNET
  if(P==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",np*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  a=(double*)realloc((void*)a,nu*nu*sizeof(double));

#ifndef MFNOSAFETYNET
  if(a==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",nu*nu*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  b=(double*)realloc((void*)b,nu*sizeof(double));

#ifndef MFNOSAFETYNET
  if(b==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",nu*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
   }
#endif

  fu=(double*)realloc((void*)fu,nu*nu*sizeof(double));

#ifndef MFNOSAFETYNET
  if(fu==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",nu*nu*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  fl=(double*)realloc((void*)fl,nu*np*sizeof(double));

#ifndef MFNOSAFETYNET
  if(fl==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",nu*np*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  for(i=1;i<nx;i++)
   {
    x=MFNV_C(u,nx*nu+np+i,e);
    xm=MFNV_C(u,nx*nu+np+i-1,e);

    for(j=0;j<nu;j++)U[j]=(MFNV_C(u,nu*i+j,e)+MFNV_C(u,nu*(i-1)+j,e))/2;
    t=(x+xm)/2;
    for(j=0;j<np;j++)P[j]=MFNV_C(u,nx*nu+j,e);

    (data->fu)(t,nu,U,np,P,U,P,fu,e);
    (data->fl)(t,nu,U,np,P,U,P,fl,e);

    for(j=0;j<nu;j++)
     {
      b[j]=MFNV_C(result,nu*(i-1)+j,e);
      for(l=0;l<nu;l++)
       {
        if(j==l)a[j+nu*l]=1.;
         else a[j+nu*l]=0.;
        a[j+nu*l]-=(x-xm)*fu[j+nu*l]/2;
        b[j]+=(x-xm)*fu[j+nu*l]*MFNV_C(result,nu*(i-1)+l,e)/2;
       }
      for(l=0;l<nu;l++)b[j]+=(x-xm)*fl[l]*MFNV_C(result,nx*nu+l,e);
     }
    MFSolveFull(nu,a,b,e);
    for(j=0;j<nu;j++)MFNVSetC(result,nu*i+j,b[j],e);
   }

  return result;
 }

double MFTPBVPTestTangent(MFImplicitMF M,MFNVector u, MFNKMatrix Phi, MFErrorHandler e)
 {
  static char RoutineName[]={"MFTPBVPTestTangent"};
  static double *A=NULL;
  static double *B=NULL;
  static double *C=NULL;
  static double *D=NULL;
  static double *r=NULL;
  static double *s=NULL;
  static double *uD=NULL;
  static double *PhiD=NULL;
  int nx,nu,np,nbc,nic,N;
  int i,j,k,l;
  int n,n1,n2,nBands,nBandsU,nBandsL,ldA;
  double error,rowsum;

  nx=MFTPBVPGetNX(M,e);
  nu=MFTPBVPGetNU(M,e);
  np=MFTPBVPGetNP(M,e);
  nbc=MFTPBVPGetNBC(M,e);
  nic=MFTPBVPGetNIC(M,e);
  k=MFIMF_K(M,e);

  r=(double*)realloc((void*)r,((nx-1)*nu+nu+np)*sizeof(double));

#ifndef MFNOSAFETYNET
  if(r==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",((nx-1)*nu+np)*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return -1.;
   }
#endif

  s=(double*)realloc((void*)s,(nbc+nic+k)*sizeof(double));

#ifndef MFNOSAFETYNET
  if(s==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",(nbc+nic+k)*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return -1.;
   }
#endif

  n=nu*nx+np;
  N=MFIMF_K(M,e);
  n1=(nx-1)*nu;
  n2=nbc+nic+k;
  nBandsL=nu-1;
  nBandsU=2*nu-1;
  ldA=2*nBandsL+nBandsU+1;
  nBands=nBandsL+nBandsU+1;

  printf("%s\n",RoutineName);fflush(stdout);

  uD=MFNV_CStar(u,e);
  PhiD=MFNKM_CStar(Phi,e);

  MFTPBVPGetJac(n,k,MFIMFGetData(M,e),uD,uD,PhiD,&A,&B,&C,&D,e);

  error=0.;
  for(l=0;l<k;l++)
   {
    for(i=0;i<n1;i++)
     {
      rowsum=0.;
      for(j=0;j<n1;j++)
       {
        if(!(j<i-nBandsL||j>i+nBandsU))
          rowsum+=A[i-j+nBands-1+ldA*j]*PhiD[j+N*l];
       }
      for(j=0;j<n2;j++)
        rowsum+=B[i+n1*j]*PhiD[n1+j+N*l];
      error+=fabs(rowsum);
      if(fabs(rowsum)>1.e-4)printf(" error in [%d,%d] is %le\n",i,l,rowsum);
     }
    for(i=0;i<n2;i++)
     {
      rowsum=0;
      for(j=0;j<n1;j++)
        rowsum+=C[i+n2*j]*PhiD[j+N*l];
      for(j=0;j<n2;j++)
        rowsum+=C[i+n2*j]*PhiD[n1+j+N*l];
      if(i<n1-k || i+k-n1-n2!=l)
       {
        error+=fabs(rowsum);
        if(fabs(rowsum)>1.e-4)printf(" error in [%d,%d] is %le\n",i,l,rowsum);
       }else{
        error+=fabs(rowsum-1.);
        if(fabs(rowsum-1)>1.e-4)printf(" error in [%d,%d] is %le\n",i,l,rowsum);
       }
     }
   }


  return error;
 }

void MFTPBVPEvaluateIntegralConstraints(MFImplicitMF M, MFNVector uv, MFNVector u0v, double *s, MFErrorHandler e)
 {
  static char RoutineName[]={"MFTPBVPEvaluateIntegralConstraints"};
  int nx,nu,np,nic,nbc,k;
  double t;
  double *u;
  double *u0;
  MFTPBVPLFUNCTION l;
  MFTPBVPMFUNCTION m;
  static double *U=NULL;
  static double *P=NULL;
  static double *Uref=NULL;
  static double *Pref=NULL;
  static double *lv=NULL;
  static double *mv=NULL;
  void *d;
  int i,j;

  d=MFIMFGetData(M,e);

  nx=((struct MFTPBVPData*)(d))->nx;
  nu=((struct MFTPBVPData*)(d))->nu;
  nbc=((struct MFTPBVPData*)(d))->nbc;
  np=((struct MFTPBVPData*)(d))->np;
  nic=((struct MFTPBVPData*)(d))->nic;
  k=((struct MFTPBVPData*)(d))->k;
  l=((struct MFTPBVPData*)(d))->l;
  m=((struct MFTPBVPData*)(d))->m;

  if(nic>0)
   {
    lv=(double*)realloc((void*)lv,nic*(nu+np)*sizeof(double));

#ifndef MFNOSAFETYNET
    if(lv==NULL)
     {
      sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",nic*(nu+np)*sizeof(double));
      MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

    mv=(double*)realloc((void*)mv,nic*(nu+np)*sizeof(double));

#ifndef MFNOSAFETYNET
    if(mv==NULL)
     {
      sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",nic*(nu+np)*sizeof(double));
      MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif
   }

  U=(double*)realloc((void*)U,nu*sizeof(double));

#ifndef MFNOSAFETYNET
  if(U==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",nu*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  P=(double*)realloc((void*)P,np*sizeof(double));

#ifndef MFNOSAFETYNET
  if(P==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",np*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  Uref=(double*)realloc((void*)Uref,nu*sizeof(double));

#ifndef MFNOSAFETYNET
  if(Uref==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",nu*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  Pref=(double*)realloc((void*)Pref,(nu+np+1)*sizeof(double));

#ifndef MFNOSAFETYNET
  if(Pref==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",(nu+np+1)*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  u=MFNV_CStar(uv,e);
  u0=MFNV_CStar(u0v,e);

  for(j=0;j<np;j++)P[j]=u[nx*nu+j];
  for(j=0;j<np;j++)Pref[j]=u0[nx*nu+j];

  m(nic,np,P,Pref,mv,e);
  for(i=0;i<nx-1;i++)
   {
    for(j=0;j<nu;j++)U[j]=.5*(u[j+nu*(i+1)]+u[j+nu*i]);
    for(j=0;j<np;j++)P[j]=u[nx*nu+j];
    t=.5*(u[nx*nu+np+i+1]+u[nx*nu+np+i]);
    for(j=0;j<nu;j++)Uref[j]=.5*(u0[j+nu*(i+1)]+u0[j+nu*i]);
    for(j=0;j<np;j++)Pref[j]=u0[nx*nu+j];
    l(nic,t,nu,U,np,P,Uref,Pref,lv,e);
    if(i==0)
     {
      m(nic,np,P,Pref,mv,e);
      for(j=0;j<nic;j++)
       s[j]=mv[j]+lv[j]*(u[nx*nu+np+1]-u[nx*nu+np+0])/2;

      for(j=0;j<nu;j++)U[j]=u[j+nu*(nx-1)];
      for(j=0;j<nu;j++)Uref[j]=u0[j+nu*(nx-1)];
      t=u[nx*nu+np+nx-1];
      l(nic,t,nu,U,np,P,Uref,Pref,lv,e);

      for(j=0;j<nic;j++)
       s[j]+=lv[j]*(u[nx*nu+np+nx]-u[nx*nu+np+nx-1])/2;
     }else{
      for(j=0;j<nic;j++)
       s[j]+=lv[j]*(u[nx*nu+np+i+1]-u[nx*nu+np+i-1])/2;
     }
   }

  return;
 }

void MFTPBVPEvaluateBoundaryConditions(MFImplicitMF M, MFNVector uv, MFNVector u0v, double *s, MFErrorHandler e)
 {
  static char RoutineName[]={"MFTPBVPEvaluateBoundaryConditions"};
  int nx,nu,np,nic,nbc,k;
  double t;
  double *u;
  double *u0;
  MFTPBVPAFUNCTION a;
  static double *UL=NULL;
  static double *UR=NULL;
  static double *P=NULL;
  static double *ULref=NULL;
  static double *URref=NULL;
  static double *Pref=NULL;
  void *d;
  int j;

  d=MFIMFGetData(M,e);

  nx=((struct MFTPBVPData*)(d))->nx;
  nu=((struct MFTPBVPData*)(d))->nu;
  nbc=((struct MFTPBVPData*)(d))->nbc;
  np=((struct MFTPBVPData*)(d))->np;
  nic=((struct MFTPBVPData*)(d))->nic;
  k=((struct MFTPBVPData*)(d))->k;
  a=((struct MFTPBVPData*)(d))->a;

  UL=(double*)realloc((void*)UL,nu*sizeof(double));

#ifndef MFNOSAFETYNET
  if(UL==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",nu*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  UR=(double*)realloc((void*)UR,nu*sizeof(double));

#ifndef MFNOSAFETYNET
  if(UR==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",nu*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  P=(double*)realloc((void*)P,np*sizeof(double));

#ifndef MFNOSAFETYNET
  if(P==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",np*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  ULref=(double*)realloc((void*)ULref,nu*sizeof(double));

#ifndef MFNOSAFETYNET
  if(ULref==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",nu*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  URref=(double*)realloc((void*)URref,nu*sizeof(double));

#ifndef MFNOSAFETYNET
  if(URref==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",nu*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  Pref=(double*)realloc((void*)Pref,(nu+np+1)*sizeof(double));

#ifndef MFNOSAFETYNET
  if(Pref==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",(nu+np+1)*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  u=MFNV_CStar(uv,e);
  u0=MFNV_CStar(u0v,e);

  for(j=0;j<nu;j++)UL[j]=.5*(u[j+nu]+u[j]);
  for(j=0;j<nu;j++)UR[j]=.5*(u[j+nu*(nx-1)]+u[j+nu*(nx-2)]);
  for(j=0;j<np;j++)P[j]=u[nx*nu+j];
  for(j=0;j<nu;j++)ULref[j]=.5*(u0[j+nu]+u0[j]);
  for(j=0;j<nu;j++)URref[j]=.5*(u0[j+nu*(nx-1)]+u0[j+nu*(nx-2)]);
  for(j=0;j<np;j++)Pref[j]=u0[nx*nu+j];
  a(nbc,nu,UL,UR,np,P,ULref,URref,Pref,s,e);

  return;
 }

int MFTPBVPAnalyze(MFImplicitMF M,MFNVector u0, MFErrorHandler e)
 {
  static char RoutineName[]={"MFTPBVPAnalyze"};
  static double *A=NULL;
  static double *B=NULL;
  static double *C=NULL;
  static double *D=NULL;
  int nx,nu,np,nbc,nic;
  int n,k;
  int result;

  nx=MFTPBVPGetNX(M,e);
  nu=MFTPBVPGetNU(M,e);
  np=MFTPBVPGetNP(M,e);
  nbc=MFTPBVPGetNBC(M,e);
  nic=MFTPBVPGetNIC(M,e);
  n=MFIMF_N(M,e);
  k=MFIMF_K(M,e);

  printf("%s, Raw k is %d, Full rank is %d\n",RoutineName,k,nu+np-nbc-nic);fflush(stdout);

  k=(nx*nu+np)-((nx-1)*nu+nbc+nic);
  MFTPBVPGetJac(n,k,MFIMFGetData(M,e),MFNV_CStar(u0,e),MFNV_CStar(u0,e),NULL,&A,&B,&C,&D,e);
  result=MFAnalyzeBordered(k,(nx-1)*nu,nbc+nic+k,nu-1,2*nu-1,A,B,C,D,e);
  printf(" Dimension of NS is %d\n",result);fflush(stdout);

  return result;
 }

int MFTPBVPTestJacobian(int nx, int nu, int np, int nbc, int nic, int n, int k,double *u0,double *Phi, void *data, MFErrorHandler e)
 {
  static char RoutineName[]={"MFTPBVPTestJacobian"};
  static double *A=NULL;
  static double *B=NULL;
  static double *C=NULL;
  static double *D=NULL;
  static double *r=NULL;
  static double *s=NULL;
  static double *r2=NULL;
  static double *s2=NULL;
  static double *u=NULL;
  static double *approxA=NULL;
  int i,j,J;
  double error,erow,delta;
  double epsilon=1.e-5;

  int verbose=0;

  approxA=(double*)realloc((void*)approxA,(nx*nu+np)*((nx-1)*nu+nbc+nic+k)*sizeof(double));

#ifndef MFNOSAFETYNET
  if(approxA==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",(nx*nu+np+1)*((nx-1)*nu+nbc+nic+k)*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0;
   }
#endif

  u=(double*)realloc((void*)u,n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(u==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0;
   }
#endif


#ifdef MFALLOWVERBOSE
  if(verbose)printf("%s nx=%d, nu=%d, np=%d, nbc=%d, nic=%d\n",RoutineName,nx,nu,np,nbc,nic);fflush(stdout);
#endif

  for(i=0;i<n;i++)u[i]=u0[i];

  error=MFTPBVPGetRes(n,k,data,u0,u,Phi,&r,&s,e);
  for(j=0;j<nx*nu+np;j++)
   {
    u[j]=u[j]+epsilon;
    MFTPBVPGetRes(n,k,data,u0,u,Phi,&r2,&s2,e);
    u[j]=u[j]-epsilon;
    error=0.;
    for(i=0;i<(nx-1)*nu;i++)
     {
      approxA[i+((nx-1)*nu+nbc+nic+k)*j]=(r2[i]-r[i])*1.e5;
     }
    for(i=0;i<nbc+nic+k;i++)
     {
      approxA[i+(nx-1)*nu+((nx-1)*nu+nbc+nic+k)*j]=(s2[i]-s[i])*1.e5;
     }
   }
  for(j=nx*nu+np;j<(nx-1)*nu+nbc+nic+k;j++)
   {
    for(i=0;i<(nx-1)*nu+nbc+nic+k;i++)
      approxA[i+((nx-1)*nu+nbc+nic+k)*j]=0.;
   }

  MFTPBVPGetJac(n,k,data,u0,u,Phi,&A,&B,&C,&D,e);

  printf("Analytic Jacobian\n");
  MFPrintBorderedBandedMatrix(stdout,(nx-1)*nu,nbc+nic+k,nu-1,2*nu-1,A,B,C,D,e);

  error=0.;
  for(i=0;i<(nx-1)*nu;i++)
   {
    erow=0.;
    for(j=0;j<(nx-1)*nu;j++)
     {
      if(j<i-(nu-1)||j>i+(2*nu-1))erow+=fabs(-approxA[i+((nx-1)*nu+nbc+nic+k)*j]);
       else erow+=fabs(A[i-j+((nu-1)+(2*nu-1)+1)-1+(2*(nu-1)+(2*nu-1)+1)*j]-approxA[i+((nx-1)*nu+nbc+nic+k)*j]);
     }
    for(j=0;j<nbc+nic+k;j++)
      erow+=fabs(B[i+(nx-1)*nu*j]-approxA[i+((nx-1)*nu+nbc+nic+k)*((nx-1)*nu+j)]);
    error+=erow;
    if(erow>10*epsilon)
     {
      printf("Large error in row %d of Jac %le (Equation %d at point %d)\n",i,erow,i%nu,i-i%nu);
      printf("    exact     approx   error\n");
      for(j=0;j<(nx-1)*nu;j++)
       {
        if(j<i-(nu-1)||j>i+(2*nu-1))printf("%3d%10.3le%10.3le%10.3le\n",j,0.,approxA[i+((nx-1)*nu+nbc+nic+k)*j],fabs(-approxA[i+((nx-1)*nu+nbc+nic+k)*j]));
         else printf("%3d%10.3le%10.3le%10.3le\n",j,
                A[i-j+((nu-1)+(2*nu-1)+1)-1+(2*(nu-1)+(2*nu-1)+1)*j],
                                                                     approxA[i+((nx-1)*nu+nbc+nic+k)*j],
           fabs(A[i-j+((nu-1)+(2*nu-1)+1)-1+(2*(nu-1)+(2*nu-1)+1)*j]-approxA[i+((nx-1)*nu+nbc+nic+k)*j]));
       }
      for(j=0;j<nbc+nic+k;j++)
        printf("%3d%10.3le%10.3le%10.3le\n",j+(nx-1)*nu,B[i+(nx-1)*nu*j],approxA[i+((nx-1)*nu+nbc+nic+k)*((nx-1)*nu+j)],fabs(B[i+(nx-1)*nu*j]-approxA[i+((nx-1)*nu+nbc+nic+k)*((nx-1)*nu+j)]));
     }
   }

  for(i=0;i<nbc+nic+k;i++)
   {
    erow=0.;
    for(j=0;j<(nx-1)*nu;j++)
      erow+=fabs(C[i+(nbc+nic+k)*j]-approxA[i+(nx-1)*nu+((nx-1)*nu+nbc+nic+k)*j]);
    for(j=0;j<nbc+nic+k;j++)
      erow+=fabs(D[i+(nbc+nic+k)*j]-approxA[i+(nx-1)*nu+((nx-1)*nu+nbc+nic+k)*(j+(nx-1)*nu)]);
    error+=erow;
/*  if(i>nbc+nic-1)*/
     {
      if(erow>10*epsilon)
       {
        printf("Large error in row %d of Jac %le",i+(nx-1)*nu,erow);
        if(i<nbc)printf(" (boundary condition %d)\n",i);
         else if(i<nbc+nic)printf(" (integral constraint %d)\n",i-nbc);
         else printf(" (arclength constraint %d)\n",i-nbc-nic);
        printf("    exact     approx   error\n");
        for(j=0;j<(nx-1)*nu;j++)
          printf("%3d%10.3le%10.3le%10.3le\n",j,C[i+(nbc+nic+k)*j],approxA[i+(nx-1)*nu+((nx-1)*nu+nbc+nic+k)*j],fabs(C[i+(nbc+nic+k)*j]-approxA[i+(nx-1)*nu+((nx-1)*nu+nbc+nic+k)*j]));
        for(j=0;j<nbc+nic+k;j++)
          printf("%3d%10.3le%10.3le%10.3le\n",j+(nx-1)*nu,
                   D[i+(nbc+nic+k)*j],
                                    approxA[i+(nx-1)*nu+((nx-1)*nu+nbc+nic+k)*(j+(nx-1)*nu)],
              fabs(D[i+(nbc+nic+k)*j]-approxA[i+(nx-1)*nu+((nx-1)*nu+nbc+nic+k)*(j+(nx-1)*nu)]));
       }
     }
   }

  printf("Error in Jacobian %le\n",error);

  return 1;
 }

int MFAnalyzeBordered(int k, int n1,int n2,int nl,int nu,double *A,double *B,double *C,double *D, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAnalyzeBordered"};
  int i,j;
  int ierr;
  static double *X=NULL;
  static double *Y=NULL;
  static int *Pivots=NULL;
  static int ldA,ldX;
  int result;

#ifdef HAVE_LAPACK

  ldA=2*nl+nu+1;
  ldX=n2+nl+1;

  X=(double*)realloc((void*)X,ldX*ldX*sizeof(double));

#ifndef MFNOSAFETYNET
  if(X==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",ldX*ldX*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0;
   }
#endif

  Y=(double*)realloc((void*)Y,(ldA+n2)*sizeof(double));

#ifndef MFNOSAFETYNET
  if(Y==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",(ldA+n2)*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0;
   }
#endif

  Pivots=(int*)realloc((void*)Pivots,(n1+n2)*sizeof(int));

#ifndef MFNOSAFETYNET
  if(Pivots==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",(n1+n2)*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0;
   }
#endif

  CALLDBOFA(A,&ldA,&n1,&nl,&nu,B,&n1,&n2,C,&n2,D,&n2,X,NULL,Pivots,&ierr);

#else
    sprintf(MFTPBVPMFErrorHandlerMsg,"The TPBVP manifold requires dgesvd from LAPACK.");
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    return 0;
#endif

  return result;
 }

int MFTangentTPBVP(int n,int k,MFNVector vu,MFNKMatrix mPhi,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFTangentTPBVP"};
  static double *A=NULL;
  static double *B=NULL;
  static double *C=NULL;
  static double *D=NULL;
  static double *X=NULL;
  static double *vr=NULL;
  static double *s=NULL;
  static double *work=NULL;
  static int    *pvt=NULL;
  int lwork;
  char jobvl;
  char jobvr;
  double t;
  int nx,nu,np,nbc,nic;
  int i,j;
  static int ldA,ldX;
  int ierr;
  int n1,n2,ml,mu;
  double *u,*Phi;
  int verbose=0;

#ifdef HAVE_LAPACK

  Phi=MFNKM_CStar(mPhi,e);
  u=MFNV_CStar(vu,e);

  nx=((struct MFTPBVPData*)(d))->nx;
  nu=((struct MFTPBVPData*)(d))->nu;
  np=((struct MFTPBVPData*)(d))->np;
  nbc=((struct MFTPBVPData*)(d))->nbc;
  nic=((struct MFTPBVPData*)(d))->nic;
  k=((struct MFTPBVPData*)(d))->k;

  n1=(nx-1)*nu;
  n2=nbc+nic+k;
  ml=nu-1;
  mu=2*nu-1;
  ldA=2*ml+mu+1;
  ldX=n2+ml+1;

  ldX=2*nu+np;

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("%s nx=%d, nu=%d, np=%d, nbc=%d, nic=%d\n",RoutineName,nx,nu,np,nbc,nic);fflush(stdout);
    printf("         n1=%d, n2=%d\n",n1,n2);fflush(stdout);
    printf("         ml=%d, mu=%d\n",ml,mu);fflush(stdout);
    printf("         u=");MFPrintNVector(stdout,vu,e);printf("\n");fflush(stdout);
   }
#endif

  i=-1;
  ierr=0;
  jobvl='N';
  jobvr='A';
  CALLDGESVD(&jobvl,&jobvr,&n2,&n2,NULL,&n2,NULL,NULL,&n2,NULL,&n2,&t,&i,&ierr);
  lwork=round(t);

  X=(double*)realloc((void*)X,ldX*ldX*sizeof(double));

#ifndef MFNOSAFETYNET
  if(X==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",ldX*ldX*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0;
   }
#endif

  vr=(double*)realloc((void*)vr,ldX*ldX*sizeof(double));

#ifndef MFNOSAFETYNET
  if(vr==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",ldX*ldX*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0;
   }
#endif

  s=(double*)realloc((void*)s,ldX*sizeof(double));

#ifndef MFNOSAFETYNET
  if(s==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",ldX*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0;
   }
#endif

  work=(double*)realloc((void*)work,lwork*sizeof(double));

#ifndef MFNOSAFETYNET
  if(work==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",lwork*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0;
   }
#endif

  pvt=(int*)realloc((void*)pvt,n*sizeof(int));

#ifndef MFNOSAFETYNET
  if(pvt==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0;
   }
#endif

  MFTPBVPGetJac(n,k,d,u,u,NULL,&A,&B,&C,&D,e);

  CALLDBOSVD(A,&ldA,&n1,&ml,&mu,B,&n1,&n2,C,&n2,D,&n2,X,s,vr,work,&lwork,pvt,&ierr);

  j=0;
  for(i=0;i<ldX;i++)
    if(fabs(s[i])<1.e-10)j++;

  if(j>k)
   {
    printf("%s failed, too many tangents (%d)\n",j);fflush(stdout);
    return 0;
   }

  for(i=0;i<k;i++)
    CALLDBONV(&i,A,&ldA,&n1,&ml,&mu,B,&n1,&n2,C,&n2,D,&n2,X,s,vr,Phi+n*i);

#ifdef MFALLOWVERBOSE
  if(0&&verbose)
   {
    for(j=0;j<nx*nu+np;j++)
     {
      if(j==0)printf("Tangent [");
        else  printf("        [");
      for(i=0;i<k;i++)
       {
        if(i>0)printf(" ");
        printf("%lf",Phi[j+n*i]);
       }
      printf("]\n");fflush(stdout);
     }
   }
#endif
  
/* Gram Schmidt */

  MFTPBVPGramSchmidtNoMat(((struct MFTPBVPData*)(d))->space,n,nx,nu,np,k,Phi,e);

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("done %s\n",RoutineName);fflush(stdout);
    printf("         u=");MFPrintNVector(stdout,vu,e);printf("\n");fflush(stdout);
   }
#endif

  return 1;
#else
    sprintf(MFTPBVPMFErrorHandlerMsg,"The TPBVP manifold requires dgesvd from LAPACK.");
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    return 0;
#endif
 }

void MFTPBVPTestTangentRaw(int n,int k,double *u,double *Phi,void *d, MFErrorHandler err)
 {
  static char RoutineName[]={"MFTPBVPTestTangentRaw"};
  static double *A=NULL;
  static double *B=NULL;
  static double *C=NULL;
  static double *D=NULL;
  int nx,nu,np,nbc,nic;
  int i,j,l;
  static int ldA,ldX;
  int ierr;
  int n1,n2,ml,mu;
  double e,error,maxerror;

  printf("%s\n",RoutineName);

  nx=((struct MFTPBVPData*)(d))->nx;
  nu=((struct MFTPBVPData*)(d))->nu;
  np=((struct MFTPBVPData*)(d))->np;
  nbc=((struct MFTPBVPData*)(d))->nbc;
  nic=((struct MFTPBVPData*)(d))->nic;
  k=((struct MFTPBVPData*)(d))->k;

  n1=(nx-1)*nu;
  n2=nbc+nic+k;
  ml=nu-1;
  mu=2*nu-1;
  ldA=2*ml+mu+1;
  ldX=n2+ml+1;

  MFTPBVPGetJac(n,k,d,u,u,Phi,&A,&B,&C,&D,err);
  for(l=0;l<k;l++)
   {
    error=0.;
    for(i=0;i<n1;i++)
     {
      if(Phi[i+n*l]!=Phi[i+n*l])printf("%s  Coordinate %d of tangent %d is NaN\n",RoutineName,i,l);
      e=0.;
      for(j=0;j<n1;j++)e+=Phi[j+n*l]*MFGetBandedMatrixElement(i,j,A,n,ml,mu,err);
      for(j=0;j<n2;j++)e+=Phi[j+n1+n*l]*B[i+n1*j];
      e=fabs(e);
      if(e>1.e-10)printf("Large error (%le)  in row %d of tangent %d\n",e,i,l);
      error+=e;
      if(i==0||e>maxerror)maxerror=e;
     }
    for(i=0;i<n2;i++)
     {
      if(Phi[i+n1+n*l]!=Phi[i+n1+n*l])printf("  Coordinate %d of tangent %d is NaN\n",n1+i,l);
      e=0.;
      for(j=0;j<n1;j++)e+=Phi[j+n*l]*C[i+n2*j];
      for(j=0;j<n2;j++)e+=Phi[j+n1+n*l]*D[i+n2*j];
      if(i-n2+k==l)e=e-1.;
      e=fabs(e);
      if(i>=n2-k)printf(" row, tangent [%d,%d] is %le\n",i-n2+k,l,err);
      if(e>1.e-10)printf("Large error (%le)  in row %d+%d of tangent %d\n",e,n1,i,l);
      error+=e;
      if(i==0||e>maxerror)maxerror=e;
     }
    printf("Average error in tangent %d is %le\n",l,fabs(error)/n1);
    printf("Maximum error in tangent %d is %le\n",l,maxerror);
   }

  return;
 }

void MFTestSolveBordered(int nbc, int nic, int n1,int n2,int ml,int mu,double *A,double *B,double *C,double *D,double *R,double *S, double *x, double *y, MFErrorHandler err)
 {
  static char RoutineName[]={"MFTestSolveBordered"};
  int i,j,l;
  int imax;
  double e,error,maxerror;

  printf("Test Solve Bordered\n");
  error=0.;
  for(i=0;i<n1;i++)
   {
    e=-R[i];
    for(j=0;j<n1;j++)e+=x[j]*MFGetBandedMatrixElement(i,j,A,n1,ml,mu,err);
    for(j=0;j<n2;j++)e+=y[j]*B[i+n1*j];
    e=fabs(e);
    error+=e;
    if(i==0||e>maxerror){maxerror=e;imax=i;}
   }
  for(i=0;i<n2;i++)
   {
    e=-S[i];
    for(j=0;j<n1;j++)e+=x[j]*C[i+n2*j];
    for(j=0;j<n2;j++)e+=y[j]*D[i+n2*j];
    e=fabs(e);
    error+=e;
    if(e>maxerror){maxerror=e;imax=i+n1;}
   }
  printf("Average error in SolveBordered is %le\n",fabs(error)/n1);
  printf("Maximum error in SolveBordered is %le, row %d\n",maxerror,imax);
  printf("Aux var is %le\n",y[n2-1]);
  printf("Diff Eq's are equations [%d,%d]\n",0,n1-1);
  printf("BC's are equations [%d,%d]\n",n1,n1+nbc-1);
  if(nic>0)printf("IC's are equations [%d,%d]\n",n1+nbc,n1+nbc+nic-1);
  printf("ALC's are equations [%d,%d]\n",n1+nbc+nic,n1+n2-1);

  return;
 }

void MFTPBVPGramSchmidtNoMat(MFNSpace space, int n, int nx, int nu, int np, int k, double *Phi, MFErrorHandler e)
 {
  static char RoutineName[]={"MFGramSchmidtNoMat"};
  double inner;
  int i,j,jj,l,o;
  int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s, n=%d, nx=%d, nu=%d, np=%d, n should be %d\n",RoutineName,n,nx,nu,np,nx*nu+np+nx);fflush(stdout);}
#endif

  o=nx*nu+np;
  for(i=0;i<k;i++)
   {

#ifdef MFALLOWVERBOSE
    if(verbose){printf("The %dth vector\n",i);fflush(stdout);}
#endif

    inner=0.;
    for(l=0;l<nu;l++)
     inner+=Phi[l+n*i]*Phi[l+n*i]*Phi[o+n*i]/2;

#ifdef MFALLOWVERBOSE
    if(verbose){printf("    dx_0=(%lf-%lf)/2=%lf, avg=%lf, Phi[%d] and Phi[%d]\n",Phi[o+1+n*i],Phi[o+n*i],(Phi[o+1+n*i]-Phi[o+n*i])/2,(Phi[o+1+n*i]+Phi[o+n*i])/2,n*i+o+1,n*i+o);fflush(stdout);}
#endif

    for(j=1;j<nx-1;j++)
     {
      for(l=0;l<nu;l++)
       inner+=Phi[l+j*nu+n*i]*Phi[l+j*nu+n*i]*Phi[o+j+n*i];

#ifdef MFALLOWVERBOSE
      if(verbose){printf("    dx_%d=(%lf-%lf)=%lf, avg=%lf, Phi[%d] and Phi[%d]\n",j,Phi[o+j+1+n*i],Phi[o+j+n*i],(Phi[o+j+1+n*i]-Phi[o+j+n*i]),(Phi[o+j+1+n*i]+Phi[o+j+n*i])/2,n*i+o+j+1,n*i+o+j);fflush(stdout);}
#endif

     }
    for(l=0;l<nu;l++)
     {
      inner+=Phi[l+(nx-1)*nu+n*i]*Phi[l+(nx-1)*nu+n*i]*Phi[o+nx-1+n*i]/2;

#ifdef MFALLOWVERBOSE
      if(verbose){printf("    dx_%d=(%lf-%lf)/2=%lf, avg=%lf, Phi[%d] and Phi[%d]\n",nx-2,Phi[o+nx+n*i],Phi[o+nx-1+n*i],(Phi[o+nx+n*i]-Phi[o+nx-1+n*i])/2,(Phi[o+nx+n*i]+Phi[o+nx-1+n*i])/2,n*i+o+nx,n*i+o+nx-1);fflush(stdout);}
#endif

      }
    for(l=0;l<np;l++)inner+=Phi[nx*nu+l+n*i]*Phi[nx*nu+l+n*i];
    inner=1./sqrt(inner);

#ifdef MFALLOWVERBOSE
    if(verbose){printf("    Norm =%lf\n",1./inner);fflush(stdout);}
#endif

    for(l=0;l<nx*nu+np;l++)Phi[l+n*i]=Phi[l+n*i]*inner;

#ifdef MFALLOWVERBOSE
    if(verbose)
     {
      printf("The %dth vector, scaled\n",i);
      for(l=0;l<nx*nu+np;l++)printf("  %d %le\n",l,Phi[l+n*i]);
      fflush(stdout);
     }
#endif

    for(j=i+1;j<k;j++)
     {
      inner=0.;
      for(l=0;l<nu;l++)
       inner+=Phi[l+n*i]*Phi[l+n*j]*Phi[o+n*i]/2;
      for(jj=1;jj<nx-1;jj++)
        for(l=0;l<nu;l++)
         inner+=Phi[l+jj*nu+n*i]*Phi[l+jj*nu+n*j]*Phi[o+jj+n*i];
      for(l=0;l<nu;l++)
       inner+=Phi[l+(nx-1)*nu+n*i]*Phi[l+(nx-1)*nu+n*j]*Phi[o+nx-1+n*i]/2;
      for(l=0;l<np;l++)inner+=Phi[nx*nu+l+n*i]*Phi[nx*nu+l+n*j];

      for(l=0;l<nx*nu+np;l++)Phi[l+n*j]=Phi[l+n*j]-inner*Phi[l+n*i];
     }
   }

  return;
 }

int MFTPBVPProjectToSave(MFNVector u, double *x, void *d, MFErrorHandler e)
 {
  struct MFTPBVPData *data;
  int i;

  data=(struct MFTPBVPData*)d;

  if(x==NULL)return data->nu+data->np;
  for(i=0;i<data->nu;i++)x[i]=MFNV_C(u,i,e);
  for(i=0;i<data->np;i++)x[data->nu+i]=MFNV_C(u,data->nx*data->nu+i,e);

  return 0;
 }

int MFTPBVPProjectToDraw(MFNVector u, double *x, void *d, MFErrorHandler e)
 {
  struct MFTPBVPData *data;
  int i,j;

  data=(struct MFTPBVPData*)d;

  if(x==NULL)return 3;

  if(data->np==2)
   {
    x[0]=0.;
    for(i=0;i<data->nx;i++)
     {
      for(j=0;j<data->nu;j++)
       {
        x[0]+=(MFNV_C(u,i*data->nu+j,e)+MFNV_C(u,(i+1)*data->nu+j,e))*(MFNV_C(u,data->nx*data->nu+data->np+i+1,e)-MFNV_C(u,data->nx*data->nu+data->np+i,e))/2.;
       }
     }
    x[0]=x[0]/data->nu;
    x[1]=MFNV_C(u,data->nx*data->nu,e);
    x[2]=MFNV_C(u,data->nx*data->nu+1,e);
   }else{
    x[0]=MFNV_C(u,data->nx*data->nu,e);
    x[1]=MFNV_C(u,data->nx*data->nu+1,e);
    x[2]=MFNV_C(u,data->nx*data->nu+2,e);
   }

  return 0;
 }

int MFTPBVPProjectForBB(MFNVector u, double *x, void *d, MFErrorHandler e)
 {
  struct MFTPBVPData *data;
  int i,j;

  data=(struct MFTPBVPData*)d;

#ifdef MFUSEPARAMETERSONLY
  if(x==NULL)return data->np;
  for(i=0;i<data->np;i++)x[i]=MFNV_C(u,data->nx*data->nu+i,e);
#else
  if(x==NULL)return data->nu+data->np;
  for(i=0;i<data->nu;i++)
   {
    x[i]=0.;
    for(j=0;j<data->nx;j++)
       x[i]+=MFNV_C(u,j+i*(data->nu),e)*MFNV_C(u,j+i*(data->nu),e);
    x[i]=x[i]/data->nx;
   }
  for(i=0;i<data->np;i++)
    x[data->nu+i]=MFNV_C(u,data->nx*data->nu+i,e);
#endif

  return 0;
 }

void MFTPBVPSetStability(MFImplicitMF M,MFNVector vU, MFNKMatrix mPhi,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFTPBVPSetStability"};
  static double *A=NULL;
  static double *B=NULL;
  static double *C=NULL;
  static double *D=NULL;
  static double *X=NULL;
  static double *X0=NULL;
  static double *r=NULL;
  static double *s=NULL;
  struct MFTPBVPData *data;
  int ldX;
  int i;
  int n,k;
  int nx,nu,np,nbc,nic;
  int result;
  int verbose=0;
  double *Phi;
  double *u;

#ifdef HAVE_LAPACK

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s, u=0x%8.8x\n",RoutineName,vU);fflush(stdout);}
#endif

  n=MFIMF_N(M,e);
  k=MFIMF_K(M,e);
  Phi=MFNKM_CStar(mPhi,e);
  u=MFNV_CStar(vU,e);

  nx=((struct MFTPBVPData*)(d))->nx;
  nu=((struct MFTPBVPData*)(d))->nu;
  np=((struct MFTPBVPData*)(d))->np;
  nbc=((struct MFTPBVPData*)(d))->nbc;
  nic=((struct MFTPBVPData*)(d))->nic;

  MFTPBVPGetJac(n,k,d,u,u,Phi,&A,&B,&C,&D,e);
  MFTPBVPGetRes(n,k,d,u,u,NULL,&r,&s,e);
  ldX=nbc+nic+k+nu;

  X=(double*)realloc((void*)X,ldX*ldX*sizeof(double));

#ifndef MFNOSAFETYNET
  if(X==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",ldX*ldX*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  X0=(double*)realloc((void*)X0,ldX*ldX*sizeof(double));

#ifndef MFNOSAFETYNET
  if(X0==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",ldX*ldX*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  MFSolveBordered((nx-1)*nu,nbc+nic+k,nu-1,2*nu-1,A,B,C,D,r,s,X,X0,e);

#ifdef MFALLOWVERBOSE
  if(0&&verbose){
   int i,j;
   printf(" X is \n");
   for(i=0;i<nbc+nic+k+nu;i++)
    {
     if(i==nu||i==nu+nbc)printf("-------------------------------------------------\n");
     j=0;
     printf(" [ %10.7lf",X[i+(nbc+nic+k+nu)*j]);
     for(j=1;j<nbc+nic+k+nu;j++)
      {
       if(j==nu||j==2*nu)printf(" |");
       printf(" %10.7lf",X[i+(nbc+nic+k+nu)*j]);
      }
     printf("]\n");fflush(stdout);
    }
  }
#endif

  if(0)
   {
    static double *A=NULL;
    static double *B=NULL;
    static double *C=NULL;
    static double *D=NULL;
    static double *Y=NULL;
    static int *Pivots=NULL;
    int nrhs;
    int ierr;
    int i,j,l;
    char trans='N';

    A=(double*)realloc((void*)A,nu* nu*sizeof(double));

#ifndef MFNOSAFETYNET
    if(A==NULL)
     {
      sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",nu*sizeof(double));
        MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

    B=(double*)realloc((void*)B,nu*(nu+np)*sizeof(double));

#ifndef MFNOSAFETYNET
    if(B==NULL)
     {
      sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",nu*(nu+np)*sizeof(double));
      MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

    C=(double*)realloc((void*)C,(nu+np)*nu*sizeof(double));

#ifndef MFNOSAFETYNET
    if(C==NULL)
     {
      sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",(nu+np)*nu*sizeof(double));
      MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

    D=(double*)realloc((void*)D,(nu+np)*(nu+np)*sizeof(double));

#ifndef MFNOSAFETYNET
    if(D==NULL)
     {
      sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",(nu+np)*(nu+np)*sizeof(double));
      MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

    Y =(double*)realloc((void*)Y,(nu+np)*(nu+np)*sizeof(double));

#ifndef MFNOSAFETYNET
    if(Y==NULL)
     {
      sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",(nu+np)*(nu+np)*sizeof(double));
      MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

    Pivots=(int*)realloc((void*)Pivots,nu*sizeof(int));

#ifndef MFNOSAFETYNET
    if(Pivots==NULL)
     {
      sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",nu*sizeof(int));
      MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

    for(i=0;i<nu;i++)
     {
      for(j=0;j<nu;j++)A[i+nu*j]=X0[i+(nbc+nic+k+nu)*j];
      for(j=0;j<nu+np;j++)B[i+nu*j]=X0[i+(nbc+nic+k+nu)*(j+nu)];
     }

    for(i=0;i<nu+np;i++)
     {
      for(j=0;j<nu;j++)C[i+(nu+np)*j]=X0[i+nu+(nbc+nic+k+nu)*j];
      for(j=0;j<nu+np;j++)D[i+(nu+np)*j]=X0[i+nu+(nbc+nic+k+nu)*(j+nu)];
     }

    CALLDGETRF(&nu,&nu,A,&nu,Pivots,&ierr);

    nrhs=nu+np;
    CALLDGETRS(&trans,&nu,&nrhs,A,&nu,Pivots,B,&nu,&ierr);

    for(i=0;i<nu+np;i++)
     {
      for(j=0;j<nu+np;j++)
       {
        Y[i+(nu+np)*j]=D[i+(nu+np)*j];
        for(l=0;l<nu;l++)
          Y[i+(nu+np)*j]-=C[i+(nu+np)*l]*B[l+nu*j];
       }
     }

#ifdef MFALLOWVERBOSE
    if(verbose){printf("parms=(%lf,%lf)\n",u[nx*nu],u[nx*nu+1]);}
#endif

    result=MFAtlasDetSign(nbc+nic+k,Y,e);

#ifdef MFALLOWVERBOSE
    if(verbose){printf(" stability index by det of Y, %d\n",result);}
#endif

   }else{
    result=MFAtlasDetSign(nu+nbc+nic+k,X0,e);

#ifdef MFALLOWVERBOSE
    if(verbose)printf(" stability index by det of X0, %d\n",result);
#endif

   }

  MFNVSetIndex(vU,result,e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("index is %d\n",result);fflush(stdout);}
#endif

  return;
#else
  sprintf(MFTPBVPMFErrorHandlerMsg,"The TPBVP manifold requires dgetrf and dgetrs from LAPACK.");
  MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
  return;
#endif
 }

/*! \fn void MFTPBVPSetEpsilon(MFImplicitMF M,double epsilon);
 * \brief Sets the tolerance on the distance between a linear approximation at the center of a chart and the 
 *        manifold.
 */
void MFTPBVPSetEpsilon(MFImplicitMF thisBVP,double e, MFErrorHandler err)
 {
  epsilon=fabs(e);
  return;
 }

int MFStopTPBVP(MFImplicitMF thisBVP,MFNVector u0,MFNKMatrix Phi0,MFNVector u1,MFNKMatrix Phi1,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFStopTPBVP"};
  int verbose=0;
  int n0,n1;
  int s;
  static double *prod=NULL;

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("In %s\n",RoutineName);fflush(stdout);
    printf("%s index 0(0x%8.8x)=%d, index 1(0x%8.8x)=%d\n",RoutineName,u0,MFNVGetIndex(u0,e),u1,MFNVGetIndex(u1,e));fflush(stdout);
   }
#endif

  n0=MFNVGetIndex(u0,e);
  n1=MFNVGetIndex(u1,e);
  if(n0==-97||n1==-97)
   {
    if(1){printf("%s, n0=%d, (index %d) n1=%d, (index %d)\n",RoutineName,n0,MFNVGetIndex2(u0,e),n1,MFNVGetIndex2(u1,e));fflush(stdout);}
    return 0;
   }

  prod=(double*)realloc((void*)prod,MFNKMatrixK(Phi0,e)*MFNKMatrixK(Phi1,e)*sizeof(double));

#ifndef MFNOSAFETYNET
  if(prod==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",MFNKMatrixK(Phi0,e)*MFNKMatrixK(Phi1,e)*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
    return 0;
   }
#endif

  MFMMMul(((struct MFTPBVPData*)d)->space,Phi0,Phi1,prod,e);
  s=MFAtlasDetSign(MFNKMatrixK(Phi0,e),prod,e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s, n0=%d, (type %d) n1=%d, (type %d) dot of TS=%d\n",RoutineName,n0,MFNVGetIndex2(u0,e),n1,MFNVGetIndex2(u1,e),s);fflush(stdout);}
#endif

  if(s<0.)n1=-n1;

  if(n0!=n1)return 1;
   else return 0;
 }

int MFSingularTPBVP(int n,int k,MFNVector vu,MFNKMatrix mPhi,MFNVector perp,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFSingularTPBVP"};
  static double *A=NULL;
  static double *B=NULL;
  static double *C=NULL;
  static double *D=NULL;
  static double *X=NULL;
  static double *vr=NULL;
  static double *s=NULL;
  static double *work=NULL;
  static int    *pvt=NULL;
  int lwork;
  char jobvl;
  char jobvr;
  double t;
  int nx,nu,np,nbc,nic;
  int i,j;
  static int ldA,ldX;
  int ierr;
  int n1,n2,ml,mu;
  double *u,*Phi,*p;
  int nullDim;
  int verbose=0;

#ifdef HAVE_LAPACK

  Phi=MFNKM_CStar(mPhi,e);
  u=MFNV_CStar(vu,e);
  if(perp!=NULL)p=MFNV_CStar(perp,e);
   else p=NULL;

  nx=((struct MFTPBVPData*)(d))->nx;
  nu=((struct MFTPBVPData*)(d))->nu;
  np=((struct MFTPBVPData*)(d))->np;
  nbc=((struct MFTPBVPData*)(d))->nbc;
  nic=((struct MFTPBVPData*)(d))->nic;
  k=((struct MFTPBVPData*)(d))->k;

  n1=(nx-1)*nu;
  n2=nbc+nic+k;
  ml=nu-1;
  mu=2*nu-1;
  ldA=2*ml+mu+1;
  ldX=n2+ml+1;

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("%s nx=%0, nu=%d, np=%d, nbc=%d, nic=%d\n",RoutineName,nx,nu,np,nbc,nic);fflush(stdout);
    printf("         n1=%d, n2=%d\n",n1,n2);fflush(stdout);
    printf("         ml=%d, mu=%d\n",ml,mu);fflush(stdout);
    printf("         u=");MFPrintNVector(stdout,vu,e);printf("\n");fflush(stdout);
   }
#endif

  MFTPBVPGetJac(n,k,d,u,u,Phi,&A,&B,&C,&D,e);
  i=-1;
  ierr=0;
  jobvl='N';
  jobvr='A';
  CALLDGESVD(&jobvl,&jobvr,&n2,&n2,NULL,&n2,NULL,NULL,&n2,NULL,&n2,&t,&i,&ierr);
  lwork=round(t);

  X=(double*)realloc((void*)X,ldX*ldX*sizeof(double));

#ifndef MFNOSAFETYNET
  if(X==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",ldX*ldX*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0;
   }
#endif

  vr=(double*)realloc((void*)vr,ldX*ldX*sizeof(double));

#ifndef MFNOSAFETYNET
  if(vr==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",ldX*ldX*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0;
   }
#endif

  s=(double*)realloc((void*)s,ldX*sizeof(double));

#ifndef MFNOSAFETYNET
  if(s==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",ldX*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0;
   }
#endif

  work=(double*)realloc((void*)work,lwork*sizeof(double));

#ifndef MFNOSAFETYNET
  if(work==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",lwork*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0;
   }
#endif

  pvt=(int*)realloc((void*)pvt,n*sizeof(int));

#ifndef MFNOSAFETYNET
  if(pvt==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0;
   }
#endif

  CALLDBOSVD(A,&ldA,&n1,&ml,&mu,B,&n1,&n2,C,&n2,D,&n2,X,s,vr,work,&lwork,pvt,&ierr);

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("Smallest Singular Value %le, p=(%lf,%lf)\n",s[ldX-1],u[nx*nu],u[nx*nu+1]);
    fflush(stdout);
   }
#endif

  nullDim=0;
  for(i=0;i<ldX;i++)
    if(fabs(s[i])<1.e-3)nullDim++;

#ifdef MFALLOWVERBOSE
  if(0&&verbose){printf("%s, dimension of null space (not including the tangent) is %d\n",RoutineName,nullDim);fflush(stdout);}
#endif

  i=0;
  if(p!=NULL)CALLDBONV(&i,A,&ldA,&n1,&ml,&mu,B,&n1,&n2,C,&n2,D,&n2,X,s,vr,p);


  return nullDim;
#else
    sprintf(MFTPBVPMFErrorHandlerMsg,"The TPBVP manifold requires dgesvd from LAPACK.");
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    return 0;
#endif
 }

void MFPrintBorderedBandedMatrixByBlock(FILE *fid, int nx,int nu,int np, int nbc, int nic,int k,double *A,double *B,double *C,double *D, MFErrorHandler e)
 {
  static char RoutineName[]={"MFPrintBorderedBandedMatrixByBlock"};
  int i,j,J;
  int ldA=0;
  int nBands,nBandsL,nBandsU;
  int n1,n2;
  double a;

  nBandsL=nu-1;
  nBandsU=2*nu-1;
  ldA=2*nBandsL+nBandsU+1;
  nBands=nBandsL+nBandsU+1;
  n1=(nx-1)*nu;
  n2=nbc+nic+k;

  fprintf(fid,"Jacobian (by block)\n");
  for(i=0;i<nx-1;i++)
   {
    for(j=0;j<nu;j++)
     {
      if(j==0)fprintf(fid,"%3d [",i);
       else   fprintf(fid,"    [");
      for(J=0;J<nu;J++)
       {
        if(J>0)fprintf(fid," ");
        a=MFGetBandedMatrixElement(j+nu*i,nu*i+J,A,n1,nBandsL,nBandsU,e);
        fprintf(fid,"%10.3le",a);
       }
      fprintf(fid," | ");
      if(i<nx-2)
       {
        for(J=0;J<nu;J++)
         {
          if(J>0)fprintf(fid," ");
          a=MFGetBandedMatrixElement(j+nu*i,nu*i+J+nu,A,n1,nBandsL,nBandsU,e);
          fprintf(fid,"%10.3le",a);
         }
       }else{
        for(J=0;J<nu;J++)
         {
          if(J>0)fprintf(fid," ");
          fprintf(fid,"%10.3le",B[j+i*nu+n1*J]);
         }
       }
      fprintf(fid," | ");
      for(J=nu;J<nu+np;J++)
       {
        if(J>0)fprintf(fid," ");
        fprintf(fid,"%10.3le",B[j+i*nu+n1*J]);
       }
      fprintf(fid,"]\n");
     }
    fprintf(fid,"    ----\n");
   }

  fprintf(fid,"    ==============\n");

  for(i=0;i<nbc;i++)
   {
    fprintf(fid,"%3d [",i+n1);
    for(j=0;j<nu;j++)
     {
      if(j>0)fprintf(fid," ");
      fprintf(fid,"%10.3le",C[i+n2*j]);
     }
    fprintf(fid," | ");
    for(j=nu;j<2*nu;j++)
     {
      if(j>0)fprintf(fid," ");
      fprintf(fid,"%10.3le",C[i+n2*j]);
     }
    fprintf(fid," | ");
    for(j=(nx-2)*nu;j<(nx-1)*nu;j++)
     {
      if(j>0)fprintf(fid," ");
      fprintf(fid,"%10.3le",C[i+n2*j]);
     }
    fprintf(fid," | ");
    for(j=0;j<nu;j++)
     {
      if(j>0)fprintf(fid," ");
      fprintf(fid,"%10.3le",D[i+n2*j]);
     }
    fprintf(fid," | ");
    for(j=nu;j<nu+np;j++)
     {
      if(j>0)fprintf(fid," ");
      fprintf(fid,"%10.3le",D[i+n2*j]);
     }
    fprintf(fid,"]\n");
   }
  fflush(fid);

  return;
 }

int MFTPBVPTransformX(int nx, int nu, int np, double *X, double **Y, MFErrorHandler e)
 {
  static char RoutineName[]={"MFTPBVPTransformX"};
  double *A=NULL;
  double *B=NULL;
  double *C=NULL;
  double *D=NULL;
  int *Pivots=NULL;
  int nrhs;
  int ierr;
  int i,j,l;
  char trans='N';

#ifdef HAVE_LAPACK

  A=(double*)realloc((void*)A,nu* nu*sizeof(double));

#ifndef MFNOSAFETYNET
  if(A==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",nu*nu*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0;
   }
#endif

  B=(double*)realloc((void*)B,nu*(nu+np)*sizeof(double));

#ifndef MFNOSAFETYNET
  if(B==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",nu*(nu+np)*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0;
   }
#endif

  C=(double*)realloc((void*)C,(nu+np)*nu*sizeof(double));

#ifndef MFNOSAFETYNET
  if(C==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",(nu+np)*nu*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0;
   }
#endif

  D=(double*)realloc((void*)D,(nu+np)*(nu+np)*sizeof(double));

#ifndef MFNOSAFETYNET
  if(D==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",(nu+np)*(nu+np)*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0;
   }
#endif

  *Y =(double*)realloc((void*)(*Y),(nu+np)*(nu+np)*sizeof(double));

#ifndef MFNOSAFETYNET
  if(*Y==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",(nu+np)*(nu+np)*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0;
   }
#endif

  Pivots=(int*)realloc((void*)Pivots,nu*sizeof(int));

#ifndef MFNOSAFETYNET
  if(Pivots==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",nu*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0;
   }
#endif

  for(i=0;i<nu;i++)
   {
    for(j=0;j<nu;j++)A[i+nu*j]=X[i+(2*nu+np)*j];
    for(j=0;j<nu+np;j++)B[i+nu*j]=X[i+(2*nu+np)*(j+nu)];
   }

  for(i=0;i<nu+np;i++)
   {
    for(j=0;j<nu;j++)C[i+(nu+np)*j]=X[i+nu+(2*nu+np)*j];
    for(j=0;j<nu+np;j++)D[i+(nu+np)*j]=X[i+nu+(2*nu+np)*(j+nu)];
   }

  CALLDGETRF(&nu,&nu,A,&nu,Pivots,&ierr);

  nrhs=nu+np;
  CALLDGETRS(&trans,&nu,&nrhs,A,&nu,Pivots,B,&nu,&ierr);

  for(i=0;i<nu+np;i++)
   {
    for(j=0;j<nu+np;j++)
     {
      (*Y)[i+(nu+np)*j]=D[i+(nu+np)*j];
      for(l=0;l<nu;l++)
        (*Y)[i+(nu+np)*j]-=C[i+(nu+np)*l]*B[l+nu*j];
     }
   }

  return 0;
#else
    sprintf(MFTPBVPMFErrorHandlerMsg,"The TPBVP manifold requires dgetrf and dgetrs from LAPACK.");
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    return 0;
#endif
 }

double det2(double a00,double a10,
            double a01,double a11, MFErrorHandler e)
 {
  double det;
 
  det=a00*a11
     -a01*a10;

  return det;
 }

double det3(double a00,double a10,double a20,
            double a01,double a11,double a21,
            double a02,double a12,double a22, MFErrorHandler e)
 {
  double det;
 
  det=a00*det2(a11,a12,a21,a22,e)
     -a01*det2(a01,a02,a21,a22,e)
     +a02*det2(a01,a02,a11,a12,e);

  return det;
 }

double det4(double a00,double a10,double a20,double a30,
            double a01,double a11,double a21,double a31,
            double a02,double a12,double a22,double a32,
            double a03,double a13,double a23,double a33, MFErrorHandler e)
 {
  double det;
/*printf(" det4:\n");
  printf(" [ %10.7lf %10.7lf %10.7lf %10.7lf ]\n",a00,a01,a02,a03);
  printf(" [ %10.7lf %10.7lf %10.7lf %10.7lf ]\n",a10,a11,a12,a13);
  printf(" [ %10.7lf %10.7lf %10.7lf %10.7lf ]\n",a20,a21,a22,a23);
  printf(" [ %10.7lf %10.7lf %10.7lf %10.7lf ]\n",a30,a31,a32,a33);*/
 
  det=a00*det3(a11,a12,a13,a21,a22,a23,a31,a32,a33,e)
     -a10*det3(a01,a02,a03,a21,a22,a23,a31,a32,a33,e)
     +a20*det3(a01,a02,a03,a11,a12,a13,a31,a32,a33,e)
     -a30*det3(a01,a02,a03,a11,a12,a13,a21,a22,a23,e);

  return det;
 }

int MFAtlasDetSign(int n,double *A, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAtlasDetSign"};
/* Since LAPACK only does partial pivoting, the determinant isn't good. This routine does guassian elimination
   to reduce the matrix */
  int i,j,l;
  int ii,jj;
  int *ipiv,*jpiv;
  double piv;
  double m;
  int det;
  int good;
  int verbose=0;
  double d;
  double t;

  static double *a =NULL;
  static double *vr=NULL;
  static double *vl=NULL;
  static double *s =NULL;
  double dr,dl;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s\n",RoutineName);fflush(stdout);}
#endif

  if(n==1 && A[0]>0)return 1;
  if(n==1 && A[0]<0)return -1;

  a =(double*)realloc((void*)a ,n*n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(a==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*n*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return -1;
   }
#endif

  vr=(double*)realloc((void*)vr,n*n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(vr==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*n*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return -1;
   }
#endif

  vl=(double*)realloc((void*)vl,n*n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(vl==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*n*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return -1;
   }
#endif

  s =(double*)realloc((void*)s ,n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(s==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return -1;
   }
#endif

#ifdef MFALLOWVERBOSE
  if(verbose){printf("Copy a\n");fflush(stdout);}
#endif

  for(i=0;i<n*n;i++)a[i]=A[i];

#ifdef MFALLOWVERBOSE
  if(verbose){printf("Call svd\n");fflush(stdout);}
#endif

  MFSVD(n,a,vl,s,vr,e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("Call EVProd(Left)\n");fflush(stdout);}
#endif

  dl=MFEVProd(n,0,vl,e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("Call EVProd(Right)\n");fflush(stdout);}
#endif

  dr=MFEVProd(n,0,vr,e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf(" dl=%lf, dr=%lf\n",dl,dr);}
#endif

  d=dl*dr;
  det=1;if(d<0)det=-1;
  if(s[n-2]<1.e-5)det=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("done %s, s[n-2]=%le\n",RoutineName,s[n-2]);fflush(stdout);}
#endif

  return det;
 }

double MFEVProd(int n,int k,double *a, MFErrorHandler e)
 {
  static char RoutineName[]={"MFEVProd"};
  int i,N;
  double rR,rI;
  double R,I;
  int info;
  char jobvl;
  char jobvr;
  int lda;
  static double *wr=NULL;
  static double *wi=NULL;
  int ldvl;
  static double *vl=NULL;
  int ldvr;
  static double *vr=NULL;
  static double *work=NULL;
  int lwork;
  double *u, *Phi;
  int verbose=0;

#ifdef HAVE_LAPACK

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s, n=%d, k=%d\n",RoutineName,n,k);fflush(stdout);}
#endif

  N=n;

  lda=n;

  wr=(double*)realloc((void*)wr,N*sizeof(double));
  if(wr==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",N*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    return -1.;
   }

  wi=(double*)realloc((void*)wi,N*sizeof(double));
  if(wi==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",N*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    return -1.;
   }

  ldvl=N;
  vl=NULL;

  ldvr=N;
  vr=NULL;

  lwork=3*N;
  work=(double*)realloc((void*)work,lwork*sizeof(double));

#ifndef MFNOSAFETYNET
  if(work==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",lwork*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return -1.;
   }
#endif

  info=0;
  jobvl='N';
  jobvr='N';
  CALLDGEEV(&jobvl,&jobvr,&N,a,&lda,wr,wi,vl,&ldvl,vr,&ldvr,work,&lwork,&info);

  rR=1.;
  rI=0.;
  for(i=0;i<N;i++)
   {
    R=rR*wr[i]-rI*wi[i];
    I=rI*wr[i]+rR*wi[i];
    rR=R;
    rI=I;
   }

  return rR;
#else
    sprintf(MFTPBVPMFErrorHandlerMsg,"The TPBVP manifold requires dgeev from LAPACK.");
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    return 0;
#endif
 }

void MFSVD(int n, double *A, double *vL, double *s, double *vR, MFErrorHandler e)
 {
  static char RoutineName[]={"MFSVD"};
  static double *work=NULL;
  int lwork;
  int ierr;
  char jobvl;
  char jobvr;
  double t;
  int nx,nu,np,nbc,nic;
  int i;

#ifdef HAVE_LAPACK

  i=-1;
  ierr=0;
  jobvl='A';
  jobvr='A';
  CALLDGESVD(&jobvl,&jobvr,&n,&n,NULL,&n,NULL,NULL,&n,NULL,&n,&t,&i,&ierr);
  lwork=round(t);
  work=(double*)realloc((void*)work,lwork*sizeof(double));

#ifndef MFNOSAFETYNET
  if(work==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",lwork*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  CALLDGESVD(&jobvl,&jobvr,&n,&n,A,&n,s,vL,&n,vR,&n,work,&lwork,&ierr);

  return;
#else
  sprintf(MFTPBVPMFErrorHandlerMsg,"The TPBVP manifold requires dgesvd from LAPACK.");
  MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
  return;
#endif
 }

void MFEV(int n, double *A, double *vL, double *sr, double *si, double *vR, MFErrorHandler e)
 {
  static char RoutineName[]={"MFSVD"};
  static double *work=NULL;
  int lwork;
  int ierr;
  char jobvl;
  char jobvr;
  double t;
  int nx,nu,np,nbc,nic;
  int i;

#ifdef HAVE_LAPACK

  i=-1;
  ierr=0;
  jobvl='V';
  jobvr='V';
  lwork=4*n;
  work=(double*)realloc((void*)work,lwork*sizeof(double));

#ifndef MFNOSAFETYNET
  if(work==NULL)
   {
    sprintf(MFTPBVPMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",lwork*sizeof(double));
    MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  CALLDGEEV(&jobvl,&jobvr,&n,A,&n,sr,si,vL,&n,vR,&n,work,&lwork,&ierr);

  return;
#else
  sprintf(MFTPBVPMFErrorHandlerMsg,"The TPBVP manifold requires dgesvd from LAPACK.");
  MFSetError(e,12,RoutineName,MFTPBVPMFErrorHandlerMsg,__LINE__,__FILE__);
  return;
#endif
 }

/*! @} */

/*! @} */

#ifdef __cplusplus
}
#endif
