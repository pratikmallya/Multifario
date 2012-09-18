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

static char *id="@(#) $Id: MFSwallow.c,v 1.4 2011/07/21 17:42:46 mhender Exp $";

static char MFSwallowMFErrorHandlerMsg[256]="";

extern double MFEpsilon;

#define MFTWOPI 6.2831853071795862320

#include <MFImplicitMF.h>
#include <MFErrorHandler.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <MFFortran.h>

#ifdef __cplusplus
 extern "C" {
#endif

static int MFProjectSwallow(int,int,MFNVector,MFNKMatrix,MFNVector,void*,int*,MFErrorHandler);
static int MFTangentSwallow(int,int,MFNVector,MFNKMatrix,void*,MFErrorHandler);
static double MFScaleSwallow(int,int,MFNVector,MFNKMatrix,void*,MFErrorHandler);
static MFImplicitMF MFReadSwallow(FILE*,MFErrorHandler);

MFNVector MFNVectorFactory(MFImplicitMF,MFErrorHandler);
MFNKMatrix MFNKMatrixFactory(MFImplicitMF,MFErrorHandler);

MFImplicitMF MFIMFCreateSwallow(MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFCreateSwallow"};
  MFImplicitMF swallow;
  MFNSpace space;

  swallow=MFIMFCreateBaseClass(4,3,"Sallow",e);

  space=MFCreateNSpace(4,e);
  MFIMFSetSpace(swallow,space,e);
  MFFreeNSpace(space,e);

  MFIMFSetProject(swallow,MFProjectSwallow,e);
  MFIMFSetTangent(swallow,MFTangentSwallow,e);
  MFIMFSetScale(swallow,MFScaleSwallow,e);

  MFIMFSetVectorFactory(swallow,MFNVectorFactory,e);
  MFIMFSetMatrixFactory(swallow,MFNKMatrixFactory,e);

  return swallow;
 }

int MFProjectSwallow(int n,int k,MFNVector vu0,MFNKMatrix mPhi,MFNVector vu,void *d,int *index, MFErrorHandler e)
 {
  static char RoutineName[]={"MFProjectSwallow"};
  double err;
  double A[16]={0.,0.,0.,0.,
                0.,0.,0.,0.,
                0.,0.,0.,0.,
                0.,0.,0.,0.};
  double b[4]={0.,0.,0.,0.};
  int ipvt[4]={0,0,0,0};
  int four=4;
  int one=1;
  int info=0;
  int itimes=0;
  char trans='N';
  double *u0,*u,*Phi;
  int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s\n",RoutineName);fflush(stdout);}
#endif

  u0=MFNV_CStar(vu0,e);
  Phi=MFNKM_CStar(mPhi,e);
  u=MFNV_CStar(vu,e);

/* u^4 + t3 u^2 + t2 u + t1 = 0 */

/* du : 4 u^3 + 2 t3 u + t2  */

  u[0]=u0[0];
  u[1]=u0[1];
  u[2]=u0[2];
  u[3]=u0[3];

  err=1.;
  itimes=0;
  while(err>1.e-10)
   {
    if(itimes>40)
     {

#ifdef MFALLOWVERBOSE
      if(verbose){printf("%s failed\n",RoutineName);fflush(stdout);}
#endif

      return 0;
     }
    A[0+4*0]=u[2]+u[0]*(2*u[3]+4*u[0]*u[0]);
    A[0+4*1]=1.;
    A[0+4*2]=u[0];
    A[0+4*3]=u[0]*u[0];
    b[0]=-(u[1]+u[0]*(u[2]+u[0]*(u[3]+u[0]*u[0])));
    A[1+4*0]=Phi[0];
    A[1+4*1]=Phi[1];
    A[1+4*2]=Phi[2];
    A[1+4*3]=Phi[3];
    b[1]=-Phi[0]*(u[0]-u0[0])-Phi[1]*(u[1]-u0[1])
         -Phi[2]*(u[2]-u0[2])-Phi[3]*(u[3]-u0[3]);
    A[2+4*0]=Phi[4];
    A[2+4*1]=Phi[5];
    A[2+4*2]=Phi[6];
    A[2+4*3]=Phi[7];
    b[2]=-Phi[4]*(u[0]-u0[0])-Phi[5]*(u[1]-u0[1])
         -Phi[6]*(u[2]-u0[2])-Phi[7]*(u[3]-u0[3]);
    A[3+4*0]=Phi[8];
    A[3+4*1]=Phi[9];
    A[3+4*2]=Phi[10];
    A[3+4*3]=Phi[11];
    b[3]=-Phi[ 8]*(u[0]-u0[0])-Phi[ 9]*(u[1]-u0[1])
         -Phi[10]*(u[2]-u0[2])-Phi[11]*(u[3]-u0[3]);

    err=sqrt(b[0]*b[0]+b[1]*b[1]+b[2]*b[2]+b[3]*b[3]);

    CALLDGETRF(&four,&four,A,&four,ipvt,&info);
    CALLDGETRS(&trans,&four,&one,A,&four,ipvt,b,&four,&info);

    u[0]=u[0]+b[0];
    u[1]=u[1]+b[1];
    u[2]=u[2]+b[2];
    u[3]=u[3]+b[3];
    itimes++;
   }
  
  *index=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("done %s, %d Newton steps\n",RoutineName,itimes);fflush(stdout);}
#endif

  return 1;
 }

int MFTangentSwallow(int n,int k,MFNVector vu,MFNKMatrix mPhi,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFTangentSwallow"};
  double A[16]={0.,0.,0.,0.,
                0.,0.,0.,0.,
                0.,0.,0.,0.,
                0.,0.,0.,0.};

  char jobvl='N';                   /* No left  eigenvectors */
  char jobvr='V';                   /*    right eigenvectors */
  int N=4;                          /* Number of columns */
  int lda=4;                        /* Leading dimension of A */
  double wr[4]={0.,0.,0.,0.};       /* Array to hold real part of eigenvalues */
  double wi[4]={0.,0.,0.,0.};       /* Array to hold imag part of eigenvalues */
  int ldvl=4;                       /* Leading dimension of vl */
  double vl[16]={0.,0.,0.,0.,       /* Array to hold left eigenvectors */
                 0.,0.,0.,0.,
                 0.,0.,0.,0.,
                 0.,0.,0.,0.};
  int ldvr=4;                       /* Leading dimension of vr */
  double vr[16]={0.,0.,0.,0.,       /* Array to hold right eigenvectors */
                 0.,0.,0.,0.,
                 0.,0.,0.,0.,
                 0.,0.,0.,0.};
  double work[16]={0.,0.,0.,0.,     /* Work array */
                   0.,0.,0.,0.,
                   0.,0.,0.,0.,
                   0.,0.,0.,0.};
  int lwork=16;                     /* Length of work array */
  int info=0;
  int i,j,l,m,mm;
  double t;
  double *u,*Phi;
  int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s\n",RoutineName);fflush(stdout);}
#endif

  Phi=MFNKM_CStar(mPhi,e);
  u=MFNV_CStar(vu,e);

  A[ 0]=u[2]+u[0]*(2*u[3]+4*u[0]*u[0]);
  A[ 1]=0.;
  A[ 2]=0.;
  A[ 3]=0.;
  A[ 4]=1.;
  A[ 5]=0.;
  A[ 6]=0.;
  A[ 7]=0.;
  A[ 8]=u[0];
  A[ 9]=0.;
  A[10]=0.;
  A[11]=0.;
  A[12]=u[0]*u[0];
  A[13]=0.;
  A[14]=0.;
  A[15]=0.;

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("A:\n");
    for(i=0;i<n;i++)
     {
      printf(" [%lf",A[i+n*0]);for(j=1;j<n;j++)printf(",%lf",A[i+n*j]);printf("]\n");fflush(stdout);
     }
    printf("\n");
   }
#endif

  info=0;
  jobvl='N';
  jobvr='V';

#ifdef MFALLOWVERBOSE
  if(verbose){printf("      call dgeev\n");fflush(stdout);}
#endif

  CALLDGEEV(&jobvl,&jobvr,&N,A,&lda,wr,wi,vl,&ldvl,vr,&ldvr,work,&lwork,&info);

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("      done dgeev\n");fflush(stdout);
    printf("\n      Null Vectors\n");fflush(stdout);
   }
#endif

  m=0;
  for(i=0;i<n;i++)
   {
    if(fabs(wr[i])+fabs(wi[i])<1.e-7)
     {
#ifdef MFALLOWVERBOSE
      if(verbose){printf("     Phi[%d]=v[%d]=(%lf",m,i,vr[0+n*i]);for(j=1;j<n;j++)printf(",%lf",vr[j+n*i]);printf(")\n");fflush(stdout);}
#endif
      for(j=0;j<n;j++)vr[j+n*m]=vr[j+n*i];
      m++;
     }
   }

#ifdef MFALLOWVERBOSE
  if(verbose){printf("\n      Gram-Schmidt\n");fflush(stdout);}
#endif

  mm=0;
  t=0.;
  for(j=0;j<n;j++)t+=vr[j]*vr[j];
  if(fabs(t)>1.e-15)
   {
    t=1./sqrt(t);
    for(j=0;j<n;j++){vr[j]=vr[j]*t;Phi[j+n*mm]=vr[j];}

#ifdef MFALLOWVERBOSE
    if(verbose){printf("\n     Phi[%d]=(%lf",mm,Phi[0+n*mm]);for(j=1;j<n;j++)printf(",%lf",Phi[j+n*mm]);printf(")\n");fflush(stdout);}
#endif

    mm++;
   }

  for(i=1;i<m;i++)
   {
    for(l=0;l<i;l++)
     {
      t=0.;
      for(j=0;j<n;j++)t+=vr[j+n*l]*vr[j+n*i];
      for(j=0;j<n;j++)vr[j+n*i]=vr[j+n*i]-t*vr[j+n*l];
     }

    t=0.;
    for(j=0;j<n;j++)t+=vr[j+n*i]*vr[j+n*i];
    if(fabs(t)>1.e-15)
     {
      t=1./sqrt(t);
      for(j=0;j<n;j++){vr[j+n*i]=vr[j+n*i]*t;Phi[j+n*mm]=vr[j+n*i];}

#ifdef MFALLOWVERBOSE
      if(verbose){printf("     Phi[%d]=(%lf",mm,Phi[0+n*mm]);for(j=1;j<n;j++)printf(",%lf",Phi[j+n*mm]);printf(")\n");fflush(stdout);}
#endif

      mm++;
     }
   }

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("\n      test Gram-Schmidt\n");
    printf("\n");
    for(i=0;i<k;i++)
     {
      printf("[");
      for(j=0;j<k;j++)
       {
        t=0.;
        for(m=0;m<n;m++)t+=Phi[m+n*i]*Phi[m+n*j];
        if(j>0)printf(" ");
        printf("%lf",t);
       }
      printf("]\n");
     }
    fflush(stdout);
    printf("done %s\n",RoutineName);fflush(stdout);
   }
#endif

  return 1;
 }

double MFScaleSwallow(int n,int k,MFNVector vu,MFNKMatrix mPhi,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFScaleSwallow"};
  double r=1.;
  double A[16]={0.,0.,0.,0.,
                0.,0.,0.,0.,
                0.,0.,0.,0.,
                0.,0.,0.,0.};
  double b[4]={0.,0.,0.,0.};
  double c[9]={0.,0.,0.,
               0.,0.,0.,
               0.,0.,0.};
  double cn;
  double tu=0.;
  double t1=0.;
  double t2=0.;
  double t3=0.;
  double su=0.;
  double s1=0.;
  double s2=0.;
  double s3=0.;
  int i,j;
  int ipvt[4]={0,0,0,0};
  int four=4;
  int one=1;
  int info=0;
  char trans='N';
  char jobvl='N';                   /* No left  eigenvectors */
  char jobvr='N';                   /*    right eigenvectors */
  int N=3;                          /* Number of columns */
  int ldc=3;                        /* Leading dimension of c */
  double wr[4]={0.,0.,0.,0.};       /* Array to hold real part of eigenvalues */
  double wi[4]={0.,0.,0.,0.};       /* Array to hold imag part of eigenvalues */
  int ldvl=2;                       /* Leading dimension of vl */
  double vl[9]={0.,0.,0.,           /* Array to hold left eigenvectors */
                0.,0.,0.,   
                0.,0.,0.};
  int ldvr=3;                       /* Leading dimension of vr */
  double vr[9]={0.,0.,0.,           /* Array to hold right eigenvectors */
                0.,0.,0.,   
                0.,0.,0.};
  double work[9]={0.,0.,0.,         /* Work array */
                  0.,0.,0.,
                  0.,0.,0.};
  int lwork=9;                     /* Length of work array */
  double *u,*Phi;

  Phi=MFNKM_CStar(mPhi,e);
  u=MFNV_CStar(vu,e);

/* u^4 + t3 u^2 + t2 u + t1 = 0 */

/* du : 4 u^3 + 2 t3 u + t2  */

  A[0+4*0]=u[2]+u[0]*(2*u[3]+4*u[0]*u[0]);
  A[0+4*1]=1.;
  A[0+4*2]=u[0];
  A[0+4*3]=u[0]*u[0];
  A[1+4*0]=Phi[0];
  A[1+4*1]=Phi[1];
  A[1+4*2]=Phi[2];
  A[1+4*3]=Phi[3];
  A[2+4*0]=Phi[4];
  A[2+4*1]=Phi[5];
  A[2+4*2]=Phi[6];
  A[2+4*3]=Phi[7];
  A[3+4*0]=Phi[8];
  A[3+4*1]=Phi[9];
  A[3+4*2]=Phi[10];
  A[3+4*3]=Phi[11];
  CALLDGETRF(&four,&four,A,&four,ipvt,&info);

  for(i=0;i<3;i++)
   {
    tu=Phi[0+4*i];
    t1=Phi[1+4*i];
    t2=Phi[2+4*i];
    t3=Phi[3+4*i];
    for(j=i;j<3;j++)
     {
      su=Phi[0+4*j];
      s1=Phi[1+4*j];
      s2=Phi[2+4*j];
      s3=Phi[3+4*j];

      b[0]=-(12*u[0]*u[0]+2*u[3])*su*tu-2*u[0]*su*t3-2*u[0]*tu*s3-su*t2-tu*s2;
      b[1]=0.;
      b[2]=0.;
      b[3]=0.;

      CALLDGETRS(&trans,&four,&one,A,&four,ipvt,b,&four,&info);

      c[i+3*j]=sqrt(b[0]*b[0]+b[1]*b[1]+b[2]*b[2]+b[3]*b[3]);
      c[j+3*i]=sqrt(b[0]*b[0]+b[1]*b[1]+b[2]*b[2]+b[3]*b[3]);
     }
   }

  info=0;
  jobvl='N';
  jobvr='N';
  CALLDGEEV(&jobvl,&jobvr,&N,c,&ldc,wr,wi,vl,&ldvl,vr,&ldvr,work,&lwork,&info);

  cn=sqrt(wr[0]*wr[0]+wi[0]*wi[0]);
  tu=sqrt(wr[1]*wr[1]+wi[1]*wi[1]);
  if(tu>cn)cn=tu;
  tu=sqrt(wr[2]*wr[2]+wi[2]*wi[2]);
  if(tu>cn)cn=tu;
  tu=sqrt(wr[3]*wr[3]+wi[3]*wi[3]);
  if(tu>cn)cn=tu;


  r=sqrt(2*MFEpsilon/cn);
  printf("ScaleSwallow, cn=%lf, %lf, %lf, %lf, max=%lf, r=%lf\n",sqrt(wr[0]*wr[0]+wi[0]*wi[0]),sqrt(wr[1]*wr[1]+wi[1]*wi[1]),sqrt(wr[2]*wr[2]+wi[2]*wi[2]),sqrt(wr[3]*wr[3]+wi[3]*wi[3]),cn,r);fflush(stdout);
  if(r<.001)r=.001;
  return r;
 }

MFImplicitMF MFReadSwallow(FILE *fid, MFErrorHandler e)
 {
  static char RoutineName[]={"MFReadSwallow"};
  MFImplicitMF swallow;

  swallow=MFIMFCreateSwallow(e);

  return swallow;
 }

#ifdef __cplusplus
}
#endif
