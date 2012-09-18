/*
 *  PROGRAM NAME:  multifario
 *
 *  (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
 *  CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
 *
 *  Please refer to the LICENSE file in the top directory
 *
 *      author: Mike Henderson mhender@watson.ibm.com
 */

static char *id="@(#) $Id: IMFFixedPt.c,v 1.11 2011/07/21 17:42:46 mhender Exp $";

static char IMFFixedPointErrorMsg[256];

#include <IMF.h>
#include <MFAtlas.h>
#include <IMFExpansion.h>
#include <IMFFixedPt.h>
#include <MFFortran.h>
#include <MFErrorHandler.h>
#include <MFPrint.h>
#include <stdio.h>
#include <math.h>

#ifdef __cplusplus
 extern "C" {
#endif

#define du(i,j) du[i+n*j]
#define dv(i,j) dv[i+n*j]
#define ddu(i,j,k) ddu[i+n*(j+m*k)]
#define dddu(i,j,k,l) dddu[i+n*(j+m*(k+m*l))]

#define up(i) up[i]
#define dup(i,j) dup[i+n*j]
#define ddup(i,j,k) ddup[i+n*(j+m*k)]

#define F(i) F[i]
#define dF(i,j) dF[i+n*j]
#define ddF(i,j,k) ddF[i+n*(j+n*k)]
#define dddF(i,j,k,l) dddF[i+n*(j+n*(k+m*l))]

#define dudFdu(i,j) dudFdu[i+m*j]
#define dvdFdu(i,j) dvdFdu[i+(n-m)*j]
#define dvdFdv(i,j) dvdFdv[i+(n-m)*j]
#define duddFdudu(i,j,k) duddFdudu[i+m*(j+m*k)]
#define dvddFdudu(i,j,k) dvddFdudu[i+(n-m)*(j+m*k)]
#define dudFddu(i,j,k) dudFddu[i+m*(j+m*k)]
#define ddudFdu(i,j,k) ddudFdu[i+m*(j+m*k)]
#define duddu(i,j,k) duddu[i+m*(j+m*k)]
#define beta(i,j,k) beta[i-1+(m-1)*(j-1+(m-1)*k)]

#define da(i,j)         da[i+m*j]
#define dda(i,j,k)     dda[i+m*(j+m*k)]
#define ddda(i,j,k,l) ddda[i+m*(j+m*(k+m*l))]

#define A0(i,j,k,p,q,r) A[i+(n-m)*(j+m*k)+(n-m)*m*m*( p+(n-m)*(q+m*r))]
#define B0(i,j,k)       b[i+(n-m)*(j+m*k)]
#define eta(i,j,k)       b[i+(n-m)*(j+m*k)]

#define A1(i,j,k,l,p,q,r,s) AA[i+(n-m)*(j+m*(k+m*l))+(n-m)*m*m*m*(p+(n-m)*(q+m*(r+m*s)))]
#define B1(i,j,k,l)         bb[i+(n-m)*(j+m*(k+m*l))]
#define round(x) ((x)>=0?(long)((x)+0.5):(long)((x)-0.5))

#ifndef MAX
#define MAX(X, Y) ((X) < (Y) ? (Y) : (X))
#endif

/*! \fn void IMFFindExpansionNearFixedPt(MFNVector U0, MFKVector P0, MFNKMatrix Du, MFNKMatrix Dv, IMFFlow F, IMFExpansion U, IMFExpansion a, MFErrorHandler e);
 *  \brief Given a fixed point in a flow and a decomposition of the phase space into two invariant linear subspaces, creates
 *            an expansion of the invariant manifold corresponding to one of the invariant subspaces.
 *
 *  \param U0 The fixed point.
 *  \param P0 The parameters of the flow of the fixed point.
 *  \param Du The first invariant linear subspace.
 *  \param Dv The second, complementary invariant linear subspace.
 *  \param F  The flow.
 *  \param U  A user provided expansion in which the expansion defining the surface is placed.
 *  \param a  A user provided expansion in which the expansion defining the complement of the surface is placed.
 *  \param e  An error handler.
 */
void IMFFindExpansionNearFixedPt(MFNVector U0, MFKVector P0, MFNKMatrix Du, MFNKMatrix Dv, IMFFlow F, IMFExpansion U, IMFExpansion a, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFFindExpansionNearFixedPt"};
/* on input:  du is a set of k orthonormal vectors in n-space, spanning the subspace */
/*            dF[i+n*p] is the Jacobian */
/*            ddF[i+n*(p+n*q)] is the second derivative */
/*            dddF[i+n*(p+n*(q+n*r))] is the third derivative */
/* on output: ddu[p+n*(i+k*j)] is the second derivative of the subspace */
/*            dddu[p+n*(i+k*(j+k*l))] is the third derivative of the subspace */
  double *A;
  double *b;
  double *AA;
  double *bb;
  int i,j,k,l,p,q,r,w;
  double dot,sum,dot1;
  double error;
  int io,tt;
  int jp,kp,lp;
  double *dudFdu;
  double *dvdFdu;
  double *dvdFdv;
  double *dvddFdudu;
  double *duddFdudu;
  double *dudFddu;
  int verboser=1;
  int verbose=0;

  int n;
  int m;
  double *u0;
  double *a0;
  double *du;
  double *dv;
  double *dF;
  double *ddF;
  double *dddF;
  double *ddu;
  double *dddu;
  double *da;
  double *dda;
  double *ddda;

  n=MFNKMatrixN(Du,e);
  m=MFNKMatrixK(Du,e);
  printf("%s, n=%d, m=%d\n",RoutineName,n,m);

  du=MFNKM_CStar(Du,e);
  dv=MFNKM_CStar(Dv,e);

  u0=MFNV_CStar(U0,e);
  a0=(double*)malloc(m*sizeof(double));

#ifndef MFNOSAFETYNET
  if(a0==NULL)
   {
    sprintf(IMFFixedPointErrorMsg,"Out of memory, trying to allocate %d bytes",m*sizeof(double));
    MFSetError(e,12,RoutineName,IMFFixedPointErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  for(i=0;i<m;i++)a0[i]=0.;

  da=(double*)malloc(m*m*sizeof(double));

#ifndef MFNOSAFETYNET
  if(da==NULL)
   {
    sprintf(IMFFixedPointErrorMsg,"Out of memory, trying to allocate %d bytes",m*m*sizeof(double));
    MFSetError(e,12,RoutineName,IMFFixedPointErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  dda=(double*)malloc(m*m*m*sizeof(double));

#ifndef MFNOSAFETYNET
  if(dda==NULL)
   {
    sprintf(IMFFixedPointErrorMsg,"Out of memory, trying to allocate %d bytes",m*m*m*sizeof(double));
    MFSetError(e,12,RoutineName,IMFFixedPointErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  ddda=(double*)malloc(m*m*m*m*sizeof(double));

#ifndef MFNOSAFETYNET
  if(ddda==NULL)
   {
    sprintf(IMFFixedPointErrorMsg,"Out of memory, trying to allocate %d bytes",m*m*m*m*sizeof(double));
    MFSetError(e,12,RoutineName,IMFFixedPointErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif


  ddu=(double*)malloc(n*m*m*sizeof(double));

#ifndef MFNOSAFETYNET
  if(ddu==NULL)
   {
    sprintf(IMFFixedPointErrorMsg,"Out of memory, trying to allocate %d bytes",n*m*m*sizeof(double));
    MFSetError(e,12,RoutineName,IMFFixedPointErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  dddu=(double*)malloc(n*m*m*m*sizeof(double));

#ifndef MFNOSAFETYNET
  if(dddu==NULL)
   {
    sprintf(IMFFixedPointErrorMsg,"Out of memory, trying to allocate %d bytes",n*m*m*m*sizeof(double));
    MFSetError(e,12,RoutineName,IMFFixedPointErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  dF=(double*)malloc(n*n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(dF==NULL)
   {
    sprintf(IMFFixedPointErrorMsg,"Out of memory, trying to allocate %d bytes",n*n*sizeof(double));
    MFSetError(e,12,RoutineName,IMFFixedPointErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  ddF=(double*)malloc(n*n*n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(ddF==NULL)
   {
    sprintf(IMFFixedPointErrorMsg,"Out of memory, trying to allocate %d bytes",n*n*n*sizeof(double));
    MFSetError(e,12,RoutineName,IMFFixedPointErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  dddF=(double*)malloc(n*n*n*n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(dddF==NULL)
   {
    sprintf(IMFFixedPointErrorMsg,"Out of memory, trying to allocate %d bytes",n*n*n*n*sizeof(double));
    MFSetError(e,12,RoutineName,IMFFixedPointErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  IMFEvaluateDerivativeOfFlow(F,U0,P0,dF,e);
  IMFEvaluateSecondDerivativeOfFlow(F,U0,P0,ddF,e);
  IMFEvaluateThirdDerivativeOfFlow(F,U0,P0,dddF,e);

/* check MFNKMatrixN(Du,e)==n */
/* check MFNKMatrixK(Dv,e)==n-m */

  printf("Expansion near Fixed Point ");
  MFPrintNVector(stdout,U0,e);
  printf("\n\n");fflush(stdout);

#ifdef MFALLOWVERBOSE
  if(verboser)
   {
    for(i=0;i<n;i++)
     {
      if(i==0)printf("F_u = ");
       else printf("      ");
      printf("(%le",dF[i]);fflush(stdout);
      for(j=1;j<n;j++)printf(",%le",dF[i+n*j]);
      printf(")\n");fflush(stdout);
     }
    printf("\n\n");fflush(stdout);

    for(j=0;j<m;j++)
     {
      printf("u_{,%d} = (%le",j,du(0,j));fflush(stdout);
      for(i=1;i<n;i++)printf(",%le",du(i,j));fflush(stdout);
      printf(")\n");fflush(stdout);
     }
    printf("\n");fflush(stdout);

    for(j=0;j<n-m;j++)
     {
      printf("v_{,%d} = (%le",j,dv(0,j));fflush(stdout);
      for(i=1;i<n;i++)printf(",%le",dv(i,j));fflush(stdout);
      printf(")\n");fflush(stdout);
     }
    printf("\n");fflush(stdout);
   }
#endif

/* Is F^i_,p,q symmetric in p and q? */

  for(i=0;i<n;i++)
   {
    for(p=0;p<n;p++)
     {
      for(q=p+1;q<n;q++)
       {
        if(fabs(ddF(i,p,q)-ddF(i,q,p))>1.e-7)
         {printf("F^%d_{,%d,%d}(%le)!=F^%d_{,%d,%d}(%le)\n",i,p,q,ddF(i,p,q),
                                                            i,q,p,ddF(i,q,p));
          fflush(stdout);}
       }
     }
   }

/* Check 1: is du o.n. */

  error=0.;
  for(i=0;i<m;i++)
   {
    for(j=0;j<m;j++)
     {
      sum=0.;
      if(i==j)sum=-1.;
      for(p=0;p<n;p++)sum+=du(p,i)*du(p,j);
      error+=fabs(sum);
      if(fabs(sum)>1.e-7)printf("  (%d,%d) of order 0, u-nrm eq. error %le\n",i,j,sum);
     }
   }
  if(error>1.e-5)printf("Total error in order 0, u-nrm eqs. %le\n",error);
   else          printf("Passed test on order 0, u-nrm eqs.\n");

  error=0.;
  for(i=0;i<n-m;i++)
   {
    for(j=0;j<n-m;j++)
     {
      sum=0.;
      if(i==j)sum=-1.;
      for(p=0;p<n;p++)sum+=dv(p,i)*dv(p,j);
      error+=fabs(sum);
      if(fabs(sum)>1.e-7)printf("  (%d,%d) of order 0, v-nrm eq. error %le\n",i,j,sum);
     }
   }
  if(error>1.e-5)printf("Total error in order 0, v-nrm eqs. %le\n",error);
   else          printf("Passed test on order 0, v-nrm eqs.\n");

  error=0.;
  for(i=0;i<n-m;i++)
   {
    for(j=0;j<m;j++)
     {
      sum=0.;
      for(p=0;p<n;p++)sum+=dv(p,i)*du(p,j);
      error+=fabs(sum);
      if(fabs(sum)>1.e-7)printf("  (%d,%d) of order 0, uv-nrm eq. error %le\n",i,j,sum);
     }
   }
  if(error>1.e-5)printf("Total error in order 0, uv-nrm eqs. %le\n",error);
   else          printf("Passed test on order 0, uv-nrm eqs.\n");

/* Solve order 1, fn eqs for da */

  for(i=0;i<m;i++)
   {
    for(j=0;j<m;j++)
     {
      da(i,j)=0;
      for(p=0;p<n;p++)
       {
        for(q=0;q<n;q++)
         {
          da(i,j)+=du(p,i)*dF(p,q)*du(q,j);
         }
       }
     }
   }

#ifdef MFALLOWVERBOSE
  if(verboser)
   {
    for(j=0;j<m;j++)
     {
      printf("a_{,%d} = (%le",j,da(0,j));fflush(stdout);
      for(i=1;i<m;i++)printf(",%le",da(i,j));fflush(stdout);
      printf(")\n");fflush(stdout);
     }
    printf("\n");fflush(stdout);
   }
#endif


/* Check 2: is dF*du = du*da ? (i.e. if span(du) closed under dF) */

  error=0.;
  for(i=0;i<n;i++)
   {
    for(j=0;j<m;j++)
     {
      sum=0.;
      for(p=0;p<n;p++)sum+=dF(i,p)*du(p,j);
      for(p=0;p<m;p++)sum-=du(i,p)*da(p,j);
      error+=fabs(sum);
      if(fabs(sum)>1.e-7)printf("  (%d,%d) of order 1, fn  eq. error=%le\n",i,j,sum);
     }
   }
  if(error>1.e-5)printf("Total error in order 1, fn  eqs. %le\n",error);
   else          printf("Passed test on order 1, fn eqs.\n");
  if(error>1.e-5)printf("     The invariant subspace defined by du is not invariant\n");

/* ------------------------------------------------------- */
/* Need udFu, vdFu, vdFv and vddFuu */

  dudFdu=(double*)malloc(m*m*sizeof(double));

#ifndef MFNOSAFETYNET
  if(dudFdu==NULL)
   {
    sprintf(IMFFixedPointErrorMsg,"Out of memory, trying to allocate %d bytes",m*m*sizeof(double));
    MFSetError(e,12,RoutineName,IMFFixedPointErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  for(i=0;i<m;i++)
   {
    for(j=0;j<m;j++)
     {
      dudFdu(i,j)=da(i,j);

#ifdef MFALLOWVERBOSE
      if(verbose){printf("u_%d,F_u u_%d = %le\n",i,j,dudFdu(i,j));fflush(stdout);}
#endif

     }
   }

  dvdFdu=(double*)malloc((n-m)*m*sizeof(double));

#ifndef MFNOSAFETYNET
  if(dvdFdu==NULL)
   {
    sprintf(IMFFixedPointErrorMsg,"Out of memory, trying to allocate %d bytes",(n-m)*m*sizeof(double));
    MFSetError(e,12,RoutineName,IMFFixedPointErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  for(i=0;i<n-m;i++)
   {
    for(j=0;j<m;j++)
     {
      dvdFdu(i,j)=0.;
      for(q=0;q<n;q++)
       {
        for(p=0;p<n;p++)
         {
          dvdFdu(i,j)+=dv(p,i)*dF(p,q)*du(q,j);
         }
       }

#ifdef MFALLOWVERBOSE
      if(verbose){printf("v_%d,F_u u_%d = %le (must be zero for invariance)\n",i,j,dvdFdu(i,j));fflush(stdout);}
#endif

     }
   }

  dvdFdv=(double*)malloc((n-m)*(n-m)*sizeof(double));

#ifndef MFNOSAFETYNET
  if(dvdFdv==NULL)
   {
    sprintf(IMFFixedPointErrorMsg,"Out of memory, trying to allocate %d bytes",(n-m)*(n-m)*sizeof(double));
    MFSetError(e,12,RoutineName,IMFFixedPointErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  for(i=0;i<n-m;i++)
   {
    for(j=0;j<n-m;j++)
     {
      dvdFdv(i,j)=0.;
      for(q=0;q<n;q++)
       {
        for(p=0;p<n;p++)
         {
          dvdFdv(i,j)+=dv(p,i)*dF(p,q)*dv(q,j);
         }
       }

#ifdef MFALLOWVERBOSE
      if(verbose){printf("v_%d,F_u v_%d = %le\n",i,j,dvdFdv(i,j));fflush(stdout);}
#endif

     }
   }

  duddFdudu=(double*)malloc(m*m*m*sizeof(double));

#ifndef MFNOSAFETYNET
  if(duddFdudu==NULL)
   {
    sprintf(IMFFixedPointErrorMsg,"Out of memory, trying to allocate %d bytes",m*m*m*sizeof(double));
    MFSetError(e,12,RoutineName,IMFFixedPointErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  for(i=0;i<m;i++)
   {
    for(j=0;j<m;j++)
     {
      for(k=0;k<m;k++)
       {
        duddFdudu(i,j,k)=0.;
        for(p=0;p<n;p++)
         {
          for(q=0;q<n;q++)
           {
            for(r=0;r<n;r++)
             {
              duddFdudu(i,j,k)+=du(p,i)*ddF(p,q,r)*du(q,j)*du(r,k);
             }
           }
         }
       }
     }
   }

  dvddFdudu=(double*)malloc((n-m)*m*m*sizeof(double));

#ifndef MFNOSAFETYNET
  if(dvddFdudu==NULL)
   {
    sprintf(IMFFixedPointErrorMsg,"Out of memory, trying to allocate %d bytes",(n-m)*m*m*sizeof(double));
    MFSetError(e,12,RoutineName,IMFFixedPointErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  for(i=0;i<n-m;i++)
   {
    for(j=0;j<m;j++)
     {
      for(k=0;k<m;k++)
       {
        dvddFdudu(i,j,k)=0.;
        for(p=0;p<n;p++)
         {
          for(q=0;q<n;q++)
           {
            for(r=0;r<n;r++)
             {
              dvddFdudu(i,j,k)+=dv(p,i)*ddF(p,q,r)*du(q,j)*du(r,k);
             }
           }
         }
/*      printf("v_%d,F_uu u_%d u_%d = %le\n",i,j,k,dvddFdudu(i,j,k));fflush(stdout);*/
       }
     }
   }

/* ------------------------------------------------------- */

  A=(double*)malloc((n-m)*m*m*(n-m)*m*m*sizeof(double));

#ifndef MFNOSAFETYNET
  if(A==NULL)
   {
    sprintf(IMFFixedPointErrorMsg,"Out of memory, trying to allocate %d bytes",(n-m)*m*m*(n-m)*sizeof(double));
    MFSetError(e,12,RoutineName,IMFFixedPointErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  b=(double*)malloc((n-m)*m*m*sizeof(double));

#ifndef MFNOSAFETYNET
  if(b==NULL)
   {
    sprintf(IMFFixedPointErrorMsg,"Out of memory, trying to allocate %d bytes",(n-m)*m*m*sizeof(double));
    MFSetError(e,12,RoutineName,IMFFixedPointErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  for(i=0;i<(n-m)*m*m*(n-m)*m*m;i++)A[i]=0;
  for(i=0;i<(n-m)*m*m;i++)b[i]=0;

  for(i=0;i<n-m;i++)
   {
    for(j=0;j<m;j++)
     {
      for(k=0;k<m;k++)
       {
        for(q=0;q<m;q++)
         {
          A0(i,j,k,i,q,j)+=dudFdu(q,k);
          A0(i,j,k,i,q,k)+=dudFdu(q,j);
         }

        for(q=0;q<n-m;q++)
          A0(i,j,k,q,j,k)-=dvdFdv(i,q);

        B0(i,j,k)=dvddFdudu(i,j,k);
       }
     }
   }

/*printf("System for second derivatives\n");fflush(stdout);
  IMFPrintFull((n-m)*m*m,A,b,e);*/

  if(!MFSolveFull((n-m)*m*m,A,b,e))exit(12);

  printf("Second derivatives\n");fflush(stdout);
  for(p=0;p<n-m;p++)
   {
    for(j=0;j<m;j++)
     {
      if(fabs(B0(p,j,j))>1.e-7)printf("   %le s_%d s_%d v_{,%d}\n",B0(p,j,j)*.5,j,j,p);
      for(k=j+1;k<m;k++)
       if(fabs(B0(p,j,k))>1.e-7)printf("   %le s_%d s_%d v_{,%d}\n",B0(p,j,k),j,k,p);
     }
   }


  for(i=0;i<n;i++)
   {
    for(j=0;j<m;j++)
     {
      for(k=0;k<m;k++)
       {
        ddu(i,j,k)=0;
        for(p=0;p<n-m;p++)
          ddu(i,j,k)+=B0(p,j,k)*dv(i,p);
       }
/*    for(k=0;k<j;k++)
        ddu(i,j,k)=ddu(i,k,j);*/
     }
   }

  if(0){
   double t;

   for(i=0;i<m;i++)
   for(j=0;j<m;j++)
     for(k=0;k<m;k++)
       for(l=0;l<m;l++)
        {
         t=0.;for(p=0;p<n;p++)t+=ddu(p,i,j)*ddu(p,k,l);
         printf(" u^p_{%d %d} u^p_{%d %d} = %le\n",i,j,k,l,t);fflush(stdout);
        }
  }

#ifdef MFALLOWVERBOSE
  if(verboser)
   {
    for(j=0;j<m;j++)
     {
      for(k=0;k<m;k++)
       {
        printf("u_{,%d,%d} = (%le",j,k,ddu(0,j,k));fflush(stdout);
        for(i=1;i<n;i++)printf(",%le",ddu(i,j,k));fflush(stdout);
        printf(")\n");fflush(stdout);
       }
     }
    printf("\n");fflush(stdout);
    for(j=0;j<m;j++)
     {
      for(k=0;k<m;k++)
       {
        printf("F_{,%d,%d} = (%le",j,k,ddF(0,j,k));fflush(stdout);
        for(i=1;i<n;i++)printf(",%le",ddF(i,j,k));fflush(stdout);
        printf(")\n");fflush(stdout);
       }
     }
    printf("\n");fflush(stdout);
   }
#endif

/* Need (du,dF ddu) */

  dudFddu=(double*)malloc(m*m*m*sizeof(double));

#ifndef MFNOSAFETYNET
  if(dudFddu==NULL)
   {
    sprintf(IMFFixedPointErrorMsg,"Out of memory, trying to allocate %d bytes",m*m*m*sizeof(double));
    MFSetError(e,12,RoutineName,IMFFixedPointErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  for(i=0;i<m;i++)
   {
    for(j=0;j<m;j++)
     {
      for(k=0;k<m;k++)
       {
        dudFddu(i,j,k)=0.;
        for(p=0;p<n;p++)
         {
          for(q=0;q<n;q++)
           {
            dudFddu(i,j,k)+=du(p,i)*dF(p,q)*ddu(q,j,k);
           }
         }
       }
     }
   }

  for(i=0;i<m;i++)
    for(j=0;j<m;j++)
      for(k=0;k<m;k++)
       {
        dda(i,j,k)=duddFdudu(i,j,k);
        for(r=0;r<n-m;r++)
          for(p=0;p<n;p++)
            for(q=0;q<n;q++)
              dda(i,j,k)+=du(p,i)*dF(p,q)*dv(q,r)*eta(r,j,k);
       }

#ifdef MFALLOWVERBOSE
  if(verboser)
   {
    for(j=0;j<m;j++)
     {
      for(k=0;k<m;k++)
       {
        printf("a_{,%d,%d} = (%le",j,k,dda(0,j,k));fflush(stdout);
        for(i=1;i<m;i++)printf(",%le",dda(i,j,k));fflush(stdout);
        printf(")\n");fflush(stdout);
       }
     }
    printf("\n");fflush(stdout);
   }
#endif

/* Check 3: order 1, nrm */

  error=0.;
  for(i=0;i<m;i++)
   {
    for(j=0;j<m;j++)
     {
      for(k=0;k<m;k++)
       {
        sum=0.;
        for(p=0;p<n;p++)sum+=du(p,i)*ddu(p,j,k);
        error+=fabs(sum);
        if(fabs(sum)>1.e-7)printf("  u^p_{,%d}u^p_{,%d,%d) of order 1, nrm eq. error=%le\n",i,j,k,sum);
        for(p=0;p<n;p++)
          if(fabs(dda(i,j,k)-dda(i,k,j))>1.e-7)printf("  order 2, dda^%d_,%d,%d!(%lf)=dda^%d_,%d,%d(%lf)\n",i,j,k,dda(i,j,k),i,k,j,dda(i,k,j));
       }
     }
   }
  if(error>1.e-5)printf("Total error in order 1, nrm eqs. %le\n",error);
   else          printf("Passed test on order 1, nrm eqs.\n");

/* Check 4: order 2, fn */

  for(i=0;i<n-m;i++)
   {
    for(j=0;j<m;j++)
     {
      for(k=0;k<m;k++)
       {
        sum=0.;

        for(p=0;p<n;p++)
          for(q=0;q<n;q++)
            for(r=0;r<n;r++)
              sum+=dv(p,i)*ddF(p,q,r)*du(q,j)*du(r,k);

        for(w=0;w<n-m;w++)
         {
          dot=0.;
          for(p=0;p<n;p++)
            for(q=0;q<n;q++)
              dot+=dv(p,i)*dF(p,q)*dv(q,w);
          dot1=0.;
          for(r=0;r<n;r++)
            dot1+=dv(r,w)*ddu(r,j,k);
          sum+=dot*dot1;
  
         }
        for(q=0;q<m;q++)
         {
          dot=0;
          for(p=0;p<n;p++)
            dot+=dv(p,i)*ddu(p,q,j);
          sum-=da(q,k)*dot;

          dot=0;
          for(p=0;p<n;p++)
            dot+=dv(p,i)*ddu(p,q,k);
          sum-=da(q,j)*dot;

          dot=0;
          for(p=0;p<n;p++)
            dot+=dv(p,i)*du(p,q);
          sum-=dda(q,j,k)*dot;
         }
        if(fabs(sum)>1.e-7)printf("  v,(%d,%d,%d) of order 2, fn eq. error=%le\n",i,j,k,sum);
       }
     }
   }


  for(i=0;i<m;i++)
   {
    for(j=0;j<m;j++)
     {
      for(k=0;k<m;k++)
       {
        sum=0.;
        for(p=0;p<n;p++)
          for(q=0;q<n;q++)
            for(r=0;r<n;r++)
              sum+=du(p,i)*ddF(p,q,r)*du(q,j)*du(r,k);

        for(p=0;p<n;p++)
          for(q=0;q<n;q++)
            sum+=du(p,i)*dF(p,q)*ddu(q,j,k);
  
        for(p=0;p<n;p++)
         {
          for(q=0;q<m;q++)
           {
            sum-=du(p,i)*ddu(p,q,j)*da(q,k);
            sum-=du(p,i)*ddu(p,q,k)*da(q,j);
            sum-=du(p,i)*du(p,q)*dda(q,j,k);
           }
         }
        if(fabs(sum)>1.e-7)printf("  u,(%d,%d,%d) of order 2, fn eq. error=%le\n",i,j,k,sum);
       }
     }
   }

  error=0.;
  for(i=0;i<n;i++)
   {
    for(j=0;j<m;j++)
     {
      for(k=0;k<m;k++)
       {
        sum=0.;
        for(p=0;p<n;p++)
          for(q=0;q<n;q++)sum+=ddF(i,p,q)*du(p,j)*du(q,k);

        for(p=0;p<n;p++)sum+=dF(i,p)*ddu(p,j,k);
        for(p=0;p<m;p++)sum-=ddu(i,p,j)*da(p,k);
        for(p=0;p<m;p++)sum-=ddu(i,p,k)*da(p,j);
        for(p=0;p<m;p++)sum-=du(i,p)*dda(p,j,k);
        error+=fabs(sum);

        if(fabs(sum)>1.e-7)printf("  (%d,%d,%d) of order 2, fn eq. error=%le\n",i,j,k,sum);

        if(fabs(ddu(i,j,k)-ddu(i,k,j))>1.e-7)printf("  order 2, ddu^%d_,%d,%d!(%lf)=ddu^%d_,%d,%d(%lf)\n",i,j,k,ddu(i,j,k),i,k,j,ddu(i,k,j));
       }
     }
   }
  if(error>1.e-5)printf("Total error in order 2, fn  eqs. %le\n",error);
   else          printf("Passed test on order 2, fn eqs.\n");


/* ------------------------------------------------------- */

  AA=(double*)malloc((n-m)*m*m*m*(n-m)*m*m*m*sizeof(double));

#ifndef MFNOSAFETYNET
  if(AA==NULL)
   {
    sprintf(IMFFixedPointErrorMsg,"Out of memory, trying to allocate %d bytes",(n-m)*m*m*m*(n-m)*m*m*m*sizeof(double));
    MFSetError(e,12,RoutineName,IMFFixedPointErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  bb=(double*)malloc((n-m)*m*m*m*sizeof(double));

#ifndef MFNOSAFETYNET
  if(bb==NULL)
   {
    sprintf(IMFFixedPointErrorMsg,"Out of memory, trying to allocate %d bytes",(n-m)*m*m*m*sizeof(double));
    MFSetError(e,12,RoutineName,IMFFixedPointErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  for(i=0;i<(n-m)*m*m*m*(n-m)*m*m*m;i++)AA[i]=0;
  for(i=0;i<(n-m)*m*m*m;i++)bb[i]=0;

  for(i=0;i<n-m;i++)
   {
    for(j=0;j<m;j++)
     {
      for(k=0;k<m;k++)
       {
        for(l=0;l<m;l++)
         {
          B1(i,j,k,l)=0.;
          for(w=0;w<n;w++)
           {
            for(p=0;p<n;p++)
             {
              for(q=0;q<n;q++)
               {
                for(r=0;r<n;r++)
                 {
                  B1(i,j,k,l)+=dv(w,i)*dddF(w,p,q,r)*du(p,j)*du(q,k)*du(r,l);
                 }
                B1(i,j,k,l)+=dv(w,i)*ddF(w,p,q)*(du(q,k)*ddu(p,j,l)+du(q,l)*ddu(p,j,k)+du(q,j)*ddu(p,k,l));
               }
             }
            for(p=0;p<m;p++)
              B1(i,j,k,l)-=dv(w,i)*ddu(w,p,j)*dda(p,k,l)+dv(w,i)*ddu(w,p,l)*dda(p,j,k)+dv(w,i)*ddu(w,p,k)*dda(p,j,l);
           }
    
          for(p=0;p<m;p++)
           {
            A1(i,j,k,l,i,p,j,k)+=dudFdu(p,l);
            A1(i,j,k,l,i,p,k,l)+=dudFdu(p,j);
            A1(i,j,k,l,i,p,j,l)+=dudFdu(p,k);
           }
          for(w=0;w<n-m;w++)
            A1(i,j,k,l,w,j,k,l)-=dvdFdv(i,w);
         }
       }
     }
   }

  if(0){printf("\nSystem for third derivatives\n");
        IMFPrintFull((n-m)*m*m*m,AA,bb,e);fflush(stdout);}

  if(!MFSolveFull((n-m)*m*m*m,AA,bb,e))
   {
    printf("The full solve for the third derivatives came back singular\n");
    fflush(stdout);
    exit(12);
   }

  if(0){
   double t;
   printf("\nSecond derivatives ddu^i_{jk}\n");
   for(i=0;i<n;i++)
   for(j=0;j<m;j++)
   for(k=0;k<m;k++)
     printf("ddu^%d_{%d %d}=%le\n",i,j,k,ddu(i,j,k));

   printf("\nConnection in the normal space eta^i_{jk}\n");
   for(i=0;i<n-m;i++)
   for(j=0;j<m;j++)
   for(k=0;k<m;k++)
    {
     t=0;for(p=0;p<n;p++)t+=dv(p,i)*ddu(p,j,k);
     printf("eta^%d_{%d %d}=%le (%le)\n",i,j,k,t,eta(i,j,k));
    }
   fflush(stdout);

   printf("\nCurvature in the normal space eta^i_{jkl}\n");
   for(i=0;i<n-m;i++)
   for(j=0;j<m;j++)
   for(k=0;k<m;k++)
   for(l=0;l<m;l++)
     printf("eta^%d_{%d %d %d}=%le\n",i,j,k,l,B1(i,j,k,l));
   fflush(stdout);
  }

  printf("\nThird derivatives\n");
  for(p=0;p<n-m;p++)
   {
    for(j=0;j<m;j++)
     {
      for(k=j;k<m;k++)
       {
        l=k;
        sum=0.;for(q=0;q<n-m;q++)sum-=eta(q,p,l)*eta(q,j,k);

        if(k==j)
          {if(fabs(sum)>1.e-7)printf("   %le s_%d s_%d s_%d u_{,%d}\n",sum/6.,j,j,j,p);
           if(fabs(B1(p,j,j,j))>1.e-7)printf("   %le s_%d s_%d s_%d v_{,%d}\n",B1(p,j,j,j)/6.,j,j,j,p);}
         else
          {if(fabs(sum)>1.e-7)printf("   %le s_%d s_%d s_%d u_{,%d}\n",sum/2.,j,k,k,p);
           if(fabs(B1(p,j,k,k))>1.e-7)printf("   %le s_%d s_%d s_%d v_{,%d}\n",B1(p,j,k,k)/2.,j,k,k,p);}

        for(l=k+1;l<m;l++)
         {
          sum=0.;for(q=0;q<n-m;q++)sum-=eta(q,p,l)*eta(q,j,k);
          if(k==j)
           {if(fabs(sum)>1.e-7)printf("   %le s_%d s_%d s_%d u_{,%d}\n",sum/2.,j,j,l,p);
            if(fabs(B1(p,j,j,l))>1.e-7)printf("   %le s_%d s_%d s_%d v_{,%d}\n",B1(p,j,k,k)/2.,j,j,l,p);}
          else
           {if(fabs(sum)>1.e-7)printf("   %le s_%d s_%d s_%d u_{,%d}\n",sum/6.,j,k,l,p);
            if(fabs(B1(p,j,k,l))>1.e-7)printf("   %le s_%d s_%d s_%d v_{,%d}\n",B1(p,j,k,l)/6.,j,k,l,p);}
         }
       }
     }
   }

#ifdef MFALLOWVERBOSE
  if(verboser)
   {
    printf("\n");fflush(stdout);
    for(j=0;j<m;j++)
     {
      for(k=0;k<m;k++)
       {
        printf("eta_{%d%d} = (%le",j,k,eta(0,j,k));fflush(stdout);
        for(i=1;i<n-m;i++)printf(",%le",eta(i,j,k));fflush(stdout);
        printf(")\n");fflush(stdout);
       }
     }
    printf("\n");fflush(stdout);
    for(j=0;j<m;j++)
     {
      for(k=0;k<m;k++)
       {
        for(l=0;l<m;l++)
         {
          printf("eta_{%d%d%d} = (%le",j,k,l,B1(0,j,k,l));fflush(stdout);
          for(i=1;i<n-m;i++)printf(",%le",B1(i,j,k,l));fflush(stdout);
          printf(")\n");fflush(stdout);
         }
       }
     }
    printf("\n");fflush(stdout);
   }
#endif

  for(i=0;i<n;i++)
   {
    for(j=0;j<m;j++)
     {
      for(k=0;k<m;k++)
       {
        for(l=0;l<m;l++)
         {
          dddu(i,j,k,l)=0.;
          for(p=0;p<n-m;p++)
            dddu(i,j,k,l)+=B1(p,j,k,l)*dv(i,p);
          for(p=0;p<n;p++)
            for(q=0;q<m;q++)
             {
              dddu(i,j,k,l)-=ddu(p,q,l)*ddu(p,j,k)*du(i,q)/3.*3.;
/*            dddu(i,j,k,l)-=ddu(p,q,j)*ddu(p,k,l)*du(i,q)/3.;
              dddu(i,j,k,l)-=ddu(p,q,k)*ddu(p,j,l)*du(i,q)/3.;*/
             }
         }
       }
     }
   }

  if(0){
   double t;

   printf("\nThird derivatives u^i_{jkl}\n");
   for(i=0;i<n;i++)
   for(j=0;j<m;j++)
   for(k=0;k<m;k++)
   for(l=0;l<m;l++)
     printf("u^%d_{%d %d %d}=%le\n",i,j,k,l,dddu(i,j,k,l));

    printf("Connection in the tangent space - G_{ijk}\n");
    for(i=0;i<m;i++)
    for(j=0;j<m;j++)
    for(k=0;k<m;k++)
     {
      t=0.;
      for(p=0;p<n;p++)
        t+=du(p,i)*ddu(p,j,k);
      printf("   G^%d_{%d %d}=%le\n",i,j,k,t);
     }

    printf("Curvature in the tangent space - G_{ijkl}\n");
    for(i=0;i<m;i++)
    for(j=0;j<m;j++)
    for(k=0;k<m;k++)
    for(l=0;l<m;l++)
     {
      t=0.;
      for(p=0;p<n;p++)
        t+=du(p,i)*dddu(p,j,k,l);
      printf("   G^%d_{%d %d %d}=%le ",i,j,k,l,t);
      t=0.;
      for(p=0;p<n-m;p++)
        t+=eta(p,i,l)*eta(p,j,k);
      printf("   (-eta^p_{%d %d} eta^p_{%d %d}=%le)\n",i,l,j,k,-t);
     }
   }

#ifdef MFALLOWVERBOSE
  if(verboser)
   {
    for(j=0;j<m;j++)
     {
      for(k=0;k<m;k++)
       {
        for(l=0;l<m;l++)
         {
          printf("u_{,%d,%d,%d} = (%le",j,k,l,dddu(0,j,k,l));fflush(stdout);
          for(i=1;i<n;i++)printf(",%le",dddu(i,j,k,l));fflush(stdout);
          printf(")\n");fflush(stdout);
         }
       }
     }
    printf("\n");fflush(stdout);
    for(j=0;j<m;j++)
     {
      for(k=0;k<m;k++)
       {
        for(l=0;l<m;l++)
         {
          printf("F_{,%d,%d,%d} = (%le",j,k,l,dddF(0,j,k,l));fflush(stdout);
          for(i=1;i<n;i++)printf(",%le",dddF(i,j,k,l));fflush(stdout);
          printf(")\n");fflush(stdout);
         }
       }
     }
    printf("\n");fflush(stdout);
   }
#endif

  for(i=0;i<m;i++)
   {
    for(j=0;j<m;j++)
     {
      for(k=0;k<m;k++)
       {
        for(l=0;l<m;l++)
         {
          ddda(i,j,k,l)=0.;
          for(p=0;p<n;p++)
           {
            for(q=0;q<n;q++)
             {
              for(r=0;r<n;r++)
               {
                for(w=0;w<n;w++)ddda(i,j,k,l)+=du(w,i)*dddF(w,p,q,r)*du(p,j)*du(q,k)*du(r,l);
                ddda(i,j,k,l)+=du(p,i)*ddF(p,q,r)*(du(r,k)*ddu(q,j,l)+du(r,l)*ddu(q,j,k)+du(r,j)*ddu(q,k,l));
               }
              ddda(i,j,k,l)+=du(p,i)*dF(p,q)*dddu(q,j,k,l);
             }
            for(q=0;q<m;q++)
             {
              ddda(i,j,k,l)-=du(p,i)*dddu(p,q,j,k)*da(q,l);
              ddda(i,j,k,l)-=du(p,i)*dddu(p,q,k,l)*da(q,j);
              ddda(i,j,k,l)-=du(p,i)*dddu(p,q,j,l)*da(q,k);
              ddda(i,j,k,l)-=du(p,i)*ddu(p,q,j)*dda(q,k,l);
              ddda(i,j,k,l)-=du(p,i)*ddu(p,q,k)*dda(q,j,l);
              ddda(i,j,k,l)-=du(p,i)*ddu(p,q,l)*dda(q,j,k);
             }
           }
         }
       }
     }
   }

#ifdef MFALLOWVERBOSE
  if(verboser)
   {
    for(j=0;j<m;j++)
     {
      for(k=0;k<m;k++)
       {
        for(l=0;l<m;l++)
         {
          printf("a_{,%d,%d,%d} = (%le",j,k,l,ddda(0,j,k,l));fflush(stdout);
          for(i=1;i<m;i++)printf(",%le",ddda(i,j,k,l));fflush(stdout);
          printf(")\n");fflush(stdout);
         }
       }
     }
    printf("\n");fflush(stdout);
   }
#endif

/* Check 5: order 2, nrm */

  error=0.;
  for(i=0;i<m;i++)
   {
    for(j=0;j<m;j++)
     {
      for(k=0;k<m;k++)
       {
        for(l=0;l<m;l++)
         {
          sum=0.;
          for(p=0;p<n;p++)sum+=du(p,i)*dddu(p,j,k,l)+ddu(p,i,l)*ddu(p,j,k);
          error+=fabs(sum);
          if(fabs(sum)>1.e-7)printf("  u^p_{,%d}u^p_{,%d,%d,%d}+u^p_{,%d,%d}u^p_{,%d,%d} of order 2, nrm eq. error=%le\n",i,j,k,l,i,l,j,k,sum);
          for(p=0;p<n;p++)
           {
            if(fabs(ddda(i,j,k,l)-ddda(i,k,j,l))>1.e-7)printf("  order 2, ddda^%d_,%d,%d,%d!(%lf)=ddda^%d_,%d,%d,%d(%lf)\n",i,j,k,l,ddda(i,j,k,l),i,k,j,l,ddda(i,k,j,l));
            if(fabs(ddda(i,j,k,l)-ddda(i,l,k,j))>1.e-7)printf("  order 2, ddda^%d_,%d,%d,%d!(%lf)=ddda^%d_,%d,%d,%d(%lf)\n",i,j,k,l,ddda(i,j,k,l),i,l,k,j,ddda(i,l,k,j));
            if(fabs(ddda(i,j,k,l)-ddda(i,j,l,k))>1.e-7)printf("  order 2, ddda^%d_,%d,%d,%d!(%lf)=ddda^%d_,%d,%d,%d(%lf)\n",i,j,k,l,ddda(i,j,k,l),i,j,l,k,ddda(i,j,l,k));
           }
         }
       }
     }
   }
  if(error>1.e-5)printf("Total error in order 2, nrm eqs. %le\n",error);
   else          printf("Passed test on order 2, nrm eqs.\n");

/* Check 6: order 3, fn */

  error=0.;
  for(i=0;i<n;i++)
   {
    for(j=0;j<m;j++)
     {
      for(k=0;k<m;k++)
       {
        for(l=0;l<m;l++)
         {
          sum=0.;
          for(p=0;p<n;p++)
            for(q=0;q<n;q++)
              for(r=0;r<n;r++)sum+=dddF(i,p,q,r)*du(p,j)*du(q,k)*du(r,l);
          for(p=0;p<n;p++)
            for(q=0;q<n;q++)sum+=ddF(i,p,q)*(ddu(p,j,l)*du(q,k)+ddu(p,j,k)*du(q,l)+du(p,j)*ddu(q,k,l));
          for(p=0;p<n;p++)sum+=dF(i,p)*dddu(p,j,k,l);
          for(p=0;p<m;p++)sum+=-dddu(i,p,j,k)*da(p,l)
                               -dddu(i,p,k,l)*da(p,j)
                               -dddu(i,p,j,l)*da(p,k)
                               -ddu(i,p,j)*dda(p,k,l)
                               -ddu(i,p,l)*dda(p,j,k)
                               -ddu(i,p,k)*dda(p,j,l)
                               -du(i,p)*ddda(p,j,k,l);
          error+=fabs(sum);
          if(fabs(sum)>1.e-7)printf("  (%d,%d,%d,%d) of order 3, fn eq. error=%le\n",i,j,k,l,sum);
          for(p=0;p<n;p++)
           {
            if(fabs(dddu(i,j,k,l)-dddu(i,k,j,l))>1.e-7)printf("  order 2, dddu^%d_,%d,%d,%d!(%lf)=dddu^%d_,%d,%d,%d(%lf)\n",i,j,k,l,dddu(i,j,k,l),i,k,j,l,dddu(i,k,j,l));
            if(fabs(dddu(i,j,k,l)-dddu(i,l,k,j))>1.e-7)printf("  order 2, dddu^%d_,%d,%d,%d!(%lf)=dddu^%d_,%d,%d,%d(%lf)\n",i,j,k,l,dddu(i,j,k,l),i,l,k,j,dddu(i,l,k,j));
            if(fabs(dddu(i,j,k,l)-dddu(i,j,l,k))>1.e-7)printf("  order 2, dddu^%d_,%d,%d,%d!(%lf)=dddu^%d_,%d,%d,%d(%lf)\n",i,j,k,l,dddu(i,j,k,l),i,j,l,k,dddu(i,j,l,k));
           }
         }
       }
     }
   }
  if(error>1.e-5)printf("Total error in order 3, fn  eqs. %le\n",error);
   else          printf("Passed test on order 3, fn eqs.\n");
  fflush(stdout);

  IMFExpansionSetDerivatives(U,u0,du,ddu,dddu,e);
  IMFExpansionSetDerivatives(a,a0,da,dda,ddda,e);

  free(a0);
  free(da);
  free(dda);
  free(ddda);
  free(ddu);
  free(dddu);

  free(A);
  free(b);
  free(AA);
  free(bb);
  free(dudFdu);
  free(dvdFdu);
  free(dvdFdv);
  free(dvddFdudu);
  free(dF);
  free(ddF);
  free(dddF);

  free(duddFdudu);
  free(dudFddu);
  printf("done %s\n",RoutineName);fflush(stdout);
  return;
 }

#ifdef HAVE_LAPACK
void CALLDGESVD(char*,char*,int*,int*,double*,int*,double*,double*,int*,double*,int*,double*,int*,int*);
#endif

/*! \fn MFNKMatrix IMFGetBasisForStableInvariantSubspace(IMFFlow F, MFNVector u, MFKVector p, MFErrorHandler e);
 *  \brief Given a hyperbolic fixed point in a flow, finds a basis for the 
 *            stable invariant linear subspace of a hyperbolic fixed point.
 *
 *  \param F The flow.
 *  \param u The fixed point. 
 *  \param p The parameters of the flow for the fixed point.
 *  \param e An error handler.
 *  \returns A new MFNKMatrix containing an orthonormal basis for the stable invariant linear subspace of the hyperbolic fixed
 *             point.
 */
MFNKMatrix IMFGetBasisForStableInvariantSubspace(IMFFlow F, MFNVector u, MFKVector p, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFGetBasisForStableInvariantSubspace"};
  int i,j,l;
  double R;
  double ev,evmax;
  int info;
  char jobvl;
  char jobvr;
  static double *a=NULL;
  int lda;
  double *wr=NULL;
  double *wi=NULL;
  int ldvl;
  double *vl=NULL;
  int ldvr;
  double *vr=NULL;
  double *work=NULL;
  int lwork;
  double dot;
  int n;
  double *dF;
  int k;
  double *du;
  MFNKMatrix Du;
  int verbose=0;

#ifdef HAVE_LAPACK

#ifdef MFALLOWVERBOSE
  if(verbose){printf("\n\n%s\n",RoutineName);}
#endif

  n=MFNV_NC(u,e);
  lda=n;
  dF=(double*)malloc(n*n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(dF==NULL)
   {
    sprintf(IMFFixedPointErrorMsg,"Out of memory, trying to allocate %d bytes",n*n*sizeof(double));
    MFSetError(e,12,RoutineName,IMFFixedPointErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  IMFEvaluateDerivativeOfFlow(F,u,p,dF,e);

  a=(double*)realloc((void*)a,n*n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(a==NULL)
   {
    sprintf(IMFFixedPointErrorMsg,"Out of memory, trying to allocate %d bytes",n*n*sizeof(double));
    MFSetError(e,12,RoutineName,IMFFixedPointErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  for(i=0;i<n*n;i++)a[i]=dF[i];

  wr=(double*)realloc((void*)wr,n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(wr==NULL)
   {
    sprintf(IMFFixedPointErrorMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(double));
    MFSetError(e,12,RoutineName,IMFFixedPointErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  wi=(double*)realloc((void*)wi,n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(wi==NULL)
   {
    sprintf(IMFFixedPointErrorMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(double));
    MFSetError(e,12,RoutineName,IMFFixedPointErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  ldvl=n;
  vl=(double*)realloc((void*)vl,n*n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(vl==NULL)
   {
    sprintf(IMFFixedPointErrorMsg,"Out of memory, trying to allocate %d bytes",n*n*sizeof(double));
    MFSetError(e,12,RoutineName,IMFFixedPointErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  ldvr=n;
  vr=(double*)realloc((void*)vr,n*n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(vr==NULL)
   {
    sprintf(IMFFixedPointErrorMsg,"Out of memory, trying to allocate %d bytes",n*n*sizeof(double));
    MFSetError(e,12,RoutineName,IMFFixedPointErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  lwork=4*n;
  work=(double*)realloc((void*)work,lwork*sizeof(double));

#ifndef MFNOSAFETYNET
  if(work==NULL)
   {
    sprintf(IMFFixedPointErrorMsg,"Out of memory, trying to allocate %d bytes",lwork*sizeof(double));
    MFSetError(e,12,RoutineName,IMFFixedPointErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  for(i=0;i<lwork;i++)work[i]=0.;

  a=(double*)realloc((void*)a,n*n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(a==NULL)
   {
    sprintf(IMFFixedPointErrorMsg,"Out of memory, trying to allocate %d bytes",n*n*sizeof(double));
    MFSetError(e,12,RoutineName,IMFFixedPointErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  for(i=0;i<n*n;i++)a[i]=dF[i];

  info=0;
  jobvl='N';
  jobvr='V';
  CALLDGEEV(&jobvl,&jobvr,&n,a,&lda,wr,wi,vl,&ldvl,vr,&ldvr,work,&lwork,&info);

  k=0;
  for(i=0;i<n;i++)
    if(wr[i]<-1.e-7)k++;

  du=(double*)malloc(n*k*sizeof(double));

#ifndef MFNOSAFETYNET
  if(du==NULL)
   {
    sprintf(IMFFixedPointErrorMsg,"Out of memory, trying to allocate %d bytes",n*k*sizeof(double));
    MFSetError(e,12,RoutineName,IMFFixedPointErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif


  k=0;
  for(i=0;i<n;i++)
   {

#ifdef MFALLOWVERBOSE
    if(verbose)
     {
      printf("Eigenvalue %d is (%lf,%lf) ",i,wr[i],wi[i]);
      printf("Eigenvector is (%lf",vr[0+n*i]);
      for(j=1;j<n;j++)printf(",%lf",vr[j+n*i]);
      printf(")\n");
      fflush(stdout);
     }
#endif

    if(wr[i]<-1.e-7)
     {
      for(j=0;j<n;j++)du[j+n*k]=vr[j+n*i];
      k++;
     }
   }

  IMFOrthonormalizeBasis(n,k,du,e);

  Du=MFCreateNKMatrixWithData(n,k,du,e);

  return Du;
#else
  sprintf(IMFFixedPointErrorMsg,"The IMFFixedPt requires dgeev from Lapack");
  MFSetError(e,12,RoutineName,IMFFixedPointErrorMsg,__LINE__,__FILE__);
  return (MFNKMatrix)NULL;
#endif
 }

/*! \fn MFNKMatrix IMFGetBasisForUnstableInvariantSubspace(IMFFlow F, MFNVector u, MFKVector p, MFErrorHandler e)
 *  \brief Given a hyperbolic fixed point in a flow, finds a basis for the 
 *            unstable invariant linear subspace of a hyperbolic fixed point.
 *
 *  \param F The flow.
 *  \param u The fixed point. 
 *  \param p The parameters of the flow for the fixed point.
 *  \param e An error handler.
 *  \returns A new MFNKMatrix containing an orthonormal basis for the unstable invariant linear subspace of the hyperbolic fixed
 *             point.
 */
MFNKMatrix IMFGetBasisForUnstableInvariantSubspace(IMFFlow F, MFNVector u, MFKVector p, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFGetBasisForUnstableInvariantSubspace"};
  int i,j,l;
  double R;
  double ev,evmax;
  int info;
  char jobvl;
  char jobvr;
  static double *a=NULL;
  int lda;
  double *wr=NULL;
  double *wi=NULL;
  int ldvl;
  double *vl=NULL;
  int ldvr;
  double *vr=NULL;
  double *work=NULL;
  int lwork;
  double dot;
  int n;
  double *dF;
  int k;
  double *du;
  MFNKMatrix Du;
  int verbose=0;

#ifdef HAVE_LAPACK

#ifdef MFALLOWVERBOSE
  if(verbose){printf("\n\n%s\n",RoutineName);}
#endif

  n=MFNV_NC(u,e);
  lda=n;
  dF=(double*)malloc(n*n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(dF==NULL)
   {
    sprintf(IMFFixedPointErrorMsg,"Out of memory, trying to allocate %d bytes",n*n*sizeof(double));
    MFSetError(e,12,RoutineName,IMFFixedPointErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  IMFEvaluateDerivativeOfFlow(F,u,p,dF,e);

  a=(double*)realloc((void*)a,n*n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(a==NULL)
   {
    sprintf(IMFFixedPointErrorMsg,"Out of memory, trying to allocate %d bytes",n*n*sizeof(double));
    MFSetError(e,12,RoutineName,IMFFixedPointErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  for(i=0;i<n*n;i++)a[i]=dF[i];

  wr=(double*)realloc((void*)wr,n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(wr==NULL)
   {
    sprintf(IMFFixedPointErrorMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(double));
    MFSetError(e,12,RoutineName,IMFFixedPointErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  wi=(double*)realloc((void*)wi,n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(wi==NULL)
   {
    sprintf(IMFFixedPointErrorMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(double));
    MFSetError(e,12,RoutineName,IMFFixedPointErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  ldvl=n;
  vl=(double*)realloc((void*)vl,n*n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(vl==NULL)
   {
    sprintf(IMFFixedPointErrorMsg,"Out of memory, trying to allocate %d bytes",n*n*sizeof(double));
    MFSetError(e,12,RoutineName,IMFFixedPointErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif


  ldvr=n;
  vr=(double*)realloc((void*)vr,n*n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(vr==NULL)
   {
    sprintf(IMFFixedPointErrorMsg,"Out of memory, trying to allocate %d bytes",n*n*sizeof(double));
    MFSetError(e,12,RoutineName,IMFFixedPointErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  lwork=4*n;
  work=(double*)realloc((void*)work,lwork*sizeof(double));

#ifndef MFNOSAFETYNET
  if(work==NULL)
   {
    sprintf(IMFFixedPointErrorMsg,"Out of memory, trying to allocate %d bytes",lwork*sizeof(double));
    MFSetError(e,12,RoutineName,IMFFixedPointErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  for(i=0;i<lwork;i++)work[i]=0.;

  info=0;
  jobvl='N';
  jobvr='V';
  CALLDGEEV(&jobvl,&jobvr,&n,a,&lda,wr,wi,vl,&ldvl,vr,&ldvr,work,&lwork,&info);

  k=0;
  for(i=0;i<n;i++)
    if(wr[i]>1.e-7)k++;

  du=(double*)malloc(n*k*sizeof(double));

#ifndef MFNOSAFETYNET
  if(du==NULL)
   {
    sprintf(IMFFixedPointErrorMsg,"Out of memory, trying to allocate %d bytes",n*k*sizeof(double));
    MFSetError(e,12,RoutineName,IMFFixedPointErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif


  k=0;
  for(i=0;i<n;i++)
   {

#ifdef MFALLOWVERBOSE
    if(verbose){printf("Eigenvalue %d is (%lf,%lf) ",i,wr[i],wi[i]);
    printf("Eigenvector is (%lf",vr[0+n*i]);
    for(j=1;j<n;j++)printf(",%lf",vr[j+n*i]);
    printf(")\n");}
#endif

    if(wr[i]>1.e-7)
     {
      for(j=0;j<n;j++)du[j+n*k]=vr[j+n*i];
      k++;
     }
   }

  IMFOrthonormalizeBasis(n,k,du,e);

  Du=MFCreateNKMatrixWithData(n,k,du,e);

  free(du);
  free(dF);
  free(a);
  free(wr);
  free(wi);
  free(vl);
  free(vr);
  free(work);

  return Du;
#else
  sprintf(IMFFixedPointErrorMsg,"The IMFFixedPt requires dgeev from Lapack");
  MFSetError(e,12,RoutineName,IMFFixedPointErrorMsg,__LINE__,__FILE__);
  return (MFNKMatrix)NULL;
#endif
 }

void IMFOrthonormalizeBasis(int n,int m,double *du, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFOrthonormalizeBasis"};
  double nrm;
  int i,j,k;

  for(j=0;j<m;j++)
   {
    for(k=0;k<j;k++)
     {
  
  /* Orthogonalize columns j and k */
  
      nrm=0;for(i=0;i<n;i++)nrm+=du(i,k)*du(i,j);
      for(i=0;i<n;i++)du(i,j)=du(i,j)-nrm*du(i,k);
     }

  /* Normalize column k */
  
    nrm=0;for(i=0;i<n;i++)nrm+=du(i,j)*du(i,j);nrm=sqrt(nrm);
          for(i=0;i<n;i++)du(i,j)=du(i,j)/nrm;
   }

  return;
 }

/*! \fn MFNKMatrix IMFGetBasisForOrthogonalComplement(MFNKMatrix Du, MFErrorHandler e);
 *  \brief Given a linear subspace, finds an orthonormal basis for its orthogonal complement.
 *
 *  \param Du An orthonormal basis a linear subspace.
 *  \param e  An error handler.
 *  \returns The orthogonal complement.
 */
MFNKMatrix IMFGetBasisForOrthogonalComplement(MFNKMatrix Du, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFGetBasisForOrthogonalComplement"};
  int i,j,m;
  int info;
  char jobvl;
  char jobvr;
  static double *a=NULL;
  int lda;
  double *s=NULL;
  int ldu;
  double *u=NULL;
  int ldv;
  double *v=NULL;
  double *work=NULL;
  int lwork;
  double sum;
  double t;
  int l;
  int n;
  int k;
  double *du;
  double *dv;
  MFNKMatrix Dv;
  int verbose=0;

#ifdef HAVE_LAPACK

  n=MFNKMatrixN(Du,e);
  k=MFNKMatrixK(Du,e);

#ifdef MFNOCONFIDENCE
  if(n==k)
   { 
    sprintf(IMFFixedPointErrorMsg,"The orthogonal complement is empty");
    MFSetError(e,12,RoutineName,IMFFixedPointErrorMsg,__LINE__,__FILE__);
    return NULL;
   }
#endif

  du=MFNKM_CStar(Du,e);
  dv=(double*)malloc(n*(n-k)*sizeof(double));

#ifndef MFNOSAFETYNET
  if(dv==NULL)
   {
    sprintf(IMFFixedPointErrorMsg,"Out of memory, trying to allocate %d bytes",n*(n-k)*sizeof(double));
    MFSetError(e,12,RoutineName,IMFFixedPointErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

#ifdef MFALLOWVERBOSE
  if(verbose){printf("\n\n%s, n=%d, k=%d\n\n",RoutineName,n,k);fflush(stdout);}
#endif

/*if(n==3&&k==2)
   {
    dv[0]=du[1]*du[5]-du[4]*du[2];
    dv[1]=du[2]*du[3]-du[5]*du[0];
    dv[2]=du[0]*du[4]-du[3]*du[1];
    return;
   }*/

  lda=n;

  a=(double*)realloc((void*)a,n*n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(a==NULL)
   {
    sprintf(IMFFixedPointErrorMsg,"Out of memory, trying to allocate %d bytes",n*n*sizeof(double));
    MFSetError(e,12,RoutineName,IMFFixedPointErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  for(i=0;i<n*n;i++)a[i]=0.;

  for(i=0;i<k;i++)
    for(j=0;j<n;j++)
      a[i+n*j]=du(j,i);

#ifdef MFALLOWVERBOSE
/*if(verbose){IMFPrintFull(n,a,NULL);}*/
#endif

  s=(double*)realloc((void*)s,n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(s==NULL)
   {
    sprintf(IMFFixedPointErrorMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(double));
    MFSetError(e,12,RoutineName,IMFFixedPointErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif


  ldu=n;
  u=(double*)realloc((void*)u,n*n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(u==NULL)
   {
    sprintf(IMFFixedPointErrorMsg,"Out of memory, trying to allocate %d bytes",n*n*sizeof(double));
    MFSetError(e,12,RoutineName,IMFFixedPointErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  ldv=n;
  v=(double*)realloc((void*)v,n*n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(v==NULL)
   {
    sprintf(IMFFixedPointErrorMsg,"Out of memory, trying to allocate %d bytes",n*n*sizeof(double));
    MFSetError(e,12,RoutineName,IMFFixedPointErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  info=0;
  jobvl='N';
  jobvr='A';
  lwork=-1;
  CALLDGESVD(&jobvl,&jobvr,&n,&n,a,&lda,s,u,&ldu,v,&ldv,&t,&lwork,&info);
  lwork=MAX(4*n,5*n-4);
  lwork=round(t);
  work=(double*)realloc((void*)work,lwork*sizeof(double));

#ifndef MFNOSAFETYNET
  if(work==NULL)
   {
    sprintf(IMFFixedPointErrorMsg,"Out of memory, trying to allocate %d bytes",lwork*sizeof(double));
    MFSetError(e,12,RoutineName,IMFFixedPointErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  for(i=0;i<lwork;i++)work[i]=0.;

  info=0;
  jobvl='N';
  jobvr='A';
  CALLDGESVD(&jobvl,&jobvr,&n,&n,a,&lda,s,u,&ldu,v,&ldv,work,&lwork,&info);

  m=0;
  for(i=0;i<n;i++)
   {

#ifdef MFALLOWVERBOSE
    if(verbose)
     {
      printf("Singular value %d is %le ",i,s[i]);
      printf("Eigenvector is (%lf",v[i+n*0]);
      for(j=1;j<n;j++)printf(",%lf",v[i+n*j]);
      printf(")\n");
      fflush(stdout);
     }
#endif

    if(i>=k)
     {
      for(j=0;j<n;j++)dv[j+n*m]=v[i+n*j];

#ifdef MFALLOWVERBOSE
      if(verbose)
       {
        printf("v_{,%d}=(%lf",m,dv[0+n*m]);
        for(j=1;j<n;j++)printf(",%lf",dv[j+n*m]);printf(")\n");
        fflush(stdout);
       }
#endif

      m++;
     }
   }

  IMFOrthonormalizeBasis(n,n-k,dv,e);

  Dv=MFCreateNKMatrixWithData(n,n-k,dv,e);
  free(dv);
  free(a);
  free(s);
  free(u);
  free(v);
  free(work);

  return Dv;
#else
  sprintf(IMFFixedPointErrorMsg,"The IMFFixedPt requires dgesvd from Lapack");
  MFSetError(e,12,RoutineName,IMFFixedPointErrorMsg,__LINE__,__FILE__);
  return (MFNKMatrix)NULL;
#endif
 }

void IMFPrintFullSchematic(int n, double *A, double *b, MFErrorHandler e)
 {
  int i,j;

  for(i=0;i<n;i++)
   {
    printf("%2d ",i);
    for(j=0;j<n;j++)
     {
      if(j>0)printf(" ");
      if(fabs(A[i+n*j])>1.e-6)printf("*");
       else printf("0");
     }
    if(b!=(double*)0)
     {
      if(fabs(b[i])>1.e-6)printf(" | *");
       else printf(" | 0");
     }
    printf("\n");
   }
  printf("\n");
 }

void IMFPrintFull(int n, double *A, double *b, MFErrorHandler e)
 {
  int i,j;

  for(i=0;i<n;i++)
   {
    printf("%2d ",i);
    for(j=0;j<n;j++)
     {
      if(j>0)printf(" ");
      printf("%le",A[i+n*j]);
     }
    if(b!=(double*)0)printf(" | %le ",b[i]);
    printf("\n");
   }
  printf("\n");
 }

#ifdef __cplusplus
}
#endif
