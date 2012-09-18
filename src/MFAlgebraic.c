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

static char *id="@(#) $Id: MFAlgebraic.c,v 1.13 2011/07/21 17:42:46 mhender Exp $";

static char MFAlgebraicMFErrorHandlerMsg[256]="";

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

#ifdef __cplusplus
 extern "C" {
#endif

int MFFullNRealEV(int,double*,MFErrorHandler);
int MFFullNPosSV(int,int,double*,MFErrorHandler);
int MFFullNZeroSV(int,int,double*,MFErrorHandler);
double MFFullProdEV(int,int,double*,MFErrorHandler);
int MFFullNPosEV(int,int,double*,MFErrorHandler);
int MFFullNUnitEV(int,int,double*,MFErrorHandler);
int  MFSolveFull(int,double*,double*,MFErrorHandler);

void MFPrintNKMatrix(FILE*,MFNKMatrix,MFErrorHandler);

static void MFFreeAlgebraicData(void*,MFErrorHandler);
static int MFProjectAlgebraic(int,int,MFNVector,MFNKMatrix,MFNVector,void*,int*,MFErrorHandler);
static int MFStopAlgebraic(MFImplicitMF,MFNVector,MFNKMatrix,MFNVector,MFNKMatrix,void*,MFErrorHandler);
static int MFTangentAlgebraic(int,int,MFNVector,MFNKMatrix,void*,MFErrorHandler);
static int MFTangentAlgebraicWithGuess(int,int,MFNVector,MFNKMatrix,MFNKMatrix,void*,MFErrorHandler);
static double MFScaleAlgebraic(int,int,MFNVector,MFNKMatrix,void*,MFErrorHandler);
static void MFWriteAlgebraicData(FILE*,void*,MFErrorHandler);
static MFImplicitMF MFReadAlgebraic(FILE*,MFErrorHandler);
static int MFSingularAlgebraic(int,int,MFNVector,MFNKMatrix,MFNVector,void*,MFErrorHandler);

static int MFAlgebraicProjectToSave(MFNVector,double*,void*,MFErrorHandler);
static int MFAlgebraicProjectToDraw(MFNVector,double*,void*,MFErrorHandler);
static int MFAlgebraicProjectForBB(MFNVector,double*,void*,MFErrorHandler);
static void MFAlgebraicSetStability(MFImplicitMF,MFNVector,MFNKMatrix,void*,MFErrorHandler);

MFNVector MFNVectorFactory(MFImplicitMF,MFErrorHandler);
MFNKMatrix MFNKMatrixFactory(MFImplicitMF,MFErrorHandler);

double MFAtlasDet(int,double*,MFErrorHandler);

struct MFAlgebraicData
 {
  int expr;

  char *variables;
  char *expression;
  int n;
  int k;
  ECFn F;
  ECFn *dF;
  ECFn *ddF;
  void (*FSub)(int*,double*,int*,double*,void*,MFErrorHandler);
  void (*dFSub)(int*,double*,int*,double*,void*,MFErrorHandler);
  void (*ddFSub)(int*,double*,int*,double*,void*,MFErrorHandler);
  MFImplicitMF thisMF;
  void *data;
 };

MFImplicitMF MFIMFCreateAlgebraicExpressionWithRadius(char *variables,char *expression,double R, MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFCreateAlgebraicExpressionWithRadius"};
  MFImplicitMF thisMF;
  MFNSpace space;
  ECFn f;
  ECFn df;
  int i,j,n,k;
  struct MFAlgebraicData *data;

  f=ECCreateFunction(variables,expression);
  n=ECFunctionM(f);
  k=ECFunctionM(f)-ECFunctionN(f);

  thisMF=MFIMFCreateBaseClass(n,k,"Algebraic",e);

  space=MFCreateNSpace(n,e);
  MFIMFSetSpace(thisMF,space,e);
  MFFreeNSpace(space,e);

  data=(struct MFAlgebraicData*)malloc(sizeof(struct MFAlgebraicData)); /*done*/

#ifndef MFNOSAFETYNET
  if(data==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct MFAlgebraicData));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(thisMF);
    return NULL;
   }
#endif

  data->variables=(char*)malloc((strlen(variables)+1)*sizeof(char));

#ifndef MFNOSAFETYNET
  if(data->variables==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",(strlen(variables)+1)*sizeof(char));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(thisMF);
    return NULL;
   }
#endif
  strcpy(data->variables,variables);

  data->expression=(char*)malloc((strlen(expression)+1)*sizeof(char));

#ifndef MFNOSAFETYNET
  if(data->variables==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",(strlen(expression)+1)*sizeof(char));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(thisMF);
    free(data->variables);
    return NULL;
   }
#endif

  strcpy(data->expression,expression);

  data->n=n;
  data->k=k;
  data->expr=1;
  data->F=f;
  data->dF=(ECFn*)malloc(n*sizeof(ECFn));
  data->thisMF=thisMF;

#ifndef MFNOSAFETYNET
  if(data->dF==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(ECFn));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(thisMF);
    free(data->variables);
    free(data->expression);
    return NULL;
   }
#endif

  data->ddF=(ECFn*)malloc(n*n*sizeof(ECFn));

#ifndef MFNOSAFETYNET
  if(data->ddF==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*n*sizeof(ECFn));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(thisMF);
    free(data->variables);
    free(data->expression);
    free(data->dF);
    return NULL;
   }
#endif

  for(i=0;i<n;i++)
   {
    df=ECCreateDerivativeOfFunction(f,i);
    data->dF[i]=df;
    for(j=0;j<n;j++)
     {
      data->ddF[i+n*j]=ECCreateDerivativeOfFunction(df,j);
     }
   }

  data->FSub=NULL;
  data->dFSub=NULL;
  data->ddFSub=NULL;
  data->thisMF=thisMF;
  data->data=NULL;

  MFIMFSetData(thisMF,(void*)data,e);
  MFIMFSetFreeData(thisMF,MFFreeAlgebraicData,e);
  MFIMFSetStop(thisMF,MFStopAlgebraic,e);
  MFIMFSetSingular(thisMF,MFSingularAlgebraic,e);
  MFIMFSetProject(thisMF,MFProjectAlgebraic,e);
  MFIMFSetTangent(thisMF,MFTangentAlgebraic,e);
  MFIMFSetTangentWithGuess(thisMF,MFTangentAlgebraicWithGuess,e);
  MFIMFSetScale(thisMF,MFScaleAlgebraic,e);
  MFIMFSetR(thisMF,R,e);
  MFIMFSetWriteData(thisMF,MFWriteAlgebraicData,e);

  MFIMFSetProjectForSave(thisMF,MFAlgebraicProjectToSave,e);
  MFIMFSetProjectForDraw(thisMF,MFAlgebraicProjectToDraw,e);
  MFIMFSetProjectForBB(thisMF,MFAlgebraicProjectForBB,e);

  MFIMFSetVectorFactory(thisMF,MFNVectorFactory,e);
  MFIMFSetMatrixFactory(thisMF,MFNKMatrixFactory,e);

  return thisMF;
 }

MFImplicitMF MFIMFCreateAlgebraicExpression(char *variables,char *expression, MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFCreateAlgebraicExpression"};
  MFImplicitMF thisMF;

  thisMF=MFIMFCreateAlgebraicExpressionWithRadius(variables,expression,.1,e);

  return thisMF;
 }

MFImplicitMF MFIMFCreateAlgebraicSubroutine(int n, int k, void (*FSub)(int*,double*,int*,double*,void*,MFErrorHandler),
                                                          void (*dFSub)(int*,double*,int*,double*,void*,MFErrorHandler),
                                                          void (*ddFSub)(int*,double*,int*,double*,void*,MFErrorHandler),
                                                          void *data,
                                                          MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFCreateAlgebraicSubroutine"};
  MFImplicitMF thisMF;

  thisMF=MFIMFCreateAlgebraicSubroutineWithRadius(n,k,FSub,dFSub,ddFSub,data,.1,e);

  return thisMF;
 }

MFImplicitMF MFIMFCreateAlgebraicSubroutineWithRadius(int n, int k,
                     void (*FSub)(int*,double*,int*,double*,void*,MFErrorHandler),
                     void (*dFSub)(int*,double*,int*,double*,void*,MFErrorHandler),
                     void (*ddFSub)(int*,double*,int*,double*,void*,MFErrorHandler),
                     void *d,
                     double R, MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFCreateAlgebraicSubroutineWithRadius"};
  MFImplicitMF thisMF;
  MFNSpace space;
  int i,j;
  struct MFAlgebraicData *data;

  thisMF=MFIMFCreateBaseClass(n,k,"Algebraic",e);

  space=MFCreateNSpace(n,e);
  MFIMFSetSpace(thisMF,space,e);
  MFFreeNSpace(space,e);

  data=(struct MFAlgebraicData*)malloc(sizeof(struct MFAlgebraicData)); /*done*/
#ifndef MFNOSAFETYNET
  if(data==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct MFAlgebraicData));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(thisMF);
    return NULL;
   }
#endif

  data->FSub=FSub;

#ifdef MFNOCONFIDENCE
  if(data->FSub==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Function (arg 3), must be provided");
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    free(thisMF);
    return NULL;
   }
#endif

  data->dFSub=dFSub;

#ifdef MFALLOWVERBOSE
  if(data->dFSub==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Jacobian (arg 4), must be provided");
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    free(thisMF);
    return NULL;
   }
#endif

  data->ddFSub=ddFSub;

#ifdef MFALLOWVERBOSE
  if(data->ddFSub==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Second Derivative (arg 5), must be provided");
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    free(thisMF);
    return NULL;
   }
#endif

  data->n=n;
  data->k=k;
  data->expr=0;
  data->thisMF=thisMF;
  data->data=d;

  MFIMFSetData(thisMF,(void*)data,e);
  MFIMFSetFreeData(thisMF,MFFreeAlgebraicData,e);
  MFIMFSetStop(thisMF,MFStopAlgebraic,e);
  MFIMFSetSingular(thisMF,MFSingularAlgebraic,e);
  MFIMFSetProject(thisMF,MFProjectAlgebraic,e);
  MFIMFSetTangent(thisMF,MFTangentAlgebraic,e);
  MFIMFSetScale(thisMF,MFScaleAlgebraic,e);
  MFIMFSetR(thisMF,R,e);
  MFIMFSetWriteData(thisMF,MFWriteAlgebraicData,e);

  MFIMFSetProjectForSave(thisMF,MFAlgebraicProjectToSave,e);
  MFIMFSetProjectForDraw(thisMF,MFAlgebraicProjectToDraw,e);
  MFIMFSetProjectForBB(thisMF,MFAlgebraicProjectForBB,e);

  MFIMFSetVectorFactory(thisMF,MFNVectorFactory,e);
  MFIMFSetMatrixFactory(thisMF,MFNKMatrixFactory,e);

  return thisMF;
 }

void MFFreeAlgebraicData(void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFFreeAlgebraicData"};
  struct MFAlgebraicData *data;
  int i,n;

  data=(struct MFAlgebraicData*)d;

  if(data->expr)
   {
    if(data->F!=NULL)ECFreeFunction(data->F);

    n=data->n;
    if(data->dF!=NULL)
     {
      for(i=0;i<n;i++)
        if((data->dF)[i]!=NULL)ECFreeFunction(data->dF[i]);
      free(data->dF);
     }

    if(data->ddF!=NULL)
     {
      for(i=0;i<n*n;i++)
        if((data->ddF)[i]!=NULL)ECFreeFunction(data->ddF[i]);
      free(data->ddF);
     }
   }

  free(data);
  return;
 }

void MFWriteAlgebraicData(FILE *fid,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFWriteAlgebraicData"};
  struct MFAlgebraicData *data;

  data=(struct MFAlgebraicData*)d;

  fprintf(fid,"%d %d %d\n",data->n,data->k,data->expression);
  if(data->expr)
   {
    fprintf(fid,"%d %d\n",strlen(data->variables),strlen(data->expression));
    fprintf(fid,"%s\n",data->variables);
    fprintf(fid,"%s\n",data->expression);
   }
  fflush(fid);
  return;
 }

MFImplicitMF MFReadAlgebraic(FILE *fid, MFErrorHandler e)
 {
  static char RoutineName[]={"MFReadAlgebraic"};
  MFImplicitMF expressionMF=NULL;
  int lnvariables=0;
  int lnexpression=0;
  char *variables;
  char *expression;
  int i;
  int n,k,flag;
  char c=' ';

  fscanf(fid,"%d %d\n",&n,&k,&flag);
  if(flag)
   {
    fscanf(fid,"%d %d\n",&lnvariables,&lnexpression);
    variables=(char*)malloc((lnvariables+1)*sizeof(char));

#ifndef MFNOSAFETYNET
    if(variables==NULL)
     {
      sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",(lnvariables+1)*sizeof(char));
      MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return NULL;
     }
#endif

    expression=(char*)malloc((lnexpression+1)*sizeof(char));

#ifndef MFNOSAFETYNET
    if(expression==NULL)
     {
      sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",(lnvariables+1)*sizeof(char));
      MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      free(variables);
      return NULL;
     }
#endif

    for(i=0;i<lnvariables;i++)fscanf(fid,"%c",variables+i);
    variables[lnvariables]=0x0;
    fscanf(fid,"%c",&c);
    for(i=0;i<lnexpression;i++)fscanf(fid,"%c",expression+i);
    fscanf(fid,"%c",&c);
    expression[lnexpression]=0x0;
   }else{
    sprintf(MFAlgebraicMFErrorHandlerMsg,"You can't read an algebraic manifold defined by subroutines");
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    return NULL;
   }

  expressionMF=MFIMFCreateAlgebraicExpression(variables,expression,e);
  free(variables);
  free(expression);

  return expressionMF;
 }

int MFProjectAlgebraic(int n,int k,MFNVector vu0,MFNKMatrix mPhi,MFNVector vu,void *d,int *index, MFErrorHandler e)
 {
  static char RoutineName[]={"MFProjectAlgebraic"};
  double err;
  double *A=NULL; /* [n*n];*/
  double *a=NULL; /* [n*n];*/
  double *b=NULL; /* [n];  */
  int itimes;
  int i,j,l,m;
  struct MFAlgebraicData *data;
  double *f=NULL;
  double *df=NULL;
  double t,dnorm,unorm;
  double *u0,*Phi,*u;
  double epsilon=1.e-7;
  int verbose=0;

  u0=MFNV_CStar(vu0,e);
  Phi=MFNKM_CStar(mPhi,e);
  u=MFNV_CStar(vu,e);

  data=(struct MFAlgebraicData*)d;

  A=(double*)malloc(n*n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(A==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*n*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0;
   }
#endif

  a=(double*)malloc(n*n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(a==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*n*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(A);
    return 0;
   }
#endif

  b=(double*)malloc(n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(b==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(A);
    free(a);
    return 0;
   }
#endif

  f=(double*)malloc((n-k)*sizeof(double));

#ifndef MFNOSAFETYNET
  if(f==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",(n-k)*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(A);
    free(a);
    free(b);
    return 0;
   }
#endif

  df=(double*)malloc(n*(n-k)*sizeof(double));

#ifndef MFNOSAFETYNET
  if(df==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*(n-k)*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(A);
    free(a);
    free(b);
    free(f);
    return 0;
   }
#endif

  for(i=0;i<n;i++)u[i]=u0[i];

  err=1.;
  itimes=0;

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("%s\n",RoutineName);fflush(stdout);
    printf("u0:(%lf",u[0]);for(l=1;l<n;l++)printf(",%lf",u[l]);printf(")\n");fflush(stdout);
   }
#endif

  while(err>epsilon)
   {
    if(data->expr)
     {
      for(j=0;j<n;j++)
       {
        ECEvaluateFunction(data->dF[j],u,f);
        for(i=0;i<n-k;i++)A[i+n*j]=f[i];
       }
      ECEvaluateFunction(data->F,u,f);
     }else{
      m=n-k;
      (data->FSub)(&n,u,&m,f,data->data,e);
      (data->dFSub)(&n,u,&m,df,data->data,e);
      for(i=0;i<m;i++)
       for(j=0;j<n;j++)A[i+n*j]=df[i+m*j];
     }

#ifdef MFALLOWVERBOSE
    if(0&&verbose){printf("f:(%lf",f[0]);for(l=1;l<n-k;l++)printf(",%lf",f[l]);printf(")\n");fflush(stdout);}
#endif

    for(i=0;i<n-k;i++)b[i]=-f[i];

    for(i=n-k;i<n;i++)b[i]=0.;
    for(j=0;j<n;j++)
     {
      for(i=0;i<n-k;i++)a[i+n*j]=A[i+n*j];
      for(i=n-k;i<n;i++)
       {
        A[i+n*j]=Phi[j+n*(i-n+k)];
        b[i]-=Phi[j+n*(i-n+k)]*(u[j]-u0[j]);
       }
     }

    for(j=0;j<n;j++)
      for(i=0;i<n;i++)
        a[i+n*j]=A[i+n*j];

#ifdef MFALLOWVERBOSE
    if(verbose)
     {
      printf("\n");
      printf("%s A:\n",RoutineName);
      for(i=0;i<n;i++)
       {
        printf(" [%lf",A[i+n*0]);for(l=1;l<n;l++)printf(",%lf",A[i+n*l]);printf("]\n");fflush(stdout);
       }
      printf("u: (%lf",u[0]);for(l=1;l<n;l++)printf(",%lf",u[l]);printf(")\n");fflush(stdout);
      printf("b: (%lf",b[0]);for(l=1;l<n;l++)printf(",%lf",b[l]);printf(")\n");fflush(stdout);
      printf("\n");
     }
#endif

    err=0.;
    for(i=0;i<n;i++)err+=b[i]*b[i];
    err=sqrt(err);

    MFSolveFull(n,A,b,e);

    t=1.;
    dnorm=fabs(b[0]);
    unorm=fabs(u[0]);
    for(i=1;i<n;i++)
     {
      if(fabs(b[i])>dnorm)dnorm=fabs(b[i]);
      if(fabs(u[i])>unorm)unorm=fabs(u[i]);
     }
    if(dnorm>.1*(unorm+1))t=.1*(unorm+1)/dnorm; /* Damping */

#ifdef MFALLOWVERBOSE
    if(verbose){printf("  %d %le %lf Du=%le, u=%le\n",itimes,err,t,dnorm,unorm);fflush(stdout);}
#endif

    for(i=0;i<n;i++)u[i]=u[i]+t*b[i];
    itimes++;
    if(itimes>10)
     {
      free(A);
      free(a);
      free(b);
      free(f);
      free(df);

#ifdef MFALLOWVERBOSE
      if(verbose){printf("done %s. Too many iterations. rc=0\n",RoutineName);fflush(stdout);}
#endif

      return 0;
     }
   }

  *index=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s, niter=%d, index is %d\n",RoutineName,itimes,*index);fflush(stdout);}
#endif

  free(A);
  free(a);
  free(b);
  free(f);
  free(df);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("done %s, rc=1\n",RoutineName);fflush(stdout);}
#endif

  return 1;
 }

int MFTangentAlgebraic(int n,int k,MFNVector vu,MFNKMatrix mPhi,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFTangentAlgebraic"};
  double *A;             /* n*n */
  char jobvl='N';                   /* No left  eigenvectors */
  char jobvr='A';                   /*    right eigenvectors */
  int lda;                        /* Leading dimension of A */
  double *s;    /* n */
  double *U;    /* n*n */
  int ldU;
  double *vr;    /* n*k */
  double *V;    /* n*n */
  int ldV;
  double *work;                     /* Work array */
  int lwork;                     /* Length of work array */
  int info=0;
  int i,j,l,m,mm;
  double t;
  struct MFAlgebraicData *data;
  double *f=NULL;
  double *df=NULL;
  double *u,*Phi;
  int verbose=0;
  int check;
  double *a;

#ifdef HAVE_LAPACK

  check=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s, k=%d, n=%d\n",RoutineName,k,n);fflush(stdout);}
#endif

  Phi=MFNKM_CStar(mPhi,e);
  u=MFNV_CStar(vu,e);

  data=(struct MFAlgebraicData*)d;

#ifdef MFALLOWVERBOSE
  if(verbose&&data->expr){printf("%s, %s\n",data->variables,data->expression);fflush(stdout);}
#endif

  A=(double*)malloc(n*n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(A==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*n*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0;
   }
#endif

  a=(double*)malloc(n*n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(a==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*n*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(A);
    return 0;
   }
#endif

  lda=n;
  ldU=n;

  U=(double*)malloc(n*n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(U==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*n*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(A);
    free(a);
    return 0;
   }
#endif

  ldV=n;
  V=(double*)malloc(n*n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(V==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*n*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(A);
    free(a);
    free(U);
    return 0;
   }
#endif

  vr=(double*)malloc(n*k*sizeof(double));

#ifndef MFNOSAFETYNET
  if(vr==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*k*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(A);
    free(a);
    free(U);
    free(V);
    return 0;
   }
#endif

  f=(double*)malloc((n-k)*sizeof(double));

#ifndef MFNOSAFETYNET
  if(f==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",(n-k)*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(A);
    free(a);
    free(U);
    free(V);
    free(vr);
    return 0;
   }
#endif

  df=(double*)malloc(n*(n-k)*sizeof(double));

#ifndef MFNOSAFETYNET
  if(df==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*(n-k)*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(A);
    free(a);
    free(U);
    free(V);
    free(vr);
    free(f);
    return 0;
   }
#endif

  s=(double*)malloc(n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(s==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(A);
    free(a);
    free(U);
    free(V);
    free(vr);
    free(f);
    free(df);
    return 0;
   }
#endif

  if(data->expr)
   {
    for(j=0;j<n;j++)
     {
      ECEvaluateFunction(data->dF[j],u,f);
      for(i=0;i<n-k;i++)A[i+n*j]=f[i];
     }
   }else{
    m=n-k;
    (data->dFSub)(&n,u,&m,df,data->data,e);
    for(i=0;i<m;i++)
      for(j=0;j<n;j++)A[i+n*j]=df[i+m*j];
   }
  for(j=0;j<n;j++)
    for(i=n-k;i<n;i++)A[i+n*j]=0.;
  for(j=0;j<n;j++)
    for(i=0;i<n;i++)a[i+n*j]=A[i+n*j];

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("u:(%lf",u[0]);for(l=1;l<n;l++)printf(",%lf",u[l]);printf(")\n");fflush(stdout);
    printf("A:\n");
    for(i=0;i<n;i++)
     {
      printf(" [%lf",A[i+n*0]);for(l=1;l<n;l++)printf(",%lf",A[i+n*l]);printf("]\n");fflush(stdout);
     }
    printf("\n");
   }
#endif


#ifdef HAVE_LAPACK
  i=-1;
  info=0;
  jobvl='N';
  jobvr='A';
  CALLDGESVD(&jobvl,&jobvr,&n,&n,A,&lda,s,U,&ldU,V,&ldV,&t,&i,&info);
  lwork=MAX(4*n,5*n-4);
  lwork=round(t);

  work=(double*)malloc(lwork*sizeof(double));

#ifndef MFNOSAFETYNET
  if(work==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",lwork*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0;
   }
#endif

  info=0;
  jobvl='N';
  jobvr='A';
  CALLDGESVD(&jobvl,&jobvr,&n,&n,A,&lda,s,U,&ldU,V,&ldV,work,&lwork,&info);
#else
  sprintf(MFAlgebraicMFErrorHandlerMsg,"Lapack is required for this routine");
  MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
  return 0;
#endif

#ifdef MFALLOWVERBOSE
  if(verbose){printf("\n      Null Vectors\n");fflush(stdout);}
#endif

  m=0;
  for(i=0;i<n;i++)
   {
    if(fabs(s[i])<1.e-7)
     {

#ifdef MFALLOWVERBOSE
      if(verbose){printf("     Phi[%d]=v[%d]=(%lf",m,i,V[0+n*i]);for(j=1;j<n;j++)printf(",%lf",V[i+n*j]);printf(")\n");fflush(stdout);}
#endif
      if(m<k)
       {
        for(j=0;j<n;j++)vr[j+n*m]=V[i+n*j];
       }
      m++;
     }
   }

  if(m>k)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Too many tangent vectors (%d when k is only %d) - Jacobian is not full rank!",m,k);
    MFSetError(e,4,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    return 0;
   }

#ifdef MFALLOWVERBOSE
  if(verbose){printf("\n      Gram-Schmidt\n");fflush(stdout);}
#endif

  t=0.;
  for(j=0;j<n;j++)t+=vr[j]*vr[j];
  if(fabs(t)>1.e-15)
   {
    t=1./sqrt(t);
    for(j=0;j<n;j++){vr[j]=vr[j]*t;Phi[j+n*0]=vr[j];}

#ifdef MFALLOWVERBOSE
    if(verbose){printf("\n     Phi[%d]=(%lf",0,Phi[0+n*0]);for(j=1;j<n;j++)printf(",%lf",Phi[j+n*0]);printf(")\n");fflush(stdout);}
#endif

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
      for(j=0;j<n;j++){vr[j+n*i]=vr[j+n*i]*t;Phi[j+n*i]=vr[j+n*i];}

#ifdef MFALLOWVERBOSE
      if(verbose){printf("     Phi[%d]=(%lf",i,Phi[0+n*i]);for(j=1;j<n;j++)printf(",%lf",Phi[j+n*i]);printf(")\n");fflush(stdout);}
#endif

     }
   }

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("\n      test Gram-Schmidt, k=%d\n",k);
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
   }
#endif

#ifdef MFALLOWVERBOSE
  if(verbose) { printf("\n      done print Q^T Q\n"); fflush(stdout); }
#endif

  if(check)
   {
    for(j=0;j<k;j++)
     {
      t=0.;for(m=0;m<n;m++)t+=a[0+n*m]*Phi[m+n*j];
      printf("     A.Phi[%d]=(%lf",j,t);
      for(i=1;i<n;i++)
       {
        t=0.;for(m=0;m<n;m++)t+=a[i+n*m]*Phi[m+n*j];
        printf(",%lf",t);
       }
      printf(")\n");fflush(stdout);
     }
   }

#ifdef MFALLOWVERBOSE
  if(verbose) { printf("\n      done test Gram-Schmidt\n"); fflush(stdout); }
#endif

  free(a);
  free(A);
  free(U);
  free(V);
  free(s);
  free(vr);
  free(work);
  free(f);
  free(df);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("done %s\n",RoutineName);fflush(stdout);}
#endif

  return 1;
#else
  sprintf(MFAlgebraicMFErrorHandlerMsg,"MFAlgebraic requires LAPACK");
  MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
  return 0;
#endif
 }

double MFScaleAlgebraic(int n,int k,MFNVector vu,MFNKMatrix mPhi,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFScaleAlgebraic"};
  double r=1.;
  double *A;
  double *b;
  double *c;
  double *u,*Phi;

  char jobvl='N';                   /* No left  eigenvectors */
  char jobvr='V';                   /*    right eigenvectors */
  double *wr;    /* k */            /* Array to hold real part of eigenvalues */
  double *wi;    /* k */            /* Array to hold imag part of eigenvalues */
  int ldvl;                       /* Leading dimension of vl */
  double *vl;   /* k*k */            /* Array to hold left eigenvectors */
  int ldvr;                       /* Leading dimension of vr */
  double *vr;   /* k*k */           /* Array to hold right eigenvectors */
  double *work;    /* k*k */        /* Work array */
  int lwork;                     /* Length of work array */
  int info=0;

  double cn;
  double tu=0.;
  double t1=0.;
  double t2=0.;
  double t3=0.;
  double su=0.;
  double s1=0.;
  double s2=0.;
  double s3=0.;
  int i,j,l;
  int ii,jj;
  int *ipvt;
  char trans='N';
  int one=1;
  struct MFAlgebraicData *data;
  double *f=NULL;
  double *df=NULL;
  double *ddf=NULL;
  int m;
  int verbose=0;

#ifdef HAVE_LAPACK

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s, k=%d, n=%d\n",RoutineName,k,n);fflush(stdout);}
#endif

  Phi=MFNKM_CStar(mPhi,e);
  u=MFNV_CStar(vu,e);

  data=(struct MFAlgebraicData*)d;

  A=(double*)malloc(n*n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(A==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*n*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0.;
   }
#endif

  b=(double*)malloc(n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(b==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(A);
    return 0.;
   }
#endif

  ipvt=(int*)malloc(n*sizeof(int));

#ifndef MFNOSAFETYNET
  if(ipvt==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(int));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(A);
    free(b);
    return 0.;
   }
#endif

  f=(double*)malloc((n-k)*sizeof(double));

#ifndef MFNOSAFETYNET
  if(f==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",(n-k)*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(A);
    free(b);
    free(ipvt);
    return 0.;
   }
#endif

  df=(double*)malloc(n*(n-k)*sizeof(double));

#ifndef MFNOSAFETYNET
  if(df==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*(n-k)*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(A);
    free(b);
    free(ipvt);
    free(f);
    return 0.;
   }
#endif

  ddf=(double*)malloc(n*n*(n-k)*sizeof(double));

#ifndef MFNOSAFETYNET
  if(df==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*n*(n-k)*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(A);
    free(b);
    free(ipvt);
    free(f);
    free(df);
    return 0.;
   }
#endif

  c=(double*)malloc(k*k*sizeof(double));

#ifndef MFNOSAFETYNET
  if(c==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",k*k*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(A);
    free(b);
    free(ipvt);
    free(f);
    free(df);
    return 0.;
   }
#endif

  wr=(double*)malloc(k*sizeof(double));

#ifndef MFNOSAFETYNET
  if(wr==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",k*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(A);
    free(b);
    free(ipvt);
    free(f);
    free(df);
    free(c);
    return 0.;
   }
#endif

  wi=(double*)malloc(k*sizeof(double));

#ifndef MFNOSAFETYNET
  if(wi==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",k*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(A);
    free(b);
    free(ipvt);
    free(f);
    free(df);
    free(c);
    free(wr);
    return 0.;
   }
#endif

  ldvl=k;
  vl=(double*)malloc(k*k*sizeof(double));

#ifndef MFNOSAFETYNET
  if(vl==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",k*k*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(A);
    free(b);
    free(ipvt);
    free(f);
    free(df);
    free(c);
    free(wr);
    free(wi);
    return 0.;
   }
#endif

  ldvr=k;
  vr=(double*)malloc(k*k*sizeof(double));

#ifndef MFNOSAFETYNET
  if(vr==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",k*k*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(A);
    free(b);
    free(ipvt);
    free(f);
    free(df);
    free(c);
    free(wr);
    free(wi);
    free(vl);
    return 0.;
   }
#endif

  lwork=4*k;
  work=(double*)malloc(lwork*sizeof(double));

#ifndef MFNOSAFETYNET
  if(work==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",lwork*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(A);
    free(b);
    free(ipvt);
    free(f);
    free(df);
    free(c);
    free(wr);
    free(wi);
    free(vl);
    free(vr);
    return 0.;
   }
#endif

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s u:(%lf",RoutineName,u[0]);for(l=1;l<n;l++)printf(",%lf",u[l]);printf(")\n");fflush(stdout);}
#endif

  if(data->expr)
   {
    for(j=0;j<n;j++)
     {
      ECEvaluateFunction(data->dF[j],u,f);
      for(i=0;i<n-k;i++)A[i+n*j]=f[i];
     }
   }else{
    m=n-k;
    (data->dFSub)(&n,u,&m,df,data->data,e);
    for(i=0;i<m;i++)
     for(j=0;j<n;j++)A[i+n*j]=df[i+m*j];
   }
  for(j=0;j<n;j++)
    for(i=n-k;i<n;i++)
      A[i+n*j]=Phi[j+n*(i-n+k)];

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("A:\n");
    for(i=0;i<n;i++)
     {
      printf(" [%lf",A[i+n*0]);for(l=1;l<n;l++)printf(",%lf",A[i+n*l]);printf("]\n");fflush(stdout);
     }
    printf("\n");

    printf("    factor\n");fflush(stdout);
   }
#endif

  CALLDGETRF(&n,&n,A,&n,ipvt,&info);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("    done factor\n");fflush(stdout);}
#endif
  if(data->expr)
   {
    for(i=0;i<n;i++)
     {
      for(j=0;j<n;j++)
       {
  if(verbose){printf("    evaluate ddF expr\n");fflush(stdout);}
        ECEvaluateFunction(data->ddF[i+n*j],u,ddf+(n-k)*(i+n*j));
       }
     }
   }else{
    m=n-k;
  if(verbose){printf("    evaluate ddF sub\n");fflush(stdout);}
    (data->ddFSub)(&n,u,&m,ddf,data->data,e);
   }
  if(verbose){printf("    done evaluate ddF\n");fflush(stdout);}

  for(i=0;i<k;i++)
   {
    for(j=i;j<k;j++)
     {

#ifdef MFALLOWVERBOSE
      if(verbose){printf("    (%d,%d) ",i,j);fflush(stdout);}
#endif

      for(l=0;l<n;l++)b[l]=0.;

#ifdef MFALLOWVERBOSE
      if(verbose)
       {
        printf("\n Phi[%d]=(%lf",i,Phi[0+n*i]);for(l=1;l<n;l++)printf(",%lf",Phi[l+n*i]);printf(")\n");fflush(stdout);
        printf(" Phi[%d]=(%lf",j,Phi[0+n*j]);for(l=1;l<n;l++)printf(",%lf",Phi[l+n*j]);printf(")\n");fflush(stdout);
       }
#endif

      for(ii=0;ii<n;ii++)
       {
        for(jj=0;jj<n;jj++)
         {
          for(l=0;l<n-k;l++)b[l]-=ddf[l+(n-k)*(ii+n*jj)]*Phi[ii+n*i]*Phi[jj+n*j];

#ifdef MFALLOWVERBOSE
          if(verbose)
           {
            printf(" dF/dx_%d/dx_%d(%lf",ii,jj,u[0]);
            for(l=1;l<n;l++)printf(",%lf",u[l]);
            printf(") = (%lf",ddf[0+(n-k)*(ii+n*jj)]*Phi[ii+n*i]*Phi[jj+n*i]);
            for(l=1;l<n-k;l++)printf(",%lf",ddf[0+(n-k)*(ii+n*jj)]*Phi[ii+n*i]*Phi[jj+n*i]);
            printf(")\n");
            fflush(stdout);
           }
#endif
         }
       }

#ifdef MFALLOWVERBOSE
      if(verbose)
       {
        printf("      backsolve ");fflush(stdout);
        printf(" rhs=(%lf",b[0]);for(l=1;l<n;l++)printf(",%lf",b[l]);printf("), info=%d\n",info);fflush(stdout);
       }
#endif

      CALLDGETRS(&trans,&n,&one,A,&n,ipvt,b,&n,&info);

#ifdef MFALLOWVERBOSE
      if(verbose)
       {
        printf(" b=(%lf",b[0]);for(l=1;l<n;l++)printf(",%lf",b[l]);printf("), info=%d\n",info);fflush(stdout);
       }
#endif

      c[i+k*j]=0.;
      for(l=0;l<n;l++)c[i+k*j]+=b[l]*b[l];
      c[i+k*j]=sqrt(c[i+k*j]);
      c[j+k*i]=c[i+k*j];

#ifdef MFALLOWVERBOSE
      if(verbose){printf("c[%d,%d]=%lf\n",i,j,c[i+k*j]);fflush(stdout);}
#endif

     }
   }

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("c:\n");
    for(i=0;i<k;i++)
     {
      printf(" [%lf",c[i+k*0]);for(l=1;l<k;l++)printf(",%lf",c[i+k*l]);printf("]\n");fflush(stdout);
     }
    printf("\n");
   }
#endif

#ifdef HAVE_LAPACK

  info=0;
  jobvl='N';
  jobvr='N';

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("    find largest eigenvalue\n");fflush(stdout);
    printf("      call dgeev\n");fflush(stdout);
   }
#endif

  CALLDGEEV(&jobvl,&jobvr,&k,c,&k,wr,wi,vl,&ldvl,vr,&ldvr,work,&lwork,&info);

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("      info from dgeev is %d\n",info);fflush(stdout);
    printf("      find max\n");fflush(stdout);
   }
#endif

#else
  sprintf(MFAlgebraicMFErrorHandlerMsg,"Lapack is required for this routine");
  MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
  return 0;
#endif


  cn=sqrt(wr[0]*wr[0]+wi[0]*wi[0]);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("ScaleAlgebraic, ev= %lf",cn);}
#endif

  l=0;
  tu=sqrt(wr[l]*wr[l]+wi[l]*wi[l]);
  cn=tu;

#ifdef MFALLOWVERBOSE
    if(verbose){printf("evs=(%lf",tu);}
#endif

  for(l=1;l<k;l++)
   {
    tu=sqrt(wr[l]*wr[l]+wi[l]*wi[l]);
    if(tu>cn)cn=tu;

#ifdef MFALLOWVERBOSE
    if(verbose){printf(",%lf",tu);}
#endif

   }

#ifdef MFALLOWVERBOSE
  if(verbose){printf(")\n");}
#endif

  r=sqrt(2*MFEpsilon/cn);

#ifdef MFALLOWVERBOSE
  if(verbose){printf(", max ev=%lf, r=%lf\n",cn,r);fflush(stdout);}
#endif

  free(A);
  free(b);
  free(c);
  free(f);
  free(df);
  free(wr);
  free(wi);
  free(vl);
  free(vr);
  free(work);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("done %s\n",RoutineName);fflush(stdout);}
#endif

  return r;
#else
  sprintf(MFAlgebraicMFErrorHandlerMsg,"MFAlgebraic requires LAPACK");
  MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
  return 0;
#endif
 }

int MFStopAlgebraic(MFImplicitMF thisMF,MFNVector u0,MFNKMatrix Phi0,MFNVector u1,MFNKMatrix Phi1,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFStopAlgebraic"};
  int verbose=0;

  if(verbose)printf("\n\n***%s\n",RoutineName);fflush(stdout);

  MFAlgebraicSetStability(thisMF,u1,Phi0,d,e);
  MFAlgebraicSetStability(thisMF,u0,Phi0,d,e);

#ifdef MFALLOWVERBOSE
  if(verbose)printf("%s index 0=%d, index 1=%d\n",RoutineName,MFNVGetIndex(u0,e),MFNVGetIndex(u1,e));fflush(stdout);
#endif

  if(MFNVGetIndex(u0,e)!=MFNVGetIndex(u1,e))return 1;
   else return 0;
 }

int MFFullNRealEV(int n,double *A, MFErrorHandler e)
 {
  static char RoutineName[]={"MFFullNRealEV"};
  int i,j;
  int lda;
  char jobvl='N';                   /* No left  eigenvectors */
  char jobvr='V';                   /*    right eigenvectors */
  static double *wr=NULL;    /* k */            /* Array to hold real part of eigenvalues */
  static double *wi=NULL;    /* k */            /* Array to hold imag part of eigenvalues */
  int ldvl;                       /* Leading dimension of vl */
  static double *vl=NULL;   /* k*k */            /* Array to hold left eigenvectors */
  int ldvr;                       /* Leading dimension of vr */
  static double *vr=NULL;   /* k*k */           /* Array to hold right eigenvectors */
  static double *work=NULL;    /* k*k */        /* Work array */
  int lwork;                     /* Length of work array */
  int info=0;
  int verbose=0;
  int result;

#ifdef HAVE_LAPACK

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    for(i=0;i<n;i++)
     {
      printf("[");
      for(j=0;j<n;j++)
       {
        if(j>0)printf(" ");
        printf("%lf",A[i+n*j]);
       }
      printf("]\n");
     }
   }
#endif

  lda=n;
  wr=(double*)realloc((void*)wr,n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(wr==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0;
   }
#endif

  wi=(double*)realloc((void*)wi,n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(wi==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(wr);
    return 0;
   }
#endif

  ldvl=n;
  vl=(double*)realloc((void*)vl,n*n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(vl==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*n*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(wr);
    free(wi);
    return 0;
   }
#endif

  ldvr=n;
  vr=(double*)realloc((void*)vr,n*n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(vr==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*n*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(wr);
    free(wi);
    free(vl);
    return 0;
   }
#endif

  lwork=4*n;
  work=(double*)realloc((void*)work,lwork*sizeof(double));

#ifndef MFNOSAFETYNET
  if(work==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",lwork*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(wr);
    free(wi);
    free(vl);
    free(vr);
    return 0;
   }
#endif

#ifdef HAVE_LAPACK


  info=0;
  jobvl='N';
  jobvr='V';

#ifdef MFALLOWVERBOSE
  if(verbose){printf("      call dgeev\n");fflush(stdout);}
#endif

  CALLDGEEV(&jobvl,&jobvr,&n,A,&lda,wr,wi,vl,&ldvl,vr,&ldvr,work,&lwork,&info);

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("      done dgeev\n");fflush(stdout);
    printf("%s\n",RoutineName);
   }
#endif

#else
  sprintf(MFAlgebraicMFErrorHandlerMsg,"Lapack is required for this routine");
  MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
  return 0;
#endif


  result=0;
  for(i=0;i<n;i++)
   {
    if(verbose){printf("   %d (%le,%le)\n",i,wr[i],wi[i]);fflush(stdout);}
    if(wr[i]>0.)result++;
   }

#ifdef MFALLOWVERBOSE
  if(verbose){printf("  number with positive real part =%d\n\n\n",result);fflush(stdout);}
#endif

  return result;
#else
  sprintf(MFAlgebraicMFErrorHandlerMsg,"MFAlgebraic requires LAPACK");
  MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
  return 0;
#endif
 }

void MFAlgebraicSetStability(MFImplicitMF M,MFNVector u, MFNKMatrix mPhi,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAlgebraicSetStability"};
  static double *A=NULL; /* [n*n];*/
  struct MFAlgebraicData *data;
  static double *f=NULL; /* [k];  */
  static double *df=NULL; /* [k];  */
  int n,k;
  int i,j,l,m;
  int result;
  int verbose=0;
  double *Phi;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s\n",RoutineName);fflush(stdout);}
#endif

  n=MFIMF_N(M,e);
  k=MFIMF_K(M,e);
  data=(struct MFAlgebraicData*)MFIMFGetData(M,e);

  A=(double*)realloc((void*)A,n*n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(A==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*n*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  f=(double*)realloc((void*)f,(n-k)*sizeof(double));

#ifndef MFNOSAFETYNET
  if(f==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",(n-k)*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(A);
    return;
   }
#endif

  df=(double*)realloc((void*)df,(n-k)*n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(df==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",(n-k)*n*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(A);
    free(f);
    return;
   }
#endif

  Phi=MFNKM_CStar(mPhi,e);

  if(data->expr)
   {
    for(j=0;j<n;j++)
     {
      ECEvaluateFunction(data->dF[j],MFNV_CStar(u,e),f);
      for(i=0;i<n-k;i++)A[i+n*j]=f[i];
     }
   }else{
    m=n-k;
    (data->dFSub)(&n,MFNV_CStar(u,e),&m,df,data->data,e);
    for(i=0;i<n-k;i++)
     for(j=0;j<n;j++)A[i+n*j]=df[i+(n-k)*j];
   }
  for(j=0;j<n;j++)
    for(i=0;i<k;i++)A[n-k+i+n*j]=Phi[j+n*i];

  result=MFFullNPosEV(n,k,A,e)%2;
  MFNVSetIndex(u,result,e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s: index is %d\n",RoutineName,result);fflush(stdout);}
#endif

  return;
 }

int MFAlgebraicProjectToSave(MFNVector u, double *x, void *d, MFErrorHandler e)
 {
  int i;
  struct MFAlgebraicData *data;

  data=(struct MFAlgebraicData*)d;

  if(x==NULL)return data->n;
  for(i=0;i<data->n;i++)x[i]=MFNV_C(u,i,e);

  return 0;
 }

int MFAlgebraicProjectToDraw(MFNVector u, double *x, void *d, MFErrorHandler e)
 {
  int i;
  struct MFAlgebraicData *data;

  data=(struct MFAlgebraicData*)d;

  if(x==NULL)return data->n;
  for(i=0;i<data->n;i++)x[i]=MFNV_C(u,i,e);

  return 0;
 }

int MFAlgebraicProjectForBB(MFNVector u, double *x, void *d, MFErrorHandler e)
 {
  int i;
  struct MFAlgebraicData *data;

  data=(struct MFAlgebraicData*)d;

  if(x==NULL)return data->n;
  for(i=0;i<data->n;i++)x[i]=MFNV_C(u,i,e);

  return 0;
 }

int MFSingularAlgebraic(int n,int k,MFNVector vu,MFNKMatrix mPhi,MFNVector perp,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFSingularAlgebraic"};
  double *A;             /* n*n */
  char jobvl='N';                   /* No left  eigenvectors */
  char jobvr='V';                   /*    right eigenvectors */
  int lda;                        /* Leading dimension of A */
  double *wr;    /* n */            /* Array to hold real part of eigenvalues */
  double *wi;    /* n */            /* Array to hold imag part of eigenvalues */
  int ldvl;                       /* Leading dimension of vl */
  double *vl;   /* n*n */            /* Array to hold left eigenvectors */
  int ldvr=4;                       /* Leading dimension of vr */
  double *vr;   /* n*n */           /* Array to hold right eigenvectors */
  double *work;    /* n*n */        /* Work array */
  int lwork;                     /* Length of work array */
  int info=0;
  int i,j,l,m,mm;
  double t;
  struct MFAlgebraicData *data;
  double *f=NULL; /* [k];  */
  double *u,*Phi,*phi;
  int verbose;
  double *s;    /* n */
  double *U;    /* n*n */
  int ldU;
  double *V;    /* n*n */
  int ldV;
  double *df=NULL;

#ifdef HAVE_LAPACK

/* Finds a NULL vector at u, perpendicular to the tangent plane Phi (assuming the Jacobian is singular) */

  verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s, k=%d, n=%d\n",RoutineName,k,n);fflush(stdout);}
#endif

  Phi=MFNKM_CStar(mPhi,e);
  u=MFNV_CStar(vu,e);
  phi=MFNV_CStar(perp,e);

  data=(struct MFAlgebraicData*)d;

  A=(double*)malloc(n*n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(A==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*n*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0;
   }
#endif

  lda=n;
  ldU=n;

  U=(double*)malloc(n*n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(U==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*n*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(A);
    return 0;
   }
#endif

  ldV=n;
  V=(double*)malloc(n*n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(V==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*n*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(A);
    free(U);
    return 0;
   }
#endif

  ldvr=n;
  vr=(double*)malloc(n*n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(vr==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*n*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(A);
    free(U);
    free(V);
    return 0;
   }
#endif

  f=(double*)malloc((n-k)*sizeof(double));

#ifndef MFNOSAFETYNET
  if(f==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",(n-k)*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(A);
    free(U);
    free(V);
    free(vr);
    return 0;
   }
#endif

  df=(double*)malloc(n*(n-k)*sizeof(double));

#ifndef MFNOSAFETYNET
  if(df==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*(n-k)*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(A);
    free(U);
    free(V);
    free(vr);
    free(f);
    return 0;
   }
#endif

  s=(double*)malloc(n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(s==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(A);
    free(U);
    free(V);
    free(vr);
    free(f);
    free(df);
    return 0;
   }
#endif

  if(data->expr)
   {

#ifdef MFALLOWVERBOSE
    if(verbose){printf("%s, %s\n",data->variables,data->expression);fflush(stdout);}
#endif

    for(j=0;j<n;j++)
     {
      ECEvaluateFunction(data->dF[j],u,f);
      for(i=0;i<n-k;i++)A[i+n*j]=f[i];
     }
   }else{

#ifdef MFALLOWVERBOSE
    if(verbose){printf("Functions\n");fflush(stdout);}
#endif

    m=n-k;
    (data->dFSub)(&n,u,&m,df,data->data,e);
    for(i=0;i<n-k;i++)
     for(j=0;j<n;j++)A[i+n*j]=df[i+(n-k)*j];
   }
  for(j=0;j<n;j++)
    for(i=0;i<k;i++)A[n-k+i+n*j]=Phi[j+n*i];

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("u:(%lf",u[0]);for(l=1;l<n;l++)printf(",%lf",u[l]);printf(")\n");fflush(stdout);
    printf("Phi:\n");MFPrintNKMatrix(stdout,mPhi,e);
    printf("A:\n");
    for(i=0;i<n;i++)
     {
      printf(" [%lf",A[i+n*0]);for(l=1;l<n;l++)printf(",%lf",A[i+n*l]);printf("]\n");fflush(stdout);
     }
    printf("\n");
   }
#endif

#ifdef HAVE_LAPACK

  i=-1;
  info=0;
  jobvl='N';
  jobvr='A';
  CALLDGESVD(&jobvl,&jobvr,&n,&n,A,&lda,s,U,&ldU,V,&ldV,&t,&i,&info);
  lwork=MAX(4*n,5*n-4);
  lwork=round(t);

  work=(double*)malloc(lwork*sizeof(double));

#ifndef MFNOSAFETYNET
  if(work==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",lwork*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0;
   }
#endif

  info=0;
  jobvl='N';
  jobvr='A';
  CALLDGESVD(&jobvl,&jobvr,&n,&n,A,&lda,s,U,&ldU,V,&ldV,work,&lwork,&info);

#else
  sprintf(MFAlgebraicMFErrorHandlerMsg,"Lapack is required for this routine");
  MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
  return 0;
#endif


#ifdef MFALLOWVERBOSE
  if(verbose){printf("\n      Null Vectors\n");fflush(stdout);}
#endif

  m=0;
  for(i=0;i<n;i++)
   {

#ifdef MFALLOWVERBOSE
    if(verbose){printf("  s[%d]=%lf,   v[%d]=(%lf",i,s[i],i,V[i+n*0]);for(j=1;j<n;j++)printf(",%lf",V[i+n*j]);printf(")\n");fflush(stdout);}
#endif

    if(fabs(s[i])<1.e-3)
     {
      for(j=0;j<n;j++)vr[j+n*m]=V[i+n*j];
      m++;
     }
   }

  if(m<1)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Jacobian is full rank!");
    MFSetError(e,4,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    free(A);
    free(U);
    free(V);
    free(s);
    free(work);
    free(f);
    free(df);
    return 0;
   }

  if(m>1)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Jacobian is rank deficient by more than one! m=%d",m);
    MFSetError(e,4,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    free(A);
    free(U);
    free(V);
    free(s);
    free(work);
    free(f);
    free(df);
    return 0;
   }

  t=0.;
  for(j=0;j<n;j++)t+=vr[j]*vr[j];
  if(fabs(t)>1.e-15)
   {
    t=1./sqrt(t);
    for(j=0;j<n;j++){vr[j]=vr[j]*t;phi[j]=vr[j];}

#ifdef MFALLOWVERBOSE
    if(verbose){printf("     phi=(%lf",phi[0]);for(j=1;j<n;j++)printf(",%lf",phi[j]);printf(")\n");fflush(stdout);}
#endif

   }

  free(A);
  free(U);
  free(V);
  free(vr);
  free(s);
  free(work);
  free(f);
  free(df);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("done %s\n",RoutineName);fflush(stdout);}
#endif

  return 1;
#else
  sprintf(MFAlgebraicMFErrorHandlerMsg,"MFAlgebraic requires LAPACK");
  MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
  return 0;
#endif
 }

int MFFullNPosSV(int n,int k, double *A, MFErrorHandler e)
 {
  static char RoutineName[]={"MFFullNPosSV"};
  int i,j;
  int lda;
  char jobvl='N';
  char jobvr='N';
  double *s=NULL;
  double *U=NULL;
  double *V=NULL;
  double *work=NULL;
  double t;
  int lwork;
  int info=0;
  int verbose=0;
  int result;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s\n",RoutineName);}
#endif

  if(0&&verbose)
   {
    for(i=0;i<n;i++)
     {
      printf("[");
      for(j=0;j<n;j++)
       {
        if(j>0)printf(" ");
        printf("%lf",A[i+n*j]);
       }
      printf("]\n");
     }
   }

  s=(double*)malloc(n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(s==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0;
   }
#endif

  U=(double*)malloc(n*n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(U==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*n*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(s);
    return 0;
   }
#endif

  V=(double*)malloc(n*n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(V==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*n*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(s);
    free(U);
    return 0;
   }
#endif

#ifdef HAVE_LAPACK

  i=-1;
  info=0;
  jobvl='N';
  jobvr='A';
  CALLDGESVD(&jobvl,&jobvr,&n,&n,A,&n,s,U,&n,V,&n,&t,&i,&info);
  lwork=MAX(4*n,5*n-4);
  lwork=round(t);
  work=(double*)malloc(lwork*sizeof(double));

#ifndef MFNOSAFETYNET
  if(work==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",lwork*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0;
   }
#endif

  info=0;
  jobvl='N';
  jobvr='N';
  CALLDGESVD(&jobvl,&jobvr,&n,&n,A,&n,s,U,&n,V,&n,work,&lwork,&info);

#else
  sprintf(MFAlgebraicMFErrorHandlerMsg,"Lapack is required for this routine");
  MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
  return 0;
#endif


  result=0;
  for(i=0;i<n;i++)
   {

#ifdef MFALLOWVERBOSE
    if(verbose)
     {
      printf("   %d %le\n",i,s[i]);
      fflush(stdout);
     }
#endif

    if(i<n-k&&s[i]>0.)result++;
   }
  if(verbose){printf("  number positive = %d\n",result);fflush(stdout);}
  free(s);
  free(U);
  free(V);
  free(work);

  return result;
 }

double MFFullProdEV(int n,int k,double *A, MFErrorHandler e)
 {
  static char RoutineName[]={"MFFullProdEV"};
  int i,j,l;
  int lda;
  char jobvl='N';                   /* No left  eigenvectors */
  char jobvr='V';                   /*    right eigenvectors */
  static double *wr=NULL;    /* k */            /* Array to hold real part of eigenvalues */
  static double *wi=NULL;    /* k */            /* Array to hold imag part of eigenvalues */
  int ldvl;                       /* Leading dimension of vl */
  static double *vl=NULL;   /* k*k */            /* Array to hold left eigenvectors */
  int ldvr;                       /* Leading dimension of vr */
  static double *vr=NULL;   /* k*k */           /* Array to hold right eigenvectors */
  static double *work=NULL;    /* k*k */        /* Work array */
  int lwork;                     /* Length of work array */
  int info=0;
  int verbose=0;
  int result;
  double resultI,resultR;
  double I,R;
  double prod;

#ifdef HAVE_LAPACK

  if(verbose)
   {
    printf("%s\n",RoutineName);
    for(i=0;i<n;i++)
     {
      printf("[");
      for(j=0;j<n;j++)
       {
        if(j>0)printf(" ");
        printf("%lf",A[i+n*j]);
       }
      printf("]\n");
     }

/*  i=0;
    printf("So x     =(%lf",A[i+n*0]/2);
    for(j=1;j<n;j++)printf(",%lf",A[i+n*j]/2);
    printf(")\n");
    i++;
    printf("   phi[1]=(%lf",A[i+n*0]);
    for(j=1;j<n;j++)printf(",%lf",A[i+n*j]);
    printf(")\n");
    i++;
    printf("   phi[0]=(%lf",A[i+n*0]);
    for(j=1;j<n;j++)printf(",%lf",A[i+n*j]);
    printf(")\n");
    printf("   2x.phi[0]=%lf\n",A[0+n*0]*A[1+n*0]+A[0+n*1]*A[1+n*1]+A[0+n*2]*A[1+n*2]);
    printf("   2x.phi[1]=%lf\n",A[0+n*0]*A[2+n*0]+A[0+n*1]*A[2+n*1]+A[0+n*2]*A[2+n*2]);
    printf("   4x.x     =%lf\n",A[0+n*0]*A[0+n*0]+A[0+n*1]*A[0+n*1]+A[0+n*2]*A[0+n*2]);
    printf("   phi[0]^2 =%lf\n",A[1+n*0]*A[1+n*0]+A[1+n*1]*A[1+n*1]+A[1+n*2]*A[1+n*2]);
    printf("   phi[1]^2 =%lf\n",A[2+n*0]*A[2+n*0]+A[2+n*1]*A[2+n*1]+A[2+n*2]*A[2+n*2]);
    printf("   phi[0].phi[1]=%lf\n",A[1+n*0]*A[2+n*0]+A[1+n*1]*A[2+n*1]+A[1+n*2]*A[2+n*2]);*/
   }

  lda=n;
  wr=(double*)realloc((void*)wr,n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(wi==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0.;
   }
#endif

  wi=(double*)realloc((void*)wi,n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(wi==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(wr);
    return 0.;
   }
#endif

  ldvl=n;
  vl=(double*)realloc((void*)vl,n*n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(vl==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*n*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(wr);
    free(wi);
    return 0.;
   }
#endif

  ldvr=n;
  vr=(double*)realloc((void*)vr,n*n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(vr==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*n*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(wr);
    free(wi);
    free(vl);
    return 0.;
   }
#endif

  lwork=4*n;
  work=(double*)realloc((void*)work,lwork*sizeof(double));

#ifndef MFNOSAFETYNET
  if(work==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",lwork*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(wr);
    free(wi);
    free(vl);
    free(vr);
    return 0.;
   }
#endif

#ifdef HAVE_LAPACK

  info=0;
  jobvl='N';
  jobvr='V';

#ifdef MFALLOWVERBOSE
  if(verbose){printf("      call dgeev\n");fflush(stdout);}
#endif

  CALLDGEEV(&jobvl,&jobvr,&n,A,&lda,wr,wi,vl,&ldvl,vr,&ldvr,work,&lwork,&info);

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("      done dgeev info=%d\n",info);
    printf("%s\n",RoutineName);
    fflush(stdout);
   }
#endif

#else
  sprintf(MFAlgebraicMFErrorHandlerMsg,"Lapack is required for this routine");
  MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
  return 0;
#endif


  resultR=1.;
  resultI=0.;
  for(i=0;i<n;i++)
   {

#ifdef MFALLOWVERBOSE
    if(verbose){printf("   %d (%le,%le)\n",i,wr[i],wi[i]);fflush(stdout);}
#endif

    R=resultR*wr[i]-resultI*wi[i];
    I=resultR*wi[i]+resultI*wr[i];
    resultR=R;
    resultI=I;
   }

#ifdef MFALLOWVERBOSE
  if(verbose){printf("  product of eigenvalues =(%lf,%lf)\n",resultR,resultI);fflush(stdout);}
#endif

  return resultR;

#else
    sprintf(MFAlgebraicMFErrorHandlerMsg,"The routine requires dgeev from LAPACK");
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    return 0.;
#endif
 }

int MFFullNPosEV(int n,int k,double *A, MFErrorHandler e)
 {
  static char RoutineName[]={"MFFullNPosEV"};
  int i,j,l;
  int lda;
  char jobvl='N';                   /* No left  eigenvectors */
  char jobvr='V';                   /*    right eigenvectors */
  static double *wr=NULL;    /* k */            /* Array to hold real part of eigenvalues */
  static double *wi=NULL;    /* k */            /* Array to hold imag part of eigenvalues */
  int ldvl;                       /* Leading dimension of vl */
  static double *vl=NULL;   /* k*k */            /* Array to hold left eigenvectors */
  int ldvr;                       /* Leading dimension of vr */
  static double *vr=NULL;   /* k*k */           /* Array to hold right eigenvectors */
  static double *work=NULL;    /* k*k */        /* Work array */
  int lwork;                     /* Length of work array */
  int info=0;
  int verbose=0;
  int result;

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("%s\n",RoutineName);
    for(i=0;i<n;i++)
     {
      printf("[");
      for(j=0;j<n;j++)
       {
        if(j>0)printf(" ");
        printf("%lf",A[i+n*j]);
       }
      printf("]\n");
     }
   }
#endif

  lda=n;
  wr=(double*)realloc((void*)wr,n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(wr==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0;
   }
#endif

  wi=(double*)realloc((void*)wi,n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(wi==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(wr);
    return 0;
   }
#endif

  ldvl=n;
  vl=(double*)realloc((void*)vl,n*n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(vl==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*n*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(wr);
    free(wi);
    return 0;
   }
#endif

  ldvr=n;
  vr=(double*)realloc((void*)vr,n*n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(vr==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*n*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(wr);
    free(wi);
    free(vl);
    return 0;
   }
#endif

  lwork=4*n;
  work=(double*)realloc((void*)work,lwork*sizeof(double));

#ifndef MFNOSAFETYNET
  if(work==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",lwork*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(wr);
    free(wi);
    free(vl);
    free(vr);
    return 0;
   }
#endif

#ifdef HAVE_LAPACK

  info=0;
  jobvl='N';
  jobvr='V';
  CALLDGEEV(&jobvl,&jobvr,&n,A,&lda,wr,wi,vl,&ldvl,vr,&ldvr,work,&lwork,&info);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s, info back from dgeev is %d\n",RoutineName,info);}
#endif

#else
  sprintf(MFAlgebraicMFErrorHandlerMsg,"Lapack is required for this routine");
  MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
  return 0;
#endif

  result=0;
  for(i=0;i<n;i++)
   {
    int j;

#ifdef MFALLOWVERBOSE
    if(verbose){printf("   %d (%le,%le)",i,wr[i],wi[i]);
                printf("          (%le",vr[0+n*i]);for(j=1;j<n;j++)printf(",%le",vr[j+n*i]);printf(")\n");fflush(stdout);}
#endif

    if(fabs(wi[i])>1.e-7)
     {
      for(j=i+1;j<n;j++)
        if(fabs(wi[j]+wi[i])<2.e-7)result+=2;
     }else if(wr[i]>0.)result++;
   }

#ifdef MFALLOWVERBOSE
  if(verbose){printf("  number of eigenvalues with positive real part = %d\n",result);fflush(stdout);}
#endif

  return result;
 }

int MFFullNUnitEV(int n,int k,double *A, MFErrorHandler e)
 {
  static char RoutineName[]={"MFFullNUnitEV"};
  int i,j,l;
  int lda;
  char jobvl='N';                   /* No left  eigenvectors */
  char jobvr='V';                   /*    right eigenvectors */
  static double *wr=NULL;    /* k */            /* Array to hold real part of eigenvalues */
  static double *wi=NULL;    /* k */            /* Array to hold imag part of eigenvalues */
  int ldvl;                       /* Leading dimension of vl */
  static double *vl=NULL;   /* k*k */            /* Array to hold left eigenvectors */
  int ldvr;                       /* Leading dimension of vr */
  static double *vr=NULL;   /* k*k */           /* Array to hold right eigenvectors */
  static double *work=NULL;    /* k*k */        /* Work array */
  int lwork;                     /* Length of work array */
  int info=0;
  int verbose=0;
  int result;

  lda=n;
  wr=(double*)realloc((void*)wr,n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(wr==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0;
   }
#endif

  wi=(double*)realloc((void*)wi,n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(wi==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(wr);
    return 0;
   }
#endif

  ldvl=n;
  vl=(double*)realloc((void*)vl,n*n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(vl==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*n*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(wr);
    free(wi);
    return 0;
   }
#endif

  ldvr=n;
  vr=(double*)realloc((void*)vr,n*n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(vr==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*n*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(wr);
    free(wi);
    free(vl);
    return 0;
   }
#endif

  lwork=4*n;
  work=(double*)realloc((void*)work,lwork*sizeof(double));

#ifndef MFNOSAFETYNET
  if(work==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",lwork*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(wr);
    free(wi);
    free(vl);
    free(vr);
    return 0;
   }
#endif

#ifdef HAVE_LAPACK

  info=0;
  jobvl='N';
  jobvr='V';
  CALLDGEEV(&jobvl,&jobvr,&n,A,&lda,wr,wi,vl,&ldvl,vr,&ldvr,work,&lwork,&info);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s\n",RoutineName);}
#endif

#else
  sprintf(MFAlgebraicMFErrorHandlerMsg,"Lapack is required for this routine");
  MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
  return 0;
#endif

  result=0;
  for(i=0;i<n;i++)
   {

#ifdef MFALLOWVERBOSE
    if(verbose){printf("   %d (%le,%le)\n",i,wr[i],wi[i]);fflush(stdout);}
#endif

    if(wr[i]*wr[i]+wi[i]*wi[i]>1.)result++;
   }

#ifdef MFALLOWVERBOSE
  if(verbose){printf("  number of eigenvalues with modulus greater than one = %d\n",result);fflush(stdout);}
#endif

  return result;
 }

int MFFullNZeroSV(int n,int k, double *A, MFErrorHandler e)
 {
  static char RoutineName[]={"MFFullNZeroSV"};
  int i,j;
  int lda;
  char jobvl='N';
  char jobvr='N';
  double *s=NULL;
  double *U=NULL;
  double *V=NULL;
  double *work=NULL;
  double t;
  int lwork;
  int info=0;
  int verbose=0;
  int result;

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("%s\n",RoutineName);
    if(0)
     {
      for(i=0;i<n;i++)
       {
        printf("[");
        for(j=0;j<n;j++)
         {
          if(j>0)printf(" ");
          printf("%lf",A[i+n*j]);
         }
        printf("]\n");
       }
     }
   }
#endif

  s=(double*)malloc(n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(s==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0;
   }
#endif

  U=(double*)malloc(n*n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(U==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*n*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(s);
    return 0;
   }
#endif

  V=(double*)malloc(n*n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(V==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*n*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(s);
    free(U);
    return 0;
   }
#endif

#ifdef HAVE_LAPACK

  i=-1;
  info=0;
  jobvl='N';
  jobvr='A';
  CALLDGESVD(&jobvl,&jobvr,&n,&n,A,&n,s,U,&n,V,&n,&t,&i,&info);
  lwork=MAX(4*n,5*n-4);
  lwork=round(t);
  work=(double*)malloc(lwork*sizeof(double));

#ifndef MFNOSAFETYNET
  if(work==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",lwork*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(s);
    free(U);
    free(V);
    return 0;
   }
#endif


  info=0;
  jobvl='N';
  jobvr='N';
  CALLDGESVD(&jobvl,&jobvr,&n,&n,A,&n,s,U,&n,V,&n,work,&lwork,&info);

  result=0;
  for(i=0;i<n;i++)
   {

#ifdef MFALLOWVERBOSE
    if(verbose){printf("   %d %le\n",i,s[i]);fflush(stdout);}
#endif

    if(fabs(s[i])<1.e-10)result++;
   }

#else
  sprintf(MFAlgebraicMFErrorHandlerMsg,"Lapack is required for this routine");
  MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
  return 0;
#endif

#ifdef MFALLOWVERBOSE
  if(verbose){printf("  number zero = %d\n",result);fflush(stdout);}
#endif

  free(s);
  free(U);
  free(V);
  free(work);

  return result;
 }

/*! @} */

/*! @} */
int MFTangentAlgebraicWithGuess(int n,int k,MFNVector vu,MFNKMatrix mPhi0,MFNKMatrix mPhi,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFTangentAlgebraicWithGuess"};
  double *A=NULL; /* [n*n];*/
  double *a=NULL; /* [n*n];*/
  double *b=NULL; /* [n];  */
  int i,j,l,m;
  int itan;
  struct MFAlgebraicData *data;
  double *f=NULL;
  double *df=NULL;
  double *Phi0,*Phi,*u;
  int rc;
  int verbose=0;

  Phi=MFNKM_CStar(mPhi,e);
  Phi0=MFNKM_CStar(mPhi0,e);
  u=MFNV_CStar(vu,e);

  data=(struct MFAlgebraicData*)d;

  A=(double*)malloc(n*n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(A==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*n*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0;
   }
#endif

  a=(double*)malloc(n*n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(a==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*n*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(A);
    return 0;
   }
#endif

  b=(double*)malloc(n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(b==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(A);
    free(a);
    return 0;
   }
#endif

  f=(double*)malloc((n-k)*sizeof(double));

#ifndef MFNOSAFETYNET
  if(f==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",(n-k)*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(A);
    free(a);
    free(b);
    return 0;
   }
#endif

  df=(double*)malloc(n*(n-k)*sizeof(double));

#ifndef MFNOSAFETYNET
  if(df==NULL)
   {
    sprintf(MFAlgebraicMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",n*(n-k)*sizeof(double));
    MFSetError(e,12,RoutineName,MFAlgebraicMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(A);
    free(a);
    free(b);
    free(f);
    return 0;
   }
#endif
  if(data->expr)
   {
    for(j=0;j<n;j++)
     {
      ECEvaluateFunction(data->dF[j],u,f);
      for(i=0;i<n-k;i++)A[i+n*j]=f[i];
     }
    ECEvaluateFunction(data->F,u,f);
   }else{
    m=n-k;
    (data->FSub)(&n,u,&m,f,data->data,e);
    (data->dFSub)(&n,u,&m,df,data->data,e);
    for(i=0;i<m;i++)
     for(j=0;j<n;j++)A[i+n*j]=df[i+m*j];
   }

#ifdef MFALLOWVERBOSE
    if(verbose)
     {
      printf("guess at Phi:\n");
      for(j=0;j<n;j++)
       {
        printf(" [%lf",Phi[j+n*0]);for(l=1;l<k;l++)printf(",%lf",Phi0[j+n*l]);printf("]\n");fflush(stdout);
       }
     }
#endif

  for(itan=0;itan<k;itan++)
   {
    for(j=0;j<n;j++)b[j]=0.;
    b[itan+n-k]=1.;

    for(j=0;j<n;j++)
     {
      for(i=0;i<n-k;i++)a[i+n*j]=A[i+n*j];
      for(i=n-k;i<n-k+itan;i++)
        a[i+n*j]=Phi [j+n*(i-n+k)];
      for(i=n-k+itan;i<n;i++)
        a[i+n*j]=Phi0[j+n*(i-n+k)];
     }

#ifdef MFALLOWVERBOSE
    if(verbose)
     {
      printf("\n");
      printf("tangent #%d\n",itan);
      printf("%s A:\n",RoutineName);
      for(i=0;i<n;i++)
       {
        printf(" [%lf",a[i+n*0]);for(l=1;l<n;l++)printf(",%lf",a[i+n*l]);printf("]  [%lf]\n",b[i]);fflush(stdout);
       }
      printf("\n");
     }
#endif

    rc=MFSolveFull(n,a,b,e);

#ifdef MFALLOWVERBOSE
    if(verbose)
     {
      for(i=0;i<n;i++)printf(" [%lf]\n",b[i]);
     }

    for(i=0;i<n;i++)Phi[i+n*itan]=b[i];
#endif
   }

  free(A);
  free(a);
  free(b);
  free(f);
  free(df);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("done %s, rc=1\n",RoutineName);fflush(stdout);}
#endif
  return 1;
 }

#ifdef __cplusplus
}
#endif
