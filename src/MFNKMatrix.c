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

static char *id="@(#) $Id: MFNKMatrix.c,v 1.9 2011/07/21 17:42:46 mhender Exp $";

static char MFNKMatrixErrorMsg[256]="";

#include <MFErrorHandler.h>
#include <MFNKMatrix.h>
#include <MFNVector.h>
#include <MFKVector.h>
#include <float.h>
#include <stdlib.h>
#include <string.h>
#include <MFPrint.h>
#include <math.h>

#define DENSE 0
#define VECTOR 1

#ifdef __cplusplus
 extern "C" {
#endif

double MFAtlasDet(int,double*,MFErrorHandler);

struct MFNKMatrixSt
 {
  int type;
  int n;
  int k;
  double *data;
  MFNVector *cols;
  int nRefs;
 };

void MFFreeNKMatrix(MFNKMatrix thisMatrix,MFErrorHandler e)
 {
  static char RoutineName[]={"MFFreeNKMatrix"};
  int i;

#ifdef MFNOCONFIDENCE
  if(thisMatrix==NULL)
   {
    sprintf(MFNKMatrixErrorMsg,"Pointer to Matrix (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFNKMatrixErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  thisMatrix->nRefs--;

  if(thisMatrix->nRefs<1)
   {
    if(thisMatrix->data!=NULL)free(thisMatrix->data);
    if(thisMatrix->cols!=NULL)
     {
      for(i=0;i<thisMatrix->k;i++)MFFreeNVector(thisMatrix->cols[i],e);
      free(thisMatrix->cols);
     }
    free(thisMatrix);
   }
  return;
 }

void MFRefMatrix(MFNKMatrix thisMatrix,MFErrorHandler e)
 {
  static char RoutineName[]={"MFRefMatrix"};

#ifdef MFNOCONFIDENCE
  if(thisMatrix==NULL)
   {
    sprintf(MFNKMatrixErrorMsg,"Pointer to Matrix (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFNKMatrixErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  thisMatrix->nRefs++;
 }

MFNKMatrix MFCreateNKMatrix(int k,MFNVector *col,MFErrorHandler e)
 {
  static char RoutineName[]={"MFCreateNKMatrix"};
  MFNKMatrix thisMatrix;
  int i,j;
  double *u;

#ifdef MFNOCONFIDENCE
  if(k<1)
   {
    sprintf(MFNKMatrixErrorMsg,"Number of Columns %d (argument 1) is Illegal. Must be positive.",k);
    MFSetError(e,12,RoutineName,MFNKMatrixErrorMsg,__LINE__,__FILE__);
    return NULL;
   }

  if(col==NULL)
   {
    sprintf(MFNKMatrixErrorMsg,"List of Columns (argument 2) is NULL.");
    MFSetError(e,12,RoutineName,MFNKMatrixErrorMsg,__LINE__,__FILE__);
    return NULL;
   }
#endif

  thisMatrix=(MFNKMatrix)malloc(sizeof(struct MFNKMatrixSt));

#ifndef MFNOSAFETYNET
  if(thisMatrix==NULL)
   {
    sprintf(MFNKMatrixErrorMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct MFNKMatrixSt));
    MFSetError(e,12,RoutineName,MFNKMatrixErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  thisMatrix->n=-1;
  thisMatrix->k=k;

  if(!strcmp(MFNVGetId(col[0],e),"DENSE"))
   {
    thisMatrix->n=MFNV_NC(col[0],e);
    thisMatrix->data=(double*)malloc(thisMatrix->n*thisMatrix->k*sizeof(double));

#ifndef MFNOSAFETYNET
    if(thisMatrix->data==NULL)
     {
      sprintf(MFNKMatrixErrorMsg,"Out of memory, trying to allocate %d bytes",thisMatrix->n*thisMatrix->k*sizeof(double));
      MFSetError(e,12,RoutineName,MFNKMatrixErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      free(thisMatrix);
      return NULL;
     }
#endif

    thisMatrix->cols=NULL;

    for(j=0;j<k;j++)
     {
      u=MFNV_CStar(col[j],e);
      for(i=0;i<thisMatrix->n;i++)
       {
        thisMatrix->data[i+thisMatrix->n*j]=u[i];
       }
     }
    thisMatrix->type=DENSE;
   }else{
    thisMatrix->data=NULL;

    thisMatrix->cols=(MFNVector*)malloc(thisMatrix->k*sizeof(MFNVector));

#ifndef MFNOSAFETYNET
    if(thisMatrix->cols==NULL)
     {
      sprintf(MFNKMatrixErrorMsg,"Out of memory, trying to allocate %d bytes",thisMatrix->k*sizeof(MFNVector));
      MFSetError(e,12,RoutineName,MFNKMatrixErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      free(thisMatrix);
      return NULL;
     }
#endif

    for(j=0;j<k;j++)
     {
      thisMatrix->cols[j]=col[j];
      MFRefNVector(col[j],e);
     }

    if(!strcmp(MFNVGetId(col[0],e),"DENSE"))thisMatrix->type=DENSE;
     else thisMatrix->type=VECTOR;
   }

  thisMatrix->nRefs=1;

  return thisMatrix;
 }

MFNKMatrix MFCreateNKMatrixWithData(int n,int k, double *data,MFErrorHandler e)
 {
  static char RoutineName[]={"MFCreateNKMatrixWithData"};
  MFNKMatrix thisMatrix;
  int i,j;

#ifdef MFNOCONFIDENCE
  if(n<1)
   {
    sprintf(MFNKMatrixErrorMsg,"Number of Rows %d (argument 1) is Illegal. Must be positive.",n);
    MFSetError(e,12,RoutineName,MFNKMatrixErrorMsg,__LINE__,__FILE__);
    return NULL;
   }

  if(k<1)
   {
    sprintf(MFNKMatrixErrorMsg,"Number of Columns %d (argument 2) is Illegal. Must be positive.",k);
    MFSetError(e,12,RoutineName,MFNKMatrixErrorMsg,__LINE__,__FILE__);
    return NULL;
   }
#endif

  thisMatrix=(MFNKMatrix)malloc(sizeof(struct MFNKMatrixSt));

#ifndef MFNOSAFETYNET
  if(thisMatrix==NULL)
   {
    sprintf(MFNKMatrixErrorMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct MFNKMatrixSt));
    MFSetError(e,12,RoutineName,MFNKMatrixErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  thisMatrix->n=n;
  thisMatrix->k=k;

  thisMatrix->data=(double*)malloc(thisMatrix->n*thisMatrix->k*sizeof(double));

#ifndef MFNOSAFETYNET
  if(thisMatrix->data==NULL)
   {
    sprintf(MFNKMatrixErrorMsg,"Out of memory, trying to allocate %d bytes",thisMatrix->n*thisMatrix->k*sizeof(double));
    MFSetError(e,12,RoutineName,MFNKMatrixErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(thisMatrix);
    return NULL;
   }
#endif

  thisMatrix->nRefs=1;

  if(data!=NULL)
   {
    for(j=0;j<k;j++)
     {
      for(i=0;i<n;i++)
       {
        thisMatrix->data[i+n*j]=data[i+n*j];
       }
     }
   }else{
    for(i=0;i<n*k;i++)thisMatrix->data[i]=0.;
   }

  thisMatrix->cols=NULL;
  thisMatrix->type=DENSE;

  return thisMatrix;
 }

void MFMVMul(MFNSpace R, MFNKMatrix thisMatrix,MFKVector s,MFNVector u,MFErrorHandler e)
 {
  static char RoutineName[]={"MFMVMul"};
  int i,j;
  MFNVector col;
  MFNVector v;
  double t;
  static int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s\n",RoutineName);fflush(stdout);}
#endif

/* Need to use Add and Scale from NSpace */
  v=MFCloneNVector(u,e);

  for(i=0;i<thisMatrix->k;i++)
   {
    if(thisMatrix->type==DENSE)
      col=MFCreateNVectorWithData(thisMatrix->n,thisMatrix->data+i*thisMatrix->n,e);
     else
      col=MFMColumn(thisMatrix,i,e);

#ifdef MFNOCONFIDENCE
    if(col==NULL)
     {
      sprintf(MFNKMatrixErrorMsg,"col %d of thisMatrix is NULL");
      MFSetError(e,12,RoutineName,MFNKMatrixErrorMsg,__LINE__,__FILE__);
      return;
     }
#endif

    t=MFKV_C(s,i,e);
    if(i==0)
     {
      MFNSpaceScale(R,t,col,u,e);
     }else{
      MFNSpaceScale(R,t,col,v,e);
      MFNSpaceAdd(R,u,v,u,e);
     }
    MFFreeNVector(col,e);
   }
  MFFreeNVector(v,e);
  return;
 }

void MFMVMulT(MFNSpace R, MFNKMatrix thisMatrix,MFNVector u,MFKVector s,MFErrorHandler e)
 {
  static char RoutineName[]={"MFMVMulT"};
  int i,j;
  double t;
  MFNVector col;
  int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s, n=%d\n",RoutineName,thisMatrix->n);fflush(stdout);}
#endif

  for(j=0;j<thisMatrix->k;j++)
   {
    col=MFMColumn(thisMatrix,j,e);
    t=MFNSpaceInner(R,col,u,e);

#ifdef MFALLOWVERBOSE
    if(verbose)
     {
      printf("   s[%d]=",j);MFPrintNVector(stdout,col,e);
      printf(".");MFPrintNVector(stdout,u,e);
      printf("=%lf\n",t);
     }
#endif

    MFKVSetC(s,j,t,e);
    MFFreeNVector(col,e);
   }

  return;
 }

double *MFNKM_CStar(MFNKMatrix thisMatrix,MFErrorHandler e)
 {
  static char RoutineName[]={"MFNKM_CStar"};

#ifdef MFNOCONFIDENCE
  if(thisMatrix==NULL)
   {
    sprintf(MFNKMatrixErrorMsg,"Pointer to Matrix (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNKMatrixErrorMsg,__LINE__,__FILE__);
    exit(12);
    return NULL;
   }

  if(thisMatrix->type!=DENSE)
   {
    sprintf(MFNKMatrixErrorMsg,"Trying to get dense arrays from non-dense Vector type %d",thisMatrix->type);
    printf("%s\n",MFNKMatrixErrorMsg);fflush(stdout);
    MFSetError(e,12,RoutineName,MFNKMatrixErrorMsg,__LINE__,__FILE__);
    return NULL;
   }
#endif

  return thisMatrix->data;
 }

MFNVector MFMColumn(MFNKMatrix thisMatrix,int j,MFErrorHandler e)
 {
  static char RoutineName[]={"MFMColumn"};
  int i;
  MFNVector c;

#ifdef MFNOCONFIDENCE
  if(thisMatrix==NULL)
   {
    sprintf(MFNKMatrixErrorMsg,"Pointer to Matrix (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNKMatrixErrorMsg,__LINE__,__FILE__);
    return NULL;
   }

  if(j>thisMatrix->k)
   {
    sprintf(MFNKMatrixErrorMsg,"Request for invalid column, %d must be in [0,%d)",j,thisMatrix->k);
    MFSetError(e,12,RoutineName,MFNKMatrixErrorMsg,__LINE__,__FILE__);
    return NULL;
   }
#endif

  switch(thisMatrix->type)
   {
    case DENSE:
     c=MFCreateWrappedNVector(thisMatrix->n,thisMatrix->data+j*thisMatrix->n,e);
     break;
    default:
     c=thisMatrix->cols[j];
     MFRefNVector(c,e);
   }

  return c;
 }

void MFMRow(MFNKMatrix thisMatrix,int i,MFKVector r,MFErrorHandler e)
 {
  static char RoutineName[]={"MFMRow"};
  int j;

  switch(thisMatrix->type)
   {
    case DENSE:
     for(j=0;j<thisMatrix->k;j++)
       MFKVSetC(r,j,thisMatrix->data[i+thisMatrix->n*j],e);
     break;
    default:
     for(j=0;j<thisMatrix->k;j++)
       MFKVSetC(r,j,MFNV_C(thisMatrix->cols[j],i,e),e);
     break;
   }

  return;
 }

int MFNKMatrixK(MFNKMatrix thisMatrix,MFErrorHandler e)
 {
  return thisMatrix->k;
 }

int MFNKMatrixN(MFNKMatrix thisMatrix,MFErrorHandler e)
 {
  return thisMatrix->n;
 }

void MFWriteNKMatrix(FILE *fid,MFNKMatrix L,MFErrorHandler e)
 {
  static char RoutineName[]={"MFWriteNKMatrix"};
  int i;

  fprintf(fid,"%s\n","NKMatrix");
  fprintf(fid,"%d %d %d %d\n",L->type,L->n,L->k,L->nRefs);
  if(L->type==DENSE)
   {
    for(i=0;i<L->n*L->k;i++)
     {
      if(i>0)fprintf(fid," ");
      fprintf(fid,"%lf",L->data[i]);
     }
    fprintf(fid,"\n");
   }else{
    for(i=0;i<L->k;i++)
      MFWriteNVector(fid,L->cols[i],e);
   }

  return;
 }

MFNKMatrix MFReadNKMatrix(FILE *fid,MFErrorHandler e)
 {
  static char RoutineName[]={"MFReadNKMatrix"};
  int i;
  MFNKMatrix L;
  char tag[100]="";

  fscanf(fid,"%s\n",tag);

#ifdef MFNOCONFIDENCE
  if(strcmp(tag,"NKMatrix"))
   {
    sprintf(MFNKMatrixErrorMsg,"Next Object is not a NKMatrix! (%s)\n",RoutineName,tag);
    MFSetError(e,12,RoutineName,MFNKMatrixErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  L=(MFNKMatrix)malloc(sizeof(struct MFNKMatrixSt));

#ifndef MFNOSAFETYNET
  if(L==NULL)
   {
    sprintf(MFNKMatrixErrorMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct MFNKMatrixSt));
    MFSetError(e,12,RoutineName,MFNKMatrixErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  fscanf(fid,"%d %d %d %d\n",&(L->type),&(L->n),&(L->k),&(L->nRefs));

  if(L->type==DENSE)
   {
    L->data=(double*)malloc(L->n*L->k*sizeof(double));

#ifndef MFNOSAFETYNET
    if(L->data==NULL)
     {
      sprintf(MFNKMatrixErrorMsg,"Out of memory, trying to allocate %d bytes",L->n*L->k*sizeof(double));
      MFSetError(e,12,RoutineName,MFNKMatrixErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return NULL;
     }
#endif

    for(i=0;i<L->n*L->k;i++)
     {
      if(i>0)fscanf(fid," ");
      fscanf(fid,"%lf",&(L->data[i]));
     }
    fscanf(fid,"\n");
    L->cols=NULL;
   }else{
    L->cols=(MFNVector*)malloc(L->k*sizeof(MFNVector));

#ifndef MFNOSAFETYNET
    if(L->cols==NULL)
     {
      sprintf(MFNKMatrixErrorMsg,"Out of memory, trying to allocate %d bytes",L->k*sizeof(MFNVector));
      MFSetError(e,12,RoutineName,MFNKMatrixErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      free(L);
      return NULL;
     }
#endif

    for(i=0;i<L->k;i++)
      L->cols[i]=MFReadNVector(fid,e);
    L->data=NULL;
   }

  return L;
 }

void MFNKMSetC(MFNKMatrix A,int i,int j,double Aij,MFErrorHandler e)
 {
  static char RoutineName[]={"MFNKMSetC"};

#ifdef MFNOCONFIDENCE
  if(A==NULL)
   {
    sprintf(MFNKMatrixErrorMsg,"Pointer to Matrix (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFNKMatrixErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(i<0)
   {
    sprintf(MFNKMatrixErrorMsg,"column index (argument 2) is less than zero");
    MFSetError(e,4,RoutineName,MFNKMatrixErrorMsg,__LINE__,__FILE__);
    return;
   }
  if(i>=A->n)
   {
    sprintf(MFNKMatrixErrorMsg,"column index (argument 2) is too large must be < %d)",A->n);
    MFSetError(e,4,RoutineName,MFNKMatrixErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(j<0)
   {
    sprintf(MFNKMatrixErrorMsg,"row index (argument 3) is less than zero");
    MFSetError(e,4,RoutineName,MFNKMatrixErrorMsg,__LINE__,__FILE__);
    return;
   }
  if(j>=A->k)
   {
    sprintf(MFNKMatrixErrorMsg,"row index (argument 3) is too large must be < %d)",A->k);
    MFSetError(e,4,RoutineName,MFNKMatrixErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  switch(A->type)
   {
    case DENSE:
     A->data[i+A->n*j]=Aij;
     break;
    default:
     MFNVSetC(A->cols[j],i,Aij,e);
     break;
   }

  return;
 }

double MFNKMGetC(MFNKMatrix A,int i,int j,MFErrorHandler e)
 {
  static char RoutineName[]={"MFNKMGetC"};
  double Aij;

#ifdef MFNOCONFIDENCE
  if(A==NULL)
   {
    sprintf(MFNKMatrixErrorMsg,"Pointer to Matrix (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFNKMatrixErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(i<0)
   {
    sprintf(MFNKMatrixErrorMsg,"column index (argument 2) is less than zero");
    MFSetError(e,4,RoutineName,MFNKMatrixErrorMsg,__LINE__,__FILE__);
    return;
   }
  if(i>=A->n)
   {
    sprintf(MFNKMatrixErrorMsg,"column index (argument 2) is too large must be < %d)",A->n);
    MFSetError(e,4,RoutineName,MFNKMatrixErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(j<0)
   {
    sprintf(MFNKMatrixErrorMsg,"row index (argument 3) is less than zero");
    MFSetError(e,4,RoutineName,MFNKMatrixErrorMsg,__LINE__,__FILE__);
    return;
   }
  if(j>=A->k)
   {
    sprintf(MFNKMatrixErrorMsg,"row index (argument 3) is too large must be < %d)",A->k);
    MFSetError(e,4,RoutineName,MFNKMatrixErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  switch(A->type)
   {
    case DENSE:
     Aij=A->data[i+A->n*j];
     break;
    default:
     Aij=MFNV_C(A->cols[j],i,e);
     break;
   }

  return Aij;
 }

void MFRefNKMatrix(MFNKMatrix thisMatrix,MFErrorHandler e)
 {
  static char RoutineName[]={"MFRefNKMatrix"};

#ifdef MFNOCONFIDENCE
  if(thisMatrix==NULL)
   {
    sprintf(MFNKMatrixErrorMsg,"Pointer to Matrix (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFNKMatrixErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  thisMatrix->nRefs++;

  return;
 }

void MFGramSchmidt(MFNSpace space,MFNKMatrix A,MFErrorHandler e)
 {
  static char RoutineName[]={"MFGramSchmidt"};
  double *a;
  int i,j,n;
  int k;
  MFNVector u;
  MFNVector v;
  MFNVector w;
  double inner;
  int verbose=0;

  n=MFNKMatrixN(A,e);
  k=MFNKMatrixK(A,e);

  switch(A->type)
   {
    case DENSE:
     a=MFNKM_CStar(A,e);
     MFGramSchmidtNoMat(n,k,a,e);
     break;
    default:
     u=MFMColumn(A,0,e);
     w=MFCloneNVector(u,e);
     MFFreeNVector(u,e);

#ifdef MFALLOWVERBOSE
     if(verbose){printf("Gram Schmidt:\n");fflush(stdout);}
#endif

     for(i=0;i<k;i++)
      {
       u=MFMColumn(A,i,e);

       inner=MFNSpaceInner(space,u,u,e);

#ifdef MFALLOWVERBOSE
       if(verbose){printf("u_%d\n",i);MFPrintNVector(stdout,u,e);printf("\n");fflush(stdout);
                   printf("<u_%d,u_%d>=%lf\n",i,i,inner);fflush(stdout);}
#endif

       inner=1./sqrt(inner);
       MFNSpaceScale(space,inner,u,u,e);

#ifdef MFALLOWVERBOSE
       if(verbose){printf("normalized u_%d\n",i);MFPrintNVector(stdout,u,e);printf("\n");fflush(stdout);}
#endif
   
       for(j=i+1;j<k;j++)
        {
         v=MFMColumn(A,j,e);

#ifdef MFALLOWVERBOSE
         if(verbose){printf("u_%d\n",j);MFPrintNVector(stdout,v,e);printf("\n");fflush(stdout);}
#endif
   
         inner=MFNSpaceInner(space,u,v,e);

#ifdef MFALLOWVERBOSE
         if(verbose){printf("<u_%d,u_%d>=%lf\n",i,j,inner);fflush(stdout);}
#endif
   
         MFNSpaceScale(space,inner,u,w,e);

#ifdef MFALLOWVERBOSE
         if(verbose){printf("<u_%d,u_%d>u_%d\n",j,i,i);MFPrintNVector(stdout,w,e);printf("\n");fflush(stdout);}
#endif
         MFNVDiff(v,w,v,e);

#ifdef MFALLOWVERBOSE
         if(verbose){printf("u_%d-<u_%d,u_%d>u_%d\n",j,j,i,i);MFPrintNVector(stdout,v,e);printf("\n");fflush(stdout);}
#endif
   
         MFFreeNVector(v,e);
        }
       MFFreeNVector(u,e);
      }
     MFFreeNVector(w,e);
   }

  return;
 }

void MFGramSchmidtNoMat(int n, int k, double *a,MFErrorHandler e)
 {
  static char RoutineName[]={"MFGramSchmidtNoMat"};
  double inner;
  int i,j,l;
  double *u;
  double *v;
  int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s\n",RoutineName);fflush(stdout);}
#endif

  for(i=0;i<k;i++)
   {
    u=a+n*i;

    inner=0.;for(l=0;l<n;l++)inner+=u[l]*u[l];
    inner=1./sqrt(inner);
    for(l=0;l<n;l++)u[l]=u[l]*inner;

    for(j=i+1;j<k;j++)
     {
      v=a+n*j;
      inner=0.;for(l=0;l<n;l++)inner+=u[l]*v[l];
      for(l=0;l<n;l++)v[l]=v[l]-inner*u[l];
     }
   }

  return;
 }

void MFNKMProjectTangentForBranchSwitch(MFNSpace space, MFNKMatrix A,MFNVector phi,MFNKMatrix Phi,MFErrorHandler e)
 {
  static char RoutineName[]={"MFNKMProjectTangentForBranchSwitch"};

  int i,j;
  int k,m;
  MFNVector u,w;
  double *p;
  double pmax;
  MFNVector pi,pj;
  int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s\n",RoutineName);fflush(stdout);
              printf("  A=");MFPrintNKMatrix(stdout,A,e);fflush(stdout);
              printf("  phi=");MFPrintNVector(stdout,phi,e);printf("\n");fflush(stdout);}
#endif

  k=MFNKMatrixK(A,e);

  p=(double*)malloc(k*sizeof(double));

#ifndef MFNOSAFETYNET
  if(p==NULL)
   {
    sprintf(MFNKMatrixErrorMsg,"Out of memory, trying to allocate %d bytes",k*sizeof(double));
    MFSetError(e,12,RoutineName,MFNKMatrixErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  j=0;
  pmax=0.;
  for(i=0;i<k;i++)
   {
    u=MFMColumn(A,i,e);
    p[i]=MFNSpaceInner(space,phi,u,e);
    if(fabs(p[i])>pmax){pmax=fabs(p[i]);j=i;}
    MFFreeNVector(u,e);

#ifdef MFALLOWVERBOSE
    if(verbose){printf("phi^* A_[%d]=%lf\n",i,p[i]);fflush(stdout);}
#endif

   }

#ifdef MFALLOWVERBOSE
  if(verbose){printf("max of phi^* A_j is for j=%d, =%lf\n",j,pmax);fflush(stdout);}
#endif


  m=0;
  u=MFMColumn(Phi,m,e);
  w=MFCloneNVector(u,e);
  MFNSpaceScale(space,1.,phi,u,e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("column 0 of Phi is phi ");MFPrintNVector(stdout,u,e);printf("\n");(stdout);}
#endif

  if(Phi->type==DENSE)MFMSetColumn(Phi,m,u,e);
  MFFreeNVector(u,e);
  m++;

  pj=MFMColumn(A,j,e);
  for(i=0;i<k;i++)
   {
    if(i!=j)
     {
      pi=MFMColumn(A,i,e);
      u=MFMColumn(Phi,m,e);

      MFNSpaceScale(space, p[i],pj,u,e);
      MFNSpaceScale(space,-p[j],pi,w,e);
      MFNSpaceAdd(space,u,w,u,e);

#ifdef MFALLOWVERBOSE
      if(verbose){printf("column %d of Phi is %lf A_%d - %lf A_%d ",m,p[i],j,p[j],i);MFPrintNVector(stdout,u,e);printf("\n");(stdout);}
#endif

      if(Phi->type==DENSE)MFMSetColumn(Phi,m,u,e);
      m++;
      MFFreeNVector(u,e);
      MFFreeNVector(pi,e);
     }
   }
  MFFreeNVector(pj,e);
  free(p);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("Before G-S:\n");MFPrintNKMatrix(stdout,Phi,e);fflush(stdout);}
#endif

  MFGramSchmidt(space,Phi,e);


#ifdef MFALLOWVERBOSE
  if(verbose){printf("After G-S:\n");MFPrintNKMatrix(stdout,Phi,e);fflush(stdout);
              printf("done %s\n",RoutineName);fflush(stdout);}
#endif

  return;
 }

MFNKMatrix MFCloneNKMatrix(MFNKMatrix A,MFErrorHandler e)
 {
  static char RoutineName[]={"MFCloneNKMatrix"};
  MFNKMatrix thisMatrix;
  int i,j;
  double *u;

#ifdef MFNOCONFIDENCE
  if(A==NULL)
   {
    sprintf(MFNKMatrixErrorMsg,"Cloning a NULL Matrix!");
    MFSetError(e,12,RoutineName,MFNKMatrixErrorMsg,__LINE__,__FILE__);
    return NULL;
   }
#endif

  thisMatrix=(MFNKMatrix)malloc(sizeof(struct MFNKMatrixSt));

#ifndef MFNOSAFETYNET
  if(thisMatrix==NULL)
   {
    sprintf(MFNKMatrixErrorMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct MFNKMatrixSt));
    MFSetError(e,12,RoutineName,MFNKMatrixErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  thisMatrix->n=A->n;
  thisMatrix->k=A->k;
  thisMatrix->type=A->type;
  if(thisMatrix->type==DENSE)
   {
    thisMatrix->data=(double*)malloc(thisMatrix->n*thisMatrix->k*sizeof(double));

#ifndef MFNOSAFETYNET
    if(thisMatrix->data==NULL)
     {
      sprintf(MFNKMatrixErrorMsg,"Out of memory, trying to allocate %d bytes",thisMatrix->n*thisMatrix->k*sizeof(double));
      MFSetError(e,12,RoutineName,MFNKMatrixErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      free(thisMatrix);
      return NULL;
     }
#endif

    thisMatrix->cols=NULL;

    for(j=0;j<thisMatrix->k;j++)
     {
      for(i=0;i<thisMatrix->n;i++)thisMatrix->data[i+thisMatrix->n*j]=A->data[i+thisMatrix->n*j];
     }
   }else{
    thisMatrix->data=NULL;

    thisMatrix->cols=(MFNVector*)malloc(thisMatrix->k*sizeof(MFNVector));

#ifndef MFNOSAFETYNET
    if(thisMatrix->cols==NULL)
     {
      sprintf(MFNKMatrixErrorMsg,"Out of memory, trying to allocate %d bytes",thisMatrix->k*sizeof(MFNVector));
      MFSetError(e,12,RoutineName,MFNKMatrixErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      free(thisMatrix);
      return NULL;
     }
#endif

    for(j=0;j<thisMatrix->k;j++)
      thisMatrix->cols[j]=MFCloneNVector(A->cols[j],e);
   }

  thisMatrix->nRefs=1;

  return thisMatrix;
 }

void MFMSetColumn(MFNKMatrix thisMatrix,int j, MFNVector u,MFErrorHandler e)
 {
  static char RoutineName[]={"MFMSetColumn"};
  int i;
  MFNVector c;

#ifdef MFNOCONFIDENCE
  if(j<0 || j >= thisMatrix->k)
   {
    sprintf(MFNKMatrixErrorMsg,"Column %d must be between 0 and %d",j,thisMatrix->k);
    MFSetError(e,12,RoutineName,MFNKMatrixErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(thisMatrix);
    return;
   }
#endif

  switch(thisMatrix->type)
   {
    case DENSE:
     for(i=0;i<thisMatrix->n;i++)
      thisMatrix->data[i+thisMatrix->n*j]=MFNV_C(u,i,e);
     break;
    case VECTOR:
      MFRefNVector(u,e);
      MFFreeNVector(thisMatrix->cols[j],e);
      thisMatrix->cols[j]=u;
     break;
    default:
     for(i=0;i<thisMatrix->n;i++)
      MFNVSetC(thisMatrix->cols[j],i,MFNV_C(u,i,e),e);
   }

  return;
 }

void MFGramSchmidtReplace(MFNSpace space,MFNKMatrix A,MFNVector phi0, MFNVector phi1,MFErrorHandler e)
 {
  static char RoutineName[]={"MFGramSchmidtReplace"};
  double *a;
  int i,j,n;
  int k;
  MFNVector u;
  MFNVector v;
  MFNVector w;
  double inner;
  int zero;
  double *B=NULL;
  int sign0,sign1;
  int verbose=0;

/* Remove vector phi0 from the basis, and replaces it with the vector phi1 */

  if(verbose){printf("%s\n",RoutineName);fflush(stdout);}

  n=MFNKMatrixN(A,e);
  k=MFNKMatrixK(A,e);

  B=(double*)malloc(k*k*sizeof(double));

#ifndef MFNOSAFETYNET
    if(B==NULL)
     {
      sprintf(MFNKMatrixErrorMsg,"Out of memory, trying to allocate %d bytes",k*k*sizeof(double));
      MFSetError(e,12,RoutineName,MFNKMatrixErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

  if(verbose)
   {
    printf("Remove vector phi0 from the basis, and replace it with the vector phi1\n");fflush(stdout);
    printf(" phi0=0x%8.8x ",phi0);MFPrintNVector(stdout,phi0,e);printf("\n");fflush(stdout);
    printf(" phi1=0x%8.8x ",phi0);MFPrintNVector(stdout,phi1,e);printf("\n");fflush(stdout);
   }

  for(i=0;i<k;i++)
   {
    u=MFMColumn(A,i,e);
    for(j=0;j<k;j++)
     {
      v=MFMColumn(A,j,e);
      B[i+k*j]=MFNSpaceInner(space,u,v,e);
      MFFreeNVector(v,e);
     }
    MFFreeNVector(u,e);
   }
  sign0=MFAtlasDet(k,B,e)>0;

  u=MFMColumn(A,0,e);
  w=MFCloneNVector(u,e);
  MFFreeNVector(u,e);

  zero=-1;
  for(i=-1;i<k;i++)
   {
    if(i>-1)u=MFMColumn(A,i,e);
     else{u=phi0;MFRefNVector(phi0,e);}

    inner=MFNSpaceInner(space,u,u,e);

#ifdef MFALLOWVERBOSE
    if(verbose){printf("u_%d\n",i);MFPrintNVector(stdout,u,e);printf("\n");fflush(stdout);
                printf("<u_%d,u_%d>=%lf\n",i,i,inner);fflush(stdout);}
#endif

    if(inner>1.e-7)
     {
      inner=1./sqrt(inner);
      MFNSpaceScale(space,inner,u,u,e);

#ifdef MFALLOWVERBOSE
      if(verbose){printf("normalized u_%d\n",i);MFPrintNVector(stdout,u,e);printf("\n");fflush(stdout);}
#endif

     }else{
      zero=i;

#ifdef MFALLOWVERBOSE
      if(verbose){printf("zero vector u_%d\n",i);MFPrintNVector(stdout,u,e);printf("\n");fflush(stdout);}
#endif

     }
 
    for(j=i+1;j<k;j++)
     {
      v=MFMColumn(A,j,e);

#ifdef MFALLOWVERBOSE
      if(verbose){printf("v=u_%d\n",j);MFPrintNVector(stdout,v,e);printf("\n");fflush(stdout);}
#endif

      inner=MFNSpaceInner(space,u,v,e);

#ifdef MFALLOWVERBOSE
      if(verbose){printf("inner=<u_%d,u_%d>=%lf\n",i,j,inner);fflush(stdout);}
#endif

      MFNSpaceScale(space,inner,u,w,e);

#ifdef MFALLOWVERBOSE
      if(verbose){printf("w=<u_%d,u_%d>u_%d\n",j,i,i);MFPrintNVector(stdout,w,e);printf("\n");fflush(stdout);}
#endif

      MFNVDiff(v,w,v,e);

#ifdef MFALLOWVERBOSE
      if(verbose){printf("v-w=u_%d-<u_%d,u_%d>u_%d\n",j,j,i,i);MFPrintNVector(stdout,v,e);printf("\n");fflush(stdout);}
#endif

      MFFreeNVector(v,e);
     }
    MFFreeNVector(u,e);

#ifdef MFALLOWVERBOSE
    if(verbose){printf("\n");fflush(stdout);}
#endif

   }
  MFFreeNVector(w,e);

  if(verbose){printf("done orthnormalizing, line %d, routine %s\n",__LINE__,RoutineName);fflush(stdout);}

  if(verbose){printf("Delete the zero column (move it to the first position) zero=%d, line %d, routine %s\n",zero,__LINE__,RoutineName);fflush(stdout);}

  if(zero>1)
   {
    for(i=zero;i>-1;i--)
     {
      u=MFMColumn(A,i-1,e);
      v=MFMColumn(A,i,e);
      MFMSetColumn(A,i,u,e);
      MFFreeNVector(u,e);
      MFFreeNVector(v,e);
     }
    if(verbose){printf("done removing the zero column line %d, routine %s\n",__LINE__,RoutineName);fflush(stdout);}
   }

/* Replace the zero column with phi1 */

  u=MFMColumn(A,0,e);
  MFFreeNVector(u,e);
  MFMSetColumn(A,0,phi1,e);

  if(verbose){printf("Gram-Schmidt on new basis line %d, routine %s\n",__LINE__,RoutineName);fflush(stdout);}
 MFGramSchmidt(space,A,e);
  if(verbose){printf("done Gram-Schmidt on new basis line %d, routine %s\n",__LINE__,RoutineName);fflush(stdout);}

  for(i=0;i<k;i++)
   {
    u=MFMColumn(A,i,e);
    for(j=0;j<k;j++)
     {
      v=MFMColumn(A,j,e);
      B[i+k*j]=MFNSpaceInner(space,u,v,e);
      MFFreeNVector(v,e);
     }
    MFFreeNVector(u,e);
   }
  sign1=MFAtlasDet(k,B,e)>0;

  if(verbose){printf("make Det the same line %d, routine %s\n",__LINE__,RoutineName);fflush(stdout);}

  if(sign0&&!sign1||sign1&&!sign0)
   {
    u=MFMColumn(A,0,e);
    MFNSpaceScale(space,-1.,u,u,e);
    MFFreeNVector(u,e);
   }

  free(B);
  if(verbose){printf("done %s\n",RoutineName);fflush(stdout);}

  return;
 }

void MFMMMul(MFNSpace R, MFNKMatrix Phi0,MFNKMatrix Phi1,double *prod,MFErrorHandler e)
 {
  static char RoutineName[]={"MFMMMul"};
  int i,j;
  MFNVector col0,col1;
  static int verbose;

/* Phi_0^T Phi_1 */

  verbose=0;
  if(verbose){printf("%s\n",RoutineName);fflush(stdout);}

#ifdef MFNOCONFIDENCE
  if(Phi0->n!=Phi1->n)
   {
    sprintf(MFNKMatrixErrorMsg,"Can't do product of incompatible matrices. A.n=%d, B.n=%d",Phi0->n,Phi1->n);
    MFSetError(e,12,RoutineName,MFNKMatrixErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(Phi0->type!=Phi1->type)
   {
    sprintf(MFNKMatrixErrorMsg,"Can't do product of incompatible matrices. A.type=%d, B.type=%d",Phi0->type,Phi1->type);
    MFSetError(e,12,RoutineName,MFNKMatrixErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  for(i=0;i<Phi0->k;i++)
   {
    if(Phi0->type==DENSE)
      col0=MFCreateNVectorWithData(Phi0->n,Phi0->data+i*Phi0->n,e);
     else
      col0=MFMColumn(Phi0,i,e);

#ifdef MFNOCONFIDENCE
    if(col0==NULL)
     {
      sprintf(MFNKMatrixErrorMsg,"col %d of A is NULL",i);
      MFSetError(e,12,RoutineName,MFNKMatrixErrorMsg,__LINE__,__FILE__);
      return;
     }
#endif

    for(j=0;j<Phi1->k;j++)
     {
      if(Phi1->type==DENSE)
        col1=MFCreateNVectorWithData(Phi1->n,Phi1->data+i*Phi1->n,e);
       else
        col1=MFMColumn(Phi1,j,e);
  
#ifdef MFNOCONFIDENCE
      if(col1==NULL)
       {
        sprintf(MFNKMatrixErrorMsg,"col %d of B is NULL",j);
        MFSetError(e,12,RoutineName,MFNKMatrixErrorMsg,__LINE__,__FILE__);
        return;
       }
#endif

      prod[i+Phi0->k*j]=MFNSpaceInner(R,col0,col1,e);
      MFFreeNVector(col1,e);
     }
    MFFreeNVector(col0,e);
   }
  return;
 }

void MFPrintNKMatrix(FILE *fid,MFNKMatrix L,MFErrorHandler e)
 {
  static char RoutineName[]={"MFPrintNKMatrix"};
  int i,n,k;
  MFKVector s;

#ifdef MFNOCONFIDENCE
  if(fid==NULL)
   {
    sprintf(MFNKMatrixErrorMsg,"fid (argument 1) is NULL.");
    MFSetError(e,12,RoutineName,MFNKMatrixErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(L==NULL)
   {
    sprintf(MFNKMatrixErrorMsg,"Matrix (argument 2) is NULL.");
    MFSetError(e,12,RoutineName,MFNKMatrixErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  n=MFNKMatrixN(L,e);
  k=MFNKMatrixK(L,e);
  s=MFCreateKVector(k,e);
  if(L->type==VECTOR)
   {
    printf("vector\n");fflush(stdout);
    printf("[\n");
    for(i=0;i<k;i++)
     {
      MFPrintNVector(fid,L->cols[i],e);
      fprintf(fid,"\n");fflush(stdout);
     }
    printf("]^T\n");
   }else if(n<16)
   {
    for(i=0;i<n;i++)
     {
      MFMRow(L,i,s,e);
      MFPrintKVector(fid,s,e);
      fprintf(fid,"\n");
     }
   }else{
    printf("dense\n");fflush(stdout);
    for(i=0;i<8;i++)
     {
      MFMRow(L,i,s,e);
      MFPrintKVector(fid,s,e);
      fprintf(fid,"\n");
     }
    fprintf(fid,"...\n");
    for(i=n-8;i<n;i++)
     {
      MFMRow(L,i,s,e);
      MFPrintKVector(fid,s,e);
      fprintf(fid,"\n");
     }
   }
  MFFreeKVector(s,e);
  return;
 }

#ifdef __cplusplus
}
#endif
