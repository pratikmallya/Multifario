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

static char *id="@(#) $Id: MFBinaryTree.c,v 1.4 2011/07/21 17:42:46 mhender Exp $";

static char MFBinaryTreeErrorMsg[256]="";

#include <MFBinaryTree.h>
#include <MFListOfCharts.h>
#include <MFAtlas.h>
#include <MFPrint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int MFBINARYTREESEQUENCENUMBER=0;
#define MFBINARYTREEMAXNUMBEROFCHARTS 5

#ifdef __cplusplus
 extern "C" {
#endif

struct MFBinaryTreeLeafSt;
typedef struct MFBinaryTreeLeafSt *MFBinaryTreeLeaf;

struct MFBinaryTreeSt {
                       int nCharts;
                       int k;
                       MFBinaryTreeLeaf root;
                      };

struct MFBinaryTreeLeafSt {
                     int sequenceNumber;
                     int nCharts;
                     int *chart;
                     double *chartR;
                     double **chartCenter;
                     double *value;
                     int d;
                     double split;
                     double *leftB;
                     double *rightB;
                     MFBinaryTreeLeaf parent;
                     MFBinaryTreeLeaf leftT;
                     MFBinaryTreeLeaf rightT;
                    };

static MFBinaryTreeLeaf MFCreateBinaryTreeLeaf(int,int,MFErrorHandler);
static void MFFreeBinaryTreeLeaf(MFBinaryTreeLeaf,MFErrorHandler);
static void MFBinaryTreeLeafAddChart(MFBinaryTreeLeaf,int,int,double,double*,MFErrorHandler);
static void MFBinaryTreeGetIntersectingCharts(MFBinaryTreeLeaf,int,int,double*,double,int*,int**,MFErrorHandler);
static void MFBinaryTreeGetNeighboringCharts(MFBinaryTreeLeaf,int,double*,double,int*,int**,MFErrorHandler);

static int nInTree=0;

MFBinaryTree MFCreateBinaryTree(int k, MFErrorHandler e)
 {
  static char RoutineName[]={"MFCreateBinaryTree"};
  int i;
  MFBinaryTree T;
  MFBinaryTreeLeaf TLeaf;
  int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s\n",RoutineName);fflush(stdout);}
#endif

#ifndef MFNOSAFETYNET
  if(k<1)
   {
    sprintf(MFBinaryTreeErrorMsg,"Dimension %d, (Argument 1) must be positive.",k);
    MFSetError(e,12,RoutineName,MFBinaryTreeErrorMsg,__LINE__,__FILE__);
    return NULL;
   }
#endif

  T=(MFBinaryTree)malloc(sizeof(struct MFBinaryTreeSt));

#ifndef MFNOSAFETYNET
  if(T==NULL)
   {
    sprintf(MFBinaryTreeErrorMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct MFBinaryTreeSt));
    MFSetError(e,12,RoutineName,MFBinaryTreeErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s, T=0x%8.8x\n",RoutineName,T);fflush(stdout);}
#endif

  T->nCharts=0;
  T->k=k;
  TLeaf=MFCreateBinaryTreeLeaf(T->k,0,e);
  TLeaf->parent=NULL;
  T->root=TLeaf;
  nInTree=0;

  return T;
 }

MFBinaryTreeLeaf MFCreateBinaryTreeLeaf(int k,int d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFCreateBinaryTreeLeaf"};
  int i;
  MFBinaryTreeLeaf TLeaf;
  double t;
  int verbose=0;

  if(d<0 || !(d<k))
   {
    sprintf(MFBinaryTreeErrorMsg,"direction %d (argument 3) is invalid. Must be positive and less than %d",d,k);
    MFSetError(e,12,RoutineName,MFBinaryTreeErrorMsg,__LINE__,__FILE__);
    return NULL;
   }

  TLeaf=(MFBinaryTreeLeaf)malloc(sizeof(struct MFBinaryTreeLeafSt));

#ifndef MFNOSAFETYNET
  if(TLeaf==NULL)
   {
    sprintf(MFBinaryTreeErrorMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct MFBinaryTreeLeafSt));
    MFSetError(e,12,RoutineName,MFBinaryTreeErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s, T=0x%8.8x\n",RoutineName,TLeaf);fflush(stdout);}
#endif

  TLeaf->sequenceNumber=MFBINARYTREESEQUENCENUMBER;
  MFBINARYTREESEQUENCENUMBER++;
  TLeaf->d=d;
  TLeaf->split=0.;
  TLeaf->nCharts=0;

  TLeaf->chart=(int*)malloc(MFBINARYTREEMAXNUMBEROFCHARTS*sizeof(int));

#ifndef MFNOSAFETYNET
  if(TLeaf->chart==NULL)
   {
    sprintf(MFBinaryTreeErrorMsg,"Out of memory, trying to allocate %d bytes",MFBINARYTREEMAXNUMBEROFCHARTS*sizeof(int));
    MFSetError(e,12,RoutineName,MFBinaryTreeErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

#ifdef MFALLOWVERBOSE
  if(verbose){printf("     T->chart=0x%8.8x\n",TLeaf->chart);fflush(stdout);}
#endif

  TLeaf->chartR=(double*)malloc(MFBINARYTREEMAXNUMBEROFCHARTS*sizeof(double));

#ifndef MFNOSAFETYNET
  if(TLeaf->chartR==NULL)
   {
    sprintf(MFBinaryTreeErrorMsg,"Out of memory, trying to allocate %d bytes",MFBINARYTREEMAXNUMBEROFCHARTS*sizeof(double));
    MFSetError(e,12,RoutineName,MFBinaryTreeErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

#ifdef MFALLOWVERBOSE
  if(verbose){printf("     T->chartR=0x%8.8x\n",TLeaf->chartR);fflush(stdout);}
#endif

  TLeaf->chartCenter=(double**)malloc(MFBINARYTREEMAXNUMBEROFCHARTS*k*sizeof(double*));

#ifndef MFNOSAFETYNET
  if(TLeaf->chartCenter==NULL)
   {
    sprintf(MFBinaryTreeErrorMsg,"Out of memory, trying to allocate %d bytes",MFBINARYTREEMAXNUMBEROFCHARTS*k*sizeof(double));
    MFSetError(e,12,RoutineName,MFBinaryTreeErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

#ifdef MFALLOWVERBOSE
  if(verbose){printf("     T->chartCenter=0x%8.8x\n",TLeaf->chartCenter);fflush(stdout);}
#endif

  TLeaf->value=(double*)malloc(MFBINARYTREEMAXNUMBEROFCHARTS*sizeof(double));

#ifndef MFNOSAFETYNET
  if(TLeaf->value==NULL)
   {
    sprintf(MFBinaryTreeErrorMsg,"Out of memory, trying to allocate %d bytes",MFBINARYTREEMAXNUMBEROFCHARTS*sizeof(double));
    MFSetError(e,12,RoutineName,MFBinaryTreeErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

#ifdef MFALLOWVERBOSE
  if(verbose){printf("     T->value=0x%8.8x\n",TLeaf->value);fflush(stdout);}
#endif

  TLeaf->leftB=(double*)malloc(k*sizeof(double));

#ifndef MFNOSAFETYNET
  if(TLeaf->leftB==NULL)
   {
    sprintf(MFBinaryTreeErrorMsg,"Out of memory, trying to allocate %d bytes",MFBINARYTREEMAXNUMBEROFCHARTS*sizeof(double));
    MFSetError(e,12,RoutineName,MFBinaryTreeErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

#ifdef MFALLOWVERBOSE
  if(verbose){printf("     T->leftB=0x%8.8x\n",TLeaf->leftB);fflush(stdout);}
#endif

  TLeaf->rightB=(double*)malloc(k*sizeof(double));

#ifndef MFNOSAFETYNET
  if(TLeaf->rightB==NULL)
   {
    sprintf(MFBinaryTreeErrorMsg,"Out of memory, trying to allocate %d bytes",MFBINARYTREEMAXNUMBEROFCHARTS*sizeof(double));
    MFSetError(e,12,RoutineName,MFBinaryTreeErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

#ifdef MFALLOWVERBOSE
  if(verbose){printf("     T->rightB=0x%8.8x\n",TLeaf->rightB);fflush(stdout);}
#endif

  for(i=0;i<k;i++)
   {
    TLeaf->leftB[i]=1.e20;
    TLeaf->rightB[i]=-1.e20;
   }

  TLeaf->parent=NULL;
  TLeaf->leftT=NULL;
  TLeaf->rightT=NULL;

  return TLeaf;
 }

void MFFreeBinaryTree(MFBinaryTree T, MFErrorHandler e)
 {
  static char RoutineName[]={"MFFreeBinaryTree"};
  int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s, T=0x%8.8x\n",RoutineName,T);fflush(stdout);}
#endif

  if(T==NULL)return;

  if(T->root!=NULL)MFFreeBinaryTreeLeaf(T->root,e);

  free(T);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("done %s, T=0x%8.8x\n",RoutineName,T);fflush(stdout);}
#endif

  return;
 }

void MFFreeBinaryTreeLeaf(MFBinaryTreeLeaf T, MFErrorHandler e)
 {
  static char RoutineName[]={"MFFreeBinaryTreeLeaf"};
  int i;
  int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("%s, T=0x%8.8x\n",RoutineName,T);fflush(stdout);
    printf("     T->chart=0x%8.8x\n",T->chart);fflush(stdout);
    printf("     T->chartR=0x%8.8x\n",T->chartR);fflush(stdout);
    printf("     T->chartCenter=0x%8.8x\n",T->chartCenter);fflush(stdout);
    printf("     T->value=0x%8.8x\n",T->value);fflush(stdout);
    printf("     T->leftB=0x%8.8x\n",T->leftB);fflush(stdout);
    printf("     T->rightB=0x%8.8x\n",T->rightB);fflush(stdout);
   }
#endif

  if(T==NULL)
   {
    return;
   }
  if(T->chart!=NULL)free(T->chart);
  T->chart=NULL;
  if(T->chartR!=NULL)free(T->chartR);
  T->chartR=NULL;
  if(T->chartCenter!=NULL)
   {
    for(i=0;i<T->nCharts;i++)
      if(T->chartCenter[i]!=NULL)free(T->chartCenter[i]);
    free(T->chartCenter);
    T->chartCenter=NULL;
   }
  if(T->value!=NULL)free(T->value);
  T->value=NULL;
  if(T->leftB!=NULL)free(T->leftB);
  T->leftB=NULL;
  if(T->rightB!=NULL)free(T->rightB);
  T->rightB=NULL;
  if(T->leftT!=NULL)MFFreeBinaryTreeLeaf(T->leftT,e);
  T->leftT=NULL;
  if(T->rightT!=NULL)MFFreeBinaryTreeLeaf(T->rightT,e);
  T->rightT=NULL;
  free(T);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("done %s, T=0x%8.8x\n",RoutineName,T);fflush(stdout);}
#endif

  return;
 }

void MFBinaryTreeAddChart(MFBinaryTree TRoot,int chartNo,double *center, double R, MFErrorHandler e)
 {
  static char RoutineName[]={"MFBinaryTreeAddChart"};
  MFBinaryTreeLeaf T;
  int k;
  int verbose=0;

  nInTree++;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s\n",RoutineName);fflush(stdout);}
#endif

  T=TRoot->root;
  k=TRoot->k;
  MFBinaryTreeLeafAddChart(T,k,chartNo,R,center,e);

  return;
 }

void MFBinaryTreeLeafAddChart(MFBinaryTreeLeaf T,int k, int chart, double R, double *center, MFErrorHandler e)
 {
  static char RoutineName[]={"MFBinaryTreeLeafAddChart"};
  int i,n;
  double t;
  int j,it;
  double *pt;
  double maxwidth,nextmaxwidth;
  double x0,x1;
  int nextdir;
  int verbose=0;

  t=center[T->d];

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("%s, %d, center (",RoutineName,T->sequenceNumber);
    for(i=0;i<k;i++){if(i>0)printf(",");printf("%lf",center[i]);}
    printf(") radius %lf, direction %d, t=%lf\n",R,T->d,t);fflush(stdout);
   }
#endif

  if(T->leftT==NULL && T->rightT==NULL &&
     T->nCharts<MFBINARYTREEMAXNUMBEROFCHARTS)
   {

#ifdef MFALLOWVERBOSE
    if(verbose){printf("  thisTree leaf (%d) is not a branch pt, and has space (%d charts/%d)\n",T->sequenceNumber,T->nCharts,MFBINARYTREEMAXNUMBEROFCHARTS);fflush(stdout);}
#endif

    if(T->nCharts>0)
     {
      for(i=0;i<T->nCharts;i++)if(T->chart[i]==chart){printf("Error, list already contains chart %d\n",i);fflush(stdout);}
      n=0;
      while(n<T->nCharts && T->value[n]<t)n++;
      for(i=T->nCharts;i>n;i--)
       {
        T->chart[i]=T->chart[i-1];
        T->chartR[i]=T->chartR[i-1];
        T->chartCenter[i]=T->chartCenter[i-1];
        T->value[i]=T->value[i-1];
       }
      T->chart[n]=chart;
      T->chartR[n]=R;
      T->chartCenter[n]=(double*)malloc(k*sizeof(double));

#ifndef MFNOSAFETYNET
      if(T->chartCenter[n]==NULL)
       {
        sprintf(MFBinaryTreeErrorMsg,"Out of memory, trying to allocate %d bytes",k*sizeof(double));
        MFSetError(e,12,RoutineName,MFBinaryTreeErrorMsg,__LINE__,__FILE__);
        MFErrorHandlerOutOfMemory(e);
        return;
       }
#endif

      for(i=0;i<k;i++)(T->chartCenter[n])[i]=center[i];
      T->value[n]=t;
     }else{
      T->chart[0]=chart;
      T->chartR[0]=R;
      T->chartCenter[0]=(double*)malloc(k*sizeof(double));

#ifndef MFNOSAFETYNET
      if(T->chartCenter[0]==NULL)
       {
        sprintf(MFBinaryTreeErrorMsg,"Out of memory, trying to allocate %d bytes",k*sizeof(double));
        MFSetError(e,12,RoutineName,MFBinaryTreeErrorMsg,__LINE__,__FILE__);
        MFErrorHandlerOutOfMemory(e);
        return;
       }
#endif

      for(i=0;i<k;i++)(T->chartCenter[0])[i]=center[i];
      T->value[0]=t;
     }
    T->nCharts++;
    for(i=0;i<k;i++)
     {
      if(T->nCharts==1)
       {
        T->leftB[i]=center[i]-R;
        T->rightB[i]=center[i]+R;
       }else{
        t=center[i]-R;
        if(t<T->leftB[i])T->leftB[i]=t;
        t=center[i]+R;
        if(t>T->rightB[i])T->rightB[i]=t;
       }
     }
   }else if(T->leftT==NULL && T->rightT==NULL)
   {

#ifdef MFALLOWVERBOSE
    if(verbose){printf("  thisTree leaf (%d) is not a branch pt, and has no space left (%d charts/%d)\n",T->sequenceNumber,T->nCharts,MFBINARYTREEMAXNUMBEROFCHARTS);fflush(stdout);}
#endif

/* recompute d and values */

#ifdef MFALLOWVERBOSE
     if(verbose){printf("  recomputing split and direction. k=%d\n",k);fflush(stdout);}
#endif

     maxwidth=-1.;
     T->split=0.;
     nextmaxwidth=0.;
     nextdir=0;
     for(i=0;i<k;i++)
      {
       x0=center[i];
       x1=center[i];
       for(j=0;j<T->nCharts;j++)
        {
         if((T->chartCenter[j])[i]>x1)x1=(T->chartCenter[j])[i];
         if((T->chartCenter[j])[i]<x0)x0=(T->chartCenter[j])[i];
        }

#ifdef MFALLOWVERBOSE
       if(verbose){printf("      direction %d, range is [%lf,%lf]\n",i,x0,x1);fflush(stdout);}
#endif

       if(x1-x0>maxwidth)
        {
         maxwidth=x1-x0;
         T->split=(x1+x0)/2;
         T->d=i;
        }else if(x1-x0>nextmaxwidth)
        {
         nextmaxwidth=x1-x0;
         nextdir=i;
        }
      }

#ifdef MFALLOWVERBOSE
    if(verbose){printf("      direction %d, range is [%lf,%lf]\n",i,x0,x1);fflush(stdout);}
#endif

/* Recompute Values */

    for(i=0;i<T->nCharts;i++)
      T->value[i]=(T->chartCenter[i])[T->d];

    T->leftT=MFCreateBinaryTreeLeaf(k,nextdir,e);
    (T->leftT)->parent=T;
    T->rightT=MFCreateBinaryTreeLeaf(k,nextdir,e);
    (T->rightT)->parent=T;

#ifdef MFALLOWVERBOSE
    if(verbose)
     {
      printf("  thisTree leaf (%d) is not a branch pt, and has no space\n",T->sequenceNumber);
      printf("     values [");
      for(i=0;i<T->nCharts;i++)
       {
        if(i>0)printf(",");
        printf("%lf",T->value[i]);
       }
      printf("], new value=%lf\n",center[T->d]);
      printf("     split is at %lf\n",T->split);
      fflush(stdout);
     }
#endif

    for(i=0;i<T->nCharts;i++)
     {
      if(T->value[i]<T->split)
       {

#ifdef MFALLOWVERBOSE
        if(verbose){printf("  move chart %d into left (%d), split=%lf, value=%lf\n",i,(T->leftT)->sequenceNumber,T->split,t);fflush(stdout);}
#endif

        MFBinaryTreeLeafAddChart(T->leftT,k,T->chart[i],T->chartR[i],T->chartCenter[i],e);
       }else if(T->value[i]>T->split){

#ifdef MFALLOWVERBOSE
        if(verbose){printf("  move chart %d into right (%d), split=%lf, value=%lf\n",i,(T->rightT)->sequenceNumber,T->split,t);fflush(stdout);}
#endif

        MFBinaryTreeLeafAddChart(T->rightT,k,T->chart[i],T->chartR[i],T->chartCenter[i],e);
       }else{
        if((T->leftT)->nCharts<(T->rightT)->nCharts)
         {

#ifdef MFALLOWVERBOSE
          if(verbose){printf("  put chart %d in left (%d), (same split value, left has fewer points) split=%lf, value=%lf\n",i,(T->leftT)->sequenceNumber,T->split,t);fflush(stdout);}
#endif

          MFBinaryTreeLeafAddChart(T->leftT,k,T->chart[i],T->chartR[i],T->chartCenter[i],e);
         }else{

#ifdef MFALLOWVERBOSE
          if(verbose){printf("  put chart %d in right (%d), (same split value, right has fewer points) split=%lf, value=%lf\n",i,(T->rightT)->sequenceNumber,T->split,t);fflush(stdout);}
#endif

          MFBinaryTreeLeafAddChart(T->rightT,k,T->chart[i],T->chartR[i],T->chartCenter[i],e);
         }
       }
     }

    t=center[T->d];
    if(t<T->split) /* if == T->split put in least full */
     {

#ifdef MFALLOWVERBOSE
      if(verbose){printf("  put new in left (%d) which has %d, split=%lf, value=%lf\n",(T->leftT)->sequenceNumber,(T->leftT)->nCharts,T->split,t);fflush(stdout);}
#endif

      MFBinaryTreeLeafAddChart(T->leftT,k,chart,R,center,e);
     }else if(t>T->split)
     {

#ifdef MFALLOWVERBOSE
      if(verbose){printf("  put new in right (%d), which has %d, split=%lf, value=%lf\n",(T->rightT)->sequenceNumber,(T->rightT)->nCharts,T->split,t);fflush(stdout);}
#endif

      MFBinaryTreeLeafAddChart(T->rightT,k,chart,R,center,e);
     }else{

#ifdef MFALLOWVERBOSE
      if(verbose){printf("  left has %d, right has %d\n",(T->leftT)->nCharts,(T->rightT)->nCharts);fflush(stdout);}
#endif

      if((T->leftT)->nCharts<(T->rightT)->nCharts)
       {

#ifdef MFALLOWVERBOSE
        if(verbose){printf("  put new in left (%d), (same split value, left has fewer points) split=%lf, value=%lf\n",(T->leftT)->sequenceNumber,T->split,t);fflush(stdout);}
#endif

        MFBinaryTreeLeafAddChart(T->leftT,k,chart,R,center,e);
       }else{

#ifdef MFALLOWVERBOSE
        if(verbose){printf("  put new in right (%d), (same split value, right has fewer points) split=%lf, value=%lf\n",(T->rightT)->sequenceNumber,T->split,t);fflush(stdout);}
#endif

        MFBinaryTreeLeafAddChart(T->rightT,k,chart,R,center,e);
       }
     }

    free(T->chart);
    T->chart=NULL;
    for(i=0;i<T->nCharts;i++)
     {
      if(T->chartCenter[i]!=NULL)free(T->chartCenter[i]);
     }
    free(T->chartCenter);
    T->chartCenter=NULL;
    free(T->chartR);
    T->chartR=NULL;
    T->nCharts=0;
   }else{

#ifdef MFALLOWVERBOSE
    if(verbose){printf("  thisTree leaf (%d) is a branch pt\n",T->sequenceNumber);fflush(stdout);}
#endif

    t=center[T->d];
    if(t<T->split) /* if == T->split put in least full */
     {

#ifdef MFALLOWVERBOSE
      if(verbose){printf("  put new in left (%d), split=%lf, value=%lf\n",(T->leftT)->sequenceNumber,T->split,t);fflush(stdout);}
#endif

      MFBinaryTreeLeafAddChart(T->leftT,k,chart,R,center,e);
     }else if(t>T->split)
     {

#ifdef MFALLOWVERBOSE
      if(verbose){printf("  put new in right (%d), split=%lf, value=%lf\n",(T->rightT)->sequenceNumber,T->split,t);fflush(stdout);}
#endif

      MFBinaryTreeLeafAddChart(T->rightT,k,chart,R,center,e);
     }else{
      if((T->rightT)->nCharts<(T->rightT)->nCharts)
       {

#ifdef MFALLOWVERBOSE
        if(verbose){printf("  put new in left (%d), (fewer points) split=%lf, value=%lf\n",(T->leftT)->sequenceNumber,T->split,t);fflush(stdout);}
#endif

        MFBinaryTreeLeafAddChart(T->leftT,k,chart,R,center,e);
       }else{

#ifdef MFALLOWVERBOSE
        if(verbose){printf("  put new in right (%d), (fewer points) split=%lf, value=%lf\n",(T->rightT)->sequenceNumber,T->split,t);fflush(stdout);}
#endif

        MFBinaryTreeLeafAddChart(T->rightT,k,chart,R,center,e);
       }
     }
   }

  for(i=0;i<k;i++)
   {
    if(T->leftT!=NULL)
     {
      if(T->leftB[i]>(T->leftT)->leftB[i])T->leftB[i]=(T->leftT)->leftB[i];
      if(T->rightB[i]<(T->leftT)->rightB[i])T->rightB[i]=(T->leftT)->rightB[i];
     }
    if(T->rightT!=NULL)
     {
      if((T->rightT)->leftB[i]<T->leftB[i])T->leftB[i]=(T->rightT)->leftB[i];
      if((T->rightT)->rightB[i]>T->rightB[i])T->rightB[i]=(T->rightT)->rightB[i];
     }
   }

  return;
 }

void MFBinaryTreeGetIntersectingCharts(MFBinaryTreeLeaf T,int k,int chart, double *center, double R,int *n,int **list, MFErrorHandler e)
 {
  static char RoutineName[]={"MFBinaryTreeGetIntersectingCharts"};
  int i,j;
  double t;
  int notin;
  int verbose=0;

  if(T==NULL)return;

#ifdef MFNOCONFIDENCE
  if(k<0)
   {
    sprintf(MFBinaryTreeErrorMsg,"k (2nd arg), is negative.");
    MFSetError(e,12,RoutineName,MFBinaryTreeErrorMsg,__LINE__,__FILE__);
    return;
   }
  if(R<0)
   {
    sprintf(MFBinaryTreeErrorMsg,"R (5th arg), is negative.");
    MFSetError(e,12,RoutineName,MFBinaryTreeErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  if(T->leftT==NULL && T->nCharts==0)return;

/* Check Bounding box */

  notin=0;

#ifdef MFALLOWVERBOSE
  if(verbose)printf("   leaf %d bb=\n",T->sequenceNumber);
#endif

#ifdef MFNOCONFIDENCE
  if(T->leftB==NULL)
   {
    sprintf(MFBinaryTreeErrorMsg,"Left subtree is NULL");
    MFSetError(e,12,RoutineName,MFBinaryTreeErrorMsg,__LINE__,__FILE__);
    return;
   }
  if(T->rightB==NULL)
   {
    sprintf(MFBinaryTreeErrorMsg,"Right subtree is NULL");
    MFSetError(e,12,RoutineName,MFBinaryTreeErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  for(i=0;i<k;i++)
   {
    t=center[i];

#ifdef MFALLOWVERBOSE
    if(verbose)printf("   %d BB %lf,%lf Disk %lf,%lf\n",i,T->leftB[i],T->rightB[i],t-R,t+R);
#endif

    if(t+R<T->leftB[i] || t-R>T->rightB[i])notin=1;
   }

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    if(notin)printf("     not in\n");
      else printf("     in\n");
   }
#endif
  if(notin)return;

/* Chart Intersects Bounding box, add charts to list, or check sub-trees */

  if(T->leftT==NULL)
   {
    *list=(int*)realloc((void*)(*list),(*n+T->nCharts)*sizeof(int));

#ifndef MFNOSAFETYNET
    if(*list==NULL)
     {
      sprintf(MFBinaryTreeErrorMsg,"Out of memory, trying to reallocate %d bytes",(*n+T->nCharts)*sizeof(int));
      MFSetError(e,12,RoutineName,MFBinaryTreeErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

    for(i=0;i<T->nCharts;i++)
     {
      if(chart>-1 && T->chart[i]!=chart)
       {
/* !!! */
        notin=0;
        for(j=0;j<k;j++)
         {
          t=center[j];
          if(t+R<(T->chartCenter[i])[j]-T->chartR[i] || t-R>(T->chartCenter[i])[j]+T->chartR[i])notin=1;
         }
        if(!notin)
/* !!! */
         {
          (*list)[*n]=T->chart[i];
          (*n)++;
         }
       }
     }
   }else{
    MFBinaryTreeGetIntersectingCharts(T->leftT,k,chart,center,R,n,list,e);
    MFBinaryTreeGetIntersectingCharts(T->rightT,k,chart,center,R,n,list,e);
   }

  return;
 }

MFListOfCharts MFCreateListOfIntersectingCharts(MFBinaryTree TRoot,int chart,double *center, double R, MFErrorHandler e)
 {
  static char RoutineName[]={"MFCreateListOfIntersectingCharts"};
  int i,n,m;
  int *charts;
  MFListOfCharts result;
  int verbose=0;

/* Start at top of tree. */
/*   if n!=0, test each member */
/*     else test bounding box of each subtree */

  n=0;
  charts=(int*)malloc(10*sizeof(int));

#ifndef MFNOSAFETYNET
  if(charts==NULL)
   {
    sprintf(MFBinaryTreeErrorMsg,"Out of memory, trying to allocate %d bytes",10*sizeof(int));
    MFSetError(e,12,RoutineName,MFBinaryTreeErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  MFBinaryTreeGetIntersectingCharts(TRoot->root,TRoot->k,chart,center,R,&n,&charts,e);
  result=MFCreateListOfCharts(n,charts,e);

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("There are %d intersecting charts\n [",n);
    for(i=0;i<n;i++)
     {
      if(i>0)printf(",");
      if(i%30==29)printf("\n  ");
      printf("%d",charts[i]);
     }
    printf("]\n");
    fflush(stdout);
   }
#endif

  return result;
 }

void MFPrintBinaryTreeLeaf(FILE *fid, int k, MFBinaryTreeLeaf T, MFErrorHandler e)
 {
  int i,j;

  if(T==NULL)
   {
    fprintf(fid,"     NULL\n");
    return;
   }
  fprintf(fid,"     sequenceNumber %d\n",T->sequenceNumber);
  fprintf(fid,"     Split is direction %d, threshold %lf\n",T->d,T->split);
  fprintf(fid,"     number of charts %d\n",T->nCharts);fflush(stdout);
  for(i=0;i<T->nCharts;i++)
   {
    fprintf(fid,"chart %d",T->chart[i]);
    fprintf(fid," Radius %lf, center  (",T->chartR[i]);
    for(j=0;j<k;j++){if(j>0)fprintf(fid,",");fprintf(fid,"%lf",(T->chartCenter[i])[j]);}
    fprintf(fid,") value %lf\n",T->value[i]);
   }
  fprintf(fid,"     Bounding box (");
  for(i=0;i<k;i++)
   {
    if(i>0)fprintf(fid,",");
    fprintf(fid,"%lf",T->leftB[i]);
   }
  fprintf(fid,")->(");
  for(i=0;i<k;i++)
   {
    if(i>0)fprintf(fid,",");
    fprintf(fid,"%lf",T->rightB[i]);
   }
  fprintf(fid,")\n");
  fprintf(fid,"LeftLeaf\n");
  MFPrintBinaryTreeLeaf(fid,k,T->leftT,e);
  fprintf(fid,"RightLeaf\n");
  MFPrintBinaryTreeLeaf(fid,k,T->rightT,e);

  return;
 }

void MFPrintBinaryTree(FILE *fid, MFBinaryTree T, MFErrorHandler e)
 {

  fprintf(fid,"Binary Tree, dimension %d\n",T->k);
  fprintf(fid,"Root\n");
  MFPrintBinaryTreeLeaf(fid,T->k,T->root,e);

  return;
 }

void MFWriteBinaryTree(FILE *fid,MFBinaryTree u, MFErrorHandler e)
 {
  static char RoutineName[]={"MFWriteBinaryTree"};
  int i;

  fprintf(fid,"%s\n","BinaryTree");
  return;
 }

MFBinaryTree MFReadBinaryTree(FILE *fid, MFErrorHandler e)
 {
  static char RoutineName[]={"MFReadBinaryTree"};
  int i;
  MFBinaryTree B;
  char tag[100]="";

  fscanf(fid,"%s\n",tag);
  if(strcmp(tag,"BinaryTree"))
   {
    sprintf(MFBinaryTreeErrorMsg,"Next Object is not a BinaryTree! (%s)\n",RoutineName,tag);
    MFSetError(e,12,RoutineName,MFBinaryTreeErrorMsg,__LINE__,__FILE__);
    return NULL;
   }
  B=NULL;
  return B;
 }

void MFRecomputeBoundingBoxes(MFBinaryTree T, int chart, double *center, double R, MFErrorHandler e)
 {
  MFBinaryTreeLeaf leaf;
  int ic;
  int i,j;

  leaf=T->root;
  while(leaf->nCharts==0)
   {
    if ( center[leaf->d]<leaf->split)
      leaf=leaf->leftT;
     else
      leaf=leaf->rightT;
   }

  for(i=0;i<leaf->nCharts;i++)
    if(leaf->chart[i]==chart)leaf->chartR[ic]=R;

  for(i=0;i<T->k;i++)
   {
    leaf->leftB[i]=center[i]-leaf->chartR[0];
    leaf->rightB[i]=center[i]+leaf->chartR[0];;
    for(j=1;j<leaf->nCharts;j++)
     {
      if(center[i]-leaf->chartR[j]<leaf->leftB[i])leaf->leftB[i]=center[i]-leaf->chartR[j];
      if(center[i]+leaf->chartR[j]>leaf->rightB[i])leaf->rightB[i]=center[i]+leaf->chartR[j];
     }
   }
  leaf=leaf->parent;

  while(leaf!=NULL)
   {
    for(i=0;i<T->k;i++)
     {
      if(leaf->leftB!=NULL)
       {
        leaf->leftB[i]=(leaf->leftT)->leftB[i];
        if(leaf->rightB!=NULL&&leaf->leftB[i]>(leaf->leftT)->leftB[i])leaf->leftB[i]=(leaf->rightT)->leftB[i];
       }else if(leaf->rightB!=NULL)
        leaf->leftB[i]=(leaf->rightT)->leftB[i];

      if(leaf->rightB!=NULL)
       {
        leaf->rightB[i]=(leaf->rightT)->rightB[i];
        if(leaf->leftB!=NULL&&leaf->rightB[i]>(leaf->rightT)->rightB[i])leaf->rightB[i]=(leaf->rightT)->rightB[i];
       }else if(leaf->leftB!=NULL)
        leaf->rightB[i]=(leaf->rightT)->rightB[i];
     }
    leaf=leaf->parent;
   }

 }

MFListOfCharts MFCreateListOfNearbyCharts(MFBinaryTree TRoot,double *center,double R, MFErrorHandler e)
 {
  static char RoutineName[]={"MFCreateListOfNearbyCharts"};
  int i,n,m;
  int *charts;
  MFListOfCharts result;
  int verbose=0;

/* Start at top of tree. */
/*   if n!=0, test each member */
/*     else test bounding box of each subtree */

  n=0;
  charts=(int*)malloc(10*sizeof(int));

#ifndef MFNOSAFETYNET
  if(charts==NULL)
   {
    sprintf(MFBinaryTreeErrorMsg,"Out of memory, trying to allocate %d bytes",10*sizeof(int));
    MFSetError(e,12,RoutineName,MFBinaryTreeErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  MFBinaryTreeGetNeighboringCharts(TRoot->root,TRoot->k,center,R,&n,&charts,e);
  result=MFCreateListOfCharts(n,charts,e);

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("There are %d intersecting charts\n [",n);
    for(i=0;i<n;i++)
     {
      if(i>0)printf(",");
      if(i%30==29)printf("\n  ");
      printf("%d",charts[i]);
     }
    printf("]\n");
    fflush(stdout);
   }
#endif

  return result;
 }

void MFBinaryTreeGetNeighboringCharts(MFBinaryTreeLeaf T,int k, double *center, double R,int *n,int **list, MFErrorHandler e)
 {
  static char RoutineName[]={"MFBinaryTreeGetIntersectingCharts"};
  int i,j;
  double t;
  int notin;
  int verbose=0;

  if(T==NULL)return;

#ifdef MFNOCONFIDENCE
  if(k<0)
   {
    sprintf(MFBinaryTreeErrorMsg,"k (2nd arg), is negative.");
    MFSetError(e,12,RoutineName,MFBinaryTreeErrorMsg,__LINE__,__FILE__);
    return;
   }
  if(R<0)
   {
    sprintf(MFBinaryTreeErrorMsg,"R (5th arg), is negative.");
    MFSetError(e,12,RoutineName,MFBinaryTreeErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  if(T->leftT==NULL && T->nCharts==0)return;

/* Check Bounding box */

  notin=0;

#ifdef MFALLOWVERBOSE
  if(verbose)printf("   leaf %d bb=\n",T->sequenceNumber);
#endif

#ifdef MFNOCONFIDENCE
  if(T->leftB==NULL)
   {
    sprintf(MFBinaryTreeErrorMsg,"Left subtree is NULL");
    MFSetError(e,12,RoutineName,MFBinaryTreeErrorMsg,__LINE__,__FILE__);
    return;
   }
  if(T->rightB==NULL)
   {
    sprintf(MFBinaryTreeErrorMsg,"Right subtree is NULL");
    MFSetError(e,12,RoutineName,MFBinaryTreeErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  for(i=0;i<k;i++)
   {
    t=center[i];

#ifdef MFALLOWVERBOSE
    if(verbose)printf("   %d BB %lf,%lf Disk %lf,%lf\n",i,T->leftB[i],T->rightB[i],t-R,t+R);
#endif

    if(t+R<T->leftB[i] || t-R>T->rightB[i])notin=1;
   }

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    if(notin)printf("     not in\n");
     else printf("     in\n");
   }
#endif

  if(notin)return;

/* Chart Intersects Bounding box, add charts to list, or check sub-trees */

  if(T->leftT==NULL)
   {
    *list=(int*)realloc((void*)(*list),(*n+T->nCharts)*sizeof(int));

#ifndef MFNOSAFETYNET
    if(*list==NULL)
     {
      sprintf(MFBinaryTreeErrorMsg,"Out of memory, trying to reallocate %d bytes",(*n+T->nCharts)*sizeof(int));
      MFSetError(e,12,RoutineName,MFBinaryTreeErrorMsg,__LINE__,__FILE__);
      return;
     }
#endif

    for(i=0;i<T->nCharts;i++)
     {
      notin=0;
      for(j=0;j<k;j++)
       {
        t=center[j];
        if(t+R<(T->chartCenter[i])[j]-T->chartR[i] || t-R>(T->chartCenter[i])[j]+T->chartR[i])notin=1;
       }
      if(!notin)
       {
        (*list)[*n]=T->chart[i];
        (*n)++;
       }
     }
   }else{
    MFBinaryTreeGetNeighboringCharts(T->leftT,k,center,R,n,list,e);
    MFBinaryTreeGetNeighboringCharts(T->rightT,k,center,R,n,list,e);
   }

  return;
 }

#ifdef __cplusplus
}
#endif
