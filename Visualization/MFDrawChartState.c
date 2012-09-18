/*
 *
 *  PROGRAM NAME:  multifario
 *
 *  (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
 *  CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
 *
 *  Please refer to the LICENSE file in the top directory
 *
 *      author: Mike Henderson mhender@watson.ibm.com
 */

static char *id="@(#) $Id: MFDrawChartState.c,v 1.3 2007/02/13 01:22:30 mhender Exp $";

static char MFChartStateErrorMsg[256]="";

#include <MFDraw.h>
#include <MFDrawChartState.h>
#include <MFChart.h>
#include <MFPolytope.h>
#include <MFErrorHandler.h>
#include <sh.h>
#include <math.h>
#include <MFErrorHandler.h>
#include <stdlib.h>

#ifdef __cplusplus
 extern "C" {
#endif

struct MFChartStateSt
 {
  int r[128],g[128],b[128];
  MFChart chart;
  int nP;
  int mP;
  int *nV;
  float **pX;
  float **pY;
  float **pZ;
  int nNP;
  int mNP;
  int *nNV;
  float **pNX;
  float **pNY;
  float **pNZ;
  float **pNNrm;
  int *pc;
  int nL;
  float *lX;
  float *lY;
  float *lZ;
  int nPlns;
  int *gInd;
  int nVert;
 };

int MFDrawNCharts=0;
int MFDrawMCharts=0;
MFChartState *MFDrawChartList=NULL;

void MFFreeChartState(MFChartState,MFErrorHandler);

MFChartState MFCreateChartState(MFChart chart, MFErrorHandler e)
 {
  static char RoutineName[]={"MFCreateChartState"};
  MFChartState thisState;
  MFPolytope P;
  int color;
  int zero=0;
  int half=128;
  int full=255;
  int i;

  thisState=(MFChartState)malloc(sizeof(struct MFChartStateSt));

#ifndef MFNOSAFETYNET
  if(thisState==NULL)
   {
    sprintf(MFChartStateErrorMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct MFChartStateSt));
    MFSetError(e,12,RoutineName,MFChartStateErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  thisState->chart=chart;
  MFRefChart(chart,e);

  thisState->nP=0;
  thisState->mP=0;
  thisState->nV=NULL;
  thisState->pX=NULL;
  thisState->pY=NULL;
  thisState->pZ=NULL;

  thisState->nNP=0;
  thisState->mNP=0;
  thisState->nNV=NULL;
  thisState->pNX=NULL;
  thisState->pNY=NULL;
  thisState->pNZ=NULL;
  thisState->pNNrm=NULL;
  thisState->pc=NULL;

  thisState->nL=0;
  thisState->lX=NULL;
  thisState->lY=NULL;
  thisState->lZ=NULL;

  thisState->r[0]=half;
  thisState->g[0]=half;
  thisState->b[0]=half;

  P=MFChartPolytope(chart,e);
  thisState->nVert=MFPolytopeNumberOfVertices(P,e);
  thisState->nPlns=MFPolytopeNumberOfFaces(P,e);
  thisState->gInd=(int*)malloc(thisState->nPlns*sizeof(int));

#ifndef MFNOSAFETYNET
  if(thisState->nPlns>0&&thisState->gInd==NULL)
   {
    sprintf(MFChartStateErrorMsg,"Out of memory, trying to allocate %d bytes",thisState->nPlns*sizeof(int));
    MFSetError(e,12,RoutineName,MFChartStateErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  for(i=0;i<thisState->nPlns;i++)
   thisState->gInd[i]=MFPolytopeFaceIndex(P,i,e);

  return(thisState);
 }

MFChartState MFAddChartState(MFChart chart, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAddChartState"};
  if(MFDrawNCharts==MFDrawMCharts)
   {
    MFDrawMCharts+=10;
    MFDrawChartList=(MFChartState*)realloc(MFDrawChartList,MFDrawMCharts*sizeof(MFChartState));

#ifndef MFNOSAFETYNET
    if(MFDrawChartList==NULL)
     {
      sprintf(MFChartStateErrorMsg,"Out of memory, trying to reallocate %d bytes",MFDrawMCharts*sizeof(MFChartState));
      MFSetError(e,12,RoutineName,MFChartStateErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return NULL;
     }
#endif

   }
  MFDrawChartList[MFDrawNCharts]=MFCreateChartState(chart,e);
  MFDrawNCharts++;

  return MFDrawChartList[MFDrawNCharts-1];
 }

int MFNChartStates(MFErrorHandler e)
 {
  static char RoutineName[]={"MFNChartStates"};

  return MFDrawNCharts;
 }

void MFClearChartState(MFChartState thisState, MFErrorHandler e)
 {
  static char RoutineName[]={"MFClearChartState"};
  MFPolytope P;
  int i;

  if(thisState->nP>0)
   {
    if(thisState->pX!=NULL)
     {
      for(i=0;i<thisState->nP;i++)
       if((thisState->pX)[i]!=NULL)free((thisState->pX)[i]);
      free(thisState->pX);thisState->pX=NULL;
     }
    if(thisState->pY!=NULL)
     {
      for(i=0;i<thisState->nP;i++)
       if((thisState->pY)[i]!=NULL)free((thisState->pY)[i]);
      free(thisState->pY);thisState->pY=NULL;
     }
    if(thisState->pZ!=NULL)
     {
      for(i=0;i<thisState->nP;i++)
       if((thisState->pZ)[i]!=NULL)free((thisState->pZ)[i]);
      free(thisState->pZ);thisState->pZ=NULL;
     }
    if(thisState->nV!=NULL){free(thisState->nV);thisState->nV=(int*)NULL;}
   }
  thisState->nP=0;
  thisState->mP=0;

  if(thisState->nNP>0)
   {
    if(thisState->pNX!=NULL)
     {
      for(i=0;i<thisState->nP;i++)
       if((thisState->pNX)[i]!=NULL)free((thisState->pNX)[i]);
      free(thisState->pNX);thisState->pNX=NULL;
     }
    if(thisState->pNY!=NULL)
     {
      for(i=0;i<thisState->nP;i++)
       if((thisState->pNY)[i]!=NULL)free((thisState->pNY)[i]);
      free(thisState->pNY);thisState->pNY=NULL;
     }
    if(thisState->pNZ!=NULL)
     {
      for(i=0;i<thisState->nP;i++)
       if((thisState->pNZ)[i]!=NULL)free((thisState->pNZ)[i]);
      free(thisState->pNZ);thisState->pNZ=NULL;
     }
    if(thisState->pNNrm!=NULL)
     {
      for(i=0;i<thisState->nNP;i++)
       if((thisState->pNNrm)[i]!=NULL)free((thisState->pNNrm)[i]);
      free(thisState->pNNrm);thisState->pNNrm=NULL;
     }
    if(thisState->nNV!=NULL){free(thisState->nNV);thisState->nNV=(int*)NULL;}
    if(thisState->pc!=NULL){free(thisState->pc);thisState->pc=(int*)NULL;}
   }
  thisState->nNP=0;
  thisState->mNP=0;

  thisState->nL=0;
  if(thisState->lX!=NULL){free(thisState->lX);thisState->lX=(float*)NULL;}
  if(thisState->lY!=NULL){free(thisState->lY);thisState->lY=(float*)NULL;}
  if(thisState->lZ!=NULL){free(thisState->lZ);thisState->lZ=(float*)NULL;}


  P=MFChartPolytope(thisState->chart,e);
  thisState->nVert=MFPolytopeNumberOfVertices(P,e);
  thisState->nPlns=MFPolytopeNumberOfFaces(P,e);
  if(thisState->gInd!=NULL)free(thisState->gInd);
  thisState->gInd=(int*)malloc(thisState->nPlns*sizeof(int));

#ifndef MFNOSAFETYNET
  if(thisState->gInd==NULL&&thisState->nPlns>0)
   {
    sprintf(MFChartStateErrorMsg,"Out of memory, trying to allocate %d bytes",thisState->nPlns*sizeof(int));
    MFSetError(e,12,RoutineName,MFChartStateErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  for(i=0;i<thisState->nPlns;i++)
   thisState->gInd[i]=MFPolytopeFaceIndex(P,i,e);

  return;
 }

void MFFreeChartState(MFChartState thisState, MFErrorHandler e)
 {
  static char RoutineName[]={"MFFreeChartState"};
  MFPolytope P;
  int i;

  if(thisState->nP>0)
   {
    if(thisState->pX!=NULL)
     {
      for(i=0;i<thisState->nP;i++)
       if((thisState->pX)[i]!=NULL)free((thisState->pX)[i]);
      free(thisState->pX);thisState->pX=NULL;
     }
    if(thisState->pY!=NULL)
     {
      for(i=0;i<thisState->nP;i++)
       if((thisState->pY)[i]!=NULL)free((thisState->pY)[i]);
      free(thisState->pY);thisState->pY=NULL;
     }
    if(thisState->pZ!=NULL)
     {
      for(i=0;i<thisState->nP;i++)
       if((thisState->pZ)[i]!=NULL)free((thisState->pZ)[i]);
      free(thisState->pZ);thisState->pZ=NULL;
     }
    if(thisState->nV!=NULL){free(thisState->nV);thisState->nV=(int*)NULL;}
   }
  thisState->nP=0;
  thisState->mP=0;

  if(thisState->nNP>0)
   {
    if(thisState->pNX!=NULL)
     {
      for(i=0;i<thisState->nP;i++)
       if((thisState->pNX)[i]!=NULL)free((thisState->pNX)[i]);
      free(thisState->pNX);thisState->pNX=NULL;
     }
    if(thisState->pNY!=NULL)
     {
      for(i=0;i<thisState->nP;i++)
       if((thisState->pNY)[i]!=NULL)free((thisState->pNY)[i]);
      free(thisState->pNY);thisState->pNY=NULL;
     }
    if(thisState->pNZ!=NULL)
     {
      for(i=0;i<thisState->nP;i++)
       if((thisState->pNZ)[i]!=NULL)free((thisState->pNZ)[i]);
      free(thisState->pNZ);thisState->pNZ=NULL;
     }
    if(thisState->pNNrm!=NULL)
     {
      for(i=0;i<thisState->nNP;i++)
       if((thisState->pNNrm)[i]!=NULL)free((thisState->pNNrm)[i]);
      free(thisState->pNNrm);thisState->pNNrm=NULL;
     }
    if(thisState->nNV!=NULL){free(thisState->nNV);thisState->nNV=(int*)NULL;}
    if(thisState->pc!=NULL){free(thisState->pc);thisState->pc=(int*)NULL;}
   }
  thisState->nNP=0;
  thisState->mNP=0;

  thisState->nL=0;
  if(thisState->lX!=NULL){free(thisState->lX);thisState->lX=(float*)NULL;}
  if(thisState->lY!=NULL){free(thisState->lY);thisState->lY=(float*)NULL;}
  if(thisState->lZ!=NULL){free(thisState->lZ);thisState->lZ=(float*)NULL;}

  thisState->nVert=0;
  thisState->nPlns=0;
  if(thisState->gInd!=NULL)free(thisState->gInd);

  free(thisState);

  return;
 }

void MFSetChartStateColor(MFChartState thisState,int c,int r,int g,int b, MFErrorHandler e)
 {
  static char RoutineName[]={"MFSetChartStateColor"};

  thisState->r[c]=r;
  thisState->g[c]=g;
  thisState->b[c]=b;
  return;
 }

void MFGetChartStateColor(MFChartState thisState,int c,int *r,int *g,int *b, MFErrorHandler e)
 {
  static char RoutineName[]={"MFGetChartStateColor"};

  *r=thisState->r[c];
  *g=thisState->g[c];
  *b=thisState->b[c];
  return;
 }

void MFChartStateAddPolygon(MFChartState thisState,int n,float *x,float *y,float *z,int c, MFErrorHandler e)
 {
  static char RoutineName[]={"MFChartStateAddPolygon"};
  int i;

  for(i=0;i<n-1;i++)
   if(0 && sqrt((x[i+1]-x[i])*(x[i+1]-x[i])
          +(y[i+1]-y[i])*(y[i+1]-y[i])
          +(z[i+1]-z[i])*(z[i+1]-z[i]))>1.)return;

  if(thisState->nP+1>thisState->mP)
   {
    thisState->mP+=1000;
    thisState->nV=(int*)realloc((thisState->nV),thisState->mP*sizeof(int));

#ifndef MFNOSAFETYNET
    if(thisState->nV==NULL)
     {
      sprintf(MFChartStateErrorMsg,"Out of memory, trying to reallocate %d bytes",thisState->mP*sizeof(int));
      MFSetError(e,12,RoutineName,MFChartStateErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

    thisState->pX=(float**)realloc((thisState->pX),thisState->mP*sizeof(float*));

#ifndef MFNOSAFETYNET
    if(thisState->pX==NULL)
     {
      sprintf(MFChartStateErrorMsg,"Out of memory, trying to reallocate %d bytes",thisState->mP*sizeof(float*));
      MFSetError(e,12,RoutineName,MFChartStateErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

    thisState->pY=(float**)realloc((thisState->pY),thisState->mP*sizeof(float*));

#ifndef MFNOSAFETYNET
    if(thisState->pY==NULL)
     {
      sprintf(MFChartStateErrorMsg,"Out of memory, trying to reallocate %d bytes",thisState->mP*sizeof(float*));
      MFSetError(e,12,RoutineName,MFChartStateErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

    thisState->pZ=(float**)realloc((thisState->pZ),thisState->mP*sizeof(float*));

#ifndef MFNOSAFETYNET
    if(thisState->pZ==NULL)
     {
      sprintf(MFChartStateErrorMsg,"Out of memory, trying to reallocate %d bytes",thisState->mP*sizeof(float*));
      MFSetError(e,12,RoutineName,MFChartStateErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

    thisState->pc=(int*)malloc(thisState->mP*sizeof(int));

#ifndef MFNOSAFETYNET
    if(thisState->pc==NULL)
     {
      sprintf(MFChartStateErrorMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(int));
      MFSetError(e,12,RoutineName,MFChartStateErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

   }

  thisState->nV[thisState->nP]=n;
  thisState->pX[thisState->nP]=(float*)malloc(n*sizeof(float));

#ifndef MFNOSAFETYNET
  if(thisState->pX[thisState->nP]==NULL)
   {
    sprintf(MFChartStateErrorMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(float));
    MFSetError(e,12,RoutineName,MFChartStateErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  thisState->pY[thisState->nP]=(float*)malloc(n*sizeof(float));

#ifndef MFNOSAFETYNET
  if(thisState->pY[thisState->nP]==NULL)
   {
    sprintf(MFChartStateErrorMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(float));
    MFSetError(e,12,RoutineName,MFChartStateErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  thisState->pZ[thisState->nP]=(float*)malloc(n*sizeof(float));

#ifndef MFNOSAFETYNET
  if(thisState->pZ[thisState->nP]==NULL)
   {
    sprintf(MFChartStateErrorMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(float));
    MFSetError(e,12,RoutineName,MFChartStateErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  for(i=0;i<n;i++)
   {
    (thisState->pX[thisState->nP])[i]=x[i];
    (thisState->pY[thisState->nP])[i]=y[i];
    (thisState->pZ[thisState->nP])[i]=z[i];
   }
  thisState->pc[thisState->nP]=c;
  thisState->nP++;
  return;
 }

void MFChartStateAddPolygonWithNormal(MFChartState thisState,int n,float *x,float *y,float *z,float *nrm, int c, MFErrorHandler e)
 {
  static char RoutineName[]={"MFChartStateAddPolygonWithNormal"};
  int i;

  for(i=0;i<n-1;i++)
   if(0 && sqrt((x[i+1]-x[i])*(x[i+1]-x[i])
          +(y[i+1]-y[i])*(y[i+1]-y[i])
          +(z[i+1]-z[i])*(z[i+1]-z[i]))>1.)return;

  if(thisState->nNP+1>thisState->mNP)
   {
    thisState->mNP+=1000;
    thisState->nNV=(int*)realloc((thisState->nNV),thisState->mNP*sizeof(int));

#ifndef MFNOSAFETYNET
    if(thisState->nNV==NULL)
     {
      sprintf(MFChartStateErrorMsg,"Out of memory, trying to reallocate %d bytes",thisState->mNP*sizeof(int));
      MFSetError(e,12,RoutineName,MFChartStateErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

    thisState->pNX=(float**)realloc((thisState->pNX),thisState->mNP*sizeof(float*));

#ifndef MFNOSAFETYNET
    if(thisState->pNX==NULL)
     {
      sprintf(MFChartStateErrorMsg,"Out of memory, trying to reallocate %d bytes",thisState->mNP*sizeof(float*));
      MFSetError(e,12,RoutineName,MFChartStateErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

    thisState->pNY=(float**)realloc((thisState->pNY),thisState->mNP*sizeof(float*));

#ifndef MFNOSAFETYNET
    if(thisState->pNY==NULL)
     {
      sprintf(MFChartStateErrorMsg,"Out of memory, trying to reallocate %d bytes",thisState->mNP*sizeof(float*));
      MFSetError(e,12,RoutineName,MFChartStateErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

    thisState->pNZ=(float**)realloc((thisState->pNZ),thisState->mNP*sizeof(float*));

#ifndef MFNOSAFETYNET
    if(thisState->pNZ==NULL)
     {
      sprintf(MFChartStateErrorMsg,"Out of memory, trying to reallocate %d bytes",thisState->mNP*sizeof(float*));
      MFSetError(e,12,RoutineName,MFChartStateErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

    thisState->pNNrm=(float**)realloc((thisState->pNNrm),thisState->mNP*sizeof(float*));

#ifndef MFNOSAFETYNET
    if(thisState->pNNrm==NULL)
     {
      sprintf(MFChartStateErrorMsg,"Out of memory, trying to reallocate %d bytes",thisState->mNP*sizeof(float*));
      MFSetError(e,12,RoutineName,MFChartStateErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

    thisState->pc=(int*)realloc((thisState->pc),thisState->mNP*sizeof(int));

#ifndef MFNOSAFETYNET
    if(thisState->pc==NULL)
     {
      sprintf(MFChartStateErrorMsg,"Out of memory, trying to reallocate %d bytes",thisState->mNP*sizeof(int));
      MFSetError(e,12,RoutineName,MFChartStateErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

   }

  thisState->nNV[thisState->nNP]=n;
  thisState->pNX[thisState->nNP]=(float*)malloc(n*sizeof(float));

#ifndef MFNOSAFETYNET
  if(thisState->pNX[thisState->nNP]==NULL)
   {
    sprintf(MFChartStateErrorMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(float));
    MFSetError(e,12,RoutineName,MFChartStateErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  thisState->pNY[thisState->nNP]=(float*)malloc(n*sizeof(float));

#ifndef MFNOSAFETYNET
  if(thisState->pNY[thisState->nNP]==NULL)
   {
    sprintf(MFChartStateErrorMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(float));
    MFSetError(e,12,RoutineName,MFChartStateErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  thisState->pNZ[thisState->nNP]=(float*)malloc(n*sizeof(float));

#ifndef MFNOSAFETYNET
  if(thisState->pNZ[thisState->nNP]==NULL)
   {
    sprintf(MFChartStateErrorMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(float));
    MFSetError(e,12,RoutineName,MFChartStateErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  thisState->pNNrm[thisState->nNP]=(float*)malloc(3*n*sizeof(float));

#ifndef MFNOSAFETYNET
  if(thisState->pNNrm[thisState->nNP]==NULL)
   {
    sprintf(MFChartStateErrorMsg,"Out of memory, trying to allocate %d bytes",3*n*sizeof(float));
    MFSetError(e,12,RoutineName,MFChartStateErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  for(i=0;i<n;i++)
   {
    (thisState->pNX[thisState->nNP])[i]=x[i];
    (thisState->pNY[thisState->nNP])[i]=y[i];
    (thisState->pNZ[thisState->nNP])[i]=z[i];
    (thisState->pNNrm[thisState->nNP])[3*i]=nrm[3*i];
    (thisState->pNNrm[thisState->nNP])[3*i+1]=nrm[3*i+1];
    (thisState->pNNrm[thisState->nNP])[3*i+1]=nrm[3*i+2];
   }
  thisState->pc[thisState->nP]=c;
  thisState->nNP++;
  return;
 }

void MFChartStateAddLine(MFChartState thisState,float x0,float y0,float z0,float x1,float y1,float z1, MFErrorHandler e)
 {
  static char RoutineName[]={"MFChartStateAddLine"};
  if(0 && sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)+(z1-z0)*(z1-z0))>1.)return;

  thisState->lX=(float*)realloc((thisState->lX),2*(thisState->nL+1)*sizeof(float));

#ifndef MFNOSAFETYNET
  if(thisState->lX==NULL)
   {
    sprintf(MFChartStateErrorMsg,"Out of memory, trying to reallocate %d bytes",2*(thisState->nL+1)*sizeof(float));
    MFSetError(e,12,RoutineName,MFChartStateErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  thisState->lY=(float*)realloc((thisState->lY),2*(thisState->nL+1)*sizeof(float));

#ifndef MFNOSAFETYNET
  if(thisState->lY==NULL)
   {
    sprintf(MFChartStateErrorMsg,"Out of memory, trying to reallocate %d bytes",2*(thisState->nL+1)*sizeof(float));
    MFSetError(e,12,RoutineName,MFChartStateErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  thisState->lZ=(float*)realloc((thisState->lZ),2*(thisState->nL+1)*sizeof(float));

#ifndef MFNOSAFETYNET
  if(thisState->lZ==NULL)
   {
    sprintf(MFChartStateErrorMsg,"Out of memory, trying to reallocate %d bytes",2*(thisState->nL+1)*sizeof(float));
    MFSetError(e,12,RoutineName,MFChartStateErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif


  thisState->lX[2*thisState->nL]=x0;
  thisState->lY[2*thisState->nL]=y0;
  thisState->lZ[2*thisState->nL]=z0;
  thisState->lX[2*thisState->nL+1]=x1;
  thisState->lY[2*thisState->nL+1]=y1;
  thisState->lZ[2*thisState->nL+1]=z1;
  thisState->nL++;
  return;
 }

MFChartState MFGetChartState(int i, MFErrorHandler e)
 {
  static char RoutineName[]={"MFGetChartState"};
  return MFDrawChartList[i];
 }

int MFChartStateChanged(MFChartState thisState, MFErrorHandler e)
 {
  static char RoutineName[]={"MFChartStateChanged"};
  MFPolytope P;
  int i;

  P=MFChartPolytope(thisState->chart,e);

  if(MFPolytopeNumberOfVertices(P,e)!=thisState->nVert)return 1;
  if(MFPolytopeNumberOfFaces(P,e)!=thisState->nPlns)return 1;

  for(i=0;i<thisState->nPlns;i++)
    if(MFPolytopeFaceIndex(P,i,e)!=thisState->gInd[i])return 1;

  return 0;
 }

void MFDrawChartFromState(MFChartState thisState, MFErrorHandler e)
 {
  static char RoutineName[]={"MFDrawChartFromState"};
  int i,n;
  int full=255;
  int zero=0;

  for(i=0;i<thisState->nP;i++)
   {
    shtric(&(thisState->r[thisState->pc[i]]),&(thisState->g[thisState->pc[i]]),&(thisState->b[thisState->pc[i]]));
    shpg(thisState->nV+i,(thisState->pX)[i],(thisState->pY)[i],(thisState->pZ)[i]);
   }

  shpgc(&full,&full,&zero);
  shpec(&full,&full,&zero);
  shbgc(&full,&full,&zero);
  shbec(&full,&full,&zero);
  for(i=0;i<thisState->nNP;i++)
   shpg(thisState->nNV+i,(thisState->pNX)[i],(thisState->pNY)[i],(thisState->pNZ)[i]);
/* shpgnrm(thisState->nNV+i,(thisState->pNX)[i],(thisState->pNY)[i],(thisState->pNZ)[i],(thisState->pNNrm)[i]);*/

  shlinc(&zero,&zero,&zero);
  for(i=0;i<thisState->nL;i++)
   shline(thisState->lX+2*i,thisState->lY+2*i,thisState->lZ+2*i,thisState->lX+2*i+1,thisState->lY+2*i+1,thisState->lZ+2*i+1);

  return;
 }

void MFFreeDrawCharts(MFErrorHandler e)
 {
  static char RoutineName[]={"MFFreeDrawCharts"};
  int i;

  for(i=0;i<MFDrawNCharts;i++)MFFreeChartState(MFDrawChartList[i],e);
  if(MFDrawChartList!=NULL)free(MFDrawChartList);
  MFDrawNCharts=0;
  MFDrawMCharts=0;
  MFDrawChartList=NULL;

  return;
 }

#ifdef __cplusplus
}
#endif
