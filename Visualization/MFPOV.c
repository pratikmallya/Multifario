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

static char *id="@(#) $Id: MFPOV.c,v 1.8 2011/07/21 17:44:29 mhender Exp $";

static char MFPOVErrorMsg[256]="";

#include <stdlib.h>
#include <math.h>
#include <MFAtlas.h>
#include <MFChart.h>
#include <MFEnumDualPolytope.h>
#include <MFEnumPolytope.h>
#include <MFPrint.h>
#include <string.h>

#ifdef __cplusplus
 extern "C" {
#endif

void MFMeshToPOVFile(char*,int,int,double*,int*,int,int*,int*,int,int*,int**,int*,MFErrorHandler);
void MFMeshToDXFile2(char*,int,int,double*,int*,int,int*,int*,int,int*,int**,int*,MFErrorHandler);
void MFMeshToDXFile3(char*,int,int,double*,int*,int,int*,int*,int,int*,int**,int*,MFErrorHandler);
void ThreeDPolytopes(MFAtlas,int,int*,int*,int*,double**,int**,int*,int*,int**,int**,int*,int*,int**,int***,int**,MFErrorHandler);
void ChartCenters(MFAtlas,int,int*,int*,int*,double**,int**,int*,int*,int**,int**,int*,int*,int**,int***,int**,MFErrorHandler);
void EdgeOnManifold(MFKVector,MFKVector,MFAtlas,int,int*,int*,int*,double**,int**,int*,int*,int**,int**,int*,int*,int**,int***,int**,MFErrorHandler);
void TriangleOnManifold(MFKVector,MFKVector,MFKVector,MFAtlas,int,int*,int*,int*,double**,int**,int*,int*,int**,int**,int*,int*,int**,int***,int**,MFErrorHandler);
void ThreeDDual(MFAtlas,int,int*,int*,int*,double**,int**,int*,int*,int**,int**,int*,int*,int**,int***,int**,MFErrorHandler);

MFChart MFAtlasChart(MFAtlas,int,MFErrorHandler);
MFPolytope MFChartPolytope(MFChart,MFErrorHandler);

void MFDrawProjectTo3d(MFChart,MFNVector,float*,float*,float*,MFErrorHandler);
int MFDrawGetData(MFChart,MFNVector,double*,MFErrorHandler);

void MFAtlasToPOV(MFAtlas A,char *name, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAtlasToPOV"};
  int mV,iV,iV0;
  double *v=NULL;
  int *vt=NULL;
  int mE,iE;
  int *ev=NULL;
  int *et=NULL;
  int mF,iF;
  int *nfv=NULL;
  int **fv=NULL;
  int *ft=NULL;
  int i,j,l,n,m;

  iV=0;
  iV0=0;
  mV=0;

  iE=0;
  mE=0;

  iF=0;
  mF=0;

  m=MFDrawGetData(MFAtlasChart(A,0,e),NULL,(double*)NULL,e)+3;

/* Have Polytope and intersection of Polytope and ball. */
/* Also Dual Polytopes (no balls) */
/* Polytope is vertices, edges and faces */
/* Dual Polytope is vertices, edges and faces */
/* Ball is center and sphere */
/* n=1, k=1 */
/* n=2, k=1 */
/* n=2, k=2 */
/* n=3, k=1 */
/* n=3, k=2 */
/* n=3, k=3 started here */

  printf("Add in chart centers\n");fflush(stdout);
  ChartCenters(A,m,&iV,&iV0,&mV,&v,&vt,&iE,&mE,&ev,&et,&iF,&mF,&nfv,&fv,&ft,e);
  printf("  There are now %d pts, %d edges %d faces\n",iV,iE,iF);fflush(stdout);
  printf("Add in Polytopes\n");fflush(stdout);
  ThreeDPolytopes   (A,m,&iV,&iV0,&mV,&v,&vt,&iE,&mE,&ev,&et,&iF,&mF,&nfv,&fv,&ft,e);
  printf("  There are now %d pts, %d edges %d faces\n",iV,iE,iF);fflush(stdout);
/*printf("Add in Dual");fflush(stdout);
  ThreeDDual        (A,m,&iV,&iV0,&mV,&v,&vt,&iE,&mE,&ev,&et,&iF,&mF,&nfv,&fv,&ft,e);
  printf("  There are now %d pts, %d edges %d faces\n",iV,iE,iF);fflush(stdout);*/
  printf("Dump it to %s\n",name);fflush(stdout);
  MFMeshToPOVFile(name,iV,m,v,vt,iE,ev,et,iF,nfv,fv,ft,e);
/*MFMeshToDXFile2(name,iV,m,v,vt,iE,ev,et,iF,nfv,fv,ft,e);*/
  if(v!=NULL)free(v);
  if(nfv!=NULL)free(nfv);
  if(fv!=NULL)
   {
    for(i=0;i<iF;i++)
     {
      if(fv[i]!=NULL)free(fv[i]);
     }
    free(fv);
   }
  return;
 }

/* Color table */

#define PI 3.1415926

double b(double v, double v0, double v1, MFErrorHandler e)
 {
  static char RoutineName[]={"b"};
  int i;
  double result;
  double t;

  while(v<v0)v+=v1-v0;
  while(v>v1)v-=v1-v0;
  i=3*(v-v0)/(v1-v0) >= 0 ? (int)(3*(v-v0)/(v1-v0)+0.5) : (int)(3*(v-v0)/(v1-v0)-0.5);
  t=3*(v-v0)/(v1-v0)-i;

  if(i==0)
    result=1.;
   else if(i==1)
    result=1-t;
   else if(i==2)
    result=0;
   else
    result=0.;

  return result;
 }

double g(double v, double v0, double v1, MFErrorHandler e)
 {
  static char RoutineName[]={"g"};
  int i;
  double result;
  double t;

  while(v<v0)v+=v1-v0;
  while(v>v1)v-=v1-v0;
  i=3*(v-v0)/(v1-v0) >= 0 ? (int)(3*(v-v0)/(v1-v0)+0.5) : (int)(3*(v-v0)/(v1-v0)-0.5);
  t=3*(v-v0)/(v1-v0)-i;

  if(i==0)
    result=t;
   else if(i==1)
    result=1;
   else if(i==2)
    result=1;
   else
    result=1-t;

  return result;
 }

double r(double v, double v0, double v1, MFErrorHandler e)
 {
  static char RoutineName[]={"r"};
  int i;
  double result;
  double t;

  while(v<v0)v+=v1-v0;
  while(v>v1)v-=v1-v0;
  i=3*(v-v0)/(v1-v0) >= 0 ? (int)(3*(v-v0)/(v1-v0)+0.5) : (int)(3*(v-v0)/(v1-v0)-0.5);
  t=3*(v-v0)/(v1-v0)-i;

  if(i==0)
    result=0;
   else if(i==1)
    result=0;
   else if(i==2)
    result=t;
   else
    result=1;

  return result;
 }

void MFMeshToPOVFile(char *name,int nV, int m, double *v, int *vt, int nE, int *ev, int *et, int nF, int *nfv, int **fv, int *ft, MFErrorHandler err)
 {
  static char RoutineName[]={"MFMeshToPOVFile"};
  int i,j,l;
  int ivt,nvt;
  int n;
  int iet,net;
  int e;
  int ift,nft;
  FILE *fid;
  char truename[1024];
  int md;

  strcpy(truename,name);
  strcat(truename,".pov");
  printf("%s, opening file %s for output\n",RoutineName,truename);fflush(stdout);
  fid=fopen(truename,"w");
  if(fid==NULL)
   {
    printf("Open failed\n");fflush(stdout);
    return;
   }

  md=3;
  if(m<3)md=m;

/* Vertices */

  if(nV>0)
   {
    nvt=vt[0];
    for(i=0;i<nV;i++)
     if(vt[i]>nvt)nvt=vt[i];
   }

  for(ivt=0;ivt<nvt+1;ivt++)
   {
    fprintf(fid,"#if(Vertices%d)\n",ivt);
    for(i=0;i<nV;i++)
     {
      if(vt[i]==ivt)
       {
        fprintf(fid,"sphere { <");
        for(j=0;j<md;j++)
         {
          if(j>0)fprintf(fid,",");
          fprintf(fid,"%lf",v[j+m*i]);fflush(fid);
         }
        fprintf(fid,">, .01 texture { Vertex%dTexture }}\n",ivt);
       }
     }
    fprintf(fid,"\n#end\n");
   }
  fflush(fid);

/*    Edges */

  if(nE>0)
   {
    net=et[0];
    for(i=0;i<nE;i++)
     if(et[i]>net)net=et[i];
   }

  for(iet=0;iet<net+1;iet++)
   {
    fprintf(fid,"\n#if(Edges%d)\n",iet);
    for(i=0;i<nE;i++)
     {
      if(et[i]==iet)
       {
        if(fabs(v[0+m*ev[2*i]]-v[0+m*ev[2*i+1]])
          +fabs(v[1+m*ev[2*i]]-v[1+m*ev[2*i+1]])
          +fabs(v[2+m*ev[2*i]]-v[2+m*ev[2*i+1]])>3.e-4)
         {
          fprintf(fid,"cylinder { <");
          e=ev[2*i];
          for(l=0;l<md;l++)
           {
            if(l>0)fprintf(fid,",");
            fprintf(fid,"%lf",v[l+m*e]);fflush(fid);
           }
          fprintf(fid,">,<");
          e=ev[2*i+1];
          for(l=0;l<md;l++)
           {
            if(l>0)fprintf(fid,",");
            fprintf(fid,"%lf",v[l+m*e]);fflush(fid);
           }
          fprintf(fid,">, .002 texture { Edge%dTexture } }\n",iet);
         }
       }
     }
    fprintf(fid,"\n#end\n");
   }

/*    Faces */

  nft=0;
  if(nF>0)
   {
    nft=ft[0];
    for(i=0;i<nF;i++)
     if(ft[i]>nft)nft=ft[i];
   }

  for(ift=0;ift<nft+1;ift++)
   {
    n=0;for(i=0;i<nF;i++)if(ft[i]==ift)n++;
    if(n>0)
     {
      fprintf(fid,"\n#if(Faces%d)\n",ift);
      fprintf(fid,"mesh{\n");
      for(i=0;i<nF;i++)
       {
        if(ft[i]==ift)
         {
          for(j=1;j<nfv[i]-1;j++)
           {
            if((fabs(v[0+m*(fv[i])[0]]-v[0+m*(fv[i])[j]])
               +fabs(v[1+m*(fv[i])[0]]-v[1+m*(fv[i])[j]])
               +fabs(v[2+m*(fv[i])[0]]-v[2+m*(fv[i])[j]])>1.e-4)
             &&
               (fabs(v[0+m*(fv[i])[0]]-v[0+m*(fv[i])[j+1]])
               +fabs(v[1+m*(fv[i])[0]]-v[1+m*(fv[i])[j+1]])
               +fabs(v[2+m*(fv[i])[0]]-v[2+m*(fv[i])[j+1]])>1.e-4)
             &&
               (fabs(v[0+m*(fv[i])[j]]-v[0+m*(fv[i])[j+1]])
               +fabs(v[1+m*(fv[i])[j]]-v[1+m*(fv[i])[j+1]])
               +fabs(v[2+m*(fv[i])[j]]-v[2+m*(fv[i])[j+1]])>1.e-4) )
             {
              if(1||m<4)
               {
                fprintf(fid,"triangle { <");
                e=(fv[i])[0];
                for(l=0;l<md;l++)
                 {
                  if(l>0)fprintf(fid,",");
                  fprintf(fid,"%lf",v[l+m*e]);fflush(fid);
                 }
                fprintf(fid,">,<");
                e=(fv[i])[j];
                for(l=0;l<md;l++)
                 {
                  if(l>0)fprintf(fid,",");
                  fprintf(fid,"%lf",v[l+m*e]);fflush(fid);
                 }
                fprintf(fid,">,<");
                e=(fv[i])[j+1];
                for(l=0;l<md;l++)
                 {
                  if(l>0)fprintf(fid,",");
                  fprintf(fid,"%lf",v[l+m*e]);fflush(fid);
                 }
                fprintf(fid,"> texture { Face%dTexture }}\n",ift);
               }else{
                fprintf(fid,"cv_triangle ( <");
                e=(fv[i])[0];
                for(l=0;l<md;l++)
                 {
                  if(l>0)fprintf(fid,",");
                  fprintf(fid,"%lf",v[l+m*e]);fflush(fid);
                 }
                fprintf(fid,">, rgb <%lf,%lf,%lf>, <",r(v[md+m*e],0.,2*PI,err),g(v[md+m*e],0.,2*PI,err),b(v[md+m*e],0.,2*PI,err));
                e=(fv[i])[j];
                for(l=0;l<md;l++)
                 {
                  if(l>0)fprintf(fid,",");
                  fprintf(fid,"%lf",v[l+m*e]);fflush(fid);
                 }
                fprintf(fid,">, rgb <%lf,%lf,%lf>, <",r(v[md+m*e],0.,2*PI,err),g(v[md+m*e],0.,2*PI,err),b(v[md+m*e],0.,2*PI,err));
                e=(fv[i])[j+1];
                for(l=0;l<md;l++)
                 {
                  if(l>0)fprintf(fid,",");
                  fprintf(fid,"%lf",v[l+m*e]);fflush(fid);
                 }
                fprintf(fid,">, rgb <%lf,%lf,%lf>",r(v[md+m*e],0.,2*PI,err),g(v[md+m*e],0.,2*PI,err),b(v[md+m*e],0.,2*PI,err));
                fprintf(fid,")\n");
               }
             }
           }
         }
       }
      fprintf(fid," }\n");
      fprintf(fid,"\n#end\n");
     }
   }

  fflush(fid);
  fclose(fid);

  return;
 }

void ThreeDPolytopes(MFAtlas A,int m, int *piV,int *piV0,int *pmV,double **pv,int **pvt,int *piE,int *pmE,int **pev,int **pet, int *piF,int *pmF,int **pnfv,int ***pfv, int **pft, MFErrorHandler e)
 {
  static char RoutineName[]={"ThreeDPolytopes"};
  int mV,iV,iV0;
  double *v;
  int *vt;
  int mE,iE;
  int *ev;
  int *et;
  int mF,iF;
  int *nfv;
  int **fv;
  int *ft;
  int i,j,l,n;
  int e0,e1;
  MFNVector u;
  int chart;
  MFEnumPolytope EP;
  int rc;
  float x,y,z;

  iV=*piV;
  iV0=*piV0;
  mV=*pmV;
  v=*pv;
  vt=*pvt;
  iE=*piE;
  mE=*pmE;
  ev=*pev;
  et=*pet;
  iF=*piF;
  mF=*pmF;
  nfv=*pnfv;
  fv=*pfv;
  ft=*pft;

  u=MFCreateNVector(MFAtlasN(A,e),e);
  for(chart=0;chart<MFAtlasNumberOfCharts(A,e);chart++)
   {
    iV0=iV;
    EP=MFEnumeratePolytope(MFChartPolytope(MFAtlasChart(A,chart,e),e),e);

/* Vertices */

    while(iV+MFEnumPolytopeNumberOfCells(EP,0,e)>=mV)
     {
      mV+=100;
      v=(double*)realloc((void*)v,m*mV*sizeof(double));

#ifndef MFNOSAFETYNET
      if(v==NULL)
       {
        sprintf(MFPOVErrorMsg,"Out of memory trying to allocate %d bytes.\n",3*mV*sizeof(double));
        MFSetError(e,12,RoutineName,MFPOVErrorMsg,__LINE__,__FILE__);
        MFErrorHandlerOutOfMemory(e);
        return;
       }
#endif

      vt=(int*)realloc((void*)vt,mV*sizeof(int));

#ifndef MFNOSAFETYNET
      if(vt==NULL)
       {
        sprintf(MFPOVErrorMsg,"Out of memory trying to allocate %d bytes.\n",mV*sizeof(int));
        MFSetError(e,12,RoutineName,MFPOVErrorMsg,__LINE__,__FILE__);
        MFErrorHandlerOutOfMemory(e);
        return;
       }
#endif

     }
    for(i=0;i<MFEnumPolytopeNumberOfCells(EP,0,e);i++)
     {
/*    rc=MFChartEvaluate(MFAtlasChart(A,chart,e),MFEnumPolytopeVertex(EP,i,e),u,e);*/
      MFChartPointInTangentSpace(MFAtlasChart(A,chart,e),MFEnumPolytopeVertex(EP,i,e),u,e);
      if(!strcmp(MFImplicitMFId(MFAtlasMF(A,e),e),"TPBVP"))
       {
        MFDrawProjectTo3d(MFAtlasChart(A,chart,e),u,&x,&y,&z,e);
        v[m*iV  ]=x;
        v[m*iV+1]=y;
        v[m*iV+2]=z;
        if(m>3)MFDrawGetData(MFAtlasChart(A,chart,e),u,v+m*iV+3,e);
       }else{
        int i;

        MFIMFProjectToDraw(MFAtlasMF(A,e),u,v+m*iV,e);
        v[m-3+m*iV]=MFChartHasBoundary(MFAtlasChart(A,chart,e),e);
        v[m-2+m*iV]=MFChartIsSingular(MFAtlasChart(A,chart,e),e);
        v[m-1+m*iV]=MFNVGetIndex(MFAtlasCenterOfChart(A,chart,e),e);
       }

      vt[iV]=0;
      iV++;
     }

/* These will need to be broken into short segments and projected onto the MF */

/* Edges */

    while(iE+MFEnumPolytopeNumberOfCells(EP,1,e)>=mE)
     {
      mE+=100;
      ev=(int*)realloc((void*)ev,2*mE*sizeof(int));

#ifndef MFNOSAFETYNET
      if(ev==NULL)
       {
        sprintf(MFPOVErrorMsg,"Out of memory trying to allocate %d bytes.\n",2*mE*sizeof(int));
        MFSetError(e,12,RoutineName,MFPOVErrorMsg,__LINE__,__FILE__);
        MFErrorHandlerOutOfMemory(e);
        return;
       }
#endif

      et=(int*)realloc((void*)et,mE*sizeof(int));

#ifndef MFNOSAFETYNET
      if(et==NULL)
       {
        sprintf(MFPOVErrorMsg,"Out of memory trying to allocate %d bytes.\n",mE*sizeof(int));
        MFSetError(e,12,RoutineName,MFPOVErrorMsg,__LINE__,__FILE__);
        MFErrorHandlerOutOfMemory(e);
        return;
       }
#endif

     }
    for(i=0;i<MFEnumPolytopeNumberOfCells(EP,1,e);i++)
     {
      ev[2*iE  ]=MFEnumPolytopeCellFace(EP,1,i,0,e)+iV0;
      ev[2*iE+1]=MFEnumPolytopeCellFace(EP,1,i,1,e)+iV0;
      et[iE]=0;
      iE++;
     }

/* These should be broken into small triangles and projected onto MF */

/* Faces */

    while(iF+MFEnumPolytopeNumberOfCells(EP,2,e)>=mF)
     {
      mF+=100;
      nfv=(int*)realloc((void*)nfv,mF*sizeof(int));

#ifndef MFNOSAFETYNET
      if(nfv==NULL)
       {
        sprintf(MFPOVErrorMsg,"Out of memory trying to allocate %d bytes.\n",mF*sizeof(int));
        MFSetError(e,12,RoutineName,MFPOVErrorMsg,__LINE__,__FILE__);
        MFErrorHandlerOutOfMemory(e);
        return;
       }
#endif

      fv=(int**)realloc((void*)fv,mF*sizeof(int*));

#ifndef MFNOSAFETYNET
      if(fv==NULL)
       {
        sprintf(MFPOVErrorMsg,"Out of memory trying to allocate %d bytes.\n",mF*sizeof(int*));
        MFSetError(e,12,RoutineName,MFPOVErrorMsg,__LINE__,__FILE__);
        MFErrorHandlerOutOfMemory(e);
        return;
       }
#endif

      for(i=mF-100;i<mF;i++){nfv[i]=0;fv[i]=NULL;}
      ft=(int*)realloc((void*)ft,mF*sizeof(int));

#ifndef MFNOSAFETYNET
      if(ft==NULL)
       {
        sprintf(MFPOVErrorMsg,"Out of memory trying to allocate %d bytes.\n",mF*sizeof(int));
        MFSetError(e,12,RoutineName,MFPOVErrorMsg,__LINE__,__FILE__);
        MFErrorHandlerOutOfMemory(e);
        return;
       }
#endif

     }

    for(i=0;i<MFEnumPolytopeNumberOfCells(EP,2,e);i++)
     {
      if(MFEnumPolytopeCellIndex(EP,2,i,0,e)>7)ft[iF]=0;
        else ft[iF]=1;
      nfv[iF]=MFEnumPolytopeNumberOfCellFaces(EP,2,i,e);
      fv[iF]=(int*)malloc(nfv[iF]*sizeof(int));

#ifndef MFNOSAFETYNET
      if(fv[iF]==NULL)
       {
        sprintf(MFPOVErrorMsg,"Out of memory trying to allocate %d bytes.\n",nfv[iF]*sizeof(int));
        MFSetError(e,12,RoutineName,MFPOVErrorMsg,__LINE__,__FILE__);
        MFErrorHandlerOutOfMemory(e);
        return;
       }
#endif

      e0=MFEnumPolytopeCellFace(EP,2,i,0,e);

      (fv[iF])[0]=MFEnumPolytopeCellFace(EP,1,e0,0,e)+iV0;
      (fv[iF])[1]=MFEnumPolytopeCellFace(EP,1,e0,1,e)+iV0;
      for(j=1;j<nfv[iF]-1;j++)
       {
        for(l=0;l<nfv[iF];l++)
         {
          e1=MFEnumPolytopeCellFace(EP,2,i,l,e);
          if(e1!=e0&&MFEnumPolytopeCellFace(EP,1,e1,0,e)+iV0==(fv[iF])[j])
           {
            (fv[iF])[j+1]=MFEnumPolytopeCellFace(EP,1,e1,1,e)+iV0;
            l=nfv[iF];
           }else if(e1!=e0&&MFEnumPolytopeCellFace(EP,1,e1,1,e)+iV0==(fv[iF])[j])
           {
            (fv[iF])[j+1]=MFEnumPolytopeCellFace(EP,1,e1,0,e)+iV0;
            l=nfv[iF];
           }
         }
        e0=e1;
       }
      iF++;
     }
    MFFreeEnumPolytope(EP,e);
   }
  MFFreeNVector(u,e);

  *piV=iV;
  *piV0=iV0;
  *pmV=mV;
  *pv=v;
  *pvt=vt;

  *piE=iE;
  *pmE=mE;
  *pev=ev;
  *pet=et;

  *piF=iF;
  *pmF=mF;
  *pnfv=nfv;
  *pfv=fv;
  *pft=ft;

  return;
 }

void ChartCenters(MFAtlas A,int m, int *piV,int *piV0,int *pmV,double **pv,int **pvt,int *piE,int *pmE,int **pev,int **pet, int *piF,int *pmF,int **pnfv,int ***pfv, int **pft, MFErrorHandler e)
 {
  static char RoutineName[]={"ChartCenters"};
  int mV,iV,iV0;
  double *v;
  int *vt;
  int mE,iE;
  int *ev;
  int *et;
  int mF,iF;
  int *nfv;
  int **fv;
  int *ft;
  MFNVector u;
  int chart;
  int j;
  float x,y,z;
  double tmp[1000];

  iV=*piV;
  iV0=*piV0;
  mV=*pmV;
  v=*pv;
  vt=*pvt;
  iE=*piE;
  mE=*pmE;
  ev=*pev;
  et=*pet;
  iF=*piF;
  mF=*pmF;
  nfv=*pnfv;
  fv=*pfv;
  ft=*pft;

  printf("%s\n",RoutineName);fflush(stdout);

  while(iV+MFAtlasNumberOfCharts(A,e)>=mV)
   {
    mV+=100;
    v=(double*)realloc((void*)v,m*mV*sizeof(double));

#ifndef MFNOSAFETYNET
      if(v==NULL)
       {
        sprintf(MFPOVErrorMsg,"Out of memory trying to allocate %d bytes.\n",3*mV*sizeof(double));
        MFSetError(e,12,RoutineName,MFPOVErrorMsg,__LINE__,__FILE__);
        MFErrorHandlerOutOfMemory(e);
        return;
       }
#endif

    vt=(int*)realloc((void*)vt,mV*sizeof(int));

#ifndef MFNOSAFETYNET
      if(vt==NULL)
       {
        sprintf(MFPOVErrorMsg,"Out of memory trying to allocate %d bytes.\n",mV*sizeof(int));
        MFSetError(e,12,RoutineName,MFPOVErrorMsg,__LINE__,__FILE__);
        MFErrorHandlerOutOfMemory(e);
        return;
       }
#endif

   }
  printf("%s, m=%d\n",RoutineName,m);
  for(chart=0;chart<MFAtlasNumberOfCharts(A,e);chart++)
   {
    u=MFChartCenter(MFAtlasChart(A,chart,e),e);

    if(!strcmp(MFImplicitMFId(MFAtlasMF(A,e),e),"TPBVP"))
     {
      MFDrawProjectTo3d(MFAtlasChart(A,chart,e),u,&x,&y,&z,e);
      v[m*iV  ]=x;
      v[m*iV+1]=y;
      v[m*iV+2]=z;
      if(m>3)MFDrawGetData(MFAtlasChart(A,chart,e),u,v+m*iV+3,e);
     }else{
      int i;

      MFIMFProjectToDraw(MFAtlasMF(A,e),u,v+m*iV,e);
      v[m-3+m*iV]=MFChartHasBoundary(MFAtlasChart(A,chart,e),e);
      v[m-2+m*iV]=MFChartIsSingular(MFAtlasChart(A,chart,e),e);
      v[m-1+m*iV]=MFNVGetIndex(MFAtlasCenterOfChart(A,chart,e),e);
     }
    vt[iV]=1;
    iV0=iV;
    iV++;
   }

  *piV=iV;
  *piV0=iV0;
  *pmV=mV;
  *pv=v;
  *pvt=vt;

  *piE=iE;
  *pmE=mE;
  *pev=ev;
  *pet=et;

  *piF=iF;
  *pmF=mF;
  *pnfv=nfv;
  *pfv=fv;
  *pft=ft;

  return;
 }

void EdgeOnManifold(MFKVector s0,MFKVector s1, MFAtlas A,int m,int *piV,int *piV0,int *pmV,double **pv,int **pvt,int *piE,int *pmE,int **pev,int **pet, int *piF,int *pmF,int **pnfv,int ***pfv, int **pft, MFErrorHandler e)
 {
  static char RoutineName[]={"EdgeOnManifold"};
  int mV,iV,iV0;
  double *v;
  int *vt;
  int mE,iE;
  int *ev;
  int *et;
  int mF,iF;
  int *nfv;
  int **fv;
  int *ft;
  MFNVector u;
  int chart;
  int j;
  float x,y,z;

  iV=*piV;
  iV0=*piV0;
  mV=*pmV;
  v=*pv;
  vt=*pvt;
  iE=*piE;
  mE=*pmE;
  ev=*pev;
  et=*pet;
  iF=*piF;
  mF=*pmF;
  nfv=*pnfv;
  fv=*pfv;
  ft=*pft;

  while(iV+MFAtlasNumberOfCharts(A,e)>=mV)
   {
    mV+=100;
    v=(double*)realloc((void*)v,m*mV*sizeof(double));

#ifndef MFNOSAFETYNET
      if(v==NULL)
       {
        sprintf(MFPOVErrorMsg,"Out of memory trying to allocate %d bytes.\n",m*mV*sizeof(double));
        MFSetError(e,12,RoutineName,MFPOVErrorMsg,__LINE__,__FILE__);
        MFErrorHandlerOutOfMemory(e);
        return;
       }
#endif

    vt=(int*)realloc((void*)vt,mV*sizeof(int));

#ifndef MFNOSAFETYNET
      if(vt==NULL)
       {
        sprintf(MFPOVErrorMsg,"Out of memory trying to allocate %d bytes.\n",mV*sizeof(int));
        MFSetError(e,12,RoutineName,MFPOVErrorMsg,__LINE__,__FILE__);
        MFErrorHandlerOutOfMemory(e);
        return;
       }
#endif

   }
  for(chart=0;chart<MFAtlasNumberOfCharts(A,e);chart++)
   {
    u=MFChartCenter(MFAtlasChart(A,chart,e),e);

    MFDrawProjectTo3d(MFAtlasChart(A,chart,e),u,&x,&y,&z,e);
    v[m*iV  ]=x;
    v[m*iV+1]=y;
    v[m*iV+2]=y;
    if(m>3)MFDrawGetData(MFAtlasChart(A,chart,e),u,v+m*iV+3,e);
    vt[iV]=1;
    iV0=iV;
    iV++;
   }

  *piV=iV;
  *piV0=iV0;
  *pmV=mV;
  *pv=v;
  *pvt=vt;

  *piE=iE;
  *pmE=mE;
  *pev=ev;
  *pet=et;

  *piF=iF;
  *pmF=mF;
  *pnfv=nfv;
  *pfv=fv;
  *pft=ft;

  return;
 }

void TriangleOnManifold(MFKVector s0, MFKVector s1, MFKVector s2, MFAtlas A, int m,int *piV,int *piV0,int *pmV,double **pv,int **pvt,int *piE,int *pmE,int **pev,int **pet, int *piF,int *pmF,int **pnfv,int ***pfv, int **pft, MFErrorHandler e)
 {
  static char RoutineName[]={"TriangleOnManifold"};
  int mV,iV,iV0;
  double *v;
  int *vt;
  int mE,iE;
  int *ev;
  int *et;
  int mF,iF;
  int *nfv;
  int **fv;
  int *ft;
  MFNVector u;
  int chart;
  int j;
  float x,y,z;

  iV=*piV;
  iV0=*piV0;
  mV=*pmV;
  v=*pv;
  vt=*pvt;
  iE=*piE;
  mE=*pmE;
  ev=*pev;
  et=*pet;
  iF=*piF;
  mF=*pmF;
  nfv=*pnfv;
  fv=*pfv;
  ft=*pft;

  while(iV+MFAtlasNumberOfCharts(A,e)>=mV)
   {
    mV+=100;
    v=(double*)realloc((void*)v,m*mV*sizeof(double));

#ifndef MFNOSAFETYNET
      if(v==NULL)
       {
        sprintf(MFPOVErrorMsg,"Out of memory trying to allocate %d bytes.\n",m*mV*sizeof(double));
        MFSetError(e,12,RoutineName,MFPOVErrorMsg,__LINE__,__FILE__);
        MFErrorHandlerOutOfMemory(e);
        return;
       }
#endif

    vt=(int*)realloc((void*)vt,mV*sizeof(int));

#ifndef MFNOSAFETYNET
      if(vt==NULL)
       {
        sprintf(MFPOVErrorMsg,"Out of memory trying to allocate %d bytes.\n",mV*sizeof(int));
        MFSetError(e,12,RoutineName,MFPOVErrorMsg,__LINE__,__FILE__);
        MFErrorHandlerOutOfMemory(e);
        return;
       }
#endif

   }
  for(chart=0;chart<MFAtlasNumberOfCharts(A,e);chart++)
   {
    u=MFChartCenter(MFAtlasChart(A,chart,e),e);

    MFDrawProjectTo3d(MFAtlasChart(A,chart,e),u,&x,&y,&z,e);
    v[m*iV  ]=x;
    v[m*iV+1]=y;
    v[m*iV+2]=z;
    if(m>3)MFDrawGetData(MFAtlasChart(A,chart,e),u,v+m*iV+3,e);
    vt[iV]=1;
    iV0=iV;
    iV++;
   }

  *piV=iV;
  *piV0=iV0;
  *pmV=mV;
  *pv=v;
  *pvt=vt;

  *piE=iE;
  *pmE=mE;
  *pev=ev;
  *pet=et;

  *piF=iF;
  *pmF=mF;
  *pnfv=nfv;
  *pfv=fv;
  *pft=ft;

  return;
 }

void ThreeDDual(MFAtlas A,int m, int *piV,int *piV0,int *pmV,double **pv,int **pvt,int *piE,int *pmE,int **pev,int **pet, int *piF,int *pmF,int **pnfv,int ***pfv, int **pft, MFErrorHandler e)
 {
  static char RoutineName[]={"ThreeDDual"};
  int mV,iV,iV0;
  double *v;
  int *vt;
  int mE,iE;
  int *ev;
  int *et;
  int mF,iF;
  int *nfv;
  int **fv;
  int *ft;
  int i,j,l,n;
  int e0,e1;
  MFNVector u;
  int chart;
  int rc;
  MFEnumDualPolytope EP;
  float x,y,z;

  iV=*piV;
  iV0=*piV0;
  mV=*pmV;
  v=*pv;
  vt=*pvt;
  iE=*piE;
  mE=*pmE;
  ev=*pev;
  et=*pet;
  iF=*piF;
  mF=*pmF;
  nfv=*pnfv;
  fv=*pfv;
  ft=*pft;

  EP=MFEnumDualOfAtlas(A,e);

/* Vertices */

  iV0=iV;
  while(iV+MFEnumDualPolytopeNumberOfCells(EP,0,e)>=mV)
   {
    mV+=100;
    v=(double*)realloc((void*)v,m*mV*sizeof(double));

#ifndef MFNOSAFETYNET
      if(v==NULL)
       {
        sprintf(MFPOVErrorMsg,"Out of memory trying to allocate %d bytes.\n",m*mV*sizeof(double));
        MFSetError(e,12,RoutineName,MFPOVErrorMsg,__LINE__,__FILE__);
        MFErrorHandlerOutOfMemory(e);
        return;
       }
#endif

    vt=(int*)realloc((void*)vt,mV*sizeof(int));

#ifndef MFNOSAFETYNET
      if(vt==NULL)
       {
        sprintf(MFPOVErrorMsg,"Out of memory trying to allocate %d bytes.\n",mV*sizeof(int));
        MFSetError(e,12,RoutineName,MFPOVErrorMsg,__LINE__,__FILE__);
        MFErrorHandlerOutOfMemory(e);
        return;
       }
#endif

   }
  for(i=0;i<MFEnumDualPolytopeNumberOfCells(EP,0,e);i++)
   {
    u=MFEnumDualPolytopeVertex(EP,i,e);
    chart=i;
    MFDrawProjectTo3d(MFAtlasChart(A,chart,e),u,&x,&y,&z,e);
    v[m*iV  ]=x;
    v[m*iV+1]=y;
    v[m*iV+2]=z;
    if(m>3)MFDrawGetData(MFAtlasChart(A,chart,e),u,v+m*iV+3,e);
    vt[iV]=2;
    iV++;
   }

/* These will need to be broken into short segments and projected onto the MF */

/* Edges */

  while(iE+MFEnumDualPolytopeNumberOfCells(EP,1,e)>=mE)
   {
    mE+=100;
    ev=(int*)realloc((void*)ev,2*mE*sizeof(int));

#ifndef MFNOSAFETYNET
    if(ev==NULL)
     {
      sprintf(MFPOVErrorMsg,"Out of memory trying to allocate %d bytes.\n",2*mE*sizeof(int));
      MFSetError(e,12,RoutineName,MFPOVErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

    et=(int*)realloc((void*)et,mE*sizeof(int));

#ifndef MFNOSAFETYNET
    if(et==NULL)
     {
      sprintf(MFPOVErrorMsg,"Out of memory trying to allocate %d bytes.\n",mE*sizeof(int));
      MFSetError(e,12,RoutineName,MFPOVErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

   }
  for(i=0;i<MFEnumDualPolytopeNumberOfCells(EP,1,e);i++)
   {
    ev[2*iE  ]=MFEnumDualPolytopeCellFace(EP,1,i,0,e)+iV0;
    ev[2*iE+1]=MFEnumDualPolytopeCellFace(EP,1,i,1,e)+iV0;
    et[iE]=2;
    iE++;
   }

/* These should be broken into small triangles and projected onto MF */

/* Faces */

  while(iF+MFEnumDualPolytopeNumberOfCells(EP,2,e)>=mF)
   {
    mF+=100;
    nfv=(int*)realloc((void*)nfv,mF*sizeof(int));

#ifndef MFNOSAFETYNET
    if(nfv==NULL)
     {
      sprintf(MFPOVErrorMsg,"Out of memory trying to allocate %d bytes.\n",mF*sizeof(int));
      MFSetError(e,12,RoutineName,MFPOVErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

    fv=(int**)realloc((void*)fv,mF*sizeof(int*));

#ifndef MFNOSAFETYNET
    if(fv==NULL)
     {
      sprintf(MFPOVErrorMsg,"Out of memory trying to allocate %d bytes.\n",mF*sizeof(int*));
      MFSetError(e,12,RoutineName,MFPOVErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

    for(i=mF-100;i<mF;i++){nfv[i]=0;fv[i]=NULL;}
    ft=(int*)realloc((void*)ft,mF*sizeof(int));

#ifndef MFNOSAFETYNET
    if(ft==NULL)
     {
      sprintf(MFPOVErrorMsg,"Out of memory trying to allocate %d bytes.\n",mF*sizeof(int));
      MFSetError(e,12,RoutineName,MFPOVErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

   }

  for(i=0;i<MFEnumDualPolytopeNumberOfCells(EP,2,e);i++)
   {
    ft[iF]=2;
    nfv[iF]=MFEnumDualPolytopeNumberOfCellFaces(EP,2,i,e);
    fv[iF]=(int*)malloc(nfv[iF]*sizeof(int));

#ifndef MFNOSAFETYNET
    if(fv[iF]==NULL)
     {
      sprintf(MFPOVErrorMsg,"Out of memory trying to allocate %d bytes.\n",nfv[iF]*sizeof(int));
      MFSetError(e,12,RoutineName,MFPOVErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

    e0=MFEnumDualPolytopeCellFace(EP,2,i,0,e);

    (fv[iF])[0]=MFEnumDualPolytopeCellFace(EP,1,e0,0,e)+iV0;
    (fv[iF])[1]=MFEnumDualPolytopeCellFace(EP,1,e0,1,e)+iV0;
    for(j=1;j<nfv[iF]-1;j++)
     {
      for(l=0;l<nfv[iF];l++)
       {
        e1=MFEnumDualPolytopeCellFace(EP,2,i,l,e);
        if(e1!=e0&&MFEnumDualPolytopeCellFace(EP,1,e1,0,e)+iV0==(fv[iF])[j])
         {
          (fv[iF])[j+1]=MFEnumDualPolytopeCellFace(EP,1,e1,1,e)+iV0;
          l=nfv[iF];
         }else if(e1!=e0&&MFEnumDualPolytopeCellFace(EP,1,e1,1,e)+iV0==(fv[iF])[j])
         {
          (fv[iF])[j+1]=MFEnumDualPolytopeCellFace(EP,1,e1,0,e)+iV0;
          l=nfv[iF];
         }
       }
      e0=e1;
     }
    iF++;
   }
  MFFreeEnumDualPolytope(EP,e);

  *piV=iV;
  *piV0=iV;
  *pmV=mV;
  *pv=v;
  *pvt=vt;

  *piE=iE;
  *pmE=mE;
  *pev=ev;
  *pet=et;

  *piF=iF;
  *pmF=mF;
  *pnfv=nfv;
  *pfv=fv;
  *pft=ft;

  return;
 }

void MFMeshToDXFile2(char *name,int nV, int m, double *v, int *vt, int nE, int *ev, int *et, int nF, int *nfv, int **fv, int *ft, MFErrorHandler e)
 {
  static char RoutineName[]={"MFMeshToDXFile2"};
  int i,j;
  int n;
  FILE *fid;
  int faces;
  int ift,nft;
  int iet,net;
  char truename[1024];
  int md;
  int k=1; /* cause I don't know what it is */
  int nC;
  int iF;
  int member;
  int nnF;

  printf("in %s\n",RoutineName);fflush(stdout);

  strcpy(truename,name);
  strcat(truename,".dx");
  fid=fopen(truename,"w");
  md=m;
  if(md>3)md=3;

  nC=0;
  for(i=0;i<nV;i++)if(vt[i]==1)nC++;

  if(nC>0)
   {
    fprintf(fid,"object \"Centers\" class array type float rank 1 shape %d items %d data follows\n",md,nC);fflush(fid);
    for(i=0;i<nV;i++)
     {
      if(vt[i]==1)
       {
        fprintf(fid,"     ");
        for(j=0;j<md;j++)
         {
          if(j>0)fprintf(fid," ");
          fprintf(fid,"%lf",v[j+m*i]);fflush(fid);
         }
        fprintf(fid,"\n");
       }
     }
   }

  fprintf(fid,"object \"Vertices\" class array type float rank 1 shape %d items %d data follows\n",md,nV);fflush(fid);
  for(i=0;i<nV;i++)
   {
    fprintf(fid,"     ");
    for(j=0;j<md;j++)
     {
      if(j>0)fprintf(fid," ");
      fprintf(fid,"%lf",v[j+m*i]);fflush(fid);
     }
    fprintf(fid,"\n");
   }
  fprintf(fid,"\n");fflush(fid);

  printf("%s, m=%d, md=%d\n",RoutineName, m,md);fflush(stdout);
  if(m>md)
   {
    printf("Dumping Data\n");fflush(stdout);
    fprintf(fid,"object \"Data\" class array type float rank 1 shape %d items %d data follows\n",m,nV);fflush(stdout);
    for(i=0;i<nV;i++)
     {
      fprintf(fid,"     ");
      for(j=0;j<m;j++)
       {
        if(j>0)fprintf(fid," ");
        if(v[j+m*i]==v[j+m*i])fprintf(fid,"%lf",v[j+m*i]);
         else fprintf(fid,"%lf",0.);
        fflush(stdout);
       }
      fprintf(fid,"\n");
     }
    fprintf(fid,"\n");fflush(stdout);
  }

/*if((md==2||md==3)&&k==1)*/
   {
    net=et[0];
    for(i=1;i<nE;i++)
      if(et[i]>net)net++;
    net++;
  
    for(iet=0;iet<net;iet++)
     {
      n=0;for(i=0;i<nE;i++)if(et[i]==iet)n++;
      fprintf(fid,"object \"lines%d\" class array type int rank 1 shape 2 items %d data follows\n",iet,n);fflush(stdout);
      for(i=0;i<nE;i++)
       {
        if(et[i]==iet)fprintf(fid," %d %d\n",ev[2*i],ev[2*i+1]);
       }
      fprintf(fid,"attribute \"ref\" string \"positions\"\n");
      fprintf(fid,"attribute \"element type\" string \"lines\"\n\n");
     }
   }
  
  nft=0;
  if(nF>0)
   {
    nft=ft[0];
    for(i=1;i<nF;i++)
      if(ft[i]>nft)nft++;
    nft++;
   }

  n=-1;
  for(i=0;i<nF;i++)n+=nfv[i]+1;

  nnF=0;
  n=-1;
  for(i=0;i<nF-1;i++)
   {
    n+=nfv[i]+1;
    if(md==2)
     {
      if((v[0+md*n]-v[0+md*(n+nfv[i+1]+1)])*(v[0+md*n]-v[0+md*(n+nfv[i+1]+1)])
        +(v[1+md*n]-v[1+md*(n+nfv[i+1]+1)])*(v[1+md*n]-v[1+md*(n+nfv[i+1]+1)])<.01)nnF++;
     }
    if(md==3)
     {
      if((v[0+md*n]-v[0+md*(n+nfv[i+1]+1)])*(v[0+md*n]-v[0+md*(n+nfv[i+1]+1)])
        +(v[1+md*n]-v[1+md*(n+nfv[i+1]+1)])*(v[1+md*n]-v[1+md*(n+nfv[i+1]+1)])
        +(v[2+md*n]-v[2+md*(n+nfv[i+1]+1)])*(v[2+md*n]-v[2+md*(n+nfv[i+1]+1)])<.01)nnF++;
     }
   }

  if(n<=nV&&nF>0)
   {
    fprintf(fid,"object \"Centers0\" class array type int rank 1 shape 2 items %d data follows\n",nnF);fflush(stdout);
    n=-1;
    for(i=0;i<nF-1;i++)
     {
      n+=nfv[i]+1;
      if(md==2)
       {
        if((v[0+md*n]-v[0+md*(n+nfv[i+1]+1)])*(v[0+md*n]-v[0+md*(n+nfv[i+1]+1)])
          +(v[1+md*n]-v[1+md*(n+nfv[i+1]+1)])*(v[1+md*n]-v[1+md*(n+nfv[i+1]+1)])<.01)fprintf(fid,"%d %d\n",n,n+nfv[i+1]+1);
       }
      if(md==3)
       {
        if((v[0+md*n]-v[0+md*(n+nfv[i+1]+1)])*(v[0+md*n]-v[0+md*(n+nfv[i+1]+1)])
          +(v[1+md*n]-v[1+md*(n+nfv[i+1]+1)])*(v[1+md*n]-v[1+md*(n+nfv[i+1]+1)])
          +(v[2+md*n]-v[2+md*(n+nfv[i+1]+1)])*(v[2+md*n]-v[2+md*(n+nfv[i+1]+1)])<.01)fprintf(fid,"%d %d\n",n,n+nfv[i+1]+1);
       }
     }
    fprintf(fid,"attribute \"ref\" string \"positions\"\n");
    fprintf(fid,"attribute \"element type\" string \"lines\"\n\n");
   }

  if(n<=nV&&nF>0)
   {
    fprintf(fid,"object \"Centers\" class array type float rank 1 shape %d items %d data follows\n",md,nF);fflush(stdout);
    n=-1;
    for(i=0;i<nF;i++)
     {
      n+=nfv[i]+1;
      for(j=0;j<md;j++)
       {
        if(j>0)fprintf(fid," ");
        fprintf(fid,"%lf",v[j+m*n]);fflush(fid);
       }
      fprintf(fid,"\n");
     }
    fprintf(fid,"\n");
   }

  printf("Number of face types is %d\n",nft);
  printf("Number of faces is %d\n",nF);

  for(ift=0;ift<nft;ift++)
   {
    printf("Doing face type %d\n",ift);
    n=0;for(i=0;i<nF;i++)if(ft[i]==ift)n++;

    printf("There are %d faces of this type\n",n);
  
    if(n>0)
     {
      n=0;for(i=0;i<nF;i++)if(ft[i]==ift)n+=nfv[i];
      fprintf(fid,"object \"Edges%d\" class array type int rank 1 shape 1 items %d data follows\n",ift,n);fflush(stdout);
      for(i=0;i<nF;i++)
       {
        if(ft[i]==ift)
         {
          for(j=0;j<nfv[i];j++)
           {
            if((fv[i])[j]<nV&&(fv[i])[j]>-1)fprintf(fid,"%d\n",(fv[i])[j]);
             else fprintf(fid,"0\n");
           }
         }
       }
      fprintf(fid,"attribute \"ref\" string \"positions\"\n\n");
    
      n=0;for(i=0;i<nF;i++)if(ft[i]==ift)n++;
      fprintf(fid,"object \"Loops%d\" class array type int rank 1 shape 1 items %d data follows\n",ift,n);fflush(stdout);
      n=0;
      for(j=0;j<nF;j++)
       {
        if(ft[j]==ift)
         {
          fprintf(fid," %d\n",n);
          n+=nfv[j];
         }
       }
      fprintf(fid,"attribute \"ref\" string \"edges\"\n\n");
    
      n=0;for(i=0;i<nF;i++)if(ft[i]==ift)n++;
      fprintf(fid,"object \"Faces%d\" class array type int rank 1 shape 1 items %d data follows\n",ift,n);fflush(stdout);
      iF=0;
      for(j=0;j<nF;j++)
       if(ft[j]==ift){fprintf(fid," %d\n",iF);iF++;}
      fprintf(fid,"attribute \"ref\" string \"loops\"\n\n");
     }
   }
  faces=1;
  for(ift=0;ift<nft;ift++)
   {
    n=0;for(i=0;i<nF;i++)if(ft[i]==ift)n++;
    if(n>0)
     {
      fprintf(fid,"object \"f%d\" class field\n",ift);
      fprintf(fid,"component \"positions\" value \"Vertices\"\n");
      if(m>md)fprintf(fid,"component \"data\" value \"Data\"\n");
      if(faces)
       {
        fprintf(fid,"component \"edges\" value \"Edges%d\"\n",ift);
        fprintf(fid,"component \"loops\" value \"Loops%d\"\n",ift);
        fprintf(fid,"component \"faces\" value \"Faces%d\"\n",ift);
        fprintf(fid,"\n");
       }
     }
   }

  n=-1;
  for(i=0;i<nF;i++)n+=nfv[i]+1;

  if(n<=nV&&nF>0)
   {
    fprintf(fid,"object \"faceVertices\" class field\n");
    fprintf(fid,"component \"positions\" value \"Vertices\"\n");
    fprintf(fid,"component \"connections\" value \"Centers0\"\n\n");
   }

  if(nC>0)
   {
    fprintf(fid,"object \"centers\" class field\n");
    fprintf(fid,"component \"positions\" value \"Centers\"\n");
    fprintf(fid,"\n");
   }

/*if((md==2||md==3)&&k==1)*/
   {
    for(iet=0;iet<net;iet++)
     {
      n=0;for(i=0;i<nE;i++)if(et[i]==iet)n++;
      if(n>0)
       {
        fprintf(fid,"object \"l%d\" class field\n",iet);
        fprintf(fid,"component \"positions\" value \"Vertices\"\n");
        fprintf(fid,"component \"connections\" value \"lines%d\"\n",iet);
        if(m>md)fprintf(fid,"component \"data\" value \"Data\"\n");
        fprintf(fid,"\n");
       }
     }
   }

  if(nC>0)
   {
    fprintf(fid,"object \"c0\" class field\n");
    fprintf(fid,"component \"positions\" value \"Vertices\"\n");
    fprintf(fid,"component \"connections\" value \"Centers0\"\n");
    if(m>md)fprintf(fid,"component \"data\" value \"Data\"\n");
    fprintf(fid,"\n");
   }

  fprintf(fid,"\n");
  member=0;
  fprintf(fid,"object \"default\" class group\n");

  for(ift=0;ift<nft;ift++)
   {
    n=0;for(i=0;i<nF;i++)if(ft[i]==ift)n++;
    if(n>0){fprintf(fid," member %d \"f%d\"\n",member,ift);member++;}
   }

  n=-1;
  for(i=0;i<nF;i++)n+=nfv[i]+1;

  if(n<=nV&&nF>0)
   {
    fprintf(fid," member %d \"faceVertices\"\n",member);member++;
   }


  for(iet=0;iet<net;iet++)
   {
    n=0;for(i=0;i<nE;i++)if(et[i]==iet)n++;
    if(n>0){fprintf(fid," member %d \"l%d\"\n",member,iet);member++;}
   }

  if(nC>0){fprintf(fid," member %d \"centers\"\n",member);member++;
           fprintf(fid," member %d \"c0\"\n",member);member++;}


  fprintf(fid,"end\n");

  fclose(fid);

  return;
 }

void MFAtlasToDX2(MFAtlas A,char *name, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAtlasToDX2"};
  int mV,iV,iV0;
  double *v=NULL;
  int *vt=NULL;
  int mE,iE;
  int *ev=NULL;
  int *et=NULL;
  int mF,iF;
  int *nfv=NULL;
  int **fv=NULL;
  int *ft=NULL;
  int i,j,l,n,m;

  printf("%s\n",RoutineName);fflush(stdout);

  iV=0;
  iV0=0;
  mV=0;

  iE=0;
  mE=0;

  iF=0;
  mF=0;

  if(!strcmp(MFImplicitMFId(MFAtlasMF(A,e),e),"TPBVP"))
    m=MFDrawGetData(MFAtlasChart(A,0,e),NULL,(double*)NULL,e)+3;
   else
    m=MFIMFProjectToDraw(MFAtlasMF(A,e),NULL,(double*)NULL,e)+3;
  printf("   m=%d\n",m);fflush(stdout);

/* Have Polytope and intersection of Polytope and ball. */
/* Also Dual Polytopes (no balls) */
/* Polytope is vertices, edges and faces */
/* Dual Polytope is vertices, edges and faces */
/* Ball is center and sphere */
/* n=1, k=1 */
/* n=2, k=1 */
/* n=2, k=2 */
/* n=3, k=1 */
/* n=3, k=2 */
/* n=3, k=3 started here */

  printf("Add in chart centers");fflush(stdout);
  ChartCenters(A,m,&iV,&iV0,&mV,&v,&vt,&iE,&mE,&ev,&et,&iF,&mF,&nfv,&fv,&ft,e);
  printf("  There are now %d pts, %d edges %d faces\n",iV,iE,iF);fflush(stdout);
  printf("Add in Polytopes");fflush(stdout);
  ThreeDPolytopes   (A,m,&iV,&iV0,&mV,&v,&vt,&iE,&mE,&ev,&et,&iF,&mF,&nfv,&fv,&ft,e);
  printf("  There are now %d pts, %d edges %d faces\n",iV,iE,iF);fflush(stdout);
/*printf("Add in Dual");fflush(stdout);
  ThreeDDual        (A,m,&iV,&iV0,&mV,&v,&vt,&iE,&mE,&ev,&et,&iF,&mF,&nfv,&fv,&ft,e);*/
  printf("  There are now %d pts, %d edges %d faces\n",iV,iE,iF);fflush(stdout);
  printf("Dump it to %s\n",name);fflush(stdout);
  MFMeshToDXFile2(name,iV,m,v,vt,iE,ev,et,iF,nfv,fv,ft,e);
  if(v!=NULL)free(v);
  if(nfv!=NULL)free(nfv);
  if(fv!=NULL)
   {
    for(i=0;i<iF;i++)
     {
      if(fv[i]!=NULL)free(fv[i]);
     }
    free(fv);
   }
  return;
 }

void MFMeshToDXFile3(char *name,int nV, int m, double *v, int *vt, int nE, int *ev, int *et, int nF, int *nfv, int **fv, int *ft, MFErrorHandler e)
 {
  static char RoutineName[]={"MFMeshToDXFile3"};
  int i,j;
  int n;
  FILE *fid;
  int faces;
  int iet,net;
  int ift,nft;
  char truename[1024];
  int md;

  strcpy(truename,name);
  strcat(truename,".dx");
  fid=fopen(truename,"w");
  md=3;

  fprintf(fid,"object \"Vertices\" class array type float rank 1 shape %d items %d data follows\n",md,nV);fflush(stdout);
  for(i=0;i<nV;i++)
   {
    fprintf(fid,"     ");
    for(j=0;j<md;j++)
     {
      if(j>0)fprintf(fid," ");
      fprintf(fid,"%lf",v[j+m*i]);fflush(stdout);
     }
    fprintf(fid,"\n");
   }
  fprintf(fid,"\n");fflush(stdout);

  printf("m=%d, md=%d\n",m,md);fflush(stdout);
  if(m>md)
   {
    printf("Dumping Data\n");fflush(stdout);
    fprintf(fid,"object \"Data\" class array type float rank 1 shape %d items %d data follows\n",m-md,nV);fflush(stdout);
    for(i=0;i<nV;i++)
     {
      fprintf(fid,"     ");
      for(j=md;j<m;j++)
       {
        if(j>md)fprintf(fid," ");
        fprintf(fid,"%lf",v[j+m*i]);fflush(stdout);
       }
      fprintf(fid,"\n");
     }
    fprintf(fid,"\n");fflush(stdout);
  }
  
  net=et[0];
  for(i=1;i<nE;i++)
    if(et[i]>net)net++;
  net++;

  printf("Number of edge types is %d\n",net);
  printf("Number of edges is %d\n",nE);

  for(iet=0;iet<net;iet++)
   {
    n=0;for(i=0;i<nE;i++)if(et[i]==iet)n++;
    fprintf(fid,"object \"Edges%d\" class array type int rank 1 shape 2 items %d data follows\n",iet,n);fflush(stdout);
    for(i=0;i<nE;i++)
     {
      if(et[i]==iet)fprintf(fid," %d %d\n",ev[2*i],ev[2*i+1]);
     }
    fprintf(fid,"attribute \"ref\" string \"positions\"\n\n");
   }
  
  nft=0;
  if(nF>0)
   {
    nft=ft[0];
    for(i=1;i<nF;i++)
      if(ft[i]>nft)nft++;
    nft++;
   }

  printf("Number of face types is %d\n",nft);
  printf("Number of faces is %d\n",nF);

  for(ift=0;ift<nft;ift++)
   {
    printf("Doing face type %d\n",ift);
    n=0;for(i=0;i<nF;i++)if(ft[i]==ift)n++;

    printf("There are %d faces of this type\n",n);
  
    if(n>0)
     {
      n=0;for(i=0;i<nF;i++)if(ft[i]==ift)n+=nfv[i]-2;
      fprintf(fid,"object \"Faces%d\" class array type int rank 1 shape 3 items %d data follows\n",ift,n);fflush(stdout);
      for(i=0;i<nF;i++)
       {
        if(ft[j]==ift)
         {
          for(j=1;j<nfv[i]-1;j++)
           {
            if((fabs(v[0+m*(fv[i])[0]]-v[0+m*(fv[i])[j]])
               +fabs(v[1+m*(fv[i])[0]]-v[1+m*(fv[i])[j]])
               +fabs(v[2+m*(fv[i])[0]]-v[2+m*(fv[i])[j]])>1.e-4)
             &&
               (fabs(v[0+m*(fv[i])[0]]-v[0+m*(fv[i])[j+1]])
               +fabs(v[1+m*(fv[i])[0]]-v[1+m*(fv[i])[j+1]])
               +fabs(v[2+m*(fv[i])[0]]-v[2+m*(fv[i])[j+1]])>1.e-4)
             &&
               (fabs(v[0+m*(fv[i])[j]]-v[0+m*(fv[i])[j+1]])
               +fabs(v[1+m*(fv[i])[j]]-v[1+m*(fv[i])[j+1]])
               +fabs(v[2+m*(fv[i])[j]]-v[2+m*(fv[i])[j+1]])>1.e-4) )
             {
               fprintf(fid," %d %d %d\n",(fv[i])[0],(fv[i])[j],(fv[i])[j+1]);
             }
           }
         }
       }
      fprintf(fid,"attribute \"ref\" string \"positions\"\n\n");
     }
   }
  
  for(iet=0;iet<net;iet++)
   {
    n=0;for(i=0;i<nE;i++)if(et[i]==iet)n++;
    if(n>0)
     {
      fprintf(fid,"object %d class field\n",iet);
      fprintf(fid,"component \"positions\" value \"Vertices\"\n");
      if(m>3)fprintf(fid,"component \"data\" value \"Data\"\n");
      fprintf(fid,"component \"edges\" value \"Edges%d\"\n",iet);
      fprintf(fid,"end\n");
     }
   }

  for(ift=0;ift<nft;ift++)
   {
    n=0;for(i=0;i<nF;i++)if(ft[i]==ift)n++;
    if(n>0)
     {
      fprintf(fid,"object %d class field\n",ift+net);
      fprintf(fid,"component \"positions\" value \"Vertices\"\n");
      if(m>3)fprintf(fid,"component \"data\" value \"Data\"\n");
      fprintf(fid,"component \"faces\" value \"Faces%d\"\n",ift);
      fprintf(fid,"end\n");
     }
   }

  fprintf(fid,"\n");
  fprintf(fid,"object \"default\" class group\n");
  for(ift=0;ift<nft;ift++)
   {
    n=0;for(i=0;i<nE;i++)if(et[i]==iet)n++;
    if(n>0)fprintf(fid," member %d %d\n",iet,iet);
   }
  for(ift=0;ift<nft;ift++)
   {
    n=0;for(i=0;i<nF;i++)if(ft[i]==ift)n++;
    if(n>0)fprintf(fid," member %d %d\n",ift+net,ift+net);
   }
  fprintf(fid,"end\n");

  fclose(fid);

  return;
 }

void MFDualAtlasToDX2(MFAtlas A,char *name, MFErrorHandler e)
 {
  static char RoutineName[]={"MFDualAtlasToPOV"};
  int mV,iV,iV0;
  double *v=NULL;
  int *vt=NULL;
  int mE,iE;
  int *ev=NULL;
  int *et=NULL;
  int mF,iF;
  int *nfv=NULL;
  int **fv=NULL;
  int *ft=NULL;
  int i,j,l,n,m;

  printf("%s\n",RoutineName);fflush(stdout);

  iV=0;
  iV0=0;
  mV=0;

  iE=0;
  mE=0;

  iF=0;
  mF=0;

  m=MFDrawGetData(MFAtlasChart(A,0,e),NULL,(double*)NULL,e)+3;

/* Have Polytope and intersection of Polytope and ball. */
/* Also Dual Polytopes (no balls) */
/* Polytope is vertices, edges and faces */
/* Dual Polytope is vertices, edges and faces */
/* Ball is center and sphere */
/* n=1, k=1 */
/* n=2, k=1 */
/* n=2, k=2 */
/* n=3, k=1 */
/* n=3, k=2 */
/* n=3, k=3 started here */

  printf("Add in chart centers");fflush(stdout);
  ChartCenters(A,m,&iV,&iV0,&mV,&v,&vt,&iE,&mE,&ev,&et,&iF,&mF,&nfv,&fv,&ft,e);
  printf("  There are now %d pts, %d edges %d faces\n",iV,iE,iF);fflush(stdout);
/*printf("Add in Polytopes");fflush(stdout);
  ThreeDPolytopes   (A,m,&iV,&iV0,&mV,&v,&vt,&iE,&mE,&ev,&et,&iF,&mF,&nfv,&fv,&ft,e);
  printf("  There are now %d pts, %d edges %d faces\n",iV,iE,iF);fflush(stdout);*/
  printf("Add in Dual");fflush(stdout);
  ThreeDDual        (A,m,&iV,&iV0,&mV,&v,&vt,&iE,&mE,&ev,&et,&iF,&mF,&nfv,&fv,&ft,e);
  printf("  There are now %d pts, %d edges %d faces\n",iV,iE,iF);fflush(stdout);
  printf("Dump it to %s\n",name);fflush(stdout);
  MFMeshToDXFile2(name,iV,m,v,vt,iE,ev,et,iF,nfv,fv,ft,e);
  if(v!=NULL)free(v);
  if(nfv!=NULL)free(nfv);
  if(fv!=NULL)
   {
    for(i=0;i<iF;i++)
     {
      if(fv[i]!=NULL)free(fv[i]);
     }
    free(fv);
   }
  return;
 }

void MFMeshToVBMFile(char *name,int k,int nV, int m, double *v, int *vt, int nE, int *ev, int *et, int nF, int *nfv, int **fv, int *ft, MFErrorHandler e)
 {
  static char RoutineName[]={"MFMeshToVBMFile"};
  int i,j;
  int n;
  FILE *fid;
  int faces;
  int ift,nft;
  int iet,net;
  char truename[1024];
  int md;

  strcpy(truename,name);
  strcat(truename,".vbm");
  fid=fopen(truename,"w");

  fprintf(fid,"# VBM_Default_Plot_Coordinates 0 1 2 3\n");fflush(fid);
  for(i=0;i<m;i++)fprintf(fid,"x%d;",i);
  fprintf(fid,"\n");fflush(fid);
  fprintf(fid,"%d\n",m);fflush(fid);
  fprintf(fid,"%d %d\n",nV,1);fflush(fid);
  fprintf(fid,"%d\n",nV);fflush(fid);
  for(i=0;i<nV;i++)
   {
    for(j=0;j<m;j++)
     {
      if(j>0)fprintf(fid," ");
      fprintf(fid,"%lf",v[j+m*i]);fflush(stdout);
     }
    fprintf(fid,"\n");
   }
  fprintf(fid,"\n");fflush(stdout);

  if(k==1)
   {
    fprintf(fid,"##Begin VBM_Connectivity\n");fflush(fid);
    for(i=0;i<nE;i++)fprintf(fid,"line 1 %d 1 %d\n",ev[2*i],ev[2*i+1]);
    fprintf(fid,"##End\n");fflush(fid);
   }

  if(k==2)
   {
    fprintf(fid,"##Begin VBM_Connectivity\n");fflush(fid);
    for(i=0;i<nF;i++)
     {
      for(j=1;j<nfv[i]-1;j++)
       fprintf(fid,"triangle 1 %d 1 %d 1 %d\n",(fv[i])[0],(fv[i])[j],(fv[i])[j+1]);
     }
    fprintf(fid,"##End\n");fflush(fid);
   }
    
  fclose(fid);

  return;
 }

#ifdef __cplusplus
}
#endif
