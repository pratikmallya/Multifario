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

static char *id="@(#) $Id: MFDrawClippedSphere.c,v 1.3 2007/02/13 01:22:30 mhender Exp $";

static char MFDrawClippedSphereErrorMsg[256]="";

#include <MFAtlas.h>
#include <MFChart.h>
#include <MFNVector.h>
#include <MFKVector.h>
#include <MFEnumPolytope.h>
#include <math.h>
#include <stdlib.h>
#include <sh.h>
#include <MFPrint.h>
#include <MFDrawChartState.h>
#include <MFDrawClippedSphere.h>

#define sqrt3 1.73205080

#ifdef __cplusplus
 extern "C" {
#endif

struct MFDrawClippedSpherePolygonListSt { 
        int nPolygons;
        int maxPolygon;
        double **polygon;
        int *nVertices;
        float **normal;
       };

typedef struct MFDrawClippedSpherePolygonListSt *MFDrawClippedSpherePolygonList;

MFKVector MFTESTPT=NULL;

MFChart MFAtlasChart(MFAtlas,int,MFErrorHandler);
int MFDrawTestInChart3d(MFChart,double,double,double,MFErrorHandler);

MFDrawClippedSpherePolygonList MFDrawCreateClippedSpherePolygonList(MFErrorHandler);
void MFDrawCheckClippedSpherePolygonListSize(MFDrawClippedSpherePolygonList,MFErrorHandler);
void MFDrawFreeClippedSpherePolygonList(MFDrawClippedSpherePolygonList,MFErrorHandler);

int MFDrawAddTriangle(MFDrawClippedSpherePolygonList,double,double,double,double,double,double,double,double,double,MFErrorHandler);
int MFDrawAddPolygonList(MFDrawClippedSpherePolygonList,int,double*,MFErrorHandler);
void MFDrawFreePolygonList(int,MFDrawClippedSpherePolygonList,MFErrorHandler);
void MFDrawCheckTriangle(int,double,double,double,double,MFDrawClippedSpherePolygonList,MFChart,double,double,double,MFErrorHandler);
void MFDrawClipTriangle(double,double,double,double,double,double,double,double,double,double,double,double,double,int*,double*,MFChart,double,MFErrorHandler);
void MFDrawEdgeBisection(double,double,double,double,double,double,double,int,double,double,double,int,MFChart,double*,double*,double*,double,MFErrorHandler);
void MFDrawProjectTo3d(MFChart,MFNVector,float*,float*,float*,MFErrorHandler);

extern MFAtlas MFMFA;

void MFDrawClippedSphere(MFChart chart, double epsilon1, double epsilon2, double epsilon3, MFChartState state, MFErrorHandler e)
 {
  static char RoutineName[]={"MFDrawClippedSphere"};
  double h;
  double xx,yy,zz;
  double xa,ya,za;
  double xb,yb,zb;
  double xc,yc,zc;
  double xma,yma,zma;
  double xmb,ymb,zmb;
  double xmc,ymc,zmc;
  double d;
  int i,j,n;
  float *x,*y,*z;
  MFDrawClippedSpherePolygonList pList;
  MFKVector s;
  MFNVector u;
  double x0,y0,z0,r;

/*   Draws the portion of a sphere which satisies test(x,y,z)>0             */
/*                                                                          */
/* (x0,y0,z0) origin of sphere                                              */
/* r          radius of sphere                                              */
/* epsilon1   size of largest triangle to draw                              */
/* epsilon2   size of smallest triangle to draw                             */
/* epsilon3   size of error in points on boundary                           */
/* test       routine to use to determine if a point is to be drawn or not. */

  x0=0.;
  y0=0.;
  z0=0.;
  r =MFChartRadius(chart,e);
  
  pList=MFDrawCreateClippedSpherePolygonList(e);

  s=MFCreateKVector(3,e);
  u=MFCreateNVector(MFChartN(chart,e),e);

/* Initial Tetrahedron */
 
  MFDrawAddTriangle(pList, x0  ,y0  ,z0+r,x0+r,y0  ,z0  ,x0  ,y0+r,z0  ,e);
  MFDrawAddTriangle(pList, x0  ,y0  ,z0+r,x0  ,y0+r,z0  ,x0-r,y0  ,z0  ,e);
  MFDrawAddTriangle(pList, x0  ,y0  ,z0+r,x0-r,y0  ,z0  ,x0  ,y0-r,z0  ,e);
  MFDrawAddTriangle(pList, x0  ,y0  ,z0+r,x0  ,y0-r,z0  ,x0+r,y0  ,z0  ,e);

  MFDrawAddTriangle(pList, x0  ,y0  ,z0-r,x0  ,y0+r,z0  ,x0+r,y0  ,z0,e);
  MFDrawAddTriangle(pList, x0  ,y0  ,z0-r,x0-r,y0  ,z0  ,x0  ,y0+r,z0,e);
  MFDrawAddTriangle(pList, x0  ,y0  ,z0-r,x0  ,y0-r,z0  ,x0-r,y0  ,z0,e);
  MFDrawAddTriangle(pList, x0  ,y0  ,z0-r,x0+r,y0  ,z0  ,x0  ,y0-r,z0,e);

/* Now subdivide until smallest side is small enough */

  h=r*sqrt(2.);
  while(h>epsilon1)
   {
    n=pList->nPolygons;
    for(i=0;i<n;i++)
     {
      xa=(pList->polygon[i])[0];
      ya=(pList->polygon[i])[1];
      za=(pList->polygon[i])[2];
      xb=(pList->polygon[i])[3];
      yb=(pList->polygon[i])[4];
      zb=(pList->polygon[i])[5];
      xc=(pList->polygon[i])[6];
      yc=(pList->polygon[i])[7];
      zc=(pList->polygon[i])[8];

/* Compute midpoints of the sides, and project out onto the sphere */

      xma=(xa+xb)*.5;
      yma=(ya+yb)*.5;
      zma=(za+zb)*.5;
      d=r/sqrt(pow(xma-x0,2)+pow(yma-y0,2)+pow(zma-z0,2));
      xma=(xma-x0)*d+x0;
      yma=(yma-y0)*d+y0;
      zma=(zma-z0)*d+z0;
      xmb=(xb+xc)*.5;
      ymb=(yb+yc)*.5;
      zmb=(zb+zc)*.5;
      d=r/sqrt(pow(xmb-x0,2)+pow(ymb-y0,2)+pow(zmb-z0,2));
      xmb=(xmb-x0)*d+x0;
      ymb=(ymb-y0)*d+y0;
      zmb=(zmb-z0)*d+z0;
      xmc=(xc+xa)*.5;
      ymc=(yc+ya)*.5;
      zmc=(zc+za)*.5;
      d=r/sqrt(pow(xmc-x0,2)+pow(ymc-y0,2)+pow(zmc-z0,2));
      xmc=(xmc-x0)*d+x0;
      ymc=(ymc-y0)*d+y0;
      zmc=(zmc-z0)*d+z0;
      (pList->polygon[i])[0]=xma;
      (pList->polygon[i])[1]=yma;
      (pList->polygon[i])[2]=zma;
      (pList->polygon[i])[3]=xmb;
      (pList->polygon[i])[4]=ymb;
      (pList->polygon[i])[5]=zmb;
      (pList->polygon[i])[6]=xmc;
      (pList->polygon[i])[7]=ymc;
      (pList->polygon[i])[8]=zmc;
      MFDrawAddTriangle(pList, xa ,ya ,za ,xma,yma,zma,xmc,ymc,zmc,e);
      MFDrawAddTriangle(pList, xb ,yb ,zb ,xmb,ymb,zmb,xma,yma,zma,e);
      MFDrawAddTriangle(pList, xc ,yc ,zc ,xmc,ymc,zmc,xmb,ymb,zmb,e);
     }
    h=h/2;
   }

/* Now have a list of triangles covering the sphere, with no edge larger than epsilon1 */
  
  n=pList->nPolygons;
  for(i=0;i<n;i++)
   {
    MFDrawCheckTriangle(i,x0,y0,z0,r,pList,chart,h,epsilon2,epsilon3,e);
   }
  
/* Compute the normals to the sphere at each point */

  for(i=0;i<pList->nPolygons;i++)
   {
    n=pList->nVertices[i];
    if(n>0)
     {
      pList->normal[i]=(float*)malloc(3*n*sizeof(float));

#ifndef MFNOSAFETYNET
      if(pList->normal[i]==NULL)
       {
        sprintf(MFDrawClippedSphereErrorMsg,"Out of memory trying to allocate %d bytes.\n",3*n*sizeof(float));
        MFSetError(e,12,RoutineName,MFDrawClippedSphereErrorMsg,__LINE__,__FILE__);
        MFErrorHandlerOutOfMemory(e);
        return;
       }
#endif

      for(j=0;j<n;j++)
       {
        (pList->normal[i])[3*j+0]=((pList->polygon[i])[3*j+0]-x0)/r;
        (pList->normal[i])[3*j+1]=((pList->polygon[i])[3*j+1]-y0)/r;
        (pList->normal[i])[3*j+2]=((pList->polygon[i])[3*j+2]-z0)/r;
       }
     }
   }

/* draw it */

  x=NULL;
  y=NULL;
  z=NULL;
  for(i=0;i<pList->nPolygons;i++)
   {
    n=pList->nVertices[i];
    if(n>0)
     {
      x=(float*)realloc((void*)x,n*sizeof(float));

#ifndef MFNOSAFETYNET
      if(x==NULL)
       {
        sprintf(MFDrawClippedSphereErrorMsg,"Out of memory trying to allocate %d bytes.\n",n*sizeof(float));
        MFSetError(e,12,RoutineName,MFDrawClippedSphereErrorMsg,__LINE__,__FILE__);
        MFErrorHandlerOutOfMemory(e);
        return;
       }
#endif

      y=(float*)realloc((void*)y,n*sizeof(float));

#ifndef MFNOSAFETYNET
      if(y==NULL)
       {
        sprintf(MFDrawClippedSphereErrorMsg,"Out of memory trying to allocate %d bytes.\n",n*sizeof(float));
        MFSetError(e,12,RoutineName,MFDrawClippedSphereErrorMsg,__LINE__,__FILE__);
        MFErrorHandlerOutOfMemory(e);
        return;
       }
#endif

      z=(float*)realloc((void*)z,n*sizeof(float));

#ifndef MFNOSAFETYNET
      if(z==NULL)
       {
        sprintf(MFDrawClippedSphereErrorMsg,"Out of memory trying to allocate %d bytes.\n",n*sizeof(float));
        MFSetError(e,12,RoutineName,MFDrawClippedSphereErrorMsg,__LINE__,__FILE__);
        MFErrorHandlerOutOfMemory(e);
        return;
       }
#endif

      for(j=0;j<n;j++)
       {
        MFKVSetC(s,0,(pList->polygon[i])[3*j  ],e);
        MFKVSetC(s,1,(pList->polygon[i])[3*j+1],e);
        MFKVSetC(s,2,(pList->polygon[i])[3*j+2],e);
        MFChartEvaluate(chart,s,u,e);
        MFDrawProjectTo3d(chart,u,x+j,y+j,z+j,e);
       }
      if(state!=NULL)MFChartStateAddPolygonWithNormal(state,n,x,y,z,pList->normal[i],1,e);
       else
        shpg(&n,x,y,z);
     }
   }
  free(x);
  free(y);
  free(z);
  MFFreeKVector(s,e);
  MFFreeNVector(u,e);

/* cleanup */

  MFDrawFreeClippedSpherePolygonList(pList,e);

  return;
 }

int MFDrawAddTriangle(MFDrawClippedSpherePolygonList pList, double x0, double y0,double z0, double x1, double y1,double z1, double x2, double y2,double z2, MFErrorHandler e)
 {
  static char RoutineName[]={"MFDrawAddTriangle"};
  MFDrawCheckClippedSpherePolygonListSize(pList,e);

  (pList->nVertices)[pList->nPolygons]=3;
  (pList->polygon)[pList->nPolygons]=(double*)malloc(3*3*sizeof(double));

#ifndef MFNOSAFETYNET
  if((pList->polygon)[pList->nPolygons]==NULL)
   {
    sprintf(MFDrawClippedSphereErrorMsg,"Out of memory trying to allocate %d bytes.\n",3*3*sizeof(double));
    MFSetError(e,12,RoutineName,MFDrawClippedSphereErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0;
   }
#endif

  ((pList->polygon)[pList->nPolygons])[0+3*0]=x0;
  ((pList->polygon)[pList->nPolygons])[1+3*0]=y0;
  ((pList->polygon)[pList->nPolygons])[2+3*0]=z0;
  ((pList->polygon)[pList->nPolygons])[0+3*1]=x1;
  ((pList->polygon)[pList->nPolygons])[1+3*1]=y1;
  ((pList->polygon)[pList->nPolygons])[2+3*1]=z1;
  ((pList->polygon)[pList->nPolygons])[0+3*2]=x2;
  ((pList->polygon)[pList->nPolygons])[1+3*2]=y2;
  ((pList->polygon)[pList->nPolygons])[2+3*2]=z2;
  (pList->nPolygons)++;
  return((pList->nPolygons)-1);
 }

void MFDrawFreePolygonList(int i, MFDrawClippedSpherePolygonList pList, MFErrorHandler e)
 {
  static char RoutineName[]={"MFDrawFreePolygonList"};
/* Delete polygon i by setting the number of vertices to -1*/

#ifdef MFNOCONFIDENCE
  if(i<0 || i>=pList->nPolygons)
   {
    sprintf(MFDrawClippedSphereErrorMsg,"Attempt to delete a non-existant polygon, (must have 0 < i(=%d) <= nPolygons(=%d))\n",i, pList->nPolygons);
    MFSetError(e,12,RoutineName,MFDrawClippedSphereErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  (pList->nVertices)[i]=-1;
  return;
 }

int MFDrawAddPolygonList(MFDrawClippedSpherePolygonList pList, int n, double* x, MFErrorHandler e)
 {
  static char RoutineName[]={"MFDrawAddPolygonList"};
  int i;
  MFDrawCheckClippedSpherePolygonListSize(pList,e);

  (pList->nVertices)[pList->nPolygons]=n;
  (pList->polygon)[pList->nPolygons]=(double*)malloc(3*n*sizeof(double));

#ifndef MFNOSAFETYNET

  if((pList->polygon)[pList->nPolygons]==NULL)
   {
    sprintf(MFDrawClippedSphereErrorMsg,"Out of memory trying to allocate %d bytes.\n",3*n*sizeof(double));
    MFSetError(e,12,RoutineName,MFDrawClippedSphereErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0;
   }
#endif

  for(i=0;i<3*n;i++)((pList->polygon)[pList->nPolygons])[i]=x[i];
  (pList->nPolygons)++;
  return((pList->nPolygons)-1);
 }

void MFDrawCheckTriangle(int i, double x0, double y0, double z0, double r, MFDrawClippedSpherePolygonList pList, MFChart chart,double size, double epsilon1,double epsilon2 , MFErrorHandler e)
 {
  static char RoutineName[]={"MFDrawCheckTriangle"};

/* Check and subdivide a triangle */

  double xa,ya,za;
  double xb,yb,zb;
  double xc,yc,zc;
  double xma,yma,zma;
  double xmb,ymb,zmb;
  double xmc,ymc,zmc;
  double d;
  int a,b,c,ma,mb,mc;
  int n=0;
  double *x=NULL;

#ifdef MFNOCONFIDENCE
  if((pList->nVertices)[i]!=3)
   {
    sprintf(MFDrawClippedSphereErrorMsg,"Attempt to delete a non-triangular polygon, (nVertices=%d)\n",(pList->nVertices)[i]);
    MFSetError(e,12,RoutineName,MFDrawClippedSphereErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

/* Get vertices, and check for all in or all out */

  xa=((pList->polygon)[i])[0];
  ya=((pList->polygon)[i])[1];
  za=((pList->polygon)[i])[2];
  a=MFDrawTestInChart3d(chart,xa,ya,za,e);

  xb=((pList->polygon)[i])[3];
  yb=((pList->polygon)[i])[4];
  zb=((pList->polygon)[i])[5];
  b=MFDrawTestInChart3d(chart,xb,yb,zb,e);

  xc=((pList->polygon)[i])[6];
  yc=((pList->polygon)[i])[7];
  zc=((pList->polygon)[i])[8];
  c=MFDrawTestInChart3d(chart,xc,yc,zc,e);

  if(a && b && c )return;
  MFDrawFreePolygonList(i,pList,e);
  if(!a && !b && !c )return;

/* Compute midpoints of the sides, and project out onto the sphere */

  xma=(xa+xb)*.5;
  yma=(ya+yb)*.5;
  zma=(za+zb)*.5;
  d=r/sqrt(pow(xma-x0,2)+pow(yma-y0,2)+pow(zma-z0,2));
  xma=(xma-x0)*d+x0;
  yma=(yma-y0)*d+y0;
  zma=(zma-z0)*d+z0;
  ma=MFDrawTestInChart3d(chart,xma,yma,zma,e);

  xmb=(xb+xc)*.5;
  ymb=(yb+yc)*.5;
  zmb=(zb+zc)*.5;
  d=r/sqrt(pow(xmb-x0,2)+pow(ymb-y0,2)+pow(zmb-z0,2));
  xmb=(xmb-x0)*d+x0;
  ymb=(ymb-y0)*d+y0;
  zmb=(zmb-z0)*d+z0;
  mb=MFDrawTestInChart3d(chart,xmb,ymb,zmb,e);

  xmc=(xc+xa)*.5;
  ymc=(yc+ya)*.5;
  zmc=(zc+za)*.5;
  d=r/sqrt(pow(xmc-x0,2)+pow(ymc-y0,2)+pow(zmc-z0,2));
  xmc=(xmc-x0)*d+x0;
  ymc=(ymc-y0)*d+y0;
  zmc=(zmc-z0)*d+z0;
  mc=MFDrawTestInChart3d(chart,xmc,ymc,zmc,e);

/* Add and check the four smaller triangles */

  x=(double*)malloc(30*sizeof(double));

#ifndef MFNOSAFETYNET
  if(x==NULL)
   {
    sprintf(MFDrawClippedSphereErrorMsg,"Out of memory trying to allocate %d bytes.\n",30*sizeof(double));
    MFSetError(e,12,RoutineName,MFDrawClippedSphereErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  if(a || ma || mc )
   {
    if(a && ma && mc )
      MFDrawAddTriangle(pList, xa ,ya ,za ,xma,yma,zma,xmc,ymc,zmc,e);
    else
     {
      if(size>epsilon1)
       {
        i=MFDrawAddTriangle(pList, xa ,ya ,za ,xma,yma,zma,xmc,ymc,zmc,e);
        MFDrawCheckTriangle(i,x0,y0,z0,r,pList,chart,size/2.,epsilon1,epsilon2,e);
       }else{
        MFDrawClipTriangle(x0,y0,z0,r, xa ,ya ,za ,xma,yma,zma,xmc,ymc,zmc,&n,x,chart,epsilon2,e);
        MFDrawAddPolygonList(pList,n,x,e);
       }
     }
   }

  if(b || mb || ma )
   {
    if(b && mb && ma )
      MFDrawAddTriangle(pList, xb ,yb ,zb ,xmb,ymb,zmb,xma,yma,zma,e);
    else
     {
      if(size>epsilon1)
       {
        i=MFDrawAddTriangle(pList, xb ,yb ,zb ,xmb,ymb,zmb,xma,yma,zma,e);
        MFDrawCheckTriangle(i,x0,y0,z0,r,pList,chart,size/2.,epsilon1,epsilon2,e);
       }else{
        MFDrawClipTriangle(x0,y0,z0,r, xb ,yb ,zb ,xmb,ymb,zmb,xma,yma,zma,&n,x,chart,epsilon2,e);
        MFDrawAddPolygonList(pList,n,x,e);
       }
     }
   }

  if(c || mc || mb )
   {
    if(c && mc && mb )
      MFDrawAddTriangle(pList, xc ,yc ,zc ,xmc,ymc,zmc,xmb,ymb,zmb,e);
    else
     {
      if(size>epsilon1)
       {
        i=MFDrawAddTriangle(pList, xc ,yc ,zc ,xmc,ymc,zmc,xmb,ymb,zmb,e);
        MFDrawCheckTriangle(i,x0,y0,z0,r,pList,chart,size/2.,epsilon1,epsilon2,e);
       }else{
        MFDrawClipTriangle(x0,y0,z0,r, xc ,yc ,zc ,xmc,ymc,zmc,xmb,ymb,zmb,&n,x,chart,epsilon2,e);
        MFDrawAddPolygonList(pList,n,x,e);
       }
     }
   }

  if(ma || mb || mc )
   {
    if(ma && mb && mc )
      MFDrawAddTriangle(pList,xma,yma,zma,xmb,ymb,zmb,xmc,ymc,zmc,e);
    else
     {
      if(size>epsilon1)
       {
        i=MFDrawAddTriangle(pList, xma,yma,zma,xmb,ymb,zmb,xmc,ymc,zmc,e);
        MFDrawCheckTriangle(i,x0,y0,z0,r,pList,chart,size/2.,epsilon1,epsilon2,e);
       }else{
        MFDrawClipTriangle(x0,y0,z0,r, xma,yma,zma,xmb,ymb,zmb,xmc,ymc,zmc,&n,x,chart,epsilon2,e);
        MFDrawAddPolygonList(pList,n,x,e);
       }
     }
   }
  free(x);
  return;
 }

void MFDrawClipTriangle(double x0, double y0, double z0, double r, double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc, int *n, double *x2, MFChart chart,double epsilon, MFErrorHandler e)
 {
  static char RoutineName[]={"MFDrawClipTriangle"};

/*    clip a triangle     */

   int i,i1,j,N;
   int in0,in1;
   double *x;
   double xx=0.;
   double yy=0.;
   double zz=0.;

   x=(double*)malloc(9*sizeof(double));

#ifndef MFNOSAFETYNET
   if(x==NULL)
    {
     sprintf(MFDrawClippedSphereErrorMsg,"Out of memory trying to allocate %d bytes.\n",9*sizeof(double));
     MFSetError(e,12,RoutineName,MFDrawClippedSphereErrorMsg,__LINE__,__FILE__);
     MFErrorHandlerOutOfMemory(e);
     return;
    }
#endif

   x[0]=xa;
   x[1]=ya;
   x[2]=za;
   x[3]=xb;
   x[4]=yb;
   x[5]=zb;
   x[6]=xc;
   x[7]=yc;
   x[8]=zc;

   in0=MFDrawTestInChart3d(chart,x[0],x[1],x[2],e);

   j=0;
   for(i=0;i<3;i++)
    {
     i1=(i+1)%3;
     in1=MFDrawTestInChart3d(chart,x[0+3*i1],x[1+3*i1],x[2+3*i1],e);

     if(in0 && in1)
      {
       x2[3*j]  =x[3*i];
       x2[3*j+1]=x[3*i+1];
       x2[3*j+2]=x[3*i+2];
       if(j==0 || fabs(x2[3*j]-x2[3*j-2])+fabs(x2[3*j+1]-x2[3*j-1])+fabs(x2[3*j+2]-x2[3*j])>1.e-5)j++;

       x2[3*j]  =x[3*i1];
       x2[3*j+1]=x[3*i1+1];
       x2[3*j+2]=x[3*i1+2];
       if(j==0 || fabs(x2[3*j]-x2[3*j-2])+fabs(x2[3*j+1]-x2[3*j-1])+fabs(x2[3*j+2]-x2[3*j])>1.e-5)j++;
      }else if(!in0 && in1)
      {
       MFDrawEdgeBisection(x0,y0,z0,r,x[3*i],x[3*i+1],x[3*i+2],in0,x[3*i1],x[3*i1+1],x[3*i1+2],in1,chart,&xx,&yy,&zz,epsilon,e);

       x2[3*j]  =xx;
       x2[3*j+1]=yy;
       x2[3*j+2]=zz;
       if(j==0 || fabs(x2[3*j]-x2[3*j-2])+fabs(x2[3*j+1]-x2[3*j-1])+fabs(x2[3*j+2]-x2[3*j])>1.e-5)j++;

       x2[3*j]  =x[3*i1];
       x2[3*j+1]=x[3*i1+1];
       x2[3*j+2]=x[3*i1+2];
       if(j==0 || fabs(x2[3*j]-x2[3*j-2])+fabs(x2[3*j+1]-x2[3*j-1])+fabs(x2[3*j+2]-x2[3*j])>1.e-5)j++;
      }else if(in0&&!in1)
      {
       MFDrawEdgeBisection(x0,y0,z0,r,x[3*i],x[3*i+1],x[3*i+2],in0,x[3*i1],x[3*i1+1],x[3*i1+2],in1,chart,&xx,&yy,&zz,epsilon,e);
       x2[3*j]  =xx;
       x2[3*j+1]=yy;
       x2[3*j+2]=zz;
       if(j==0 || fabs(x2[3*j]-x2[3*j-2])+fabs(x2[3*j+1]-x2[3*j-1])+fabs(x2[3*j+2]-x2[3*j])>1.e-5)j++;
      }
     in0=in1;
    }

   if(j>0 && fabs(x2[3*j]-x2[3*j-2])+fabs(x2[3*j+1]-x2[3*j-1])+fabs(x2[3*j+2]-x2[3*j])<=1.e-5)j--;
   *n=j;

   free(x);
   return;
 }

void MFDrawEdgeBisection(double xo,double yo,double zo,double r,double x0,double y0,double z0,int l0,double x1,double y1,double z1,int l1, MFChart chart, double *x, double *y,double *z, double epsilon, MFErrorHandler e)
 {
  static char RoutineName[]={"MFDrawEdgeBisection"};
  int l;
  double d;

  epsilon=3*epsilon;
  if(l0)
   {
    *x=x0;
    *y=y0;
    *z=z0;
   }else{
    *x=x1;
    *y=y1;
    *z=z1;
   }

  while(fabs(x1-x0)+fabs(y1-y0)+fabs(z1-z0)>epsilon)
   {
    *x=.5*(x0+x1);
    *y=.5*(y0+y1);
    *z=.5*(z0+z1);
    d=r/sqrt(pow(*x-xo,2)+pow(*y-yo,2)+pow(*z-zo,2));
    *x=(*x-xo)*d+xo;
    *y=(*y-yo)*d+yo;
    *z=(*z-zo)*d+zo;
    l=MFDrawTestInChart3d(chart,*x,*y,*z,e);

    if(l0&&!l || !l0&&l)
     {
      x1=*x;
      y1=*y;
      z1=*z;
      l1=l;
     }else{
      x0=*x;
      y0=*y;
      z0=*z;
      l0=l;
     }
   }

  return;
 }

MFDrawClippedSpherePolygonList MFDrawCreateClippedSpherePolygonList(MFErrorHandler e)
 {
  static char RoutineName[]={"MFDrawCreateClippedSpherePolygonList"};
  int i;
  MFDrawClippedSpherePolygonList pList;

  pList=(struct MFDrawClippedSpherePolygonListSt*)malloc(sizeof(struct MFDrawClippedSpherePolygonListSt));

#ifndef MFNOSAFETYNET
  if(pList==NULL)
   {
    sprintf(MFDrawClippedSphereErrorMsg,"Out of memory trying to allocate %d bytes.\n",sizeof(struct MFDrawClippedSpherePolygonListSt));
    MFSetError(e,12,RoutineName,MFDrawClippedSphereErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  pList->maxPolygon=100;
  pList->polygon=(double**)malloc(pList->maxPolygon*sizeof(double*));

#ifndef MFNOSAFETYNET
  if(pList->polygon==NULL)
   {
    sprintf(MFDrawClippedSphereErrorMsg,"Out of memory trying to allocate %d bytes.\n",pList->maxPolygon*sizeof(double*));
    MFSetError(e,12,RoutineName,MFDrawClippedSphereErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  pList->nVertices=(int*)malloc(pList->maxPolygon*sizeof(int));

#ifndef MFNOSAFETYNET
  if(pList->nVertices==NULL)
   {
    sprintf(MFDrawClippedSphereErrorMsg,"Out of memory trying to allocate %d bytes.\n",pList->maxPolygon*sizeof(int));
    MFSetError(e,12,RoutineName,MFDrawClippedSphereErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  for(i=0;i<pList->maxPolygon;i++)pList->nVertices[i]=-1;
  pList->nPolygons=0;
  pList->normal=(float**)malloc(pList->maxPolygon*sizeof(float*));

#ifndef MFNOSAFETYNET
  if(pList->normal==NULL)
   {
    sprintf(MFDrawClippedSphereErrorMsg,"Out of memory trying to allocate %d bytes.\n",pList->maxPolygon*sizeof(float*));
    MFSetError(e,12,RoutineName,MFDrawClippedSphereErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  return(pList);
 }

void MFDrawCheckClippedSpherePolygonListSize(MFDrawClippedSpherePolygonList pList, MFErrorHandler e)
 {
  static char RoutineName[]={"MFDrawCheckClippedSpherePolygonListSize"};
  int i;

  if(pList->nPolygons>=pList->maxPolygon)
   {
    pList->maxPolygon+=100;
    pList->polygon=(double**)realloc(pList->polygon,3*(pList->maxPolygon)*sizeof(double*));

#ifndef MFNOSAFETYNET
    if(pList->polygon==NULL)
     {
      sprintf(MFDrawClippedSphereErrorMsg,"Out of memory trying to allocate %d bytes.\n",3*pList->maxPolygon*sizeof(double*));
      MFSetError(e,12,RoutineName,MFDrawClippedSphereErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

    pList->normal=(float**)realloc(pList->normal,3*(pList->maxPolygon)*sizeof(float*));

#ifndef MFNOSAFETYNET
    if(pList->normal==NULL)
     {
      sprintf(MFDrawClippedSphereErrorMsg,"Out of memory trying to allocate %d bytes.\n",3*pList->maxPolygon*sizeof(float*));
      MFSetError(e,12,RoutineName,MFDrawClippedSphereErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

    pList->nVertices=(int*)realloc(pList->nVertices,(pList->maxPolygon)*sizeof(int));

#ifndef MFNOSAFETYNET
    if(pList->nVertices==NULL)
     {
      sprintf(MFDrawClippedSphereErrorMsg,"Out of memory trying to allocate %d bytes.\n",pList->maxPolygon*sizeof(int));
      MFSetError(e,12,RoutineName,MFDrawClippedSphereErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif
  
    for(i=(pList->maxPolygon)-100;i<(pList->maxPolygon);i++)(pList->nVertices)[i]=-1;
   }
 }

void MFDrawFreeClippedSpherePolygonList(MFDrawClippedSpherePolygonList pList, MFErrorHandler e)
 {
  static char RoutineName[]={"MFDrawFreeClippedSpherePolygonList"};
  int i;
   
  for(i=0;i<pList->nPolygons;i++)
   {
    free(pList->polygon[i]);
    if(pList->nVertices[i]>0)free(pList->normal[i]);
   }
  free(pList->polygon);
  free(pList->normal);
  free(pList->nVertices);
  free(pList);
 
  return;
 }

int MFDrawTestInChart3d(MFChart chart,double x,double y,double z, MFErrorHandler e)
 {
  static char RoutineName[]={"MFDrawTestInChart3d"};

  if(MFTESTPT==NULL)MFTESTPT=MFCreateKVector(MFChartK(chart,e),e);
  MFKVSetC(MFTESTPT,0,x,e);
  MFKVSetC(MFTESTPT,1,y,e);
  MFKVSetC(MFTESTPT,2,z,e);
  return MFPolytopeInterior(MFChartPolytope(chart,e),MFTESTPT,e);
 }

void MFTestClippedSphere(MFChart chart, double epsilon1, double epsilon2, double epsilon3, MFErrorHandler e)
 {
  static char RoutineName[]={"MFTestClippedSphere"};
  double h;
  double xx,yy,zz;
  double xa,ya,za;
  double xb,yb,zb;
  double xc,yc,zc;
  double xma,yma,zma;
  double xmb,ymb,zmb;
  double xmc,ymc,zmc;
  double d;
  int i,j,n;
  float *x=NULL;
  float *y=NULL;
  float *z=NULL;
  MFDrawClippedSpherePolygonList pList;
  MFKVector s;
  MFNVector u;
  double x0,y0,z0,r;
  int jj;

/*   Draws the portion of a sphere which satisies test(x,y,z)>0             */
/*                                                                          */
/* (x0,y0,z0) origin of sphere                                              */
/* r          radius of sphere                                              */
/* epsilon1   size of largest triangle to draw                              */
/* epsilon2   size of smallest triangle to draw                             */
/* epsilon3   size of error in points on boundary                           */
/* test       routine to use to determine if a point is to be drawn or not. */

  x0=0.;
  y0=0.;
  z0=0.;
  r =MFChartRadius(chart,e);
  
  pList=MFDrawCreateClippedSpherePolygonList(e);

  s=MFCreateKVector(3,e);
  u=MFCreateNVector(MFChartN(chart,e),e);

/* Initial Tetrahedron */
 
  MFDrawAddTriangle(pList, x0  ,y0  ,z0+r,x0+r,y0  ,z0  ,x0  ,y0+r,z0  ,e);
  MFDrawAddTriangle(pList, x0  ,y0  ,z0+r,x0  ,y0+r,z0  ,x0-r,y0  ,z0  ,e);
  MFDrawAddTriangle(pList, x0  ,y0  ,z0+r,x0-r,y0  ,z0  ,x0  ,y0-r,z0  ,e);
  MFDrawAddTriangle(pList, x0  ,y0  ,z0+r,x0  ,y0-r,z0  ,x0+r,y0  ,z0  ,e);

  MFDrawAddTriangle(pList, x0  ,y0  ,z0-r,x0  ,y0+r,z0  ,x0+r,y0  ,z0,e);
  MFDrawAddTriangle(pList, x0  ,y0  ,z0-r,x0-r,y0  ,z0  ,x0  ,y0+r,z0,e);
  MFDrawAddTriangle(pList, x0  ,y0  ,z0-r,x0  ,y0-r,z0  ,x0-r,y0  ,z0,e);
  MFDrawAddTriangle(pList, x0  ,y0  ,z0-r,x0+r,y0  ,z0  ,x0  ,y0-r,z0,e);

/* Now subdivide until smallest side is small enough */

  h=r*sqrt(2.);
  while(h>epsilon1)
   {
    n=pList->nPolygons;
    for(i=0;i<n;i++)
     {
      xa=(pList->polygon[i])[0];
      ya=(pList->polygon[i])[1];
      za=(pList->polygon[i])[2];
      xb=(pList->polygon[i])[3];
      yb=(pList->polygon[i])[4];
      zb=(pList->polygon[i])[5];
      xc=(pList->polygon[i])[6];
      yc=(pList->polygon[i])[7];
      zc=(pList->polygon[i])[8];

/* Compute midpoints of the sides, and project out onto the sphere */

      xma=(xa+xb)*.5;
      yma=(ya+yb)*.5;
      zma=(za+zb)*.5;
      d=r/sqrt(pow(xma-x0,2)+pow(yma-y0,2)+pow(zma-z0,2));
      xma=(xma-x0)*d+x0;
      yma=(yma-y0)*d+y0;
      zma=(zma-z0)*d+z0;
      xmb=(xb+xc)*.5;
      ymb=(yb+yc)*.5;
      zmb=(zb+zc)*.5;
      d=r/sqrt(pow(xmb-x0,2)+pow(ymb-y0,2)+pow(zmb-z0,2));
      xmb=(xmb-x0)*d+x0;
      ymb=(ymb-y0)*d+y0;
      zmb=(zmb-z0)*d+z0;
      xmc=(xc+xa)*.5;
      ymc=(yc+ya)*.5;
      zmc=(zc+za)*.5;
      d=r/sqrt(pow(xmc-x0,2)+pow(ymc-y0,2)+pow(zmc-z0,2));
      xmc=(xmc-x0)*d+x0;
      ymc=(ymc-y0)*d+y0;
      zmc=(zmc-z0)*d+z0;
      (pList->polygon[i])[0]=xma;
      (pList->polygon[i])[1]=yma;
      (pList->polygon[i])[2]=zma;
      (pList->polygon[i])[3]=xmb;
      (pList->polygon[i])[4]=ymb;
      (pList->polygon[i])[5]=zmb;
      (pList->polygon[i])[6]=xmc;
      (pList->polygon[i])[7]=ymc;
      (pList->polygon[i])[8]=zmc;
      MFDrawAddTriangle(pList, xa ,ya ,za ,xma,yma,zma,xmc,ymc,zmc,e);
      MFDrawAddTriangle(pList, xb ,yb ,zb ,xmb,ymb,zmb,xma,yma,zma,e);
      MFDrawAddTriangle(pList, xc ,yc ,zc ,xmc,ymc,zmc,xmb,ymb,zmb,e);
     }
    h=h/2;
   }

/* Now have a list of triangles covering the sphere, with no edge larger than epsilon1 */
  
  n=pList->nPolygons;
  for(i=0;i<n;i++)
    MFDrawCheckTriangle(i,x0,y0,z0,r,pList,chart,h,epsilon2,epsilon3,e);

/* now have a list of triangles which are on the boundary of the chart, and inside the polytope */
/* Are these points inside any other chart? */

  for(i=0;i<pList->nPolygons;i++)
   {
    n=pList->nVertices[i];
    if(n>0)
     {
      for(j=0;j<n;j++)
       {
        MFKVSetC(s,0,(pList->polygon[i])[3*j  ],e);
        MFKVSetC(s,1,(pList->polygon[i])[3*j+1],e);
        MFKVSetC(s,2,(pList->polygon[i])[3*j+2],e);
        MFChartEvaluate(chart,s,u,e);
        for(jj=0;jj<MFAtlasNumberOfCharts(MFMFA,e);jj++)
         {
          if(MFAtlasChart(MFMFA,j,e)!=chart)
           {
            MFChartProjectIntoTangentSpace(MFAtlasChart(MFMFA,j,e),u,s,e);
            if(MFChartInterior(MFAtlasChart(MFMFA,j,e),s,e))
              printf("ERROR: here's a point that is supposedly on the boundary, but which is in chart %d\n",j);fflush(stdout);
           }
         }
       }
     }
   }
  free(x);
  free(y);
  free(z);
  MFFreeKVector(s,e);
  MFFreeNVector(u,e);

/* cleanup */

  MFDrawFreeClippedSpherePolygonList(pList,e);

  return;
 }

#ifdef __cplusplus
}
#endif
