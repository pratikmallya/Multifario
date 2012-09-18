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

static char *id="@(#) $Id: MFPrint.c,v 1.3 2011/07/21 17:42:46 mhender Exp $";

static char MFPrintErrorMsg[256]="";

#include <MFPrint.h>
#include <MFErrorHandler.h>
#include <MFKVector.h>
#include <MFAtlasFriends.h>

#define DENSE 0
#define VECTOR 1
#define LOCA 2

#ifdef __cplusplus
 extern "C" {
#endif

void MFPrintAtlas(FILE *fid,MFAtlas A,MFErrorHandler e)
 {
  static char RoutineName[]={"MFPrintAtlas"};
  int i,j;

#ifdef MFNOCONFIDENCE
  if(fid==NULL)
   {
    sprintf(MFPrintErrorMsg,"fid (argument 1) is NULL.");
    MFSetError(e,12,RoutineName,MFPrintErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(A==NULL)
   {
    sprintf(MFPrintErrorMsg,"Atlas (argument 2) is NULL.");
    MFSetError(e,12,RoutineName,MFPrintErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  fprintf(fid,"     Atlas: \n\n");
  fprintf(fid,"       Base space dimension: %d \n",MFAtlasK(A,e));
  fprintf(fid,"       Target space dimension: %d \n",MFAtlasN(A,e));

  fprintf(fid,"\n       %d Half Spaces \n",MFAtlasNumberOfHalfSpaces(A,e));
  for(i=0;i<MFAtlasNumberOfHalfSpaces(A,e);i++)
   {
    fprintf(fid,"        %d (%d), [%d,%d]",i,MFAtlasHalfSpaceIndex(A,i,e),MFAtlasHalfSpaceLeftChart(A,i,e),MFAtlasHalfSpaceRightChart(A,i,e));
    fprintf(fid,"  x.");
    MFPrintKVector(fid,MFAtlasHalfSpaceNormal(A,i,e),e);
    fprintf(fid,"-%lf<0\n",MFAtlasHalfSpaceOrigin(A,i,e));
   }

  fprintf(fid,"\n       %d Charts \n",MFAtlasNumberOfCharts(A,e));
  for(i=0;i<MFAtlasNumberOfCharts(A,e);i++)
   {
    fprintf(fid,"         Chart %d\n",i);
    fprintf(fid,"           Center: ");
    MFPrintNVector(fid,MFChartCenter(MFAtlasChart(A,i,e),e),e);
    fprintf(fid,"           Radius %lf\n",MFAtlasChartRadius(A,i,e),e);
    fprintf(fid,"           Basis for tangent space:\n");
    MFPrintNKMatrix(fid,MFChartTangentSpace(MFAtlasChart(A,i,e),e),e);
    fprintf(fid,"           Polytope:\n");
    MFPrintPolytope(fid,MFChartPolytope(MFAtlasChart(A,i,e),e),e);
   }

  return;
 }

void MFPrintKVector(FILE *fid, MFKVector s, MFErrorHandler e)
 {
  static char RoutineName[]={"MFPrintKVector"};
  int i,n;

#ifdef MFNOCONFIDENCE
  if(fid==NULL)
   {
    sprintf(MFPrintErrorMsg,"fid (argument 1) is NULL.");
    MFSetError(e,12,RoutineName,MFPrintErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(s==NULL)
   {
    sprintf(MFPrintErrorMsg,"Vector (argument 2) is NULL.");
    MFSetError(e,12,RoutineName,MFPrintErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  n=MFKV_NC(s,e);
  fprintf(fid,"(");
  for(i=0;i<n;i++)
   {
    if(i>0)fprintf(fid,",");
    fprintf(fid,"%lf",MFKV_C(s,i,e));
   }
  fprintf(fid,")");

  return;
 }

void MFPrintChart(FILE *fid,MFChart chart, MFErrorHandler e)
 {
  static char RoutineName[]={"MFPrintChart"};
  MFNVector phi;
  int i;

#ifdef MFNOCONFIDENCE
  if(fid==NULL)
   {
    sprintf(MFPrintErrorMsg,"fid (argument 1) is NULL.");
    MFSetError(e,12,RoutineName,MFPrintErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(chart==NULL)
   {
    sprintf(MFPrintErrorMsg,"Chart (argument 2) is NULL.");
    MFSetError(e,12,RoutineName,MFPrintErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  fprintf(fid,"Chart center ");
  MFPrintNVector(fid,MFChartCenter(chart,e),e);
  fprintf(fid,"\n");
  fprintf(fid,"      R %lf\n",MFChartRadius(chart,e));
  fprintf(fid,"      Basis for TS ");
  for(i=0;i<MFChartK(chart,e);i++)
   {
    if(i>0)fprintf(fid,"                   ");
    phi=MFMColumn(MFChartTangentSpace(chart,e),i,e);
    MFPrintNVector(fid,phi,e);fprintf(fid,"\n");
    MFFreeNVector(phi,e);
   }
  fprintf(fid,"      Polytope\n");MFPrintPolytope(fid,MFChartPolytope(chart,e),e);
  fflush(fid);
  return;
 }

void MFPrintPolytope(FILE *file,MFPolytope P, MFErrorHandler e)
 {
  static char RoutineName[]={"MFPrintPolytope"};
  int i,j;
  MFKVector v;

#ifdef MFNOCONFIDENCE
  if(file==NULL)
   {
    sprintf(MFPrintErrorMsg,"fid (argument 1) is NULL.");
    MFSetError(e,12,RoutineName,MFPrintErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(P==NULL)
   {
    sprintf(MFPrintErrorMsg,"Polytope (argument 2) is NULL.");
    MFSetError(e,12,RoutineName,MFPrintErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  if(P==NULL)return;
  v=MFCreateKVector(MFPolytopeDimension(P,e),e);

  fprintf(file," Polytope 0x%8.8x\n",P);
  fprintf(file,"               dimension %d\n",MFPolytopeDimension(P,e));

  fprintf(file,"               %d vertices\n",MFPolytopeNumberOfVertices(P,e));
  for(i=0;i<MFPolytopeNumberOfVertices(P,e);i++)
   {
    fprintf(file,"                   %2d (",i);fflush(file);
    MFPolytopeVertex(P,i,v,e);
    MFPrintKVector(file,v,e);
    fprintf(file,") indices %d (",MFPolytopeNumberOfVertexIndices(P,i,e));fflush(file);
    for(j=0;j<MFPolytopeNumberOfVertexIndices(P,i,e);j++)
     {
      if(j>0)fprintf(file,",");
      fprintf(file,"%4d",MFPolytopeVertexIndex(P,i,j,e));fflush(file);
     }
      
    fprintf(file,")");
    fprintf(file," r=%lf",MFPolytopeRadiusOfVertex(P,i,e));
    fprintf(file,"\n");
    fflush(file);
   }

  fprintf(file,"               %d faces\n",MFPolytopeNumberOfFaces(P,e));
  fflush(file);
  for(i=0;i<MFPolytopeNumberOfFaces(P,e);i++)
   {
    fprintf(file,"                   %2d, index %4d contains %d vertices ",i,MFPolytopeFaceIndex(P,i,e),MFPolytopeNumberOfVerticesOnFace(P,i,e));
    if(MFPolytopeFaceNormal(P,i,e)!=NULL)
     {
/*    fprintf(file,"                          ");*/
      MFPrintKVector(file,MFPolytopeFaceNormal(P,i,e),e);
      fprintf(file,".s<%lf    0x%8.8x\n",MFPolytopeFaceOrigin(P,i,e),MFPolytopeFaceNormal(P,i,e));
      fflush(file);
     }
   }

  MFFreeKVector(v,e);
  return;
 }

void MFPrintPolytopeTerse(FILE *file,MFPolytope P, MFErrorHandler e)
 {
  static char RoutineName[]={"MFPrintPolytopeTerse"};
  int i,j;

#ifdef MFNOCONFIDENCE
  if(file==NULL)
   {
    sprintf(MFPrintErrorMsg,"fid (argument 1) is NULL.");
    MFSetError(e,12,RoutineName,MFPrintErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(P==NULL)
   {
    sprintf(MFPrintErrorMsg,"Polytope (argument 2) is NULL.");
    MFSetError(e,12,RoutineName,MFPrintErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  fprintf(file,"    Dimension %d\n",MFPolytopeDimension(P,e));

  fprintf(file,"    Vertices\n",MFPolytopeNumberOfVertices(P,e));
  for(i=0;i<MFPolytopeNumberOfVertices(P,e);i++)
   {
    fprintf(file,"         %3d [",i);fflush(file);
    for(j=0;j<MFPolytopeNumberOfVertexIndices(P,i,e);j++)
     {
      if(j>0)fprintf(file,",");
      fprintf(file,"%3d",MFPolytopeVertexIndex(P,i,j,e));fflush(file);
     }
      
    fprintf(file,"]");
    fprintf(file,"\n");
    fflush(file);
   }

  fprintf(file,"    Faces\n",MFPolytopeNumberOfFaces(P,e));
  fflush(file);
  for(i=0;i<MFPolytopeNumberOfFaces(P,e);i++)
   {
    fprintf(file,"     %3d %3d vertices\n",MFPolytopeFaceIndex(P,i,e),MFPolytopeNumberOfVerticesOnFace(P,i,e));
   }

  fflush(file);
  return;
 }

void MFPrintAtlasFaceList(FILE *fid,MFAtlas A, MFErrorHandler e)
 {
  static char RoutineName[]={"MFPrintAtlasFaceList"};
  int i,j;

#ifdef MFNOCONFIDENCE
  if(fid==NULL)
   {
    sprintf(MFPrintErrorMsg,"fid (argument 1) is NULL.");
    MFSetError(e,12,RoutineName,MFPrintErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(A==NULL)
   {
    sprintf(MFPrintErrorMsg,"Atlas (argument 2) is NULL.");
    MFSetError(e,12,RoutineName,MFPrintErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  fprintf(fid,"\n     %3d Half Spaces \n",MFAtlasNumberOfHalfSpaces(A,e));
  for(i=0;i<MFAtlasNumberOfHalfSpaces(A,e);i++)
   {
    if(MFAtlasIsHalfSpaceHyperCube(A,i,e))
     {
      fprintf(fid,"     %3d, Hypercube face\n",i);
     }else{
      fprintf(fid,"     %3d, [%3d,%3d]\n",i,MFAtlasHalfSpaceLeftChart(A,i,e),MFAtlasHalfSpaceRightChart(A,i,e));
     }
   }

  return;
 }

#ifdef __cplusplus
}
#endif
