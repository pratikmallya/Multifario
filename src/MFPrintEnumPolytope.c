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

static char *id="@(#) $Id: MFPrintEnumPolytope.c,v 1.2 2007/02/13 01:22:34 mhender Exp $";

static char MFPrintErrorMsg[256]="";

#include <MFPrint.h>
#include <MFErrorHandler.h>
#include <MFKVector.h>
#include <MFAtlasFriends.h>
#include <MFEnumPolytope.h>
#include <MFEnumDualPolytope.h>

#define DENSE 0
#define VECTOR 1
#define LOCA 2

#ifdef __cplusplus
 extern "C" {
#endif

void MFPrintEnumPolytope(FILE *fout,MFEnumPolytope EP, MFErrorHandler e)
 {
  static char RoutineName[]={"MFPrintEnumPolytope"};
  int k,d;
  int i,n;
  int j,m;
  MFKVector v;

#ifdef MFNOCONFIDENCE
  if(fout==NULL)
   {
    sprintf(MFPrintErrorMsg,"fid (argument 1) is NULL.");
    MFSetError(e,12,RoutineName,MFPrintErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(EP==NULL)
   {
    sprintf(MFPrintErrorMsg,"Enumerated Polytope (argument 2) is NULL.");
    MFSetError(e,12,RoutineName,MFPrintErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  k=MFEnumPolytopeDimension(EP,e);
  for(d=0;d<k+1;d++)
   {
    fprintf(fout," %d cells\n",d);
    n=MFEnumPolytopeNumberOfCells(EP,d,e);
    for(i=0;i<n;i++)
     {
      fprintf(fout,"  %2d  ",i);fflush(fout);
      fprintf(fout,"[");fflush(fout);
      m=MFEnumPolytopeNumberOfCellIndices(EP,d,i,e);
      for(j=0;j<m;j++)
       {
        if(j>0)fprintf(fout," ");
        fprintf(fout,"%2d",MFEnumPolytopeCellIndex(EP,d,i,j,e));fflush(fout);
       }
      fprintf(fout,"]   (");
      m=MFEnumPolytopeNumberOfCellFaces(EP,d,i,e);
      for(j=0;j<m;j++)
       {
        if(j>0)fprintf(fout,",");
        fprintf(fout,"%2d",MFEnumPolytopeCellFace(EP,d,i,j,e));fflush(fout);
       }
      fprintf(fout,") (");fflush(fout);
      m=MFEnumPolytopeNumberOfFaceCells(EP,d,i,e);
      for(j=0;j<m;j++)
       {
        if(j>0)fprintf(fout,",");
        fprintf(fout,"%2d",MFEnumPolytopeFaceCell(EP,d,i,j,e));fflush(fout);
       }
      fprintf(fout,") ");fflush(fout);
      v=MFEnumPolytopeVertex(EP,i,e);
      if(d==0&&v!=NULL)MFPrintKVector(stdout,v,e);
      fprintf(fout,"\n");fflush(fout);
     }
   }

  return;
 }

void MFPrintEnumDualPolytope(FILE *fout,MFEnumDualPolytope EP, MFErrorHandler e)
 {
  static char RoutineName[]={"MFPrintEnumDualPolytope"};
  int k,d;
  int i,n;
  int j,m;
  MFNVector v;

#ifdef MFNOCONFIDENCE
  if(fout==NULL)
   {
    sprintf(MFPrintErrorMsg,"fid (argument 1) is NULL.");
    MFSetError(e,12,RoutineName,MFPrintErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(EP==NULL)
   {
    sprintf(MFPrintErrorMsg,"Enumerated Dual Polytope (argument 2) is NULL.");
    MFSetError(e,12,RoutineName,MFPrintErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  k=MFEnumDualPolytopeDimension(EP,e);
  for(d=0;d<k+1;d++)
   {
    fprintf(fout," %d cells\n",d);
    n=MFEnumDualPolytopeNumberOfCells(EP,d,e);
    for(i=0;i<n;i++)
     {
      fprintf(fout,"  %2d  ",i);fflush(fout);
      fprintf(fout,"[");fflush(fout);
      m=MFEnumDualPolytopeNumberOfCellIndices(EP,d,i,e);
      for(j=0;j<m;j++)
       {
        if(j>0)fprintf(fout," ");
        fprintf(fout,"%2d",MFEnumDualPolytopeCellIndex(EP,d,i,j,e));fflush(fout);
       }
      fprintf(fout,"]   (");
      m=MFEnumDualPolytopeNumberOfCellFaces(EP,d,i,e);
      for(j=0;j<m;j++)
       {
        if(j>0)fprintf(fout,",");
        fprintf(fout,"%2d",MFEnumDualPolytopeCellFace(EP,d,i,j,e));fflush(fout);
       }
      fprintf(fout,") (");fflush(fout);
      m=MFEnumDualPolytopeNumberOfFaceCells(EP,d,i,e);
      for(j=0;j<m;j++)
       {
        if(j>0)fprintf(fout,",");
        fprintf(fout,"%2d",MFEnumDualPolytopeFaceCell(EP,d,i,j,e));fflush(fout);
       }
      fprintf(fout,") ");fflush(fout);
      v=MFEnumDualPolytopeVertex(EP,i,e);
      if(d==0&&v!=NULL)MFPrintNVector(stdout,v,e);
      fprintf(fout,"\n");fflush(fout);
     }
   }

  return;
 }

#ifdef __cplusplus
}
#endif
