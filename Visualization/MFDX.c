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

static char *id="@(#) $Id: MFDX.c,v 1.3 2007/02/13 01:22:30 mhender Exp $";

static char MFDXErrorMsg[256]="";

#include <stdlib.h>
#include <MFAtlas.h>
#include <MFChart.h>
#include <MFEnumDualPolytope.h>
#include <MFEnumPolytope.h>
#include <MFPrint.h>
#include <MFErrorHandler.h>

#ifdef __cplusplus
 extern "C" {
#endif

void MFMeshToDXFile(char*,int,int,double*,int,int*,int**,MFErrorHandler);

MFChart MFAtlasChart(MFAtlas,int,MFErrorHandler);
MFPolytope MFChartPolytope(MFChart,MFErrorHandler);

void MFDualPolytopeToDXFile(char *name,MFEnumDualPolytope P, MFErrorHandler e)
 {
  static char RoutineName[]={"MFDualPolytopeToDXFile"};
  int i,nV;
  int j,nT,nE,k;
  MFNVector u;
  int iV[24]={0};
  int r,f,nF;
  int n,m;
  FILE *fid;
  int edges,faces,vols;

  int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s\n",RoutineName);fflush(stdout);}
#endif

  fid=fopen(name,"w");

  nV=MFEnumDualPolytopeNumberOfVertices(P,e);
  u=MFEnumDualPolytopeVertex(P,0,e);
  m=MFNV_NC(u,e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf(" This Polytope has %d vertices of dimension %d\n",nV,m);fflush(stdout);}
#endif

  fprintf(fid,"object \"Vertices\" class array type float rank 1 shape %d items %d data follows\n",m,nV);fflush(stdout);
  for(i=0;i<nV;i++)
   {
    u=MFEnumDualPolytopeVertex(P,i,e);
    fprintf(fid,"     ");
    for(j=0;j<m;j++)
     {
      if(j>0)fprintf(fid," ");
      fprintf(fid,"%lf",MFNV_C(u,j,e));fflush(stdout);
     }
    fprintf(fid,"\n");
   }
  fprintf(fid,"\n");fflush(stdout);

  m=MFEnumDualPolytopeDimension(P,e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf(" The Polytope itself is of dimension %d\n",m);fflush(stdout);}
#endif

  if(m>0)nE=MFEnumDualPolytopeNumberOfCells(P,1,e);
   else nE=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf(" There are %d 1-cells\n",nE);fflush(stdout);}
#endif

  n=0;
  for(i=0;i<nE;i++)
   {
    nF=MFEnumDualPolytopeNumberOfFaceCells(P,1,i,e);
    if(nF==0)n++;
   }

#ifdef MFALLOWVERBOSE
  if(verbose){printf("   %d of them have no adjacent 2-cells\n",n);fflush(stdout);}
#endif

  if(n>0)
   {
    edges=1;
    fprintf(fid,"object \"Edges\" class array type int rank 1 shape 2 items %d data follows\n",n);fflush(stdout);
    for(i=0;i<nE;i++)
     {
      nF=MFEnumDualPolytopeNumberOfFaceCells(P,1,i,e);
      if(nF==0)
       {
        nF=MFEnumDualPolytopeNumberOfCellFaces(P,1,i,e);
        for(j=0;j<nF;j++)
         {
          if(j>0)fprintf(fid," ");
          r=MFEnumDualPolytopeCellFace(P,1,i,j,e);
          fprintf(fid,"%d",r);
         }
        fprintf(fid,"\n");fflush(stdout);
       }
     }
    fprintf(fid,"attribute \"element type\" string \"lines\"\n");
    fprintf(fid,"attribute \"ref\" string \"positions\"\n\n");
   }else{
    edges=0;
   }



  if(m>1)nT=MFEnumDualPolytopeNumberOfCells(P,2,e);
   else nT=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf(" There are %d 2-cells\n",nT);fflush(stdout);}
#endif

  n=0;
  for(i=0;i<nT;i++)
   {
    nF=MFEnumDualPolytopeNumberOfCellFaces(P,2,i,e);
    nV=0;
    for(j=0;j<nF;j++)
     {
      r=MFEnumDualPolytopeCellFace(P,2,i,j,e);
      iV[nV]=MFEnumDualPolytopeCellFace(P,1,r,0,e);
      nV++;
      iV[nV]=MFEnumDualPolytopeCellFace(P,1,r,1,e);
      nV++;
     }
    for(j=0;j+1<nV;j++)
     {
      for(k=j+1;k<nV;k++)
       {
        if(iV[j]>iV[k])
         {
          r=iV[j];
          iV[j]=iV[k];
          iV[k]=r;
         }
       }
     }
    k=0;
    for(j=0;j<nV;j++)
     {
      if(j==0||iV[j]!=iV[j-1])k++;
     }
    nF=MFEnumDualPolytopeNumberOfFaceCells(P,2,i,e);
    if(nF==0&&k==3)n++;
    if(0){fprintf(stderr,"Face %d has %d face cells, %d vertices raw and %d vertices collapsed\n",i,nF,nV,k);fflush(stderr);}
   }

#ifdef MFALLOWVERBOSE
  if(verbose){printf("   %d of them have no adjacent 3-cells\n",n);fflush(stdout);}
#endif

  if(n>0)
   {
    faces=1;
    fprintf(fid,"object \"Faces\" class array type int rank 1 shape 3 items %d data follows\n",n);fflush(stdout);
    for(i=0;i<nT;i++)
     {
      nF=MFEnumDualPolytopeNumberOfFaceCells(P,2,i,e);
      if(nF==0)
       {
        nF=MFEnumDualPolytopeNumberOfCellFaces(P,2,i,e);
        nV=0;
        for(j=0;j<nF;j++)
         {
          r=MFEnumDualPolytopeCellFace(P,2,i,j,e);
          iV[nV]=MFEnumDualPolytopeCellFace(P,1,r,0,e);
          nV++;
          iV[nV]=MFEnumDualPolytopeCellFace(P,1,r,1,e);
          nV++;
         }
        for(j=0;j+1<nV;j++)
         {
          for(k=j+1;k<nV;k++)
           {
            if(iV[j]>iV[k])
             {
              r=iV[j];
              iV[j]=iV[k];
              iV[k]=r;
             }
           }
         }
        k=0;
        for(j=0;j<nV;j++)
         {
          if(j==0||iV[j]!=iV[j-1])k++;
         }
        if(k==3)
         {
          for(j=0;j<nV;j++)
           {
            if(j==0||iV[j]!=iV[j-1])
             {
              if(j>0)fprintf(fid," ");
              fprintf(fid,"%d",iV[j]);
             }
           }
          fprintf(fid,"\n");fflush(stdout);
         }
       }
     }
    fprintf(fid,"attribute \"element type\" string \"triangles\"\n");
    fprintf(fid,"attribute \"ref\" string \"positions\"\n\n");
   }else{
    faces=0;
   }

  if(m>2)nT=MFEnumDualPolytopeNumberOfCells(P,3,e);
   else nT=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf(" There are %d 3-cells\n",nT);fflush(stdout);}
#endif

  n=0;
  for(i=0;i<nT;i++)
   {
    nF=MFEnumDualPolytopeNumberOfCellFaces(P,3,i,e);
    nV=0;
    for(j=0;j<nF;j++)
     {
      f=MFEnumDualPolytopeCellFace(P,3,i,j,e);
      nE=MFEnumDualPolytopeNumberOfCellFaces(P,2,f,e);
      for(k=0;k<nE;k++)
       {
        r=MFEnumDualPolytopeCellFace(P,2,f,k,e);
        iV[nV]=MFEnumDualPolytopeCellFace(P,1,r,0,e);
        nV++;
        iV[nV]=MFEnumDualPolytopeCellFace(P,1,r,1,e);
        nV++;
       }
     }
/* Sort */
    for(j=0;j+1<nV;j++)
     {
      for(k=j+1;k<nV;k++)
       {
        if(iV[j]>iV[k])
         {
          r=iV[j];
          iV[j]=iV[k];
          iV[k]=r;
         }
       }
     }
    k=0;
    for(j=0;j<nV;j++)
     {
      if(j==0||iV[j]!=iV[j-1])k++;
     }
    if(k==4)n++;
   }

  if(n>0)
   {
    vols=1;
    fprintf(fid,"object \"Volumes\" class array type int rank 1 shape 4 items %d data follows\n",n);fflush(stdout);
    for(i=0;i<nT;i++)
     {
      nF=MFEnumDualPolytopeNumberOfCellFaces(P,3,i,e);
      nV=0;
      for(j=0;j<nF;j++)
       {
        f=MFEnumDualPolytopeCellFace(P,3,i,j,e);
        nE=MFEnumDualPolytopeNumberOfCellFaces(P,2,f,e);
        for(k=0;k<nE;k++)
         {
          r=MFEnumDualPolytopeCellFace(P,2,f,k,e);
          iV[nV]=MFEnumDualPolytopeCellFace(P,1,r,0,e);
          nV++;
          iV[nV]=MFEnumDualPolytopeCellFace(P,1,r,1,e);
          nV++;
         }
       }
/* Sort */
      for(j=0;j+1<nV;j++)
       {
        for(k=j+1;k<nV;k++)
         {
          if(iV[j]>iV[k])
           {
            r=iV[j];
            iV[j]=iV[k];
            iV[k]=r;
           }
         }
       }
      k=0;
      for(j=0;j<nV;j++)
       {
        if(j==0||iV[j]!=iV[j-1])k++;
       }
      if(k==4)
       {
        for(j=0;j<nV;j++)
         {
          if(j==0||iV[j]!=iV[j-1])
           {
            if(j>0)fprintf(fid," ");
            fprintf(fid,"%d",iV[j]);
           }
         }
        fprintf(fid,"\n");fflush(stdout);
       }
     }
    fprintf(fid,"attribute \"element type\" string \"tetrahedra\"\n");
    fprintf(fid,"attribute \"ref\" string \"positions\"\n\n");
   }else{
    vols=0;
   }

  if(edges)
   {
    fprintf(fid,"object \"Lines\" class field\n");
    fprintf(fid,"component \"positions\" value \"Vertices\"\n");
    fprintf(fid,"component \"connections\" value \"Edges\"\n");
    fprintf(fid,"\n");
   }

  if(faces)
   {
    fprintf(fid,"object \"Triangles\" class field\n");
    fprintf(fid,"component \"positions\" value \"Vertices\"\n");
    fprintf(fid,"component \"connections\" value \"Faces\"\n");
    fprintf(fid,"\n");
   }
  
  if(vols)
   {
    fprintf(fid,"object \"Tetrahedra\" class field\n");
    fprintf(fid,"component \"positions\" value \"Vertices\"\n");
    fprintf(fid,"component \"connections\" value \"Volumes\"\n");
    fprintf(fid,"\n");
   }

  fprintf(fid,"object \"All\" class group\n");
  i=0;
  if(edges){fprintf(fid,"member \"lines\" value \"Lines\"\n",i);i++;}
  if(faces){fprintf(fid,"member \"triangles\" value \"Triangles\"\n",i);i++;};
  if(vols){fprintf(fid,"member  \"tetrahedra\" value \"Tetrahedra\"\n",i);i++;};
  fprintf(fid,"\nDefault \"All\"\n");
  fprintf(fid,"\nend\n");

  fclose(fid);

  return;
 }

void MFAtlasToDX(MFAtlas A,char *name, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAtlasToDX"};
  int chart;
  int mV,iV,iV0;
  double *v=NULL;
  int mF,iF;
  int *nfv=NULL;
  int **fv=NULL;
  int e0,e1;
  int i,j,l,n;
  MFNVector u;
  MFEnumPolytope EP;
  int rc;

  nfv=NULL;
  fv=NULL;
  v=NULL;

  printf("%s\n",RoutineName);fflush(stdout);

  iV=0;
  mV=0;

  iF=0;
  mF=0;

  printf("Loop over all charts\n");fflush(stdout);
  u=MFCreateNVector(MFAtlasN(A,e),e);
  for(chart=0;chart<MFAtlasNumberOfCharts(A,e);chart++)
   {
    iV0=iV;
/*  printf("\n  Chart %d -- iV0=%d\n\n",chart,iV0);fflush(stdout);*/
    EP=MFEnumeratePolytope(MFChartPolytope(MFAtlasChart(A,chart,e),e),e);
/* Vertices */

    if(iV+MFEnumPolytopeNumberOfCells(EP,0,e)>=mV)
     {
      mV+=100;
      v=(double*)realloc((void*)v,3*mV*sizeof(double));

#ifndef MFNOSAFETYNET
      if(v==NULL)
       {
        sprintf(MFDXErrorMsg,"Out of memory trying to allocate %d bytes.\n",3*mV*sizeof(float));
        MFSetError(e,12,RoutineName,MFDXErrorMsg,__LINE__,__FILE__);
        MFErrorHandlerOutOfMemory(e);
        return;
       }
#endif

     }
    for(i=0;i<MFEnumPolytopeNumberOfCells(EP,0,e);i++)
     {
      rc=MFChartEvaluate(MFAtlasChart(A,chart,e),MFEnumPolytopeVertex(EP,i,e),u,e);
      v[3*iV  ]=MFNV_C(u,0,e);
      v[3*iV+1]=MFNV_C(u,1,e);
      v[3*iV+2]=MFNV_C(u,2,e);
      iV++;
     }

    if(iF+MFEnumPolytopeNumberOfCells(EP,2,e)>=mF)
     {
      mF+=100;
      nfv=(int*)realloc((void*)nfv,mF*sizeof(int));

#ifndef MFNOSAFETYNET
      if(nfv==NULL)
       {
        sprintf(MFDXErrorMsg,"Out of memory trying to allocate %d bytes.\n",mF*sizeof(float));
        MFSetError(e,12,RoutineName,MFDXErrorMsg,__LINE__,__FILE__);
        MFErrorHandlerOutOfMemory(e);
        return;
       }
#endif

      fv=(int**)realloc((void*)fv,mF*sizeof(int*));

#ifndef MFNOSAFETYNET
      if(fv==NULL)
       {
        sprintf(MFDXErrorMsg,"Out of memory trying to allocate %d bytes.\n",mF*sizeof(float));
        MFSetError(e,12,RoutineName,MFDXErrorMsg,__LINE__,__FILE__);
        MFErrorHandlerOutOfMemory(e);
        return;
       }
#endif

      for(i=mF-100;i<mF;i++){nfv[i]=0;fv[i]=NULL;}
     }
    for(i=0;i<MFEnumPolytopeNumberOfCells(EP,2,e);i++)
     {
/*    printf("     Face %d, %d edges\n",i,MFEnumPolytopeNumberOfCellFaces(EP,2,i));fflush(stdout);*/
      if(MFEnumPolytopeCellIndex(EP,2,i,0,e)>7)
       {
        nfv[iF]=MFEnumPolytopeNumberOfCellFaces(EP,2,i,e);
        fv[iF]=(int*)malloc(nfv[iF]*sizeof(int));

#ifndef MFNOSAFETYNET
        if(fv[iF]==NULL)
         {
          sprintf(MFDXErrorMsg,"Out of memory trying to allocate %d bytes.\n",nfv[iF]*sizeof(float));
          MFSetError(e,12,RoutineName,MFDXErrorMsg,__LINE__,__FILE__);
          MFErrorHandlerOutOfMemory(e);
          return;
         }
#endif

        e0=MFEnumPolytopeCellFace(EP,2,i,0,e);

/*      printf("     [");
        for(j=0;j<nfv[iF];j++)
         {
          if(j>0)printf(" ");
          printf("%d {%d,%d}",MFEnumPolytopeCellFace(EP,2,i,j,e),MFEnumPolytopeCellFace(EP,1,MFEnumPolytopeCellFace(EP,2,i,j,e),0,e),MFEnumPolytopeCellFace(EP,1,MFEnumPolytopeCellFace(EP,2,i,j,e),1,e)+iV0);
         }
        printf("]\n");fflush(stdout);*/

        (fv[iF])[0]=MFEnumPolytopeCellFace(EP,1,e0,0,e)+iV0;
        (fv[iF])[1]=MFEnumPolytopeCellFace(EP,1,e0,1,e)+iV0;
/*      printf("        e0=%d [%d,%d]\n",e0,MFEnumPolytopeCellFace(EP,1,e0,0,e)+iV0,MFEnumPolytopeCellFace(EP,1,e0,1,e)+iV0);fflush(stdout);*/
        for(j=1;j<nfv[iF]-1;j++)
         {
          for(l=0;l<nfv[iF];l++)
           {
            e1=MFEnumPolytopeCellFace(EP,2,i,l,e);
            if(e1!=e0&&MFEnumPolytopeCellFace(EP,1,e1,0,e)+iV0==(fv[iF])[j])
             {
              (fv[iF])[j+1]=MFEnumPolytopeCellFace(EP,1,e1,1,e)+iV0;
              l=nfv[iF];
/*          printf("        next edge e1=%d [%d,%d]\n",e1,MFEnumPolytopeCellFace(EP,1,e1,0,e),MFEnumPolytopeCellFace(EP,1,e1,1,e));fflush(stdout);*/
             }else if(e1!=e0&&MFEnumPolytopeCellFace(EP,1,e1,1,e)+iV0==(fv[iF])[j])
             {
              (fv[iF])[j+1]=MFEnumPolytopeCellFace(EP,1,e1,0,e)+iV0;
              l=nfv[iF];
/*            printf("        next edge e1=%d [%d,%d]\n",e1,MFEnumPolytopeCellFace(EP,1,e1,1,e),MFEnumPolytopeCellFace(EP,1,e1,0,e));fflush(stdout);*/
             }
           }
          e0=e1;
         }
/*      printf("Done this face  [");fflush(stdout);
        for(j=0;j<nfv[iF];j++)
         {
          if(j>0)printf(" ");
          printf("%d",(fv[iF])[j]);
         }
        printf("]\n");fflush(stdout);*/
        iF++;
       }
     }
    MFFreeEnumPolytope(EP,e);
   }

  MFMeshToDXFile(name,iV,3,v,iF,nfv,fv,e);
  MFFreeNVector(u,e);
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

void MFMeshToDXFile(char *name,int nV, int m, double *v, int nF, int *nfv, int **fv, MFErrorHandler e)
 {
  static char RoutineName[]={"MFMeshToDXFile"};
  int i,j;
  int n;
  FILE *fid;
  int faces;

  fid=fopen(name,"w");

  fprintf(fid,"object \"Vertices\" class array type float rank 1 shape %d items %d data follows\n",m,nV);fflush(stdout);
  for(i=0;i<nV;i++)
   {
    fprintf(fid,"     ");
    for(j=0;j<m;j++)
     {
      if(j>0)fprintf(fid," ");
      fprintf(fid,"%lf",v[j+m*i]);fflush(stdout);
     }
    fprintf(fid,"\n");
   }
  fprintf(fid,"\n");fflush(stdout);

  n=0;
  for(i=0;i<nF;i++)n+=nfv[i];

  faces=0;
  if(n>0)
   {
    faces=1;
    fprintf(fid,"object \"Edges\" class array type int rank 1 shape 1 items %d data follows\n",n);fflush(stdout);
    for(i=0;i<nF;i++)
      for(j=0;j<nfv[i];j++)fprintf(fid,"%d\n",(fv[i])[j]);
    fprintf(fid,"attribute \"element type\" string \"edges\"\n");
    fprintf(fid,"attribute \"ref\" string \"positions\"\n\n");
  
    fprintf(fid,"object \"Loops\" class array type int rank 1 shape 1 items %d data follows\n",nF);fflush(stdout);
    n=0;
    for(i=0;i<nF;i++)
     {
      fprintf(fid," %d\n",n);
      n+=nfv[i];
     }
    fprintf(fid,"attribute \"element type\" string \"loops\"\n");
    fprintf(fid,"attribute \"ref\" string \"edges\"\n\n");
  
    fprintf(fid,"object \"Faces\" class array type int rank 1 shape 1 items %d data follows\n",nF);fflush(stdout);
    for(i=0;i<nF;i++)
      fprintf(fid," %d\n",i);
    fprintf(fid,"attribute \"element type\" string \"faces\"\n");
    fprintf(fid,"attribute \"ref\" string \"loops\"\n\n");
   }


  fprintf(fid,"object \"default\" class field\n");
  i=0;
  fprintf(fid,"component \"positions\" value \"Vertices\"\n",i);i++;
  if(faces)
   {
    fprintf(fid,"component \"edges\" value \"Edges\"\n",i);i++;
    fprintf(fid,"component \"loops\" value \"Loops\"\n",i);i++;
    fprintf(fid,"component \"faces\" value \"Faces\"\n",i);i++;
   }
  fprintf(fid,"\nend\n");

  fclose(fid);

  return;
 }

#ifdef __cplusplus
}
#endif
