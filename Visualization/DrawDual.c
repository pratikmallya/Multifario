/* 
    @(#)DrawDual.c	1.3
    02/04/19 16:20:17
   
    PROGRAM NAME:  multifario

    (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
    CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.

    Please refer to the LICENSE file in the top directory

*/
/*      author: Mike Henderson mhender@watson.ibm.com */

#include <MFAtlas.h>
#include <MFImplicitMF.h>
#include <MFNRegion.h>
#include <MFNVector.h>
#include <MFDraw.h>
#include <MFEnumPolytope.h>
#include <MFEnumDualPolytope.h>
#include <MFDX.h>
#include <MFErrorHandler.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include <sh.h>

#ifdef __cplusplus
 extern "C" {
#endif

int main(int argc, char *argv[])
 {
  float xmin,xmax,dist,alpha,beta;
  int gray;
  MFAtlas S;
  MFEnumDualPolytope P;
  MFErrorHandler e;
  FILE *fid;
  char name[1024]="";

  if(argc<2){printf("Usage %s filename\n\n",argv[0]);fflush(stdout);return 8;}

  strcpy(name,argv[1]);
  strcat(name,".atlas");

  fid=fopen(name,"r");
  if(fid==(FILE*)NULL)
   {
    printf("Error, could not open file %s, %s\n",name,strerror(errno));
    fflush(stdout);
    return 12;
   }

  e=MFCreateErrorHandler();

  printf("Reading Atlas %s\n",argv[1]);fflush(stdout);
  S=MFReadAtlas(fid,e);
  fclose(fid);
  printf("Done reading Atlas, %d charts\n",MFAtlasNumberOfCharts(S,e));fflush(stdout);

  shSetOutputResolution(2048,2048);
  shSetOutputFormat("tiff");
  strcpy(name,argv[1]);
  strcat(name,"Dual");
  shSetOutputFilename(name);
  MFDrawInitializeFromFile(argv[1],e);

  printf("EnumerateDual of Atlas\n");fflush(stdout);
  P=MFEnumDualOfAtlas(S,e);
  printf("Draw Enumerated Dual of Atlas\n");fflush(stdout);
  MFDrawEnumDualPolytope(P,e);MFDrawDisplay(e);

  MFFreeEnumDualPolytope(P,e);
  MFFreeAtlas(S,e);
  MFDrawClose(e);

  MFFreeErrorHandler(e);

  return(0);
 }

#ifdef __cplusplus
}
#endif
