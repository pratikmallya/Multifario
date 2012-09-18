/* 
    @(#)DualToDX.c	1.4
    02/04/19 16:20:28
   
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
#include <MFPrint.h>
#include <MFDraw.h>
#include <MFEnumPolytope.h>
#include <MFEnumDualPolytope.h>
#include <MFDX.h>
#include <MFErrorHandler.h>
#include <math.h>
#include <errno.h>
#include <string.h>

#ifdef __cplusplus
 extern "C" {
#endif

#define N 200

void MFDualAtlasToDX2(MFAtlas,char*,MFErrorHandler);

int main(int argc, char *argv[])
 {
  int gray;
  MFAtlas S;
  MFEnumDualPolytope P;
  FILE *fid;
  char name[1024]="";
  MFErrorHandler e;

  if(argc<2){printf("Usage %s filename\n",argv[0]);fflush(stdout);return 8;}

  printf("%s %s\n",argv[0],argv[1]);

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
  printf("Done reading Atlas\n");fflush(stdout);

  strcpy(name,argv[1]);
  strcat(name,"Dual");

  MFDualAtlasToDX2(S,name,e);

/*printf("EnumerateDual of Atlas\n");fflush(stdout);
  P=MFEnumDualOfAtlas(S,e);
  printf("Dump Enumerated Dual of Atlas as DX file\n");fflush(stdout);
  MFDualPolytopeToDXFile(name,P,e);
  MFFreeEnumDualPolytope(P,e);*/

  MFFreeAtlas(S,e);

  MFFreeErrorHandler(e);

  return(0);
 }

#ifdef __cplusplus
}
#endif
