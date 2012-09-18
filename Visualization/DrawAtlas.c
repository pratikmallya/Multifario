/* 
    @(#)DrawAtlas.c	1.7
    03/02/26 16:38:10

    PROGRAM NAME:  multifario

    (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
    CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.

    Please refer to the LICENSE file in the top directory
   
*/
/*      author: Mike Henderson mhender@watson.ibm.com */

#include <multifarioConfig.h>
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
#include <sh.h>
#ifdef HAVE_GETOPT_H
#include <getopt.h>
#endif

#ifdef __cplusplus
 extern "C" {
#endif

int main(int argc, char *argv[])
 {
  float xmin,xmax,dist,alpha,beta;
  int gray;
  MFAtlas S;
  FILE *fid;
  char name[1024]="";
  int npixels;
  char format[1024]="";
  int c;
  MFErrorHandler e;

  if(argc<2){printf("Usage %s filename\n",argv[0]);fflush(stdout);return 8;}

/*printf("%s %s\n",argv[0],argv[1]);*/

  npixels=1024;
  format[0]=0x0;
  while((c=getopt(argc,argv,"n:f:"))!=EOF)
   {
    switch(c)
     {
      case 'n':
       sscanf(optarg,"%d",&npixels);
       break;
      case 'f':
       sscanf(optarg,"%s",format);
       break;
     }
   }


  printf("Options: npixels %d\n",npixels);
  printf("         atlas  %s\n",argv[optind]);fflush(stdout);
  printf("         format  %s\n",format);fflush(stdout);

  if(argc<2)
   {
    fprintf(stderr,"%s, no atlas file specified\n",argv[0]);
    return 0;
   }
  strcpy(name,argv[optind]);
  strcat(name,".atlas");
  printf("argc=%d, argv[0]=%s, argv[1]=%s\n",argc,argv[0],argv[optind]);fflush(stdout);
  fid=fopen(name,"r");
  if(fid==(FILE*)NULL)
   {
    printf("Error, could not open file %s, %s\n",name,strerror(errno));
    fflush(stdout);
    return 12;
   }

  e=MFCreateErrorHandler();

  printf("Reading Atlas %s\n",argv[optind]);fflush(stdout);
  S=MFReadAtlas(fid,e);
  fclose(fid);
  printf("Done reading Atlas, %d charts\n",MFAtlasNumberOfCharts(S,e));fflush(stdout);

  shSetOutputResolution(npixels,npixels);
  if(format[0]!=0x0)shSetOutputFormat(format);

  shSetOutputFilename(argv[optind]);
  MFDrawInitializeFromFile(argv[optind],e);
  MFDrawAtlasOnce(S,e);
  MFDrawDisplay(e);

  MFFreeAtlas(S,e);
  MFDrawClose(e);
  MFFreeErrorHandler(e);

  return(0);
 }

#ifdef __cplusplus
}
#endif
