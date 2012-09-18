/* 
 *  %W%
 *  %D% %T%
 * 
 *  author: Mike Henderson mhender@watson.ibm.com
 */

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

void MFAtlasToDX2(MFAtlas,char*,MFErrorHandler);

int main(int argc, char *argv[])
 {
  int gray;
  MFAtlas S;
  MFEnumPolytope P;
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

  printf("Dump Atlas as DX file\n");fflush(stdout);
  strcpy(name,argv[1]);
/*strcat(name,".dx");*/
/*MFAtlasToDX(S,name,e);*/
  MFAtlasToDX2(S,name,e);

  MFFreeAtlas(S,e);
  MFFreeErrorHandler(e);

  return(0);
 }

#ifdef __cplusplus
}
#endif
