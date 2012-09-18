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

#define N 200
#ifdef SHSOFT
extern int shIMax;
extern int shJMax;
#else
int shIMax;
int shJMax;
#endif

#ifdef __cplusplus
 extern "C" {
#endif

int main(int argc, char *argv[])
 {
  int gray;
  MFAtlas S;
  MFEnumPolytope P;
  MFErrorHandler e;
  FILE *fid,*fout;
  char name[1024]="";

  if(argc<2)
   {
    printf("Usage:\n     %s filename\n Opens file \"filename.atlas\"s\n",argv[0]);
    fflush(stdout);
    return 12;
   }

  strcpy(name,argv[1]);
  strcat(name,".atlas");

  fid=fopen(name,"r");
  if(fid==(FILE*)NULL)
   {
    printf("Error, could not open file %s, %s\n",name,strerror(errno));
    fflush(stdout);
    return 12;
   }

  strcpy(name,argv[1]);
  strcat(name,".centers");
  fout=fopen(name,"w");
  if(fout==(FILE*)NULL)
   {
    printf("Error, could not open file %s for output, %s\n",name,strerror(errno));
    fflush(stdout);
    return 12;
   }

  printf("Reading Atlas %s\n",argv[1]);fflush(stdout);

  e=MFCreateErrorHandler();

  S=MFReadAtlas(fid,e);
  fclose(fid);
  printf("Done reading Atlas\n");fflush(stdout);

  printf("Dump Atlas as centers file\n");fflush(stdout);
  
  MFWriteChartCenters(fout,S,e);
  fclose(fout);

  MFFreeAtlas(S,e);

  MFFreeErrorHandler(e);

  return(0);
 }

#ifdef __cplusplus
}
#endif
