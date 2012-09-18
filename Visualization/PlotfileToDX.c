#include <stdio.h>
#include <MFAtlas.h>
#include <MFDraw.h>
#include <sh.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifdef __cplusplus
 extern "C" {
#endif

void MFPlotfileToDX(FILE*,char*);

int main(int argc, char *argv[])
 {
  FILE *fid;
  char file[4096];

  if(argc<2){printf("Usage %s filename\n",argv[0]);fflush(stdout);return 8;}

  strcpy(file,argv[1]);
  strcat(file,".plotfile");
  fid=fopen(file,"r");
  if(fid==(FILE*)NULL){printf("Error opening file %s\n",file);fflush(stdout);return 8;}

  MFPlotfileToDX(fid,argv[1]);

  fclose(fid);

  return 0;
 }

#ifdef __cplusplus
}
#endif
