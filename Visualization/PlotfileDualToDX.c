#include <stdio.h>
#include <MFAtlas.h>
#include <MFDraw.h>
#include <sh.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <MF2dCellComplex.h>

#ifdef __cplusplus
 extern "C" {
#endif

void MFPlotfileDualToDX(FILE *fid,char *name);

int main(int argc, char *argv[])
 {
  FILE *fid;
  char file[4096];
  MF2dCellComplex C;
  MFErrorHandler e;

  if(argc<2){printf("Usage %s filename\n",argv[0]);fflush(stdout);return 8;}

  e=MFCreateErrorHandler();

  strcpy(file,argv[1]);
  strcat(file,".plotfile");
  fid=fopen(file,"r");
  if(fid==(FILE*)NULL){printf("Error opening file %s\n",file);fflush(stdout);return 8;}

  strcpy(file,argv[1]);
  strcat(file,"Dual");
/*MFPlotfileDualToDX(fid,file);*/


  C=MFPlotfileDual(fid,e);
  strcpy(file,argv[1]);
  strcat(file,"Dual");
  MFOutput2dCellComplexAsDX(C,file,e);
  MFFree2dCellComplex(C,e);

  fclose(fid);
  MFFreeErrorHandler(e);

  return 0;
 }

#ifdef __cplusplus
}
#endif
