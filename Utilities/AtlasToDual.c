/* 
    @(#)ComputeDual.c	1.3
    02/04/19 14:47:04
   
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
#include <MFEnumPolytope.h>
#include <MFEnumDualPolytope.h>
#include <MFDX.h>
#include <math.h>

#ifdef __cplusplus
 extern "C" {
#endif

int main(int argc, char *argv[])
 {
  int i,n;
  int kkk;
  float xmin,xmax,dist,alpha,beta;
  int gray;
  MFAtlas S;
  MFEnumDualPolytope P;
  FILE *fid;
  MFErrorHandler e;

  if(argc<2)
   {
    printf("Usage:\n     %s filename\n Opens file \"filename.atlas\"s\n",argv[0]);
    fflush(stdout);
    return 12;
   }

  e=MFCreateErrorHandler();

  fid=fopen(argv[1],"r");
  S=MFReadAtlas(fid,e);
  fclose(fid);

  P=MFEnumDualOfAtlas(S,e);
  MFDualPolytopeToDXFile("Pendula.dx",P,e);
  MFFreeEnumDualPolytope(P,e);

  MFFreeAtlas(S,e);

  MFFreeErrorHandler(e);

  return(0);
 }

#ifdef __cplusplus
}
#endif
