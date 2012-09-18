/* 
    @(#)shputps.c	1.11
    02/04/19 16:40:23
   
    PROGRAM NAME:  multifario

    (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
    CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.

    Please refer to the LICENSE file in the top directory

*/
/*      author: Mike Henderson mhender@watson.ibm.com */

#include <shInternal.h>
#include <stdio.h>

#ifdef __cplusplus
 extern "C" {
#endif

void shputps(char *file,int *ln,int *m,int *n,unsigned char *imageR,unsigned char *imageG,unsigned char *imageB,int fln)
 {
  char title[255]="";
  FILE *fid;
  int i,j;
  int i0,j0,i1,j1;
  int xll,yll,xur,yur;
  float w;

  i0=*m;
  j0=*n;
  i1=-1;
  j1=-1;
  for(j=0;j<*n;j++)
   {
    for(i=0;i<*m;i++)
     {
      if(shRedBuffer[i+j*shIMax]!=255||shGreenBuffer[i+j*shIMax]!=255||shBlueBuffer[i+j*shIMax]!=255)
       {
        if(j<=j0)j0=j;
        if(i<=i0)i0=i;
        if(j>=j1)j1=j;
        if(i>=i1)i1=i;
       }
     }
   }
  w=(*n)*504./(*m);
  xll=i0*w/(*m) >= 0 ? (int)(i0*w/(*m)+0.5) : (int)(i0*w/(*m)-0.5);
  yll=j0*w/(*m) >= 0 ? (int)(j0*w/(*m)+0.5) : (int)(j0*w/(*m)-0.5);
  xur=i1*w/(*m) >= 0 ? (int)(i1*w/(*m)+0.5) : (int)(i1*w/(*m)-0.5);
  yur=j1*w/(*m) >= 0 ? (int)(j1*w/(*m)+0.5) : (int)(j1*w/(*m)-0.5);
      
  strncpy(title,file,*ln);
  title[*ln]=0x0;

  printf("Writing image to the postscript file %s\n",file);fflush(stdout);

  fid=fopen(file,"w");
  if(fid==(FILE*)NULL)
   {
    fprintf(stderr,"shputps: Problem opening file -->%s<-- for writing.\n",file);
    fflush(stderr);
    return;
   }
  fprintf(fid,"%s\n","%!PS-Adobe-3.0");
  fprintf(fid,"%s\n","%%Creator: SH by Michael E. Henderson, 03/17/93");
  fprintf(fid,"%s%s\n","%%Title: ",title);
  fprintf(fid,"%s %d %d %d %d\n","%%BoundingBox:",xll+50,yll+50,xur+52,yur+52);
  fprintf(fid,"%s\n","%%EndComments");
  fprintf(fid,"%s\n","%%BeginProlog");
  fprintf(fid,"%s\n","%%BeginResource: ");
  fprintf(fid,"%s\n"," ");
  fprintf(fid,"%s\n","% see if we have the ""colorimage"" operator.");
  fprintf(fid,"%s\n","% define one if we don""t");
  fprintf(fid,"%s\n","/colorimage where   % do we know about ""colorimage""?");
  fprintf(fid,"%s\n","  { pop }           % yes: pop off the ""dict"" returned");
  fprintf(fid,"%s\n","  {                 % no:  define one");
  fprintf(fid,"%s\n","    /str1 1 string def");
  fprintf(fid,"%s\n","    /str3 3 string def");
  fprintf(fid,"%s\n","    /colorimage");
  fprintf(fid,"%s\n","      { pop pop     % pop off ""false"", ""3"" operands");
  fprintf(fid,"%s\n","        pop         % pop off old ""readhexstring"" proc");
  fprintf(fid,"%s\n","                    % and define a new one for ""image""");
  fprintf(fid,"%s\n","        { currentfile str3 readhexstring pop pop");
  fprintf(fid,"%s\n","          str1 0    % for the ""put"" below");
  fprintf(fid,"%s\n","          str3 0 get 20 mul    % Red");
  fprintf(fid,"%s\n","          str3 1 get 32 mul    % Green");
  fprintf(fid,"%s\n","          str3 2 get 12 mul    % Blue");
  fprintf(fid,"%s\n","          add add 64 idiv      % I = .5G + .31R + .18B");
  fprintf(fid,"%s\n","          put str1  % str1 = intensity of R,G,B triplet");
  fprintf(fid,"%s\n","        } image");
  fprintf(fid,"%s\n","      } def         % end of colorimage def");
  fprintf(fid,"%s\n","  } ifelse          % end of ""false"" case");
  fprintf(fid,"%s\n"," ");
  fprintf(fid,"%s\n","%%EndResource");
  fprintf(fid,"%s\n","%%EndProlog");
  fprintf(fid,"%s\n"," ");
  fprintf(fid,"%s\n","%%Page: 1 1");
  fprintf(fid,"%s\n"," ");
  fprintf(fid,"%s\n","%%BeginPageSetup");
  fprintf(fid,"%s\n"," ");
  fprintf(fid," newpath\n");
  fprintf(fid,"    %d %d moveto\n",xll+50,yll+50);
  fprintf(fid,"    %d %d lineto\n",xur+52,yll+50);
  fprintf(fid,"    %d %d lineto\n",xur+52,yur+52);
  fprintf(fid,"    %d %d lineto\n",xll+50,yur+52);
  fprintf(fid,"    closepath clip\n");
  fprintf(fid,"%s\n","% Image size (1/72 inch coords,");
  fprintf(fid,"%d %f %s\n",504,w," scale");

  fprintf(fid,"%s\n"," ");
  fprintf(fid,"%s\n","%%EndPageSetup");
  fprintf(fid,"%s\n"," ");
  fprintf(fid,"%f %f %s\n",.1,.1," translate");
  fprintf(fid,"%s %d %s\n","/readstr ",3*(*m)," string def");
  fprintf(fid,"%d %d %d %s %d %d %d %d %d %d %s\n",*m,*n,8,"[ ", *m, 0,0,*n,0,0, " ]");
  fprintf(fid,"%s\n","{currentfile readstr readhexstring pop}");
  fprintf(fid,"%s\n","false 3 colorimage");

  for(j=0;j<*n;j++)
   {
    for(i=0;i<*m;i++)
     {
      fprintf(fid,"%2.2x%2.2x%2.2x",shRedBuffer[i+j*shIMax],shGreenBuffer[i+j*shIMax],shBlueBuffer[i+j*shIMax]);
      if(i%20==19)fprintf(fid,"\n");
     }
    if((*m-1)%20!=19)fprintf(fid,"\n");
   }

  fprintf(fid,"%s\n"," ");
  fprintf(fid,"%f %f %s\n",1./504.,1./w," scale");
  shdoStrings(fid,shIMax/512.,shJMax/w);

  fprintf(fid,"%s\n"," showpage");
  fclose(fid);

  return;
 }

#ifdef __cplusplus
}
#endif
