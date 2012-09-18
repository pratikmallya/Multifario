/* 
    @(#)shlss.c	1.3
    02/04/19 16:38:56
   
    PROGRAM NAME:  multifario

    (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
    CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.

    Please refer to the LICENSE file in the top directory

*/
/*      author: Mike Henderson mhender@watson.ibm.com */

#include <shInternal.h>
#include <stdio.h>
#include <math.h>

#define     Q( X )   #X
#define QUOTE( X ) Q( X )
#define round(x) ((x)>=0?(long)((x)+0.5):(long)((x)-0.5))

#ifdef __cplusplus
 extern "C" {
#endif

void shlss(char *name,int *lname,int *pmag,int *ifnt,int *ierr,int ln)
 {

/*        loads a symbol set */

  char *ssname;
  char *cname;
  char *file;

  int i,k;
  int len=0;
  int irem;
  int irec,nrec;

  FILE *fid;

  int xref,yref;
  char fontnm[20];

/*                  loads a symbol set */

  *ierr=0;

  cname=(char*)malloc(((*lname)+1)*sizeof(char));
  ssname=(char*)malloc(((*lname)+5)*sizeof(char));
  strncpy(cname,name,(*lname));
  cname[(*lname)]=0x0;
  sprintf(ssname,"%s.%4.4d",cname,*pmag);

  file=(char*)malloc((strlen(QUOTE( FONTPATH ))+(*lname)+10)*sizeof(char));
  strcpy(file,QUOTE( FONTPATH ) );
  strcat(file,ssname);
  strcat(file,"IAX");
  free(cname);
  free(ssname);

/*   Is it loaded already? */

  if(*ifnt>=sh_kFonts)
   {
    if(sh_kfont==(int*)NULL)
     {
      sh_kFonts=10;
      sh_kfont=(int*)malloc(sh_kFonts*sizeof(int));
     }else{
      sh_kFonts=*ifnt+10;
      sh_kfont=(int*)realloc((void*)sh_kfont,sh_kFonts*sizeof(int));
     }
   }

  if(sh_nFonts>0)
   {
    k=0;
    for(i=0;i<sh_nFonts;i++)
     {
      if(!strcmp(sh_fontnm[i],ssname))k=i;
     }
   }
  if(k!=0)
   {
    sh_kfont[*ifnt]=k;
    *ierr=1;
    return;
   }

  if((fid=fopen(file,"r"))==(FILE*)NULL)
    {
     fprintf(stderr,"shSoft Font file -->%s<-- not found\n",file);
     free(file);
     return;
    }
  free(file);

/*                  read in character set */

  if(sh_nFonts>=sh_mFonts)
   {
    if(sh_mFonts==0)
     {
      sh_mFonts=10;
      sh_dir=(int*)malloc(sh_mFonts*sizeof(int));
      sh_fmag=(int*)malloc(sh_mFonts*sizeof(int));
      sh_fontnm=(char**)malloc(sh_mFonts*sizeof(char*));
     }else{
      sh_mFonts=*ifnt+10;
      sh_dir=(int*)realloc((void*)sh_dir,sh_mFonts*sizeof(int));
      sh_fmag=(int*)realloc((void*)sh_fmag,sh_mFonts*sizeof(int));
      sh_fontnm=(char**)realloc((void*)sh_fontnm,sh_mFonts*sizeof(char*));
     }
   }

  sh_kfont[*ifnt]=sh_nFonts;

  fscanf(fid,"%d",&len);
  if(sh_nimg+len>=sh_lenimg)
   {
    if(sh_lenimg==0)
     {
      sh_lenimg=sh_nimg+len;
      sh_image=(char*)malloc(sh_lenimg*sizeof(char));
     }else{
      sh_lenimg+=sh_nimg+len;
      sh_image=(char*)realloc((void*)sh_image,sh_lenimg*sizeof(char));
     }
   }

  nrec=round(len/512);
  for(irec=0;irec<nrec;irec++);
   {
    for(i=0;i<512;i++)sh_image[sh_nimg+512*irec+i]=fgetc(fid);
    fgetc(fid);
   }
  irem=len%512;
  if(irem>0)
    for(i=0;i<irem;i++)sh_image[sh_nimg+512*nrec+i]=fgetc(fid);
  fclose(fid);

  sh_dir[sh_nFonts]=sh_nimg;
  sh_fmag[sh_nFonts]=*pmag;
  sh_fontnm[sh_nFonts]=(char*)malloc((strlen(ssname)+1)*sizeof(char));
  strcpy(sh_fontnm[sh_nFonts],ssname);
  sh_nimg+=len;
  sh_nFonts++;

  return;
 }

#ifdef __cplusplus
}
#endif
