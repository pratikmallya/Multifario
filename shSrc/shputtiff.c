/* 
 *  %W%
 *  %D% %T%
 * 
 *  author: Mike Henderson mhender@watson.ibm.com
 */
#include <shInternal.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <multifarioConfig.h>

#ifdef HAVE_LIBTIFF
#include <tiffio.h>
#include <tiff.h>
#endif

#ifdef __cplusplus
 extern "C" {
#endif

void shputtiff(char *file, int nrows,int ncols,unsigned char *r,unsigned char *g,unsigned char *b)
 {
#ifdef HAVE_LIBTIFF
  TIFF *out = NULL;
  
  char *fn;
  char *buf;
  unsigned int i,j;
  unsigned red=0;
  unsigned green=1;
  unsigned blue=2;

  int rc;
  int ln;
  int separate;

  separate=0;

  ln=nrows*ncols;

  fn=(char *)malloc((strlen(file)+1)*sizeof(char));
  strcpy(fn,file);

/*printf("  puttiff for ==>%s<==\n",fn);
  printf("    Image is %dx%d\n",nrows,ncols);*/

  out = TIFFOpen(fn, "w");
  free(fn);
  if (out == NULL)return;

  TIFFSetField(out, TIFFTAG_IMAGEWIDTH, (unsigned long) ncols);
  TIFFSetField(out, TIFFTAG_IMAGELENGTH, (unsigned long) nrows);
  TIFFSetField(out, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
  TIFFSetField(out, TIFFTAG_SAMPLESPERPIXEL, 3);
  TIFFSetField(out, TIFFTAG_BITSPERSAMPLE, 8,8,8);
  TIFFSetField(out, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
  TIFFSetField(out, TIFFTAG_ROWSPERSTRIP, 1);
/*TIFFSetField(out, TIFFTAG_STRIPBYTECOUNTS, ncols);*/
  
  TIFFSetField(out, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);

  if(separate)
   {
    TIFFSetField(out, TIFFTAG_PLANARCONFIG, PLANARCONFIG_SEPARATE);
    buf=(char *)malloc(ncols*sizeof(char));
   }else{
    TIFFSetField(out, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    buf=(char *)malloc(3*ncols*sizeof(char));
   }
  TIFFSetField(out, TIFFTAG_SOFTWARE, "puttiff");

  for(i=0;i<nrows;i++)
   {
    if(separate)
     {
      for(j=0;j<ncols;j++)buf[j]=r[j+i*ncols];
/*    rc=TIFFWriteScanline(out, (unsigned char*)buf, (unsigned)i, red);*/
      rc=TIFFWriteScanline(out, (unsigned char*)buf, (unsigned)(nrows-i-1), red);

      for(j=0;j<ncols;j++)buf[j]=g[j+i*ncols];
/*    rc=TIFFWriteScanline(out, (unsigned char*)buf, (unsigned)i, green);*/
      rc=TIFFWriteScanline(out, (unsigned char*)buf, (unsigned)(nrows-i-1), green);

      for(j=0;j<ncols;j++)buf[j]=b[j+i*ncols];
/*    rc=TIFFWriteScanline(out, (unsigned char*)buf, (unsigned)i, blue);*/
      rc=TIFFWriteScanline(out, (unsigned char*)buf, (unsigned)(nrows-i-1), blue);
     }else{
      for(j=0;j<ncols;j++)buf[3*j  ]=r[j+i*ncols];
      for(j=0;j<ncols;j++)buf[3*j+1]=g[j+i*ncols];
      for(j=0;j<ncols;j++)buf[3*j+2]=b[j+i*ncols];
/*    rc=TIFFWriteScanline(out, (unsigned char*)buf, (unsigned)i, 0);*/
      rc=TIFFWriteScanline(out, (unsigned char*)buf, (unsigned)(nrows-i-1), 0);
     }
   }
  free(buf);

  TIFFWriteDirectory(out);

  TIFFClose(out);
#endif
  return;
 }

#ifdef __cplusplus
}
#endif
