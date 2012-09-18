/* 
    @(#)shdopg.c	1.5
    02/04/19 16:38:04
   
    PROGRAM NAME:  Manifold

    (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
    CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.

    Please refer to the LICENSE file in the top directory

*/
/*      author: Mike Henderson mhender@watson.ibm.com */

#include <shInternal.h>
#include <math.h>
#include <stdio.h>

int SHmaxpol=0;
float *SHx=(float*)NULL;
float *SHy=(float*)NULL;
float *SHz=(float*)NULL;
int   *SHi=(int*)NULL;
int   *SHj=(int*)NULL;
int   *SHk0=(int*)NULL;
int   *SHk1=(int*)NULL;
int   *SHke=(int*)NULL;
float *SHo=(float*)NULL;
float *SHs=(float*)NULL;
float *SHxo=(float*)NULL;
float *SHyo=(float*)NULL;
float *SHzo=(float*)NULL;
float *SHno=(float*)NULL;
float *SHxs=(float*)NULL;
float *SHys=(float*)NULL;
float *SHzs=(float*)NULL;
float *SHns=(float*)NULL;
int   *SHi0=(int*)NULL;
float *SHx0=(float*)NULL;
float *SHy0=(float*)NULL;
float *SHz0=(float*)NULL;
float *SHn0=(float*)NULL;
float *SHs0=(float*)NULL;
int   *SHiflg=(int*)NULL;

/* #define VERBOSE 1*/

void shdopg(int *nt,float *xt,float *yt,float *zt,float *n)
 {

/*   shade and z-buffer a polygon defined by the vertices xt,yt,zt.  */
/*    with specified normal                                          */

  float nn[3]={1.,0.,0.};
  float t;
  float de0,de1;

  float xx,yy,zz;
  float xxx=0.;
  float yyy=0.;
  float zzz=0.;
  float zbuf=0.;

  int r=0;
  int g=0;
  int b=0;
  int nv,imax,jmax,imin,jmin;
  int m,mi,mj;
  int mleft,mright;
  int kt;
  int nedges;
  float an;
  int imask=0;
  int ipix,jpix;

/*  Project the vertices into Perspective coordinates */
/*  Transform into raster coordinates                 */
/*   SHk0[m] points to one point of each edge,          */
/*   SHk1[m] points to the other                        */

#ifdef VERBOSE
  printf("  shdopg n=%d\n",*nt);
  for(m=0;m<*nt;m++)
   {
    printf("    %d (%f,%f,%f)  (%f,%f,%f)\n",m,xt[m],yt[m],zt[m],n[3*m],n[1+3*m],n[2+3*m]);
    fflush(stdout);
   }
#endif

  if(*nt>SHmaxpol)
   {
    SHmaxpol+=1000;
    SHx=(float*)realloc((void*)SHx,SHmaxpol*sizeof(float));
    if(SHx==(float*)NULL){fprintf(stderr,"shdopg: out of memory!");fflush(stderr);return;}
    SHy=(float*)realloc((void*)SHy,SHmaxpol*sizeof(float));
    if(SHy==(float*)NULL){fprintf(stderr,"shdopg: out of memory!");fflush(stderr);return;}
    SHz=(float*)realloc((void*)SHz,SHmaxpol*sizeof(float));
    if(SHz==(float*)NULL){fprintf(stderr,"shdopg: out of memory!");fflush(stderr);return;}
    SHi=(int*)realloc((void*)SHi,SHmaxpol*sizeof(int));
    if(SHi==(int*)NULL){fprintf(stderr,"shdopg: out of memory!");fflush(stderr);return;}
    SHj=(int*)realloc((void*)SHj,SHmaxpol*sizeof(int));
    if(SHj==(int*)NULL){fprintf(stderr,"shdopg: out of memory!");fflush(stderr);return;}
    SHk0=(int*)realloc((void*)SHk0,SHmaxpol*sizeof(int));
    if(SHk0==(int*)NULL){fprintf(stderr,"shdopg: out of memory!");fflush(stderr);return;}
    SHk1=(int*)realloc((void*)SHk1,SHmaxpol*sizeof(int));
    if(SHk1==(int*)NULL){fprintf(stderr,"shdopg: out of memory!");fflush(stderr);return;}
    SHke=(int*)realloc((void*)SHke,SHmaxpol*sizeof(int));
    if(SHke==(int*)NULL){fprintf(stderr,"shdopg: out of memory!");fflush(stderr);return;}
    SHo=(float*)realloc((void*)SHo,SHmaxpol*sizeof(float));
    if(SHo==(float*)NULL){fprintf(stderr,"shdopg: out of memory!");fflush(stderr);return;}
    SHs=(float*)realloc((void*)SHs,SHmaxpol*sizeof(float));
    if(SHs==(float*)NULL){fprintf(stderr,"shdopg: out of memory!");fflush(stderr);return;}
    SHxo=(float*)realloc((void*)SHxo,SHmaxpol*sizeof(float));
    if(SHxo==(float*)NULL){fprintf(stderr,"shdopg: out of memory!");fflush(stderr);return;}
    SHyo=(float*)realloc((void*)SHyo,SHmaxpol*sizeof(float));
    if(SHyo==(float*)NULL){fprintf(stderr,"shdopg: out of memory!");fflush(stderr);return;}
    SHzo=(float*)realloc((void*)SHzo,SHmaxpol*sizeof(float));
    if(SHzo==(float*)NULL){fprintf(stderr,"shdopg: out of memory!");fflush(stderr);return;}
    SHno=(float*)realloc((void*)SHno,(2+3*SHmaxpol)*sizeof(float));
    if(SHno==(float*)NULL){fprintf(stderr,"shdopg: out of memory!");fflush(stderr);return;}
    SHxs=(float*)realloc((void*)SHxs,SHmaxpol*sizeof(float));
    if(SHxs==(float*)NULL){fprintf(stderr,"shdopg: out of memory!");fflush(stderr);return;}
    SHys=(float*)realloc((void*)SHys,SHmaxpol*sizeof(float));
    if(SHys==(float*)NULL){fprintf(stderr,"shdopg: out of memory!");fflush(stderr);return;}
    SHzs=(float*)realloc((void*)SHzs,SHmaxpol*sizeof(float));
    if(SHzs==(float*)NULL){fprintf(stderr,"shdopg: out of memory!");fflush(stderr);return;}
    SHns=(float*)realloc((void*)SHns,(2+3*SHmaxpol)*sizeof(float));
    if(SHns==(float*)NULL){fprintf(stderr,"shdopg: out of memory!");fflush(stderr);return;}
    SHi0=(int*)realloc((void*)SHi0,SHmaxpol*sizeof(int));
    if(SHi0==(int*)NULL){fprintf(stderr,"shdopg: out of memory!");fflush(stderr);return;}
    SHx0=(float*)realloc((void*)SHx0,SHmaxpol*sizeof(float));
    if(SHx0==(float*)NULL){fprintf(stderr,"shdopg: out of memory!");fflush(stderr);return;}
    SHy0=(float*)realloc((void*)SHy0,SHmaxpol*sizeof(float));
    if(SHy0==(float*)NULL){fprintf(stderr,"shdopg: out of memory!");fflush(stderr);return;}
    SHz0=(float*)realloc((void*)SHz0,SHmaxpol*sizeof(float));
    if(SHz0==(float*)NULL){fprintf(stderr,"shdopg: out of memory!");fflush(stderr);return;}
    SHn0=(float*)realloc((void*)SHn0,(2+3*SHmaxpol)*sizeof(float));
    if(SHn0==(float*)NULL){fprintf(stderr,"shdopg: out of memory!");fflush(stderr);return;}
    SHs0=(float*)realloc((void*)SHs0,SHmaxpol*sizeof(float));
    if(SHs0==(float*)NULL){fprintf(stderr,"shdopg: out of memory!");fflush(stderr);return;}
    SHiflg=(int*)realloc((void*)SHiflg,SHmaxpol*sizeof(int));
    if(SHiflg==(int*)NULL){fprintf(stderr,"shdopg: out of memory!");fflush(stderr);return;}
   }

  nv=*nt;
  if(fabs(xt[0]-xt[nv-1])+fabs(yt[0]-yt[nv-1])+fabs(zt[0]-zt[nv-1])<3.e-7)nv--;

  imax=0;
  jmax=0;
  imin=shIMax;
  jmin=shJMax;
  for(m=0;m<nv;m++)
   {
    shpers(xt+m,yt+m,zt+m,SHx+m,SHy+m,SHz+m);
    SHi[m]=SHx[m]*shMax+.5;
    SHj[m]=SHy[m]*shMax+.5;
    SHk0[m]=m;
    SHk1[m]=m+1;
    if(m==nv-1)SHk1[m]=0;
    if(SHi[m]>imax)imax=SHi[m];
    if(SHj[m]>jmax)jmax=SHj[m];
    if(SHi[m]<imin)imin=SHi[m];
    if(SHj[m]<jmin)jmin=SHj[m];
   }

#ifdef VERBOSE
  printf("  shdopg n=%d\n",nv);
  for(m=0;m<nv;m++)
   {
    printf("    %d (%d,%d)\n",m,SHi[m],SHj[m]);
    fflush(stdout);
   }
#endif

/*  Quick clip if polygon is not in image */

  if((imax<0||imin>=shIMax||jmax<0||jmin>=shJMax))
   {
#ifdef VERBOSE
    printf("  done shdopg -- polygon out of range (%d,%d)->(%d,%d), (0,0)->(%d,%d)\n",imin,jmin,imax,jmax,shIMax,shJMax); 
#endif
    return;
   }

  if(jmax==jmin||nv<3)
   {
    shpl(&nv,xt,yt,zt);
#ifdef VERBOSE
    printf("  done shdopg -- too few vertices, or horiSHzontal line, nv=%d, j in (%d,%d)\n",nv,jmin,jmax); 
#endif
    return;
   }

/*  flip edges so that SHy[SHk0[m]]<SHy[SHk1[m]] */

  for(m=0;m<nv;m++)
   {
    if(SHj[SHk0[m]]>SHj[SHk1[m]])
     {
      kt=SHk0[m];
      SHk0[m]=SHk1[m];
      SHk1[m]=kt;
     }
   }

/*  sort edges by increasing SHy[SHk0[m]] */

  for(mi=0;mi<nv-1;mi++)
   {
    for(mj=mi+1;mj<nv;mj++)
     {
      if(SHj[SHk0[mj]]<SHj[SHk0[mi]])
       {
        kt=SHk0[mi];
        SHk0[mi]=SHk0[mj];
        SHk0[mj]=kt;
        kt=SHk1[mi];
        SHk1[mi]=SHk1[mj];
        SHk1[mj]=kt;
       }
     }
   }

#ifdef VERBOSE
  printf("  shdopg edges\n");
  for(m=0;m<nv;m++)
   {
    printf("    %d (%d,%d)  (%f,%f,%f)->(%f,%f,%f)\n",m,SHk0[m],SHk1[m],SHx[SHk0[m]],SHy[SHk0[m]],SHz[SHk0[m]],SHx[SHk1[m]],SHy[SHk1[m]],SHz[SHk1[m]]);
    fflush(stdout);
   }
#endif

/*  Normalize the normal vector at each vertex */

  for(m=0;m<nv;m++)
   {
    an=sqrt(n[3*m]*n[3*m]+n[1+3*m]*n[1+3*m]+n[2+3*m]*n[2+3*m]);
    if(fabs(an)<1.e-6)
     {
      shpl(&nv,xt,yt,zt);
#ifdef VERBOSE
      printf("  done shdopg -- flat line m=%d, n=(%f,%f,%f), |n|=%f\n",m,n[3*m],n[1+3*m],n[2+3*m],an);
#endif
      return;
     }else{
      an=1./an;
      n[3*m]=n[3*m]*an;
      n[1+3*m]=n[1+3*m]*an;
      n[2+3*m]=n[2+3*m]*an;
     }

/*  Compute the equations of the edges of the polygon */


/*   (x,y,z)=(SHxo[m],SHyo[m],SHzo[m])+j*(SHxs[m],SHys[m],SHzs[m]) */

  if(SHj[SHk1[m]]==SHj[SHk0[m]])
   {
    SHs[m]=0.;
    SHxs[m]=0.;
    SHys[m]=0.;
    SHzs[m]=0.;
    SHns[3*m]=0;
    SHns[1+3*m]=0;
    SHns[2+3*m]=0;
   }else{
    t=1./(float)(SHj[SHk1[m]]-SHj[SHk0[m]]);
    SHs[m]=     (SHi[SHk1[m]]-SHi[SHk0[m]])*t;
    SHxs[m]=    (SHx[SHk1[m]]-SHx[SHk0[m]])*t;
    SHys[m]=    (SHy[SHk1[m]]-SHy[SHk0[m]])*t;
    SHzs[m]=    (SHz[SHk1[m]]-SHz[SHk0[m]])*t;
    SHns[3*m]=  (n[  3*SHk1[m]]-n[  3*SHk0[m]])*t;
    SHns[1+3*m]=(n[1+3*SHk1[m]]-n[1+3*SHk0[m]])*t;
    SHns[2+3*m]=(n[2+3*SHk1[m]]-n[2+3*SHk0[m]])*t;
   }
  SHo[m]=     SHi[SHk0[m]]-SHj[SHk0[m]]*SHs[m];
  SHxo[m]=    SHx[SHk0[m]]-SHj[SHk0[m]]*SHxs[m];
  SHyo[m]=    SHy[SHk0[m]]-SHj[SHk0[m]]*SHys[m];
  SHzo[m]=    SHz[SHk0[m]]-SHj[SHk0[m]]*SHzs[m];
  SHno[3*m]=  n[  3*SHk0[m]]-SHj[SHk0[m]]*SHns[3*m];
  SHno[1+3*m]=n[1+3*SHk0[m]]-SHj[SHk0[m]]*SHns[1+3*m];
  SHno[2+3*m]=n[2+3*SHk0[m]]-SHj[SHk0[m]]*SHns[2+3*m];
 }

#ifdef VERBOSE
  printf("  shdopg equations\n");
  for(m=0;m<nv;m++)
   {
    printf("    %4.4d s %f,  SHxs (%f,%f,%f), SHns (%f,%f,%f)\n",m,SHs[m],SHxs[m],SHys[m],SHzs[m],SHns[3*m],SHns[1+3*m],SHns[2+3*m]);
    printf("         o %f,  SHxo (%f,%f,%f), SHno (%f,%f,%f)\n",SHo[m],SHxo[m],SHyo[m],SHzo[m],SHno[3*m],SHno[1+3*m],SHno[2+3*m]);
    fflush(stdout);
   }
#endif

/*  Now fill the polygon, bottom to top */

 for(jpix=jmin;jpix<=jmax;jpix++)
  {

/*  Find the range of edges in k which cross this value of jpix */

    mleft=0;
    mright=nv-1;
    for(m=0;m<nv;m++)
     {
      if(SHj[SHk1[mleft]]<jpix)mleft++;
      if(SHj[SHk0[mright]]>jpix)mright--;
     }

/*    Compute the intersection points on each edge */

#ifdef VERBOSE
    printf("scan line %d, left %d, right %d\n",jpix,mleft,mright);
    fflush(stdout);
#endif

    nedges=0;
    for(m=mleft;m<=mright;m++)
     {
      if(SHj[SHk1[m]]>=jpix)
       {
        if(SHj[SHk1[m]]==jpix&&SHj[SHk0[m]]==jpix)
         {
          SHi0[nedges]= SHi[SHk0[m]];
          SHiflg[nedges]=3;
          SHs0[nedges]=fabs( SHi[SHk1[m]] - SHi[SHk0[m]] )+.5;
          SHx0[nedges]=SHx[SHk0[m]];
          SHy0[nedges]=SHy[SHk0[m]];
          SHz0[nedges]=SHz[SHk0[m]];
          SHn0[3*nedges]=n[3*SHk0[m]];
          SHn0[1+3*nedges]=n[1+3*SHk0[m]];
          SHn0[2+3*nedges]=n[2+3*SHk0[m]];
          SHke[nedges]=nedges;
          nedges++;
  
          SHi0[nedges]= SHi[SHk1[m]];
          SHiflg[nedges]=3;
          SHs0[nedges]=fabs( SHi[SHk1[m]] - SHi[SHk0[m]] )+.5;
          SHx0[nedges]=SHx[SHk1[m]];
          SHy0[nedges]=SHy[SHk1[m]];
          SHz0[nedges]=SHz[SHk1[m]];
          SHn0[3*nedges]=n[3*SHk1[m]];
          SHn0[1+3*nedges]=n[1+3*SHk1[m]];
          SHn0[2+3*nedges]=n[2+3*SHk1[m]];
          SHke[nedges]=nedges;
          nedges++;
         }else{
          SHiflg[nedges]=0;
          if(SHj[SHk1[m]]==jpix)SHiflg[nedges]=1;
          if(SHj[SHk0[m]]==jpix)SHiflg[nedges]=2;
          if(SHo[m]+jpix*SHs[m]>=0.)
            SHi0[nedges]= SHo[m]+jpix* SHs[m]+.5;
           else
            SHi0[nedges]= SHo[m]+jpix* SHs[m]-.5;
          if(SHj[SHk1[m]]==SHj[SHk0[m]])
           {
            fprintf(stderr," divide by zero in shpgnrm, SHj[%d]=SHj[%d]=%d, jpix=%d\n",SHk1[m],SHk0[m],SHj[SHk0[m]],jpix);
            abort();
           }
          SHs0[nedges]=fabs((SHi[SHk1[m]]-SHi[SHk0[m]])/(float)(SHj[SHk1[m]]-SHj[SHk0[m]]));
          SHx0[nedges]=SHxo[m]+jpix*SHxs[m];
          SHy0[nedges]=SHyo[m]+jpix*SHys[m];
          SHz0[nedges]=SHzo[m]+jpix*SHzs[m];
          SHn0[3*nedges]=SHno[3*m]+jpix*SHns[3*m];
          SHn0[1+3*nedges]=SHno[1+3*m]+jpix*SHns[1+3*m];
          SHn0[2+3*nedges]=SHno[2+3*m]+jpix*SHns[2+3*m];
          SHke[nedges]=nedges;
          nedges++;
         }
       }
     }

/*    Sort SHke by increasing SHi0 (intersection points left to right.) */

    for(mi=0;mi<nedges;mi++)
     {
      for(mj=mi+1;mj<nedges;mj++)
       {
        if(SHi0[SHke[mj]]<SHi0[SHke[mi]])
         {
          kt=SHke[mi];
          SHke[mi]=SHke[mj];
          SHke[mj]=kt;
         }
        if(SHi0[SHke[mj]]==SHi0[SHke[mi]]&&SHiflg[SHke[mj]]<SHiflg[SHke[mi]])
         {
          kt=SHke[mi];
          SHke[mi]=SHke[mj];
          SHke[mj]=kt;
         }
       }
     }

/*    Remove duplicates */

    mi=0;
    while(mi<nedges-1)
     {
      if(SHi0[SHke[mi]]==SHi0[SHke[mi+1]]&&SHiflg[SHke[mi]]!=0&&SHiflg[SHke[mi+1]]!=0&&SHiflg[SHke[mi]]!=SHiflg[SHke[mi+1]])
       {
        if(SHs0[SHke[mi+1]]>SHs0[SHke[mi]])SHs0[SHke[mi]]=SHs0[SHke[mi+1]];
        SHs0[SHke[mi+1]]=SHs0[SHke[mi]];
        for(mj=mi+1;mj<nedges;mj++)SHke[mj]=SHke[mj+1];
        nedges--;
       }else{
        mi++;
       }
     }

/*   Fill the Slabs */

    for(m=0;m<nedges;m+=2)
     {
      de0=.5+SHs0[SHke[m]];
      de1=.5+SHs0[SHke[m+1]];
      for(ipix=SHi0[SHke[m]];ipix<=SHi0[SHke[m+1]];ipix++)
       {
        if(ipix>-1&&ipix<shIMax&&jpix>-1&&jpix<shJMax)
         {

/*   Compute perspective coordinates of the pixel being filled (xx,yy,zz) */

          if(SHi0[SHke[m]]==SHi0[SHke[m+1]])
           {
            xx=SHx0[SHke[m]];
            yy=SHy0[SHke[m]];
            zz=SHz0[SHke[m]];
            nn[0]=SHn0[3*SHke[m]];
            nn[1]=SHn0[1+3*SHke[m]];
            nn[2]=SHn0[2+3*SHke[m]];
           }else{
            t=(float)(ipix-SHi0[SHke[m]])/(float)(SHi0[SHke[m+1]]-SHi0[SHke[m]]);
            xx=SHx0[SHke[m]]+t*(SHx0[SHke[m+1]]-SHx0[SHke[m]]);
            yy=SHy0[SHke[m]]+t*(SHy0[SHke[m+1]]-SHy0[SHke[m]]);
            zz=SHz0[SHke[m]]+t*(SHz0[SHke[m+1]]-SHz0[SHke[m]]);
            nn[0]=SHn0[3*SHke[m]]+t*(SHn0[3*SHke[m+1]]-SHn0[3*SHke[m]]);
            nn[1]=SHn0[1+3*SHke[m]]+t*(SHn0[1+3*SHke[m+1]]-SHn0[1+3*SHke[m]]);
            nn[2]=SHn0[2+3*SHke[m]]+t*(SHn0[2+3*SHke[m+1]]-SHn0[2+3*SHke[m]]);
           }
          an=sqrt(nn[0]*nn[0]+nn[1]*nn[1]+nn[2]*nn[2]);
          if(fabs(an)<1.e-7)
           {
            nn[0]=SHn0[3*SHke[m]];
            nn[1]=SHn0[1+3*SHke[m]];
            nn[2]=SHn0[2+3*SHke[m]];
            an=sqrt(nn[0]*nn[0]+nn[1]*nn[1]+nn[2]*nn[2]);
            nn[0]=nn[0]*an;
            nn[1]=nn[1]*an;
            nn[2]=nn[2]*an;
           }else{
            an=1./an;
            nn[0]=nn[0]*an;
            nn[1]=nn[1]*an;
            nn[2]=nn[2]*an;
           }

          shunpers(&xx,&yy,&zz,&xxx,&yyy,&zzz);

/*    Check against the Z-buffer */

          shgetz(&ipix,&jpix,&zbuf);
#ifdef VERBOSE
          printf("shgetz %d,%d %f %f\n",ipix,jpix,zbuf,zz);fflush(stdout);
#endif
          if(zbuf>zz)
           {

/*      If pixel is visible, shade it and put it in the image */

            shcolor(&xxx,&yyy,&zzz,nn,&(sh_rd[0]),&(sh_gd[0]),&(sh_bd[0]),&(sh_rd[2]),&(sh_gd[2]),&(sh_bd[2]),&r,&g,&b);
#ifdef VERBOSE
            printf("shcolor, x-( %f %f %f ) n-( %f %f %f ) fc-(%d %d %d ) bc-(%d %d %d) c %d %d %d\n",xxx,yyy,zzz,nn[0],nn[1],nn[2],sh_rd[0],sh_gd[0],sh_bd[0],sh_rd[2],sh_gd[2],sh_bd[2],r,g,b);
#endif

            if(shMask)
             {
              sh3dmask(&xxx,&yyy,&zzz,nn,&imask);
              if(imask==1)shsetp(&ipix,&jpix,&r,&g,&b,&zz);
             }else{
              shsetp(&ipix,&jpix,&r,&g,&b,&zz);
             }
           }
         }
       }
     }
   }

#ifdef VERBOSE
  printf("  done shdopg -- normal end\n");
#endif
  return;
 }
