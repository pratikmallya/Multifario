/* 
    @(#)shclippg.c	1.3
    02/04/19 16:37:41
   
    PROGRAM NAME:  multifario

    (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
    CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.

    Please refer to the LICENSE file in the top directory

*/
/*      author: Mike Henderson mhender@watson.ibm.com */

#include <shInternal.h>
#include <math.h>
#include <stdio.h>

#ifdef __cplusplus
 extern "C" {
#endif

#define maxpol 1000

void shclippg(int *n0,float *x0,float *y0,float *z0,float *nrm0,int *ic0,int *n1,float **x1,float **y1,float **z1,float **nrm1,int **ic1)
 {

/*     clip a polygon */

  float *x2,*y2,*z2;
  float *nrm2;
  int *ic2;
  int in0,in1;
  int m1,m2;
  int i,j,ip;
  float direc;
  float t;
  int i0;
  int verbose;

  verbose=0;

  m2=maxpol;
  x2=(float*)malloc(m2*sizeof(float));
  y2=(float*)malloc(m2*sizeof(float));
  z2=(float*)malloc(m2*sizeof(float));
  nrm2=(float*)malloc(3*m2*sizeof(float));
  ic2=(int*)malloc(m2*sizeof(int));

  m1=(*n0)+1;
  if(*x1==(float*)NULL)*x1=(float*)malloc(m1*sizeof(float));
  if(*y1==(float*)NULL)*y1=(float*)malloc(m1*sizeof(float));
  if(*z1==(float*)NULL)*z1=(float*)malloc(m1*sizeof(float));
  if(*nrm1==(float*)NULL)*nrm1=(float*)malloc(3*m1*sizeof(float));
  if(*ic1==(int*)NULL)*ic1=(int*)malloc(m1*sizeof(int));

/* Clipping Planes */

/*    NPLNS     Number of clipping planes */
/*    PLNO[2+3*M) point on plane M */
/*    PLNN[2+3*M) normal to plane M */
/*    IPLN(M)   side of plane M to clip (IPLN*PLNN<0 ==> clipped) */
/*    oper(M)   whether clipping is ored or anded. */


/* copy polygon 0 into 1, removing near duplicates. */

  *n1=0;
  (*x1)[*n1]=x0[0];
  (*y1)[*n1]=y0[0];
  (*z1)[*n1]=z0[0];
  (*nrm1)[0+3*(*n1)]=nrm0[0];
  (*nrm1)[1+3*(*n1)]=nrm0[1];
  (*nrm1)[2+3*(*n1)]=nrm0[2];
  (*ic1)[*n1]=ic0[0];
  (*n1)++;
  for(i=1;i<*n0;i++)
   {
    if(fabs(x0[i]-x0[i-1])+fabs(y0[i]-y0[i-1])+fabs(z0[i]-z0[i-1])>1.e-6)
     {
      (*x1)[*n1]=x0[i];
      (*y1)[*n1]=y0[i];
      (*z1)[*n1]=z0[i];
      (*nrm1)[  3*(*n1)]=nrm0[  3*i];
      (*nrm1)[1+3*(*n1)]=nrm0[1+3*i];
      (*nrm1)[2+3*(*n1)]=nrm0[2+3*i];
      (*ic1)[(*n1)]=ic0[i];
      (*n1)++;
     }
   }

/* Duplicate first point */

  (*x1)[*n1]=x0[0];
  (*y1)[*n1]=y0[0];
  (*z1)[*n1]=z0[0];
  (*nrm1)[  3*(*n1)]=nrm0[  3*0];
  (*nrm1)[1+3*(*n1)]=nrm0[1+3*0];
  (*nrm1)[2+3*(*n1)]=nrm0[2+3*0];
  (*ic1)[(*n1)]=ic0[0];
  (*n1)++;

/* Loop over the clipping planes */

/*    in0=this point visible, in1=next point visible */

  if(sh_nplns==0)
   {
    free(x2);
    free(y2);
    free(z2);
    free(nrm2);
    free(ic2);
    return;
   }

  for(ip=0;ip<sh_nplns;ip++)
   {
    if(verbose)printf("plane %d is ((x,y,z)-(%f,%f,%f)).(%f,%f,%f)*%d>0\n",ip,sh_plno[  3*ip],sh_plno[1+3*ip],sh_plno[2+3*ip],sh_plnn[  3*ip],sh_plnn[1+3*ip],sh_plnn[2+3*ip],sh_ipln[ip]);
    in0=0;
    i0=0;
    for(i=0;i<*n1;i++)
     {
      if(verbose){printf("vertex %d(%d) (%f,%f,%f) is ",i,*n1,(*x1)[i],(*y1)[i],(*z1)[i]);fflush(stdout);}
      direc=sh_plnn[  3*ip]*((*x1)[i]-sh_plno[  3*ip])+sh_plnn[1+3*ip]*((*y1)[i]-sh_plno[1+3*ip])+sh_plnn[2+3*ip]*((*z1)[i]-sh_plno[2+3*ip]);
      in1=sh_ipln[ip]*direc>0;
      if(verbose){if(in1)printf("in\n");else printf("out\n");}
      if(!in0)i0=i;
      in0=in0||in1;
     }
    if(verbose)fflush(stdout);
    if(!in0)
     {
      if(verbose){printf("Polygon completely outside\n\n");fflush(stdout);}
      *n1=0;
      free(x2);
      free(y2);
      free(z2);
      free(nrm2);
      free(ic2);
      return;
     }

    if(i0>0)i0--;

    if(verbose){printf("Starting at vertex %d\n",i0);fflush(stdout);}
#define JUNK
#ifdef JUNK
    j=0;
    direc=sh_plnn[  3*ip]*((*x1)[i0]-sh_plno[  3*ip])+sh_plnn[1+3*ip]*((*y1)[i0]-sh_plno[1+3*ip])+sh_plnn[2+3*ip]*((*z1)[i0]-sh_plno[2+3*ip]);
    in0=sh_ipln[ip]*direc>0;
    i=i0;
    while(i+1<*n1)
     {
      if(verbose){printf("Vertex %d(%d) (%f,%f,%f) next is %d (%f,%f,%f)\n",i,*n1,(*x1)[i],(*y1)[i],(*z1)[i],i+1,(*x1)[i+1],(*y1)[i+1],(*z1)[i+1]);fflush(stdout);}
      direc=sh_plnn[  3*ip]*((*x1)[i+1]-sh_plno[  3*ip])
           +sh_plnn[1+3*ip]*((*y1)[i+1]-sh_plno[1+3*ip])
           +sh_plnn[2+3*ip]*((*z1)[i+1]-sh_plno[2+3*ip]);
      in1=sh_ipln[ip]*direc>0;
      if(verbose)
       {
        if(in0){printf("   %d is in\n",i);fflush(stdout);}
        if(!in0){printf("   %d is out\n",i);fflush(stdout);}
        if(in1){printf("   %d is in\n",i+1);fflush(stdout);}
        if(!in1){printf("   %d is out\n",i+1);fflush(stdout);}
       }
      if(!in0&&!in1)
       {
        if(verbose){printf("   !in0&&!in1\n");fflush(stdout);}
        in0=in1;
        i++;
       }else if(in0&&in1)
       {
        if(verbose){printf("   in0&&in1\n");fflush(stdout);}
        in0=in1;
        if(j+1>m2)
         {
          m2+=5;
          x2=(float*)realloc((void*)x2,m2*sizeof(float));
          y2=(float*)realloc((void*)y2,m2*sizeof(float));
          z2=(float*)realloc((void*)z2,m2*sizeof(float));
          nrm2=(float*)realloc((void*)nrm2,3*m2*sizeof(float));
          ic2=(int*)realloc((void*)ic2,m2*sizeof(int));
         }
        x2[j]=(*x1)[i];
        y2[j]=(*y1)[i];
        z2[j]=(*z1)[i];
        nrm2[  3*j]=(*nrm1)[  3*i];
        nrm2[1+3*j]=(*nrm1)[1+3*i];
        nrm2[2+3*j]=(*nrm1)[2+3*i];
        ic2[j]=(*ic1)[i];
        j++;

        i++;
       }else if(!in0&&in1)
       {
        if(verbose){printf("   !in0&&in1\n");fflush(stdout);}
        t=-( sh_plnn[  3*ip]*((*x1)[i]-sh_plno[  3*ip])
            +sh_plnn[1+3*ip]*((*y1)[i]-sh_plno[1+3*ip])
            +sh_plnn[2+3*ip]*((*z1)[i]-sh_plno[2+3*ip]) )
          /( sh_plnn[  3*ip]*((*x1)[i+1]-(*x1)[i])
            +sh_plnn[1+3*ip]*((*y1)[i+1]-(*y1)[i])
            +sh_plnn[2+3*ip]*((*z1)[i+1]-(*z1)[i]) );
        (*x1)[i]=(*x1)[i]+t*((*x1)[i+1]-(*x1)[i]);
        (*y1)[i]=(*y1)[i]+t*((*y1)[i+1]-(*y1)[i]);
        (*z1)[i]=(*z1)[i]+t*((*z1)[i+1]-(*z1)[i]);
        (*nrm1)[  3*i]=(*nrm1)[  3*i]+t*((*nrm1)[  3*(i+1)]-(*nrm1)[  3*i]);
        (*nrm1)[1+3*i]=(*nrm1)[1+3*i]+t*((*nrm1)[1+3*(i+1)]-(*nrm1)[1+3*i]);
        (*nrm1)[2+3*i]=(*nrm1)[2+3*i]+t*((*nrm1)[2+3*(i+1)]-(*nrm1)[2+3*i]);
        (*ic1)[i]=-1;
        if(j+1>m2)
         {
          m2+=5;
          x2=(float*)realloc((void*)x2,m2*sizeof(float));
          y2=(float*)realloc((void*)y2,m2*sizeof(float));
          z2=(float*)realloc((void*)z2,m2*sizeof(float));
          nrm2=(float*)realloc((void*)nrm2,3*m2*sizeof(float));
          ic2=(int*)realloc((void*)ic2,m2*sizeof(int));
         }
        x2[j]=(*x1)[i];
        y2[j]=(*y1)[i];
        z2[j]=(*z1)[i];
        nrm2[  3*j]=(*nrm1)[  3*i];
        nrm2[1+3*j]=(*nrm1)[1+3*i];
        nrm2[2+3*j]=(*nrm1)[2+3*i];
        ic2[j]=(*ic1)[i];
        j++;
        in0=1;
        i++;
       }else if(in0&&!in1)
       {
        if(verbose){printf("   in0&&!in1\n");fflush(stdout);}
        t=-( sh_plnn[  3*ip]*((*x1)[i]-sh_plno[  3*ip])
            +sh_plnn[1+3*ip]*((*y1)[i]-sh_plno[1+3*ip])
            +sh_plnn[2+3*ip]*((*z1)[i]-sh_plno[2+3*ip]) )
          /( sh_plnn[  3*ip]*((*x1)[i+1]-(*x1)[i])
            +sh_plnn[1+3*ip]*((*y1)[i+1]-(*y1)[i])
            +sh_plnn[2+3*ip]*((*z1)[i+1]-(*z1)[i]) );
        if(j+2>m2)
         {
          m2+=5;
          x2=(float*)realloc((void*)x2,m2*sizeof(float));
          y2=(float*)realloc((void*)y2,m2*sizeof(float));
          z2=(float*)realloc((void*)z2,m2*sizeof(float));
          nrm2=(float*)realloc((void*)nrm2,3*m2*sizeof(float));
          ic2=(int*)realloc((void*)ic2,m2*sizeof(int));
         }
        x2[j]=(*x1)[i];
        y2[j]=(*y1)[i];
        z2[j]=(*z1)[i];
        nrm2[  3*j]=(*nrm1)[  3*i];
        nrm2[1+3*j]=(*nrm1)[1+3*i];
        nrm2[2+3*j]=(*nrm1)[2+3*i];
        ic2[j]=(*ic1)[i];
        j++;

        (*x1)[i]=(*x1)[i]+t*((*x1)[i+1]-(*x1)[i]);
        (*y1)[i]=(*y1)[i]+t*((*y1)[i+1]-(*y1)[i]);
        (*z1)[i]=(*z1)[i]+t*((*z1)[i+1]-(*z1)[i]);
        (*nrm1)[  3*i]=(*nrm1)[  3*i]+t*((*nrm1)[  3*(i+1)]-(*nrm1)[  3*i]);
        (*nrm1)[1+3*i]=(*nrm1)[1+3*i]+t*((*nrm1)[1+3*(i+1)]-(*nrm1)[1+3*i]);
        (*nrm1)[2+3*i]=(*nrm1)[2+3*i]+t*((*nrm1)[2+3*(i+1)]-(*nrm1)[2+3*i]);
        (*ic1)[i]=1;

        x2[j]=(*x1)[i];
        y2[j]=(*y1)[i];
        z2[j]=(*z1)[i];
        nrm2[  3*j]=(*nrm1)[  3*i];
        nrm2[1+3*j]=(*nrm1)[1+3*i];
        nrm2[2+3*j]=(*nrm1)[2+3*i];
        ic2[j]=(*ic1)[i];
        j++;

        in0=0;
       }
     }                             /* End of while(i+1<*n1)       */
    if(verbose){printf("done plane %d, Caught %d vertices\n",ip,j);fflush(stdout);} 

/*  Copy back from 2 into 1 */

    if(j+1>m1)
     {
      m1=j+3;
      *x1=(float*)realloc((void*)(*x1),m1*sizeof(float));
      *y1=(float*)realloc((void*)(*y1),m1*sizeof(float));
      *z1=(float*)realloc((void*)(*z1),m1*sizeof(float));
      *nrm1=(float*)realloc((void*)(*nrm1),3*m1*sizeof(float));
      *ic1=(int*)realloc((void*)(*ic1),m1*sizeof(int));
     }

    for(i=0;i<j;i++)
     {
      (*x1)[i]=x2[i];
      (*y1)[i]=y2[i];
      (*z1)[i]=z2[i];
      (*nrm1)[  3*i]=nrm2[  3*i];
      (*nrm1)[1+3*i]=nrm2[1+3*i];
      (*nrm1)[2+3*i]=nrm2[2+3*i];
      (*ic1)[i]=ic2[i];
     }
    (*x1)[j]=x2[0];
    (*y1)[j]=y2[0];
    (*z1)[j]=z2[0];
    (*nrm1)[  3*j]=nrm2[  3*0];
    (*nrm1)[1+3*j]=nrm2[1+3*0];
    (*nrm1)[2+3*j]=nrm2[2+3*0];
    (*ic1)[j]=ic2[0];
    (*n1)=j+1;
#endif
   }                               /* End of for(ip=0;ip<sh_nplns;ip++) */

  if(verbose){printf("Part of Polygon remains unclipped, %d vertices\n\n",j);fflush(stdout);} 

  free(x2);
  free(y2);
  free(z2);
  free(nrm2);
  free(ic2);

  return;
 }

#ifdef __cplusplus
}
#endif
