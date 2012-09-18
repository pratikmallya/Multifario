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

int main(int argc, char *argv[])
 {
  FILE *fin;
  FILE *fout;
  char file[4096];

  int i;
  int branch,nBranches,iBranch;
  double epsl,epsu,epss;
  double ds,dsmin,dsmax;
  double rl,rr,ul,ur;
  int ndim,ips,irs,ilp;
  int ntst,ncol,iad,isp;
  int isw,iplt,nbc,nint;
  int nmx,npr,mxbf,iid;
  int itmx,itnw,nwtn,jac,nuzr;
  int ipar;
  int iact;
  int nPts,mPts;
  int object;
  int point,type,label;
  double *par=NULL,*unorm=NULL;
  char c;

  if(argc<2){printf("Usage %s filename\n\n     where filename is the auto.7 file to be converted.",argv[0]);fflush(stdout);return 8;}

  strcpy(file,argv[1]);
  fin=fopen(file,"r");
  if(fin==(FILE*)NULL){printf("Error opening input file %s\n",file);fflush(stdout);return 8;}

  strcpy(file,argv[1]);
  strcat(file,".dx");
  fout=fopen(file,"w");
  if(fout==(FILE*)NULL){printf("Error opening output file %s\n",file);fflush(stdout);return 8;}

  fscanf(fin," %d %le %le %le %le\n",&branch,&rl,&rr,&ul,&ur);
  fscanf(fin," %d EPSL=%le EPSU =%le EPSS =%le\n",&branch,&epsl,&epsu,&epss);
  fscanf(fin," %d DS =%le DSMIN=%le DSMAX=%le\n",&branch,&ds,&dsmin,&dsmax);
  fscanf(fin," %d NDIM=%d IPS =%d IRS =%d ILP =%d\n",&branch,&ndim,&ips,&irs,&ilp);
  fscanf(fin," %d NTST=%d NCOL=%d IAD =%d ISP =%d\n",&branch,&ntst,&ncol,&iad,&isp);
  fscanf(fin," %d ISW =%d IPLT=%d NBC =%d NINT=%d\n",&branch,&isw,&iplt,&nbc,&nint);
  fscanf(fin," %d NMX=%d NPR =%d MXBF=%d IID =%d\n",&branch,&nmx,&npr,&mxbf,&iid);
  fscanf(fin," %d ITMX=%d ITNW =%d NWTN=%d JAC=%d NUZR=%d\n",&branch,&itmx,&itnw,&nwtn,&jac,&nuzr);
  fscanf(fin," %d User-specified parameter: %d\n",&branch,&ipar);
  fscanf(fin," %d Active continuation parameter: %d\n",&branch,&iact);
  while((c=fgetc(fin))!=(int)'\n');
  while((c=fgetc(fin))!=(int)'\n');

  nBranches=1;
  iBranch=1;
  object=0;
  nPts=0;
  mPts=0;
  while(1)
   {
    while(mPts<=nPts)
     {
      mPts+=100;
      par=(double*)realloc((void*)par,mPts*sizeof(double));
      unorm=(double*)realloc((void*)unorm,mPts*sizeof(double));
     }

    fscanf(fin," %d %d %d %d %le %le",&branch,&point,&type,&label,par+nPts,unorm+nPts);
    printf("%d %d %d %d %le %le",branch,point,type,label,par[nPts],unorm[nPts]);fflush(stdout);
    while((c=fgetc(fin))!=(int)'\n' && !feof(fin)){printf("%c",c);fflush(stdout);}
    if(feof(fin)){printf("\n");fflush(stdout);}
     else {printf("%c",c);fflush(stdout);}

    if(feof(fin))goto finish;
    nPts++;

    if(branch!=iBranch)
     {
      fprintf(fout,"object %d class array type float rank 1 shape 2 items %d data follows\n",object,nPts);
      for(i=0;i<nPts;i++)
       {
        fprintf(fout,"         %lf   %lf\n",par[i],unorm[i]);
       }
      fprintf(fout,"\n");
      object++;

      fprintf(fout,"object %d class gridconnections counts %d\n",object,nPts-1);
      fprintf(fout,"\n");
      object++;

      fprintf(fout,"object \"branch%d\" class field\n",iBranch);
      fprintf(fout,"component \"positions\" value %d\n",object-2);
      fprintf(fout,"component \"connections\" value %d\n",object-1);
      fprintf(fout,"\n");

      nPts=0;
      iBranch=branch;
      nBranches++;
      object++;
     }
   }

finish:
  if(nPts>0)
   {
    fprintf(fout,"object %d class array type float rank 1 shape 2 items %d data follows\n",object,nPts);
    for(i=0;i<nPts;i++)
      fprintf(fout,"         %lf   %lf\n",par[i],unorm[i]);
    fprintf(fout,"\n");
    object++;

    fprintf(fout,"object %d class gridconnections counts %d\n",object,nPts-1);
    fprintf(fout,"\n");
    object++;
 
    fprintf(fout,"object \"branch%d\" class field\n",iBranch);
    fprintf(fout,"component \"positions\" value %d\n",object-2);
    fprintf(fout,"component \"connections\" value %d\n",object-1);
    fprintf(fout,"\n");
    nBranches++;
   }

  fprintf(fout,"object \"default\" class group\n");
  for(i=1;i<nBranches;i++)
    fprintf(fout,"   member %d \"branch%d\"\n",i-1,i);
  fprintf(fout,"end\n");
  fclose(fin);
  fclose(fout);

  return 0;
 }

#ifdef __cplusplus
}
#endif
