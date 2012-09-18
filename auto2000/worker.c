#include "auto_f2c.h"
#include "auto_c.h"
#include "auto_types.h"

/* The memory for these are taken care of in autobv_ and autoae_ */
extern struct {
  integer irtn;
  integer *nrtn;
} global_rotations;

#ifdef DEBUG
#include <stdio.h>
#include <stdarg.h>
void log_message(char *fmt, ...)
{
  va_list argp;
  FILE *log_file;
  log_file=fopen("/tmp/redrod_log","a");
  va_start(argp, fmt);
  vfprintf(log_file, fmt, argp);
  va_end(argp);
  fflush(log_file);
  fclose(log_file);
}
#endif

#ifdef MPI
double ***local_storage_aa=NULL;
double ***local_storage_bb=NULL;
double ***local_storage_cc=NULL;
double **local_storage_dd=NULL;
double **local_storage_fa=NULL;
double *local_storage_fc=NULL;

int mpi_worker() {
  MPI_Status stat;
  int my_rank;
  int message_type;
  integer funi_icni_params[5]={0,0,0,0,0};
  integer setup_common_blocks=0;
  MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
  
  while(1) {
    MPI_Recv(&message_type,1,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
    switch(message_type){
    case AUTO_MPI_KILL_MESSAGE:  /*The kill message*/
      MPI_Finalize();
      exit(0);
      break;
    case AUTO_MPI_INIT_MESSAGE:
      MPI_Bcast(funi_icni_params,5,MPI_LONG,0,MPI_COMM_WORLD);
      /* After we INIT we want setubv to initialize the common blocks
	 for the functions is autlib3.c and autlib5.c */
      setup_common_blocks = 1;
      break;
    case AUTO_MPI_SETUBV_MESSAGE:  /*The setubv message is not done yet*/
      mpi_setubv_worker(funi_icni_params,setup_common_blocks);
      /* We only want to initialize the blocks once for each
	 outer loop of main() */
      setup_common_blocks = 0;
      break;
    case AUTO_MPI_CONPAR_MESSAGE: /*The setubv message*/
      mpi_conpar_worker();
      break;
    default:
      fprintf(stderr,"Unknown message recieved: %d\n",message_type);
      break;
    }
  }
  return 0;
}
#endif

#ifdef MPI
int mpi_conpar_worker() {
  MPI_Status stat;
  int my_rank;
  integer params[5];
  
  integer nov, nra, nca;
  doublereal ***a;
  integer ncb;
  doublereal ***b;
  integer nrc;
  doublereal ***c, **d;
  integer *irf, *icf;
  integer na;

  /* find out which worker I am */
  MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
  /* input scalars */
  MPI_Recv(&na ,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);

  MPI_Bcast(params        ,5,MPI_LONG,0,MPI_COMM_WORLD);
  nov = params[0];
  nra = params[1];
  nca = params[2];
  ncb = params[3];
  nrc = params[4];

  /* input/output arrays */

  irf = (integer *)malloc(sizeof(integer)*na*nra);
  icf = (integer *)malloc(sizeof(integer)*na*nca);

  a=local_storage_aa;
  b=local_storage_bb;
  c=local_storage_cc;
  d=local_storage_dd;

  conpar_process(nov, na, nra, nca, a, ncb, b, nrc, c, d, irf, icf);

  MPI_Reduce(d[0],NULL,ncb*nrc,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Gatherv(a[0][0],nca*nra*na,MPI_DOUBLE,NULL,NULL,NULL,MPI_DOUBLE,
	      0,MPI_COMM_WORLD);
  MPI_Gatherv(b[0][0],ncb*nra*na,MPI_DOUBLE,NULL,NULL,NULL,MPI_DOUBLE,
	      0,MPI_COMM_WORLD);
  MPI_Gatherv(c[0][0],nca*nrc*na,MPI_DOUBLE,NULL,NULL,NULL,MPI_DOUBLE,
	      0,MPI_COMM_WORLD);
  MPI_Gatherv(irf,na*nra,MPI_LONG,NULL,NULL,NULL,MPI_LONG,0,MPI_COMM_WORLD);
  MPI_Gatherv(icf,na*nca,MPI_LONG,NULL,NULL,NULL,MPI_LONG,0,MPI_COMM_WORLD);
  
  /*free arrays*/
  free(irf);
  free(icf); 
  return 1;
}
#endif

#ifdef MPI
int mpi_setubv_worker(integer *funi_icni_params,integer setup_common_blocks) {
  MPI_Status stat;
  int my_rank,comm_size;
  int i,j;
  integer na, ndim, ncol, nint, ncb, nrc, nra, nca, ndxloc, loop_offset;
#ifdef MANIFOLD
  integer params[9];
  integer nalc;
#else
  integer params[8];
#endif
  iap_type *iap;
  rap_type *rap;
  doublereal *par;
  integer *icp;
  doublereal **ups, **uoldps, **wp, **wt, *wi, **udotps, **upoldp, *thu,
    *thl, *rldot, ***aa, ***bb, ***cc, **dd, **fa, *fc, *dtm;
  FUNI_TYPE((*funi));
  ICNI_TYPE((*icni));

  MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
  MPI_Comm_size(MPI_COMM_WORLD,&comm_size);
  /* input scalars */
  MPI_Recv(&na,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  MPI_Recv(&loop_offset,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);

#ifdef MANIFOLD
  MPI_Bcast(params          ,9,MPI_LONG,0,MPI_COMM_WORLD);
  nalc=params[8];
#else
  MPI_Bcast(params          ,8,MPI_LONG,0,MPI_COMM_WORLD);
#endif
  ndim=params[0];
  ncol=params[1];
  nint=params[2];
  ncb=params[3];
  nrc=params[4];
  nra=params[5];
  nca=params[6];
  ndxloc=params[7];

  iap=(iap_type *)malloc(sizeof(iap_type));
  rap=(rap_type *)malloc(sizeof(rap_type));
  par=(double *)malloc(sizeof(double)*NPARX2);
  icp=(integer *)malloc(sizeof(integer)*NPARX);
  ups=dmatrix(ndxloc, ndim*ncol);
  uoldps=dmatrix(ndxloc, ndim*ncol);
  wp=dmatrix(ncol + 1, ncol);
  wt=dmatrix(ncol + 1, ncol);
  wi=(double *)malloc(sizeof(double)*(ncol + 1));
  udotps=dmatrix(ndxloc, ndim*ncol);
  upoldp=dmatrix(ndxloc, ndim*ncol);
  thu=(double *)malloc(sizeof(double)*ndim*8);
  thl=(double *)malloc(sizeof(double)*ncb);
  rldot=(double *)malloc(sizeof(double)*NPARX);

  {
    int position=0;
    void *buffer;
    int bufsize;
    int size_int,size_double;
    int niap,nrap;
    /* Here we compute the number of elements in the iap and rap structures.
       Since each of the structures is homogeneous we just divide the total
       size by the size of the individual elements.*/
    niap = sizeof(iap_type)/sizeof(integer);
    nrap = sizeof(rap_type)/sizeof(doublereal);

    MPI_Pack_size(niap+NPARX,MPI_LONG,MPI_COMM_WORLD,&size_int);
    MPI_Pack_size(nrap+NPARX2+
		  ndxloc*ndim*ncol+
		  ndxloc*ndim*ncol+
		  (ncol + 1)*ncol+
		  (ncol + 1)*ncol+
		  (ncol + 1)+
		  ndxloc*ndim*ncol+
		  ndxloc*ndim*ncol+
		  ndim*8+
		  ncb+
		  NPARX,
		  MPI_DOUBLE,MPI_COMM_WORLD,&size_double);

    bufsize = size_int + size_double;
    buffer=malloc((unsigned)bufsize);

    MPI_Bcast(buffer  ,bufsize,MPI_PACKED,0,MPI_COMM_WORLD);

    MPI_Unpack(buffer,bufsize,&position,iap     ,niap,MPI_LONG,MPI_COMM_WORLD);
    MPI_Unpack(buffer,bufsize,&position,rap     ,nrap,MPI_DOUBLE,MPI_COMM_WORLD);
    /***********************************/
    MPI_Unpack(buffer,bufsize,&position,par    ,NPARX2,MPI_DOUBLE,MPI_COMM_WORLD);
    MPI_Unpack(buffer,bufsize,&position,icp    ,NPARX,MPI_LONG,MPI_COMM_WORLD);
    MPI_Unpack(buffer,bufsize,&position,ups[0]   ,ndxloc*ndim*ncol,MPI_DOUBLE,MPI_COMM_WORLD);
    MPI_Unpack(buffer,bufsize,&position,uoldps[0],ndxloc*ndim*ncol,MPI_DOUBLE,MPI_COMM_WORLD);
    MPI_Unpack(buffer,bufsize,&position,wp[0]  ,(ncol + 1)*ncol,MPI_DOUBLE,MPI_COMM_WORLD);
    MPI_Unpack(buffer,bufsize,&position,wt[0]  ,(ncol + 1)*ncol,MPI_DOUBLE,MPI_COMM_WORLD);
    MPI_Unpack(buffer,bufsize,&position,wi     ,(ncol + 1),MPI_DOUBLE,MPI_COMM_WORLD);
    MPI_Unpack(buffer,bufsize,&position,udotps[0],ndxloc*ndim*ncol,MPI_DOUBLE,MPI_COMM_WORLD);
    MPI_Unpack(buffer,bufsize,&position,upoldp[0],ndxloc*ndim*ncol,MPI_DOUBLE,MPI_COMM_WORLD);

    MPI_Unpack(buffer,bufsize,&position,thu    ,ndim*8,MPI_DOUBLE,MPI_COMM_WORLD);
    MPI_Unpack(buffer,bufsize,&position,thl    ,ncb,MPI_DOUBLE,MPI_COMM_WORLD);
    MPI_Unpack(buffer,bufsize,&position,rldot  ,NPARX,MPI_DOUBLE,MPI_COMM_WORLD);

  }

  if(setup_common_blocks) {
    /* At this point the iap structure is set up, so we can allocate some
       arrays which are used in autlib3.c and autlib5.c */
    allocate_global_memory(*(iap));
    
    /* This function sets up the blrtn global.  Some of the routines
       pointed to by icni and bcni use this structure */
    setrtn(iap->ndm,iap->ntst,ups,par);
  } 

  dtm=(double *)malloc(sizeof(double)*na);
  MPI_Scatterv(NULL,NULL,NULL,MPI_DOUBLE,
	      dtm,na,MPI_DOUBLE,
	      0,MPI_COMM_WORLD);

  /* output arrays */
  if(local_storage_aa==NULL)
    local_storage_aa=dmatrix_3d(na, nra, nca);
  aa = local_storage_aa;

  if(local_storage_bb==NULL)
    local_storage_bb=dmatrix_3d(na, nra, ncb);
  bb = local_storage_bb;

  if(local_storage_cc==NULL)
    local_storage_cc=dmatrix_3d(na, nrc, nca);
  cc = local_storage_cc;

  if(local_storage_dd==NULL)
    local_storage_dd=dmatrix(nrc, ncb);
  dd = local_storage_dd;

  if(local_storage_fa==NULL)
    local_storage_fa=dmatrix(na, ndim * ncol);
  fa = local_storage_fa;

  if(local_storage_fc==NULL && nint > 0)
    local_storage_fc=(doublereal *)malloc(sizeof(doublereal)*nint);
  fc = local_storage_fc;

  /* 
     A little explanation of what is going on here
     is in order I believe.  This array is
     created by a summation across all workers,
     hence it needs a summation in the master.
  */
  
  /* We sum into dd, which is a local variable initialized to
     0.0. We then sum our part with the masters part
     in the master. */

  /*zero pseudo-arclength part of array, rest is done in setubv() */
  for(i=nint;i<nrc;i++)
    for (j=0; j<ncb;j++)
      dd[i][j]=0.0;

  /* figure out what funi and icni are from
     the iap array.  This is originally done 
     in autlib1.c.  We do it here, since I
     don't know how to pass function pointers
     through MPI in a possibly heterogeneous 
     environment :-) */
  {
    function_list list;
    iap->ips  = funi_icni_params[0];
    iap->irs  = funi_icni_params[1];
    iap->isw  = funi_icni_params[2];
    iap->itp  = funi_icni_params[3];
    iap->nfpr = funi_icni_params[4];
    if(set_funi_and_icni(iap,&list)) {
      MPI_Finalize();
      exit(0);
    }
    funi = list.bvlist.funi;
    icni = list.bvlist.icni;
  }

  /*this call uses the loop_offset variable since up and uoldps
    and sent by the MPI version in their entirety, but
    aa, bb, cc, and fa have been shifted.  The loop_offset
    variable contains the original value of loop_start and shifts
    ups etc too*/
  setubv_make_aa_bb_cc_dd(ndim, na, ncol, nint,
#ifdef MANIFOLD
			  nalc,
#endif
			  ncb, nrc, nra, funi, icni, iap, rap, par, icp,
			  aa, bb, cc, dd, fa, fc, ups+loop_offset,
			  uoldps+loop_offset, udotps+loop_offset,
			  upoldp+loop_offset, dtm, thu, wi, wp, wt);

  MPI_Gatherv(fa[0],ndim*ncol*na,MPI_DOUBLE,NULL,NULL,NULL,MPI_DOUBLE,
	      0,MPI_COMM_WORLD);
  MPI_Reduce(fc,NULL,nint,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

  /*free input arrays*/
  free(iap);
  free(rap);
  free(par);
  free(icp);
  free_dmatrix(ups);
  free_dmatrix(uoldps);
  free(dtm);
  free(wp);
  free(wt);
  free_dmatrix(udotps);
  free_dmatrix(upoldp);
  free(wi);
  free(thu);
  free(thl);
  free(rldot);
  return 1;
}

int set_funi_and_icni(iap_type *iap,function_list *list) {
  set_function_pointers(*iap,list);

  if (list->type != AUTOBV) {
    printf("Illegal problem type in set_funi_and_icni\n");
    exit(1);
  }
  return 0;
}


#endif



