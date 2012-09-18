#include "auto_f2c.h"
#include "auto_c.h"
#include "auto_types.h"

int 
setubv_make_aa_bb_cc_dd(integer ndim, integer na, integer ncol, integer nint, 
			integer nalc,
			integer ncb, integer nrc, integer nra, 
			FUNI_TYPE((*funi)), ICNI_TYPE((*icni)),
			iap_type *iap, rap_type *rap, doublereal *par, integer *icp,
			doublereal ***aa, doublereal ***bb, doublereal ***cc, doublereal **dd,
			doublereal **fa, doublereal *fc,
			doublereal **ups, doublereal **uoldps, doublereal **udotps, doublereal **upoldp,
			doublereal *dtm, doublereal *thu, doublereal *wi, doublereal **wp, doublereal **wt)
{
  /* System generated locals */
  integer dicd_dim1, dfdu_dim1, dfdp_dim1;
  
  /* Local variables */
  integer i, j, k, l, m;
  integer k1, l1;
  integer i1,j1;

  integer ib, ic;
  integer ic1;
  doublereal ddt;
  integer ndxloc;

  doublereal *dicd, *ficd, *dfdp, *dfdu, *uold;
  doublereal *f;
  doublereal *u, *wploc;
  doublereal *uic, *uio, *prm, *uid, *uip;

#ifdef USAGE
  struct rusage *setubv_make_aa_bb_cc_usage,*fa_usage;
  usage_start(&setubv_make_aa_bb_cc_usage);
#endif

  if (nint > 0) {
      dicd = (doublereal *)malloc(sizeof(doublereal)*nint*(ndim + NPARX));
      ficd = (doublereal *)malloc(sizeof(doublereal)*nint);
  }
  else
      ficd = dicd = NULL;

  ndxloc=iap->ntst+1;
  
  dfdp = (doublereal *)malloc(sizeof(doublereal)*ndim*NPARX);
  dfdu = (doublereal *)malloc(sizeof(doublereal)*ndim*ndim);
  uold = (doublereal *)malloc(sizeof(doublereal)*ndim);
  f    = (doublereal *)malloc(sizeof(doublereal)*ndim);
  u    = (doublereal *)malloc(sizeof(doublereal)*ndim);
  wploc= (doublereal *)malloc(sizeof(doublereal)*(ncol+1));
  uic  = (doublereal *)malloc(sizeof(doublereal)*ndim);
  uio  = (doublereal *)malloc(sizeof(doublereal)*ndim);
  prm  = (doublereal *)malloc(sizeof(doublereal)*NPARX);
  uid  = (doublereal *)malloc(sizeof(doublereal)*ndim);
  uip  = (doublereal *)malloc(sizeof(doublereal)*ndim);

  dicd_dim1 = nint;
  dfdu_dim1 = ndim;
  dfdp_dim1 = ndim;

  /* Generate AA, BB and DD: */
  
  /* Initialize to zero. */
  for (i = 0; i < nint; ++i) {
    for (k = 0; k < ncb; ++k) {
      dd[i][k] = 0.;
    }
    fc[i] = 0;
  }

  /*      Partition the mesh intervals */
  /*j will be replaced with 0 and na*/
  for (j = 0; j < na; ++j) {
    doublereal *up = ups[j];
    doublereal *up1 = ups[j + 1];
    doublereal *uoldp = uoldps[j];
    doublereal *uoldp1 = uoldps[j + 1];
      
    ddt = 1. / dtm[j];
    for (ic = 0; ic < ncol; ++ic) {
      for (k = 0; k < ndim; ++k) {
	u[k] = wt[ncol][ic] * up1[k];
	uold[k] = wt[ncol][ic] * uoldp1[k];
	for (l = 0; l < ncol; ++l) {
	  l1 = l * ndim + k;
	  u[k] += wt[l][ic] * up[l1];
	  uold[k] += wt[l][ic] * uoldp[l1];
	}
      }

      for (i = 0; i < NPARX; ++i) {
	prm[i] = par[i];
      }
      /*  
	  Ok this is a little wierd, so hold tight.  This function
	  is actually a pointer to a wrapper function, which eventually
	  calls the user defined func_.  Which wrapper is used
	  depends on what kind of problem it is.
      */
      (*(funi))(iap, rap, ndim, u, uold, icp, prm, 2, f, dfdu, dfdp);
      /* transpose dfdu for optimal access */
      {
      integer ii, jj;
      doublereal tmp;
      for (ii = 0; ii < ndim; ++ii) {
        for (jj = 0; jj < ii; ++jj) {
          tmp = dfdu[ii + jj * ndim];
          dfdu[ii + jj * ndim] = dfdu[jj + ii * ndim];
          dfdu[jj + ii * ndim] = tmp;
        }
      }
      ic1 = ic * ndim;
      for (ib = 0; ib < ncol + 1; ++ib) {
	wploc[ib] = ddt * wp[ib][ic];
      }
      for (i = 0; i < ndim; ++i) {
	double *aa_offset = aa[j][ic1 + i];
	double *dfdu_offset = &ARRAY2D(dfdu, 0, i);
	for (ib = 0; ib < ncol + 1; ++ib) {
	  double wt_tmp = -wt[ib][ic];
	  for (k = 0; k < ndim; ++k) {
	    aa_offset[k] = wt_tmp * dfdu_offset[k];
	  }
	  aa_offset[i] += wploc[ib];
	  aa_offset += ndim;
	}
	for (k = 0; k < ncb; ++k) {
	  bb[j][ic1 + i][k] = -ARRAY2D(dfdp, i, icp[k]);
	}
	fa[j][ic1 + i] = f[i] - wploc[ncol] * up1[i];
	for (k = 0; k < ncol; ++k) {
	  k1 = k * ndim + i;
	  fa[j][ic1 + i] -= wploc[k] * up[k1];
	}
      }
      }
    }
  }

  /* generate CC and DD; the boundary conditions are not
     done parallelly */
  
  /*     Integral constraints : */
  if (nint > 0)
   {
    for (j = 0; j < na; ++j)
     {
      int jp1 = j + 1;
      for (k = 0; k <= ncol; ++k)
       {
	for (i = 0; i < ndim; ++i)
         {
	  i1 = k * ndim + i;
	  j1 = j;
	  if (k == ncol)i1 = i;
	  if (k == ncol)j1 = jp1;
	  uic[i] = ups[j1][i1];
	  uio[i] = uoldps[j1][i1];
	  uid[i] = udotps[j1][i1];
	  uip[i] = upoldp[j1][i1];
	}
	
	(*(icni))(iap, rap, ndim, par, icp, nint, uic, uio, uid, uip, ficd, 2, dicd);
	
	for (m = 0; m < nint; ++m)
         {
	  k1 = k * ndim;
	  for (i = 0; i < ndim; ++i)
           {
             cc[j][m][k1+i] = dtm[j] * wi[k ] * ARRAY2D(dicd, m, i);
           }
	  fc[m] -= dtm[j] * wi[k] * ficd[m];
	  for (i = 0; i < ncb; ++i) dd[m][i] += dtm[j] * wi[k] * ARRAY2D(dicd, m, ndim + icp[i]);
	}
      }
    }
  }

  /*     Pseudo-arclength equation : */

  for (j = 0; j < na; ++j)
   {
    for (m = 0; m < nalc; ++m)
     {
      for (i = 0; i < ndim; ++i)
       {
        for (k = 0; k < ncol; ++k)
         {
	  k1 = k * ndim + i;
          cc[j][nint+m][k1] = dtm[j] * thu[i] * wi[k] * udotps[j + m * ndxloc][k1];
         }
        cc[j][nint+m][nra + i] = dtm[j] * thu[i] * wi[ncol] * udotps[j + 1 + m*ndxloc][i];
       }
     }
   }

  free(dicd );
  free(ficd );
  free(dfdp );
  free(dfdu );
  free(uold );
  free(f    );
  free(u    );
  free(wploc);
  free(uic  );
  free(uio  );
  free(prm  );
  free(uid  );
  free(uip  );

#ifdef USAGE
  usage_end(setubv_make_aa_bb_cc_usage,"in setubv worker");
#endif
  return 0;
}

#ifdef PTHREADS
static void *setubv_threads_make_aa_bb_cc_dd(void * arg)
{
  setubv_parallel_arglist *p =  (setubv_parallel_arglist *)arg;
  /* In the shared memory case we sum into a local
     variable our contribution, and then sum
     into shared memory at the end */
  p->dd = dmatrix(p->nint, p->ncb);
  p->fc = malloc(p->nint * sizeof(doublereal));
  setubv_make_aa_bb_cc_dd(p->ndim, p->na, p->ncol, p->nint, 
	p->nalc,
	p->ncb, p->nrc, p->nra, p->funi, p->icni, p->iap,
	p->rap, p->par, p->icp, p->aa, p->bb, p->cc, p->dd,
	p->fa, p->fc, p->ups, p->uoldps, p->udotps, p->upoldp,
	p->dtm, p->thu, p->wi, p->wp, p->wt);
  return NULL;
}

/* Fill in a setubv_parallel_arglist for the individual variables */
static int setubv_threads_wrapper(integer ndim, integer na, integer ncol, integer nint,
				  integer nalc,
				  integer ncb, integer nrc, integer nra,
				  FUNI_TYPE((*funi)), ICNI_TYPE((*icni)), iap_type *iap, rap_type *rap, doublereal *par, 
				  integer *icp, doublereal ***aa, doublereal ***bb, 
				  doublereal ***cc, doublereal **dd, doublereal **fa, doublereal *fc, doublereal **ups, 
				  doublereal **uoldps, doublereal **udotps, 
				  doublereal **upoldp, doublereal *dtm, 
				  doublereal *thu, doublereal *wi, doublereal **wp, doublereal **wt)
{
  setubv_parallel_arglist *send_data;
  int i,j,k;
  pthread_t *th;
  void * retval;
  pthread_attr_t attr;
  int retcode;

#ifdef USAGE
  struct timeval *pthreads_create,*pthreads_join,*pthreads_all;
  time_start(&pthreads_create);
  time_start(&pthreads_all);
#endif

  th = (pthread_t *)malloc(sizeof(pthread_t)*global_num_procs);
  send_data = (setubv_parallel_arglist *)malloc(sizeof(setubv_parallel_arglist)*global_num_procs);
  pthread_attr_init(&attr);
  pthread_attr_setscope(&attr,PTHREAD_SCOPE_SYSTEM);

  for(i=0;i<global_num_procs;i++) {
    integer loop_start;
    setubv_parallel_arglist *p = &send_data[i];
    loop_start = (i*na)/global_num_procs;
    p->ndim   = ndim;
    p->na     = ((i+1)*na)/global_num_procs - loop_start;
    p->ncol   = ncol;
    p->nint   = nint;
    p->nalc   = nalc;
    p->ncb    = ncb;
    p->nrc    = nrc;
    p->nra    = nra;
    p->funi   = funi;
    p->icni   = icni;
    p->iap    = iap;
    p->rap    = rap;
    p->par    = par;
    p->icp    = icp;
    p->aa     = aa + loop_start;
    p->bb     = bb + loop_start;
    p->cc     = cc + loop_start;
    p->dd     = NULL;
    p->fa     = fa + loop_start;
    p->fc     = NULL;
    p->ups    = ups + loop_start;
    p->uoldps = uoldps + loop_start;
    p->udotps = udotps + loop_start;
    p->upoldp = upoldp + loop_start;
    p->dtm    = dtm + loop_start;
    p->wp     = wp;
    p->wt     = wt;
    p->wi     = wi;
    p->thu    = thu;
    retcode = pthread_create(&th[i], &attr, setubv_threads_make_aa_bb_cc_dd, (void *)p);
    if (retcode != 0) fprintf(stderr, "create %d failed %d\n", i, retcode);
  }
#ifdef USAGE
  time_end(pthreads_create,"setubv pthreads create",fp9);
  time_start(&pthreads_join);
#endif
  for(i=0;i<global_num_procs;i++) {
    setubv_parallel_arglist *p = &send_data[i];
    retcode = pthread_join(th[i], &retval);
    if (retcode != 0) fprintf(stderr, "join %d failed %d\n", i, retcode);
    /* This is where we sum into the global copy of the dd
       array, in the shared memory case */
    for(j=0;j<nint;j++) {
      fc[j] += p->fc[j];
      for (k=0; k<ncb;k++)
        dd[j][k] += p->dd[j][k];
    }
    free(p->fc);
    free_dmatrix(p->dd);
  }
  free(send_data);
  free(th);
#ifdef USAGE
  time_end(pthreads_join,"setubv pthreads join",fp9);
  time_end(pthreads_all,"setubv pthreads all",fp9);
#endif

  return 0;
}
#endif

#ifdef MPI
int 
setubv_mpi_wrapper(integer ndim, integer na, integer ncol, integer nint,
       integer nalc,
       integer ncb, integer nrc, integer nra, integer nca,
       integer ndxloc, iap_type *iap, rap_type *rap, doublereal *par, integer *icp,
       doublereal **fa, doublereal *fc,
       doublereal *rldot, doublereal **ups, doublereal **uoldps, doublereal **udotps, doublereal **upoldp,
       doublereal *dtm, doublereal *thl, doublereal *thu,
       doublereal *wi, doublereal **wp, doublereal **wt)
{
  integer loop_start,loop_end,local_na;
  int i,comm_size;
  int *fa_counts,*fa_displacements;
  int *dtm_counts,*dtm_displacements;

  MPI_Comm_size(MPI_COMM_WORLD,&comm_size);
  fa_counts=(int *)malloc(sizeof(int)*comm_size);
  fa_displacements=(int *)malloc(sizeof(int)*comm_size);
  dtm_counts=(int *)malloc(sizeof(int)*comm_size);
  dtm_displacements=(int *)malloc(sizeof(int)*comm_size);
  fa_counts[0] = 0;
  fa_displacements[0] = 0;
  dtm_counts[0] = 0;
  dtm_displacements[0] = 0;
  
  for(i=1;i<comm_size;i++){
    
    /*Send message to get worker into setubv mode*/
    {
      int message=AUTO_MPI_SETUBV_MESSAGE;
      MPI_Send(&message,1,MPI_INT,i,0,MPI_COMM_WORLD);
    }
    loop_start = ((i-1)*na)/(comm_size - 1);
    loop_end = (i*na)/(comm_size - 1);
    fa_counts[i] = ndim*ncol*(loop_end-loop_start);
    fa_displacements[i] = ndim*ncol*loop_start;
    dtm_counts[i] = (loop_end-loop_start);
    dtm_displacements[i] = (loop_start);

    local_na = loop_end-loop_start;
    MPI_Send(&local_na ,1,MPI_LONG,i,0,MPI_COMM_WORLD);
    MPI_Send(&loop_start    ,1,MPI_LONG,i,0,MPI_COMM_WORLD);
  }

  {
    integer params[11];
    params[0]=ndim;
    params[1]=ncol;
    params[2]=nint;
    params[3]=ncb;
    params[4]=nrc;
    params[5]=nra;
    params[6]=nca;
    params[7]=ndxloc;
    params[8]=nalc;
    MPI_Bcast(params     ,9,MPI_LONG,0,MPI_COMM_WORLD);
  }    

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

    MPI_Pack(iap    ,niap,MPI_LONG,buffer,bufsize,&position,MPI_COMM_WORLD);
    MPI_Pack(rap    ,nrap,MPI_DOUBLE,buffer,bufsize,&position,MPI_COMM_WORLD);
    /**********************************************/
    MPI_Pack(par    ,NPARX2,MPI_DOUBLE,buffer,bufsize,&position,MPI_COMM_WORLD);
    MPI_Pack(icp    ,NPARX,MPI_LONG,buffer,bufsize,&position,MPI_COMM_WORLD);
    MPI_Pack(ups[0] ,ndxloc*ndim*ncol,MPI_DOUBLE,buffer,bufsize,&position,MPI_COMM_WORLD);
    MPI_Pack(uoldps[0],ndxloc*ndim*ncol,MPI_DOUBLE,buffer,bufsize,&position,MPI_COMM_WORLD);
    MPI_Pack(wp[0]  ,(ncol + 1)*ncol,MPI_DOUBLE,buffer,bufsize,&position,MPI_COMM_WORLD);
    MPI_Pack(wt[0]  ,(ncol + 1)*ncol,MPI_DOUBLE,buffer,bufsize,&position,MPI_COMM_WORLD);
    MPI_Pack(wi     ,(ncol + 1),MPI_DOUBLE,buffer,bufsize,&position,MPI_COMM_WORLD);
    MPI_Pack(udotps[0],ndxloc*ndim*ncol,MPI_DOUBLE,buffer,bufsize,&position,MPI_COMM_WORLD);
    MPI_Pack(upoldp[0],ndxloc*ndim*ncol,MPI_DOUBLE,buffer,bufsize,&position,MPI_COMM_WORLD);

    MPI_Pack(thu    ,ndim*8,MPI_DOUBLE,buffer,bufsize,&position,MPI_COMM_WORLD);
    MPI_Pack(thl    ,ncb,MPI_DOUBLE,buffer,bufsize,&position,MPI_COMM_WORLD);
    MPI_Pack(rldot  ,NPARX,MPI_DOUBLE,buffer,bufsize,&position,MPI_COMM_WORLD);
    
    MPI_Bcast(buffer     ,position,MPI_PACKED,0,MPI_COMM_WORLD);
  }

  MPI_Scatterv(dtm        ,dtm_counts,dtm_displacements,MPI_DOUBLE,
	       NULL,0,MPI_DOUBLE,
	       0,MPI_COMM_WORLD);

  /* Worker runs here */

  MPI_Gatherv(NULL,0,MPI_DOUBLE,
	      fa[0],fa_counts,fa_displacements,MPI_DOUBLE,
	      0,MPI_COMM_WORLD);
  {
    /*I create a temporary send buffer for the MPI_Reduce
      command.  This is because there isn't an
      asymmetric version (like MPI_Scatterv).*/
    double *fctemp = malloc(nint*sizeof(doublereal));
    for(i=0;i<nint;i++)
      fctemp[i]=fc[i];
    MPI_Reduce(fctemp,fc,nint,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    free(fctemp);
  }
  
  return 0;
}
#endif

static void setubv_pseudo_arclength(integer ndim, integer nalc, integer ncb, integer nrc, iap_type *iap, doublereal *rds, doublereal **dd, doublereal *fc, doublereal *rlcur, doublereal *rlold, doublereal *rldot, doublereal **udotps, doublereal **dups, doublereal *dtm, doublereal *thl, doublereal *thu);

int setubv(integer ndim, integer ips, integer na, integer ncol, integer nint, integer nalc, integer ncb, integer nrc, integer nra, integer nca, FUNI_TYPE((*funi)), ICNI_TYPE((*icni)), integer ndxloc, iap_type *iap, rap_type *rap, doublereal *par, integer *icp, doublereal *rds, doublereal ***aa, doublereal ***bb, doublereal ***cc, doublereal **dd, doublereal **fa, doublereal *fc, doublereal *rlcur, doublereal *rlold, doublereal *rldot, doublereal **ups, doublereal **uoldps, doublereal **udotps, doublereal **upoldp, doublereal **dups, doublereal *dtm, doublereal *thl, doublereal *thu)
 {
  doublereal *wi, **wp, **wt;
  int i,j;
  
#ifdef USAGE
  struct rusage *initialization_usage,*fc_usage,*parallel_overhead_usage;
  usage_start(&initialization_usage);
#endif
  wi   = (doublereal *)malloc(sizeof(doublereal)*(ncol+1) );
  wp   = dmatrix(ncol+1, ncol);
  wt   = dmatrix(ncol+1, ncol);

  wint(ncol + 1, wi);
  genwts(ncol, ncol + 1, wt, wp);
  
#ifdef USAGE
  usage_end(initialization_usage,"setubv initialization");
#endif

#ifdef USAGE
  usage_start(&fc_usage);
#endif

/* dups should be precomputed! (it doesn't seem to be) */

  for(i=0;i<iap->ntst;i++)
   for(j=0;j<iap->ndim*iap->ncol;j++)
    dups[i][j]=ups[i][j]-uoldps[i][j];
  for(j=0;j<iap->ncol;j++)
   dups[iap->ntst][j]=ups[iap->ntst][j]-uoldps[iap->ntst][j];

  setubv_pseudo_arclength(ndim, nalc, ncb, nrc, iap, rds, dd, fc, rlcur, rlold, rldot, udotps, dups, dtm, thl, thu);

#ifdef USAGE
  usage_end(fc_usage,"setubv make fc");
#endif

  switch(global_setubv_type) {

#ifdef PTHREADS
  case SETUBV_PTHREADS:
    {
      setubv_threads_wrapper(ndim, na, ncol, nint,
			     nalc,
			     ncb, nrc, nra, funi, icni, iap, rap, par, icp,
			     aa, bb, cc, dd, fa,
			     fc, ups, uoldps, udotps, upoldp, dtm,
			     thu, wi, wp, wt);
      break;
    }
#endif

#ifdef MPI
    case SETUBV_MPI:
      if(global_verbose_flag)
	printf("Setubv MPI start\n");
      setubv_mpi_wrapper(ndim, na, ncol, nint, nalc, ncb, nrc, nra, nca, ndxloc, iap, rap, par, icp, fa, fc, rldot, ups, uoldps, udotps, upoldp, dtm, thl, thu, wi, wp, wt);
      if(global_verbose_flag)
	printf("Setubv MPI end\n");
      break;
#endif

    default:
      setubv_make_aa_bb_cc_dd(ndim, na, ncol, nint, nalc, ncb, nrc, nra, funi, icni, iap, rap, par, icp, aa, bb, cc, dd, fa, fc, ups, uoldps, udotps, upoldp, dtm, thu, wi, wp, wt);
      break;
    }

  free(wi   );
  free_dmatrix(wp);
  free_dmatrix(wt);
  return 0;
}

void setubv_make_boundary(integer ndim, integer na, integer nbc,
       integer ncb, integer nra, BCNI_TYPE((*bcni)),
       iap_type *iap, rap_type *rap, doublereal *par,
       integer *icp, doublereal ***ccbc, doublereal **dd, doublereal *fc,
       doublereal *rlcur, doublereal *rlold,
       doublereal **ups, doublereal **uoldps, doublereal **dups)
{
  integer i,j,k;
  integer dbc_dim1 = nbc;

  doublereal *dbc  = (doublereal *)malloc(sizeof(doublereal)*(nbc)*(2*ndim + NPARX));
  doublereal *fbc  = (doublereal *)malloc(sizeof(doublereal)*(nbc));
  doublereal *ubc0 = (doublereal *)malloc(sizeof(doublereal)*ndim);
  doublereal *ubc1 = (doublereal *)malloc(sizeof(doublereal)*ndim);
  
  /* Set constants. */
  for (i = 0; i < ncb; ++i) {
    par[icp[i]] = rlcur[i];
  }
  
  /*     ** Time evolution computations (parabolic systems) */
  if (iap->ips == 14 || iap->ips == 16) {
    rap->tivp = rlold[0];
  } 

  /* Boundary condition part of FC */
  if (nbc > 0) {
    for (i = 0; i < ndim; ++i) {
      ubc0[i] = ups[0][i];
      ubc1[i] = ups[na][i];
    }
    
    (*(bcni))(iap, rap, ndim, par, icp, nbc, ubc0, ubc1, fbc, 2, dbc);
    for (i = 0; i < nbc; ++i) {
      fc[i] = -fbc[i];
      for (k = 0; k < ndim; ++k) {
	/*NOTE!!
	  This needs to split up.  Only the first processor does the first part
	  and only the last processors does the last part.
          (I leave this non-parallel for now since
           a) it doesn't play well with HomCont
           b) there is almost nothing to be gained -- Bart)
        */
	ccbc[0][i][k] = ARRAY2D(dbc, i, k);
	ccbc[1][i][k] = ARRAY2D(dbc ,i , ndim + k);
      }

      for (k = 0; k < ncb; ++k) {
	dd[i][k] = ARRAY2D(dbc, i, (ndim *2) + icp[k]);
      }
    }
    /*       Save difference : */
    for (j = 0; j < na + 1; ++j) {
      for (i = 0; i < nra; ++i) {
	dups[j][i] = ups[j][i] - uoldps[j][i];
      }
    }
  }

  free(dbc);
  free(fbc);
  free(ubc0);
  free(ubc1);
}

static void setubv_pseudo_arclength(integer ndim, integer nalc, integer ncb, integer nrc, iap_type *iap, doublereal *rds, doublereal **dd, doublereal *fc, doublereal *rlcur, doublereal *rlold, doublereal *rldot, doublereal **udotps, doublereal **dups, doublereal *dtm, doublereal *thl, doublereal *thu)
 {
  doublereal rlsum;
  integer i,j,k;
  integer m;
  int nint;
  int ntst;
  int ndxloc;
  int ncol;
  int nfpr;

  nint=iap->nint;
  ntst=iap->ntst;
  ncol=iap->ncol;
  nfpr=iap->nfpr;
  ndxloc=ntst+1;

  if(0){printf("setubv:\n");fflush(stdout);}

  for (m = 0; m < nalc; ++m)

  for (m = 0; m < nalc; ++m)
   {
    for (i = 0; i < ncb; ++i)
     {
      dd[nint+m][i] = thl[i] * rldot[i+m*NPARX];
     }

    rlsum = 0.;
    for (i = 0; i < ncb; ++i)rlsum += thl[i] * (rlcur[i] - rlold[i]) * rldot[i+m*NPARX];

    fc[nint+m] = rds[m] - rinpr(iap, ndim, udotps+m*ndxloc, dups, dtm, thu) - rlsum;
   }

  return;
}
