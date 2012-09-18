#include "auto_f2c.h"
#include "auto_c.h"
#include "auto_types.h"

/*This is the process function.  It can be called on a SMP using shared memory,
  or wrapped inside another routine for message passing*/

int conpar_process(integer nov, integer na, integer nra, integer nca, doublereal ***a, integer ncb, doublereal ***b, integer nrc, doublereal ***c, doublereal **d, integer *irf, integer *icf)
{
  /* Local variables */
  integer ipiv, jpiv, itmp;
  doublereal tpiv;
  integer i, j, l, k1, m2, ic, ir;
  doublereal rm;
  doublereal piv;
  int *amaxima;

#ifdef USAGE
  struct rusage *conpar_process_usage;
  usage_start(&conpar_process_usage);
#endif

  amaxima = malloc(nra*sizeof(*amaxima));

  /* Note that the summation of the adjacent overlapped part of C */
  /* is delayed until REDUCE, in order to merge it with other communications.*/
  /* NA is the local NTST. */
  
  /* Condensation of parameters (Elimination of local variables). */
  m2 = nca - nov;

  for (i = 0;i < na; i++) {
    doublereal **a_i = a[i];
    doublereal **b_i = b[i];
    integer *icf_i = &icf[i * nca];
    integer *irf_i = &irf[i * nra];
    for (j = 0; j < nra; ++j) {
      integer one = 1;
      irf_i[j] = j;
      amaxima[j] = nov + auto_idamax(m2 - nov, &a_i[j][nov], &one) - 1;
    }
    for (j = 0; j < nca; ++j) {
      icf_i[j] = j;
    }
    for (ic = nov; ic < m2; ++ic) {
      int irp = ic - nov;
      /*	     **Search for pivot (Complete pivoting) */
      piv = 0.0;
      ipiv = irp;
      jpiv = ic;
      for (k1 = irp; k1 < nra; ++k1) {
	int row = irf_i[k1];
	tpiv = fabs(a_i[row][icf_i[amaxima[row]]]);
	if (piv < tpiv) {
	  piv = tpiv;
	  ipiv = k1;
	  jpiv = amaxima[row];
	}
      }
      /*	     **Move indices */
      itmp = icf_i[ic];
      icf_i[ic] = icf_i[jpiv];
      icf_i[jpiv] = itmp;
      itmp = irf_i[irp];
      irf_i[irp] = irf_i[ipiv];
      irf_i[ipiv] = itmp;
      {
	int icf_ic_i = icf_i[ic];
	int irf_irp_i = irf_i[irp];
	doublereal *a_offset2 = a_i[irf_irp_i];
	doublereal *b_offset2 = b_i[irf_irp_i];
	/*	     **End of pivoting; elimination starts here */
	int icp1 = ic + 1;
	piv = a_offset2[icf_ic_i];
	for (ir = irp + 1; ir < nra; ++ir) {
	  int irf_ir_i = irf_i[ir];
          doublereal *a_offset1 = a_i[irf_ir_i];
	  doublereal *b_offset1 = b_i[irf_ir_i];
	  rm = a_offset1[icf_ic_i] / piv;
	  a_offset1[icf_ic_i] = rm;
	  if (rm != (double)0.) {
	    for (l = 0; l < nov; ++l) {
	      a_offset1[l] -= rm * a_offset2[l];
	    }
	    for (l = icp1; l < nca; ++l) {
	      int icf_l_i = icf_i[l];
	      a_offset1[icf_l_i] -= rm * a_offset2[icf_l_i];
	    }
	    for (l = 0; l < ncb; ++l) {
	      b_offset1[l] -= rm * b_offset2[l];
	    }
	  }
	  if (rm != (double)0. || amaxima[irf_ir_i] == jpiv) {
	    doublereal ppiv = 0;
	    integer jppiv = icp1;
	    /* recalculate absolute maximum for current row */
	    for (l = icp1; l < m2; ++l) {
	      tpiv = fabs(a_offset1[icf_i[l]]);
	      if (ppiv < tpiv) {
		ppiv = tpiv;
		jppiv = l;
	      }
	    }
	    amaxima[irf_ir_i] = jppiv;
	  } else if (amaxima[irf_ir_i] == ic) {
	    amaxima[irf_ir_i] = jpiv;
	  }
	}
	for (ir = 0; ir < nrc; ++ir) {
	  doublereal *c_offset1 = c[i][ir];
	  rm = c_offset1[icf_ic_i] / piv;
	  c_offset1[icf_ic_i]=rm;
	  if (rm != (double)0.) {
	    doublereal *d_offset1 = d[ir];
	    for (l = 0; l < nov; ++l) {
	      c_offset1[l] -= rm * a_offset2[l];
	    }
	    for (l = icp1; l < nca; ++l) {
	      int icf_l_i = icf_i[l];
	      c_offset1[icf_l_i] -= rm * a_offset2[icf_l_i];
	    }
	    for (l = 0; l < ncb; ++l) {
              d_offset1[l] -= rm * b_offset2[l];
	    }
	  }
	}
      }
    }
  }
  free(amaxima);
#ifdef USAGE
  usage_end(conpar_process_usage,"in conpar worker");
#endif
  return 0;
}

/*This is the process function.  It is meant to be called either
  on a SMP using shared memory, or wrapped inside another
  routine for message passing*/

#ifdef PTHREADS

static void *conpar_threads_process(void * arg)
{
  conpar_parallel_arglist *p = arg;
  integer i,j;
  /* A little explanation of what is going on here
     is in order I believe.  This array is
     created by a summation across all workers. */
  p->d = dmatrix(p->nrc, p->ncb);
  /* In the shared memory case we sum into a local
     variable our contribution, and then sum
     into shared memory at the end */
  for(i=0;i<p->nrc;i++)
    for (j=0; j<p->ncb;j++)
      p->d[i][j]=0.0;
  conpar_process(p->nov, p->na, p->nra, p->nca, p->a, p->ncb, p->b,
		 p->nrc, p->c, p->d, p->irf, p->icf);
  return NULL;
}

static int conpar_threads_wrapper(integer nov, integer na, integer nra, integer nca, doublereal ***a, integer ncb, doublereal ***b, integer nrc, doublereal ***c, doublereal **d, integer *irf, integer *icf)

{
  conpar_parallel_arglist *data;
  int i,j,k;
  pthread_t *th;
  void * retval;
  pthread_attr_t attr;
  int retcode;
#ifdef USAGE
  struct timeval *pthreads_wait;
  time_start(&pthreads_wait);
#endif
  
  data = (conpar_parallel_arglist *)malloc(sizeof(conpar_parallel_arglist)*global_num_procs);
  th = (pthread_t *)malloc(sizeof(pthread_t)*global_num_procs);

  for(i=0;i<global_num_procs;i++) {
    integer loop_start, loop_end;
    
    /*start and end of the computed loop*/
    loop_start = (i*na)/global_num_procs;
    loop_end = ((i+1)*na)/global_num_procs;
    
    /*3D Arrays*/ 
    data[i].a = a + loop_start;
    data[i].b = b + loop_start;
    data[i].c = c + loop_start;
    
    /*2D Arrays*/
    data[i].irf = irf + loop_start*nra;
    data[i].icf = icf + loop_start*nca;
    
    /*Scalars*/
    data[i].nrc = nrc;
    data[i].ncb = ncb;
    data[i].nov = nov;
    data[i].nra = nra;
    data[i].nca = nca;
    
    data[i].na = loop_end - loop_start;
    
  }
  pthread_attr_init(&attr);
  pthread_attr_setscope(&attr,PTHREAD_SCOPE_SYSTEM);
  for(i=0;i<global_num_procs;i++) {
    retcode = pthread_create(&th[i], &attr, conpar_threads_process, (void *) &data[i]);
    if (retcode != 0) fprintf(stderr, "create %d failed %d\n", i, retcode);
  }
  
  for(i=0;i<global_num_procs;i++) {
    retcode = pthread_join(th[i], &retval);
    if (retcode != 0) fprintf(stderr, "join %d failed %d\n", i, retcode);
    /* This is were we sum into the global copy of the d array */  
    for(j=0;j<nrc;j++)
      for (k=0; k<ncb;k++)
	d[j][k] += data[i].d[j][k];
    free_dmatrix(data[i].d);
  }  

  free(data);
  free(th);
#ifdef USAGE
  time_end(pthreads_wait,"conpar pthreads wrapper",fp9);
#endif
  return 0;
}
#endif

#ifdef MPI
int 
conpar_mpi_wrapper(integer nov, integer na, integer nra, 
		   integer nca, doublereal ***a, integer ncb, 
		   doublereal ***b, integer nrc, 
		   doublereal ***c, doublereal **d, integer *irf, integer *icf)

{
    integer loop_start,loop_end;
    integer loop_end_tmp;
    int i,j,comm_size;
    int *a_counts,*a_displacements;
    int *b_counts,*b_displacements;
    int *c_counts,*c_displacements;
    int *irf_counts,*irf_displacements;
    int *icf_counts,*icf_displacements;


    MPI_Comm_size(MPI_COMM_WORLD,&comm_size);
    a_counts=(int *)malloc(sizeof(int)*comm_size);
    a_displacements=(int *)malloc(sizeof(int)*comm_size);
    b_counts=(int *)malloc(sizeof(int)*comm_size);
    b_displacements=(int *)malloc(sizeof(int)*comm_size);
    c_counts=(int *)malloc(sizeof(int)*comm_size);
    c_displacements=(int *)malloc(sizeof(int)*comm_size);
    irf_counts=(int *)malloc(sizeof(int)*comm_size);
    irf_displacements=(int *)malloc(sizeof(int)*comm_size);
    icf_counts=(int *)malloc(sizeof(int)*comm_size);
    icf_displacements=(int *)malloc(sizeof(int)*comm_size);
    a_counts[0] = 0;
    a_displacements[0] = 0;
    b_counts[0] = 0;
    b_displacements[0] = 0;
    c_counts[0] = 0;
    c_displacements[0] = 0;
    irf_counts[0] = 0;
    irf_displacements[0] = 0;
    icf_counts[0] = 0;
    icf_displacements[0] = 0;

    for(i=1;i<comm_size;i++){
      
      /*Send message to get worker into conpar mode*/
      {
	int message=AUTO_MPI_CONPAR_MESSAGE;
	MPI_Send(&message,1,MPI_INT,i,0,MPI_COMM_WORLD);
      }
      loop_start = ((i-1)*na)/(comm_size - 1);
      loop_end = ((i)*na)/(comm_size - 1);
      a_counts[i] = nca*nra*(loop_end-loop_start);
      a_displacements[i] = nca*nra*loop_start;
      b_counts[i] = ncb*nra*(loop_end-loop_start);
      b_displacements[i] = ncb*nra*loop_start;
      c_counts[i] = nca*nrc*(loop_end-loop_start);
      c_displacements[i] = nca*nrc*loop_start;
      irf_counts[i] = nra*(loop_end-loop_start);
      irf_displacements[i] = nra*loop_start;
      icf_counts[i] = nca*(loop_end-loop_start);
      icf_displacements[i] = nca*loop_start;
      loop_end_tmp = loop_end-loop_start;
      MPI_Send(&loop_end_tmp   ,1,MPI_LONG,i,0,MPI_COMM_WORLD);
    }
    {
      integer params[5];
      params[0]=nov;
      params[1]=nra;
      params[2]=nca;
      params[3]=ncb;
      params[4]=nrc;

      
      MPI_Bcast(params        ,5,MPI_LONG,0,MPI_COMM_WORLD);
    }

    /* Worker is running now */

    {
      /*I create a temporary send buffer for the MPI_Reduce
	command.  This is because there isn't an
	asymmetric version (like MPI_Scatterv).*/
      double **dtemp;
      dtemp = dmatrix(nrc,ncb);
      for(i=0;i<nrc;i++)
        for(j=0;j<ncb;j++)
          dtemp[i][j]=d[i][j];
      MPI_Reduce(dtemp[0],d[0],ncb*nrc,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      free_dmatrix(dtemp);
    }
    MPI_Gatherv(NULL,0,MPI_DOUBLE,
		a[0][0],a_counts,a_displacements,MPI_DOUBLE,
		0,MPI_COMM_WORLD);
    MPI_Gatherv(NULL,0,MPI_DOUBLE,
		b[0][0],b_counts,b_displacements,MPI_DOUBLE,
		0,MPI_COMM_WORLD);
    MPI_Gatherv(NULL,0,MPI_DOUBLE,
		c[0][0],c_counts,c_displacements,MPI_DOUBLE,
		0,MPI_COMM_WORLD);
    MPI_Gatherv(NULL,0,MPI_LONG,
		irf,irf_counts,irf_displacements,MPI_LONG,
		0,MPI_COMM_WORLD);
    MPI_Gatherv(NULL,0,MPI_LONG,
		icf,icf_counts,icf_displacements,MPI_LONG,
		0,MPI_COMM_WORLD);
    return 0;
}
#endif

int 
conpar(integer nov, integer na, integer nra, integer nca, doublereal ***a, integer ncb, doublereal ***b, integer nrc, doublereal ***c, doublereal **d, integer *irf, integer *icf)
{
  integer nex;

  nex = nca - (nov << 1);
  if (nex == 0) {
    return 0;
  }

  switch(global_conpar_type) {
#ifdef PTHREADS
  case CONPAR_PTHREADS:
    conpar_threads_wrapper(nov, na, nra, nca, a, 
			    ncb, b, nrc, c, d, irf, icf);
    break;
#endif
#ifdef MPI
  case CONPAR_MPI:
    if(global_verbose_flag)
      printf("MPI conpar start\n");
    conpar_mpi_wrapper(nov, na, nra, nca, a, 
			ncb, b, nrc, c, d,irf, icf);
    if(global_verbose_flag)
      printf("MPI conpar end\n");
    break;
#endif
  default:
    conpar_process(nov, na, nra, nca, a, ncb, b, nrc, c, d, irf, icf);
    break;
  }
  return 0;
}
