/*
 *  PROGRAM NAME:  multifario
 *
 *  (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
 *  CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
 *
 *  Please refer to the LICENSE file in the top directory
 *
 *      author: Mike Henderson mhender@watson.ibm.com
 *      date:   July 15, 2002 Modified MFAUTO
 */

static char *id="@(#) $Id: MFAUTOIMF.c,v 1.14 2011/07/21 17:42:46 mhender Exp $";

#include <multifarioConfig.h>

static char MFAUTOMFErrorHandlerMsg[256]="";

#include <MFImplicitMF.h>
#include <MFErrorHandler.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <MFAUTO.h>
#include <MFNSpace.h>
#include <auto_c.h>

#ifdef __cplusplus
 extern "C" {
#endif
int testTangent(doublereal **udot, doublereal *rdot,iap_type *iap, rap_type *rap, doublereal *par, integer *icp,
                FUNI_TYPE((*funi)), BCNI_TYPE ((*bcni)), ICNI_TYPE((*icni)), doublereal *rds,
                integer nllv, doublereal *rlcur, doublereal *rlold, doublereal *rldot, integer ndxloc,
                doublereal **ups, doublereal **dups, doublereal **uoldps,
                doublereal **udotps, doublereal **upoldp, doublereal *dtm,doublereal **fa, doublereal *fc,
                doublereal **p0, doublereal **p1, doublereal *thl, doublereal *thu);

extern user_function_list user;    /* This is where AUTO gets the things to evaluate */

int adaptK(iap_type *iap, integer nold, integer ncold, integer nnew, integer ncnew, doublereal *tm, doublereal *dtm, integer ndxloc, doublereal **ups, doublereal ***vps);

static int MFProjectAUTO(int,int,MFNVector,MFNKMatrix,MFNVector,void*,int*,MFErrorHandler);
static int MFTangentAUTOWithGuess(int,int,MFNVector,MFNKMatrix,MFNKMatrix,void*,MFErrorHandler);
static int MFTangentAUTOOriginal(int,int,MFNVector,MFNKMatrix,MFNKMatrix,void*,MFErrorHandler);
static double MFScaleAUTO(int,int,MFNVector,MFNKMatrix,void*,MFErrorHandler);
static int MFAUTOProjectToSave(MFNVector,double*,void*,MFErrorHandler);
static int MFAUTOProjectToDraw(MFNVector,double*,void*,MFErrorHandler);
static int MFAUTOProjectForBB(MFNVector,double*,void*,MFErrorHandler);
static int MFAUTOStop(MFImplicitMF,MFNVector,MFNKMatrix,MFNVector,MFNKMatrix,void*,MFErrorHandler);
static void MFFreeDataAUTO(void*,MFErrorHandler);

void MFAUTOBVNVCopyDataValues(MFNVector,MFNVector,MFErrorHandler);
void MFAUTOBVNVSetT(MFNVector,doublereal*,MFErrorHandler);
void MFAUTOBVNVSetType(MFNVector,char*,MFErrorHandler);
void MFAUTOBVNVSetThu(MFNVector,doublereal*,MFErrorHandler);

doublereal **MFAUTOBVNVGetU(MFNVector,MFErrorHandler);
doublereal *MFAUTOBVNVGetT(MFNVector,MFErrorHandler);
doublereal *MFAUTOBVNVGetDt(MFNVector,MFErrorHandler);
doublereal *MFAUTOBVNVGetPar(MFNVector,MFErrorHandler);
doublereal MFAUTOBVNVGetR(MFNVector,MFErrorHandler);
int MFAUTOBVNVGetNit(MFNVector,MFErrorHandler);
void MFAUTOBVNVSetNit(MFNVector,int,MFErrorHandler);
int MFAUTOBVNVGetNtst(MFNVector,MFErrorHandler);
int MFAUTOBVNVGetNcol(MFNVector,MFErrorHandler);
int MFAUTOBVNVGetNdim(MFNVector,MFErrorHandler);
int MFAUTOBVNVGetNpar(MFNVector,MFErrorHandler);

void MFAUTOBVNVSetNUZ(MFNVector,int,MFErrorHandler);
void MFAUTOBVNVSetUZBV(MFNVector,int,doublereal,MFErrorHandler);
void MFAUTOBVNVSetLPBV(MFNVector,doublereal,MFErrorHandler);
void MFAUTOBVNVSetBPBV(MFNVector,doublereal,MFErrorHandler);
void MFAUTOBVNVSetSPBV(MFNVector,doublereal,MFErrorHandler);
void MFAUTOBVNVSetEV(MFNVector,doublecomplex*,MFErrorHandler);
void MFAUTOBVNVSetP0(MFNVector,doublereal**,MFErrorHandler);
void MFAUTOBVNVSetP1(MFNVector,doublereal**,MFErrorHandler);

int MFAUTOBVNVGetNUZ(MFNVector,MFErrorHandler);
doublereal MFAUTOBVNVGetUZBV(MFNVector,int,MFErrorHandler);
doublereal MFAUTOBVNVGetLPBV(MFNVector,MFErrorHandler);
doublereal MFAUTOBVNVGetBPBV(MFNVector,MFErrorHandler);
doublereal MFAUTOBVNVGetSPBV(MFNVector,MFErrorHandler);
doublecomplex *MFAUTOBVNVGetEV(MFNVector,MFErrorHandler);
doublereal** MFAUTOBVNVGetP0(MFNVector,MFErrorHandler);
doublereal** MFAUTOBVNVGetP1(MFNVector,MFErrorHandler);

void MFPrintAUTOBVNVectorFull(FILE*,MFNVector,MFErrorHandler);

double MFAtlasDet(int,double*,MFErrorHandler);

MFNVector MFAUTOVectorFactory(MFImplicitMF,MFErrorHandler);
MFNKMatrix MFAUTOMatrixFactory(MFImplicitMF,MFErrorHandler);

doublereal *MFAUTONSpaceGetThl(MFNSpace,MFErrorHandler);
doublereal *MFAUTONSpaceGetThu(MFNSpace,MFErrorHandler);

void MFPrintAUTOBVNVectorFull(FILE*,MFNVector,MFErrorHandler);

void MFAUTOSetLimitPointStop(MFNVector,MFErrorHandler);
int MFAUTOTestLimitPointStop(MFImplicitMF,MFNVector,MFNKMatrix,MFNVector,MFNKMatrix,MFErrorHandler);
void MFAUTOSetBifurcationPointStop(MFNVector,MFErrorHandler);
int MFAUTOTestBifurcationPointStop(MFImplicitMF,MFNVector,MFNKMatrix,MFNVector,MFNKMatrix,MFErrorHandler);
void MFAUTOSetSpecialPointStop(MFNVector,MFErrorHandler);
int MFAUTOTestSpecialPointStop(MFImplicitMF,MFNVector,MFNKMatrix,MFNVector,MFNKMatrix,MFErrorHandler);

double MFAUTOBVGetRealParameter(MFImplicitMF,char*,MFErrorHandler);

extern struct {
  integer irtn;
  integer *nrtn;
} global_rotations;

struct MFAUTOData 
  {
   MFNSpace space;
   iap_type *iap;
   rap_type *rap;
   long *icp;

   MFfunc_type func;
   MFbcnd_type bcnd;
   MFicnd_type icnd;
   MFpvls_type pvls;

   FUNI_TYPE((*funi));
   BCNI_TYPE((*bcni));
   ICNI_TYPE((*icni));
   PVLI_TYPE_BVP((*pvli));
   STPNT_TYPE_BVP((*stpnt));

/* "global memory" */

   double **ups;
   double **uoldps;
   double **dups;
   double **upoldp;
   double **udotps;
   double *rds;
   double **fa;
   double *fc;
   double *rldot;
   double *rlcur;
   double *rlold;
   double *rhs;
   double ds0;

/* "testing for singularity" */

   int nStops;
   int mStops;
   MFAUTOSetStopDataRtn *setStopData;
   MFAUTOTestStopDataRtn *testStopData;
   int nUserZeros;
   int mUserZeros;
   integer *userZeroParm;
   double *userZeroValue;
   int skipFirstFloquetMult;
  };

struct MFAUTOTPBVPSt
  {
   integer k;
   integer ndim;
   integer nbc;
   integer nic;
   integer npar;
   integer nfpr;
   integer jac;
   integer ntst;
   integer ncol;
   integer nicp;
   integer *icp;

   MFfunc_type func;
   MFbcnd_type bcnd;
   MFicnd_type icnd;
   MFpvls_type pvls;

   long nrefs;
  };

extern FILE *fp3;
extern FILE *fp7;
extern FILE *fp9;
extern FILE *fp12;
extern int global_conpar_type;
extern int global_setubv_type;
extern int global_reduce_type;
extern int global_num_procs;
extern int global_verbose_flag;

/*!    \fn MFImplicitMF MFCreateAUTOBV(MFAUTOTPBVP tpbvp, MFNSpace space, MFErrorHandler e);
 *     \brief Creates an implicit representation of the solution manifold of a two point boundary value problem.
 *
 *     \param tpbvp A two point boundary value problem of the type used by AUTO.
 *     \param space The space in which the solution manifold lives.
 *     \param e     An error handler.
 *     \returns The solution manifold. 
 */
MFImplicitMF MFCreateAUTOBV(MFAUTOTPBVP tpbvp, MFNSpace space, MFErrorHandler e)
 {
  static char RoutineName[]={"MFCreateAUTOBV"};
  MFImplicitMF Auto;

#ifdef HAVE_AUTO

  struct MFAUTOData *data;
  int ndim,ncol,ntst,npar;
  int nint,n;
  function_list *list;
  iap_type *iap;
  rap_type *rap;
  integer iad;
  integer jac;
  integer nicp;
  integer i;
  integer *icp;
  doublereal rl0;
  doublereal rl1;
  doublereal a0;
  doublereal a1;
  integer ips;
  integer isw;
  integer irs;
  integer ilp;
  integer isp;
  integer iplt;
  integer npr;
  integer iid;
  integer nmx;
  integer mxbf;
  integer itmx;
  integer itnw;
  integer nwtn;
  doublereal epsl;
  doublereal epsu;
  doublereal epss;
  doublereal ds;
  doublereal dsmin;
  doublereal dsmax;
  integer iads;
  integer nuzr;
  integer *iuz;
  doublereal *vuz;
  integer nalc;
  integer nbc ;
  integer nfpr;
  int verbose=0;
  doublereal amp;
  doublereal det;
  doublereal tivp;
  doublereal fldf;
  doublereal hbff;
  doublereal biff;
  doublereal spbf;

  if(fp3==NULL)fp3=fopen("auto.3","w");
  if(fp7==NULL)fp7=fopen("auto.7","w");
  if(fp9==NULL)fp9=fopen("auto.9","w");
  if(fp12==NULL)fp12=fopen("auto.12","w");

  iap=malloc(sizeof(iap_type));

#ifndef MFNOSAFETYNET
  if(iap==NULL)
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Out of memory trying to allocate %d bytes.",sizeof(iap_type));
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  rap=malloc(sizeof(rap_type));

#ifndef MFNOSAFETYNET
  if(rap==NULL)
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Out of memory trying to allocate %d bytes.",sizeof(rap_type));
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  list=malloc(sizeof(function_list));

#ifndef MFNOSAFETYNET
  if(list==NULL)
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Out of memory trying to allocate %d bytes.",sizeof(function_list));
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

/* Read constants from file and allocate thu. */
/*  init(&iap, &rap, par, icp, thl, thu, &iuz, &vuz, &eof,k);*/

/* Set things normally read from file */

/* Problem stuff */

  ndim=MFAUTOTPBVPGetNDIM(tpbvp,e);
  ntst=MFAUTOTPBVPGetNTST(tpbvp,e);
  ncol=MFAUTOTPBVPGetNCOL(tpbvp,e);
  npar=MFAUTOTPBVPGetNPAR(tpbvp,e);
  nalc=MFAUTOTPBVPGetK(tpbvp,e);
  nbc =MFAUTOTPBVPGetNBC(tpbvp,e);
  nint=MFAUTOTPBVPGetNIC(tpbvp,e);

  iad=3;                                    /* 0=fixed mesh, >0 number of steps before adapt */
  jac=MFAUTOTPBVPGetJAC(tpbvp,e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s\n",RoutineName);fflush(stdout);}
#endif

  icp=malloc(2*NPARX*sizeof(integer));

#ifndef MFNOSAFETYNET
  if(icp==NULL)
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Out of memory trying to allocate %d bytes.",2*NPARX*sizeof(integer));
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  nicp=tpbvp->nicp;
  for(i=0;i<nicp;i++)icp[i]=(MFAUTOTPBVPGetICP(tpbvp,e))[i];  /* List of available parameters (permutation) */

#ifdef MFALLOWVERBOSE
  if(verbose){printf("icp = [");for(i=0;i<nicp;i++){if(i>0)printf(",");printf("%d",icp[i]);}printf("]\n");fflush(stdout);}
#endif

/* Omega (not used) */

  rl0=-1.e6;
  rl1=1.e6;
  a0=-1.e6;
  a1=1.e6;

/* ** Boundary value problems. */

  ips=4;                 /* Type of Problem */

                         /*    ips= 1    Stationary Solutions of ODE's w/ Hopf bifs */
                         /*    ips=-1;   Fixed Points of maps w/ Hopf bifs */
                         /*    ips=-2;   Time integration with Euler */
                         /*    ips= 2;   Computation of Periodic Solutions */
                         /*    ips= 4;   Boundary Value Problem */
                         /*    ips= 5;   Algebraic Optimization Problem*/
                         /*    ips= 7;   Boundary Value Problem with Floquet Multipliers */
                         /*    ips= 9;   For detection and continuation of homoclinic bifurcations */
                         /*    ips=11;   Spatially uniform solutions of a system of parabolic PDEs */
                         /*    ips=12;   Continuation of traveling waave solution of a system of parabolic PDEs */
                         /*    ips=14;   Time evolution of traveling wave solutions to a system of parabolic PDEs */
                         /*    ips=15;   Optimization of periodic solutions */
                         /*    ips=16;   Similar to 14 but with boundary conditions */
                         /*    ips=17;   Stationary Solutions of parabolic systems */

  isw=1;                 /* If restart use bifurcation point (branch switch) */
  irs=0;                 /* Restart 0=New problem, >0 label of point */

/* Detection of Bifurcations */

  ilp=0;                /* 0=no detection of limit points, >0=locat limit points */
  isp=1;                /* 0=no detection, 1,2,3 = detection (for TPBVP) */

  iplt=0;    /* For plotting projection 0=L_2 norm */
  npr=20;    /* Interval for printing out point */
  iid= 0;    /* Output level */

  nmx=100;   /* Maximum number of steps to take along a branch */
  mxbf=5;    /* Maximum number of biifurcating branches to trace out */

/* For the projection */

  itmx=8;    /* Maximum iterations to locate a bifurcation point */
  itnw=5;    /* Maximum iterations for projection */
  nwtn=3;    /* Freeze Jacobian after nwtn iterations */

/* Tolerances */

  epsl=1.e-4;  /* Relative tolerance for PAR in Newton */
  epsu=1.e-4;  /* Relative tolerance for U in Newton */
  epss=1.e-4;  /* Relative tolerance for arclength in detecting bifurcations */

  if(epsl<0.0)
   {
    printf("Warning : EPSL less then 0.0, will use absolute value instead.");
    epsl=fabs(epsl);
   }
  if(epsu<0.0)
   {
    printf("Warning : EPSU less then 0.0, will use absolute value instead.");
    epsu=fabs(epsu);
   }
  if(epss < 0.0)
   {
    printf("Warning : EPSS less then 0.0, will use absolute value instead.");
    epss=fabs(epss);
   }

/* Stepsize */

  ds=0.01;     /* Stepsize */
  dsmin=0.001; /* Minimum Stepsize */
  dsmax=2.0;   /* Maximum Stepsize */
  iads=1;      /* Adapt the stepsize after this many steps */

/* User zeroes */
/*       These should be added to the continuation method as "stop" functions */

  nuzr=0;
  iuz=NULL;
  vuz=NULL;
/*
  ds=MFAUTOBVGetRealParameter(tpbvp,"ds",e);
  dsmin=MFAUTOBVGetRealParameter(tpbvp,"dsmin",e);
  dsmax=MFAUTOBVGetRealParameter(tpbvp,"dsmax",e);
  amp=MFAUTOBVGetRealParameter(tpbvp,"amp",e);
  epsl=MFAUTOBVGetRealParameter(tpbvp,"epsl",e);
  epsu=MFAUTOBVGetRealParameter(tpbvp,"epsu",e);
  epss=MFAUTOBVGetRealParameter(tpbvp,"epss",e);
  det=MFAUTOBVGetRealParameter(tpbvp,"det",e);
  tivp=MFAUTOBVGetRealParameter(tpbvp,"tivp",e);
  fldf=MFAUTOBVGetRealParameter(tpbvp,"fldf",e);
  hbff=MFAUTOBVGetRealParameter(tpbvp,"hbff",e);
  biff=MFAUTOBVGetRealParameter(tpbvp,"biff",e);
  spbf=MFAUTOBVGetRealParameter(tpbvp,"spbf",e);
 */
  if(dsmin<0.0)
   {
    printf("Warning : DSMIN less then 0.0, will use absolute value instead.");
    dsmin = fabs(dsmin);
   }

  if(dsmax<0.0)
   {
    printf("Warning : DSMAX less then 0.0, will use absolute value instead.");
    dsmax = fabs(dsmax);
   }

/* Copy them into iap and rap */

  iap->ndim = ndim;
  iap->ips = ips;
  iap->irs = irs;
  iap->ilp = ilp;
  iap->ntst = ntst;
  iap->ncol = ncol;
  iap->iad = iad;
  iap->iads = iads;
  iap->isp = isp;
  iap->isw = isw;
  iap->iplt = iplt;
  iap->nbc = nbc;
  iap->nint = nint;
  iap->nalc = nalc;
  iap->nmx = nmx;
  iap->nuzr = nuzr;
  iap->npr = npr;
  iap->mxbf = mxbf;
  iap->iid = iid;
  iap->itmx = itmx;
  iap->itnw = itnw;
  iap->nwtn = nwtn;
  iap->jac = jac;

  iap->ndm = ndim;
  iap->nbc0 = 1; if(nbc!=0)iap->nbc0 = nbc;
  iap->nnt0 = 1; if(nint!=0)iap->nnt0 = nint;
  iap->iuzr = 1;
  iap->itp = 0;
  iap->itpst = 0;
  iap->ibr = 1;
  iap->nit = 0;
  iap->ntot = 0;
  iap->nins = 0;
  iap->istop = 0;
  iap->nbif = 0;
  iap->ipos = 1;
  iap->lab = 0;
  iap->nicp = nicp;

  iap->mynode = 0;
  iap->numnodes = 1;
  iap->parallel_flag = 0;

  rap->ds = dsmax;
  rap->dsmin = dsmin;
  rap->dsmax = dsmax;
  rap->dsold = ds;

  rap->rl0 = rl0;
  rap->rl1 = rl1;
  rap->a0 = a0;
  rap->a1 = a1;

  rap->amp = amp;
  rap->epsl = epsl;
  rap->epsu = epsu;
  rap->epss = epss;
  rap->det = det;
  rap->tivp = tivp;
  rap->fldf = fldf;
  rap->hbff = hbff;
  rap->biff = biff;
  rap->spbf = spbf;

/* End of init() */

/* Analyze iap (type of problem) and set functions in list */
/*  set_function_pointers(iap,&list); */

/*  This bit depends on the problem type . */
/*  ** Boundary value problems. */

  list->type           = AUTOBV;
  list->bvlist.funi  = funi;
  list->bvlist.bcni  = bcni;
  list->bvlist.icni  = icni;
  list->bvlist.stpnt = stpnub;
  list->bvlist.pvli  = pvlsbv;

/*  end of set_function_pointers(); */

/*init1(&iap, &rap, icp, par);*/

  if(iap->isw==0)iap->isw=1;

  if(rap->ds==0.)rap->ds=(double).1;
  if(rap->dsmin==0.)rap->dsmin=fabs(rap->ds)*1e-4;
  rap->ds=HMACH1*rap->ds;
  rap->dsmin/=HMACH1;
  rap->dsmax=HMACH1*rap->dsmax;

/*  This bit depends on the problem type . */
/*        ** Boundary value problems */

  nfpr=iap->nbc+iap->nint-iap->ndim+iap->nalc;
  iap->nfpr=nfpr;

/* end of init1();*/

/*chdim(&iap);   This just makes sure that nfpr<NPARX */

/* These were from StartOfCNRL */

  rap->dsold=rap->ds;
  iap->isp=abs(iap->isp);
  iap->nit=0;
  iap->ntot=0;
  iap->istop=0;

/* Now have all the details of the problem */

  n=(iap->ntst+1)*iap->ndim*iap->ncol+iap->nfpr;

  nbc=iap->nbc;
  nint=iap->nint;
  ntst=iap->ntst;
  ncol=iap->ncol;
  ndim=iap->ndim;

  data=malloc(sizeof(struct MFAUTOData));

#ifndef MFNOSAFETYNET
  if(data==NULL)
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Out of memory trying to allocate %d bytes.",sizeof(struct MFAUTOData));
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  data->space=space;MFRefNSpace(space,e);
  data->iap=iap;
  data->rap=rap;
  data->icp=icp;
  data->funi=list->bvlist.funi;
  data->bcni=list->bvlist.bcni;
  data->icni=list->bvlist.icni;
  data->pvli=list->bvlist.pvli;
  data->stpnt=list->bvlist.stpnt;

  data->ups=dmatrix(ntst+1,ndim*ncol);
#ifndef MFNOSAFETYNET
  if(data->ups==NULL)
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Out of memory trying to allocate %d bytes.",(ntst+1)*ndim*ncol*sizeof(double));
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  data->uoldps=dmatrix(ntst+1,ndim*ncol);
#ifndef MFNOSAFETYNET
  if(data->uoldps==NULL)
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Out of memory trying to allocate %d bytes.",(ntst+1)*ndim*ncol*sizeof(double));
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  data->dups=dmatrix(ntst+1,ndim*ncol);
#ifndef MFNOSAFETYNET
  if(data->dups==NULL)
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Out of memory trying to allocate %d bytes.",(ntst+1)*ndim*ncol*sizeof(double));
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif


  data->upoldp=dmatrix(ntst+1,ndim*ncol);

#ifndef MFNOSAFETYNET
  if(data->upoldp==NULL)
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Out of memory trying to allocate %d bytes.",(ntst+1)*ndim*ncol*sizeof(double));
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  data->udotps=dmatrix(NPARX*(ntst+1),ndim*ncol);

#ifndef MFNOSAFETYNET
  if(data->udotps==NULL)
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Out of memory trying to allocate %d bytes.",NPARX*(ntst+1)*ndim*ncol*sizeof(double));
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  data->rds=malloc(NPARX*sizeof(double));

#ifndef MFNOSAFETYNET
  if(data->rds==NULL)
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Out of memory trying to allocate %d bytes.",NPARX*sizeof(double));
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  data->rldot=malloc(10*NPARX*sizeof(double));

#ifndef MFNOSAFETYNET
  if(data->rldot==NULL)
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Out of memory trying to allocate %d bytes.",NPARX*sizeof(double));
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  data->rlcur=malloc(NPARX*sizeof(double));

#ifndef MFNOSAFETYNET
  if(data->rlcur==NULL)
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Out of memory trying to allocate %d bytes.",NPARX*sizeof(double));
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  data->rlold=malloc(NPARX*sizeof(double));

#ifndef MFNOSAFETYNET
  if(data->rlold==NULL)
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Out of memory trying to allocate %d bytes.",NPARX*sizeof(double));
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif
  data->rhs=malloc(nalc*sizeof(double));

#ifndef MFNOSAFETYNET
  if(data->rhs==NULL)
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Out of memory trying to allocate %d bytes.",(ndim+nalc)*sizeof(double));
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  data->fa=dmatrix(ntst+1,ncol*ndim);

#ifndef MFNOSAFETYNET
  if(data->fa==NULL)
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Out of memory trying to allocate %d bytes.",(ntst+1)*ndim*ncol*sizeof(double));
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  data->fc=malloc((ndim+iap->nbc+iap->nint+iap->nalc)*sizeof(double));

#ifndef MFNOSAFETYNET
  if(data->fc==NULL)
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Out of memory trying to allocate %d bytes.",(ndim+nbc+nint+nalc)*sizeof(double));
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  data->func=MFAUTOTPBVPGetFUNC(tpbvp,e);
  data->bcnd=MFAUTOTPBVPGetBCND(tpbvp,e);
  data->icnd=MFAUTOTPBVPGetICND(tpbvp,e);
  data->pvls=MFAUTOTPBVPGetPVLS(tpbvp,e);

  data->nStops=0;
  data->mStops=0;
  data->setStopData=NULL;
  data->testStopData=NULL;
  data->nUserZeros=0;
  data->mUserZeros=0;
  data->userZeroParm=NULL;
  data->userZeroValue=NULL;
  data->skipFirstFloquetMult=0;
  data->ds0=data->rap->ds;

/*  allocate_global_memory(iap);*/

  if(global_rotations.nrtn!=NULL)free(global_rotations.nrtn);

  global_rotations.nrtn = malloc(sizeof(integer)*(iap->ndim));
  global_rotations.irtn=0;
  for(i=0;i<iap->ndim;i++)global_rotations.nrtn[i]=0;

#ifndef MFNOSAFETYNET
  if(global_rotations.nrtn==NULL)
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Out of memory trying to allocate %d bytes.",sizeof(integer)*(iap->nbc));
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

/*  end of allocate_global_memory();*/

  Auto=MFIMFCreateBaseClass(n,nalc,"AUTOBV",e);

  MFIMFSetSpace(Auto,space,e);
  MFIMFSetData(Auto,(void*)data,e);
  MFIMFSetFreeData(Auto,MFFreeDataAUTO,e);
  MFIMFSetProject(Auto,MFProjectAUTO,e);
  MFIMFSetTangentWithGuess(Auto,MFTangentAUTOOriginal,e);  /* **** */
  MFIMFSetTangentWithGuess(Auto,MFTangentAUTOWithGuess,e);
  MFIMFSetScale(Auto,MFScaleAUTO,e);
  MFIMFSetStop(Auto,MFAUTOStop,e);
  MFIMFSetProjectForSave(Auto,MFAUTOProjectToSave,e);
  MFIMFSetProjectForDraw(Auto,MFAUTOProjectToDraw,e);
  MFIMFSetProjectForBB(Auto,MFAUTOProjectForBB,e);
  MFIMFSetRMin(Auto,rap->dsmin,e);
  MFIMFSetVectorFactory(Auto,MFAUTOVectorFactory,e);
  MFIMFSetMatrixFactory(Auto,MFAUTOMatrixFactory,e);

  free(list);

#else
  Auto=NULL;
  sprintf(MFAUTOMFErrorHandlerMsg,"%s called, but AUTO was not found when multifario was installed.",RoutineName);
  MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
#endif

  return Auto;
 }

/*!    \fn MFImplicitMF MFCreateAUTOPeriodicOrbit(MFAUTOTPBVP tpbvp, MFNSpace space, MFerrorHandler e);
 *     \brief Creates an implicit representation of the solution manifold of a Periodic Orbit
 *
 *     \param tpbvp A two point boundary value problem of the type used by AUTO.
 *     \param space The space in which the solution manifold lives.
 *     \param e     An error handler.
 *     \returns The solution manifold. 
 */
MFImplicitMF MFCreateAUTOPeriodicOrbit(MFAUTOTPBVP tpbvp, MFNSpace space, MFErrorHandler e)
 {
  static char RoutineName[]={"MFCreateAUTOPeriodicOrbit"};
  MFImplicitMF Auto;

#ifdef HAVE_AUTO

  struct MFAUTOData *data;
  int ndim,ncol,ntst,npar;
  int nint,n;
  function_list *list;
  iap_type *iap;
  rap_type *rap;
  integer iad;
  integer jac;
  integer nicp;
  integer i;
  integer *icp;
  doublereal rl0;
  doublereal rl1;
  doublereal a0;
  doublereal a1;
  integer ips;
  integer isw;
  integer irs;
  integer ilp;
  integer isp;
  integer iplt;
  integer npr;
  integer iid;
  integer nmx;
  integer mxbf;
  integer itmx;
  integer itnw;
  integer nwtn;
  doublereal epsl;
  doublereal epsu;
  doublereal epss;
  doublereal ds;
  doublereal dsmin;
  doublereal dsmax;
  integer iads;
  integer nuzr;
  integer *iuz;
  doublereal *vuz;
  integer nalc;
  integer nbc ;
  integer nfpr;
  int verbose=0;
  doublereal amp;
  doublereal det;
  doublereal tivp;
  doublereal fldf;
  doublereal hbff;
  doublereal biff;
  doublereal spbf;

  if(fp3==NULL)fp3=fopen("auto.3","w");
  if(fp7==NULL)fp7=fopen("auto.7","w");
  if(fp9==NULL)fp9=fopen("auto.9","w");
  if(fp12==NULL)fp12=fopen("auto.12","w");

  iap=malloc(sizeof(iap_type));

#ifndef MFNOSAFETYNET
  if(iap==NULL)
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Out of memory trying to allocate %d bytes.",sizeof(iap_type));
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  rap=malloc(sizeof(rap_type));

#ifndef MFNOSAFETYNET
  if(rap==NULL)
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Out of memory trying to allocate %d bytes.",sizeof(rap_type));
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  list=malloc(sizeof(function_list));

#ifndef MFNOSAFETYNET
  if(list==NULL)
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Out of memory trying to allocate %d bytes.",sizeof(function_list));
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

/* Read constants from file and allocate thu. */
/*  init(&iap, &rap, par, icp, thl, thu, &iuz, &vuz, &eof,k);*/

/* Set things normally read from file */

/* Problem stuff */

  ndim=MFAUTOTPBVPGetNDIM(tpbvp,e);
  ntst=MFAUTOTPBVPGetNTST(tpbvp,e);
  ncol=MFAUTOTPBVPGetNCOL(tpbvp,e);
  npar=MFAUTOTPBVPGetNPAR(tpbvp,e);
  nalc=MFAUTOTPBVPGetK(tpbvp,e);
  nbc =MFAUTOTPBVPGetNBC(tpbvp,e);
  nint=MFAUTOTPBVPGetNIC(tpbvp,e);

  iad=3;                                    /* 0=fixed mesh, >0 number of steps before adapt */
  jac=MFAUTOTPBVPGetJAC(tpbvp,e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s\n",RoutineName);fflush(stdout);}
#endif

  icp=malloc(2*NPARX*sizeof(integer));

#ifndef MFNOSAFETYNET
  if(icp==NULL)
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Out of memory trying to allocate %d bytes.",2*NPARX*sizeof(integer));
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  nicp=tpbvp->nicp;
  for(i=0;i<nicp;i++)icp[i]=(MFAUTOTPBVPGetICP(tpbvp,e))[i];  /* List of available parameters (permutation) */

#ifdef MFALLOWVERBOSE
  if(verbose){printf("icp = [");for(i=0;i<nicp;i++){if(i>0)printf(",");printf("%d",icp[i]);}printf("]\n");fflush(stdout);}
#endif

/* Omega (not used) */

  rl0=-1.e6;
  rl1=1.e6;
  a0=-1.e6;
  a1=1.e6;

/* ** Boundary value problems. */

#if 0

  } else if (ips == -2) {
    /*        ** Time integration */
#ifndef MANIFOLD
    nfpr = 1;
#else
    nfpr = iap->nalc;
#endif
    isp = 0;
    ilp = 0;
    icp = realloc(icp, nfpr * sizeof(*icp));
    icp[0] = 13;

  } else if (ips == 2 && abs(isw) == 1) {
    /*        ** Periodic Solutions */
    nbc = ndim;
    nint = 1;
#ifndef MANIFOLD
    nfpr = nbc + nint - ndim + 1;
#else
    nfpr = nbc + nint - ndim + iap->nalc;
#endif
    icp = realloc(icp, nfpr * sizeof(*icp));
    /*        **ISW=1 when starting from a HB */
    if (itp == 3 || abs(itp) / 10 == 3) {
      isw = 1;
    }
    if (nicp == 1) {
      /*          **Variable period */
      icp[1] = 10;
    }

  } else if (iap.ips == 2 && abs(iap.isw) != 2) {
    /*  ** Periodic solutions */
    if (iap.itp != 3 && abs(iap.itp / 10) != 3) {
      if (iap.irs > 0) {
        data->type         = AUTOBV;
        data->bvlist.funi  = fnps;
        data->bvlist.bcni  = bcps;
        data->bvlist.icni  = icps;
        data->bvlist.stpnt = stpnbv;
        data->bvlist.pvli  = pvlsbv;
        if (iap.isw == -1 && iap.itp == 7) {
          data->bvlist.stpnt = stpnpdble;
        }
      } else {
        data->type         = AUTOBV;
        data->bvlist.funi  = fnps;
        data->bvlist.bcni  = bcps;
        data->bvlist.icni  = icps;
        data->bvlist.stpnt = stpnub;
        data->bvlist.pvli  = pvlsbv;
      }
    } else {
      data->type           = AUTOBV;
      data->bvlist.funi  = fnps;
      data->bvlist.bcni  = bcps;
      data->bvlist.icni  = icps;
      data->bvlist.stpnt = stpnps;
      data->bvlist.pvli  = pvlsbv;
    }
#endif

  ips=2;                 /* Type of Problem */
                         /*    ips= 0;   Algebraic */

                         /*    ips= 1    Stationary Solutions of ODE's w/ Hopf bifs */
                         /*    ips=-1;   Fixed Points of maps w/ Hopf bifs */
                         /*    ips=-2;   Time integration with Euler */
                         /*    ips= 2;   Computation of Periodic Solutions */
                         /*    ips= 4;   Boundary Value Problem */
                         /*    ips= 5;   Algebraic Optimization Problem*/
                         /*    ips= 7;   Boundary Value Problem with Floquet Multipliers */
                         /*    ips= 9;   For detection and continuation of homoclinic bifurcations */
                         /*    ips=11;   Spatially uniform solutions of a system of parabolic PDEs */
                         /*    ips=12;   Continuation of traveling waave solution of a system of parabolic PDEs */
                         /*    ips=14;   Time evolution of traveling wave solutions to a system of parabolic PDEs */
                         /*    ips=15;   Optimization of periodic solutions */
                         /*    ips=16;   Similar to 14 but with boundary conditions */
                         /*    ips=17;   Stationary Solutions of parabolic systems */

  isw=1;                 /* If restart use bifurcation point (branch switch) */
  irs=0;                 /* Restart 0=New problem, >0 label of point */

  nbc = ndim;
  nint = 1;
  nfpr = nbc + nint - ndim + nalc;

  tpbvp->nbc=nbc;
  tpbvp->nic=nint;
  tpbvp->nfpr=nfpr;

  iad=3;                                    /* 0=fixed mesh, >0 number of steps before adapt */
  jac=MFAUTOTPBVPGetJAC(tpbvp,e);

/* Detection of Bifurcations */

  ilp=0;                /* 0=no detection of limit points, >0=locat limit points */
  isp=1;                /* 0=no detection, 1,2,3 = detection (for TPBVP) */

  iplt=0;    /* For plotting projection 0=L_2 norm */
  npr=20;    /* Interval for printing out point */
  iid= 0;    /* Output level */

  nmx=100;   /* Maximum number of steps to take along a branch */
  mxbf=5;    /* Maximum number of biifurcating branches to trace out */

/* For the projection */

  itmx=8;    /* Maximum iterations to locate a bifurcation point */
  itnw=5;    /* Maximum iterations for projection */
  nwtn=3;    /* Freeze Jacobian after nwtn iterations */

/* Tolerances */

  epsl=1.e-4;  /* Relative tolerance for PAR in Newton */
  epsu=1.e-4;  /* Relative tolerance for U in Newton */
  epss=1.e-4;  /* Relative tolerance for arclength in detecting bifurcations */

  if(epsl<0.0)
   {
    printf("Warning : EPSL less then 0.0, will use absolute value instead.");
    epsl=fabs(epsl);
   }
  if(epsu<0.0)
   {
    printf("Warning : EPSU less then 0.0, will use absolute value instead.");
    epsu=fabs(epsu);
   }
  if(epss < 0.0)
   {
    printf("Warning : EPSS less then 0.0, will use absolute value instead.");
    epss=fabs(epss);
   }

/* Stepsize */

  ds=0.01;     /* Stepsize */
  dsmin=0.001; /* Minimum Stepsize */
  dsmax=2.0;   /* Maximum Stepsize */
  iads=1;      /* Adapt the stepsize after this many steps */

/* User zeroes */
/*       These should be added to the continuation method as "stop" functions */

  nuzr=0;
  iuz=NULL;
  vuz=NULL;
/*
  ds=MFAUTOBVGetRealParameter(tpbvp,"ds",e);
  dsmin=MFAUTOBVGetRealParameter(tpbvp,"dsmin",e);
  dsmax=MFAUTOBVGetRealParameter(tpbvp,"dsmax",e);
  amp=MFAUTOBVGetRealParameter(tpbvp,"amp",e);
  epsl=MFAUTOBVGetRealParameter(tpbvp,"epsl",e);
  epsu=MFAUTOBVGetRealParameter(tpbvp,"epsu",e);
  epss=MFAUTOBVGetRealParameter(tpbvp,"epss",e);
  det=MFAUTOBVGetRealParameter(tpbvp,"det",e);
  tivp=MFAUTOBVGetRealParameter(tpbvp,"tivp",e);
  fldf=MFAUTOBVGetRealParameter(tpbvp,"fldf",e);
  hbff=MFAUTOBVGetRealParameter(tpbvp,"hbff",e);
  biff=MFAUTOBVGetRealParameter(tpbvp,"biff",e);
  spbf=MFAUTOBVGetRealParameter(tpbvp,"spbf",e);
 */
  if(dsmin<0.0)
   {
    printf("Warning : DSMIN less then 0.0, will use absolute value instead.");
    dsmin = fabs(dsmin);
   }

  if(dsmax<0.0)
   {
    printf("Warning : DSMAX less then 0.0, will use absolute value instead.");
    dsmax = fabs(dsmax);
   }

/* Copy them into iap and rap */

  iap->ndim = ndim;
  iap->ips = ips;
  iap->irs = irs;
  iap->ilp = ilp;
  iap->ntst = ntst;
  iap->ncol = ncol;
  iap->iad = iad;
  iap->iads = iads;
  iap->isp = isp;
  iap->isw = isw;
  iap->iplt = iplt;
  iap->nbc = nbc;
  iap->nint = nint;
  iap->nalc = nalc;
  iap->nmx = nmx;
  iap->nuzr = nuzr;
  iap->npr = npr;
  iap->mxbf = mxbf;
  iap->iid = iid;
  iap->itmx = itmx;
  iap->itnw = itnw;
  iap->nwtn = nwtn;
  iap->jac = jac;

  iap->ndm = ndim;
  iap->nbc0 = 1; if(nbc!=0)iap->nbc0 = nbc;
  iap->nnt0 = 1; if(nint!=0)iap->nnt0 = nint;
  iap->iuzr = 1;
  iap->itp = 0;
  iap->itpst = 0;
  iap->ibr = 1;
  iap->nit = 0;
  iap->ntot = 0;
  iap->nins = 0;
  iap->istop = 0;
  iap->nbif = 0;
  iap->ipos = 1;
  iap->lab = 0;
  iap->nicp = nicp;

  iap->mynode = 0;
  iap->numnodes = 1;
  iap->parallel_flag = 0;

  rap->ds = dsmax;
  rap->dsmin = dsmin;
  rap->dsmax = dsmax;
  rap->dsold = ds;

  rap->rl0 = rl0;
  rap->rl1 = rl1;
  rap->a0 = a0;
  rap->a1 = a1;

  rap->amp = amp;
  rap->epsl = epsl;
  rap->epsu = epsu;
  rap->epss = epss;
  rap->det = det;
  rap->tivp = tivp;
  rap->fldf = fldf;
  rap->hbff = hbff;
  rap->biff = biff;
  rap->spbf = spbf;

/* End of init() */

/* Analyze iap (type of problem) and set functions in list */
/*  set_function_pointers(iap,&list); */

/*  This bit depends on the problem type . */

/*  ** Periodic solutions */
  if (iap->itp != 3 && abs(iap->itp / 10) != 3)
   {
    if (iap->irs > 0)
     {
      list->type         = AUTOBV;
      list->bvlist.funi  = fnps;
      list->bvlist.bcni  = bcps;
      list->bvlist.icni  = icps;
      list->bvlist.stpnt = stpnbv;
      list->bvlist.pvli  = pvlsbv;
      if (iap->isw == -1 && iap->itp == 7)
       {
        list->bvlist.stpnt = stpnpdble;
       }
     }else{
      list->type         = AUTOBV;
      list->bvlist.funi  = fnps;
      list->bvlist.bcni  = bcps;
      list->bvlist.icni  = icps;
      list->bvlist.stpnt = stpnub;
      list->bvlist.pvli  = pvlsbv;
     }
   }else{
    list->type           = AUTOBV;
    list->bvlist.funi  = fnps;
    list->bvlist.bcni  = bcps;
    list->bvlist.icni  = icps;
    list->bvlist.stpnt = stpnps;
    list->bvlist.pvli  = pvlsbv;
   }

/*  end of set_function_pointers(); */

/*init1(&iap, &rap, icp, par);*/

  if(iap->isw==0)iap->isw=1;

  if(rap->ds==0.)rap->ds=(double).1;
  if(rap->dsmin==0.)rap->dsmin=fabs(rap->ds)*1e-4;
  rap->ds=HMACH1*rap->ds;
  rap->dsmin/=HMACH1;
  rap->dsmax=HMACH1*rap->dsmax;

/*  This bit depends on the problem type . */
/*        ** Boundary value problems */

  nfpr=iap->nbc+iap->nint-iap->ndim+iap->nalc;
  iap->nfpr=nfpr;

/* end of init1();*/

/*chdim(&iap);   This just makes sure that nfpr<NPARX */

/* These were from StartOfCNRL */

  rap->dsold=rap->ds;
  iap->isp=abs(iap->isp);
  iap->nit=0;
  iap->ntot=0;
  iap->istop=0;

/* Now have all the details of the problem */

  n=(iap->ntst+1)*iap->ndim*iap->ncol+iap->nfpr;

  nbc=iap->nbc;
  nint=iap->nint;
  ntst=iap->ntst;
  ncol=iap->ncol;
  ndim=iap->ndim;

  data=malloc(sizeof(struct MFAUTOData));

#ifndef MFNOSAFETYNET
  if(data==NULL)
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Out of memory trying to allocate %d bytes.",sizeof(struct MFAUTOData));
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  data->space=space;MFRefNSpace(space,e);
  data->iap=iap;
  data->rap=rap;
  data->icp=icp;
  data->funi=list->bvlist.funi;
  data->bcni=list->bvlist.bcni;
  data->icni=list->bvlist.icni;
  data->pvli=list->bvlist.pvli;
  data->stpnt=list->bvlist.stpnt;

  data->ups=dmatrix(ntst+1,ndim*ncol);
#ifndef MFNOSAFETYNET
  if(data->ups==NULL)
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Out of memory trying to allocate %d bytes.",(ntst+1)*ndim*ncol*sizeof(double));
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  data->uoldps=dmatrix(ntst+1,ndim*ncol);
#ifndef MFNOSAFETYNET
  if(data->uoldps==NULL)
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Out of memory trying to allocate %d bytes.",(ntst+1)*ndim*ncol*sizeof(double));
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  data->dups=dmatrix(ntst+1,ndim*ncol);
#ifndef MFNOSAFETYNET
  if(data->dups==NULL)
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Out of memory trying to allocate %d bytes.",(ntst+1)*ndim*ncol*sizeof(double));
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif


  data->upoldp=dmatrix(ntst+1,ndim*ncol);

#ifndef MFNOSAFETYNET
  if(data->upoldp==NULL)
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Out of memory trying to allocate %d bytes.",(ntst+1)*ndim*ncol*sizeof(double));
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  data->udotps=dmatrix(NPARX*(ntst+1),ndim*ncol);

#ifndef MFNOSAFETYNET
  if(data->udotps==NULL)
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Out of memory trying to allocate %d bytes.",NPARX*(ntst+1)*ndim*ncol*sizeof(double));
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  data->rds=malloc(NPARX*sizeof(double));

#ifndef MFNOSAFETYNET
  if(data->rds==NULL)
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Out of memory trying to allocate %d bytes.",NPARX*sizeof(double));
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  data->rldot=malloc(10*NPARX*sizeof(double));

#ifndef MFNOSAFETYNET
  if(data->rldot==NULL)
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Out of memory trying to allocate %d bytes.",NPARX*sizeof(double));
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  data->rlcur=malloc(NPARX*sizeof(double));

#ifndef MFNOSAFETYNET
  if(data->rlcur==NULL)
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Out of memory trying to allocate %d bytes.",NPARX*sizeof(double));
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  data->rlold=malloc(NPARX*sizeof(double));

#ifndef MFNOSAFETYNET
  if(data->rlold==NULL)
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Out of memory trying to allocate %d bytes.",NPARX*sizeof(double));
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif
  data->rhs=malloc(nalc*sizeof(double));

#ifndef MFNOSAFETYNET
  if(data->rhs==NULL)
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Out of memory trying to allocate %d bytes.",(ndim+nalc)*sizeof(double));
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  data->fa=dmatrix(ntst+1,ncol*ndim);

#ifndef MFNOSAFETYNET
  if(data->fa==NULL)
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Out of memory trying to allocate %d bytes.",(ntst+1)*ndim*ncol*sizeof(double));
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  data->fc=malloc((ndim+iap->nbc+iap->nint+iap->nalc)*sizeof(double));

#ifndef MFNOSAFETYNET
  if(data->fc==NULL)
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Out of memory trying to allocate %d bytes.",(ndim+nbc+nint+nalc)*sizeof(double));
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  data->func=MFAUTOTPBVPGetFUNC(tpbvp,e);
  data->bcnd=MFAUTOTPBVPGetBCND(tpbvp,e);
  data->icnd=MFAUTOTPBVPGetICND(tpbvp,e);
  data->pvls=MFAUTOTPBVPGetPVLS(tpbvp,e);

  data->nStops=0;
  data->mStops=0;
  data->setStopData=NULL;
  data->testStopData=NULL;
  data->nUserZeros=0;
  data->mUserZeros=0;
  data->userZeroParm=NULL;
  data->userZeroValue=NULL;
  data->skipFirstFloquetMult=0;
  data->ds0=data->rap->ds;

/*  allocate_global_memory(iap);*/

  if(global_rotations.nrtn!=NULL)free(global_rotations.nrtn);
  global_rotations.nrtn = malloc(sizeof(integer)*(iap->ndim));
  global_rotations.irtn=0;
  for(i=0;i<iap->ndim;i++)global_rotations.nrtn[i]=0;

#ifndef MFNOSAFETYNET
  if(global_rotations.nrtn==NULL)
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Out of memory trying to allocate %d bytes.",sizeof(integer)*(iap->nbc));
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

/*  end of allocate_global_memory();*/

  Auto=MFIMFCreateBaseClass(n,nalc,"AUTOBV",e);

  MFIMFSetSpace(Auto,space,e);
  MFIMFSetData(Auto,(void*)data,e);
  MFIMFSetFreeData(Auto,MFFreeDataAUTO,e);
  MFIMFSetProject(Auto,MFProjectAUTO,e);
  MFIMFSetTangentWithGuess(Auto,MFTangentAUTOOriginal,e);  /* **** */
  MFIMFSetTangentWithGuess(Auto,MFTangentAUTOWithGuess,e);
  MFIMFSetScale(Auto,MFScaleAUTO,e);
  MFIMFSetStop(Auto,MFAUTOStop,e);
  MFIMFSetProjectForSave(Auto,MFAUTOProjectToSave,e);
  MFIMFSetProjectForDraw(Auto,MFAUTOProjectToDraw,e);
  MFIMFSetProjectForBB(Auto,MFAUTOProjectForBB,e);
  MFIMFSetRMin(Auto,rap->dsmin,e);
  MFIMFSetVectorFactory(Auto,MFAUTOVectorFactory,e);
  MFIMFSetMatrixFactory(Auto,MFAUTOMatrixFactory,e);

  free(list);

#else
  Auto=NULL;
  sprintf(MFAUTOMFErrorHandlerMsg,"%s called, but AUTO was not found when multifario was installed.",RoutineName);
  MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
#endif

  return Auto;
 }

int MFAUTOProjectToSave(MFNVector u, double *x, void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOProjectToSave"};

#ifndef HAVE_AUTO
  sprintf(MFAUTOMFErrorHandlerMsg,"%s called, but AUTO was not found when multifario was installed.",RoutineName);
  MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
#endif

  return MFAUTOProjectToDraw(u,x,d,e);
 }

int MFAUTOProjectForBB(MFNVector u, double *x, void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOProjectForBB"};

#ifndef HAVE_AUTO
  sprintf(MFAUTOMFErrorHandlerMsg,"%s called, but AUTO was not found when multifario was installed.",RoutineName);
  MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
#endif

  return MFAUTOProjectToDraw(u,x,d,e);
 }

int MFProjectAUTO(int n,int k,MFNVector u0,MFNKMatrix Phi,MFNVector u,void *d,int *stab, MFErrorHandler e)
 {
  static char RoutineName[]={"MFProjectAUTO"};

#ifdef HAVE_AUTO
  int i,j,l;
  long iuz;
  double *tm=NULL;
  double *dtm=NULL;
  double **ups=NULL;
  double **uoldps=NULL;
  double **udotps=NULL;
  double **upoldp=NULL;
  double **dups=NULL;
  double *rldot=NULL;
  double *rlcur=NULL;
  double *rlold=NULL;
  double **fa=NULL;
  double *fc=NULL;
  double *par=NULL;
  MFNVector du=NULL;
  double *rds=NULL;
  struct MFAUTOData *data;
  long ndim,ntst,ncol,npar,nfpr,ntot;
  long ndxloc;
  doublereal **p0=NULL;
  doublereal **p1=NULL;
  doublereal *thu=NULL;
  doublereal *thl=NULL;
  double nuz,uzbv,spbv,bpbv,lpbv;
  long chng;
  doublecomplex *ev=NULL;
  int ii,jj;
  double t;
  int ifst;
  double umx,dumx,rdumx;
  int nrow,done;
  double au,adu;
  double rdrl,adrl;
  double delref,delmax;

  int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("Enter %s\n\n",RoutineName);fflush(stdout);}
#endif

  data=(struct MFAUTOData*)d;

  ((user_function_list*)&user)->func=data->func;
  ((user_function_list*)&user)->stpnt=NULL;
  ((user_function_list*)&user)->bcnd=data->bcnd;
  ((user_function_list*)&user)->icnd=data->icnd;
  ((user_function_list*)&user)->fopt=NULL;
  ((user_function_list*)&user)->pvls=data->pvls;

  ntst=MFAUTOBVNVGetNtst(u0,e);
  ncol=MFAUTOBVNVGetNcol(u0,e);
  ndim=MFAUTOBVNVGetNdim(u0,e);
  npar=MFAUTOBVNVGetNpar(u0,e);
  nfpr = data->iap->nfpr;
  ndxloc=ntst+1;
  ntot = data->iap->ntot;

  MFAUTOBVNVCopyDataValues(u0,u,e);

  p0=MFAUTOBVNVGetP0(u,e);
  p1=MFAUTOBVNVGetP1(u,e);
  ev=MFAUTOBVNVGetEV(u,e);

  setrtn(ndim,ntst,MFAUTOBVNVGetU(u0,e),MFAUTOBVNVGetPar(u0,e));

  dups=data->dups;                 /* Temporary */
  upoldp=data->upoldp;             /* Temporary (used to store f(ups)) */

  rds=data->rds;

  fa=data->fa;
  fc=data->fc;

  thl=MFAUTONSpaceGetThl(data->space,e);
  thu=MFAUTONSpaceGetThu(data->space,e);

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("ups (predicted solution)\n");
    MFPrintAUTOBVNVectorFull(stdout,u0,e);printf("\n");fflush(stdout);
    printf("\n");
   }
#endif

/* Copy of Predicted solution to hold new solution (ups,rlcur) */

  ups=MFAUTOBVNVGetU(u,e);
/*
 *for(i=0;i<ndxloc;i++)
 *  for(j=0;j<ndim*ncol;j++)
 *   ups[i][j]=(MFAUTOBVNVGetU(u0,e))[i][j];
 */

  rlcur=data->rlcur;
  for(i=0;i<nfpr;i++)
   {
    MFAUTOBVNVGetPar(u,e)[data->icp[i]]=MFAUTOBVNVGetPar(u0,e)[data->icp[i]];
    rlcur[i]=MFAUTOBVNVGetPar(u0,e)[data->icp[i]];
   }

/*  The "previous" solution is (uoldps,rlold). Since rds=0 for us, this is the same as the predicted solution u0    */

  uoldps=MFAUTOBVNVGetU(u0,e);
  rlold=data->rlold;
  for(i=0;i<nfpr;i++)rlold[i]=MFAUTOBVNVGetPar(u0,e)[data->icp[i]];

/* Tangent at Previous: (udotps,rldot) */

  udotps=data->udotps;
  rldot=data->rldot;
  for(i=0;i<k;i++)
   {
    du=MFMColumn(Phi,i,e);

#ifdef MFALLOWVERBOSE
  if(verbose)
     {
      printf("tangent %d)\n",i);
      MFPrintAUTOBVNVectorFull(stdout,du,e);printf("\n");fflush(stdout);
      printf("\n");fflush(stdout);
     }
#endif

    for(j=0;j<ntst+1;j++)
      for(l=0;l<ndim*ncol;l++)
        udotps[j+i*(ntst+1)][l]=(MFAUTOBVNVGetU(du,e))[j][l];

    for(j=0;j<nfpr;j++)rldot[j+NPARX*i]=(MFAUTOBVNVGetPar(du,e))[data->icp[j]];
    MFFreeNVector(du,e);
   }

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    for(i=0;i<k;i++)
     {
      du=MFMColumn(Phi,i,e);
      printf("tangent #%d (udotps))\n",i);MFPrintAUTOBVNVectorFull(stdout,du,e);printf("\n");printf("\n");fflush(stdout);
      MFFreeNVector(du,e);
     }
   }
#endif

/* Adapt the mesh to the previous solution. */

  data->iap->itp = 0;

  tm=MFAUTOBVNVGetT(u,e);
  dtm=MFAUTOBVNVGetDt(u,e);
  par=MFAUTOBVNVGetPar(u,e);
  
  for(i=0;i<k;i++)rds[i]=0.;

/*
 *   Solve the boundary value problem
 */

  (data->iap)->istop=0;

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
     MFNVector diff;
     MFKVector s;
 
     diff=MFCloneNVector(u,e);
     s=MFCreateKVector(k,e);

/* Check: is Phi^T.(u-u0) = 0?*/

     MFNSpaceDirection(data->space,u,u0,diff,e);
     MFMVMulT(data->space,Phi,diff,s,e);

     printf("Initial: the projection of u-u0 onto Phi is ");MFPrintKVector(stdout,s,e);printf(" (should be zero!)\n");fflush(stdout);

     MFFreeNVector(diff,e);
     MFFreeKVector(s,e);
    }
#endif

  nrow=ndim*ncol;
  for(data->iap->nit=1;data->iap->nit<=data->iap->itnw;(data->iap->nit)++)
   {
    ifst=0;
    if(data->iap->nit<=data->iap->nwtn)ifst=1;

/*
 *   Stores U-prime (derivative with respect to T) in (upoldp,rldot).
 */
    stupbv(data->iap,data->rap,par,data->icp,data->funi,rlcur,rlold,rldot,ndxloc,ups,uoldps,upoldp);

/*
 *   (dups,rlcur-rlold) is the difference 
 */

    for(i=0;i<ntst;i++)
      for(j=0;j<ndim*ncol;j++)
        dups[i][j]=ups[i][j]-uoldps[i][j];
    for(j=0;j<ndim;j++)
      dups[ntst][j]=ups[ntst][j]-uoldps[ntst][j];

#ifdef MFALLOWVERBOSE
    if(verbose)
     {
      MFNVector diff;
      MFKVector s;
      double delta;
 
      diff=MFCloneNVector(u0,e);
      s=MFCreateKVector(k,e);

      MFNSpaceDirection(data->space,u,u0,diff,e);
      printf("Project u-u0 onto Phi (left is tangent, right is (u-u0) Line %d, File %s\n",__LINE__,__FILE__);
      MFMVMulT(data->space,Phi,diff,s,e);
      printf("\n------The projection of u-u0 onto Phi is ");MFPrintKVector(stdout,s,e);printf(" (should be zero!)\n");fflush(stdout);

/* Check: is Phi^T.(u-u0) = 0?*/

      for(i=0;i<ntst;i++)
        for(j=0;j<ndim*ncol;j++)
          MFAUTOBVNVGetU(diff,e)[i][j]=dups[i][j];
      for(j=0;j<ndim;j++)
        MFAUTOBVNVGetU(diff,e)[ntst][j]=dups[ntst][j];
      for(j=0;j<nfpr;j++)(MFAUTOBVNVGetPar(diff,e))[data->icp[j]]=rlcur[j]-rlold[j];
      printf("rlcur[0]=%lf, rlold[0]=%lf\n",rlcur[0],rlold[0]);fflush(stdout);

      printf("Find the norm of u-u0\n");
      delta=MFNSpaceInner(data->space,diff,diff,e);
      printf("\n------The norm of u-u0 is %lf\n");fflush(stdout);
 
      MFFreeNVector(diff,e);
      MFFreeKVector(s,e);
     }
#endif

    solvbv(ifst,data->iap,data->rap,par,data->icp,funi,bcni,icni,rds,0,rlcur,rlold,rldot,ndxloc,ups,dups,uoldps,udotps,upoldp,dtm,fa,fc,p0,p1,thl,thu,1,NULL);

    dumx=0.;
    umx=0.;
    for(j=0;j<ntst;++j)
     {
      for(i=0;i<nrow;++i)
       {
	adu=fabs(fa[j][i]);
	if(adu>dumx)dumx = adu;
	au=fabs(ups[j][i]);
	if(au>umx)umx=au;
	ups[j][i]+=fa[j][i];
       }
     }

    for(i=0;i<ndim;++i)ups[ntst][i]+=fc[i];

    for(i=0;i<nfpr;++i)rlcur[i]+=fc[ndim+i];
    for(i=0;i<nfpr;++i)(MFAUTOBVNVGetPar(u,e))[data->icp[i]]=rlcur[i];

    done=1;
    rdrl=0.;
    for(i=0;i<nfpr;++i)
     {
      adrl=fabs(fc[ndim+i])/(fabs(rlcur[i])+1.);
      if(adrl>data->rap->epsl)done=0;
      if(adrl>rdrl)rdrl=adrl;
     }
    rdumx=dumx/(umx+1.);

    {
     MFNVector diff;

     diff=MFCloneNVector(u0,e);
     MFNSpaceDirection(data->space,u,u0,diff,e);
     for(i=0;i<k;i++)
      {
       du=MFMColumn(Phi,i,e);
       MFFreeNVector(du,e);
      }
    }

    if(done&&rdumx<data->rap->epsu)goto endloop;
    if(data->iap->nit==1)
     {
      delref=max(rdrl,rdumx)*20;
     }else{
      delmax=max(rdrl,rdumx);
      if(delmax>delref)
       {
        (data->iap)->istop=1;
        goto endloop;
       }
     }
   }

endloop:
  for(i=0;i<nfpr;++i)(MFAUTOBVNVGetPar(u,e))[data->icp[i]]=rlcur[i];

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
   MFNVector diff;
   MFKVector s;
   double delta;
 
   diff=MFCloneNVector(u0,e);
   s=MFCreateKVector(k,e);

/* Check: is Phi^T.(u-u0) = 0?*/

   MFNSpaceDirection(data->space,u,u0,diff,e);
   MFMVMulT(data->space,Phi,diff,s,e);
   printf("Find delta Line %d, File %s\n",__LINE__,__FILE__);
   delta=MFNSpaceInner(data->space,diff,diff,e);

   printf("Final: the projection of u-u0 onto Phi is ");MFPrintKVector(stdout,s,e);printf(" (should be zero!)\n");fflush(stdout);

   MFFreeNVector(diff,e);
   MFFreeKVector(s,e);

    printf("%s, %d iterations, istop=%d\n",RoutineName,data->iap->nit,data->iap->istop);fflush(stdout);
    printf("Solution\n");fflush(stdout);
    MFPrintAUTOBVNVectorFull(stdout,u,e);
    printf("\n\n\n");fflush(stdout);
   }
#endif

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    for(i=0;i<k;i++)
     {
      du=MFMColumn(Phi,i,e);
      printf("udotps[%d] (tangent of chart)\n",i);
      MFPrintAUTOBVNVectorFull(stdout,du,e);printf("\n");fflush(stdout);
      MFFreeNVector(du,e);
      printf("\n");
     }
   }
#endif

  MFAUTOBVNVSetNit(u,(data->iap)->nit,e);
  MFAUTOBVNVSetR(u,MFAUTOBVNVGetR(u0,e),e);

  if((data->iap)->istop==1)return 0;

  for(i=0;i<data->nUserZeros;i++)
    MFAUTOBVNVSetUZBV(u,i,par[data->userZeroParm[i]]-data->userZeroValue[i],e);
  for(i=0;i<data->nStops;i++)(*data->setStopData)(u,e);

  *stab=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("done %s, rc=%d\n",RoutineName,(data->iap)->istop);fflush(stdout);}
#endif

  return 1;
#else
  sprintf(MFAUTOMFErrorHandlerMsg,"%s called, but AUTO was not found when multifario was installed.",RoutineName);
  MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
  return 0;
#endif
 }

int MFTangentAUTOWithGuess(int n,int k,MFNVector u,MFNKMatrix Phi0,MFNKMatrix Phi,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFTangentAUTOWithGuess"};

#ifdef HAVE_AUTO

  long nllv;
  double *par;
  double *tm;
  double *dtm;
  double **ups;
  double **uoldps;
  double **udotps;
  double **upoldp;
  double **dups;
  double *rldot;
  double *rlcur;
  double *rlold;
  double **fa;
  double *fc;
  MFNVector c0;
  MFNVector c;
  double *rds;
  struct MFAUTOData *data;
  long ndim,ntst,ncol,npar;
  long ntot,ndxloc;
  doublereal **p0;
  doublereal **p1;
  doublereal *thu;
  doublereal *thl;
  int i,j,l;
  double t;
  long ifst=1;
  double *rhs;
  int ii,jj;
  int nbc,nint,nalc,nfpr;
  int verbose=0;
  int test=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s\n\n",RoutineName);fflush(stdout);}
#endif

  data=(struct MFAUTOData*)d;

  ((user_function_list*)&user)->func=data->func;
  ((user_function_list*)&user)->stpnt=NULL;
  ((user_function_list*)&user)->bcnd=data->bcnd;
  ((user_function_list*)&user)->icnd=data->icnd;
  ((user_function_list*)&user)->fopt=NULL;
  ((user_function_list*)&user)->pvls=data->pvls;

  tm=MFAUTOBVNVGetT(u,e);
  dtm=MFAUTOBVNVGetDt(u,e);

  ntst=MFAUTOBVNVGetNtst(u,e);
  ncol=MFAUTOBVNVGetNcol(u,e);
  ndim=MFAUTOBVNVGetNdim(u,e);
  npar=MFAUTOBVNVGetNpar(u,e);
  nfpr=MFAUTOBVNVGetNfpr(u,e);
  nbc =data->iap->nbc;
  nint=data->iap->nint;
  nalc=data->iap->nalc;
  ndxloc=ntst+1;
  ntot=(ntst+1)*ndim*ncol;

/* Current */

  upoldp=data->upoldp;    /* Temporary */
  dups=data->dups;        /* Temporary */

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    for(i=0;i<k;i++)
     {
      if(verbose)
       {
        c0=MFMColumn(Phi0,i,e);
        printf("udotps[%d] (reference tangent)\n",i);
        MFPrintAUTOBVNVectorFull(stdout,c0,e);
        printf("\n\n");fflush(stdout);
        MFFreeNVector(c0,e);
       }
     }
    printf("chart center:\n");
    MFPrintAUTOBVNVectorFull(stdout,u,e);
    printf("\n\n");fflush(stdout);
   }
#endif

/* This loop finds each tangent */

  for(i=0;i<k;i++)
   {

/* Previous solution */

    uoldps=data->uoldps;
    rlold=data->rlold;
    for(j=0;j<ntst+1;j++)
     for(l=0;l<ndim*ncol;l++)
       uoldps[j][l]=MFAUTOBVNVGetU(u,e)[j][l];
    for(ii=0;ii<nfpr;ii++)rlold[ii]=MFAUTOBVNVGetPar(u,e)[data->icp[ii]];

/* Tangents */

    udotps=data->udotps;
    rldot=data->rldot;
    for(ii=0;ii<k;ii++)
     {
      jj=i+ii+1;
      if(jj>k-1)jj=jj-k;

      c0=MFMColumn(Phi0,jj,e);
  
      for(j=0;j<ntst+1;j++)
       {
        for(l=0;l<ndim*ncol;l++)
         {
          if(ii<k)udotps[j+ii*(ntst+1)][l]=(MFAUTOBVNVGetU(c0,e))[j][l];
           else udotps[j+ii*(ntst+1)][l]=0.;
         }
       }

      for(j=0;j<nfpr;j++)
       {
        if(ii<k)
          rldot[j+NPARX*ii]=(MFAUTOBVNVGetPar(c0,e))[data->icp[j]];
         else
          rldot[j+NPARX*ii]=0.;
       }

      MFFreeNVector(c0,e);
     }

/* stepsize (not needed as useDefaultRHS=false) */

    rds=data->rds;
    for(ii=0;ii<k;ii++)rds[i]=0.;

/* where the result is returned */

    fa=data->fa;
    fc=data->fc;

    for(j=0;j<ntst+1;j++)
     for(l=0;l<ndim*ncol;l++)
       fa[j][l]=0.;
    for(j=0;j<nbc+nint+nfpr;j++)fc[j]=0.;

/* The previous solution */

    rlcur=data->rlcur;
    ups=MFAUTOBVNVGetU(u,e);
    par=MFAUTOBVNVGetPar(u,e);
    for(ii=0;ii<k;ii++)rlcur[ii]=par[data->icp[ii]];

    thl=MFAUTONSpaceGetThl(data->space,e);
    thu=MFAUTONSpaceGetThu(data->space,e);

    nllv=0;
    ifst=1;
    rhs=data->rhs;

    p0=MFAUTOBVNVGetP0(u,e);
    p1=MFAUTOBVNVGetP1(u,e);
/*
 *   Stores U-prime (derivative with respect to T) in (upoldp,rldot).
 */
    stupbv(data->iap,data->rap,par,data->icp,data->funi,rlcur,rlold,rldot,ndxloc,ups,uoldps,upoldp);

/*
 *   (dups,rlcur-rlold) is the difference 
 */
    c=MFMColumn(Phi,i,e);

    for(l=0;l<ntst;l++)
      for(j=0;j<ndim*ncol;j++)
        dups[l][j]=ups[l][j]-uoldps[l][j];
    for(j=0;j<ndim;j++)
      dups[ntst][j]=ups[ntst][j]-uoldps[ntst][j];

    for(j=0;j<k;j++)rhs[j]=0.;
    rhs[k-1]=1.;

#ifdef MFALLOWVERBOSE
    if(verbose)
     {
      printf("rhs=\n");fflush(stdout);
      for(j=0;j<nalc;j++)printf("     %d %lf\n",j,rhs[j]);
      printf("\n");fflush(stdout);
      fflush(stdout);
     }
#endif
    for(j=0;j<k;j++){rlcur[j]=par[data->icp[j]];rlold[j]=rlcur[j];}

    solvbv(ifst,data->iap,data->rap,par,data->icp,data->funi,data->bcni,data->icni,rds,nllv,rlcur,rlold,rldot,ndxloc,ups,dups,uoldps,udotps,upoldp,dtm,fa,fc,p0,p1,thl,thu,0,rhs);

    for(ii=0;ii<ntst;ii++)
     for(jj=0;jj<ndim*ncol;jj++)MFAUTOBVNVGetU(c,e)[ii][jj]=fa[ii][jj];
    for(j=0;j<ndim;j++)MFAUTOBVNVGetU(c,e)[ntst][j]=fc[j];
    for(j=0;j<k;j++)   MFAUTOBVNVGetPar(c,e)[data->icp[j]]=fc[j+ndim];

/*  testTangent(MFAUTOBVNVGetU(c,e),MFAUTOBVNVGetPar(c,e),data->iap,data->rap,par,data->icp,data->funi,data->bcni,data->icni,rds,nllv,rlcur,rlold,rldot,ndxloc,ups,dups,uoldps,udotps,upoldp,dtm,fa,fc,p0,p1,thl,thu);
 */

#ifdef MFALLOWVERBOSE
    if(verbose)
     {
      printf("\n%s, OUTPUT from solvbv for tangent %d\n",RoutineName,i);
      MFPrintAUTOBVNVectorFull(stdout,c,e);
      printf("\n\n");fflush(stdout);
      printf("fa=\n");fflush(stdout);
      for(j=0;j<ntst;++j)
        for(jj=0;jj<ncol;++jj)
         {
          printf(" %d %d",j,jj);
          for(ii=0;ii<ndim;++ii){printf("  %14.7lf",fa[j][ii+ndim*jj]);fflush(stdout);}
          printf("\n");fflush(stdout);
         }
      printf("\n");fflush(stdout);
      printf("fc=\n");fflush(stdout);
      for(j=0;j<nbc+nint+nalc;j++)printf("     %d %14.7lf\n",j,fc[j]);
      printf("\n");fflush(stdout);
      fflush(stdout);
     }
#endif

    MFFreeNVector(c,e);
    ifst=1;
/*  ifst=0;*/
   }

  MFGramSchmidt(data->space,Phi,e);

#ifdef MFALLOWVERBOSE
  if(verbose)
  {
   for(i=0;i<k;i++)
    {
     c=MFMColumn(Phi,i,e);
     printf("tangent[%d]\n",i);
     MFPrintAUTOBVNVectorFull(stdout,c,e);
     printf("\n\n");fflush(stdout);
     MFFreeNVector(c,e);
    }
   printf("done %s\n",RoutineName);fflush(stdout);
  }
#endif

 {
  doublereal *tm;
  doublereal *dtm;
  doublereal **up;
  doublereal ***PhiU;
  doublereal **PhiP;
  MFNVector phi;

  tm=MFAUTOBVNVGetT(u,e);
  dtm=MFAUTOBVNVGetDt(u,e);
  up=MFAUTOBVNVGetU(u,e);
  PhiU=(doublereal***)malloc(((data->iap)->nalc)*sizeof(doublereal**));
  PhiP=(doublereal**)malloc(((data->iap)->nalc)*sizeof(doublereal*));
  for(i=0;i<(data->iap)->nalc;i++)
   {
    phi=MFMColumn(Phi,i,e);
    PhiU[i]=MFAUTOBVNVGetU(phi,e);
    PhiP[i]=MFAUTOBVNVGetPar(phi,e);
    MFFreeNVector(phi,e);
   }

#ifdef MFALLOWVERBOSE
  if(test)
   {
    integer ncol,ntst,ndim;
    int j,idim,itst,icol;
    double e,error;
    doublereal *f,*df,*dfdp;
    doublereal *du;
    doublereal *dp;
    doublereal *t;
    integer it,tt,i0;
    doublereal **wt;
    doublereal **wp;
    integer jcol,jpar,jdim;

    ncol=(data->iap)->ncol;
    ntst=(data->iap)->ntst;
    ndim=(data->iap)->ndim;

    du=(doublereal*)malloc(ndim*(ncol*ntst+1)*sizeof(doublereal));
    dp=(doublereal*)malloc(k*sizeof(doublereal));
    t =(doublereal*)malloc((ncol*ntst+1)*sizeof(doublereal));

    f=(doublereal*)malloc(ndim*sizeof(doublereal));
    df=(doublereal*)malloc(ndim*ndim*sizeof(doublereal));
    dfdp=(doublereal*)malloc(ndim*k*sizeof(doublereal));

    wt=(doublereal**)malloc((ncol+1)*sizeof(doublereal*));
    wp=(doublereal**)malloc((ncol+1)*sizeof(doublereal*));
    for(icol=0;icol<ncol+1;icol++)
     {
      wt[icol]=(doublereal*)malloc(ncol*sizeof(doublereal));
      wp[icol]=(doublereal*)malloc(ncol*sizeof(doublereal));
     }

    genwts(ncol, ncol + 1, wt, wp);

    for(i=0;i<k;i++)
     {

      for(itst=0;itst<ntst;itst++)
       {
        for(icol=0;icol<ncol;icol++)
         {
          for(idim=0;idim<ndim;idim++)
           {
            du[idim+ndim*(icol+ncol*itst)]=PhiU[i][itst][idim+ndim*icol];
           }
          t[icol+ncol*itst]=tm[itst]+icol*(tm[itst+1]-tm[itst])/ncol;
         }
       }
      for(idim=0;idim<ndim;idim++)
        du[idim+ndim*ncol*ntst]=PhiU[i][ntst][idim];
      t[+ncol*ntst]=tm[ntst];

      for(j=0;j<k;j++)dp[j]=PhiP[i][data->icp[j]];

      printf("t: %lf",t[0]);
      for(j=1;j<ncol*ntst+1;j++)printf(",%lf",t[j]);
      printf("\n");fflush(stdout);

      printf("du: %lf",du[0]);
      for(j=1;j<ndim*(ncol*ntst+1);j++)printf(",%lf",du[j]);
      printf("\n");fflush(stdout);

      error=0.;
      for(itst=0;itst<ntst;itst++)
       {
        for(icol=1;icol<ncol;icol++)
         {
          for(idim=1;idim<ndim;idim++)
           {
            it=idim+ndim*(icol+ncol*itst);
            i0=     ndim*(icol+ncol*itst);
 
            e=du[it]-du[i0];
            printf("  e=(%10.7lf)-(%10.7lf)                                      s:%d, j:%d, i: %d\n",du[it],du[i0],idim,itst,icol);fflush(stdout);
            for(jcol=0;jcol<ncol;jcol++)
             {
              (data->funi)(data->iap,data->rap,ndim,ups[itst]+ndim*jcol,uoldps[itst]+ndim*jcol,data->icp,par,1,f,df,dfdp);
              for(jdim=0;jdim<ndim;jdim++)
               {
                e-=(t[icol+ncol*itst]-t[ncol*itst])*wt[icol][jcol]*df[idim+ndim*jdim]*du[jdim+ndim*(jcol+ncol*itst)];
                printf("                    -(%10.7lf)*(%10.7lf)*(%10.7lf)*(%10.7lf)\n",t[icol+ncol*itst]-t[ncol*itst],wt[icol][jcol],df[idim+ndim*jdim],du[jdim+ndim*(jcol+ncol*itst)]);fflush(stdout);
               }
              for(jpar=0;jpar<k;jpar++)
               {
                e-=(t[icol+ncol*itst]-t[ncol*itst])*wt[icol][jcol]*dfdp[idim+ndim*jdim]*dp[jpar];
                printf("                    -(%10.7lf)*(%10.7lf)*(%10.7lf)*(%10.7lf)",t[icol+ncol*itst]-t[ncol*itst],wt[icol][jcol],dfdp[jpar],dp[jpar]);fflush(stdout);
                if(jpar<k-1)printf("\n");
               }
             }
            printf("=%10.7lf\n",e);fflush(stdout);
            error+=e*e;
           }
         }
        error+=e*e;
       }
      error=sqrt(error)/ndim/(ncol*ntst+1);
      printf("Approximate error in tangent %d = %14.7le\n",i,error);fflush(stdout);
     }

    free(du);
    free(dp);

    for(icol=0;icol<ncol+1;icol++)
     {
      free(wt[icol]);
      free(wp[icol]);
     }
    free(wt);
    free(wp);

    free(f);
    free(df);
    free(dfdp);
   }
#endif


  adaptK(data->iap, MFAUTOBVNVGetNtst(u,e),  MFAUTOBVNVGetNcol(u,e), MFAUTOBVNVGetNtst(u,e), MFAUTOBVNVGetNcol(u,e), tm, dtm,  MFAUTOBVNVGetNtst(u,e)+1, up, PhiU);

  MFAUTOBVNVSetT(u,tm,e);
  for(i=0;i<MFAUTOBVNVGetNfpr(u,e);i++)
   {
    phi=MFMColumn(Phi,i,e);
    MFAUTOBVNVSetT(phi,tm,e);
    MFFreeNVector(phi,e);
   }

  free(PhiU);
  free(PhiP);
 }

  return 1;

#else
  sprintf(MFAUTOMFErrorHandlerMsg,"%s called, but AUTO was not found when multifario was installed.",RoutineName);
  MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
#endif

  return 0;
 }

double MFScaleAUTO(int n,int k,MFNVector u,MFNKMatrix Phi,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFScaleAUTO"};
  struct MFAUTOData *data;
  double R;
  int nit;
  double *tm;
  double *dtm;
  double **up;
  double ***Phip;
  MFNVector phi;
  int i;
  int verbose=0;

  data=(struct MFAUTOData*)d;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s\n",RoutineName);fflush(stdout);}
#endif

  nit=MFAUTOBVNVGetNit(u,e);
  if(nit<0)
   {
    R=data->rap->ds;
   }else if(nit<100)
   {
    R=MFAUTOBVNVGetR(u,e);

#ifdef MFALLOWVERBOSE
    if(verbose){printf("   original R=%le, nit=%d\n",R,nit);fflush(stdout);}
#endif

  /* Adapt the stepsize along the branch. */

    data->iap->nit=nit;
    if(data->iap->iads!=0)adptds(data->iap, data->rap, &R);

/*  if(MFAUTOBVNVGetNit(u,e)>3)R=.5*R;
    if(MFAUTOBVNVGetNit(u,e)<3)R=1.8*R;*/
   }else{
    R=data->rap->dsmin;
   }

#ifdef MFALLOWVERBOSE
  if(verbose){printf("   new R=%lf\n",R);fflush(stdout);}
#endif

#ifdef MFALLOWVERBOSE
  if(verbose&&R<data->rap->dsmin){printf("   R<dsmin=%lf\n",data->rap->dsmin);fflush(stdout);}
  if(verbose&&R>data->rap->dsmax){printf("   R>dsmax=%lf\n",data->rap->dsmax);fflush(stdout);}
#endif

  if(R>data->rap->dsmax)R=data->rap->dsmax;
  if(R<data->rap->dsmin)R=data->rap->dsmin;

  MFAUTOBVNVSetR(u,R);

  return R;
 }

int MFAUTOStop(MFImplicitMF this,MFNVector u0, MFNKMatrix Phi0 ,MFNVector u1, MFNKMatrix mPhi1,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOStop"};
  int i;
  int result;
  int nuz;
  struct MFAUTOData *data;
  double uzbv0,uzbv1;
  int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s\n",RoutineName);fflush(stdout);}
#endif

  data=(struct MFAUTOData*)d;

/* Set AUTO to this problem */

  ((user_function_list*)&user)->func=data->func;
  ((user_function_list*)&user)->stpnt=NULL;
  ((user_function_list*)&user)->bcnd=data->bcnd;
  ((user_function_list*)&user)->icnd=data->icnd;
  ((user_function_list*)&user)->fopt=NULL;
  ((user_function_list*)&user)->pvls=data->pvls;

  result=FALSE;
  for(i=0;i<data->nStops;i++)
   {
    result=result||(*data->testStopData)(this,u0,Phi0,u1,mPhi1,e);
    if(result)
     {
      MFAUTOBVNVSetR(u1,data->ds0);
      return;
     }
   }

/* User Zeroes */

#ifdef MFALLOWVERBOSE
  if(verbose){printf("  nUserZeros=%d\n",data->nUserZeros);fflush(stdout);}
#endif

  for(i=0;i<data->nUserZeros;i++)
   {
    uzbv0=MFAUTOBVNVGetUZBV(u0,i,e);
    uzbv1=MFAUTOBVNVGetUZBV(u1,i,e);

#ifdef MFALLOWVERBOSE
    if(verbose){printf("  User Zero %d, u0=%lf, u1=%lf\n",i,uzbv0,uzbv1);fflush(stdout);}
#endif

    result=result||uzbv0*uzbv1<0.;
    if(result)
     {
      MFAUTOBVNVSetType(u1,"UZ",e);
      MFAUTOBVNVSetR(u1,data->ds0,e);
      return;
     }
   }

#ifdef MFALLOWVERBOSE
  if(verbose){printf("done %s, result=%d\n",RoutineName,result);fflush(stdout);}
#endif

  return result;
 }


int MFAUTOProjectToDraw(MFNVector u, double *x, void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOProjectToDraw"};
  struct MFAUTOData *data;
  double *dtm;
  double **p0;
  double **p1;
  double **ups;
  double *par;
  double *thl;
  long ndim,ntst,ndxloc;
  int n1,n2,nt;
  int i;

  data=(struct MFAUTOData*)d;

  if(x==NULL)return data->iap->nfpr+1;

  ((user_function_list*)&user)->func=data->func;
  ((user_function_list*)&user)->stpnt=NULL;
  ((user_function_list*)&user)->bcnd=data->bcnd;
  ((user_function_list*)&user)->icnd=data->icnd;
  ((user_function_list*)&user)->fopt=NULL;
  ((user_function_list*)&user)->pvls=data->pvls;

  dtm=MFAUTOBVNVGetDt(u,e);
  p0=MFAUTOBVNVGetP0(u,e);
  p1=MFAUTOBVNVGetP1(u,e);
  par=MFAUTOBVNVGetPar(u,e);
  ups=MFAUTOBVNVGetU(u,e);
  ntst=MFAUTOBVNVGetNtst(u,e);
  ndim=MFAUTOBVNVGetNdim(u,e);
  ndxloc=ntst+1;

#ifdef MFNOCONFIDENCE
  if(data->pvli==NULL)
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"pvli is NULL.");
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    return;
   }
#endif

  (*(data->pvli))(data->iap, data->rap, data->icp, MFAUTOBVNVGetDt(u,e), ndxloc, ups, &ndim, p0,p1,par);

  for(i=0;i<data->iap->nfpr;i++)x[i]=par[data->icp[i]];
  x[data->iap->nfpr]=sqrt(rinpr(data->iap, ndim, ups, ups, dtm, MFAUTONSpaceGetThu(data->space,e)));

  return 0;
 }

iap_type *MFAUTOIMFGetIAP(MFImplicitMF M, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOIMFGetIAP(MFImplicitMF"};
  struct MFAUTOData *data;

#ifdef MFNOCONFIDENCE
  if(strcmp(MFImplicitMFId(M,e),"AUTOBV"))
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Manifold is not an AUTO problem.");
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    return;
   }
#endif

  data=(struct MFAUTOData*)MFIMFGetData(M,e);

  return data->iap;
 }

rap_type *MFAUTOIMFGetRAP(MFImplicitMF M, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOIMFGetRAP(MFImplicitMF"};
  struct MFAUTOData *data;

#ifdef MFNOCONFIDENCE
  if(strcmp(MFImplicitMFId(M,e),"AUTOBV"))
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Manifold is not an AUTO problem.");
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    return;
   }
#endif

  data=(struct MFAUTOData*)MFIMFGetData(M,e);

  return data->rap;
 }

integer *MFAUTOIMFGetICP(MFImplicitMF M, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOIMFGetICP(MFImplicitMF"};
  struct MFAUTOData *data;

#ifdef MFNOCONFIDENCE
  if(strcmp(MFImplicitMFId(M,e),"AUTOBV"))
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Manifold is not an AUTO problem.");
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    return;
   }
#endif

  data=(struct MFAUTOData*)MFIMFGetData(M,e);

  return data->icp;
 }

integer *MFAUTOIMFGetIUZ(MFImplicitMF M, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOIMFGetIUZ(MFImplicitMF"};
  struct MFAUTOData *data;

#ifdef MFNOCONFIDENCE
  if(strcmp(MFImplicitMFId(M,e),"AUTOBV"))
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Manifold is not an AUTO problem.");
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    return;
   }
#endif

  data=(struct MFAUTOData*)MFIMFGetData(M,e);

  return data->userZeroParm;
 }

doublereal *MFAUTOIMFGetVUZ(MFImplicitMF M, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOIMFGetVUZ(MFImplicitMF"};
  struct MFAUTOData *data;

#ifdef MFNOCONFIDENCE
  if(strcmp(MFImplicitMFId(M,e),"AUTOBV"))
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Manifold is not an AUTO problem.");
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    return;
   }
#endif

  data=(struct MFAUTOData*)MFIMFGetData(M,e);

  return data->userZeroValue;
 }

MFfunc_type MFAUTOIMFGetF(MFImplicitMF M, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOIMFGetF(MFImplicitMF"};
  struct MFAUTOData *data;

#ifdef MFNOCONFIDENCE
  if(strcmp(MFImplicitMFId(M,e),"AUTOBV"))
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Manifold is not an AUTO problem.");
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    return;
   }
#endif

  data=(struct MFAUTOData*)MFIMFGetData(M,e);

  return data->func;
 }

MFbcnd_type MFAUTOIMFGetBC(MFImplicitMF M, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOIMFGetBC(MFImplicitMF"};
  struct MFAUTOData *data;

#ifdef MFNOCONFIDENCE
  if(strcmp(MFImplicitMFId(M,e),"AUTOBV"))
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Manifold is not an AUTO problem.");
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    return;
   }
#endif

  data=(struct MFAUTOData*)MFIMFGetData(M,e);

  return data->bcnd;
 }

MFicnd_type MFAUTOIMFGetIC(MFImplicitMF M, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOIMFGetIC(MFImplicitMF"};
  struct MFAUTOData *data;

#ifdef MFNOCONFIDENCE
  if(strcmp(MFImplicitMFId(M,e),"AUTOBV"))
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Manifold is not an AUTO problem.");
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    return;
   }
#endif

  data=(struct MFAUTOData*)MFIMFGetData(M,e);

  return data->icnd;
 }

MFpvls_type MFAUTOIMFGetPV(MFImplicitMF M, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOIMFGetPV(MFImplicitMF"};
  struct MFAUTOData *data;

#ifdef MFNOCONFIDENCE
  if(strcmp(MFImplicitMFId(M,e),"AUTOBV"))
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Manifold is not an AUTO problem.");
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    return;
   }
#endif

  data=(struct MFAUTOData*)MFIMFGetData(M,e);

  return data->pvls;
 }

/* --------------------------------------------------------------------------------------------------*/

void MFFreeDataAUTO(void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFFreeDataAUTO"};
  struct MFAUTOData *data;

  data=(struct MFAUTOData*)d;

  MFFreeNSpace(data->space,e);
  free(data->rap);
  free(data->iap);
  free(data->icp);
  free_dmatrix(data->ups);
  free_dmatrix(data->uoldps);
  free_dmatrix(data->dups);
  free_dmatrix(data->upoldp);
  free_dmatrix(data->udotps);
  free(data->rds);
  free_dmatrix(data->fa);
  free(data->fc);
  free(data->rldot);
  free(data->rlcur);
  free(data->rlold);
  free(data->rhs);
  if(data->setStopData!=NULL)free(data->setStopData);
  if(data->testStopData!=NULL)free(data->testStopData);
  if(data->userZeroParm!=NULL)free(data->userZeroParm);
  if(data->userZeroValue!=NULL)free(data->userZeroValue);
  free(data);
 }

/*!    \fn MFAUTOTPBVP MFCreateAUTOTPBVP(integer k,integer ndim,MFfunc_type func,integer jac,
 *                            integer nbc,MFbcnd_type bcnd,
 *                            integer nic,MFicnd_type icnd,
 *                            integer npar,integer nicp, integer *icp,integer ntst,integer ncol,
 *                            MFpvls_type pvls, MFErrorHandler e);
 *     \brief Creates a representation of the type of two point boundary value problem used by AUTO.
 *
 *     \param k The dimension of the solution manifold.
 *     \param ndim The dimension of the phase space.
 *     \param func A routine to evaluate a flow field defined over the phase space.
 *     \param jac TRUE if the derivatives are also available through calls to func. 
 *     \param nbc The number of boundary conditions applied to a trajectory in phase space.
 *     \param bcnd A routine to evaluate the residual of the boundary conditions.
 *     \param nic The number of integral conditions applied to a trajectory in phase space.
 *     \param icnd A routine to evaluate the kernal of the integral conditions.
 *     \param npar The total number of parameters on which the flow field depends.
 *     \param nicp The number of parameters which AUTO can use to follow paths of special solutions (e.g. bifurcations).
 *     \param icp an array of length at least nicp which lists the available parameters.
 *     \param ntst The number of mesh intervals to use to discretize the trajectory.
 *     \param ncol The number of collocation points to use on wach mesh interval.
 *     \param pvls OK, I haven't figured out what this one does. I suspect I don't have to require it for BVPs
 *     \param e    An error handler.
 *     \returns A representation of the boundary value problem.
 */
MFAUTOTPBVP MFCreateAUTOTPBVP(integer k,integer ndim,MFfunc_type func,integer jac,
                              integer nbc,MFbcnd_type bcnd,
                              integer nic,MFicnd_type icnd,
                              integer npar,integer nicp, integer *icp,integer ntst,integer ncol,
                              MFpvls_type pvls, MFErrorHandler e)
 {
  static char RoutineName[]={"MFCreateAUTOTPBVP"};
  MFAUTOTPBVP tpbvp;
  int i;

  tpbvp=malloc(sizeof(struct MFAUTOTPBVPSt));

#ifndef MFNOSAFETYNET
  if(tpbvp==NULL)
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Out of memory trying to allocate %d bytes.",sizeof(struct MFAUTOTPBVPSt));
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  tpbvp->k=k;
  tpbvp->ndim=ndim;
  tpbvp->nbc=nbc;
  tpbvp->nic=nic;
  tpbvp->npar=npar;
  tpbvp->nfpr=nbc+nic-ndim+k;
  tpbvp->jac=jac;
  tpbvp->ntst=ntst;
  tpbvp->ncol=ncol;
  tpbvp->nicp=nicp;
  tpbvp->icp=icp;
  for(i=0;i<npar;i++)tpbvp->icp[i]=icp[i];

  tpbvp->func=func;
  tpbvp->bcnd=bcnd;
  tpbvp->icnd=icnd;
  tpbvp->pvls=pvls;

  tpbvp->nrefs=1;

  return tpbvp;
 }

/*! \fn void MFWriteAUTOTPBVP(FILE* fid,MFAUTOTPBVP tpbvp, MFErrorHandler e);
 *  \brief Writes a boundary value problem to a file.
 *
 *  \param fid The file to write to.
 *  \param e     An error handler.
 *  \param tpbvp The boundary value problem being queried.
 */
void MFWriteAUTOTPBVP(FILE* fid,MFAUTOTPBVP tpbvp, MFErrorHandler e);

/*! \fn MFAUTOTPBVP MFReadAUTOTPBVP(FILE* fid, MFErrorHandler e);
 *  \brief Reads a boundary value problem from a file.
 *
 *  \param fid The file to write to.
 *  \param e     An error handler.
 *  \returns tpbvp The boundary value problem.
 */
MFAUTOTPBVP MFReadAUTOTPBVP(FILE* fid, MFErrorHandler e);

/*! \fn void MFRefAUTOTPBVP(MFAUTOTPBVP tpbvp, MFErrorHandler e);
 *  \brief Adds a reference to a boundary value problem.
 *
 *  \param tpbvp The boundary value problem being referenced.
 *  \param e     An error handler.
 *  \sa ReferenceCounting MFFreeAUTOTPBVP
 */
void MFRefAUTOTPBVP(MFAUTOTPBVP tpbvp, MFErrorHandler e)
 {
  static char RoutineName[]={"MFRefAUTOTPBVP"};
  tpbvp->nrefs++;
 }

/*! \fn void MFFreeAUTOTPBVP(MFAUTOTPBVP tpbvp, MFErrorHandler e);
 *  \brief Frees a reference to a boundary value problem, and deletes the chart if there are no references left.
 *
 *  \param tpbvp The boundary value problem being unreferenced.
 *  \param e     An error handler.
 *  \sa ReferenceCounting MFRefAUTOTPBVP
 */
void MFFreeAUTOTPBVP(MFAUTOTPBVP tpbvp, MFErrorHandler e)
 {
  static char RoutineName[]={"MFFreeAUTOTPBVP"};
  tpbvp->nrefs--;
  if(tpbvp->nrefs<=0)
   {
    free(tpbvp);
   }

  return;
 }

integer MFAUTOTPBVPGetJAC(MFAUTOTPBVP tpbvp, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOTPBVPGetJAC"};
  return tpbvp->jac;
 }

integer MFAUTOTPBVPGetNBC(MFAUTOTPBVP tpbvp, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOTPBVPGetNBC"};
  return tpbvp->nbc;
 }

integer MFAUTOTPBVPGetNIC(MFAUTOTPBVP tpbvp, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOTPBVPGetNIC"};
  return tpbvp->nic;
 }

integer MFAUTOTPBVPGetNDIM(MFAUTOTPBVP tpbvp, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOTPBVPGetNDIM"};
  return tpbvp->ndim;
 }

integer MFAUTOTPBVPGetNTST(MFAUTOTPBVP tpbvp, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOTPBVPGetNTST"};
  return tpbvp->ntst;
 }

integer MFAUTOTPBVPGetNCOL(MFAUTOTPBVP tpbvp, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOTPBVPGetNCOL"};
  return tpbvp->ncol;
 }

integer MFAUTOTPBVPGetNFPR(MFAUTOTPBVP tpbvp, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOTPBVPGetNFPR"};
  struct MFAUTOData *data;

  return tpbvp->nfpr;
 }

integer MFAUTOTPBVPGetNPAR(MFAUTOTPBVP tpbvp, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOTPBVPGetNPAR"};
  return tpbvp->npar;
 }

integer *MFAUTOTPBVPGetICP(MFAUTOTPBVP tpbvp, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOTPBVPGetICP"};
  return tpbvp->icp;
 }

MFfunc_type MFAUTOTPBVPGetFUNC(MFAUTOTPBVP tpbvp, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOTPBVPGetFUNC"};
  return tpbvp->func;
 }

MFbcnd_type MFAUTOTPBVPGetBCND(MFAUTOTPBVP tpbvp, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOTPBVPGetBCND"};
  return tpbvp->bcnd;
 }

MFicnd_type MFAUTOTPBVPGetICND(MFAUTOTPBVP tpbvp, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOTPBVPGetICND"};
  return tpbvp->icnd;
 }

MFpvls_type MFAUTOTPBVPGetPVLS(MFAUTOTPBVP tpbvp, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOTPBVPGetPVLS"};
  return tpbvp->pvls;
 }

integer MFAUTOTPBVPGetK(MFAUTOTPBVP tpbvp, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOTPBVPGetK"};
  return tpbvp->k;
 }

/*!    \fn int MFAUTOGetStartPoint(MFImplicitMF M,MFAUTOTPBVP tpbvp, MFstpnt_type stpnt,doublereal *p0,MFNVector *u0,MFNKMatrix *Phi0, MFErrorHandler e);
 *     \brief Creates a starting point for AUTO.
 *
 *     \param M The solution manifold.
 *     \param tpbvp The boundary value problem.
 *     \param stpnt The routine which evalutes the starting point at each point in time [0,1].
 *     \param p0 An array containing the parameter values of the starting point.
 *     \param u0 (Output) A place to put the start point.
 *     \param Phi0 (Output) A place to put the starting tangent, or NULL if the tangent is not wanted.
 *  \param e     An error handler.
 *     \returns TRUE if sucessful.
 */
int MFAUTOGetStartPoint(MFImplicitMF M,MFAUTOTPBVP tpbvp, MFstpnt_type stpnt,doublereal *p0,MFNVector *u0,MFNKMatrix *Phi0, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOGetStartPoint"};
  integer ndim,ntst,npar,nicp,nfpr;
  integer ndxloc;
  doublereal **ups;
  doublereal *par;
  doublecomplex *ev;
  integer ncol;
  integer *icp;
  struct MFAUTOData *data;
  doublereal *thl;
  doublereal *thu;

  doublereal *rlcur;
  doublereal *rldot;
  doublereal **udotps;
  doublereal **upoldp;
  integer nodir;
  integer rc;
  int verbose=0;

  int i,j,n,k;
  MFNVector u1;
  MFNKMatrix Phi1;
  MFNVector col;

#ifdef MFNOCONFIDENCE
  if(strcmp(MFImplicitMFId(M,e),"AUTOBV"))
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Parameter 1 must be an AUTOBV");
    MFSetError(e,4,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    return 0;
   }
#endif

  data=(struct MFAUTOData*)MFIMFGetData(M,e);

  ndim=(data->iap)->ndim;
  ntst=(data->iap)->ntst;
  ncol=(data->iap)->ncol;
  nfpr=(data->iap)->nfpr;
  npar=MFIMF_K(M,e);
  nicp=(data->iap)->nicp;
  icp=data->icp;

  ndxloc=ntst+1;

  thl=MFAUTONSpaceGetThl(data->space,e);
  thu=MFAUTONSpaceGetThu(data->space,e);

  u1=MFIMFVectorFactory(M,e);
  MFAUTOBVNVSetR(u1,data->ds0,e);
  MFAUTOBVNVSetNit(u1,3,e);

  *Phi0=NULL;
  ups=MFAUTOBVNVGetU(u1,e);
  upoldp=data->upoldp;

#ifndef MFNOSAFETYNET
  if(upoldp==NULL)
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Out of memory trying to allocate %d bytes.",(ntst+1)*ndim*ncol*sizeof(double));
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0;
   }
#endif

  udotps=data->udotps;

#ifndef MFNOSAFETYNET
  if(udotps==NULL)
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Out of memory trying to allocate %d bytes.",(ntst+1)*ndim*ncol*sizeof(double));
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0;
   }
#endif

  par=MFAUTOBVNVGetPar(u1,e);
  for(i=0;i<npar;i++)par[i]=p0[i];

  rlcur=data->rlcur;

#ifndef MFNOSAFETYNET
  if(rlcur==NULL)
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Out of memory trying to allocate %d bytes.",NPARX*sizeof(doublereal));
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0;
   }
#endif

  rldot=data->rldot;

#ifndef MFNOSAFETYNET
  if(rldot==NULL)
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Out of memory trying to allocate %d bytes.",NPARX*sizeof(doublereal));
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0;
   }
#endif

  for(i=0;i<nfpr;i++)
   {
    rlcur[i]=par[icp[i]];
    rldot[i]=0.;
   }

  ev=malloc(ndim*sizeof(doublecomplex));

#ifndef MFNOSAFETYNET
  if(ev==NULL)
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Out of memory trying to allocate %d bytes.",ndim*sizeof(doublereal));
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0;
   }
#endif

  for(i=0;i<ndim;i++){ev[i].r=0.;ev[i].i=0.;}
  MFAUTOBVNVSetEV(u1,ev,e);

  ((user_function_list*)&user)->func=tpbvp->func;
  ((user_function_list*)&user)->stpnt=stpnt;
  ((user_function_list*)&user)->bcnd=tpbvp->bcnd;
  ((user_function_list*)&user)->icnd=tpbvp->icnd;
  ((user_function_list*)&user)->fopt=NULL;
  ((user_function_list*)&user)->pvls=tpbvp->pvls;

  (*data->stpnt)(data->iap,data->rap,par,data->icp,&ntst,&ncol,rlcur,rldot,ndxloc,ups,udotps,upoldp,MFAUTOBVNVGetT(u1,e),MFAUTOBVNVGetDt(u1,e),&nodir,thl,thu);
  MFAUTOBVNVSetType(u1,"EP",e);

  n=MFIMF_N(M,e);
  k=MFIMF_K(M,e);
  *u0=MFCloneNVector(u1,e);
  MFAUTOBVNVSetR(*u0,data->ds0,e);
  MFAUTOBVNVSetNit(*u0,3,e);
  *Phi0=MFIMFMatrixFactory(M,e);

  rc=MFIMFProject(M,u1,*Phi0,*u0,e);

#ifndef MFNOSAFETYNET
  if(!rc)
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Initial guess from stpnt failed to converge to the manifold");
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    goto freeAndReturn;
   }
#endif

  Phi1=MFIMFTangentSpaceWithGuess(M,*u0,*Phi0,e);
  MFFreeNKMatrix(*Phi0,e);
  *Phi0=Phi1;

#ifdef MFALLOWVERBOSE
  if(verbose)
  {
   for(i=0;i<k;i++)
    {
     col=MFMColumn(*Phi0,i,e);
     printf("Tangent[%d]\n",i);
     MFPrintAUTOBVNVectorFull(stdout,col,e);
     printf("\n\n");fflush(stdout);
     MFFreeNVector(col,e);
    }
  }
#endif

#ifndef MFNOSAFETYNET
  if(!rc)
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"The projected initial guess from stpnt failed to get tangent space");
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    goto freeAndReturn;
   }
#endif

freeAndReturn:

  free(u1);

  return rc;
 }

void MFAUTOAddStop(MFImplicitMF M, MFAUTOSetStopDataRtn setStop, MFAUTOTestStopDataRtn testStop, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOAddStop"};
  struct MFAUTOData *data;

#ifdef MFNOCONFIDENCE
  if(strcmp(MFImplicitMFId(M,e),"AUTOBV"))
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Parameter 1 must be an AUTOBV");
    MFSetError(e,4,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    return;
   }
#endif

  data=(struct MFAUTOData*)MFIMFGetData(M,e);

  if(data->nStops>=data->mStops)
   {
    data->mStops+=10;
    data->setStopData=realloc(data->setStopData,data->mStops*sizeof(MFAUTOSetStopDataRtn*));

#ifndef MFNOSAFETYNET
    if(data->setStopData==NULL)
     {
      sprintf(MFAUTOMFErrorHandlerMsg,"Out of memory trying to allocate %d bytes.",(data->mStops)*sizeof(MFAUTOSetStopDataRtn*));
      MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

    data->testStopData=realloc(data->testStopData,data->mStops*sizeof(MFAUTOTestStopDataRtn*));

#ifndef MFNOSAFETYNET
    if(data->testStopData==NULL)
     {
      sprintf(MFAUTOMFErrorHandlerMsg,"Out of memory trying to allocate %d bytes.",(data->mStops)*sizeof(MFAUTOSetStopDataRtn*));
      MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

   }
  data->setStopData[data->nStops]=setStop;
  data->testStopData[data->nStops]=testStop;
  data->nStops++;

  return;
 }

/*!    \fn void MFAUTOAddUserZero(MFImplicitMF M,int parm, double value, MFErrorHandler e);
 *     \brief Adds a linear function of a single parameter whose level set indicates a curve on the solution that
 *               needs to be resolved.
 *
 *     \param M The solution manifold.
 *     \param parm The number of the parameter (in the list of all parameters, not just the ones in icp).
 *     \param value The level at which the level set is being tagged.
 *  \param e     An error handler.
 */
void MFAUTOAddUserZero(MFImplicitMF M,int parm, double value, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOAddUserZero"};
  struct MFAUTOData *data;

#ifdef MFNOCONFIDENCE
  if(strcmp(MFImplicitMFId(M,e),"AUTOBV"))
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Parameter 1 must be an AUTOBV");
    MFSetError(e,4,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    return;
   }
#endif

  data=(struct MFAUTOData*)MFIMFGetData(M,e);

  if(data->nUserZeros>=data->mUserZeros)
   {
    data->mUserZeros+=10;
    data->userZeroParm=realloc(data->userZeroParm,data->mUserZeros*sizeof(integer*));

#ifndef MFNOSAFETYNET
    if(data->userZeroParm==NULL)
     {
      sprintf(MFAUTOMFErrorHandlerMsg,"Out of memory trying to allocate %d bytes.",(data->mUserZeros)*sizeof(integer*));
      MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

    data->userZeroValue=realloc(data->userZeroValue,data->mUserZeros*sizeof(double*));

#ifndef MFNOSAFETYNET
    if(data->userZeroValue==NULL)
     {
      sprintf(MFAUTOMFErrorHandlerMsg,"Out of memory trying to allocate %d bytes.",(data->mUserZeros)*sizeof(double*));
      MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

   }
  data->userZeroParm[data->nUserZeros]=parm;
  data->userZeroValue[data->nUserZeros]=value;
  data->nUserZeros++;
  data->iap->nuzr = data->nUserZeros;

  return;
 }

/*!    \fn void MFAUTODetectLimitPoints(MFImplicitMF M, MFErrorHandler e);
 *     \brief Indicates that limit points are to be found (places the partial derivative of the solution manifold with respect
 *             to a single parameter is zero). The default is not to locate these submanifolds.
 *
 *     \param M The solution manifold.
 *  \param e     An error handler.
 */
void MFAUTODetectLimitPoints(MFImplicitMF M, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTODetectLimitPoints"};
  struct MFAUTOData *data;
  iap_type *iap;

#ifdef MFNOCONFIDENCE
  if(strcmp(MFImplicitMFId(M,e),"AUTOBV"))
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Parameter 1 must be an AUTOBV");
    MFSetError(e,4,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    return;
   }
#endif

  data=(struct MFAUTOData*)MFIMFGetData(M,e);
  iap=data->iap;

  MFAUTOAddStop(M,MFAUTOSetLimitPointStop,MFAUTOTestLimitPointStop,e);
  iap->ilp=1;

  return;
 }

/*!    \fn void MFAUTODetectBifurcationPoints(MFImplicitMF M, MFErrorHandler e);
 *     \brief Indicates that bifurcation points are to be found (places the Jacobian of the boundary value problem has 
 *               a null space with odd dimension (usually 1). The default is not to monitor these.
 *
 *     \param M The solution manifold.
 *  \param e     An error handler.
 */
void MFAUTODetectBifurcationPoints(MFImplicitMF M, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTODetectBifurcationPoints"};

  MFAUTOAddStop(M,MFAUTOSetBifurcationPointStop,MFAUTOTestBifurcationPointStop,e);
  return;
 }

/*!    \fn void MFAUTODetectSpecialPoints(MFImplicitMF M,int skip, MFErrorHandler e);
 *     \brief Indicates that special bifurcation points are to be found (for example, Hopf bifurcations).
 *               The default is not to detect.
 *
 *     \param M The solution manifold.
 *     \param skip Indicates that skip multipliers are ignored. 
 *  \param e     An error handler.
 */
void MFAUTODetectSpecialPoints(MFImplicitMF M,int skip, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTODetectSpecialPoints"};
  struct MFAUTOData *data;

  MFAUTOAddStop(M,MFAUTOSetSpecialPointStop,MFAUTOTestSpecialPointStop,e);

  data=(struct MFAUTOData*)MFIMFGetData(M,e);
  data->skipFirstFloquetMult=skip;
  return;
 }

void MFAUTOSetLimitPointStop(MFNVector u, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOSetLimitPointStop"};
  return;
 }

void MFAUTOSetBifurcationPointStop(MFNVector u, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOSetBifurcationPointStop"};
  int i,j;
  int ndim;
  doublereal **pp;
  doublereal **p1;
  doublereal det;

/* Bifurcation Points (sign of the determinant) */

  ndim=MFAUTOBVNVGetNdim(u,e);
  p1=MFAUTOBVNVGetP1(u,e);
  pp=dmatrix(ndim,ndim);

#ifndef MFNOSAFETYNET
  if(pp==NULL)
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Out of memory trying to allocate %d bytes.",ndim*ndim*sizeof(double*));
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  for (i = 0; i < ndim; ++i) {
    for (j = 0; j < ndim; ++j) {
      pp[i][j] = p1[j][i];
    }
  }
  ge(ndim, ndim, *pp, 0, 0, NULL, 0, (doublereal*)NULL, &det);
  free_dmatrix(pp);

/* Have to worry about the sign! (see MFAlgebraic) */

  MFAUTOBVNVSetBPBV(u,det,e);

  return;
 }

void MFAUTOSetSpecialPointStop(MFNVector u, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOSetSpecialPointStop"};

  int i,j;
  int ndim;
  int iid;
  doublereal **p0; 
  doublereal **p1;
  doublereal det;
  doublecomplex *ev;
  doublereal a,b;
  doublereal amin,azm1;
  int loc;

/* Special Points (number of Floquet multipliers inside the unit circle */

  ndim=MFAUTOBVNVGetNdim(u,e);
  ev=MFAUTOBVNVGetEV(u,e);
  p0=MFAUTOBVNVGetP0(u,e);
  p1=MFAUTOBVNVGetP1(u,e);

/* this is basically fnspbv */

/*  Compute the Floquet multipliers */

  iid=0;
  flowkm(ndim, p0, p1, iid, ev);

/* Find the multiplier closest to z=1. */

  amin = 1.;
  loc=0;
  for (j = 0; j < ndim; ++j) {
    a=ev[j].r - 1.;
    b=ev[j].i;
    azm1 = sqrt(a*a+b*b);
    if (j==0 || azm1 <= amin) {
      amin = azm1;
      loc = j;
    }
  }
  if (loc != 0) {
    a=ev[loc].r;
    b=ev[loc].i;
    ev[loc].r = ev[0].r;
    ev[loc].i = ev[0].i;
    ev[0].r=a;
    ev[0].i=b;
  }

  /* Order the remaining Floquet multipliers by distance from |z|=1. */

  if(ndim>=3)
   {
    for(i=1;i<ndim-1;i++)
     {
      amin=1.;
      loc=0;
      for(j=i;j<ndim;j++)
       {
        azm1=z_abs(&ev[j])-1.;
        azm1=fabs(azm1);
        if(j==i||azm1<=amin)
         {
          amin=azm1;
          loc=j;
         }
       }
      if(loc!=i)
       {
        a=ev[loc].r;
        b=ev[loc].i;
        ev[loc].r=ev[i].r;
        ev[loc].i=ev[i].i;
        ev[i].r=a;
        ev[i].i=b;
       }
     }
   }
  return;
 }

int MFAUTOTestLimitPointStop(MFImplicitMF M,MFNVector u0,MFNKMatrix Phi0,MFNVector u1,MFNKMatrix Phi1, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOTestLimitPointStop"};
/* Have two tangents, du0=Phi0 Phi0^T(u1-u0) and du1=Phi1 Phi1^T(u1-u0)
                      du0=(dv0,dp0) and du1=(dv1,dp1)

                if du0.du1>0 and dp0.dp1<0 Limit point passed
                if du0.du1<0 and dp0.dp1>0 Limit point passed
*/
  double lpbv0,lpbv1;

  MFNVector du;
  MFNVector t0;
  MFNVector t1;
  MFKVector s0;
  MFKVector s1;
  int k;
  double dot;
  double ds;
  int result;
  int verbose=0;
  int i;
  MFNVector c;
  MFNVector c0,c1;

  k=MFIMF_K(M,e);

#ifdef MFALLOWVERBOSE
  if(0&&verbose)
   {
    printf(" %s\n",RoutineName);
    printf(" u0:\n");MFPrintAUTOBVNVectorFull(stdout,u0,e);printf("\n");fflush(stdout);
    for(i=0;i<k;i++)
     {
      c=MFMColumn(Phi0,i,e);
      printf(" Phi0[%d]:\n",i);MFPrintAUTOBVNVectorFull(stdout,c,e);printf("\n");fflush(stdout);
      MFFreeNVector(c,e);
     }
    printf(" u1:\n");MFPrintAUTOBVNVectorFull(stdout,u1,e);printf("\n");fflush(stdout);
    for(i=0;i<k;i++)
     {
      c=MFMColumn(Phi1,i,e);
      printf(" Phi1[%d]:\n",i);MFPrintAUTOBVNVectorFull(stdout,c,e);printf("\n");fflush(stdout);
      MFFreeNVector(c,e);
     }
    fflush(stdout);
   }
#endif

  if(MFNVGetIndex2(u0,e)<0)return 0;
  if(MFNVGetIndex2(u1,e)<0)return 0;

  c0=MFMColumn(Phi0,0,e);
  c1=MFMColumn(Phi1,0,e);
  printf("Find dot product of two tangent vectors Line %d, File %s\n",__LINE__,__FILE__);
  dot=MFNSpaceInner(MFIMFNSpace(M,e),c0,c1,e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("               <phi0,phi1>=%lf\n",dot);fflush(stdout);}
#endif

  dot=MFAUTOBVNVParmDot(c0,c1,e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("               r0dot*r1dot=%lf\n",dot);fflush(stdout);}
#endif

  MFFreeNVector(c0,e);
  MFFreeNVector(c1,e);

  du=MFCloneNVector(u0,e);
  s0=MFCreateKVector(k,e);
  t0=MFCloneNVector(u0,e);
  s1=MFCreateKVector(k,e);
  t1=MFCloneNVector(u0,e);
  MFNSpaceDirection(MFIMFNSpace(M,e),u1,u0,du,e);
  MFMVMulT(MFIMFNSpace(M,e),Phi0,du,s0,e);
  ds=MFKVNorm(s0,e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("        s0=");MFPrintKVector(stdout,s0);printf("  |s0|=%lf\n",ds);fflush(stdout);}
#endif

  MFKVScale(1./ds,s0,e);
  MFMVMulT(MFIMFNSpace(M,e),Phi1,du,s1,e);
  ds=MFKVNorm(s1,e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("        s1=");MFPrintKVector(stdout,s1);printf("  |s1|=%lf\n",ds);fflush(stdout);}
#endif

  MFKVScale(1./ds,s1,e);
  MFMVMul(MFIMFNSpace(M,e),Phi0,s0,t0,e);
  MFMVMul(MFIMFNSpace(M,e),Phi1,s1,t1,e);
  dot=MFAUTOBVNVParmDot(t0,t1,e);
  result=dot<0.;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("  Limit Point, r0dot*r1dot=%lf\n",dot);fflush(stdout);}
#endif

  printf("Find dot product of two tangent vectors Line %d, File %s\n",__LINE__,__FILE__);
  dot=MFNSpaceInner(MFIMFNSpace(M,e),t0,t1,e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("               <phi0,phi1>=%lf\n",dot);fflush(stdout);}
#endif

  MFFreeNVector(du,e);
  MFFreeNVector(t0,e);
  MFFreeNVector(t1,e);
  MFFreeKVector(s0,e);
  MFFreeKVector(s1,e);

  if(result)
    MFAUTOBVNVSetType(u1,"LP",e);

  return result;
 }

int MFAUTOTestBifurcationPointStop(MFImplicitMF M,MFNVector u0,MFNKMatrix Phi0,MFNVector u1,MFNKMatrix Phi1, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOTestBifurcationPointStop"};
  double bpbv0,bpbv1;
  int verbose=0;
  int det;
  int result;
  doublereal *P;
  doublereal D;
  int i,j,k;
  MFNVector coli,colj;

  bpbv0=MFAUTOBVNVGetBPBV(u0,e);
  bpbv1=MFAUTOBVNVGetBPBV(u1,e);
  k=MFIMF_K(M,e);

  P=malloc(k*k*sizeof(doublereal));

#ifndef MFNOSAFETYNET
  if(P==NULL)
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Out of memory trying to allocate %d bytes.",k*k*sizeof(doublereal));
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0;
   }
#endif

  for(i=0;i<k;i++)
   {
    coli=MFMColumn(Phi0,i,e);
    for(j=0;j<k;j++)
     {
      colj=MFMColumn(Phi1,j,e);
  printf("Find dot product of two tangent vectors Line %d, File %s\n",__LINE__,__FILE__);
      P[i+k*j]=MFNSpaceInner(MFIMFNSpace(M,e),coli,colj,e);
      MFFreeNVector(colj,e);
     }
    MFFreeNVector(coli,e);
   }

  det=1;
  D=MFAtlasDet(k,P,e);
  if(D<0)det=-1;

  result=bpbv1*bpbv0*det<0.;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("  Bifurcation, u0=%lf, u1=%lf |Phi0^TPhi_1|=%d\n",bpbv0,bpbv1,D);fflush(stdout);}
#endif

  free(P);
  if(result)MFAUTOBVNVSetType(u1,"BP",e);

  return 0;
 }

int MFAUTOTestSpecialPointStop(MFImplicitMF M,MFNVector u0,MFNKMatrix Phi0,MFNVector u1,MFNKMatrix Phi1, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOTestSpecialPointStop"};
  double spbv0,spbv1;
  int i;
  int verbose=0;
  int result;
  doublecomplex *ev0;
  doublecomplex *ev1;
  doublereal a0,b0;
  doublereal a1,b1;
  doublereal n0,n1;
  integer nins0,nins1;
  integer ndim;

  ndim=MFAUTOBVNVGetNdim(u0,e);
  ev0=MFAUTOBVNVGetEV(u0,e);
  ev1=MFAUTOBVNVGetEV(u1,e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("  Special Point, u0=%lf, u1=%lf\n",spbv0,spbv1);fflush(stdout);}
#endif

  /* Print error message if the Floquet multiplier at z=1 is inaccurate. */
  /* (ISP is set to negative and detection of bifurations is discontinued) */

  a0=ev0[0].r-1.;
  b0=ev0[0].i;
  n1=sqrt(a0*a0+b0*b0);
  if(n0>(double).05||n1>(double).05)
   {
    fprintf(stderr,"NOTE:Multiplier inaccurate\n");
    for(i=0;i<ndim;i++)
      fprintf(stderr,"        Multiplier %3li %14.6E %14.6E,  %14.6E %14.6E,\n",i,ev0[i].r,ev0[i].i,ev1[i].r,ev1[i],i);
    return 0;
   }

  if(n0<(double).01&&n1<(double).01)
   {
    fprintf(stderr,"NOTE:Multiplier accurate again\n");	
    for(i=0;i<ndim;i++)
      fprintf(stderr,"        Multiplier %3li %14.6E %14.6E,  %14.6E %14.6E,\n",i,ev0[i].r,ev0[i].i,ev1[i].r,ev1[i],i);
   }

  /* Count the number of Floquet multipliers inside the unit circle. */

  if(ndim==1)
   {
    return 0;
   }else{
    nins0=1;
    for(i=1;i<ndim;i++)
      if(z_abs(&ev0[i])<=1.)nins0++;
    nins1=1;
    for(i=1;i<ndim;i++)
      if(z_abs(&ev1[i])<=1.)nins1++;
    if(0)
      if(d_imag(&ev0[1])==0.&&ev0[1].r>0.||d_imag(&ev1[1])==0.&&ev1[1].r>0.)return 0;
   }

  return nins0!=nins1;
 }

/*!    \fn MFImplicitMF MFCreateAUTOPeriodicSolution(MFAUTOTPBVP tpbvp, MFNSpace space, MFErrorHandler e);
 *     \brief Creates an implicit representation of the solution manifold of periodic motions (closed curves in phase space).
 *
 *     \param tpbvp A two point boundary value problem of the type used by AUTO.
 *     \param space The space in which the solution manifold lives.
 *  \param e     An error handler.
 *     \returns The solution manifold. 
 */
MFImplicitMF MFCreateAUTOPeriodicSolution(MFAUTOTPBVP tpbvp, MFNSpace space, MFErrorHandler e)
 {
  static char RoutineName[]={"MFCreateAUTOPeriodicSolution"};
  MFImplicitMF Auto;
  struct MFAUTOData *data;
  int ndim,ncol,ntst,npar;
  int nbc,nint,n;
  function_list *list;
  iap_type *iap;
  rap_type *rap;
  integer iad;
  integer jac;
  integer nicp;
  integer i;
  integer *icp;
  doublereal rl0;
  doublereal rl1;
  doublereal a0;
  doublereal a1;
  integer ips;
  integer isw;
  integer irs;
  integer ilp;
  integer isp;
  integer iplt;
  integer npr;
  integer iid;
  integer nmx;
  integer mxbf;
  integer itmx;
  integer itnw;
  integer nwtn;
  integer itp;
  doublereal epsl;
  doublereal epsu;
  doublereal epss;
  doublereal ds;
  doublereal dsmin;
  doublereal dsmax;
  integer iads;
  integer nuzr;
  integer *iuz;
  doublereal *vuz;
  integer nalc;

  if(fp3==NULL)fp3=fopen("junk3.auto","rw");
  if(fp7==NULL)fp7=fopen("junk7.auto","rw");
  if(fp9==NULL)fp9=fopen("junk9.auto","rw");
  if(fp12==NULL)fp12=fopen("junk12.auto","rw");

  iap=malloc(sizeof(iap_type));

#ifndef MFNOSAFETYNET
  if(iap==NULL)
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Out of memory trying to allocate %d bytes.",sizeof(iap_type));
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  rap=malloc(sizeof(rap_type));

#ifndef MFNOSAFETYNET
  if(rap==NULL)
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Out of memory trying to allocate %d bytes.",sizeof(rap_type));
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  list=malloc(sizeof(function_list));

#ifndef MFNOSAFETYNET
  if(list==NULL)
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Out of memory trying to allocate %d bytes.",sizeof(function_list));
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

/* Read constants from file and allocate thu. */
/*  init(&iap, &rap, par, icp, thl, thu, &iuz, &vuz, &eof,k);*/

/* Set things normally read from file */

/* Problem stuff */

  ndim=MFAUTOTPBVPGetNDIM(tpbvp,e);
  nbc=MFAUTOTPBVPGetNBC(tpbvp,e);
  if(nbc!=0)abort();
  nbc=ndim;
  nint=MFAUTOTPBVPGetNIC(tpbvp,e);
  if(nint!=0)abort();
  nint=1;
  ntst=MFAUTOTPBVPGetNTST(tpbvp,e);
  ncol=MFAUTOTPBVPGetNCOL(tpbvp,e);
  npar=MFAUTOTPBVPGetNPAR(tpbvp,e);
  nalc=MFAUTOTPBVPGetK(tpbvp,e);
  iad=3;                                    /* 0=fixed mesh, >0 number of steps before adapt */
  jac=MFAUTOTPBVPGetJAC(tpbvp,e);
  nicp=tpbvp->nicp;                         /* Number of available parameters */
  icp=malloc(2*NPARX*sizeof(integer));

#ifndef MFNOSAFETYNET
  if(icp==NULL)
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Out of memory trying to allocate %d bytes.",2*NPARX*sizeof(integer));
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  for(i=0;i<nicp;i++)icp[i]=(MFAUTOTPBVPGetICP(tpbvp,e))[i];  /* List of available parameters (permutation) */
  printf("icp = [");for(i=0;i<nicp;i++){if(i>0)printf(",");printf("%d",icp[i]);}printf("]\n");fflush(stdout);


/* Omega (not used) */

  rl0=-1.e6;
  rl1=1.e6;
  a0=-1.e6;
  a1=1.e6;

/* ** Periodic Solutions. */

  ips=2;                 /* Type of Problem */
  isw=0;                 /* If restart use bifurcation point (branch switch) */
  irs=0;                 /* Restart 0=New problem, >0 label of point */
  itp=0;                 /* Restart 0=New problem, >0 label of point */

/* Detection of Bifurcations */

  ilp=0;                /* 0=no detection of limit points, >0=locat limit points */
  isp=1;                /* 0=no detection, 1,2,3 = detection (for TPBVP) */

  iplt=0;    /* For plotting projection 0=L_2 norm */
  npr=20;    /* Interval for printing out point */
  iid= 0;    /* Output level 2 prints iterations */

  nmx=100;   /* Maximum number of steps to take along a branch */
  mxbf=5;    /* Maximum number of biifurcating branches to trace out */

/* For the projection */

  itmx=8;    /* Maximum iterations to locate a bifurcation point */
  itnw=5;    /* Maximum iterations for projection */
  nwtn=3;    /* Freeze Jacobian after nwtn iterations */

/* Tolerances */

  epsl=1.e-4;  /* Relative tolerance for PAR in Newton */
  epsu=1.e-4;  /* Relative tolerance for U in Newton */
  epss=1.e-4;  /* Relative tolerance for arclength in detecting bifurcations */

  if(epsl<0.0)
   {
    printf("Warning : EPSL less then 0.0, will use absolute value instead.");
    epsl=fabs(epsl);
   }
  if(epsu<0.0)
   {
    printf("Warning : EPSU less then 0.0, will use absolute value instead.");
    epsu=fabs(epsu);
   }
  if(epss < 0.0)
   {
    printf("Warning : EPSS less then 0.0, will use absolute value instead.");
    epss=fabs(epss);
   }

/* Stepsize */

  ds=0.01;     /* Stepsize */
  dsmin=0.001; /* Minimum Stepsize */
  dsmax=2.0;   /* Maximum Stepsize */
  iads=1;      /* Adapt the stepsize after this many steps */

  if(dsmin<0.0)
   {
    printf("Warning : DSMIN less then 0.0, will use absolute value instead.");
    dsmin = fabs(dsmin);
   }

  if(dsmax<0.0)
   {
    printf("Warning : DSMAX less then 0.0, will use absolute value instead.");
    dsmax = fabs(dsmax);
   }

/* User zeroes */
/*       These should be added to the continuation method as "stop" functions */

  nuzr=0;
  iuz=NULL;
  vuz=NULL;

/* Copy them into iap and rap */

  iap->ndim = ndim;
  iap->ips = ips;
  iap->irs = irs;
  iap->ilp = ilp;
  iap->ntst = ntst;
  iap->ncol = ncol;
  iap->iad = iad;
  iap->iads = iads;
  iap->isp = isp;
  iap->isw = isw;
  iap->iplt = iplt;
  iap->nbc = nbc;
  iap->nint = nint;
  iap->nalc = nalc;
  iap->nmx = nmx;
  iap->nuzr = nuzr;
  iap->npr = npr;
  iap->mxbf = mxbf;
  iap->iid = iid;
  printf("Setting iid to %d\n",iap->iid);fflush(stdout);
  iap->itmx = itmx;
  iap->itnw = itnw;
  iap->nwtn = nwtn;
  iap->jac = jac;

  iap->ndm = ndim;
  iap->nbc0 = 1; if(nbc!=0)iap->nbc0 = nbc;
  iap->nnt0 = 1; if(nint!=0)iap->nnt0 = nint;
  iap->iuzr = 1;
  iap->itp = itp;
  iap->itpst = 0;
  iap->nfpr=nbc+nint-ndim+nalc;
  iap->ibr = 1;
  iap->nit = 0;
  iap->ntot = 0;
  iap->nins = 0;
  iap->istop = 0;
  iap->nbif = 0;
  iap->ipos = 1;
  iap->lab = 0;
  iap->nicp = nicp;

  iap->mynode = 0;
  iap->numnodes = 1;
  iap->parallel_flag = 0;

  rap->ds = ds;
  rap->dsmin = dsmin;
  rap->dsmax = dsmax;
  rap->dsold = ds;

  rap->rl0 = rl0;
  rap->rl1 = rl1;
  rap->a0 = a0;
  rap->a1 = a1;

  rap->amp = 0.;
  rap->epsl = epsl;
  rap->epsu = epsu;
  rap->epss = epss;
  rap->det = 0.;
  rap->tivp = 0.;
  rap->fldf = 0.;
  rap->hbff = 0.;
  rap->biff = 0.;
  rap->spbf = 0.;

/* End of init() */

/* Analyze iap (type of problem) and set functions in list */
/*  set_function_pointers(iap,&list); */

/*  This bit depends on the problem type . */
/*  ** Boundary value problems. */

  list->type           = AUTOBV;
  list->bvlist.funi  = fnps;
  list->bvlist.bcni  = bcps;
  list->bvlist.icni  = icps;
  list->bvlist.stpnt = stpnub;
  list->bvlist.pvli  = pvlsbv;

/*  end of set_function_pointers(); */

/*init1(&iap, &rap, icp, par);*/

  if(iap->isw==0)iap->isw=1;

  if(rap->ds==0.)rap->ds=(double).1;
  if(rap->dsmin==0.)rap->dsmin=fabs(rap->ds)*1e-4;
  rap->ds=HMACH1*rap->ds;
  rap->dsmin/=HMACH1;
  rap->dsmax=HMACH1*rap->dsmax;

/*  This bit depends on the problem type . */
/*        ** Boundary value problems */

  printf("ndim = %d\n",iap->ndim);fflush(stdout);
  printf("ntst = %d\n",iap->ntst);fflush(stdout);
  printf("ncol = %d\n",iap->ncol);fflush(stdout);
  printf("nbc = %d\n",iap->nbc);fflush(stdout);
  printf("nic = %d\n",iap->nint);fflush(stdout);
  printf("nalc = %d\n",iap->nalc);fflush(stdout);
  iap->nfpr=iap->nbc+iap->nint-iap->ndim+iap->nalc;

/* end of init1();*/

/*chdim(&iap);   This just makes sure that nfpr<NPARX */

/* These were from StartOfCNRL */

  rap->dsold=rap->ds;
  iap->isp=abs(iap->isp);
  iap->nit=0;
  iap->ntot=0;
  iap->istop=0;

/* Now have all the details of the problem */

  n=(iap->ntst+1)*iap->ndim*iap->ncol+iap->nfpr;

  nbc=iap->nbc;
  nint=iap->nint;
  ntst=iap->ntst;
  ncol=iap->ncol;
  ndim=iap->ndim;

  data=malloc(sizeof(struct MFAUTOData));

#ifndef MFNOSAFETYNET
  if(data==NULL)
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Out of memory trying to allocate %d bytes.",sizeof(struct MFAUTOData));
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  data->space=space;MFRefNSpace(space,e);
  data->iap=iap;
  data->rap=rap;
  data->icp=icp;
  data->funi=list->bvlist.funi;
  data->bcni=list->bvlist.bcni;
  data->icni=list->bvlist.icni;
  data->pvli=list->bvlist.pvli;
  data->stpnt=list->bvlist.stpnt;
  data->dups=dmatrix(ntst+1,ndim*ncol);

#ifndef MFNOSAFETYNET
  if(data->dups==NULL)
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Out of memory trying to allocate %d bytes.",(ntst+1)*ndim*ncol*sizeof(double));
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  data->upoldp=dmatrix(ntst+1,ndim*ncol);

#ifndef MFNOSAFETYNET
  if(data->upoldp==NULL)
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Out of memory trying to allocate %d bytes.",(ntst+1)*ndim*ncol*sizeof(double));
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  data->udotps=dmatrix((iap->nfpr)*(ntst+1),ndim*ncol);

#ifndef MFNOSAFETYNET
  if(data->udotps==NULL)
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Out of memory trying to allocate %d bytes.",(ntst+1)*ndim*ncol*sizeof(double));
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif


  data->rds=malloc((iap->nfpr)*sizeof(double));

#ifndef MFNOSAFETYNET
  if(data->rds==NULL)
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Out of memory trying to allocate %d bytes.",(iap->nfpr)*sizeof(double));
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  data->rldot=malloc((iap->nfpr)*NPARX*sizeof(double));

#ifndef MFNOSAFETYNET
  if(data->rldot==NULL)
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Out of memory trying to allocate %d bytes.",(iap->nfpr)*NPARX*sizeof(double));
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  data->rlcur=malloc(NPARX*sizeof(double));

#ifndef MFNOSAFETYNET
  if(data->rlcur==NULL)
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Out of memory trying to allocate %d bytes.",NPARX*sizeof(double));
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  data->rlold=malloc(NPARX*sizeof(double));

#ifndef MFNOSAFETYNET
  if(data->rlold==NULL)
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Out of memory trying to allocate %d bytes.",NPARX*sizeof(double));
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  data->rhs=malloc((iap->nfpr)*sizeof(double));

#ifndef MFNOSAFETYNET
  if(data->rhs==NULL)
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Out of memory trying to allocate %d bytes.",(iap->nfpr)*sizeof(double));
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  data->fa=dmatrix(ntst+1,ncol*ndim);

#ifndef MFNOSAFETYNET
  if(data->fa==NULL)
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Out of memory trying to allocate %d bytes.",(ntst+1)*ncol*ndim*sizeof(double));
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  data->fc=malloc((nbc+nint+iap->nalc)*sizeof(double));

#ifndef MFNOSAFETYNET
  if(data->fc==NULL)
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Out of memory trying to allocate %d bytes.",(nbc+nint+iap->nalc)*sizeof(double));
    MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  data->func=MFAUTOTPBVPGetFUNC(tpbvp,e);
  data->bcnd=MFAUTOTPBVPGetBCND(tpbvp,e);
  data->icnd=MFAUTOTPBVPGetICND(tpbvp,e);
  data->pvls=MFAUTOTPBVPGetPVLS(tpbvp,e);

  data->nStops=0;
  data->mStops=0;
  data->setStopData=NULL;
  data->testStopData=NULL;
  data->nUserZeros=0;
  data->mUserZeros=0;
  data->skipFirstFloquetMult=0;
  data->ds0=data->rap->ds;

/*  allocate_global_memory(iap);*/

  if(global_rotations.nrtn!=NULL)free(global_rotations.nrtn);
  if(iap->nbc>0)
   {
    global_rotations.nrtn = malloc(sizeof(integer)*(iap->nbc));

#ifndef MFNOSAFETYNET
    if(global_rotations.nrtn==NULL)
     {
      sprintf(MFAUTOMFErrorHandlerMsg,"Out of memory trying to allocate %d bytes.",sizeof(integer)*(iap->nbc));
      MFSetError(e,12,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return NULL;
     }
#endif

  }else
    global_rotations.nrtn = NULL;

/*  end of allocate_global_memory();*/

  Auto=MFIMFCreateBaseClass(n,iap->nalc,"AUTOBV",e);
  MFIMFSetSpace(Auto,space,e);
  MFIMFSetData(Auto,(void*)data,e);
  MFIMFSetFreeData(Auto,MFFreeDataAUTO,e);
  MFIMFSetProject(Auto,MFProjectAUTO,e);
  MFIMFSetTangentWithGuess(Auto,MFTangentAUTOWithGuess,e);
  MFIMFSetScale(Auto,MFScaleAUTO,e);
  MFIMFSetStop(Auto,MFAUTOStop,e);
  MFIMFSetProjectForSave(Auto,MFAUTOProjectToSave,e);
  MFIMFSetProjectForDraw(Auto,MFAUTOProjectToDraw,e);
  MFIMFSetProjectForBB(Auto,MFAUTOProjectForBB,e);
  MFIMFSetRMin(Auto,rap->dsmin,e);
  MFIMFSetVectorFactory(Auto,MFAUTOVectorFactory,e);
  MFIMFSetMatrixFactory(Auto,MFAUTOMatrixFactory,e);

  free(list);

  return Auto;
 }

MFNVector MFAUTOVectorFactory(MFImplicitMF this, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOVectorFactory"};
  struct MFAUTOData *data;
  MFNVector result;
  int i;
  int ntst,ncol;

  data=(struct MFAUTOData*)MFIMFGetData(this,e);
  ntst=data->iap->ntst;
  ncol=data->iap->ncol;

  if(!strcmp(MFImplicitMFId(this,e),"AUTOBV"))
   {
    result=MFCreateAUTOBVNVector(ntst,ncol,(data->iap)->ndim,(data->iap)->nicp,(data->iap)->nuzr,(data->iap)->nfpr,(data->iap)->nicp,data->icp,e);
    MFAUTOBVNVSetThu(result,MFAUTONSpaceGetThu(data->space,e),e);
   }else
    result=NULL;

  return result;
 }

MFNKMatrix MFAUTOMatrixFactory(MFImplicitMF this, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOMatrixFactory"};
  MFNKMatrix A;
  struct MFAUTOData *data;
  int i,j,k;
  MFNVector *col;
  long *icp;
  doublereal *par;

  k=MFIMF_K(this,e);
  col=malloc(k*sizeof(MFNVector));

#ifndef MFNOSAFETYNET
  if(col==NULL)
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Pointer to AUTOMF (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  data=(struct MFAUTOData*)MFIMFGetData(this,e);
  icp=data->icp;

  if(!strcmp(MFImplicitMFId(this,e),"AUTOBV"))
   {
    for(i=0;i<k;i++)
     {
      col[i]=MFAUTOVectorFactory(this,e);
      par=MFAUTOBVNVGetPar(col[i],e);
      for(j=0;j<k;j++)par[j]=0.;
      par[icp[i]]=1.;
     }
    A=MFCreateNKMatrix(k,col,e);
    for(i=0;i<k;i++)MFFreeNVector(col[i],e);
    free(col);
    return A;
   }else{
    sprintf(MFAUTOMFErrorHandlerMsg,"For AUTOMF the ImplicitMF must be an AUTOBV (argument 1)");
    MFSetError(e,4,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    return NULL;
   }
 }

/*!    \fn int MFAUTOBVSetIntegerParameter(MFImplicitMF M, char *parameterName, int value, MFErrorHandler e);
 *     \brief Allows the user to set AUTO's integer parameters. These are usually read from a file. Instead, default values
 *            are used when the AUTOBV is created, and the user may change them.
 *
 * Legal integer parameter names. See AUTO documentation for descriptions.
 *  <ul>
 *    <li>                  ips
 *    <li>                  irs
 *    <li>                  ilp
 *    <li>                  ntst
 *    <li>                  ncol
 *    <li>                  iad
 *    <li>                  iads
 *    <li>                  isp
 *    <li>                  isw
 *    <li>                  iplt
 *    <li>                  nbc
 *    <li>                  nint
 *    <li>                  nalc
 *    <li>                  nmx
 *    <li>                  nuzr
 *    <li>                  npr
 *    <li>                  mxbf
 *    <li>                  iid
 *    <li>                  itmx
 *    <li>                  itnw
 *    <li>                  nwtn
 *    <li>                  jac
 *    <li>                  iuzr
 *    <li>                  itp
 *    <li>                  itpst
 *    <li>                  ibr
 *    <li>                  nit
 *    <li>                  ntot
 *    <li>                  nins
 *    <li>                  istop
 *    <li>                  nbif
 *    <li>                  ipos
 *    <li>                  lab;
 *    <li>                  mynode
 *    <li>                  numnodes
 *    <li>                  parallel_flag
 *  </ul>
 *     \param M An AUTOBV solution manifold.
 *     \param parameterName Which parameter to set.
 *     \param value The new value.
 *  \param e     An error handler.
 *     \returns FALSE if the parameter name does not match a parameter.
 */
int MFAUTOBVSetIntegerParameter(MFImplicitMF M, char *parameterName, int value, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOBVSetIntegerParameer"};
  struct MFAUTOData *data;
  iap_type *iap;

  if(strcmp(MFImplicitMFId(M,e),"AUTOBV"))
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Parameter 1 must be an AUTOBV");
    MFSetError(e,4,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    return 0;
   }

  data=(struct MFAUTOData*)MFIMFGetData(M,e);
  iap=data->iap;

  if(!strcmp(parameterName,"ips"))iap->ips = value;
   else if(!strcmp(parameterName,"irs"))iap->irs = value;
   else if(!strcmp(parameterName,"ilp"))iap->ilp = value;
   else if(!strcmp(parameterName,"ntst"))iap->ntst = value;
   else if(!strcmp(parameterName,"ncol"))iap->ncol = value;
   else if(!strcmp(parameterName,"iad"))iap->iad = value;
   else if(!strcmp(parameterName,"iads"))iap->iads = value;
   else if(!strcmp(parameterName,"isp"))iap->isp = value;
   else if(!strcmp(parameterName,"isw"))iap->isw = value;
   else if(!strcmp(parameterName,"iplt"))iap->iplt = value;
   else if(!strcmp(parameterName,"nbc"))iap->nbc = value;
   else if(!strcmp(parameterName,"nint"))iap->nint = value;
   else if(!strcmp(parameterName,"nalc"))iap->nalc = value;
   else if(!strcmp(parameterName,"nmx"))iap->nmx = value;
   else if(!strcmp(parameterName,"nuzr"))iap->nuzr = value;
   else if(!strcmp(parameterName,"npr"))iap->npr = value;
   else if(!strcmp(parameterName,"mxbf"))iap->mxbf = value;
   else if(!strcmp(parameterName,"iid"))iap->iid = value;
   else if(!strcmp(parameterName,"itmx"))iap->itmx = value;
   else if(!strcmp(parameterName,"itnw"))iap->itnw = value;
   else if(!strcmp(parameterName,"nwtn"))iap->nwtn = value;
   else if(!strcmp(parameterName,"jac"))iap->jac = value;
   else if(!strcmp(parameterName,"iuzr"))iap->iuzr = value;
   else if(!strcmp(parameterName,"itp"))iap->itp = value;
   else if(!strcmp(parameterName,"itpst"))iap->itpst = value;
   else if(!strcmp(parameterName,"ibr"))iap->ibr = value;
   else if(!strcmp(parameterName,"nit"))iap->nit = value;
   else if(!strcmp(parameterName,"ntot"))iap->ntot = value;
   else if(!strcmp(parameterName,"nins"))iap->nins = value;
   else if(!strcmp(parameterName,"istop"))iap->istop = value;
   else if(!strcmp(parameterName,"nbif"))iap->nbif = value;
   else if(!strcmp(parameterName,"ipos"))iap->ipos = value;
   else if(!strcmp(parameterName,"lab"))iap->lab = value;
   else if(!strcmp(parameterName,"mynode"))iap->mynode = value;
   else if(!strcmp(parameterName,"numnodes"))iap->numnodes = value;
   else if(!strcmp(parameterName,"parallel_flag"))iap->parallel_flag = value;
   else return 0;

  return 1;
 }

/*!    \fn int MFAUTOBVSetRealParameter(MFImplicitMF M, char *parameterName, double value, MFErrorHandler e);
 *     \brief Allows the user to set AUTO's real valued parameters. These are usually read from a file. Instead, default values
 *            are used when the AUTOBV is created, and the user may change them.
 *
 * Legal real parameter names. See AUTO
 *    <ul>
 *     <li>                 ds
 *     <li>                 dsmin
 *     <li>                 dsmax
 *     <li>                 amp
 *     <li>                 epsl
 *     <li>                 epsu
 *     <li>                 epss
 *     <li>                 det
 *     <li>                 tivp
 *     <li>                 fldf
 *     <li>                 hbff
 *     <li>                 biff
 *     <li>                 spbf
 *    </ul>
 *     \param M An AUTOBV solution manifold.
 *     \param parameterName Which parameter to set.
 *     \param value The new value.
 *  \param e     An error handler.
 *     \returns FALSE if the parameter name does not match a parameter.
 */
int MFAUTOBVSetRealParameter(MFImplicitMF M, char *parameterName, double value, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOBVSetRealParameer"};
  struct MFAUTOData *data;
  rap_type *rap;

  if(strcmp(MFImplicitMFId(M,e),"AUTOBV"))
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Parameter 1 must be an AUTOBV");
    MFSetError(e,4,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    return 0;
   }

  data=(struct MFAUTOData*)MFIMFGetData(M,e);
  rap=data->rap;

  if(!strcmp(parameterName,"ds"))rap->ds = value;
   else if(!strcmp(parameterName,"dsmin"))rap->dsmin = value;
   else if(!strcmp(parameterName,"dsmax"))rap->dsmax = value;
   else if(!strcmp(parameterName,"amp"))rap->amp =value;
   else if(!strcmp(parameterName,"epsl"))rap->epsl =value;
   else if(!strcmp(parameterName,"epsu"))rap->epsu =value;
   else if(!strcmp(parameterName,"epss"))rap->epss =value;
   else if(!strcmp(parameterName,"det"))rap->det =value;
   else if(!strcmp(parameterName,"tivp"))rap->tivp =value;
   else if(!strcmp(parameterName,"fldf"))rap->fldf =value;
   else if(!strcmp(parameterName,"hbff"))rap->hbff =value;
   else if(!strcmp(parameterName,"biff"))rap->biff =value;
   else if(!strcmp(parameterName,"spbf"))rap->spbf =value;
   else return 0;

  return 1;
 }

/*!    \fn int MFAUTOBVGetIntegerParameter(MFImplicitMF M, char *parameterName, MFErrorHandler e);
 *     \brief Allows the user to set AUTO's integer parameters. These are usually read from a file. Instead, default values
 *            are used when the AUTOBV is created, and the user may change them.
 *
 * Legal integer parameter names. See AUTO documentation for descriptions.
 *  <ul>
 *    <li>                  ips
 *    <li>                  irs
 *    <li>                  ilp
 *    <li>                  ntst
 *    <li>                  ncol
 *    <li>                  iad
 *    <li>                  iads
 *    <li>                  isp
 *    <li>                  isw
 *    <li>                  iplt
 *    <li>                  nbc
 *    <li>                  nint
 *    <li>                  nalc
 *    <li>                  nmx
 *    <li>                  nuzr
 *    <li>                  npr
 *    <li>                  mxbf
 *    <li>                  iid
 *    <li>                  itmx
 *    <li>                  itnw
 *    <li>                  nwtn
 *    <li>                  jac
 *    <li>                  iuzr
 *    <li>                  itp
 *    <li>                  itpst
 *    <li>                  ibr
 *    <li>                  nit
 *    <li>                  ntot
 *    <li>                  nins
 *    <li>                  istop
 *    <li>                  nbif
 *    <li>                  ipos
 *    <li>                  lab;
 *    <li>                  mynode
 *    <li>                  numnodes
 *    <li>                  parallel_flag
 *  </ul>
 *     \param M An AUTOBV solution manifold.
 *     \param parameterName Which parameter value to retreive. A warning is issued if the parameter name does not match a parameter.
 *  \param e     An error handler.
 *     \returns The current value of the parameter.
 */
int MFAUTOBVGetIntegerParameter(MFImplicitMF M, char *parameterName, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOBVGetIntegerParameer"};
  struct MFAUTOData *data;
  iap_type *iap;

  if(strcmp(MFImplicitMFId(M,e),"AUTOBV"))
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Parameter 1 must be an AUTOBV");
    MFSetError(e,4,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    return 0;
   }

  data=(struct MFAUTOData*)MFIMFGetData(M,e);
  iap=data->iap;

  if(!strcmp(parameterName,"ips"))return iap->ips;
   else if(!strcmp(parameterName,"irs"))return iap->irs;
   else if(!strcmp(parameterName,"ilp"))return iap->ilp;
   else if(!strcmp(parameterName,"ntst"))return iap->ntst;
   else if(!strcmp(parameterName,"ncol"))return iap->ncol;
   else if(!strcmp(parameterName,"iad"))return iap->iad;
   else if(!strcmp(parameterName,"iads"))return iap->iads;
   else if(!strcmp(parameterName,"isp"))return iap->isp;
   else if(!strcmp(parameterName,"isw"))return iap->isw;
   else if(!strcmp(parameterName,"iplt"))return iap->iplt;
   else if(!strcmp(parameterName,"nbc"))return iap->nbc;
   else if(!strcmp(parameterName,"nint"))return iap->nint;
   else if(!strcmp(parameterName,"nalc"))return iap->nalc;
   else if(!strcmp(parameterName,"nmx"))return iap->nmx;
   else if(!strcmp(parameterName,"nuzr"))return iap->nuzr;
   else if(!strcmp(parameterName,"npr"))return iap->npr;
   else if(!strcmp(parameterName,"mxbf"))return iap->mxbf;
   else if(!strcmp(parameterName,"iid"))return iap->iid;
   else if(!strcmp(parameterName,"itmx"))return iap->itmx;
   else if(!strcmp(parameterName,"itnw"))return iap->itnw;
   else if(!strcmp(parameterName,"nwtn"))return iap->nwtn;
   else if(!strcmp(parameterName,"jac"))return iap->jac;
   else if(!strcmp(parameterName,"iuzr"))return iap->iuzr;
   else if(!strcmp(parameterName,"itp"))return iap->itp;
   else if(!strcmp(parameterName,"itpst"))return iap->itpst;
   else if(!strcmp(parameterName,"ibr"))return iap->ibr;
   else if(!strcmp(parameterName,"nit"))return iap->nit;
   else if(!strcmp(parameterName,"ntot"))return iap->ntot;
   else if(!strcmp(parameterName,"nins"))return iap->nins;
   else if(!strcmp(parameterName,"istop"))return iap->istop;
   else if(!strcmp(parameterName,"nbif"))return iap->nbif;
   else if(!strcmp(parameterName,"ipos"))return iap->ipos;
   else if(!strcmp(parameterName,"lab"))return iap->lab;
   else if(!strcmp(parameterName,"mynode"))return iap->mynode;
   else if(!strcmp(parameterName,"numnodes"))return iap->numnodes;
   else if(!strcmp(parameterName,"parallel_flag"))return iap->parallel_flag;
   else{
    sprintf(MFAUTOMFErrorHandlerMsg,"Parameter %s is not a real valued parameter",parameterName);
    MFSetError(e,4,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    return 0;
   }

  return 0;
 }

/*!    \fn double MFAUTOBVGetRealParameter(MFImplicitMF M, char *parameterName, MFErrorHandler e);
 *     \brief Allows the user to set AUTO's real valued parameters. These are usually read from a file. Instead, default values
 *            are used when the AUTOBV is created, and the user may change them.
 *
 * Legal real parameter names. See AUTO
 *    <ul>
 *     <li>                 ds
 *     <li>                 dsmin
 *     <li>                 dsmax
 *     <li>                 amp
 *     <li>                 epsl
 *     <li>                 epsu
 *     <li>                 epss
 *     <li>                 det
 *     <li>                 tivp
 *     <li>                 fldf
 *     <li>                 hbff
 *     <li>                 biff
 *     <li>                 spbf
 *    </ul>
 *     \param M An AUTOBV solution manifold.
 *     \param parameterName Which parameter value to retreive. A warning is issued if the parameter name does not match a parameter.
 *  \param e     An error handler.
 *     \returns The current value of the parameter.
 */
double MFAUTOBVGetRealParameter(MFImplicitMF M, char *parameterName, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTOBVGetRealParameer"};
  struct MFAUTOData *data;
  rap_type *rap;

  if(strcmp(MFImplicitMFId(M,e),"AUTOBV"))
   {
    sprintf(MFAUTOMFErrorHandlerMsg,"Parameter 1 must be an AUTOBV");
    MFSetError(e,4,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    return 0;
   }

  data=(struct MFAUTOData*)MFIMFGetData(M,e);
  rap=data->rap;

  if(!strcmp(parameterName,"ds"))return rap->ds;
   else if(!strcmp(parameterName,"dsmin"))return rap->dsmin;
   else if(!strcmp(parameterName,"dsmax"))return rap->dsmax;
   else if(!strcmp(parameterName,"amp"))return rap->amp;
   else if(!strcmp(parameterName,"epsl"))return rap->epsl;
   else if(!strcmp(parameterName,"epsu"))return rap->epsu;
   else if(!strcmp(parameterName,"epss"))return rap->epss;
   else if(!strcmp(parameterName,"det"))return rap->det;
   else if(!strcmp(parameterName,"tivp"))return rap->tivp;
   else if(!strcmp(parameterName,"fldf"))return rap->fldf;
   else if(!strcmp(parameterName,"hbff"))return rap->hbff;
   else if(!strcmp(parameterName,"biff"))return rap->biff;
   else if(!strcmp(parameterName,"spbf"))return rap->spbf;
   else{
    sprintf(MFAUTOMFErrorHandlerMsg,"Parameter %s is not a real valued parameter",parameterName);
    MFSetError(e,4,RoutineName,MFAUTOMFErrorHandlerMsg,__LINE__,__FILE__);
    return 0.;
   }

  return 0.;
 }

int testTangent(doublereal **udot, doublereal *rdot,iap_type *iap, rap_type *rap, doublereal *par, integer *icp,
                FUNI_TYPE((*funi)), BCNI_TYPE ((*bcni)), ICNI_TYPE((*icni)), doublereal *rds,
                integer nllv, doublereal *rlcur, doublereal *rlold, doublereal *rldot, integer ndxloc,
                doublereal **ups, doublereal **dups, doublereal **uoldps,
                doublereal **udotps, doublereal **upoldp, doublereal *dtm,doublereal **fa, doublereal *fc,
                doublereal **p0, doublereal **p1, doublereal *thl, doublereal *thu)
 {
  integer ndim,ntst,ncol,nbc,nint,nalc,nfpr;
  integer na,nra,ncb,nrc,nca;
  int i,j,k,ii,jj,kk;

  doublereal ***aa;
  doublereal ***bb;
  doublereal ***cc;
  doublereal **dd;
  doublereal *wi;
  doublereal **wp;
  doublereal **wt;

  /* System generated locals */
  integer dicd_dim1, dfdu_dim1, dfdp_dim1;
  
  /* Local variables */
  integer l, m;

  integer ib, ic;
  doublereal ddt;

  doublereal *dicd, *ficd, *dfdp, *dfdu, *uold;
  doublereal *f;
  doublereal *u, *wploc;
  doublereal *uic, *uio, *prm, *uid, *uip;

  doublereal *up;
  doublereal *up1;
  doublereal *uoldp;
  doublereal *uoldp1;

  doublereal tmp;

  double *aa_offset;
  double *dfdu_offset;
  double wt_tmp;
  int jp1;

  integer dbc_dim1;

  doublereal *dbc;
  doublereal *fbc;
  doublereal *ubc0;
  doublereal *ubc1;

  ntst=iap->ntst;
  ncol = iap->ncol;
  nbc =iap->nbc ;
  ndim=iap->ndim;
  nint=iap->nint;
  nalc=iap->nalc;
  nfpr=iap->nfpr;
  ncb =nbc+nint+nalc-ndim;
  nrc =ncb;
  nra =ndim*ncol;
  nca =ndim*(ncol+1);
  na=ntst;

  printf("ncol=%d, nbc=%d, nint=%d, nalc=%d, ncb=%d, nrc=%d, nra=%d, ntst=%d\n",ncol,nbc,nint,nalc,ncb,nrc,nra,ntst);fflush(stdout);

  wi   = (doublereal *)malloc(sizeof(doublereal)*(ncol+1) );
  wp   = dmatrix(ncol+1, ncol);
  wt   = dmatrix(ncol+1, ncol);

  if(nint>0)
   {
    dicd = (doublereal *)malloc(sizeof(doublereal)*nint*(ndim + NPARX));
    ficd = (doublereal *)malloc(sizeof(doublereal)*nint);
   }else{
    ficd = dicd = NULL;
   }

  ndxloc=iap->ntst+1;
  nbc =iap->nbc;
  ntst=na;
  nfpr =iap->nfpr;

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

  dbc  = (doublereal *)malloc(sizeof(doublereal)*(nbc)*(2*ndim + NPARX));
  fbc  = (doublereal *)malloc(sizeof(doublereal)*(nbc));
  ubc0 = (doublereal *)malloc(sizeof(doublereal)*ndim);
  ubc1 = (doublereal *)malloc(sizeof(doublereal)*ndim);

  aa=dmatrix_3d(ntst         ,ndim*ncol    ,ndim*(ncol+1));
  bb=dmatrix_3d(ntst         ,ndim*ncol    ,nbc+nint+nalc-ndim);
  cc=dmatrix_3d(ntst         ,nbc+nint+nalc,ndim*(ncol+1));
  dd=dmatrix   (nbc+nint+nalc              ,ncb);

  dicd_dim1 = nint;
  dfdu_dim1 = ndim;
  dfdp_dim1 = ndim;
  dbc_dim1 = nbc;

  /* Generate AA, BB and DD: */
  
  /* Initialize to zero. */

  for(i=0;i<nint+nalc;++i)
    for(k=0;k<ncb;++k)
      dd[i][k] = 0.;

  for(j=0;j<ntst;++j)
   {
    up = ups[j];
    up1 = ups[j + 1];
    uoldp = uoldps[j];
    uoldp1 = uoldps[j + 1];
      
    ddt = 1. / dtm[j];
    for(ic=0;ic<ncol;++ic)
     {
      for(k=0;k<ndim;++k)
       {
	u[k]=wt[ncol][ic]*up1[k];
	uold[k]=wt[ncol][ic]*uoldp1[k];
	for(l=0;l<ncol;++l)
         {
	  u[k]+=wt[l][ic]*up[k+ndim*l];
	  uold[k]+=wt[l][ic]*uoldp[k+ndim*l];
 	 }
       }

      for(i=0;i<NPARX;++i)prm[i]=par[i];

      (*(funi))(iap, rap, ndim, u, uold, icp, prm, 2, f, dfdu, dfdp);

      /* transpose dfdu for optimal access */

      for(ii=0;ii<ndim;++ii)
       {
        for(jj=0;jj<ii;++jj)
         {
          tmp = dfdu[ii + jj * ndim];
          dfdu[ii+jj*ndim]=dfdu[jj+ii*ndim];
          dfdu[jj+ii*ndim]=tmp;
         }
       }

      for(ib=0;ib<ncol+1;++ib)wploc[ib]=ddt*wp[ib][ic];

      for(i=0;i<ndim;++i)
       {
        aa_offset=aa[j][i+ndim*ic];
        dfdu_offset=&ARRAY2D(dfdu,0,i);
        for(ib=0;ib<ncol+1;++ib)
         {
          wt_tmp=-wt[ib][ic];
          for(k=0;k<ndim;++k)aa_offset[k]=wt_tmp*dfdu_offset[k];
          aa_offset[i]+=wploc[ib];
          aa_offset+=ndim;
         }
        for (k = 0; k < ncb; ++k)
         {
          bb[j][i+ndim*ic][k] = -ARRAY2D(dfdp, i, icp[k]);
         }
       }
     }
   }

  /* generate CC and DD; */

  /* Boundary conditions */

  if(nbc>0)
   {
    for(i=0;i<ndim;++i)
     {
      ubc0[i] = ups[0][i];
      ubc1[i] = ups[na][i];
     }

    (*(bcni))(iap, rap, ndim, par, icp, nbc, ubc0, ubc1, fbc, 2, dbc);
    for(i=0;i<nbc;++i)
     {
      fc[i]=0.;
      for(k=0;k<ndim;++k)
       {
        fc[i]+=ARRAY2D(dbc,i,k)*udot[0][k+ndim*0];
        fc[i]+=ARRAY2D(dbc,i,ndim+k)*udot[ntst][k+ndim*0];
       }

      for(k=0;k<ncb;++k)
       {
        fc[i]+=ARRAY2D(dbc,i,(ndim*2)+icp[k])*rdot[k];
       }
     }
   }
  
  /*     Integral constraints : */

  if(nint>0)
   {
    for(m=0;m<nint;++m)fc[nbc+m]=0.;

    for(j=0;j<ntst;++j)
     {
      for(i=0;i<ncol;++i)
       {
	for(k=0;k<ndim;++k)
         {
	  uic[k]=ups   [j][k+ndim*i];
	  uio[k]=uoldps[j][k+ndim*i];
	  uid[k]=udotps[j][k+ndim*i];
	  uip[k]=upoldp[j][k+ndim*i];
	 }
	
	(*(icni))(iap, rap, ndim, par, icp, nint, uic, uio, uid, uip, ficd, 2, dicd);
	
	for(m=0;m<nint;++m)
         {
	  for(k=0;k<ndim;++k)
           {
            fc[nbc+m]+=dtm[j]*wi[i]*ARRAY2D(dicd,m,k)*udot[j][i+ndim*k];
           }
	  for(k=0;k<nfpr;++k)
           {
            fc[nbc+m]+=dtm[j]*wi[i]*ARRAY2D(dicd,m,ndim+icp[k])*rdot[k];
           }
	 }
       }
     }
    j=ntst;
    i=0;
    for(k=0;k<ndim;++k)
     {
      uic[k]=ups   [j][k+ndim*i];
      uio[k]=uoldps[j][k+ndim*i];
      uid[k]=udotps[j][k+ndim*i];
      uip[k]=upoldp[j][k+ndim*i];
     }
	
    (*(icni))(iap, rap, ndim, par, icp, nint, uic, uio, uid, uip, ficd, 2, dicd);

    for(m=0;m<nint;++m)
     {
      for(k=0;k<ndim;++k)
       {
        fc[nbc+m]+=dtm[j]*wi[i]*ARRAY2D(dicd,m,k)*udot[j][i+ndim*k];
       }
      for(k=0;k<nfpr;++k)
       {
        fc[nbc+m]+=dtm[j]*wi[i]*ARRAY2D(dicd,m,ndim+icp[k])*rdot[k];
       }
     }
   }
	
  /*     Pseudo-arclength equation : */

  for(m=0;m<nalc;++m)
   {
    fc[nbc+nint+m]=0.;
    for(j=0;j<ntst;++j)
     {
      for(i=0;i<ncol;++i)
       {
        for(k=0;k<ndim;++k)
         {
          fc[nbc+nint+m]+=dtm[j]*thu[k]*wi[i]*udotps[j+m*(ntst+1)][k+ndim*i]*udot[j][k+ndim*i];
         }
       }
      for(k=0;k<ndim;++k)
       {
        fc[nbc+nint+m]+=dtm[ntst-1]*thu[k]*wi[ncol]*udotps[ntst+m*(ntst+1)][k]*udot[ntst][k];
       }
     }
    for(i=0;i<nfpr;++i)
      fc[nbc+nint+m]+=thl[i]*rdot[i];
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

  for(i=0;i<=ntst;i++)
   for(j=0;j<ndim*ncol;j++)fa[i][j]=0.;

  for(i=0;i<ntst;i++)
   {
    for(m=0;m<ncol*ndim;m++)
     {
      for(j=0;j<ncol*ndim;j++)
       {
        fa[i][m]+=aa[i][m][j]*udot[i][j];
       }
      for(j=0;j<ncb;j++)
        fa[i][m]+=bb[i][m][j]*rdot[j];
     }
   }

  return 0;
 }

static int MFTangentAUTOOriginal(int n,int k,MFNVector u,MFNKMatrix Phi,MFNKMatrix Phi0,void *d,MFErrorHandler e)
 {
  static char RoutineName[]={"MFTangentAUTOOriginal"};
  long nllv;
  doublereal *par;
  double *tm;
  double *dtm;
  doublereal **ups;
  doublereal **uoldps;
  doublereal **udotps;
  doublereal **upoldp;
  doublereal **dups;
  double *rldot;
  double *rlcur;
  double *rlold;
  doublereal **fa;
  double *fc;
  MFNVector du;
  double *rds;
  struct MFAUTOData *data;
  long ndim,ntst,ncol,npar;
  long ntot,ndxloc;
  doublereal **p0;
  doublereal **p1;
  doublereal *thu;
  doublereal *thl;
  int i,j,l;
  MFNVector col;
  long ifst=1;
  doublereal *rhs;
  int ii,jj;
  int nbc,nint,nalc,nfpr;

  data=(struct MFAUTOData*)d;

  par=MFAUTOBVNVGetPar(u,e);

  tm=MFAUTOBVNVGetT(u,e);
  dtm=MFAUTOBVNVGetDt(u,e);

/* Previous: */

  uoldps=MFAUTOBVNVGetU(u,e);
  rlold=(double*)malloc(data->iap->nfpr*sizeof(double));if(rlold==(double*)NULL)abort();
  for(i=0;i<data->iap->nfpr;i++)rlold[i]=par[data->icp[i]];

  udotps=data->udotps;
  rldot=data->rldot;
  rlcur=data->rlcur;

  ntst=MFAUTOBVNVGetNtst(u,e);
  ncol=MFAUTOBVNVGetNcol(u,e);
  ndim=MFAUTOBVNVGetNdim(u,e);
  npar=MFAUTOBVNVGetNpar(u,e);
  ndxloc=ntst+1;

  ntot=(ntst+1)*ndim*ncol;
  ndxloc=ntst+1;

  nbc=data->iap->nbc;
  nint=data->iap->nint;
  nalc=data->iap->nalc;
  nfpr=data->iap->nfpr;

  for(i=0;i<k;i++)
   {
    du = MFMColumn(Phi,i,e);
    for(j=0;j<ntst+1;j++)
     {
      for(l=0;l<ndim*ncol;l++)
       {
        udotps[j+i*(ntst+1)][l]=0.;
       }
     }
    for(j=0;j<npar;j++)rldot[j+npar*i]=0.;
    MFFreeNVector(du,e);
   }

/* Temporaries */

  dups=data->dups;

  upoldp=data->upoldp;

  rds=data->rds;
  for(i=0;i<k;i++)rds[i]=0.;

  fa=data->fa;
  fc=data->fc;

  for(j=0;j<ntst+1;j++)
   for(l=0;l<ndim*ncol;l++)
    {
     fa[j][l]=0.;
     dups[j][l]=0.;
     upoldp[j][l]=0.;
    }
  for(j=0;j<nfpr;j++)
   {
    fc[j]=0.;
   }

  thl=MFAUTONSpaceGetThl(data->space,e);
  thu=MFAUTONSpaceGetThu(data->space,e);

  nllv=k;
  ifst=1;
  rhs=data->rhs;
  for(i=0;i<k;i++)
   {
    for(j=0;j<k;j++)rhs[j]=0.;
    rhs[k-i-1]=1.;
    col=MFMColumn(Phi,i,e);
    MFAUTOBVNVCopyDataValues(u,col,e);
    ups=MFAUTOBVNVGetU(col,e);
    par=MFAUTOBVNVGetPar(col,e);
    p0=MFAUTOBVNVGetP0(col,e);
    p1=MFAUTOBVNVGetP1(col,e);

    solvbv( ifst,data->iap,data->rap,par,data->icp,data->funi,data->bcni,data->icni,rds, nllv,rlcur,rlold,rldot, ndxloc,ups,dups,uoldps,udotps,upoldp,dtm,fa,fc,p0,p1,thl,thu,0,rhs);

    for(ii=0;ii<ntst;ii++)
     for(jj=0;jj<ndim*ncol;jj++)ups[ii][jj]=fa[ii][jj];
    for(j=0;j<ndim;j++)ups[ntst][j]=fc[j];
    for(j=0;j<nfpr;j++)rlcur[j]=fc[ndim+j];
    for(j=0;j<nfpr;j++)par[data->icp[j]]=rlcur[j];

    for(ii=0;ii<ntst;ii++)
     for(jj=0;jj<ndim*ncol;jj++)ups[ii][jj]=0.;
    for(j=0;j<ndim;j++)ups[ntst][j]=0.;
/*  for(j=0;j<nfpr;j++)rlcur[j]=0.;
    rlcur[i]=1.;
    for(j=0;j<nfpr;j++)par[data->icp[j]]=rlcur[j];*/

    MFFreeNVector(col,e);
    ifst=0;
   }

  MFGramSchmidt(data->space,Phi,e);

  return;
 }

#ifdef __cplusplus
}
#endif
