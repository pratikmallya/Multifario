/*
 *  PROGRAM NAME:  multifario
 *
 *  (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
 *  CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
 *
 *  Please refer to the LICENSE file in the top directory
 *
 *      author: Mike Henderson mhender@watson.ibm.com
 *      date:   August 1, 2006
 */

static char *id="@(#) $Id: MFAUTO.c,v 1.9 2009/04/23 18:15:55 mhender Exp $";

#include <multifarioConfig.h>

#ifdef HAVE_AUTO

static char MFAUTOErrorMsg[256];

#include <MFAtlas.h>
#include <MFContinuationMethod.h>
#include <MFImplicitMF.h>
#include <MFNRegion.h>
#include <MFNSpace.h>
#include <MFNVector.h>
#include <MFPrint.h>
#include <stdio.h>
#include <MFAUTO.h>

#ifdef __cplusplus
 extern "C" {
#endif

extern FILE *fp3;
extern FILE *fp7;
extern FILE *fp9;
extern FILE *fp12;
extern int global_conpar_type;
extern int global_setubv_type;
extern int global_reduce_type;
extern int global_num_procs;
extern int global_verbose_flag;

extern user_function_list user;    /* This is where AUTO gets the things to evaluate */

extern struct {
  integer irtn;
  integer *nrtn;
} global_rotations;

void MFAtlasPageOutAllCharts(MFAtlas,MFErrorHandler);
doublereal *MFAUTONSpaceGetThl(MFNSpace,MFErrorHandler);
doublereal *MFAUTONSpaceGetThu(MFNSpace,MFErrorHandler);

int MFStpntFromAUTOBVNVector(integer,doublereal,doublereal*,doublereal*);
static MFNVector MFStpntFromAUTOBVNVectorU0;

static MFErrorHandler MFStpntFromAUTOBVNVectorError;

struct MFAUTOParmBlockST
 {
  int verbose;
  MFstpnt_type stpnt;
 };

typedef struct MFAUTOParmBlockST *MFAUTOParmBlock;

static void MFAUTORunAtlas(MFContinuationMethod CM,MFAtlas S,MFImplicitMF M, MFNRegion Omega, int m, MFNVector *u0, MFNKMatrix *Phi0,MFErrorHandler);
static void MFFreeAUTOParmBlock(void*,MFErrorHandler);

/*!  \fn MFContinuationMethod MFCreateAUTOsMethod(MFErrorHandler e);
 *   \brief Creates a continuation method which is implemented as a call to AUTO's toplevel.
 *
 *   \param A place to return errors.
 *   \returns A new continuation method.
 */
MFContinuationMethod MFCreateAUTOsMethod(MFErrorHandler e)
 {
  static char RoutineName[]={"MFCreateAUTOsMethod"};
  MFContinuationMethod result;
  MFAUTOParmBlock parms;

  parms=malloc(sizeof(struct MFAUTOParmBlockST));

#ifndef MFNOSAFETYNET
  if(parms==NULL)
   {
    sprintf(MFAUTOErrorMsg,"Out of space trying to allocate %d bytes.",sizeof(struct MFAUTOParmBlockST));
    MFSetError(e,12,RoutineName,MFAUTOErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    MFFreeContinuationMethod(result,e);
    return NULL;
   }
#endif

  result=MFCreateContinuationMethodBaseClass("AUTO",e);
  MFContinuationMethodSetParmBlock(result,parms,e);
  MFContinuationMethodSetFreeParmBlock(result,MFFreeAUTOParmBlock,e);
  MFContinuationMethodSetExtendAtlasMultipleWithTangents(result,MFAUTORunAtlas,e);

  return result;
 }

void MFFreeAUTOParmBlock(void *p, MFErrorHandler e)
 { 
  free((MFAUTOParmBlock)p);

  return;
 }

void MFAUTORunAtlas(MFContinuationMethod CM,MFAtlas S,MFImplicitMF M, MFNRegion Omega, int m, MFNVector *u0, MFNKMatrix *Phi0,MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTORunAtlas"};

  function_list *list;
  iap_type *iap;
  rap_type *rap;
  integer *icp;
  integer *icu;
  integer *iuz;
  doublereal *vuz;
  doublereal *par;
  doublereal *thl;
  doublereal *thu;
  int i;

#ifdef MFNOCONFIDENCE
  if(strcmp(MFImplicitMFId(M,e),"AUTOBV"))
   {
    sprintf(MFAUTOErrorMsg,"Manifold is not an AUTO problem.");
    MFSetError(e,12,RoutineName,MFAUTOErrorMsg,__LINE__,__FILE__);
    return;
   }

  for(i=0;i<m;i++)
   {
    if(u0[i]==NULL)
     {
      sprintf(MFAUTOErrorMsg,"Initial Guess u0[%d] is NULL.",m);
      MFSetError(e,12,RoutineName,MFAUTOErrorMsg,__LINE__,__FILE__);
      return;
     }
   }
#endif

  par=MFAUTOBVNVGetPar(u0[0],e);

  if(fp3==NULL)fp3=fopen("auto.3","w");
  if(fp7==NULL)fp7=fopen("auto.7","w");
  if(fp9==NULL)fp9=stderr;
  if(fp12==NULL)fp12=fopen("auto.12","w");

  iap=MFAUTOIMFGetIAP(M,e);
  rap=MFAUTOIMFGetRAP(M,e);
  icp=MFAUTOIMFGetICP(M,e);
  iuz=MFAUTOIMFGetIUZ(M,e);
  vuz=MFAUTOIMFGetVUZ(M,e);

  thl=MFAUTONSpaceGetThl(MFIMFNSpace(M,e),e);
  thu=MFAUTONSpaceGetThu(MFIMFNSpace(M,e),e);

  if(iap->isw==0)iap->isw=1;

  if(rap->ds==0.)rap->ds=(double).1;
  if(rap->dsmin==0.)rap->dsmin=fabs(rap->ds)*1e-4;
  rap->ds=HMACH1*rap->ds;
  rap->dsmin/=HMACH1;
  rap->dsmax=HMACH1*rap->dsmax;

  rap->dsold=rap->ds;
  iap->isp=abs(iap->isp);
  iap->nit=0;
  iap->ntot=0;
  iap->istop=0;

/* Now have all the details of the problem */

  if(global_rotations.nrtn!=NULL)free(global_rotations.nrtn);
  if(iap->ndm>0)
   {
    global_rotations.nrtn = malloc(sizeof(integer)*(iap->ndm));

#ifndef MFNOSAFETYNET
   if(global_rotations.nrtn==NULL)
    {
     sprintf(MFAUTOErrorMsg,"Out of space trying to allocate %d bytes.",sizeof(integer)*(iap->nbc));
     MFSetError(e,12,RoutineName,MFAUTOErrorMsg,__LINE__,__FILE__);
     MFErrorHandlerOutOfMemory(e);
     return;
    }
#endif

   }else
    global_rotations.nrtn = NULL;

  icu=icp; /* Selects "output parameters" */

  ((user_function_list*)&user)->func=MFAUTOIMFGetF(M,e);
  ((user_function_list*)&user)->stpnt=MFStpntFromAUTOBVNVector;
    MFStpntFromAUTOBVNVectorU0=u0[0];
     MFRefNVector(u0[0],e);
    MFStpntFromAUTOBVNVectorError=e;
     MFRefError(e);
  ((user_function_list*)&user)->bcnd=MFAUTOIMFGetBC(M,e);
  ((user_function_list*)&user)->icnd=MFAUTOIMFGetIC(M,e);
  ((user_function_list*)&user)->fopt=NULL;
  ((user_function_list*)&user)->pvls=MFAUTOIMFGetPV(M,e);

  autobv(iap, rap, par, icp, icu, funi, bcni, icni, stpnub, pvlsbv, thl, thu, iuz, vuz);

  MFFreeNVector(u0[0],e);

  return;
 }

int MFStpntFromAUTOBVNVector(integer ndim, doublereal t, doublereal *u, doublereal *par)
 {
  int ntst,npar;
  int i,it;
  double *p0;
  double **u0;
  double *t0;

  npar=MFAUTOBVNVGetNpar(MFStpntFromAUTOBVNVectorU0,MFStpntFromAUTOBVNVectorError);
  p0=MFAUTOBVNVGetPar(MFStpntFromAUTOBVNVectorU0,MFStpntFromAUTOBVNVectorError);
  for(i=0;i<npar;i++)par[i]=p0[i];

  u0=MFAUTOBVNVGetU(MFStpntFromAUTOBVNVectorU0,MFStpntFromAUTOBVNVectorError);
  t0=MFAUTOBVNVGetT(MFStpntFromAUTOBVNVectorU0,MFStpntFromAUTOBVNVectorError);

  ntst=MFAUTOBVNVGetNtst(MFStpntFromAUTOBVNVectorU0,MFStpntFromAUTOBVNVectorError);

  it=0;while(it<ntst+1 && t0[it]<=t)it++;

  if(it<ntst)
   {
    for(i=0;i<ndim;i++)u[i]=((t-t0[it])*u0[it][i]+(t0[it+1]-t)*u0[it+1][i])/(t0[it+1]-t0[it]);
   }else{
    for(i=0;i<ndim;i++)u[i]=u0[it-1][i];
   }

  return 1;
 }

#else

int MFThereIsNoAUTO_MFAUTO()
 {
  return 0;
 }

#endif

#ifdef __cplusplus
}
#endif
