/*    @(#)MFAUTO.h	1.4
      02/07/24 09:35:38
 
      Author: Mike Henderson
      Date:   July 18, 2002
*/
#ifndef __MFAUTO_H__
#define __MFAUTO_H__

#include <multifarioConfig.h>

#include <MFErrorHandler.h>
#include <MFAtlas.h>
#include <MFNRegion.h>
#include <MFNVector.h>
#include <MFNSpace.h>
#include <MFImplicitMF.h>

#include <auto_f2c.h>
#include <auto_c.h>

#ifdef __cplusplus
 extern "C" {
#endif

doublereal **MFAUTOBVNVGetU(MFNVector,MFErrorHandler);
doublereal *MFAUTOBVNVGetT(MFNVector,MFErrorHandler);
doublereal *MFAUTOBVNVGetDt(MFNVector,MFErrorHandler);
doublereal *MFAUTOBVNVGetPar(MFNVector,MFErrorHandler);
int MFAUTOBVNVGetNtst(MFNVector,MFErrorHandler);
int MFAUTOBVNVGetNcol(MFNVector,MFErrorHandler);
int MFAUTOBVNVGetNdim(MFNVector,MFErrorHandler);
int MFAUTOBVNVGetNpar(MFNVector,MFErrorHandler);
int MFAUTOBVNVGetNfpr(MFNVector,MFErrorHandler);
long *MFAUTOBVNVGetICP(MFNVector,MFErrorHandler);
doublereal MFAUTOBVNVParmDot(MFNVector,MFNVector,MFErrorHandler);

struct MFAUTOTPBVPSt;
typedef struct MFAUTOTPBVPSt *MFAUTOTPBVP;

typedef int (*MFfunc_type)(integer,const doublereal*,const integer*,const doublereal*,integer,doublereal*,doublereal*,doublereal*);
typedef int (*MFbcnd_type)(integer,const doublereal*,const integer*,integer,const doublereal*,const doublereal*,integer,doublereal*,doublereal*);
typedef int (*MFicnd_type)(integer,const doublereal*,const integer*,integer,const doublereal*,const doublereal*,const doublereal*,const doublereal*,integer,doublereal*,doublereal*);
typedef int (*MFstpnt_type)(integer,doublereal,doublereal*,doublereal*);
typedef int (*MFpvls_type)(integer,const void*,doublereal*);

MFAUTOTPBVP MFCreateAUTOTPBVP(integer,integer,MFfunc_type,integer,
                              integer,MFbcnd_type,
                              integer,MFicnd_type,
                              integer,integer,integer*,integer,integer,MFpvls_type,MFErrorHandler);


int MFAUTOSetIntegerParameter(MFAUTOTPBVP,char*,int,MFErrorHandler);
int MFAUTOSetRealParameter(MFAUTOTPBVP,char*,double,MFErrorHandler);
int MFAUTOGetIntegerParameter(MFAUTOTPBVP,char*,MFErrorHandler);
double MFAUTOGetRealParameter(MFAUTOTPBVP,char*,MFErrorHandler);
MFImplicitMF MFCreateAUTOBV(MFAUTOTPBVP,MFNSpace,MFErrorHandler);
MFImplicitMF MFCreateAUTOPeriodicOrbit(MFAUTOTPBVP,MFNSpace,MFErrorHandler);

integer MFAUTOTPBVPGetK(MFAUTOTPBVP,MFErrorHandler);
integer MFAUTOTPBVPGetNDIM(MFAUTOTPBVP,MFErrorHandler);
integer MFAUTOTPBVPGetNTST(MFAUTOTPBVP,MFErrorHandler);
integer MFAUTOTPBVPGetNCOL(MFAUTOTPBVP,MFErrorHandler);
integer MFAUTOTPBVPGetNPAR(MFAUTOTPBVP,MFErrorHandler);
integer MFAUTOTPBVPGetNFPR(MFAUTOTPBVP,MFErrorHandler);
integer *MFAUTOTPBVPGetICP(MFAUTOTPBVP,MFErrorHandler);
void MFRefAUTOTPBVP(MFAUTOTPBVP,MFErrorHandler);
void MFFreeAUTOTPBVP(MFAUTOTPBVP,MFErrorHandler);
integer MFAUTOTPBVPGetJAC(MFAUTOTPBVP,MFErrorHandler);
integer MFAUTOTPBVPGetNBC(MFAUTOTPBVP,MFErrorHandler);
integer MFAUTOTPBVPGetNIC(MFAUTOTPBVP,MFErrorHandler);

MFfunc_type MFAUTOTPBVPGetFUNC(MFAUTOTPBVP,MFErrorHandler);
MFbcnd_type MFAUTOTPBVPGetBCND(MFAUTOTPBVP,MFErrorHandler);
MFicnd_type MFAUTOTPBVPGetICND(MFAUTOTPBVP,MFErrorHandler);
MFpvls_type MFAUTOTPBVPGetPVLS(MFAUTOTPBVP,MFErrorHandler);

MFNSpace MFCreateAUTONSpace(struct MFAUTOTPBVPSt*,doublereal*,doublereal*,MFErrorHandler);
MFImplicitMF MFCreateAUTOBV(struct MFAUTOTPBVPSt*,MFNSpace,MFErrorHandler);
MFNRegion MFNRegionCreateAUTO(MFNSpace,int,long*,double*,double*,double,double,MFErrorHandler);
MFNVector MFCreateAUTOBVNVector(int,int,int,int,int,int,int,long*,MFErrorHandler);

int MFAUTOGetStartPoint(MFImplicitMF,MFAUTOTPBVP,MFstpnt_type,doublereal*,MFNVector*,MFNKMatrix*,MFErrorHandler);

typedef void (*MFAUTOSetStopDataRtn)(MFNVector,MFErrorHandler);
typedef int (*MFAUTOTestStopDataRtn)(MFImplicitMF,MFNVector,MFNKMatrix,MFNVector,MFNKMatrix,MFErrorHandler);

void MFAUTOAddStop(MFImplicitMF,MFAUTOSetStopDataRtn,MFAUTOTestStopDataRtn,MFErrorHandler);
void MFAUTOAddUserZero(MFImplicitMF,int,double,MFErrorHandler);
void MFAUTODetectLimitPoints(MFImplicitMF,MFErrorHandler);
void MFAUTODetectBifurcationPoints(MFImplicitMF,MFErrorHandler);
void MFAUTODetectSpecialPoints(MFImplicitMF,int,MFErrorHandler);

MFContinuationMethod MFCreateAUTOsMethod(MFErrorHandler);

iap_type *MFAUTOIMFGetIAP(MFImplicitMF,MFErrorHandler);
rap_type *MFAUTOIMFGetRAP(MFImplicitMF,MFErrorHandler);
integer *MFAUTOIMFGetICP(MFImplicitMF,MFErrorHandler);
integer *MFAUTOIMFGetIUZ(MFImplicitMF,MFErrorHandler);
doublereal *MFAUTOIMFGetVUZ(MFImplicitMF,MFErrorHandler);
MFfunc_type MFAUTOIMFGetF(MFImplicitMF,MFErrorHandler);
MFbcnd_type MFAUTOIMFGetBC(MFImplicitMF,MFErrorHandler);
MFicnd_type MFAUTOIMFGetIC(MFImplicitMF,MFErrorHandler);
MFpvls_type MFAUTOIMFGetPV(MFImplicitMF,MFErrorHandler);

#ifdef __cplusplus
}
#endif

#endif
