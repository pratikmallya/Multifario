#ifndef __MFCONTINUATIONMETHOD__
#define __MFCONTINUATIONMETHOD__
#include <MFImplicitMF.h>
#include <MFNRegion.h>
#include <MFNVector.h>
#include <MFNKMatrix.h>
#include <MFErrorHandler.h>

/*! \defgroup MFContinuationMethod */

/*! \addtogroup MFContinuationMethod
 *  @{
 */

/*! \class MFContinuationMethod MFContinuationMethod.h MFContinuationMethod.h
 *  \brief A continuation algorithm for computing a manifold. This is a base class which is passed to the Compute and Extend
 *         Atlas routines in MFAtlas. It provides the routine that these routines actually call. In addtion each method will
 *         have it's own parameters that control the algorith.
 */
struct MFContinuationMethodST;
typedef struct MFContinuationMethodST *MFContinuationMethod;

#ifdef __cplusplus
 extern "C" {
#endif

MFContinuationMethod MFCreateContinuationMethodBaseClass(char *type,MFErrorHandler);
void MFFreeContinuationMethod(MFContinuationMethod,MFErrorHandler);
void MFRefContinuationMethod(MFContinuationMethod,MFErrorHandler);
char *MFContinuationMethodGetType(MFContinuationMethod,MFErrorHandler);
void *MFContinuationMethodGetParmBlock(MFContinuationMethod,MFErrorHandler);
void MFContinuationMethodSetParmBlock(MFContinuationMethod, void *parms,MFErrorHandler);
void MFContinuationMethodSetFreeParmBlock(MFContinuationMethod, void (*FreeParmBlock)(void*,MFErrorHandler),MFErrorHandler);
void MFContinuationMethodSetExtendAtlasMultipleWithTangents(MFContinuationMethod, void (*ExtendAtlasMultipleWithTangents)(struct MFContinuationMethodST*,MFAtlas,MFImplicitMF,MFNRegion,int,MFNVector*,MFNKMatrix*,MFErrorHandler),MFErrorHandler);
void MFContinuationMethodSetCloseAtlas(MFContinuationMethod, void (*CloseAtlas)(struct MFContinuationMethodST*,MFAtlas,MFErrorHandler),MFErrorHandler);
void MFContinuationMethodSetFlushAtlas(MFContinuationMethod, void (*FlushAtlas)(struct MFContinuationMethodST*,MFAtlas,MFErrorHandler),MFErrorHandler);

MFAtlas MFComputeAtlas(MFContinuationMethod H, MFImplicitMF M, MFNRegion Omega, MFNVector u0,MFErrorHandler);
MFAtlas MFComputeAtlasWithTangent(MFContinuationMethod H, MFImplicitMF M, MFNRegion Omega, MFNVector u0, MFNKMatrix Phi0,MFErrorHandler);
MFAtlas MFComputeAtlasMultiple(MFContinuationMethod H, MFImplicitMF M, MFNRegion Omega, int m, MFNVector *u0,MFErrorHandler);
MFAtlas MFComputeAtlasMultipleWithTangents(MFContinuationMethod H, MFImplicitMF M, MFNRegion Omega, int m, MFNVector *u0, MFNKMatrix *Phi0,MFErrorHandler);
void MFExtendAtlas(MFAtlas S, MFContinuationMethod H, MFImplicitMF M, MFNRegion Omega, MFNVector u0,MFErrorHandler);

/*! \fn void MFContinuationMethodExtendAtlasMultiple(MFContinuationMethod algorithm,MFAtlas A,MFImplicitMF M,MFNRegion Omega,int n,MFNVector *u0,MFNKMatrix *Tan0,MFErrorHandler e);
 *  \brief This is the basic continuation routine. It is called by the routines for the MFAtlas, though it can be called directly.
 *
 *  \param algorithm The Algorithm.
 *  \param A The atlas.
 *  \param M The implicitly defined manifold.
 *  \param Omega A region for the computation.
 *  \param n The number of initial points in the arrays u0 and Tan0.
 *  \param u0 A list of initial points on the manifold.
 *  \param Tan0 A list of bases for the tangent spaces of the manifold at the initial points.
 *  \param e A place to return errors.
 */
void MFExtendAtlasMultiple(MFAtlas S, MFContinuationMethod H, MFImplicitMF M, MFNRegion Omega, int m, MFNVector *u0,MFErrorHandler);

void MFExtendAtlasWithTangent(MFAtlas S, MFContinuationMethod H, MFImplicitMF M, MFNRegion Omega, MFNVector u0,MFNKMatrix Phi0,MFErrorHandler);
void MFExtendAtlasMultipleWithTangents(MFAtlas S, MFContinuationMethod H, MFImplicitMF M, MFNRegion Omega, int m, MFNVector *u0, MFNKMatrix *Phi0,MFErrorHandler);

/*! \fn void MFContinuationMethodFlushAtlas(MFContinuationMethod algorithm,MFAtlas A,MFErrorHandler e);
 *  \brief Flush an atlas. It allows the algorithm to finish up any pending steps (like writing output files), however the
 *         algorithm may be invoked again without starting from scratch.
 *
 *  \param algorithm The Algorithm.
 *  \param A The Atlas.
 *  \param e A place to return errors.
 */
void MFFlushAtlas(MFContinuationMethod H, MFAtlas S,MFErrorHandler);

/*! \fn void MFContinuationMethodCloseAtlas(MFContinuationMethod algorithm,MFAtlas A,MFErrorHandler e);
 *  \brief Closes an atlas. It allows the algorithm to finish all of its operations (like closing output files).
 *
 *  \param algorithm The Algorithm.
 *  \param A The Atlas.
 *  \param e A place to return errors.
 */
void MFCloseAtlas(MFContinuationMethod H, MFAtlas S,MFErrorHandler);

MFContinuationMethod MFCreateRheinboldtsMethod(MFErrorHandler);
int MFRheinboldtSetIntegerParameter(MFContinuationMethod,char*,int,MFErrorHandler);
int MFRheinboldtGetIntegerParameter(MFContinuationMethod,char*,MFErrorHandler);
int MFRheinboldtSetRealParameter(MFContinuationMethod,char*,double,MFErrorHandler);
double MFRheinboldtGetRealParameter(MFContinuationMethod,char*,MFErrorHandler);

MFContinuationMethod MFCreateAllgowerSchmidtsMethod(MFErrorHandler);
int MFAllgowerSchmidtSetIntegerParameter(MFContinuationMethod,char*,int,MFErrorHandler);
int MFAllgowerSchmidtGetIntegerParameter(MFContinuationMethod,char*,MFErrorHandler);
int MFAllgowerSchmidtSetRealParameter(MFContinuationMethod,char*,double,MFErrorHandler);
double MFAllgowerSchmidtGetRealParameter(MFContinuationMethod,char*,MFErrorHandler);

#ifdef __cplusplus
}
#endif

/*! @} */

#endif
