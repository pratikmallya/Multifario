#ifndef __IMF_H__
#define __IMF_H__
#include <MFAtlas.h>
#include <IMFFlow.h>
#include <MFErrorHandler.h>

#ifdef __cplusplus
 extern "C" {
#endif

int  MFSolveFull(int,double*,double*,MFErrorHandler);
int  IMFSolveFull(int,double*,double*,MFErrorHandler);
void IMFPrintFull(int,double*,double*,MFErrorHandler);
void IMFPrintFullSchematic(int,double*,double*,MFErrorHandler);

MFAtlas IMFComputeStableInvariantManifold2(IMFFlow,char*,MFNVector,MFKVector,MFNRegion,double,double,double,int,int,double,MFErrorHandler);
MFAtlas IMFComputeUnstableInvariantManifold2(IMFFlow,char*,MFNVector,MFKVector,MFNRegion,double,double,double,int,int,double,MFErrorHandler);
MFAtlas IMFComputeInvariantManifold(IMFFlow,MFKVector,char*,MFAtlas,MFNRegion,double,double,double,int,int,double,MFErrorHandler);

/*! \defgroup InvariantManifolds */

/*! \addtogroup InvariantManifolds
 *  @{
 */

/*! \fn MFAtlas IMFComputeStableInvariantManifold(IMFFlow F,char *name,MFNVector u0,MFKVector p0,MFNRegion Omega,double eps,double dt,double tmax,int maxInterp,int maxCharts,double R0,MFErrorHandlere)
 * \brief Computes the stable invariant manifold of a fixed point.
 *
 * \param F The flow.
 * \param name A name to used for paging files and so forth.
 * \param u0 The fixed point.
 * \param p0 The parameters of the flow for the fixed point.
 * \param Omega A region to bound the computation.
 * \param eps The tolerance on the size of the quadratic terms over a balls. Controls the stepsize.
 * \param dt The initial timestep to use along a fat trajectory.
 * \param tmax The upper limit on trajectory length.
 * \param maxInterp The upper limit on the number of interpolation performed.
 * \param maxCharts The upper limit on the number of charts in the atlas.
 * \param R0 The radius for the initial ball about the fixedpoint that serves as initial surface for the manifold.
 * \param e A place to return errors.
 * \returns A new Atlas
 */
MFAtlas IMFComputeStableInvariantManifold(IMFFlow,char*,MFNVector,MFKVector,MFNRegion,double,double,double,int,int,double,MFErrorHandler);

/*! \fn MFAtlas IMFComputeUntableInvariantManifold(IMFFlow F,char *name,MFNVector u0,MFKVector p0,MFNRegion Omega,double eps,double dt,double tmax,int maxInterp,int maxCharts,double R0,MFErrorHandlere)
 * \brief Computes the unstable invariant manifold of a fixed point.
 *
 * \param F The flow.
 * \param name A name to used for paging files and so forth.
 * \param u0 The fixed point.
 * \param p0 The parameters of the flow for the fixed point.
 * \param Omega A region to bound the computation.
 * \param eps The tolerance on the size of the quadratic terms over a balls. Controls the stepsize.
 * \param dt The initial timestep to use along a fat trajectory.
 * \param tmax The upper limit on trajectory length.
 * \param maxInterp The upper limit on the number of interpolation performed.
 * \param maxCharts The upper limit on the number of charts in the atlas.
 * \param R0 The radius for the initial ball about the fixedpoint that serves as initial surface for the manifold.
 * \param e A place to return errors.
 * \returns A new Atlas
 */
MFAtlas IMFComputeUnstableInvariantManifold(IMFFlow,char*,MFNVector,MFKVector,MFNRegion,double,double,double,int,int,double,MFErrorHandler);

/*! @} */

#ifdef __cplusplus
 }
#endif

#endif
