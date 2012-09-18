#ifndef __IMFINTERPOLATION_H__
#define __IMFINTERPOLATION_H__

#include <MFAtlas.h>
#include <IMFFlow.h>
#include <MFErrorHandler.h>

#ifdef __cplusplus
 extern "C" {
#endif

/*! \defgroup InvariantManifolds */
/*! \defgroup IMFInterpolation */

/*! \addtogroup InvariantManifolds
 *  @{
 */

/*! \addtogroup IMFInterpolation
 *  @{
 */

/*! \class IMFInterpolation IMFInterpolation.h IMFInterpolation.h
 *  \brief The routine which locates an interpolation point for a new fat trajectory.
 */

/*! \fn MFNVector IMFGetInterpolationPoint(MFAtlas A,IMFFlow F, MFKVector p0, MFAtlas c, double tmax, MFNRegion Omega, MFErrorHandler e)
 *  \brief Finds a point on the boundary of the current charts in the atlas at which the flow is outward.
 *
 *  \param A The atlas.
 *  \param F The flow.
 *  \param p0 The parameters for the flow.
 *  \param c The manifold of initial conditions. Points on this have first priority.
 *  \param tmax The largest time allowed on the manifold.
 *  \param Omega A region which limits the computation.
 *  \param e A place to return errors.
 *  \returns A new point satisfying the requirements.
 */
MFNVector IMFGetInterpolationPoint(MFAtlas,IMFFlow,MFKVector,MFAtlas,double,MFNRegion,MFErrorHandler);

/*! @} */

/*!  @} */

#ifdef __cplusplus
}
#endif

#endif
