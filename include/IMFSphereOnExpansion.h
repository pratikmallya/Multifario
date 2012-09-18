#ifndef __IMFSPHEREONEXPANSION_H__
#define __IMFSPHEREONEXPANSION_H__

#include <MFAtlas.h>
#include <IMFExpansion.h>
#include <MFErrorHandler.h>

/*! \defgroup IMFSphereOnExpansion */

#ifdef __cplusplus
 extern "C" {
#endif

/*! \addtogroup MFImplicitMF
 *  @{
 */

/*! \addtogroup IMFSphereOnExpansion
 *  @{
 */

/*! \class IMFSphereOnExpansion IMFSphereOnExpansion.h IMFSphereOnExpansion.h
 *  \brief A manifold which is the intersection of a ball about a point and an expansion about the same point.
 *         This manifold is used as a manifold of initial conditions for computing stable and unstable manifolds
 *         of hyperbolic fixed points.
 */

/*! \fn MFImplicitMF IMFCreateSphereOnExpansion(IMFExpansion E, IMFFlow F, MFKVector p0,double eps, double R, double r,MFErrorHandler e);
 *  \brief Creates a manifold which is the intersection of a sphere and the surface defined by an expansion.
 *
 *  \param E The expansion.
 *  \param F The flow.
 *  \param p0 The parameters for the flow.
 *  \param eps The tolerance allowed for the error on the conditions that the manifold lie on the expansion.
 *  \param R The raidus of the sphere.
 *  \param r The default scale associated with the new manifold (MFIMFSetR).
 *  \param e A place to return errors.
 *  \returns A new MFImplicitMF.
 */
MFImplicitMF IMFCreateSphereOnExpansion(IMFExpansion,IMFFlow,MFKVector,double,double,double,MFErrorHandler);

IMFExpansion IMFSphereOnExpansionGetLocal(MFAtlas,MFNVector,MFErrorHandler);

/*! @} */

/*! @} */

#ifdef __cplusplus
}
#endif

#endif
