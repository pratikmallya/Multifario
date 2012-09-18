#ifndef __IMFFIXEDPT_H__
#define __IMFFIXEDPT_H__
#include <MFAtlas.h>
#include <IMFExpansion.h>
#include <IMFFlow.h>

#ifdef __cplusplus
 extern "C" {
#endif

/*! \defgroup InvariantManifolds */
/*! \defgroup IMFFixedPoint */

/*! \addtogroup InvariantManifolds
 *  @{
 */

/*! \addtogroup IMFFixedPoint
 *  @{
 */

/*! \class IMFFixedPoint , IMFFixedPt.h IMFFixedPt.h
 *  \brief A group of utilities for computing stable and unstable manifolds of a hyperbolic fixed point.
 *
 */

/*! \fn void IMFFindExpansionNearFixedPt(MFNVector U0, MFKVector P0, MFNKMatrix Du, MFNKMatrix Dv, IMFFlow F, IMFExpansion U, IMFExpansion a,MFErrorHandler e);
 *  \brief Given a fixed point in a flow and a decomposition of the phase space into two invariant linear subspaces, creates
 *            an expansion of the invariant manifold corresponding to one of the invariant subspaces.
 *
 *  \param U0 The fixed point.
 *  \param P0 The parameters of the flow of the fixed point.
 *  \param Du The first invariant linear subspace.
 *  \param Dv The second, complementary invariant linear subspace.
 *  \param F  The flow.
 *  \param U  A user provided expansion in which the expansion defining the surface is placed.
 *  \param a  A user provided expansion in which the expansion defining the complement of the surface is placed.
 *  \param e A place to return errors.
 */
void IMFFindExpansionNearFixedPt(MFNVector,MFKVector,MFNKMatrix,MFNKMatrix,IMFFlow,IMFExpansion,IMFExpansion,MFErrorHandler);

/*! \fn MFNKMatrix IMFGetBasisForUnstableInvariantSubspace(IMFFlow F, MFNVector u, MFKVector p,MFErrorHandler e);
 *  \brief Given a hyperbolic fixed point in a flow, finds a basis for the 
 *            unstable invariant linear subspace of a hyperbolic fixed point.
 *
 *  \param F The flow.
 *  \param u The fixed point. 
 *  \param p The parameters of the flow for the fixed point.
 *  \param e A place to return errors.
 *  \returns A new MFNKMatrix containing an orthonormal basis for the unstable invariant linear subspace of the hyperbolic fixed
 *             point.
 */
MFNKMatrix IMFGetBasisForUnstableInvariantSubspace(IMFFlow,MFNVector,MFKVector,MFErrorHandler);

/*! \fn MFNKMatrix IMFGetBasisForStableInvariantSubspace(IMFFlow F, MFNVector u, MFKVector p,MFErrorHandler e);
 *  \brief Given a hyperbolic fixed point in a flow, finds a basis for the 
 *            stable invariant linear subspace of a hyperbolic fixed point.
 *
 *  \param F The flow.
 *  \param u The fixed point. 
 *  \param p The parameters of the flow for the fixed point.
 *  \param e A place to return errors.
 *  \returns A new MFNKMatrix containing an orthonormal basis for the stable invariant linear subspace of the hyperbolic fixed
 *             point.
 */
MFNKMatrix IMFGetBasisForStableInvariantSubspace(IMFFlow,MFNVector,MFKVector,MFErrorHandler);

/*! \fn void IMFOrthonormalizeBasis(int n,int m,double *du, MFErrorHandler e);
 *  \brief Given a basis, uses Gram-Schmidt to orthonormalize the basis.
 *
 *  \param n  The length of a basis vector.
 *  \param k  The number of basis vectors.
 *  \param du A matrix (column order) containing the basis vectors.
 *  \param e A place to return errors.
 *  \returns The orthogonal complement.
 */
void IMFOrthonormalizeBasis(int,int,double*,MFErrorHandler);

/*! \fn MFNKMatrix IMFGetBasisForOrthogonalComplement(MFNKMatrix Du,MFErrorHandler e);
 *  \brief Given a linear subspace, finds an orthonormal basis for its orthogonal complement.
 *
 *  \param Du An orthonormal basis a linear subspace.
 *  \param e A place to return errors.
 *  \returns The orthogonal complement.
 */
MFNKMatrix IMFGetBasisForOrthogonalComplement(MFNKMatrix,MFErrorHandler);

/*! @} */

/*! @} */

#ifdef __cplusplus
}
#endif

#endif
