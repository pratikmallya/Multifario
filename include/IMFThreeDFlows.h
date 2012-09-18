#ifndef __IMFTHREEDFLOWS_H__
#define __IMFTHREEDFLOWS_H__

#include <MFErrorHandler.h>

#ifdef __cplusplus
 extern "C" {
#endif

/*! \defgroup Flows */
/*! \defgroup Source         */
/*! \defgroup SaddleFocus    */
/*! \defgroup SaddleCenter   */
/*! \defgroup Hyperbolic     */
/*! \defgroup Lorenz         */
/*! \defgroup StandardLorenz */

/*! \addtogroup Flows
 *  @{
 */

/*! \addtogroup Source
 *  @{
 */

/*! \fn IMFFlow IMFCreateThreeDSourceFlow(MFErrorHandler e);
 *  \brief A flow in three dimensional phase space with a source at the origin
 *
 *  \returns A new flow.
 */
IMFFlow IMFCreateThreeDSourceFlow(MFErrorHandler);

/*! @} */

/*! \addtogroup Hyperbolic
 *  @{
 */

/*! \fn IMFFlow IMFCreateThreeDHyperbolicFlow(MFErrorHandler e);
 *  \brief A flow in three dimensional phase space with a hyperbolic fixed point at the origin with the xy-plane
 *         the unstable invariant manifold and the z-axis the stable manifold.
 *
 *  \param e A place to return errors.
 *  \returns A new flow.
 */
IMFFlow IMFCreateThreeDHyperbolicFlow(MFErrorHandler);

/*! @} */

/*! \addtogroup SaddleFocus
 *  @{
 */

/*! \fn IMFFlow IMFCreateThreeDSaddleFocusFlow(MFErrorHandler e);
 *  \brief A flow in three dimensional phase space with a hyperbolic fixed point at the origin with the xy-plane
 *         the unstable invariant manifold with complex conjugate eigenvalues, and the z-axis the stable manifold.
 *         The stable eigenvalues are determined by the parameter eps, and are 1+i*eps and 1-i*eps.
 *
 *  \param e A place to return errors.
 *  \returns A new flow.
 */
IMFFlow IMFCreateThreeDSaddleFocusFlow(MFErrorHandler);

/*! @} */

/*! \addtogroup SaddleCenter
 *  @{
 */

/*! \fn IMFFlow IMFCreateThreeDSaddleCenterFlow(MFErrorHandler e);
 *  \brief A flow in three dimensional phase space with a fixed point at the origin with the xy-plane
 *         its center invariant manifold and the z-axis the stable manifold.
 *
 *  \param e A place to return errors.
 *  \returns A new flow.
 */
IMFFlow IMFCreateThreeDSaddleCenterFlow(MFErrorHandler);

/*! @} */

/*! \addtogroup Lorenz
 *  @{
 */

/*! \fn IMFFlow IMFCreateLorenzFlow(MFErrorHandler e);
 *  \brief The Lorenz flow. The three parameers are sigmma, rho and beta.
 *
 *  \param e A place to return errors.
 *  \returns A new flow.
 */
IMFFlow IMFCreateLorenzFlow(MFErrorHandler);

/*! @} */

/*! \addtogroup StandardLorenz
 *  @{
 */

/*! \fn IMFFlow IMFCreateStandardLorenzFlow(MFErrorHandler e);
 *  \brief The "Standard" Lorenz flow.  sigma=10.; rho=28.; beta=8./3.;
 *
 *  \param e A place to return errors.
 *  \returns A new flow.
 */
IMFFlow IMFCreateStandardLorenzFlow(MFErrorHandler);

/*! @} */

/*! @} */

#ifdef __cplusplus
}
#endif

#endif
