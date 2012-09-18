#ifndef __IMFExpansionPt__
#define __IMFExpansionPt__

#include <MFAtlas.h>
#include <MFErrorHandler.h>

#ifdef __cplusplus
 extern "C" {
#endif

/*! \addtogroup MFNVector
 *  @{
 */

/*! \addtogroup IMFExpansionNVector
 *  @{
 */

/*! \class IMFExpansionNVector IMFExpansionPt.h IMFExpansionPt.h 
 *  \brief An MFNVector type for points on a fat trajectory.
 *
 *  An IMFExpansionNVector is a point whose coordinates are the coefficients (constant, linear, quadratic) in a
 *  Taylor series for the surface containing a trajectory -- given by the constant term. Various information
 *  tags the point, including a link to the preceeding point along the fat trajectory, type information (an end point,
 *  interpolated or on the manifold of initial conditions, or a point interior to the trajectory).
 *  The remaining tags are used to determine if two charts overlap.
 */

/*! \fn MFNVector IMFCreateExpansionNVector(IMFExpansion E,double t, MFNVector sigma, int prevChart, int type,MFErrorHandler e);
 *  \brief Creates an MFNVector which has the coefficients of the expansion as coordinates. This is used for points
 *            on a fat trajectory.
 *
 *  \param E The expansion.
 *  \param t The value of the time coordinate to associate with the MFNVector (used for determined overlap).
 *  \param sigma The point (in the embedding space) of the initial "front" point to associate with the MFNVector (used for determined overlap).
 *  \param prevChart An identifying number carried with the MFNVector to indicate the preceding point on a fat trajectory.
 *  \param e A place to return errors.
 *  \param type Another identifier used to indicate endpoints of a fat trajectory and whether they are interpolation points
 *          or points on the manifold of initial conditions.
 */
MFNVector IMFCreateExpansionNVector(IMFExpansion E,double t, MFNVector sigma, int prevChart, int type,MFErrorHandler e);

/*! \fn IMFExpansion IMFExpansionNVGetE(MFNVector this,MFErrorHandler e);
 *  \brief Extracts the expansion from a point on a fat trajectory.
 *
 *  \param thisExpansion The MFNVector (must be an IMFExpansionNVector).
 *  \param e A place to return errors.
 *  \returns The expression.
 */
IMFExpansion IMFExpansionNVGetE(MFNVector this,MFErrorHandler e);

/*! \fn double IMFExpansionNVGetT(MFNVector this,MFErrorHandler e);
 *  \brief Extracts the time coordinate along the invariant manifold.
 *
 *  \param thisExpansion The MFNVector (must be an MFExpansionNVector).
 *  \param e A place to return errors.
 *  \returns The time value.
 */
double IMFExpansionNVGetT(MFNVector this,MFErrorHandler e);

/*! \fn MFNVector IMFExpansionNVGetSigma(MFNVector thisExpansion);
 *  \brief Extracts the sigma coordinate from an IMFExpansionNVector
 *
 *  \param thisExpansion The MFNVector (must be an IMFExpansionNVector).
 *  \param e A place to return errors.
 *  \returns The value of sigma.
 */
MFNVector IMFExpansionNVGetSigma(MFNVector this,MFErrorHandler e);

/*! \fn int IMFExpansionNVGetType(MFNVector this,MFErrorHandler e);
 *  \brief Extracts the type of an IMFExpansionNVector
 *
 *  \param this The MFNVector (must be an IMFExpansionNVector).
 *  \param e A place to return errors.
 *  \returns The type.
 */
int IMFExpansionNVGetType(MFNVector this,MFErrorHandler e);

/*! \fn MFKVector IMFExpansionNVGetS0(MFNVector this,MFErrorHandler e);
 *  \brief Extracts the manifold coordinates (the manifold of inital conditions) associated with the IMFExpansionNVector
 *
 *  \param this The MFNVector (must be an IMFExpansionNVector,MFErrorHandler e);
 *  \param e A place to return errors.
 *  \returns The point.
 */
MFKVector IMFExpansionNVGetS0(MFNVector this,MFErrorHandler e);

/*! \fn int IMFExpansionNVGetChart0(MFNVector this,MFErrorHandler e);
 *  \brief Extracts the chart (on the manifold of inital conditions) associated with the IMFExpansionNVector
 *
 *  \param this The MFNVector (must be an IMFExpansionNVector).
 *  \param e A place to return errors.
 *  \returns The chart number.
 */
int IMFExpansionNVGetChart0(MFNVector this,MFErrorHandler e);

/*! \fn int IMFExpansionNVGetPrev(MFNVector this,MFErrorHandler e);
 *  \brief Extracts the chart that lies on the fat before the one whose center is the IMFExpansionNVector.
 *
 *  \param this The MFNVector (must be an IMFExpansionNVector).
 *  \param e A place to return errors.
 *  \returns The nukmber of the preceeding chart.
 */
int IMFExpansionNVGetPrev(MFNVector this,MFErrorHandler e);

/*! \fn void IMFExpansionNVSetS0(MFNVector this, MFKVector S, MFErrorHandler e);
 *  \brief Sets the manifold coordinates (the manifold of inital conditions) associated with the IMFExpansionNVector
 *
 *  \param this The MFNVector (must be an IMFExpansionNVector).
 *  \param S0 The initial condition.
 *  \param e A place to return errors.
 *  \returns The nukmber of the preceeding chart.
 */
void IMFExpansionNVSetS0(MFNVector thisExpansion, MFKVector S, MFErrorHandler e);

/*! \fn void IMFExpansionNVSetChart0(MFNVector this, int chart, MFErrorHandler e);
 *  \brief Sets the manifold coordinates (the manifold of inital conditions) associated with the IMFExpansionNVector
 *
 *  \param this The MFNVector (must be an IMFExpansionNVector).
 *  \param chart The number of the intial chart (in the Atlas).
 *  \param e A place to return errors.
 *  \returns The nukmber of the preceeding chart.
 */
void IMFExpansionNVSetChart0(MFNVector thisExpansion, int chart0, MFErrorHandler e);

/*! @} */

/*! @} */

#ifdef __cplusplus
 }
#endif

#endif
