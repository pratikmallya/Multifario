#ifndef __IMFINTEGRATEFAT_H__
#define __IMFINTEGRATEFAT_H__

#include <MFAtlas.h>
#include <MFChart.h>
#include <IMFFlow.h>
#include <MFErrorHandler.h>

#ifdef __cplusplus
 extern "C" {
#endif

/*! \defgroup InvariantManifolds */

/*! \addtogroup InvariantManifolds
 *  @{
 */

/*! \addtogroup IMFFatTrajectory
 *  @{
 */

/*! \class IMFFatTrajectory IMFIntegrateFat.h IMFIntegrateFat.h
 *  \brief A fat trajectory is a trajectory in a flow with a local expansion attached that also evolves.
 */

/*! \fn void IMFExtendAtlasAlongFatTraj(MFAtlas I, IMFFlow F, char *name, MFNVector u0, MFKVector p0, double dt, double ti, double tf,double epsilon,MFNRegion Omega, int maxSteps, int maxCharts,double Rmax,int nSkip, MFErrorHandler e);
 *  \brief Integrates a fat trajectory and adds the points to an atlas.
 *
 * \param I The atlas to which new charts will be added.
 * \param F The flow.
 * \param name The name of the problem (used to page results to disk).
 * \param u0 The initial point on the fat trajectory.
 * \param p0 The flow parameters.
 * \param dt The initial time step.
 * \param ti The initial time.
 * \param tf The upper limit on time.
 * \param epsilon The maximum allowed size of the second order terms within a ball. Determines the stepsize.
 * \param Omega A region to limit the computation.
 * \param maxSteps The maximum number of steps to take along the fat trajectory.
 * \param maxCharts The maximum number of charts allowed in the atlas.
 * \param Rmax The largest radius of a ball.
 * \param nSkip The number of steps to make without adding points to the atlas.
 */
void IMFExtendAtlasAlongFatTraj(MFAtlas,IMFFlow,char*,MFNVector,MFKVector,double,double,double,double,MFNRegion,int,int,double,int, MFErrorHandler);

/*! \fn MFChart IMFStepAlongFatTraj(MFAtlas I, IMFFlow F, MFNVector u0, MFKVector p0, double R0, double dt, double epsilon,MFNRegion Omega, double Rmax, MFErrorHandler e)
 *  \brief Integrates one step along a fat trajectory.
 *
 * \param I The atlas to which new charts will be added.
 * \param F The flow.
 * \param u0 The initial point on the fat trajectory.
 * \param p0 The flow parameters.
 * \param R0 The distance requested between initial and final points of the step.
 * \param dt The initial time step.
 * \param epsilon The maximum allowed size of the second order terms within a ball. Determines the stepsize.
 * \param Omega A region to limit the computation.
 * \param Rmax The largest radius of a ball.
 */
MFChart IMFStepAlongFatTraj(MFAtlas,IMFFlow,MFNVector,MFKVector,double,double,double,MFNRegion,double,MFErrorHandler);

/*! @} */

/*! @} */

#ifdef __cplusplus
}
#endif

#endif
