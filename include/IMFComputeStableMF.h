/*! \defgroup InvariantManifolds */

/*! \addtogroup InvariantManifolds
 *  @{
 */

/*! \class IMFComputeInvariantManifold IMFComputeStableMF.h IMFComputeStableMF.h
 *  \brief A group of routines to compute invariant manifolds.
 *
 */

/*! \fn MFAtlas IMFComputeStableInvariantManifold(IMFFlow F,char *name, MFNVector u0, MFKVector p0, MFNRegion Omega, double eps, double dt, double tmax, int maxInterp, int maxCharts, double R0,MFErrorHandler e);
 * \brief Computes an atlas of charts that cover the stable manifold of a hyerbolic fixed point.
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
MFAtlas IMFComputeStableInvariantManifold(IMFFlow F,char *name, MFNVector u0, MFKVector p0, MFNRegion Omega, double eps, double dt, double tmax, int maxInterp, int maxCharts, double R0,MFErrorHandler e);

/*! \fn MFAtlas IMFComputeUnstableInvariantManifold(IMFFlow F,char *name, MFNVector u0, MFKVector p0, MFNRegion Omega, double eps, double dt, double tmax, int maxInterp, int maxCharts, double R0,MFErrorHandler e);
 * \brief Computes an atlas of charts that cover the unstable manifold of a hyerbolic fixed point.
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
MFAtlas IMFComputeUnstableInvariantManifold(IMFFlow F,char *name, MFNVector u0, MFKVector p0, MFNRegion Omega, double eps, double dt, double tmax, int maxInterp, int maxCharts, double R0,MFErrorHandler e);

/*! \fn MFAtlas IMFComputeInvariantManifold(IMFFlow F,MFKVector p0,char *name, MFAtlas c,MFNRegion Omega, double eps, double dt, double tmax, int maxInterp, int maxCharts, double Rmax,MFErrorHandler e);
 * \brief Computes an atlas of charts that cover the image of a manifold under a flow. A streamsurface.
 *
 * \param F The flow.
 * \param name A name to used for paging files and so forth.
 * \param p0 The parameters of the flow for the fixed point.
 * \param c The manifold of initial conditions.
 * \param Omega A region to bound the computation.
 * \param eps The tolerance on the size of the quadratic terms over a balls. Controls the stepsize.
 * \param dt The initial timestep to use along a fat trajectory.
 * \param tmax The upper limit on trajectory length.
 * \param maxInterp The upper limit on the number of interpolation performed.
 * \param maxCharts The upper limit on the number of charts in the atlas.
 * \param Rmax An upper limit to impose on the radius of the balls along the fat trajectories.
 * \param e A place to return errors.
 * \returns A new Atlas
 */
MFAtlas IMFComputeInvariantManifold(IMFFlow F,MFKVector p0,char *name, MFAtlas c,MFNRegion Omega, double eps, double dt, double tmax, int maxInterp, int maxCharts, double Rmax,MFErrorHandler e);

/*! @} */
