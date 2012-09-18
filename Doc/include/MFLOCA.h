/*! \defgroup MFLOCANSpace */
/*! \defgroup MFLOCANVector */
/*! \defgroup MFLOCAIMF */
/*! \defgroup MFLOCANRegion */

/*! \addtogroup MFNSpace
 *  @{
 */

/*! \addtogroup MFLOCANSpace
 *  @{
 */

/*! \fn MFNSpace MFCreateLOCANSpace(void *data,MFErrorHandler e,MFErrorHandler e);
 *  \brief Creates a NSpace for LOCA
 *
 *  \param data A LOCA object.
 *  \param e A place to return errors.
 *  \returns A new MFNSpace.
 */
MFNSpace MFCreateLOCANSpace(void *data,MFErrorHandler e,MFErrorHandler e);

/*! @} */

/*! @} */

/*! \addtogroup MFNVector
 *  @{
 */

/*! \addtogroup MFLOCANVector
 *  @{
 */

/*! \fn MFNVector MFCreateLOCANVectorWithData(int nx,double *x,int np,double *p,MFErrorHandler e,MFErrorHandler e);
 *  \brief Creates an implicitly defined manifold for LOCA
 *
 *  \param nx The number of unknowns.
 *  \param x An array of values for the unknowns.
 *  \param np The number of parameters.
 *  \param p An array of values for the parameters.
 *  \param e A place to return errors.
 *  \returns A new MFNVector.
 */
MFNVector MFCreateLOCANVectorWithData(int nx,double *,int np,double *p,MFErrorHandler e,MFErrorHandler e);

/*! \fn MFNVector MFCreateLOCANVector(int nx,int np,MFErrorHandler e,MFErrorHandler e);
 *  \brief Creates an implicitly defined manifold for LOCA
 *
 *  \param nx The number of unknowns.
 *  \param np The number of parameters.
 *  \param e A place to return errors.
 *  \returns A new MFNVector.
 */
MFNVector MFCreateLOCANVector(int nx,int np,MFErrorHandler e,MFErrorHandler e);

/*! @} */

/*! @} */

/*! \addtogroup MFNRegion
 *  @{
 */

/*! \addtogroup MFLOCANRegion
 *  @{
 */

/*! \fn MFNRegion MFNRegionCreateLOCA(void *data,MFErrorHandler e,MFErrorHandler e);
 *  \brief Creates a NSpace for LOCA
 *
 *  \param data A LOCA (struct con_struct *)object.
 *  \param e A place to return errors.
 *  \returns A new MFNRegion.
 */
MFNRegion MFNRegionCreateLOCA(void *data,MFErrorHandler e,MFErrorHandler e);

/*! @} */

/*! @} */

/*! \addtogroup MFImplicitMF
 *  @{
 */

/*! \addtogroup MFLOCAIMF
 *  @{
 */

/*! \fn MFImplicitMF MFIMFCreateLOCA(int np,void *data,MFErrorHandler e,MFErrorHandler e);
 *  \brief Creates an implicitly defined manifold for LOCA
 *
 *  \param np The number of parameters.
 *  \param data A LOCA (struct con_struct *)object.
 *  \param e A place to return errors.
 *  \returns A new MFImplicitMF.
 */
MFImplicitMF MFIMFCreateLOCA(int np,void *data,MFErrorHandler e,MFErrorHandler e);

/*! @} */

/*! @} */
