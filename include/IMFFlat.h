/*! \defgroup IMFFlat */

/*! \addtogroup MFImplicitMF
 *  @{
 */

/*! \addtogroup IMFFlat
 *  @{
 */

/*! \fn MFImplicitMF IMFCreateFlat(int n, int k,MFErrorHandler e);
 *  \brief A flat k-manifold embedded in Euclidean n-space
 *  
 *  \param n The dimension of the embedding space.
 *  \param k The dimension of the manifold.
 *  \param e A place to return errors.
 *  \returns A new MFImplicitMF.
 */
MFImplicitMF IMFCreateFlat(int n, int k,MFErrorHandler e);

/*!    @} */

/*!    @} */

