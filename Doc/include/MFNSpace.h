/*! \defgroup MFNSpace */

/*! \addtogroup MFNSpace
 *  @{
 */

/*! \class MFNSpace MFNSpace.h MFNSpace.h
 *  \brief An embedding space (a metric space), which provides inner products and distances.
 *
 *  An MFNSpace represents an embedding space. It is a base class -- the user cannot create
 *  base class spaces directly. An implicitly defined manifold will create an appropriate space,
 *  and the user can query an atlas or an implicitly defined manifold to find the space. Otherwise,
 *  constructors for some common spaces are provided.
 */

/*! \fn MFNSpace MFCreateNSpace(int n,MFErrorHandler e);
 *  \brief Creates a Euclidean n-space.
 *
 *  \param n The dimension of the space.
 *  \param e A place to return errors.
 *  \returns A new n-space.
 */
MFNSpace MFCreateNSpace(int n,MFErrorHandler e);

/*! \fn MFNSpace MFCreateTPBVPNSpace(int nx,int nu,int np, MFErrorHandler e);;
 *  \brief Creates a space for the solutions of an MFTPBVP (two point boundary value problem).
 *
 *  \param nx The number of mesh intervals.
 *  \param nu The number of functions defined on the mesh.
 *  \param np The number of scalar parameters.
 *  \param e A place to return errors.
 *  \returns A new n-space.
 */
MFNSpace MFCreateTPBVPNSpace(int nx,int nu,int np, MFErrorHandler e);;

/*! \fn double MFNSpaceInner(MFNSpace space,MFNVector a,MFNVector b, MFErrorHandler e);;
 *  \brief Calculates the inner product of two vectors using the inner product of the space.
 *
 *  \param space The space.
 *  \param a The first n-vector.
 *  \param b The second n-vector.
 *  \param e A place to return errors.
 *  \returns <a,b>.
 */
double MFNSpaceInner(MFNSpace space,MFNVector a,MFNVector b, MFErrorHandler e);;

/*! \fn double MFNSpaceDistance(MFNSpace space,MFNVector a,MFNVector b, MFErrorHandler e);;
 *  \brief Calculates the distance between two vectors using the metric of the space.
 *
 *  \param space The space.
 *  \param a The origin.
 *  \param b The destination.
 *  \param e A place to return errors.
 *  \returns dist(a,b).
 */
double MFNSpaceDistance(MFNSpace space,MFNVector a,MFNVector b, MFErrorHandler e);;

/*! \fn double MFNSpaceTangentDistance(MFNSpace space,MFNKMatrix A,MFNKMatrix B, MFErrorHandler e);;
 *  \brief Calculates the "distance" between two "Tangent Spaces". Really just a measure of how different the span of the
 *         bases are.
 *
 *  \param space The space.
 *  \param A The origin.
 *  \param B The destination.
 *  \param e A place to return errors.
 *  \returns dist(A,B).
 */
double MFNSpaceTangentDistance(MFNSpace space,MFNKMatrix A,MFNKMatrix B, MFErrorHandler e);;

/*! \fn void MFNSpaceDirection(MFNSpace space,MFNVector a,MFNVector b,MFNVector c, MFErrorHandler e);;
 *  \brief Calculates the tangent to the geodesic between two points.
 *
 *  \param space The space.
 *  \param a The origin.
 *  \param b The destination.
 *  \param c The direction (infinitesimal) from a to b..
 *  \param e A place to return errors.
 */
void MFNSpaceDirection(MFNSpace space,MFNVector a,MFNVector b,MFNVector c, MFErrorHandler e);;

/*! \fn void MFNSpaceAdd(MFNSpace space,MFNVector a,MFNVector b,MFNVector c, MFErrorHandler e);;
 *  \brief Calculates the sum of two vectors in the embedding space.
 *
 *  \param space The space.
 *  \param a The first n-vector.
 *  \param b The second n-vector.
 *  \param c The sum a+b.
 *  \param e A place to return errors.
 */
void MFNSpaceAdd(MFNSpace space,MFNVector a,MFNVector b,MFNVector c, MFErrorHandler e);;

/*! \fn void MFNSpaceScale(MFNSpace space,double scalar,MFNVector a,MFNVector b, MFErrorHandler e);;
 *  \brief Calculates the product of a scalar and a vector in the embedding space.
 *
 *  \param space The space.
 *  \param scalar The scalar.
 *  \param a The n-vector.
 *  \param b The product scalar*a.
 *  \param e A place to return errors.
 */
void MFNSpaceScale(MFNSpace space,double scalar,MFNVector a,MFNVector b, MFErrorHandler e);;

/*! \fn void *MFNSpaceGetData(MFNSpace space, MFErrorHandler e);;
 *  \brief Returns the internal representation of an embedding space. Be careful!
 *
 *  \param space The space.
 *  \param e A place to return errors.
 *  \returns A pointer to the data block.
 */
void *MFNSpaceGetData(MFNSpace space, MFErrorHandler e);;

/*! \fn char *MFNSpaceGetId(MFNSpace space, MFErrorHandler e);;
 *  \brief Query the type an embedding space. Don't free the string, and don't change it.
 *
 *  \param space The space.
 *  \param e A place to return errors.
 *  \returns The character string with the ID.
 */
char *MFNSpaceGetId(MFNSpace space, MFErrorHandler e);;

/*! \fn void MFRefNSpace(MFNSpace n-space, MFErrorHandler e);
 *  \brief Adds a reference to the n-space.
 *
 *  \param space The n-space being referenced.
 *  \param e A place to return errors.
 *  \sa ReferenceCounting MFFreeNSpace
 */
void MFRefNSpace(MFNSpace n-space,MFErrorHandler e);

/*! \fn void MFFreeNSpace(MFNSpace space, MFErrorHandler e);
 *  \brief Frees a reference to the n-space, and deletes the n-space if there are no references left.
 *
 *  \param space The n-space being unreferenced.
 *  \param e A place to return errors.
 *  \sa ReferenceCounting MFRefNSpace
 */
void MFFreeNSpace(MFNSpace space,MFErrorHandler e);

/*! \fn void MFWriteNSpace(FILE* fid,MFNSpace space, MFErrorHandler e);
 *  \brief Writes a n-space to a file.
 *
 *  \param fid The file to write to.
 *  \param space The n-space being queried.
 *  \param e A place to return errors.
 */
void MFWriteNSpace(FILE* fid,MFNSpace space,MFErrorHandler e);

/*! \fn MFNSpace MFReadNSpace(FILE* fid,MFAtlas A, MFErrorHandler e);
 *  \brief Reads a n-space from a file.
 *
 *  \param fid The file to write to.
 *  \param A The Atlas for the n-space.
 *  \param e A place to return errors.
 *  \returns n-space The n-space.
 */
MFNSpace MFReadNSpace(FILE* fid,MFAtlas A,MFErrorHandler e);

/*! @} */
