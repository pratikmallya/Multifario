/*! \defgroup MFNRegion */
/*! \defgroup MFCubeNRegion */
/*! \defgroup MFRectangleNRegion */
/*! \defgroup MFDodecahedronMinusIcosahedronNRegion */
/*! \defgroup MF3dPolygonalNRegion */
/*! \defgroup MF3dEdgeNRegion */
/*! \defgroup MF3dPolyhedralNRegion */
/*! \defgroup MFTPBVPNRegion */

/*! \addtogroup MFNRegion
 *  @{
 */

/*! \class MFNRegion  MFNRegion.h MFNRegion.h
 *  \brief A subset of N-space.
 *
 *  An MFNRegion represents a subset of n-dimensional space. Given an MFNVector the MFNRegion
 *  provides a way of testing for "inside" and "outside". It is a base class, with constructors
 *  for some common regions.
 */

/*! \addtogroup MFCubeNRegion
 *  @{
 */

/*! \fn MFNRegion MFNRegionCreateCube(double x0,double y0,double z0,double x1,double y1,double z1,MFErrorHandler e);
 *  \brief Creates a cube in 3-space with edges parallel to the coordinate axes.
 *
 *  \param x0 The lower limit of the first coordinate.
 *  \param y0 The lower limit of the second coordinate.
 *  \param z0 The lower limit of the third coordinate.
 *  \param x1 The upper limit of the first coordinate.
 *  \param y1 The upper limit of the second coordinate.
 *  \param z1 The upper limit of the third coordinate.
 *  \param e A place to return errors.
 *  \return A new region.
 */
MFNRegion MFNRegionCreateCube(double x0,double y0,double z0,double x1,double y1,double z1,MFErrorHandler e);

/*! \fn MFNRegion MFNRegionCreateRectangle(double x0,double y0,double x1,double y1,MFErrorHandler e);
 *  \brief Creates a rectangle in 2-space with edges parallel to the coordinate axes.
 *
 *  \param x0 The lower limit of the first coordinate.
 *  \param y0 The lower limit of the second coordinate.
 *  \param x1 The upper limit of the first coordinate.
 *  \param y1 The upper limit of the second coordinate.
 *  \param e A place to return errors.
 *  \return A new region.
 */

/*! @} */

/*! \addtogroup MFRectangleNRegion
 *  @{
 */

MFNRegion MFNRegionCreateRectangle(double x0,double y0,double x1,double y1,MFErrorHandler e);

/*! \fn MFNRegion MFNRegionCreateHyperCube(int n,double R,MFErrorHandler e);
 *  \brief Creates a cube in n-space centered at the origin with edges parallel to the coordinate axes and length 2*R.
 *
 *  \param n The dimension of the cube.
 *  \param R Half of the edge length of the cube.
 *  \param e A place to return errors.
 *  \return A new region.
 */

/*! @} */

/*! \addtogroup MFCubeNRegion
 *  @{
 */

MFNRegion MFNRegionCreateHyperCube(int n,double R,MFErrorHandler e);

/*! \fn MFNRegion MFNRegionCreateHyperCubeByCorners(int n,MFNVector min,MFNVector max,MFErrorHandler e);
 *  \brief Creates a cube in n-space with edges parallel to the coordinate axes and corners max and min.
 *
 *  \param n The dimension of the cube.
 *  \param min The coordinates of min define the lower limits on each coordinate value.
 *  \param max The coordinates of max define the upper limits on each coordinate value.
 *  \param e A place to return errors.
 *  \return A new region.
 */
MFNRegion MFNRegionCreateHyperCubeByCorners(int n,MFNVector min,MFNVector max,MFErrorHandler e);

/*! \fn MFNRegion MFNRegionCreateDodecahedronMinusIcosahedron(MFErrorHandler e);
 *  \brief Creates a region in 333pace which lies between a Dodecahedron and a rotated Iscosahedron. It was used in 
 *         one of my projects, and you're welcome to it (it is not easy to mesh).
 *
 *  \param e A place to return errors.
 *  \return A new region.
 */

/*! @} */

/*! \addtogroup MFDodecahedronMinusIcosahedronNRegion
 *  @{
 */

MFNRegion MFNRegionCreateDodecahedronMinusIcosahedron(MFErrorHandler e);

/*! \fn int MFNRegionInterior(MFNRegion Omega,MFNVector u,MFErrorHandler e);
 *  \brief Determine if a point is inside the region.
 *
 *  \param Omega The region.
 *  \param u The point.
 *  \param e A place to return errors.
 *  \return TRUE if u is in Omega.
 */

/*! @} */

int MFNRegionInterior(MFNRegion Omega,MFNVector u,MFErrorHandler e);

/*! \fn void MFFreeNRegion(MFNRegion Omega,MFErrorHandler e);
 *  \brief Removes a reference to a region, and if the number of references goes to zero, deletes it.
 *
 *  \param Omega The region.
 *  \param e A place to return errors.
 */
void MFFreeNRegion(MFNRegion Omega,MFErrorHandler e);

/*! \fn void MFRefNRegion(MFNRegion omega, MFErrorHandler e);
 *  \brief Adds a references to the region.
 *
 *  \param omega The region being referenced.
 *  \param e A place to return errors.
 *  \sa ReferenceCounting MFFreeNRegion
 */
void MFRefNRegion(MFNRegion omega,MFErrorHandler e);

/*! \fn void MFWriteNRegion(FILE* fid,MFNRegion omega, MFErrorHandler e);
 *  \brief Writes a region to a file.
 *
 *  \param fid The file to write to.
 *  \param omega The region being queried.
 *  \param e A place to return errors.
 */
void MFWriteNRegion(FILE* fid,MFNRegion omega,MFErrorHandler e);

/*! \fn MFNRegion MFReadNRegion(FILE* fid, MFErrorHandler e);
 *  \brief Reads a region from a file.
 *
 *  \param fid The file to write to.
 *  \param e A place to return errors.
 *  \returns The region.
 */
MFNRegion MFReadNRegion(FILE* fid,MFErrorHandler e);

/*! \addtogroup MF3dPolygonalNRegion
 *  @{
 */

/*! \fn MFNRegion MFNRegionCreatePolygonal3dRegion(int nV,double *vertices,MFErrorHandler e);
 *  \brief Creates a planar polygonal region in 3-space with given vertices.
 *
 *  \param nV The number of vertices.
 *  \param vertices The vertices.
 *  \param e A place to return errors.
 *  \return A new region.
 */
MFNRegion MFNRegionCreatePolygonal3dRegion(int nV,double *vertices,MFErrorHandler e);

/*! @} */

/*! \addtogroup MF3dEdgeNRegion
 *  @{
 */

/*! \fn MFNRegion MFNRegionCreateEdge3dRegion(double *v0,double *v1,MFErrorHandler e);
 *  \brief Creates a 1 dimensional edge in 3-space.
 *
 *  \param v0 One endpoint of the edge v0[3].
 *  \param v1 The other endpoint v1[3].
 *  \param e A place to return errors.
 *  \return A new region.
 */
MFNRegion MFNRegionCreateEdge3dRegion(double *v0,double *v1,MFErrorHandler e);

/*! @} */

/*! \addtogroup MF3dPolyhedralNRegion
 *  @{
 */

/*! \fn MFNRegion MFNRegionCreatePolyhedral3dRegion(int nv,double *v,int nf,int *nFaceVertices,int **faceVertices,MFErrorHandler e);
 *  \brief Creates a region in 3-space which is the interior of a polyhedron.
 *
 *  \param nv The number of vertices.
 *  \param v An array containing the 3*nv coordinates of the vertices v[0+3*iv],v[1+3*iv],v[2+3*iv].
 *  \param nf The number of faces.
 *  \param nFaceVertices The number of vertices on each face, an array of length at least nf.
 *  \param faceVertices The vertex numbers of the vertices on each face faceVertices[if][ifv].
 *  \param e A place to return errors.
 *  \return A new region.
 */
MFNRegion MFNRegionCreatePolyhedral3dRegion(int nv,double *v,int nf,int *nFaceVertices,int **faceVertices,MFErrorHandler e);

/*! @} */

/*! \addtogroup MFCSGBallsNRegion
 *  @{
 */

/*! \fn MFNRegion MFNRegionCreateCSGBalls(int n,int nb,double *x0,double *R,int *dir, MFErrorHandler e);
 *  \brief Creates a region which is a simple set union of the interior or exterior of spherical balls.
 *
 *  \param n The dimension of the balls.
 *  \param nb The number of balls.
 *  \param x0 The center of the balls. Stored x0[0]=c0[0],..,x0[n-1]=c0[n-1], x0[n]=c1[0],...
 *  \param R The radii of the balls.
 *  \param dir if positive the interior is added, if negative the exterior.
 *  \param e A place to return errors.
 *  \return A new region.
 */
MFNRegion MFNRegionCreateCSGBalls(int n,int nb,double *x0,double *R,int *dir, MFErrorHandler e);

/*! @} */

/*! \addtogroup MFTPBVPNRegion
 *  @{
 */

/*! \fn MFNRegion MFNRegionCreateTPBVP(int nx,int nu,int np,double *p0,double *p1,double u0, double u1, MFErrorHandler e);
 *  \brief Creates a region which is a used for two point boundary value problems.
 *
 *  \param nx The number of mesh points.
 *  \param nu The dimension of the function.
 *  \param np The number of parameters.
 *  \param *p0 Lower bounds on the parameters.
 *  \param *p1 Upper bounds on the parameters.
 *  \param u0 Lower bound on the norm of the function.
 *  \param u1 Upper bound on the norm of the function.
 *  \param e A place to return errors.
 *  \return A new region.
 */
MFNRegion MFNRegionCreateTPBVP(int nx,int nu,int np,double *p0,double *p1,double u0, double u1, MFErrorHandler e);

/*! @} */


/*! @} */
