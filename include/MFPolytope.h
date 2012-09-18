/* 
    @(#)MFPolytope.h	1.6
    02/11/11 15:04:50
   
    PROGRAM NAME:  Manifold
   
    (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
    CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
   
    Please refer to the LICENSE file in the top directory

*/
/*      author: Mike Henderson mhender@watson.ibm.com */

#ifndef __NEWPOLYTOPE_H__
#define __NEWPOLYTOPE_H__
#include <MFBase.h>

#include <stdio.h>
#include <MFKVector.h>
#include <MFErrorHandler.h>

/*! \defgroup MFPolytope */

/*! \addtogroup MFChart
 *
 *  @{
 */

/*! \addtogroup MFPolytope
 *
 *  @{
 */

/*! \class MFPolytope 
 *  \brief A polytope, which represents a neighborhood on a manifold.
 */
struct MFPolytopeSt;
typedef struct MFPolytopeSt *MFPolytope;

#ifdef __cplusplus
 extern "C" {
#endif

/*! \fn MFPolytope MFCreateHyperCubeAtOrigin(int d,double halfedge, MFErrorHandler e);
 *  \brief Creates a d-dimensional hypercube whose faces are orthogonal to one coorindate axis, whose
 *         edges are length 2*halfedge, and whose center is at the origin.
 *
 *  \param d The dimension.
 *  \param halfedge Half the length of an edge.
 *  \param e A place to report errors.
 *  \returns A new polytope.
 */
MFPolytope MFCreateHyperCubeAtOrigin(int,double,MFErrorHandler);

/*! \fn MFPolytope MFCreateSimplexAtOrigin(int d,double side, MFErrorHandler e);
 *  \brief Creates a d-dimensional simplex. The last d-k coordinates of the kth vertex are zero.
 *
 *  \param d The dimension.
 *  \param side The length of an edge.
 *  \param e A place to report errors.
 *  \returns A new polytope.
 */
MFPolytope MFCreateSimplexAtOrigin(int,double,MFErrorHandler);

/*! \fn void MFSubtractHalfSpaceFromPolytope(MFPolytope P,int side,MFKVector n,double o, MFErrorHandler e);
 *  \brief Removes a half space from a polytope. The positive (side>0) side is the set of points x that
 *                       satisfy x.n-o>0.
 *
 *  \param P The original polytope (will be modified).
 *  \param side Which side, positive or negative of the hyperplane.
 *  \param n The normal of the hyperplane.
 *  \param o The distance from the origin along the normal to the hyperplane.
 *  \param e A place to report errors.
 */
void MFSubtractHalfSpaceFromPolytope(MFPolytope,int,MFKVector,double,MFErrorHandler);

/*! \fn void MFClipPolytope(MFPolytope P,int index,double *din,int mark, MFErrorHandler e);
 *  \brief Clips a polytope. index is the label to use for the new face, din[nv] is an array indicating which
 *             side of the plane each vertex lies. mark is simply an integer carried along with the face. multifario
 *             uses it to determine the polytop which lies opposite the face.
 *
 *  \param P The original polytope (will be modified).
 *  \param index The index used for this new face.
 *  \param din An array of length at least the number of vertices. Vertices with positive din are kept, negative are removed.
 *  \param mark A number carried along on the face.
 *  \param e A place to report errors.
 */
void MFClipPolytope(MFPolytope P,int,double*,int,MFErrorHandler);

/*! \fn void MFRefPolytope(MFPolytope polytope, MFErrorHandler e);
 *  \brief Adds a reference to a polytope.
 *
 *  \param polytope The polytope being referenced.
 *  \param e A place to report errors.
 *  \sa ReferenceCounting MFFreePolytope
 */
void MFRefPolytope(MFPolytope,MFErrorHandler);

/*! \fn void MFFreePolytope(MFPolytope polytope, MFErrorHandler e);
 *  \brief Frees a reference to a polytope, and deletes the polytope if there are no references left.
 *
 *  \param polytope The polytope being unreferenced.
 *  \param e A place to report errors.
 *  \sa ReferenceCounting MFRefPolytope
 */
void MFFreePolytope(MFPolytope,MFErrorHandler);

/*! \fn int MFPolytopeDimension(MFPolytope P, MFErrorHandler e);
 *  \brief returns the embedding dimension of the polytope.
 *
 *  \param P The polytope.
 *  \param e A place to report errors.
 *  \returns The dimension.
 */
int MFPolytopeDimension(MFPolytope,MFErrorHandler);

/*! \fn int MFPolytopeNumberOfVertices(MFPolytope P, MFErrorHandler e);
 *  \brief returns the number of vertices of the polytope.
 *
 *  \param P The polytope.
 *  \param e A place to report errors.
 *  \returns The dimension.
 */
int MFPolytopeNumberOfVertices(MFPolytope,MFErrorHandler);

/*! \fn void MFPolytopeVertex(MFPolytope P,int i,MFKVector v, MFErrorHandler e);
 *  \brief returns the embedding dimension of the polytope.
 *
 *  \param P The polytope.
 *  \param i Which vertex.
 *  \param v A place to put the coordinates of the vertex.
 *  \param e A place to report errors.
 */
void MFPolytopeVertex(MFPolytope,int,MFKVector,MFErrorHandler);

/*! \fn int MFPolytopeVertexIndex(MFPolytope P,int v,int i, MFErrorHandler e);
 *  \brief Returns an index (face number) of a vertex.
 *
 *  \param P The polytope.
 *  \param v Which vertex.
 *  \param i Which index.
 *  \param e A place to report errors.
 *  \returns The index.
 */
int MFPolytopeVertexIndex(MFPolytope,int,int,MFErrorHandler);

/*! \fn int MFPolytopeNumberOfVertexIndices(MFPolytope P,int v, MFErrorHandler e);
 *  \brief Returns the number of indices of a vertex (how many faces it lies on).
 *
 *  \param P The polytope.
 *  \param v Which vertex.
 *  \param e A place to report errors.
 *  \returns The number of vertices.
 */
int MFPolytopeNumberOfVertexIndices(MFPolytope,int,MFErrorHandler);

/*! \fn double MFPolytopeRadiusOfVertex(MFPolytope P,int v, MFErrorHandler e);
 *  \brief Returns the distance between a vertex and the origin.
 *
 *  \param P The polytope.
 *  \param v Which vertex.
 *  \param e A place to report errors.
 *  \returns The distance from the origin to the vertex.
 */
double MFPolytopeRadiusOfVertex(MFPolytope,int,MFErrorHandler);

/*! \fn double MFPolytopeLargestRadiusOfVertex(MFPolytope P, MFErrorHandler e);
 *  \brief Returns the largest distance between a vertex and the origin over all the vertices.
 *
 *  \param P The polytope.
 *  \param e A place to report errors.
 *  \returns The greatest distance from the origin to a vertex.
 */
double MFPolytopeLargestRadiusOfVertex(MFPolytope,MFErrorHandler);

/*! \fn int MFPolytopeNumberOfFaces(MFPolytope P, MFErrorHandler e);
 *  \brief Returns the number of faces of the polytope.
 *
 *  \param P The polytope.
 *  \param e A place to report errors.
 *  \returns The number of faces.
 */
int MFPolytopeNumberOfFaces(MFPolytope,MFErrorHandler);

/*! \fn int MFPolytopeFaceIndex(MFPolytope P,int f, MFErrorHandler e);
 *  \brief Returns the index of a face.
 *
 *  \param P The polytope.
 *  \param f Which face.
 *  \param e A place to report errors.
 *  \returns The index of the face.
 */
int MFPolytopeFaceIndex(MFPolytope,int,MFErrorHandler);

/*! \fn MFKVector MFPolytopeFaceNormal(MFPolytope P,int f, MFErrorHandler e);
 *  \brief Returns the nroaml to a face.
 *
 *  \param P The polytope.
 *  \param f Which face.
 *  \param e A place to report errors.
 *  \returns The normal to a face.
 */
MFKVector MFPolytopeFaceNormal(MFPolytope,int,MFErrorHandler);

/*! \fn double MFPolytopeFaceOrigin(MFPolytope P,int f, MFErrorHandler e);
 *  \brief Returns the origin of a face (the signed distance from the origin along the normal to the face).
 *
 *  \param P The polytope.
 *  \param f Which face.
 *  \param e A place to report errors.
 *  \returns The origin of the face.
 */
double MFPolytopeFaceOrigin(MFPolytope,int,MFErrorHandler);

/*! \fn int MFPolytopeNumberOfVerticesOnFace(MFPolytope P,int f, MFErrorHandler e);
 *  \brief Returns the number of vertices that lie on a face.
 *
 *  \param P The polytope.
 *  \param f Which face.
 *  \param e A place to report errors.
 *  \returns The number of vertices on the face.
 */
int MFPolytopeNumberOfVerticesOnFace(MFPolytope,int,MFErrorHandler);

/*! \fn int MFPolytopeInterior(MFPolytope P,MFKVector x, MFErrorHandler e);
 *  \brief Returns TRUE if x is inside the polytope.
 *
 *  \param P The polytope.
 *  \param x The point.
 *  \param e A place to report errors.
 *  \returns TRUE if the point is interior, FALSE otherwise.
 */
int MFPolytopeInterior(MFPolytope,MFKVector,MFErrorHandler);

/*! \fn MFPolytope MFPolytopeRemoveSmallEdges(MFPolytope P,double epsilon, MFErrorHandler e);
 *  \brief Returns a new polytope which has no edges shorter than epsilon. This is done by merging the endpoints of such edges.
 *
 *  \param P The polytope.
 *  \param epsilon The smallest edge allowed..
 *  \param e A place to report errors.
 *  \returns A new polytope with no short edges.
 */
MFPolytope MFPolytopeRemoveSmallEdges(MFPolytope,double,MFErrorHandler);

/*! \fn void MFWritePolytope(FILE* fid,MFPolytope polytope, MFErrorHandler e);
 *  \brief Writes a polytope to a file.
 *
 *  \param fid The file to write to.
 *  \param polytope The polytope being queried.
 *  \param e A place to report errors.
 */
void MFWritePolytope(FILE*,MFPolytope,MFErrorHandler);

void MFWritePolytope(FILE* fid,MFPolytope polytope, MFErrorHandler e);

/*! \fn MFPolytope MFReadPolytope(FILE* fid,MFAtlas A, MFErrorHandler e);
 *  \brief Reads a polytope from a file.
 *
 *  \param fid The file to write to.
 *  \param A The Atlas for the polytope.
 *  \param e A place to report errors.
 *  \returns polytope The polytope.
 */
MFPolytope MFReadPolytope(FILE*,MFErrorHandler);

/*! \fn int MFPolytopeClosestFace(MFPolytope P,MFKVector x,double *d, MFErrorHandler e);
 *  \brief Finds the closest face to a point.
 *
 *  \param P The polytope.
 *  \param x The point.
 *  \param d A pointer to double to set to the distance.
 *  \param e A place to report errors.
 *  \returns The face.
 */
int MFPolytopeClosestFace(MFPolytope,MFKVector,double*,MFErrorHandler);

/*! \fn void MFPolytopeSetVertexMark(MFPolytope P,int v,int mark, MFErrorHandler e);
 *  \brief Sets the mark associated with a vertex.
 *
 *  \param P The polytope.
 *  \param v The vertex.
 *  \param mark The new vertex mark.
 *  \param e A place to report errors.
 */
void MFPolytopeSetVertexMark(MFPolytope,int,int,MFErrorHandler);

/*! \fn int MFPolytopeGetVertexMark(MFPolytope P,int v,int mark, MFErrorHandler e);
 *  \brief Sets the mark associated with a vertex.
 *
 *  \param P The polytope.
 *  \param v The vertex.
 *  \param mark The new vertex mark.
 *  \param e A place to report errors.
 *  \returns The current vertex mark.
 */
int MFPolytopeGetVertexMark(MFPolytope,int,MFErrorHandler);

#ifdef __cplusplus
}
#endif

#endif

/*! @} */

/*! @} */
