/* 
 *  @(#)MFChart.h	1.9
 *  03/01/02 10:14:58
 * 
 *  PROGRAM NAME:  multifario
 * 
 *  (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
 *  CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
 * 
 *  Please refer to the LICENSE file in the top directory
 *
 *
 *      author: Mike Henderson mhender@watson.ibm.com
 *      date:   February 18, 1999
 */

#ifndef __MFCHART_H__
#define __MFCHART_H__
#include <MFBase.h>

#include <MFNVector.h>
#include <MFKVector.h>
#include <MFNKMatrix.h>
#include <MFImplicitMF.h>
#include <MFPolytope.h>
#include <MFErrorHandler.h>
#include <MFNRegion.h>
#include <stdio.h>

/*! \addtogroup MFAtlas
 *  @{
 */

/*! \defgroup MFChart */

/*! \addtogroup MFChart
 *  @{
 */


/*! \class MFChart MFChart.h MFChart.h
 *  \brief A chart, which represents a neighborhood on a manifold.
 *
 *  An MFChart represents an individual chart in the atlas of a manifold.
 *  The domain of the chart is a spherical ball in a k-dimensional Euclidean space, where k is the dimension
 *  of the manifold.
 *  The chart mapping maps from the domain of the chart onto the manifold, which is embedded in an n-dimensional
 *  Euclidean space, called the embedding space.
 *  In addition, each chart has a k-dimensional polyhedron (or polytope) which approximately partitions
 *  the manifold into disjoint pieces, and allows the boundary of the manifold to be determined.
 */
struct MFChartSt;
typedef struct MFChartSt *MFChart;

#ifdef __cplusplus
 extern "C" {
#endif

/*! \fn MFChart MFCreateChart(MFImplicitMF M,MFNVector u,MFNKMatrix Phi,double R, MFErrorHandler e);
 *  \brief Creates a chart.
 *
 *  \param M the implicitly defined manifold this chart describes.
 *  \param u the center of the new chart (an MFNVector).
 *  \param Phi a basis for the tangent space of the chart (an MFNKMatrix).
 *  \param R the radius of the chart.
 *  \param e A place to return errors.
 *  \return A new chart.
 */
MFChart MFCreateChart(MFImplicitMF,MFNVector,MFNKMatrix,double,MFErrorHandler);

/*! \fn MFChart MFCreateChartWIthCubeSize(MFImplicitMF M,MFNVector u,MFNKMatrix Phi,double R, double L, MFErrorHandler e);
 *  \brief Creates a chart.
 *
 *  \param M the implicitly defined manifold this chart describes.
 *  \param u the center of the new chart (an MFNVector).
 *  \param Phi a basis for the tangent space of the chart (an MFNKMatrix).
 *  \param R the radius of the chart.
 *  \param L the length of the edges of the initial cube around the chart.
 *  \param e A place to return errors.
 *  \return A new chart.
 */
MFChart MFCreateChartWithCubeSize(MFImplicitMF M,MFNVector u,MFNKMatrix TS, double R, double L, MFErrorHandler e);

/*! \fn void MFSubtractHalfSpaceFromChart(MFChart chart,int i,MFKVector n,double o, MFErrorHandler e);
 *  \brief Subtracts a half space from a chart, i.e. the polyhedron associated with the chart is intersected with the set x.n>o
 *
 *  \param chart the chart being modified.
 *  \param i a label to be associated with the boundary of the half space.
 *  \param n the normal to the bounding plane of the halfspace (an MFKVector s).
 *  \param o the origin of the bounding plane of the half space.
 *  \param e A place to return errors.
 */
void MFSubtractHalfSpaceFromChart(MFChart,int,MFKVector,double,MFErrorHandler);

MFChart MFCreateBoundaryChart(MFImplicitMF,MFNVector,MFNKMatrix,int,int,double,MFErrorHandler);

/*! \fn void MFRefChart(MFChart chart, MFErrorHandler e);
 *  \brief Adds a reference to the chart.
 *
 *  \param chart The chart being referenced.
 *  \param e A place to return errors.
 *  \sa ReferenceCounting MFFreeChart
 */
void MFRefChart(MFChart,MFErrorHandler);

/*! \fn void MFFreeChart(MFChart chart, MFErrorHandler e);
 *  \brief Frees a reference to the chart, and deletes the chart if there are no references left.
 *
 *  \param chart The chart being unreferenced.
 *  \param e A place to return errors.
 *  \sa ReferenceCounting MFRefChart
 */
void MFFreeChart(MFChart,MFErrorHandler);

/*! \fn MFPolytope MFChartPolytope(MFChart chart, MFErrorHandler e);
 *  \brief Returns the polyhedron (polytope) for the chart.
 *
 *  \param chart The chart being accessed.
 *  \param e A place to return errors.
 *  \returns The polyhedron of the chart.
 */
MFPolytope MFChartPolytope(MFChart,MFErrorHandler);

/*! \fn MFNVector MFChartCenter(MFChart chart, MFErrorHandler e);
 *  \brief Returns the center of the chart.
 *
 *  \param chart The chart being accessed.
 *  \param e A place to return errors.
 *  \returns The center of the chart.
 */
MFNVector MFChartCenter(MFChart,MFErrorHandler);

/*! \fn double MFChartRadius(MFChart chart, MFErrorHandler e);
 *  \brief Returns the radius of the domain of the chart.
 *
 *  \param chart The chart being accessed.
 *  \param e A place to return errors.
 *  \returns The radius of the domain of the chart.
 */
double MFChartRadius(MFChart,MFErrorHandler);

/*! \fn MFNKMatrix MFChartTangentSpace(MFChart chart, MFErrorHandler e);
 *  \brief Returns the nxk matrix whose columns are a basis for the tangent space of the chart.
 *
 *  \param chart The chart being accessed.
 *  \param e A place to return errors.
 *  \returns The basis of the tangent space of the chart.
 */
MFNKMatrix MFChartTangentSpace(MFChart,MFErrorHandler);

/*! \fn int MFChartEvaluate(MFChart chart,MFKVector s,MFNVector u, MFErrorHandler e);
 *  \brief Projects a point s in the domain of the chart, onto the manifold.
 *
 *  \param chart The chart defining the projection.
 *  \param s The k-vector defining a point in the tangent space.
 *  \param u (output) The point in the embedding space corresponding to s.
 *  \param e A place to return errors.
 *  \returns TRUE if the projection was succesful.
 */
int MFChartEvaluate(MFChart,MFKVector,MFNVector,MFErrorHandler);

/*! \fn int MFChartInterior(MFChart chart,MFKVector s, MFErrorHandler e);
 *  \brief Tests to see if a point s is in the domain of the chart (the spherical ball, not the polyhedron).
 *
 *  \param chart The chart being queried.
 *  \param s The k-vector defining a point in the tangent space.
 *  \param e A place to return errors.
 *  \returns TRUE if s is in the domain of the chart.
 */
int MFChartInterior(MFChart,MFKVector,MFErrorHandler);

/*! \fn int MFChartHasBoundary(MFChart chart, MFErrorHandler e);
 *  \brief Determines if the chart contributes to the boundary of the atlas of charts.
 *
 *  \param chart The chart.
 *  \param e A place to return errors.
 *  \returns Returns TRUE if the chart contributes to the boundary of the atlas of charts.
 */
int MFChartHasBoundary(MFChart,MFErrorHandler);

/*! \fn int MFChartK(MFChart chart, MFErrorHandler e);
 *  \brief Returns the dimension of the manifold, which is also the dimension of the domain of the chart mapping.
 *
 *  \param chart The chart being queried.
 *  \param e A place to return errors.
 */
int MFChartK(MFChart,MFErrorHandler);

/*! \fn int MFChartN(MFChart chart, MFErrorHandler e);
 *  \brief Returns the dimension of the embedding space of the manifold, which is the dimension of the range of the chart mapping.
 *
 *  \param chart The chart being queried.
 *  \param e A place to return errors.
 */
int MFChartN(MFChart,MFErrorHandler);

/*! \fn void MFChartProjectIntoTangentSpace(MFChart chart,MFNVector u,MFKVector s, MFErrorHandler e);
 *  \brief Projects a point in the embedding space onto the tangent space of the chart.
 *
 *  \param chart The chart.
 *  \param u The n-vector defining a point in the embedding space.
 *  \param s (output) The k-vector of the point in the tangent space.
 *  \param e A place to return errors.
 */
void MFChartProjectIntoTangentSpace(MFChart,MFNVector,MFKVector,MFErrorHandler);

/*! \fn void MFChartProjectVectorIntoTangentSpace(MFChart chart,MFNVector u,MFKVector s, MFErrorHandler e);
 *  \brief Projects a point in the embedding space onto the tangent space of the chart.
 *
 *  \param chart The chart.
 *  \param u The n-vector defining a point in the embedding space, relative to the center of the chart.
 *  \param s (output) The k-vector of the point in the tangent space.
 *  \param e A place to return errors.
 */
void MFChartProjectVectorIntoTangentSpace(MFChart,MFNVector,MFKVector,MFErrorHandler);

/*! \fn void MFChartPointInTangentSpace(MFChart chart,MFKVector s,MFNVector u, MFErrorHandler e);
 *  \brief Finds the point in the embedding space which lies on the linear approximation of the chart mapping.
 *
 *  \param chart The chart.
 *  \param s The k-vector of the point in the tangent space.
 *  \param u (output) The n-vector defining a point in the embedding space which lies on the linear approximation of the manifold
 *       at the center of the chart.
 *  \param e A place to return errors.
 */
void MFChartPointInTangentSpace(MFChart,MFKVector,MFNVector,MFErrorHandler);

/*! \fn void MFWriteChart(FILE* fid,MFChart chart, MFErrorHandler e);
 *  \brief Writes a chart to a file.
 *
 *  \param fid The file to write to.
 *  \param chart The chart being queried.
 *  \param e A place to return errors.
 */
void MFWriteChart(FILE*,MFChart,MFErrorHandler);

struct MFAtlasSt;

/*! \fn MFChart MFReadChart(FILE* fid,MFAtlas A, MFErrorHandler e);
 *  \brief Reads a chart from a file.
 *
 *  \param fid The file to write to.
 *  \param A The Atlas for the chart.
 *  \param e A place to return errors.
 *  \returns chart The chart.
 */
MFChart MFReadChart(FILE*,struct MFAtlasSt*,MFErrorHandler);

int MFChartGetPositionInBoundaryList(MFChart,MFErrorHandler);
void MFChartSetPositionInBoundaryList(MFChart,int,MFErrorHandler);

/*! \fn MFImplicitMF MFChartGetManifold(MFChart chart, MFErrorHandler e);
 *  \brief Returns the implicitly defined manifold which the chart represents.
 *
 *  \param chart The chart being queried.
 *  \param e A place to return errors.
 *  \returns Returns the implicitly defined manifold which the chart represents.
 */
MFImplicitMF MFChartGetManifold(MFChart,MFErrorHandler);

/*! \fn int MFChartIsSingular(MFChart chart, MFErrorHandler e);
 *  \brief Returns TRUE if the chart has been marked as singular (i.e. the chart mapping may not be one to one.
 *
 *  \param chart The chart being queried.
 *  \param e A place to return errors.
 *  \returns Returns TRUE if the chart has been marked as singular (i.e. the chart mapping may not be one to one.
 */
int MFChartIsSingular(MFChart,MFErrorHandler);

void MFChartSetSingular(MFChart,MFErrorHandler);
void MFChartSetNonSingular(MFChart,MFErrorHandler);
void MFChartSetNonSingular(MFChart,MFErrorHandler);

void MFChartClean(MFChart,MFErrorHandler);
int MFChartPageOut(MFChart,FILE*,int,MFErrorHandler);

/*! \fn int MFChartPaged(MFChart chart, MFErrorHandler e);
 *  \brief Returns TRUE if the chart has been paged to file.
 *
 *  \param chart The chart being queried.
 *  \param e A place to return errors.
 *  \returns Returns TRUE if the chart has been paged to file.
 */
int MFChartPaged(MFChart,MFErrorHandler);

/*! \fn int MFChartNearBoundary(MFChart chart, MFErrorHandler e);
 *  \brief Returns TRUE if the chart has been marked as near the boundary of the atlas.
 *
 *  \param chart The chart being queried.
 *  \param e A place to return errors.
 *  \returns TRUE if the chart has been marked as near the boundary of the atlas.
 */
int MFChartNearBoundary(MFChart,MFErrorHandler);

void MFChartSetNearBoundary(MFChart,int,MFErrorHandler);

/*! \fn int MFChartNumberOfFaces(MFChart chart, MFErrorHandler e);
 *  \brief Returns the number of faces of the chart polyhedron (i.e. the number of neighbors of the chart.
 *
 *  \param chart The chart being queried.
 *  \param e A place to return errors.
 *  \returns The number of faces.
 */
int MFChartNumberOfFaces(MFChart,MFErrorHandler);

/*! \fn double MFChartSuggestedRadius(MFChart chart, MFErrorHandler e);
 *  \brief Returns the hint at the radius of the chart domain.
 *
 *  \param chart The chart being queried.
 *  \param e A place to return errors.
 *  \returns The suggest chart radius.
 */
double MFChartSuggestedRadius(MFChart,MFErrorHandler);

void MFSetSuggestedChartRadius(MFChart,double,MFErrorHandler);

/*! \fn int MFChartReferenceNumber(MFChart chart, MFErrorHandler e);
 *  \brief Returns the reference number identifying the chart.
 *
 *  \param chart The chart being queried.
 *  \param e A place to return errors.
 *  \returns The reference number identifying the chart.
 */
int MFChartReferenceNumber(MFChart,MFErrorHandler);

void MFChartSetReferenceNumber(MFChart,int,MFErrorHandler);

/*! \fn double MFVolumeOfChart(MFChart chart,MFNRegion Omega,MFErrorHandler e);
 *  \brief Finds the volume of the polyhedra of a chart.
 *
 *  \param chart The chart.
 *  \param Omega A region to limit the calculation.
 *  \param e A place to return errors.
 *  \returns The volume.
 */
double MFVolumeOfChart(MFChart chart,MFNRegion Omega,MFErrorHandler e);

MFChart MFChartGetNextChart(MFChart thisChart, MFErrorHandler e);
MFChart MFChartGetPrevChart(MFChart thisChart, MFErrorHandler e);

#ifdef __cplusplus
}
#endif

/*! @} */

/*! @} */

#endif
