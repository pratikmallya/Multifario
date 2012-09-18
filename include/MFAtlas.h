/* 
    %W%
    %D% %T%

    PROGRAM NAME:  Manifold
   
    (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
    CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
   
    Please refer to the LICENSE file in the top directory

*/
/*      author: Mike Henderson mhender@watson.ibm.com */
/*      date:   February 22, 1999                     */

#ifndef __MFATLAS_H__
#define __MFATLAS_H__

#include <stdlib.h>
#include <MFBase.h>
#include <stdio.h>

#include <MFKVector.h>
#include <MFNVector.h>
#include <MFNRegion.h>
#include <MFNKMatrix.h>
#include <MFChart.h>
#include <MFImplicitMF.h>
#include <MFErrorHandler.h>

/*! \mainpage Multifario

<table width="90%" border=0>
<tr> <td><img src="http://multifario.sourceforge.net/LorenzU0Icon.jpg"></td> <td><img src="http://multifario.sourceforge.net/Torus1400Icon.jpg"></td> </tr>
</table>

<strong>Multifario</strong> is a set of subroutines and data structures for computing
manifolds that occur in Dynamical Systems. Fixed points, periodic orbits,
heteroclinic and homoclinic connections, and many types of singular motions
have been formulated as solutions of systems of parameterized algebraic
or two point boundary value problems. (See the set of objects which can
be computed with <a href="http://sourceforge.net/projects/auto2000/">AUTO</a>.)

<p>The code is structured so that the user can provide their own solvers,
but two stand-alone solvers are provided. Both use the Lapack linear equation
solvers (available through <a href="http://www.netlib.org/">Netlib</a>).
The first solver is for algebraic systems, the second is for two point
boundary value problems, using a fixed mesh and Keller's Box scheme.

<p>Interfaces to <a href="http://sourceforge.net/projects/auto2000/">AUTO</a>,
(a full package for algebraic and two point boundary value problems), and
<a href="http://www.cs.sandia.gov/projects/loca/main.html">LOCA</a>
(Sandia's continuation code) are available. Note: the AUTO interface is
not yet quite a full implementation --
<p>The output of the code is a "plotfile" and/or "centerfile", or an "atlas".
The plotfile is a list of polyhedra which cover the manifold, whose vertices
are have been projected via a routine the user may provide. The centerfile
is a list of the points on the manifold at the center of each polyhedron.
The atlas is a full dump of the data structure which the algorithm uses.
<p>Basic conversion utilities are included for converting plotfiles to
<a href="http://www.opendx.org/">DataExplorer</a>
format, and to <a href="http://www.povray.org/">POV-RAY</a> files. In addition
a simple z-buffer code is included which will render a plotfile to a Tiff
file, if&nbsp; the <a href="http://www.libtiff.org/">LIBTIFF</a> library
is available, or Postscript.

<p>I've added the support for computing invariant
manifolds. This is the implementation of the algorithm described in the 
paper

<p><center><a href=http://epubs.siam.org/sam-bin/dbq/article/60289>
"Computing Invariant Manifolds by Integrating Fat Trajectories"</a>, SIADS v4(4), 2005 pages 832-882.</center> 
<p>The packaging is a bit primitive yet, and documentation is still being written.

<p align=right> <a href="http://www.research.ibm.com/people/h/henderson/">Mike Henderson</a>

<h2>Documentation</h2>

<p>The algorithm is described in detail the paper
<p>Henderson, Michael E., "Multiple Parameter Continuation: Computing
Implicitly Defined k-manifolds", IJBC v12(3), pages 451-76.</p><p>
I've put together a quick <a href="http://multifario.sourceforge.net/AlgorithmOverview.html">overview
of the algorithm</a>, and a full <a href="http://www.fields.utoronto.ca/audio/01-02/dynamsys/henderson/">presentation</a>
on the algorithm is available online (thanks to the <a href="http://www.fields.utoronto.ca/">Fields
Institute</a>). This was at the Fields Institute in Toronto in December,
2001 as part of their <a href="http://www.fields.utoronto.ca/programs/scientific/01-02/numerical/dynamsys/">Workshop
on Computational Challenges in Dynamical Systems.</a>
<p>An application, where the algorithm is used to compute configurations
of Clamped Elastica are computed may be found in the paper
<br>&nbsp;
<p><center>Michael E. Henderson and Sebastion Neukirch, <a href="http://www.lps.ens.fr/~neukirch/classif_clamped_elas_part2.html">"Classification
of the spatial equilibria of the clamped elastica"</a>, IJBC 14(4): 1223-39, 2003.</center>

<p> A TeX file with a description of how to use the code, and how to code your
own problems is included in the file Doc/MF.tex (<a href="http://multifario.sourceforge.net/MF.pdf">PDF
version</a>).

<h2> Examples</h2>

<p>&nbsp;
<table BORDER COLS=2 WIDTH="80%" >
<tr>
<td VALIGN=TOP>Computing a 1-manifold embedded in 2-space:
<p>&nbsp;&nbsp;&nbsp;&nbsp;<a href="http://multifario.sourceforge.net/Circle.html">
Circle</a> - 1.3Mb animated gif

<p>Computing a 2-manifold embedded in 2-space:
<p>&nbsp;&nbsp;&nbsp;&nbsp; <a href="http://multifario.sourceforge.net/Plane.html">Plane</a>
- 1.9Mb animated gif
<p>and 2-manifolds embedded in 3-space:
<p>&nbsp;&nbsp;&nbsp;&nbsp; <a href="http://multifario.sourceforge.net/Sphere.html">Sphere</a>
- 3.9Mb animated gif&nbsp; (<a href="http://multifario.sourceforge.net/Sphere.mov">.mov
version</a>, 10.6Mb)
<br>&nbsp;&nbsp;&nbsp;&nbsp; <a href="http://multifario.sourceforge.net/Torus.html">Torus</a>

- 4.9Mb animated gif (<a href="http://multifario.sourceforge.net/Torus.mov">.mov
version</a>, 32.5Mb)
<p>Computing a 3-manifold embedded in 3-space:
<p>&nbsp;&nbsp;&nbsp;&nbsp; <a href="http://multifario.sourceforge.net/Cube.html">Cube</a>
- 8.6 Mb&nbsp; animated gif</td>

<td VALIGN=TOP>Computing a manifold of stationary points
<p><a href="http://multifario.sourceforge.net/ClampedElastica.html">Equilibrium
of clamped elastica</a></p><p>
Computing a manifold of periodic motions
<p><a href="http://multifario.sourceforge.net/CoupledPendula.html">Periodic
Motions of Coupled Pendula</a></p><p>

Computing Invariant Manifolds
<p><a href="http://multifario.sourceforge.net/Lorenz.html">Invariant
Manifolds in the Lorenz System</a></p><p>
</td>
</tr>
</table>
The raw output of the examples can be seen in these images:
<center><table>
<tr>
<td>ComputeLine</td>

<td><a href="http://multifario.sourceforge.net/Circle.jpg">ComputeCircle</a></td>

<td><a href="http://multifario.sourceforge.net/Transcritical.jpg">ComputeTranscritical</a></td>

<td><a href="http://multifario.sourceforge.net/Domokos.jpg">ComputeDomokos</a> (<a href="http://multifario.sourceforge.net/Domokos2.jpg">further</a>)</td>
</tr>

<tr>
<td><a href="http://multifario.sourceforge.net/Plane.jpg">ComputePlane</a></td>

<td><a href="http://multifario.sourceforge.net/Sphere.jpg">ComputeSphere</a></td>

<td><a href="http://multifario.sourceforge.net/ComplexCusp.jpg">ComputeCusp</a></td>

<td><a href="http://multifario.sourceforge.net/Rod.jpg">ComputeRod</a></td>
</tr>

<tr>
<td><a href="http://multifario.sourceforge.net/PlaneClip.jpg">ComputePlaneClip</a></td>

<td><a href="http://multifario.sourceforge.net/SphereSub.jpg">ComputeSphereSub</a></td>

<td><a href="http://multifario.sourceforge.net/TaylorA.jpg">ComputeTaylorA</a></td>

<td></td>
</tr>

<tr>
<td><a href="http://multifario.sourceforge.net/3Space.jpg">Compute3Space</a></td>

<td><a href="http://multifario.sourceforge.net/Torus.jpg">ComputeTorus</a></td>

<td><a href="http://multifario.sourceforge.net/Taylor24.jpg">ComputeTaylor24</a></td>

<td></td>
</tr>

<tr>
<td>Compute4Space</td>

<td><a href="http://multifario.sourceforge.net/GenusTwo.jpg">ComputeGenusTwo</a></td>

<td></td>

<td></td>
</tr>

<tr>
<td></td>

<td><a href="http://multifario.sourceforge.net/SpherePacking3.jpg">ComputeSpherePacking</a></td>

<td></td>

<td></td>
</tr>
</table></center>
</p>

<h2>
Downloading</h2>

<p>Multifario is now available through the SourceForge code repository.
It can be downloaded as a <a href="http://multifario.sourceforge.net/Tarballs">gzipped
tar file</a>, or via <a href="http://multifario.sourceforge.net/faqs.html#ObtainSrcCode">SVN</a>.
<p>The package is released under the Common Public License (CPL); see
<a href="http://www.opensource.org/licenses/cpl.html">http://www.opensource.org/licenses/cpl.html</a>

or the enclosed LICENSE file.&nbsp; As stated in the license, you can use
the provided code without charge, but if you modify the code and you distribute
your modifications in some way (e.g. in an executable), you must make those
modifications public.</p>

<h2>
Installation</h2>

<p>To install
<p>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; download the source code.
<br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; edit the file share/config.site
to give local lib and include dirs for Lapack. (Last two lines in the file).
<br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; run the configure script
-- "./configure"
<br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; make the libraries, utiliies
and examples -- "make"</p>

 */
/*! \defgroup MFAtlas */

/*! \addtogroup MFAtlas
 *  @{
 */

#ifdef __cplusplus
 extern "C" {
#endif

/*! \class MFAtlas MFAtlas.h MFAtlas.h
 *  \brief An atlas (or list) of charts, which represents a manifold.
 *
 *  An MFAtlas represents a manifold by way of a list of charts.
 *  Each chart represents a non-empty neighborhood on the manifold. The charts must all have the same domain (a k-dimensional
 *  Euclidean, and range (an n-dimensional Euclidean space). The dimension of the range is the dimension of the manifold
 *  and the dimension of the domain is the dimension of the manifold.
 *
 *  Incompletely computed manifolds are allowed, but in principle every point on the manifold should be contained in the
 *  image of at least one of the charts in the atlas. In addition, if two charts overlap, then all of the points in the
 *  overlap must map to the same point in the embedding space.
 */

struct MFAtlasSt;
typedef struct MFAtlasSt *MFAtlas;
#include <MFContinuationMethod.h>

/*! \fn MFAtlas MFCreateAtlas(MFImplicitMF M,MFErrorHandler e);
 *  \brief Creates an atlas with no charts.
 *
 *  \param M the implicitly defined manifold this atlas describes.
 *  \param e A place to return errors.
 *  \return A new atlas.
 */
MFAtlas MFCreateAtlas(MFImplicitMF,MFErrorHandler);

/*! \fn void MFFreeAtlas(MFAtlas A,MFErrorHandler e);
 *  \brief Releases an atlas. This subtracts one from the reference count, and if the count is zero, frees the storage
 *         associated with the atlas.
 *
 *  \param A The atlas being released.
 *  \param e A place to return errors.
 */
void MFFreeAtlas(MFAtlas,MFErrorHandler);

/*! \fn void MFRefAtlas(MFAtlas A, MFErrorHandler e);
 *  \brief Adds a reference to the atlas.
 *
 *  \param A The Atlas being referenced.
 *  \param e A place to return errors.
 *  \sa ReferenceCounting MFFreeAtlas
 */
void MFRefAtlas(MFAtlas,MFErrorHandler);

/*! \fn int MFAtlasAddChartWithApproxTS(MFAtlas A,MFNVector u,MFNKMatrix Phi,double delta,double R,MFErrorHandler e);
 *  \brief Computes a neighborhood of a point u on the manifold, adds a chart to the list of charts in an atlas and updates
 *         the polyhedra of the charts. An approximation of the tangent space is passed in and used when the new chart is
 *         created. 
 *
 *  \param A The atlas being modified.
 *  \param u The center of the chart to be created and added.
 *  \param Phi An approximation to the tangent space of the new chart.
 *  \param delta The radius of the new chart is chosen so that in the domain of the chart the distance between the linear
 *         approximation of the manifold at the center of the chart and the manifold is less than delta.
 *  \param R A radius to use for the new chart if one can't be computed.
 *  \param e A place to return errors.
 *  \returns TRUE if the operation was sucessful.
 */
int MFAtlasAddChartWithApproxTS(MFAtlas,MFNVector,MFNKMatrix,double,double,MFErrorHandler);

/*! \fn int MFAtlasK(MFAtlas A,MFErrorHandler e);
 *  \brief Returns the dimension of the manifold.
 *
 *  \param A The atlas being queried.
 *  \param e A place to return errors.
 *  \returns The dimension of the manifold.
 */
int MFAtlasK(MFAtlas,MFErrorHandler);

/*! \fn int MFAtlasN(MFAtlas A,MFErrorHandler e);
 *  \brief Returns the dimension of the embedding space of the manifold.
 *
 *  \param A The atlas being queried.
 *  \param e A place to return errors.
 *  \returns The embedding dimension of the manifold.
 */
int MFAtlasN(MFAtlas,MFErrorHandler);

/*! \fn MFImplicitMF MFAtlasMF(MFAtlas A,MFErrorHandler e);
 *  \brief Returns the implicitly defined manifold which the atlas represents.
 *
 *  \param A The atlas being queried.
 *  \param e A place to return errors.
 *  \returns The implicit defined manifold.
 */
MFImplicitMF MFAtlasMF(MFAtlas,MFErrorHandler);

/*! \fn int MFAtlasNumberOfCharts(MFAtlas A,MFErrorHandler e);
 *  \brief Returns the number of charts in the atlas.
 *
 *  \param A The atlas being queried.
 *  \param e A place to return errors.
 *  \returns The number of charts.
 */
int MFAtlasNumberOfCharts(MFAtlas,MFErrorHandler);

/*! \fn double MFAtlasChartRadius(MFAtlas A,int chart,MFErrorHandler e);
 *  \brief Returns the radius of the specified chart.
 *
 *  \param A The atlas being queried.
 *  \param chart The chart being queried (between 0 and MFAtlasNumberOfCharts(A)-1).
 *  \param e A place to return errors.
 *  \returns The radius.
 */
double MFAtlasChartRadius(MFAtlas,int,MFErrorHandler);

/*! \fn MFNVector MFAtlasChartCenter(MFAtlas A,int chart,MFErrorHandler e);
 *  \brief Returns the center of the specified chart.
 *
 *  \param A The atlas being queried.
 *  \param chart The chart being queried (between 0 and MFAtlasNumberOfCharts(A)-1).
 *  \param e A place to return errors.
 *  \returns The center of the chart.
 */
MFNVector MFAtlasChartCenter(MFAtlas,int,MFErrorHandler);

MFNVector MFAtlasCenterOfChart(MFAtlas,int,MFErrorHandler);

/*! \fn MFNKMatrix MFAtlasChartTangentSpace(MFAtlas A,int chart,MFErrorHandler e);
 *  \brief Returns the tangent space of the specified chart.
 *
 *  \param A The atlas being queried.
 *  \param chart The chart being queried (between 0 and MFAtlasNumberOfCharts(A)-1).
 *  \param e A place to return errors.
 *  \returns The tangent space of the chart.
 */
MFNKMatrix MFAtlasChartTangentSpace(MFAtlas,int,MFErrorHandler);

/*! \fn int MFAtlasIsPointInChart(MFAtlas A,int chart,MFKVector s,MFErrorHandler e);
 *  \brief Tests whether a k-dimensional point is in the domain of a particular chart in the atlas.
 *
 *  \param A The atlas being queried.
 *  \param chart The chart being queried (between 0 and MFAtlasNumberOfCharts(A)-1).
 *  \param s The k-dimensional point.
 *  \param e A place to return errors.
 *  \returns TRUE if the point is in the domain.
 */
int MFAtlasIsPointInChart(MFAtlas,int,MFKVector,MFErrorHandler);

int MFAtlasIsChartNearBoundary(MFAtlas,int,MFErrorHandler);

/*! \fn void MFAtlasEvaluateChart(MFAtlas A,int chart,MFKVector s,MFNVector u,MFErrorHandler e);
 *  \brief Evaluates the chart mapping of the specified chart.
 *
 *  \param A The atlas being queried.
 *  \param chart The chart being queried (between 0 and MFAtlasNumberOfCharts(A)-1).
 *  \param s The k-dimensional point.
 *  \param u The image of s.
 *  \param e A place to return errors.
 */
void MFAtlasEvaluateChart(MFAtlas,int,MFKVector,MFNVector,MFErrorHandler);

/*! \fn int MFAtlasNumberOfChartNeighbors(MFAtlas A,int chart,MFErrorHandler e);
 *  \brief Returns the number of charts in the atlas that are neighbors of the specified chart (the number of faces in the 
 *         chart's polyhedron).
 *
 *  \param A The atlas being queried.
 *  \param chart The chart being queried (between 0 and MFAtlasNumberOfCharts(A)-1).
 *  \param e A place to return errors.
 *  \returns The number of neighbors.
 */
int MFAtlasNumberOfChartNeighbors(MFAtlas,int,MFErrorHandler);

/*! \fn int MFAtlasChartNeighbor(MFAtlas A,int chart,int neighbor,MFErrorHandler e);
 *  \brief Returns the index of the chart in the atlas of a neighbor of the specified chart.
 *
 *  \param A The atlas being queried.
 *  \param chart The chart being queried (between 0 and MFAtlasNumberOfCharts(A)-1).
 *  \param neighbor The number of the neighbor being requested (between 0 and MFAtlasNumberOfChartNeighbors(A,chart)-1).
 *  \param e A place to return errors.
 *  \returns The neighbor (between 0 and MFAtlasNumberOfCharts(A)-1).
 */
int MFAtlasChartNeighbor(MFAtlas,int,int,MFErrorHandler);

/*! \fn int MFAtlasNumberOfChartsWithBoundary(MFAtlas A,MFErrorHandler e);
 *  \brief Returns the number of charts that are currently on the  boundary of the union of the chart images.
 *
 *  \param A The atlas being queried.
 *  \param e A place to return errors.
 *  \returns The number of boundary charts.
 */
int MFAtlasNumberOfChartsWithBoundary(MFAtlas,MFErrorHandler);

/*! \fn int MFAtlasChartWithBoundary(MFAtlas A,int chart,MFErrorHandler e);
 *  \brief Returns the index of a chart in the boundary list of the atlas.
 *
 *  \param A The atlas being queried.
 *  \param chart The chart being queried (between 0 and MFAtlasNumberOfChartsWithBoundary(A)-1).
 *  \param e A place to return errors.
 *  \returns The boundary chart.
 */
int MFAtlasChartWithBoundary(MFAtlas,int,MFErrorHandler);

/*! \fn int MFAtlasPointOnBoundaryInsideRegion(MFAtlas A,MFNRegion Omega,MFNVector u,MFNKMatrix *Phi,double *R,MFErrorHandler e);
 *  \brief Creates a point very near the current boundary of the manifold and a neighborhood of the point.
 *
 *  \param A The atlas being queried.
 *  \param Omega The point is required to lie within Omega.
 *  \param u The center of the new neighborhood (if the return code is TRUE).
 *  \param Phi If the return code is TRUE, the tangent space at the center.
 *  \param R If the return code is TRUE, an estimate of the radius of the new neigborhood.
 *  \param e A place to return errors.
 *  \returns TRUE if a boundary point exists in Omega.
 */
int MFAtlasPointOnBoundaryInsideRegion(MFAtlas,MFNRegion,MFNVector,MFNKMatrix*,double*,MFErrorHandler);

/*! \fn MFAtlas MFReadAtlas(FILE *fid,MFErrorHandler e);
 *  \brief Reads an Atlas from a file.
 *
 *  \param fid the file containing an atlas.
 *  \param e A place to return errors.
 *  \returns The atlas.
 */
MFAtlas MFReadAtlas(FILE*,MFErrorHandler);

/*! \fn void MFWriteAtlas(FILE *fid,MFAtlas A,MFErrorHandler e);
 *  \brief Writes an Atlas to a file.
 *
 *  \param fid the file to write to an atlas.
 *  \param A The atlas.
 *  \param e A place to return errors.
 */
void MFWriteAtlas(FILE*,MFAtlas,MFErrorHandler);

/*! \fn void MFWriteChartCenters(FILE *fid,MFAtlas A,MFErrorHandler e);
 *  \brief Writes the centers, tangent space bases, and radii of the charts in an Atlas to a file.
 *
 *  \param fid the file to write the centers into.
 *  \param A The atlas.
 *  \param e A place to return errors.
 */
void MFWriteChartCenters(FILE*,MFAtlas,MFErrorHandler);

/*! \fn int MFReadChartCenters(FILE *fid,MFNVector** centerList,MFNKMatrix** PhiList ,double **RList,MFErrorHandler e);
 *  \brief Reads the centers, tangent space bases, and radii of the charts in a file.
 *
 *  \param fid the file containing an atlas.
 *  \param centerList A pointer to an array of MFNVectors for the centers
 *  \param PhiList A pointer to an array of MFNKMatrices for the tangent space bases.
 *  \param RList A pointer to an array of double for the radii.
 *  \param e A place to return errors.
 *  \returns The number of centers, tangent space bases and radii read.
 */
int MFReadChartCenters(FILE*,MFNVector**,MFNKMatrix**,double**,MFErrorHandler);

/*! \fn void MFAtlasSetEpsilon(MFAtlas A,double epsilon,MFErrorHandler e);
 *  \brief Sets the internal parameter epsilon, which is the tolerance allowed between the tangent space and the manifold.
 *
 *  \param A The atlas.
 *  \param epsilon The tolerance.
 *  \param e A place to return errors.
 */
void MFAtlasSetEpsilon(MFAtlas,double,MFErrorHandler);

/*! \fn void MFAtlasSetRMin(MFAtlas A,double Rmin,MFErrorHandler e);
 *  \brief Sets the internal parameter Rmin, which is the minimum radius allowed for a chart.
 *
 *  \param A The atlas.
 *  \param Rmin The minimum radius.
 *  \param e A place to return errors.
 */
void MFAtlasSetRMin(MFAtlas,double,MFErrorHandler);

/*! \fn void MFAtlasSetDotMin(MFAtlas A,double dotmin,MFErrorHandler e);
 *  \brief Sets the internal parameter dotmin, which is the smallest dot product between neighboring tangent spaces. (Generalized
 *         to more than one dimensional manifolds. The dot product for a flat manifold will be 1.
 *
 *  \param A The atlas.
 *  \param dotmin The minimum allowed dot product.
 *  \param e A place to return errors.
 */
void MFAtlasSetDotMin(MFAtlas,double,MFErrorHandler);

/*! \fn void MFAtlasSetVerbose(MFAtlas A,int verbose,MFErrorHandler e);
 *  \brief Controls the output to stdout from operations on the atlas. A value of 0 is quiet.
 *
 *  \param A The atlas.
 *  \param verbose The level of verbosity.
 *  \param e A place to return errors.
 */
void MFAtlasSetVerbose(MFAtlas,int,MFErrorHandler);

void MFAtlasReduceChartRadius(MFAtlas,int,double,MFErrorHandler);

struct MFChartSt;

/*! \fn MFAtlas MFComputeAtlas(MFContinuationMethod algorithm,MFImplicitMF M,MFNRegion Omega,MFNVector u0,MFErrorHandler e);
 *  \brief Creates an atlas which covers the connected component of the implicitly defined manifold M which is inside
 *   the region Omega, and contains the initial point u0.
 *
 *  \param algorithm An algorithm to use for the computation.
 *  \param M An implicitly defined manifold.
 *  \param Omega The region.
 *  \param u0 A point on M.
 *  \param e A place to return errors.
 *  \returns A new atlas.
 */
MFAtlas MFComputeAtlas(MFContinuationMethod,MFImplicitMF,MFNRegion,MFNVector,MFErrorHandler);

/*! \fn MFAtlas MFComputeAtlasWithTangent(MFContinuationMethod algorithm,MFImplicitMF M,MFNRegion Omega,MFNVector u0,MFNKMatrix Phi0,MFErrorHandler e);
 *  \brief Creates an atlas which covers the connected component of the implicitly defined manifold M which is inside
 *   the region Omega, and contains the initial point u0.
 *
 *  \param algorithm An algorithm to use for the computation.
 *  \param M An implicitly defined manifold.
 *  \param Omega The region.
 *  \param u0 A point on M.
 *  \param Phi0 The tangent space of M at u.
 *  \param e A place to return errors.
 *  \returns A new atlas.
 */
MFAtlas MFComputeAtlasWithTangent(MFContinuationMethod,MFImplicitMF,MFNRegion,MFNVector,MFNKMatrix,MFErrorHandler);

/*! \fn MFAtlas MFComputeAtlasMultiple(MFContinuationMethod algorithm,MFImplicitMF M,MFNRegion Omega,int n,MFNVector *u0,MFErrorHandler e);
 *  \brief Creates an atlas which covers the connected component of the implicitly defined manifold M which is inside
 *   the region Omega, and contains the n initial points in the array u0[].
 *
 *  \param algorithm An algorithm to use for the computation.
 *  \param M An implicitly defined manifold.
 *  \param Omega The region.
 *  \param n The number of initial points on M.
 *  \param u0 A list of points on M.
 *  \param e A place to return errors.
 *  \returns A new atlas.
 */
MFAtlas MFComputeAtlasMultiple(MFContinuationMethod,MFImplicitMF,MFNRegion,int,MFNVector*,MFErrorHandler);

/*! \fn MFAtlas MFComputeAtlasMultipleWithTangents(MFContinuationMethod algorithm,MFImplicitMF M,MFNRegion Omega,int n,MFNVector *u0,MFNKMatrix *Phi0,MFErrorHandler e);
 *  \brief Creates an atlas which covers the connected component of the implicitly defined manifold M which is inside
 *   the region Omega, and contains the n initial points in the array u0[].
 *
 *  \param algorithm An algorithm to use for the computation.
 *  \param M An implicitly defined manifold.
 *  \param Omega The region.
 *  \param n The number of initial points on M.
 *  \param u0 A list of points on M.
 *  \param Phi0 A list of tangent spaces of M at the initial points.
 *  \param e A place to return errors.
 *  \returns A new atlas.
 */
MFAtlas MFComputeAtlasMultipleWithTangents(MFContinuationMethod,MFImplicitMF,MFNRegion,int,MFNVector*,MFNKMatrix*,MFErrorHandler);

/*! \fn MFAtlas MFAnimateAtlas(MFImplicitMF M, MFNRegion Omega, MFNVector u0, double epsilon, int kmax, int verbose, char *stub, int DrawSimp, int DrawAfter, int DrawEvery,MFErrorHandler e);
 *  \brief Creates an animation of the computation of an atlas. That is, it writes images of the Atlas after adding a chart.
 *         The viewpoint etc are read from a .view file.
 *
 *  \param M An implicitly defined manifold.
 *  \param Omega The region.
 *  \param u0 An initial point on M.
 *  \param epsilon The tolerance for the distance between the manifold and a linear approximation at the chart centers. Used to estimate R.
 *  \param kmax The largest number of charts to add.
 *  \param verbose Controls the verbosity on stdout. A value of 0 is quiet.
 *  \param stub A prefix to use for the files.
 *  \param DrawSimp TRUE if the dual is to be drawn prefix to use for the files.
 *  \param DrawAfter skips the first DrawAfter-1 frames in the animation.
 *  \param DrawEvery Creates an image every time DrawEvery charts have been added.
 *  \param e A place to return errors.
 *  \returns A new atlas.
 */
MFAtlas MFAnimateAtlas(MFImplicitMF,MFNRegion,MFNVector,double,int,int,char*,int,int,int,MFErrorHandler);

MFAtlas MFComputeInvariantMap(MFImplicitMF,MFNRegion,MFNVector,MFNKMatrix,double,double,int,int,MFErrorHandler);

double MFVolumeOfChart(struct MFChartSt*,MFNRegion,MFErrorHandler);

/*! \fn double MFVolumeOfAtlas(MFAtlas A,MFNRegion Omega,MFErrorHandler e);
 *  \brief Computes the volume of a manifold. This is approximated by adding the volumes of the convex chart polyhedra.
 *
 *  \param A An atlas representing the manifold.
 *  \param Omega The calculation will include the part of the manifold in this region.
 *  \param e A place to return errors.
 *  \returns The volume.
 */
double MFVolumeOfAtlas(MFAtlas,MFNRegion,MFErrorHandler);

/*! \fn int MFAtlasAddChartWithAll(MFAtlas A,MFNVector u,MFNKMatrix Phi,double R,MFErrorHandler e);
 *  \brief Adds a chart to the Atlas.A The chart is created with the center, tangents space and radius given.
 *
 *  \param A The atlas being modified.
 *  \param u The center of the chart to be created and added.
 *  \param Phi A basis for the tangent space of the manifold at u.
 *  \param R The radius of the domain of the chart.
 *  \param e A place to return errors.
 *  \returns TRUE if the operation was sucessful.
 */
int MFAtlasAddChartWithAll(MFAtlas,MFNVector,MFNKMatrix,double,MFErrorHandler);

/*! \fn int MFAtlasAddChart(MFAtlas A,MFNVector u,MFErrorHandler e);
 *  \brief Computes a neighborhood of a point u on the manifold, adds a chart to the list of charts in an atlas and updates
 *         the polyhedra of the charts.
 *
 *  \param A The atlas being modified.
 *  \param u The center of the chart to be created and added.
 *  \param e A place to return errors.
 *  \returns TRUE if the operation was sucessful.
 */
int MFAtlasAddChart(MFAtlas,MFChart,MFErrorHandler);

/*! \fn int MFAtlasPointOnBoundaryWOProject(MFAtlas A,MFNRegion Omega,MFKVector s,MFErrorHandler e);
 *  \brief Finds a point near the boundary of the atlas, on an existing chart, but not projected onto the manifold.
 *
 *  \param A The atlas being modified.
 *  \param Omega A region within which the boundary point must lie.
 *  \param s The point in the chart's domain which gives the required point.
 *  \param e A place to return errors.
 *  \returns The number of the chart.
 */
int MFAtlasPointOnBoundaryWOProject(MFAtlas,MFNRegion,MFKVector,MFErrorHandler);

int MFComputePointsOn3dEdge(double,double**,double,double,double,double,double,double,MFErrorHandler);
int MFComputePointsOn3dPolygon(double,int,double*,int,double*,double*,double**,double**,MFErrorHandler);

/*! \fn MFNVector MFAtlasGetPointOnBoundaryChart(MFAtlas A,MFNRegion Omega,int chart,double t0,MFErrorHandler e);
 *  \brief Refines a point near the boundary by bisection.
 *
 *  \param A The atlas being modified.
 *  \param Omega A region within which the boundary point must lie.
 *  \param chart A boundary chart.
 *  \param t0 A distance which generates a point outside Omega.
 *  \param e A place to return errors.
 *  \returns The point on the manifold.
 */
MFNVector MFAtlasGetPointOnBoundaryChart(MFAtlas,MFNRegion,int,double,MFErrorHandler);

/*! \fn int MFAtlasGetSingularChartWithBoundary(MFAtlas A,MFNRegion Omega,MFErrorHandler e);
 *  \brief Get a singular chart that lies on the boundary.
 *
 *  \param A The atlas.
 *  \param Omega A region within which the boundary point must lie.
 *  \param e A place to return errors.
 *  \returns The number of the chart.
 */
int MFAtlasGetSingularChartWithBoundary(MFAtlas,MFNRegion,MFErrorHandler);

/*! \fn void MFExtendAtlas(MFAtlas A,MFContinuationMethod algorithm,MFImplicitMF M,MFNRegion Omega,MFNVector u0,MFErrorHandler e);
 *  \brief Extend a manifold.
 *
 *  \param A The atlas.
 *  \param algorithm A continuation algorithm.
 *  \param M The implicitly defined manifold.
 *  \param Omega A region for the computation.
 *  \param u0 An initial point on the manifold.
 *  \param e A place to return errors.
 */
void MFExtendAtlas(MFAtlas,MFContinuationMethod,MFImplicitMF,MFNRegion,MFNVector,MFErrorHandler);

/*! \fn void MFExtendAtlasMultiple(MFAtlas A,MFContinuationMethod algorithm,MFImplicitMF M,MFNRegion Omega,int n,MFNVector *u0,MFErrorHandler e);
 *  \brief Extend a manifold starting with a list of initial points.
 *
 *  \param A The atlas.
 *  \param algorithm A continuation algorithm.
 *  \param M The implicitly defined manifold.
 *  \param Omega A region for the computation.
 *  \param n The number of initial points in the array u0.
 *  \param u0 An array of initial points on the manifold.
 *  \param e A place to return errors.
 */
void MFExtendAtlasMultiple(MFAtlas,MFContinuationMethod,MFImplicitMF,MFNRegion,int,MFNVector*,MFErrorHandler);

/*! \fn void MFExtendAtlasWithTangent(MFAtlas A,MFContinuationMethod algorithm,MFImplicitMF M,MFNRegion Omega,MFNVector u0,MFNKMatrix Tan0,MFErrorHandler e);
 *  \brief Extend a manifold, with one initial solution and tangent space.
 *
 *  \param A The atlas.
 *  \param algorithm A continuation algorithm.
 *  \param M The implicitly defined manifold.
 *  \param Omega A region for the computation.
 *  \param u0 An initial point on the manifold.
 *  \param Tan0 A basis for the tangent space of the manifold at the initial point.
 *  \param e A place to return errors.
 */
void MFExtendAtlasWithTangent(MFAtlas,MFContinuationMethod,MFImplicitMF,MFNRegion,MFNVector,MFNKMatrix,MFErrorHandler);

/*! \fn void MFExtendAtlasMultipleWithTangents(MFAtlas A,MFContinuationMethod algorithm,MFImplicitMF M,MFNRegion Omega,int n,MFNVector *u0,MFNKMatrix *Tan0,MFErrorHandler e);
 *  \brief Extend a manifold starting with a list of initial points and tangent spaces.
 *
 *  \param A The atlas.
 *  \param algorithm A continuation algorithm.
 *  \param M The implicitly defined manifold.
 *  \param Omega A region for the computation.
 *  \param n The number of initial points in the arrays u0 and Tan0.
 *  \param u0 A list of initial points on the manifold.
 *  \param Tan0 A list of bases for the tangent spaces of the manifold at the initial points.
 *  \param e A place to return errors.
 */
void MFExtendAtlasMultipleWithTangents(MFAtlas,MFContinuationMethod,MFImplicitMF,MFNRegion,int,MFNVector*,MFNKMatrix*,MFErrorHandler);

/*! \fn void MFCloseAtlas(MFContinuationMethod algorithm,MFAtlas A,MFErrorHandler e);
 *  \brief Finish any processing for the atlas.
 *
 *  \param algorithm A continuation algorithm.
 *  \param A The atlas.
 *  \param e A place to return errors.
 */
void MFCloseAtlas(MFContinuationMethod,MFAtlas,MFErrorHandler);

/*! \fn void MFFlushAtlas(MFContinuationMethod algorithm,MFAtlas A,MFErrorHandler e);
 *  \brief Finish any i/o for the atlas.
 *
 *  \param algorithm A continuation algorithm.
 *  \param A The atlas.
 *  \param e A place to return errors.
 */
void MFFlushAtlas(MFContinuationMethod,MFAtlas,MFErrorHandler);

/*! \fn int MFAtlasIsChartSingular(MFAtlas A,int chart,MFErrorHandler e);
 *  \brief Queries the specified chart to see if it is marked singular.
 *
 *  \param A The atlas being queried.
 *  \param chart The chart being queried (between 0 and MFAtlasNumberOfCharts(A)-1).
 *  \param e A place to return errors.
 *  \returns TRUE If the chart being queried is marked as singular.
 */
int MFAtlasIsChartSingular(MFAtlas,int,MFErrorHandler);

/*! \fn double MFAtlasChartSuggestedRadius(MFAtlas A,int chart,MFErrorHandler e);
 *  \brief Returns the suggested radius of the specified chart.
 *
 *  \param A The atlas being queried.
 *  \param chart The chart being queried (between 0 and MFAtlasNumberOfCharts(A)-1).
 *  \param e A place to return errors.
 *  \returns The suggested chart radius.
 */
double MFAtlasChartSuggestedRadius(MFAtlas,int,MFErrorHandler);

/*! \fn void MFAtlasAddClipF(MFAtlas A,double (*side)(MFNVector),MFErrorHandler e);
 *  \brief Add a clipping plane to the atlas. This serves as a way to clean up plots. The chart polyhedra are clipped against
 *         the function "side". Linear interpolation is used on the value of "side" at the vertices, and the positive part
 *         is kept.
 *
 *  \param A The atlas.
 *  \param side A function giving the distance from the clipping surface.
 *  \param e A place to return errors.
 */
void MFAtlasAddClipF(MFAtlas,double (*)(MFNVector,MFErrorHandler),MFErrorHandler);

/*! \fn void MFAtlasClearClipF(MFAtlas A,MFErrorHandler e);
 *  \brief Removes all clipping planes.
 *
 *  \param A The atlas.
 *  \param e A place to return errors.
 */
void MFAtlasClearClipF(MFAtlas,MFErrorHandler);

/*! \fn void MFAtlasSetExpFactor(MFAtlas A,double factor,MFErrorHandler e);
 *  \brief Sets the parameter controlling the output to the plotfile. Each vertex of a polygon is scaled by this factor, so that
 *         any gaps between polyhedral faces are closed.
 *
 *  \param A The atlas.
 *  \param factor The expansion factor.
 *  \param e A place to return errors.
 */
void MFAtlasSetExpFactor(MFAtlas,double,MFErrorHandler);

/*! \fn void MFAtlasGetExpFactor(MFAtlas A,MFErrorHandler e);
 *  \brief Gets the current value of the expansion factor controlling the output to the plotfile.
 *
 *  \param A The atlas.
 *  \param e A place to return errors.
 *  \returns The expansion factor.
 */
double MFAtlasGetExpFactor(MFAtlas,MFErrorHandler);


/*! \fn void MFAtlasNSpace(MFAtlas A,MFErrorHandler e);
 *  \brief Gets the NSpace in which the Atlas is embedded.
 *
 *  \param A The atlas.
 *  \param e A place to return errors.
 *  \returns The NSpace.
 */
MFNSpace MFAtlasNSpace(MFAtlas,MFErrorHandler);

/*! \fn MFChart MFAtlasChart(MFAtlas A, int chart, MFErrorHandler e);
 *  \brief Gets the specified chart.
 *
 *  \param A The atlas.
 *  \param i The number of the chart.
 *  \param e A place to return errors.
 *  \returns The chart.
 */
MFChart MFAtlasChart(MFAtlas A, int chart, MFErrorHandler e);

void MFAtlasClosePlotfile(MFAtlas,MFErrorHandler);
void MFAtlasCloseCenterfile(MFAtlas,MFErrorHandler);

void MFUpdateNeighborsReferenceMarks(MFAtlas,int,int,MFErrorHandler);
int MFAtlasAddChartToList(MFAtlas,MFChart,MFErrorHandler);

void MFAtlasReleaseSingularCharts(MFAtlas,int,int,MFErrorHandler);

int MFAtlasAddChartWithCenter(MFAtlas,MFNVector,MFErrorHandler);
int MFAtlasAddChart(MFAtlas,MFChart,MFErrorHandler);
void MFAtlasReduceChartRadius(MFAtlas,int,double,MFErrorHandler);

#ifdef __cplusplus
}
#endif

#endif

/*! @} */

/*! @} */
