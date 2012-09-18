#ifndef __IMFEXPANSION_H__
#define __IMFEXPANSION_H__
#include <MFAtlas.h>
#include <MFErrorHandler.h>

#ifdef __cplusplus
 extern "C" {
#endif

/*! \defgroup InvariantManifolds */
/*! \defgroup IMFExpansion */

/*! \addtogroup InvariantManifolds
 *  @{
 */

/*! \addtogroup IMFExpansion
 *  @{
 */

/*! \class IMFExpansion IMFExpansion.h IMFExpansion.h
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
struct IMFExpansionSt;
typedef struct IMFExpansionSt *IMFExpansion;

/*!  \fn IMFExpansion IMFCreateExpansion(int n, int k,MFErrorHandler e);
 *   \brief Creates a Taylor series expansion to third order for a mapping from IR^k to IR^n. This is the default
 *            constructor, and returns an expansion of the function which is always zero.
 *
 *   \param n The dimension of the range of the expansion.
 *   \param k The dimension of the domain of the expansion.
 *   \param e A place to return errors.
 *   \returns A third order expansion of the zero valued funciton.
 */
IMFExpansion IMFCreateExpansion(int,int,MFErrorHandler);

/*! \fn void IMFFreeExpansion(IMFExpansion E,MFErrorHandler e);
 *  \brief Frees a reference to an expansion, and deletes the expansion if there are no references left.
 *
 *  \param E The expansion being unreferenced.
 *  \param e A place to return errors.
 *  \sa ReferenceCounting MFRefExpansion()
 */
void IMFFreeExpansion(IMFExpansion,MFErrorHandler);

/*! \fn void IMFRefExpansion(IMFExpansion E,MFErrorHandler e);
 *  \brief Adds a reference to a expansion.
 *
 *  \param E The expansion being referenced.
 *  \param e A place to return errors.
 *  \sa ReferenceCounting MFFreeExpansion()
 */
void IMFRefExpansion(IMFExpansion,MFErrorHandler);

/*! \fn void IMFEvaluateExpansion(IMFExpansion E,MFKVector s,MFNVector u,MFErrorHandler e);
 *  \brief Evaluates an expansion.
 *
 *  \param E An expansion.
 *  \param s A point in the domain of the expansion.
 *  \param u Provided by the user, a place to put the result.
 *   \param e A place to return errors.
 */
void IMFEvaluateExpansion(IMFExpansion,MFKVector,MFNVector,MFErrorHandler);

/*! \fn void IMFEvaluateExpansionDirectionalDerivative(IMFExpansion E,MFKVector s,MFKVector ds,MFNVector du,MFErrorHandler e);
 *  \brief Evaluates an expansion.
 *
 *  \param E An expansion.
 *  \param s A point in the domain of the expansion.
 *  \param ds A direction in which to find the partial derivative.
 *  \param du Provided by the user, a place to put the result.
 *   \param e A place to return errors.
 */
void IMFEvaluateExpansionDirectionalDerivative(IMFExpansion,MFKVector,MFKVector,MFNVector,MFErrorHandler);

/*! \fn void IMFEvaluateExpansionSecondDirectionalDerivative(IMFExpansion E,MFKVector s,MFKVector ds0,MFKVector ds1,MFNVector ddu,MFErrorHandler e);
 *  \brief Evaluates an expansion.
 *
 *  \param E An expansion.
 *  \param s A point in the domain of the expansion.
 *  \param ds0 The first direction.
 *  \param ds1 The second direction.
 *  \param ddu Provided by the user, a place to put the result.
 *  \param e A place to return errors.
 */
void IMFEvaluateExpansionSecondDirectionalDerivative(IMFExpansion,MFKVector,MFKVector,MFKVector,MFNVector,MFErrorHandler);

/*! \fn void IMFEvaluateExpansionDerivative(IMFExpansion E,MFKVector s,int d,MFNVector du,MFErrorHandler e);
 *  \brief Evaluates the derivative of an expansion along a coordinate direction.
 *
 *  \param E An expansion.
 *  \param s A point in the domain of the expansion.
 *  \param d A coordinate axis in which to find the partial derivative.
 *  \param du Provided by the user, a place to put the result.
 *  \param e A place to return errors.
 */
void IMFEvaluateExpansionDerivative(IMFExpansion,MFKVector,int,MFNVector,MFErrorHandler);

/*! \fn void IMFEvaluateExpansionSecondDerivative(IMFExpansion E,MFKVector s,int d0,int d1,MFNVector ddu,MFErrorHandler e);
 *  \brief Evaluates the second derivative of an expansion along coordinate directions.
 *
 *  \param E An expansion.
 *  \param s A point in the domain of the expansion.
 *  \param d0 A coordinate axis in which to find the partial derivative.
 *  \param d1 A second coordinate axis in which to find the partial derivative.
 *  \param ddu Provided by the user, a place to put the result.
 *  \param e A place to return errors.
 */
void IMFEvaluateExpansionSecondDerivative(IMFExpansion,MFKVector,int,int,MFNVector,MFErrorHandler);

/*! \fn void IMFEvaluateExpansionThirdDerivative(IMFExpansion E,MFKVector s,int d0,int d1,int d2,MFNVector dddu,MFErrorHandler e);
 *  \brief Evaluates the third derivative of an expansion along coordinate directions.
 *
 *  \param E An expansion.
 *  \param s A point in the domain of the expansion.
 *  \param d0 A coordinate axis in which to find the partial derivative.
 *  \param d1 A second coordinate axis in which to find the partial derivative.
 *  \param d2 A third coordinate axis in which to find the partial derivative.
 *  \param dddu Provided by the user, a place to put the result.
 *  \param e A place to return errors.
 */
void IMFEvaluateExpansionThirdDerivative(IMFExpansion,MFKVector,int,int,int,MFNVector,MFErrorHandler);

void IMFEvaluateExpansionDirectionalDerivativeE(IMFExpansion,double*,double*,double*,MFErrorHandler);

/*! \fn int IMFExpansionN(IMFExpansion E,MFErrorHandler e);
 *  \brief Returns the dimension of the range of an expansion.
 *
 *  \param E An expansion.
 *  \param e A place to return errors.
 *  \returns The dimension of the range of the expansion.
 */
int IMFExpansionN(IMFExpansion,MFErrorHandler);

/*! \fn int IMFExpansionK(IMFExpansion E,MFErrorHandler e);
 *  \brief Returns the dimension of the domain of an expansion.
 *
 *  \param E An expansion.
 *  \param e A place to return errors.
 *  \returns The dimension of the domain of the expansion.
 */
int IMFExpansionK(IMFExpansion,MFErrorHandler);

/*! \fn int IMFExpansionOrder(IMFExpansion E,MFErrorHandler e);
 *  \brief Returns the order of an expansion (the degree of the highest non-zero term in the Taylor series).
 *
 *  \param E An expansion.
 *  \param e A place to return errors.
 *  \returns The order of the expansion.
 */
int IMFExpansionOrder(IMFExpansion,MFErrorHandler);

double *IMFExpansionU(IMFExpansion,MFErrorHandler);
double *IMFExpansionDu(IMFExpansion,MFErrorHandler);
double *IMFExpansionDDu(IMFExpansion,MFErrorHandler);
double *IMFExpansionDDDu(IMFExpansion,MFErrorHandler);

void IMFExpnSetDerivatives(IMFExpansion,double*,double*,double*,double*,MFErrorHandler);
char *IMFC(double,char*,int*);

/*! \fn void IMFPrintExpansion(FILE *fid,IMFExpansion E,MFErrorHandler e);
 *  \brief Prints a \"pretty\" version of the expansion.
 *
 *  \param fid The file.
 *  \param E The expansion.
 *  \param e A place to return errors.
 */
void IMFPrintExpansion(FILE*,IMFExpansion,MFErrorHandler);

/*! \fn IMFExpansion IMFCloneExpansion(IMFExpansion E,MFErrorHandler e);
 *  \brief Creates a deep copy of an expansion.
 *
 *  \param E The expansion.
 *  \param e A place to return errors.
 *  \returns A clone.
 */
IMFExpansion IMFCloneExpansion(IMFExpansion,MFErrorHandler);

/*! \fn IMFExpansion IMFRectify(IMFExpansion E,MFErrorHandler e);
 *  \brief Changes the metric of the expansion to the identity. 
 *
 *  \param E The expansion.
 *  \param e A place to return errors.
 *  \returns A new expansion which has been rectified (the metric changed to the identity).
 */
IMFExpansion IMFRectifyExpansion(IMFExpansion,MFErrorHandler);

struct IMFFlowSt;

/*! \fn IMFExpansion IMFInflateExpansionWithFlow(IMFExpansion E,IMFFlow f, MFKVector p,MFErrorHandler e);
 *  \brief Creates the product of an expansion with a trjectory of a flow.
 *
 *  \param E The starting expansion.
 *  \param f The flow. (At the origin of the expansion should be transverse to the surface defined by the expansion).
 *  \param p The flow's parameter values for the trajectory.
 *  \param e A place to return errors.
 *  \returns A new expansion.
 */
IMFExpansion IMFInflateExpansionWithFlow(IMFExpansion,struct IMFFlowSt*,MFKVector,MFErrorHandler);

/*! \fn double IMFExpansionR(IMFExpansion E, double epsilon,MFErrorHandler e);
 *  \brief Estimates a ball such that within the ball the second order terms are bounded in absolute value by epsilon.
 *
 *  \param E The expansion.
 *  \param epsilon The required bound.
 *  \param e A place to return errors.
 *  \returns The radius of the ball.
 */
double IMFExpansionR(IMFExpansion,double,MFErrorHandler);

/*! \fn MFNKMatrix IMFExpansionTS(IMFExpansion E,MFErrorHandler e);
 *  \brief Returns a new MFNKMatrix with an orthonormal basis for the tangent space at the origin of the expansion.
 *
 *  \param E The expansion.
 *  \param e A place to return errors.
 *  \returns The basis in an MFNKMatrix.
 */
MFNKMatrix IMFExpansionTS(IMFExpansion,MFErrorHandler);

double *IMFExpansionData(IMFExpansion,MFErrorHandler);
int IMFExpansionDataLn(IMFExpansion,MFErrorHandler);

void IMFExpansionSetDerivatives(IMFExpansion,double*,double*,double*,double*,MFErrorHandler);

/*! \fn void IMFWriteExpansion(FILE* fid,IMFExpansion E,MFErrorHandler e);
 *  \brief Writes an expansion to a file.
 *
 *  \param fid The file to write to.
 *  \param E The expansion being queried.
 *  \param e A place to return errors.
 */
void IMFWriteExpansion(FILE* fid,IMFExpansion E,MFErrorHandler e);

/*! \fn IMFExpansion IMFReadExpansion(FILE* fid,MFErrorHandler e);
 *  \brief Reads an expansion from a file.
 *
 *  \param fid The file to read from.
 *  \param e A place to return errors.
 *  \returns The expansion.
 */
IMFExpansion IMFReadExpansion(FILE* fid,MFErrorHandler e);

/*! @} */

/*! @} */

#ifdef __cplusplus
}
#endif

#endif
