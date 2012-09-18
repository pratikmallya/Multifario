/* 
    @(#)MFImplicitMF.h	1.21
    03/07/24 11:20:52
   
    PROGRAM NAME:  Manifold
   
    (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
    CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
   
    Please refer to the LICENSE file in the top directory

*/
/*      author: Mike Henderson mhender@watson.ibm.com */
/*      date:   February 19, 1999                     */
/*              October 6, 2004   Added ProjectFromCenter */

#ifndef __MFIMPLICITMF_H__
#define __MFIMPLICITMF_H__
#include <MFBase.h>

#include <MFKVector.h>
#include <MFNVector.h>
#include <MFNKMatrix.h>
#include <MFNSpace.h>
#include <IMFFlow.h>
#include <MFErrorHandler.h>

/*! \defgroup MFImplicitMF */

/*! \defgroup MFPlaneMF */
/*! \defgroup MFNSpaceMF */
/*! \defgroup MFEdgeIn3SpaceMF */
/*! \defgroup MFPolygonIn3SpaceMF */
/*! \defgroup MFAlgebraicMF */
/*! \defgroup MFCircleMF */
/*! \defgroup MFFlatMF */
/*! \defgroup MFTorusMF */
/*! \defgroup MFTPBVPMF */
/*! \defgroup MFSphereMF */
/*! \defgroup MFCoupledPendula */
/*! \defgroup MFSwallowTail */
/*! \defgroup MFForcedOscillator */
/*! \defgroup MFSpline */

/*! \addtogroup MFImplicitMF
 *  @{
 */
/*! \class MFImplicitMF 
 *  \brief A manifold, defined by an equation F(u)=0.
 *
 *  An MFImplicitMF represents a manifold which is described as the inverse image of zero. That's just a fancy
 *  way of saying that points "u" on the manifold satisfy an equation F(u)=0.
 *
 *  MFImplicitMF is a base class, so constructors for particular types of IMF's have to be called. These register
 *  routines that perform the operations that MFImplicitMF requires.
 */

struct MFImplicitMFSt;
typedef struct MFImplicitMFSt *MFImplicitMF;

#ifdef __cplusplus
 extern "C" {
#endif

/*! \addtogroup MFSpline
 *  @{
 */

/*! \fn MFImplicitMF MFIMFCreatePeriodicSpline(int npts, int n, double **x, MFErrorHandler e);
 *  \brief Create a manifold which is a closed curve defined by a spline
 *
 *  \param npts The number of knots on the spline. The first point is not duplicated as the last.
 *  \param n    The dimension of the space in which the spline lives.
 *  \param x    An array of arrays of the data defining the n splines. Dimensions are x[n][npts].
 *  \returns    An implicitly defined manifold.
 */
MFImplicitMF MFIMFCreatePeriodicSpline(int npts, int n, double **x, MFErrorHandler e);

/*! \fn MFImplicitMF MFIMFCreateSpline(int npts, int n, double **x, MFErrorHandler e)
 *  \brief Create a manifold which is a spline curve
 *
 *  \param npts The number of knots on the spline.
 *  \param n    The dimension of the space in which the spline lives.
 *  \param x    An array of arrays of the data defining the n splines. Dimensions are x[n][npts].
 *  \returns    An implicitly defined manifold.
 */

/*! @} */


/*! \addtogroup MFAlgebraicMF
 *  @{
 */

/*! \fn MFImplicitMF MFIMFCreateAlgebraicExpressionWithRadius(char *vars,char *expression, double R,MFErrorHandler e);
 *  \brief Create an algebraic system of nonlinear equations from an expression in a string. This uses the ExpCmp
 *         utility which is pert of multifario.
 *
 *  \param vars A string containing a list of the variables used in the expression. e.g. "[x,y,z]".
 *  \param expression A string containing an expression for the function value in terms of the variables in the first
 *                          argeument. e.g. "[sin(x),1+z*y]". 
 *  \param R The radius.
 *  \param e A place to return errors.
 *  \returns An implicitly defined manifold.
 */
MFImplicitMF MFIMFCreateAlgebraicExpressionWithRadius(char*,char*,double,MFErrorHandler);

/*! \fn MFImplicitMF MFIMFCreateAlgebraicExpression(char *vars,char *expression,MFErrorHandler e);
 *  \brief Create an algebraic system of nonlinear equations from an expression in a string. This uses the ExpCmp
 *         utility which is pert of multifario.
 *
 *  \param vars A string containing a list of the variables used in the expression. e.g. "[x,y,z]".
 *  \param expression A string containing an expression for the function value in terms of the variables in the first
 *                          argeument. e.g. "[sin(x),1+z*y]". 
 *  \param e A place to return errors.
 *  \returns An implicitly defined manifold.
 */
MFImplicitMF MFIMFCreateAlgebraicExpression(char*,char*,MFErrorHandler);

/*! \fn MFImplicitMF MFIMFCreateAlgebraicSubroutine(int n,int k,void (*F)(int*,double*,int*,double*,void*,MFErrorHandler), void (*dF)(int*,double*,int*,double*,void*,MFErrorHandler)), void (*ddF)(int*,double*,int*,double*,void*,MFErrorHandler),void *data,MFErrorHandler e);
 *  \brief This constructor creates a nonlinear set of equations defined by the subroutine F.
 *
 *  \param n The embedding dimension (the dimension of the range of F).
 *  \param k The dimension of the solution manifold (the dimension of the range of F is n+k).
 *  \param F A subroutine that evaluates F.
 *  \param dF A subroutine that evaluates the derivative of F (if NULL differenceing will be used).
 *  \param ddF A subroutine that evaluates the second derivative of F (if NULL differenceing will be used).
 *  \param data A pointer which is passed to each routine (F, dF, and ddF) which may contain any necessary user data.
 *  \param e An error handler.
 *  \returns An implicitly defined manifold.
 */
MFImplicitMF MFIMFCreateAlgebraicSubroutine(int,int,
                              void (*)(int*,double*,int*,double*,void*,MFErrorHandler),
                              void (*)(int*,double*,int*,double*,void*,MFErrorHandler),
                              void (*)(int*,double*,int*,double*,void*,MFErrorHandler),
                              void *,
                              MFErrorHandler);

/*! \fn MFImplicitMF MFIMFCreateAlgebraicSubroutineWithRadius(int n,int k,void (*F)(int*,double*,int*,void*,double*,void*,MFErrorHandler*), void (*dF)(int*,double*,int*,void*,double*,void*,MFErrorHandler*), void (*ddF)(int*,double*,int*,void*,double*,void*,MFErrorHandler*), void *data, double R, MFErrorHandler e);
 *  \brief This constructor is exactly the same as calling MFIMFCreateAlgebraicSiubroutine and then MFIMFSetR.
 *
 *  \param n The embedding dimension (the dimension of the range of F).
 *  \param k The dimension of the solution manifold (the dimension of the range of F is n+k).
 *  \param F A subroutine that evaluates F.
 *  \param dF A subroutine that evaluates the derivative of F (if NULL differenceing will be used).
 *  \param ddF A subroutine that evaluates the second derivative of F (if NULL differenceing will be used).
 *  \param data A pointer which is passed to each routine (F, dF, and ddF) which may contain any necessary user data.
 *  \param R The radius.
 *  \param e An error handler.
 *  \returns An implicitly defined manifold.
 */
MFImplicitMF MFIMFCreateAlgebraicSubroutineWithRadius(int,int,void (*)(int*,double*,int*,double*,void*,MFErrorHandler),
                                                              void (*)(int*,double*,int*,double*,void*,MFErrorHandler),
                                                              void (*)(int*,double*,int*,double*,void*,MFErrorHandler),
                                                              void *,
                                                              double,
                                                              MFErrorHandler);

/*! @} */

/*! \addtogroup MFTorusMF
 *  @{
 */

/*! \fn MFImplicitMF MFIMFCreateTorus(double x,double y,double z,double RI,double RO,MFErrorHandler e);
 *  \brief Creates an implicitly defined manifold for a 2-torus embedded in 3-space, centered at (x,y,z) with
 *         inner radius RI and outer radius RO. If RO were 0. (don't! the eqs are singular then), the torus
 *         would be a circle of radius RI in the xy-plane. You are free to set RI<RO, but it won't look like a torus.
 *
 *  \param x The x coordinate of the center of the torus
 *  \param y The y coordinate of the center of the torus
 *  \param z The z coordinate of the center of the torus
 *  \param RI The inner radius.
 *  \param RO The outer radius.
 *  \param e A place to return errors.
 *  \returns An implicitly defined manifold.
 */
MFImplicitMF MFIMFCreateTorus(double,double,double,double,double,MFErrorHandler);

/*! @} */

/*! \addtogroup MFSphereMF
 *  @{
 */

/*! \fn MFImplicitMF MFIMFCreateSphere(double x,double y,double z,double R,MFErrorHandler e);
 *  \brief Creates an implicitly defined manifold for a 2-sphere embedded in 3-space, centered at (x,y,z) with
 *         radius R.
 *
 *  \param x The x coordinate of the center of the sphere
 *  \param y The y coordinate of the center of the sphere
 *  \param z The z coordinate of the center of the sphere
 *  \param R The radius.
 *  \param e A place to return errors.
 *  \returns An implicitly defined manifold.
 */
MFImplicitMF MFIMFCreateSphere(double,double,double,double,MFErrorHandler);

/*! @} */

/*! \addtogroup MFCircleMF
 *  @{
 */

/*! \fn MFImplicitMF MFIMFCreateCircle(double x,double y,double R,MFErrorHandler e);
 *  \brief Creates an implicitly defined manifold for a circle embedded in the plane, centered at (x,y) with
 *         radius R.
 *
 *  \param x The x coordinate of the center of the circle
 *  \param y The y coordinate of the center of the circle
 *  \param R The radius.
 *  \param e A place to return errors.
 *  \returns An implicitly defined manifold.
 */
MFImplicitMF MFIMFCreateCircle(double,double,double,MFErrorHandler);

/*! @} */

/*! \addtogroup MFPlaneMF
 *  @{
 */

/*! \fn MFImplicitMF MFIMFCreatePlane(MFErrorHandler e);
 *  \brief Creates the plane. (not a very exciting manifold I'll grant you).
 *
 *  \param e A place to return errors.
 *  \returns An implicitly defined manifold.
 */
MFImplicitMF MFIMFCreatePlane(MFErrorHandler);

/*! @} */

/*! \addtogroup MFNSpaceMF
 *  @{
 */

/*! \fn MFImplicitMF MFIMFCreateNSpaceWithRadius(int k,double R,MFErrorHandler e);
 *  \brief Creates Euclidean k-space embedded in k-space. This is the same as MFIMFCreateNSpace(k,e) and MFSetR(.,R,e)
 *
 *  \param   k The dimension of the manifold and the embedding space.
 *  \param   R The radius.
 *  \param e A place to return errors.
 *  \returns An implicitly defined manifold.
 */
MFImplicitMF MFIMFCreateNSpaceWithRadius(int,double,MFErrorHandler);

/*! \fn MFImplicitMF MFIMFCreateNSpace(int k,MFErrorHandler e);

 *  \brief Creates Euclidean k-space embedded in k-space.
 *
 *  \param k The dimension of the space.
 *  \param e A place to return errors.
 *  \returns An implicitly defined manifold.
 */
MFImplicitMF MFIMFCreateNSpace(int,MFErrorHandler);

/*! @} */

MFImplicitMF MFIMFCreateCSTR(double,double,double,double,double,double,double,double,double,double,double,double,MFErrorHandler);

/*! \addtogroup MFFlatMF
 *  @{
 */

/*! \fn MFImplicitMF MFIMFCreateFlat(int n, int k, double *o,double *v,MFErrorHandler e);
 *  \brief Creates a flat k-space embedded in n-space. o is an array of length n, and v is an array of length n*k giving 
 *         the tangent plane. The jth basis vector is v[i+n*j]. The manifold is o[i]+sum_j v[i+n*j].s[j].
 *
 *  \param n The dimension of the embedding space.
 *  \param k The dimension of the manifold.
 *  \param o A point on the manifold.
 *  \param v The (constant) tangent space.  The jth basis vector is v[i+n*j].
 *  \param e A place to return errors.
 *  \returns An implicitly defined manifold.
 */
MFImplicitMF MFIMFCreateFlat(int,int,double*,double*,MFErrorHandler);

/*! @} */

/*! \addtogroup MFPolygonIn3SpaceMF
 *  @{
 */

/*! \fn MFImplicitMF MFIMFCreatePolygonIn3SpaceWithRadius(int nv,double *v,double R,MFErrorHandler e);
 *  \brief Creates a Euclidean plane which contains the given polygon (which is assumed to be flat).
 *
 *  \param nv The number of vertices.
 *  \param v  The vertices stored as (v[0+3*iv],v[1+3*iv],v[2+3*iv]).
 *  \param R  A radius to use for charts on the manifold.
 *  \param e A place to return errors.
 *  \returns An implicitly defined manifold.
 */
MFImplicitMF MFIMFCreatePolygonIn3SpaceWithRadius(int,double*,double,MFErrorHandler);

/*! @} */

/*! \addtogroup MFEdgeIn3SpaceMF
 *  @{
 *
/
/*! \fn MFImplicitMF MFIMFCreateEdgeIn3SpaceWithRadius(double *o,double *d,double R,MFErrorHandler e);
 *  \brief Creates a Euclidean line which contains the segment between the points l and r (left and right).
 *
 *  \param o An array of length at least 3, with the coordinates of the left end point of the interval.
 *  \param d An array of length at least 3, with the coordinates of the right end point of the interval.
 *  \param R  A radius to use for charts on the manifold.
 *  \param e A place to return errors.
 *  \returns An implicitly defined manifold.
 */
MFImplicitMF MFIMFCreateEdgeIn3SpaceWithRadius(double*,double*,double,MFErrorHandler);

/*! @} */

/*! \addtogroup MFTPBVPMF
 *  @{
 */

/*!
 * \brief A function defining the right hand side of a two point boundary value problem.
 */
typedef void (*MFTPBVPFFUNCTION)(double,int,double*,int,double*,double*,double*,double*,MFErrorHandler);

/*!
 * \brief A function defining the boundary value equations of a two point boundary value problem.
 */
typedef void (*MFTPBVPAFUNCTION)(int,int,double*,double*,int,double*,double*,double*,double*,double*,MFErrorHandler);

/*!
 * \brief A function defining the integral piece of an integral constraint of a two point boundary value problem.
 */
typedef void (*MFTPBVPLFUNCTION)(int,double,int,double*,int,double*,double*,double*,double*,MFErrorHandler);

/*!
 * \brief A function defining the non-vector piece of an integral constraint of a two point boundary value problem.
 */
typedef void (*MFTPBVPMFUNCTION)(int,int,double*,double*,double*,MFErrorHandler);

/*! \fn MFImplicitMF MFIMFCreateTPBVP(int k, int nx,int nu,int np, MFTPBVPFFUNCTION f, MFTPBVPFFUNCTION fu, MFTPBVPFFUNCTION fl, int nbc, MFTPBVPAFUNCTION a, MFTPBVPAFUNCTION au, MFTPBVPAFUNCTION al, int nic, MFTPBVPLFUNCTION l, MFTPBVPLFUNCTION lu, MFTPBVPLFUNCTION ll, MFTPBVPMFUNCTION m, MFTPBVPMFUNCTION ml,MFErrorHandler e);
 *  \brief Creates a manifold which is the solution manifold of a two point boundary value problem with integral constraints.
 *         Keller's second order box scheme is used.
 *
 *  \param k The number of degrees of freedom (the dimension of the solution manifold).
 *  \param nx The number of mesh intervals to use in the discretization.
 *  \param nu The number of functions defined on the mesh.
 *  \param np The number of scalar parameters.
 *  \param f  The right hand siade of the ODE's u'=f(u,p).
 *  \param fu The first derivatives of the right hand side with respect to u.
 *  \param fl The first derivatives of the right hand side with respect to the parameters l.
 *  \param nbc The number of boundary conditions.
 *  \param a  The function which defines the boundary conditions a(u(0),u(1),p)=0
 *  \param au The derivative of the boundary conditions with respect to u.
 *  \param al The derivative of the boundary conditions with respect to the parameters l.
 *  \param nic The number of integral conditions
 *  \param l The function which defines the integral part of the integral conditions int_0^1 l(u(t),p) dt + m(l)=0
 *  \param lu The derivative of the integral part of the integral conditions with respect to u.
 *  \param ll The derivative of the integral part of the integral conditions with respect to the parameters l.
 *  \param m The function which defines the non-integral part of the integral conditions int_0^1 l(u(t),p) dt + m(l)=0
 *  \param ml The derivative of the non-integral part of the integral conditions with respect to the parameters l.
 *  \param e A place to return errors.
 *  \returns An implicitly defined manifold.
 */
MFImplicitMF MFIMFCreateTPBVP(int k, int nx,int nu,int np, MFTPBVPFFUNCTION f, MFTPBVPFFUNCTION fu, MFTPBVPFFUNCTION fl, int nbc, MFTPBVPAFUNCTION a, MFTPBVPAFUNCTION au, MFTPBVPAFUNCTION al, int nic, MFTPBVPLFUNCTION l, MFTPBVPLFUNCTION lu, MFTPBVPLFUNCTION ll, MFTPBVPMFUNCTION m, MFTPBVPMFUNCTION ml,MFErrorHandler);

/*! \fn void MFTPBVPSetEpsilon(MFImplicitMF M,double epsilon,MFErrorHandler e);
 * \brief Sets the tolerance on the distance between a linear approximation at the center of a chart and the 
 *        manifold.
 *  \param M An MFTPBVPMF
 *  \param epsilon The value for the parameter epsilon.
 *  \param e A place to return errors.
 */
void MFTPBVPSetEpsilon(MFImplicitMF,double,MFErrorHandler);

/*! \fn MFNVector MFTPBVPIntegrateForInitialSolution(MFImplicitMF M,double *u0,double *p,double *x,MFErrorHandler e);
 *  \brief Solves an initial value problem in place of a MFTPBVPMF, and returns the solution. This may be useful in
 *         constructing initial guesses.
 *
 *  \param M An MFTPBVPMF
 *  \param u0 An array of length nu with the initial condition.
 *  \param p  An array of length nl with the parameter values.
 *  \param x  The nx+1 mesh points on [0,1].
 *  \param e A place to return errors.
 *  \returns A solution (u(t),p).
 */
MFNVector MFTPBVPIntegrateForInitialSolution(MFImplicitMF,double*,double*,double*,MFErrorHandler);

/*! \fn MFNVector MFTPBVPIntegrateForTangent(MFImplicitMF M,MFNVector u,double *du0,double *dp,MFErrorHandler e);
 *  \brief Solves an initial value problem of the linearization of a MFTPBVPMF, and returns the solution. This may be useful in
 *         constructing initial approximations of the columns of the basis for the tangent space.
 *
 *  \param M An MFTPBVPMF
 *  \param u A solution (u(t),p) which is the "point" at which the variational equations are written.
 *  \param du0 The initial perturbation (at x=0).
 *  \param dp The perturbation of the parameters.
 *  \param e A place to return errors.
 *  \returns A solution (du(t),dp) that might be used as a basis vector for the tangent space.
 */
MFNVector MFTPBVPIntegrateForTangent(MFImplicitMF,MFNVector,double*,double*,MFErrorHandler);

/*! @} */

/*! \addtogroup MFCoupledPendula
 *  @{
 */

/*! \fn MFImplicitMF MFIMFCreatePendula(int nt,double kappa,double gamma,int windingno, MFErrorHandler e)
 *  \brief Creates a manifold which is the motion of a pair of linearly coupled pendula.
 *
 *  \param nt The number of time steps to use.
 *  \param kappa The coupling constant.
 *  \param gamma The damping parameter.
 *  \param windingno The number of 2 pi intervals the average angle increases before the motion repeats.
 *  \param e A place to return errors.
 *  \returns The manifold.
 */
MFImplicitMF MFIMFCreatePendula(int,double,double,int,MFErrorHandler);

/*! @} */

/*! \addtogroup MFSwallowTail
 *  @{
 */

/*! \fn MFImplicitMF MFIMFCreateSwallow(MFErrorHandler e);
 *  \brief Creates a manifold which is the swallowtail catastrophe
 *
 *  \param e A place to return errors.
 *  \returns The swallowtail manifold.
 */
MFImplicitMF MFIMFCreateSwallow(MFErrorHandler);

/*! @} */

/*! \addtogroup MFForcedOscillator
 *  @{
 */

/*! \fn MFImplicitMF MFIMFCreateForcedOscillator(int nt,int windingno,MFErrorHandler e)
 *  \brief Creates a manifold which is the motion of a forced oscillator
 *
 *  \param nt The number of time steps to use.
 *  \param windingno The number of 2 pi intervals the angle increases before the motion repeats.
 *  \param e A place to return errors.
 *  \returns The manifold.
 */
MFImplicitMF MFIMFCreateForcedOscillator(int,int,MFErrorHandler);

/*! @} */

/*! \fn char *MFImplicitMFId(MFImplicitMF M,MFErrorHandler e);
 *  \brief Returns a character string which assigns a type to an ImplicitMF. Do not change or free the string!
 *
 *  \param M The implicitly defined manifold.
 *  \param e A place to return errors.
 *  \returns The type of the ImplicitMF.
 */
char *MFImplicitMFId(MFImplicitMF,MFErrorHandler);

/*! \fn int MFIMF_N(MFImplicitMF M,MFErrorHandler e);
 *  \brief Gets the embedding space dimension of an ImplicitMF.
 *
 *  \param M The implicitly defined manifold.
 *  \param e A place to return errors.
 *  \returns The dimension of the embedding space (n).
 */
int MFIMF_N(MFImplicitMF,MFErrorHandler);

/*! \fn int MFIMF_K(MFImplicitMF M,MFErrorHandler e);
 *  \brief Gets the dimension of an ImplicitMF.
 *
 *  \param M The implicitly defined manifold.
 *  \param e A place to return errors.
 *  \returns The dimension of the manifold (k).
 */
int MFIMF_K(MFImplicitMF,MFErrorHandler);

/*! \fn MFNSpace MFIMFNSpace(MFImplicitMF M,MFErrorHandler e);
 *  \brief Gets the space in which the manifold is embedded.
 *
 *  \param M The implicitly defined manifold.
 *  \param e A place to return errors.
 *  \returns The embedding space.
 */
MFNSpace MFIMFNSpace(MFImplicitMF,MFErrorHandler);

/*! \fn int MFIMFProjectFromCenter(MFImplicitMF M,MFNVector u0,MFNKMatrix Phi0,MFKVector s,MFNVector u,MFErrorHandler e);
 *  \brief This is the generalizatio of the prediction and correction operations in pseudo-arclength continuation.
 *         The prediction is u0+Phi0.s, and the correction is from this point to M orthogonal to Phi0.
 *
 *  \param M The implicitly defined manifold.
 *  \param u0 The point to be projected.
 *  \param Phi0 The tangent space.
 *  \param s The point in the tangent space that is to be projected.
 *  \param u The projected point, which lies on M and Phi0^T(u-u0)-0.
 *  \param e A place to return errors.
 *  \returns TRUE if the projection was sucessful
 */
int MFIMFProjectFromCenter(MFImplicitMF,MFNVector,MFNKMatrix,MFKVector,MFNVector,MFErrorHandler);

/*! \fn int MFIMFProject(MFImplicitMF M,MFNVector u0 ,MFNKMatrix Phi0,MFNVector u,MFErrorHandler e);
 *  \brief Projects a point u0 onto an ImplicitMF orthogonal to a tangent space Phi0
 *
 *  \param M The implicitly defined manifold.
 *  \param u0 The point to be projected.
 *  \param Phi0 The tangent space.
 *  \param u The projected point, which lies on M and Phi0^T(u-u0)-0.
 *  \param e A place to return errors.
 *  \returns TRUE if the projection was sucessful
 */
int MFIMFProject(MFImplicitMF,MFNVector,MFNKMatrix,MFNVector,MFErrorHandler);

/*! \fn MFNKMatrix MFIMFTangentSpace(MFImplicitMF M,MFNVector u,MFErrorHandler e);
 *  \brief Computes and returns an orthonormal basis for the tangent space of M at the point u on M.
 *
 *  \param M The implicitly defined manifold.
 *  \param u The point at which the tangent space is needed.
 *  \param e A place to return errors.
 *  \returns An orthonormal basis for the tangent space of M at u on M.
 */
MFNKMatrix MFIMFTangentSpace(MFImplicitMF,MFNVector,MFErrorHandler);

/*! \fn MFNKMatrix MFIMFTangentSpaceWithGuess(MFImplicitMF M,MFNVector u,MFNKMatrix guess,MFErrorHandler e);
 *  \brief Computes and returns an orthonormal basis for the tangent space of M at the point u on M.
 *
 *  \param M The implicitly defined manifold.
 *  \param u The point at which the tangent space is needed.
 *  \param guess An approximate tangent space. The new basis will align roughly with the basis in guess.
 *  \param e A place to return errors.
 *  \returns An orthonormal basis for the tangent space of M at u on M.
 */
MFNKMatrix MFIMFTangentSpaceWithGuess(MFImplicitMF,MFNVector,MFNKMatrix,MFErrorHandler);

/*! \fn double MFIMFScale(MFImplicitMF M,MFNVector u,MFNKMatrix Tan,MFErrorHandler e);
 *  \brief Estimates the radius of a spherical ball in the tangent space of M at u on M.
 *
 *  \param M The implicitly defined manifold.
 *  \param u The point at which the tangent space is needed.
 *  \param Tan The tangent space of M at u.
 *  \param e A place to return errors.
 *  \returns The radius.
 */
double MFIMFScale(MFImplicitMF,MFNVector,MFNKMatrix,MFErrorHandler);

/*! \fn int MFIMFProjectToSave(MFImplicitMF M,MFNVector u,double *x,MFErrorHandler e);
 *  \brief Projects a point so that it can be saved to file (used for the centerfile).
 *
 *  The convention is that if x is NULL, the routine returns the dimension of x, so that the user can allocate storage.
 *
 *  \param M The implicitly defined manifold.
 *  \param u The point that is to be projected, or (double*)NULL if the dimension is needed.
 *  \param x The implicitly defined manifold.
 *  \param e A place to return errors.
 *  \returns The dimension of the projected point if x is (double*)NULL.
 */
int MFIMFProjectToSave(MFImplicitMF,MFNVector,double*,MFErrorHandler);

/*! \fn int MFIMFProjectToDraw(MFImplicitMF M,MFNVector u,double *x,MFErrorHandler e);
 *  \brief Projects a point so that it can be drawn (used for the plotfile).
 *
 *  The convention is that if x is NULL, the routine returns the dimension of x, so that the user can allocate storage.
 *
 *  \param M The implicitly defined manifold.
 *  \param u The point that is to be projected, or (double*)NULL if the dimension is needed.
 *  \param x The implicitly defined manifold.
 *  \param e A place to return errors.
 *  \returns The dimension of the projected point if x is (double*)NULL.
 */
int MFIMFProjectToDraw(MFImplicitMF,MFNVector,double*,MFErrorHandler);

/*! \fn int MFIMFProjectToBB(MFImplicitMF M,MFNVector u,double *x,MFErrorHandler e);
 *  \brief Projects a point so that it can be placed in the hierarchical bounding boxes that are used to speed up
 *         the location of neighboring charts.
 *
 *  The convention is that if x is NULL, the routine returns the dimension of x, so that the user can allocate storage.
 *
 *  \param M The implicitly defined manifold.
 *  \param u The point that is to be projected, or (double*)NULL if the dimension is needed.
 *  \param x The implicitly defined manifold.
 *  \param e A place to return errors.
 *  \returns The dimension of the projected point if x is (double*)NULL.
 */
int MFIMFProjectToBB(MFImplicitMF,MFNVector,double*,MFErrorHandler);

/*! \fn void *MFIMFGetData(MFImplicitMF M,MFErrorHandler e);
 *  \brief Returns a pointer to the internal data used by the manifold. Be very careful using this. You need to know what you're
 *         doing.
 *
 *  \param M The implicitly defined manifold.
 *  \param e A place to return errors.
 *  \returns A pointer to the internal data.
 */
void *MFIMFGetData(MFImplicitMF,MFErrorHandler);

/*! \fn int (*MFIMFGetProjectForSave(MFImplicitMF M,MFErrorHandler e))(MFNVector,double*,void*,MFErrorHandler);
 *  \brief Returns a pointer to the routine used by the manifold to project for saving to file.
 *
 *  \param M The implicitly defined manifold.
 *  \param e A place to return errors.
 *  \returns A pointer to the routine.
 */
int (*MFIMFGetProjectForSave(MFImplicitMF M,MFErrorHandler e))(MFNVector,double*,void*,MFErrorHandler);

/*! \fn int (*MFIMFGetProjectForDraw(MFImplicitMF M,MFErrorHandler e))(MFNVector,double*,void*,MFErrorHandler);
 *  \brief Returns a pointer to the routine used by the manifold to project for saving to the plotfile.
 *
 *  \param M The implicitly defined manifold.
 *  \param e A place to return errors.
 *  \returns A pointer to the routine.
 */
int (*MFIMFGetProjectForDraw(MFImplicitMF M,MFErrorHandler e))(MFNVector,double*,void*,MFErrorHandler);

/*! \fn int (*MFIMFGetProjectForBB(MFImplicitMF M,MFErrorHandler e))(MFNVector,double*,void*,MFErrorHandler);
 *  \brief Returns a pointer to the routine used by the manifold to project for the boundaing boxes.
 *
 *  \param M The implicitly defined manifold.
 *  \param e A place to return errors.
 *  \returns A pointer to the routine.
 */
int (*MFIMFGetProjectForBB  (MFImplicitMF M,MFErrorHandler e))(MFNVector,double*,void*,MFErrorHandler);

/*! \fn double MFIMFGetR(MFImplicitMF M,MFErrorHandler e);
 *  \brief Returns the radius currently associated with a manifold. It provides a guess if no other information about the
 *          scale is available.
 *
 *  \param M The implicitly defined manifold.
 *  \param e A place to return errors.
 *  \returns The radius.
 */
double MFIMFGetR(MFImplicitMF,MFErrorHandler);

/*! \fn void MFIMFSetR(MFImplicitMF M,double R,MFErrorHandler e);
 *  \brief Associates a radius with a manifold. It provides a guess if no other information about the scale is available.
 *
 *  \param M The implicitly defined manifold.
 *  \param R The radius.
 *  \param e A place to return errors.
 */
void MFIMFSetR(MFImplicitMF,double,MFErrorHandler);

/*! \fn double MFIMFGetRMin(MFImplicitMF M,MFErrorHandler e);
 *  \brief Returns the current value of the radius that the user has provided as a smallest radius for the manifold.
 *
 *  \param M The implicitly defined manifold.
 *  \param e A place to return errors.
 *  \returns The radius.
 */
double MFIMFGetRMin(MFImplicitMF,MFErrorHandler);

/*! \fn void MFIMFSetRMin(MFImplicitMF M,double R,MFErrorHandler e);
 *  \brief Allows the user to set a smallest radius for a manifold. A continuation method may ignore it, but it is there.
 *
 *  \param M The implicitly defined manifold.
 *  \param R The radius.
 *  \param e A place to return errors.
 */
void MFIMFSetRMin(MFImplicitMF,double,MFErrorHandler);

/*! \fn MFNVector MFIMFVectorFactory(MFImplicitMF M,MFErrorHandler e);
 *  \brief This is a factory to create an NVector that is the approriate type for a point on this manifold. All allocation
 *         of vectors in multifario is by cloning, so the type of this vector is important. NVectors and ImplicitMF's are
 *         both base classes, so the user has no other way (besides the documentation) of knowing the type of NVector to
 *         use as an initial point.
 *
 *  \param M The implicitly defined manifold.
 *  \param e A place to return errors.
 *  \returns A clean and shiny new NVector of the right type for computations of this manifold.
 */
MFNVector MFIMFVectorFactory(MFImplicitMF,MFErrorHandler);

/*! \fn MFNKMatrix MFIMFMatrixFactory(MFImplicitMF M,MFErrorHandler e);
 *  \brief This is a factory to create an NKMatrix that is the approriate type for a point on this manifold. Unless
 *         the manifold uses a dense array of doubles the columns of the matrix will be NVectors of the appropriate type.
 *  \param M The implicitly defined manifold.
 *  \param e A place to return errors.
 *  \returns A clean and shiny new NKMatrix of the right type for the tangent space of this manifold.
 */
MFNKMatrix MFIMFMatrixFactory(MFImplicitMF,MFErrorHandler);

/*! \fn void MFRefImplicitMF(MFImplicitMF M,MFErrorHandler e);
 *  \brief Adds a reference to an MFImplicitMF.
 *
 *  \param M The an ImplicitMF being referenced.
 *  \param e A place to return errors.
 *  \sa ReferenceCounting MFFreeImplicitMF
 */
void MFRefImplicitMF(MFImplicitMF,MFErrorHandler);

/*! \fn void MFFreeImplicitMF(MFImplicitMF M,MFErrorHandler e);
 *  \brief Frees a reference to the ImplicitMF, and deletes the ImplicitMF if there are no references left.
 *
 *  \param M The ImplicitMF being unreferenced.
 *  \param e A place to return errors.
 *  \sa ReferenceCounting MFRefImplicitMF
 */
void MFFreeImplicitMF(MFImplicitMF,MFErrorHandler);

/*! \fn void MFWriteImplicitMF(FILE* fid,MFImplicitMF M,MFErrorHandler e);
 *  \brief Writes a ImplicitMF to a file.
 *
 *  \param fid The file to write to.
 *  \param M The ImplicitMF being written.
 *  \param e A place to return errors.
 */
void MFWriteImplicitMF(FILE*,MFImplicitMF,MFErrorHandler);

/*! \fn MFImplicitMF MFReadImplicitMF(FILE* fid,MFErrorHandler e);
 *  \brief Reads a ImplicitMF from a file.
 *
 *  \param fid The file to read from.
 *  \param e A place to return errors.
 */
MFImplicitMF MFReadImplicitMF(FILE*,MFErrorHandler);

/* These are new */

void MFIMFEvaluate(MFImplicitMF,MFNVector,MFNVector,MFErrorHandler);
void MFIMFApplyJacobian(MFImplicitMF,MFNVector,MFNKMatrix,MFNKMatrix,MFErrorHandler);
void MFIMFApplySecDer(MFImplicitMF,MFNVector,MFNVector,MFNVector,MFNVector,MFErrorHandler);
void MFIMFSetSpace(MFImplicitMF,MFNSpace,MFErrorHandler);

int MFIMFStop(MFImplicitMF,MFNVector,MFNKMatrix,MFNVector,MFNKMatrix,MFErrorHandler);

MFImplicitMF MFIMFCreateBaseClass(int,int,char*,MFErrorHandler);
void MFIMFSetData(MFImplicitMF,void*,MFErrorHandler);

void MFIMFSetWriteData(MFImplicitMF,void (*)(FILE*,void*,MFErrorHandler),MFErrorHandler);
void MFIMFSetFreeData(MFImplicitMF,void (*)(void*,MFErrorHandler),MFErrorHandler);
void MFIMFSetProjectFromCenter(MFImplicitMF,int (*)(int,int,MFNVector,MFNKMatrix, MFKVector,MFNVector,void*,int*,MFErrorHandler),MFErrorHandler);
void MFIMFSetProject(MFImplicitMF,int (*)(int,int,MFNVector,MFNKMatrix,MFNVector,void*,int*,MFErrorHandler),MFErrorHandler);
void MFIMFSetTangent(MFImplicitMF,int (*)(int,int,MFNVector,MFNKMatrix,void*,MFErrorHandler),MFErrorHandler);
void MFIMFSetTangentWithGuess(MFImplicitMF,int (*)(int,int,MFNVector,MFNKMatrix,MFNKMatrix,void*,MFErrorHandler),MFErrorHandler);
void MFIMFSetEvaluate(MFImplicitMF,void (*)(int,MFNVector,MFNVector,void*,MFErrorHandler),MFErrorHandler);
void MFIMFSetApplyJacobian(MFImplicitMF,void (*)(int,int,MFNVector,MFNKMatrix,MFNKMatrix,void*,MFErrorHandler),MFErrorHandler);
void MFIMFSetApplySecDer(MFImplicitMF,void (*)(int,int,MFNVector,MFNVector,MFNVector,MFNVector,void*,MFErrorHandler),MFErrorHandler);
void MFIMFSetScale(MFImplicitMF,double (*)(int,int,MFNVector,MFNKMatrix,void*,MFErrorHandler),MFErrorHandler);
void MFIMFSetStop(MFImplicitMF,int (*)(MFImplicitMF,MFNVector,MFNKMatrix,MFNVector,MFNKMatrix,void*,MFErrorHandler),MFErrorHandler);
void MFIMFSetProjectForSave(MFImplicitMF,int (*)(MFNVector,double*,void*,MFErrorHandler),MFErrorHandler);
void MFIMFSetProjectForDraw(MFImplicitMF,int (*)(MFNVector,double*,void*,MFErrorHandler),MFErrorHandler);
void MFIMFSetProjectForBB(MFImplicitMF,int (*)(MFNVector,double*,void*,MFErrorHandler),MFErrorHandler);
void MFIMFSetVectorFactory(MFImplicitMF,MFNVector (*)(MFImplicitMF,MFErrorHandler),MFErrorHandler);
void MFIMFSetMatrixFactory(MFImplicitMF,MFNKMatrix (*)(MFImplicitMF,MFErrorHandler),MFErrorHandler);

void MFIMFSetSingular(MFImplicitMF,int (*)(int,int,MFNVector,MFNKMatrix,MFNVector,void*,MFErrorHandler),MFErrorHandler);
int MFIMFSingular(MFImplicitMF,MFNVector,MFNKMatrix,MFNVector,MFErrorHandler);

void MFIMFSetSetStability(MFImplicitMF,void (*)(MFImplicitMF,MFNVector,MFNKMatrix,void*,MFErrorHandler),MFErrorHandler);
void MFIMFSetStability(MFImplicitMF,MFNVector,MFNKMatrix,MFErrorHandler);

MFImplicitMF MFIMFCreateOrthogonalFoliation(IMFFlow,MFErrorHandler);

#ifdef __cplusplus
}
#endif

/*! @} */

#endif
