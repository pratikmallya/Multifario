#ifndef __IMFFLOW_H__
#define __IMFFLOW_H__

#include <MFNVector.h>
#include <MFNKMatrix.h>
#include <MFErrorHandler.h>
#include <stdio.h>

#ifdef __cplusplus
 extern "C" {
#endif

/*! \defgroup Flows */
/*! \defgroup FatFlow */
/*! \defgroup BackwardsFatFlow */

/*! \addtogroup Flows
 *  @{
 */

/*! \class IMFFlow 
 *  \brief A flow is a way of defining a continuous dynamical system. It gives the time derivative of the state as a 
 *                 function of the state. If the flow depends on time it is call non-autonomous. An
 *                 IMFFlow represents an autonomous flow, which is independant of time.
 *
 *  A flow is a way of defining a continuous dynamical system. It gives the time derivative of the state as a 
 *  function of the state. If the flow depends on time it is call non-autonomous. An IMFFlow
 * represents an autonomous flow, which is independant of time. If the state variable is u, the flow is
 *
 *                        u'(t)=F(u(t))
 *
 *  where l is a vector of parameters.
 */
struct IMFFlowSt;
typedef struct IMFFlowSt *IMFFlow;

/*! typedef void (*MFFlowFunction)(double*,double*,double*,void*,MFErrorHandler e);
 *    /brief A convenience for the signature of the functions which evaluate the flow and its derivatives.
 */
typedef void (*MFFlowFunction)(double*,double*,double*,void*,MFErrorHandler);

/*! typedef void (*MFFlowFreeData)(void*,MFErrorHandler e);
 *    /brief A convenience for the signature of the "FreeData" routine of a flow.
 */
typedef void (*MFFlowFreeData)(void*,MFErrorHandler);

/*! \fn IMFFlow IMFCreateFlow(int nu, int np, MFFlowFunction F, MFFlowFunction dF, MFFlowFunction dFdp, MFFlowFunction ddF, MFFlowFunction dddF,void *data, MFFlowFreeData freeData,MFErrorHandler e);
 * \brief Creates a parameterized flow.
 *
 * \param nu The dimension of the phase space.
 * \param np The dimension of the parameter space.
 * \param F The flow. Assigns a vector to each point in phase space.
 * \param dF The derivative of the flow with respect to the phase space variables (i.e. the Jacobian).
 * \param dFdp The derivative of the flow with respect to the parameter space variables.
 * \param ddF The second derivative of the flow with respect to the phase space variables.
 * \param dddF The thired derivative of the flow with respect to the phase space variables.
 * \param data A parameter block that will be passed to F and it's derivatives.
 * \param freeData A routine that may be provided to release the parameter block when the flow is deleted.
 *  \param e A place to return errors.
 * \returns A flow.
 */
IMFFlow IMFCreateFlow(int,int,MFFlowFunction,MFFlowFunction,MFFlowFunction,MFFlowFunction,MFFlowFunction,void*,MFFlowFreeData,MFErrorHandler);

/*! \addtogroup FatFlow
 *  @{
 */

/*! \fn IMFFlow IMFCreateFatFlow(IMFFlow F,int k, MFErrorHandler e)
 * \brief Creates a fat parameterized flow. A fat flow is a flow for the point in phase space, an orthonormal basis for a
 *                           "tangent space" (it doesn't have to be the tangent space of anything, it's just a linear subspace),
 *                            and a quadratic in the normal space (orthogonal to the tangent space).
 *
 * \param F The flow that will be fattened.
 * \param k The dimension of the tangent space.
 * \returns The fattened flow.
 */
IMFFlow IMFCreateFatFlow(IMFFlow,int,MFErrorHandler);

/*! @} */

/*! \addtogroup FatFlow
 *  @{
 */

/*! \fn IMFFlow IMFCreateBackwardFatFlow(IMFFlow F,int k,MFErrorHandler e);
 * \brief Creates a fat parameterized flow, based on the reverse time flow of a given flow F. This is the same as
 *                                  IMFCreateFatFlow(IMFCreateBackwardFlow(F),k,MFErrorHandler e);
 *
 * \param F The flow that will be fattened.
 * \param k The dimension of the tangent space.
 *  \param e A place to return errors.
 * \returns The fattened flow.
 */
IMFFlow IMFCreateBackwardFatFlow(IMFFlow,int,MFErrorHandler);

/*! @} */

/*! \addtogroup FatFlow
 *  @{
 */

IMFFlow IMFCreateBackwardFlow(IMFFlow,MFErrorHandler);

/*! @} */

/*! \fn void IMFFreeFlow(IMFFlow flow,MFErrorHandler e);
 *  \brief Frees a reference to the flow, and deletes the flow if there are no references left.
 *
 *  \param flow The flow being unreferenced.
 *  \param e A place to return errors.
 *  \sa ReferenceCounting IMFRefFlow
 */
void IMFFreeFlow(IMFFlow,MFErrorHandler);

/*! \fn void IMFRefFlow(IMFFlow flow,MFErrorHandler e);
 *  \brief Adds a reference to the flow.
 *
 *  \param flow The flow being referenced.
 *  \param e A place to return errors.
 *  \sa ReferenceCounting IMFFreeFlow
 */
void IMFRefFlow(IMFFlow,MFErrorHandler);

/*! \fn void IMFEvaluateFlow(IMFFlow F, MFNVector vu, MFKVector vp, double *f,MFErrorHandler e);
 * \brief Evaluates the flow direction at a point in phase space.
 *
 * \param F The flow.
 * \param vu A point in phase space.
 * \param vp A point in parameter space.
 * \param f An array of length at least the dimension of the phase space to hold the flow vector.
 *  \param e A place to return errors.
 */
void IMFEvaluateFlow(IMFFlow,MFNVector,MFKVector,double*,MFErrorHandler);

/*! \fn void IMFEvaluateDerivativeOfFlow(IMFFlow F, MFNVector vu, MFKVector vp, double *df,MFErrorHandler e);
 * \brief Evaluates the derivative of the flow at a point in phase space (the Jacobian).
 *
 * \param F The flow.
 * \param vu A point in phase space.
 * \param vp A point in parameter space.
 * \param df An array of length at least the square of the dimension of the phase space to hold the Jacobian. It is
 *           filled with entries in column order. That is dF^i_j = df[i+nu*j], where nu is the dimension of the phase space.
 *  \param e A place to return errors.
 */
void IMFEvaluateDerivativeOfFlow(IMFFlow,MFNVector,MFKVector,double*,MFErrorHandler);

/*! \fn void IMFEvaluateParameterDerivativeOfFlow(IMFFlow F, MFNVector vu, MFKVector vp,double *dp,MFErrorHandler e);
 * \brief Evaluates the derivative of the flow at a point in phase space with respect to the parameter space variables (dFdp).
 *
 * \param F The flow.
 * \param vu A point in phase space.
 * \param vp A point in parameter space.
 * \param dp An array of length at least the dimension of the phase space times the dimension of the parameter space,
 *           to hold the derivatives. It is
 *           filled with entries in column order. That is dF^i_j = dp[i+nu*j], where nu is the dimension
 *           of the phase space.
 *  \param e A place to return errors.
 */
void IMFEvaluateParameterDerivativeOfFlow(IMFFlow,MFNVector,MFKVector,double*,MFErrorHandler);

/*! \fn void IMFEvaluateSecondDerivativeOfFlow(IMFFlow F, MFNVector vu, MFKVector vp,double *ddf,MFErrorHandler e);
 * \brief Evaluates the second derivative of the flow at a point in phase space (dFdudu).
 *
 * \param F The flow.
 * \param vu A point in phase space.
 * \param vp A point in parameter space.
 * \param ddf An array of length at least the cube of the dimension of the phase space to hold the derivatives. It is
 *           filled with entries in column order. That is dF^i_jk = df[i+nu*(j+nu*k)], where nu is the dimension of the phase space.
 *  \param e A place to return errors.
 */
void IMFEvaluateSecondDerivativeOfFlow(IMFFlow,MFNVector,MFKVector,double*,MFErrorHandler);

/*! \fn void IMFEvaluateThirdDerivativeOfFlow(IMFFlow F, MFNVector vu, MFKVector vp,double *ddf,MFErrorHandler e);
 * \brief Evaluates the third derivative of the flow at a point in phase space (dFdududu).
 *
 * \param F The flow.
 * \param vu A point in phase space.
 * \param vp A point in parameter space.
 * \param dddf An array of length at least the cube of the dimension of the phase space to hold the derivatives. It is
 *           filled with entries in column order. That is dF^i_jkl = df[i+nu*(j+nu*(k+nu*l))], where nu is the dimension
 *           of the phase space.
 *  \param e A place to return errors.
 */
void IMFEvaluateThirdDerivativeOfFlow(IMFFlow,MFNVector,MFKVector,double*,MFErrorHandler);

/*! \fn int IMFFlowNU(IMFFlow F,MFErrorHandler e);
 * \brief Returns the dimension of the phase space of the flow F.
 *
 * \param F The flow.
 *  \param e A place to return errors.
 * \returns The dimension of the phase space.
 */
int IMFFlowNU(IMFFlow,MFErrorHandler);

/*! \fn int IMFFlowNP(IMFFlow F,MFErrorHandler e);
 * \brief Returns the dimension of the parameter space of the flow F.
 *
 * \param F The flow.
 *  \param e A place to return errors.
 * \returns The dimension of the parameter space.
 */
int IMFFlowNP(IMFFlow,MFErrorHandler);

void *IMFFlowData(IMFFlow,MFErrorHandler);
MFFlowFreeData IMFFreeData(IMFFlow,MFErrorHandler);

/*! \fn double IMFFlowR(IMFFlow F, double eps, MFNVector u, MFKVector p, MFNKMatrix mPhi,MFErrorHandler e);
 * \brief Estimates the radius of a negihborhood of a point in phase space within which the flow F
 *        does not change more than eps.
 *
 * \param F A flow.
 * \param eps The bound on the change of F allowed within the neighborhood.
 * \param u A point in phase space.
 * \param p A point in parameter space.
 * \param mPhi An orthonormal basis for a subspace in phase space. The bound on the change in F need only hold in this
 *             subspace (a spherical ball in the subspace centered at the point in phase space).
 *  \param e A place to return errors.
 * \returns The size of a neighborhood in which the bound applies.
 */
double IMFFlowR(IMFFlow,double,MFNVector,MFKVector,MFNKMatrix,MFErrorHandler);

/*! \fn void IMFWriteFlow(FILE* fid,IMFFlow flow,MFErrorHandler e);
 *  \brief Writes a flow to a file.
 *
 *  \param fid The file to write to.
 *  \param flow The flow being queried.
 *  \param e A place to return errors.
 */
void IMFWriteFlow(FILE* fid,IMFFlow flow,MFErrorHandler e);

/*! \fn IMFFlow IMFReadFlow(FILE* fid,MFErrorHandler e);
 *  \brief Reads a flow from a file.
 *
 *  \param fid The file to write to.
 *  \param e A place to return errors.
 *  \returns flow The flow.
 */
IMFFlow IMFReadFlow(FILE*,MFErrorHandler);

/*! \fn void *IMFFlowData(IMFFlow F,MFErrorHandler e);
 * \brief Returns the data block that was provided when F was constructed. BE VERY CAREFUL WITH THIS.
 *
 * \param F The flow.
 *  \param e A place to return errors.
 * \returns The data block which belongs to the flow.
 */
void *IMFFlowData(IMFFlow F,MFErrorHandler e);

/*! \fn MFFlowFreeData IMFFlowFreeData(IMFFlow F,MFErrorHandler e);
 * \brief Returns the routine to free the data block that was provided when F was constructed.
 *
 * \param F The flow.
 *  \param e A place to return errors.
 * \returns The routine to free the data block.
 */
MFFlowFreeData IMFFlowFreeData(IMFFlow F,MFErrorHandler e);

/*! \fn void MFFlowSetNVectorFactory(IMFFlow this,MFNVector (*factory)(IMFFlow,MFErrorHandler), MFErrorHandler e);
 *  \brief This sets the factory to create an NVector that is the appropriate type for a point on the manifold.
 *
 *  \param this The flow.
 *  \param factory The routine to use as the factory.
 *  \param e A place to handle errors.
 */
void MFFlowSetNVectorFactory(IMFFlow this,MFNVector (*factory)(IMFFlow,MFErrorHandler), MFErrorHandler e);

/*! \fn void MFFlowSetKVectorFactory(IMFFlow this,MFKVector (*factory)(IMFFlow,MFErrorHandler), MFErrorHandler e);
 *  \brief This sets the factory to create an KVector that is the appropriate type for a point on the manifold.
 *
 *  \param this The flow.
 *  \param factory The routine to use as the factory.
 *  \param e A place to handle errors.
 */
void MFFlowSetKVectorFactory(IMFFlow this,MFKVector (*factory)(IMFFlow,MFErrorHandler), MFErrorHandler e);

/*! \fn MFNVector MFFlowNVectorFactory(IMFFlow F, MFError e);
 *  \brief This is a factory to create an NVector that is the approriate type for a point on thisIDM manifold. All allocation
 *         of vectors in multifario is by cloning, so the type of thisIDM vector is important. NVectors and ImplicitMF's are
 *         both base classes, so the user has no other way (besides the documentation) of knowing the type of NVector to
 *         use as an initial point.
 *
 *  \param F The flow.
 *  \param e A place to handle errors.
 *  \returns A clean and shiny new NVector of the right type for computations of thisIDM manifold.
 */
MFNVector MFFlowNVectorFactory(IMFFlow F, MFErrorHandler e);

/*! \fn MFKVector MFFlowKVectorFactory(IMFFlow F, MFError e);
 *  \brief This is a factory to create an NVector that is the approriate type for a point on thisIDM manifold. All allocation
 *         of vectors in multifario is by cloning, so the type of thisIDM vector is important. NVectors and ImplicitMF's are
 *         both base classes, so the user has no other way (besides the documentation) of knowing the type of NVector to
 *         use as an initial point.
 *
 *  \param F The flow.
 *  \param e A place to handle errors.
 *  \returns A clean and shiny new NVector of the right type for computations of thisIDM manifold.
 */
MFKVector MFFlowKVectorFactory(IMFFlow F, MFErrorHandler e);

/*! @} */

#ifdef __cplusplus
}
#endif

#endif
