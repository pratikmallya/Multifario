/* 
    @(#)MFKVector.h	1.3
    02/04/19 14:41:38
   
    PROGRAM NAME:  Manifold
   
    (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
    CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
   
    Please refer to the LICENSE file in the top directory

*/
/*      author: Mike Henderson mhender@watson.ibm.com */
/*      date:   November 11, 1997                     */
/*              February 2, 1999   Ported to C        */

#ifndef __MFKVECTOR_H__
#define __MFKVECTOR_H__
#include <MFBase.h>

#include <stdio.h>
#include <MFErrorHandler.h>

/*! \defgroup MFKVector */

/*! \addtogroup MFAtlas
 *  @{
 */

/*! \addtogroup MFKVector
 *  @{
 */

/*! \class MFKVector MFKVector.h MFKVector.h
 *  \brief A euclidean vector with that same dimension as a manifold.
 *
 *  An MFKVector represents a point in the domain of a chart.
 */
struct MFKVectorSt;
typedef struct MFKVectorSt *MFKVector;

#ifdef __cplusplus
 extern "C" {
#endif

/*! \fn MFKVector MFCreateKVector(int k,MFErrorHandler e);
 *  \brief Creates a k-dimensional vector.
 *
 *  \param k The dimension.
 *  \param e A place to return errors.
 *  \return A new k-vector.
 */
MFKVector MFCreateKVector(int,MFErrorHandler);

/*! \fn MFKVector MFCreateKVectorWithData(int k,double *c,MFErrorHandler e);
 *  \brief Creates a k-dimensional vector.
 *
 *  \param k The dimension.
 *  \param c A list of k doubles to be the coordinates of the new vector. Values are copied, so the user should free the
 *         list of doubles when appropriate. The copy is deleted when the k-vector is deleted.
 *  \param e A place to return errors.
 *  \return A new k-vector.
 */
MFKVector MFCreateKVectorWithData(int,double*,MFErrorHandler);

/*! \fn void MFRefKVector(MFKVector s,MFErrorHandler e);
 *  \brief Adds a reference to a k-dimensional vector.
 *
 *  \param s The k-vector.
 *  \param e A place to return errors.
 */
void MFRefKVector(MFKVector,MFErrorHandler);

/*! \fn void MFFreeKVector(MFKVector s,MFErrorHandler e);
 *  \brief Removes a reference to a k-dimensional vector, and deletes it if the number of references becomes zero.
 *
 *  \param s The k-vector.
 *  \param e A place to return errors.
 */
void MFFreeKVector(MFKVector,MFErrorHandler);

/*! \fn int MFKV_NC(MFKVector s,MFErrorHandler e);
 *  \brief Returns k, the doimension of the k-vector.
 *
 *  \param s The k-vector.
 *  \param e A place to return errors.
 */
int MFKV_NC(MFKVector,MFErrorHandler);

/*! \fn double MFKV_C(MFKVector s,int i,MFErrorHandler e);
 *  \brief Returns a coordinate value of a k-vector.
 *
 *  \param s The k-vector.
 *  \param i Which coordinate (btween 0 and k-1).
 *  \param e A place to return errors.
 *  \returns The requested coordinate value.
 */
double MFKV_C(MFKVector,int,MFErrorHandler);

/*! \fn void MFKVSetC(MFKVector s,int i,double c,MFErrorHandler e);
 *  \brief Changes a coordinate value of a k-vector.
 *
 *  \param s The k-vector.
 *  \param i Which coordinate (btween 0 and k-1).
 *  \param c The new coordinate value.
 *  \param e A place to return errors.
 */
void MFKVSetC(MFKVector,int,double,MFErrorHandler);

/*! \fn double MFKVDot(MFKVector s0,MFKVector s1,MFErrorHandler e);
 *  \brief Calculates the Euclidean inner product of two k-vectors.
 *
 *  \param s0 The first k-vector.
 *  \param s1 The second k-vector.
 *  \param e A place to return errors.
 *  \returns The inner product <s0,s1>.
 */
double MFKVDot(MFKVector,MFKVector,MFErrorHandler);

/*! \fn double MFKVNorm(MFKVector s,MFErrorHandler e);
 *  \brief Calculates the Euclidean norm of a k-vector.
 *
 *  \param s The k-vector.
 *  \param e A place to return errors.
 *  \returns The norm: sqrt(<s,s>)
 */
double MFKVNorm(MFKVector,MFErrorHandler);

/*! \fn void MFKVScale(double a,MFKVector s,MFErrorHandler e);
 *  \brief Multiply a k-vector by a scalar. The result is placed back in s.
 *
 *  \param a The scalar.
 *  \param s The k-vector.
 *  \param e A place to return errors.
 */
void MFKVScale(double,MFKVector,MFErrorHandler);

/*! \fn void MFKVScaleMul(double a,MFKVector s,MFKVector t,MFErrorHandler e);
 *  \brief Multiply a k-vector by a scalar.
 *
 *  \param a The scalar.
 *  \param s The k-vector.
 *  \param t The product: t=a*s.
 *  \param e A place to return errors.
 */
void MFKVScaleMul(double,MFKVector,MFKVector,MFErrorHandler);

/*! \fn void MFKVDiff(MFKVector a,MFKVector b,MFKVector c,MFErrorHandler e);
 *  \brief Find the difference of two k-vectors. The result is put in c.
 *
 *  \param a The first k-vector.
 *  \param b The second k-vector.
 *  \param c The difference: c=a-b.
 *  \param e A place to return errors.
 */
void MFKVDiff(MFKVector,MFKVector,MFKVector,MFErrorHandler);

/*! \fn void MFKVAdd(MFKVector a,MFKVector b,MFKVector c,MFErrorHandler e);
 *  \brief Find the sum of two k-vectors. The result is put in c.
 *
 *  \param a The first k-vector.
 *  \param b The second k-vector.
 *  \param c The sum: c=a+b.
 *  \param e A place to return errors.
 */
void MFKVAdd(MFKVector,MFKVector,MFKVector,MFErrorHandler);

/*! \fn void MFWriteKVector(FILE *fid,MFKVector s,MFErrorHandler e);
 *  \brief Write a k-vector to a file.
 *
 *  \param fid The file descriptor.
 *  \param s The k-vector.
 *  \param e A place to return errors.
 */
void MFWriteKVector(FILE*,MFKVector,MFErrorHandler);

/*! \fn MFKVector MFReadKVector(FILE *fid,MFErrorHandler e);
 *  \brief Read a k-vector from a file.
 *
 *  \param fid The file descriptor.
 *  \param e A place to return errors.
 *  \returns The k-vector.
 */
MFKVector MFReadKVector(FILE*,MFErrorHandler);

double *MFKV_CStar(MFKVector,MFErrorHandler);

#ifdef __cplusplus
}
#endif

/*! @} */

/*! @} */

#endif
