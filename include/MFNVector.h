/* 
    @(#)MFNVector.h	1.14
    03/07/07 22:24:11
   
    PROGRAM NAME:  Manifold
   
    (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
    CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
   
    Please refer to the LICENSE file in the top directory

*/
/*      author: Mike Henderson mhender@watson.ibm.com */
/*      date:   November 11, 1997                     */
/*              February 2, 1999   Ported to C        */

#ifndef __MFNVECTOR_H__
#define __MFNVECTOR_H__
#include <MFBase.h>

#include <stdio.h>
#include <MFErrorHandler.h>

/*! \defgroup MFNVector */
/*! \defgroup MFWrappedNVector */
/*! \defgroup MFDenseNVector */

/*! \addtogroup MFNVector
 *  @{
 */

/*! \class MFNVector MFNVector.h MFNVector.h
 *  \brief A vector that lies in the embedding space of a matrix.
 *
 *  An MFNVector represents a vector in the embedding space of a manifold. It is a base class. That is,
 *  The domain of the chart is a spherical ball in a k-dimensional Euclidean space, where k is the dimension
 *  some Implicitly Defined Manifolds store information in the NVector and so the data structure for an
 *  NVector is not the same for all NVectors. For example, LOCA can distribute a solution across a number of 
 *  processors.
 *
 *  The user generally should not create a base class NVector directly. A constructor routine should exist if
 *  the type of NVector needed is known. Otherwise, use the Vector Factory in the implicitly defined manifold.
 *  This routine will return the type of NVector which is appropriate to that manifold.
 */
struct MFNVectorSt;
typedef struct MFNVectorSt *MFNVector;

#ifdef __cplusplus
 extern "C" {
#endif

/*! \fn MFNVector MFCloneNVector(MFNVector u,MFErrorHandler e);
 *  \brief Creates a copy of an MFNVector. The copy is the same type and has the same coordinates. Changing the coordinates
 *         of the copy cwwill not change the original's coordinates.
 *
 *  \param u An MFNVector.
 *  \param e A place to return errors.
 *  \return A clone of u.
 */
MFNVector MFCloneNVector(MFNVector,MFErrorHandler);


/*! \addtogroup MFDenseNVector
 *  @{
 */

/*! \fn MFNVector MFCreateNVector(int n,MFErrorHandler e);
 *  \brief Creates an n-vector with zero coordinates values, stored in an array.
 *
 *  \param n The number of coordinates of the MFNVector.
 *  \param e A place to return errors.
 *  \return A new MFNVector.
 */
MFNVector MFCreateNVector(int,MFErrorHandler);

/*! \fn MFNVector MFCreateNVectorWithData(int n,double *c,MFErrorHandler e);
 *  \brief Creates an n-vector with given coordinates values, stored in an array.
 *
 *  \param n The number of coordinates of the MFNVector.
 *  \param c An array of the coordinates values. The array is copied, so changing c will not change the coordinates values.
 *  \param e A place to return errors.
 *  \return A new MFNVector.
 */
MFNVector MFCreateNVectorWithData(int,double*,MFErrorHandler);

/*! @} */

/*! \addtogroup MFWrappedNVector
 *  @{
 */

/*! \fn MFNVector MFCreateWrappedNVector(int n,double *c,MFErrorHandler e);
 *  \brief Creates an n-vector with given coordinates values, stored in an array.
 *
 *  \param n The number of coordinates of the MFNVector.
 *  \param c An array of the coordinates values. The array is NOT copied, so changing c will change the coordinates values. The array should not be deleted by the user. 
 *  \param e A place to return errors.
 *  \return A new MFNVector.
 */
MFNVector MFCreateWrappedNVector(int,double*,MFErrorHandler);

/*! \fn void *MFNVectorGetData(MFNVector u,MFErrorHandler e);
 *  \brief Returns a pointer to the internal data of a vector. The user needs to know how it is stored.
 *
 *  \param u The MFNVector.
 *  \param e A place to return errors.
 *  \returns A pointer to the data.
 */
void *MFNVectorGetData(MFNVector,MFErrorHandler);

/*! @} */

/*! \fn void MFRefNVector(MFNVector u,MFErrorHandler e);
 *  \brief Adds a reference to the MFNVector.
 *
 *  \param u The MFNVector to reference.
 *  \param e A place to return errors.
 */
void MFRefNVector(MFNVector,MFErrorHandler);

/*! \fn void MFFreeNVector(MFNVector u,MFErrorHandler e);
 *  \brief Removes a reference to the MFNVector and if the number of references becomes zero deletes the MFNVector.
 *
 *  \param u The MFNVector to free.
 *  \param e A place to return errors.
 */
void MFFreeNVector(MFNVector,MFErrorHandler);

/*! \fn int MFNV_NC(MFNVector u,MFErrorHandler e);
 *  \brief Returns the dimension of a MFNVector.
 *
 *  \param u The MFNVector.
 *  \param e A place to return errors.
 *  \returns The dimension of u.
 */
int MFNV_NC(MFNVector,MFErrorHandler);

/*! \fn double MFNV_C(MFNVector u,int i,MFErrorHandler e);
 *  \brief Returns the value of the ith coordinate of a MFNVector. This is inefficient, check the type and get the data.
 *
 *  \param u The MFNVector.
 *  \param i The coordinate requested.
 *  \param e A place to return errors.
 *  \returns The coordinate value.
 */
double MFNV_C(MFNVector,int,MFErrorHandler);

/*! \fn void MFNVSetC(MFNVector u,int i ,double c,MFErrorHandler e);
 *  \brief Changes the value of the ith coordinate of a MFNVector. This is inefficient, check the type and get the data.
 *
 *  \param u The MFNVector.
 *  \param i The coordinate to be changed.
 *  \param c The new coordinate value.
 *  \param e A place to return errors.
 */
void MFNVSetC(MFNVector,int,double,MFErrorHandler);

/*! \fn void MFNVAdd(MFNVector a,MFNVector b,MFNVector c,MFErrorHandler e);
 *  \brief Adds a and b coordinatewise and puts the result in c. It is better to use the NSpace.
 *
 *  \param a The first vector.
 *  \param b The first vector.
 *  \param c a+b
 *  \param e A place to return errors.
 */
void MFNVAdd(MFNVector,MFNVector,MFNVector,MFErrorHandler);

/*! \fn void MFNVDiff(MFNVector a,MFNVector b,MFNVector c,MFErrorHandler e);
 *  \brief Takes the difference of a and b coordinatewise and puts the result in c. It is better to use the NSpace.
 *
 *  \param a The first vector.
 *  \param b The first vector.
 *  \param c a-b
 *  \param e A place to return errors.
 */
void MFNVDiff(MFNVector,MFNVector,MFNVector,MFErrorHandler);

char *MFNVGetType(MFNVector,MFErrorHandler);

int MFNVGetIndex(MFNVector,MFErrorHandler);
void MFNVSetIndex(MFNVector,int,MFErrorHandler);

int MFNVGetIndex2(MFNVector,MFErrorHandler);
void MFNVSetIndex2(MFNVector,int,MFErrorHandler);

/*! \fn void MFWriteNVector(FILE *fid,MFNVector u,MFErrorHandler e);
 *  \brief Writes a vector to a file.
 *
 *  \param fid The file.
 *  \param u The vector.
 *  \param e A place to return errors.
 */
void MFWriteNVector(FILE*,MFNVector,MFErrorHandler);

/*! \fn void MFReadNVector(FILE *fid,MFErrorHandler e);
 *  \brief Read a vector from a file.
 *
 *  \param fid The file.
 *  \param e A place to return errors.
 *  \returns The vector.
 */
MFNVector MFReadNVector(FILE*,MFErrorHandler);

/*! \fn char *MFNVGetId(MFNVector u,MFErrorHandler e);
 *  \brief Get the type of a vector. This is a string which is set when the vector was created. Don't delete it!
 *
 *  \param u The vector.
 *  \param e A place to return errors.
 *  \returns The type.
 */
char *MFNVGetId(MFNVector,MFErrorHandler);

double *MFNV_CStar(MFNVector,MFErrorHandler);

/* These are new or weren't exposed */

MFNVector MFCreateNVectorBaseClass(char*,MFErrorHandler);
void MFNVectorSetPrint(MFNVector,void (*)(FILE*,void*,MFErrorHandler),MFErrorHandler);
void MFNVectorSetData(MFNVector,void*,MFErrorHandler);
void MFNVectorSetFreeData(MFNVector,void (*)(void*,MFErrorHandler),MFErrorHandler);
void MFNVectorSetWriteData(MFNVector,void (*)(FILE*,void*,MFErrorHandler),MFErrorHandler);
void MFNVectorSetClone(MFNVector,MFNVector (*)(void*,MFErrorHandler),MFErrorHandler);
void MFNVectorSetGetNC(MFNVector,int (*)(void*,MFErrorHandler),MFErrorHandler);
void MFNVectorSetGetC(MFNVector,double (*)(int,void*,MFErrorHandler),MFErrorHandler);
void MFNVectorSetSetC(MFNVector,void (*)(int,double,void*,MFErrorHandler),MFErrorHandler);
void MFNVectorSetAdd(MFNVector,void (*)(void*,void*,void*,MFErrorHandler),MFErrorHandler);
void MFNVectorSetDiff(MFNVector,void (*)(void*,void*,void*,MFErrorHandler),MFErrorHandler);

#ifdef __cplusplus
}
#endif

/*! @} */

#endif
