#ifndef __MFERROR_H__
#define __MFERROR_H__
/*
    (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
    CORPORATION 06/21/2005.  ALL RIGHTS RESERVED.

    Please refer to the LICENSE file in the top directory
*/
/*      author: Mike Henderson mhender@watson.ibm.com */
/*      version: %W% %D% %T% */
/*      date:   June 21, 2005                     */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

/*! \defgroup MFAtlas */
/*! \defgroup MFErrorHandler */

/*! \addtogroup MFAtlas
 *  @{
 */

/*! \addtogroup MFErrorHandler
 *  @{
 */

/*! \class MFErrorHandler 
 *  \brief An error handler, which can be changed how errors are handled.
 *
 *  An error handler receives a call when an error occurs, and receives information about the severity
 *  of the error and where it occured.
 */
struct MFErrorHandlerSt;
typedef struct MFErrorHandlerSt *MFErrorHandler;

#ifdef __cplusplus
 extern "C" {
#endif

/*! \fn MFErrorHandler MFCreateErrorHandler();
 *  \brief Creates a new error handler.
 *
 *  \returns The new error handler.
 */
MFErrorHandler MFCreateErrorHandler();

/*! \fn void MFSetError(MFErrorHandler e,int severity,char *routine,char *message,int line,char *filename);
 *  \brief Reports an error to an error handler.
 *
 *  \param e The error handler.
 *  \param severity The severity of the error.
 *  \param routine The routine in which the error occured.
 *  \param message A message describing the error.
 *  \param line The line in the source code at which the error was issued (use the __LINE__ macro).
 *  \param filename The name of the source file in which the error was issued (use the __FILE__ macro).
 */
void MFSetError(MFErrorHandler,int,char*,char*,int,char*);

/*! \fn void MFRefErrorHandler(MFErrorHandler e);
 *  \brief Adds a reference to an Error Handler.
 *
 *  \param e The e being referenced.
 *  \sa ReferenceCounting MFFreeErrorHandler
 */
void MFRefErrorHandler(MFErrorHandler);

/*! \fn void MFFreeErrorHandler(MFErrorHandler e);
 *  \brief Frees a reference to a e, and deletes the e if there are no references left.
 *
 *  \param e The error handler being unreferenced.
 *  \sa ReferenceCounting MFRefErrorHandler
 */
void MFFreeErrorHandler(MFErrorHandler);

/*! \fn int MFGetNErrors(MFErrorHandler e);
 *  \brief Gets the number of errors that have occuered up to this point.
 *
 *  \param e The error handler.
 *  \returns The number of errors.
 */
int MFGetNErrors(MFErrorHandler);

/*! \fn int MFErrorHandlerGetSev(MFErrorHandler e,int errorNumber);
 *  \brief Gets the severity of a particular error. The errorNumber must be between 0 and MFGetNErrors(e)-1.
 *
 *  \param e The error handler.
 *  \param errorNumber Which error.
 *  \returns The severity of the error.
 */
int MFErrorHandlerGetSev(MFErrorHandler,int);

/*! \fn char *MFErrorHandlerGetRoutine(MFErrorHandler e,int errorNumber);
 *  \brief Gets the name of the routine where a particular error occured. The errorNumber must be between 0 and MFGetNErrors(e)-1.
 *
 *  \param e The error handler.
 *  \param errorNumber Which error.
 *  \returns The name of the routine.
 */
char *MFErrorHandlerGetRoutine(MFErrorHandler,int);

/*! \fn char *MFErrorHandlerGetMsg(MFErrorHandler e,int errorNumber);
 *  \brief Gets the message associated with a particular error. The errorNumber must be between 0 and MFGetNErrors(e)-1.
 *
 *  \param e The error handler.
 *  \param errorNumber Which error.
 *  \returns The message.
 */
char *MFErrorHandlerGetMsg(MFErrorHandler,int);

/*! \fn int MFErrorHandlerGetLine(MFErrorHandler e,int errorNumber);
 *  \brief Gets the line in the source file where a particular error was issued. The errorNumber must be between 0 and MFGetNErrors(e)-1.
 *
 *  \param e The error handler.
 *  \param errorNumber Which error.
 *  \returns The line number.
 */
int MFErrorHandlerGetLine(MFErrorHandler,int);

/*! \fn char *MFErrorHandlerGetFile(MFErrorHandler e,int errorNumber);
 *  \brief Gets the name of the source file where a particular error was issued. The errorNumber must be between 0 and MFGetNErrors(e)-1.
 *
 *  \param e The error handler.
 *  \param errorNumber Which error.
 *  \returns The name of the source file.
 */
char *MFErrorHandlerGetFile(MFErrorHandler,int);

/*! \fn int MFIsError(MFErrorHandler e);
 *  \brief Returns TRUE (1) if any error has been reported to the error handler.
 *
 *  \param e The error handler.
 *  \returns TRUE if an error has occured.
 */
int MFIsError(MFErrorHandler);

/*! \fn void MFClearErrors(MFErrorHandler e);
 *  \brief Clears the list of errors.
 *
 *  \param e The error handler.
 */
void MFClearErrors(MFErrorHandler);

/*! \fn void MFErrorHandlerOutOfMemory(MFErrorHandler e);
 *  \brief A special routine to report to the error handler that no space is available to malloc and it's ilk.
 *
 *  \param e The error handler.
 */
void MFErrorHandlerOutOfMemory(MFErrorHandler);

/*! @} */

/*! @} */

#ifdef __cplusplus
}
#endif

#endif
