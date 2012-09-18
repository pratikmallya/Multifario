#ifndef __INVMFNSPACE_H__
#define __INVMFNSPACE_H__

#include <MFNSpace.h>
#include <MFErrorHandler.h>

#ifdef __cplusplus
 extern "C" {
#endif

/*! \defgroup MFNSpace */
/*! \defgroup InvMFNSpace */

/*! \addtogroup MFNSpace
 *  @{
 */

/*! \addtogroup InvMFNSpace
 *  @{
 */

/*! \fn MFNSpace InvMFCreateNSpace(int n,MFErrorHandler e);
 *  \brief Creates an inner product space for the intersection of a sphere and expansion.
 *
 *  \param n The dimension of the space.
 *  \param e A place to return errors.
 *  \returns A new MFNSpace
 */
MFNSpace InvMFCreateNSpace(int, MFErrorHandler);

/*! @} */

/*! @} */

#ifdef __cplusplus
}
#endif

#endif
