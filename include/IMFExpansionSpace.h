#ifndef __IMFEXPANSIONSPACE_H__
#define __IMFEXPANSIONSPACE_H__
#include <IMFExpansion.h>
#include <MFNSpace.h>
#include <MFImplicitMF.h>

#ifdef __cplusplus
 extern "C" {
#endif
#include <MFErrorHandler.h>

/*! \defgroup MFNSpace */
/*! \defgroup IMFExpansionSpace */

/*! \addtogroup MFNSpace
 *  @{
 */

/*! \addtogroup MFExpansionSpace
 *  @{
 */

/*! \fn MFNSpace InvMFCreateNSpace(int n,MFErrorHandler e);
 *  \brief Creates an inner product space for the intersection of a sphere and expansion.
 *
 *  \param n The dimension of the space.
 *  \param e A place to return errors.
 *  \returns A new MFNSpace
 */
MFNSpace IMFCreateExpansionSpace(int,MFErrorHandler);

/*! @} */

/*! @} */

#ifdef __cplusplus
}
#endif

#endif
