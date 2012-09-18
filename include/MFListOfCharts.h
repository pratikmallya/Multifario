/* 
    @(#)MFListOfCharts.h	1.2
    02/04/19 14:41:46
   
    PROGRAM NAME:  Manifold
   
    (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
    CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
   
    Please refer to the LICENSE file in the top directory

*/
/*      author: Mike Henderson mhender@watson.ibm.com */

#ifndef __MFLISTOFCHARTS_H__
#define __MFLISTOFCHARTS_H__
#include <MFBase.h>
#include <MFErrorHandler.h>

#ifdef __cplusplus
 extern "C" {
#endif

/*! \defgroup MFListOfCharts */

/*! \addtogroup MFAtlas
 *  @{
 */

/*! \addtogroup MFListOfCharts
 *  @{
 */

/*! \class MFListOFCharts MFListOfCharts.h MFListOfCharts.h
 *  \brief Like it says, a list of charts. A return type for the binary tree search operation.
 */
struct MFListOfChartsSt;
typedef struct MFListOfChartsSt *MFListOfCharts;

/*! \fn MFListOfCharts MFCreateListOfCharts(int n,int *chartNumbers,MFErrorHandler e);
 *  \brief Creates a list of charts from a list of chart numbers.
 *
 *  \param n The number of charts in the list.
 *  \param chartNumbers An integer array of length at least n containing the identifying numbers, usually the
 *            position in the list of charts in an atlas.
 *  \param e A place to return errors.
 *  \returns The list.
 */
MFListOfCharts MFCreateListOfCharts(int,int*,MFErrorHandler);

/*! \fn int MFNumberOfIntersectingCharts(MFListOfCharts list,MFErrorHandler e);
 *  \brief Returns the number of charts in the list.
 *
 *  \param list The list.
 *  \param e A place to return errors.
 *  \returns The number of charts in the list.
 */
int MFNumberOfIntersectingCharts(MFListOfCharts,MFErrorHandler);

/*! \fn int MFIntersectingChart(MFListOfCharts list,int i,MFErrorHandler e);
 *  \brief Returns the identifying number of the indicted chart in the list.
 *
 *  \param list The list.
 *  \param i The element being queried.
 *  \param e A place to return errors.
 *  \returns The identifying number of the chart.
 */
int MFIntersectingChart(MFListOfCharts,int,MFErrorHandler);

/*! \fn void MFRefListOfCharts(MFListOfCharts list,MFErrorHandler e);
 *  \brief Adds a reference to the chart.
 *
 *  \param list The list being referenced.
 *  \param e A place to return errors.
 *  \sa ReferenceCounting MFFreeListOFCharts
 */
void MFRefListOfIntersectingCharts(MFListOfCharts,MFErrorHandler);

/*! \fn void MFFreeListOfCharts(MFListOfCharts list,MFErrorHandler e);
 *  \brief Frees a reference to the chart, and deletes the chart if there are no references left.
 *
 *  \param list The list being unreferenced.
 *  \param e A place to return errors.
 *  \sa ReferenceCounting MFRefListOFCharts
 */
void MFFreeListOfIntersectingCharts(MFListOfCharts,MFErrorHandler);

/*! @} */

/*! @} */

#ifdef __cplusplus
}
#endif

#endif
