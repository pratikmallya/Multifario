/* 
    %W%
    %D% %T%

    PROGRAM NAME:  Manifold
   
    (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
    CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
   
    Please refer to the LICENSE file in the top directory

*/
/*      author: Mike Henderson mhender@watson.ibm.com */
/*      date:   February 12, 2007                     */

#ifndef __MFMULTIFARIOMETHOD_H__
#define __MFMULTIFARIOMETHOD_H__

#include <MFContinuationMethod.h>
#include <MFErrorHandler.h>

#ifdef __cplusplus
 extern "C" {
#endif

/*! \defgroup MFContinuationMethod */
/*! \defgroup MFMultifariosMethod */

/*! \addtogroup MFContinuationMethod
 *  @{
 */

/*! \addtogroup MFMultfiariosMethod
 *  @{
 */

/*! \fn MFContinuationMethod MFCreateMultfiariosMethod(MFErrorHandler e);
 *  \brief Creates the continuation algorithm in Multfiario 2001.
 *
 *  \param e A place to return errors.
 *  \return A new algorithm.
 */
MFContinuationMethod MFCreateMultifariosMethod(MFErrorHandler);

/*! \fn void MFFreeMultfiariosMethod(MFContinuationMethod algorithm,MFErrorHandler e);
 *  \brief Removes a reference to an algorithm, and deletes it if the reference count goes to zero.
 *
 *  \param algorithm The Algorithm.
 *  \param e A place to return errors.
 */
void MFFreeMultfiariosMethod(MFContinuationMethod algorithm,MFErrorHandler e);

/*! \fn void MFMultfiarioSetFilename(MFContinuationMethod algorithm,char *name,MFErrorHandler e);
 *  \brief When the atlas computed with this algorithm is close the atlas is written to the file name.atlas.
 *
 *  \param algorithm The Algorithm.
 *  \param name The name for the atlas file saved when the atlas is closed.
 *  \param e A place to return errors.
 */
void MFMultifarioSetFilename(MFContinuationMethod,char*,MFErrorHandler);

/*! \fn char *MFMultfiarioGetFilename(MFContinuationMethod algorithm,MFErrorHandler e);
 *  \brief Queries the file name to which the atlas is saved when the atlas is closed (name.atlas).
 *
 *  \param algorithm The Algorithm.
 *  \param e A place to return errors.
 *  \returns The name for the atlas file saved when the atlas is closed.
 */
char *MFMultifarioGetFilename(MFContinuationMethod,MFErrorHandler);

/*! \fn void MFMultfiarioAddClipF(MFContinuationMethod algorithm,double (*side)(MFNVector, MFErrorHandler),MFErrorHandler e);
 *  \brief Add a clipping plane to the computed atlas. This serves as a way to clean up plots. The chart polyhedra are clipped
 *         against the function "side". Linear interpolation is used on the value of "side" at the vertices, and the positive part
 *         is kept.
 *
 *  \param algorithm The Algorithm.
 *  \param side A function giving the distance from the clipping surface.
 *  \param e A place to return errors.
 */
void MFMultifarioAddClipF(MFContinuationMethod,double (*)(MFNVector,MFErrorHandler),MFErrorHandler);

/*! \fn void MFMultfiarioClearClipF(MFContinuationMethod algorithm,MFErrorHandler e);
 *  \brief Removes all clipping planes.
 *
 *  \param algorithm The Algorithm.
 *  \param e A place to return errors.
 */
void MFMultifarioClearClipF(MFContinuationMethod,MFErrorHandler);

/*!    \fn int MFMultfiarioSetIntegerParameter(MFContinuationMethod algorithm, char *parameterName, int value, MFErrorHandler e);
 *     \brief Allows the user to set integer parameters.
 *
 * Legal integer parameter names and their default values.
 *  <ul>
 *    <li>                   verbose                      default:    0
 *    <li>                   maxCharts                    default:    1
 *    <li>                   page                         default:    1
 *    <li>                   pageEvery                    default: 1000
 *    <li>                   useBB                        default:    1
 *    <li>                   dumpToPlotFile               default:    1
 *    <li>                   dumpToCenterFile             default:    1
 *    <li>                   dumpToRestartFile            default:    1
 *    <li>                   dumpToRestartFileEvery       default:  100
 *    <li>                   checkPoint                   default:    0
 *    <li>                   checkPointEvery              default:  100
 *    <li>                   branchSwitch                 default:    1
 *  </ul>
 *
 *  \param algorithm The Algorithm.
 *  \param parameterName Which parameter to set.
 *  \param value The new value.
 *  \param e A place to return errors.
 *  \returns FALSE if the parameter name does not match a parameter.
 */
int MFMultifarioSetIntegerParameter(MFContinuationMethod,char*,int,MFErrorHandler);

/*!    \fn int MFMultfiarioSetRealParameter(MFContinuationMethod algorithm, char *parameterName, double value, MFErrorHandler e);
 *     \brief Allows the user to set real valued parameters.
 *
 * Legal real parameter names and their default values.
 *   <ul>
 *     <li>                   minR                        default:    0.01
 *     <li>                   maxR                        default:    1.
 *     <li>                   epsilon                     default:    1.e-7
 *     <li>                   dotmin                      default:    0.2
 *   </ul>
 *
 *  \param algorithm The Algorithm.
 *  \param parameterName Which parameter to set.
 *  \param value The new value.
 *  \param e A place to return errors.
 *  \returns FALSE if the parameter name does not match a parameter.
 */
int MFMultifarioSetIntegerParameter(MFContinuationMethod,char*,int,MFErrorHandler);

/*!    \fn int MFMultfiarioGetIntegerParameter(MFContinuationMethod algorithm, char *parameterName,MFErrorHandler e);
 *     \brief Allows the user to set integer parameters. These are usually read from a file. Instead, default values
 *            are used when the Multfiario is created, and the user may change them.
 *
 * Legal integer parameter names.
 *  <ul>
 *    <li>                   verbose                      1 if output wanted, 0 OW
 *    <li>                   maxCharts                    The maximum number of charts to add.
 *    <li>                   page                         Whether or not to page interior charts.
 *    <li>                   pageEvery                    The interval (number of charts) between checks for paging.
 *    <li>                   useBB                        1 if bounding boxes are to be used.
 *    <li>                   dumpToPlotFile               1 if a plotfile, containing the chart polygons is wanted.
 *    <li>                   dumpToCenterFile             1 is a centerfile, with the centers of the charts, is wanted
 *    <li>                   dumpToRestartFile            1 is a restart file is wanted.
 *    <li>                   dumpToRestartFileEvery       The interval (number of charts) at which a restart file is written.
 *    <li>                   checkPoint                   1 is a checkpoint file is wanted.
 *    <li>                   checkPointEvery              The interval (number of charts) at which a checkpoint file is written.
 *    <li>                   branchSwitch                 1 if the algorithm is to branch switch, 0 OW.
 *  </ul>
 *
 *  \param algorithm A Multfiario continuation method.
 *  \param parameterName Which parameter value to retreive. A warning is issued if the parameter name does not match a parameter.
 *  \param e A place to return errors.
 *  \returns The current value of the parameter.
 */
int MFMultifarioGetIntegerParameter(MFContinuationMethod,char*,MFErrorHandler);

/*!    \fn double MFMultfiarioGetRealParameter(MFContinuationMethod algorithm, char *parameterName,MFErrorHandler e);
 *     \brief Allows the user to set real valued parameters. These are usually read from a file. Instead, default values
 *            are used when the Multfiario is created, and the user may change them.
 *
 * Legal real parameter names.
 *    <ul>
 *     <li>                   minR      The smallest chart radius to allow before giving up.
 *     <li>                   maxR      The largest chart radius to allow.
 *     <li>                   epsilon   The maximum distance allowed between the tangent space and manifold. (controls stepsize, MFErrorHandler e);
 *     <li>                   dotmin    The smallest dot product allowed between adjacent tangent spaces (controls stepsize). A
 *                                              value of 1 means that the tangent spaces must be parallel.
 *    </ul>
 *
 *  \param algorithm The Algorithm.
 *  \param parameterName Which parameter value to retreive. A warning is issued if the parameter name does not match a parameter.
 *  \param e A place to return errors.
 *  \returns The current value of the parameter.
 */
double MFMultifarioGetRealParameter(MFContinuationMethod,char*,MFErrorHandler);

void MFMultifariosMethodClearClipF(MFContinuationMethod,MFErrorHandler);
int MFMultifarioSetRealParameter(MFContinuationMethod,char*,double,MFErrorHandler);
void MFMultifarioSetFilename(MFContinuationMethod,char*,MFErrorHandler);
char *MFMultifarioGetFilename(MFContinuationMethod,MFErrorHandler);
int MFMultifarioGetIntegerParameter(MFContinuationMethod,char*,MFErrorHandler);
int MFMultifarioSetRealParameter(MFContinuationMethod,char*,double,MFErrorHandler);
double MFMultifarioGetRealParameter(MFContinuationMethod,char*,MFErrorHandler);

/*! @} */

/*! @} */

#ifdef __cplusplus
}
#endif

#endif
