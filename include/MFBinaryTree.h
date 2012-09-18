/* 
    @(#)MFBinaryTree.h	1.5
    02/07/26 09:10:35
   
    PROGRAM NAME:  Manifold
    
    (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
    CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
    
    Please refer to the LICENSE file in the top directory 

*/
/*      author: Mike Henderson mhender@watson.ibm.com */

#ifndef __MFBINARYTREE_H__
#define __MFBINARYTREE_H__
#include <MFBase.h>
#include <MFChart.h>
#include <MFListOfCharts.h>

#ifdef __cplusplus
 extern "C" {
#endif
#include <MFErrorHandler.h>
/*! \defgroup MFBinaryTree */

/*! \addtogroup MFAtlas
 *  @{
 */

/*! \addtogroup MFBinaryTree
 *  @{
 */

/*! \class MFBinaryTree 
 *  \brief A set of hierarchical bounding boxes for storing spherical balls.
 *
 *  An MFBinaryTree reduces the quadratic complexity of finding the balls which are within a distance R of 
 *  a point. At first balls are added to the root. When the number of balls exceeds a parameter value the
 *  root is split into two boxes, and balls are moved into the two subboxes according to which they lie in.
 *  Each node (or leaf, MFErrorHandler e); has a bounding box, and as balls are added to the box or it's children the bounding
 *  box is updated.
 */
struct MFBinaryTreeSt;
typedef struct MFBinaryTreeSt *MFBinaryTree;

/*! \fn MFBinaryTree MFCreateBinaryTree(int d,MFErrorHandler e);
 *  \brief Creates a binary tree in dimension d.
 *
 *  \param d The dimension
 *  \param e A place to return errors.
 *  \returns A new BinaryTree
 */
MFBinaryTree MFCreateBinaryTree(int,MFErrorHandler);

/*! \fn void MFFreeBinaryTree(MFBinaryTree tree, MFErrorHandler e);
 *  \brief Frees a reference to the tree, and deletes the tree if there are no references left.
 *
 *  \param tree The tree being unreferenced.
 *  \param e A place to return errors.
 *  \sa ReferenceCounting MFRefBinaryTree
 */
void MFFreeBinaryTree(MFBinaryTree tree,MFErrorHandler e);

/*! \fn void MFBinaryTreeAddChart(MFBinaryTree tree,int seq,double *x,double R,MFErrorHandler e);
 *  \brief adds a ball to a binary tree.
 *
 *  \param tree The BinaryTree
 *  \param seq A number that stays with the ball to identify it.
 *  \param x A double array of length at least d (the dimension of the binary tree, MFErrorHandler e); with the center of the ball.
 *  \param e A place to return errors.
 *  \param R The center of the ball.
 */
void MFBinaryTreeAddChart(MFBinaryTree,int,double*,double,MFErrorHandler);

/*! \fn MFListOfCharts MFCreateListOfIntersectingCharts(MFBinaryTree tree,int seq,double *x,double R,MFErrorHandler e);
 *  \brief Adds a ball to the tree and gets a list of charts (the sequence numbers, MFErrorHandler e); that overlap a ball.
 *
 *  \param tree The BinaryTree
 *  \param seq A number that identifies the ball.
 *  \param x A double array of length at least d (the dimension of the binary tree, MFErrorHandler e); with the center of the ball.
 *  \param R The center of the ball.
 *  \param e A place to return errors.
 *  \returns A list of balls.
 */
MFListOfCharts MFCreateListOfIntersectingCharts(MFBinaryTree,int,double*,double,MFErrorHandler);

/*! \fn MFListOfCharts MFCreateListOfNearbyCharts(MFBinaryTree tree,double *x,double R,MFErrorHandler e);
 *  \brief Gets a list of charts (the sequence numbers, MFErrorHandler e); that overlap a ball without adding the ball.
 *
 *  \param tree The BinaryTree
 *  \param x A double array of length at least d (the dimension of the binary tree, MFErrorHandler e); with the center of the ball.
 *  \param R The center of the ball.
 *  \param e A place to return errors.
 *  \returns A list of balls.
 */
MFListOfCharts MFCreateListOfNearbyCharts(MFBinaryTree,double*,double,MFErrorHandler);

/*! \fn void MFWriteBinaryTree(FILE* fid,MFBinaryTree tree, MFErrorHandler e);
 *  \brief Writes a tree to a file.
 *
 *  \param fid The file to write to.
 *  \param tree The tree being queried.
 *  \param e A place to return errors.
 */
void MFWriteBinaryTree(FILE*,MFBinaryTree,MFErrorHandler);

/*! \fn MFBinaryTree MFReadBinaryTree(FILE* fid,MFAtlas A, MFErrorHandler e);
 *  \brief Reads a tree from a file.
 *
 *  \param fid The file to write to.
 *  \param A The Atlas for the tree.
 *  \param e A place to return errors.
 *  \returns tree The tree.
 */
MFBinaryTree MFReadBinaryTree(FILE*,MFErrorHandler);

/*! \fn void MFRecomputeBoundingBoxes(MFBinaryTree tree,int seq,double *x,double R,MFErrorHandler e);
 *  \brief Updates the boundaing boxes as if the ball were being added, but does not add the ball.
 *
 *  \param tree The BinaryTree
 *  \param seq A number that identifies the ball.
 *  \param x A double array of length at least d (the dimension of the binary tree, MFErrorHandler e); with the center of the ball.
 *  \param R The center of the ball.
 *  \param e A place to return errors.
 */
void MFRecomputeBoundingBoxes(MFBinaryTree,int,double*,double,MFErrorHandler);

/*! \fn void MFRefBinaryTree(MFBinaryTree tree, MFErrorHandler e);
 *  \brief Adds a reference to the tree.
 *
 *  \param tree The tree being referenced.
 *  \param e A place to return errors.
 *  \sa ReferenceCounting MFFreeBinaryTree
 */
void MFRefBinaryTree(MFBinaryTree tree,MFErrorHandler e);


#ifdef __cplusplus
}
#endif

#endif

/*! @} */

/*! @} */
