/*! \defgroup MFNKMatrix */

/*! \addtogroup MFAtlas
 *  @{
 */

/*! \addtogroup MFNKMatrix
 *  @{
 */

/*! \class MFNKMatrix  MFNKMatrix.h MFNKMatrix.h
 *  \brief An n by k matrix, mainly used to store a basis for the tangent space of a k-dimensional manifold
 *  embedded in an n- dimensional space.
 *
 *  There are several types of MFNKMatrices. The user should use the Matrix factory routine in an implicitly defined
 *  manifold so that the correct type for that manifold is used.
 */

/*! \fn MFNKMatrix MFCreateNKMatrix(int k,MFNVector *columns,MFErrorHandler e);
 *  \brief Creates an MFNKMatrix which is stored as k MFNVectors.
 *
 *  \param k The number of columns.
 *  \param columns A list of the columns.
 *  \param e A place to return errors.
 *  \returns A new MFNKMatrix.
 */
MFNKMatrix MFCreateNKMatrix(int k,MFNVector *columns,MFErrorHandler e);

/*! \fn MFNKMatrix MFCreateNKMatrixWithData(int n,int k,double *entries,MFErrorHandler e);
 *  \brief Creates an MFNKMatrix which is stored as an array of doubles.
 *
 *  \param n The number of rows.
 *  \param k The number of columns.
 *  \param entries A list of n*k doubles which are copied into the entries of the new matrix. Entry (i,j) is entries[i+n*k].
 *  \param e A place to return errors.
 *  \returns A new MFNKMatrix.
 */
MFNKMatrix MFCreateNKMatrixWithData(int n,int k,double *entries,MFErrorHandler e);

/*! \fn int MFNKMatrixK(MFNKMatrix A,MFErrorHandler e);
 *  \brief Get the number of columns in the basis.
 *
 *  \param A The matrix containing the basis.
 *  \param e A place to return errors.
 *  \returns The number of columns, k.
 */
int MFNKMatrixK(MFNKMatrix A,MFErrorHandler e);

/*! \fn int MFNKMatrixN(MFNKMatrix A,MFErrorHandler e);
 *  \brief Get the number of coordinate values in the basis vectors.
 *
 *  \param A The matrix containing the basis.
 *  \param e A place to return errors.
 *  \returns The number of coordinate values in a column, n.
 */
int MFNKMatrixN(MFNKMatrix A,MFErrorHandler e);

/*! \fn MFNVector MFMColumn(MFNKMatrix A,int col,MFErrorHandler e);
 *  \brief Returns a column of the matrix A. A reference has been added to the n-vector, so the user must "free" it.
 *
 *  \param A The matrix containing the basis.
 *  \param col  Which column
 *  \param e A place to return errors.
 *  \returns The column.
 */
MFNVector MFMColumn(MFNKMatrix A,int col,MFErrorHandler e);

/*! \fn void MFMRow(MFNKMatrix A,int row,MFKVector s,MFErrorHandler e);
 *  \brief Copies a row of the matrix A into a k-vector.
 *
 *  \param A The matrix containing the basis.
 *  \param row  Which row
 *  \param s The row.
 *  \param e A place to return errors.
 */
void MFMRow(MFNKMatrix A,int row,MFKVector s,MFErrorHandler e);

/*! \fn void MFNKMSetC(MFNKMatrix A,int row,int col,double value,MFErrorHandler e);
 *  \brief Changes a single element of a basis. This is slow if many elements are to be changed. Get the column and downcast.
 *
 *  \param A The matrix containing the basis.
 *  \param row  Which row
 *  \param col  Which column
 *  \param value The new element value.
 *  \param e A place to return errors.
 */
void MFNKMSetC(MFNKMatrix A,int row,int col,double value,MFErrorHandler e);

/*! \fn void MFMVMul(MFNSpace space,MFNKMatrix A,MFKVector x,MFNVector b,MFErrorHandler e);
 *  \brief Matrix vector multiply.
 *
 *  \param space The space to use to add the columns.
 *  \param A The matrix containing the basis.
 *  \param x The vector.
 *  \param b The product Ax=b.
 *  \param e A place to return errors.
 */
void MFMVMul(MFNSpace space,MFNKMatrix A,MFKVector x,MFNVector b,MFErrorHandler e);

/*! \fn void MFMVMulT(MFNSpace space,MFNKMatrix A,MFNVector x,MFKVector b,MFErrorHandler e);
 *  \brief Matrix transpose vector multiply.
 *
 *  \param space The space to use to add the columns.
 *  \param A The matrix containing the basis.
 *  \param x The vector.
 *  \param b The product A^Tx=b.
 *  \param e A place to return errors.
 */
void MFMVMulT(MFNSpace space,MFNKMatrix A,MFNVector x,MFKVector b,MFErrorHandler e);

/*! \fn void MFGramSchmidt(MFNSpace space,MFNKMatrix A,MFErrorHandler e);
 *  \brief Orthonormalize the basis vectors.
 *
 *  \param space The space to use for norms and inner products.
 *  \param A The matrix containing the basis.
 *  \param e A place to return errors.
 */
void MFGramSchmidt(MFNSpace space,MFNKMatrix A,MFErrorHandler e);

/*! \fn void MFGramSchmidtNoMat(int n,int k,double *A,MFErrorHandler e);
 *  \brief A "raw" orthonormalization of the basis vectors stored as the columns of a dense array A[i+n*j]. A Euclidean space
 *         is assumed.
 *
 *  \param n The number of rows in A.
 *  \param k The number of columns in A.
 *  \param A The matrix containing the basis. A is length n*k, and stored by column A[i+n*j].
 *  \param e A place to return errors.
 */
void MFGramSchmidtNoMat(int n,int k,double *A,MFErrorHandler e);

/*! \fn void MFNKMProjectTangentForBranchSwitch(MFNSpace space,MFNKMatrix A,MFNVector x, MFNKMatrix B,MFErrorHandler e);
 *  \brief Projects the basis onto the vector x. A basis for the orthogonal complement of x in the span of the columns, plus 
 *         the unit vector in the x direction. This is used in the parallel search branch switching.
 *
 *  \param space The space to use for inner products.
 *  \param A The matrix containing the basis.
 *  \param x The vector defining the projection.
 *  \param B The projected basis.
 *  \param e A place to return errors.
 */
void MFNKMProjectTangentForBranchSwitch(MFNSpace space,MFNKMatrix A,MFNVector x, MFNKMatrix B,MFErrorHandler e);

/*! \fn MFNKMatrix MFCloneNKMatrix(MFNKMatrix A,MFErrorHandler e);
 *  \brief Clones (deep copy) a matrix.
 *
 *  \param A The matrix.
 *  \param e A place to return errors.
 *  \returns The clone.
 */
MFNKMatrix MFCloneNKMatrix(MFNKMatrix A,MFErrorHandler e);

/*! \fn void MFMSetColumn(MFNKMatrix A,int col,MFNVector v,MFErrorHandler e);
 *  \brief Replaces a column of a matrix.
 *
 *  \param A The matrix.
 *  \param col Which column.
 *  \param v The new column to use.
 *  \param e A place to return errors.
 */
void MFMSetColumn(MFNKMatrix A,int col,MFNVector v,MFErrorHandler e);

/*! \fn void MFRefNKMatrix(MFNKMatrix A, MFErrorHandler e);
 *  \brief Adds a reference to the matrix.
 *
 *  \param A The matrix being referenced.
 *  \param e A place to return errors.
 *  \sa ReferenceCounting MFFreeNKMatrix
 */
void MFRefNKMatrix(MFNKMatrix A,MFErrorHandler e);

/*! \fn void MFFreeNKMatrix(MFNKMatrix A, MFErrorHandler e);
 *  \brief Frees a reference to the A, and deletes the A if there are no references left.
 *
 *  \param A The matrix being unreferenced.
 *  \param e A place to return errors.
 *  \sa ReferenceCounting MFRefNKMatrix
 */
void MFFreeNKMatrix(MFNKMatrix A,MFErrorHandler e);

/*! \fn void MFWriteNKMatrix(FILE* fid,MFNKMatrix A, MFErrorHandler e);
 *  \brief Writes a A to a file.
 *
 *  \param fid The file to write to.
 *  \param A The matrix being queried.
 *  \param e A place to return errors.
 */
void MFWriteNKMatrix(FILE* fid,MFNKMatrix A,MFErrorHandler e);

/*! \fn MFNKMatrix MFReadNKMatrix(FILE* fid,MFErrorHandler e);
 *  \brief Reads a A from a file.
 *
 *  \param fid The file to write to.
 *  \param e A place to return errors.
 *  \returns The matrix.
 */
MFNKMatrix MFReadNKMatrix(FILE* fid,MFErrorHandler e);

/*! @} */

/*! @} */
