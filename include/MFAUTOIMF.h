/*! \defgroup MFAUTOBV */
/*! \defgroup MFAUTOTPBVP */

/*! \addtogroup MFImplicitMF
 *  @{
 */

/*! \addtogroup MFAUTOBV
 *  @{
 */

/*! \class MFAUTOBV MFAUTO.h MFAUTO.h
 *  \brief A solution manifold of a boundary value problem posed in AUTO's representation.
 *
 */

/*! typedef int (*MFfunc_type)(integer,const doublereal*,const integer*,const doublereal*,integer,doublereal*,doublereal*,doublereal*,MFErrorHandler e);
 *
 *    The signature of the routine that AUTO requires to evaluate the flow.
 */
typedef int (*MFfunc_type)(integer,const doublereal*,const integer*,const doublereal*,integer,doublereal*,doublereal*,doublereal*,MFErrorHandler e);

/*! typedef int (*MFbcnd_type)(integer,const doublereal*,const integer*,integer,const doublereal*,const doublereal*,integer,doublereal*,doublereal*,MFErrorHandler e);
 *
 *    The signature of the routine that AUTO requires to evaluate the boundary conditions.
 */
typedef int (*MFbcnd_type)(integer,const doublereal*,const integer*,integer,const doublereal*,const doublereal*,integer,doublereal*,doublereal*,MFErrorHandler e);

/*! typedef int (*MFicnd_type)(integer,const doublereal*,const integer*,integer,const doublereal*,const doublereal*,const doublereal*,const doublereal*,integer,doublereal*,doublereal*,MFErrorHandler e);
 *
 *    The signature of the routine that AUTO requires to evaluate integral conditions.
 */
typedef int (*MFicnd_type)(integer,const doublereal*,const integer*,integer,const doublereal*,const doublereal*,const doublereal*,const doublereal*,integer,doublereal*,doublereal*,MFErrorHandler e);

/*! typedef int (*MFstpnt_type)(integer,doublereal,doublereal*,doublereal*,MFErrorHandler e);
 *
 *    The signature of the routine that AUTO requires to evaluate an initial solution over t in [0,1].
 */
typedef int (*MFstpnt_type)(integer,doublereal,doublereal*,doublereal*,MFErrorHandler e);

/*! typedef int (*MFpvls_type)(integer,const void*,doublereal*,MFErrorHandler e);
 *
 *    The signature of the routine that AUTO requires to project solutions of algebraic problems.
 */
typedef int (*MFpvls_type)(integer,const void*,doublereal*,MFErrorHandler e);

/*!    \fn MFImplicitMF MFCreateAUTOBV(MFAUTOTPBVP tpbvp, MFNSpace space,MFErrorHandler e);
 *     \brief Creates an implicit representation of the solution manifold of a two point boundary value problem.
 *
 *     \param tpbvp A two point boundary value problem of the type used by AUTO.
 *     \param space The space in which the solution manifold lives.
 *  \param e A place to return errors.
 *     \returns The solution manifold. 
 */
MFImplicitMF MFCreateAUTOBV(MFAUTOTPBVP tpbvp, MFNSpace space,MFErrorHandler e);

/*!    \fn int MFAUTOGetStartPoint(MFImplicitMF M,MFAUTOTPBVP tpbvp, MFstpnt_type stpnt,doublereal *p0,MFNVector *u0,MFNKMatrix *Phi0,MFErrorHandler e);
 *     \brief Creates a starting point for AUTO.
 *
 *     \param M The solution manifold.
 *     \param tpbvp The boundary value problem.
 *     \param stpnt The routine which evalutes the starting point at each point in time [0,1].
 *     \param p0 An array containing the parameter values of the starting point.
 *     \param u0 (Output) A place to put the start point.
 *     \param Phi0 (Output) A place to put the starting tangent, or NULL if the tangent is not wanted.
 *  \param e A place to return errors.
 *     \returns TRUE if sucessful.
 */
int MFAUTOGetStartPoint(MFImplicitMF M,MFAUTOTPBVP tpbvp, MFstpnt_type stpnt,doublereal *p0,MFNVector *u0,MFNKMatrix *Phi0,MFErrorHandler e);

/*!    \fn void MFAUTOAddUserZero(MFImplicitMF M,int parm, double value,MFErrorHandler e);
 *     \brief Adds a linear function of a single parameter whose level set indicates a curve on the solution that
 *               needs to be resolved.
 *
 *     \param M The solution manifold.
 *     \param parm The number of the parameter (in the list of all parameters, not just the ones in icp).
 *     \param value The level at which the level set is being tagged.
 *  \param e A place to return errors.
 */
void MFAUTOAddUserZero(MFImplicitMF M,int parm, double value,MFErrorHandler e);

/*!    \fn void MFAUTODetectLimitPoints(MFImplicitMF M,MFErrorHandler e);
 *     \brief Indicates that limit points are to be found (places the partial derivative of the solution manifold with respect
 *             to a single parameter is zero). The default is not to locate these submanifolds.
 *
 *     \param M The solution manifold.
 *  \param e A place to return errors.
 */
void MFAUTODetectLimitPoints(MFImplicitMF M,MFErrorHandler e);

/*!    \fn void MFAUTODetectBifurcationPoints(MFImplicitMF M,MFErrorHandler e);
 *     \brief Indicates that bifurcation points are to be found (places the Jacobian of the boundary value problem has 
 *               a null space with odd dimension (usually 1). The default is not to monitor these.
 *
 *     \param M The solution manifold.
 *  \param e A place to return errors.
 */
void MFAUTODetectBifurcationPoints(MFImplicitMF M,MFErrorHandler e);

/*!   \fn void MFAUTODetectSpecialPoints(MFImplicitMF M,int skip, MFErrorHandler e);
 *     \brief Indicates that special bifurcation points are to be found (for example, Hopf bifurcations).
 *               The default is not to detect.
 *
 *     \param M The solution manifold.
 *     \param skip Indicates that skip multipliers are ignored. 
 *  \param e A place to return errors.
 */
void MFAUTODetectSpecialPoints(MFImplicitMF M,int skip, MFErrorHandler e);

/*!    \fn MFImplicitMF MFCreateAUTOPeriodicSolution(MFAUTOTPBVP tpbvp, MFNSpace space,MFErrorHandler e);
 *     \brief Creates an implicit representation of the solution manifold of periodic motions (closed curves in phase space).
 *
 *     \param tpbvp A two point boundary value problem of the type used by AUTO.
 *     \param space The space in which the solution manifold lives.
 *  \param e A place to return errors.
 *     \returns The solution manifold. 
 */
MFImplicitMF MFCreateAUTOPeriodicSolution(MFAUTOTPBVP tpbvp, MFNSpace space,MFErrorHandler e);

/*!    \fn int MFAUTOBVSetIntegerParameter(MFImplicitMF M, char *parameterName, int value,MFErrorHandler e);
 *     \brief Allows the user to set AUTO's integer parameters. These are usually read from a file. Instead, default values
 *            are used when the AUTOBV is created, and the user may change them.
 *
 * Legal integer parameter names. See AUTO documentation for descriptions.
 *  <ul>
 *    <li>                  ips
 *    <li>                  irs
 *    <li>                  ilp
 *    <li>                  ntst
 *    <li>                  ncol
 *    <li>                  iad
 *    <li>                  iads
 *    <li>                  isp
 *    <li>                  isw
 *    <li>                  iplt
 *    <li>                  nbc
 *    <li>                  nint
 *    <li>                  nalc
 *    <li>                  nmx
 *    <li>                  nuzr
 *    <li>                  npr
 *    <li>                  mxbf
 *    <li>                  iid
 *    <li>                  itmx
 *    <li>                  itnw
 *    <li>                  nwtn
 *    <li>                  jac
 *    <li>                  iuzr
 *    <li>                  itp
 *    <li>                  itpst
 *    <li>                  ibr
 *    <li>                  nit
 *    <li>                  ntot
 *    <li>                  nins
 *    <li>                  istop
 *    <li>                  nbif
 *    <li>                  ipos
 *    <li>                  lab;
 *    <li>                  mynode
 *    <li>                  numnodes
 *    <li>                  parallel_flag
 *  </ul>
 *     \param M An AUTOBV solution manifold.
 *     \param parameterName Which parameter to set.
 *     \param value The new value.
 *  \param e A place to return errors.
 *     \returns FALSE if the parameter name does not match a parameter.
 */
int MFAUTOBVSetIntegerParameter(MFImplicitMF M, char *parameterName, int value,MFErrorHandler e);

/*!    \fn int MFAUTOBVSetRealParameter(MFImplicitMF M, char *parameterName, double value,MFErrorHandler e);
 *     \brief Allows the user to set AUTO's real valued parameters. These are usually read from a file. Instead, default values
 *            are used when the AUTOBV is created, and the user may change them.
 *
 * Legal real parameter names. See AUTO documentation for descriptions.
 *    <ul>
 *     <li>                 ds
 *     <li>                 dsmin
 *     <li>                 dsmax
 *     <li>                 amp
 *     <li>                 epsl
 *     <li>                 epsu
 *     <li>                 epss
 *     <li>                 det
 *     <li>                 tivp
 *     <li>                 fldf
 *     <li>                 hbff
 *     <li>                 biff
 *     <li>                 spbf
 *    </ul>
 *     \param M An AUTOBV solution manifold.
 *     \param parameterName Which parameter to set.
 *     \param value The new value.
 *  \param e A place to return errors.
 *     \returns FALSE if the parameter name does not match a parameter.
 */
int MFAUTOBVSetRealParameter(MFImplicitMF M, char *parameterName, double value,MFErrorHandler e);

/*!    \fn int MFAUTOBVGetIntegerParameter(MFImplicitMF M, char *parameterName,MFErrorHandler e);
 *     \brief Allows the user to set AUTO's integer parameters. These are usually read from a file. Instead, default values
 *            are used when the AUTOBV is created, and the user may change them.
 *
 * Legal integer parameter names. See AUTO documentation for descriptions.
 *  <ul>
 *    <li>                  ips
 *    <li>                  irs
 *    <li>                  ilp
 *    <li>                  ntst
 *    <li>                  ncol
 *    <li>                  iad
 *    <li>                  iads
 *    <li>                  isp
 *    <li>                  isw
 *    <li>                  iplt
 *    <li>                  nbc
 *    <li>                  nint
 *    <li>                  nalc
 *    <li>                  nmx
 *    <li>                  nuzr
 *    <li>                  npr
 *    <li>                  mxbf
 *    <li>                  iid
 *    <li>                  itmx
 *    <li>                  itnw
 *    <li>                  nwtn
 *    <li>                  jac
 *    <li>                  iuzr
 *    <li>                  itp
 *    <li>                  itpst
 *    <li>                  ibr
 *    <li>                  nit
 *    <li>                  ntot
 *    <li>                  nins
 *    <li>                  istop
 *    <li>                  nbif
 *    <li>                  ipos
 *    <li>                  lab;
 *    <li>                  mynode
 *    <li>                  numnodes
 *    <li>                  parallel_flag
 *  </ul>
 *     \param M An AUTOBV solution manifold.
 *     \param parameterName Which parameter value to retreive. A warning is issued if the parameter name does not match a parameter.
 *  \param e A place to return errors.
 *     \returns The current value of the parameter.
 */
int MFAUTOBVGetIntegerParameter(MFImplicitMF M, char *parameterName,MFErrorHandler e);

/*!    \fn double MFAUTOBVGetRealParameter(MFImplicitMF M, char *parameterName,MFErrorHandler e);
 *     \brief Allows the user to set AUTO's real valued parameters. These are usually read from a file. Instead, default values
 *            are used when the AUTOBV is created, and the user may change them.
 *
 * Legal real parameter names. See AUTO documentation for descriptions.
 *    <ul>
 *     <li>                 ds
 *     <li>                 dsmin
 *     <li>                 dsmax
 *     <li>                 amp
 *     <li>                 epsl
 *     <li>                 epsu
 *     <li>                 epss
 *     <li>                 det
 *     <li>                 tivp
 *     <li>                 fldf
 *     <li>                 hbff
 *     <li>                 biff
 *     <li>                 spbf
 *    </ul>
 *     \param M An AUTOBV solution manifold.
 *     \param parameterName Which parameter value to retreive. A warning is issued if the parameter name does not match a parameter.
 *  \param e A place to return errors.
 *     \returns The current value of the parameter.
 */
double MFAUTOBVGetRealParameter(MFImplicitMF M, char *parameterName,MFErrorHandler e);

/*!    @} */

/*!    @} */

/*! \addtogroup MFAUTOTPBVP
 *  @{
 */

/*! \class MFAUTOTPBVP MFAUTO.h MFAUTO.h
 *  \brief A boundary value problem posed in AUTO's representation.
 *
 */

/*!    \fn MFAUTOTPBVP MFCreateAUTOTPBVP(integer k,integer ndim,MFfunc_type func,integer jac,
 *                            integer nbc,MFbcnd_type bcnd,
 *                            integer nic,MFicnd_type icnd,
 *                            integer npar,integer nicp, integer *icp,integer ntst,integer ncol,
 *                            MFpvls_type pvls,MFErrorHandler e);
 *     \brief Creates a representation of the type of two point boundary value problem used by AUTO.
 *
 *     \param k The dimension of the solution manifold.
 *     \param ndim The dimension of the phase space.
 *     \param func A routine to evaluate a flow field defined over the phase space.
 *     \param jac TRUE if the derivatives are also available through calls to func. 
 *     \param nbc The number of boundary conditions applied to a trajectory in phase space.
 *     \param bcnd A routine to evaluate the residual of the boundary conditions.
 *     \param nic The number of integral conditions applied to a trajectory in phase space.
 *     \param icnd A routine to evaluate the kernal of the integral conditions.
 *     \param npar The total number of parameters on which the flow field depends.
 *     \param nicp The number of parameters which AUTO can use to follow paths of special solutions (e.g. bifurcations).
 *     \param icp an array of length at least nicp which lists the available parameters.
 *     \param ntst The number of mesh intervals to use to discretize the trajectory.
 *     \param ncol The number of collocation points to use on wach mesh interval.
 *     \param pvls OK, I haven't figured out what this one does. I suspect I don't have to require it for BVPs
 *  \param e A place to return errors.
 *     \returns A representation of the boundary value problem.
 */
MFAUTOTPBVP MFCreateAUTOTPBVP(integer k,integer ndim,MFfunc_type func,integer jac,
                              integer nbc,MFbcnd_type bcnd,
                              integer nic,MFicnd_type icnd,
                              integer npar,integer nicp, integer *icp,integer ntst,integer ncol,
                              MFpvls_type pvls,MFErrorHandler e);

/*! \fn void MFWriteAUTOTPBVP(FILE* fid,MFAUTOTPBVP tpbvp,MFErrorHandler e);
 *  \brief Writes a boundary value problem to a file.
 *
 *  \param fid The file to write to.
 *  \param tpbvp The boundary value problem being queried.
 *  \param e A place to return errors.
 */
void MFWriteAUTOTPBVP(FILE* fid,MFAUTOTPBVP tpbvp,MFErrorHandler e);


/*! \fn MFAUTOTPBVP MFReadAUTOTPBVP(FILE* fid,MFErrorHandler e);
 *  \brief Reads a boundary value problem from a file.
 *
 *  \param fid The file to write to.
 *  \param e A place to return errors.
 *  \returns tpbvp The boundary value problem.
 */
MFAUTOTPBVP MFReadAUTOTPBVP(FILE* fid,MFErrorHandler e);


/*! \fn void MFRefAUTOTPBVP(MFAUTOTPBVP tpbvp,MFErrorHandler e);
 *  \brief Adds a reference to a boundary value problem.
 *
 *  \param tpbvp The boundary value problem being referenced.
 *  \param e A place to return errors.
 *  \sa ReferenceCounting MFFreeAUTOTPBVP
 */
void MFRefAUTOTPBVP(MFAUTOTPBVP tpbvp,MFErrorHandler e);

/*! \fn void MFFreeAUTOTPBVP(MFAUTOTPBVP tpbvp,MFErrorHandler e);
 *  \brief Frees a reference to a boundary value problem, and deletes the chart if there are no references left.
 *
 *  \param tpbvp The boundary value problem being unreferenced.
 *  \param e A place to return errors.
 *  \sa ReferenceCounting MFRefAUTOTPBVP
 */
void MFFreeAUTOTPBVP(MFAUTOTPBVP tpbvp,MFErrorHandler e);

/*!    @} */
