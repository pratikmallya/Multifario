#include <multifarioConfig.h>

#ifdef NEED_BLAS_WRAPPERS

void F77_FUNC(daxpy,DAXPY)()
 {
  FC_FUNC(daxpy,DAXPY)();
  return;
 }

void F77_FUNC(dscal,DSCAL)(int *n,double *da,double *dx,int *incx)                                       
 {
  FC_FUNC(dscal,DSCAL)(n,da,dx,incx)                                       
  return;
 }

void F77_FUNC(daxpy,DAXPY)(int *n,double *da,double *dx,int *incx, double *dy, int *incy)
 {
  FC_FUNC(daxpy,DAXPY)(n,da,dx,incx, dy, incy);
  return;
 }

void F77_FUNC(dgetrf,DGETRF)(int *m,int *n,double *a,int *lda, int *ipiv, int *info);
 {
  FC_FUNC(dgetrf,DGETRF)(m,n,a,lda, ipiv, info);
  return;
 }

void F77_FUNC(dlaswp,DLASWP)(int *n, double *a, int *lda, int *k1, int *k2, int *ipiv,int *incx)
 {
  FC_FUNC(dlaswp,DLASWP)(n, a, lda, k1, k2, ipiv,incx);
  return;
 }

void F77_FUNC(dtrsm,DTRSM)(char *side, char *uplo, char *transa, char *diag, int *n, int *m, double *a, int *lda, double *b, int *ldb)
 {
  FC_FUNC(dtrsm,DTRSM)(side, uplo, transa, diag, n, m, a, lda, b, ldb);
  return;
 }

void F77_FUNC(dgesvd,DGESVD)(char *jobu, char *jobvt, int *m, int *n, double *A, int *lda, double *s,double *u, int *ldu, double *vt, int *ldvt, double *work, int *lwork, int *info)
 {
  FC_FUNC(dgesvd,DGESVD)(jobu, jobvt, m, n, A, lda, s,u, ldu, vt, ldvt, work, lwork, info);
  return;
 }
#else

static int dummy()
 {
  return 0;
 }

#endif
