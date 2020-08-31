#pragma once

#include "_f2c.h"

// C D S Z A

// C
//
#pragma region CBLAS1 -- Single precision complex BLAS routines
//
//set(CBLAS1 icamax.c scabs1.c scasum.c scnrm2.c
//	caxpy.c ccopy.c	cdotc.c cdotu.c csscal.c crotg.c cscal.c cswap.c csrot.c)

CBLAS_API
integer icamax_(integer* n, complex* cx, integer* incx);

CBLAS_API
doublereal scabs1_(complex* z__);

CBLAS_API
doublereal scasum_(integer* n, complex* cx, integer* incx);

CBLAS_API
doublereal scnrm2_(integer* n, complex* x, integer* incx);


/* Subroutine */
CBLAS_API
int caxpy_(integer* n, complex* ca,
	complex* cx, integer* incx, complex* cy, integer* incy);

/* Subroutine */
CBLAS_API
int ccopy_(integer* n,
	complex* cx, integer* incx, complex* cy, integer* incy);

/* Complex */
CBLAS_API
VOID cdotc_(complex* ret_val, integer* n,
	complex* cx, integer* incx, complex* cy, integer* incy);

/* Complex */
CBLAS_API
VOID cdotu_(complex* ret_val, integer* n,
	complex* cx, integer* incx, complex* cy, integer* incy);

/* Subroutine */
CBLAS_API
int csscal_(integer* n, real* sa, complex* cx, integer* incx);

/* Subroutine */
CBLAS_API
int crotg_(complex* ca, complex* cb, real* c__, complex* s);

/* Subroutine */
CBLAS_API
int cscal_(integer* n, complex* ca, complex* cx, integer* incx);

/* Subroutine */
CBLAS_API
int cswap_(integer* n, complex* cx, integer* incx, complex* cy, integer* incy);

/* Subroutine */
CBLAS_API
int csrot_(integer* n,
	complex* cx, integer* incx,
	complex* cy, integer* incy, real* c__, real* s);


#pragma endregion

#pragma region CBLAS2 -- Single precision complex BLAS2 routines
//
//set(CBLAS2 cgbmv.c cgemv.c cgerc.c cgeru.c 
//	chbmv.c chemv.c cher.c cher2.c chpmv.c chpr.c chpr2.c
//	ctbmv.c ctpmv.c ctrmv.c ctbsv.c ctpsv.c ctrsv.c)

/* Subroutine */
CBLAS_API
int cgbmv_(char* trans, integer* m, integer* n, integer* kl,
	integer* ku, complex* alpha, complex* a, integer* lda, complex* x,
	integer* incx, complex* beta, complex* y, integer* incy);

/* Subroutine */
CBLAS_API
int cgemv_(char* trans, integer* m, integer* n, complex*
	alpha, complex* a, integer* lda, complex* x, integer* incx, complex*
	beta, complex* y, integer* incy);

/* Subroutine */
CBLAS_API
int cgerc_(integer* m, integer* n, complex* alpha, complex*
	x, integer* incx, complex* y, integer* incy, complex* a, integer* lda);

/* Subroutine */
CBLAS_API
int cgeru_(integer* m, integer* n, complex* alpha, complex*
	x, integer* incx, complex* y, integer* incy, complex* a, integer* lda);

/* Subroutine */
CBLAS_API
int chbmv_(char* uplo, integer* n, integer* k, complex*
	alpha, complex* a, integer* lda, complex* x, integer* incx, complex*
	beta, complex* y, integer* incy);

/* Subroutine */
CBLAS_API
int cher_(char* uplo, integer* n, real* alpha, complex* x,
	integer* incx, complex* a, integer* lda);

/* Subroutine */
CBLAS_API
int cher2_(char* uplo, integer* n, complex* alpha, complex*
	x, integer* incx, complex* y, integer* incy, complex* a, integer* lda);

/* Subroutine */
CBLAS_API
int chemv_(char* uplo, integer* n, complex* alpha, complex*
	a, integer* lda, complex* x, integer* incx, complex* beta, complex* y,
	integer* incy);

/* Subroutine */
CBLAS_API
int chpmv_(char* uplo, integer* n, complex* alpha, complex*
	ap, complex* x, integer* incx, complex* beta, complex* y, integer*
	incy);

/* Subroutine */
CBLAS_API
int chpr_(char* uplo, integer* n, real* alpha, complex* x,
	integer* incx, complex* ap);

/* Subroutine */
CBLAS_API
int chpr2_(char* uplo, integer* n, complex* alpha, complex*
	x, integer* incx, complex* y, integer* incy, complex* ap);

/* Subroutine */
CBLAS_API
int ctbmv_(char* uplo, char* trans, char* diag, integer* n,
	integer* k, complex* a, integer* lda, complex* x, integer* incx);

/* Subroutine */
CBLAS_API
int ctpmv_(char* uplo, char* trans, char* diag, integer* n,
	complex* ap, complex* x, integer* incx);

/* Subroutine */
CBLAS_API
int ctrmv_(char* uplo, char* trans, char* diag, integer* n,
	complex* a, integer* lda, complex* x, integer* incx);

/* Subroutine */
CBLAS_API
int ctbsv_(char* uplo, char* trans, char* diag, integer* n,
	integer* k, complex* a, integer* lda, complex* x, integer* incx);

/* Subroutine */
CBLAS_API
int ctpsv_(char* uplo, char* trans, char* diag, integer* n,
	complex* ap, complex* x, integer* incx);

/* Subroutine */
CBLAS_API
int ctrsv_(char* uplo, char* trans, char* diag, integer* n,
	complex* a, integer* lda, complex* x, integer* incx);




#pragma endregion

#pragma region CBLAS3 -- Single precision complex BLAS3 routines
//
//set(CBLAS3 cgemm.c chemm.c cher2k.c cherk.c
//	csymm.c csyr2k.c csyrk.c ctrmm.c ctrsm.c)

/* Subroutine */
CBLAS_API
int cgemm_(char* transa, char* transb, integer* m, integer*
	n, integer* k, complex* alpha, complex* a, integer* lda, complex* b,
	integer* ldb, complex* beta, complex* c__, integer* ldc);

/* Subroutine */
CBLAS_API
int chemm_(char* side, char* uplo, integer* m, integer* n,
	complex* alpha, complex* a, integer* lda, complex* b, integer* ldb,
	complex* beta, complex* c__, integer* ldc);

/* Subroutine */
CBLAS_API
int cher2k_(char* uplo, char* trans, integer* n, integer* k,
	complex* alpha, complex* a, integer* lda, complex* b, integer* ldb,
	real* beta, complex* c__, integer* ldc);

/* Subroutine */
CBLAS_API
int cherk_(char* uplo, char* trans, integer* n, integer* k,
	real* alpha, complex* a, integer* lda, real* beta, complex* c__,
	integer* ldc);

/* Subroutine */
CBLAS_API
int csymm_(char* side, char* uplo, integer* m, integer* n,
	complex* alpha, complex* a, integer* lda, complex* b, integer* ldb,
	complex* beta, complex* c__, integer* ldc);

/* Subroutine */
CBLAS_API
int csyr2k_(char* uplo, char* trans, integer* n, integer* k,
	complex* alpha, complex* a, integer* lda, complex* b, integer* ldb,
	complex* beta, complex* c__, integer* ldc);

/* Subroutine */
CBLAS_API
int csyrk_(char* uplo, char* trans, integer* n, integer* k,
	complex* alpha, complex* a, integer* lda, complex* beta, complex* c__,
	integer* ldc);

/* Subroutine */
CBLAS_API
int ctrmm_(char* side, char* uplo, char* transa, char* diag,
	integer* m, integer* n, complex* alpha, complex* a, integer* lda,
	complex* b, integer* ldb);

/* Subroutine */
CBLAS_API
int ctrsm_(char* side, char* uplo, char* transa, char* diag,
	integer* m, integer* n, complex* alpha, complex* a, integer* lda,
	complex* b, integer* ldb);


#pragma endregion



// D
//
#pragma region DBLAS1 -- Double precision real BLAS routines
//
//set(DBLAS1 idamax.c dasum.c daxpy.c dcopy.c ddot.c dnrm2.c
//	drot.c drotg.c dscal.c dsdot.c dswap.c drotm.c drotmg.c)

CBLAS_API
integer idamax_(integer* n, doublereal* dx, integer* incx);

CBLAS_API
doublereal dasum_(integer* n, doublereal* dx, integer* incx);

/* Subroutine */
CBLAS_API
int daxpy_(integer* n, doublereal* da,
	doublereal* dx, integer* incx,
	doublereal* dy, integer* incy);


/* Subroutine */
CBLAS_API
int dcopy_(integer* n,
	doublereal* dx, integer* incx,
	doublereal* dy, integer* incy);

CBLAS_API
doublereal ddot_(integer* n,
	doublereal* dx, integer* incx,
	doublereal* dy, integer* incy);


CBLAS_API
doublereal dnrm2_(integer* n,
	doublereal* x, integer* incx);


/* Subroutine */
CBLAS_API
int drot_(integer* n,
	doublereal* dx, integer* incx,
	doublereal* dy, integer* incy,
	doublereal* c__, doublereal* s);


/* Subroutine */
CBLAS_API
int drotg_(doublereal* da, doublereal* db, doublereal* c__, doublereal* s);

/* Subroutine */
CBLAS_API
int dscal_(integer* n, doublereal* da,
	doublereal* dx, integer* incx);

CBLAS_API
doublereal dsdot_(integer* n,
	real* sx, integer* incx,
	real* sy, integer* incy);


/* Subroutine */
CBLAS_API
int dswap_(integer* n,
	doublereal* dx, integer* incx,
	doublereal* dy, integer* incy);

/* Subroutine */
CBLAS_API
int drotm_(integer* n,
	doublereal* dx, integer* incx,
	doublereal* dy, integer* incy, doublereal* dparam);

/* Subroutine */
CBLAS_API
int drotmg_(doublereal* dd1, doublereal* dd2,
	doublereal* dx1, doublereal* dy1, doublereal* dparam);

#pragma endregion

#pragma region DBLAS2 -- Double precision real BLAS2 routines
//
//set(DBLAS2 dgbmv.c dgemv.c dger.c 
//	dsbmv.c dspmv.c dspr.c dspr2.c dsymv.c dsyr.c  dsyr2.c
//	dtbmv.c dtbsv.c dtpmv.c dtpsv.c dtrmv.c dtrsv.c	)

/* Subroutine */
CBLAS_API
int dgbmv_(char* trans, integer* m, integer* n, integer* kl,
	integer* ku, doublereal* alpha, doublereal* a, integer* lda,
	doublereal* x, integer* incx, doublereal* beta, doublereal* y,
	integer* incy);

/* Subroutine */
CBLAS_API
int dgemv_(char* trans, integer* m, integer* n, doublereal*
	alpha, doublereal* a, integer* lda, doublereal* x, integer* incx,
	doublereal* beta, doublereal* y, integer* incy);

/* Subroutine */
CBLAS_API
int dger_(integer* m, integer* n, doublereal* alpha,
	doublereal* x, integer* incx, doublereal* y, integer* incy,
	doublereal* a, integer* lda);

/* Subroutine */
CBLAS_API
int dsbmv_(char* uplo, integer* n, integer* k, doublereal*
	alpha, doublereal* a, integer* lda, doublereal* x, integer* incx,
	doublereal* beta, doublereal* y, integer* incy);

/* Subroutine */
CBLAS_API
int dspmv_(char* uplo, integer* n, doublereal* alpha,
	doublereal* ap, doublereal* x, integer* incx, doublereal* beta,
	doublereal* y, integer* incy);

/* Subroutine */
CBLAS_API
int dspr_(char* uplo, integer* n, doublereal* alpha,
	doublereal* x, integer* incx, doublereal* ap);

/* Subroutine */
CBLAS_API
int dspr2_(char* uplo, integer* n, doublereal* alpha,
	doublereal* x, integer* incx, doublereal* y, integer* incy,
	doublereal* ap);

/* Subroutine */
CBLAS_API
int dsymv_(char* uplo, integer* n, doublereal* alpha,
	doublereal* a, integer* lda, doublereal* x, integer* incx, doublereal
	* beta, doublereal* y, integer* incy);

/* Subroutine */
CBLAS_API
int dsyr_(char* uplo, integer* n, doublereal* alpha,
	doublereal* x, integer* incx, doublereal* a, integer* lda);

/* Subroutine */
CBLAS_API
int dsyr2_(char* uplo, integer* n, doublereal* alpha,
	doublereal* x, integer* incx, doublereal* y, integer* incy,
	doublereal* a, integer* lda);

/* Subroutine */
CBLAS_API
int dtbmv_(char* uplo, char* trans, char* diag, integer* n,
	integer* k, doublereal* a, integer* lda, doublereal* x, integer* incx);

/* Subroutine */
CBLAS_API
int dtbsv_(char* uplo, char* trans, char* diag, integer* n,
	integer* k, doublereal* a, integer* lda, doublereal* x, integer* incx);

/* Subroutine */
CBLAS_API
int dtpmv_(char* uplo, char* trans, char* diag, integer* n,
	doublereal* ap, doublereal* x, integer* incx);

/* Subroutine */
CBLAS_API
int dtpsv_(char* uplo, char* trans, char* diag, integer* n,
	doublereal* ap, doublereal* x, integer* incx);

/* Subroutine */
CBLAS_API
int dtrmv_(char* uplo, char* trans, char* diag, integer* n,
	doublereal* a, integer* lda, doublereal* x, integer* incx);

/* Subroutine */
CBLAS_API
int dtrsv_(char* uplo, char* trans, char* diag, integer* n,
	doublereal* a, integer* lda, doublereal* x, integer* incx);


#pragma endregion

#pragma region DBLAS3 -- Double precision real BLAS3 routines
//
//set(DBLAS3 dgemm.c dsymm.c dsyr2k.c dsyrk.c dtrmm.c dtrsm.c)


/* Subroutine */
CBLAS_API
int dgemm_(char* transa, char* transb, integer* m, integer*
	n, integer* k, doublereal* alpha, doublereal* a, integer* lda,
	doublereal* b, integer* ldb, doublereal* beta, doublereal* c__,
	integer* ldc);

/* Subroutine */
CBLAS_API
int dsymm_(char* side, char* uplo, integer* m, integer* n,
	doublereal* alpha, doublereal* a, integer* lda, doublereal* b,
	integer* ldb, doublereal* beta, doublereal* c__, integer* ldc);

/* Subroutine */
CBLAS_API
int dsyr2k_(char* uplo, char* trans, integer* n, integer* k,
	doublereal* alpha, doublereal* a, integer* lda, doublereal* b,
	integer* ldb, doublereal* beta, doublereal* c__, integer* ldc);

/* Subroutine */
CBLAS_API
int dsyrk_(char* uplo, char* trans, integer* n, integer* k,
	doublereal* alpha, doublereal* a, integer* lda, doublereal* beta,
	doublereal* c__, integer* ldc);

/* Subroutine */
CBLAS_API
int dtrmm_(char* side, char* uplo, char* transa, char* diag,
	integer* m, integer* n, doublereal* alpha, doublereal* a, integer*
	lda, doublereal* b, integer* ldb);

/* Subroutine */
CBLAS_API
int dtrsm_(char* side, char* uplo, char* transa, char* diag,
	integer* m, integer* n, doublereal* alpha, doublereal* a, integer*
	lda, doublereal* b, integer* ldb);


#pragma endregion



// S
//
#pragma region SBLAS1 -- Single precision real BLAS routines
//
//set(SBLAS1 isamax.c sasum.c saxpy.c scopy.c sdot.c snrm2.c
//	srot.c srotg.c sscal.c sswap.c sdsdot.c srotm.c srotmg.c)

CBLAS_API
integer isamax_(integer* n, real* sx, integer* incx);

CBLAS_API
doublereal sasum_(integer* n, real* sx, integer* incx);

/* Subroutine */
CBLAS_API
int saxpy_(integer* n, real* sa,
	real* sx, integer* incx, real* sy, integer* incy);

/* Subroutine */
CBLAS_API
int scopy_(integer* n,
	real* sx, integer* incx, real* sy, integer* incy);

CBLAS_API
doublereal sdot_(integer* n,
	real* sx, integer* incx, real* sy, integer* incy);

CBLAS_API
doublereal snrm2_(integer* n, real* x, integer* incx);


/* Subroutine */
CBLAS_API
int srot_(integer* n, real* sx, integer* incx,
	real* sy, integer* incy, real* c__, real* s);

/* Subroutine */
CBLAS_API
int srotg_(real* sa, real* sb, real* c__, real* s);

/* Subroutine */
CBLAS_API
int sscal_(integer* n, real* sa, real* sx, integer* incx);

/* Subroutine */
CBLAS_API
int sswap_(integer* n, real* sx, integer* incx, real* sy, integer* incy);

CBLAS_API
doublereal sdsdot_(integer* n, real* sb,
	real* sx, integer* incx, real* sy,	integer* incy);

/* Subroutine */
CBLAS_API
int srotm_(integer* n, real* sx, integer* incx, real* sy, integer* incy, real* sparam);

/* Subroutine */
CBLAS_API
int srotmg_(real* sd1, real* sd2, real* sx1, real* sy1, real* sparam);



#pragma endregion

#pragma region SBLAS2 -- Single precision real BLAS2 routines
//
//set(SBLAS2
//	sgbmv.c sgemv.c sger.c 
//	ssbmv.c sspmv.c sspr.c sspr2.c ssymv.c ssyr.c ssyr2.c
//	stbmv.c stbsv.c stpmv.c stpsv.c strmv.c strsv.c)

/* Subroutine */
CBLAS_API
int sgbmv_(char* trans, integer* m, integer* n, integer* kl,
	integer* ku, real* alpha, real* a, integer* lda, real* x, integer*
	incx, real* beta, real* y, integer* incy);

/* Subroutine */
CBLAS_API
int sgemv_(char* trans, integer* m, integer* n, real* alpha,
	real* a, integer* lda, real* x, integer* incx, real* beta, real* y,
	integer* incy);

/* Subroutine */
CBLAS_API
int sger_(integer* m, integer* n, real* alpha, real* x,
	integer* incx, real* y, integer* incy, real* a, integer* lda);

/* Subroutine */
CBLAS_API
int ssbmv_(char* uplo, integer* n, integer* k, real* alpha,
	real* a, integer* lda, real* x, integer* incx, real* beta, real* y,
	integer* incy);

/* Subroutine */
CBLAS_API
int sspmv_(char* uplo, integer* n, real* alpha, real* ap,
	real* x, integer* incx, real* beta, real* y, integer* incy);

/* Subroutine */
CBLAS_API
int sspr_(char* uplo, integer* n, real* alpha, real* x,
	integer* incx, real* ap);

/* Subroutine */
CBLAS_API
int sspr2_(char* uplo, integer* n, real* alpha, real* x,
	integer* incx, real* y, integer* incy, real* ap);

/* Subroutine */
CBLAS_API
int ssymv_(char* uplo, integer* n, real* alpha, real* a,
	integer* lda, real* x, integer* incx, real* beta, real* y, integer*
	incy);

/* Subroutine */
CBLAS_API
int ssyr_(char* uplo, integer* n, real* alpha, real* x,
	integer* incx, real* a, integer* lda);

/* Subroutine */
CBLAS_API
int ssyr2_(char* uplo, integer* n, real* alpha, real* x,
	integer* incx, real* y, integer* incy, real* a, integer* lda);

/* Subroutine */
CBLAS_API
int stbmv_(char* uplo, char* trans, char* diag, integer* n,
	integer* k, real* a, integer* lda, real* x, integer* incx);

/* Subroutine */
CBLAS_API
int stbsv_(char* uplo, char* trans, char* diag, integer* n,
	integer* k, real* a, integer* lda, real* x, integer* incx);

/* Subroutine */
CBLAS_API
int stpmv_(char* uplo, char* trans, char* diag, integer* n,
	real* ap, real* x, integer* incx);

/* Subroutine */
CBLAS_API
int stpsv_(char* uplo, char* trans, char* diag, integer* n,
	real* ap, real* x, integer* incx);

/* Subroutine */
CBLAS_API
int strmv_(char* uplo, char* trans, char* diag, integer* n,
	real* a, integer* lda, real* x, integer* incx);

/* Subroutine */
CBLAS_API
int strsv_(char* uplo, char* trans, char* diag, integer* n,
	real* a, integer* lda, real* x, integer* incx);


#pragma endregion

#pragma region SBLAS3 -- Single precision real BLAS3 routines
//
//set(SBLAS3 sgemm.c ssymm.c ssyr2k.c ssyrk.c strmm.c strsm.c )

/* Subroutine */
CBLAS_API
int sgemm_(char* transa, char* transb, integer* m, integer*
	n, integer* k, real* alpha, real* a, integer* lda, real* b, integer*
	ldb, real* beta, real* c__, integer* ldc);

/* Subroutine */
CBLAS_API
int ssymm_(char* side, char* uplo, integer* m, integer* n,
	real* alpha, real* a, integer* lda, real* b, integer* ldb, real* beta,
	real* c__, integer* ldc);

/* Subroutine */
CBLAS_API
int ssyr2k_(char* uplo, char* trans, integer* n, integer* k,
	real* alpha, real* a, integer* lda, real* b, integer* ldb, real* beta,
	real* c__, integer* ldc);

/* Subroutine */
CBLAS_API
int ssyrk_(char* uplo, char* trans, integer* n, integer* k,
	real* alpha, real* a, integer* lda, real* beta, real* c__, integer*
	ldc);

/* Subroutine */
CBLAS_API
int strmm_(char* side, char* uplo, char* transa, char* diag,
	integer* m, integer* n, real* alpha, real* a, integer* lda, real* b,
	integer* ldb);

/* Subroutine */
CBLAS_API
int strsm_(char* side, char* uplo, char* transa, char* diag,
	integer* m, integer* n, real* alpha, real* a, integer* lda, real* b,
	integer* ldb);





#pragma endregion


// Z
//
#pragma region ZBLAS1 -- Double precision complex BLAS routines
//
//set(ZBLAS1 izamax.c dcabs1.c dzasum.c dznrm2.c
//	zaxpy.c zcopy.c	zdotc.c zdotu.c zdrot.c zdscal.c zrotg.c zscal.c zswap.c)

CBLAS_API
integer izamax_(integer* n, doublecomplex* zx, integer* incx);

CBLAS_API
doublereal dcabs1_(doublecomplex* z__);

CBLAS_API
doublereal dzasum_(integer* n, doublecomplex* zx, integer* incx);

CBLAS_API
doublereal dznrm2_(integer* n, doublecomplex* x, integer* incx);

/* Subroutine */
CBLAS_API
int zaxpy_(integer* n, doublecomplex* za, doublecomplex* zx,
	integer* incx, doublecomplex* zy, integer* incy);

/* Subroutine */
CBLAS_API
int zcopy_(integer* n, doublecomplex* zx, integer* incx,
	doublecomplex* zy, integer* incy);

/* Double Complex */
CBLAS_API
VOID zdotc_(doublecomplex* ret_val, integer* n,
	doublecomplex* zx, integer* incx, doublecomplex* zy, integer* incy);

/* Double Complex */
CBLAS_API
VOID zdotu_(doublecomplex* ret_val, integer* n,
	doublecomplex* zx, integer* incx, doublecomplex* zy, integer* incy);

/* Subroutine */
CBLAS_API
int zdrot_(integer* n, doublecomplex* cx, integer* incx,
	doublecomplex* cy, integer* incy, doublereal* c__, doublereal* s);

/* Subroutine */
CBLAS_API
int zdscal_(integer* n, doublereal* da, doublecomplex* zx,
	integer* incx);

/* Subroutine */
CBLAS_API
int zrotg_(doublecomplex* ca, doublecomplex* cb, doublereal*
	c__, doublecomplex* s);

/* Subroutine */
CBLAS_API
int zscal_(integer* n, doublecomplex* za, doublecomplex* zx,
	integer* incx);

/* Subroutine */
CBLAS_API
int zswap_(integer* n, doublecomplex* zx, integer* incx,
	doublecomplex* zy, integer* incy);


#pragma endregion

#pragma region ZBLAS2 -- Double precision complex BLAS2 routines
//
//set(ZBLAS2
//	zgbmv.c zgemv.c zgerc.c zgeru.c
//	zhbmv.c zhemv.c zher.c zher2.c  zhpmv.c zhpr.c zhpr2.c
//	ztbmv.c ztbsv.c ztpmv.c ztpsv.c ztrmv.c ztrsv.c)

/* Subroutine */
CBLAS_API
int zgbmv_(char* trans, integer* m, integer* n, integer* kl,
	integer* ku, doublecomplex* alpha, doublecomplex* a, integer* lda,
	doublecomplex* x, integer* incx, doublecomplex* beta, doublecomplex*
	y, integer* incy);

/* Subroutine */
CBLAS_API
int zgemv_(char* trans, integer* m, integer* n,
	doublecomplex* alpha, doublecomplex* a, integer* lda, doublecomplex*
	x, integer* incx, doublecomplex* beta, doublecomplex* y, integer*
	incy);

/* Subroutine */
CBLAS_API
int zgerc_(integer* m, integer* n, doublecomplex* alpha,
	doublecomplex* x, integer* incx, doublecomplex* y, integer* incy,
	doublecomplex* a, integer* lda);

/* Subroutine */
CBLAS_API
int zgeru_(integer* m, integer* n, doublecomplex* alpha,
	doublecomplex* x, integer* incx, doublecomplex* y, integer* incy,
	doublecomplex* a, integer* lda);

/* Subroutine */
CBLAS_API
int zhbmv_(char* uplo, integer* n, integer* k, doublecomplex
	* alpha, doublecomplex* a, integer* lda, doublecomplex* x, integer*
	incx, doublecomplex* beta, doublecomplex* y, integer* incy);

/* Subroutine */
CBLAS_API
int zhemv_(char* uplo, integer* n, doublecomplex* alpha,
	doublecomplex* a, integer* lda, doublecomplex* x, integer* incx,
	doublecomplex* beta, doublecomplex* y, integer* incy);

/* Subroutine */
CBLAS_API
int zher_(char* uplo, integer* n, doublereal* alpha,
	doublecomplex* x, integer* incx, doublecomplex* a, integer* lda);

/* Subroutine */
CBLAS_API
int zher2_(char* uplo, integer* n, doublecomplex* alpha,
	doublecomplex* x, integer* incx, doublecomplex* y, integer* incy,
	doublecomplex* a, integer* lda);

/* Subroutine */
CBLAS_API
int zhpmv_(char* uplo, integer* n, doublecomplex* alpha,
	doublecomplex* ap, doublecomplex* x, integer* incx, doublecomplex*
	beta, doublecomplex* y, integer* incy);

/* Subroutine */
CBLAS_API
int zhpr_(char* uplo, integer* n, doublereal* alpha,
	doublecomplex* x, integer* incx, doublecomplex* ap);

/* Subroutine */
CBLAS_API
int zhpr2_(char* uplo, integer* n, doublecomplex* alpha,
	doublecomplex* x, integer* incx, doublecomplex* y, integer* incy,
	doublecomplex* ap);

/* Subroutine */
CBLAS_API
int ztbmv_(char* uplo, char* trans, char* diag, integer* n,
	integer* k, doublecomplex* a, integer* lda, doublecomplex* x, integer
	* incx);

/* Subroutine */
CBLAS_API
int ztbsv_(char* uplo, char* trans, char* diag, integer* n,
	integer* k, doublecomplex* a, integer* lda, doublecomplex* x, integer
	* incx);

/* Subroutine */
CBLAS_API
int ztpmv_(char* uplo, char* trans, char* diag, integer* n,
	doublecomplex* ap, doublecomplex* x, integer* incx);

/* Subroutine */
CBLAS_API
int ztpsv_(char* uplo, char* trans, char* diag, integer* n,
	doublecomplex* ap, doublecomplex* x, integer* incx);

/* Subroutine */
CBLAS_API
int ztrmv_(char* uplo, char* trans, char* diag, integer* n,
	doublecomplex* a, integer* lda, doublecomplex* x, integer* incx);

/* Subroutine */
CBLAS_API
int ztrsv_(char* uplo, char* trans, char* diag, integer* n,
	doublecomplex* a, integer* lda, doublecomplex* x, integer* incx);


#pragma endregion

#pragma region ZBLAS3 -- Double precision complex BLAS3 routines
//
//set(ZBLAS3
//	zgemm.c zhemm.c zher2k.c zherk.c
//	zsymm.c zsyr2k.c zsyrk.c ztrmm.c ztrsm.c)

/* Subroutine */
CBLAS_API
int zgemm_(char* transa, char* transb, integer* m, integer*
	n, integer* k, doublecomplex* alpha, doublecomplex* a, integer* lda,
	doublecomplex* b, integer* ldb, doublecomplex* beta, doublecomplex*
	c__, integer* ldc);

/* Subroutine */
CBLAS_API
int zhemm_(char* side, char* uplo, integer* m, integer* n,
	doublecomplex* alpha, doublecomplex* a, integer* lda, doublecomplex*
	b, integer* ldb, doublecomplex* beta, doublecomplex* c__, integer*
	ldc);

/* Subroutine */
CBLAS_API
int zher2k_(char* uplo, char* trans, integer* n, integer* k,
	doublecomplex* alpha, doublecomplex* a, integer* lda, doublecomplex*
	b, integer* ldb, doublereal* beta, doublecomplex* c__, integer* ldc);

/* Subroutine */
CBLAS_API
int zherk_(char* uplo, char* trans, integer* n, integer* k,
	doublereal* alpha, doublecomplex* a, integer* lda, doublereal* beta,
	doublecomplex* c__, integer* ldc);

/* Subroutine */
CBLAS_API
int zsymm_(char* side, char* uplo, integer* m, integer* n,
	doublecomplex* alpha, doublecomplex* a, integer* lda, doublecomplex*
	b, integer* ldb, doublecomplex* beta, doublecomplex* c__, integer*
	ldc);

/* Subroutine */
CBLAS_API
int zsyr2k_(char* uplo, char* trans, integer* n, integer* k,
	doublecomplex* alpha, doublecomplex* a, integer* lda, doublecomplex*
	b, integer* ldb, doublecomplex* beta, doublecomplex* c__, integer*
	ldc);

/* Subroutine */
CBLAS_API
int zsyrk_(char* uplo, char* trans, integer* n, integer* k,
	doublecomplex* alpha, doublecomplex* a, integer* lda, doublecomplex*
	beta, doublecomplex* c__, integer* ldc);

/* Subroutine */
CBLAS_API
int ztrmm_(char* side, char* uplo, char* transa, char* diag,
	integer* m, integer* n, doublecomplex* alpha, doublecomplex* a,
	integer* lda, doublecomplex* b, integer* ldb);

/* Subroutine */
CBLAS_API
int ztrsm_(char* side, char* uplo, char* transa, char* diag,
	integer* m, integer* n, doublecomplex* alpha, doublecomplex* a,
	integer* lda, doublecomplex* b, integer* ldb);


#pragma endregion

// Auxiliary
//
#pragma region ALLBLAS -- Auxiliary routines for Level 2 and 3 BLAS
//
//set(ALLBLAS  lsame.c xerbla.c xerbla_array.c)

CBLAS_API
logical lsame_(char* ca, char* cb);

/* Subroutine */
CBLAS_API
int xerbla_(char* srname, integer * info);

/* Subroutine */
CBLAS_API
int xerbla_array__(char* srname_array__, integer * srname_len__,
	integer * info, ftnlen srname_array_len);


#pragma endregion

