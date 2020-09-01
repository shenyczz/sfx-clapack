#pragma once

#include "_f2c.h"
#include "cblas.h"


CLAPACK_API int fn_clapack(void);

// variants - 变体
/// Auxiliary - 辅助

//	ala		- ALLAUX -- Auxiliary routines called from all precisions
///	aladz	- DZLAUX -- Auxiliary routines called from all precisions but only from routines using extra precision.
//	alasc	- SCLAUX -- Auxiliary routines called from both REALand COMPLEX

//	cla		- CLASRC -- Single precision complex LAPACK routines
///	cxla	- CXLASRC -- Single precision complex LAPACK routines using extra precision.

//	dla		- DLASRC -- Double precision real LAPACK routines
///	dxla	- DXLASRC -- Double precision real LAPACK routines using extra precision.

//	sla		- SLASRC -- Single precision real LAPACK routines
///	sxla	- SXLASRC -- Single precision real LAPACK routines using extra precision.

//	zla		- ZLASRC -- Double precision complex LAPACK routines
///	zxla	- ZXLASRC -- Double precision complex LAPACK routines using extra precision.



#pragma region ALLAUX -- Auxiliary routines called from all precisions (ala)

/* Character */
CLAPACK_API
VOID chla_transtype__(char* ret_val, ftnlen ret_val_len,
	integer* trans);

CLAPACK_API
integer ieeeck_(integer* ispec, real* zero, real* one);

CLAPACK_API
integer iladiag_(char* diag);

CLAPACK_API
integer ilaenv_(integer* ispec, char* name__, char* opts, integer* n1,
	integer* n2, integer* n3, integer* n4);

CLAPACK_API
integer ilaprec_(char* prec);

CLAPACK_API
integer ilatrans_(char* trans);

CLAPACK_API
integer ilauplo_(char* uplo);

/* Subroutine */
CLAPACK_API
int ilaver_(integer* vers_major__, integer* vers_minor__,
	integer* vers_patch__);

CLAPACK_API
integer iparmq_(integer* ispec, char* name__, char* opts, integer* n, integer
	* ilo, integer* ihi, integer* lwork);

CLAPACK_API
logical lsamen_(integer* n, char* ca, char* cb);


#pragma endregion



#pragma region CLASRC -- Single precision complex LAPACK routines (cla)

/* Subroutine */
CLAPACK_API
int cbdsqr_(char* uplo, integer* n, integer* ncvt, integer*
	nru, integer* ncc, real* d__, real* e, complex* vt, integer* ldvt,
	complex* u, integer* ldu, complex* c__, integer* ldc, real* rwork,
	integer* info);

/* Subroutine */
CLAPACK_API
int cgbbrd_(char* vect, integer* m, integer* n, integer* ncc,
	integer* kl, integer* ku, complex* ab, integer* ldab, real* d__,
	real* e, complex* q, integer* ldq, complex* pt, integer* ldpt,
	complex* c__, integer* ldc, complex* work, real* rwork, integer* info);

/* Subroutine */
CLAPACK_API
int cgbcon_(char* norm, integer* n, integer* kl, integer* ku,
	complex* ab, integer* ldab, integer* ipiv, real* anorm, real* rcond,
	complex* work, real* rwork, integer* info);

/* Subroutine */
CLAPACK_API
int cgbequ_(integer* m, integer* n, integer* kl, integer* ku,
	complex* ab, integer* ldab, real* r__, real* c__, real* rowcnd, real
	* colcnd, real* amax, integer* info);

/* Subroutine */
CLAPACK_API
int cgbequb_(integer* m, integer* n, integer* kl, integer*
	ku, complex* ab, integer* ldab, real* r__, real* c__, real* rowcnd,
	real* colcnd, real* amax, integer* info);

/* Subroutine */
CLAPACK_API
int cgbrfs_(char* trans, integer* n, integer* kl, integer*
	ku, integer* nrhs, complex* ab, integer* ldab, complex* afb, integer*
	ldafb, integer* ipiv, complex* b, integer* ldb, complex* x, integer*
	ldx, real* ferr, real* berr, complex* work, real* rwork, integer*
	info);


/* Subroutine */
CLAPACK_API
int cgbsv_(integer* n, integer* kl, integer* ku, integer*
	nrhs, complex* ab, integer* ldab, integer* ipiv, complex* b, integer*
	ldb, integer* info);

/* Subroutine */
CLAPACK_API
int cgbsvx_(char* fact, char* trans, integer* n, integer* kl,
	integer* ku, integer* nrhs, complex* ab, integer* ldab, complex* afb,
	integer* ldafb, integer* ipiv, char* equed, real* r__, real* c__,
	complex* b, integer* ldb, complex* x, integer* ldx, real* rcond, real
	* ferr, real* berr, complex* work, real* rwork, integer* info);

/* Subroutine */
CLAPACK_API
int cgbtf2_(integer* m, integer* n, integer* kl, integer* ku,
	complex* ab, integer* ldab, integer* ipiv, integer* info);

/* Subroutine */
CLAPACK_API
int cgbtrf_(integer* m, integer* n, integer* kl, integer* ku,
	complex* ab, integer* ldab, integer* ipiv, integer* info);

/* Subroutine */
CLAPACK_API
int cgbtrs_(char* trans, integer* n, integer* kl, integer*
	ku, integer* nrhs, complex* ab, integer* ldab, integer* ipiv, complex
	* b, integer* ldb, integer* info);

/* Subroutine */
CLAPACK_API
int cgebak_(char* job, char* side, integer* n, integer* ilo,
	integer* ihi, real* scale, integer* m, complex* v, integer* ldv,
	integer* info);

/* Subroutine */
CLAPACK_API
int cgebal_(char* job, integer* n, complex* a, integer* lda,
	integer* ilo, integer* ihi, real* scale, integer* info);

/* Subroutine */
CLAPACK_API
int cgebd2_(integer* m, integer* n, complex* a, integer* lda,
	real* d__, real* e, complex* tauq, complex* taup, complex* work,
	integer* info);

/* Subroutine */
CLAPACK_API
int cgebrd_(integer* m, integer* n, complex* a, integer* lda,
	real* d__, real* e, complex* tauq, complex* taup, complex* work,
	integer* lwork, integer* info);

/* Subroutine */
CLAPACK_API
int cgecon_(char* norm, integer* n, complex* a, integer* lda,
	real* anorm, real* rcond, complex* work, real* rwork, integer* info);

/* Subroutine */
CLAPACK_API
int cgeequ_(integer* m, integer* n, complex* a, integer* lda,
	real* r__, real* c__, real* rowcnd, real* colcnd, real* amax,
	integer* info);

/* Subroutine */
CLAPACK_API
int cgeequb_(integer* m, integer* n, complex* a, integer*
	lda, real* r__, real* c__, real* rowcnd, real* colcnd, real* amax,
	integer* info);

/* Subroutine */
CLAPACK_API
int cgees_(char* jobvs, char* sort, L_fp select, integer* n,
	complex* a, integer* lda, integer* sdim, complex* w, complex* vs,
	integer* ldvs, complex* work, integer* lwork, real* rwork, logical*
	bwork, integer* info);

/* Subroutine */
CLAPACK_API
int cgeesx_(char* jobvs, char* sort, L_fp select, char*
	sense, integer* n, complex* a, integer* lda, integer* sdim, complex*
	w, complex* vs, integer* ldvs, real* rconde, real* rcondv, complex*
	work, integer* lwork, real* rwork, logical* bwork, integer* info);

/* Subroutine */
CLAPACK_API
int cgeev_(char* jobvl, char* jobvr, integer* n, complex* a,
	integer* lda, complex* w, complex* vl, integer* ldvl, complex* vr,
	integer* ldvr, complex* work, integer* lwork, real* rwork, integer*
	info);

/* Subroutine */
CLAPACK_API
int cgeevx_(char* balanc, char* jobvl, char* jobvr, char*
	sense, integer* n, complex* a, integer* lda, complex* w, complex* vl,
	integer* ldvl, complex* vr, integer* ldvr, integer* ilo, integer* ihi,
	real* scale, real* abnrm, real* rconde, real* rcondv, complex* work,
	integer* lwork, real* rwork, integer* info);

/* Subroutine */
CLAPACK_API
int cgegs_(char* jobvsl, char* jobvsr, integer* n, complex*
	a, integer* lda, complex* b, integer* ldb, complex* alpha, complex*
	beta, complex* vsl, integer* ldvsl, complex* vsr, integer* ldvsr,
	complex* work, integer* lwork, real* rwork, integer* info);

/* Subroutine */
CLAPACK_API
int cgegv_(char* jobvl, char* jobvr, integer* n, complex* a,
	integer* lda, complex* b, integer* ldb, complex* alpha, complex* beta,
	complex* vl, integer* ldvl, complex* vr, integer* ldvr, complex*
	work, integer* lwork, real* rwork, integer* info);

/* Subroutine */
CLAPACK_API
int cgehd2_(integer* n, integer* ilo, integer* ihi, complex*
	a, integer* lda, complex* tau, complex* work, integer* info);

/* Subroutine */
CLAPACK_API
int cgehrd_(integer* n, integer* ilo, integer* ihi, complex*
	a, integer* lda, complex* tau, complex* work, integer* lwork, integer
	* info);

/* Subroutine */
CLAPACK_API
int cgelq2_(integer* m, integer* n, complex* a, integer* lda,
	complex* tau, complex* work, integer* info);

/* Subroutine */
CLAPACK_API
int cgelqf_(integer* m, integer* n, complex* a, integer* lda,
	complex* tau, complex* work, integer* lwork, integer* info);

/* Subroutine */
CLAPACK_API
int cgels_(char* trans, integer* m, integer* n, integer*
	nrhs, complex* a, integer* lda, complex* b, integer* ldb, complex*
	work, integer* lwork, integer* info);

/* Subroutine */
CLAPACK_API
int cgelsd_(integer* m, integer* n, integer* nrhs, complex*
	a, integer* lda, complex* b, integer* ldb, real* s, real* rcond,
	integer* rank, complex* work, integer* lwork, real* rwork, integer*
	iwork, integer* info);

/* Subroutine */
CLAPACK_API
int cgelss_(integer* m, integer* n, integer* nrhs, complex*
	a, integer* lda, complex* b, integer* ldb, real* s, real* rcond,
	integer* rank, complex* work, integer* lwork, real* rwork, integer*
	info);

/* Subroutine */
CLAPACK_API
int cgelsx_(integer* m, integer* n, integer* nrhs, complex*
	a, integer* lda, complex* b, integer* ldb, integer* jpvt, real* rcond,
	integer* rank, complex* work, real* rwork, integer* info);

/* Subroutine */
CLAPACK_API
int cgelsy_(integer* m, integer* n, integer* nrhs, complex*
	a, integer* lda, complex* b, integer* ldb, integer* jpvt, real* rcond,
	integer* rank, complex* work, integer* lwork, real* rwork, integer*
	info);

/* Subroutine */
CLAPACK_API
int cgeql2_(integer* m, integer* n, complex* a, integer* lda,
	complex* tau, complex* work, integer* info);

/* Subroutine */
CLAPACK_API
int cgeqlf_(integer* m, integer* n, complex* a, integer* lda,
	complex* tau, complex* work, integer* lwork, integer* info);

/* Subroutine */
CLAPACK_API
int cgeqp3_(integer* m, integer* n, complex* a, integer* lda,
	integer* jpvt, complex* tau, complex* work, integer* lwork, real*
	rwork, integer* info);

/* Subroutine */
CLAPACK_API
int cgeqpf_(integer* m, integer* n, complex* a, integer* lda,
	integer* jpvt, complex* tau, complex* work, real* rwork, integer*
	info);

/* Subroutine */
CLAPACK_API
int cgeqr2_(integer* m, integer* n, complex* a, integer* lda,
	complex* tau, complex* work, integer* info);

/* Subroutine */
CLAPACK_API
int cgeqrf_(integer* m, integer* n, complex* a, integer* lda,
	complex* tau, complex* work, integer* lwork, integer* info);

/* Subroutine */
CLAPACK_API
int cgerfs_(char* trans, integer* n, integer* nrhs, complex*
	a, integer* lda, complex* af, integer* ldaf, integer* ipiv, complex*
	b, integer* ldb, complex* x, integer* ldx, real* ferr, real* berr,
	complex* work, real* rwork, integer* info);

/* Subroutine */
CLAPACK_API
int cgerq2_(integer* m, integer* n, complex* a, integer* lda,
	complex* tau, complex* work, integer* info);

/* Subroutine */
CLAPACK_API
int cgerqf_(integer* m, integer* n, complex* a, integer* lda,
	complex* tau, complex* work, integer* lwork, integer* info);

/* Subroutine */
CLAPACK_API
int cgesc2_(integer* n, complex* a, integer* lda, complex*
	rhs, integer* ipiv, integer* jpiv, real* scale);

/* Subroutine */
CLAPACK_API
int cgesdd_(char* jobz, integer* m, integer* n, complex* a,
	integer* lda, real* s, complex* u, integer* ldu, complex* vt, integer
	* ldvt, complex* work, integer* lwork, real* rwork, integer* iwork,
	integer* info);

/* Subroutine */
CLAPACK_API
int cgesv_(integer* n, integer* nrhs, complex* a, integer*
	lda, integer* ipiv, complex* b, integer* ldb, integer* info);

/* Subroutine */
CLAPACK_API
int cgesvd_(char* jobu, char* jobvt, integer* m, integer* n,
	complex* a, integer* lda, real* s, complex* u, integer* ldu, complex*
	vt, integer* ldvt, complex* work, integer* lwork, real* rwork,
	integer* info);

/* Subroutine */
CLAPACK_API
int cgesvx_(char* fact, char* trans, integer* n, integer*
	nrhs, complex* a, integer* lda, complex* af, integer* ldaf, integer*
	ipiv, char* equed, real* r__, real* c__, complex* b, integer* ldb,
	complex* x, integer* ldx, real* rcond, real* ferr, real* berr,
	complex* work, real* rwork, integer* info);

/* Subroutine */
CLAPACK_API
int cgetc2_(integer* n, complex* a, integer* lda, integer*
	ipiv, integer* jpiv, integer* info);

/* Subroutine */
CLAPACK_API
int cgetf2_(integer* m, integer* n, complex* a, integer* lda,
	integer* ipiv, integer* info);

/* Subroutine */
CLAPACK_API
int cgetrf_(integer* m, integer* n, complex* a, integer* lda,
	integer* ipiv, integer* info);

/* Subroutine */
CLAPACK_API
int cgetri_(integer* n, complex* a, integer* lda, integer*
	ipiv, complex* work, integer* lwork, integer* info);

/* Subroutine */
CLAPACK_API
int cgetrs_(char* trans, integer* n, integer* nrhs, complex*
	a, integer* lda, integer* ipiv, complex* b, integer* ldb, integer*
	info);

/* Subroutine */
CLAPACK_API
int cggbak_(char* job, char* side, integer* n, integer* ilo,
	integer* ihi, real* lscale, real* rscale, integer* m, complex* v,
	integer* ldv, integer* info);

/* Subroutine */
CLAPACK_API
int cggbal_(char* job, integer* n, complex* a, integer* lda,
	complex* b, integer* ldb, integer* ilo, integer* ihi, real* lscale,
	real* rscale, real* work, integer* info);

/* Subroutine */
CLAPACK_API
int cgges_(char* jobvsl, char* jobvsr, char* sort, L_fp
	selctg, integer* n, complex* a, integer* lda, complex* b, integer*
	ldb, integer* sdim, complex* alpha, complex* beta, complex* vsl,
	integer* ldvsl, complex* vsr, integer* ldvsr, complex* work, integer*
	lwork, real* rwork, logical* bwork, integer* info);

/* Subroutine */
CLAPACK_API
int cggesx_(char* jobvsl, char* jobvsr, char* sort, L_fp
	selctg, char* sense, integer* n, complex* a, integer* lda, complex* b,
	integer* ldb, integer* sdim, complex* alpha, complex* beta, complex*
	vsl, integer* ldvsl, complex* vsr, integer* ldvsr, real* rconde, real
	* rcondv, complex* work, integer* lwork, real* rwork, integer* iwork,
	integer* liwork, logical* bwork, integer* info);

/* Subroutine */
CLAPACK_API
int cggev_(char* jobvl, char* jobvr, integer* n, complex* a,
	integer* lda, complex* b, integer* ldb, complex* alpha, complex* beta,
	complex* vl, integer* ldvl, complex* vr, integer* ldvr, complex*
	work, integer* lwork, real* rwork, integer* info);

/* Subroutine */
CLAPACK_API
int cggevx_(char* balanc, char* jobvl, char* jobvr, char*
	sense, integer* n, complex* a, integer* lda, complex* b, integer* ldb,
	complex* alpha, complex* beta, complex* vl, integer* ldvl, complex*
	vr, integer* ldvr, integer* ilo, integer* ihi, real* lscale, real*
	rscale, real* abnrm, real* bbnrm, real* rconde, real* rcondv, complex
	* work, integer* lwork, real* rwork, integer* iwork, logical* bwork,
	integer* info);

/* Subroutine */
CLAPACK_API
int cggglm_(integer* n, integer* m, integer* p, complex* a,
	integer* lda, complex* b, integer* ldb, complex* d__, complex* x,
	complex* y, complex* work, integer* lwork, integer* info);

/* Subroutine */
CLAPACK_API
int cgghrd_(char* compq, char* compz, integer* n, integer*
	ilo, integer* ihi, complex* a, integer* lda, complex* b, integer* ldb,
	complex* q, integer* ldq, complex* z__, integer* ldz, integer* info);

/* Subroutine */
CLAPACK_API
int cgglse_(integer* m, integer* n, integer* p, complex* a,
	integer* lda, complex* b, integer* ldb, complex* c__, complex* d__,
	complex* x, complex* work, integer* lwork, integer* info);

/* Subroutine */
CLAPACK_API
int cggqrf_(integer* n, integer* m, integer* p, complex* a,
	integer* lda, complex* taua, complex* b, integer* ldb, complex* taub,
	complex* work, integer* lwork, integer* info);

/* Subroutine */
CLAPACK_API
int cggrqf_(integer* m, integer* p, integer* n, complex* a,
	integer* lda, complex* taua, complex* b, integer* ldb, complex* taub,
	complex* work, integer* lwork, integer* info);

/* Subroutine */
CLAPACK_API
int cggsvd_(char* jobu, char* jobv, char* jobq, integer* m,
	integer* n, integer* p, integer* k, integer* l, complex* a, integer*
	lda, complex* b, integer* ldb, real* alpha, real* beta, complex* u,
	integer* ldu, complex* v, integer* ldv, complex* q, integer* ldq,
	complex* work, real* rwork, integer* iwork, integer* info);

/* Subroutine */
CLAPACK_API
int cggsvp_(char* jobu, char* jobv, char* jobq, integer* m,
	integer* p, integer* n, complex* a, integer* lda, complex* b, integer
	* ldb, real* tola, real* tolb, integer* k, integer* l, complex* u,
	integer* ldu, complex* v, integer* ldv, complex* q, integer* ldq,
	integer* iwork, real* rwork, complex* tau, complex* work, integer* info);

/* Subroutine */
CLAPACK_API
int cgtcon_(char* norm, integer* n, complex* dl, complex*
	d__, complex* du, complex* du2, integer* ipiv, real* anorm, real*
	rcond, complex* work, integer* info);

/* Subroutine */
CLAPACK_API
int cgtrfs_(char* trans, integer* n, integer* nrhs, complex*
	dl, complex* d__, complex* du, complex* dlf, complex* df, complex*
	duf, complex* du2, integer* ipiv, complex* b, integer* ldb, complex*
	x, integer* ldx, real* ferr, real* berr, complex* work, real* rwork,
	integer* info);

/* Subroutine */
CLAPACK_API
int cgtsv_(integer* n, integer* nrhs, complex* dl, complex*
	d__, complex* du, complex* b, integer* ldb, integer* info);

/* Subroutine */
CLAPACK_API
int cgtsvx_(char* fact, char* trans, integer* n, integer*
	nrhs, complex* dl, complex* d__, complex* du, complex* dlf, complex*
	df, complex* duf, complex* du2, integer* ipiv, complex* b, integer*
	ldb, complex* x, integer* ldx, real* rcond, real* ferr, real* berr,
	complex* work, real* rwork, integer* info);

/* Subroutine */
CLAPACK_API
int cgttrf_(integer* n, complex* dl, complex* d__, complex*
	du, complex* du2, integer* ipiv, integer* info);

/* Subroutine */
CLAPACK_API
int cgttrs_(char* trans, integer* n, integer* nrhs, complex*
	dl, complex* d__, complex* du, complex* du2, integer* ipiv, complex*
	b, integer* ldb, integer* info);

/* Subroutine */
CLAPACK_API
int cgtts2_(integer* itrans, integer* n, integer* nrhs,
	complex* dl, complex* d__, complex* du, complex* du2, integer* ipiv,
	complex* b, integer* ldb);

/* Subroutine */
CLAPACK_API
int chbev_(char* jobz, char* uplo, integer* n, integer* kd,
	complex* ab, integer* ldab, real* w, complex* z__, integer* ldz,
	complex* work, real* rwork, integer* info);

/* Subroutine */
CLAPACK_API
int chbevd_(char* jobz, char* uplo, integer* n, integer* kd,
	complex* ab, integer* ldab, real* w, complex* z__, integer* ldz,
	complex* work, integer* lwork, real* rwork, integer* lrwork, integer*
	iwork, integer* liwork, integer* info);

/* Subroutine */
CLAPACK_API
int chbevx_(char* jobz, char* range, char* uplo, integer* n,
	integer* kd, complex* ab, integer* ldab, complex* q, integer* ldq,
	real* vl, real* vu, integer* il, integer* iu, real* abstol, integer*
	m, real* w, complex* z__, integer* ldz, complex* work, real* rwork,
	integer* iwork, integer* ifail, integer* info);

/* Subroutine */
CLAPACK_API
int chbgst_(char* vect, char* uplo, integer* n, integer* ka,
	integer* kb, complex* ab, integer* ldab, complex* bb, integer* ldbb,
	complex* x, integer* ldx, complex* work, real* rwork, integer* info);







#pragma endregion

#pragma region CXLASRC -- Single precision complex LAPACK routines using extra precision (cxla)
//#       
//set(CXLASRC     cgesvxx.c cgerfsx.c cla_gerfsx_extended.c cla_geamv.c
//	cla_gercond_c.c cla_gercond_x.c cla_rpvgrw.c
//	csysvxx.c csyrfsx.c cla_syrfsx_extended.c cla_syamv.c
//	cla_syrcond_c.c cla_syrcond_x.c cla_syrpvgrw.c
//	cposvxx.c cporfsx.c cla_porfsx_extended.c
//	cla_porcond_c.c cla_porcond_x.c cla_porpvgrw.c
//	cgbsvxx.c cgbrfsx.c cla_gbrfsx_extended.c cla_gbamv.c
//	cla_gbrcond_c.c cla_gbrcond_x.c cla_gbrpvgrw.c
//	chesvxx.c cherfsx.c cla_herfsx_extended.c cla_heamv.c
//	cla_hercond_c.c cla_hercond_x.c cla_herpvgrw.c
//	cla_lin_berr.c clarscl2.c clascl2.c cla_wwaddw.c)

#pragma endregion



#pragma region DZLAUX -- Auxiliary routines called from both DOUBLE PRECISION COMPLEX * 16 (aladz)

/* Subroutine */
CLAPACK_API
int dbdsdc_(char* uplo, char* compq, integer* n, doublereal*
	d__, doublereal* e, doublereal* u, integer* ldu, doublereal* vt,
	integer* ldvt, doublereal* q, integer* iq, doublereal* work, integer*
	iwork, integer* info);

/* Subroutine */
CLAPACK_API
int dbdsqr_(char* uplo, integer* n, integer* ncvt, integer*
	nru, integer* ncc, doublereal* d__, doublereal* e, doublereal* vt,
	integer* ldvt, doublereal* u, integer* ldu, doublereal* c__, integer*
	ldc, doublereal* work, integer* info);

/* Subroutine */
CLAPACK_API
int ddisna_(char* job, integer* m, integer* n, doublereal*
	d__, doublereal* sep, integer* info);

CLAPACK_API
logical disnan_(doublereal* din);

/* Subroutine */
CLAPACK_API
int dlabad_(doublereal* small, doublereal* large);

/* Subroutine */
CLAPACK_API
int dlacpy_(char* uplo, integer* m, integer* n, doublereal*
	a, integer* lda, doublereal* b, integer* ldb);

/* Subroutine */
CLAPACK_API
int dladiv_(doublereal* a, doublereal* b, doublereal* c__,
	doublereal* d__, doublereal* p, doublereal* q);

/* Subroutine */
CLAPACK_API
int dlae2_(doublereal* a, doublereal* b, doublereal* c__,
	doublereal* rt1, doublereal* rt2);

/* Subroutine */
CLAPACK_API
int dlaebz_(integer* ijob, integer* nitmax, integer* n,
	integer* mmax, integer* minp, integer* nbmin, doublereal* abstol,
	doublereal* reltol, doublereal* pivmin, doublereal* d__, doublereal*
	e, doublereal* e2, integer* nval, doublereal* ab, doublereal* c__,
	integer* mout, integer* nab, doublereal* work, integer* iwork,
	integer* info);

/* Subroutine */
CLAPACK_API
int dlaed0_(integer* icompq, integer* qsiz, integer* n,
	doublereal* d__, doublereal* e, doublereal* q, integer* ldq,
	doublereal* qstore, integer* ldqs, doublereal* work, integer* iwork,
	integer* info);

/* Subroutine */
CLAPACK_API
int dlaed1_(integer* n, doublereal* d__, doublereal* q,
	integer* ldq, integer* indxq, doublereal* rho, integer* cutpnt,
	doublereal* work, integer* iwork, integer* info);

/* Subroutine */
CLAPACK_API
int dlaed2_(integer* k, integer* n, integer* n1, doublereal*
	d__, doublereal* q, integer* ldq, integer* indxq, doublereal* rho,
	doublereal* z__, doublereal* dlamda, doublereal* w, doublereal* q2,
	integer* indx, integer* indxc, integer* indxp, integer* coltyp,
	integer* info);

/* Subroutine */
CLAPACK_API
int dlaed3_(integer* k, integer* n, integer* n1, doublereal*
	d__, doublereal* q, integer* ldq, doublereal* rho, doublereal* dlamda,
	doublereal* q2, integer* indx, integer* ctot, doublereal* w,
	doublereal* s, integer* info);

/* Subroutine */
CLAPACK_API
int dlaed4_(integer* n, integer* i__, doublereal* d__,
	doublereal* z__, doublereal* delta, doublereal* rho, doublereal* dlam,
	integer* info);

/* Subroutine */
CLAPACK_API
int dlaed5_(integer* i__, doublereal* d__, doublereal* z__,
	doublereal* delta, doublereal* rho, doublereal* dlam);

/* Subroutine */
CLAPACK_API
int dlaed6_(integer* kniter, logical* orgati, doublereal*
	rho, doublereal* d__, doublereal* z__, doublereal* finit, doublereal*
	tau, integer* info);

/* Subroutine */
CLAPACK_API
int dlaed7_(integer* icompq, integer* n, integer* qsiz,
	integer* tlvls, integer* curlvl, integer* curpbm, doublereal* d__,
	doublereal* q, integer* ldq, integer* indxq, doublereal* rho, integer
	* cutpnt, doublereal* qstore, integer* qptr, integer* prmptr, integer*
	perm, integer* givptr, integer* givcol, doublereal* givnum,
	doublereal* work, integer* iwork, integer* info);

/* Subroutine */
CLAPACK_API
int dlaed8_(integer* icompq, integer* k, integer* n, integer
	* qsiz, doublereal* d__, doublereal* q, integer* ldq, integer* indxq,
	doublereal* rho, integer* cutpnt, doublereal* z__, doublereal* dlamda,
	doublereal* q2, integer* ldq2, doublereal* w, integer* perm, integer
	* givptr, integer* givcol, doublereal* givnum, integer* indxp, integer
	* indx, integer* info);

/* Subroutine */
CLAPACK_API
int dlaed9_(integer* k, integer* kstart, integer* kstop,
	integer* n, doublereal* d__, doublereal* q, integer* ldq, doublereal*
	rho, doublereal* dlamda, doublereal* w, doublereal* s, integer* lds,
	integer* info);

/* Subroutine */
CLAPACK_API
int dlaeda_(integer* n, integer* tlvls, integer* curlvl,
	integer* curpbm, integer* prmptr, integer* perm, integer* givptr,
	integer* givcol, doublereal* givnum, doublereal* q, integer* qptr,
	doublereal* z__, doublereal* ztemp, integer* info);

/* Subroutine */
CLAPACK_API
int dlaev2_(doublereal* a, doublereal* b, doublereal* c__,
	doublereal* rt1, doublereal* rt2, doublereal* cs1, doublereal* sn1);

/* Subroutine */
CLAPACK_API
int dlagtf_(integer* n, doublereal* a, doublereal* lambda,
	doublereal* b, doublereal* c__, doublereal* tol, doublereal* d__,
	integer* in, integer* info);

/* Subroutine */
CLAPACK_API
int dlagts_(integer* job, integer* n, doublereal* a,
	doublereal* b, doublereal* c__, doublereal* d__, integer* in,
	doublereal* y, doublereal* tol, integer* info);

CLAPACK_API
logical dlaisnan_(doublereal* din1, doublereal* din2);

CLAPACK_API
doublereal dlamch_(char* cmach);

/* Subroutine */
CLAPACK_API
int dlamrg_(integer* n1, integer* n2, doublereal* a, integer
	* dtrd1, integer* dtrd2, integer* index);

CLAPACK_API
integer dlaneg_(integer* n, doublereal* d__, doublereal* lld, doublereal*
	sigma, doublereal* pivmin, integer* r__);

CLAPACK_API
doublereal dlanst_(char* norm, integer* n, doublereal* d__, doublereal* e);

CLAPACK_API
doublereal dlapy2_(doublereal* x, doublereal* y);

CLAPACK_API
doublereal dlapy3_(doublereal* x, doublereal* y, doublereal* z__);

/* Subroutine */
CLAPACK_API
int dlarnv_(integer* idist, integer* iseed, integer* n,
	doublereal* x);

/* Subroutine */
CLAPACK_API
int dlarra_(integer* n, doublereal* d__, doublereal* e,
	doublereal* e2, doublereal* spltol, doublereal* tnrm, integer* nsplit,
	integer* isplit, integer* info);

/* Subroutine */
CLAPACK_API
int dlarrb_(integer* n, doublereal* d__, doublereal* lld,
	integer* ifirst, integer* ilast, doublereal* rtol1, doublereal* rtol2,
	integer* offset, doublereal* w, doublereal* wgap, doublereal* werr,
	doublereal* work, integer* iwork, doublereal* pivmin, doublereal*
	spdiam, integer* twist, integer* info);

/* Subroutine */
CLAPACK_API
int dlarrc_(char* jobt, integer* n, doublereal* vl,
	doublereal* vu, doublereal* d__, doublereal* e, doublereal* pivmin,
	integer* eigcnt, integer* lcnt, integer* rcnt, integer* info);

/* Subroutine */
CLAPACK_API
int dlarrd_(char* range, char* order, integer* n, doublereal
	* vl, doublereal* vu, integer* il, integer* iu, doublereal* gers,
	doublereal* reltol, doublereal* d__, doublereal* e, doublereal* e2,
	doublereal* pivmin, integer* nsplit, integer* isplit, integer* m,
	doublereal* w, doublereal* werr, doublereal* wl, doublereal* wu,
	integer* iblock, integer* indexw, doublereal* work, integer* iwork,
	integer* info);

/* Subroutine */
CLAPACK_API
int dlarre_(char* range, integer* n, doublereal* vl,
	doublereal* vu, integer* il, integer* iu, doublereal* d__, doublereal
	* e, doublereal* e2, doublereal* rtol1, doublereal* rtol2, doublereal*
	spltol, integer* nsplit, integer* isplit, integer* m, doublereal* w,
	doublereal* werr, doublereal* wgap, integer* iblock, integer* indexw,
	doublereal* gers, doublereal* pivmin, doublereal* work, integer*
	iwork, integer* info);

/* Subroutine */
CLAPACK_API
int dlarrf_(integer* n, doublereal* d__, doublereal* l,
	doublereal* ld, integer* clstrt, integer* clend, doublereal* w,
	doublereal* wgap, doublereal* werr, doublereal* spdiam, doublereal*
	clgapl, doublereal* clgapr, doublereal* pivmin, doublereal* sigma,
	doublereal* dplus, doublereal* lplus, doublereal* work, integer* info);

/* Subroutine */
CLAPACK_API
int dlarrj_(integer* n, doublereal* d__, doublereal* e2,
	integer* ifirst, integer* ilast, doublereal* rtol, integer* offset,
	doublereal* w, doublereal* werr, doublereal* work, integer* iwork,
	doublereal* pivmin, doublereal* spdiam, integer* info);

/* Subroutine */
CLAPACK_API
int dlarrk_(integer* n, integer* iw, doublereal* gl,
	doublereal* gu, doublereal* d__, doublereal* e2, doublereal* pivmin,
	doublereal* reltol, doublereal* w, doublereal* werr, integer* info);

/* Subroutine */
CLAPACK_API
int dlarrr_(integer* n, doublereal* d__, doublereal* e,
	integer* info);

/* Subroutine */
CLAPACK_API
int dlartg_(doublereal* f, doublereal* g, doublereal* cs,
	doublereal* sn, doublereal* r__);

/* Subroutine */
CLAPACK_API
int dlaruv_(integer* iseed, integer* n, doublereal* x);

/* Subroutine */
CLAPACK_API
int dlas2_(doublereal* f, doublereal* g, doublereal* h__,
	doublereal* ssmin, doublereal* ssmax);

/* Subroutine */
CLAPACK_API
int dlascl_(char* type__, integer* kl, integer* ku,
	doublereal* cfrom, doublereal* cto, integer* m, integer* n,
	doublereal* a, integer* lda, integer* info);

/* Subroutine */
CLAPACK_API
int dlasd0_(integer* n, integer* sqre, doublereal* d__,
	doublereal* e, doublereal* u, integer* ldu, doublereal* vt, integer*
	ldvt, integer* smlsiz, integer* iwork, doublereal* work, integer*
	info);

/* Subroutine */
CLAPACK_API
int dlasd1_(integer* nl, integer* nr, integer* sqre,
	doublereal* d__, doublereal* alpha, doublereal* beta, doublereal* u,
	integer* ldu, doublereal* vt, integer* ldvt, integer* idxq, integer*
	iwork, doublereal* work, integer* info);

/* Subroutine */
CLAPACK_API
int dlasd2_(integer* nl, integer* nr, integer* sqre, integer
	* k, doublereal* d__, doublereal* z__, doublereal* alpha, doublereal*
	beta, doublereal* u, integer* ldu, doublereal* vt, integer* ldvt,
	doublereal* dsigma, doublereal* u2, integer* ldu2, doublereal* vt2,
	integer* ldvt2, integer* idxp, integer* idx, integer* idxc, integer*
	idxq, integer* coltyp, integer* info);

/* Subroutine */
CLAPACK_API
int dlasd3_(integer* nl, integer* nr, integer* sqre, integer
	* k, doublereal* d__, doublereal* q, integer* ldq, doublereal* dsigma,
	doublereal* u, integer* ldu, doublereal* u2, integer* ldu2,
	doublereal* vt, integer* ldvt, doublereal* vt2, integer* ldvt2,
	integer* idxc, integer* ctot, doublereal* z__, integer* info);

/* Subroutine */
CLAPACK_API
int dlasd4_(integer* n, integer* i__, doublereal* d__,
	doublereal* z__, doublereal* delta, doublereal* rho, doublereal*
	sigma, doublereal* work, integer* info);

/* Subroutine */
CLAPACK_API
int dlasd5_(integer* i__, doublereal* d__, doublereal* z__,
	doublereal* delta, doublereal* rho, doublereal* dsigma, doublereal*
	work);

/* Subroutine */
CLAPACK_API
int dlasd6_(integer* icompq, integer* nl, integer* nr,
	integer* sqre, doublereal* d__, doublereal* vf, doublereal* vl,
	doublereal* alpha, doublereal* beta, integer* idxq, integer* perm,
	integer* givptr, integer* givcol, integer* ldgcol, doublereal* givnum,
	integer* ldgnum, doublereal* poles, doublereal* difl, doublereal*
	difr, doublereal* z__, integer* k, doublereal* c__, doublereal* s,
	doublereal* work, integer* iwork, integer* info);

/* Subroutine */
CLAPACK_API
int dlasd7_(integer* icompq, integer* nl, integer* nr,
	integer* sqre, integer* k, doublereal* d__, doublereal* z__,
	doublereal* zw, doublereal* vf, doublereal* vfw, doublereal* vl,
	doublereal* vlw, doublereal* alpha, doublereal* beta, doublereal*
	dsigma, integer* idx, integer* idxp, integer* idxq, integer* perm,
	integer* givptr, integer* givcol, integer* ldgcol, doublereal* givnum,
	integer* ldgnum, doublereal* c__, doublereal* s, integer* info);

/* Subroutine */
CLAPACK_API
int dlasd8_(integer* icompq, integer* k, doublereal* d__,
	doublereal* z__, doublereal* vf, doublereal* vl, doublereal* difl,
	doublereal* difr, integer* lddifr, doublereal* dsigma, doublereal*
	work, integer* info);

/* Subroutine */
CLAPACK_API
int dlasda_(integer* icompq, integer* smlsiz, integer* n,
	integer* sqre, doublereal* d__, doublereal* e, doublereal* u, integer
	* ldu, doublereal* vt, integer* k, doublereal* difl, doublereal* difr,
	doublereal* z__, doublereal* poles, integer* givptr, integer* givcol,
	integer* ldgcol, integer* perm, doublereal* givnum, doublereal* c__,
	doublereal* s, doublereal* work, integer* iwork, integer* info);

/* Subroutine */
CLAPACK_API
int dlasdq_(char* uplo, integer* sqre, integer* n, integer*
	ncvt, integer* nru, integer* ncc, doublereal* d__, doublereal* e,
	doublereal* vt, integer* ldvt, doublereal* u, integer* ldu,
	doublereal* c__, integer* ldc, doublereal* work, integer* info);

/* Subroutine */
CLAPACK_API
int dlasdt_(integer* n, integer* lvl, integer* nd, integer*
	inode, integer* ndiml, integer* ndimr, integer* msub);

/* Subroutine */
CLAPACK_API
int dlaset_(char* uplo, integer* m, integer* n, doublereal*
	alpha, doublereal* beta, doublereal* a, integer* lda);

/* Subroutine */
CLAPACK_API
int dlasq1_(integer* n, doublereal* d__, doublereal* e,
	doublereal* work, integer* info);

/* Subroutine */
CLAPACK_API
int dlasq2_(integer* n, doublereal* z__, integer* info);

/* Subroutine */
CLAPACK_API
int dlasq3_(integer* i0, integer* n0, doublereal* z__,
	integer* pp, doublereal* dmin__, doublereal* sigma, doublereal* desig,
	doublereal* qmax, integer* nfail, integer* iter, integer* ndiv,
	logical* ieee, integer* ttype, doublereal* dmin1, doublereal* dmin2,
	doublereal* dn, doublereal* dn1, doublereal* dn2, doublereal* g,
	doublereal* tau);

/* Subroutine */
CLAPACK_API
int dlasq4_(integer* i0, integer* n0, doublereal* z__,
	integer* pp, integer* n0in, doublereal* dmin__, doublereal* dmin1,
	doublereal* dmin2, doublereal* dn, doublereal* dn1, doublereal* dn2,
	doublereal* tau, integer* ttype, doublereal* g);

/* Subroutine */
CLAPACK_API
int dlasq5_(integer* i0, integer* n0, doublereal* z__,
	integer* pp, doublereal* tau, doublereal* dmin__, doublereal* dmin1,
	doublereal* dmin2, doublereal* dn, doublereal* dnm1, doublereal* dnm2,
	logical* ieee);

/* Subroutine */
CLAPACK_API
int dlasq6_(integer* i0, integer* n0, doublereal* z__,
	integer* pp, doublereal* dmin__, doublereal* dmin1, doublereal* dmin2,
	doublereal* dn, doublereal* dnm1, doublereal* dnm2);

/* Subroutine */
CLAPACK_API
int dlasr_(char* side, char* pivot, char* direct, integer* m,
	integer* n, doublereal* c__, doublereal* s, doublereal* a, integer*
	lda);

/* Subroutine */
CLAPACK_API
int dlasrt_(char* id, integer* n, doublereal* d__, integer*
	info);

/* Subroutine */
CLAPACK_API
int dlassq_(integer* n, doublereal* x, integer* incx,
	doublereal* scale, doublereal* sumsq);

/* Subroutine */
CLAPACK_API
int dlasv2_(doublereal* f, doublereal* g, doublereal* h__,
	doublereal* ssmin, doublereal* ssmax, doublereal* snr, doublereal*
	csr, doublereal* snl, doublereal* csl);

CLAPACK_API
integer dmaxloc_(doublereal* a, integer* dimm);

/* Subroutine */
CLAPACK_API
int dpttrf_(integer* n, doublereal* d__, doublereal* e,
	integer* info);

/* Subroutine */
CLAPACK_API
int dstebz_(char* range, char* order, integer* n, doublereal
	* vl, doublereal* vu, integer* il, integer* iu, doublereal* abstol,
	doublereal* d__, doublereal* e, integer* m, integer* nsplit,
	doublereal* w, integer* iblock, integer* isplit, doublereal* work,
	integer* iwork, integer* info);

/* Subroutine */
CLAPACK_API
int dstedc_(char* compz, integer* n, doublereal* d__,
	doublereal* e, doublereal* z__, integer* ldz, doublereal* work,
	integer* lwork, integer* iwork, integer* liwork, integer* info);

/* Subroutine */
CLAPACK_API
int dsteqr_(char* compz, integer* n, doublereal* d__,
	doublereal* e, doublereal* z__, integer* ldz, doublereal* work,
	integer* info);

/* Subroutine */
CLAPACK_API
int dsterf_(integer* n, doublereal* d__, doublereal* e,
	integer* info);

CLAPACK_API
doublereal dsecnd_();


#pragma endregion

#pragma region DLASRC -- Double precision real LAPACK routines (dla)

/* Subroutine */
CLAPACK_API
int dgbbrd_(char* vect, integer* m, integer* n, integer* ncc,
	integer* kl, integer* ku, doublereal* ab, integer* ldab, doublereal*
	d__, doublereal* e, doublereal* q, integer* ldq, doublereal* pt,
	integer* ldpt, doublereal* c__, integer* ldc, doublereal* work,
	integer* info);

/* Subroutine */
CLAPACK_API
int dgbcon_(char* norm, integer* n, integer* kl, integer* ku,
	doublereal* ab, integer* ldab, integer* ipiv, doublereal* anorm,
	doublereal* rcond, doublereal* work, integer* iwork, integer* info);

/* Subroutine */
CLAPACK_API
int dgbequ_(integer* m, integer* n, integer* kl, integer* ku,
	doublereal* ab, integer* ldab, doublereal* r__, doublereal* c__,
	doublereal* rowcnd, doublereal* colcnd, doublereal* amax, integer*
	info);

/* Subroutine */
CLAPACK_API
int dgbequb_(integer* m, integer* n, integer* kl, integer*
	ku, doublereal* ab, integer* ldab, doublereal* r__, doublereal* c__,
	doublereal* rowcnd, doublereal* colcnd, doublereal* amax, integer*
	info);

/* Subroutine */
CLAPACK_API
int dgbrfs_(char* trans, integer* n, integer* kl, integer*
	ku, integer* nrhs, doublereal* ab, integer* ldab, doublereal* afb,
	integer* ldafb, integer* ipiv, doublereal* b, integer* ldb,
	doublereal* x, integer* ldx, doublereal* ferr, doublereal* berr,
	doublereal* work, integer* iwork, integer* info);

/* Subroutine */
CLAPACK_API
int dgbsv_(integer* n, integer* kl, integer* ku, integer*
	nrhs, doublereal* ab, integer* ldab, integer* ipiv, doublereal* b,
	integer* ldb, integer* info);

/* Subroutine */
CLAPACK_API
int dgbsvx_(char* fact, char* trans, integer* n, integer* kl,
	integer* ku, integer* nrhs, doublereal* ab, integer* ldab,
	doublereal* afb, integer* ldafb, integer* ipiv, char* equed,
	doublereal* r__, doublereal* c__, doublereal* b, integer* ldb,
	doublereal* x, integer* ldx, doublereal* rcond, doublereal* ferr,
	doublereal* berr, doublereal* work, integer* iwork, integer* info);

/* Subroutine */
CLAPACK_API
int dgbtf2_(integer* m, integer* n, integer* kl, integer* ku,
	doublereal* ab, integer* ldab, integer* ipiv, integer* info);

/* Subroutine */
CLAPACK_API
int dgbtrf_(integer* m, integer* n, integer* kl, integer* ku,
	doublereal* ab, integer* ldab, integer* ipiv, integer* info);

/* Subroutine */
CLAPACK_API
int dgbtrs_(char* trans, integer* n, integer* kl, integer*
	ku, integer* nrhs, doublereal* ab, integer* ldab, integer* ipiv,
	doublereal* b, integer* ldb, integer* info);

/* Subroutine */
CLAPACK_API
int dgebak_(char* job, char* side, integer* n, integer* ilo,
	integer* ihi, doublereal* scale, integer* m, doublereal* v, integer*
	ldv, integer* info);

/* Subroutine */
CLAPACK_API
int dgebal_(char* job, integer* n, doublereal* a, integer*
	lda, integer* ilo, integer* ihi, doublereal* scale, integer* info);

/* Subroutine */
CLAPACK_API
int dgebd2_(integer* m, integer* n, doublereal* a, integer*
	lda, doublereal* d__, doublereal* e, doublereal* tauq, doublereal*
	taup, doublereal* work, integer* info);

/* Subroutine */
CLAPACK_API
int dgebrd_(integer* m, integer* n, doublereal* a, integer*
	lda, doublereal* d__, doublereal* e, doublereal* tauq, doublereal*
	taup, doublereal* work, integer* lwork, integer* info);

/* Subroutine */
CLAPACK_API
int dgecon_(char* norm, integer* n, doublereal* a, integer*
	lda, doublereal* anorm, doublereal* rcond, doublereal* work, integer*
	iwork, integer* info);

/* Subroutine */
CLAPACK_API
int dgeequ_(integer* m, integer* n, doublereal* a, integer*
	lda, doublereal* r__, doublereal* c__, doublereal* rowcnd, doublereal
	* colcnd, doublereal* amax, integer* info);

/* Subroutine */
CLAPACK_API
int dgeequb_(integer* m, integer* n, doublereal* a, integer*
	lda, doublereal* r__, doublereal* c__, doublereal* rowcnd, doublereal
	* colcnd, doublereal* amax, integer* info);

/* Subroutine */
CLAPACK_API
int dgees_(char* jobvs, char* sort, L_fp select, integer* n,
	doublereal* a, integer* lda, integer* sdim, doublereal* wr,
	doublereal* wi, doublereal* vs, integer* ldvs, doublereal* work,
	integer* lwork, logical* bwork, integer* info);

/* Subroutine */
CLAPACK_API
int dgeesx_(char* jobvs, char* sort, L_fp select, char*
	sense, integer* n, doublereal* a, integer* lda, integer* sdim,
	doublereal* wr, doublereal* wi, doublereal* vs, integer* ldvs,
	doublereal* rconde, doublereal* rcondv, doublereal* work, integer*
	lwork, integer* iwork, integer* liwork, logical* bwork, integer* info);

/* Subroutine */
CLAPACK_API
int dgeev_(char* jobvl, char* jobvr, integer* n, doublereal*
	a, integer* lda, doublereal* wr, doublereal* wi, doublereal* vl,
	integer* ldvl, doublereal* vr, integer* ldvr, doublereal* work,
	integer* lwork, integer* info);

/* Subroutine */
CLAPACK_API
int dgeevx_(char* balanc, char* jobvl, char* jobvr, char*
	sense, integer* n, doublereal* a, integer* lda, doublereal* wr,
	doublereal* wi, doublereal* vl, integer* ldvl, doublereal* vr,
	integer* ldvr, integer* ilo, integer* ihi, doublereal* scale,
	doublereal* abnrm, doublereal* rconde, doublereal* rcondv, doublereal
	* work, integer* lwork, integer* iwork, integer* info);

/* Subroutine */
CLAPACK_API
int dgegs_(char* jobvsl, char* jobvsr, integer* n,
	doublereal* a, integer* lda, doublereal* b, integer* ldb, doublereal*
	alphar, doublereal* alphai, doublereal* beta, doublereal* vsl,
	integer* ldvsl, doublereal* vsr, integer* ldvsr, doublereal* work,
	integer* lwork, integer* info);

/* Subroutine */
CLAPACK_API
int dgegv_(char* jobvl, char* jobvr, integer* n, doublereal*
	a, integer* lda, doublereal* b, integer* ldb, doublereal* alphar,
	doublereal* alphai, doublereal* beta, doublereal* vl, integer* ldvl,
	doublereal* vr, integer* ldvr, doublereal* work, integer* lwork,
	integer* info);

/* Subroutine */
CLAPACK_API
int dgehd2_(integer* n, integer* ilo, integer* ihi,
	doublereal* a, integer* lda, doublereal* tau, doublereal* work,
	integer* info);

/* Subroutine */
CLAPACK_API
int dgehrd_(integer* n, integer* ilo, integer* ihi,
	doublereal* a, integer* lda, doublereal* tau, doublereal* work,
	integer* lwork, integer* info);

/* Subroutine */
CLAPACK_API
int dgejsv_(char* joba, char* jobu, char* jobv, char* jobr,
	char* jobt, char* jobp, integer* m, integer* n, doublereal* a,
	integer* lda, doublereal* sva, doublereal* u, integer* ldu,
	doublereal* v, integer* ldv, doublereal* work, integer* lwork,
	integer* iwork, integer* info);

/* Subroutine */
CLAPACK_API
int dgelq2_(integer* m, integer* n, doublereal* a, integer*
	lda, doublereal* tau, doublereal* work, integer* info);

/* Subroutine */
CLAPACK_API
int dgelqf_(integer* m, integer* n, doublereal* a, integer*
	lda, doublereal* tau, doublereal* work, integer* lwork, integer* info);

/* Subroutine */
CLAPACK_API
int dgels_(char* trans, integer* m, integer* n, integer*
	nrhs, doublereal* a, integer* lda, doublereal* b, integer* ldb,
	doublereal* work, integer* lwork, integer* info);

/* Subroutine */
CLAPACK_API
int dgelsd_(integer* m, integer* n, integer* nrhs,
	doublereal* a, integer* lda, doublereal* b, integer* ldb, doublereal*
	s, doublereal* rcond, integer* rank, doublereal* work, integer* lwork,
	integer* iwork, integer* info);

/* Subroutine */
CLAPACK_API
int dgelss_(integer* m, integer* n, integer* nrhs,
	doublereal* a, integer* lda, doublereal* b, integer* ldb, doublereal*
	s, doublereal* rcond, integer* rank, doublereal* work, integer* lwork,
	integer* info);

/* Subroutine */
CLAPACK_API
int dgelsx_(integer* m, integer* n, integer* nrhs,
	doublereal* a, integer* lda, doublereal* b, integer* ldb, integer*
	jpvt, doublereal* rcond, integer* rank, doublereal* work, integer*
	info);

/* Subroutine */
CLAPACK_API
int dgelsy_(integer* m, integer* n, integer* nrhs,
	doublereal* a, integer* lda, doublereal* b, integer* ldb, integer*
	jpvt, doublereal* rcond, integer* rank, doublereal* work, integer*
	lwork, integer* info);

/* Subroutine */
CLAPACK_API
int dgeql2_(integer* m, integer* n, doublereal* a, integer*
	lda, doublereal* tau, doublereal* work, integer* info);

/* Subroutine */
CLAPACK_API
int dgeqlf_(integer* m, integer* n, doublereal* a, integer*
	lda, doublereal* tau, doublereal* work, integer* lwork, integer* info);

/* Subroutine */
CLAPACK_API
int dgeqp3_(integer* m, integer* n, doublereal* a, integer*
	lda, integer* jpvt, doublereal* tau, doublereal* work, integer* lwork,
	integer* info);

/* Subroutine */
CLAPACK_API
int dgeqpf_(integer* m, integer* n, doublereal* a, integer*
	lda, integer* jpvt, doublereal* tau, doublereal* work, integer* info);

/* Subroutine */
CLAPACK_API
int dgeqr2_(integer* m, integer* n, doublereal* a, integer*
	lda, doublereal* tau, doublereal* work, integer* info);

/* Subroutine */
CLAPACK_API
int dgeqrf_(integer* m, integer* n, doublereal* a, integer*
	lda, doublereal* tau, doublereal* work, integer* lwork, integer* info);

/* Subroutine */
CLAPACK_API
int dgerfs_(char* trans, integer* n, integer* nrhs,
	doublereal* a, integer* lda, doublereal* af, integer* ldaf, integer*
	ipiv, doublereal* b, integer* ldb, doublereal* x, integer* ldx,
	doublereal* ferr, doublereal* berr, doublereal* work, integer* iwork,
	integer* info);

/* Subroutine */
CLAPACK_API
int dgerq2_(integer* m, integer* n, doublereal* a, integer*
	lda, doublereal* tau, doublereal* work, integer* info);

/* Subroutine */
CLAPACK_API
int dgerqf_(integer* m, integer* n, doublereal* a, integer*
	lda, doublereal* tau, doublereal* work, integer* lwork, integer* info);

/* Subroutine */
CLAPACK_API
int dgesc2_(integer* n, doublereal* a, integer* lda,
	doublereal* rhs, integer* ipiv, integer* jpiv, doublereal* scale);

/* Subroutine */
CLAPACK_API
int dgesdd_(char* jobz, integer* m, integer* n, doublereal*
	a, integer* lda, doublereal* s, doublereal* u, integer* ldu,
	doublereal* vt, integer* ldvt, doublereal* work, integer* lwork,
	integer* iwork, integer* info);

/* Subroutine */
CLAPACK_API
int dgesv_(integer* n, integer* nrhs, doublereal* a, integer
	* lda, integer* ipiv, doublereal* b, integer* ldb, integer* info);

/* Subroutine */
CLAPACK_API
int dgesvd_(char* jobu, char* jobvt, integer* m, integer* n,
	doublereal* a, integer* lda, doublereal* s, doublereal* u, integer*
	ldu, doublereal* vt, integer* ldvt, doublereal* work, integer* lwork,
	integer* info);

/* Subroutine */
CLAPACK_API
int dgesvj_(char* joba, char* jobu, char* jobv, integer* m,
	integer* n, doublereal* a, integer* lda, doublereal* sva, integer* mv,
	doublereal* v, integer* ldv, doublereal* work, integer* lwork,
	integer* info);

/* Subroutine */
CLAPACK_API
int dgesvx_(char* fact, char* trans, integer* n, integer*
	nrhs, doublereal* a, integer* lda, doublereal* af, integer* ldaf,
	integer* ipiv, char* equed, doublereal* r__, doublereal* c__,
	doublereal* b, integer* ldb, doublereal* x, integer* ldx, doublereal*
	rcond, doublereal* ferr, doublereal* berr, doublereal* work, integer*
	iwork, integer* info);

/* Subroutine */
CLAPACK_API
int dgetc2_(integer* n, doublereal* a, integer* lda, integer
	* ipiv, integer* jpiv, integer* info);

/* Subroutine */
CLAPACK_API
int dgetf2_(integer* m, integer* n, doublereal* a, integer*
	lda, integer* ipiv, integer* info);

/* Subroutine */
CLAPACK_API
int dgetrf_(integer* m, integer* n, doublereal* a, integer*
	lda, integer* ipiv, integer* info);

/* Subroutine */
CLAPACK_API
int dgetri_(integer* n, doublereal* a, integer* lda, integer
	* ipiv, doublereal* work, integer* lwork, integer* info);

/* Subroutine */
CLAPACK_API
int dgetrs_(char* trans, integer* n, integer* nrhs,
	doublereal* a, integer* lda, integer* ipiv, doublereal* b, integer*
	ldb, integer* info);

/* Subroutine */
CLAPACK_API
int dggbak_(char* job, char* side, integer* n, integer* ilo,
	integer* ihi, doublereal* lscale, doublereal* rscale, integer* m,
	doublereal* v, integer* ldv, integer* info);

/* Subroutine */
CLAPACK_API
int dggbal_(char* job, integer* n, doublereal* a, integer*
	lda, doublereal* b, integer* ldb, integer* ilo, integer* ihi,
	doublereal* lscale, doublereal* rscale, doublereal* work, integer*
	info);

/* Subroutine */
CLAPACK_API
int dgges_(char* jobvsl, char* jobvsr, char* sort, L_fp
	selctg, integer* n, doublereal* a, integer* lda, doublereal* b,
	integer* ldb, integer* sdim, doublereal* alphar, doublereal* alphai,
	doublereal* beta, doublereal* vsl, integer* ldvsl, doublereal* vsr,
	integer* ldvsr, doublereal* work, integer* lwork, logical* bwork,
	integer* info);

/* Subroutine */
CLAPACK_API
int dggesx_(char* jobvsl, char* jobvsr, char* sort, L_fp
	selctg, char* sense, integer* n, doublereal* a, integer* lda,
	doublereal* b, integer* ldb, integer* sdim, doublereal* alphar,
	doublereal* alphai, doublereal* beta, doublereal* vsl, integer* ldvsl,
	doublereal* vsr, integer* ldvsr, doublereal* rconde, doublereal*
	rcondv, doublereal* work, integer* lwork, integer* iwork, integer*
	liwork, logical* bwork, integer* info);

/* Subroutine */
CLAPACK_API
int dggev_(char* jobvl, char* jobvr, integer* n, doublereal*
	a, integer* lda, doublereal* b, integer* ldb, doublereal* alphar,
	doublereal* alphai, doublereal* beta, doublereal* vl, integer* ldvl,
	doublereal* vr, integer* ldvr, doublereal* work, integer* lwork,
	integer* info);

/* Subroutine */
CLAPACK_API
int dggevx_(char* balanc, char* jobvl, char* jobvr, char*
	sense, integer* n, doublereal* a, integer* lda, doublereal* b,
	integer* ldb, doublereal* alphar, doublereal* alphai, doublereal*
	beta, doublereal* vl, integer* ldvl, doublereal* vr, integer* ldvr,
	integer* ilo, integer* ihi, doublereal* lscale, doublereal* rscale,
	doublereal* abnrm, doublereal* bbnrm, doublereal* rconde, doublereal*
	rcondv, doublereal* work, integer* lwork, integer* iwork, logical*
	bwork, integer* info);

/* Subroutine */
CLAPACK_API
int dggglm_(integer* n, integer* m, integer* p, doublereal*
	a, integer* lda, doublereal* b, integer* ldb, doublereal* d__,
	doublereal* x, doublereal* y, doublereal* work, integer* lwork,
	integer* info);

/* Subroutine */
CLAPACK_API
int dgghrd_(char* compq, char* compz, integer* n, integer*
	ilo, integer* ihi, doublereal* a, integer* lda, doublereal* b,
	integer* ldb, doublereal* q, integer* ldq, doublereal* z__, integer*
	ldz, integer* info);

/* Subroutine */
CLAPACK_API
int dgglse_(integer* m, integer* n, integer* p, doublereal*
	a, integer* lda, doublereal* b, integer* ldb, doublereal* c__,
	doublereal* d__, doublereal* x, doublereal* work, integer* lwork,
	integer* info);

/* Subroutine */
CLAPACK_API
int dggqrf_(integer* n, integer* m, integer* p, doublereal*
	a, integer* lda, doublereal* taua, doublereal* b, integer* ldb,
	doublereal* taub, doublereal* work, integer* lwork, integer* info);

/* Subroutine */
CLAPACK_API
int dggrqf_(integer* m, integer* p, integer* n, doublereal*
	a, integer* lda, doublereal* taua, doublereal* b, integer* ldb,
	doublereal* taub, doublereal* work, integer* lwork, integer* info);

/* Subroutine */
CLAPACK_API
int dggsvd_(char* jobu, char* jobv, char* jobq, integer* m,
	integer* n, integer* p, integer* k, integer* l, doublereal* a,
	integer* lda, doublereal* b, integer* ldb, doublereal* alpha,
	doublereal* beta, doublereal* u, integer* ldu, doublereal* v, integer
	* ldv, doublereal* q, integer* ldq, doublereal* work, integer* iwork,
	integer* info);

/* Subroutine */
CLAPACK_API
int dggsvp_(char* jobu, char* jobv, char* jobq, integer* m,
	integer* p, integer* n, doublereal* a, integer* lda, doublereal* b,
	integer* ldb, doublereal* tola, doublereal* tolb, integer* k, integer
	* l, doublereal* u, integer* ldu, doublereal* v, integer* ldv,
	doublereal* q, integer* ldq, integer* iwork, doublereal* tau,
	doublereal* work, integer* info);

/* Subroutine */
CLAPACK_API
int dgsvj0_(char* jobv, integer* m, integer* n, doublereal*
	a, integer* lda, doublereal* d__, doublereal* sva, integer* mv,
	doublereal* v, integer* ldv, doublereal* eps, doublereal* sfmin,
	doublereal* tol, integer* nsweep, doublereal* work, integer* lwork,
	integer* info);

/* Subroutine */
CLAPACK_API
int dgsvj1_(char* jobv, integer* m, integer* n, integer* n1,
	doublereal* a, integer* lda, doublereal* d__, doublereal* sva,
	integer* mv, doublereal* v, integer* ldv, doublereal* eps, doublereal
	* sfmin, doublereal* tol, integer* nsweep, doublereal* work, integer*
	lwork, integer* info);

/* Subroutine */
CLAPACK_API
int dgtcon_(char* norm, integer* n, doublereal* dl,
	doublereal* d__, doublereal* du, doublereal* du2, integer* ipiv,
	doublereal* anorm, doublereal* rcond, doublereal* work, integer*
	iwork, integer* info);

/* Subroutine */
CLAPACK_API
int dgtrfs_(char* trans, integer* n, integer* nrhs,
	doublereal* dl, doublereal* d__, doublereal* du, doublereal* dlf,
	doublereal* df, doublereal* duf, doublereal* du2, integer* ipiv,
	doublereal* b, integer* ldb, doublereal* x, integer* ldx, doublereal*
	ferr, doublereal* berr, doublereal* work, integer* iwork, integer*
	info);

/* Subroutine */
CLAPACK_API
int dgtsv_(integer* n, integer* nrhs, doublereal* dl,
	doublereal* d__, doublereal* du, doublereal* b, integer* ldb, integer
	* info);

/* Subroutine */
CLAPACK_API
int dgtsvx_(char* fact, char* trans, integer* n, integer*
	nrhs, doublereal* dl, doublereal* d__, doublereal* du, doublereal*
	dlf, doublereal* df, doublereal* duf, doublereal* du2, integer* ipiv,
	doublereal* b, integer* ldb, doublereal* x, integer* ldx, doublereal*
	rcond, doublereal* ferr, doublereal* berr, doublereal* work, integer*
	iwork, integer* info);

/* Subroutine */
CLAPACK_API
int dgttrf_(integer* n, doublereal* dl, doublereal* d__,
	doublereal* du, doublereal* du2, integer* ipiv, integer* info);

/* Subroutine */
CLAPACK_API
int dgttrs_(char* trans, integer* n, integer* nrhs,
	doublereal* dl, doublereal* d__, doublereal* du, doublereal* du2,
	integer* ipiv, doublereal* b, integer* ldb, integer* info);

/* Subroutine */
CLAPACK_API
int dgtts2_(integer* itrans, integer* n, integer* nrhs,
	doublereal* dl, doublereal* d__, doublereal* du, doublereal* du2,
	integer* ipiv, doublereal* b, integer* ldb);

/* Subroutine */
CLAPACK_API
int dhgeqz_(char* job, char* compq, char* compz, integer* n,
	integer* ilo, integer* ihi, doublereal* h__, integer* ldh, doublereal
	* t, integer* ldt, doublereal* alphar, doublereal* alphai, doublereal*
	beta, doublereal* q, integer* ldq, doublereal* z__, integer* ldz,
	doublereal* work, integer* lwork, integer* info);

/* Subroutine */
CLAPACK_API
int dhsein_(char* side, char* eigsrc, char* initv, logical*
	select, integer* n, doublereal* h__, integer* ldh, doublereal* wr,
	doublereal* wi, doublereal* vl, integer* ldvl, doublereal* vr,
	integer* ldvr, integer* mm, integer* m, doublereal* work, integer*
	ifaill, integer* ifailr, integer* info);

/* Subroutine */
CLAPACK_API
int dhseqr_(char* job, char* compz, integer* n, integer* ilo,
	integer* ihi, doublereal* h__, integer* ldh, doublereal* wr,
	doublereal* wi, doublereal* z__, integer* ldz, doublereal* work,
	integer* lwork, integer* info);

/* Subroutine */
CLAPACK_API
int dlabrd_(integer* m, integer* n, integer* nb, doublereal*
	a, integer* lda, doublereal* d__, doublereal* e, doublereal* tauq,
	doublereal* taup, doublereal* x, integer* ldx, doublereal* y, integer
	* ldy);

/* Subroutine */
CLAPACK_API
int dlacn2_(integer* n, doublereal* v, doublereal* x,
	integer* isgn, doublereal* est, integer* kase, integer* isave);

/* Subroutine */
CLAPACK_API
int dlacon_(integer* n, doublereal* v, doublereal* x,
	integer* isgn, doublereal* est, integer* kase);

/* Subroutine */
CLAPACK_API
int dlaein_(logical* rightv, logical* noinit, integer* n,
	doublereal* h__, integer* ldh, doublereal* wr, doublereal* wi,
	doublereal* vr, doublereal* vi, doublereal* b, integer* ldb,
	doublereal* work, doublereal* eps3, doublereal* smlnum, doublereal*
	bignum, integer* info);

/* Subroutine */
CLAPACK_API
int dlaexc_(logical* wantq, integer* n, doublereal* t,
	integer* ldt, doublereal* q, integer* ldq, integer* j1, integer* n1,
	integer* n2, doublereal* work, integer* info);

/* Subroutine */
CLAPACK_API
int dlag2_(doublereal* a, integer* lda, doublereal* b,
	integer* ldb, doublereal* safmin, doublereal* scale1, doublereal*
	scale2, doublereal* wr1, doublereal* wr2, doublereal* wi);

/* Subroutine */
CLAPACK_API
int dlag2s_(integer* m, integer* n, doublereal* a, integer*
	lda, real* sa, integer* ldsa, integer* info);

/* Subroutine */
CLAPACK_API
int dlags2_(logical* upper, doublereal* a1, doublereal* a2,
	doublereal* a3, doublereal* b1, doublereal* b2, doublereal* b3,
	doublereal* csu, doublereal* snu, doublereal* csv, doublereal* snv,
	doublereal* csq, doublereal* snq);

/* Subroutine */
CLAPACK_API
int dlagtm_(char* trans, integer* n, integer* nrhs,
	doublereal* alpha, doublereal* dl, doublereal* d__, doublereal* du,
	doublereal* x, integer* ldx, doublereal* beta, doublereal* b, integer
	* ldb);

/* Subroutine */
CLAPACK_API
int dlagv2_(doublereal* a, integer* lda, doublereal* b,
	integer* ldb, doublereal* alphar, doublereal* alphai, doublereal*
	beta, doublereal* csl, doublereal* snl, doublereal* csr, doublereal*
	snr);

/* Subroutine */
CLAPACK_API
int dlahqr_(logical* wantt, logical* wantz, integer* n,
	integer* ilo, integer* ihi, doublereal* h__, integer* ldh, doublereal
	* wr, doublereal* wi, integer* iloz, integer* ihiz, doublereal* z__,
	integer* ldz, integer* info);

/* Subroutine */
CLAPACK_API
int dlahr2_(integer* n, integer* k, integer* nb, doublereal*
	a, integer* lda, doublereal* tau, doublereal* t, integer* ldt,
	doublereal* y, integer* ldy);

/* Subroutine */
CLAPACK_API
int dlahrd_(integer* n, integer* k, integer* nb, doublereal*
	a, integer* lda, doublereal* tau, doublereal* t, integer* ldt,
	doublereal* y, integer* ldy);

/* Subroutine */
CLAPACK_API
int dlaic1_(integer* job, integer* j, doublereal* x,
	doublereal* sest, doublereal* w, doublereal* gamma, doublereal*
	sestpr, doublereal* s, doublereal* c__);

/* Subroutine */
CLAPACK_API
int dlaln2_(logical* ltrans, integer* na, integer* nw,
	doublereal* smin, doublereal* ca, doublereal* a, integer* lda,
	doublereal* d1, doublereal* d2, doublereal* b, integer* ldb,
	doublereal* wr, doublereal* wi, doublereal* x, integer* ldx,
	doublereal* scale, doublereal* xnorm, integer* info);

/* Subroutine */
CLAPACK_API
int dlals0_(integer* icompq, integer* nl, integer* nr,
	integer* sqre, integer* nrhs, doublereal* b, integer* ldb, doublereal
	* bx, integer* ldbx, integer* perm, integer* givptr, integer* givcol,
	integer* ldgcol, doublereal* givnum, integer* ldgnum, doublereal*
	poles, doublereal* difl, doublereal* difr, doublereal* z__, integer*
	k, doublereal* c__, doublereal* s, doublereal* work, integer* info);

/* Subroutine */
CLAPACK_API
int dlalsa_(integer* icompq, integer* smlsiz, integer* n,
	integer* nrhs, doublereal* b, integer* ldb, doublereal* bx, integer*
	ldbx, doublereal* u, integer* ldu, doublereal* vt, integer* k,
	doublereal* difl, doublereal* difr, doublereal* z__, doublereal*
	poles, integer* givptr, integer* givcol, integer* ldgcol, integer*
	perm, doublereal* givnum, doublereal* c__, doublereal* s, doublereal*
	work, integer* iwork, integer* info);

/* Subroutine */
CLAPACK_API
int dlalsd_(char* uplo, integer* smlsiz, integer* n, integer
	* nrhs, doublereal* d__, doublereal* e, doublereal* b, integer* ldb,
	doublereal* rcond, integer* rank, doublereal* work, integer* iwork,
	integer* info);

CLAPACK_API
doublereal dlangb_(char* norm, integer* n, integer* kl, integer* ku,
	doublereal* ab, integer* ldab, doublereal* work);

CLAPACK_API
doublereal dlange_(char* norm, integer* m, integer* n, doublereal* a, integer
	* lda, doublereal* work);

CLAPACK_API
doublereal dlangt_(char* norm, integer* n, doublereal* dl, doublereal* d__,
	doublereal* du);

CLAPACK_API
doublereal dlanhs_(char* norm, integer* n, doublereal* a, integer* lda,
	doublereal* work);

CLAPACK_API
doublereal dlansb_(char* norm, char* uplo, integer* n, integer* k, doublereal
	* ab, integer* ldab, doublereal* work);

CLAPACK_API
doublereal dlansf_(char* norm, char* transr, char* uplo, integer* n,
	doublereal* a, doublereal* work);

CLAPACK_API
doublereal dlansp_(char* norm, char* uplo, integer* n, doublereal* ap,
	doublereal* work);

CLAPACK_API
doublereal dlansy_(char* norm, char* uplo, integer* n, doublereal* a, integer
	* lda, doublereal* work);

CLAPACK_API
doublereal dlantb_(char* norm, char* uplo, char* diag, integer* n, integer* k,
	doublereal* ab, integer* ldab, doublereal* work);

CLAPACK_API
doublereal dlantp_(char* norm, char* uplo, char* diag, integer* n, doublereal
	* ap, doublereal* work);

CLAPACK_API
doublereal dlantr_(char* norm, char* uplo, char* diag, integer* m, integer* n,
	doublereal* a, integer* lda, doublereal* work);

/* Subroutine */
CLAPACK_API
int dlanv2_(doublereal* a, doublereal* b, doublereal* c__,
	doublereal* d__, doublereal* rt1r, doublereal* rt1i, doublereal* rt2r,
	doublereal* rt2i, doublereal* cs, doublereal* sn);

/* Subroutine */
CLAPACK_API
int dlapll_(integer* n, doublereal* x, integer* incx,
	doublereal* y, integer* incy, doublereal* ssmin);

/* Subroutine */
CLAPACK_API
int dlapmt_(logical* forwrd, integer* m, integer* n,
	doublereal* x, integer* ldx, integer* k);

/* Subroutine */
CLAPACK_API
int dlaqgb_(integer* m, integer* n, integer* kl, integer* ku,
	doublereal* ab, integer* ldab, doublereal* r__, doublereal* c__,
	doublereal* rowcnd, doublereal* colcnd, doublereal* amax, char* equed);

/* Subroutine */
CLAPACK_API
int dlaqge_(integer* m, integer* n, doublereal* a, integer*
	lda, doublereal* r__, doublereal* c__, doublereal* rowcnd, doublereal
	* colcnd, doublereal* amax, char* equed);

/* Subroutine */
CLAPACK_API
int dlaqp2_(integer* m, integer* n, integer* offset,
	doublereal* a, integer* lda, integer* jpvt, doublereal* tau,
	doublereal* vn1, doublereal* vn2, doublereal* work);

/* Subroutine */
CLAPACK_API
int dlaqps_(integer* m, integer* n, integer* offset, integer
	* nb, integer* kb, doublereal* a, integer* lda, integer* jpvt,
	doublereal* tau, doublereal* vn1, doublereal* vn2, doublereal* auxv,
	doublereal* f, integer* ldf);

/* Subroutine */
CLAPACK_API
int dlaqr0_(logical* wantt, logical* wantz, integer* n,
	integer* ilo, integer* ihi, doublereal* h__, integer* ldh, doublereal
	* wr, doublereal* wi, integer* iloz, integer* ihiz, doublereal* z__,
	integer* ldz, doublereal* work, integer* lwork, integer* info);

/* Subroutine */
CLAPACK_API
int dlaqr1_(integer* n, doublereal* h__, integer* ldh,
	doublereal* sr1, doublereal* si1, doublereal* sr2, doublereal* si2,
	doublereal* v);

/* Subroutine */
CLAPACK_API
int dlaqr2_(logical* wantt, logical* wantz, integer* n,
	integer* ktop, integer* kbot, integer* nw, doublereal* h__, integer*
	ldh, integer* iloz, integer* ihiz, doublereal* z__, integer* ldz,
	integer* ns, integer* nd, doublereal* sr, doublereal* si, doublereal*
	v, integer* ldv, integer* nh, doublereal* t, integer* ldt, integer*
	nv, doublereal* wv, integer* ldwv, doublereal* work, integer* lwork);

/* Subroutine */
CLAPACK_API
int dlaqr3_(logical* wantt, logical* wantz, integer* n,
	integer* ktop, integer* kbot, integer* nw, doublereal* h__, integer*
	ldh, integer* iloz, integer* ihiz, doublereal* z__, integer* ldz,
	integer* ns, integer* nd, doublereal* sr, doublereal* si, doublereal*
	v, integer* ldv, integer* nh, doublereal* t, integer* ldt, integer*
	nv, doublereal* wv, integer* ldwv, doublereal* work, integer* lwork);

/* Subroutine */CLAPACK_API
int dlaqr4_(logical* wantt, logical* wantz, integer* n,
	integer* ilo, integer* ihi, doublereal* h__, integer* ldh, doublereal
	* wr, doublereal* wi, integer* iloz, integer* ihiz, doublereal* z__,
	integer* ldz, doublereal* work, integer* lwork, integer* info);

/* Subroutine */
CLAPACK_API
int dlaqr5_(logical* wantt, logical* wantz, integer* kacc22,
	integer* n, integer* ktop, integer* kbot, integer* nshfts, doublereal
	* sr, doublereal* si, doublereal* h__, integer* ldh, integer* iloz,
	integer* ihiz, doublereal* z__, integer* ldz, doublereal* v, integer*
	ldv, doublereal* u, integer* ldu, integer* nv, doublereal* wv,
	integer* ldwv, integer* nh, doublereal* wh, integer* ldwh);

/* Subroutine */
CLAPACK_API
int dlaqsb_(char* uplo, integer* n, integer* kd, doublereal*
	ab, integer* ldab, doublereal* s, doublereal* scond, doublereal* amax,
	char* equed);

/* Subroutine */
CLAPACK_API
int dlaqsp_(char* uplo, integer* n, doublereal* ap,
	doublereal* s, doublereal* scond, doublereal* amax, char* equed);

/* Subroutine */
CLAPACK_API
int dlaqsy_(char* uplo, integer* n, doublereal* a, integer*
	lda, doublereal* s, doublereal* scond, doublereal* amax, char* equed);

/* Subroutine */
CLAPACK_API
int dlaqtr_(logical* ltran, logical* lreal, integer* n,
	doublereal* t, integer* ldt, doublereal* b, doublereal* w, doublereal
	* scale, doublereal* x, doublereal* work, integer* info);

/* Subroutine */
CLAPACK_API
int dlar1v_(integer* n, integer* b1, integer* bn, doublereal
	* lambda, doublereal* d__, doublereal* l, doublereal* ld, doublereal*
	lld, doublereal* pivmin, doublereal* gaptol, doublereal* z__, logical
	* wantnc, integer* negcnt, doublereal* ztz, doublereal* mingma,
	integer* r__, integer* isuppz, doublereal* nrminv, doublereal* resid,
	doublereal* rqcorr, doublereal* work);

/* Subroutine */
CLAPACK_API
int dlar2v_(integer* n, doublereal* x, doublereal* y,
	doublereal* z__, integer* incx, doublereal* c__, doublereal* s,
	integer* incc);

/* Subroutine */
CLAPACK_API
int dlarf_(char* side, integer* m, integer* n, doublereal* v,
	integer* incv, doublereal* tau, doublereal* c__, integer* ldc,
	doublereal* work);

/* Subroutine */
CLAPACK_API
int dlarfb_(char* side, char* trans, char* direct, char*
	storev, integer* m, integer* n, integer* k, doublereal* v, integer*
	ldv, doublereal* t, integer* ldt, doublereal* c__, integer* ldc,
	doublereal* work, integer* ldwork);

/* Subroutine */
CLAPACK_API
int dlarfg_(integer* n, doublereal* alpha, doublereal* x,
	integer* incx, doublereal* tau);

/* Subroutine */
CLAPACK_API
int dlarfp_(integer* n, doublereal* alpha, doublereal* x,
	integer* incx, doublereal* tau);

/* Subroutine */
CLAPACK_API
int dlarft_(char* direct, char* storev, integer* n, integer*
	k, doublereal* v, integer* ldv, doublereal* tau, doublereal* t,
	integer* ldt);

/* Subroutine */
CLAPACK_API
int dlarfx_(char* side, integer* m, integer* n, doublereal*
	v, doublereal* tau, doublereal* c__, integer* ldc, doublereal* work);

/* Subroutine */
CLAPACK_API
int dlargv_(integer* n, doublereal* x, integer* incx,
	doublereal* y, integer* incy, doublereal* c__, integer* incc);

/* Subroutine */
CLAPACK_API
int dlarrv_(integer* n, doublereal* vl, doublereal* vu,
	doublereal* d__, doublereal* l, doublereal* pivmin, integer* isplit,
	integer* m, integer* dol, integer* dou, doublereal* minrgp,
	doublereal* rtol1, doublereal* rtol2, doublereal* w, doublereal* werr,
	doublereal* wgap, integer* iblock, integer* indexw, doublereal* gers,
	doublereal* z__, integer* ldz, integer* isuppz, doublereal* work,
	integer* iwork, integer* info);

/* Subroutine */
CLAPACK_API
int dlartv_(integer* n, doublereal* x, integer* incx,
	doublereal* y, integer* incy, doublereal* c__, doublereal* s, integer
	* incc);

/* Subroutine */
CLAPACK_API
int dlarz_(char* side, integer* m, integer* n, integer* l,
	doublereal* v, integer* incv, doublereal* tau, doublereal* c__,
	integer* ldc, doublereal* work);

/* Subroutine */
CLAPACK_API
int dlarzb_(char* side, char* trans, char* direct, char*
	storev, integer* m, integer* n, integer* k, integer* l, doublereal* v,
	integer* ldv, doublereal* t, integer* ldt, doublereal* c__, integer*
	ldc, doublereal* work, integer* ldwork);

/* Subroutine */
CLAPACK_API
int dlarzt_(char* direct, char* storev, integer* n, integer*
	k, doublereal* v, integer* ldv, doublereal* tau, doublereal* t,
	integer* ldt);

/* Subroutine */
CLAPACK_API
int dlaswp_(integer* n, doublereal* a, integer* lda, integer
	* k1, integer* k2, integer* ipiv, integer* incx);

/* Subroutine */
CLAPACK_API
int dlasy2_(logical* ltranl, logical* ltranr, integer* isgn,
	integer* n1, integer* n2, doublereal* tl, integer* ldtl, doublereal*
	tr, integer* ldtr, doublereal* b, integer* ldb, doublereal* scale,
	doublereal* x, integer* ldx, doublereal* xnorm, integer* info);

/* Subroutine */
CLAPACK_API
int dlasyf_(char* uplo, integer* n, integer* nb, integer* kb,
	doublereal* a, integer* lda, integer* ipiv, doublereal* w, integer*
	ldw, integer* info);

/* Subroutine */
CLAPACK_API
int dlat2s_(char* uplo, integer* n, doublereal* a, integer*
	lda, real* sa, integer* ldsa, integer* info);

/* Subroutine */
CLAPACK_API
int dlatbs_(char* uplo, char* trans, char* diag, char*
	normin, integer* n, integer* kd, doublereal* ab, integer* ldab,
	doublereal* x, doublereal* scale, doublereal* cnorm, integer* info);

/* Subroutine */
CLAPACK_API
int dlatdf_(integer* ijob, integer* n, doublereal* z__,
	integer* ldz, doublereal* rhs, doublereal* rdsum, doublereal* rdscal,
	integer* ipiv, integer* jpiv);

/* Subroutine */
CLAPACK_API
int dlatps_(char* uplo, char* trans, char* diag, char*
	normin, integer* n, doublereal* ap, doublereal* x, doublereal* scale,
	doublereal* cnorm, integer* info);

/* Subroutine */
CLAPACK_API
int dlatrd_(char* uplo, integer* n, integer* nb, doublereal*
	a, integer* lda, doublereal* e, doublereal* tau, doublereal* w,
	integer* ldw);

/* Subroutine */
CLAPACK_API
int dlatrs_(char* uplo, char* trans, char* diag, char*
	normin, integer* n, doublereal* a, integer* lda, doublereal* x,
	doublereal* scale, doublereal* cnorm, integer* info);

/* Subroutine */
CLAPACK_API
int dlatrz_(integer* m, integer* n, integer* l, doublereal*
	a, integer* lda, doublereal* tau, doublereal* work);

/* Subroutine */
CLAPACK_API
int dlatzm_(char* side, integer* m, integer* n, doublereal*
	v, integer* incv, doublereal* tau, doublereal* c1, doublereal* c2,
	integer* ldc, doublereal* work);

/* Subroutine */
CLAPACK_API
int dlauu2_(char* uplo, integer* n, doublereal* a, integer*
	lda, integer* info);

/* Subroutine */
CLAPACK_API
int dlauum_(char* uplo, integer* n, doublereal* a, integer*
	lda, integer* info);

/* Subroutine */
CLAPACK_API
int dopgtr_(char* uplo, integer* n, doublereal* ap,
	doublereal* tau, doublereal* q, integer* ldq, doublereal* work,
	integer* info);

/* Subroutine */
CLAPACK_API
int dopmtr_(char* side, char* uplo, char* trans, integer* m,
	integer* n, doublereal* ap, doublereal* tau, doublereal* c__, integer
	* ldc, doublereal* work, integer* info);

/* Subroutine */
CLAPACK_API
int dorg2l_(integer* m, integer* n, integer* k, doublereal*
	a, integer* lda, doublereal* tau, doublereal* work, integer* info);

/* Subroutine */
CLAPACK_API
int dorg2r_(integer* m, integer* n, integer* k, doublereal*
	a, integer* lda, doublereal* tau, doublereal* work, integer* info);

/* Subroutine */
CLAPACK_API
int dorgbr_(char* vect, integer* m, integer* n, integer* k,
	doublereal* a, integer* lda, doublereal* tau, doublereal* work,
	integer* lwork, integer* info);

/* Subroutine */
CLAPACK_API
int dorghr_(integer* n, integer* ilo, integer* ihi,
	doublereal* a, integer* lda, doublereal* tau, doublereal* work,
	integer* lwork, integer* info);

/* Subroutine */
CLAPACK_API
int dorgl2_(integer* m, integer* n, integer* k, doublereal*
	a, integer* lda, doublereal* tau, doublereal* work, integer* info);

/* Subroutine */
CLAPACK_API
int dorglq_(integer* m, integer* n, integer* k, doublereal*
	a, integer* lda, doublereal* tau, doublereal* work, integer* lwork,
	integer* info);

/* Subroutine */
CLAPACK_API
int dorgql_(integer* m, integer* n, integer* k, doublereal*
	a, integer* lda, doublereal* tau, doublereal* work, integer* lwork,
	integer* info);

/* Subroutine */
CLAPACK_API
int dorgqr_(integer* m, integer* n, integer* k, doublereal*
	a, integer* lda, doublereal* tau, doublereal* work, integer* lwork,
	integer* info);

/* Subroutine */
CLAPACK_API
int dorgr2_(integer* m, integer* n, integer* k, doublereal*
	a, integer* lda, doublereal* tau, doublereal* work, integer* info);

/* Subroutine */
CLAPACK_API
int dorgrq_(integer* m, integer* n, integer* k, doublereal*
	a, integer* lda, doublereal* tau, doublereal* work, integer* lwork,
	integer* info);

/* Subroutine */
CLAPACK_API
int dorgtr_(char* uplo, integer* n, doublereal* a, integer*
	lda, doublereal* tau, doublereal* work, integer* lwork, integer* info);

/* Subroutine */
CLAPACK_API
int dorm2l_(char* side, char* trans, integer* m, integer* n,
	integer* k, doublereal* a, integer* lda, doublereal* tau, doublereal*
	c__, integer* ldc, doublereal* work, integer* info);

/* Subroutine */
CLAPACK_API
int dorm2r_(char* side, char* trans, integer* m, integer* n,
	integer* k, doublereal* a, integer* lda, doublereal* tau, doublereal*
	c__, integer* ldc, doublereal* work, integer* info);

/* Subroutine */
CLAPACK_API
int dormbr_(char* vect, char* side, char* trans, integer* m,
	integer* n, integer* k, doublereal* a, integer* lda, doublereal* tau,
	doublereal* c__, integer* ldc, doublereal* work, integer* lwork,
	integer* info);

/* Subroutine */
CLAPACK_API
int dormhr_(char* side, char* trans, integer* m, integer* n,
	integer* ilo, integer* ihi, doublereal* a, integer* lda, doublereal*
	tau, doublereal* c__, integer* ldc, doublereal* work, integer* lwork,
	integer* info);

/* Subroutine */
CLAPACK_API
int dorml2_(char* side, char* trans, integer* m, integer* n,
	integer* k, doublereal* a, integer* lda, doublereal* tau, doublereal*
	c__, integer* ldc, doublereal* work, integer* info);

/* Subroutine */
CLAPACK_API
int dormlq_(char* side, char* trans, integer* m, integer* n,
	integer* k, doublereal* a, integer* lda, doublereal* tau, doublereal*
	c__, integer* ldc, doublereal* work, integer* lwork, integer* info);

/* Subroutine */
CLAPACK_API
int dormql_(char* side, char* trans, integer* m, integer* n,
	integer* k, doublereal* a, integer* lda, doublereal* tau, doublereal*
	c__, integer* ldc, doublereal* work, integer* lwork, integer* info);

/* Subroutine */
CLAPACK_API
int dormqr_(char* side, char* trans, integer* m, integer* n,
	integer* k, doublereal* a, integer* lda, doublereal* tau, doublereal*
	c__, integer* ldc, doublereal* work, integer* lwork, integer* info);

/* Subroutine */
CLAPACK_API
int dormr2_(char* side, char* trans, integer* m, integer* n,
	integer* k, doublereal* a, integer* lda, doublereal* tau, doublereal*
	c__, integer* ldc, doublereal* work, integer* info);

/* Subroutine */
CLAPACK_API
int dormr3_(char* side, char* trans, integer* m, integer* n,
	integer* k, integer* l, doublereal* a, integer* lda, doublereal* tau,
	doublereal* c__, integer* ldc, doublereal* work, integer* info);

/* Subroutine */
CLAPACK_API
int dormrq_(char* side, char* trans, integer* m, integer* n,
	integer* k, doublereal* a, integer* lda, doublereal* tau, doublereal*
	c__, integer* ldc, doublereal* work, integer* lwork, integer* info);

/* Subroutine */
CLAPACK_API
int dormrz_(char* side, char* trans, integer* m, integer* n,
	integer* k, integer* l, doublereal* a, integer* lda, doublereal* tau,
	doublereal* c__, integer* ldc, doublereal* work, integer* lwork,
	integer* info);

/* Subroutine */
CLAPACK_API
int dormtr_(char* side, char* uplo, char* trans, integer* m,
	integer* n, doublereal* a, integer* lda, doublereal* tau, doublereal*
	c__, integer* ldc, doublereal* work, integer* lwork, integer* info);

/* Subroutine */
CLAPACK_API
int dpbcon_(char* uplo, integer* n, integer* kd, doublereal*
	ab, integer* ldab, doublereal* anorm, doublereal* rcond, doublereal*
	work, integer* iwork, integer* info);

/* Subroutine */
CLAPACK_API
int dpbequ_(char* uplo, integer* n, integer* kd, doublereal*
	ab, integer* ldab, doublereal* s, doublereal* scond, doublereal* amax,
	integer* info);

/* Subroutine */
CLAPACK_API
int dpbrfs_(char* uplo, integer* n, integer* kd, integer*
	nrhs, doublereal* ab, integer* ldab, doublereal* afb, integer* ldafb,
	doublereal* b, integer* ldb, doublereal* x, integer* ldx, doublereal*
	ferr, doublereal* berr, doublereal* work, integer* iwork, integer*
	info);

/* Subroutine */
CLAPACK_API
int dpbstf_(char* uplo, integer* n, integer* kd, doublereal*
	ab, integer* ldab, integer* info);

/* Subroutine */
CLAPACK_API
int dpbsv_(char* uplo, integer* n, integer* kd, integer*
	nrhs, doublereal* ab, integer* ldab, doublereal* b, integer* ldb,
	integer* info);

/* Subroutine */
CLAPACK_API
int dpbsvx_(char* fact, char* uplo, integer* n, integer* kd,
	integer* nrhs, doublereal* ab, integer* ldab, doublereal* afb,
	integer* ldafb, char* equed, doublereal* s, doublereal* b, integer*
	ldb, doublereal* x, integer* ldx, doublereal* rcond, doublereal* ferr,
	doublereal* berr, doublereal* work, integer* iwork, integer* info);

/* Subroutine */
CLAPACK_API
int dpbtf2_(char* uplo, integer* n, integer* kd, doublereal*
	ab, integer* ldab, integer* info);

/* Subroutine */
CLAPACK_API
int dpbtrf_(char* uplo, integer* n, integer* kd, doublereal*
	ab, integer* ldab, integer* info);

/* Subroutine */
CLAPACK_API
int dpbtrs_(char* uplo, integer* n, integer* kd, integer*
	nrhs, doublereal* ab, integer* ldab, doublereal* b, integer* ldb,
	integer* info);

/* Subroutine */
CLAPACK_API
int dpftrf_(char* transr, char* uplo, integer* n, doublereal
	* a, integer* info);

/* Subroutine */
CLAPACK_API
int dpftri_(char* transr, char* uplo, integer* n, doublereal
	* a, integer* info);

/* Subroutine */
CLAPACK_API
int dpftrs_(char* transr, char* uplo, integer* n, integer*
	nrhs, doublereal* a, doublereal* b, integer* ldb, integer* info);

/* Subroutine */
CLAPACK_API
int dpocon_(char* uplo, integer* n, doublereal* a, integer*
	lda, doublereal* anorm, doublereal* rcond, doublereal* work, integer*
	iwork, integer* info);

/* Subroutine */
CLAPACK_API
int dpoequ_(integer* n, doublereal* a, integer* lda,
	doublereal* s, doublereal* scond, doublereal* amax, integer* info);

/* Subroutine */
CLAPACK_API
int dpoequb_(integer* n, doublereal* a, integer* lda,
	doublereal* s, doublereal* scond, doublereal* amax, integer* info);

/* Subroutine */
CLAPACK_API
int dporfs_(char* uplo, integer* n, integer* nrhs,
	doublereal* a, integer* lda, doublereal* af, integer* ldaf,
	doublereal* b, integer* ldb, doublereal* x, integer* ldx, doublereal*
	ferr, doublereal* berr, doublereal* work, integer* iwork, integer*
	info);

/* Subroutine */
CLAPACK_API
int dposv_(char* uplo, integer* n, integer* nrhs, doublereal
	* a, integer* lda, doublereal* b, integer* ldb, integer* info);

/* Subroutine */
CLAPACK_API
int dposvx_(char* fact, char* uplo, integer* n, integer*
	nrhs, doublereal* a, integer* lda, doublereal* af, integer* ldaf,
	char* equed, doublereal* s, doublereal* b, integer* ldb, doublereal*
	x, integer* ldx, doublereal* rcond, doublereal* ferr, doublereal*
	berr, doublereal* work, integer* iwork, integer* info);

/* Subroutine */
CLAPACK_API
int dpotf2_(char* uplo, integer* n, doublereal* a, integer*
	lda, integer* info);

/* Subroutine */
CLAPACK_API
int dpotrf_(char* uplo, integer* n, doublereal* a, integer*
	lda, integer* info);

/* Subroutine */
CLAPACK_API
int dpotri_(char* uplo, integer* n, doublereal* a, integer*
	lda, integer* info);

/* Subroutine */
CLAPACK_API
int dpotrs_(char* uplo, integer* n, integer* nrhs,
	doublereal* a, integer* lda, doublereal* b, integer* ldb, integer*
	info);

/* Subroutine */
CLAPACK_API
int dppcon_(char* uplo, integer* n, doublereal* ap,
	doublereal* anorm, doublereal* rcond, doublereal* work, integer*
	iwork, integer* info);

/* Subroutine */
CLAPACK_API
int dppequ_(char* uplo, integer* n, doublereal* ap,
	doublereal* s, doublereal* scond, doublereal* amax, integer* info);

/* Subroutine */
CLAPACK_API
int dpprfs_(char* uplo, integer* n, integer* nrhs,
	doublereal* ap, doublereal* afp, doublereal* b, integer* ldb,
	doublereal* x, integer* ldx, doublereal* ferr, doublereal* berr,
	doublereal* work, integer* iwork, integer* info);

/* Subroutine */
CLAPACK_API
int dppsv_(char* uplo, integer* n, integer* nrhs, doublereal
	* ap, doublereal* b, integer* ldb, integer* info);

/* Subroutine */
CLAPACK_API
int dppsvx_(char* fact, char* uplo, integer* n, integer*
	nrhs, doublereal* ap, doublereal* afp, char* equed, doublereal* s,
	doublereal* b, integer* ldb, doublereal* x, integer* ldx, doublereal*
	rcond, doublereal* ferr, doublereal* berr, doublereal* work, integer*
	iwork, integer* info);

/* Subroutine */
CLAPACK_API
int dpptrf_(char* uplo, integer* n, doublereal* ap, integer*
	info);

/* Subroutine */
CLAPACK_API
int dpptri_(char* uplo, integer* n, doublereal* ap, integer*
	info);

/* Subroutine */
CLAPACK_API
int dpptrs_(char* uplo, integer* n, integer* nrhs,
	doublereal* ap, doublereal* b, integer* ldb, integer* info);

/* Subroutine */
CLAPACK_API
int dpstf2_(char* uplo, integer* n, doublereal* a, integer*
	lda, integer* piv, integer* rank, doublereal* tol, doublereal* work,
	integer* info);

/* Subroutine */
CLAPACK_API
int dpstrf_(char* uplo, integer* n, doublereal* a, integer*
	lda, integer* piv, integer* rank, doublereal* tol, doublereal* work,
	integer* info);

/* Subroutine */
CLAPACK_API
int dptcon_(integer* n, doublereal* d__, doublereal* e,
	doublereal* anorm, doublereal* rcond, doublereal* work, integer* info);

/* Subroutine */
CLAPACK_API
int dpteqr_(char* compz, integer* n, doublereal* d__,
	doublereal* e, doublereal* z__, integer* ldz, doublereal* work,
	integer* info);

/* Subroutine */
CLAPACK_API
int dptrfs_(integer* n, integer* nrhs, doublereal* d__,
	doublereal* e, doublereal* df, doublereal* ef, doublereal* b, integer
	* ldb, doublereal* x, integer* ldx, doublereal* ferr, doublereal* berr,
	doublereal* work, integer* info);

/* Subroutine */
CLAPACK_API
int dptsv_(integer* n, integer* nrhs, doublereal* d__,
	doublereal* e, doublereal* b, integer* ldb, integer* info);

/* Subroutine */
CLAPACK_API
int dptsvx_(char* fact, integer* n, integer* nrhs,
	doublereal* d__, doublereal* e, doublereal* df, doublereal* ef,
	doublereal* b, integer* ldb, doublereal* x, integer* ldx, doublereal*
	rcond, doublereal* ferr, doublereal* berr, doublereal* work, integer*
	info);

/* Subroutine */
CLAPACK_API
int dpttrs_(integer* n, integer* nrhs, doublereal* d__,
	doublereal* e, doublereal* b, integer* ldb, integer* info);

/* Subroutine */
CLAPACK_API
int dptts2_(integer* n, integer* nrhs, doublereal* d__,
	doublereal* e, doublereal* b, integer* ldb);

/* Subroutine */
CLAPACK_API
int drscl_(integer* n, doublereal* sa, doublereal* sx,
	integer* incx);

/* Subroutine */
CLAPACK_API
int dsbev_(char* jobz, char* uplo, integer* n, integer* kd,
	doublereal* ab, integer* ldab, doublereal* w, doublereal* z__,
	integer* ldz, doublereal* work, integer* info);

/* Subroutine */
CLAPACK_API
int dsbevd_(char* jobz, char* uplo, integer* n, integer* kd,
	doublereal* ab, integer* ldab, doublereal* w, doublereal* z__,
	integer* ldz, doublereal* work, integer* lwork, integer* iwork,
	integer* liwork, integer* info);

/* Subroutine */
CLAPACK_API
int dsbevx_(char* jobz, char* range, char* uplo, integer* n,
	integer* kd, doublereal* ab, integer* ldab, doublereal* q, integer*
	ldq, doublereal* vl, doublereal* vu, integer* il, integer* iu,
	doublereal* abstol, integer* m, doublereal* w, doublereal* z__,
	integer* ldz, doublereal* work, integer* iwork, integer* ifail,
	integer* info);

/* Subroutine */
CLAPACK_API
int dsbgst_(char* vect, char* uplo, integer* n, integer* ka,
	integer* kb, doublereal* ab, integer* ldab, doublereal* bb, integer*
	ldbb, doublereal* x, integer* ldx, doublereal* work, integer* info);

/* Subroutine */
CLAPACK_API
int dsbgv_(char* jobz, char* uplo, integer* n, integer* ka,
	integer* kb, doublereal* ab, integer* ldab, doublereal* bb, integer*
	ldbb, doublereal* w, doublereal* z__, integer* ldz, doublereal* work,
	integer* info);

/* Subroutine */
CLAPACK_API
int dsbgvd_(char* jobz, char* uplo, integer* n, integer* ka,
	integer* kb, doublereal* ab, integer* ldab, doublereal* bb, integer*
	ldbb, doublereal* w, doublereal* z__, integer* ldz, doublereal* work,
	integer* lwork, integer* iwork, integer* liwork, integer* info);

/* Subroutine */
CLAPACK_API
int dsbgvx_(char* jobz, char* range, char* uplo, integer* n,
	integer* ka, integer* kb, doublereal* ab, integer* ldab, doublereal*
	bb, integer* ldbb, doublereal* q, integer* ldq, doublereal* vl,
	doublereal* vu, integer* il, integer* iu, doublereal* abstol, integer
	* m, doublereal* w, doublereal* z__, integer* ldz, doublereal* work,
	integer* iwork, integer* ifail, integer* info);

/* Subroutine */
CLAPACK_API
int dsbtrd_(char* vect, char* uplo, integer* n, integer* kd,
	doublereal* ab, integer* ldab, doublereal* d__, doublereal* e,
	doublereal* q, integer* ldq, doublereal* work, integer* info);

/* Subroutine */
CLAPACK_API
int dsfrk_(char* transr, char* uplo, char* trans, integer* n,
	integer* k, doublereal* alpha, doublereal* a, integer* lda,
	doublereal* beta, doublereal* c__);

/* Subroutine */
CLAPACK_API int dsgesv_(integer* n, integer* nrhs, doublereal* a,
	integer* lda, integer* ipiv, doublereal* b, integer* ldb, doublereal*
	x, integer* ldx, doublereal* work, real* swork, integer* iter,
	integer* info);

/* Subroutine */
CLAPACK_API
int dspcon_(char* uplo, integer* n, doublereal* ap, integer*
	ipiv, doublereal* anorm, doublereal* rcond, doublereal* work, integer
	* iwork, integer* info);

/* Subroutine */
CLAPACK_API
int dspev_(char* jobz, char* uplo, integer* n, doublereal*
	ap, doublereal* w, doublereal* z__, integer* ldz, doublereal* work,
	integer* info);

/* Subroutine */
CLAPACK_API int dspevd_(char* jobz, char* uplo, integer* n, doublereal*
	ap, doublereal* w, doublereal* z__, integer* ldz, doublereal* work,
	integer* lwork, integer* iwork, integer* liwork, integer* info);

/* Subroutine */
CLAPACK_API
int dspevx_(char* jobz, char* range, char* uplo, integer* n,
	doublereal* ap, doublereal* vl, doublereal* vu, integer* il, integer*
	iu, doublereal* abstol, integer* m, doublereal* w, doublereal* z__,
	integer* ldz, doublereal* work, integer* iwork, integer* ifail,
	integer* info);

/* Subroutine */
CLAPACK_API
int dspgst_(integer* itype, char* uplo, integer* n,
	doublereal* ap, doublereal* bp, integer* info);

/* Subroutine */
CLAPACK_API
int dspgv_(integer* itype, char* jobz, char* uplo, integer*
	n, doublereal* ap, doublereal* bp, doublereal* w, doublereal* z__,
	integer* ldz, doublereal* work, integer* info);

/* Subroutine */
CLAPACK_API
int dspgvd_(integer* itype, char* jobz, char* uplo, integer*
	n, doublereal* ap, doublereal* bp, doublereal* w, doublereal* z__,
	integer* ldz, doublereal* work, integer* lwork, integer* iwork,
	integer* liwork, integer* info);

/* Subroutine */
CLAPACK_API
int dspgvx_(integer* itype, char* jobz, char* range, char*
	uplo, integer* n, doublereal* ap, doublereal* bp, doublereal* vl,
	doublereal* vu, integer* il, integer* iu, doublereal* abstol, integer
	* m, doublereal* w, doublereal* z__, integer* ldz, doublereal* work,
	integer* iwork, integer* ifail, integer* info);

/* Subroutine */
CLAPACK_API
int dsposv_(char* uplo, integer* n, integer* nrhs,
	doublereal* a, integer* lda, doublereal* b, integer* ldb, doublereal*
	x, integer* ldx, doublereal* work, real* swork, integer* iter,
	integer* info);

/* Subroutine */
CLAPACK_API
int dsprfs_(char* uplo, integer* n, integer* nrhs,
	doublereal* ap, doublereal* afp, integer* ipiv, doublereal* b,
	integer* ldb, doublereal* x, integer* ldx, doublereal* ferr,
	doublereal* berr, doublereal* work, integer* iwork, integer* info);

/* Subroutine */
CLAPACK_API
int dspsv_(char* uplo, integer* n, integer* nrhs, doublereal
	* ap, integer* ipiv, doublereal* b, integer* ldb, integer* info);

/* Subroutine */
CLAPACK_API
int dspsvx_(char* fact, char* uplo, integer* n, integer*
	nrhs, doublereal* ap, doublereal* afp, integer* ipiv, doublereal* b,
	integer* ldb, doublereal* x, integer* ldx, doublereal* rcond,
	doublereal* ferr, doublereal* berr, doublereal* work, integer* iwork,
	integer* info);

/* Subroutine */
CLAPACK_API
int dsptrd_(char* uplo, integer* n, doublereal* ap,
	doublereal* d__, doublereal* e, doublereal* tau, integer* info);

/* Subroutine */
CLAPACK_API
int dsptrf_(char* uplo, integer* n, doublereal* ap, integer*
	ipiv, integer* info);

/* Subroutine */
CLAPACK_API
int dsptri_(char* uplo, integer* n, doublereal* ap, integer*
	ipiv, doublereal* work, integer* info);

/* Subroutine */
CLAPACK_API
int dsptrs_(char* uplo, integer* n, integer* nrhs,
	doublereal* ap, integer* ipiv, doublereal* b, integer* ldb, integer*
	info);

/* Subroutine */
CLAPACK_API
int dstegr_(char* jobz, char* range, integer* n, doublereal*
	d__, doublereal* e, doublereal* vl, doublereal* vu, integer* il,
	integer* iu, doublereal* abstol, integer* m, doublereal* w,
	doublereal* z__, integer* ldz, integer* isuppz, doublereal* work,
	integer* lwork, integer* iwork, integer* liwork, integer* info);

/* Subroutine */
CLAPACK_API
int dstein_(integer* n, doublereal* d__, doublereal* e,
	integer* m, doublereal* w, integer* iblock, integer* isplit,
	doublereal* z__, integer* ldz, doublereal* work, integer* iwork,
	integer* ifail, integer* info);

/* Subroutine */
CLAPACK_API
int dstemr_(char* jobz, char* range, integer* n, doublereal*
	d__, doublereal* e, doublereal* vl, doublereal* vu, integer* il,
	integer* iu, integer* m, doublereal* w, doublereal* z__, integer* ldz,
	integer* nzc, integer* isuppz, logical* tryrac, doublereal* work,
	integer* lwork, integer* iwork, integer* liwork, integer* info);

/* Subroutine */
CLAPACK_API
int dstev_(char* jobz, integer* n, doublereal* d__,
	doublereal* e, doublereal* z__, integer* ldz, doublereal* work,
	integer* info);

/* Subroutine */
CLAPACK_API
int dstevd_(char* jobz, integer* n, doublereal* d__,
	doublereal* e, doublereal* z__, integer* ldz, doublereal* work,
	integer* lwork, integer* iwork, integer* liwork, integer* info);

/* Subroutine */
CLAPACK_API
int dstevr_(char* jobz, char* range, integer* n, doublereal*
	d__, doublereal* e, doublereal* vl, doublereal* vu, integer* il,
	integer* iu, doublereal* abstol, integer* m, doublereal* w,
	doublereal* z__, integer* ldz, integer* isuppz, doublereal* work,
	integer* lwork, integer* iwork, integer* liwork, integer* info);

/* Subroutine */
CLAPACK_API
int dstevx_(char* jobz, char* range, integer* n, doublereal*
	d__, doublereal* e, doublereal* vl, doublereal* vu, integer* il,
	integer* iu, doublereal* abstol, integer* m, doublereal* w,
	doublereal* z__, integer* ldz, doublereal* work, integer* iwork,
	integer* ifail, integer* info);

/* Subroutine */
CLAPACK_API
int dsycon_(char* uplo, integer* n, doublereal* a, integer*
	lda, integer* ipiv, doublereal* anorm, doublereal* rcond, doublereal*
	work, integer* iwork, integer* info);

/* Subroutine */
CLAPACK_API
int dsyequb_(char* uplo, integer* n, doublereal* a, integer*
	lda, doublereal* s, doublereal* scond, doublereal* amax, doublereal*
	work, integer* info);

/* Subroutine */
CLAPACK_API
int dsyev_(char* jobz, char* uplo, integer* n, doublereal* a,
	integer* lda, doublereal* w, doublereal* work, integer* lwork,
	integer* info);

/* Subroutine */
CLAPACK_API
int dsyevd_(char* jobz, char* uplo, integer* n, doublereal*
	a, integer* lda, doublereal* w, doublereal* work, integer* lwork,
	integer* iwork, integer* liwork, integer* info);

/* Subroutine */
CLAPACK_API
int dsyevr_(char* jobz, char* range, char* uplo, integer* n,
	doublereal* a, integer* lda, doublereal* vl, doublereal* vu, integer*
	il, integer* iu, doublereal* abstol, integer* m, doublereal* w,
	doublereal* z__, integer* ldz, integer* isuppz, doublereal* work,
	integer* lwork, integer* iwork, integer* liwork, integer* info);

/* Subroutine */
CLAPACK_API
int dsyevx_(char* jobz, char* range, char* uplo, integer* n,
	doublereal* a, integer* lda, doublereal* vl, doublereal* vu, integer*
	il, integer* iu, doublereal* abstol, integer* m, doublereal* w,
	doublereal* z__, integer* ldz, doublereal* work, integer* lwork,
	integer* iwork, integer* ifail, integer* info);

/* Subroutine */
CLAPACK_API
int dsygs2_(integer* itype, char* uplo, integer* n,
	doublereal* a, integer* lda, doublereal* b, integer* ldb, integer*
	info);

/* Subroutine */
CLAPACK_API
int dsygst_(integer* itype, char* uplo, integer* n,
	doublereal* a, integer* lda, doublereal* b, integer* ldb, integer*
	info);

/* Subroutine */
CLAPACK_API
int dsygv_(integer* itype, char* jobz, char* uplo, integer*
	n, doublereal* a, integer* lda, doublereal* b, integer* ldb,
	doublereal* w, doublereal* work, integer* lwork, integer* info);

/* Subroutine */
CLAPACK_API
int dsygvd_(integer* itype, char* jobz, char* uplo, integer*
	n, doublereal* a, integer* lda, doublereal* b, integer* ldb,
	doublereal* w, doublereal* work, integer* lwork, integer* iwork,
	integer* liwork, integer* info);

/* Subroutine */
CLAPACK_API
int dsygvx_(integer* itype, char* jobz, char* range, char*
	uplo, integer* n, doublereal* a, integer* lda, doublereal* b, integer
	* ldb, doublereal* vl, doublereal* vu, integer* il, integer* iu,
	doublereal* abstol, integer* m, doublereal* w, doublereal* z__,
	integer* ldz, doublereal* work, integer* lwork, integer* iwork,
	integer* ifail, integer* info);

/* Subroutine */
CLAPACK_API
int dsyrfs_(char* uplo, integer* n, integer* nrhs,
	doublereal* a, integer* lda, doublereal* af, integer* ldaf, integer*
	ipiv, doublereal* b, integer* ldb, doublereal* x, integer* ldx,
	doublereal* ferr, doublereal* berr, doublereal* work, integer* iwork,
	integer* info);

/* Subroutine */
CLAPACK_API
int dsysv_(char* uplo, integer* n, integer* nrhs, doublereal
	* a, integer* lda, integer* ipiv, doublereal* b, integer* ldb,
	doublereal* work, integer* lwork, integer* info);

/* Subroutine */
CLAPACK_API
int dsysvx_(char* fact, char* uplo, integer* n, integer*
	nrhs, doublereal* a, integer* lda, doublereal* af, integer* ldaf,
	integer* ipiv, doublereal* b, integer* ldb, doublereal* x, integer*
	ldx, doublereal* rcond, doublereal* ferr, doublereal* berr,
	doublereal* work, integer* lwork, integer* iwork, integer* info);

/* Subroutine */
CLAPACK_API
int dsytd2_(char* uplo, integer* n, doublereal* a, integer*
	lda, doublereal* d__, doublereal* e, doublereal* tau, integer* info);

/* Subroutine */
CLAPACK_API
int dsytf2_(char* uplo, integer* n, doublereal* a, integer*
	lda, integer* ipiv, integer* info);

/* Subroutine */
CLAPACK_API
int dsytrd_(char* uplo, integer* n, doublereal* a, integer*
	lda, doublereal* d__, doublereal* e, doublereal* tau, doublereal*
	work, integer* lwork, integer* info);

/* Subroutine */
CLAPACK_API
int dsytrf_(char* uplo, integer* n, doublereal* a, integer*
	lda, integer* ipiv, doublereal* work, integer* lwork, integer* info);

/* Subroutine */
CLAPACK_API
int dsytri_(char* uplo, integer* n, doublereal* a, integer*
	lda, integer* ipiv, doublereal* work, integer* info);

/* Subroutine */
CLAPACK_API
int dsytrs_(char* uplo, integer* n, integer* nrhs,
	doublereal* a, integer* lda, integer* ipiv, doublereal* b, integer*
	ldb, integer* info);

/* Subroutine */
CLAPACK_API
int dtbcon_(char* norm, char* uplo, char* diag, integer* n,
	integer* kd, doublereal* ab, integer* ldab, doublereal* rcond,
	doublereal* work, integer* iwork, integer* info);

/* Subroutine */
CLAPACK_API
int dtbrfs_(char* uplo, char* trans, char* diag, integer* n,
	integer* kd, integer* nrhs, doublereal* ab, integer* ldab, doublereal
	* b, integer* ldb, doublereal* x, integer* ldx, doublereal* ferr,
	doublereal* berr, doublereal* work, integer* iwork, integer* info);

/* Subroutine */
CLAPACK_API
int dtbtrs_(char* uplo, char* trans, char* diag, integer* n,
	integer* kd, integer* nrhs, doublereal* ab, integer* ldab, doublereal
	* b, integer* ldb, integer* info);

/* Subroutine */
CLAPACK_API
int dtfsm_(char* transr, char* side, char* uplo, char* trans,
	char* diag, integer* m, integer* n, doublereal* alpha, doublereal* a,
	doublereal* b, integer* ldb);

/* Subroutine */
CLAPACK_API
int dtftri_(char* transr, char* uplo, char* diag, integer* n,
	doublereal* a, integer* info);

/* Subroutine */
CLAPACK_API int dtfttp_(char* transr, char* uplo, integer* n, doublereal
	* arf, doublereal* ap, integer* info);

/* Subroutine */
CLAPACK_API
int dtfttr_(char* transr, char* uplo, integer* n, doublereal
	* arf, doublereal* a, integer* lda, integer* info);

/* Subroutine */
CLAPACK_API
int dtgevc_(char* side, char* howmny, logical* select,
	integer* n, doublereal* s, integer* lds, doublereal* p, integer* ldp,
	doublereal* vl, integer* ldvl, doublereal* vr, integer* ldvr, integer
	* mm, integer* m, doublereal* work, integer* info);

/* Subroutine */
CLAPACK_API int dtgex2_(logical* wantq, logical* wantz, integer* n,
	doublereal* a, integer* lda, doublereal* b, integer* ldb, doublereal*
	q, integer* ldq, doublereal* z__, integer* ldz, integer* j1, integer*
	n1, integer* n2, doublereal* work, integer* lwork, integer* info);

/* Subroutine */
CLAPACK_API
int dtgexc_(logical* wantq, logical* wantz, integer* n,
	doublereal* a, integer* lda, doublereal* b, integer* ldb, doublereal*
	q, integer* ldq, doublereal* z__, integer* ldz, integer* ifst,
	integer* ilst, doublereal* work, integer* lwork, integer* info);

/* Subroutine */
CLAPACK_API
int dtgsen_(integer* ijob, logical* wantq, logical* wantz,
	logical* select, integer* n, doublereal* a, integer* lda, doublereal*
	b, integer* ldb, doublereal* alphar, doublereal* alphai, doublereal*
	beta, doublereal* q, integer* ldq, doublereal* z__, integer* ldz,
	integer* m, doublereal* pl, doublereal* pr, doublereal* dif,
	doublereal* work, integer* lwork, integer* iwork, integer* liwork,
	integer* info);

/* Subroutine */
CLAPACK_API
int dtgsja_(char* jobu, char* jobv, char* jobq, integer* m,
	integer* p, integer* n, integer* k, integer* l, doublereal* a,
	integer* lda, doublereal* b, integer* ldb, doublereal* tola,
	doublereal* tolb, doublereal* alpha, doublereal* beta, doublereal* u,
	integer* ldu, doublereal* v, integer* ldv, doublereal* q, integer*
	ldq, doublereal* work, integer* ncycle, integer* info);

/* Subroutine */
CLAPACK_API
int dtgsna_(char* job, char* howmny, logical* select,
	integer* n, doublereal* a, integer* lda, doublereal* b, integer* ldb,
	doublereal* vl, integer* ldvl, doublereal* vr, integer* ldvr,
	doublereal* s, doublereal* dif, integer* mm, integer* m, doublereal*
	work, integer* lwork, integer* iwork, integer* info);

/* Subroutine */
CLAPACK_API
int dtgsy2_(char* trans, integer* ijob, integer* m, integer*
	n, doublereal* a, integer* lda, doublereal* b, integer* ldb,
	doublereal* c__, integer* ldc, doublereal* d__, integer* ldd,
	doublereal* e, integer* lde, doublereal* f, integer* ldf, doublereal*
	scale, doublereal* rdsum, doublereal* rdscal, integer* iwork, integer
	* pq, integer* info);

/* Subroutine */
CLAPACK_API
int dtgsyl_(char* trans, integer* ijob, integer* m, integer*
	n, doublereal* a, integer* lda, doublereal* b, integer* ldb,
	doublereal* c__, integer* ldc, doublereal* d__, integer* ldd,
	doublereal* e, integer* lde, doublereal* f, integer* ldf, doublereal*
	scale, doublereal* dif, doublereal* work, integer* lwork, integer*
	iwork, integer* info);

/* Subroutine */
CLAPACK_API
int dtpcon_(char* norm, char* uplo, char* diag, integer* n,
	doublereal* ap, doublereal* rcond, doublereal* work, integer* iwork,
	integer* info);

/* Subroutine */
CLAPACK_API
int dtprfs_(char* uplo, char* trans, char* diag, integer* n,
	integer* nrhs, doublereal* ap, doublereal* b, integer* ldb,
	doublereal* x, integer* ldx, doublereal* ferr, doublereal* berr,
	doublereal* work, integer* iwork, integer* info);

/* Subroutine */
CLAPACK_API
int dtptri_(char* uplo, char* diag, integer* n, doublereal*
	ap, integer* info);

/* Subroutine */
CLAPACK_API int dtptrs_(char* uplo, char* trans, char* diag, integer* n,
	integer* nrhs, doublereal* ap, doublereal* b, integer* ldb, integer*
	info);

/* Subroutine */
CLAPACK_API
int dtpttf_(char* transr, char* uplo, integer* n, doublereal
	* ap, doublereal* arf, integer* info);

/* Subroutine */
CLAPACK_API
int dtpttr_(char* uplo, integer* n, doublereal* ap,
	doublereal* a, integer* lda, integer* info);

/* Subroutine */
CLAPACK_API
int dtrcon_(char* norm, char* uplo, char* diag, integer* n,
	doublereal* a, integer* lda, doublereal* rcond, doublereal* work,
	integer* iwork, integer* info);

/* Subroutine */
CLAPACK_API
int dtrevc_(char* side, char* howmny, logical* select,
	integer* n, doublereal* t, integer* ldt, doublereal* vl, integer*
	ldvl, doublereal* vr, integer* ldvr, integer* mm, integer* m,
	doublereal* work, integer* info);

/* Subroutine */
CLAPACK_API
int dtrexc_(char* compq, integer* n, doublereal* t, integer*
	ldt, doublereal* q, integer* ldq, integer* ifst, integer* ilst,
	doublereal* work, integer* info);

/* Subroutine */
CLAPACK_API
int dtrrfs_(char* uplo, char* trans, char* diag, integer* n,
	integer* nrhs, doublereal* a, integer* lda, doublereal* b, integer*
	ldb, doublereal* x, integer* ldx, doublereal* ferr, doublereal* berr,
	doublereal* work, integer* iwork, integer* info);

/* Subroutine */
CLAPACK_API
int dtrsen_(char* job, char* compq, logical* select, integer
	* n, doublereal* t, integer* ldt, doublereal* q, integer* ldq,
	doublereal* wr, doublereal* wi, integer* m, doublereal* s, doublereal
	* sep, doublereal* work, integer* lwork, integer* iwork, integer*
	liwork, integer* info);

/* Subroutine */
CLAPACK_API
int dtrsna_(char* job, char* howmny, logical* select,
	integer* n, doublereal* t, integer* ldt, doublereal* vl, integer*
	ldvl, doublereal* vr, integer* ldvr, doublereal* s, doublereal* sep,
	integer* mm, integer* m, doublereal* work, integer* ldwork, integer*
	iwork, integer* info);

/* Subroutine */
CLAPACK_API
int dtrsyl_(char* trana, char* tranb, integer* isgn, integer
	* m, integer* n, doublereal* a, integer* lda, doublereal* b, integer*
	ldb, doublereal* c__, integer* ldc, doublereal* scale, integer* info);

/* Subroutine */
CLAPACK_API
int dtrti2_(char* uplo, char* diag, integer* n, doublereal*
	a, integer* lda, integer* info);

/* Subroutine */
CLAPACK_API
int dtrtri_(char* uplo, char* diag, integer* n, doublereal*
	a, integer* lda, integer* info);

/* Subroutine */
CLAPACK_API
int dtrtrs_(char* uplo, char* trans, char* diag, integer* n,
	integer* nrhs, doublereal* a, integer* lda, doublereal* b, integer*
	ldb, integer* info);

/* Subroutine */
CLAPACK_API
int dtrttf_(char* transr, char* uplo, integer* n, doublereal
	* a, integer* lda, doublereal* arf, integer* info);

/* Subroutine */
CLAPACK_API
int dtrttp_(char* uplo, integer* n, doublereal* a, integer*
	lda, doublereal* ap, integer* info);

/* Subroutine */
CLAPACK_API
int dtzrqf_(integer* m, integer* n, doublereal* a, integer*
	lda, doublereal* tau, integer* info);

/* Subroutine */
CLAPACK_API
int dtzrzf_(integer* m, integer* n, doublereal* a, integer*
	lda, doublereal* tau, doublereal* work, integer* lwork, integer* info);




#pragma endregion

#pragma region DXLASRC -- Double precision real LAPACK routines using extra precision (dxla)
//#       
//set(DXLASRC dgesvxx.c dgerfsx.c dla_gerfsx_extended.c dla_geamv.c
//	dla_gercond.c dla_rpvgrw.c dsysvxx.c dsyrfsx.c
//	dla_syrfsx_extended.c dla_syamv.c dla_syrcond.c dla_syrpvgrw.c
//	dposvxx.c dporfsx.c dla_porfsx_extended.c dla_porcond.c
//	dla_porpvgrw.c dgbsvxx.c dgbrfsx.c dla_gbrfsx_extended.c
//	dla_gbamv.c dla_gbrcond.c dla_gbrpvgrw.c dla_lin_berr.c dlarscl2.c
//	dlascl2.c dla_wwaddw.c)
//

#pragma endregion



#pragma region SCLAUX -- Auxiliary routines called from both REALand COMPLEX (alasc)


/* Subroutine */
CLAPACK_API
int sbdsdc_(char* uplo, char* compq, integer* n, real* d__,
	real* e, real* u, integer* ldu, real* vt, integer* ldvt, real* q,
	integer* iq, real* work, integer* iwork, integer* info);

/* Subroutine */
CLAPACK_API
int sbdsqr_(char* uplo, integer* n, integer* ncvt, integer*
	nru, integer* ncc, real* d__, real* e, real* vt, integer* ldvt, real*
	u, integer* ldu, real* c__, integer* ldc, real* work, integer* info);

/* Subroutine */
CLAPACK_API
int sdisna_(char* job, integer* m, integer* n, real* d__,
	real* sep, integer* info);

CLAPACK_API
doublereal second_();

CLAPACK_API
logical sisnan_(real* sin__);

/* Subroutine */
CLAPACK_API
int slabad_(real* small, real* large);

/* Subroutine */
CLAPACK_API
int slacpy_(char* uplo, integer* m, integer* n, real* a,
	integer* lda, real* b, integer* ldb);

/* Subroutine */
CLAPACK_API
int sladiv_(real* a, real* b, real* c__, real* d__, real* p,
	real* q);

/* Subroutine */
CLAPACK_API
int slae2_(real* a, real* b, real* c__, real* rt1, real* rt2);

/* Subroutine */
CLAPACK_API
int slaebz_(integer* ijob, integer* nitmax, integer* n,
	integer* mmax, integer* minp, integer* nbmin, real* abstol, real*
	reltol, real* pivmin, real* d__, real* e, real* e2, integer* nval,
	real* ab, real* c__, integer* mout, integer* nab, real* work, integer
	* iwork, integer* info);

/* Subroutine */
CLAPACK_API
int slaed0_(integer* icompq, integer* qsiz, integer* n, real
	* d__, real* e, real* q, integer* ldq, real* qstore, integer* ldqs,
	real* work, integer* iwork, integer* info);

/* Subroutine */
CLAPACK_API
int slaed1_(integer* n, real* d__, real* q, integer* ldq,
	integer* indxq, real* rho, integer* cutpnt, real* work, integer*
	iwork, integer* info);

/* Subroutine */
CLAPACK_API
int slaed2_(integer* k, integer* n, integer* n1, real* d__,
	real* q, integer* ldq, integer* indxq, real* rho, real* z__, real*
	dlamda, real* w, real* q2, integer* indx, integer* indxc, integer*
	indxp, integer* coltyp, integer* info);

/* Subroutine */
CLAPACK_API
int slaed3_(integer* k, integer* n, integer* n1, real* d__,
	real* q, integer* ldq, real* rho, real* dlamda, real* q2, integer*
	indx, integer* ctot, real* w, real* s, integer* info);

/* Subroutine */
CLAPACK_API
int slaed4_(integer* n, integer* i__, real* d__, real* z__,
	real* delta, real* rho, real* dlam, integer* info);

/* Subroutine */
CLAPACK_API
int slaed5_(integer* i__, real* d__, real* z__, real* delta,
	real* rho, real* dlam);

/* Subroutine */
CLAPACK_API
int slaed6_(integer* kniter, logical* orgati, real* rho,
	real* d__, real* z__, real* finit, real* tau, integer* info);

/* Subroutine */
CLAPACK_API
int slaed7_(integer* icompq, integer* n, integer* qsiz,
	integer* tlvls, integer* curlvl, integer* curpbm, real* d__, real* q,
	integer* ldq, integer* indxq, real* rho, integer* cutpnt, real*
	qstore, integer* qptr, integer* prmptr, integer* perm, integer*
	givptr, integer* givcol, real* givnum, real* work, integer* iwork,
	integer* info);

/* Subroutine */
CLAPACK_API
int slaed8_(integer* icompq, integer* k, integer* n, integer
	* qsiz, real* d__, real* q, integer* ldq, integer* indxq, real* rho,
	integer* cutpnt, real* z__, real* dlamda, real* q2, integer* ldq2,
	real* w, integer* perm, integer* givptr, integer* givcol, real*
	givnum, integer* indxp, integer* indx, integer* info);

/* Subroutine */
CLAPACK_API
int slaed9_(integer* k, integer* kstart, integer* kstop,
	integer* n, real* d__, real* q, integer* ldq, real* rho, real* dlamda,
	real* w, real* s, integer* lds, integer* info);

/* Subroutine */
CLAPACK_API
int slaeda_(integer* n, integer* tlvls, integer* curlvl,
	integer* curpbm, integer* prmptr, integer* perm, integer* givptr,
	integer* givcol, real* givnum, real* q, integer* qptr, real* z__,
	real* ztemp, integer* info);

/* Subroutine */
CLAPACK_API
int slaev2_(real* a, real* b, real* c__, real* rt1, real*
	rt2, real* cs1, real* sn1);

/* Subroutine */
CLAPACK_API
int slagtf_(integer* n, real* a, real* lambda, real* b, real
	* c__, real* tol, real* d__, integer* in, integer* info);

/* Subroutine */
CLAPACK_API
int slagts_(integer* job, integer* n, real* a, real* b, real
	* c__, real* d__, integer* in, real* y, real* tol, integer* info);

CLAPACK_API
logical slaisnan_(real* sin1, real* sin2);

CLAPACK_API
doublereal slamch_(char* cmach);

/* Subroutine */
CLAPACK_API
int slamrg_(integer* n1, integer* n2, real* a, integer*
	strd1, integer* strd2, integer* index);

CLAPACK_API
integer slaneg_(integer* n, real* d__, real* lld, real* sigma, real* pivmin,
	integer* r__);

CLAPACK_API
doublereal slanst_(char* norm, integer* n, real* d__, real* e);

CLAPACK_API
doublereal slapy2_(real* x, real* y);

CLAPACK_API
doublereal slapy3_(real* x, real* y, real* z__);

/* Subroutine */
CLAPACK_API
int slarnv_(integer* idist, integer* iseed, integer* n, real* x);

/* Subroutine */
CLAPACK_API
int slarra_(integer* n, real* d__, real* e, real* e2, real*
	spltol, real* tnrm, integer* nsplit, integer* isplit, integer* info);

/* Subroutine */
CLAPACK_API
int slarrb_(integer* n, real* d__, real* lld, integer*
	ifirst, integer* ilast, real* rtol1, real* rtol2, integer* offset,
	real* w, real* wgap, real* werr, real* work, integer* iwork, real*
	pivmin, real* spdiam, integer* twist, integer* info);

/* Subroutine */
CLAPACK_API
int slarrc_(char* jobt, integer* n, real* vl, real* vu, real
	* d__, real* e, real* pivmin, integer* eigcnt, integer* lcnt, integer*
	rcnt, integer* info);

/* Subroutine */
CLAPACK_API
int slarrd_(char* range, char* order, integer* n, real* vl,
	real* vu, integer* il, integer* iu, real* gers, real* reltol, real*
	d__, real* e, real* e2, real* pivmin, integer* nsplit, integer*
	isplit, integer* m, real* w, real* werr, real* wl, real* wu, integer*
	iblock, integer* indexw, real* work, integer* iwork, integer* info);

/* Subroutine */
CLAPACK_API
int slarre_(char* range, integer* n, real* vl, real* vu,
	integer* il, integer* iu, real* d__, real* e, real* e2, real* rtol1,
	real* rtol2, real* spltol, integer* nsplit, integer* isplit, integer*
	m, real* w, real* werr, real* wgap, integer* iblock, integer* indexw,
	real* gers, real* pivmin, real* work, integer* iwork, integer* info);

/* Subroutine */
CLAPACK_API
int slarrf_(integer* n, real* d__, real* l, real* ld,
	integer* clstrt, integer* clend, real* w, real* wgap, real* werr,
	real* spdiam, real* clgapl, real* clgapr, real* pivmin, real* sigma,
	real* dplus, real* lplus, real* work, integer* info);

/* Subroutine */
CLAPACK_API
int slarrj_(integer* n, real* d__, real* e2, integer* ifirst,
	integer* ilast, real* rtol, integer* offset, real* w, real* werr,
	real* work, integer* iwork, real* pivmin, real* spdiam, integer* info);

/* Subroutine */
CLAPACK_API
int slarrk_(integer* n, integer* iw, real* gl, real* gu,
	real* d__, real* e2, real* pivmin, real* reltol, real* w, real* werr,
	integer* info);

/* Subroutine */
CLAPACK_API
int slarrr_(integer* n, real* d__, real* e, integer* info);

/* Subroutine */
CLAPACK_API
int slartg_(real* f, real* g, real* cs, real* sn, real* r__);

/* Subroutine */
CLAPACK_API
int slaruv_(integer* iseed, integer* n, real* x);

/* Subroutine */
CLAPACK_API
int slas2_(real* f, real* g, real* h__, real* ssmin, real* ssmax);

/* Subroutine */
CLAPACK_API
int slascl_(char* type__, integer* kl, integer* ku, real*
	cfrom, real* cto, integer* m, integer* n, real* a, integer* lda,
	integer* info);

/* Subroutine */
CLAPACK_API
int slasd0_(integer* n, integer* sqre, real* d__, real* e,
	real* u, integer* ldu, real* vt, integer* ldvt, integer* smlsiz,
	integer* iwork, real* work, integer* info);

/* Subroutine */
CLAPACK_API
int slasd1_(integer* nl, integer* nr, integer* sqre, real*
	d__, real* alpha, real* beta, real* u, integer* ldu, real* vt,
	integer* ldvt, integer* idxq, integer* iwork, real* work, integer*
	info);

/* Subroutine */
CLAPACK_API
int slasd2_(integer* nl, integer* nr, integer* sqre, integer
	* k, real* d__, real* z__, real* alpha, real* beta, real* u, integer*
	ldu, real* vt, integer* ldvt, real* dsigma, real* u2, integer* ldu2,
	real* vt2, integer* ldvt2, integer* idxp, integer* idx, integer* idxc,
	integer* idxq, integer* coltyp, integer* info);

/* Subroutine */
CLAPACK_API
int slasd3_(integer* nl, integer* nr, integer* sqre, integer
	* k, real* d__, real* q, integer* ldq, real* dsigma, real* u, integer*
	ldu, real* u2, integer* ldu2, real* vt, integer* ldvt, real* vt2,
	integer* ldvt2, integer* idxc, integer* ctot, real* z__, integer*
	info);

/* Subroutine */
CLAPACK_API
int slasd4_(integer* n, integer* i__, real* d__, real* z__,
	real* delta, real* rho, real* sigma, real* work, integer* info);

/* Subroutine */
CLAPACK_API
int slasd5_(integer* i__, real* d__, real* z__, real* delta,
	real* rho, real* dsigma, real* work);

/* Subroutine */
CLAPACK_API
int slasd6_(integer* icompq, integer* nl, integer* nr,
	integer* sqre, real* d__, real* vf, real* vl, real* alpha, real* beta,
	integer* idxq, integer* perm, integer* givptr, integer* givcol,
	integer* ldgcol, real* givnum, integer* ldgnum, real* poles, real*
	difl, real* difr, real* z__, integer* k, real* c__, real* s, real*
	work, integer* iwork, integer* info);

/* Subroutine */
CLAPACK_API
int slasd7_(integer* icompq, integer* nl, integer* nr,
	integer* sqre, integer* k, real* d__, real* z__, real* zw, real* vf,
	real* vfw, real* vl, real* vlw, real* alpha, real* beta, real* dsigma,
	integer* idx, integer* idxp, integer* idxq, integer* perm, integer*
	givptr, integer* givcol, integer* ldgcol, real* givnum, integer*
	ldgnum, real* c__, real* s, integer* info);

/* Subroutine */
CLAPACK_API
int slasd8_(integer* icompq, integer* k, real* d__, real*
	z__, real* vf, real* vl, real* difl, real* difr, integer* lddifr,
	real* dsigma, real* work, integer* info);

/* Subroutine */
CLAPACK_API
int slasda_(integer* icompq, integer* smlsiz, integer* n,
	integer* sqre, real* d__, real* e, real* u, integer* ldu, real* vt,
	integer* k, real* difl, real* difr, real* z__, real* poles, integer*
	givptr, integer* givcol, integer* ldgcol, integer* perm, real* givnum,
	real* c__, real* s, real* work, integer* iwork, integer* info);

/* Subroutine */
CLAPACK_API
int slasdq_(char* uplo, integer* sqre, integer* n, integer*
	ncvt, integer* nru, integer* ncc, real* d__, real* e, real* vt,
	integer* ldvt, real* u, integer* ldu, real* c__, integer* ldc, real*
	work, integer* info);

/* Subroutine */
CLAPACK_API
int slasdt_(integer* n, integer* lvl, integer* nd, integer*
	inode, integer* ndiml, integer* ndimr, integer* msub);

/* Subroutine */
CLAPACK_API
int slaset_(char* uplo, integer* m, integer* n, real* alpha,
	real* beta, real* a, integer* lda);

/* Subroutine */
CLAPACK_API
int slasq1_(integer* n, real* d__, real* e, real* work,
	integer* info);

/* Subroutine */
CLAPACK_API
int slasq2_(integer* n, real* z__, integer* info);

/* Subroutine */
CLAPACK_API
int slasq3_(integer* i0, integer* n0, real* z__, integer* pp,
	real* dmin__, real* sigma, real* desig, real* qmax, integer* nfail,
	integer* iter, integer* ndiv, logical* ieee, integer* ttype, real*
	dmin1, real* dmin2, real* dn, real* dn1, real* dn2, real* g, real*
	tau);

/* Subroutine */
CLAPACK_API
int slasq4_(integer* i0, integer* n0, real* z__, integer* pp,
	integer* n0in, real* dmin__, real* dmin1, real* dmin2, real* dn,
	real* dn1, real* dn2, real* tau, integer* ttype, real* g);

/* Subroutine */
CLAPACK_API
int slasq5_(integer* i0, integer* n0, real* z__, integer* pp,
	real* tau, real* dmin__, real* dmin1, real* dmin2, real* dn, real*
	dnm1, real* dnm2, logical* ieee);

/* Subroutine */
CLAPACK_API
int slasq6_(integer* i0, integer* n0, real* z__, integer* pp,
	real* dmin__, real* dmin1, real* dmin2, real* dn, real* dnm1, real*
	dnm2);

/* Subroutine */
CLAPACK_API
int slasr_(char* side, char* pivot, char* direct, integer* m,
	integer* n, real* c__, real* s, real* a, integer* lda);

/* Subroutine */
CLAPACK_API
int slasrt_(char* id, integer* n, real* d__, integer* info);

/* Subroutine */
CLAPACK_API
int slassq_(integer* n, real* x, integer* incx, real* scale,
	real* sumsq);

/* Subroutine */
CLAPACK_API
int slasv2_(real* f, real* g, real* h__, real* ssmin, real*
	ssmax, real* snr, real* csr, real* snl, real* csl);

CLAPACK_API
integer smaxloc_(real* a, integer* dimm);

/* Subroutine */
CLAPACK_API
int spttrf_(integer* n, real* d__, real* e, integer* info);

/* Subroutine */
CLAPACK_API
int sstebz_(char* range, char* order, integer* n, real* vl,
	real* vu, integer* il, integer* iu, real* abstol, real* d__, real* e,
	integer* m, integer* nsplit, real* w, integer* iblock, integer*
	isplit, real* work, integer* iwork, integer* info);

/* Subroutine */
CLAPACK_API
int sstedc_(char* compz, integer* n, real* d__, real* e,
	real* z__, integer* ldz, real* work, integer* lwork, integer* iwork,
	integer* liwork, integer* info);

/* Subroutine */
CLAPACK_API
int ssteqr_(char* compz, integer* n, real* d__, real* e,
	real* z__, integer* ldz, real* work, integer* info);

/* Subroutine */
CLAPACK_API
int ssterf_(integer* n, real* d__, real* e, integer* info);





#pragma endregion

#pragma region SLASRC -- Single precision real LAPACK routines (sla)
//#       
//set(SLASRC
//	sgbbrd.c sgbcon.c sgbequ.c sgbrfs.c sgbsv.c
//	sgbsvx.c sgbtf2.c sgbtrf.c sgbtrs.c sgebak.c sgebal.c sgebd2.c
//	sgebrd.c sgecon.c sgeequ.c sgees.c  sgeesx.c sgeev.c  sgeevx.c
//	sgegs.c  sgegv.c  sgehd2.c sgehrd.c sgelq2.c sgelqf.c
//	sgels.c  sgelsd.c sgelss.c sgelsx.c sgelsy.c sgeql2.c sgeqlf.c
//	sgeqp3.c sgeqpf.c sgeqr2.c sgeqrf.c sgerfs.c sgerq2.c sgerqf.c
//	sgesc2.c sgesdd.c sgesv.c  sgesvd.c sgesvx.c sgetc2.c sgetf2.c
//	sgetrf.c sgetri.c
//	sgetrs.c sggbak.c sggbal.c sgges.c  sggesx.c sggev.c  sggevx.c
//	sggglm.c sgghrd.c sgglse.c sggqrf.c
//	sggrqf.c sggsvd.c sggsvp.c sgtcon.c sgtrfs.c sgtsv.c
//	sgtsvx.c sgttrf.c sgttrs.c sgtts2.c shgeqz.c
//	shsein.c shseqr.c slabrd.c slacon.c slacn2.c
//	slaein.c slaexc.c slag2.c  slags2.c slagtm.c slagv2.c slahqr.c
//	slahrd.c slahr2.c slaic1.c slaln2.c slals0.c slalsa.c slalsd.c
//	slangb.c slange.c slangt.c slanhs.c slansb.c slansp.c
//	slansy.c slantb.c slantp.c slantr.c slanv2.c
//	slapll.c slapmt.c
//	slaqgb.c slaqge.c slaqp2.c slaqps.c slaqsb.c slaqsp.c slaqsy.c
//	slaqr0.c slaqr1.c slaqr2.c slaqr3.c slaqr4.c slaqr5.c
//	slaqtr.c slar1v.c slar2v.c ilaslr.c ilaslc.c
//	slarf.c  slarfb.c slarfg.c slarft.c slarfx.c slargv.c
//	slarrv.c slartv.c slarfp.c
//	slarz.c  slarzb.c slarzt.c slaswp.c slasy2.c slasyf.c
//	slatbs.c slatdf.c slatps.c slatrd.c slatrs.c slatrz.c slatzm.c
//	slauu2.c slauum.c sopgtr.c sopmtr.c sorg2l.c sorg2r.c
//	sorgbr.c sorghr.c sorgl2.c sorglq.c sorgql.c sorgqr.c sorgr2.c
//	sorgrq.c sorgtr.c sorm2l.c sorm2r.c
//	sormbr.c sormhr.c sorml2.c sormlq.c sormql.c sormqr.c sormr2.c
//	sormr3.c sormrq.c sormrz.c sormtr.c spbcon.c spbequ.c spbrfs.c
//	spbstf.c spbsv.c  spbsvx.c
//	spbtf2.c spbtrf.c spbtrs.c spocon.c spoequ.c sporfs.c sposv.c
//	sposvx.c spotf2.c spotrf.c spotri.c spotrs.c spstrf.c spstf2.c
//	sppcon.c sppequ.c
//	spprfs.c sppsv.c  sppsvx.c spptrf.c spptri.c spptrs.c sptcon.c
//	spteqr.c sptrfs.c sptsv.c  sptsvx.c spttrs.c sptts2.c srscl.c
//	ssbev.c  ssbevd.c ssbevx.c ssbgst.c ssbgv.c  ssbgvd.c ssbgvx.c
//	ssbtrd.c sspcon.c sspev.c  sspevd.c sspevx.c sspgst.c
//	sspgv.c  sspgvd.c sspgvx.c ssprfs.c sspsv.c  sspsvx.c ssptrd.c
//	ssptrf.c ssptri.c ssptrs.c sstegr.c sstein.c sstev.c  sstevd.c sstevr.c
//	sstevx.c ssycon.c ssyev.c  ssyevd.c ssyevr.c ssyevx.c ssygs2.c
//	ssygst.c ssygv.c  ssygvd.c ssygvx.c ssyrfs.c ssysv.c  ssysvx.c
//	ssytd2.c ssytf2.c ssytrd.c ssytrf.c ssytri.c ssytrs.c stbcon.c
//	stbrfs.c stbtrs.c stgevc.c stgex2.c stgexc.c stgsen.c
//	stgsja.c stgsna.c stgsy2.c stgsyl.c stpcon.c stprfs.c stptri.c
//	stptrs.c
//	strcon.c strevc.c strexc.c strrfs.c strsen.c strsna.c strsyl.c
//	strti2.c strtri.c strtrs.c stzrqf.c stzrzf.c sstemr.c
//	slansf.c spftrf.c spftri.c spftrs.c ssfrk.c stfsm.c stftri.c stfttp.c
//	stfttr.c stpttf.c stpttr.c strttf.c strttp.c
//	sgejsv.c  sgesvj.c  sgsvj0.c  sgsvj1.c
//	sgeequb.c ssyequb.c spoequb.c sgbequb.c)

#pragma endregion

#pragma region SXLASRC -- Single precision real LAPACK routines using extra precision (sxla)
//#       
//set(SXLASRC sgesvxx.c sgerfsx.c sla_gerfsx_extended.c sla_geamv.c
//	sla_gercond.c sla_rpvgrw.c ssysvxx.c ssyrfsx.c
//	sla_syrfsx_extended.c sla_syamv.c sla_syrcond.c sla_syrpvgrw.c
//	sposvxx.c sporfsx.c sla_porfsx_extended.c sla_porcond.c
//	sla_porpvgrw.c sgbsvxx.c sgbrfsx.c sla_gbrfsx_extended.c
//	sla_gbamv.c sla_gbrcond.c sla_gbrpvgrw.c sla_lin_berr.c slarscl2.c
//	slascl2.c sla_wwaddw.c)


#pragma endregion



#pragma region ZLASRC -- Double precision complex LAPACK routines (zla)
//#       
//set(ZLASRC
//	zbdsqr.c zgbbrd.c zgbcon.c zgbequ.c zgbrfs.c zgbsv.c  zgbsvx.c
//	zgbtf2.c zgbtrf.c zgbtrs.c zgebak.c zgebal.c zgebd2.c zgebrd.c
//	zgecon.c zgeequ.c zgees.c  zgeesx.c zgeev.c  zgeevx.c
//	zgegs.c  zgegv.c  zgehd2.c zgehrd.c zgelq2.c zgelqf.c
//	zgels.c  zgelsd.c zgelss.c zgelsx.c zgelsy.c zgeql2.c zgeqlf.c zgeqp3.c
//	zgeqpf.c zgeqr2.c zgeqrf.c zgerfs.c zgerq2.c zgerqf.c
//	zgesc2.c zgesdd.c zgesv.c  zgesvd.c zgesvx.c zgetc2.c zgetf2.c zgetrf.c
//	zgetri.c zgetrs.c
//	zggbak.c zggbal.c zgges.c  zggesx.c zggev.c  zggevx.c zggglm.c
//	zgghrd.c zgglse.c zggqrf.c zggrqf.c
//	zggsvd.c zggsvp.c
//	zgtcon.c zgtrfs.c zgtsv.c  zgtsvx.c zgttrf.c zgttrs.c zgtts2.c zhbev.c
//	zhbevd.c zhbevx.c zhbgst.c zhbgv.c  zhbgvd.c zhbgvx.c zhbtrd.c
//	zhecon.c zheev.c  zheevd.c zheevr.c zheevx.c zhegs2.c zhegst.c
//	zhegv.c  zhegvd.c zhegvx.c zherfs.c zhesv.c  zhesvx.c zhetd2.c
//	zhetf2.c zhetrd.c
//	zhetrf.c zhetri.c zhetrs.c zhgeqz.c zhpcon.c zhpev.c  zhpevd.c
//	zhpevx.c zhpgst.c zhpgv.c  zhpgvd.c zhpgvx.c zhprfs.c zhpsv.c
//	zhpsvx.c
//	zhptrd.c zhptrf.c zhptri.c zhptrs.c zhsein.c zhseqr.c zlabrd.c
//	zlacgv.c zlacon.c zlacn2.c zlacp2.c zlacpy.c zlacrm.c zlacrt.c zladiv.c
//	zlaed0.c zlaed7.c zlaed8.c
//	zlaein.c zlaesy.c zlaev2.c zlags2.c zlagtm.c
//	zlahef.c zlahqr.c
//	zlahrd.c zlahr2.c zlaic1.c zlals0.c zlalsa.c zlalsd.c zlangb.c zlange.c
//	zlangt.c zlanhb.c
//	zlanhe.c
//	zlanhp.c zlanhs.c zlanht.c zlansb.c zlansp.c zlansy.c zlantb.c
//	zlantp.c zlantr.c zlapll.c zlapmt.c zlaqgb.c zlaqge.c
//	zlaqhb.c zlaqhe.c zlaqhp.c zlaqp2.c zlaqps.c zlaqsb.c
//	zlaqr0.c zlaqr1.c zlaqr2.c zlaqr3.c zlaqr4.c zlaqr5.c
//	zlaqsp.c zlaqsy.c zlar1v.c zlar2v.c ilazlr.c ilazlc.c
//	zlarcm.c zlarf.c  zlarfb.c
//	zlarfg.c zlarft.c zlarfp.c
//	zlarfx.c zlargv.c zlarnv.c zlarrv.c zlartg.c zlartv.c
//	zlarz.c  zlarzb.c zlarzt.c zlascl.c zlaset.c zlasr.c
//	zlassq.c zlaswp.c zlasyf.c
//	zlatbs.c zlatdf.c zlatps.c zlatrd.c zlatrs.c zlatrz.c zlatzm.c zlauu2.c
//	zlauum.c zpbcon.c zpbequ.c zpbrfs.c zpbstf.c zpbsv.c
//	zpbsvx.c zpbtf2.c zpbtrf.c zpbtrs.c zpocon.c zpoequ.c zporfs.c
//	zposv.c  zposvx.c zpotf2.c zpotrf.c zpotri.c zpotrs.c zpstrf.c zpstf2.c
//	zppcon.c zppequ.c zpprfs.c zppsv.c  zppsvx.c zpptrf.c zpptri.c zpptrs.c
//	zptcon.c zpteqr.c zptrfs.c zptsv.c  zptsvx.c zpttrf.c zpttrs.c zptts2.c
//	zrot.c   zspcon.c zspmv.c  zspr.c   zsprfs.c zspsv.c
//	zspsvx.c zsptrf.c zsptri.c zsptrs.c zdrscl.c zstedc.c
//	zstegr.c zstein.c zsteqr.c zsycon.c zsymv.c
//	zsyr.c   zsyrfs.c zsysv.c  zsysvx.c zsytf2.c zsytrf.c zsytri.c
//	zsytrs.c ztbcon.c ztbrfs.c ztbtrs.c ztgevc.c ztgex2.c
//	ztgexc.c ztgsen.c ztgsja.c ztgsna.c ztgsy2.c ztgsyl.c ztpcon.c
//	ztprfs.c ztptri.c
//	ztptrs.c ztrcon.c ztrevc.c ztrexc.c ztrrfs.c ztrsen.c ztrsna.c
//	ztrsyl.c ztrti2.c ztrtri.c ztrtrs.c ztzrqf.c ztzrzf.c zung2l.c
//	zung2r.c zungbr.c zunghr.c zungl2.c zunglq.c zungql.c zungqr.c zungr2.c
//	zungrq.c zungtr.c zunm2l.c zunm2r.c zunmbr.c zunmhr.c zunml2.c
//	zunmlq.c zunmql.c zunmqr.c zunmr2.c zunmr3.c zunmrq.c zunmrz.c
//	zunmtr.c zupgtr.c
//	zupmtr.c izmax1.c dzsum1.c zstemr.c
//	zcgesv.c zcposv.c zlag2c.c clag2z.c zlat2c.c
//	zhfrk.c ztfttp.c zlanhf.c zpftrf.c zpftri.c zpftrs.c ztfsm.c ztftri.c
//	ztfttr.c ztpttf.c ztpttr.c ztrttf.c ztrttp.c
//	zgeequb.c zgbequb.c zsyequb.c zpoequb.c zheequb.c)


doublereal dzsum1_(integer* n, doublecomplex* cx, integer* incx);

#pragma endregion

#pragma region ZXLASRC -- Double precision complex LAPACK routines using extra precision (zxla)
//#       
//set(ZXLASRC  zgesvxx.c zgerfsx.c zla_gerfsx_extended.c zla_geamv.c
//	zla_gercond_c.c zla_gercond_x.c zla_rpvgrw.c zsysvxx.c zsyrfsx.c
//	zla_syrfsx_extended.c zla_syamv.c zla_syrcond_c.c zla_syrcond_x.c
//	zla_syrpvgrw.c zposvxx.c zporfsx.c zla_porfsx_extended.c
//	zla_porcond_c.c zla_porcond_x.c zla_porpvgrw.c zgbsvxx.c zgbrfsx.c
//	zla_gbrfsx_extended.c zla_gbamv.c zla_gbrcond_c.c zla_gbrcond_x.c
//	zla_gbrpvgrw.c zhesvxx.c zherfsx.c zla_herfsx_extended.c
//	zla_heamv.c zla_hercond_c.c zla_hercond_x.c zla_herpvgrw.c
//	zla_lin_berr.c zlarscl2.c zlascl2.c zla_wwaddw.c)

#pragma endregion




#pragma region IXX


CLAPACK_API
integer icmax1_(integer* n, complex* cx, integer* incx);

CLAPACK_API
integer ilaclc_(integer* m, integer* n, complex* a, integer* lda);

CLAPACK_API
integer ilaclr_(integer* m, integer* n, complex* a, integer* lda);

CLAPACK_API
integer iladlc_(integer* m, integer* n, doublereal* a, integer* lda);

CLAPACK_API
integer iladlr_(integer* m, integer* n, doublereal* a, integer* lda);

CLAPACK_API
integer ilaslc_(integer* m, integer* n, real* a, integer* lda);

CLAPACK_API
integer ilaslr_(integer* m, integer* n, real* a, integer* lda);

CLAPACK_API
integer ilazlc_(integer* m, integer* n, doublecomplex* a, integer* lda);

CLAPACK_API
integer ilazlr_(integer* m, integer* n, doublecomplex* a, integer* lda);

CLAPACK_API
integer izmax1_(integer* n, doublecomplex* cx, integer* incx);


/* Subroutine */
CLAPACK_API
int dlamc1_(integer* beta, integer* t, logical* rnd, logical
	* ieee1);

/* Subroutine */
CLAPACK_API
int slamc1_(integer* beta, integer* t, logical* rnd, logical
	* ieee1);

/* Subroutine */
CLAPACK_API
int slamc2_(integer* beta, integer* t, logical* rnd, real*
	eps, integer* emin, real* rmin, integer* emax, real* rmax);

CLAPACK_API
doublereal slamc3_(real* a, real* b);

/* Subroutine */
CLAPACK_API
int slamc4_(integer* emin, real* start, integer* base);

/* Subroutine */
CLAPACK_API
int slamc5_(integer* beta, integer* p, integer* emin,
	logical* ieee, integer* emax, real* rmax);

/* Subroutine */
CLAPACK_API
int dlamc1_(integer* beta, integer* t, logical* rnd, logical
	* ieee1);

/* Subroutine */
CLAPACK_API
int dlamc2_(integer* beta, integer* t, logical* rnd,
	doublereal* eps, integer* emin, doublereal* rmin, integer* emax,
	doublereal* rmax);

CLAPACK_API
doublereal dlamc3_(doublereal* a, doublereal* b);

/* Subroutine */
CLAPACK_API
int dlamc4_(integer* emin, doublereal* start, integer* base);

/* Subroutine */
CLAPACK_API
int dlamc5_(integer* beta, integer* p, integer* emin,
	logical* ieee, integer* emax, doublereal* rmax);

CLAPACK_API
integer ilaenv_(integer* ispec, char* name__, char* opts, integer* n1,
	integer* n2, integer* n3, integer* n4);

#pragma endregion




#pragma region UNDEFINE


/*
CLAPACK_API
int cgbrfsx_(char* trans, char* equed, integer* n, integer*
	kl, integer* ku, integer* nrhs, complex* ab, integer* ldab, complex*
	afb, integer* ldafb, integer* ipiv, real* r__, real* c__, complex* b,
	integer* ldb, complex* x, integer* ldx, real* rcond, real* berr,
	integer* n_err_bnds__, real* err_bnds_norm__, real* err_bnds_comp__,
	integer* nparams, real* params, complex* work, real* rwork, integer*
	info);


CLAPACK_API
int cgbsvxx_(char* fact, char* trans, integer* n, integer*
	kl, integer* ku, integer* nrhs, complex* ab, integer* ldab, complex*
	afb, integer* ldafb, integer* ipiv, char* equed, real* r__, real* c__,
	complex* b, integer* ldb, complex* x, integer* ldx, real* rcond,
	real* rpvgrw, real* berr, integer* n_err_bnds__, real*
	err_bnds_norm__, real* err_bnds_comp__, integer* nparams, real*
	params, complex* work, real* rwork, integer* info);


CLAPACK_API
int cgerfsx_(char* trans, char* equed, integer* n, integer*
	nrhs, complex* a, integer* lda, complex* af, integer* ldaf, integer*
	ipiv, real* r__, real* c__, complex* b, integer* ldb, complex* x,
	integer* ldx, real* rcond, real* berr, integer* n_err_bnds__, real*
	err_bnds_norm__, real* err_bnds_comp__, integer* nparams, real*
	params, complex* work, real* rwork, integer* info);


CLAPACK_API
int cgesvxx_(char* fact, char* trans, integer* n, integer*
	nrhs, complex* a, integer* lda, complex* af, integer* ldaf, integer*
	ipiv, char* equed, real* r__, real* c__, complex* b, integer* ldb,
	complex* x, integer* ldx, real* rcond, real* rpvgrw, real* berr,
	integer* n_err_bnds__, real* err_bnds_norm__, real* err_bnds_comp__,
	integer* nparams, real* params, complex* work, real* rwork, integer*
	info);


CLAPACK_API
int dgbrfsx_(char* trans, char* equed, integer* n, integer*
	kl, integer* ku, integer* nrhs, doublereal* ab, integer* ldab,
	doublereal* afb, integer* ldafb, integer* ipiv, doublereal* r__,
	doublereal* c__, doublereal* b, integer* ldb, doublereal* x, integer*
	ldx, doublereal* rcond, doublereal* berr, integer* n_err_bnds__,
	doublereal* err_bnds_norm__, doublereal* err_bnds_comp__, integer*
	nparams, doublereal* params, doublereal* work, integer* iwork,
	integer* info);


CLAPACK_API
int dgbsvxx_(char* fact, char* trans, integer* n, integer*
	kl, integer* ku, integer* nrhs, doublereal* ab, integer* ldab,
	doublereal* afb, integer* ldafb, integer* ipiv, char* equed,
	doublereal* r__, doublereal* c__, doublereal* b, integer* ldb,
	doublereal* x, integer* ldx, doublereal* rcond, doublereal* rpvgrw,
	doublereal* berr, integer* n_err_bnds__, doublereal* err_bnds_norm__,
	doublereal* err_bnds_comp__, integer* nparams, doublereal* params,
	doublereal* work, integer* iwork, integer* info);


CLAPACK_API
int dgerfsx_(char* trans, char* equed, integer* n, integer*
	nrhs, doublereal* a, integer* lda, doublereal* af, integer* ldaf,
	integer* ipiv, doublereal* r__, doublereal* c__, doublereal* b,
	integer* ldb, doublereal* x, integer* ldx, doublereal* rcond,
	doublereal* berr, integer* n_err_bnds__, doublereal* err_bnds_norm__,
	doublereal* err_bnds_comp__, integer* nparams, doublereal* params,
	doublereal* work, integer* iwork, integer* info);


CLAPACK_API
int dgesvxx_(char* fact, char* trans, integer* n, integer*
	nrhs, doublereal* a, integer* lda, doublereal* af, integer* ldaf,
	integer* ipiv, char* equed, doublereal* r__, doublereal* c__,
	doublereal* b, integer* ldb, doublereal* x, integer* ldx, doublereal*
	rcond, doublereal* rpvgrw, doublereal* berr, integer* n_err_bnds__,
	doublereal* err_bnds_norm__, doublereal* err_bnds_comp__, integer*
	nparams, doublereal* params, doublereal* work, integer* iwork,
	integer* info);


CLAPACK_API
int dla_gbamv__(integer* trans, integer* m, integer* n,
	integer* kl, integer* ku, doublereal* alpha, doublereal* ab, integer*
	ldab, doublereal* x, integer* incx, doublereal* beta, doublereal* y,
	integer* incy);

doublereal dla_gbrcond__(char* trans, integer* n, integer* kl, integer* ku,
	doublereal* ab, integer* ldab, doublereal* afb, integer* ldafb,
	integer* ipiv, integer* cmode, doublereal* c__, integer* info,
	doublereal* work, integer* iwork, ftnlen trans_len);


CLAPACK_API
int dla_gbrfsx_extended__(integer* prec_type__, integer*
	trans_type__, integer* n, integer* kl, integer* ku, integer* nrhs,
	doublereal* ab, integer* ldab, doublereal* afb, integer* ldafb,
	integer* ipiv, logical* colequ, doublereal* c__, doublereal* b,
	integer* ldb, doublereal* y, integer* ldy, doublereal* berr_out__,
	integer* n_norms__, doublereal* errs_n__, doublereal* errs_c__,
	doublereal* res, doublereal* ayb, doublereal* dy, doublereal*
	y_tail__, doublereal* rcond, integer* ithresh, doublereal* rthresh,
	doublereal* dz_ub__, logical* ignore_cwise__, integer* info);

CLAPACK_API
doublereal dla_gbrpvgrw__(integer* n, integer* kl, integer* ku, integer*
	ncols, doublereal* ab, integer* ldab, doublereal* afb, integer* ldafb);


CLAPACK_API
int dla_geamv__(integer* trans, integer* m, integer* n,
	doublereal* alpha, doublereal* a, integer* lda, doublereal* x,
	integer* incx, doublereal* beta, doublereal* y, integer* incy);

CLAPACK_API
doublereal dla_gercond__(char* trans, integer* n, doublereal* a, integer* lda,
	doublereal* af, integer* ldaf, integer* ipiv, integer* cmode,
	doublereal* c__, integer* info, doublereal* work, integer* iwork,
	ftnlen trans_len);


CLAPACK_API
int dla_gerfsx_extended__(integer* prec_type__, integer*
	trans_type__, integer* n, integer* nrhs, doublereal* a, integer* lda,
	doublereal* af, integer* ldaf, integer* ipiv, logical* colequ,
	doublereal* c__, doublereal* b, integer* ldb, doublereal* y, integer*
	ldy, doublereal* berr_out__, integer* n_norms__, doublereal* errs_n__,
	doublereal* errs_c__, doublereal* res, doublereal* ayb, doublereal*
	dy, doublereal* y_tail__, doublereal* rcond, integer* ithresh,
	doublereal* rthresh, doublereal* dz_ub__, logical* ignore_cwise__,
	integer* info);


CLAPACK_API
int dla_lin_berr__(integer* n, integer* nz, integer* nrhs,
	doublereal* res, doublereal* ayb, doublereal* berr);

CLAPACK_API
doublereal dla_porcond__(char* uplo, integer* n, doublereal* a, integer* lda,
	doublereal* af, integer* ldaf, integer* cmode, doublereal* c__,
	integer* info, doublereal* work, integer* iwork, ftnlen uplo_len);


CLAPACK_API
int dla_porfsx_extended__(integer* prec_type__, char* uplo,
	integer* n, integer* nrhs, doublereal* a, integer* lda, doublereal*
	af, integer* ldaf, logical* colequ, doublereal* c__, doublereal* b,
	integer* ldb, doublereal* y, integer* ldy, doublereal* berr_out__,
	integer* n_norms__, doublereal* errs_n__, doublereal* errs_c__,
	doublereal* res, doublereal* ayb, doublereal* dy, doublereal*
	y_tail__, doublereal* rcond, integer* ithresh, doublereal* rthresh,
	doublereal* dz_ub__, logical* ignore_cwise__, integer* info, ftnlen
	uplo_len);

CLAPACK_API
doublereal dla_porpvgrw__(char* uplo, integer* ncols, doublereal* a, integer*
	lda, doublereal* af, integer* ldaf, doublereal* work, ftnlen uplo_len);

CLAPACK_API
doublereal dla_rpvgrw__(integer* n, integer* ncols, doublereal* a, integer*
	lda, doublereal* af, integer* ldaf);

CLAPACK_API
int dla_syamv__(integer* uplo, integer* n, doublereal* alpha,
	doublereal* a, integer* lda, doublereal* x, integer* incx,
	doublereal* beta, doublereal* y, integer* incy);

CLAPACK_API
doublereal dla_syrcond__(char* uplo, integer* n, doublereal* a, integer* lda,
	doublereal* af, integer* ldaf, integer* ipiv, integer* cmode,
	doublereal* c__, integer* info, doublereal* work, integer* iwork,
	ftnlen uplo_len);

CLAPACK_API
int dla_syrfsx_extended__(integer* prec_type__, char* uplo,
	integer* n, integer* nrhs, doublereal* a, integer* lda, doublereal*
	af, integer* ldaf, integer* ipiv, logical* colequ, doublereal* c__,
	doublereal* b, integer* ldb, doublereal* y, integer* ldy, doublereal*
	berr_out__, integer* n_norms__, doublereal* errs_n__, doublereal*
	errs_c__, doublereal* res, doublereal* ayb, doublereal* dy,
	doublereal* y_tail__, doublereal* rcond, integer* ithresh, doublereal
	* rthresh, doublereal* dz_ub__, logical* ignore_cwise__, integer* info,
	ftnlen uplo_len);

CLAPACK_API
doublereal dla_syrpvgrw__(char* uplo, integer* n, integer* info, doublereal*
	a, integer* lda, doublereal* af, integer* ldaf, integer* ipiv,
	doublereal* work, ftnlen uplo_len);

CLAPACK_API
int dla_wwaddw__(integer* n, doublereal* x, doublereal* y,
	doublereal* w);

CLAPACK_API
int dlarscl2_(integer* m, integer* n, doublereal* d__,
	doublereal* x, integer* ldx);

CLAPACK_API
int dlascl2_(integer* m, integer* n, doublereal* d__,
	doublereal* x, integer* ldx);

CLAPACK_API
int dporfsx_(char* uplo, char* equed, integer* n, integer*
	nrhs, doublereal* a, integer* lda, doublereal* af, integer* ldaf,
	doublereal* s, doublereal* b, integer* ldb, doublereal* x, integer*
	ldx, doublereal* rcond, doublereal* berr, integer* n_err_bnds__,
	doublereal* err_bnds_norm__, doublereal* err_bnds_comp__, integer*
	nparams, doublereal* params, doublereal* work, integer* iwork,
	integer* info);

CLAPACK_API
int dposvxx_(char* fact, char* uplo, integer* n, integer*
	nrhs, doublereal* a, integer* lda, doublereal* af, integer* ldaf,
	char* equed, doublereal* s, doublereal* b, integer* ldb, doublereal*
	x, integer* ldx, doublereal* rcond, doublereal* rpvgrw, doublereal*
	berr, integer* n_err_bnds__, doublereal* err_bnds_norm__, doublereal*
	err_bnds_comp__, integer* nparams, doublereal* params, doublereal*
	work, integer* iwork, integer* info);

CLAPACK_API
int dsyrfsx_(char* uplo, char* equed, integer* n, integer*
	nrhs, doublereal* a, integer* lda, doublereal* af, integer* ldaf,
	integer* ipiv, doublereal* s, doublereal* b, integer* ldb, doublereal
	* x, integer* ldx, doublereal* rcond, doublereal* berr, integer*
	n_err_bnds__, doublereal* err_bnds_norm__, doublereal*
	err_bnds_comp__, integer* nparams, doublereal* params, doublereal*
	work, integer* iwork, integer* info);

CLAPACK_API
int dsysvxx_(char* fact, char* uplo, integer* n, integer*
	nrhs, doublereal* a, integer* lda, doublereal* af, integer* ldaf,
	integer* ipiv, char* equed, doublereal* s, doublereal* b, integer*
	ldb, doublereal* x, integer* ldx, doublereal* rcond, doublereal*
	rpvgrw, doublereal* berr, integer* n_err_bnds__, doublereal*
	err_bnds_norm__, doublereal* err_bnds_comp__, integer* nparams,
	doublereal* params, doublereal* work, integer* iwork, integer* info);
*/

#pragma endregion







