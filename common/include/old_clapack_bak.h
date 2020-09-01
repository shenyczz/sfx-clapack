/* header file for clapack 3.2.1 */

#ifndef __CLAPACK_H
#define __CLAPACK_H

#include "old_f2c.h"


/*************************************************************
* 
*                        LAPACK
* 
**************************************************************/

#pragma region CLA


/* Subroutine */ int chbgv_(char *jobz, char *uplo, integer *n, integer *ka, 
	integer *kb, complex *ab, integer *ldab, complex *bb, integer *ldbb, 
	real *w, complex *z__, integer *ldz, complex *work, real *rwork, 
	integer *info);

/* Subroutine */ int chbgvd_(char *jobz, char *uplo, integer *n, integer *ka, 
	integer *kb, complex *ab, integer *ldab, complex *bb, integer *ldbb, 
	real *w, complex *z__, integer *ldz, complex *work, integer *lwork, 
	real *rwork, integer *lrwork, integer *iwork, integer *liwork, 
	integer *info);

/* Subroutine */ int chbgvx_(char *jobz, char *range, char *uplo, integer *n, 
	integer *ka, integer *kb, complex *ab, integer *ldab, complex *bb, 
	integer *ldbb, complex *q, integer *ldq, real *vl, real *vu, integer *
	il, integer *iu, real *abstol, integer *m, real *w, complex *z__, 
	integer *ldz, complex *work, real *rwork, integer *iwork, integer *
	ifail, integer *info);

/* Subroutine */ int chbtrd_(char *vect, char *uplo, integer *n, integer *kd, 
	complex *ab, integer *ldab, real *d__, real *e, complex *q, integer *
	ldq, complex *work, integer *info);

/* Subroutine */ int checon_(char *uplo, integer *n, complex *a, integer *lda, 
	 integer *ipiv, real *anorm, real *rcond, complex *work, integer *
	info);

/* Subroutine */ int cheequb_(char *uplo, integer *n, complex *a, integer *
	lda, real *s, real *scond, real *amax, complex *work, integer *info);

/* Subroutine */ int cheev_(char *jobz, char *uplo, integer *n, complex *a, 
	integer *lda, real *w, complex *work, integer *lwork, real *rwork, 
	integer *info);

/* Subroutine */ int cheevd_(char *jobz, char *uplo, integer *n, complex *a, 
	integer *lda, real *w, complex *work, integer *lwork, real *rwork, 
	integer *lrwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int cheevr_(char *jobz, char *range, char *uplo, integer *n, 
	complex *a, integer *lda, real *vl, real *vu, integer *il, integer *
	iu, real *abstol, integer *m, real *w, complex *z__, integer *ldz, 
	integer *isuppz, complex *work, integer *lwork, real *rwork, integer *
	lrwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int cheevx_(char *jobz, char *range, char *uplo, integer *n, 
	complex *a, integer *lda, real *vl, real *vu, integer *il, integer *
	iu, real *abstol, integer *m, real *w, complex *z__, integer *ldz, 
	complex *work, integer *lwork, real *rwork, integer *iwork, integer *
	ifail, integer *info);

/* Subroutine */ int chegs2_(integer *itype, char *uplo, integer *n, complex *
	a, integer *lda, complex *b, integer *ldb, integer *info);

/* Subroutine */ int chegst_(integer *itype, char *uplo, integer *n, complex *
	a, integer *lda, complex *b, integer *ldb, integer *info);

/* Subroutine */ int chegv_(integer *itype, char *jobz, char *uplo, integer *
	n, complex *a, integer *lda, complex *b, integer *ldb, real *w, 
	complex *work, integer *lwork, real *rwork, integer *info);

/* Subroutine */ int chegvd_(integer *itype, char *jobz, char *uplo, integer *
	n, complex *a, integer *lda, complex *b, integer *ldb, real *w, 
	complex *work, integer *lwork, real *rwork, integer *lrwork, integer *
	iwork, integer *liwork, integer *info);

/* Subroutine */ int chegvx_(integer *itype, char *jobz, char *range, char *
	uplo, integer *n, complex *a, integer *lda, complex *b, integer *ldb, 
	real *vl, real *vu, integer *il, integer *iu, real *abstol, integer *
	m, real *w, complex *z__, integer *ldz, complex *work, integer *lwork, 
	 real *rwork, integer *iwork, integer *ifail, integer *info);

/* Subroutine */ int cherfs_(char *uplo, integer *n, integer *nrhs, complex *
	a, integer *lda, complex *af, integer *ldaf, integer *ipiv, complex *
	b, integer *ldb, complex *x, integer *ldx, real *ferr, real *berr, 
	complex *work, real *rwork, integer *info);

/* Subroutine */ int cherfsx_(char *uplo, char *equed, integer *n, integer *
	nrhs, complex *a, integer *lda, complex *af, integer *ldaf, integer *
	ipiv, real *s, complex *b, integer *ldb, complex *x, integer *ldx, 
	real *rcond, real *berr, integer *n_err_bnds__, real *err_bnds_norm__, 
	 real *err_bnds_comp__, integer *nparams, real *params, complex *work, 
	 real *rwork, integer *info);

/* Subroutine */ int chesv_(char *uplo, integer *n, integer *nrhs, complex *a, 
	 integer *lda, integer *ipiv, complex *b, integer *ldb, complex *work, 
	 integer *lwork, integer *info);

/* Subroutine */ int chesvx_(char *fact, char *uplo, integer *n, integer *
	nrhs, complex *a, integer *lda, complex *af, integer *ldaf, integer *
	ipiv, complex *b, integer *ldb, complex *x, integer *ldx, real *rcond, 
	 real *ferr, real *berr, complex *work, integer *lwork, real *rwork, 
	integer *info);

/* Subroutine */ int chesvxx_(char *fact, char *uplo, integer *n, integer *
	nrhs, complex *a, integer *lda, complex *af, integer *ldaf, integer *
	ipiv, char *equed, real *s, complex *b, integer *ldb, complex *x, 
	integer *ldx, real *rcond, real *rpvgrw, real *berr, integer *
	n_err_bnds__, real *err_bnds_norm__, real *err_bnds_comp__, integer *
	nparams, real *params, complex *work, real *rwork, integer *info);

/* Subroutine */ int chetd2_(char *uplo, integer *n, complex *a, integer *lda, 
	 real *d__, real *e, complex *tau, integer *info);

/* Subroutine */ int chetf2_(char *uplo, integer *n, complex *a, integer *lda, 
	 integer *ipiv, integer *info);

/* Subroutine */ int chetrd_(char *uplo, integer *n, complex *a, integer *lda, 
	 real *d__, real *e, complex *tau, complex *work, integer *lwork, 
	integer *info);

/* Subroutine */ int chetrf_(char *uplo, integer *n, complex *a, integer *lda, 
	 integer *ipiv, complex *work, integer *lwork, integer *info);

/* Subroutine */ int chetri_(char *uplo, integer *n, complex *a, integer *lda, 
	 integer *ipiv, complex *work, integer *info);

/* Subroutine */ int chetrs_(char *uplo, integer *n, integer *nrhs, complex *
	a, integer *lda, integer *ipiv, complex *b, integer *ldb, integer *
	info);

/* Subroutine */ int chfrk_(char *transr, char *uplo, char *trans, integer *n, 
	 integer *k, real *alpha, complex *a, integer *lda, real *beta, 
	complex *c__);

/* Subroutine */ int chgeqz_(char *job, char *compq, char *compz, integer *n, 
	integer *ilo, integer *ihi, complex *h__, integer *ldh, complex *t, 
	integer *ldt, complex *alpha, complex *beta, complex *q, integer *ldq, 
	 complex *z__, integer *ldz, complex *work, integer *lwork, real *
	rwork, integer *info);

/* Subroutine */ int chpcon_(char *uplo, integer *n, complex *ap, integer *
	ipiv, real *anorm, real *rcond, complex *work, integer *info);

/* Subroutine */ int chpev_(char *jobz, char *uplo, integer *n, complex *ap, 
	real *w, complex *z__, integer *ldz, complex *work, real *rwork, 
	integer *info);

/* Subroutine */ int chpevd_(char *jobz, char *uplo, integer *n, complex *ap, 
	real *w, complex *z__, integer *ldz, complex *work, integer *lwork, 
	real *rwork, integer *lrwork, integer *iwork, integer *liwork, 
	integer *info);

/* Subroutine */ int chpevx_(char *jobz, char *range, char *uplo, integer *n, 
	complex *ap, real *vl, real *vu, integer *il, integer *iu, real *
	abstol, integer *m, real *w, complex *z__, integer *ldz, complex *
	work, real *rwork, integer *iwork, integer *ifail, integer *info);

/* Subroutine */ int chpgst_(integer *itype, char *uplo, integer *n, complex *
	ap, complex *bp, integer *info);

/* Subroutine */ int chpgv_(integer *itype, char *jobz, char *uplo, integer *
	n, complex *ap, complex *bp, real *w, complex *z__, integer *ldz, 
	complex *work, real *rwork, integer *info);

/* Subroutine */ int chpgvd_(integer *itype, char *jobz, char *uplo, integer *
	n, complex *ap, complex *bp, real *w, complex *z__, integer *ldz, 
	complex *work, integer *lwork, real *rwork, integer *lrwork, integer *
	iwork, integer *liwork, integer *info);

/* Subroutine */ int chpgvx_(integer *itype, char *jobz, char *range, char *
	uplo, integer *n, complex *ap, complex *bp, real *vl, real *vu, 
	integer *il, integer *iu, real *abstol, integer *m, real *w, complex *
	z__, integer *ldz, complex *work, real *rwork, integer *iwork, 
	integer *ifail, integer *info);

/* Subroutine */ int chprfs_(char *uplo, integer *n, integer *nrhs, complex *
	ap, complex *afp, integer *ipiv, complex *b, integer *ldb, complex *x, 
	 integer *ldx, real *ferr, real *berr, complex *work, real *rwork, 
	integer *info);

/* Subroutine */ int chpsv_(char *uplo, integer *n, integer *nrhs, complex *
	ap, integer *ipiv, complex *b, integer *ldb, integer *info);

/* Subroutine */ int chpsvx_(char *fact, char *uplo, integer *n, integer *
	nrhs, complex *ap, complex *afp, integer *ipiv, complex *b, integer *
	ldb, complex *x, integer *ldx, real *rcond, real *ferr, real *berr, 
	complex *work, real *rwork, integer *info);

/* Subroutine */ int chptrd_(char *uplo, integer *n, complex *ap, real *d__, 
	real *e, complex *tau, integer *info);

/* Subroutine */ int chptrf_(char *uplo, integer *n, complex *ap, integer *
	ipiv, integer *info);

/* Subroutine */ int chptri_(char *uplo, integer *n, complex *ap, integer *
	ipiv, complex *work, integer *info);

/* Subroutine */ int chptrs_(char *uplo, integer *n, integer *nrhs, complex *
	ap, integer *ipiv, complex *b, integer *ldb, integer *info);

/* Subroutine */ int chsein_(char *side, char *eigsrc, char *initv, logical *
	select, integer *n, complex *h__, integer *ldh, complex *w, complex *
	vl, integer *ldvl, complex *vr, integer *ldvr, integer *mm, integer *
	m, complex *work, real *rwork, integer *ifaill, integer *ifailr, 
	integer *info);

/* Subroutine */ int chseqr_(char *job, char *compz, integer *n, integer *ilo, 
	 integer *ihi, complex *h__, integer *ldh, complex *w, complex *z__, 
	integer *ldz, complex *work, integer *lwork, integer *info);

/* Subroutine */ int cla_gbamv__(integer *trans, integer *m, integer *n, 
	integer *kl, integer *ku, real *alpha, complex *ab, integer *ldab, 
	complex *x, integer *incx, real *beta, real *y, integer *incy);

doublereal cla_gbrcond_c__(char *trans, integer *n, integer *kl, integer *ku, 
	complex *ab, integer *ldab, complex *afb, integer *ldafb, integer *
	ipiv, real *c__, logical *capply, integer *info, complex *work, real *
	rwork, ftnlen trans_len);

doublereal cla_gbrcond_x__(char *trans, integer *n, integer *kl, integer *ku, 
	complex *ab, integer *ldab, complex *afb, integer *ldafb, integer *
	ipiv, complex *x, integer *info, complex *work, real *rwork, ftnlen 
	trans_len);

/* Subroutine */ int cla_gbrfsx_extended__(integer *prec_type__, integer *
	trans_type__, integer *n, integer *kl, integer *ku, integer *nrhs, 
	complex *ab, integer *ldab, complex *afb, integer *ldafb, integer *
	ipiv, logical *colequ, real *c__, complex *b, integer *ldb, complex *
	y, integer *ldy, real *berr_out__, integer *n_norms__, real *errs_n__,
	 real *errs_c__, complex *res, real *ayb, complex *dy, complex *
	y_tail__, real *rcond, integer *ithresh, real *rthresh, real *dz_ub__,
	 logical *ignore_cwise__, integer *info);

doublereal cla_gbrpvgrw__(integer *n, integer *kl, integer *ku, integer *
	ncols, complex *ab, integer *ldab, complex *afb, integer *ldafb);

/* Subroutine */ int cla_geamv__(integer *trans, integer *m, integer *n, real 
	*alpha, complex *a, integer *lda, complex *x, integer *incx, real *
	beta, real *y, integer *incy);

doublereal cla_gercond_c__(char *trans, integer *n, complex *a, integer *lda, 
	complex *af, integer *ldaf, integer *ipiv, real *c__, logical *capply,
	 integer *info, complex *work, real *rwork, ftnlen trans_len);

doublereal cla_gercond_x__(char *trans, integer *n, complex *a, integer *lda, 
	complex *af, integer *ldaf, integer *ipiv, complex *x, integer *info, 
	complex *work, real *rwork, ftnlen trans_len);

/* Subroutine */ int cla_gerfsx_extended__(integer *prec_type__, integer *
	trans_type__, integer *n, integer *nrhs, complex *a, integer *lda, 
	complex *af, integer *ldaf, integer *ipiv, logical *colequ, real *c__,
	 complex *b, integer *ldb, complex *y, integer *ldy, real *berr_out__,
	 integer *n_norms__, real *errs_n__, real *errs_c__, complex *res, 
	real *ayb, complex *dy, complex *y_tail__, real *rcond, integer *
	ithresh, real *rthresh, real *dz_ub__, logical *ignore_cwise__, 
	integer *info);

/* Subroutine */ int cla_heamv__(integer *uplo, integer *n, real *alpha, 
	complex *a, integer *lda, complex *x, integer *incx, real *beta, real 
	*y, integer *incy);

doublereal cla_hercond_c__(char *uplo, integer *n, complex *a, integer *lda, 
	complex *af, integer *ldaf, integer *ipiv, real *c__, logical *capply,
	 integer *info, complex *work, real *rwork, ftnlen uplo_len);

doublereal cla_hercond_x__(char *uplo, integer *n, complex *a, integer *lda, 
	complex *af, integer *ldaf, integer *ipiv, complex *x, integer *info, 
	complex *work, real *rwork, ftnlen uplo_len);

/* Subroutine */ int cla_herfsx_extended__(integer *prec_type__, char *uplo, 
	integer *n, integer *nrhs, complex *a, integer *lda, complex *af, 
	integer *ldaf, integer *ipiv, logical *colequ, real *c__, complex *b, 
	integer *ldb, complex *y, integer *ldy, real *berr_out__, integer *
	n_norms__, real *errs_n__, real *errs_c__, complex *res, real *ayb, 
	complex *dy, complex *y_tail__, real *rcond, integer *ithresh, real *
	rthresh, real *dz_ub__, logical *ignore_cwise__, integer *info, 
	ftnlen uplo_len);

doublereal cla_herpvgrw__(char *uplo, integer *n, integer *info, complex *a, 
	integer *lda, complex *af, integer *ldaf, integer *ipiv, real *work, 
	ftnlen uplo_len);

/* Subroutine */ int cla_lin_berr__(integer *n, integer *nz, integer *nrhs, 
	complex *res, real *ayb, real *berr);

doublereal cla_porcond_c__(char *uplo, integer *n, complex *a, integer *lda, 
	complex *af, integer *ldaf, real *c__, logical *capply, integer *info,
	 complex *work, real *rwork, ftnlen uplo_len);

doublereal cla_porcond_x__(char *uplo, integer *n, complex *a, integer *lda, 
	complex *af, integer *ldaf, complex *x, integer *info, complex *work, 
	real *rwork, ftnlen uplo_len);

/* Subroutine */ int cla_porfsx_extended__(integer *prec_type__, char *uplo, 
	integer *n, integer *nrhs, complex *a, integer *lda, complex *af, 
	integer *ldaf, logical *colequ, real *c__, complex *b, integer *ldb, 
	complex *y, integer *ldy, real *berr_out__, integer *n_norms__, real *
	errs_n__, real *errs_c__, complex *res, real *ayb, complex *dy, 
	complex *y_tail__, real *rcond, integer *ithresh, real *rthresh, real 
	*dz_ub__, logical *ignore_cwise__, integer *info, ftnlen uplo_len);

doublereal cla_porpvgrw__(char *uplo, integer *ncols, complex *a, integer *
	lda, complex *af, integer *ldaf, real *work, ftnlen uplo_len);

doublereal cla_rpvgrw__(integer *n, integer *ncols, complex *a, integer *lda, 
	complex *af, integer *ldaf);

/* Subroutine */ int cla_syamv__(integer *uplo, integer *n, real *alpha, 
	complex *a, integer *lda, complex *x, integer *incx, real *beta, real 
	*y, integer *incy);

doublereal cla_syrcond_c__(char *uplo, integer *n, complex *a, integer *lda, 
	complex *af, integer *ldaf, integer *ipiv, real *c__, logical *capply,
	 integer *info, complex *work, real *rwork, ftnlen uplo_len);

doublereal cla_syrcond_x__(char *uplo, integer *n, complex *a, integer *lda, 
	complex *af, integer *ldaf, integer *ipiv, complex *x, integer *info, 
	complex *work, real *rwork, ftnlen uplo_len);

/* Subroutine */ int cla_syrfsx_extended__(integer *prec_type__, char *uplo, 
	integer *n, integer *nrhs, complex *a, integer *lda, complex *af, 
	integer *ldaf, integer *ipiv, logical *colequ, real *c__, complex *b, 
	integer *ldb, complex *y, integer *ldy, real *berr_out__, integer *
	n_norms__, real *errs_n__, real *errs_c__, complex *res, real *ayb, 
	complex *dy, complex *y_tail__, real *rcond, integer *ithresh, real *
	rthresh, real *dz_ub__, logical *ignore_cwise__, integer *info, 
	ftnlen uplo_len);

doublereal cla_syrpvgrw__(char *uplo, integer *n, integer *info, complex *a, 
	integer *lda, complex *af, integer *ldaf, integer *ipiv, real *work, 
	ftnlen uplo_len);

/* Subroutine */ int cla_wwaddw__(integer *n, complex *x, complex *y, complex 
	*w);

/* Subroutine */ int clabrd_(integer *m, integer *n, integer *nb, complex *a, 
	integer *lda, real *d__, real *e, complex *tauq, complex *taup, 
	complex *x, integer *ldx, complex *y, integer *ldy);

/* Subroutine */ int clacgv_(integer *n, complex *x, integer *incx);

/* Subroutine */ int clacn2_(integer *n, complex *v, complex *x, real *est, 
	integer *kase, integer *isave);

/* Subroutine */ int clacon_(integer *n, complex *v, complex *x, real *est, 
	integer *kase);

/* Subroutine */ int clacp2_(char *uplo, integer *m, integer *n, real *a, 
	integer *lda, complex *b, integer *ldb);

/* Subroutine */ int clacpy_(char *uplo, integer *m, integer *n, complex *a, 
	integer *lda, complex *b, integer *ldb);

/* Subroutine */ int clacrm_(integer *m, integer *n, complex *a, integer *lda, 
	 real *b, integer *ldb, complex *c__, integer *ldc, real *rwork);

/* Subroutine */ int clacrt_(integer *n, complex *cx, integer *incx, complex *
	cy, integer *incy, complex *c__, complex *s);

/* Complex */ VOID cladiv_(complex * ret_val, complex *x, complex *y);

/* Subroutine */ int claed0_(integer *qsiz, integer *n, real *d__, real *e, 
	complex *q, integer *ldq, complex *qstore, integer *ldqs, real *rwork, 
	 integer *iwork, integer *info);

/* Subroutine */ int claed7_(integer *n, integer *cutpnt, integer *qsiz, 
	integer *tlvls, integer *curlvl, integer *curpbm, real *d__, complex *
	q, integer *ldq, real *rho, integer *indxq, real *qstore, integer *
	qptr, integer *prmptr, integer *perm, integer *givptr, integer *
	givcol, real *givnum, complex *work, real *rwork, integer *iwork, 
	integer *info);

/* Subroutine */ int claed8_(integer *k, integer *n, integer *qsiz, complex *
	q, integer *ldq, real *d__, real *rho, integer *cutpnt, real *z__, 
	real *dlamda, complex *q2, integer *ldq2, real *w, integer *indxp, 
	integer *indx, integer *indxq, integer *perm, integer *givptr, 
	integer *givcol, real *givnum, integer *info);

/* Subroutine */ int claein_(logical *rightv, logical *noinit, integer *n, 
	complex *h__, integer *ldh, complex *w, complex *v, complex *b, 
	integer *ldb, real *rwork, real *eps3, real *smlnum, integer *info);

/* Subroutine */ int claesy_(complex *a, complex *b, complex *c__, complex *
	rt1, complex *rt2, complex *evscal, complex *cs1, complex *sn1);

/* Subroutine */ int claev2_(complex *a, complex *b, complex *c__, real *rt1, 
	real *rt2, real *cs1, complex *sn1);

/* Subroutine */ int clag2z_(integer *m, integer *n, complex *sa, integer *
	ldsa, doublecomplex *a, integer *lda, integer *info);

/* Subroutine */ int clags2_(logical *upper, real *a1, complex *a2, real *a3, 
	real *b1, complex *b2, real *b3, real *csu, complex *snu, real *csv, 
	complex *snv, real *csq, complex *snq);

/* Subroutine */ int clagtm_(char *trans, integer *n, integer *nrhs, real *
	alpha, complex *dl, complex *d__, complex *du, complex *x, integer *
	ldx, real *beta, complex *b, integer *ldb);

/* Subroutine */ int clahef_(char *uplo, integer *n, integer *nb, integer *kb, 
	 complex *a, integer *lda, integer *ipiv, complex *w, integer *ldw, 
	integer *info);

/* Subroutine */ int clahqr_(logical *wantt, logical *wantz, integer *n, 
	integer *ilo, integer *ihi, complex *h__, integer *ldh, complex *w, 
	integer *iloz, integer *ihiz, complex *z__, integer *ldz, integer *
	info);

/* Subroutine */ int clahr2_(integer *n, integer *k, integer *nb, complex *a, 
	integer *lda, complex *tau, complex *t, integer *ldt, complex *y, 
	integer *ldy);

/* Subroutine */ int clahrd_(integer *n, integer *k, integer *nb, complex *a, 
	integer *lda, complex *tau, complex *t, integer *ldt, complex *y, 
	integer *ldy);

/* Subroutine */ int claic1_(integer *job, integer *j, complex *x, real *sest, 
	 complex *w, complex *gamma, real *sestpr, complex *s, complex *c__);

/* Subroutine */ int clals0_(integer *icompq, integer *nl, integer *nr, 
	integer *sqre, integer *nrhs, complex *b, integer *ldb, complex *bx, 
	integer *ldbx, integer *perm, integer *givptr, integer *givcol, 
	integer *ldgcol, real *givnum, integer *ldgnum, real *poles, real *
	difl, real *difr, real *z__, integer *k, real *c__, real *s, real *
	rwork, integer *info);

/* Subroutine */ int clalsa_(integer *icompq, integer *smlsiz, integer *n, 
	integer *nrhs, complex *b, integer *ldb, complex *bx, integer *ldbx, 
	real *u, integer *ldu, real *vt, integer *k, real *difl, real *difr, 
	real *z__, real *poles, integer *givptr, integer *givcol, integer *
	ldgcol, integer *perm, real *givnum, real *c__, real *s, real *rwork, 
	integer *iwork, integer *info);

/* Subroutine */ int clalsd_(char *uplo, integer *smlsiz, integer *n, integer 
	*nrhs, real *d__, real *e, complex *b, integer *ldb, real *rcond, 
	integer *rank, complex *work, real *rwork, integer *iwork, integer *
	info);

doublereal clangb_(char *norm, integer *n, integer *kl, integer *ku, complex *
	ab, integer *ldab, real *work);

doublereal clange_(char *norm, integer *m, integer *n, complex *a, integer *
	lda, real *work);

doublereal clangt_(char *norm, integer *n, complex *dl, complex *d__, complex 
	*du);

doublereal clanhb_(char *norm, char *uplo, integer *n, integer *k, complex *
	ab, integer *ldab, real *work);

doublereal clanhe_(char *norm, char *uplo, integer *n, complex *a, integer *
	lda, real *work);

doublereal clanhf_(char *norm, char *transr, char *uplo, integer *n, complex *
	a, real *work);

doublereal clanhp_(char *norm, char *uplo, integer *n, complex *ap, real *
	work);

doublereal clanhs_(char *norm, integer *n, complex *a, integer *lda, real *
	work);

doublereal clanht_(char *norm, integer *n, real *d__, complex *e);

doublereal clansb_(char *norm, char *uplo, integer *n, integer *k, complex *
	ab, integer *ldab, real *work);

doublereal clansp_(char *norm, char *uplo, integer *n, complex *ap, real *
	work);

doublereal clansy_(char *norm, char *uplo, integer *n, complex *a, integer *
	lda, real *work);

doublereal clantb_(char *norm, char *uplo, char *diag, integer *n, integer *k, 
	 complex *ab, integer *ldab, real *work);

doublereal clantp_(char *norm, char *uplo, char *diag, integer *n, complex *
	ap, real *work);

doublereal clantr_(char *norm, char *uplo, char *diag, integer *m, integer *n, 
	 complex *a, integer *lda, real *work);

/* Subroutine */ int clapll_(integer *n, complex *x, integer *incx, complex *
	y, integer *incy, real *ssmin);

/* Subroutine */ int clapmt_(logical *forwrd, integer *m, integer *n, complex 
	*x, integer *ldx, integer *k);

/* Subroutine */ int claqgb_(integer *m, integer *n, integer *kl, integer *ku, 
	 complex *ab, integer *ldab, real *r__, real *c__, real *rowcnd, real 
	*colcnd, real *amax, char *equed);

/* Subroutine */ int claqge_(integer *m, integer *n, complex *a, integer *lda, 
	 real *r__, real *c__, real *rowcnd, real *colcnd, real *amax, char *
	equed);

/* Subroutine */ int claqhb_(char *uplo, integer *n, integer *kd, complex *ab, 
	 integer *ldab, real *s, real *scond, real *amax, char *equed);

/* Subroutine */ int claqhe_(char *uplo, integer *n, complex *a, integer *lda, 
	 real *s, real *scond, real *amax, char *equed);

/* Subroutine */ int claqhp_(char *uplo, integer *n, complex *ap, real *s, 
	real *scond, real *amax, char *equed);

/* Subroutine */ int claqp2_(integer *m, integer *n, integer *offset, complex 
	*a, integer *lda, integer *jpvt, complex *tau, real *vn1, real *vn2, 
	complex *work);

/* Subroutine */ int claqps_(integer *m, integer *n, integer *offset, integer 
	*nb, integer *kb, complex *a, integer *lda, integer *jpvt, complex *
	tau, real *vn1, real *vn2, complex *auxv, complex *f, integer *ldf);

/* Subroutine */ int claqr0_(logical *wantt, logical *wantz, integer *n, 
	integer *ilo, integer *ihi, complex *h__, integer *ldh, complex *w, 
	integer *iloz, integer *ihiz, complex *z__, integer *ldz, complex *
	work, integer *lwork, integer *info);

/* Subroutine */ int claqr1_(integer *n, complex *h__, integer *ldh, complex *
	s1, complex *s2, complex *v);

/* Subroutine */ int claqr2_(logical *wantt, logical *wantz, integer *n, 
	integer *ktop, integer *kbot, integer *nw, complex *h__, integer *ldh, 
	 integer *iloz, integer *ihiz, complex *z__, integer *ldz, integer *
	ns, integer *nd, complex *sh, complex *v, integer *ldv, integer *nh, 
	complex *t, integer *ldt, integer *nv, complex *wv, integer *ldwv, 
	complex *work, integer *lwork);

/* Subroutine */ int claqr3_(logical *wantt, logical *wantz, integer *n, 
	integer *ktop, integer *kbot, integer *nw, complex *h__, integer *ldh, 
	 integer *iloz, integer *ihiz, complex *z__, integer *ldz, integer *
	ns, integer *nd, complex *sh, complex *v, integer *ldv, integer *nh, 
	complex *t, integer *ldt, integer *nv, complex *wv, integer *ldwv, 
	complex *work, integer *lwork);

/* Subroutine */ int claqr4_(logical *wantt, logical *wantz, integer *n, 
	integer *ilo, integer *ihi, complex *h__, integer *ldh, complex *w, 
	integer *iloz, integer *ihiz, complex *z__, integer *ldz, complex *
	work, integer *lwork, integer *info);

/* Subroutine */ int claqr5_(logical *wantt, logical *wantz, integer *kacc22, 
	integer *n, integer *ktop, integer *kbot, integer *nshfts, complex *s, 
	 complex *h__, integer *ldh, integer *iloz, integer *ihiz, complex *
	z__, integer *ldz, complex *v, integer *ldv, complex *u, integer *ldu, 
	 integer *nv, complex *wv, integer *ldwv, integer *nh, complex *wh, 
	integer *ldwh);

/* Subroutine */ int claqsb_(char *uplo, integer *n, integer *kd, complex *ab, 
	 integer *ldab, real *s, real *scond, real *amax, char *equed);

/* Subroutine */ int claqsp_(char *uplo, integer *n, complex *ap, real *s, 
	real *scond, real *amax, char *equed);

/* Subroutine */ int claqsy_(char *uplo, integer *n, complex *a, integer *lda, 
	 real *s, real *scond, real *amax, char *equed);

/* Subroutine */ int clar1v_(integer *n, integer *b1, integer *bn, real *
	lambda, real *d__, real *l, real *ld, real *lld, real *pivmin, real *
	gaptol, complex *z__, logical *wantnc, integer *negcnt, real *ztz, 
	real *mingma, integer *r__, integer *isuppz, real *nrminv, real *
	resid, real *rqcorr, real *work);

/* Subroutine */ int clar2v_(integer *n, complex *x, complex *y, complex *z__, 
	 integer *incx, real *c__, complex *s, integer *incc);

/* Subroutine */ int clarcm_(integer *m, integer *n, real *a, integer *lda, 
	complex *b, integer *ldb, complex *c__, integer *ldc, real *rwork);

/* Subroutine */ int clarf_(char *side, integer *m, integer *n, complex *v, 
	integer *incv, complex *tau, complex *c__, integer *ldc, complex *
	work);

/* Subroutine */ int clarfb_(char *side, char *trans, char *direct, char *
	storev, integer *m, integer *n, integer *k, complex *v, integer *ldv, 
	complex *t, integer *ldt, complex *c__, integer *ldc, complex *work, 
	integer *ldwork);

/* Subroutine */ int clarfg_(integer *n, complex *alpha, complex *x, integer *
	incx, complex *tau);

/* Subroutine */ int clarfp_(integer *n, complex *alpha, complex *x, integer *
	incx, complex *tau);

/* Subroutine */ int clarft_(char *direct, char *storev, integer *n, integer *
	k, complex *v, integer *ldv, complex *tau, complex *t, integer *ldt);

/* Subroutine */ int clarfx_(char *side, integer *m, integer *n, complex *v, 
	complex *tau, complex *c__, integer *ldc, complex *work);

/* Subroutine */ int clargv_(integer *n, complex *x, integer *incx, complex *
	y, integer *incy, real *c__, integer *incc);

/* Subroutine */ int clarnv_(integer *idist, integer *iseed, integer *n, 
	complex *x);

/* Subroutine */ int clarrv_(integer *n, real *vl, real *vu, real *d__, real *
	l, real *pivmin, integer *isplit, integer *m, integer *dol, integer *
	dou, real *minrgp, real *rtol1, real *rtol2, real *w, real *werr, 
	real *wgap, integer *iblock, integer *indexw, real *gers, complex *
	z__, integer *ldz, integer *isuppz, real *work, integer *iwork, 
	integer *info);

/* Subroutine */ int clarscl2_(integer *m, integer *n, real *d__, complex *x, 
	integer *ldx);

/* Subroutine */ int clartg_(complex *f, complex *g, real *cs, complex *sn, 
	complex *r__);

/* Subroutine */ int clartv_(integer *n, complex *x, integer *incx, complex *
	y, integer *incy, real *c__, complex *s, integer *incc);

/* Subroutine */ int clarz_(char *side, integer *m, integer *n, integer *l, 
	complex *v, integer *incv, complex *tau, complex *c__, integer *ldc, 
	complex *work);

/* Subroutine */ int clarzb_(char *side, char *trans, char *direct, char *
	storev, integer *m, integer *n, integer *k, integer *l, complex *v, 
	integer *ldv, complex *t, integer *ldt, complex *c__, integer *ldc, 
	complex *work, integer *ldwork);

/* Subroutine */ int clarzt_(char *direct, char *storev, integer *n, integer *
	k, complex *v, integer *ldv, complex *tau, complex *t, integer *ldt);

/* Subroutine */ int clascl_(char *type__, integer *kl, integer *ku, real *
	cfrom, real *cto, integer *m, integer *n, complex *a, integer *lda, 
	integer *info);

/* Subroutine */ int clascl2_(integer *m, integer *n, real *d__, complex *x, 
	integer *ldx);

/* Subroutine */ int claset_(char *uplo, integer *m, integer *n, complex *
	alpha, complex *beta, complex *a, integer *lda);

/* Subroutine */ int clasr_(char *side, char *pivot, char *direct, integer *m, 
	 integer *n, real *c__, real *s, complex *a, integer *lda);

/* Subroutine */ int classq_(integer *n, complex *x, integer *incx, real *
	scale, real *sumsq);

/* Subroutine */ int claswp_(integer *n, complex *a, integer *lda, integer *
	k1, integer *k2, integer *ipiv, integer *incx);

/* Subroutine */ int clasyf_(char *uplo, integer *n, integer *nb, integer *kb, 
	 complex *a, integer *lda, integer *ipiv, complex *w, integer *ldw, 
	integer *info);

/* Subroutine */ int clatbs_(char *uplo, char *trans, char *diag, char *
	normin, integer *n, integer *kd, complex *ab, integer *ldab, complex *
	x, real *scale, real *cnorm, integer *info);

/* Subroutine */ int clatdf_(integer *ijob, integer *n, complex *z__, integer 
	*ldz, complex *rhs, real *rdsum, real *rdscal, integer *ipiv, integer 
	*jpiv);

/* Subroutine */ int clatps_(char *uplo, char *trans, char *diag, char *
	normin, integer *n, complex *ap, complex *x, real *scale, real *cnorm, 
	 integer *info);

/* Subroutine */ int clatrd_(char *uplo, integer *n, integer *nb, complex *a, 
	integer *lda, real *e, complex *tau, complex *w, integer *ldw);

/* Subroutine */ int clatrs_(char *uplo, char *trans, char *diag, char *
	normin, integer *n, complex *a, integer *lda, complex *x, real *scale, 
	 real *cnorm, integer *info);

/* Subroutine */ int clatrz_(integer *m, integer *n, integer *l, complex *a, 
	integer *lda, complex *tau, complex *work);

/* Subroutine */ int clatzm_(char *side, integer *m, integer *n, complex *v, 
	integer *incv, complex *tau, complex *c1, complex *c2, integer *ldc, 
	complex *work);

/* Subroutine */ int clauu2_(char *uplo, integer *n, complex *a, integer *lda, 
	 integer *info);

/* Subroutine */ int clauum_(char *uplo, integer *n, complex *a, integer *lda, 
	 integer *info);

/* Subroutine */ int cpbcon_(char *uplo, integer *n, integer *kd, complex *ab, 
	 integer *ldab, real *anorm, real *rcond, complex *work, real *rwork, 
	integer *info);

/* Subroutine */ int cpbequ_(char *uplo, integer *n, integer *kd, complex *ab, 
	 integer *ldab, real *s, real *scond, real *amax, integer *info);

/* Subroutine */ int cpbrfs_(char *uplo, integer *n, integer *kd, integer *
	nrhs, complex *ab, integer *ldab, complex *afb, integer *ldafb, 
	complex *b, integer *ldb, complex *x, integer *ldx, real *ferr, real *
	berr, complex *work, real *rwork, integer *info);

/* Subroutine */ int cpbstf_(char *uplo, integer *n, integer *kd, complex *ab, 
	 integer *ldab, integer *info);

/* Subroutine */ int cpbsv_(char *uplo, integer *n, integer *kd, integer *
	nrhs, complex *ab, integer *ldab, complex *b, integer *ldb, integer *
	info);

/* Subroutine */ int cpbsvx_(char *fact, char *uplo, integer *n, integer *kd, 
	integer *nrhs, complex *ab, integer *ldab, complex *afb, integer *
	ldafb, char *equed, real *s, complex *b, integer *ldb, complex *x, 
	integer *ldx, real *rcond, real *ferr, real *berr, complex *work, 
	real *rwork, integer *info);

/* Subroutine */ int cpbtf2_(char *uplo, integer *n, integer *kd, complex *ab, 
	 integer *ldab, integer *info);

/* Subroutine */ int cpbtrf_(char *uplo, integer *n, integer *kd, complex *ab, 
	 integer *ldab, integer *info);

/* Subroutine */ int cpbtrs_(char *uplo, integer *n, integer *kd, integer *
	nrhs, complex *ab, integer *ldab, complex *b, integer *ldb, integer *
	info);

/* Subroutine */ int cpftrf_(char *transr, char *uplo, integer *n, complex *a, 
	 integer *info);

/* Subroutine */ int cpftri_(char *transr, char *uplo, integer *n, complex *a, 
	 integer *info);

/* Subroutine */ int cpftrs_(char *transr, char *uplo, integer *n, integer *
	nrhs, complex *a, complex *b, integer *ldb, integer *info);

/* Subroutine */ int cpocon_(char *uplo, integer *n, complex *a, integer *lda, 
	 real *anorm, real *rcond, complex *work, real *rwork, integer *info);

/* Subroutine */ int cpoequ_(integer *n, complex *a, integer *lda, real *s, 
	real *scond, real *amax, integer *info);

/* Subroutine */ int cpoequb_(integer *n, complex *a, integer *lda, real *s, 
	real *scond, real *amax, integer *info);

/* Subroutine */ int cporfs_(char *uplo, integer *n, integer *nrhs, complex *
	a, integer *lda, complex *af, integer *ldaf, complex *b, integer *ldb, 
	 complex *x, integer *ldx, real *ferr, real *berr, complex *work, 
	real *rwork, integer *info);

/* Subroutine */ int cporfsx_(char *uplo, char *equed, integer *n, integer *
	nrhs, complex *a, integer *lda, complex *af, integer *ldaf, real *s, 
	complex *b, integer *ldb, complex *x, integer *ldx, real *rcond, real 
	*berr, integer *n_err_bnds__, real *err_bnds_norm__, real *
	err_bnds_comp__, integer *nparams, real *params, complex *work, real *
	rwork, integer *info);

/* Subroutine */ int cposv_(char *uplo, integer *n, integer *nrhs, complex *a, 
	 integer *lda, complex *b, integer *ldb, integer *info);

/* Subroutine */ int cposvx_(char *fact, char *uplo, integer *n, integer *
	nrhs, complex *a, integer *lda, complex *af, integer *ldaf, char *
	equed, real *s, complex *b, integer *ldb, complex *x, integer *ldx, 
	real *rcond, real *ferr, real *berr, complex *work, real *rwork, 
	integer *info);

/* Subroutine */ int cposvxx_(char *fact, char *uplo, integer *n, integer *
	nrhs, complex *a, integer *lda, complex *af, integer *ldaf, char *
	equed, real *s, complex *b, integer *ldb, complex *x, integer *ldx, 
	real *rcond, real *rpvgrw, real *berr, integer *n_err_bnds__, real *
	err_bnds_norm__, real *err_bnds_comp__, integer *nparams, real *
	params, complex *work, real *rwork, integer *info);

/* Subroutine */ int cpotf2_(char *uplo, integer *n, complex *a, integer *lda, 
	 integer *info);

/* Subroutine */ int cpotrf_(char *uplo, integer *n, complex *a, integer *lda, 
	 integer *info);

/* Subroutine */ int cpotri_(char *uplo, integer *n, complex *a, integer *lda, 
	 integer *info);

/* Subroutine */ int cpotrs_(char *uplo, integer *n, integer *nrhs, complex *
	a, integer *lda, complex *b, integer *ldb, integer *info);

/* Subroutine */ int cppcon_(char *uplo, integer *n, complex *ap, real *anorm, 
	 real *rcond, complex *work, real *rwork, integer *info);

/* Subroutine */ int cppequ_(char *uplo, integer *n, complex *ap, real *s, 
	real *scond, real *amax, integer *info);

/* Subroutine */ int cpprfs_(char *uplo, integer *n, integer *nrhs, complex *
	ap, complex *afp, complex *b, integer *ldb, complex *x, integer *ldx, 
	real *ferr, real *berr, complex *work, real *rwork, integer *info);

/* Subroutine */ int cppsv_(char *uplo, integer *n, integer *nrhs, complex *
	ap, complex *b, integer *ldb, integer *info);

/* Subroutine */ int cppsvx_(char *fact, char *uplo, integer *n, integer *
	nrhs, complex *ap, complex *afp, char *equed, real *s, complex *b, 
	integer *ldb, complex *x, integer *ldx, real *rcond, real *ferr, real 
	*berr, complex *work, real *rwork, integer *info);

/* Subroutine */ int cpptrf_(char *uplo, integer *n, complex *ap, integer *
	info);

/* Subroutine */ int cpptri_(char *uplo, integer *n, complex *ap, integer *
	info);

/* Subroutine */ int cpptrs_(char *uplo, integer *n, integer *nrhs, complex *
	ap, complex *b, integer *ldb, integer *info);

/* Subroutine */ int cpstf2_(char *uplo, integer *n, complex *a, integer *lda, 
	 integer *piv, integer *rank, real *tol, real *work, integer *info);

/* Subroutine */ int cpstrf_(char *uplo, integer *n, complex *a, integer *lda, 
	 integer *piv, integer *rank, real *tol, real *work, integer *info);

/* Subroutine */ int cptcon_(integer *n, real *d__, complex *e, real *anorm, 
	real *rcond, real *rwork, integer *info);

/* Subroutine */ int cpteqr_(char *compz, integer *n, real *d__, real *e, 
	complex *z__, integer *ldz, real *work, integer *info);

/* Subroutine */ int cptrfs_(char *uplo, integer *n, integer *nrhs, real *d__, 
	 complex *e, real *df, complex *ef, complex *b, integer *ldb, complex 
	*x, integer *ldx, real *ferr, real *berr, complex *work, real *rwork, 
	integer *info);

/* Subroutine */ int cptsv_(integer *n, integer *nrhs, real *d__, complex *e, 
	complex *b, integer *ldb, integer *info);

/* Subroutine */ int cptsvx_(char *fact, integer *n, integer *nrhs, real *d__, 
	 complex *e, real *df, complex *ef, complex *b, integer *ldb, complex 
	*x, integer *ldx, real *rcond, real *ferr, real *berr, complex *work, 
	real *rwork, integer *info);

/* Subroutine */ int cpttrf_(integer *n, real *d__, complex *e, integer *info);

/* Subroutine */ int cpttrs_(char *uplo, integer *n, integer *nrhs, real *d__, 
	 complex *e, complex *b, integer *ldb, integer *info);

/* Subroutine */ int cptts2_(integer *iuplo, integer *n, integer *nrhs, real *
	d__, complex *e, complex *b, integer *ldb);

/* Subroutine */ int crot_(integer *n, complex *cx, integer *incx, complex *
	cy, integer *incy, real *c__, complex *s);

/* Subroutine */ int cspcon_(char *uplo, integer *n, complex *ap, integer *
	ipiv, real *anorm, real *rcond, complex *work, integer *info);

/* Subroutine */ int cspmv_(char *uplo, integer *n, complex *alpha, complex *
	ap, complex *x, integer *incx, complex *beta, complex *y, integer *
	incy);

/* Subroutine */ int cspr_(char *uplo, integer *n, complex *alpha, complex *x, 
	 integer *incx, complex *ap);

/* Subroutine */ int csprfs_(char *uplo, integer *n, integer *nrhs, complex *
	ap, complex *afp, integer *ipiv, complex *b, integer *ldb, complex *x, 
	 integer *ldx, real *ferr, real *berr, complex *work, real *rwork, 
	integer *info);

/* Subroutine */ int cspsv_(char *uplo, integer *n, integer *nrhs, complex *
	ap, integer *ipiv, complex *b, integer *ldb, integer *info);

/* Subroutine */ int cspsvx_(char *fact, char *uplo, integer *n, integer *
	nrhs, complex *ap, complex *afp, integer *ipiv, complex *b, integer *
	ldb, complex *x, integer *ldx, real *rcond, real *ferr, real *berr, 
	complex *work, real *rwork, integer *info);

/* Subroutine */ int csptrf_(char *uplo, integer *n, complex *ap, integer *
	ipiv, integer *info);

/* Subroutine */ int csptri_(char *uplo, integer *n, complex *ap, integer *
	ipiv, complex *work, integer *info);

/* Subroutine */ int csptrs_(char *uplo, integer *n, integer *nrhs, complex *
	ap, integer *ipiv, complex *b, integer *ldb, integer *info);

/* Subroutine */ int csrscl_(integer *n, real *sa, complex *sx, integer *incx);

/* Subroutine */ int cstedc_(char *compz, integer *n, real *d__, real *e, 
	complex *z__, integer *ldz, complex *work, integer *lwork, real *
	rwork, integer *lrwork, integer *iwork, integer *liwork, integer *
	info);

/* Subroutine */ int cstegr_(char *jobz, char *range, integer *n, real *d__, 
	real *e, real *vl, real *vu, integer *il, integer *iu, real *abstol, 
	integer *m, real *w, complex *z__, integer *ldz, integer *isuppz, 
	real *work, integer *lwork, integer *iwork, integer *liwork, integer *
	info);

/* Subroutine */ int cstein_(integer *n, real *d__, real *e, integer *m, real 
	*w, integer *iblock, integer *isplit, complex *z__, integer *ldz, 
	real *work, integer *iwork, integer *ifail, integer *info);

/* Subroutine */ int cstemr_(char *jobz, char *range, integer *n, real *d__, 
	real *e, real *vl, real *vu, integer *il, integer *iu, integer *m, 
	real *w, complex *z__, integer *ldz, integer *nzc, integer *isuppz, 
	logical *tryrac, real *work, integer *lwork, integer *iwork, integer *
	liwork, integer *info);

/* Subroutine */ int csteqr_(char *compz, integer *n, real *d__, real *e, 
	complex *z__, integer *ldz, real *work, integer *info);

/* Subroutine */ int csycon_(char *uplo, integer *n, complex *a, integer *lda, 
	 integer *ipiv, real *anorm, real *rcond, complex *work, integer *
	info);

/* Subroutine */ int csyequb_(char *uplo, integer *n, complex *a, integer *
	lda, real *s, real *scond, real *amax, complex *work, integer *info);

/* Subroutine */ int csymv_(char *uplo, integer *n, complex *alpha, complex *
	a, integer *lda, complex *x, integer *incx, complex *beta, complex *y, 
	 integer *incy);

/* Subroutine */ int csyr_(char *uplo, integer *n, complex *alpha, complex *x, 
	 integer *incx, complex *a, integer *lda);

/* Subroutine */ int csyrfs_(char *uplo, integer *n, integer *nrhs, complex *
	a, integer *lda, complex *af, integer *ldaf, integer *ipiv, complex *
	b, integer *ldb, complex *x, integer *ldx, real *ferr, real *berr, 
	complex *work, real *rwork, integer *info);

/* Subroutine */ int csyrfsx_(char *uplo, char *equed, integer *n, integer *
	nrhs, complex *a, integer *lda, complex *af, integer *ldaf, integer *
	ipiv, real *s, complex *b, integer *ldb, complex *x, integer *ldx, 
	real *rcond, real *berr, integer *n_err_bnds__, real *err_bnds_norm__, 
	 real *err_bnds_comp__, integer *nparams, real *params, complex *work, 
	 real *rwork, integer *info);

/* Subroutine */ int csysv_(char *uplo, integer *n, integer *nrhs, complex *a, 
	 integer *lda, integer *ipiv, complex *b, integer *ldb, complex *work, 
	 integer *lwork, integer *info);

/* Subroutine */ int csysvx_(char *fact, char *uplo, integer *n, integer *
	nrhs, complex *a, integer *lda, complex *af, integer *ldaf, integer *
	ipiv, complex *b, integer *ldb, complex *x, integer *ldx, real *rcond, 
	 real *ferr, real *berr, complex *work, integer *lwork, real *rwork, 
	integer *info);

/* Subroutine */ int csysvxx_(char *fact, char *uplo, integer *n, integer *
	nrhs, complex *a, integer *lda, complex *af, integer *ldaf, integer *
	ipiv, char *equed, real *s, complex *b, integer *ldb, complex *x, 
	integer *ldx, real *rcond, real *rpvgrw, real *berr, integer *
	n_err_bnds__, real *err_bnds_norm__, real *err_bnds_comp__, integer *
	nparams, real *params, complex *work, real *rwork, integer *info);

/* Subroutine */ int csytf2_(char *uplo, integer *n, complex *a, integer *lda, 
	 integer *ipiv, integer *info);

/* Subroutine */ int csytrf_(char *uplo, integer *n, complex *a, integer *lda, 
	 integer *ipiv, complex *work, integer *lwork, integer *info);

/* Subroutine */ int csytri_(char *uplo, integer *n, complex *a, integer *lda, 
	 integer *ipiv, complex *work, integer *info);

/* Subroutine */ int csytrs_(char *uplo, integer *n, integer *nrhs, complex *
	a, integer *lda, integer *ipiv, complex *b, integer *ldb, integer *
	info);

/* Subroutine */ int ctbcon_(char *norm, char *uplo, char *diag, integer *n, 
	integer *kd, complex *ab, integer *ldab, real *rcond, complex *work, 
	real *rwork, integer *info);

/* Subroutine */ int ctbrfs_(char *uplo, char *trans, char *diag, integer *n, 
	integer *kd, integer *nrhs, complex *ab, integer *ldab, complex *b, 
	integer *ldb, complex *x, integer *ldx, real *ferr, real *berr, 
	complex *work, real *rwork, integer *info);

/* Subroutine */ int ctbtrs_(char *uplo, char *trans, char *diag, integer *n, 
	integer *kd, integer *nrhs, complex *ab, integer *ldab, complex *b, 
	integer *ldb, integer *info);

/* Subroutine */ int ctfsm_(char *transr, char *side, char *uplo, char *trans, 
	 char *diag, integer *m, integer *n, complex *alpha, complex *a, 
	complex *b, integer *ldb);

/* Subroutine */ int ctftri_(char *transr, char *uplo, char *diag, integer *n, 
	 complex *a, integer *info);

/* Subroutine */ int ctfttp_(char *transr, char *uplo, integer *n, complex *
	arf, complex *ap, integer *info);

/* Subroutine */ int ctfttr_(char *transr, char *uplo, integer *n, complex *
	arf, complex *a, integer *lda, integer *info);

/* Subroutine */ int ctgevc_(char *side, char *howmny, logical *select, 
	integer *n, complex *s, integer *lds, complex *p, integer *ldp, 
	complex *vl, integer *ldvl, complex *vr, integer *ldvr, integer *mm, 
	integer *m, complex *work, real *rwork, integer *info);

/* Subroutine */ int ctgex2_(logical *wantq, logical *wantz, integer *n, 
	complex *a, integer *lda, complex *b, integer *ldb, complex *q, 
	integer *ldq, complex *z__, integer *ldz, integer *j1, integer *info);

/* Subroutine */ int ctgexc_(logical *wantq, logical *wantz, integer *n, 
	complex *a, integer *lda, complex *b, integer *ldb, complex *q, 
	integer *ldq, complex *z__, integer *ldz, integer *ifst, integer *
	ilst, integer *info);

/* Subroutine */ int ctgsen_(integer *ijob, logical *wantq, logical *wantz, 
	logical *select, integer *n, complex *a, integer *lda, complex *b, 
	integer *ldb, complex *alpha, complex *beta, complex *q, integer *ldq, 
	 complex *z__, integer *ldz, integer *m, real *pl, real *pr, real *
	dif, complex *work, integer *lwork, integer *iwork, integer *liwork, 
	integer *info);

/* Subroutine */ int ctgsja_(char *jobu, char *jobv, char *jobq, integer *m, 
	integer *p, integer *n, integer *k, integer *l, complex *a, integer *
	lda, complex *b, integer *ldb, real *tola, real *tolb, real *alpha, 
	real *beta, complex *u, integer *ldu, complex *v, integer *ldv, 
	complex *q, integer *ldq, complex *work, integer *ncycle, integer *
	info);

/* Subroutine */ int ctgsna_(char *job, char *howmny, logical *select, 
	integer *n, complex *a, integer *lda, complex *b, integer *ldb, 
	complex *vl, integer *ldvl, complex *vr, integer *ldvr, real *s, real 
	*dif, integer *mm, integer *m, complex *work, integer *lwork, integer 
	*iwork, integer *info);

/* Subroutine */ int ctgsy2_(char *trans, integer *ijob, integer *m, integer *
	n, complex *a, integer *lda, complex *b, integer *ldb, complex *c__, 
	integer *ldc, complex *d__, integer *ldd, complex *e, integer *lde, 
	complex *f, integer *ldf, real *scale, real *rdsum, real *rdscal, 
	integer *info);

/* Subroutine */ int ctgsyl_(char *trans, integer *ijob, integer *m, integer *
	n, complex *a, integer *lda, complex *b, integer *ldb, complex *c__, 
	integer *ldc, complex *d__, integer *ldd, complex *e, integer *lde, 
	complex *f, integer *ldf, real *scale, real *dif, complex *work, 
	integer *lwork, integer *iwork, integer *info);

/* Subroutine */ int ctpcon_(char *norm, char *uplo, char *diag, integer *n, 
	complex *ap, real *rcond, complex *work, real *rwork, integer *info);

/* Subroutine */ int ctprfs_(char *uplo, char *trans, char *diag, integer *n, 
	integer *nrhs, complex *ap, complex *b, integer *ldb, complex *x, 
	integer *ldx, real *ferr, real *berr, complex *work, real *rwork, 
	integer *info);

/* Subroutine */ int ctptri_(char *uplo, char *diag, integer *n, complex *ap, 
	integer *info);

/* Subroutine */ int ctptrs_(char *uplo, char *trans, char *diag, integer *n, 
	integer *nrhs, complex *ap, complex *b, integer *ldb, integer *info);

/* Subroutine */ int ctpttf_(char *transr, char *uplo, integer *n, complex *
	ap, complex *arf, integer *info);

/* Subroutine */ int ctpttr_(char *uplo, integer *n, complex *ap, complex *a, 
	integer *lda, integer *info);

/* Subroutine */ int ctrcon_(char *norm, char *uplo, char *diag, integer *n, 
	complex *a, integer *lda, real *rcond, complex *work, real *rwork, 
	integer *info);

/* Subroutine */ int ctrevc_(char *side, char *howmny, logical *select, 
	integer *n, complex *t, integer *ldt, complex *vl, integer *ldvl, 
	complex *vr, integer *ldvr, integer *mm, integer *m, complex *work, 
	real *rwork, integer *info);

/* Subroutine */ int ctrexc_(char *compq, integer *n, complex *t, integer *
	ldt, complex *q, integer *ldq, integer *ifst, integer *ilst, integer *
	info);

/* Subroutine */ int ctrrfs_(char *uplo, char *trans, char *diag, integer *n, 
	integer *nrhs, complex *a, integer *lda, complex *b, integer *ldb, 
	complex *x, integer *ldx, real *ferr, real *berr, complex *work, real 
	*rwork, integer *info);

/* Subroutine */ int ctrsen_(char *job, char *compq, logical *select, integer 
	*n, complex *t, integer *ldt, complex *q, integer *ldq, complex *w, 
	integer *m, real *s, real *sep, complex *work, integer *lwork, 
	integer *info);

/* Subroutine */ int ctrsna_(char *job, char *howmny, logical *select, 
	integer *n, complex *t, integer *ldt, complex *vl, integer *ldvl, 
	complex *vr, integer *ldvr, real *s, real *sep, integer *mm, integer *
	m, complex *work, integer *ldwork, real *rwork, integer *info);

/* Subroutine */ int ctrsyl_(char *trana, char *tranb, integer *isgn, integer 
	*m, integer *n, complex *a, integer *lda, complex *b, integer *ldb, 
	complex *c__, integer *ldc, real *scale, integer *info);

/* Subroutine */ int ctrti2_(char *uplo, char *diag, integer *n, complex *a, 
	integer *lda, integer *info);

/* Subroutine */ int ctrtri_(char *uplo, char *diag, integer *n, complex *a, 
	integer *lda, integer *info);

/* Subroutine */ int ctrtrs_(char *uplo, char *trans, char *diag, integer *n, 
	integer *nrhs, complex *a, integer *lda, complex *b, integer *ldb, 
	integer *info);

/* Subroutine */ int ctrttf_(char *transr, char *uplo, integer *n, complex *a, 
	 integer *lda, complex *arf, integer *info);

/* Subroutine */ int ctrttp_(char *uplo, integer *n, complex *a, integer *lda, 
	 complex *ap, integer *info);

/* Subroutine */ int ctzrqf_(integer *m, integer *n, complex *a, integer *lda, 
	 complex *tau, integer *info);

/* Subroutine */ int ctzrzf_(integer *m, integer *n, complex *a, integer *lda, 
	 complex *tau, complex *work, integer *lwork, integer *info);

/* Subroutine */ int cung2l_(integer *m, integer *n, integer *k, complex *a, 
	integer *lda, complex *tau, complex *work, integer *info);

/* Subroutine */ int cung2r_(integer *m, integer *n, integer *k, complex *a, 
	integer *lda, complex *tau, complex *work, integer *info);

/* Subroutine */ int cungbr_(char *vect, integer *m, integer *n, integer *k, 
	complex *a, integer *lda, complex *tau, complex *work, integer *lwork, 
	 integer *info);

/* Subroutine */ int cunghr_(integer *n, integer *ilo, integer *ihi, complex *
	a, integer *lda, complex *tau, complex *work, integer *lwork, integer 
	*info);

/* Subroutine */ int cungl2_(integer *m, integer *n, integer *k, complex *a, 
	integer *lda, complex *tau, complex *work, integer *info);

/* Subroutine */ int cunglq_(integer *m, integer *n, integer *k, complex *a, 
	integer *lda, complex *tau, complex *work, integer *lwork, integer *
	info);

/* Subroutine */ int cungql_(integer *m, integer *n, integer *k, complex *a, 
	integer *lda, complex *tau, complex *work, integer *lwork, integer *
	info);

/* Subroutine */ int cungqr_(integer *m, integer *n, integer *k, complex *a, 
	integer *lda, complex *tau, complex *work, integer *lwork, integer *
	info);

/* Subroutine */ int cungr2_(integer *m, integer *n, integer *k, complex *a, 
	integer *lda, complex *tau, complex *work, integer *info);

/* Subroutine */ int cungrq_(integer *m, integer *n, integer *k, complex *a, 
	integer *lda, complex *tau, complex *work, integer *lwork, integer *
	info);

/* Subroutine */ int cungtr_(char *uplo, integer *n, complex *a, integer *lda, 
	 complex *tau, complex *work, integer *lwork, integer *info);

/* Subroutine */ int cunm2l_(char *side, char *trans, integer *m, integer *n, 
	integer *k, complex *a, integer *lda, complex *tau, complex *c__, 
	integer *ldc, complex *work, integer *info);

/* Subroutine */ int cunm2r_(char *side, char *trans, integer *m, integer *n, 
	integer *k, complex *a, integer *lda, complex *tau, complex *c__, 
	integer *ldc, complex *work, integer *info);

/* Subroutine */ int cunmbr_(char *vect, char *side, char *trans, integer *m, 
	integer *n, integer *k, complex *a, integer *lda, complex *tau, 
	complex *c__, integer *ldc, complex *work, integer *lwork, integer *
	info);

/* Subroutine */ int cunmhr_(char *side, char *trans, integer *m, integer *n, 
	integer *ilo, integer *ihi, complex *a, integer *lda, complex *tau, 
	complex *c__, integer *ldc, complex *work, integer *lwork, integer *
	info);

/* Subroutine */ int cunml2_(char *side, char *trans, integer *m, integer *n, 
	integer *k, complex *a, integer *lda, complex *tau, complex *c__, 
	integer *ldc, complex *work, integer *info);

/* Subroutine */ int cunmlq_(char *side, char *trans, integer *m, integer *n, 
	integer *k, complex *a, integer *lda, complex *tau, complex *c__, 
	integer *ldc, complex *work, integer *lwork, integer *info);

/* Subroutine */ int cunmql_(char *side, char *trans, integer *m, integer *n, 
	integer *k, complex *a, integer *lda, complex *tau, complex *c__, 
	integer *ldc, complex *work, integer *lwork, integer *info);

/* Subroutine */ int cunmqr_(char *side, char *trans, integer *m, integer *n, 
	integer *k, complex *a, integer *lda, complex *tau, complex *c__, 
	integer *ldc, complex *work, integer *lwork, integer *info);

/* Subroutine */ int cunmr2_(char *side, char *trans, integer *m, integer *n, 
	integer *k, complex *a, integer *lda, complex *tau, complex *c__, 
	integer *ldc, complex *work, integer *info);

/* Subroutine */ int cunmr3_(char *side, char *trans, integer *m, integer *n, 
	integer *k, integer *l, complex *a, integer *lda, complex *tau, 
	complex *c__, integer *ldc, complex *work, integer *info);

/* Subroutine */ int cunmrq_(char *side, char *trans, integer *m, integer *n, 
	integer *k, complex *a, integer *lda, complex *tau, complex *c__, 
	integer *ldc, complex *work, integer *lwork, integer *info);

/* Subroutine */ int cunmrz_(char *side, char *trans, integer *m, integer *n, 
	integer *k, integer *l, complex *a, integer *lda, complex *tau, 
	complex *c__, integer *ldc, complex *work, integer *lwork, integer *
	info);

/* Subroutine */ int cunmtr_(char *side, char *uplo, char *trans, integer *m, 
	integer *n, complex *a, integer *lda, complex *tau, complex *c__, 
	integer *ldc, complex *work, integer *lwork, integer *info);

/* Subroutine */ int cupgtr_(char *uplo, integer *n, complex *ap, complex *
	tau, complex *q, integer *ldq, complex *work, integer *info);

/* Subroutine */ int cupmtr_(char *side, char *uplo, char *trans, integer *m, 
	integer *n, complex *ap, complex *tau, complex *c__, integer *ldc, 
	complex *work, integer *info);

#pragma endregion



#pragma region SLA

doublereal scsum1_(integer *n, complex *cx, integer *incx);

/* Subroutine */ int sgbbrd_(char *vect, integer *m, integer *n, integer *ncc, 
	 integer *kl, integer *ku, real *ab, integer *ldab, real *d__, real *
	e, real *q, integer *ldq, real *pt, integer *ldpt, real *c__, integer 
	*ldc, real *work, integer *info);

/* Subroutine */ int sgbcon_(char *norm, integer *n, integer *kl, integer *ku, 
	 real *ab, integer *ldab, integer *ipiv, real *anorm, real *rcond, 
	real *work, integer *iwork, integer *info);

/* Subroutine */ int sgbequ_(integer *m, integer *n, integer *kl, integer *ku, 
	 real *ab, integer *ldab, real *r__, real *c__, real *rowcnd, real *
	colcnd, real *amax, integer *info);

/* Subroutine */ int sgbequb_(integer *m, integer *n, integer *kl, integer *
	ku, real *ab, integer *ldab, real *r__, real *c__, real *rowcnd, real 
	*colcnd, real *amax, integer *info);

/* Subroutine */ int sgbrfs_(char *trans, integer *n, integer *kl, integer *
	ku, integer *nrhs, real *ab, integer *ldab, real *afb, integer *ldafb, 
	 integer *ipiv, real *b, integer *ldb, real *x, integer *ldx, real *
	ferr, real *berr, real *work, integer *iwork, integer *info);

/* Subroutine */ int sgbrfsx_(char *trans, char *equed, integer *n, integer *
	kl, integer *ku, integer *nrhs, real *ab, integer *ldab, real *afb, 
	integer *ldafb, integer *ipiv, real *r__, real *c__, real *b, integer 
	*ldb, real *x, integer *ldx, real *rcond, real *berr, integer *
	n_err_bnds__, real *err_bnds_norm__, real *err_bnds_comp__, integer *
	nparams, real *params, real *work, integer *iwork, integer *info);

/* Subroutine */ int sgbsv_(integer *n, integer *kl, integer *ku, integer *
	nrhs, real *ab, integer *ldab, integer *ipiv, real *b, integer *ldb, 
	integer *info);

/* Subroutine */ int sgbsvx_(char *fact, char *trans, integer *n, integer *kl, 
	 integer *ku, integer *nrhs, real *ab, integer *ldab, real *afb, 
	integer *ldafb, integer *ipiv, char *equed, real *r__, real *c__, 
	real *b, integer *ldb, real *x, integer *ldx, real *rcond, real *ferr, 
	 real *berr, real *work, integer *iwork, integer *info);

/* Subroutine */ int sgbsvxx_(char *fact, char *trans, integer *n, integer *
	kl, integer *ku, integer *nrhs, real *ab, integer *ldab, real *afb, 
	integer *ldafb, integer *ipiv, char *equed, real *r__, real *c__, 
	real *b, integer *ldb, real *x, integer *ldx, real *rcond, real *
	rpvgrw, real *berr, integer *n_err_bnds__, real *err_bnds_norm__, 
	real *err_bnds_comp__, integer *nparams, real *params, real *work, 
	integer *iwork, integer *info);

/* Subroutine */ int sgbtf2_(integer *m, integer *n, integer *kl, integer *ku, 
	 real *ab, integer *ldab, integer *ipiv, integer *info);

/* Subroutine */ int sgbtrf_(integer *m, integer *n, integer *kl, integer *ku, 
	 real *ab, integer *ldab, integer *ipiv, integer *info);

/* Subroutine */ int sgbtrs_(char *trans, integer *n, integer *kl, integer *
	ku, integer *nrhs, real *ab, integer *ldab, integer *ipiv, real *b, 
	integer *ldb, integer *info);

/* Subroutine */ int sgebak_(char *job, char *side, integer *n, integer *ilo, 
	integer *ihi, real *scale, integer *m, real *v, integer *ldv, integer 
	*info);

/* Subroutine */ int sgebal_(char *job, integer *n, real *a, integer *lda, 
	integer *ilo, integer *ihi, real *scale, integer *info);

/* Subroutine */ int sgebd2_(integer *m, integer *n, real *a, integer *lda, 
	real *d__, real *e, real *tauq, real *taup, real *work, integer *info);

/* Subroutine */ int sgebrd_(integer *m, integer *n, real *a, integer *lda, 
	real *d__, real *e, real *tauq, real *taup, real *work, integer *
	lwork, integer *info);

/* Subroutine */ int sgecon_(char *norm, integer *n, real *a, integer *lda, 
	real *anorm, real *rcond, real *work, integer *iwork, integer *info);

/* Subroutine */ int sgeequ_(integer *m, integer *n, real *a, integer *lda, 
	real *r__, real *c__, real *rowcnd, real *colcnd, real *amax, integer 
	*info);

/* Subroutine */ int sgeequb_(integer *m, integer *n, real *a, integer *lda, 
	real *r__, real *c__, real *rowcnd, real *colcnd, real *amax, integer 
	*info);

/* Subroutine */ int sgees_(char *jobvs, char *sort, L_fp select, integer *n, 
	real *a, integer *lda, integer *sdim, real *wr, real *wi, real *vs, 
	integer *ldvs, real *work, integer *lwork, logical *bwork, integer *
	info);

/* Subroutine */ int sgeesx_(char *jobvs, char *sort, L_fp select, char *
	sense, integer *n, real *a, integer *lda, integer *sdim, real *wr, 
	real *wi, real *vs, integer *ldvs, real *rconde, real *rcondv, real *
	work, integer *lwork, integer *iwork, integer *liwork, logical *bwork, 
	 integer *info);

/* Subroutine */ int sgeev_(char *jobvl, char *jobvr, integer *n, real *a, 
	integer *lda, real *wr, real *wi, real *vl, integer *ldvl, real *vr, 
	integer *ldvr, real *work, integer *lwork, integer *info);

/* Subroutine */ int sgeevx_(char *balanc, char *jobvl, char *jobvr, char *
	sense, integer *n, real *a, integer *lda, real *wr, real *wi, real *
	vl, integer *ldvl, real *vr, integer *ldvr, integer *ilo, integer *
	ihi, real *scale, real *abnrm, real *rconde, real *rcondv, real *work, 
	 integer *lwork, integer *iwork, integer *info);

/* Subroutine */ int sgegs_(char *jobvsl, char *jobvsr, integer *n, real *a, 
	integer *lda, real *b, integer *ldb, real *alphar, real *alphai, real 
	*beta, real *vsl, integer *ldvsl, real *vsr, integer *ldvsr, real *
	work, integer *lwork, integer *info);

/* Subroutine */ int sgegv_(char *jobvl, char *jobvr, integer *n, real *a, 
	integer *lda, real *b, integer *ldb, real *alphar, real *alphai, real 
	*beta, real *vl, integer *ldvl, real *vr, integer *ldvr, real *work, 
	integer *lwork, integer *info);

/* Subroutine */ int sgehd2_(integer *n, integer *ilo, integer *ihi, real *a, 
	integer *lda, real *tau, real *work, integer *info);

/* Subroutine */ int sgehrd_(integer *n, integer *ilo, integer *ihi, real *a, 
	integer *lda, real *tau, real *work, integer *lwork, integer *info);

/* Subroutine */ int sgejsv_(char *joba, char *jobu, char *jobv, char *jobr, 
	char *jobt, char *jobp, integer *m, integer *n, real *a, integer *lda, 
	 real *sva, real *u, integer *ldu, real *v, integer *ldv, real *work, 
	integer *lwork, integer *iwork, integer *info);

/* Subroutine */ int sgelq2_(integer *m, integer *n, real *a, integer *lda, 
	real *tau, real *work, integer *info);

/* Subroutine */ int sgelqf_(integer *m, integer *n, real *a, integer *lda, 
	real *tau, real *work, integer *lwork, integer *info);

/* Subroutine */ int sgels_(char *trans, integer *m, integer *n, integer *
	nrhs, real *a, integer *lda, real *b, integer *ldb, real *work, 
	integer *lwork, integer *info);

/* Subroutine */ int sgelsd_(integer *m, integer *n, integer *nrhs, real *a, 
	integer *lda, real *b, integer *ldb, real *s, real *rcond, integer *
	rank, real *work, integer *lwork, integer *iwork, integer *info);

/* Subroutine */ int sgelss_(integer *m, integer *n, integer *nrhs, real *a, 
	integer *lda, real *b, integer *ldb, real *s, real *rcond, integer *
	rank, real *work, integer *lwork, integer *info);

/* Subroutine */ int sgelsx_(integer *m, integer *n, integer *nrhs, real *a, 
	integer *lda, real *b, integer *ldb, integer *jpvt, real *rcond, 
	integer *rank, real *work, integer *info);

/* Subroutine */ int sgelsy_(integer *m, integer *n, integer *nrhs, real *a, 
	integer *lda, real *b, integer *ldb, integer *jpvt, real *rcond, 
	integer *rank, real *work, integer *lwork, integer *info);

/* Subroutine */ int sgeql2_(integer *m, integer *n, real *a, integer *lda, 
	real *tau, real *work, integer *info);

/* Subroutine */ int sgeqlf_(integer *m, integer *n, real *a, integer *lda, 
	real *tau, real *work, integer *lwork, integer *info);

/* Subroutine */ int sgeqp3_(integer *m, integer *n, real *a, integer *lda, 
	integer *jpvt, real *tau, real *work, integer *lwork, integer *info);

/* Subroutine */ int sgeqpf_(integer *m, integer *n, real *a, integer *lda, 
	integer *jpvt, real *tau, real *work, integer *info);

/* Subroutine */ int sgeqr2_(integer *m, integer *n, real *a, integer *lda, 
	real *tau, real *work, integer *info);

/* Subroutine */ int sgeqrf_(integer *m, integer *n, real *a, integer *lda, 
	real *tau, real *work, integer *lwork, integer *info);

/* Subroutine */ int sgerfs_(char *trans, integer *n, integer *nrhs, real *a, 
	integer *lda, real *af, integer *ldaf, integer *ipiv, real *b, 
	integer *ldb, real *x, integer *ldx, real *ferr, real *berr, real *
	work, integer *iwork, integer *info);

/* Subroutine */ int sgerfsx_(char *trans, char *equed, integer *n, integer *
	nrhs, real *a, integer *lda, real *af, integer *ldaf, integer *ipiv, 
	real *r__, real *c__, real *b, integer *ldb, real *x, integer *ldx, 
	real *rcond, real *berr, integer *n_err_bnds__, real *err_bnds_norm__, 
	 real *err_bnds_comp__, integer *nparams, real *params, real *work, 
	integer *iwork, integer *info);

/* Subroutine */ int sgerq2_(integer *m, integer *n, real *a, integer *lda, 
	real *tau, real *work, integer *info);

/* Subroutine */ int sgerqf_(integer *m, integer *n, real *a, integer *lda, 
	real *tau, real *work, integer *lwork, integer *info);

/* Subroutine */ int sgesc2_(integer *n, real *a, integer *lda, real *rhs, 
	integer *ipiv, integer *jpiv, real *scale);

/* Subroutine */ int sgesdd_(char *jobz, integer *m, integer *n, real *a, 
	integer *lda, real *s, real *u, integer *ldu, real *vt, integer *ldvt, 
	 real *work, integer *lwork, integer *iwork, integer *info);

/* Subroutine */ int sgesv_(integer *n, integer *nrhs, real *a, integer *lda, 
	integer *ipiv, real *b, integer *ldb, integer *info);

/* Subroutine */ int sgesvd_(char *jobu, char *jobvt, integer *m, integer *n, 
	real *a, integer *lda, real *s, real *u, integer *ldu, real *vt, 
	integer *ldvt, real *work, integer *lwork, integer *info);

/* Subroutine */ int sgesvj_(char *joba, char *jobu, char *jobv, integer *m, 
	integer *n, real *a, integer *lda, real *sva, integer *mv, real *v, 
	integer *ldv, real *work, integer *lwork, integer *info);

/* Subroutine */ int sgesvx_(char *fact, char *trans, integer *n, integer *
	nrhs, real *a, integer *lda, real *af, integer *ldaf, integer *ipiv, 
	char *equed, real *r__, real *c__, real *b, integer *ldb, real *x, 
	integer *ldx, real *rcond, real *ferr, real *berr, real *work, 
	integer *iwork, integer *info);

/* Subroutine */ int sgesvxx_(char *fact, char *trans, integer *n, integer *
	nrhs, real *a, integer *lda, real *af, integer *ldaf, integer *ipiv, 
	char *equed, real *r__, real *c__, real *b, integer *ldb, real *x, 
	integer *ldx, real *rcond, real *rpvgrw, real *berr, integer *
	n_err_bnds__, real *err_bnds_norm__, real *err_bnds_comp__, integer *
	nparams, real *params, real *work, integer *iwork, integer *info);

/* Subroutine */ int sgetc2_(integer *n, real *a, integer *lda, integer *ipiv, 
	 integer *jpiv, integer *info);

/* Subroutine */ int sgetf2_(integer *m, integer *n, real *a, integer *lda, 
	integer *ipiv, integer *info);

/* Subroutine */ int sgetrf_(integer *m, integer *n, real *a, integer *lda, 
	integer *ipiv, integer *info);

/* Subroutine */ int sgetri_(integer *n, real *a, integer *lda, integer *ipiv, 
	 real *work, integer *lwork, integer *info);

/* Subroutine */ int sgetrs_(char *trans, integer *n, integer *nrhs, real *a, 
	integer *lda, integer *ipiv, real *b, integer *ldb, integer *info);

/* Subroutine */ int sggbak_(char *job, char *side, integer *n, integer *ilo, 
	integer *ihi, real *lscale, real *rscale, integer *m, real *v, 
	integer *ldv, integer *info);

/* Subroutine */ int sggbal_(char *job, integer *n, real *a, integer *lda, 
	real *b, integer *ldb, integer *ilo, integer *ihi, real *lscale, real 
	*rscale, real *work, integer *info);

/* Subroutine */ int sgges_(char *jobvsl, char *jobvsr, char *sort, L_fp 
	selctg, integer *n, real *a, integer *lda, real *b, integer *ldb, 
	integer *sdim, real *alphar, real *alphai, real *beta, real *vsl, 
	integer *ldvsl, real *vsr, integer *ldvsr, real *work, integer *lwork, 
	 logical *bwork, integer *info);

/* Subroutine */ int sggesx_(char *jobvsl, char *jobvsr, char *sort, L_fp 
	selctg, char *sense, integer *n, real *a, integer *lda, real *b, 
	integer *ldb, integer *sdim, real *alphar, real *alphai, real *beta, 
	real *vsl, integer *ldvsl, real *vsr, integer *ldvsr, real *rconde, 
	real *rcondv, real *work, integer *lwork, integer *iwork, integer *
	liwork, logical *bwork, integer *info);

/* Subroutine */ int sggev_(char *jobvl, char *jobvr, integer *n, real *a, 
	integer *lda, real *b, integer *ldb, real *alphar, real *alphai, real 
	*beta, real *vl, integer *ldvl, real *vr, integer *ldvr, real *work, 
	integer *lwork, integer *info);

/* Subroutine */ int sggevx_(char *balanc, char *jobvl, char *jobvr, char *
	sense, integer *n, real *a, integer *lda, real *b, integer *ldb, real 
	*alphar, real *alphai, real *beta, real *vl, integer *ldvl, real *vr, 
	integer *ldvr, integer *ilo, integer *ihi, real *lscale, real *rscale, 
	 real *abnrm, real *bbnrm, real *rconde, real *rcondv, real *work, 
	integer *lwork, integer *iwork, logical *bwork, integer *info);

/* Subroutine */ int sggglm_(integer *n, integer *m, integer *p, real *a, 
	integer *lda, real *b, integer *ldb, real *d__, real *x, real *y, 
	real *work, integer *lwork, integer *info);

/* Subroutine */ int sgghrd_(char *compq, char *compz, integer *n, integer *
	ilo, integer *ihi, real *a, integer *lda, real *b, integer *ldb, real 
	*q, integer *ldq, real *z__, integer *ldz, integer *info);

/* Subroutine */ int sgglse_(integer *m, integer *n, integer *p, real *a, 
	integer *lda, real *b, integer *ldb, real *c__, real *d__, real *x, 
	real *work, integer *lwork, integer *info);

/* Subroutine */ int sggqrf_(integer *n, integer *m, integer *p, real *a, 
	integer *lda, real *taua, real *b, integer *ldb, real *taub, real *
	work, integer *lwork, integer *info);

/* Subroutine */ int sggrqf_(integer *m, integer *p, integer *n, real *a, 
	integer *lda, real *taua, real *b, integer *ldb, real *taub, real *
	work, integer *lwork, integer *info);

/* Subroutine */ int sggsvd_(char *jobu, char *jobv, char *jobq, integer *m, 
	integer *n, integer *p, integer *k, integer *l, real *a, integer *lda, 
	 real *b, integer *ldb, real *alpha, real *beta, real *u, integer *
	ldu, real *v, integer *ldv, real *q, integer *ldq, real *work, 
	integer *iwork, integer *info);

/* Subroutine */ int sggsvp_(char *jobu, char *jobv, char *jobq, integer *m, 
	integer *p, integer *n, real *a, integer *lda, real *b, integer *ldb, 
	real *tola, real *tolb, integer *k, integer *l, real *u, integer *ldu, 
	 real *v, integer *ldv, real *q, integer *ldq, integer *iwork, real *
	tau, real *work, integer *info);

/* Subroutine */ int sgsvj0_(char *jobv, integer *m, integer *n, real *a, 
	integer *lda, real *d__, real *sva, integer *mv, real *v, integer *
	ldv, real *eps, real *sfmin, real *tol, integer *nsweep, real *work, 
	integer *lwork, integer *info);

/* Subroutine */ int sgsvj1_(char *jobv, integer *m, integer *n, integer *n1, 
	real *a, integer *lda, real *d__, real *sva, integer *mv, real *v, 
	integer *ldv, real *eps, real *sfmin, real *tol, integer *nsweep, 
	real *work, integer *lwork, integer *info);

/* Subroutine */ int sgtcon_(char *norm, integer *n, real *dl, real *d__, 
	real *du, real *du2, integer *ipiv, real *anorm, real *rcond, real *
	work, integer *iwork, integer *info);

/* Subroutine */ int sgtrfs_(char *trans, integer *n, integer *nrhs, real *dl, 
	 real *d__, real *du, real *dlf, real *df, real *duf, real *du2, 
	integer *ipiv, real *b, integer *ldb, real *x, integer *ldx, real *
	ferr, real *berr, real *work, integer *iwork, integer *info);

/* Subroutine */ int sgtsv_(integer *n, integer *nrhs, real *dl, real *d__, 
	real *du, real *b, integer *ldb, integer *info);

/* Subroutine */ int sgtsvx_(char *fact, char *trans, integer *n, integer *
	nrhs, real *dl, real *d__, real *du, real *dlf, real *df, real *duf, 
	real *du2, integer *ipiv, real *b, integer *ldb, real *x, integer *
	ldx, real *rcond, real *ferr, real *berr, real *work, integer *iwork, 
	integer *info);

/* Subroutine */ int sgttrf_(integer *n, real *dl, real *d__, real *du, real *
	du2, integer *ipiv, integer *info);

/* Subroutine */ int sgttrs_(char *trans, integer *n, integer *nrhs, real *dl, 
	 real *d__, real *du, real *du2, integer *ipiv, real *b, integer *ldb, 
	 integer *info);

/* Subroutine */ int sgtts2_(integer *itrans, integer *n, integer *nrhs, real 
	*dl, real *d__, real *du, real *du2, integer *ipiv, real *b, integer *
	ldb);

/* Subroutine */ int shgeqz_(char *job, char *compq, char *compz, integer *n, 
	integer *ilo, integer *ihi, real *h__, integer *ldh, real *t, integer 
	*ldt, real *alphar, real *alphai, real *beta, real *q, integer *ldq, 
	real *z__, integer *ldz, real *work, integer *lwork, integer *info);

/* Subroutine */ int shsein_(char *side, char *eigsrc, char *initv, logical *
	select, integer *n, real *h__, integer *ldh, real *wr, real *wi, real 
	*vl, integer *ldvl, real *vr, integer *ldvr, integer *mm, integer *m, 
	real *work, integer *ifaill, integer *ifailr, integer *info);

/* Subroutine */ int shseqr_(char *job, char *compz, integer *n, integer *ilo, 
	 integer *ihi, real *h__, integer *ldh, real *wr, real *wi, real *z__, 
	 integer *ldz, real *work, integer *lwork, integer *info);

/* Subroutine */ int sla_gbamv__(integer *trans, integer *m, integer *n, 
	integer *kl, integer *ku, real *alpha, real *ab, integer *ldab, real *
	x, integer *incx, real *beta, real *y, integer *incy);

doublereal sla_gbrcond__(char *trans, integer *n, integer *kl, integer *ku, 
	real *ab, integer *ldab, real *afb, integer *ldafb, integer *ipiv, 
	integer *cmode, real *c__, integer *info, real *work, integer *iwork, 
	ftnlen trans_len);

/* Subroutine */ int sla_gbrfsx_extended__(integer *prec_type__, integer *
	trans_type__, integer *n, integer *kl, integer *ku, integer *nrhs, 
	real *ab, integer *ldab, real *afb, integer *ldafb, integer *ipiv, 
	logical *colequ, real *c__, real *b, integer *ldb, real *y, integer *
	ldy, real *berr_out__, integer *n_norms__, real *errs_n__, real *
	errs_c__, real *res, real *ayb, real *dy, real *y_tail__, real *rcond,
	 integer *ithresh, real *rthresh, real *dz_ub__, logical *
	ignore_cwise__, integer *info);

doublereal sla_gbrpvgrw__(integer *n, integer *kl, integer *ku, integer *
	ncols, real *ab, integer *ldab, real *afb, integer *ldafb);

/* Subroutine */ int sla_geamv__(integer *trans, integer *m, integer *n, real 
	*alpha, real *a, integer *lda, real *x, integer *incx, real *beta, 
	real *y, integer *incy);

doublereal sla_gercond__(char *trans, integer *n, real *a, integer *lda, real 
	*af, integer *ldaf, integer *ipiv, integer *cmode, real *c__, integer 
	*info, real *work, integer *iwork, ftnlen trans_len);

/* Subroutine */ int sla_gerfsx_extended__(integer *prec_type__, integer *
	trans_type__, integer *n, integer *nrhs, real *a, integer *lda, real *
	af, integer *ldaf, integer *ipiv, logical *colequ, real *c__, real *b,
	 integer *ldb, real *y, integer *ldy, real *berr_out__, integer *
	n_norms__, real *errs_n__, real *errs_c__, real *res, real *ayb, real 
	*dy, real *y_tail__, real *rcond, integer *ithresh, real *rthresh, 
	real *dz_ub__, logical *ignore_cwise__, integer *info);

/* Subroutine */ int sla_lin_berr__(integer *n, integer *nz, integer *nrhs, 
	real *res, real *ayb, real *berr);

doublereal sla_porcond__(char *uplo, integer *n, real *a, integer *lda, real *
	af, integer *ldaf, integer *cmode, real *c__, integer *info, real *
	work, integer *iwork, ftnlen uplo_len);

/* Subroutine */ int sla_porfsx_extended__(integer *prec_type__, char *uplo, 
	integer *n, integer *nrhs, real *a, integer *lda, real *af, integer *
	ldaf, logical *colequ, real *c__, real *b, integer *ldb, real *y, 
	integer *ldy, real *berr_out__, integer *n_norms__, real *errs_n__, 
	real *errs_c__, real *res, real *ayb, real *dy, real *y_tail__, real *
	rcond, integer *ithresh, real *rthresh, real *dz_ub__, logical *
	ignore_cwise__, integer *info, ftnlen uplo_len);

doublereal sla_porpvgrw__(char *uplo, integer *ncols, real *a, integer *lda, 
	real *af, integer *ldaf, real *work, ftnlen uplo_len);

doublereal sla_rpvgrw__(integer *n, integer *ncols, real *a, integer *lda, 
	real *af, integer *ldaf);

/* Subroutine */ int sla_syamv__(integer *uplo, integer *n, real *alpha, real 
	*a, integer *lda, real *x, integer *incx, real *beta, real *y, 
	integer *incy);

doublereal sla_syrcond__(char *uplo, integer *n, real *a, integer *lda, real *
	af, integer *ldaf, integer *ipiv, integer *cmode, real *c__, integer *
	info, real *work, integer *iwork, ftnlen uplo_len);

/* Subroutine */ int sla_syrfsx_extended__(integer *prec_type__, char *uplo, 
	integer *n, integer *nrhs, real *a, integer *lda, real *af, integer *
	ldaf, integer *ipiv, logical *colequ, real *c__, real *b, integer *
	ldb, real *y, integer *ldy, real *berr_out__, integer *n_norms__, 
	real *errs_n__, real *errs_c__, real *res, real *ayb, real *dy, real *
	y_tail__, real *rcond, integer *ithresh, real *rthresh, real *dz_ub__,
	 logical *ignore_cwise__, integer *info, ftnlen uplo_len);

doublereal sla_syrpvgrw__(char *uplo, integer *n, integer *info, real *a, 
	integer *lda, real *af, integer *ldaf, integer *ipiv, real *work, 
	ftnlen uplo_len);

/* Subroutine */ int sla_wwaddw__(integer *n, real *x, real *y, real *w);


/* Subroutine */ int slabrd_(integer *m, integer *n, integer *nb, real *a, 
	integer *lda, real *d__, real *e, real *tauq, real *taup, real *x, 
	integer *ldx, real *y, integer *ldy);

/* Subroutine */ int slacn2_(integer *n, real *v, real *x, integer *isgn, 
	real *est, integer *kase, integer *isave);

/* Subroutine */ int slacon_(integer *n, real *v, real *x, integer *isgn, 
	real *est, integer *kase);

/* Subroutine */ int slaein_(logical *rightv, logical *noinit, integer *n, 
	real *h__, integer *ldh, real *wr, real *wi, real *vr, real *vi, real 
	*b, integer *ldb, real *work, real *eps3, real *smlnum, real *bignum, 
	integer *info);

/* Subroutine */ int slaexc_(logical *wantq, integer *n, real *t, integer *
	ldt, real *q, integer *ldq, integer *j1, integer *n1, integer *n2, 
	real *work, integer *info);

/* Subroutine */ int slag2_(real *a, integer *lda, real *b, integer *ldb, 
	real *safmin, real *scale1, real *scale2, real *wr1, real *wr2, real *
	wi);

/* Subroutine */ int slag2d_(integer *m, integer *n, real *sa, integer *ldsa, 
	doublereal *a, integer *lda, integer *info);

/* Subroutine */ int slags2_(logical *upper, real *a1, real *a2, real *a3, 
	real *b1, real *b2, real *b3, real *csu, real *snu, real *csv, real *
	snv, real *csq, real *snq);

/* Subroutine */ int slagtm_(char *trans, integer *n, integer *nrhs, real *
	alpha, real *dl, real *d__, real *du, real *x, integer *ldx, real *
	beta, real *b, integer *ldb);

/* Subroutine */ int slagv2_(real *a, integer *lda, real *b, integer *ldb, 
	real *alphar, real *alphai, real *beta, real *csl, real *snl, real *
	csr, real *snr);

/* Subroutine */ int slahqr_(logical *wantt, logical *wantz, integer *n, 
	integer *ilo, integer *ihi, real *h__, integer *ldh, real *wr, real *
	wi, integer *iloz, integer *ihiz, real *z__, integer *ldz, integer *
	info);

/* Subroutine */ int slahr2_(integer *n, integer *k, integer *nb, real *a, 
	integer *lda, real *tau, real *t, integer *ldt, real *y, integer *ldy);

/* Subroutine */ int slahrd_(integer *n, integer *k, integer *nb, real *a, 
	integer *lda, real *tau, real *t, integer *ldt, real *y, integer *ldy);

/* Subroutine */ int slaic1_(integer *job, integer *j, real *x, real *sest, 
	real *w, real *gamma, real *sestpr, real *s, real *c__);

/* Subroutine */ int slaln2_(logical *ltrans, integer *na, integer *nw, real *
	smin, real *ca, real *a, integer *lda, real *d1, real *d2, real *b, 
	integer *ldb, real *wr, real *wi, real *x, integer *ldx, real *scale, 
	real *xnorm, integer *info);

/* Subroutine */ int slals0_(integer *icompq, integer *nl, integer *nr, 
	integer *sqre, integer *nrhs, real *b, integer *ldb, real *bx, 
	integer *ldbx, integer *perm, integer *givptr, integer *givcol, 
	integer *ldgcol, real *givnum, integer *ldgnum, real *poles, real *
	difl, real *difr, real *z__, integer *k, real *c__, real *s, real *
	work, integer *info);

/* Subroutine */ int slalsa_(integer *icompq, integer *smlsiz, integer *n, 
	integer *nrhs, real *b, integer *ldb, real *bx, integer *ldbx, real *
	u, integer *ldu, real *vt, integer *k, real *difl, real *difr, real *
	z__, real *poles, integer *givptr, integer *givcol, integer *ldgcol, 
	integer *perm, real *givnum, real *c__, real *s, real *work, integer *
	iwork, integer *info);

/* Subroutine */ int slalsd_(char *uplo, integer *smlsiz, integer *n, integer 
	*nrhs, real *d__, real *e, real *b, integer *ldb, real *rcond, 
	integer *rank, real *work, integer *iwork, integer *info);

doublereal slangb_(char *norm, integer *n, integer *kl, integer *ku, real *ab, 
	 integer *ldab, real *work);

doublereal slange_(char *norm, integer *m, integer *n, real *a, integer *lda, 
	real *work);

doublereal slangt_(char *norm, integer *n, real *dl, real *d__, real *du);

doublereal slanhs_(char *norm, integer *n, real *a, integer *lda, real *work);

doublereal slansb_(char *norm, char *uplo, integer *n, integer *k, real *ab, 
	integer *ldab, real *work);

doublereal slansf_(char *norm, char *transr, char *uplo, integer *n, real *a, 
	real *work);

doublereal slansp_(char *norm, char *uplo, integer *n, real *ap, real *work);

doublereal slansy_(char *norm, char *uplo, integer *n, real *a, integer *lda, 
	real *work);

doublereal slantb_(char *norm, char *uplo, char *diag, integer *n, integer *k, 
	 real *ab, integer *ldab, real *work);

doublereal slantp_(char *norm, char *uplo, char *diag, integer *n, real *ap, 
	real *work);

doublereal slantr_(char *norm, char *uplo, char *diag, integer *m, integer *n, 
	 real *a, integer *lda, real *work);

/* Subroutine */ int slanv2_(real *a, real *b, real *c__, real *d__, real *
	rt1r, real *rt1i, real *rt2r, real *rt2i, real *cs, real *sn);

/* Subroutine */ int slapll_(integer *n, real *x, integer *incx, real *y, 
	integer *incy, real *ssmin);

/* Subroutine */ int slapmt_(logical *forwrd, integer *m, integer *n, real *x, 
	 integer *ldx, integer *k);

/* Subroutine */ int slaqgb_(integer *m, integer *n, integer *kl, integer *ku, 
	 real *ab, integer *ldab, real *r__, real *c__, real *rowcnd, real *
	colcnd, real *amax, char *equed);

/* Subroutine */ int slaqge_(integer *m, integer *n, real *a, integer *lda, 
	real *r__, real *c__, real *rowcnd, real *colcnd, real *amax, char *
	equed);

/* Subroutine */ int slaqp2_(integer *m, integer *n, integer *offset, real *a, 
	 integer *lda, integer *jpvt, real *tau, real *vn1, real *vn2, real *
	work);

/* Subroutine */ int slaqps_(integer *m, integer *n, integer *offset, integer 
	*nb, integer *kb, real *a, integer *lda, integer *jpvt, real *tau, 
	real *vn1, real *vn2, real *auxv, real *f, integer *ldf);

/* Subroutine */ int slaqr0_(logical *wantt, logical *wantz, integer *n, 
	integer *ilo, integer *ihi, real *h__, integer *ldh, real *wr, real *
	wi, integer *iloz, integer *ihiz, real *z__, integer *ldz, real *work, 
	 integer *lwork, integer *info);

/* Subroutine */ int slaqr1_(integer *n, real *h__, integer *ldh, real *sr1, 
	real *si1, real *sr2, real *si2, real *v);

/* Subroutine */ int slaqr2_(logical *wantt, logical *wantz, integer *n, 
	integer *ktop, integer *kbot, integer *nw, real *h__, integer *ldh, 
	integer *iloz, integer *ihiz, real *z__, integer *ldz, integer *ns, 
	integer *nd, real *sr, real *si, real *v, integer *ldv, integer *nh, 
	real *t, integer *ldt, integer *nv, real *wv, integer *ldwv, real *
	work, integer *lwork);

/* Subroutine */ int slaqr3_(logical *wantt, logical *wantz, integer *n, 
	integer *ktop, integer *kbot, integer *nw, real *h__, integer *ldh, 
	integer *iloz, integer *ihiz, real *z__, integer *ldz, integer *ns, 
	integer *nd, real *sr, real *si, real *v, integer *ldv, integer *nh, 
	real *t, integer *ldt, integer *nv, real *wv, integer *ldwv, real *
	work, integer *lwork);

/* Subroutine */ int slaqr4_(logical *wantt, logical *wantz, integer *n, 
	integer *ilo, integer *ihi, real *h__, integer *ldh, real *wr, real *
	wi, integer *iloz, integer *ihiz, real *z__, integer *ldz, real *work, 
	 integer *lwork, integer *info);

/* Subroutine */ int slaqr5_(logical *wantt, logical *wantz, integer *kacc22, 
	integer *n, integer *ktop, integer *kbot, integer *nshfts, real *sr, 
	real *si, real *h__, integer *ldh, integer *iloz, integer *ihiz, real 
	*z__, integer *ldz, real *v, integer *ldv, real *u, integer *ldu, 
	integer *nv, real *wv, integer *ldwv, integer *nh, real *wh, integer *
	ldwh);

/* Subroutine */ int slaqsb_(char *uplo, integer *n, integer *kd, real *ab, 
	integer *ldab, real *s, real *scond, real *amax, char *equed);

/* Subroutine */ int slaqsp_(char *uplo, integer *n, real *ap, real *s, real *
	scond, real *amax, char *equed);

/* Subroutine */ int slaqsy_(char *uplo, integer *n, real *a, integer *lda, 
	real *s, real *scond, real *amax, char *equed);

/* Subroutine */ int slaqtr_(logical *ltran, logical *lreal, integer *n, real 
	*t, integer *ldt, real *b, real *w, real *scale, real *x, real *work, 
	integer *info);

/* Subroutine */ int slar1v_(integer *n, integer *b1, integer *bn, real *
	lambda, real *d__, real *l, real *ld, real *lld, real *pivmin, real *
	gaptol, real *z__, logical *wantnc, integer *negcnt, real *ztz, real *
	mingma, integer *r__, integer *isuppz, real *nrminv, real *resid, 
	real *rqcorr, real *work);

/* Subroutine */ int slar2v_(integer *n, real *x, real *y, real *z__, integer 
	*incx, real *c__, real *s, integer *incc);

/* Subroutine */ int slarf_(char *side, integer *m, integer *n, real *v, 
	integer *incv, real *tau, real *c__, integer *ldc, real *work);

/* Subroutine */ int slarfb_(char *side, char *trans, char *direct, char *
	storev, integer *m, integer *n, integer *k, real *v, integer *ldv, 
	real *t, integer *ldt, real *c__, integer *ldc, real *work, integer *
	ldwork);

/* Subroutine */ int slarfg_(integer *n, real *alpha, real *x, integer *incx, 
	real *tau);

/* Subroutine */ int slarfp_(integer *n, real *alpha, real *x, integer *incx, 
	real *tau);

/* Subroutine */ int slarft_(char *direct, char *storev, integer *n, integer *
	k, real *v, integer *ldv, real *tau, real *t, integer *ldt);

/* Subroutine */ int slarfx_(char *side, integer *m, integer *n, real *v, 
	real *tau, real *c__, integer *ldc, real *work);

/* Subroutine */ int slargv_(integer *n, real *x, integer *incx, real *y, 
	integer *incy, real *c__, integer *incc);

/* Subroutine */ int slarrv_(integer *n, real *vl, real *vu, real *d__, real *
	l, real *pivmin, integer *isplit, integer *m, integer *dol, integer *
	dou, real *minrgp, real *rtol1, real *rtol2, real *w, real *werr, 
	real *wgap, integer *iblock, integer *indexw, real *gers, real *z__, 
	integer *ldz, integer *isuppz, real *work, integer *iwork, integer *
	info);

/* Subroutine */ int slarscl2_(integer *m, integer *n, real *d__, real *x, 
	integer *ldx);

/* Subroutine */ int slartv_(integer *n, real *x, integer *incx, real *y, 
	integer *incy, real *c__, real *s, integer *incc);

/* Subroutine */ int slarz_(char *side, integer *m, integer *n, integer *l, 
	real *v, integer *incv, real *tau, real *c__, integer *ldc, real *
	work);

/* Subroutine */ int slarzb_(char *side, char *trans, char *direct, char *
	storev, integer *m, integer *n, integer *k, integer *l, real *v, 
	integer *ldv, real *t, integer *ldt, real *c__, integer *ldc, real *
	work, integer *ldwork);

/* Subroutine */ int slarzt_(char *direct, char *storev, integer *n, integer *
	k, real *v, integer *ldv, real *tau, real *t, integer *ldt);

/* Subroutine */ int slascl2_(integer *m, integer *n, real *d__, real *x, 
	integer *ldx);

/* Subroutine */ int slaswp_(integer *n, real *a, integer *lda, integer *k1, 
	integer *k2, integer *ipiv, integer *incx);

/* Subroutine */ int slasy2_(logical *ltranl, logical *ltranr, integer *isgn, 
	integer *n1, integer *n2, real *tl, integer *ldtl, real *tr, integer *
	ldtr, real *b, integer *ldb, real *scale, real *x, integer *ldx, real 
	*xnorm, integer *info);

/* Subroutine */ int slasyf_(char *uplo, integer *n, integer *nb, integer *kb, 
	 real *a, integer *lda, integer *ipiv, real *w, integer *ldw, integer 
	*info);

/* Subroutine */ int slatbs_(char *uplo, char *trans, char *diag, char *
	normin, integer *n, integer *kd, real *ab, integer *ldab, real *x, 
	real *scale, real *cnorm, integer *info);

/* Subroutine */ int slatdf_(integer *ijob, integer *n, real *z__, integer *
	ldz, real *rhs, real *rdsum, real *rdscal, integer *ipiv, integer *
	jpiv);

/* Subroutine */ int slatps_(char *uplo, char *trans, char *diag, char *
	normin, integer *n, real *ap, real *x, real *scale, real *cnorm, 
	integer *info);

/* Subroutine */ int slatrd_(char *uplo, integer *n, integer *nb, real *a, 
	integer *lda, real *e, real *tau, real *w, integer *ldw);

/* Subroutine */ int slatrs_(char *uplo, char *trans, char *diag, char *
	normin, integer *n, real *a, integer *lda, real *x, real *scale, real 
	*cnorm, integer *info);

/* Subroutine */ int slatrz_(integer *m, integer *n, integer *l, real *a, 
	integer *lda, real *tau, real *work);

/* Subroutine */ int slatzm_(char *side, integer *m, integer *n, real *v, 
	integer *incv, real *tau, real *c1, real *c2, integer *ldc, real *
	work);

/* Subroutine */ int slauu2_(char *uplo, integer *n, real *a, integer *lda, 
	integer *info);

/* Subroutine */ int slauum_(char *uplo, integer *n, real *a, integer *lda, 
	integer *info);

/* Subroutine */ int sopgtr_(char *uplo, integer *n, real *ap, real *tau, 
	real *q, integer *ldq, real *work, integer *info);

/* Subroutine */ int sopmtr_(char *side, char *uplo, char *trans, integer *m, 
	integer *n, real *ap, real *tau, real *c__, integer *ldc, real *work, 
	integer *info);

/* Subroutine */ int sorg2l_(integer *m, integer *n, integer *k, real *a, 
	integer *lda, real *tau, real *work, integer *info);

/* Subroutine */ int sorg2r_(integer *m, integer *n, integer *k, real *a, 
	integer *lda, real *tau, real *work, integer *info);

/* Subroutine */ int sorgbr_(char *vect, integer *m, integer *n, integer *k, 
	real *a, integer *lda, real *tau, real *work, integer *lwork, integer 
	*info);

/* Subroutine */ int sorghr_(integer *n, integer *ilo, integer *ihi, real *a, 
	integer *lda, real *tau, real *work, integer *lwork, integer *info);

/* Subroutine */ int sorgl2_(integer *m, integer *n, integer *k, real *a, 
	integer *lda, real *tau, real *work, integer *info);

/* Subroutine */ int sorglq_(integer *m, integer *n, integer *k, real *a, 
	integer *lda, real *tau, real *work, integer *lwork, integer *info);

/* Subroutine */ int sorgql_(integer *m, integer *n, integer *k, real *a, 
	integer *lda, real *tau, real *work, integer *lwork, integer *info);

/* Subroutine */ int sorgqr_(integer *m, integer *n, integer *k, real *a, 
	integer *lda, real *tau, real *work, integer *lwork, integer *info);

/* Subroutine */ int sorgr2_(integer *m, integer *n, integer *k, real *a, 
	integer *lda, real *tau, real *work, integer *info);

/* Subroutine */ int sorgrq_(integer *m, integer *n, integer *k, real *a, 
	integer *lda, real *tau, real *work, integer *lwork, integer *info);

/* Subroutine */ int sorgtr_(char *uplo, integer *n, real *a, integer *lda, 
	real *tau, real *work, integer *lwork, integer *info);

/* Subroutine */ int sorm2l_(char *side, char *trans, integer *m, integer *n, 
	integer *k, real *a, integer *lda, real *tau, real *c__, integer *ldc, 
	 real *work, integer *info);

/* Subroutine */ int sorm2r_(char *side, char *trans, integer *m, integer *n, 
	integer *k, real *a, integer *lda, real *tau, real *c__, integer *ldc, 
	 real *work, integer *info);

/* Subroutine */ int sormbr_(char *vect, char *side, char *trans, integer *m, 
	integer *n, integer *k, real *a, integer *lda, real *tau, real *c__, 
	integer *ldc, real *work, integer *lwork, integer *info);

/* Subroutine */ int sormhr_(char *side, char *trans, integer *m, integer *n, 
	integer *ilo, integer *ihi, real *a, integer *lda, real *tau, real *
	c__, integer *ldc, real *work, integer *lwork, integer *info);

/* Subroutine */ int sorml2_(char *side, char *trans, integer *m, integer *n, 
	integer *k, real *a, integer *lda, real *tau, real *c__, integer *ldc, 
	 real *work, integer *info);

/* Subroutine */ int sormlq_(char *side, char *trans, integer *m, integer *n, 
	integer *k, real *a, integer *lda, real *tau, real *c__, integer *ldc, 
	 real *work, integer *lwork, integer *info);

/* Subroutine */ int sormql_(char *side, char *trans, integer *m, integer *n, 
	integer *k, real *a, integer *lda, real *tau, real *c__, integer *ldc, 
	 real *work, integer *lwork, integer *info);

/* Subroutine */ int sormqr_(char *side, char *trans, integer *m, integer *n, 
	integer *k, real *a, integer *lda, real *tau, real *c__, integer *ldc, 
	 real *work, integer *lwork, integer *info);

/* Subroutine */ int sormr2_(char *side, char *trans, integer *m, integer *n, 
	integer *k, real *a, integer *lda, real *tau, real *c__, integer *ldc, 
	 real *work, integer *info);

/* Subroutine */ int sormr3_(char *side, char *trans, integer *m, integer *n, 
	integer *k, integer *l, real *a, integer *lda, real *tau, real *c__, 
	integer *ldc, real *work, integer *info);

/* Subroutine */ int sormrq_(char *side, char *trans, integer *m, integer *n, 
	integer *k, real *a, integer *lda, real *tau, real *c__, integer *ldc, 
	 real *work, integer *lwork, integer *info);

/* Subroutine */ int sormrz_(char *side, char *trans, integer *m, integer *n, 
	integer *k, integer *l, real *a, integer *lda, real *tau, real *c__, 
	integer *ldc, real *work, integer *lwork, integer *info);

/* Subroutine */ int sormtr_(char *side, char *uplo, char *trans, integer *m, 
	integer *n, real *a, integer *lda, real *tau, real *c__, integer *ldc, 
	 real *work, integer *lwork, integer *info);

/* Subroutine */ int spbcon_(char *uplo, integer *n, integer *kd, real *ab, 
	integer *ldab, real *anorm, real *rcond, real *work, integer *iwork, 
	integer *info);

/* Subroutine */ int spbequ_(char *uplo, integer *n, integer *kd, real *ab, 
	integer *ldab, real *s, real *scond, real *amax, integer *info);

/* Subroutine */ int spbrfs_(char *uplo, integer *n, integer *kd, integer *
	nrhs, real *ab, integer *ldab, real *afb, integer *ldafb, real *b, 
	integer *ldb, real *x, integer *ldx, real *ferr, real *berr, real *
	work, integer *iwork, integer *info);

/* Subroutine */ int spbstf_(char *uplo, integer *n, integer *kd, real *ab, 
	integer *ldab, integer *info);

/* Subroutine */ int spbsv_(char *uplo, integer *n, integer *kd, integer *
	nrhs, real *ab, integer *ldab, real *b, integer *ldb, integer *info);

/* Subroutine */ int spbsvx_(char *fact, char *uplo, integer *n, integer *kd, 
	integer *nrhs, real *ab, integer *ldab, real *afb, integer *ldafb, 
	char *equed, real *s, real *b, integer *ldb, real *x, integer *ldx, 
	real *rcond, real *ferr, real *berr, real *work, integer *iwork, 
	integer *info);

/* Subroutine */ int spbtf2_(char *uplo, integer *n, integer *kd, real *ab, 
	integer *ldab, integer *info);

/* Subroutine */ int spbtrf_(char *uplo, integer *n, integer *kd, real *ab, 
	integer *ldab, integer *info);

/* Subroutine */ int spbtrs_(char *uplo, integer *n, integer *kd, integer *
	nrhs, real *ab, integer *ldab, real *b, integer *ldb, integer *info);

/* Subroutine */ int spftrf_(char *transr, char *uplo, integer *n, real *a, 
	integer *info);

/* Subroutine */ int spftri_(char *transr, char *uplo, integer *n, real *a, 
	integer *info);

/* Subroutine */ int spftrs_(char *transr, char *uplo, integer *n, integer *
	nrhs, real *a, real *b, integer *ldb, integer *info);

/* Subroutine */ int spocon_(char *uplo, integer *n, real *a, integer *lda, 
	real *anorm, real *rcond, real *work, integer *iwork, integer *info);

/* Subroutine */ int spoequ_(integer *n, real *a, integer *lda, real *s, real 
	*scond, real *amax, integer *info);

/* Subroutine */ int spoequb_(integer *n, real *a, integer *lda, real *s, 
	real *scond, real *amax, integer *info);

/* Subroutine */ int sporfs_(char *uplo, integer *n, integer *nrhs, real *a, 
	integer *lda, real *af, integer *ldaf, real *b, integer *ldb, real *x, 
	 integer *ldx, real *ferr, real *berr, real *work, integer *iwork, 
	integer *info);

/* Subroutine */ int sporfsx_(char *uplo, char *equed, integer *n, integer *
	nrhs, real *a, integer *lda, real *af, integer *ldaf, real *s, real *
	b, integer *ldb, real *x, integer *ldx, real *rcond, real *berr, 
	integer *n_err_bnds__, real *err_bnds_norm__, real *err_bnds_comp__, 
	integer *nparams, real *params, real *work, integer *iwork, integer *
	info);

/* Subroutine */ int sposv_(char *uplo, integer *n, integer *nrhs, real *a, 
	integer *lda, real *b, integer *ldb, integer *info);

/* Subroutine */ int sposvx_(char *fact, char *uplo, integer *n, integer *
	nrhs, real *a, integer *lda, real *af, integer *ldaf, char *equed, 
	real *s, real *b, integer *ldb, real *x, integer *ldx, real *rcond, 
	real *ferr, real *berr, real *work, integer *iwork, integer *info);

/* Subroutine */ int sposvxx_(char *fact, char *uplo, integer *n, integer *
	nrhs, real *a, integer *lda, real *af, integer *ldaf, char *equed, 
	real *s, real *b, integer *ldb, real *x, integer *ldx, real *rcond, 
	real *rpvgrw, real *berr, integer *n_err_bnds__, real *
	err_bnds_norm__, real *err_bnds_comp__, integer *nparams, real *
	params, real *work, integer *iwork, integer *info);

/* Subroutine */ int spotf2_(char *uplo, integer *n, real *a, integer *lda, 
	integer *info);

/* Subroutine */ int spotrf_(char *uplo, integer *n, real *a, integer *lda, 
	integer *info);

/* Subroutine */ int spotri_(char *uplo, integer *n, real *a, integer *lda, 
	integer *info);

/* Subroutine */ int spotrs_(char *uplo, integer *n, integer *nrhs, real *a, 
	integer *lda, real *b, integer *ldb, integer *info);

/* Subroutine */ int sppcon_(char *uplo, integer *n, real *ap, real *anorm, 
	real *rcond, real *work, integer *iwork, integer *info);

/* Subroutine */ int sppequ_(char *uplo, integer *n, real *ap, real *s, real *
	scond, real *amax, integer *info);

/* Subroutine */ int spprfs_(char *uplo, integer *n, integer *nrhs, real *ap, 
	real *afp, real *b, integer *ldb, real *x, integer *ldx, real *ferr, 
	real *berr, real *work, integer *iwork, integer *info);

/* Subroutine */ int sppsv_(char *uplo, integer *n, integer *nrhs, real *ap, 
	real *b, integer *ldb, integer *info);

/* Subroutine */ int sppsvx_(char *fact, char *uplo, integer *n, integer *
	nrhs, real *ap, real *afp, char *equed, real *s, real *b, integer *
	ldb, real *x, integer *ldx, real *rcond, real *ferr, real *berr, real 
	*work, integer *iwork, integer *info);

/* Subroutine */ int spptrf_(char *uplo, integer *n, real *ap, integer *info);

/* Subroutine */ int spptri_(char *uplo, integer *n, real *ap, integer *info);

/* Subroutine */ int spptrs_(char *uplo, integer *n, integer *nrhs, real *ap, 
	real *b, integer *ldb, integer *info);

/* Subroutine */ int spstf2_(char *uplo, integer *n, real *a, integer *lda, 
	integer *piv, integer *rank, real *tol, real *work, integer *info);

/* Subroutine */ int spstrf_(char *uplo, integer *n, real *a, integer *lda, 
	integer *piv, integer *rank, real *tol, real *work, integer *info);

/* Subroutine */ int sptcon_(integer *n, real *d__, real *e, real *anorm, 
	real *rcond, real *work, integer *info);

/* Subroutine */ int spteqr_(char *compz, integer *n, real *d__, real *e, 
	real *z__, integer *ldz, real *work, integer *info);

/* Subroutine */ int sptrfs_(integer *n, integer *nrhs, real *d__, real *e, 
	real *df, real *ef, real *b, integer *ldb, real *x, integer *ldx, 
	real *ferr, real *berr, real *work, integer *info);

/* Subroutine */ int sptsv_(integer *n, integer *nrhs, real *d__, real *e, 
	real *b, integer *ldb, integer *info);

/* Subroutine */ int sptsvx_(char *fact, integer *n, integer *nrhs, real *d__, 
	 real *e, real *df, real *ef, real *b, integer *ldb, real *x, integer 
	*ldx, real *rcond, real *ferr, real *berr, real *work, integer *info);

/* Subroutine */ int spttrs_(integer *n, integer *nrhs, real *d__, real *e, 
	real *b, integer *ldb, integer *info);

/* Subroutine */ int sptts2_(integer *n, integer *nrhs, real *d__, real *e, 
	real *b, integer *ldb);

/* Subroutine */ int srscl_(integer *n, real *sa, real *sx, integer *incx);

/* Subroutine */ int ssbev_(char *jobz, char *uplo, integer *n, integer *kd, 
	real *ab, integer *ldab, real *w, real *z__, integer *ldz, real *work, 
	 integer *info);

/* Subroutine */ int ssbevd_(char *jobz, char *uplo, integer *n, integer *kd, 
	real *ab, integer *ldab, real *w, real *z__, integer *ldz, real *work, 
	 integer *lwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int ssbevx_(char *jobz, char *range, char *uplo, integer *n, 
	integer *kd, real *ab, integer *ldab, real *q, integer *ldq, real *vl, 
	 real *vu, integer *il, integer *iu, real *abstol, integer *m, real *
	w, real *z__, integer *ldz, real *work, integer *iwork, integer *
	ifail, integer *info);

/* Subroutine */ int ssbgst_(char *vect, char *uplo, integer *n, integer *ka, 
	integer *kb, real *ab, integer *ldab, real *bb, integer *ldbb, real *
	x, integer *ldx, real *work, integer *info);

/* Subroutine */ int ssbgv_(char *jobz, char *uplo, integer *n, integer *ka, 
	integer *kb, real *ab, integer *ldab, real *bb, integer *ldbb, real *
	w, real *z__, integer *ldz, real *work, integer *info);

/* Subroutine */ int ssbgvd_(char *jobz, char *uplo, integer *n, integer *ka, 
	integer *kb, real *ab, integer *ldab, real *bb, integer *ldbb, real *
	w, real *z__, integer *ldz, real *work, integer *lwork, integer *
	iwork, integer *liwork, integer *info);

/* Subroutine */ int ssbgvx_(char *jobz, char *range, char *uplo, integer *n, 
	integer *ka, integer *kb, real *ab, integer *ldab, real *bb, integer *
	ldbb, real *q, integer *ldq, real *vl, real *vu, integer *il, integer 
	*iu, real *abstol, integer *m, real *w, real *z__, integer *ldz, real 
	*work, integer *iwork, integer *ifail, integer *info);

/* Subroutine */ int ssbtrd_(char *vect, char *uplo, integer *n, integer *kd, 
	real *ab, integer *ldab, real *d__, real *e, real *q, integer *ldq, 
	real *work, integer *info);

/* Subroutine */ int ssfrk_(char *transr, char *uplo, char *trans, integer *n, 
	 integer *k, real *alpha, real *a, integer *lda, real *beta, real *
	c__);

/* Subroutine */ int sspcon_(char *uplo, integer *n, real *ap, integer *ipiv, 
	real *anorm, real *rcond, real *work, integer *iwork, integer *info);

/* Subroutine */ int sspev_(char *jobz, char *uplo, integer *n, real *ap, 
	real *w, real *z__, integer *ldz, real *work, integer *info);

/* Subroutine */ int sspevd_(char *jobz, char *uplo, integer *n, real *ap, 
	real *w, real *z__, integer *ldz, real *work, integer *lwork, integer 
	*iwork, integer *liwork, integer *info);

/* Subroutine */ int sspevx_(char *jobz, char *range, char *uplo, integer *n, 
	real *ap, real *vl, real *vu, integer *il, integer *iu, real *abstol, 
	integer *m, real *w, real *z__, integer *ldz, real *work, integer *
	iwork, integer *ifail, integer *info);

/* Subroutine */ int sspgst_(integer *itype, char *uplo, integer *n, real *ap, 
	 real *bp, integer *info);

/* Subroutine */ int sspgv_(integer *itype, char *jobz, char *uplo, integer *
	n, real *ap, real *bp, real *w, real *z__, integer *ldz, real *work, 
	integer *info);

/* Subroutine */ int sspgvd_(integer *itype, char *jobz, char *uplo, integer *
	n, real *ap, real *bp, real *w, real *z__, integer *ldz, real *work, 
	integer *lwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int sspgvx_(integer *itype, char *jobz, char *range, char *
	uplo, integer *n, real *ap, real *bp, real *vl, real *vu, integer *il, 
	 integer *iu, real *abstol, integer *m, real *w, real *z__, integer *
	ldz, real *work, integer *iwork, integer *ifail, integer *info);

/* Subroutine */ int ssprfs_(char *uplo, integer *n, integer *nrhs, real *ap, 
	real *afp, integer *ipiv, real *b, integer *ldb, real *x, integer *
	ldx, real *ferr, real *berr, real *work, integer *iwork, integer *
	info);

/* Subroutine */ int sspsv_(char *uplo, integer *n, integer *nrhs, real *ap, 
	integer *ipiv, real *b, integer *ldb, integer *info);

/* Subroutine */ int sspsvx_(char *fact, char *uplo, integer *n, integer *
	nrhs, real *ap, real *afp, integer *ipiv, real *b, integer *ldb, real 
	*x, integer *ldx, real *rcond, real *ferr, real *berr, real *work, 
	integer *iwork, integer *info);

/* Subroutine */ int ssptrd_(char *uplo, integer *n, real *ap, real *d__, 
	real *e, real *tau, integer *info);

/* Subroutine */ int ssptrf_(char *uplo, integer *n, real *ap, integer *ipiv, 
	integer *info);

/* Subroutine */ int ssptri_(char *uplo, integer *n, real *ap, integer *ipiv, 
	real *work, integer *info);

/* Subroutine */ int ssptrs_(char *uplo, integer *n, integer *nrhs, real *ap, 
	integer *ipiv, real *b, integer *ldb, integer *info);

/* Subroutine */ int sstegr_(char *jobz, char *range, integer *n, real *d__, 
	real *e, real *vl, real *vu, integer *il, integer *iu, real *abstol, 
	integer *m, real *w, real *z__, integer *ldz, integer *isuppz, real *
	work, integer *lwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int sstein_(integer *n, real *d__, real *e, integer *m, real 
	*w, integer *iblock, integer *isplit, real *z__, integer *ldz, real *
	work, integer *iwork, integer *ifail, integer *info);

/* Subroutine */ int sstemr_(char *jobz, char *range, integer *n, real *d__, 
	real *e, real *vl, real *vu, integer *il, integer *iu, integer *m, 
	real *w, real *z__, integer *ldz, integer *nzc, integer *isuppz, 
	logical *tryrac, real *work, integer *lwork, integer *iwork, integer *
	liwork, integer *info);

/* Subroutine */ int sstev_(char *jobz, integer *n, real *d__, real *e, real *
	z__, integer *ldz, real *work, integer *info);

/* Subroutine */ int sstevd_(char *jobz, integer *n, real *d__, real *e, real 
	*z__, integer *ldz, real *work, integer *lwork, integer *iwork, 
	integer *liwork, integer *info);

/* Subroutine */ int sstevr_(char *jobz, char *range, integer *n, real *d__, 
	real *e, real *vl, real *vu, integer *il, integer *iu, real *abstol, 
	integer *m, real *w, real *z__, integer *ldz, integer *isuppz, real *
	work, integer *lwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int sstevx_(char *jobz, char *range, integer *n, real *d__, 
	real *e, real *vl, real *vu, integer *il, integer *iu, real *abstol, 
	integer *m, real *w, real *z__, integer *ldz, real *work, integer *
	iwork, integer *ifail, integer *info);

/* Subroutine */ int ssycon_(char *uplo, integer *n, real *a, integer *lda, 
	integer *ipiv, real *anorm, real *rcond, real *work, integer *iwork, 
	integer *info);

/* Subroutine */ int ssyequb_(char *uplo, integer *n, real *a, integer *lda, 
	real *s, real *scond, real *amax, real *work, integer *info);

/* Subroutine */ int ssyev_(char *jobz, char *uplo, integer *n, real *a, 
	integer *lda, real *w, real *work, integer *lwork, integer *info);

/* Subroutine */ int ssyevd_(char *jobz, char *uplo, integer *n, real *a, 
	integer *lda, real *w, real *work, integer *lwork, integer *iwork, 
	integer *liwork, integer *info);

/* Subroutine */ int ssyevr_(char *jobz, char *range, char *uplo, integer *n, 
	real *a, integer *lda, real *vl, real *vu, integer *il, integer *iu, 
	real *abstol, integer *m, real *w, real *z__, integer *ldz, integer *
	isuppz, real *work, integer *lwork, integer *iwork, integer *liwork, 
	integer *info);

/* Subroutine */ int ssyevx_(char *jobz, char *range, char *uplo, integer *n, 
	real *a, integer *lda, real *vl, real *vu, integer *il, integer *iu, 
	real *abstol, integer *m, real *w, real *z__, integer *ldz, real *
	work, integer *lwork, integer *iwork, integer *ifail, integer *info);

/* Subroutine */ int ssygs2_(integer *itype, char *uplo, integer *n, real *a, 
	integer *lda, real *b, integer *ldb, integer *info);

/* Subroutine */ int ssygst_(integer *itype, char *uplo, integer *n, real *a, 
	integer *lda, real *b, integer *ldb, integer *info);

/* Subroutine */ int ssygv_(integer *itype, char *jobz, char *uplo, integer *
	n, real *a, integer *lda, real *b, integer *ldb, real *w, real *work, 
	integer *lwork, integer *info);

/* Subroutine */ int ssygvd_(integer *itype, char *jobz, char *uplo, integer *
	n, real *a, integer *lda, real *b, integer *ldb, real *w, real *work, 
	integer *lwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int ssygvx_(integer *itype, char *jobz, char *range, char *
	uplo, integer *n, real *a, integer *lda, real *b, integer *ldb, real *
	vl, real *vu, integer *il, integer *iu, real *abstol, integer *m, 
	real *w, real *z__, integer *ldz, real *work, integer *lwork, integer 
	*iwork, integer *ifail, integer *info);

/* Subroutine */ int ssyrfs_(char *uplo, integer *n, integer *nrhs, real *a, 
	integer *lda, real *af, integer *ldaf, integer *ipiv, real *b, 
	integer *ldb, real *x, integer *ldx, real *ferr, real *berr, real *
	work, integer *iwork, integer *info);

/* Subroutine */ int ssyrfsx_(char *uplo, char *equed, integer *n, integer *
	nrhs, real *a, integer *lda, real *af, integer *ldaf, integer *ipiv, 
	real *s, real *b, integer *ldb, real *x, integer *ldx, real *rcond, 
	real *berr, integer *n_err_bnds__, real *err_bnds_norm__, real *
	err_bnds_comp__, integer *nparams, real *params, real *work, integer *
	iwork, integer *info);

/* Subroutine */ int ssysv_(char *uplo, integer *n, integer *nrhs, real *a, 
	integer *lda, integer *ipiv, real *b, integer *ldb, real *work, 
	integer *lwork, integer *info);

/* Subroutine */ int ssysvx_(char *fact, char *uplo, integer *n, integer *
	nrhs, real *a, integer *lda, real *af, integer *ldaf, integer *ipiv, 
	real *b, integer *ldb, real *x, integer *ldx, real *rcond, real *ferr, 
	 real *berr, real *work, integer *lwork, integer *iwork, integer *
	info);

/* Subroutine */ int ssysvxx_(char *fact, char *uplo, integer *n, integer *
	nrhs, real *a, integer *lda, real *af, integer *ldaf, integer *ipiv, 
	char *equed, real *s, real *b, integer *ldb, real *x, integer *ldx, 
	real *rcond, real *rpvgrw, real *berr, integer *n_err_bnds__, real *
	err_bnds_norm__, real *err_bnds_comp__, integer *nparams, real *
	params, real *work, integer *iwork, integer *info);

/* Subroutine */ int ssytd2_(char *uplo, integer *n, real *a, integer *lda, 
	real *d__, real *e, real *tau, integer *info);

/* Subroutine */ int ssytf2_(char *uplo, integer *n, real *a, integer *lda, 
	integer *ipiv, integer *info);

/* Subroutine */ int ssytrd_(char *uplo, integer *n, real *a, integer *lda, 
	real *d__, real *e, real *tau, real *work, integer *lwork, integer *
	info);

/* Subroutine */ int ssytrf_(char *uplo, integer *n, real *a, integer *lda, 
	integer *ipiv, real *work, integer *lwork, integer *info);

/* Subroutine */ int ssytri_(char *uplo, integer *n, real *a, integer *lda, 
	integer *ipiv, real *work, integer *info);

/* Subroutine */ int ssytrs_(char *uplo, integer *n, integer *nrhs, real *a, 
	integer *lda, integer *ipiv, real *b, integer *ldb, integer *info);

/* Subroutine */ int stbcon_(char *norm, char *uplo, char *diag, integer *n, 
	integer *kd, real *ab, integer *ldab, real *rcond, real *work, 
	integer *iwork, integer *info);

/* Subroutine */ int stbrfs_(char *uplo, char *trans, char *diag, integer *n, 
	integer *kd, integer *nrhs, real *ab, integer *ldab, real *b, integer 
	*ldb, real *x, integer *ldx, real *ferr, real *berr, real *work, 
	integer *iwork, integer *info);

/* Subroutine */ int stbtrs_(char *uplo, char *trans, char *diag, integer *n, 
	integer *kd, integer *nrhs, real *ab, integer *ldab, real *b, integer 
	*ldb, integer *info);

/* Subroutine */ int stfsm_(char *transr, char *side, char *uplo, char *trans, 
	 char *diag, integer *m, integer *n, real *alpha, real *a, real *b, 
	integer *ldb);

/* Subroutine */ int stftri_(char *transr, char *uplo, char *diag, integer *n, 
	 real *a, integer *info);

/* Subroutine */ int stfttp_(char *transr, char *uplo, integer *n, real *arf, 
	real *ap, integer *info);

/* Subroutine */ int stfttr_(char *transr, char *uplo, integer *n, real *arf, 
	real *a, integer *lda, integer *info);

/* Subroutine */ int stgevc_(char *side, char *howmny, logical *select, 
	integer *n, real *s, integer *lds, real *p, integer *ldp, real *vl, 
	integer *ldvl, real *vr, integer *ldvr, integer *mm, integer *m, real 
	*work, integer *info);

/* Subroutine */ int stgex2_(logical *wantq, logical *wantz, integer *n, real 
	*a, integer *lda, real *b, integer *ldb, real *q, integer *ldq, real *
	z__, integer *ldz, integer *j1, integer *n1, integer *n2, real *work, 
	integer *lwork, integer *info);

/* Subroutine */ int stgexc_(logical *wantq, logical *wantz, integer *n, real 
	*a, integer *lda, real *b, integer *ldb, real *q, integer *ldq, real *
	z__, integer *ldz, integer *ifst, integer *ilst, real *work, integer *
	lwork, integer *info);

/* Subroutine */ int stgsen_(integer *ijob, logical *wantq, logical *wantz, 
	logical *select, integer *n, real *a, integer *lda, real *b, integer *
	ldb, real *alphar, real *alphai, real *beta, real *q, integer *ldq, 
	real *z__, integer *ldz, integer *m, real *pl, real *pr, real *dif, 
	real *work, integer *lwork, integer *iwork, integer *liwork, integer *
	info);

/* Subroutine */ int stgsja_(char *jobu, char *jobv, char *jobq, integer *m, 
	integer *p, integer *n, integer *k, integer *l, real *a, integer *lda, 
	 real *b, integer *ldb, real *tola, real *tolb, real *alpha, real *
	beta, real *u, integer *ldu, real *v, integer *ldv, real *q, integer *
	ldq, real *work, integer *ncycle, integer *info);

/* Subroutine */ int stgsna_(char *job, char *howmny, logical *select, 
	integer *n, real *a, integer *lda, real *b, integer *ldb, real *vl, 
	integer *ldvl, real *vr, integer *ldvr, real *s, real *dif, integer *
	mm, integer *m, real *work, integer *lwork, integer *iwork, integer *
	info);

/* Subroutine */ int stgsy2_(char *trans, integer *ijob, integer *m, integer *
	n, real *a, integer *lda, real *b, integer *ldb, real *c__, integer *
	ldc, real *d__, integer *ldd, real *e, integer *lde, real *f, integer 
	*ldf, real *scale, real *rdsum, real *rdscal, integer *iwork, integer 
	*pq, integer *info);

/* Subroutine */ int stgsyl_(char *trans, integer *ijob, integer *m, integer *
	n, real *a, integer *lda, real *b, integer *ldb, real *c__, integer *
	ldc, real *d__, integer *ldd, real *e, integer *lde, real *f, integer 
	*ldf, real *scale, real *dif, real *work, integer *lwork, integer *
	iwork, integer *info);

/* Subroutine */ int stpcon_(char *norm, char *uplo, char *diag, integer *n, 
	real *ap, real *rcond, real *work, integer *iwork, integer *info);

/* Subroutine */ int stprfs_(char *uplo, char *trans, char *diag, integer *n, 
	integer *nrhs, real *ap, real *b, integer *ldb, real *x, integer *ldx, 
	 real *ferr, real *berr, real *work, integer *iwork, integer *info);

/* Subroutine */ int stptri_(char *uplo, char *diag, integer *n, real *ap, 
	integer *info);

/* Subroutine */ int stptrs_(char *uplo, char *trans, char *diag, integer *n, 
	integer *nrhs, real *ap, real *b, integer *ldb, integer *info);

/* Subroutine */ int stpttf_(char *transr, char *uplo, integer *n, real *ap, 
	real *arf, integer *info);

/* Subroutine */ int stpttr_(char *uplo, integer *n, real *ap, real *a, 
	integer *lda, integer *info);

/* Subroutine */ int strcon_(char *norm, char *uplo, char *diag, integer *n, 
	real *a, integer *lda, real *rcond, real *work, integer *iwork, 
	integer *info);

/* Subroutine */ int strevc_(char *side, char *howmny, logical *select, 
	integer *n, real *t, integer *ldt, real *vl, integer *ldvl, real *vr, 
	integer *ldvr, integer *mm, integer *m, real *work, integer *info);

/* Subroutine */ int strexc_(char *compq, integer *n, real *t, integer *ldt, 
	real *q, integer *ldq, integer *ifst, integer *ilst, real *work, 
	integer *info);

/* Subroutine */ int strrfs_(char *uplo, char *trans, char *diag, integer *n, 
	integer *nrhs, real *a, integer *lda, real *b, integer *ldb, real *x, 
	integer *ldx, real *ferr, real *berr, real *work, integer *iwork, 
	integer *info);

/* Subroutine */ int strsen_(char *job, char *compq, logical *select, integer 
	*n, real *t, integer *ldt, real *q, integer *ldq, real *wr, real *wi, 
	integer *m, real *s, real *sep, real *work, integer *lwork, integer *
	iwork, integer *liwork, integer *info);

/* Subroutine */ int strsna_(char *job, char *howmny, logical *select, 
	integer *n, real *t, integer *ldt, real *vl, integer *ldvl, real *vr, 
	integer *ldvr, real *s, real *sep, integer *mm, integer *m, real *
	work, integer *ldwork, integer *iwork, integer *info);

/* Subroutine */ int strsyl_(char *trana, char *tranb, integer *isgn, integer 
	*m, integer *n, real *a, integer *lda, real *b, integer *ldb, real *
	c__, integer *ldc, real *scale, integer *info);

/* Subroutine */ int strti2_(char *uplo, char *diag, integer *n, real *a, 
	integer *lda, integer *info);

/* Subroutine */ int strtri_(char *uplo, char *diag, integer *n, real *a, 
	integer *lda, integer *info);

/* Subroutine */ int strtrs_(char *uplo, char *trans, char *diag, integer *n, 
	integer *nrhs, real *a, integer *lda, real *b, integer *ldb, integer *
	info);

/* Subroutine */ int strttf_(char *transr, char *uplo, integer *n, real *a, 
	integer *lda, real *arf, integer *info);

/* Subroutine */ int strttp_(char *uplo, integer *n, real *a, integer *lda, 
	real *ap, integer *info);

/* Subroutine */ int stzrqf_(integer *m, integer *n, real *a, integer *lda, 
	real *tau, integer *info);

/* Subroutine */ int stzrzf_(integer *m, integer *n, real *a, integer *lda, 
	real *tau, real *work, integer *lwork, integer *info);

#pragma endregion


#pragma region ZLA


/* Subroutine */ int zbdsqr_(char *uplo, integer *n, integer *ncvt, integer *
	nru, integer *ncc, doublereal *d__, doublereal *e, doublecomplex *vt, 
	integer *ldvt, doublecomplex *u, integer *ldu, doublecomplex *c__, 
	integer *ldc, doublereal *rwork, integer *info);

/* Subroutine */ int zcgesv_(integer *n, integer *nrhs, doublecomplex *a, 
	integer *lda, integer *ipiv, doublecomplex *b, integer *ldb, 
	doublecomplex *x, integer *ldx, doublecomplex *work, complex *swork, 
	doublereal *rwork, integer *iter, integer *info);

/* Subroutine */ int zcposv_(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	doublecomplex *x, integer *ldx, doublecomplex *work, complex *swork, 
	doublereal *rwork, integer *iter, integer *info);

/* Subroutine */ int zdrscl_(integer *n, doublereal *sa, doublecomplex *sx, 
	integer *incx);

/* Subroutine */ int zgbbrd_(char *vect, integer *m, integer *n, integer *ncc, 
	 integer *kl, integer *ku, doublecomplex *ab, integer *ldab, 
	doublereal *d__, doublereal *e, doublecomplex *q, integer *ldq, 
	doublecomplex *pt, integer *ldpt, doublecomplex *c__, integer *ldc, 
	doublecomplex *work, doublereal *rwork, integer *info);

/* Subroutine */ int zgbcon_(char *norm, integer *n, integer *kl, integer *ku, 
	 doublecomplex *ab, integer *ldab, integer *ipiv, doublereal *anorm, 
	doublereal *rcond, doublecomplex *work, doublereal *rwork, integer *
	info);

/* Subroutine */ int zgbequ_(integer *m, integer *n, integer *kl, integer *ku, 
	 doublecomplex *ab, integer *ldab, doublereal *r__, doublereal *c__, 
	doublereal *rowcnd, doublereal *colcnd, doublereal *amax, integer *
	info);

/* Subroutine */ int zgbequb_(integer *m, integer *n, integer *kl, integer *
	ku, doublecomplex *ab, integer *ldab, doublereal *r__, doublereal *
	c__, doublereal *rowcnd, doublereal *colcnd, doublereal *amax, 
	integer *info);

/* Subroutine */ int zgbrfs_(char *trans, integer *n, integer *kl, integer *
	ku, integer *nrhs, doublecomplex *ab, integer *ldab, doublecomplex *
	afb, integer *ldafb, integer *ipiv, doublecomplex *b, integer *ldb, 
	doublecomplex *x, integer *ldx, doublereal *ferr, doublereal *berr, 
	doublecomplex *work, doublereal *rwork, integer *info);

/* Subroutine */ int zgbrfsx_(char *trans, char *equed, integer *n, integer *
	kl, integer *ku, integer *nrhs, doublecomplex *ab, integer *ldab, 
	doublecomplex *afb, integer *ldafb, integer *ipiv, doublereal *r__, 
	doublereal *c__, doublecomplex *b, integer *ldb, doublecomplex *x, 
	integer *ldx, doublereal *rcond, doublereal *berr, integer *
	n_err_bnds__, doublereal *err_bnds_norm__, doublereal *
	err_bnds_comp__, integer *nparams, doublereal *params, doublecomplex *
	work, doublereal *rwork, integer *info);

/* Subroutine */ int zgbsv_(integer *n, integer *kl, integer *ku, integer *
	nrhs, doublecomplex *ab, integer *ldab, integer *ipiv, doublecomplex *
	b, integer *ldb, integer *info);

/* Subroutine */ int zgbsvx_(char *fact, char *trans, integer *n, integer *kl, 
	 integer *ku, integer *nrhs, doublecomplex *ab, integer *ldab, 
	doublecomplex *afb, integer *ldafb, integer *ipiv, char *equed, 
	doublereal *r__, doublereal *c__, doublecomplex *b, integer *ldb, 
	doublecomplex *x, integer *ldx, doublereal *rcond, doublereal *ferr, 
	doublereal *berr, doublecomplex *work, doublereal *rwork, integer *
	info);

/* Subroutine */ int zgbsvxx_(char *fact, char *trans, integer *n, integer *
	kl, integer *ku, integer *nrhs, doublecomplex *ab, integer *ldab, 
	doublecomplex *afb, integer *ldafb, integer *ipiv, char *equed, 
	doublereal *r__, doublereal *c__, doublecomplex *b, integer *ldb, 
	doublecomplex *x, integer *ldx, doublereal *rcond, doublereal *rpvgrw, 
	 doublereal *berr, integer *n_err_bnds__, doublereal *err_bnds_norm__, 
	 doublereal *err_bnds_comp__, integer *nparams, doublereal *params, 
	doublecomplex *work, doublereal *rwork, integer *info);

/* Subroutine */ int zgbtf2_(integer *m, integer *n, integer *kl, integer *ku, 
	 doublecomplex *ab, integer *ldab, integer *ipiv, integer *info);

/* Subroutine */ int zgbtrf_(integer *m, integer *n, integer *kl, integer *ku, 
	 doublecomplex *ab, integer *ldab, integer *ipiv, integer *info);

/* Subroutine */ int zgbtrs_(char *trans, integer *n, integer *kl, integer *
	ku, integer *nrhs, doublecomplex *ab, integer *ldab, integer *ipiv, 
	doublecomplex *b, integer *ldb, integer *info);

/* Subroutine */ int zgebak_(char *job, char *side, integer *n, integer *ilo, 
	integer *ihi, doublereal *scale, integer *m, doublecomplex *v, 
	integer *ldv, integer *info);

/* Subroutine */ int zgebal_(char *job, integer *n, doublecomplex *a, integer 
	*lda, integer *ilo, integer *ihi, doublereal *scale, integer *info);

/* Subroutine */ int zgebd2_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublereal *d__, doublereal *e, doublecomplex *tauq, 
	doublecomplex *taup, doublecomplex *work, integer *info);

/* Subroutine */ int zgebrd_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublereal *d__, doublereal *e, doublecomplex *tauq, 
	doublecomplex *taup, doublecomplex *work, integer *lwork, integer *
	info);

/* Subroutine */ int zgecon_(char *norm, integer *n, doublecomplex *a, 
	integer *lda, doublereal *anorm, doublereal *rcond, doublecomplex *
	work, doublereal *rwork, integer *info);

/* Subroutine */ int zgeequ_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublereal *r__, doublereal *c__, doublereal *rowcnd, 
	doublereal *colcnd, doublereal *amax, integer *info);

/* Subroutine */ int zgeequb_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublereal *r__, doublereal *c__, doublereal *rowcnd, 
	doublereal *colcnd, doublereal *amax, integer *info);

/* Subroutine */ int zgees_(char *jobvs, char *sort, L_fp select, integer *n, 
	doublecomplex *a, integer *lda, integer *sdim, doublecomplex *w, 
	doublecomplex *vs, integer *ldvs, doublecomplex *work, integer *lwork, 
	 doublereal *rwork, logical *bwork, integer *info);

/* Subroutine */ int zgeesx_(char *jobvs, char *sort, L_fp select, char *
	sense, integer *n, doublecomplex *a, integer *lda, integer *sdim, 
	doublecomplex *w, doublecomplex *vs, integer *ldvs, doublereal *
	rconde, doublereal *rcondv, doublecomplex *work, integer *lwork, 
	doublereal *rwork, logical *bwork, integer *info);

/* Subroutine */ int zgeev_(char *jobvl, char *jobvr, integer *n, 
	doublecomplex *a, integer *lda, doublecomplex *w, doublecomplex *vl, 
	integer *ldvl, doublecomplex *vr, integer *ldvr, doublecomplex *work, 
	integer *lwork, doublereal *rwork, integer *info);

/* Subroutine */ int zgeevx_(char *balanc, char *jobvl, char *jobvr, char *
	sense, integer *n, doublecomplex *a, integer *lda, doublecomplex *w, 
	doublecomplex *vl, integer *ldvl, doublecomplex *vr, integer *ldvr, 
	integer *ilo, integer *ihi, doublereal *scale, doublereal *abnrm, 
	doublereal *rconde, doublereal *rcondv, doublecomplex *work, integer *
	lwork, doublereal *rwork, integer *info);

/* Subroutine */ int zgegs_(char *jobvsl, char *jobvsr, integer *n, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	doublecomplex *alpha, doublecomplex *beta, doublecomplex *vsl, 
	integer *ldvsl, doublecomplex *vsr, integer *ldvsr, doublecomplex *
	work, integer *lwork, doublereal *rwork, integer *info);

/* Subroutine */ int zgegv_(char *jobvl, char *jobvr, integer *n, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	doublecomplex *alpha, doublecomplex *beta, doublecomplex *vl, integer 
	*ldvl, doublecomplex *vr, integer *ldvr, doublecomplex *work, integer 
	*lwork, doublereal *rwork, integer *info);

/* Subroutine */ int zgehd2_(integer *n, integer *ilo, integer *ihi, 
	doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *
	work, integer *info);

/* Subroutine */ int zgehrd_(integer *n, integer *ilo, integer *ihi, 
	doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *
	work, integer *lwork, integer *info);

/* Subroutine */ int zgelq2_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublecomplex *tau, doublecomplex *work, integer *info);

/* Subroutine */ int zgelqf_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublecomplex *tau, doublecomplex *work, integer *lwork, 
	 integer *info);

/* Subroutine */ int zgels_(char *trans, integer *m, integer *n, integer *
	nrhs, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	doublecomplex *work, integer *lwork, integer *info);

/* Subroutine */ int zgelsd_(integer *m, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	doublereal *s, doublereal *rcond, integer *rank, doublecomplex *work, 
	integer *lwork, doublereal *rwork, integer *iwork, integer *info);

/* Subroutine */ int zgelss_(integer *m, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	doublereal *s, doublereal *rcond, integer *rank, doublecomplex *work, 
	integer *lwork, doublereal *rwork, integer *info);

/* Subroutine */ int zgelsx_(integer *m, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	integer *jpvt, doublereal *rcond, integer *rank, doublecomplex *work, 
	doublereal *rwork, integer *info);

/* Subroutine */ int zgelsy_(integer *m, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	integer *jpvt, doublereal *rcond, integer *rank, doublecomplex *work, 
	integer *lwork, doublereal *rwork, integer *info);

/* Subroutine */ int zgeql2_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublecomplex *tau, doublecomplex *work, integer *info);

/* Subroutine */ int zgeqlf_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublecomplex *tau, doublecomplex *work, integer *lwork, 
	 integer *info);

/* Subroutine */ int zgeqp3_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, integer *jpvt, doublecomplex *tau, doublecomplex *work, 
	integer *lwork, doublereal *rwork, integer *info);

/* Subroutine */ int zgeqpf_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, integer *jpvt, doublecomplex *tau, doublecomplex *work, 
	doublereal *rwork, integer *info);

/* Subroutine */ int zgeqr2_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublecomplex *tau, doublecomplex *work, integer *info);

/* Subroutine */ int zgeqrf_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublecomplex *tau, doublecomplex *work, integer *lwork, 
	 integer *info);

/* Subroutine */ int zgerfs_(char *trans, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, doublecomplex *af, integer *ldaf, 
	integer *ipiv, doublecomplex *b, integer *ldb, doublecomplex *x, 
	integer *ldx, doublereal *ferr, doublereal *berr, doublecomplex *work, 
	 doublereal *rwork, integer *info);

/* Subroutine */ int zgerfsx_(char *trans, char *equed, integer *n, integer *
	nrhs, doublecomplex *a, integer *lda, doublecomplex *af, integer *
	ldaf, integer *ipiv, doublereal *r__, doublereal *c__, doublecomplex *
	b, integer *ldb, doublecomplex *x, integer *ldx, doublereal *rcond, 
	doublereal *berr, integer *n_err_bnds__, doublereal *err_bnds_norm__, 
	doublereal *err_bnds_comp__, integer *nparams, doublereal *params, 
	doublecomplex *work, doublereal *rwork, integer *info);

/* Subroutine */ int zgerq2_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublecomplex *tau, doublecomplex *work, integer *info);

/* Subroutine */ int zgerqf_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublecomplex *tau, doublecomplex *work, integer *lwork, 
	 integer *info);

/* Subroutine */ int zgesc2_(integer *n, doublecomplex *a, integer *lda, 
	doublecomplex *rhs, integer *ipiv, integer *jpiv, doublereal *scale);

/* Subroutine */ int zgesdd_(char *jobz, integer *m, integer *n, 
	doublecomplex *a, integer *lda, doublereal *s, doublecomplex *u, 
	integer *ldu, doublecomplex *vt, integer *ldvt, doublecomplex *work, 
	integer *lwork, doublereal *rwork, integer *iwork, integer *info);

/* Subroutine */ int zgesv_(integer *n, integer *nrhs, doublecomplex *a, 
	integer *lda, integer *ipiv, doublecomplex *b, integer *ldb, integer *
	info);

/* Subroutine */ int zgesvd_(char *jobu, char *jobvt, integer *m, integer *n, 
	doublecomplex *a, integer *lda, doublereal *s, doublecomplex *u, 
	integer *ldu, doublecomplex *vt, integer *ldvt, doublecomplex *work, 
	integer *lwork, doublereal *rwork, integer *info);

/* Subroutine */ int zgesvx_(char *fact, char *trans, integer *n, integer *
	nrhs, doublecomplex *a, integer *lda, doublecomplex *af, integer *
	ldaf, integer *ipiv, char *equed, doublereal *r__, doublereal *c__, 
	doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx, 
	doublereal *rcond, doublereal *ferr, doublereal *berr, doublecomplex *
	work, doublereal *rwork, integer *info);

/* Subroutine */ int zgesvxx_(char *fact, char *trans, integer *n, integer *
	nrhs, doublecomplex *a, integer *lda, doublecomplex *af, integer *
	ldaf, integer *ipiv, char *equed, doublereal *r__, doublereal *c__, 
	doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx, 
	doublereal *rcond, doublereal *rpvgrw, doublereal *berr, integer *
	n_err_bnds__, doublereal *err_bnds_norm__, doublereal *
	err_bnds_comp__, integer *nparams, doublereal *params, doublecomplex *
	work, doublereal *rwork, integer *info);

/* Subroutine */ int zgetc2_(integer *n, doublecomplex *a, integer *lda, 
	integer *ipiv, integer *jpiv, integer *info);

/* Subroutine */ int zgetf2_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, integer *ipiv, integer *info);

/* Subroutine */ int zgetrf_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, integer *ipiv, integer *info);

/* Subroutine */ int zgetri_(integer *n, doublecomplex *a, integer *lda, 
	integer *ipiv, doublecomplex *work, integer *lwork, integer *info);

/* Subroutine */ int zgetrs_(char *trans, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, integer *ipiv, doublecomplex *b, 
	integer *ldb, integer *info);

/* Subroutine */ int zggbak_(char *job, char *side, integer *n, integer *ilo, 
	integer *ihi, doublereal *lscale, doublereal *rscale, integer *m, 
	doublecomplex *v, integer *ldv, integer *info);

/* Subroutine */ int zggbal_(char *job, integer *n, doublecomplex *a, integer 
	*lda, doublecomplex *b, integer *ldb, integer *ilo, integer *ihi, 
	doublereal *lscale, doublereal *rscale, doublereal *work, integer *
	info);

/* Subroutine */ int zgges_(char *jobvsl, char *jobvsr, char *sort, L_fp 
	selctg, integer *n, doublecomplex *a, integer *lda, doublecomplex *b, 
	integer *ldb, integer *sdim, doublecomplex *alpha, doublecomplex *
	beta, doublecomplex *vsl, integer *ldvsl, doublecomplex *vsr, integer 
	*ldvsr, doublecomplex *work, integer *lwork, doublereal *rwork, 
	logical *bwork, integer *info);

/* Subroutine */ int zggesx_(char *jobvsl, char *jobvsr, char *sort, L_fp 
	selctg, char *sense, integer *n, doublecomplex *a, integer *lda, 
	doublecomplex *b, integer *ldb, integer *sdim, doublecomplex *alpha, 
	doublecomplex *beta, doublecomplex *vsl, integer *ldvsl, 
	doublecomplex *vsr, integer *ldvsr, doublereal *rconde, doublereal *
	rcondv, doublecomplex *work, integer *lwork, doublereal *rwork, 
	integer *iwork, integer *liwork, logical *bwork, integer *info);

/* Subroutine */ int zggev_(char *jobvl, char *jobvr, integer *n, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	doublecomplex *alpha, doublecomplex *beta, doublecomplex *vl, integer 
	*ldvl, doublecomplex *vr, integer *ldvr, doublecomplex *work, integer 
	*lwork, doublereal *rwork, integer *info);

/* Subroutine */ int zggevx_(char *balanc, char *jobvl, char *jobvr, char *
	sense, integer *n, doublecomplex *a, integer *lda, doublecomplex *b, 
	integer *ldb, doublecomplex *alpha, doublecomplex *beta, 
	doublecomplex *vl, integer *ldvl, doublecomplex *vr, integer *ldvr, 
	integer *ilo, integer *ihi, doublereal *lscale, doublereal *rscale, 
	doublereal *abnrm, doublereal *bbnrm, doublereal *rconde, doublereal *
	rcondv, doublecomplex *work, integer *lwork, doublereal *rwork, 
	integer *iwork, logical *bwork, integer *info);

/* Subroutine */ int zggglm_(integer *n, integer *m, integer *p, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	doublecomplex *d__, doublecomplex *x, doublecomplex *y, doublecomplex 
	*work, integer *lwork, integer *info);

/* Subroutine */ int zgghrd_(char *compq, char *compz, integer *n, integer *
	ilo, integer *ihi, doublecomplex *a, integer *lda, doublecomplex *b, 
	integer *ldb, doublecomplex *q, integer *ldq, doublecomplex *z__, 
	integer *ldz, integer *info);

/* Subroutine */ int zgglse_(integer *m, integer *n, integer *p, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	doublecomplex *c__, doublecomplex *d__, doublecomplex *x, 
	doublecomplex *work, integer *lwork, integer *info);

/* Subroutine */ int zggqrf_(integer *n, integer *m, integer *p, 
	doublecomplex *a, integer *lda, doublecomplex *taua, doublecomplex *b, 
	 integer *ldb, doublecomplex *taub, doublecomplex *work, integer *
	lwork, integer *info);

/* Subroutine */ int zggrqf_(integer *m, integer *p, integer *n, 
	doublecomplex *a, integer *lda, doublecomplex *taua, doublecomplex *b, 
	 integer *ldb, doublecomplex *taub, doublecomplex *work, integer *
	lwork, integer *info);

/* Subroutine */ int zggsvd_(char *jobu, char *jobv, char *jobq, integer *m, 
	integer *n, integer *p, integer *k, integer *l, doublecomplex *a, 
	integer *lda, doublecomplex *b, integer *ldb, doublereal *alpha, 
	doublereal *beta, doublecomplex *u, integer *ldu, doublecomplex *v, 
	integer *ldv, doublecomplex *q, integer *ldq, doublecomplex *work, 
	doublereal *rwork, integer *iwork, integer *info);

/* Subroutine */ int zggsvp_(char *jobu, char *jobv, char *jobq, integer *m, 
	integer *p, integer *n, doublecomplex *a, integer *lda, doublecomplex 
	*b, integer *ldb, doublereal *tola, doublereal *tolb, integer *k, 
	integer *l, doublecomplex *u, integer *ldu, doublecomplex *v, integer 
	*ldv, doublecomplex *q, integer *ldq, integer *iwork, doublereal *
	rwork, doublecomplex *tau, doublecomplex *work, integer *info);

/* Subroutine */ int zgtcon_(char *norm, integer *n, doublecomplex *dl, 
	doublecomplex *d__, doublecomplex *du, doublecomplex *du2, integer *
	ipiv, doublereal *anorm, doublereal *rcond, doublecomplex *work, 
	integer *info);

/* Subroutine */ int zgtrfs_(char *trans, integer *n, integer *nrhs, 
	doublecomplex *dl, doublecomplex *d__, doublecomplex *du, 
	doublecomplex *dlf, doublecomplex *df, doublecomplex *duf, 
	doublecomplex *du2, integer *ipiv, doublecomplex *b, integer *ldb, 
	doublecomplex *x, integer *ldx, doublereal *ferr, doublereal *berr, 
	doublecomplex *work, doublereal *rwork, integer *info);

/* Subroutine */ int zgtsv_(integer *n, integer *nrhs, doublecomplex *dl, 
	doublecomplex *d__, doublecomplex *du, doublecomplex *b, integer *ldb, 
	 integer *info);

/* Subroutine */ int zgtsvx_(char *fact, char *trans, integer *n, integer *
	nrhs, doublecomplex *dl, doublecomplex *d__, doublecomplex *du, 
	doublecomplex *dlf, doublecomplex *df, doublecomplex *duf, 
	doublecomplex *du2, integer *ipiv, doublecomplex *b, integer *ldb, 
	doublecomplex *x, integer *ldx, doublereal *rcond, doublereal *ferr, 
	doublereal *berr, doublecomplex *work, doublereal *rwork, integer *
	info);

/* Subroutine */ int zgttrf_(integer *n, doublecomplex *dl, doublecomplex *
	d__, doublecomplex *du, doublecomplex *du2, integer *ipiv, integer *
	info);

/* Subroutine */ int zgttrs_(char *trans, integer *n, integer *nrhs, 
	doublecomplex *dl, doublecomplex *d__, doublecomplex *du, 
	doublecomplex *du2, integer *ipiv, doublecomplex *b, integer *ldb, 
	integer *info);

/* Subroutine */ int zgtts2_(integer *itrans, integer *n, integer *nrhs, 
	doublecomplex *dl, doublecomplex *d__, doublecomplex *du, 
	doublecomplex *du2, integer *ipiv, doublecomplex *b, integer *ldb);

/* Subroutine */ int zhbev_(char *jobz, char *uplo, integer *n, integer *kd, 
	doublecomplex *ab, integer *ldab, doublereal *w, doublecomplex *z__, 
	integer *ldz, doublecomplex *work, doublereal *rwork, integer *info);

/* Subroutine */ int zhbevd_(char *jobz, char *uplo, integer *n, integer *kd, 
	doublecomplex *ab, integer *ldab, doublereal *w, doublecomplex *z__, 
	integer *ldz, doublecomplex *work, integer *lwork, doublereal *rwork, 
	integer *lrwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int zhbevx_(char *jobz, char *range, char *uplo, integer *n, 
	integer *kd, doublecomplex *ab, integer *ldab, doublecomplex *q, 
	integer *ldq, doublereal *vl, doublereal *vu, integer *il, integer *
	iu, doublereal *abstol, integer *m, doublereal *w, doublecomplex *z__, 
	 integer *ldz, doublecomplex *work, doublereal *rwork, integer *iwork, 
	 integer *ifail, integer *info);

/* Subroutine */ int zhbgst_(char *vect, char *uplo, integer *n, integer *ka, 
	integer *kb, doublecomplex *ab, integer *ldab, doublecomplex *bb, 
	integer *ldbb, doublecomplex *x, integer *ldx, doublecomplex *work, 
	doublereal *rwork, integer *info);

/* Subroutine */ int zhbgv_(char *jobz, char *uplo, integer *n, integer *ka, 
	integer *kb, doublecomplex *ab, integer *ldab, doublecomplex *bb, 
	integer *ldbb, doublereal *w, doublecomplex *z__, integer *ldz, 
	doublecomplex *work, doublereal *rwork, integer *info);

/* Subroutine */ int zhbgvd_(char *jobz, char *uplo, integer *n, integer *ka, 
	integer *kb, doublecomplex *ab, integer *ldab, doublecomplex *bb, 
	integer *ldbb, doublereal *w, doublecomplex *z__, integer *ldz, 
	doublecomplex *work, integer *lwork, doublereal *rwork, integer *
	lrwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int zhbgvx_(char *jobz, char *range, char *uplo, integer *n, 
	integer *ka, integer *kb, doublecomplex *ab, integer *ldab, 
	doublecomplex *bb, integer *ldbb, doublecomplex *q, integer *ldq, 
	doublereal *vl, doublereal *vu, integer *il, integer *iu, doublereal *
	abstol, integer *m, doublereal *w, doublecomplex *z__, integer *ldz, 
	doublecomplex *work, doublereal *rwork, integer *iwork, integer *
	ifail, integer *info);

/* Subroutine */ int zhbtrd_(char *vect, char *uplo, integer *n, integer *kd, 
	doublecomplex *ab, integer *ldab, doublereal *d__, doublereal *e, 
	doublecomplex *q, integer *ldq, doublecomplex *work, integer *info);

/* Subroutine */ int zhecon_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *ipiv, doublereal *anorm, doublereal *rcond, 
	doublecomplex *work, integer *info);

/* Subroutine */ int zheequb_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, doublereal *s, doublereal *scond, doublereal *amax, 
	doublecomplex *work, integer *info);

/* Subroutine */ int zheev_(char *jobz, char *uplo, integer *n, doublecomplex 
	*a, integer *lda, doublereal *w, doublecomplex *work, integer *lwork, 
	doublereal *rwork, integer *info);

/* Subroutine */ int zheevd_(char *jobz, char *uplo, integer *n, 
	doublecomplex *a, integer *lda, doublereal *w, doublecomplex *work, 
	integer *lwork, doublereal *rwork, integer *lrwork, integer *iwork, 
	integer *liwork, integer *info);

/* Subroutine */ int zheevr_(char *jobz, char *range, char *uplo, integer *n, 
	doublecomplex *a, integer *lda, doublereal *vl, doublereal *vu, 
	integer *il, integer *iu, doublereal *abstol, integer *m, doublereal *
	w, doublecomplex *z__, integer *ldz, integer *isuppz, doublecomplex *
	work, integer *lwork, doublereal *rwork, integer *lrwork, integer *
	iwork, integer *liwork, integer *info);

/* Subroutine */ int zheevx_(char *jobz, char *range, char *uplo, integer *n, 
	doublecomplex *a, integer *lda, doublereal *vl, doublereal *vu, 
	integer *il, integer *iu, doublereal *abstol, integer *m, doublereal *
	w, doublecomplex *z__, integer *ldz, doublecomplex *work, integer *
	lwork, doublereal *rwork, integer *iwork, integer *ifail, integer *
	info);

/* Subroutine */ int zhegs2_(integer *itype, char *uplo, integer *n, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	integer *info);

/* Subroutine */ int zhegst_(integer *itype, char *uplo, integer *n, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	integer *info);

/* Subroutine */ int zhegv_(integer *itype, char *jobz, char *uplo, integer *
	n, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	doublereal *w, doublecomplex *work, integer *lwork, doublereal *rwork, 
	 integer *info);

/* Subroutine */ int zhegvd_(integer *itype, char *jobz, char *uplo, integer *
	n, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	doublereal *w, doublecomplex *work, integer *lwork, doublereal *rwork, 
	 integer *lrwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int zhegvx_(integer *itype, char *jobz, char *range, char *
	uplo, integer *n, doublecomplex *a, integer *lda, doublecomplex *b, 
	integer *ldb, doublereal *vl, doublereal *vu, integer *il, integer *
	iu, doublereal *abstol, integer *m, doublereal *w, doublecomplex *z__, 
	 integer *ldz, doublecomplex *work, integer *lwork, doublereal *rwork, 
	 integer *iwork, integer *ifail, integer *info);

/* Subroutine */ int zherfs_(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, doublecomplex *af, integer *ldaf, 
	integer *ipiv, doublecomplex *b, integer *ldb, doublecomplex *x, 
	integer *ldx, doublereal *ferr, doublereal *berr, doublecomplex *work, 
	 doublereal *rwork, integer *info);

/* Subroutine */ int zherfsx_(char *uplo, char *equed, integer *n, integer *
	nrhs, doublecomplex *a, integer *lda, doublecomplex *af, integer *
	ldaf, integer *ipiv, doublereal *s, doublecomplex *b, integer *ldb, 
	doublecomplex *x, integer *ldx, doublereal *rcond, doublereal *berr, 
	integer *n_err_bnds__, doublereal *err_bnds_norm__, doublereal *
	err_bnds_comp__, integer *nparams, doublereal *params, doublecomplex *
	work, doublereal *rwork, integer *info);

/* Subroutine */ int zhesv_(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, integer *ipiv, doublecomplex *b, 
	integer *ldb, doublecomplex *work, integer *lwork, integer *info);

/* Subroutine */ int zhesvx_(char *fact, char *uplo, integer *n, integer *
	nrhs, doublecomplex *a, integer *lda, doublecomplex *af, integer *
	ldaf, integer *ipiv, doublecomplex *b, integer *ldb, doublecomplex *x, 
	 integer *ldx, doublereal *rcond, doublereal *ferr, doublereal *berr, 
	doublecomplex *work, integer *lwork, doublereal *rwork, integer *info);

/* Subroutine */ int zhesvxx_(char *fact, char *uplo, integer *n, integer *
	nrhs, doublecomplex *a, integer *lda, doublecomplex *af, integer *
	ldaf, integer *ipiv, char *equed, doublereal *s, doublecomplex *b, 
	integer *ldb, doublecomplex *x, integer *ldx, doublereal *rcond, 
	doublereal *rpvgrw, doublereal *berr, integer *n_err_bnds__, 
	doublereal *err_bnds_norm__, doublereal *err_bnds_comp__, integer *
	nparams, doublereal *params, doublecomplex *work, doublereal *rwork, 
	integer *info);

/* Subroutine */ int zhetd2_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, doublereal *d__, doublereal *e, doublecomplex *tau, 
	integer *info);

/* Subroutine */ int zhetf2_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *ipiv, integer *info);

/* Subroutine */ int zhetrd_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, doublereal *d__, doublereal *e, doublecomplex *tau, 
	doublecomplex *work, integer *lwork, integer *info);

/* Subroutine */ int zhetrf_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *ipiv, doublecomplex *work, integer *lwork, 
	integer *info);

/* Subroutine */ int zhetri_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *ipiv, doublecomplex *work, integer *info);

/* Subroutine */ int zhetrs_(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, integer *ipiv, doublecomplex *b, 
	integer *ldb, integer *info);

/* Subroutine */ int zhfrk_(char *transr, char *uplo, char *trans, integer *n, 
	 integer *k, doublereal *alpha, doublecomplex *a, integer *lda, 
	doublereal *beta, doublecomplex *c__);

/* Subroutine */ int zhgeqz_(char *job, char *compq, char *compz, integer *n, 
	integer *ilo, integer *ihi, doublecomplex *h__, integer *ldh, 
	doublecomplex *t, integer *ldt, doublecomplex *alpha, doublecomplex *
	beta, doublecomplex *q, integer *ldq, doublecomplex *z__, integer *
	ldz, doublecomplex *work, integer *lwork, doublereal *rwork, integer *
	info);

/* Subroutine */ int zhpcon_(char *uplo, integer *n, doublecomplex *ap, 
	integer *ipiv, doublereal *anorm, doublereal *rcond, doublecomplex *
	work, integer *info);

/* Subroutine */ int zhpev_(char *jobz, char *uplo, integer *n, doublecomplex 
	*ap, doublereal *w, doublecomplex *z__, integer *ldz, doublecomplex *
	work, doublereal *rwork, integer *info);

/* Subroutine */ int zhpevd_(char *jobz, char *uplo, integer *n, 
	doublecomplex *ap, doublereal *w, doublecomplex *z__, integer *ldz, 
	doublecomplex *work, integer *lwork, doublereal *rwork, integer *
	lrwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int zhpevx_(char *jobz, char *range, char *uplo, integer *n, 
	doublecomplex *ap, doublereal *vl, doublereal *vu, integer *il, 
	integer *iu, doublereal *abstol, integer *m, doublereal *w, 
	doublecomplex *z__, integer *ldz, doublecomplex *work, doublereal *
	rwork, integer *iwork, integer *ifail, integer *info);

/* Subroutine */ int zhpgst_(integer *itype, char *uplo, integer *n, 
	doublecomplex *ap, doublecomplex *bp, integer *info);

/* Subroutine */ int zhpgv_(integer *itype, char *jobz, char *uplo, integer *
	n, doublecomplex *ap, doublecomplex *bp, doublereal *w, doublecomplex 
	*z__, integer *ldz, doublecomplex *work, doublereal *rwork, integer *
	info);

/* Subroutine */ int zhpgvd_(integer *itype, char *jobz, char *uplo, integer *
	n, doublecomplex *ap, doublecomplex *bp, doublereal *w, doublecomplex 
	*z__, integer *ldz, doublecomplex *work, integer *lwork, doublereal *
	rwork, integer *lrwork, integer *iwork, integer *liwork, integer *
	info);

/* Subroutine */ int zhpgvx_(integer *itype, char *jobz, char *range, char *
	uplo, integer *n, doublecomplex *ap, doublecomplex *bp, doublereal *
	vl, doublereal *vu, integer *il, integer *iu, doublereal *abstol, 
	integer *m, doublereal *w, doublecomplex *z__, integer *ldz, 
	doublecomplex *work, doublereal *rwork, integer *iwork, integer *
	ifail, integer *info);

/* Subroutine */ int zhprfs_(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *ap, doublecomplex *afp, integer *ipiv, doublecomplex *
	b, integer *ldb, doublecomplex *x, integer *ldx, doublereal *ferr, 
	doublereal *berr, doublecomplex *work, doublereal *rwork, integer *
	info);

/* Subroutine */ int zhpsv_(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *ap, integer *ipiv, doublecomplex *b, integer *ldb, 
	integer *info);

/* Subroutine */ int zhpsvx_(char *fact, char *uplo, integer *n, integer *
	nrhs, doublecomplex *ap, doublecomplex *afp, integer *ipiv, 
	doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx, 
	doublereal *rcond, doublereal *ferr, doublereal *berr, doublecomplex *
	work, doublereal *rwork, integer *info);

/* Subroutine */ int zhptrd_(char *uplo, integer *n, doublecomplex *ap, 
	doublereal *d__, doublereal *e, doublecomplex *tau, integer *info);

/* Subroutine */ int zhptrf_(char *uplo, integer *n, doublecomplex *ap, 
	integer *ipiv, integer *info);

/* Subroutine */ int zhptri_(char *uplo, integer *n, doublecomplex *ap, 
	integer *ipiv, doublecomplex *work, integer *info);

/* Subroutine */ int zhptrs_(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *ap, integer *ipiv, doublecomplex *b, integer *ldb, 
	integer *info);

/* Subroutine */ int zhsein_(char *side, char *eigsrc, char *initv, logical *
	select, integer *n, doublecomplex *h__, integer *ldh, doublecomplex *
	w, doublecomplex *vl, integer *ldvl, doublecomplex *vr, integer *ldvr, 
	 integer *mm, integer *m, doublecomplex *work, doublereal *rwork, 
	integer *ifaill, integer *ifailr, integer *info);

/* Subroutine */ int zhseqr_(char *job, char *compz, integer *n, integer *ilo, 
	 integer *ihi, doublecomplex *h__, integer *ldh, doublecomplex *w, 
	doublecomplex *z__, integer *ldz, doublecomplex *work, integer *lwork, 
	 integer *info);

/* Subroutine */ int zla_gbamv__(integer *trans, integer *m, integer *n, 
	integer *kl, integer *ku, doublereal *alpha, doublecomplex *ab, 
	integer *ldab, doublecomplex *x, integer *incx, doublereal *beta, 
	doublereal *y, integer *incy);

doublereal zla_gbrcond_c__(char *trans, integer *n, integer *kl, integer *ku, 
	doublecomplex *ab, integer *ldab, doublecomplex *afb, integer *ldafb, 
	integer *ipiv, doublereal *c__, logical *capply, integer *info, 
	doublecomplex *work, doublereal *rwork, ftnlen trans_len);

doublereal zla_gbrcond_x__(char *trans, integer *n, integer *kl, integer *ku, 
	doublecomplex *ab, integer *ldab, doublecomplex *afb, integer *ldafb, 
	integer *ipiv, doublecomplex *x, integer *info, doublecomplex *work, 
	doublereal *rwork, ftnlen trans_len);

/* Subroutine */ int zla_gbrfsx_extended__(integer *prec_type__, integer *
	trans_type__, integer *n, integer *kl, integer *ku, integer *nrhs, 
	doublecomplex *ab, integer *ldab, doublecomplex *afb, integer *ldafb, 
	integer *ipiv, logical *colequ, doublereal *c__, doublecomplex *b, 
	integer *ldb, doublecomplex *y, integer *ldy, doublereal *berr_out__, 
	integer *n_norms__, doublereal *errs_n__, doublereal *errs_c__, 
	doublecomplex *res, doublereal *ayb, doublecomplex *dy, doublecomplex 
	*y_tail__, doublereal *rcond, integer *ithresh, doublereal *rthresh, 
	doublereal *dz_ub__, logical *ignore_cwise__, integer *info);

doublereal zla_gbrpvgrw__(integer *n, integer *kl, integer *ku, integer *
	ncols, doublecomplex *ab, integer *ldab, doublecomplex *afb, integer *
	ldafb);

/* Subroutine */ int zla_geamv__(integer *trans, integer *m, integer *n, 
	doublereal *alpha, doublecomplex *a, integer *lda, doublecomplex *x, 
	integer *incx, doublereal *beta, doublereal *y, integer *incy);

doublereal zla_gercond_c__(char *trans, integer *n, doublecomplex *a, integer 
	*lda, doublecomplex *af, integer *ldaf, integer *ipiv, doublereal *
	c__, logical *capply, integer *info, doublecomplex *work, doublereal *
	rwork, ftnlen trans_len);

doublereal zla_gercond_x__(char *trans, integer *n, doublecomplex *a, integer 
	*lda, doublecomplex *af, integer *ldaf, integer *ipiv, doublecomplex *
	x, integer *info, doublecomplex *work, doublereal *rwork, ftnlen 
	trans_len);

/* Subroutine */ int zla_gerfsx_extended__(integer *prec_type__, integer *
	trans_type__, integer *n, integer *nrhs, doublecomplex *a, integer *
	lda, doublecomplex *af, integer *ldaf, integer *ipiv, logical *colequ,
	 doublereal *c__, doublecomplex *b, integer *ldb, doublecomplex *y, 
	integer *ldy, doublereal *berr_out__, integer *n_norms__, doublereal *
	errs_n__, doublereal *errs_c__, doublecomplex *res, doublereal *ayb, 
	doublecomplex *dy, doublecomplex *y_tail__, doublereal *rcond, 
	integer *ithresh, doublereal *rthresh, doublereal *dz_ub__, logical *
	ignore_cwise__, integer *info);

/* Subroutine */ int zla_heamv__(integer *uplo, integer *n, doublereal *alpha,
	 doublecomplex *a, integer *lda, doublecomplex *x, integer *incx, 
	doublereal *beta, doublereal *y, integer *incy);

doublereal zla_hercond_c__(char *uplo, integer *n, doublecomplex *a, integer *
	lda, doublecomplex *af, integer *ldaf, integer *ipiv, doublereal *c__,
	 logical *capply, integer *info, doublecomplex *work, doublereal *
	rwork, ftnlen uplo_len);

doublereal zla_hercond_x__(char *uplo, integer *n, doublecomplex *a, integer *
	lda, doublecomplex *af, integer *ldaf, integer *ipiv, doublecomplex *
	x, integer *info, doublecomplex *work, doublereal *rwork, ftnlen 
	uplo_len);

/* Subroutine */ int zla_herfsx_extended__(integer *prec_type__, char *uplo, 
	integer *n, integer *nrhs, doublecomplex *a, integer *lda, 
	doublecomplex *af, integer *ldaf, integer *ipiv, logical *colequ, 
	doublereal *c__, doublecomplex *b, integer *ldb, doublecomplex *y, 
	integer *ldy, doublereal *berr_out__, integer *n_norms__, doublereal *
	errs_n__, doublereal *errs_c__, doublecomplex *res, doublereal *ayb, 
	doublecomplex *dy, doublecomplex *y_tail__, doublereal *rcond, 
	integer *ithresh, doublereal *rthresh, doublereal *dz_ub__, logical *
	ignore_cwise__, integer *info, ftnlen uplo_len);

doublereal zla_herpvgrw__(char *uplo, integer *n, integer *info, 
	doublecomplex *a, integer *lda, doublecomplex *af, integer *ldaf, 
	integer *ipiv, doublereal *work, ftnlen uplo_len);

/* Subroutine */ int zla_lin_berr__(integer *n, integer *nz, integer *nrhs, 
	doublecomplex *res, doublereal *ayb, doublereal *berr);

doublereal zla_porcond_c__(char *uplo, integer *n, doublecomplex *a, integer *
	lda, doublecomplex *af, integer *ldaf, doublereal *c__, logical *
	capply, integer *info, doublecomplex *work, doublereal *rwork, ftnlen 
	uplo_len);

doublereal zla_porcond_x__(char *uplo, integer *n, doublecomplex *a, integer *
	lda, doublecomplex *af, integer *ldaf, doublecomplex *x, integer *
	info, doublecomplex *work, doublereal *rwork, ftnlen uplo_len);

/* Subroutine */ int zla_porfsx_extended__(integer *prec_type__, char *uplo, 
	integer *n, integer *nrhs, doublecomplex *a, integer *lda, 
	doublecomplex *af, integer *ldaf, logical *colequ, doublereal *c__, 
	doublecomplex *b, integer *ldb, doublecomplex *y, integer *ldy, 
	doublereal *berr_out__, integer *n_norms__, doublereal *errs_n__, 
	doublereal *errs_c__, doublecomplex *res, doublereal *ayb, 
	doublecomplex *dy, doublecomplex *y_tail__, doublereal *rcond, 
	integer *ithresh, doublereal *rthresh, doublereal *dz_ub__, logical *
	ignore_cwise__, integer *info, ftnlen uplo_len);

doublereal zla_porpvgrw__(char *uplo, integer *ncols, doublecomplex *a, 
	integer *lda, doublecomplex *af, integer *ldaf, doublereal *work, 
	ftnlen uplo_len);

doublereal zla_rpvgrw__(integer *n, integer *ncols, doublecomplex *a, integer 
	*lda, doublecomplex *af, integer *ldaf);

/* Subroutine */ int zla_syamv__(integer *uplo, integer *n, doublereal *alpha,
	 doublecomplex *a, integer *lda, doublecomplex *x, integer *incx, 
	doublereal *beta, doublereal *y, integer *incy);

doublereal zla_syrcond_c__(char *uplo, integer *n, doublecomplex *a, integer *
	lda, doublecomplex *af, integer *ldaf, integer *ipiv, doublereal *c__,
	 logical *capply, integer *info, doublecomplex *work, doublereal *
	rwork, ftnlen uplo_len);

doublereal zla_syrcond_x__(char *uplo, integer *n, doublecomplex *a, integer *
	lda, doublecomplex *af, integer *ldaf, integer *ipiv, doublecomplex *
	x, integer *info, doublecomplex *work, doublereal *rwork, ftnlen 
	uplo_len);

/* Subroutine */ int zla_syrfsx_extended__(integer *prec_type__, char *uplo, 
	integer *n, integer *nrhs, doublecomplex *a, integer *lda, 
	doublecomplex *af, integer *ldaf, integer *ipiv, logical *colequ, 
	doublereal *c__, doublecomplex *b, integer *ldb, doublecomplex *y, 
	integer *ldy, doublereal *berr_out__, integer *n_norms__, doublereal *
	errs_n__, doublereal *errs_c__, doublecomplex *res, doublereal *ayb, 
	doublecomplex *dy, doublecomplex *y_tail__, doublereal *rcond, 
	integer *ithresh, doublereal *rthresh, doublereal *dz_ub__, logical *
	ignore_cwise__, integer *info, ftnlen uplo_len);

doublereal zla_syrpvgrw__(char *uplo, integer *n, integer *info, 
	doublecomplex *a, integer *lda, doublecomplex *af, integer *ldaf, 
	integer *ipiv, doublereal *work, ftnlen uplo_len);

/* Subroutine */ int zla_wwaddw__(integer *n, doublecomplex *x, doublecomplex 
	*y, doublecomplex *w);

/* Subroutine */ int zlabrd_(integer *m, integer *n, integer *nb, 
	doublecomplex *a, integer *lda, doublereal *d__, doublereal *e, 
	doublecomplex *tauq, doublecomplex *taup, doublecomplex *x, integer *
	ldx, doublecomplex *y, integer *ldy);

/* Subroutine */ int zlacgv_(integer *n, doublecomplex *x, integer *incx);

/* Subroutine */ int zlacn2_(integer *n, doublecomplex *v, doublecomplex *x, 
	doublereal *est, integer *kase, integer *isave);

/* Subroutine */ int zlacon_(integer *n, doublecomplex *v, doublecomplex *x, 
	doublereal *est, integer *kase);

/* Subroutine */ int zlacp2_(char *uplo, integer *m, integer *n, doublereal *
	a, integer *lda, doublecomplex *b, integer *ldb);

/* Subroutine */ int zlacpy_(char *uplo, integer *m, integer *n, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb);

/* Subroutine */ int zlacrm_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublereal *b, integer *ldb, doublecomplex *c__, 
	integer *ldc, doublereal *rwork);

/* Subroutine */ int zlacrt_(integer *n, doublecomplex *cx, integer *incx, 
	doublecomplex *cy, integer *incy, doublecomplex *c__, doublecomplex *
	s);

/* Double Complex */ VOID zladiv_(doublecomplex * ret_val, doublecomplex *x, 
	doublecomplex *y);

/* Subroutine */ int zlaed0_(integer *qsiz, integer *n, doublereal *d__, 
	doublereal *e, doublecomplex *q, integer *ldq, doublecomplex *qstore, 
	integer *ldqs, doublereal *rwork, integer *iwork, integer *info);

/* Subroutine */ int zlaed7_(integer *n, integer *cutpnt, integer *qsiz, 
	integer *tlvls, integer *curlvl, integer *curpbm, doublereal *d__, 
	doublecomplex *q, integer *ldq, doublereal *rho, integer *indxq, 
	doublereal *qstore, integer *qptr, integer *prmptr, integer *perm, 
	integer *givptr, integer *givcol, doublereal *givnum, doublecomplex *
	work, doublereal *rwork, integer *iwork, integer *info);

/* Subroutine */ int zlaed8_(integer *k, integer *n, integer *qsiz, 
	doublecomplex *q, integer *ldq, doublereal *d__, doublereal *rho, 
	integer *cutpnt, doublereal *z__, doublereal *dlamda, doublecomplex *
	q2, integer *ldq2, doublereal *w, integer *indxp, integer *indx, 
	integer *indxq, integer *perm, integer *givptr, integer *givcol, 
	doublereal *givnum, integer *info);

/* Subroutine */ int zlaein_(logical *rightv, logical *noinit, integer *n, 
	doublecomplex *h__, integer *ldh, doublecomplex *w, doublecomplex *v, 
	doublecomplex *b, integer *ldb, doublereal *rwork, doublereal *eps3, 
	doublereal *smlnum, integer *info);

/* Subroutine */ int zlaesy_(doublecomplex *a, doublecomplex *b, 
	doublecomplex *c__, doublecomplex *rt1, doublecomplex *rt2, 
	doublecomplex *evscal, doublecomplex *cs1, doublecomplex *sn1);

/* Subroutine */ int zlaev2_(doublecomplex *a, doublecomplex *b, 
	doublecomplex *c__, doublereal *rt1, doublereal *rt2, doublereal *cs1, 
	 doublecomplex *sn1);

/* Subroutine */ int zlag2c_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, complex *sa, integer *ldsa, integer *info);

/* Subroutine */ int zlags2_(logical *upper, doublereal *a1, doublecomplex *
	a2, doublereal *a3, doublereal *b1, doublecomplex *b2, doublereal *b3, 
	 doublereal *csu, doublecomplex *snu, doublereal *csv, doublecomplex *
	snv, doublereal *csq, doublecomplex *snq);

/* Subroutine */ int zlagtm_(char *trans, integer *n, integer *nrhs, 
	doublereal *alpha, doublecomplex *dl, doublecomplex *d__, 
	doublecomplex *du, doublecomplex *x, integer *ldx, doublereal *beta, 
	doublecomplex *b, integer *ldb);

/* Subroutine */ int zlahef_(char *uplo, integer *n, integer *nb, integer *kb, 
	 doublecomplex *a, integer *lda, integer *ipiv, doublecomplex *w, 
	integer *ldw, integer *info);

/* Subroutine */ int zlahqr_(logical *wantt, logical *wantz, integer *n, 
	integer *ilo, integer *ihi, doublecomplex *h__, integer *ldh, 
	doublecomplex *w, integer *iloz, integer *ihiz, doublecomplex *z__, 
	integer *ldz, integer *info);

/* Subroutine */ int zlahr2_(integer *n, integer *k, integer *nb, 
	doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *t, 
	integer *ldt, doublecomplex *y, integer *ldy);

/* Subroutine */ int zlahrd_(integer *n, integer *k, integer *nb, 
	doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *t, 
	integer *ldt, doublecomplex *y, integer *ldy);

/* Subroutine */ int zlaic1_(integer *job, integer *j, doublecomplex *x, 
	doublereal *sest, doublecomplex *w, doublecomplex *gamma, doublereal *
	sestpr, doublecomplex *s, doublecomplex *c__);

/* Subroutine */ int zlals0_(integer *icompq, integer *nl, integer *nr, 
	integer *sqre, integer *nrhs, doublecomplex *b, integer *ldb, 
	doublecomplex *bx, integer *ldbx, integer *perm, integer *givptr, 
	integer *givcol, integer *ldgcol, doublereal *givnum, integer *ldgnum, 
	 doublereal *poles, doublereal *difl, doublereal *difr, doublereal *
	z__, integer *k, doublereal *c__, doublereal *s, doublereal *rwork, 
	integer *info);

/* Subroutine */ int zlalsa_(integer *icompq, integer *smlsiz, integer *n, 
	integer *nrhs, doublecomplex *b, integer *ldb, doublecomplex *bx, 
	integer *ldbx, doublereal *u, integer *ldu, doublereal *vt, integer *
	k, doublereal *difl, doublereal *difr, doublereal *z__, doublereal *
	poles, integer *givptr, integer *givcol, integer *ldgcol, integer *
	perm, doublereal *givnum, doublereal *c__, doublereal *s, doublereal *
	rwork, integer *iwork, integer *info);

/* Subroutine */ int zlalsd_(char *uplo, integer *smlsiz, integer *n, integer 
	*nrhs, doublereal *d__, doublereal *e, doublecomplex *b, integer *ldb, 
	 doublereal *rcond, integer *rank, doublecomplex *work, doublereal *
	rwork, integer *iwork, integer *info);

doublereal zlangb_(char *norm, integer *n, integer *kl, integer *ku, 
	doublecomplex *ab, integer *ldab, doublereal *work);

doublereal zlange_(char *norm, integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublereal *work);

doublereal zlangt_(char *norm, integer *n, doublecomplex *dl, doublecomplex *
	d__, doublecomplex *du);

doublereal zlanhb_(char *norm, char *uplo, integer *n, integer *k, 
	doublecomplex *ab, integer *ldab, doublereal *work);

doublereal zlanhe_(char *norm, char *uplo, integer *n, doublecomplex *a, 
	integer *lda, doublereal *work);

doublereal zlanhf_(char *norm, char *transr, char *uplo, integer *n, 
	doublecomplex *a, doublereal *work);

doublereal zlanhp_(char *norm, char *uplo, integer *n, doublecomplex *ap, 
	doublereal *work);

doublereal zlanhs_(char *norm, integer *n, doublecomplex *a, integer *lda, 
	doublereal *work);

doublereal zlanht_(char *norm, integer *n, doublereal *d__, doublecomplex *e);

doublereal zlansb_(char *norm, char *uplo, integer *n, integer *k, 
	doublecomplex *ab, integer *ldab, doublereal *work);

doublereal zlansp_(char *norm, char *uplo, integer *n, doublecomplex *ap, 
	doublereal *work);

doublereal zlansy_(char *norm, char *uplo, integer *n, doublecomplex *a, 
	integer *lda, doublereal *work);

doublereal zlantb_(char *norm, char *uplo, char *diag, integer *n, integer *k, 
	 doublecomplex *ab, integer *ldab, doublereal *work);

doublereal zlantp_(char *norm, char *uplo, char *diag, integer *n, 
	doublecomplex *ap, doublereal *work);

doublereal zlantr_(char *norm, char *uplo, char *diag, integer *m, integer *n, 
	 doublecomplex *a, integer *lda, doublereal *work);

/* Subroutine */ int zlapll_(integer *n, doublecomplex *x, integer *incx, 
	doublecomplex *y, integer *incy, doublereal *ssmin);

/* Subroutine */ int zlapmt_(logical *forwrd, integer *m, integer *n, 
	doublecomplex *x, integer *ldx, integer *k);

/* Subroutine */ int zlaqgb_(integer *m, integer *n, integer *kl, integer *ku, 
	 doublecomplex *ab, integer *ldab, doublereal *r__, doublereal *c__, 
	doublereal *rowcnd, doublereal *colcnd, doublereal *amax, char *equed);

/* Subroutine */ int zlaqge_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublereal *r__, doublereal *c__, doublereal *rowcnd, 
	doublereal *colcnd, doublereal *amax, char *equed);

/* Subroutine */ int zlaqhb_(char *uplo, integer *n, integer *kd, 
	doublecomplex *ab, integer *ldab, doublereal *s, doublereal *scond, 
	doublereal *amax, char *equed);

/* Subroutine */ int zlaqhe_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, doublereal *s, doublereal *scond, doublereal *amax, 
	char *equed);

/* Subroutine */ int zlaqhp_(char *uplo, integer *n, doublecomplex *ap, 
	doublereal *s, doublereal *scond, doublereal *amax, char *equed);

/* Subroutine */ int zlaqp2_(integer *m, integer *n, integer *offset, 
	doublecomplex *a, integer *lda, integer *jpvt, doublecomplex *tau, 
	doublereal *vn1, doublereal *vn2, doublecomplex *work);

/* Subroutine */ int zlaqps_(integer *m, integer *n, integer *offset, integer 
	*nb, integer *kb, doublecomplex *a, integer *lda, integer *jpvt, 
	doublecomplex *tau, doublereal *vn1, doublereal *vn2, doublecomplex *
	auxv, doublecomplex *f, integer *ldf);

/* Subroutine */ int zlaqr0_(logical *wantt, logical *wantz, integer *n, 
	integer *ilo, integer *ihi, doublecomplex *h__, integer *ldh, 
	doublecomplex *w, integer *iloz, integer *ihiz, doublecomplex *z__, 
	integer *ldz, doublecomplex *work, integer *lwork, integer *info);

/* Subroutine */ int zlaqr1_(integer *n, doublecomplex *h__, integer *ldh, 
	doublecomplex *s1, doublecomplex *s2, doublecomplex *v);

/* Subroutine */ int zlaqr2_(logical *wantt, logical *wantz, integer *n, 
	integer *ktop, integer *kbot, integer *nw, doublecomplex *h__, 
	integer *ldh, integer *iloz, integer *ihiz, doublecomplex *z__, 
	integer *ldz, integer *ns, integer *nd, doublecomplex *sh, 
	doublecomplex *v, integer *ldv, integer *nh, doublecomplex *t, 
	integer *ldt, integer *nv, doublecomplex *wv, integer *ldwv, 
	doublecomplex *work, integer *lwork);

/* Subroutine */ int zlaqr3_(logical *wantt, logical *wantz, integer *n, 
	integer *ktop, integer *kbot, integer *nw, doublecomplex *h__, 
	integer *ldh, integer *iloz, integer *ihiz, doublecomplex *z__, 
	integer *ldz, integer *ns, integer *nd, doublecomplex *sh, 
	doublecomplex *v, integer *ldv, integer *nh, doublecomplex *t, 
	integer *ldt, integer *nv, doublecomplex *wv, integer *ldwv, 
	doublecomplex *work, integer *lwork);

/* Subroutine */ int zlaqr4_(logical *wantt, logical *wantz, integer *n, 
	integer *ilo, integer *ihi, doublecomplex *h__, integer *ldh, 
	doublecomplex *w, integer *iloz, integer *ihiz, doublecomplex *z__, 
	integer *ldz, doublecomplex *work, integer *lwork, integer *info);

/* Subroutine */ int zlaqr5_(logical *wantt, logical *wantz, integer *kacc22, 
	integer *n, integer *ktop, integer *kbot, integer *nshfts, 
	doublecomplex *s, doublecomplex *h__, integer *ldh, integer *iloz, 
	integer *ihiz, doublecomplex *z__, integer *ldz, doublecomplex *v, 
	integer *ldv, doublecomplex *u, integer *ldu, integer *nv, 
	doublecomplex *wv, integer *ldwv, integer *nh, doublecomplex *wh, 
	integer *ldwh);

/* Subroutine */ int zlaqsb_(char *uplo, integer *n, integer *kd, 
	doublecomplex *ab, integer *ldab, doublereal *s, doublereal *scond, 
	doublereal *amax, char *equed);

/* Subroutine */ int zlaqsp_(char *uplo, integer *n, doublecomplex *ap, 
	doublereal *s, doublereal *scond, doublereal *amax, char *equed);

/* Subroutine */ int zlaqsy_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, doublereal *s, doublereal *scond, doublereal *amax, 
	char *equed);

/* Subroutine */ int zlar1v_(integer *n, integer *b1, integer *bn, doublereal 
	*lambda, doublereal *d__, doublereal *l, doublereal *ld, doublereal *
	lld, doublereal *pivmin, doublereal *gaptol, doublecomplex *z__, 
	logical *wantnc, integer *negcnt, doublereal *ztz, doublereal *mingma, 
	 integer *r__, integer *isuppz, doublereal *nrminv, doublereal *resid, 
	 doublereal *rqcorr, doublereal *work);

/* Subroutine */ int zlar2v_(integer *n, doublecomplex *x, doublecomplex *y, 
	doublecomplex *z__, integer *incx, doublereal *c__, doublecomplex *s, 
	integer *incc);

/* Subroutine */ int zlarcm_(integer *m, integer *n, doublereal *a, integer *
	lda, doublecomplex *b, integer *ldb, doublecomplex *c__, integer *ldc, 
	 doublereal *rwork);

/* Subroutine */ int zlarf_(char *side, integer *m, integer *n, doublecomplex 
	*v, integer *incv, doublecomplex *tau, doublecomplex *c__, integer *
	ldc, doublecomplex *work);

/* Subroutine */ int zlarfb_(char *side, char *trans, char *direct, char *
	storev, integer *m, integer *n, integer *k, doublecomplex *v, integer 
	*ldv, doublecomplex *t, integer *ldt, doublecomplex *c__, integer *
	ldc, doublecomplex *work, integer *ldwork);

/* Subroutine */ int zlarfg_(integer *n, doublecomplex *alpha, doublecomplex *
	x, integer *incx, doublecomplex *tau);

/* Subroutine */ int zlarfp_(integer *n, doublecomplex *alpha, doublecomplex *
	x, integer *incx, doublecomplex *tau);

/* Subroutine */ int zlarft_(char *direct, char *storev, integer *n, integer *
	k, doublecomplex *v, integer *ldv, doublecomplex *tau, doublecomplex *
	t, integer *ldt);

/* Subroutine */ int zlarfx_(char *side, integer *m, integer *n, 
	doublecomplex *v, doublecomplex *tau, doublecomplex *c__, integer *
	ldc, doublecomplex *work);

/* Subroutine */ int zlargv_(integer *n, doublecomplex *x, integer *incx, 
	doublecomplex *y, integer *incy, doublereal *c__, integer *incc);

/* Subroutine */ int zlarnv_(integer *idist, integer *iseed, integer *n, 
	doublecomplex *x);

/* Subroutine */ int zlarrv_(integer *n, doublereal *vl, doublereal *vu, 
	doublereal *d__, doublereal *l, doublereal *pivmin, integer *isplit, 
	integer *m, integer *dol, integer *dou, doublereal *minrgp, 
	doublereal *rtol1, doublereal *rtol2, doublereal *w, doublereal *werr, 
	 doublereal *wgap, integer *iblock, integer *indexw, doublereal *gers, 
	 doublecomplex *z__, integer *ldz, integer *isuppz, doublereal *work, 
	integer *iwork, integer *info);

/* Subroutine */ int zlarscl2_(integer *m, integer *n, doublereal *d__, 
	doublecomplex *x, integer *ldx);

/* Subroutine */ int zlartg_(doublecomplex *f, doublecomplex *g, doublereal *
	cs, doublecomplex *sn, doublecomplex *r__);

/* Subroutine */ int zlartv_(integer *n, doublecomplex *x, integer *incx, 
	doublecomplex *y, integer *incy, doublereal *c__, doublecomplex *s, 
	integer *incc);

/* Subroutine */ int zlarz_(char *side, integer *m, integer *n, integer *l, 
	doublecomplex *v, integer *incv, doublecomplex *tau, doublecomplex *
	c__, integer *ldc, doublecomplex *work);

/* Subroutine */ int zlarzb_(char *side, char *trans, char *direct, char *
	storev, integer *m, integer *n, integer *k, integer *l, doublecomplex 
	*v, integer *ldv, doublecomplex *t, integer *ldt, doublecomplex *c__, 
	integer *ldc, doublecomplex *work, integer *ldwork);

/* Subroutine */ int zlarzt_(char *direct, char *storev, integer *n, integer *
	k, doublecomplex *v, integer *ldv, doublecomplex *tau, doublecomplex *
	t, integer *ldt);

/* Subroutine */ int zlascl_(char *type__, integer *kl, integer *ku, 
	doublereal *cfrom, doublereal *cto, integer *m, integer *n, 
	doublecomplex *a, integer *lda, integer *info);

/* Subroutine */ int zlascl2_(integer *m, integer *n, doublereal *d__, 
	doublecomplex *x, integer *ldx);

/* Subroutine */ int zlaset_(char *uplo, integer *m, integer *n, 
	doublecomplex *alpha, doublecomplex *beta, doublecomplex *a, integer *
	lda);

/* Subroutine */ int zlasr_(char *side, char *pivot, char *direct, integer *m, 
	 integer *n, doublereal *c__, doublereal *s, doublecomplex *a, 
	integer *lda);

/* Subroutine */ int zlassq_(integer *n, doublecomplex *x, integer *incx, 
	doublereal *scale, doublereal *sumsq);

/* Subroutine */ int zlaswp_(integer *n, doublecomplex *a, integer *lda, 
	integer *k1, integer *k2, integer *ipiv, integer *incx);

/* Subroutine */ int zlasyf_(char *uplo, integer *n, integer *nb, integer *kb, 
	 doublecomplex *a, integer *lda, integer *ipiv, doublecomplex *w, 
	integer *ldw, integer *info);

/* Subroutine */ int zlat2c_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, complex *sa, integer *ldsa, integer *info);

/* Subroutine */ int zlatbs_(char *uplo, char *trans, char *diag, char *
	normin, integer *n, integer *kd, doublecomplex *ab, integer *ldab, 
	doublecomplex *x, doublereal *scale, doublereal *cnorm, integer *info);

/* Subroutine */ int zlatdf_(integer *ijob, integer *n, doublecomplex *z__, 
	integer *ldz, doublecomplex *rhs, doublereal *rdsum, doublereal *
	rdscal, integer *ipiv, integer *jpiv);

/* Subroutine */ int zlatps_(char *uplo, char *trans, char *diag, char *
	normin, integer *n, doublecomplex *ap, doublecomplex *x, doublereal *
	scale, doublereal *cnorm, integer *info);

/* Subroutine */ int zlatrd_(char *uplo, integer *n, integer *nb, 
	doublecomplex *a, integer *lda, doublereal *e, doublecomplex *tau, 
	doublecomplex *w, integer *ldw);

/* Subroutine */ int zlatrs_(char *uplo, char *trans, char *diag, char *
	normin, integer *n, doublecomplex *a, integer *lda, doublecomplex *x, 
	doublereal *scale, doublereal *cnorm, integer *info);

/* Subroutine */ int zlatrz_(integer *m, integer *n, integer *l, 
	doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *
	work);

/* Subroutine */ int zlatzm_(char *side, integer *m, integer *n, 
	doublecomplex *v, integer *incv, doublecomplex *tau, doublecomplex *
	c1, doublecomplex *c2, integer *ldc, doublecomplex *work);

/* Subroutine */ int zlauu2_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *info);

/* Subroutine */ int zlauum_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *info);

/* Subroutine */ int zpbcon_(char *uplo, integer *n, integer *kd, 
	doublecomplex *ab, integer *ldab, doublereal *anorm, doublereal *
	rcond, doublecomplex *work, doublereal *rwork, integer *info);

/* Subroutine */ int zpbequ_(char *uplo, integer *n, integer *kd, 
	doublecomplex *ab, integer *ldab, doublereal *s, doublereal *scond, 
	doublereal *amax, integer *info);

/* Subroutine */ int zpbrfs_(char *uplo, integer *n, integer *kd, integer *
	nrhs, doublecomplex *ab, integer *ldab, doublecomplex *afb, integer *
	ldafb, doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx, 
	 doublereal *ferr, doublereal *berr, doublecomplex *work, doublereal *
	rwork, integer *info);

/* Subroutine */ int zpbstf_(char *uplo, integer *n, integer *kd, 
	doublecomplex *ab, integer *ldab, integer *info);

/* Subroutine */ int zpbsv_(char *uplo, integer *n, integer *kd, integer *
	nrhs, doublecomplex *ab, integer *ldab, doublecomplex *b, integer *
	ldb, integer *info);

/* Subroutine */ int zpbsvx_(char *fact, char *uplo, integer *n, integer *kd, 
	integer *nrhs, doublecomplex *ab, integer *ldab, doublecomplex *afb, 
	integer *ldafb, char *equed, doublereal *s, doublecomplex *b, integer 
	*ldb, doublecomplex *x, integer *ldx, doublereal *rcond, doublereal *
	ferr, doublereal *berr, doublecomplex *work, doublereal *rwork, 
	integer *info);

/* Subroutine */ int zpbtf2_(char *uplo, integer *n, integer *kd, 
	doublecomplex *ab, integer *ldab, integer *info);

/* Subroutine */ int zpbtrf_(char *uplo, integer *n, integer *kd, 
	doublecomplex *ab, integer *ldab, integer *info);

/* Subroutine */ int zpbtrs_(char *uplo, integer *n, integer *kd, integer *
	nrhs, doublecomplex *ab, integer *ldab, doublecomplex *b, integer *
	ldb, integer *info);

/* Subroutine */ int zpftrf_(char *transr, char *uplo, integer *n, 
	doublecomplex *a, integer *info);

/* Subroutine */ int zpftri_(char *transr, char *uplo, integer *n, 
	doublecomplex *a, integer *info);

/* Subroutine */ int zpftrs_(char *transr, char *uplo, integer *n, integer *
	nrhs, doublecomplex *a, doublecomplex *b, integer *ldb, integer *info);

/* Subroutine */ int zpocon_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, doublereal *anorm, doublereal *rcond, doublecomplex *
	work, doublereal *rwork, integer *info);

/* Subroutine */ int zpoequ_(integer *n, doublecomplex *a, integer *lda, 
	doublereal *s, doublereal *scond, doublereal *amax, integer *info);

/* Subroutine */ int zpoequb_(integer *n, doublecomplex *a, integer *lda, 
	doublereal *s, doublereal *scond, doublereal *amax, integer *info);

/* Subroutine */ int zporfs_(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, doublecomplex *af, integer *ldaf, 
	doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx, 
	doublereal *ferr, doublereal *berr, doublecomplex *work, doublereal *
	rwork, integer *info);

/* Subroutine */ int zporfsx_(char *uplo, char *equed, integer *n, integer *
	nrhs, doublecomplex *a, integer *lda, doublecomplex *af, integer *
	ldaf, doublereal *s, doublecomplex *b, integer *ldb, doublecomplex *x, 
	 integer *ldx, doublereal *rcond, doublereal *berr, integer *
	n_err_bnds__, doublereal *err_bnds_norm__, doublereal *
	err_bnds_comp__, integer *nparams, doublereal *params, doublecomplex *
	work, doublereal *rwork, integer *info);

/* Subroutine */ int zposv_(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	integer *info);

/* Subroutine */ int zposvx_(char *fact, char *uplo, integer *n, integer *
	nrhs, doublecomplex *a, integer *lda, doublecomplex *af, integer *
	ldaf, char *equed, doublereal *s, doublecomplex *b, integer *ldb, 
	doublecomplex *x, integer *ldx, doublereal *rcond, doublereal *ferr, 
	doublereal *berr, doublecomplex *work, doublereal *rwork, integer *
	info);

/* Subroutine */ int zposvxx_(char *fact, char *uplo, integer *n, integer *
	nrhs, doublecomplex *a, integer *lda, doublecomplex *af, integer *
	ldaf, char *equed, doublereal *s, doublecomplex *b, integer *ldb, 
	doublecomplex *x, integer *ldx, doublereal *rcond, doublereal *rpvgrw, 
	 doublereal *berr, integer *n_err_bnds__, doublereal *err_bnds_norm__, 
	 doublereal *err_bnds_comp__, integer *nparams, doublereal *params, 
	doublecomplex *work, doublereal *rwork, integer *info);

/* Subroutine */ int zpotf2_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *info);

/* Subroutine */ int zpotrf_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *info);

/* Subroutine */ int zpotri_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *info);

/* Subroutine */ int zpotrs_(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	integer *info);

/* Subroutine */ int zppcon_(char *uplo, integer *n, doublecomplex *ap, 
	doublereal *anorm, doublereal *rcond, doublecomplex *work, doublereal 
	*rwork, integer *info);

/* Subroutine */ int zppequ_(char *uplo, integer *n, doublecomplex *ap, 
	doublereal *s, doublereal *scond, doublereal *amax, integer *info);

/* Subroutine */ int zpprfs_(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *ap, doublecomplex *afp, doublecomplex *b, integer *ldb, 
	 doublecomplex *x, integer *ldx, doublereal *ferr, doublereal *berr, 
	doublecomplex *work, doublereal *rwork, integer *info);

/* Subroutine */ int zppsv_(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *ap, doublecomplex *b, integer *ldb, integer *info);

/* Subroutine */ int zppsvx_(char *fact, char *uplo, integer *n, integer *
	nrhs, doublecomplex *ap, doublecomplex *afp, char *equed, doublereal *
	s, doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx, 
	doublereal *rcond, doublereal *ferr, doublereal *berr, doublecomplex *
	work, doublereal *rwork, integer *info);

/* Subroutine */ int zpptrf_(char *uplo, integer *n, doublecomplex *ap, 
	integer *info);

/* Subroutine */ int zpptri_(char *uplo, integer *n, doublecomplex *ap, 
	integer *info);

/* Subroutine */ int zpptrs_(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *ap, doublecomplex *b, integer *ldb, integer *info);

/* Subroutine */ int zpstf2_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *piv, integer *rank, doublereal *tol, 
	doublereal *work, integer *info);

/* Subroutine */ int zpstrf_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *piv, integer *rank, doublereal *tol, 
	doublereal *work, integer *info);

/* Subroutine */ int zptcon_(integer *n, doublereal *d__, doublecomplex *e, 
	doublereal *anorm, doublereal *rcond, doublereal *rwork, integer *
	info);

/* Subroutine */ int zpteqr_(char *compz, integer *n, doublereal *d__, 
	doublereal *e, doublecomplex *z__, integer *ldz, doublereal *work, 
	integer *info);

/* Subroutine */ int zptrfs_(char *uplo, integer *n, integer *nrhs, 
	doublereal *d__, doublecomplex *e, doublereal *df, doublecomplex *ef, 
	doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx, 
	doublereal *ferr, doublereal *berr, doublecomplex *work, doublereal *
	rwork, integer *info);

/* Subroutine */ int zptsv_(integer *n, integer *nrhs, doublereal *d__, 
	doublecomplex *e, doublecomplex *b, integer *ldb, integer *info);

/* Subroutine */ int zptsvx_(char *fact, integer *n, integer *nrhs, 
	doublereal *d__, doublecomplex *e, doublereal *df, doublecomplex *ef, 
	doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx, 
	doublereal *rcond, doublereal *ferr, doublereal *berr, doublecomplex *
	work, doublereal *rwork, integer *info);

/* Subroutine */ int zpttrf_(integer *n, doublereal *d__, doublecomplex *e, 
	integer *info);

/* Subroutine */ int zpttrs_(char *uplo, integer *n, integer *nrhs, 
	doublereal *d__, doublecomplex *e, doublecomplex *b, integer *ldb, 
	integer *info);

/* Subroutine */ int zptts2_(integer *iuplo, integer *n, integer *nrhs, 
	doublereal *d__, doublecomplex *e, doublecomplex *b, integer *ldb);

/* Subroutine */ int zrot_(integer *n, doublecomplex *cx, integer *incx, 
	doublecomplex *cy, integer *incy, doublereal *c__, doublecomplex *s);

/* Subroutine */ int zspcon_(char *uplo, integer *n, doublecomplex *ap, 
	integer *ipiv, doublereal *anorm, doublereal *rcond, doublecomplex *
	work, integer *info);

/* Subroutine */ int zspmv_(char *uplo, integer *n, doublecomplex *alpha, 
	doublecomplex *ap, doublecomplex *x, integer *incx, doublecomplex *
	beta, doublecomplex *y, integer *incy);

/* Subroutine */ int zspr_(char *uplo, integer *n, doublecomplex *alpha, 
	doublecomplex *x, integer *incx, doublecomplex *ap);

/* Subroutine */ int zsprfs_(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *ap, doublecomplex *afp, integer *ipiv, doublecomplex *
	b, integer *ldb, doublecomplex *x, integer *ldx, doublereal *ferr, 
	doublereal *berr, doublecomplex *work, doublereal *rwork, integer *
	info);

/* Subroutine */ int zspsv_(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *ap, integer *ipiv, doublecomplex *b, integer *ldb, 
	integer *info);

/* Subroutine */ int zspsvx_(char *fact, char *uplo, integer *n, integer *
	nrhs, doublecomplex *ap, doublecomplex *afp, integer *ipiv, 
	doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx, 
	doublereal *rcond, doublereal *ferr, doublereal *berr, doublecomplex *
	work, doublereal *rwork, integer *info);

/* Subroutine */ int zsptrf_(char *uplo, integer *n, doublecomplex *ap, 
	integer *ipiv, integer *info);

/* Subroutine */ int zsptri_(char *uplo, integer *n, doublecomplex *ap, 
	integer *ipiv, doublecomplex *work, integer *info);

/* Subroutine */ int zsptrs_(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *ap, integer *ipiv, doublecomplex *b, integer *ldb, 
	integer *info);

/* Subroutine */ int zstedc_(char *compz, integer *n, doublereal *d__, 
	doublereal *e, doublecomplex *z__, integer *ldz, doublecomplex *work, 
	integer *lwork, doublereal *rwork, integer *lrwork, integer *iwork, 
	integer *liwork, integer *info);

/* Subroutine */ int zstegr_(char *jobz, char *range, integer *n, doublereal *
	d__, doublereal *e, doublereal *vl, doublereal *vu, integer *il, 
	integer *iu, doublereal *abstol, integer *m, doublereal *w, 
	doublecomplex *z__, integer *ldz, integer *isuppz, doublereal *work, 
	integer *lwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int zstein_(integer *n, doublereal *d__, doublereal *e, 
	integer *m, doublereal *w, integer *iblock, integer *isplit, 
	doublecomplex *z__, integer *ldz, doublereal *work, integer *iwork, 
	integer *ifail, integer *info);

/* Subroutine */ int zstemr_(char *jobz, char *range, integer *n, doublereal *
	d__, doublereal *e, doublereal *vl, doublereal *vu, integer *il, 
	integer *iu, integer *m, doublereal *w, doublecomplex *z__, integer *
	ldz, integer *nzc, integer *isuppz, logical *tryrac, doublereal *work, 
	 integer *lwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ int zsteqr_(char *compz, integer *n, doublereal *d__, 
	doublereal *e, doublecomplex *z__, integer *ldz, doublereal *work, 
	integer *info);

/* Subroutine */ int zsycon_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *ipiv, doublereal *anorm, doublereal *rcond, 
	doublecomplex *work, integer *info);

/* Subroutine */ int zsyequb_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, doublereal *s, doublereal *scond, doublereal *amax, 
	doublecomplex *work, integer *info);

/* Subroutine */ int zsymv_(char *uplo, integer *n, doublecomplex *alpha, 
	doublecomplex *a, integer *lda, doublecomplex *x, integer *incx, 
	doublecomplex *beta, doublecomplex *y, integer *incy);

/* Subroutine */ int zsyr_(char *uplo, integer *n, doublecomplex *alpha, 
	doublecomplex *x, integer *incx, doublecomplex *a, integer *lda);

/* Subroutine */ int zsyrfs_(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, doublecomplex *af, integer *ldaf, 
	integer *ipiv, doublecomplex *b, integer *ldb, doublecomplex *x, 
	integer *ldx, doublereal *ferr, doublereal *berr, doublecomplex *work, 
	 doublereal *rwork, integer *info);

/* Subroutine */ int zsyrfsx_(char *uplo, char *equed, integer *n, integer *
	nrhs, doublecomplex *a, integer *lda, doublecomplex *af, integer *
	ldaf, integer *ipiv, doublereal *s, doublecomplex *b, integer *ldb, 
	doublecomplex *x, integer *ldx, doublereal *rcond, doublereal *berr, 
	integer *n_err_bnds__, doublereal *err_bnds_norm__, doublereal *
	err_bnds_comp__, integer *nparams, doublereal *params, doublecomplex *
	work, doublereal *rwork, integer *info);

/* Subroutine */ int zsysv_(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, integer *ipiv, doublecomplex *b, 
	integer *ldb, doublecomplex *work, integer *lwork, integer *info);

/* Subroutine */ int zsysvx_(char *fact, char *uplo, integer *n, integer *
	nrhs, doublecomplex *a, integer *lda, doublecomplex *af, integer *
	ldaf, integer *ipiv, doublecomplex *b, integer *ldb, doublecomplex *x, 
	 integer *ldx, doublereal *rcond, doublereal *ferr, doublereal *berr, 
	doublecomplex *work, integer *lwork, doublereal *rwork, integer *info);

/* Subroutine */ int zsysvxx_(char *fact, char *uplo, integer *n, integer *
	nrhs, doublecomplex *a, integer *lda, doublecomplex *af, integer *
	ldaf, integer *ipiv, char *equed, doublereal *s, doublecomplex *b, 
	integer *ldb, doublecomplex *x, integer *ldx, doublereal *rcond, 
	doublereal *rpvgrw, doublereal *berr, integer *n_err_bnds__, 
	doublereal *err_bnds_norm__, doublereal *err_bnds_comp__, integer *
	nparams, doublereal *params, doublecomplex *work, doublereal *rwork, 
	integer *info);

/* Subroutine */ int zsytf2_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *ipiv, integer *info);

/* Subroutine */ int zsytrf_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *ipiv, doublecomplex *work, integer *lwork, 
	integer *info);

/* Subroutine */ int zsytri_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, integer *ipiv, doublecomplex *work, integer *info);

/* Subroutine */ int zsytrs_(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, integer *ipiv, doublecomplex *b, 
	integer *ldb, integer *info);

/* Subroutine */ int ztbcon_(char *norm, char *uplo, char *diag, integer *n, 
	integer *kd, doublecomplex *ab, integer *ldab, doublereal *rcond, 
	doublecomplex *work, doublereal *rwork, integer *info);

/* Subroutine */ int ztbrfs_(char *uplo, char *trans, char *diag, integer *n, 
	integer *kd, integer *nrhs, doublecomplex *ab, integer *ldab, 
	doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx, 
	doublereal *ferr, doublereal *berr, doublecomplex *work, doublereal *
	rwork, integer *info);

/* Subroutine */ int ztbtrs_(char *uplo, char *trans, char *diag, integer *n, 
	integer *kd, integer *nrhs, doublecomplex *ab, integer *ldab, 
	doublecomplex *b, integer *ldb, integer *info);

/* Subroutine */ int ztfsm_(char *transr, char *side, char *uplo, char *trans, 
	 char *diag, integer *m, integer *n, doublecomplex *alpha, 
	doublecomplex *a, doublecomplex *b, integer *ldb);

/* Subroutine */ int ztftri_(char *transr, char *uplo, char *diag, integer *n, 
	 doublecomplex *a, integer *info);

/* Subroutine */ int ztfttp_(char *transr, char *uplo, integer *n, 
	doublecomplex *arf, doublecomplex *ap, integer *info);

/* Subroutine */ int ztfttr_(char *transr, char *uplo, integer *n, 
	doublecomplex *arf, doublecomplex *a, integer *lda, integer *info);

/* Subroutine */ int ztgevc_(char *side, char *howmny, logical *select, 
	integer *n, doublecomplex *s, integer *lds, doublecomplex *p, integer 
	*ldp, doublecomplex *vl, integer *ldvl, doublecomplex *vr, integer *
	ldvr, integer *mm, integer *m, doublecomplex *work, doublereal *rwork, 
	 integer *info);

/* Subroutine */ int ztgex2_(logical *wantq, logical *wantz, integer *n, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	doublecomplex *q, integer *ldq, doublecomplex *z__, integer *ldz, 
	integer *j1, integer *info);

/* Subroutine */ int ztgexc_(logical *wantq, logical *wantz, integer *n, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	doublecomplex *q, integer *ldq, doublecomplex *z__, integer *ldz, 
	integer *ifst, integer *ilst, integer *info);

/* Subroutine */ int ztgsen_(integer *ijob, logical *wantq, logical *wantz, 
	logical *select, integer *n, doublecomplex *a, integer *lda, 
	doublecomplex *b, integer *ldb, doublecomplex *alpha, doublecomplex *
	beta, doublecomplex *q, integer *ldq, doublecomplex *z__, integer *
	ldz, integer *m, doublereal *pl, doublereal *pr, doublereal *dif, 
	doublecomplex *work, integer *lwork, integer *iwork, integer *liwork, 
	integer *info);

/* Subroutine */ int ztgsja_(char *jobu, char *jobv, char *jobq, integer *m, 
	integer *p, integer *n, integer *k, integer *l, doublecomplex *a, 
	integer *lda, doublecomplex *b, integer *ldb, doublereal *tola, 
	doublereal *tolb, doublereal *alpha, doublereal *beta, doublecomplex *
	u, integer *ldu, doublecomplex *v, integer *ldv, doublecomplex *q, 
	integer *ldq, doublecomplex *work, integer *ncycle, integer *info);

/* Subroutine */ int ztgsna_(char *job, char *howmny, logical *select, 
	integer *n, doublecomplex *a, integer *lda, doublecomplex *b, integer 
	*ldb, doublecomplex *vl, integer *ldvl, doublecomplex *vr, integer *
	ldvr, doublereal *s, doublereal *dif, integer *mm, integer *m, 
	doublecomplex *work, integer *lwork, integer *iwork, integer *info);

/* Subroutine */ int ztgsy2_(char *trans, integer *ijob, integer *m, integer *
	n, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	doublecomplex *c__, integer *ldc, doublecomplex *d__, integer *ldd, 
	doublecomplex *e, integer *lde, doublecomplex *f, integer *ldf, 
	doublereal *scale, doublereal *rdsum, doublereal *rdscal, integer *
	info);

/* Subroutine */ int ztgsyl_(char *trans, integer *ijob, integer *m, integer *
	n, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	doublecomplex *c__, integer *ldc, doublecomplex *d__, integer *ldd, 
	doublecomplex *e, integer *lde, doublecomplex *f, integer *ldf, 
	doublereal *scale, doublereal *dif, doublecomplex *work, integer *
	lwork, integer *iwork, integer *info);

/* Subroutine */ int ztpcon_(char *norm, char *uplo, char *diag, integer *n, 
	doublecomplex *ap, doublereal *rcond, doublecomplex *work, doublereal 
	*rwork, integer *info);

/* Subroutine */ int ztprfs_(char *uplo, char *trans, char *diag, integer *n, 
	integer *nrhs, doublecomplex *ap, doublecomplex *b, integer *ldb, 
	doublecomplex *x, integer *ldx, doublereal *ferr, doublereal *berr, 
	doublecomplex *work, doublereal *rwork, integer *info);

/* Subroutine */ int ztptri_(char *uplo, char *diag, integer *n, 
	doublecomplex *ap, integer *info);

/* Subroutine */ int ztptrs_(char *uplo, char *trans, char *diag, integer *n, 
	integer *nrhs, doublecomplex *ap, doublecomplex *b, integer *ldb, 
	integer *info);

/* Subroutine */ int ztpttf_(char *transr, char *uplo, integer *n, 
	doublecomplex *ap, doublecomplex *arf, integer *info);

/* Subroutine */ int ztpttr_(char *uplo, integer *n, doublecomplex *ap, 
	doublecomplex *a, integer *lda, integer *info);

/* Subroutine */ int ztrcon_(char *norm, char *uplo, char *diag, integer *n, 
	doublecomplex *a, integer *lda, doublereal *rcond, doublecomplex *
	work, doublereal *rwork, integer *info);

/* Subroutine */ int ztrevc_(char *side, char *howmny, logical *select, 
	integer *n, doublecomplex *t, integer *ldt, doublecomplex *vl, 
	integer *ldvl, doublecomplex *vr, integer *ldvr, integer *mm, integer 
	*m, doublecomplex *work, doublereal *rwork, integer *info);

/* Subroutine */ int ztrexc_(char *compq, integer *n, doublecomplex *t, 
	integer *ldt, doublecomplex *q, integer *ldq, integer *ifst, integer *
	ilst, integer *info);

/* Subroutine */ int ztrrfs_(char *uplo, char *trans, char *diag, integer *n, 
	integer *nrhs, doublecomplex *a, integer *lda, doublecomplex *b, 
	integer *ldb, doublecomplex *x, integer *ldx, doublereal *ferr, 
	doublereal *berr, doublecomplex *work, doublereal *rwork, integer *
	info);

/* Subroutine */ int ztrsen_(char *job, char *compq, logical *select, integer 
	*n, doublecomplex *t, integer *ldt, doublecomplex *q, integer *ldq, 
	doublecomplex *w, integer *m, doublereal *s, doublereal *sep, 
	doublecomplex *work, integer *lwork, integer *info);

/* Subroutine */ int ztrsna_(char *job, char *howmny, logical *select, 
	integer *n, doublecomplex *t, integer *ldt, doublecomplex *vl, 
	integer *ldvl, doublecomplex *vr, integer *ldvr, doublereal *s, 
	doublereal *sep, integer *mm, integer *m, doublecomplex *work, 
	integer *ldwork, doublereal *rwork, integer *info);

/* Subroutine */ int ztrsyl_(char *trana, char *tranb, integer *isgn, integer 
	*m, integer *n, doublecomplex *a, integer *lda, doublecomplex *b, 
	integer *ldb, doublecomplex *c__, integer *ldc, doublereal *scale, 
	integer *info);

/* Subroutine */ int ztrti2_(char *uplo, char *diag, integer *n, 
	doublecomplex *a, integer *lda, integer *info);

/* Subroutine */ int ztrtri_(char *uplo, char *diag, integer *n, 
	doublecomplex *a, integer *lda, integer *info);

/* Subroutine */ int ztrtrs_(char *uplo, char *trans, char *diag, integer *n, 
	integer *nrhs, doublecomplex *a, integer *lda, doublecomplex *b, 
	integer *ldb, integer *info);

/* Subroutine */ int ztrttf_(char *transr, char *uplo, integer *n, 
	doublecomplex *a, integer *lda, doublecomplex *arf, integer *info);

/* Subroutine */ int ztrttp_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, doublecomplex *ap, integer *info);

/* Subroutine */ int ztzrqf_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublecomplex *tau, integer *info);

/* Subroutine */ int ztzrzf_(integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublecomplex *tau, doublecomplex *work, integer *lwork, 
	 integer *info);

/* Subroutine */ int zung2l_(integer *m, integer *n, integer *k, 
	doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *
	work, integer *info);

/* Subroutine */ int zung2r_(integer *m, integer *n, integer *k, 
	doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *
	work, integer *info);

/* Subroutine */ int zungbr_(char *vect, integer *m, integer *n, integer *k, 
	doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *
	work, integer *lwork, integer *info);

/* Subroutine */ int zunghr_(integer *n, integer *ilo, integer *ihi, 
	doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *
	work, integer *lwork, integer *info);

/* Subroutine */ int zungl2_(integer *m, integer *n, integer *k, 
	doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *
	work, integer *info);

/* Subroutine */ int zunglq_(integer *m, integer *n, integer *k, 
	doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *
	work, integer *lwork, integer *info);

/* Subroutine */ int zungql_(integer *m, integer *n, integer *k, 
	doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *
	work, integer *lwork, integer *info);

/* Subroutine */ int zungqr_(integer *m, integer *n, integer *k, 
	doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *
	work, integer *lwork, integer *info);

/* Subroutine */ int zungr2_(integer *m, integer *n, integer *k, 
	doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *
	work, integer *info);

/* Subroutine */ int zungrq_(integer *m, integer *n, integer *k, 
	doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *
	work, integer *lwork, integer *info);

/* Subroutine */ int zungtr_(char *uplo, integer *n, doublecomplex *a, 
	integer *lda, doublecomplex *tau, doublecomplex *work, integer *lwork, 
	 integer *info);

/* Subroutine */ int zunm2l_(char *side, char *trans, integer *m, integer *n, 
	integer *k, doublecomplex *a, integer *lda, doublecomplex *tau, 
	doublecomplex *c__, integer *ldc, doublecomplex *work, integer *info);

/* Subroutine */ int zunm2r_(char *side, char *trans, integer *m, integer *n, 
	integer *k, doublecomplex *a, integer *lda, doublecomplex *tau, 
	doublecomplex *c__, integer *ldc, doublecomplex *work, integer *info);

/* Subroutine */ int zunmbr_(char *vect, char *side, char *trans, integer *m, 
	integer *n, integer *k, doublecomplex *a, integer *lda, doublecomplex 
	*tau, doublecomplex *c__, integer *ldc, doublecomplex *work, integer *
	lwork, integer *info);

/* Subroutine */ int zunmhr_(char *side, char *trans, integer *m, integer *n, 
	integer *ilo, integer *ihi, doublecomplex *a, integer *lda, 
	doublecomplex *tau, doublecomplex *c__, integer *ldc, doublecomplex *
	work, integer *lwork, integer *info);

/* Subroutine */ int zunml2_(char *side, char *trans, integer *m, integer *n, 
	integer *k, doublecomplex *a, integer *lda, doublecomplex *tau, 
	doublecomplex *c__, integer *ldc, doublecomplex *work, integer *info);

/* Subroutine */ int zunmlq_(char *side, char *trans, integer *m, integer *n, 
	integer *k, doublecomplex *a, integer *lda, doublecomplex *tau, 
	doublecomplex *c__, integer *ldc, doublecomplex *work, integer *lwork, 
	 integer *info);

/* Subroutine */ int zunmql_(char *side, char *trans, integer *m, integer *n, 
	integer *k, doublecomplex *a, integer *lda, doublecomplex *tau, 
	doublecomplex *c__, integer *ldc, doublecomplex *work, integer *lwork, 
	 integer *info);

/* Subroutine */ int zunmqr_(char *side, char *trans, integer *m, integer *n, 
	integer *k, doublecomplex *a, integer *lda, doublecomplex *tau, 
	doublecomplex *c__, integer *ldc, doublecomplex *work, integer *lwork, 
	 integer *info);

/* Subroutine */ int zunmr2_(char *side, char *trans, integer *m, integer *n, 
	integer *k, doublecomplex *a, integer *lda, doublecomplex *tau, 
	doublecomplex *c__, integer *ldc, doublecomplex *work, integer *info);

/* Subroutine */ int zunmr3_(char *side, char *trans, integer *m, integer *n, 
	integer *k, integer *l, doublecomplex *a, integer *lda, doublecomplex 
	*tau, doublecomplex *c__, integer *ldc, doublecomplex *work, integer *
	info);

/* Subroutine */ int zunmrq_(char *side, char *trans, integer *m, integer *n, 
	integer *k, doublecomplex *a, integer *lda, doublecomplex *tau, 
	doublecomplex *c__, integer *ldc, doublecomplex *work, integer *lwork, 
	 integer *info);

/* Subroutine */ int zunmrz_(char *side, char *trans, integer *m, integer *n, 
	integer *k, integer *l, doublecomplex *a, integer *lda, doublecomplex 
	*tau, doublecomplex *c__, integer *ldc, doublecomplex *work, integer *
	lwork, integer *info);

/* Subroutine */ int zunmtr_(char *side, char *uplo, char *trans, integer *m, 
	integer *n, doublecomplex *a, integer *lda, doublecomplex *tau, 
	doublecomplex *c__, integer *ldc, doublecomplex *work, integer *lwork, 
	 integer *info);

/* Subroutine */ int zupgtr_(char *uplo, integer *n, doublecomplex *ap, 
	doublecomplex *tau, doublecomplex *q, integer *ldq, doublecomplex *
	work, integer *info);

/* Subroutine */ int zupmtr_(char *side, char *uplo, char *trans, integer *m, 
	integer *n, doublecomplex *ap, doublecomplex *tau, doublecomplex *c__, 
	 integer *ldc, doublecomplex *work, integer *info);

#pragma endregion



#endif /* __CLAPACK_H */
