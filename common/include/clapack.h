#pragma once

#include "_f2c.h"
#include "cblas.h"


CLAPACK_API int fn_clapack(void);

//  variants  - 变体
/// Auxiliary - 辅助

//	ala		- ALLAUX  -- Auxiliary routines called from all precisions
///	aladz	- DZLAUX  -- Auxiliary routines called from all precisions but only from routines using extra precision.
//	alasc	- SCLAUX  -- Auxiliary routines called from both REALand COMPLEX

//	cla		- CLASRC  -- Single precision complex LAPACK routines
///	clax	- CXLASRC -- Single precision complex LAPACK routines using extra precision.

//	dla		- DLASRC  -- Double precision real LAPACK routines
///	dlax	- DXLASRC -- Double precision real LAPACK routines using extra precision.

//	sla		- SLASRC  -- Single precision real LAPACK routines
///	slax	- SXLASRC -- Single precision real LAPACK routines using extra precision.

//	zla		- ZLASRC  -- Double precision complex LAPACK routines
///	zlax	- ZXLASRC -- Double precision complex LAPACK routines using extra precision.




// A
#pragma region ALLAUX -- Auxiliary routines called from all precisions (ala) - ok

CLAPACK_API
VOID chla_transtype__(char* ret_val, ftnlen ret_val_len, integer* trans);

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


CLAPACK_API
int ilaver_(integer* vers_major__, integer* vers_minor__, integer* vers_patch__);

CLAPACK_API
integer iparmq_(integer* ispec, char* name__, char* opts, integer* n, integer
	* ilo, integer* ihi, integer* lwork);

CLAPACK_API
logical lsamen_(integer* n, char* ca, char* cb);

//}}@@@ finished.

#pragma endregion

#pragma region DZLAUX -- Auxiliary routines called from both DOUBLE PRECISION COMPLEX * 16 (aladz) - ok


CLAPACK_API
int dbdsdc_(char* uplo, char* compq, integer* n, doublereal*
	d__, doublereal* e, doublereal* u, integer* ldu, doublereal* vt,
	integer* ldvt, doublereal* q, integer* iq, doublereal* work, integer*
	iwork, integer* info);


CLAPACK_API
int dbdsqr_(char* uplo, integer* n, integer* ncvt, integer*
	nru, integer* ncc, doublereal* d__, doublereal* e, doublereal* vt,
	integer* ldvt, doublereal* u, integer* ldu, doublereal* c__, integer*
	ldc, doublereal* work, integer* info);


CLAPACK_API
int ddisna_(char* job, integer* m, integer* n, doublereal*
	d__, doublereal* sep, integer* info);

CLAPACK_API
logical disnan_(doublereal* din);


CLAPACK_API
int dlabad_(doublereal* small, doublereal* large);


CLAPACK_API
int dlacpy_(char* uplo, integer* m, integer* n, doublereal*
	a, integer* lda, doublereal* b, integer* ldb);


CLAPACK_API
int dladiv_(doublereal* a, doublereal* b, doublereal* c__,
	doublereal* d__, doublereal* p, doublereal* q);


CLAPACK_API
int dlae2_(doublereal* a, doublereal* b, doublereal* c__,
	doublereal* rt1, doublereal* rt2);


CLAPACK_API
int dlaebz_(integer* ijob, integer* nitmax, integer* n,
	integer* mmax, integer* minp, integer* nbmin, doublereal* abstol,
	doublereal* reltol, doublereal* pivmin, doublereal* d__, doublereal*
	e, doublereal* e2, integer* nval, doublereal* ab, doublereal* c__,
	integer* mout, integer* nab, doublereal* work, integer* iwork,
	integer* info);


CLAPACK_API
int dlaed0_(integer* icompq, integer* qsiz, integer* n,
	doublereal* d__, doublereal* e, doublereal* q, integer* ldq,
	doublereal* qstore, integer* ldqs, doublereal* work, integer* iwork,
	integer* info);


CLAPACK_API
int dlaed1_(integer* n, doublereal* d__, doublereal* q,
	integer* ldq, integer* indxq, doublereal* rho, integer* cutpnt,
	doublereal* work, integer* iwork, integer* info);


CLAPACK_API
int dlaed2_(integer* k, integer* n, integer* n1, doublereal*
	d__, doublereal* q, integer* ldq, integer* indxq, doublereal* rho,
	doublereal* z__, doublereal* dlamda, doublereal* w, doublereal* q2,
	integer* indx, integer* indxc, integer* indxp, integer* coltyp,
	integer* info);


CLAPACK_API
int dlaed3_(integer* k, integer* n, integer* n1, doublereal*
	d__, doublereal* q, integer* ldq, doublereal* rho, doublereal* dlamda,
	doublereal* q2, integer* indx, integer* ctot, doublereal* w,
	doublereal* s, integer* info);


CLAPACK_API
int dlaed4_(integer* n, integer* i__, doublereal* d__,
	doublereal* z__, doublereal* delta, doublereal* rho, doublereal* dlam,
	integer* info);


CLAPACK_API
int dlaed5_(integer* i__, doublereal* d__, doublereal* z__,
	doublereal* delta, doublereal* rho, doublereal* dlam);


CLAPACK_API
int dlaed6_(integer* kniter, logical* orgati, doublereal*
	rho, doublereal* d__, doublereal* z__, doublereal* finit, doublereal*
	tau, integer* info);


CLAPACK_API
int dlaed7_(integer* icompq, integer* n, integer* qsiz,
	integer* tlvls, integer* curlvl, integer* curpbm, doublereal* d__,
	doublereal* q, integer* ldq, integer* indxq, doublereal* rho, integer
	* cutpnt, doublereal* qstore, integer* qptr, integer* prmptr, integer*
	perm, integer* givptr, integer* givcol, doublereal* givnum,
	doublereal* work, integer* iwork, integer* info);


CLAPACK_API
int dlaed8_(integer* icompq, integer* k, integer* n, integer
	* qsiz, doublereal* d__, doublereal* q, integer* ldq, integer* indxq,
	doublereal* rho, integer* cutpnt, doublereal* z__, doublereal* dlamda,
	doublereal* q2, integer* ldq2, doublereal* w, integer* perm, integer
	* givptr, integer* givcol, doublereal* givnum, integer* indxp, integer
	* indx, integer* info);


CLAPACK_API
int dlaed9_(integer* k, integer* kstart, integer* kstop,
	integer* n, doublereal* d__, doublereal* q, integer* ldq, doublereal*
	rho, doublereal* dlamda, doublereal* w, doublereal* s, integer* lds,
	integer* info);


CLAPACK_API
int dlaeda_(integer* n, integer* tlvls, integer* curlvl,
	integer* curpbm, integer* prmptr, integer* perm, integer* givptr,
	integer* givcol, doublereal* givnum, doublereal* q, integer* qptr,
	doublereal* z__, doublereal* ztemp, integer* info);


CLAPACK_API
int dlaev2_(doublereal* a, doublereal* b, doublereal* c__,
	doublereal* rt1, doublereal* rt2, doublereal* cs1, doublereal* sn1);


CLAPACK_API
int dlagtf_(integer* n, doublereal* a, doublereal* lambda,
	doublereal* b, doublereal* c__, doublereal* tol, doublereal* d__,
	integer* in, integer* info);


CLAPACK_API
int dlagts_(integer* job, integer* n, doublereal* a,
	doublereal* b, doublereal* c__, doublereal* d__, integer* in,
	doublereal* y, doublereal* tol, integer* info);

CLAPACK_API
logical dlaisnan_(doublereal* din1, doublereal* din2);

CLAPACK_API
doublereal dlamch_(char* cmach);

CLAPACK_API
int dlamc1_(integer* beta, integer* t, logical* rnd, logical* ieee1);

CLAPACK_API
int dlamc2_(integer* beta, integer* t, logical* rnd,
	doublereal* eps, integer* emin, doublereal* rmin, integer* emax,
	doublereal* rmax);

CLAPACK_API
doublereal dlamc3_(doublereal* a, doublereal* b);

CLAPACK_API
int dlamc4_(integer* emin, doublereal* start, integer* base);

CLAPACK_API
int dlamc5_(integer* beta, integer* p, integer* emin,
	logical* ieee, integer* emax, doublereal* rmax);

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

CLAPACK_API
int dlarnv_(integer* idist, integer* iseed, integer* n,
	doublereal* x);

CLAPACK_API
int dlarra_(integer* n, doublereal* d__, doublereal* e,
	doublereal* e2, doublereal* spltol, doublereal* tnrm, integer* nsplit,
	integer* isplit, integer* info);

CLAPACK_API
int dlarrb_(integer* n, doublereal* d__, doublereal* lld,
	integer* ifirst, integer* ilast, doublereal* rtol1, doublereal* rtol2,
	integer* offset, doublereal* w, doublereal* wgap, doublereal* werr,
	doublereal* work, integer* iwork, doublereal* pivmin, doublereal*
	spdiam, integer* twist, integer* info);

CLAPACK_API
int dlarrc_(char* jobt, integer* n, doublereal* vl,
	doublereal* vu, doublereal* d__, doublereal* e, doublereal* pivmin,
	integer* eigcnt, integer* lcnt, integer* rcnt, integer* info);

CLAPACK_API
int dlarrd_(char* range, char* order, integer* n, doublereal
	* vl, doublereal* vu, integer* il, integer* iu, doublereal* gers,
	doublereal* reltol, doublereal* d__, doublereal* e, doublereal* e2,
	doublereal* pivmin, integer* nsplit, integer* isplit, integer* m,
	doublereal* w, doublereal* werr, doublereal* wl, doublereal* wu,
	integer* iblock, integer* indexw, doublereal* work, integer* iwork,
	integer* info);

CLAPACK_API
int dlarre_(char* range, integer* n, doublereal* vl,
	doublereal* vu, integer* il, integer* iu, doublereal* d__, doublereal
	* e, doublereal* e2, doublereal* rtol1, doublereal* rtol2, doublereal*
	spltol, integer* nsplit, integer* isplit, integer* m, doublereal* w,
	doublereal* werr, doublereal* wgap, integer* iblock, integer* indexw,
	doublereal* gers, doublereal* pivmin, doublereal* work, integer*
	iwork, integer* info);

CLAPACK_API
int dlarrf_(integer* n, doublereal* d__, doublereal* l,
	doublereal* ld, integer* clstrt, integer* clend, doublereal* w,
	doublereal* wgap, doublereal* werr, doublereal* spdiam, doublereal*
	clgapl, doublereal* clgapr, doublereal* pivmin, doublereal* sigma,
	doublereal* dplus, doublereal* lplus, doublereal* work, integer* info);

CLAPACK_API
int dlarrj_(integer* n, doublereal* d__, doublereal* e2,
	integer* ifirst, integer* ilast, doublereal* rtol, integer* offset,
	doublereal* w, doublereal* werr, doublereal* work, integer* iwork,
	doublereal* pivmin, doublereal* spdiam, integer* info);

CLAPACK_API
int dlarrk_(integer* n, integer* iw, doublereal* gl,
	doublereal* gu, doublereal* d__, doublereal* e2, doublereal* pivmin,
	doublereal* reltol, doublereal* w, doublereal* werr, integer* info);


CLAPACK_API
int dlarrr_(integer* n, doublereal* d__, doublereal* e,
	integer* info);


CLAPACK_API
int dlartg_(doublereal* f, doublereal* g, doublereal* cs,
	doublereal* sn, doublereal* r__);


CLAPACK_API
int dlaruv_(integer* iseed, integer* n, doublereal* x);


CLAPACK_API
int dlas2_(doublereal* f, doublereal* g, doublereal* h__,
	doublereal* ssmin, doublereal* ssmax);


CLAPACK_API
int dlascl_(char* type__, integer* kl, integer* ku,
	doublereal* cfrom, doublereal* cto, integer* m, integer* n,
	doublereal* a, integer* lda, integer* info);


CLAPACK_API
int dlasd0_(integer* n, integer* sqre, doublereal* d__,
	doublereal* e, doublereal* u, integer* ldu, doublereal* vt, integer*
	ldvt, integer* smlsiz, integer* iwork, doublereal* work, integer*
	info);


CLAPACK_API
int dlasd1_(integer* nl, integer* nr, integer* sqre,
	doublereal* d__, doublereal* alpha, doublereal* beta, doublereal* u,
	integer* ldu, doublereal* vt, integer* ldvt, integer* idxq, integer*
	iwork, doublereal* work, integer* info);


CLAPACK_API
int dlasd2_(integer* nl, integer* nr, integer* sqre, integer
	* k, doublereal* d__, doublereal* z__, doublereal* alpha, doublereal*
	beta, doublereal* u, integer* ldu, doublereal* vt, integer* ldvt,
	doublereal* dsigma, doublereal* u2, integer* ldu2, doublereal* vt2,
	integer* ldvt2, integer* idxp, integer* idx, integer* idxc, integer*
	idxq, integer* coltyp, integer* info);


CLAPACK_API
int dlasd3_(integer* nl, integer* nr, integer* sqre, integer
	* k, doublereal* d__, doublereal* q, integer* ldq, doublereal* dsigma,
	doublereal* u, integer* ldu, doublereal* u2, integer* ldu2,
	doublereal* vt, integer* ldvt, doublereal* vt2, integer* ldvt2,
	integer* idxc, integer* ctot, doublereal* z__, integer* info);


CLAPACK_API
int dlasd4_(integer* n, integer* i__, doublereal* d__,
	doublereal* z__, doublereal* delta, doublereal* rho, doublereal*
	sigma, doublereal* work, integer* info);


CLAPACK_API
int dlasd5_(integer* i__, doublereal* d__, doublereal* z__,
	doublereal* delta, doublereal* rho, doublereal* dsigma, doublereal*
	work);


CLAPACK_API
int dlasd6_(integer* icompq, integer* nl, integer* nr,
	integer* sqre, doublereal* d__, doublereal* vf, doublereal* vl,
	doublereal* alpha, doublereal* beta, integer* idxq, integer* perm,
	integer* givptr, integer* givcol, integer* ldgcol, doublereal* givnum,
	integer* ldgnum, doublereal* poles, doublereal* difl, doublereal*
	difr, doublereal* z__, integer* k, doublereal* c__, doublereal* s,
	doublereal* work, integer* iwork, integer* info);


CLAPACK_API
int dlasd7_(integer* icompq, integer* nl, integer* nr,
	integer* sqre, integer* k, doublereal* d__, doublereal* z__,
	doublereal* zw, doublereal* vf, doublereal* vfw, doublereal* vl,
	doublereal* vlw, doublereal* alpha, doublereal* beta, doublereal*
	dsigma, integer* idx, integer* idxp, integer* idxq, integer* perm,
	integer* givptr, integer* givcol, integer* ldgcol, doublereal* givnum,
	integer* ldgnum, doublereal* c__, doublereal* s, integer* info);


CLAPACK_API
int dlasd8_(integer* icompq, integer* k, doublereal* d__,
	doublereal* z__, doublereal* vf, doublereal* vl, doublereal* difl,
	doublereal* difr, integer* lddifr, doublereal* dsigma, doublereal*
	work, integer* info);


CLAPACK_API
int dlasda_(integer* icompq, integer* smlsiz, integer* n,
	integer* sqre, doublereal* d__, doublereal* e, doublereal* u, integer
	* ldu, doublereal* vt, integer* k, doublereal* difl, doublereal* difr,
	doublereal* z__, doublereal* poles, integer* givptr, integer* givcol,
	integer* ldgcol, integer* perm, doublereal* givnum, doublereal* c__,
	doublereal* s, doublereal* work, integer* iwork, integer* info);


CLAPACK_API
int dlasdq_(char* uplo, integer* sqre, integer* n, integer*
	ncvt, integer* nru, integer* ncc, doublereal* d__, doublereal* e,
	doublereal* vt, integer* ldvt, doublereal* u, integer* ldu,
	doublereal* c__, integer* ldc, doublereal* work, integer* info);


CLAPACK_API
int dlasdt_(integer* n, integer* lvl, integer* nd, integer*
	inode, integer* ndiml, integer* ndimr, integer* msub);


CLAPACK_API
int dlaset_(char* uplo, integer* m, integer* n, doublereal*
	alpha, doublereal* beta, doublereal* a, integer* lda);


CLAPACK_API
int dlasq1_(integer* n, doublereal* d__, doublereal* e,
	doublereal* work, integer* info);


CLAPACK_API
int dlasq2_(integer* n, doublereal* z__, integer* info);


CLAPACK_API
int dlasq3_(integer* i0, integer* n0, doublereal* z__,
	integer* pp, doublereal* dmin__, doublereal* sigma, doublereal* desig,
	doublereal* qmax, integer* nfail, integer* iter, integer* ndiv,
	logical* ieee, integer* ttype, doublereal* dmin1, doublereal* dmin2,
	doublereal* dn, doublereal* dn1, doublereal* dn2, doublereal* g,
	doublereal* tau);


CLAPACK_API
int dlasq4_(integer* i0, integer* n0, doublereal* z__,
	integer* pp, integer* n0in, doublereal* dmin__, doublereal* dmin1,
	doublereal* dmin2, doublereal* dn, doublereal* dn1, doublereal* dn2,
	doublereal* tau, integer* ttype, doublereal* g);


CLAPACK_API
int dlasq5_(integer* i0, integer* n0, doublereal* z__,
	integer* pp, doublereal* tau, doublereal* dmin__, doublereal* dmin1,
	doublereal* dmin2, doublereal* dn, doublereal* dnm1, doublereal* dnm2,
	logical* ieee);


CLAPACK_API
int dlasq6_(integer* i0, integer* n0, doublereal* z__,
	integer* pp, doublereal* dmin__, doublereal* dmin1, doublereal* dmin2,
	doublereal* dn, doublereal* dnm1, doublereal* dnm2);


CLAPACK_API
int dlasr_(char* side, char* pivot, char* direct, integer* m,
	integer* n, doublereal* c__, doublereal* s, doublereal* a, integer*
	lda);


CLAPACK_API
int dlasrt_(char* id, integer* n, doublereal* d__, integer*
	info);


CLAPACK_API
int dlassq_(integer* n, doublereal* x, integer* incx,
	doublereal* scale, doublereal* sumsq);


CLAPACK_API
int dlasv2_(doublereal* f, doublereal* g, doublereal* h__,
	doublereal* ssmin, doublereal* ssmax, doublereal* snr, doublereal*
	csr, doublereal* snl, doublereal* csl);

CLAPACK_API
integer dmaxloc_(doublereal* a, integer* dimm);


CLAPACK_API
int dpttrf_(integer* n, doublereal* d__, doublereal* e,
	integer* info);


CLAPACK_API
int dstebz_(char* range, char* order, integer* n, doublereal
	* vl, doublereal* vu, integer* il, integer* iu, doublereal* abstol,
	doublereal* d__, doublereal* e, integer* m, integer* nsplit,
	doublereal* w, integer* iblock, integer* isplit, doublereal* work,
	integer* iwork, integer* info);


CLAPACK_API
int dstedc_(char* compz, integer* n, doublereal* d__,
	doublereal* e, doublereal* z__, integer* ldz, doublereal* work,
	integer* lwork, integer* iwork, integer* liwork, integer* info);


CLAPACK_API
int dsteqr_(char* compz, integer* n, doublereal* d__,
	doublereal* e, doublereal* z__, integer* ldz, doublereal* work,
	integer* info);


CLAPACK_API
int dsterf_(integer* n, doublereal* d__, doublereal* e,
	integer* info);


CLAPACK_API
doublereal dsecnd_();

//}}@@@ finished

#pragma endregion

#pragma region SCLAUX -- Auxiliary routines called from both REALand COMPLEX (alasc) - ok


CLAPACK_API
int sbdsdc_(char* uplo, char* compq, integer* n, real* d__,
	real* e, real* u, integer* ldu, real* vt, integer* ldvt, real* q,
	integer* iq, real* work, integer* iwork, integer* info);

CLAPACK_API
int sbdsqr_(char* uplo, integer* n, integer* ncvt, integer*
	nru, integer* ncc, real* d__, real* e, real* vt, integer* ldvt, real*
	u, integer* ldu, real* c__, integer* ldc, real* work, integer* info);

CLAPACK_API
int sdisna_(char* job, integer* m, integer* n, real* d__,
	real* sep, integer* info);

CLAPACK_API
doublereal second_();

CLAPACK_API
logical sisnan_(real* sin__);

CLAPACK_API
int slabad_(real* small, real* large);

CLAPACK_API
int slacpy_(char* uplo, integer* m, integer* n, real* a,
	integer* lda, real* b, integer* ldb);

CLAPACK_API
int sladiv_(real* a, real* b, real* c__, real* d__, real* p,
	real* q);

CLAPACK_API
int slae2_(real* a, real* b, real* c__, real* rt1, real* rt2);

CLAPACK_API
int slaebz_(integer* ijob, integer* nitmax, integer* n,
	integer* mmax, integer* minp, integer* nbmin, real* abstol, real*
	reltol, real* pivmin, real* d__, real* e, real* e2, integer* nval,
	real* ab, real* c__, integer* mout, integer* nab, real* work, integer
	* iwork, integer* info);

CLAPACK_API
int slaed0_(integer* icompq, integer* qsiz, integer* n, real
	* d__, real* e, real* q, integer* ldq, real* qstore, integer* ldqs,
	real* work, integer* iwork, integer* info);

CLAPACK_API
int slaed1_(integer* n, real* d__, real* q, integer* ldq,
	integer* indxq, real* rho, integer* cutpnt, real* work, integer*
	iwork, integer* info);

CLAPACK_API
int slaed2_(integer* k, integer* n, integer* n1, real* d__,
	real* q, integer* ldq, integer* indxq, real* rho, real* z__, real*
	dlamda, real* w, real* q2, integer* indx, integer* indxc, integer*
	indxp, integer* coltyp, integer* info);

CLAPACK_API
int slaed3_(integer* k, integer* n, integer* n1, real* d__,
	real* q, integer* ldq, real* rho, real* dlamda, real* q2, integer*
	indx, integer* ctot, real* w, real* s, integer* info);

CLAPACK_API
int slaed4_(integer* n, integer* i__, real* d__, real* z__,
	real* delta, real* rho, real* dlam, integer* info);

CLAPACK_API
int slaed5_(integer* i__, real* d__, real* z__, real* delta,
	real* rho, real* dlam);


CLAPACK_API
int slaed6_(integer* kniter, logical* orgati, real* rho,
	real* d__, real* z__, real* finit, real* tau, integer* info);

CLAPACK_API
int slaed7_(integer* icompq, integer* n, integer* qsiz,
	integer* tlvls, integer* curlvl, integer* curpbm, real* d__, real* q,
	integer* ldq, integer* indxq, real* rho, integer* cutpnt, real*
	qstore, integer* qptr, integer* prmptr, integer* perm, integer*
	givptr, integer* givcol, real* givnum, real* work, integer* iwork,
	integer* info);

CLAPACK_API
int slaed8_(integer* icompq, integer* k, integer* n, integer
	* qsiz, real* d__, real* q, integer* ldq, integer* indxq, real* rho,
	integer* cutpnt, real* z__, real* dlamda, real* q2, integer* ldq2,
	real* w, integer* perm, integer* givptr, integer* givcol, real*
	givnum, integer* indxp, integer* indx, integer* info);

CLAPACK_API
int slaed9_(integer* k, integer* kstart, integer* kstop,
	integer* n, real* d__, real* q, integer* ldq, real* rho, real* dlamda,
	real* w, real* s, integer* lds, integer* info);

CLAPACK_API
int slaeda_(integer* n, integer* tlvls, integer* curlvl,
	integer* curpbm, integer* prmptr, integer* perm, integer* givptr,
	integer* givcol, real* givnum, real* q, integer* qptr, real* z__,
	real* ztemp, integer* info);

CLAPACK_API
int slaev2_(real* a, real* b, real* c__, real* rt1, real*
	rt2, real* cs1, real* sn1);

CLAPACK_API
int slagtf_(integer* n, real* a, real* lambda, real* b, real
	* c__, real* tol, real* d__, integer* in, integer* info);

CLAPACK_API
int slagts_(integer* job, integer* n, real* a, real* b, real
	* c__, real* d__, integer* in, real* y, real* tol, integer* info);

CLAPACK_API
logical slaisnan_(real* sin1, real* sin2);

CLAPACK_API
doublereal slamch_(char* cmach);

CLAPACK_API
int slamc1_(integer* beta, integer* t, logical* rnd, logical* ieee1);

CLAPACK_API
int slamc2_(integer* beta, integer* t, logical* rnd, real* eps,
	integer* emin, real* rmin, integer* emax, real* rmax);

CLAPACK_API
doublereal slamc3_(real* a, real* b);

CLAPACK_API
int slamc4_(integer* emin, real* start, integer* base);

CLAPACK_API
int slamc5_(integer* beta, integer* p, integer* emin,
	logical* ieee, integer* emax, real* rmax);

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

CLAPACK_API
int slarnv_(integer* idist, integer* iseed, integer* n, real* x);

CLAPACK_API
int slarra_(integer* n, real* d__, real* e, real* e2, real*
	spltol, real* tnrm, integer* nsplit, integer* isplit, integer* info);

CLAPACK_API
int slarrb_(integer* n, real* d__, real* lld, integer*
	ifirst, integer* ilast, real* rtol1, real* rtol2, integer* offset,
	real* w, real* wgap, real* werr, real* work, integer* iwork, real*
	pivmin, real* spdiam, integer* twist, integer* info);

CLAPACK_API
int slarrc_(char* jobt, integer* n, real* vl, real* vu, real
	* d__, real* e, real* pivmin, integer* eigcnt, integer* lcnt, integer*
	rcnt, integer* info);

CLAPACK_API
int slarrd_(char* range, char* order, integer* n, real* vl,
	real* vu, integer* il, integer* iu, real* gers, real* reltol, real*
	d__, real* e, real* e2, real* pivmin, integer* nsplit, integer*
	isplit, integer* m, real* w, real* werr, real* wl, real* wu, integer*
	iblock, integer* indexw, real* work, integer* iwork, integer* info);

CLAPACK_API
int slarre_(char* range, integer* n, real* vl, real* vu,
	integer* il, integer* iu, real* d__, real* e, real* e2, real* rtol1,
	real* rtol2, real* spltol, integer* nsplit, integer* isplit, integer*
	m, real* w, real* werr, real* wgap, integer* iblock, integer* indexw,
	real* gers, real* pivmin, real* work, integer* iwork, integer* info);

CLAPACK_API
int slarrf_(integer* n, real* d__, real* l, real* ld,
	integer* clstrt, integer* clend, real* w, real* wgap, real* werr,
	real* spdiam, real* clgapl, real* clgapr, real* pivmin, real* sigma,
	real* dplus, real* lplus, real* work, integer* info);

CLAPACK_API
int slarrj_(integer* n, real* d__, real* e2, integer* ifirst,
	integer* ilast, real* rtol, integer* offset, real* w, real* werr,
	real* work, integer* iwork, real* pivmin, real* spdiam, integer* info);

CLAPACK_API
int slarrk_(integer* n, integer* iw, real* gl, real* gu,
	real* d__, real* e2, real* pivmin, real* reltol, real* w, real* werr,
	integer* info);

CLAPACK_API
int slarrr_(integer* n, real* d__, real* e, integer* info);

CLAPACK_API
int slartg_(real* f, real* g, real* cs, real* sn, real* r__);

CLAPACK_API
int slaruv_(integer* iseed, integer* n, real* x);

CLAPACK_API
int slas2_(real* f, real* g, real* h__, real* ssmin, real* ssmax);

CLAPACK_API
int slascl_(char* type__, integer* kl, integer* ku, real*
	cfrom, real* cto, integer* m, integer* n, real* a, integer* lda,
	integer* info);

CLAPACK_API
int slasd0_(integer* n, integer* sqre, real* d__, real* e,
	real* u, integer* ldu, real* vt, integer* ldvt, integer* smlsiz,
	integer* iwork, real* work, integer* info);

CLAPACK_API
int slasd1_(integer* nl, integer* nr, integer* sqre, real*
	d__, real* alpha, real* beta, real* u, integer* ldu, real* vt,
	integer* ldvt, integer* idxq, integer* iwork, real* work, integer*
	info);

CLAPACK_API
int slasd2_(integer* nl, integer* nr, integer* sqre, integer
	* k, real* d__, real* z__, real* alpha, real* beta, real* u, integer*
	ldu, real* vt, integer* ldvt, real* dsigma, real* u2, integer* ldu2,
	real* vt2, integer* ldvt2, integer* idxp, integer* idx, integer* idxc,
	integer* idxq, integer* coltyp, integer* info);

CLAPACK_API
int slasd3_(integer* nl, integer* nr, integer* sqre, integer
	* k, real* d__, real* q, integer* ldq, real* dsigma, real* u, integer*
	ldu, real* u2, integer* ldu2, real* vt, integer* ldvt, real* vt2,
	integer* ldvt2, integer* idxc, integer* ctot, real* z__, integer*
	info);

CLAPACK_API
int slasd4_(integer* n, integer* i__, real* d__, real* z__,
	real* delta, real* rho, real* sigma, real* work, integer* info);

CLAPACK_API
int slasd5_(integer* i__, real* d__, real* z__, real* delta,
	real* rho, real* dsigma, real* work);

CLAPACK_API
int slasd6_(integer* icompq, integer* nl, integer* nr,
	integer* sqre, real* d__, real* vf, real* vl, real* alpha, real* beta,
	integer* idxq, integer* perm, integer* givptr, integer* givcol,
	integer* ldgcol, real* givnum, integer* ldgnum, real* poles, real*
	difl, real* difr, real* z__, integer* k, real* c__, real* s, real*
	work, integer* iwork, integer* info);

CLAPACK_API
int slasd7_(integer* icompq, integer* nl, integer* nr,
	integer* sqre, integer* k, real* d__, real* z__, real* zw, real* vf,
	real* vfw, real* vl, real* vlw, real* alpha, real* beta, real* dsigma,
	integer* idx, integer* idxp, integer* idxq, integer* perm, integer*
	givptr, integer* givcol, integer* ldgcol, real* givnum, integer*
	ldgnum, real* c__, real* s, integer* info);

CLAPACK_API
int slasd8_(integer* icompq, integer* k, real* d__, real*
	z__, real* vf, real* vl, real* difl, real* difr, integer* lddifr,
	real* dsigma, real* work, integer* info);

CLAPACK_API
int slasda_(integer* icompq, integer* smlsiz, integer* n,
	integer* sqre, real* d__, real* e, real* u, integer* ldu, real* vt,
	integer* k, real* difl, real* difr, real* z__, real* poles, integer*
	givptr, integer* givcol, integer* ldgcol, integer* perm, real* givnum,
	real* c__, real* s, real* work, integer* iwork, integer* info);

CLAPACK_API
int slasdq_(char* uplo, integer* sqre, integer* n, integer*
	ncvt, integer* nru, integer* ncc, real* d__, real* e, real* vt,
	integer* ldvt, real* u, integer* ldu, real* c__, integer* ldc, real*
	work, integer* info);

CLAPACK_API
int slasdt_(integer* n, integer* lvl, integer* nd, integer*
	inode, integer* ndiml, integer* ndimr, integer* msub);

CLAPACK_API
int slaset_(char* uplo, integer* m, integer* n, real* alpha,
	real* beta, real* a, integer* lda);

CLAPACK_API
int slasq1_(integer* n, real* d__, real* e, real* work,
	integer* info);

CLAPACK_API
int slasq2_(integer* n, real* z__, integer* info);

CLAPACK_API
int slasq3_(integer* i0, integer* n0, real* z__, integer* pp,
	real* dmin__, real* sigma, real* desig, real* qmax, integer* nfail,
	integer* iter, integer* ndiv, logical* ieee, integer* ttype, real*
	dmin1, real* dmin2, real* dn, real* dn1, real* dn2, real* g, real*
	tau);

CLAPACK_API
int slasq4_(integer* i0, integer* n0, real* z__, integer* pp,
	integer* n0in, real* dmin__, real* dmin1, real* dmin2, real* dn,
	real* dn1, real* dn2, real* tau, integer* ttype, real* g);

CLAPACK_API
int slasq5_(integer* i0, integer* n0, real* z__, integer* pp,
	real* tau, real* dmin__, real* dmin1, real* dmin2, real* dn, real*
	dnm1, real* dnm2, logical* ieee);

CLAPACK_API
int slasq6_(integer* i0, integer* n0, real* z__, integer* pp,
	real* dmin__, real* dmin1, real* dmin2, real* dn, real* dnm1, real*
	dnm2);

CLAPACK_API
int slasr_(char* side, char* pivot, char* direct, integer* m,
	integer* n, real* c__, real* s, real* a, integer* lda);

CLAPACK_API
int slasrt_(char* id, integer* n, real* d__, integer* info);

CLAPACK_API
int slassq_(integer* n, real* x, integer* incx, real* scale,
	real* sumsq);

CLAPACK_API
int slasv2_(real* f, real* g, real* h__, real* ssmin, real*
	ssmax, real* snr, real* csr, real* snl, real* csl);

CLAPACK_API
integer smaxloc_(real* a, integer* dimm);

CLAPACK_API
int spttrf_(integer* n, real* d__, real* e, integer* info);

CLAPACK_API
int sstebz_(char* range, char* order, integer* n, real* vl,
	real* vu, integer* il, integer* iu, real* abstol, real* d__, real* e,
	integer* m, integer* nsplit, real* w, integer* iblock, integer*
	isplit, real* work, integer* iwork, integer* info);

CLAPACK_API
int sstedc_(char* compz, integer* n, real* d__, real* e,
	real* z__, integer* ldz, real* work, integer* lwork, integer* iwork,
	integer* liwork, integer* info);

CLAPACK_API
int ssteqr_(char* compz, integer* n, real* d__, real* e,
	real* z__, integer* ldz, real* work, integer* info);

CLAPACK_API
int ssterf_(integer* n, real* d__, real* e, integer* info);

CLAPACK_API
doublereal second_();

//}}@@@ finished

#pragma endregion



// C
#pragma region CLASRC -- Single precision complex LAPACK routines (cla) - ok


CLAPACK_API
int cbdsqr_(char* uplo, integer* n, integer* ncvt, integer*
	nru, integer* ncc, real* d__, real* e, complex* vt, integer* ldvt,
	complex* u, integer* ldu, complex* c__, integer* ldc, real* rwork,
	integer* info);

CLAPACK_API
int cgbbrd_(char* vect, integer* m, integer* n, integer* ncc,
	integer* kl, integer* ku, complex* ab, integer* ldab, real* d__,
	real* e, complex* q, integer* ldq, complex* pt, integer* ldpt,
	complex* c__, integer* ldc, complex* work, real* rwork, integer* info);

CLAPACK_API
int cgbcon_(char* norm, integer* n, integer* kl, integer* ku,
	complex* ab, integer* ldab, integer* ipiv, real* anorm, real* rcond,
	complex* work, real* rwork, integer* info);

CLAPACK_API
int cgbequ_(integer* m, integer* n, integer* kl, integer* ku,
	complex* ab, integer* ldab, real* r__, real* c__, real* rowcnd, real
	* colcnd, real* amax, integer* info);

CLAPACK_API
int cgbequb_(integer* m, integer* n, integer* kl, integer*
	ku, complex* ab, integer* ldab, real* r__, real* c__, real* rowcnd,
	real* colcnd, real* amax, integer* info);

CLAPACK_API
int cgbrfs_(char* trans, integer* n, integer* kl, integer*
	ku, integer* nrhs, complex* ab, integer* ldab, complex* afb, integer*
	ldafb, integer* ipiv, complex* b, integer* ldb, complex* x, integer*
	ldx, real* ferr, real* berr, complex* work, real* rwork, integer*
	info);

CLAPACK_API
int cgbsv_(integer* n, integer* kl, integer* ku, integer*
	nrhs, complex* ab, integer* ldab, integer* ipiv, complex* b, integer*
	ldb, integer* info);

CLAPACK_API
int cgbsvx_(char* fact, char* trans, integer* n, integer* kl,
	integer* ku, integer* nrhs, complex* ab, integer* ldab, complex* afb,
	integer* ldafb, integer* ipiv, char* equed, real* r__, real* c__,
	complex* b, integer* ldb, complex* x, integer* ldx, real* rcond, real
	* ferr, real* berr, complex* work, real* rwork, integer* info);

CLAPACK_API
int cgbtf2_(integer* m, integer* n, integer* kl, integer* ku,
	complex* ab, integer* ldab, integer* ipiv, integer* info);

CLAPACK_API
int cgbtrf_(integer* m, integer* n, integer* kl, integer* ku,
	complex* ab, integer* ldab, integer* ipiv, integer* info);

CLAPACK_API
int cgbtrs_(char* trans, integer* n, integer* kl, integer*
	ku, integer* nrhs, complex* ab, integer* ldab, integer* ipiv, complex
	* b, integer* ldb, integer* info);

CLAPACK_API
int cgebak_(char* job, char* side, integer* n, integer* ilo,
	integer* ihi, real* scale, integer* m, complex* v, integer* ldv,
	integer* info);

CLAPACK_API
int cgebal_(char* job, integer* n, complex* a, integer* lda,
	integer* ilo, integer* ihi, real* scale, integer* info);

CLAPACK_API
int cgebd2_(integer* m, integer* n, complex* a, integer* lda,
	real* d__, real* e, complex* tauq, complex* taup, complex* work,
	integer* info);

CLAPACK_API
int cgebrd_(integer* m, integer* n, complex* a, integer* lda,
	real* d__, real* e, complex* tauq, complex* taup, complex* work,
	integer* lwork, integer* info);

CLAPACK_API
int cgecon_(char* norm, integer* n, complex* a, integer* lda,
	real* anorm, real* rcond, complex* work, real* rwork, integer* info);

CLAPACK_API
int cgeequ_(integer* m, integer* n, complex* a, integer* lda,
	real* r__, real* c__, real* rowcnd, real* colcnd, real* amax,
	integer* info);

CLAPACK_API
int cgeequb_(integer* m, integer* n, complex* a, integer*
	lda, real* r__, real* c__, real* rowcnd, real* colcnd, real* amax,
	integer* info);

CLAPACK_API
int cgees_(char* jobvs, char* sort, L_fp select, integer* n,
	complex* a, integer* lda, integer* sdim, complex* w, complex* vs,
	integer* ldvs, complex* work, integer* lwork, real* rwork, logical*
	bwork, integer* info);

CLAPACK_API
int cgeesx_(char* jobvs, char* sort, L_fp select, char*
	sense, integer* n, complex* a, integer* lda, integer* sdim, complex*
	w, complex* vs, integer* ldvs, real* rconde, real* rcondv, complex*
	work, integer* lwork, real* rwork, logical* bwork, integer* info);


CLAPACK_API
int cgeev_(char* jobvl, char* jobvr, integer* n, complex* a,
	integer* lda, complex* w, complex* vl, integer* ldvl, complex* vr,
	integer* ldvr, complex* work, integer* lwork, real* rwork, integer*
	info);


CLAPACK_API
int cgeevx_(char* balanc, char* jobvl, char* jobvr, char*
	sense, integer* n, complex* a, integer* lda, complex* w, complex* vl,
	integer* ldvl, complex* vr, integer* ldvr, integer* ilo, integer* ihi,
	real* scale, real* abnrm, real* rconde, real* rcondv, complex* work,
	integer* lwork, real* rwork, integer* info);


CLAPACK_API
int cgegs_(char* jobvsl, char* jobvsr, integer* n, complex*
	a, integer* lda, complex* b, integer* ldb, complex* alpha, complex*
	beta, complex* vsl, integer* ldvsl, complex* vsr, integer* ldvsr,
	complex* work, integer* lwork, real* rwork, integer* info);


CLAPACK_API
int cgegv_(char* jobvl, char* jobvr, integer* n, complex* a,
	integer* lda, complex* b, integer* ldb, complex* alpha, complex* beta,
	complex* vl, integer* ldvl, complex* vr, integer* ldvr, complex*
	work, integer* lwork, real* rwork, integer* info);


CLAPACK_API
int cgehd2_(integer* n, integer* ilo, integer* ihi, complex*
	a, integer* lda, complex* tau, complex* work, integer* info);


CLAPACK_API
int cgehrd_(integer* n, integer* ilo, integer* ihi, complex*
	a, integer* lda, complex* tau, complex* work, integer* lwork, integer
	* info);


CLAPACK_API
int cgelq2_(integer* m, integer* n, complex* a, integer* lda,
	complex* tau, complex* work, integer* info);


CLAPACK_API
int cgelqf_(integer* m, integer* n, complex* a, integer* lda,
	complex* tau, complex* work, integer* lwork, integer* info);


CLAPACK_API
int cgels_(char* trans, integer* m, integer* n, integer*
	nrhs, complex* a, integer* lda, complex* b, integer* ldb, complex*
	work, integer* lwork, integer* info);


CLAPACK_API
int cgelsd_(integer* m, integer* n, integer* nrhs, complex*
	a, integer* lda, complex* b, integer* ldb, real* s, real* rcond,
	integer* rank, complex* work, integer* lwork, real* rwork, integer*
	iwork, integer* info);


CLAPACK_API
int cgelss_(integer* m, integer* n, integer* nrhs, complex*
	a, integer* lda, complex* b, integer* ldb, real* s, real* rcond,
	integer* rank, complex* work, integer* lwork, real* rwork, integer*
	info);


CLAPACK_API
int cgelsx_(integer* m, integer* n, integer* nrhs, complex*
	a, integer* lda, complex* b, integer* ldb, integer* jpvt, real* rcond,
	integer* rank, complex* work, real* rwork, integer* info);


CLAPACK_API
int cgelsy_(integer* m, integer* n, integer* nrhs, complex*
	a, integer* lda, complex* b, integer* ldb, integer* jpvt, real* rcond,
	integer* rank, complex* work, integer* lwork, real* rwork, integer*
	info);


CLAPACK_API
int cgeql2_(integer* m, integer* n, complex* a, integer* lda,
	complex* tau, complex* work, integer* info);


CLAPACK_API
int cgeqlf_(integer* m, integer* n, complex* a, integer* lda,
	complex* tau, complex* work, integer* lwork, integer* info);


CLAPACK_API
int cgeqp3_(integer* m, integer* n, complex* a, integer* lda,
	integer* jpvt, complex* tau, complex* work, integer* lwork, real*
	rwork, integer* info);


CLAPACK_API
int cgeqpf_(integer* m, integer* n, complex* a, integer* lda,
	integer* jpvt, complex* tau, complex* work, real* rwork, integer*
	info);


CLAPACK_API
int cgeqr2_(integer* m, integer* n, complex* a, integer* lda,
	complex* tau, complex* work, integer* info);


CLAPACK_API
int cgeqrf_(integer* m, integer* n, complex* a, integer* lda,
	complex* tau, complex* work, integer* lwork, integer* info);


CLAPACK_API
int cgerfs_(char* trans, integer* n, integer* nrhs, complex*
	a, integer* lda, complex* af, integer* ldaf, integer* ipiv, complex*
	b, integer* ldb, complex* x, integer* ldx, real* ferr, real* berr,
	complex* work, real* rwork, integer* info);


CLAPACK_API
int cgerq2_(integer* m, integer* n, complex* a, integer* lda,
	complex* tau, complex* work, integer* info);


CLAPACK_API
int cgerqf_(integer* m, integer* n, complex* a, integer* lda,
	complex* tau, complex* work, integer* lwork, integer* info);


CLAPACK_API
int cgesc2_(integer* n, complex* a, integer* lda, complex*
	rhs, integer* ipiv, integer* jpiv, real* scale);


CLAPACK_API
int cgesdd_(char* jobz, integer* m, integer* n, complex* a,
	integer* lda, real* s, complex* u, integer* ldu, complex* vt, integer
	* ldvt, complex* work, integer* lwork, real* rwork, integer* iwork,
	integer* info);


CLAPACK_API
int cgesv_(integer* n, integer* nrhs, complex* a, integer*
	lda, integer* ipiv, complex* b, integer* ldb, integer* info);


CLAPACK_API
int cgesvd_(char* jobu, char* jobvt, integer* m, integer* n,
	complex* a, integer* lda, real* s, complex* u, integer* ldu, complex*
	vt, integer* ldvt, complex* work, integer* lwork, real* rwork,
	integer* info);


CLAPACK_API
int cgesvx_(char* fact, char* trans, integer* n, integer*
	nrhs, complex* a, integer* lda, complex* af, integer* ldaf, integer*
	ipiv, char* equed, real* r__, real* c__, complex* b, integer* ldb,
	complex* x, integer* ldx, real* rcond, real* ferr, real* berr,
	complex* work, real* rwork, integer* info);


CLAPACK_API
int cgetc2_(integer* n, complex* a, integer* lda, integer*
	ipiv, integer* jpiv, integer* info);


CLAPACK_API
int cgetf2_(integer* m, integer* n, complex* a, integer* lda,
	integer* ipiv, integer* info);


CLAPACK_API
int cgetrf_(integer* m, integer* n, complex* a, integer* lda,
	integer* ipiv, integer* info);


CLAPACK_API
int cgetri_(integer* n, complex* a, integer* lda, integer*
	ipiv, complex* work, integer* lwork, integer* info);


CLAPACK_API
int cgetrs_(char* trans, integer* n, integer* nrhs, complex*
	a, integer* lda, integer* ipiv, complex* b, integer* ldb, integer*
	info);


CLAPACK_API
int cggbak_(char* job, char* side, integer* n, integer* ilo,
	integer* ihi, real* lscale, real* rscale, integer* m, complex* v,
	integer* ldv, integer* info);


CLAPACK_API
int cggbal_(char* job, integer* n, complex* a, integer* lda,
	complex* b, integer* ldb, integer* ilo, integer* ihi, real* lscale,
	real* rscale, real* work, integer* info);


CLAPACK_API
int cgges_(char* jobvsl, char* jobvsr, char* sort, L_fp
	selctg, integer* n, complex* a, integer* lda, complex* b, integer*
	ldb, integer* sdim, complex* alpha, complex* beta, complex* vsl,
	integer* ldvsl, complex* vsr, integer* ldvsr, complex* work, integer*
	lwork, real* rwork, logical* bwork, integer* info);


CLAPACK_API
int cggesx_(char* jobvsl, char* jobvsr, char* sort, L_fp
	selctg, char* sense, integer* n, complex* a, integer* lda, complex* b,
	integer* ldb, integer* sdim, complex* alpha, complex* beta, complex*
	vsl, integer* ldvsl, complex* vsr, integer* ldvsr, real* rconde, real
	* rcondv, complex* work, integer* lwork, real* rwork, integer* iwork,
	integer* liwork, logical* bwork, integer* info);


CLAPACK_API
int cggev_(char* jobvl, char* jobvr, integer* n, complex* a,
	integer* lda, complex* b, integer* ldb, complex* alpha, complex* beta,
	complex* vl, integer* ldvl, complex* vr, integer* ldvr, complex*
	work, integer* lwork, real* rwork, integer* info);


CLAPACK_API
int cggevx_(char* balanc, char* jobvl, char* jobvr, char*
	sense, integer* n, complex* a, integer* lda, complex* b, integer* ldb,
	complex* alpha, complex* beta, complex* vl, integer* ldvl, complex*
	vr, integer* ldvr, integer* ilo, integer* ihi, real* lscale, real*
	rscale, real* abnrm, real* bbnrm, real* rconde, real* rcondv, complex
	* work, integer* lwork, real* rwork, integer* iwork, logical* bwork,
	integer* info);


CLAPACK_API
int cggglm_(integer* n, integer* m, integer* p, complex* a,
	integer* lda, complex* b, integer* ldb, complex* d__, complex* x,
	complex* y, complex* work, integer* lwork, integer* info);


CLAPACK_API
int cgghrd_(char* compq, char* compz, integer* n, integer*
	ilo, integer* ihi, complex* a, integer* lda, complex* b, integer* ldb,
	complex* q, integer* ldq, complex* z__, integer* ldz, integer* info);


CLAPACK_API
int cgglse_(integer* m, integer* n, integer* p, complex* a,
	integer* lda, complex* b, integer* ldb, complex* c__, complex* d__,
	complex* x, complex* work, integer* lwork, integer* info);


CLAPACK_API
int cggqrf_(integer* n, integer* m, integer* p, complex* a,
	integer* lda, complex* taua, complex* b, integer* ldb, complex* taub,
	complex* work, integer* lwork, integer* info);


CLAPACK_API
int cggrqf_(integer* m, integer* p, integer* n, complex* a,
	integer* lda, complex* taua, complex* b, integer* ldb, complex* taub,
	complex* work, integer* lwork, integer* info);


CLAPACK_API
int cggsvd_(char* jobu, char* jobv, char* jobq, integer* m,
	integer* n, integer* p, integer* k, integer* l, complex* a, integer*
	lda, complex* b, integer* ldb, real* alpha, real* beta, complex* u,
	integer* ldu, complex* v, integer* ldv, complex* q, integer* ldq,
	complex* work, real* rwork, integer* iwork, integer* info);


CLAPACK_API
int cggsvp_(char* jobu, char* jobv, char* jobq, integer* m,
	integer* p, integer* n, complex* a, integer* lda, complex* b, integer
	* ldb, real* tola, real* tolb, integer* k, integer* l, complex* u,
	integer* ldu, complex* v, integer* ldv, complex* q, integer* ldq,
	integer* iwork, real* rwork, complex* tau, complex* work, integer* info);


CLAPACK_API
int cgtcon_(char* norm, integer* n, complex* dl, complex*
	d__, complex* du, complex* du2, integer* ipiv, real* anorm, real*
	rcond, complex* work, integer* info);


CLAPACK_API
int cgtrfs_(char* trans, integer* n, integer* nrhs, complex*
	dl, complex* d__, complex* du, complex* dlf, complex* df, complex*
	duf, complex* du2, integer* ipiv, complex* b, integer* ldb, complex*
	x, integer* ldx, real* ferr, real* berr, complex* work, real* rwork,
	integer* info);


CLAPACK_API
int cgtsv_(integer* n, integer* nrhs, complex* dl, complex*
	d__, complex* du, complex* b, integer* ldb, integer* info);


CLAPACK_API
int cgtsvx_(char* fact, char* trans, integer* n, integer*
	nrhs, complex* dl, complex* d__, complex* du, complex* dlf, complex*
	df, complex* duf, complex* du2, integer* ipiv, complex* b, integer*
	ldb, complex* x, integer* ldx, real* rcond, real* ferr, real* berr,
	complex* work, real* rwork, integer* info);


CLAPACK_API
int cgttrf_(integer* n, complex* dl, complex* d__, complex*
	du, complex* du2, integer* ipiv, integer* info);


CLAPACK_API
int cgttrs_(char* trans, integer* n, integer* nrhs, complex*
	dl, complex* d__, complex* du, complex* du2, integer* ipiv, complex*
	b, integer* ldb, integer* info);


CLAPACK_API
int cgtts2_(integer* itrans, integer* n, integer* nrhs,
	complex* dl, complex* d__, complex* du, complex* du2, integer* ipiv,
	complex* b, integer* ldb);


CLAPACK_API
int chbev_(char* jobz, char* uplo, integer* n, integer* kd,
	complex* ab, integer* ldab, real* w, complex* z__, integer* ldz,
	complex* work, real* rwork, integer* info);


CLAPACK_API
int chbevd_(char* jobz, char* uplo, integer* n, integer* kd,
	complex* ab, integer* ldab, real* w, complex* z__, integer* ldz,
	complex* work, integer* lwork, real* rwork, integer* lrwork, integer*
	iwork, integer* liwork, integer* info);


CLAPACK_API
int chbevx_(char* jobz, char* range, char* uplo, integer* n,
	integer* kd, complex* ab, integer* ldab, complex* q, integer* ldq,
	real* vl, real* vu, integer* il, integer* iu, real* abstol, integer*
	m, real* w, complex* z__, integer* ldz, complex* work, real* rwork,
	integer* iwork, integer* ifail, integer* info);


CLAPACK_API
int chbgst_(char* vect, char* uplo, integer* n, integer* ka,
	integer* kb, complex* ab, integer* ldab, complex* bb, integer* ldbb,
	complex* x, integer* ldx, complex* work, real* rwork, integer* info);



CLAPACK_API
int chbgv_(char* jobz, char* uplo, integer* n, integer* ka,
	integer* kb, complex* ab, integer* ldab, complex* bb, integer* ldbb,
	real* w, complex* z__, integer* ldz, complex* work, real* rwork,
	integer* info);


CLAPACK_API
int chbgvd_(char* jobz, char* uplo, integer* n, integer* ka,
	integer* kb, complex* ab, integer* ldab, complex* bb, integer* ldbb,
	real* w, complex* z__, integer* ldz, complex* work, integer* lwork,
	real* rwork, integer* lrwork, integer* iwork, integer* liwork,
	integer* info);


CLAPACK_API
int chbgvx_(char* jobz, char* range, char* uplo, integer* n,
	integer* ka, integer* kb, complex* ab, integer* ldab, complex* bb,
	integer* ldbb, complex* q, integer* ldq, real* vl, real* vu, integer*
	il, integer* iu, real* abstol, integer* m, real* w, complex* z__,
	integer* ldz, complex* work, real* rwork, integer* iwork, integer*
	ifail, integer* info);



CLAPACK_API
int chbtrd_(char* vect, char* uplo, integer* n, integer* kd,
	complex* ab, integer* ldab, real* d__, real* e, complex* q, integer*
	ldq, complex* work, integer* info);


CLAPACK_API
int checon_(char* uplo, integer* n, complex* a, integer* lda,
	integer* ipiv, real* anorm, real* rcond, complex* work, integer*
	info);


CLAPACK_API
int cheequb_(char* uplo, integer* n, complex* a, integer*
	lda, real* s, real* scond, real* amax, complex* work, integer* info);


CLAPACK_API
int cheev_(char* jobz, char* uplo, integer* n, complex* a,
	integer* lda, real* w, complex* work, integer* lwork, real* rwork,
	integer* info);


CLAPACK_API
int cheevd_(char* jobz, char* uplo, integer* n, complex* a,
	integer* lda, real* w, complex* work, integer* lwork, real* rwork,
	integer* lrwork, integer* iwork, integer* liwork, integer* info);


CLAPACK_API
int cheevr_(char* jobz, char* range, char* uplo, integer* n,
	complex* a, integer* lda, real* vl, real* vu, integer* il, integer*
	iu, real* abstol, integer* m, real* w, complex* z__, integer* ldz,
	integer* isuppz, complex* work, integer* lwork, real* rwork, integer*
	lrwork, integer* iwork, integer* liwork, integer* info);

CLAPACK_API
int cheevx_(char* jobz, char* range, char* uplo, integer* n,
	complex* a, integer* lda, real* vl, real* vu, integer* il, integer*
	iu, real* abstol, integer* m, real* w, complex* z__, integer* ldz,
	complex* work, integer* lwork, real* rwork, integer* iwork, integer*
	ifail, integer* info);


CLAPACK_API
int chegs2_(integer* itype, char* uplo, integer* n, complex*
	a, integer* lda, complex* b, integer* ldb, integer* info);


CLAPACK_API
int chegst_(integer* itype, char* uplo, integer* n, complex*
	a, integer* lda, complex* b, integer* ldb, integer* info);


CLAPACK_API
int chegv_(integer* itype, char* jobz, char* uplo, integer*
	n, complex* a, integer* lda, complex* b, integer* ldb, real* w,
	complex* work, integer* lwork, real* rwork, integer* info);


CLAPACK_API
int chegvd_(integer* itype, char* jobz, char* uplo, integer*
	n, complex* a, integer* lda, complex* b, integer* ldb, real* w,
	complex* work, integer* lwork, real* rwork, integer* lrwork, integer*
	iwork, integer* liwork, integer* info);


CLAPACK_API
int chegvx_(integer* itype, char* jobz, char* range, char*
	uplo, integer* n, complex* a, integer* lda, complex* b, integer* ldb,
	real* vl, real* vu, integer* il, integer* iu, real* abstol, integer*
	m, real* w, complex* z__, integer* ldz, complex* work, integer* lwork,
	real* rwork, integer* iwork, integer* ifail, integer* info);


CLAPACK_API
int cherfs_(char* uplo, integer* n, integer* nrhs, complex*
	a, integer* lda, complex* af, integer* ldaf, integer* ipiv, complex*
	b, integer* ldb, complex* x, integer* ldx, real* ferr, real* berr,
	complex* work, real* rwork, integer* info);


CLAPACK_API
int chesv_(char* uplo, integer* n, integer* nrhs, complex* a,
	integer* lda, integer* ipiv, complex* b, integer* ldb, complex* work,
	integer* lwork, integer* info);


CLAPACK_API
int chesvx_(char* fact, char* uplo, integer* n, integer*
	nrhs, complex* a, integer* lda, complex* af, integer* ldaf, integer*
	ipiv, complex* b, integer* ldb, complex* x, integer* ldx, real* rcond,
	real* ferr, real* berr, complex* work, integer* lwork, real* rwork,
	integer* info);

CLAPACK_API
int chetd2_(char* uplo, integer* n, complex* a, integer* lda,
	real* d__, real* e, complex* tau, integer* info);

CLAPACK_API
int chetf2_(char* uplo, integer* n, complex* a, integer* lda,
	integer* ipiv, integer* info);

CLAPACK_API
int chetrd_(char* uplo, integer* n, complex* a, integer* lda,
	real* d__, real* e, complex* tau, complex* work, integer* lwork,
	integer* info);

CLAPACK_API
int chetrf_(char* uplo, integer* n, complex* a, integer* lda,
	integer* ipiv, complex* work, integer* lwork, integer* info);

CLAPACK_API
int chetri_(char* uplo, integer* n, complex* a, integer* lda,
	integer* ipiv, complex* work, integer* info);

CLAPACK_API
int chetrs_(char* uplo, integer* n, integer* nrhs, complex*
	a, integer* lda, integer* ipiv, complex* b, integer* ldb, integer*
	info);

CLAPACK_API
int chfrk_(char* transr, char* uplo, char* trans, integer* n,
	integer* k, real* alpha, complex* a, integer* lda, real* beta,
	complex* c__);

CLAPACK_API
int chgeqz_(char* job, char* compq, char* compz, integer* n,
	integer* ilo, integer* ihi, complex* h__, integer* ldh, complex* t,
	integer* ldt, complex* alpha, complex* beta, complex* q, integer* ldq,
	complex* z__, integer* ldz, complex* work, integer* lwork, real*
	rwork, integer* info);

CLAPACK_API int chpcon_(char* uplo, integer* n, complex* ap, integer*
	ipiv, real* anorm, real* rcond, complex* work, integer* info);

CLAPACK_API int chpev_(char* jobz, char* uplo, integer* n, complex* ap,
	real* w, complex* z__, integer* ldz, complex* work, real* rwork,
	integer* info);

CLAPACK_API int chpevd_(char* jobz, char* uplo, integer* n, complex* ap,
	real* w, complex* z__, integer* ldz, complex* work, integer* lwork,
	real* rwork, integer* lrwork, integer* iwork, integer* liwork,
	integer* info);

CLAPACK_API int chpevx_(char* jobz, char* range, char* uplo, integer* n,
	complex* ap, real* vl, real* vu, integer* il, integer* iu, real*
	abstol, integer* m, real* w, complex* z__, integer* ldz, complex*
	work, real* rwork, integer* iwork, integer* ifail, integer* info);

CLAPACK_API int chpgst_(integer* itype, char* uplo, integer* n, complex*
	ap, complex* bp, integer* info);

CLAPACK_API int chpgv_(integer* itype, char* jobz, char* uplo, integer*
	n, complex* ap, complex* bp, real* w, complex* z__, integer* ldz,
	complex* work, real* rwork, integer* info);

CLAPACK_API int chpgvd_(integer* itype, char* jobz, char* uplo, integer*
	n, complex* ap, complex* bp, real* w, complex* z__, integer* ldz,
	complex* work, integer* lwork, real* rwork, integer* lrwork, integer*
	iwork, integer* liwork, integer* info);

CLAPACK_API int chpgvx_(integer* itype, char* jobz, char* range, char*
	uplo, integer* n, complex* ap, complex* bp, real* vl, real* vu,
	integer* il, integer* iu, real* abstol, integer* m, real* w, complex*
	z__, integer* ldz, complex* work, real* rwork, integer* iwork,
	integer* ifail, integer* info);

CLAPACK_API int chprfs_(char* uplo, integer* n, integer* nrhs, complex*
	ap, complex* afp, integer* ipiv, complex* b, integer* ldb, complex* x,
	integer* ldx, real* ferr, real* berr, complex* work, real* rwork,
	integer* info);

CLAPACK_API int chpsv_(char* uplo, integer* n, integer* nrhs, complex*
	ap, integer* ipiv, complex* b, integer* ldb, integer* info);

CLAPACK_API int chpsvx_(char* fact, char* uplo, integer* n, integer*
	nrhs, complex* ap, complex* afp, integer* ipiv, complex* b, integer*
	ldb, complex* x, integer* ldx, real* rcond, real* ferr, real* berr,
	complex* work, real* rwork, integer* info);

CLAPACK_API int chptrd_(char* uplo, integer* n, complex* ap, real* d__,
	real* e, complex* tau, integer* info);

CLAPACK_API int chptrf_(char* uplo, integer* n, complex* ap, integer*
	ipiv, integer* info);

CLAPACK_API int chptri_(char* uplo, integer* n, complex* ap, integer*
	ipiv, complex* work, integer* info);

CLAPACK_API int chptrs_(char* uplo, integer* n, integer* nrhs, complex*
	ap, integer* ipiv, complex* b, integer* ldb, integer* info);

CLAPACK_API int chsein_(char* side, char* eigsrc, char* initv, logical*
	select, integer* n, complex* h__, integer* ldh, complex* w, complex*
	vl, integer* ldvl, complex* vr, integer* ldvr, integer* mm, integer*
	m, complex* work, real* rwork, integer* ifaill, integer* ifailr,
	integer* info);

CLAPACK_API int chseqr_(char* job, char* compz, integer* n, integer* ilo,
	integer* ihi, complex* h__, integer* ldh, complex* w, complex* z__,
	integer* ldz, complex* work, integer* lwork, integer* info);

CLAPACK_API int clabrd_(integer* m, integer* n, integer* nb, complex* a,
	integer* lda, real* d__, real* e, complex* tauq, complex* taup,
	complex* x, integer* ldx, complex* y, integer* ldy);

CLAPACK_API int clacgv_(integer* n, complex* x, integer* incx);

CLAPACK_API int clacn2_(integer* n, complex* v, complex* x, real* est,
	integer* kase, integer* isave);

CLAPACK_API int clacon_(integer* n, complex* v, complex* x, real* est,
	integer* kase);

CLAPACK_API int clacp2_(char* uplo, integer* m, integer* n, real* a,
	integer* lda, complex* b, integer* ldb);

CLAPACK_API int clacpy_(char* uplo, integer* m, integer* n, complex* a,
	integer* lda, complex* b, integer* ldb);

CLAPACK_API int clacrm_(integer* m, integer* n, complex* a, integer* lda,
	real* b, integer* ldb, complex* c__, integer* ldc, real* rwork);

CLAPACK_API int clacrt_(integer* n, complex* cx, integer* incx, complex*
	cy, integer* incy, complex* c__, complex* s);

/* Complex */ VOID cladiv_(complex* ret_val, complex* x, complex* y);

CLAPACK_API int claed0_(integer* qsiz, integer* n, real* d__, real* e,
	complex* q, integer* ldq, complex* qstore, integer* ldqs, real* rwork,
	integer* iwork, integer* info);

CLAPACK_API int claed7_(integer* n, integer* cutpnt, integer* qsiz,
	integer* tlvls, integer* curlvl, integer* curpbm, real* d__, complex*
	q, integer* ldq, real* rho, integer* indxq, real* qstore, integer*
	qptr, integer* prmptr, integer* perm, integer* givptr, integer*
	givcol, real* givnum, complex* work, real* rwork, integer* iwork,
	integer* info);

CLAPACK_API int claed8_(integer* k, integer* n, integer* qsiz, complex*
	q, integer* ldq, real* d__, real* rho, integer* cutpnt, real* z__,
	real* dlamda, complex* q2, integer* ldq2, real* w, integer* indxp,
	integer* indx, integer* indxq, integer* perm, integer* givptr,
	integer* givcol, real* givnum, integer* info);

CLAPACK_API int claein_(logical* rightv, logical* noinit, integer* n,
	complex* h__, integer* ldh, complex* w, complex* v, complex* b,
	integer* ldb, real* rwork, real* eps3, real* smlnum, integer* info);

CLAPACK_API int claesy_(complex* a, complex* b, complex* c__, complex*
	rt1, complex* rt2, complex* evscal, complex* cs1, complex* sn1);

CLAPACK_API int claev2_(complex* a, complex* b, complex* c__, real* rt1,
	real* rt2, real* cs1, complex* sn1);

CLAPACK_API int clag2z_(integer* m, integer* n, complex* sa, integer*
	ldsa, doublecomplex* a, integer* lda, integer* info);

CLAPACK_API int clags2_(logical* upper, real* a1, complex* a2, real* a3,
	real* b1, complex* b2, real* b3, real* csu, complex* snu, real* csv,
	complex* snv, real* csq, complex* snq);

CLAPACK_API int clagtm_(char* trans, integer* n, integer* nrhs, real*
	alpha, complex* dl, complex* d__, complex* du, complex* x, integer*
	ldx, real* beta, complex* b, integer* ldb);

CLAPACK_API int clahef_(char* uplo, integer* n, integer* nb, integer* kb,
	complex* a, integer* lda, integer* ipiv, complex* w, integer* ldw,
	integer* info);

CLAPACK_API int clahqr_(logical* wantt, logical* wantz, integer* n,
	integer* ilo, integer* ihi, complex* h__, integer* ldh, complex* w,
	integer* iloz, integer* ihiz, complex* z__, integer* ldz, integer*
	info);

CLAPACK_API int clahr2_(integer* n, integer* k, integer* nb, complex* a,
	integer* lda, complex* tau, complex* t, integer* ldt, complex* y,
	integer* ldy);

CLAPACK_API int clahrd_(integer* n, integer* k, integer* nb, complex* a,
	integer* lda, complex* tau, complex* t, integer* ldt, complex* y,
	integer* ldy);

CLAPACK_API int claic1_(integer* job, integer* j, complex* x, real* sest,
	complex* w, complex* gamma, real* sestpr, complex* s, complex* c__);

CLAPACK_API int clals0_(integer* icompq, integer* nl, integer* nr,
	integer* sqre, integer* nrhs, complex* b, integer* ldb, complex* bx,
	integer* ldbx, integer* perm, integer* givptr, integer* givcol,
	integer* ldgcol, real* givnum, integer* ldgnum, real* poles, real*
	difl, real* difr, real* z__, integer* k, real* c__, real* s, real*
	rwork, integer* info);

CLAPACK_API int clalsa_(integer* icompq, integer* smlsiz, integer* n,
	integer* nrhs, complex* b, integer* ldb, complex* bx, integer* ldbx,
	real* u, integer* ldu, real* vt, integer* k, real* difl, real* difr,
	real* z__, real* poles, integer* givptr, integer* givcol, integer*
	ldgcol, integer* perm, real* givnum, real* c__, real* s, real* rwork,
	integer* iwork, integer* info);

CLAPACK_API int clalsd_(char* uplo, integer* smlsiz, integer* n, integer
	* nrhs, real* d__, real* e, complex* b, integer* ldb, real* rcond,
	integer* rank, complex* work, real* rwork, integer* iwork, integer*
	info);

CLAPACK_API
doublereal clangb_(char* norm, integer* n, integer* kl, integer* ku, complex*
	ab, integer* ldab, real* work);

CLAPACK_API
doublereal clange_(char* norm, integer* m, integer* n, complex* a, integer*
	lda, real* work);

CLAPACK_API
doublereal clangt_(char* norm, integer* n, complex* dl, complex* d__, complex
	* du);

CLAPACK_API
doublereal clanhb_(char* norm, char* uplo, integer* n, integer* k, complex*
	ab, integer* ldab, real* work);

CLAPACK_API
doublereal clanhe_(char* norm, char* uplo, integer* n, complex* a, integer*
	lda, real* work);

CLAPACK_API
doublereal clanhf_(char* norm, char* transr, char* uplo, integer* n, complex*
	a, real* work);

CLAPACK_API
doublereal clanhp_(char* norm, char* uplo, integer* n, complex* ap, real*
	work);

CLAPACK_API
doublereal clanhs_(char* norm, integer* n, complex* a, integer* lda, real*
	work);

CLAPACK_API
doublereal clanht_(char* norm, integer* n, real* d__, complex* e);

CLAPACK_API
doublereal clansb_(char* norm, char* uplo, integer* n, integer* k, complex*
	ab, integer* ldab, real* work);

CLAPACK_API
doublereal clansp_(char* norm, char* uplo, integer* n, complex* ap, real*
	work);

CLAPACK_API
doublereal clansy_(char* norm, char* uplo, integer* n, complex* a, integer*
	lda, real* work);

CLAPACK_API
doublereal clantb_(char* norm, char* uplo, char* diag, integer* n, integer* k,
	complex* ab, integer* ldab, real* work);

CLAPACK_API
doublereal clantp_(char* norm, char* uplo, char* diag, integer* n, complex*
	ap, real* work);

CLAPACK_API
doublereal clantr_(char* norm, char* uplo, char* diag, integer* m, integer* n,
	complex* a, integer* lda, real* work);

CLAPACK_API int clapll_(integer* n, complex* x, integer* incx, complex*
	y, integer* incy, real* ssmin);

CLAPACK_API int clapmt_(logical* forwrd, integer* m, integer* n, complex
	* x, integer* ldx, integer* k);

CLAPACK_API int claqgb_(integer* m, integer* n, integer* kl, integer* ku,
	complex* ab, integer* ldab, real* r__, real* c__, real* rowcnd, real
	* colcnd, real* amax, char* equed);

CLAPACK_API int claqge_(integer* m, integer* n, complex* a, integer* lda,
	real* r__, real* c__, real* rowcnd, real* colcnd, real* amax, char*
	equed);

CLAPACK_API int claqhb_(char* uplo, integer* n, integer* kd, complex* ab,
	integer* ldab, real* s, real* scond, real* amax, char* equed);

CLAPACK_API int claqhe_(char* uplo, integer* n, complex* a, integer* lda,
	real* s, real* scond, real* amax, char* equed);

CLAPACK_API int claqhp_(char* uplo, integer* n, complex* ap, real* s,
	real* scond, real* amax, char* equed);

CLAPACK_API int claqp2_(integer* m, integer* n, integer* offset, complex
	* a, integer* lda, integer* jpvt, complex* tau, real* vn1, real* vn2,
	complex* work);

CLAPACK_API int claqps_(integer* m, integer* n, integer* offset, integer
	* nb, integer* kb, complex* a, integer* lda, integer* jpvt, complex*
	tau, real* vn1, real* vn2, complex* auxv, complex* f, integer* ldf);

CLAPACK_API int claqr0_(logical* wantt, logical* wantz, integer* n,
	integer* ilo, integer* ihi, complex* h__, integer* ldh, complex* w,
	integer* iloz, integer* ihiz, complex* z__, integer* ldz, complex*
	work, integer* lwork, integer* info);

CLAPACK_API int claqr1_(integer* n, complex* h__, integer* ldh, complex*
	s1, complex* s2, complex* v);

CLAPACK_API int claqr2_(logical* wantt, logical* wantz, integer* n,
	integer* ktop, integer* kbot, integer* nw, complex* h__, integer* ldh,
	integer* iloz, integer* ihiz, complex* z__, integer* ldz, integer*
	ns, integer* nd, complex* sh, complex* v, integer* ldv, integer* nh,
	complex* t, integer* ldt, integer* nv, complex* wv, integer* ldwv,
	complex* work, integer* lwork);

CLAPACK_API int claqr3_(logical* wantt, logical* wantz, integer* n,
	integer* ktop, integer* kbot, integer* nw, complex* h__, integer* ldh,
	integer* iloz, integer* ihiz, complex* z__, integer* ldz, integer*
	ns, integer* nd, complex* sh, complex* v, integer* ldv, integer* nh,
	complex* t, integer* ldt, integer* nv, complex* wv, integer* ldwv,
	complex* work, integer* lwork);

CLAPACK_API int claqr4_(logical* wantt, logical* wantz, integer* n,
	integer* ilo, integer* ihi, complex* h__, integer* ldh, complex* w,
	integer* iloz, integer* ihiz, complex* z__, integer* ldz, complex*
	work, integer* lwork, integer* info);

CLAPACK_API int claqr5_(logical* wantt, logical* wantz, integer* kacc22,
	integer* n, integer* ktop, integer* kbot, integer* nshfts, complex* s,
	complex* h__, integer* ldh, integer* iloz, integer* ihiz, complex*
	z__, integer* ldz, complex* v, integer* ldv, complex* u, integer* ldu,
	integer* nv, complex* wv, integer* ldwv, integer* nh, complex* wh,
	integer* ldwh);

CLAPACK_API int claqsb_(char* uplo, integer* n, integer* kd, complex* ab,
	integer* ldab, real* s, real* scond, real* amax, char* equed);

CLAPACK_API int claqsp_(char* uplo, integer* n, complex* ap, real* s,
	real* scond, real* amax, char* equed);

CLAPACK_API int claqsy_(char* uplo, integer* n, complex* a, integer* lda,
	real* s, real* scond, real* amax, char* equed);

CLAPACK_API int clar1v_(integer* n, integer* b1, integer* bn, real*
	lambda, real* d__, real* l, real* ld, real* lld, real* pivmin, real*
	gaptol, complex* z__, logical* wantnc, integer* negcnt, real* ztz,
	real* mingma, integer* r__, integer* isuppz, real* nrminv, real*
	resid, real* rqcorr, real* work);

CLAPACK_API int clar2v_(integer* n, complex* x, complex* y, complex* z__,
	integer* incx, real* c__, complex* s, integer* incc);

CLAPACK_API int clarcm_(integer* m, integer* n, real* a, integer* lda,
	complex* b, integer* ldb, complex* c__, integer* ldc, real* rwork);

CLAPACK_API int clarf_(char* side, integer* m, integer* n, complex* v,
	integer* incv, complex* tau, complex* c__, integer* ldc, complex*
	work);

CLAPACK_API int clarfb_(char* side, char* trans, char* direct, char*
	storev, integer* m, integer* n, integer* k, complex* v, integer* ldv,
	complex* t, integer* ldt, complex* c__, integer* ldc, complex* work,
	integer* ldwork);

CLAPACK_API int clarfg_(integer* n, complex* alpha, complex* x, integer*
	incx, complex* tau);

CLAPACK_API int clarfp_(integer* n, complex* alpha, complex* x, integer*
	incx, complex* tau);

CLAPACK_API int clarft_(char* direct, char* storev, integer* n, integer*
	k, complex* v, integer* ldv, complex* tau, complex* t, integer* ldt);

CLAPACK_API int clarfx_(char* side, integer* m, integer* n, complex* v,
	complex* tau, complex* c__, integer* ldc, complex* work);

CLAPACK_API int clargv_(integer* n, complex* x, integer* incx, complex*
	y, integer* incy, real* c__, integer* incc);

CLAPACK_API int clarnv_(integer* idist, integer* iseed, integer* n,
	complex* x);

CLAPACK_API int clarrv_(integer* n, real* vl, real* vu, real* d__, real*
	l, real* pivmin, integer* isplit, integer* m, integer* dol, integer*
	dou, real* minrgp, real* rtol1, real* rtol2, real* w, real* werr,
	real* wgap, integer* iblock, integer* indexw, real* gers, complex*
	z__, integer* ldz, integer* isuppz, real* work, integer* iwork,
	integer* info);

CLAPACK_API int clartg_(complex* f, complex* g, real* cs, complex* sn,
	complex* r__);

CLAPACK_API int clartv_(integer* n, complex* x, integer* incx, complex*
	y, integer* incy, real* c__, complex* s, integer* incc);

CLAPACK_API int clarz_(char* side, integer* m, integer* n, integer* l,
	complex* v, integer* incv, complex* tau, complex* c__, integer* ldc,
	complex* work);

CLAPACK_API int clarzb_(char* side, char* trans, char* direct, char*
	storev, integer* m, integer* n, integer* k, integer* l, complex* v,
	integer* ldv, complex* t, integer* ldt, complex* c__, integer* ldc,
	complex* work, integer* ldwork);

CLAPACK_API int clarzt_(char* direct, char* storev, integer* n, integer*
	k, complex* v, integer* ldv, complex* tau, complex* t, integer* ldt);

CLAPACK_API int clascl_(char* type__, integer* kl, integer* ku, real*
	cfrom, real* cto, integer* m, integer* n, complex* a, integer* lda,
	integer* info);

CLAPACK_API int claset_(char* uplo, integer* m, integer* n, complex*
	alpha, complex* beta, complex* a, integer* lda);

CLAPACK_API int clasr_(char* side, char* pivot, char* direct, integer* m,
	integer* n, real* c__, real* s, complex* a, integer* lda);

CLAPACK_API int classq_(integer* n, complex* x, integer* incx, real*
	scale, real* sumsq);

CLAPACK_API int claswp_(integer* n, complex* a, integer* lda, integer*
	k1, integer* k2, integer* ipiv, integer* incx);

CLAPACK_API int clasyf_(char* uplo, integer* n, integer* nb, integer* kb,
	complex* a, integer* lda, integer* ipiv, complex* w, integer* ldw,
	integer* info);

CLAPACK_API int clatbs_(char* uplo, char* trans, char* diag, char*
	normin, integer* n, integer* kd, complex* ab, integer* ldab, complex*
	x, real* scale, real* cnorm, integer* info);

CLAPACK_API int clatdf_(integer* ijob, integer* n, complex* z__, integer
	* ldz, complex* rhs, real* rdsum, real* rdscal, integer* ipiv, integer
	* jpiv);

CLAPACK_API int clatps_(char* uplo, char* trans, char* diag, char*
	normin, integer* n, complex* ap, complex* x, real* scale, real* cnorm,
	integer* info);

CLAPACK_API int clatrd_(char* uplo, integer* n, integer* nb, complex* a,
	integer* lda, real* e, complex* tau, complex* w, integer* ldw);

CLAPACK_API int clatrs_(char* uplo, char* trans, char* diag, char*
	normin, integer* n, complex* a, integer* lda, complex* x, real* scale,
	real* cnorm, integer* info);

CLAPACK_API int clatrz_(integer* m, integer* n, integer* l, complex* a,
	integer* lda, complex* tau, complex* work);

CLAPACK_API int clatzm_(char* side, integer* m, integer* n, complex* v,
	integer* incv, complex* tau, complex* c1, complex* c2, integer* ldc,
	complex* work);

CLAPACK_API int clauu2_(char* uplo, integer* n, complex* a, integer* lda,
	integer* info);

CLAPACK_API int clauum_(char* uplo, integer* n, complex* a, integer* lda,
	integer* info);

CLAPACK_API int cpbcon_(char* uplo, integer* n, integer* kd, complex* ab,
	integer* ldab, real* anorm, real* rcond, complex* work, real* rwork,
	integer* info);

CLAPACK_API int cpbequ_(char* uplo, integer* n, integer* kd, complex* ab,
	integer* ldab, real* s, real* scond, real* amax, integer* info);

CLAPACK_API int cpbrfs_(char* uplo, integer* n, integer* kd, integer*
	nrhs, complex* ab, integer* ldab, complex* afb, integer* ldafb,
	complex* b, integer* ldb, complex* x, integer* ldx, real* ferr, real*
	berr, complex* work, real* rwork, integer* info);

CLAPACK_API int cpbstf_(char* uplo, integer* n, integer* kd, complex* ab,
	integer* ldab, integer* info);

CLAPACK_API int cpbsv_(char* uplo, integer* n, integer* kd, integer*
	nrhs, complex* ab, integer* ldab, complex* b, integer* ldb, integer*
	info);

CLAPACK_API int cpbsvx_(char* fact, char* uplo, integer* n, integer* kd,
	integer* nrhs, complex* ab, integer* ldab, complex* afb, integer*
	ldafb, char* equed, real* s, complex* b, integer* ldb, complex* x,
	integer* ldx, real* rcond, real* ferr, real* berr, complex* work,
	real* rwork, integer* info);

CLAPACK_API int cpbtf2_(char* uplo, integer* n, integer* kd, complex* ab,
	integer* ldab, integer* info);

CLAPACK_API int cpbtrf_(char* uplo, integer* n, integer* kd, complex* ab,
	integer* ldab, integer* info);

CLAPACK_API int cpbtrs_(char* uplo, integer* n, integer* kd, integer*
	nrhs, complex* ab, integer* ldab, complex* b, integer* ldb, integer*
	info);

CLAPACK_API int cpftrf_(char* transr, char* uplo, integer* n, complex* a,
	integer* info);

CLAPACK_API int cpftri_(char* transr, char* uplo, integer* n, complex* a,
	integer* info);

CLAPACK_API int cpftrs_(char* transr, char* uplo, integer* n, integer*
	nrhs, complex* a, complex* b, integer* ldb, integer* info);

CLAPACK_API int cpocon_(char* uplo, integer* n, complex* a, integer* lda,
	real* anorm, real* rcond, complex* work, real* rwork, integer* info);

CLAPACK_API int cpoequ_(integer* n, complex* a, integer* lda, real* s,
	real* scond, real* amax, integer* info);

CLAPACK_API int cpoequb_(integer* n, complex* a, integer* lda, real* s,
	real* scond, real* amax, integer* info);

CLAPACK_API int cporfs_(char* uplo, integer* n, integer* nrhs, complex*
	a, integer* lda, complex* af, integer* ldaf, complex* b, integer* ldb,
	complex* x, integer* ldx, real* ferr, real* berr, complex* work,
	real* rwork, integer* info);

CLAPACK_API int cposv_(char* uplo, integer* n, integer* nrhs, complex* a,
	integer* lda, complex* b, integer* ldb, integer* info);

CLAPACK_API int cposvx_(char* fact, char* uplo, integer* n, integer*
	nrhs, complex* a, integer* lda, complex* af, integer* ldaf, char*
	equed, real* s, complex* b, integer* ldb, complex* x, integer* ldx,
	real* rcond, real* ferr, real* berr, complex* work, real* rwork,
	integer* info);

CLAPACK_API int cpotf2_(char* uplo, integer* n, complex* a, integer* lda,
	integer* info);

CLAPACK_API int cpotrf_(char* uplo, integer* n, complex* a, integer* lda,
	integer* info);

CLAPACK_API int cpotri_(char* uplo, integer* n, complex* a, integer* lda,
	integer* info);

CLAPACK_API int cpotrs_(char* uplo, integer* n, integer* nrhs, complex*
	a, integer* lda, complex* b, integer* ldb, integer* info);

CLAPACK_API int cppcon_(char* uplo, integer* n, complex* ap, real* anorm,
	real* rcond, complex* work, real* rwork, integer* info);

CLAPACK_API int cppequ_(char* uplo, integer* n, complex* ap, real* s,
	real* scond, real* amax, integer* info);

CLAPACK_API int cpprfs_(char* uplo, integer* n, integer* nrhs, complex*
	ap, complex* afp, complex* b, integer* ldb, complex* x, integer* ldx,
	real* ferr, real* berr, complex* work, real* rwork, integer* info);

CLAPACK_API int cppsv_(char* uplo, integer* n, integer* nrhs, complex*
	ap, complex* b, integer* ldb, integer* info);

CLAPACK_API int cppsvx_(char* fact, char* uplo, integer* n, integer*
	nrhs, complex* ap, complex* afp, char* equed, real* s, complex* b,
	integer* ldb, complex* x, integer* ldx, real* rcond, real* ferr, real
	* berr, complex* work, real* rwork, integer* info);

CLAPACK_API int cpptrf_(char* uplo, integer* n, complex* ap, integer*
	info);

CLAPACK_API int cpptri_(char* uplo, integer* n, complex* ap, integer*
	info);

CLAPACK_API int cpptrs_(char* uplo, integer* n, integer* nrhs, complex*
	ap, complex* b, integer* ldb, integer* info);

CLAPACK_API int cpstf2_(char* uplo, integer* n, complex* a, integer* lda,
	integer* piv, integer* rank, real* tol, real* work, integer* info);

CLAPACK_API int cpstrf_(char* uplo, integer* n, complex* a, integer* lda,
	integer* piv, integer* rank, real* tol, real* work, integer* info);

CLAPACK_API int cptcon_(integer* n, real* d__, complex* e, real* anorm,
	real* rcond, real* rwork, integer* info);

CLAPACK_API int cpteqr_(char* compz, integer* n, real* d__, real* e,
	complex* z__, integer* ldz, real* work, integer* info);

CLAPACK_API int cptrfs_(char* uplo, integer* n, integer* nrhs, real* d__,
	complex* e, real* df, complex* ef, complex* b, integer* ldb, complex
	* x, integer* ldx, real* ferr, real* berr, complex* work, real* rwork,
	integer* info);

CLAPACK_API int cptsv_(integer* n, integer* nrhs, real* d__, complex* e,
	complex* b, integer* ldb, integer* info);

CLAPACK_API int cptsvx_(char* fact, integer* n, integer* nrhs, real* d__,
	complex* e, real* df, complex* ef, complex* b, integer* ldb, complex
	* x, integer* ldx, real* rcond, real* ferr, real* berr, complex* work,
	real* rwork, integer* info);

CLAPACK_API int cpttrf_(integer* n, real* d__, complex* e, integer* info);

CLAPACK_API int cpttrs_(char* uplo, integer* n, integer* nrhs, real* d__,
	complex* e, complex* b, integer* ldb, integer* info);

CLAPACK_API int cptts2_(integer* iuplo, integer* n, integer* nrhs, real*
	d__, complex* e, complex* b, integer* ldb);

CLAPACK_API int crot_(integer* n, complex* cx, integer* incx, complex*
	cy, integer* incy, real* c__, complex* s);

CLAPACK_API int cspcon_(char* uplo, integer* n, complex* ap, integer*
	ipiv, real* anorm, real* rcond, complex* work, integer* info);

CLAPACK_API int cspmv_(char* uplo, integer* n, complex* alpha, complex*
	ap, complex* x, integer* incx, complex* beta, complex* y, integer*
	incy);

CLAPACK_API int cspr_(char* uplo, integer* n, complex* alpha, complex* x,
	integer* incx, complex* ap);

CLAPACK_API int csprfs_(char* uplo, integer* n, integer* nrhs, complex*
	ap, complex* afp, integer* ipiv, complex* b, integer* ldb, complex* x,
	integer* ldx, real* ferr, real* berr, complex* work, real* rwork,
	integer* info);

CLAPACK_API int cspsv_(char* uplo, integer* n, integer* nrhs, complex*
	ap, integer* ipiv, complex* b, integer* ldb, integer* info);

CLAPACK_API int cspsvx_(char* fact, char* uplo, integer* n, integer*
	nrhs, complex* ap, complex* afp, integer* ipiv, complex* b, integer*
	ldb, complex* x, integer* ldx, real* rcond, real* ferr, real* berr,
	complex* work, real* rwork, integer* info);

CLAPACK_API int csptrf_(char* uplo, integer* n, complex* ap, integer*
	ipiv, integer* info);

CLAPACK_API int csptri_(char* uplo, integer* n, complex* ap, integer*
	ipiv, complex* work, integer* info);

CLAPACK_API int csptrs_(char* uplo, integer* n, integer* nrhs, complex*
	ap, integer* ipiv, complex* b, integer* ldb, integer* info);

CLAPACK_API int csrscl_(integer* n, real* sa, complex* sx, integer* incx);

CLAPACK_API int cstedc_(char* compz, integer* n, real* d__, real* e,
	complex* z__, integer* ldz, complex* work, integer* lwork, real*
	rwork, integer* lrwork, integer* iwork, integer* liwork, integer*
	info);

CLAPACK_API int cstegr_(char* jobz, char* range, integer* n, real* d__,
	real* e, real* vl, real* vu, integer* il, integer* iu, real* abstol,
	integer* m, real* w, complex* z__, integer* ldz, integer* isuppz,
	real* work, integer* lwork, integer* iwork, integer* liwork, integer*
	info);

CLAPACK_API int cstein_(integer* n, real* d__, real* e, integer* m, real
	* w, integer* iblock, integer* isplit, complex* z__, integer* ldz,
	real* work, integer* iwork, integer* ifail, integer* info);

CLAPACK_API int cstemr_(char* jobz, char* range, integer* n, real* d__,
	real* e, real* vl, real* vu, integer* il, integer* iu, integer* m,
	real* w, complex* z__, integer* ldz, integer* nzc, integer* isuppz,
	logical* tryrac, real* work, integer* lwork, integer* iwork, integer*
	liwork, integer* info);

CLAPACK_API int csteqr_(char* compz, integer* n, real* d__, real* e,
	complex* z__, integer* ldz, real* work, integer* info);

CLAPACK_API int csycon_(char* uplo, integer* n, complex* a, integer* lda,
	integer* ipiv, real* anorm, real* rcond, complex* work, integer*
	info);

CLAPACK_API int csyequb_(char* uplo, integer* n, complex* a, integer*
	lda, real* s, real* scond, real* amax, complex* work, integer* info);

CLAPACK_API int csymv_(char* uplo, integer* n, complex* alpha, complex*
	a, integer* lda, complex* x, integer* incx, complex* beta, complex* y,
	integer* incy);

CLAPACK_API int csyr_(char* uplo, integer* n, complex* alpha, complex* x,
	integer* incx, complex* a, integer* lda);

CLAPACK_API int csyrfs_(char* uplo, integer* n, integer* nrhs, complex*
	a, integer* lda, complex* af, integer* ldaf, integer* ipiv, complex*
	b, integer* ldb, complex* x, integer* ldx, real* ferr, real* berr,
	complex* work, real* rwork, integer* info);

CLAPACK_API int csyrfsx_(char* uplo, char* equed, integer* n, integer*
	nrhs, complex* a, integer* lda, complex* af, integer* ldaf, integer*
	ipiv, real* s, complex* b, integer* ldb, complex* x, integer* ldx,
	real* rcond, real* berr, integer* n_err_bnds__, real* err_bnds_norm__,
	real* err_bnds_comp__, integer* nparams, real* params, complex* work,
	real* rwork, integer* info);

CLAPACK_API int csysv_(char* uplo, integer* n, integer* nrhs, complex* a,
	integer* lda, integer* ipiv, complex* b, integer* ldb, complex* work,
	integer* lwork, integer* info);

CLAPACK_API int csysvx_(char* fact, char* uplo, integer* n, integer*
	nrhs, complex* a, integer* lda, complex* af, integer* ldaf, integer*
	ipiv, complex* b, integer* ldb, complex* x, integer* ldx, real* rcond,
	real* ferr, real* berr, complex* work, integer* lwork, real* rwork,
	integer* info);

CLAPACK_API int csytf2_(char* uplo, integer* n, complex* a, integer* lda,
	integer* ipiv, integer* info);

CLAPACK_API int csytrf_(char* uplo, integer* n, complex* a, integer* lda,
	integer* ipiv, complex* work, integer* lwork, integer* info);

CLAPACK_API int csytri_(char* uplo, integer* n, complex* a, integer* lda,
	integer* ipiv, complex* work, integer* info);

CLAPACK_API int csytrs_(char* uplo, integer* n, integer* nrhs, complex*
	a, integer* lda, integer* ipiv, complex* b, integer* ldb, integer*
	info);

CLAPACK_API int ctbcon_(char* norm, char* uplo, char* diag, integer* n,
	integer* kd, complex* ab, integer* ldab, real* rcond, complex* work,
	real* rwork, integer* info);

CLAPACK_API int ctbrfs_(char* uplo, char* trans, char* diag, integer* n,
	integer* kd, integer* nrhs, complex* ab, integer* ldab, complex* b,
	integer* ldb, complex* x, integer* ldx, real* ferr, real* berr,
	complex* work, real* rwork, integer* info);

CLAPACK_API int ctbtrs_(char* uplo, char* trans, char* diag, integer* n,
	integer* kd, integer* nrhs, complex* ab, integer* ldab, complex* b,
	integer* ldb, integer* info);

CLAPACK_API int ctfsm_(char* transr, char* side, char* uplo, char* trans,
	char* diag, integer* m, integer* n, complex* alpha, complex* a,
	complex* b, integer* ldb);

CLAPACK_API int ctftri_(char* transr, char* uplo, char* diag, integer* n,
	complex* a, integer* info);

CLAPACK_API int ctfttp_(char* transr, char* uplo, integer* n, complex*
	arf, complex* ap, integer* info);

CLAPACK_API int ctfttr_(char* transr, char* uplo, integer* n, complex*
	arf, complex* a, integer* lda, integer* info);

CLAPACK_API int ctgevc_(char* side, char* howmny, logical* select,
	integer* n, complex* s, integer* lds, complex* p, integer* ldp,
	complex* vl, integer* ldvl, complex* vr, integer* ldvr, integer* mm,
	integer* m, complex* work, real* rwork, integer* info);

CLAPACK_API int ctgex2_(logical* wantq, logical* wantz, integer* n,
	complex* a, integer* lda, complex* b, integer* ldb, complex* q,
	integer* ldq, complex* z__, integer* ldz, integer* j1, integer* info);

CLAPACK_API int ctgexc_(logical* wantq, logical* wantz, integer* n,
	complex* a, integer* lda, complex* b, integer* ldb, complex* q,
	integer* ldq, complex* z__, integer* ldz, integer* ifst, integer*
	ilst, integer* info);

CLAPACK_API int ctgsen_(integer* ijob, logical* wantq, logical* wantz,
	logical* select, integer* n, complex* a, integer* lda, complex* b,
	integer* ldb, complex* alpha, complex* beta, complex* q, integer* ldq,
	complex* z__, integer* ldz, integer* m, real* pl, real* pr, real*
	dif, complex* work, integer* lwork, integer* iwork, integer* liwork,
	integer* info);

CLAPACK_API int ctgsja_(char* jobu, char* jobv, char* jobq, integer* m,
	integer* p, integer* n, integer* k, integer* l, complex* a, integer*
	lda, complex* b, integer* ldb, real* tola, real* tolb, real* alpha,
	real* beta, complex* u, integer* ldu, complex* v, integer* ldv,
	complex* q, integer* ldq, complex* work, integer* ncycle, integer*
	info);

CLAPACK_API int ctgsna_(char* job, char* howmny, logical* select,
	integer* n, complex* a, integer* lda, complex* b, integer* ldb,
	complex* vl, integer* ldvl, complex* vr, integer* ldvr, real* s, real
	* dif, integer* mm, integer* m, complex* work, integer* lwork, integer
	* iwork, integer* info);

CLAPACK_API int ctgsy2_(char* trans, integer* ijob, integer* m, integer*
	n, complex* a, integer* lda, complex* b, integer* ldb, complex* c__,
	integer* ldc, complex* d__, integer* ldd, complex* e, integer* lde,
	complex* f, integer* ldf, real* scale, real* rdsum, real* rdscal,
	integer* info);

CLAPACK_API int ctgsyl_(char* trans, integer* ijob, integer* m, integer*
	n, complex* a, integer* lda, complex* b, integer* ldb, complex* c__,
	integer* ldc, complex* d__, integer* ldd, complex* e, integer* lde,
	complex* f, integer* ldf, real* scale, real* dif, complex* work,
	integer* lwork, integer* iwork, integer* info);

CLAPACK_API int ctpcon_(char* norm, char* uplo, char* diag, integer* n,
	complex* ap, real* rcond, complex* work, real* rwork, integer* info);

CLAPACK_API int ctprfs_(char* uplo, char* trans, char* diag, integer* n,
	integer* nrhs, complex* ap, complex* b, integer* ldb, complex* x,
	integer* ldx, real* ferr, real* berr, complex* work, real* rwork,
	integer* info);

CLAPACK_API int ctptri_(char* uplo, char* diag, integer* n, complex* ap,
	integer* info);

CLAPACK_API int ctptrs_(char* uplo, char* trans, char* diag, integer* n,
	integer* nrhs, complex* ap, complex* b, integer* ldb, integer* info);

CLAPACK_API int ctpttf_(char* transr, char* uplo, integer* n, complex*
	ap, complex* arf, integer* info);

CLAPACK_API int ctpttr_(char* uplo, integer* n, complex* ap, complex* a,
	integer* lda, integer* info);

CLAPACK_API int ctrcon_(char* norm, char* uplo, char* diag, integer* n,
	complex* a, integer* lda, real* rcond, complex* work, real* rwork,
	integer* info);

CLAPACK_API int ctrevc_(char* side, char* howmny, logical* select,
	integer* n, complex* t, integer* ldt, complex* vl, integer* ldvl,
	complex* vr, integer* ldvr, integer* mm, integer* m, complex* work,
	real* rwork, integer* info);

CLAPACK_API int ctrexc_(char* compq, integer* n, complex* t, integer*
	ldt, complex* q, integer* ldq, integer* ifst, integer* ilst, integer*
	info);

CLAPACK_API int ctrrfs_(char* uplo, char* trans, char* diag, integer* n,
	integer* nrhs, complex* a, integer* lda, complex* b, integer* ldb,
	complex* x, integer* ldx, real* ferr, real* berr, complex* work, real
	* rwork, integer* info);

CLAPACK_API int ctrsen_(char* job, char* compq, logical* select, integer
	* n, complex* t, integer* ldt, complex* q, integer* ldq, complex* w,
	integer* m, real* s, real* sep, complex* work, integer* lwork,
	integer* info);

CLAPACK_API int ctrsna_(char* job, char* howmny, logical* select,
	integer* n, complex* t, integer* ldt, complex* vl, integer* ldvl,
	complex* vr, integer* ldvr, real* s, real* sep, integer* mm, integer*
	m, complex* work, integer* ldwork, real* rwork, integer* info);

CLAPACK_API int ctrsyl_(char* trana, char* tranb, integer* isgn, integer
	* m, integer* n, complex* a, integer* lda, complex* b, integer* ldb,
	complex* c__, integer* ldc, real* scale, integer* info);

CLAPACK_API int ctrti2_(char* uplo, char* diag, integer* n, complex* a,
	integer* lda, integer* info);

CLAPACK_API int ctrtri_(char* uplo, char* diag, integer* n, complex* a,
	integer* lda, integer* info);

CLAPACK_API int ctrtrs_(char* uplo, char* trans, char* diag, integer* n,
	integer* nrhs, complex* a, integer* lda, complex* b, integer* ldb,
	integer* info);

CLAPACK_API int ctrttf_(char* transr, char* uplo, integer* n, complex* a,
	integer* lda, complex* arf, integer* info);

CLAPACK_API int ctrttp_(char* uplo, integer* n, complex* a, integer* lda,
	complex* ap, integer* info);

CLAPACK_API int ctzrqf_(integer* m, integer* n, complex* a, integer* lda,
	complex* tau, integer* info);

CLAPACK_API int ctzrzf_(integer* m, integer* n, complex* a, integer* lda,
	complex* tau, complex* work, integer* lwork, integer* info);

CLAPACK_API int cung2l_(integer* m, integer* n, integer* k, complex* a,
	integer* lda, complex* tau, complex* work, integer* info);

CLAPACK_API int cung2r_(integer* m, integer* n, integer* k, complex* a,
	integer* lda, complex* tau, complex* work, integer* info);

CLAPACK_API int cungbr_(char* vect, integer* m, integer* n, integer* k,
	complex* a, integer* lda, complex* tau, complex* work, integer* lwork,
	integer* info);

CLAPACK_API int cunghr_(integer* n, integer* ilo, integer* ihi, complex*
	a, integer* lda, complex* tau, complex* work, integer* lwork, integer
	* info);

CLAPACK_API int cungl2_(integer* m, integer* n, integer* k, complex* a,
	integer* lda, complex* tau, complex* work, integer* info);

CLAPACK_API int cunglq_(integer* m, integer* n, integer* k, complex* a,
	integer* lda, complex* tau, complex* work, integer* lwork, integer*
	info);

CLAPACK_API int cungql_(integer* m, integer* n, integer* k, complex* a,
	integer* lda, complex* tau, complex* work, integer* lwork, integer*
	info);

CLAPACK_API int cungqr_(integer* m, integer* n, integer* k, complex* a,
	integer* lda, complex* tau, complex* work, integer* lwork, integer*
	info);

CLAPACK_API int cungr2_(integer* m, integer* n, integer* k, complex* a,
	integer* lda, complex* tau, complex* work, integer* info);

CLAPACK_API int cungrq_(integer* m, integer* n, integer* k, complex* a,
	integer* lda, complex* tau, complex* work, integer* lwork, integer*
	info);

CLAPACK_API int cungtr_(char* uplo, integer* n, complex* a, integer* lda,
	complex* tau, complex* work, integer* lwork, integer* info);

CLAPACK_API int cunm2l_(char* side, char* trans, integer* m, integer* n,
	integer* k, complex* a, integer* lda, complex* tau, complex* c__,
	integer* ldc, complex* work, integer* info);

CLAPACK_API int cunm2r_(char* side, char* trans, integer* m, integer* n,
	integer* k, complex* a, integer* lda, complex* tau, complex* c__,
	integer* ldc, complex* work, integer* info);

CLAPACK_API int cunmbr_(char* vect, char* side, char* trans, integer* m,
	integer* n, integer* k, complex* a, integer* lda, complex* tau,
	complex* c__, integer* ldc, complex* work, integer* lwork, integer*
	info);

CLAPACK_API int cunmhr_(char* side, char* trans, integer* m, integer* n,
	integer* ilo, integer* ihi, complex* a, integer* lda, complex* tau,
	complex* c__, integer* ldc, complex* work, integer* lwork, integer*
	info);

CLAPACK_API int cunml2_(char* side, char* trans, integer* m, integer* n,
	integer* k, complex* a, integer* lda, complex* tau, complex* c__,
	integer* ldc, complex* work, integer* info);

CLAPACK_API int cunmlq_(char* side, char* trans, integer* m, integer* n,
	integer* k, complex* a, integer* lda, complex* tau, complex* c__,
	integer* ldc, complex* work, integer* lwork, integer* info);

CLAPACK_API int cunmql_(char* side, char* trans, integer* m, integer* n,
	integer* k, complex* a, integer* lda, complex* tau, complex* c__,
	integer* ldc, complex* work, integer* lwork, integer* info);

CLAPACK_API int cunmqr_(char* side, char* trans, integer* m, integer* n,
	integer* k, complex* a, integer* lda, complex* tau, complex* c__,
	integer* ldc, complex* work, integer* lwork, integer* info);

CLAPACK_API int cunmr2_(char* side, char* trans, integer* m, integer* n,
	integer* k, complex* a, integer* lda, complex* tau, complex* c__,
	integer* ldc, complex* work, integer* info);

CLAPACK_API int cunmr3_(char* side, char* trans, integer* m, integer* n,
	integer* k, integer* l, complex* a, integer* lda, complex* tau,
	complex* c__, integer* ldc, complex* work, integer* info);

CLAPACK_API int cunmrq_(char* side, char* trans, integer* m, integer* n,
	integer* k, complex* a, integer* lda, complex* tau, complex* c__,
	integer* ldc, complex* work, integer* lwork, integer* info);

CLAPACK_API int cunmrz_(char* side, char* trans, integer* m, integer* n,
	integer* k, integer* l, complex* a, integer* lda, complex* tau,
	complex* c__, integer* ldc, complex* work, integer* lwork, integer*
	info);

CLAPACK_API int cunmtr_(char* side, char* uplo, char* trans, integer* m,
	integer* n, complex* a, integer* lda, complex* tau, complex* c__,
	integer* ldc, complex* work, integer* lwork, integer* info);

CLAPACK_API int cupgtr_(char* uplo, integer* n, complex* ap, complex*
	tau, complex* q, integer* ldq, complex* work, integer* info);

CLAPACK_API int cupmtr_(char* side, char* uplo, char* trans, integer* m,
	integer* n, complex* ap, complex* tau, complex* c__, integer* ldc,
	complex* work, integer* info);

CLAPACK_API
integer icmax1_(integer* n, complex* cx, integer* incx);

CLAPACK_API
integer ilaclc_(integer* m, integer* n, complex* a, integer* lda);

CLAPACK_API
integer ilaclr_(integer* m, integer* n, complex* a, integer* lda);


//}} finished

#pragma endregion

#pragma region CXLASRC -- Single precision complex LAPACK routines using extra precision (clax) - ok
//#

/*
	cgbsvxx.c cgbrfsx.c
	cgerfsx.c cgesvxx.c
	chesvxx.c cherfsx.c
	cposvxx.c cporfsx.c
	csysvxx.c csyrfsx.c
	
	cla_gbamv.c
	cla_gbrcond_c.c
	cla_gbrcond_x.c
	cla_gbrpvgrw.c
	cla_gbrfsx_extended.c
	
	cla_geamv.c
	cla_gercond_c.c
	cla_gercond_x.c
	cla_gerfsx_extended.c 
	
	cla_herfsx_extended.c
	cla_heamv.c
	cla_hercond_c.c
	cla_hercond_x.c
	cla_herpvgrw.c
	
	cla_lin_berr.c
	
	cla_rpvgrw.c
	
	cla_syamv.c
	cla_syrcond_c.c
	cla_syrcond_x.c
	cla_syrfsx_extended.c
	cla_syrpvgrw.c
	
	cla_porfsx_extended.c
	cla_porcond_c.c
	cla_porcond_x.c
	cla_porpvgrw.c
	
	cla_wwaddw.c
*/

CLAPACK_API int clarscl2_(integer* m, integer* n, real* d__, complex* x,
	integer* ldx);

CLAPACK_API int clascl2_(integer* m, integer* n, real* d__, complex* x,
	integer* ldx);



/*
int cgbrfsx_(char* trans, char* equed, integer* n, integer*
	kl, integer* ku, integer* nrhs, complex* ab, integer* ldab, complex*
	afb, integer* ldafb, integer* ipiv, real* r__, real* c__, complex* b,
	integer* ldb, complex* x, integer* ldx, real* rcond, real* berr,
	integer* n_err_bnds__, real* err_bnds_norm__, real* err_bnds_comp__,
	integer* nparams, real* params, complex* work, real* rwork, integer*
	info);

int cgbsvxx_(char* fact, char* trans, integer* n, integer*
	kl, integer* ku, integer* nrhs, complex* ab, integer* ldab, complex*
	afb, integer* ldafb, integer* ipiv, char* equed, real* r__, real* c__,
	complex* b, integer* ldb, complex* x, integer* ldx, real* rcond,
	real* rpvgrw, real* berr, integer* n_err_bnds__, real*
	err_bnds_norm__, real* err_bnds_comp__, integer* nparams, real*
	params, complex* work, real* rwork, integer* info);

int cgerfsx_(char* trans, char* equed, integer* n, integer*
	nrhs, complex* a, integer* lda, complex* af, integer* ldaf, integer*
	ipiv, real* r__, real* c__, complex* b, integer* ldb, complex* x,
	integer* ldx, real* rcond, real* berr, integer* n_err_bnds__, real*
	err_bnds_norm__, real* err_bnds_comp__, integer* nparams, real*
	params, complex* work, real* rwork, integer* info);

int cgesvxx_(char* fact, char* trans, integer* n, integer*
	nrhs, complex* a, integer* lda, complex* af, integer* ldaf, integer*
	ipiv, char* equed, real* r__, real* c__, complex* b, integer* ldb,
	complex* x, integer* ldx, real* rcond, real* rpvgrw, real* berr,
	integer* n_err_bnds__, real* err_bnds_norm__, real* err_bnds_comp__,
	integer* nparams, real* params, complex* work, real* rwork, integer*
	info);

int cherfsx_(char* uplo, char* equed, integer* n, integer*
	nrhs, complex* a, integer* lda, complex* af, integer* ldaf, integer*
	ipiv, real* s, complex* b, integer* ldb, complex* x, integer* ldx,
	real* rcond, real* berr, integer* n_err_bnds__, real* err_bnds_norm__,
	real* err_bnds_comp__, integer* nparams, real* params, complex* work,
	real* rwork, integer* info);


int chesvxx_(char* fact, char* uplo, integer* n, integer*
	nrhs, complex* a, integer* lda, complex* af, integer* ldaf, integer*
	ipiv, char* equed, real* s, complex* b, integer* ldb, complex* x,
	integer* ldx, real* rcond, real* rpvgrw, real* berr, integer*
	n_err_bnds__, real* err_bnds_norm__, real* err_bnds_comp__, integer*
	nparams, real* params, complex* work, real* rwork, integer* info);


int cla_gbamv__(integer *trans, integer *m, integer *n,
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

int cla_gbrfsx_extended__(integer *prec_type__, integer *
	trans_type__, integer *n, integer *kl, integer *ku, integer *nrhs,
	complex *ab, integer *ldab, complex *afb, integer *ldafb, integer *
	ipiv, logical *colequ, real *c__, complex *b, integer *ldb, complex *
	y, integer *ldy, real *berr_out__, integer *n_norms__, real *errs_n__,
	 real *errs_c__, complex *res, real *ayb, complex *dy, complex *
	y_tail__, real *rcond, integer *ithresh, real *rthresh, real *dz_ub__,
	 logical *ignore_cwise__, integer *info);

doublereal cla_gbrpvgrw__(integer *n, integer *kl, integer *ku, integer *
	ncols, complex *ab, integer *ldab, complex *afb, integer *ldafb);

int cla_geamv__(integer *trans, integer *m, integer *n, real
	*alpha, complex *a, integer *lda, complex *x, integer *incx, real *
	beta, real *y, integer *incy);

doublereal cla_gercond_c__(char *trans, integer *n, complex *a, integer *lda,
	complex *af, integer *ldaf, integer *ipiv, real *c__, logical *capply,
	 integer *info, complex *work, real *rwork, ftnlen trans_len);

doublereal cla_gercond_x__(char *trans, integer *n, complex *a, integer *lda,
	complex *af, integer *ldaf, integer *ipiv, complex *x, integer *info,
	complex *work, real *rwork, ftnlen trans_len);

int cla_gerfsx_extended__(integer *prec_type__, integer *
	trans_type__, integer *n, integer *nrhs, complex *a, integer *lda,
	complex *af, integer *ldaf, integer *ipiv, logical *colequ, real *c__,
	 complex *b, integer *ldb, complex *y, integer *ldy, real *berr_out__,
	 integer *n_norms__, real *errs_n__, real *errs_c__, complex *res,
	real *ayb, complex *dy, complex *y_tail__, real *rcond, integer *
	ithresh, real *rthresh, real *dz_ub__, logical *ignore_cwise__,
	integer *info);

int cla_heamv__(integer *uplo, integer *n, real *alpha,
	complex *a, integer *lda, complex *x, integer *incx, real *beta, real
	*y, integer *incy);

doublereal cla_hercond_c__(char *uplo, integer *n, complex *a, integer *lda,
	complex *af, integer *ldaf, integer *ipiv, real *c__, logical *capply,
	 integer *info, complex *work, real *rwork, ftnlen uplo_len);

doublereal cla_hercond_x__(char *uplo, integer *n, complex *a, integer *lda,
	complex *af, integer *ldaf, integer *ipiv, complex *x, integer *info,
	complex *work, real *rwork, ftnlen uplo_len);

int cla_herfsx_extended__(integer *prec_type__, char *uplo,
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

int cla_lin_berr__(integer *n, integer *nz, integer *nrhs,
	complex *res, real *ayb, real *berr);

doublereal cla_porcond_c__(char *uplo, integer *n, complex *a, integer *lda,
	complex *af, integer *ldaf, real *c__, logical *capply, integer *info,
	 complex *work, real *rwork, ftnlen uplo_len);

doublereal cla_porcond_x__(char *uplo, integer *n, complex *a, integer *lda,
	complex *af, integer *ldaf, complex *x, integer *info, complex *work,
	real *rwork, ftnlen uplo_len);

int cla_porfsx_extended__(integer *prec_type__, char *uplo,
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

int cla_syamv__(integer *uplo, integer *n, real *alpha,
	complex *a, integer *lda, complex *x, integer *incx, real *beta, real
	*y, integer *incy);

doublereal cla_syrcond_c__(char *uplo, integer *n, complex *a, integer *lda,
	complex *af, integer *ldaf, integer *ipiv, real *c__, logical *capply,
	 integer *info, complex *work, real *rwork, ftnlen uplo_len);

doublereal cla_syrcond_x__(char *uplo, integer *n, complex *a, integer *lda,
	complex *af, integer *ldaf, integer *ipiv, complex *x, integer *info,
	complex *work, real *rwork, ftnlen uplo_len);

int cla_syrfsx_extended__(integer *prec_type__, char *uplo,
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

int cla_wwaddw__(integer *n, complex *x, complex *y, complex *w);


int cporfsx_(char *uplo, char *equed, integer *n, integer *
	nrhs, complex *a, integer *lda, complex *af, integer *ldaf, real *s,
	complex *b, integer *ldb, complex *x, integer *ldx, real *rcond, real
	*berr, integer *n_err_bnds__, real *err_bnds_norm__, real *
	err_bnds_comp__, integer *nparams, real *params, complex *work, real *
	rwork, integer *info);

int cposvxx_(char *fact, char *uplo, integer *n, integer *
	nrhs, complex *a, integer *lda, complex *af, integer *ldaf, char *
	equed, real *s, complex *b, integer *ldb, complex *x, integer *ldx,
	real *rcond, real *rpvgrw, real *berr, integer *n_err_bnds__, real *
	err_bnds_norm__, real *err_bnds_comp__, integer *nparams, real *
	params, complex *work, real *rwork, integer *info);

int csysvxx_(char *fact, char *uplo, integer *n, integer *
	nrhs, complex *a, integer *lda, complex *af, integer *ldaf, integer *
	ipiv, char *equed, real *s, complex *b, integer *ldb, complex *x,
	integer *ldx, real *rcond, real *rpvgrw, real *berr, integer *
	n_err_bnds__, real *err_bnds_norm__, real *err_bnds_comp__, integer *
	nparams, real *params, complex *work, real *rwork, integer *info);






*/

#pragma endregion


// D
#pragma region DLASRC -- Double precision real LAPACK routines (dla) - ok


CLAPACK_API
int dgbbrd_(char* vect, integer* m, integer* n, integer* ncc,
	integer* kl, integer* ku, doublereal* ab, integer* ldab, doublereal*
	d__, doublereal* e, doublereal* q, integer* ldq, doublereal* pt,
	integer* ldpt, doublereal* c__, integer* ldc, doublereal* work,
	integer* info);


CLAPACK_API
int dgbcon_(char* norm, integer* n, integer* kl, integer* ku,
	doublereal* ab, integer* ldab, integer* ipiv, doublereal* anorm,
	doublereal* rcond, doublereal* work, integer* iwork, integer* info);


CLAPACK_API
int dgbequ_(integer* m, integer* n, integer* kl, integer* ku,
	doublereal* ab, integer* ldab, doublereal* r__, doublereal* c__,
	doublereal* rowcnd, doublereal* colcnd, doublereal* amax, integer*
	info);


CLAPACK_API
int dgbequb_(integer* m, integer* n, integer* kl, integer*
	ku, doublereal* ab, integer* ldab, doublereal* r__, doublereal* c__,
	doublereal* rowcnd, doublereal* colcnd, doublereal* amax, integer*
	info);


CLAPACK_API
int dgbrfs_(char* trans, integer* n, integer* kl, integer*
	ku, integer* nrhs, doublereal* ab, integer* ldab, doublereal* afb,
	integer* ldafb, integer* ipiv, doublereal* b, integer* ldb,
	doublereal* x, integer* ldx, doublereal* ferr, doublereal* berr,
	doublereal* work, integer* iwork, integer* info);


CLAPACK_API
int dgbsv_(integer* n, integer* kl, integer* ku, integer*
	nrhs, doublereal* ab, integer* ldab, integer* ipiv, doublereal* b,
	integer* ldb, integer* info);


CLAPACK_API
int dgbsvx_(char* fact, char* trans, integer* n, integer* kl,
	integer* ku, integer* nrhs, doublereal* ab, integer* ldab,
	doublereal* afb, integer* ldafb, integer* ipiv, char* equed,
	doublereal* r__, doublereal* c__, doublereal* b, integer* ldb,
	doublereal* x, integer* ldx, doublereal* rcond, doublereal* ferr,
	doublereal* berr, doublereal* work, integer* iwork, integer* info);


CLAPACK_API
int dgbtf2_(integer* m, integer* n, integer* kl, integer* ku,
	doublereal* ab, integer* ldab, integer* ipiv, integer* info);


CLAPACK_API
int dgbtrf_(integer* m, integer* n, integer* kl, integer* ku,
	doublereal* ab, integer* ldab, integer* ipiv, integer* info);


CLAPACK_API
int dgbtrs_(char* trans, integer* n, integer* kl, integer*
	ku, integer* nrhs, doublereal* ab, integer* ldab, integer* ipiv,
	doublereal* b, integer* ldb, integer* info);


CLAPACK_API
int dgebak_(char* job, char* side, integer* n, integer* ilo,
	integer* ihi, doublereal* scale, integer* m, doublereal* v, integer*
	ldv, integer* info);


CLAPACK_API
int dgebal_(char* job, integer* n, doublereal* a, integer*
	lda, integer* ilo, integer* ihi, doublereal* scale, integer* info);


CLAPACK_API
int dgebd2_(integer* m, integer* n, doublereal* a, integer*
	lda, doublereal* d__, doublereal* e, doublereal* tauq, doublereal*
	taup, doublereal* work, integer* info);


CLAPACK_API
int dgebrd_(integer* m, integer* n, doublereal* a, integer*
	lda, doublereal* d__, doublereal* e, doublereal* tauq, doublereal*
	taup, doublereal* work, integer* lwork, integer* info);


CLAPACK_API
int dgecon_(char* norm, integer* n, doublereal* a, integer*
	lda, doublereal* anorm, doublereal* rcond, doublereal* work, integer*
	iwork, integer* info);


CLAPACK_API
int dgeequ_(integer* m, integer* n, doublereal* a, integer*
	lda, doublereal* r__, doublereal* c__, doublereal* rowcnd, doublereal
	* colcnd, doublereal* amax, integer* info);


CLAPACK_API
int dgeequb_(integer* m, integer* n, doublereal* a, integer*
	lda, doublereal* r__, doublereal* c__, doublereal* rowcnd, doublereal
	* colcnd, doublereal* amax, integer* info);


CLAPACK_API
int dgees_(char* jobvs, char* sort, L_fp select, integer* n,
	doublereal* a, integer* lda, integer* sdim, doublereal* wr,
	doublereal* wi, doublereal* vs, integer* ldvs, doublereal* work,
	integer* lwork, logical* bwork, integer* info);


CLAPACK_API
int dgeesx_(char* jobvs, char* sort, L_fp select, char*
	sense, integer* n, doublereal* a, integer* lda, integer* sdim,
	doublereal* wr, doublereal* wi, doublereal* vs, integer* ldvs,
	doublereal* rconde, doublereal* rcondv, doublereal* work, integer*
	lwork, integer* iwork, integer* liwork, logical* bwork, integer* info);


CLAPACK_API
int dgeev_(char* jobvl, char* jobvr, integer* n, doublereal*
	a, integer* lda, doublereal* wr, doublereal* wi, doublereal* vl,
	integer* ldvl, doublereal* vr, integer* ldvr, doublereal* work,
	integer* lwork, integer* info);


CLAPACK_API
int dgeevx_(char* balanc, char* jobvl, char* jobvr, char*
	sense, integer* n, doublereal* a, integer* lda, doublereal* wr,
	doublereal* wi, doublereal* vl, integer* ldvl, doublereal* vr,
	integer* ldvr, integer* ilo, integer* ihi, doublereal* scale,
	doublereal* abnrm, doublereal* rconde, doublereal* rcondv, doublereal
	* work, integer* lwork, integer* iwork, integer* info);


CLAPACK_API
int dgegs_(char* jobvsl, char* jobvsr, integer* n,
	doublereal* a, integer* lda, doublereal* b, integer* ldb, doublereal*
	alphar, doublereal* alphai, doublereal* beta, doublereal* vsl,
	integer* ldvsl, doublereal* vsr, integer* ldvsr, doublereal* work,
	integer* lwork, integer* info);


CLAPACK_API
int dgegv_(char* jobvl, char* jobvr, integer* n, doublereal*
	a, integer* lda, doublereal* b, integer* ldb, doublereal* alphar,
	doublereal* alphai, doublereal* beta, doublereal* vl, integer* ldvl,
	doublereal* vr, integer* ldvr, doublereal* work, integer* lwork,
	integer* info);


CLAPACK_API
int dgehd2_(integer* n, integer* ilo, integer* ihi,
	doublereal* a, integer* lda, doublereal* tau, doublereal* work,
	integer* info);


CLAPACK_API
int dgehrd_(integer* n, integer* ilo, integer* ihi,
	doublereal* a, integer* lda, doublereal* tau, doublereal* work,
	integer* lwork, integer* info);


CLAPACK_API
int dgejsv_(char* joba, char* jobu, char* jobv, char* jobr,
	char* jobt, char* jobp, integer* m, integer* n, doublereal* a,
	integer* lda, doublereal* sva, doublereal* u, integer* ldu,
	doublereal* v, integer* ldv, doublereal* work, integer* lwork,
	integer* iwork, integer* info);


CLAPACK_API
int dgelq2_(integer* m, integer* n, doublereal* a, integer*
	lda, doublereal* tau, doublereal* work, integer* info);


CLAPACK_API
int dgelqf_(integer* m, integer* n, doublereal* a, integer*
	lda, doublereal* tau, doublereal* work, integer* lwork, integer* info);


CLAPACK_API
int dgels_(char* trans, integer* m, integer* n, integer*
	nrhs, doublereal* a, integer* lda, doublereal* b, integer* ldb,
	doublereal* work, integer* lwork, integer* info);


CLAPACK_API
int dgelsd_(integer* m, integer* n, integer* nrhs,
	doublereal* a, integer* lda, doublereal* b, integer* ldb, doublereal*
	s, doublereal* rcond, integer* rank, doublereal* work, integer* lwork,
	integer* iwork, integer* info);


CLAPACK_API
int dgelss_(integer* m, integer* n, integer* nrhs,
	doublereal* a, integer* lda, doublereal* b, integer* ldb, doublereal*
	s, doublereal* rcond, integer* rank, doublereal* work, integer* lwork,
	integer* info);


CLAPACK_API
int dgelsx_(integer* m, integer* n, integer* nrhs,
	doublereal* a, integer* lda, doublereal* b, integer* ldb, integer*
	jpvt, doublereal* rcond, integer* rank, doublereal* work, integer*
	info);


CLAPACK_API
int dgelsy_(integer* m, integer* n, integer* nrhs,
	doublereal* a, integer* lda, doublereal* b, integer* ldb, integer*
	jpvt, doublereal* rcond, integer* rank, doublereal* work, integer*
	lwork, integer* info);


CLAPACK_API
int dgeql2_(integer* m, integer* n, doublereal* a, integer*
	lda, doublereal* tau, doublereal* work, integer* info);


CLAPACK_API
int dgeqlf_(integer* m, integer* n, doublereal* a, integer*
	lda, doublereal* tau, doublereal* work, integer* lwork, integer* info);


CLAPACK_API
int dgeqp3_(integer* m, integer* n, doublereal* a, integer*
	lda, integer* jpvt, doublereal* tau, doublereal* work, integer* lwork,
	integer* info);


CLAPACK_API
int dgeqpf_(integer* m, integer* n, doublereal* a, integer*
	lda, integer* jpvt, doublereal* tau, doublereal* work, integer* info);


CLAPACK_API
int dgeqr2_(integer* m, integer* n, doublereal* a, integer*
	lda, doublereal* tau, doublereal* work, integer* info);


CLAPACK_API
int dgeqrf_(integer* m, integer* n, doublereal* a, integer*
	lda, doublereal* tau, doublereal* work, integer* lwork, integer* info);


CLAPACK_API
int dgerfs_(char* trans, integer* n, integer* nrhs,
	doublereal* a, integer* lda, doublereal* af, integer* ldaf, integer*
	ipiv, doublereal* b, integer* ldb, doublereal* x, integer* ldx,
	doublereal* ferr, doublereal* berr, doublereal* work, integer* iwork,
	integer* info);


CLAPACK_API
int dgerq2_(integer* m, integer* n, doublereal* a, integer*
	lda, doublereal* tau, doublereal* work, integer* info);


CLAPACK_API
int dgerqf_(integer* m, integer* n, doublereal* a, integer*
	lda, doublereal* tau, doublereal* work, integer* lwork, integer* info);


CLAPACK_API
int dgesc2_(integer* n, doublereal* a, integer* lda,
	doublereal* rhs, integer* ipiv, integer* jpiv, doublereal* scale);


CLAPACK_API
int dgesdd_(char* jobz, integer* m, integer* n, doublereal*
	a, integer* lda, doublereal* s, doublereal* u, integer* ldu,
	doublereal* vt, integer* ldvt, doublereal* work, integer* lwork,
	integer* iwork, integer* info);


CLAPACK_API
int dgesv_(integer* n, integer* nrhs, doublereal* a, integer
	* lda, integer* ipiv, doublereal* b, integer* ldb, integer* info);


CLAPACK_API
int dgesvd_(char* jobu, char* jobvt, integer* m, integer* n,
	doublereal* a, integer* lda, doublereal* s, doublereal* u, integer*
	ldu, doublereal* vt, integer* ldvt, doublereal* work, integer* lwork,
	integer* info);


CLAPACK_API
int dgesvj_(char* joba, char* jobu, char* jobv, integer* m,
	integer* n, doublereal* a, integer* lda, doublereal* sva, integer* mv,
	doublereal* v, integer* ldv, doublereal* work, integer* lwork,
	integer* info);


CLAPACK_API
int dgesvx_(char* fact, char* trans, integer* n, integer*
	nrhs, doublereal* a, integer* lda, doublereal* af, integer* ldaf,
	integer* ipiv, char* equed, doublereal* r__, doublereal* c__,
	doublereal* b, integer* ldb, doublereal* x, integer* ldx, doublereal*
	rcond, doublereal* ferr, doublereal* berr, doublereal* work, integer*
	iwork, integer* info);


CLAPACK_API
int dgetc2_(integer* n, doublereal* a, integer* lda, integer
	* ipiv, integer* jpiv, integer* info);


CLAPACK_API
int dgetf2_(integer* m, integer* n, doublereal* a, integer*
	lda, integer* ipiv, integer* info);


CLAPACK_API
int dgetrf_(integer* m, integer* n, doublereal* a, integer*
	lda, integer* ipiv, integer* info);


CLAPACK_API
int dgetri_(integer* n, doublereal* a, integer* lda, integer
	* ipiv, doublereal* work, integer* lwork, integer* info);


CLAPACK_API
int dgetrs_(char* trans, integer* n, integer* nrhs,
	doublereal* a, integer* lda, integer* ipiv, doublereal* b, integer*
	ldb, integer* info);


CLAPACK_API
int dggbak_(char* job, char* side, integer* n, integer* ilo,
	integer* ihi, doublereal* lscale, doublereal* rscale, integer* m,
	doublereal* v, integer* ldv, integer* info);


CLAPACK_API
int dggbal_(char* job, integer* n, doublereal* a, integer*
	lda, doublereal* b, integer* ldb, integer* ilo, integer* ihi,
	doublereal* lscale, doublereal* rscale, doublereal* work, integer*
	info);


CLAPACK_API
int dgges_(char* jobvsl, char* jobvsr, char* sort, L_fp
	selctg, integer* n, doublereal* a, integer* lda, doublereal* b,
	integer* ldb, integer* sdim, doublereal* alphar, doublereal* alphai,
	doublereal* beta, doublereal* vsl, integer* ldvsl, doublereal* vsr,
	integer* ldvsr, doublereal* work, integer* lwork, logical* bwork,
	integer* info);


CLAPACK_API
int dggesx_(char* jobvsl, char* jobvsr, char* sort, L_fp
	selctg, char* sense, integer* n, doublereal* a, integer* lda,
	doublereal* b, integer* ldb, integer* sdim, doublereal* alphar,
	doublereal* alphai, doublereal* beta, doublereal* vsl, integer* ldvsl,
	doublereal* vsr, integer* ldvsr, doublereal* rconde, doublereal*
	rcondv, doublereal* work, integer* lwork, integer* iwork, integer*
	liwork, logical* bwork, integer* info);


CLAPACK_API
int dggev_(char* jobvl, char* jobvr, integer* n, doublereal*
	a, integer* lda, doublereal* b, integer* ldb, doublereal* alphar,
	doublereal* alphai, doublereal* beta, doublereal* vl, integer* ldvl,
	doublereal* vr, integer* ldvr, doublereal* work, integer* lwork,
	integer* info);


CLAPACK_API
int dggevx_(char* balanc, char* jobvl, char* jobvr, char*
	sense, integer* n, doublereal* a, integer* lda, doublereal* b,
	integer* ldb, doublereal* alphar, doublereal* alphai, doublereal*
	beta, doublereal* vl, integer* ldvl, doublereal* vr, integer* ldvr,
	integer* ilo, integer* ihi, doublereal* lscale, doublereal* rscale,
	doublereal* abnrm, doublereal* bbnrm, doublereal* rconde, doublereal*
	rcondv, doublereal* work, integer* lwork, integer* iwork, logical*
	bwork, integer* info);


CLAPACK_API
int dggglm_(integer* n, integer* m, integer* p, doublereal*
	a, integer* lda, doublereal* b, integer* ldb, doublereal* d__,
	doublereal* x, doublereal* y, doublereal* work, integer* lwork,
	integer* info);


CLAPACK_API
int dgghrd_(char* compq, char* compz, integer* n, integer*
	ilo, integer* ihi, doublereal* a, integer* lda, doublereal* b,
	integer* ldb, doublereal* q, integer* ldq, doublereal* z__, integer*
	ldz, integer* info);


CLAPACK_API
int dgglse_(integer* m, integer* n, integer* p, doublereal*
	a, integer* lda, doublereal* b, integer* ldb, doublereal* c__,
	doublereal* d__, doublereal* x, doublereal* work, integer* lwork,
	integer* info);


CLAPACK_API
int dggqrf_(integer* n, integer* m, integer* p, doublereal*
	a, integer* lda, doublereal* taua, doublereal* b, integer* ldb,
	doublereal* taub, doublereal* work, integer* lwork, integer* info);


CLAPACK_API
int dggrqf_(integer* m, integer* p, integer* n, doublereal*
	a, integer* lda, doublereal* taua, doublereal* b, integer* ldb,
	doublereal* taub, doublereal* work, integer* lwork, integer* info);


CLAPACK_API
int dggsvd_(char* jobu, char* jobv, char* jobq, integer* m,
	integer* n, integer* p, integer* k, integer* l, doublereal* a,
	integer* lda, doublereal* b, integer* ldb, doublereal* alpha,
	doublereal* beta, doublereal* u, integer* ldu, doublereal* v, integer
	* ldv, doublereal* q, integer* ldq, doublereal* work, integer* iwork,
	integer* info);


CLAPACK_API
int dggsvp_(char* jobu, char* jobv, char* jobq, integer* m,
	integer* p, integer* n, doublereal* a, integer* lda, doublereal* b,
	integer* ldb, doublereal* tola, doublereal* tolb, integer* k, integer
	* l, doublereal* u, integer* ldu, doublereal* v, integer* ldv,
	doublereal* q, integer* ldq, integer* iwork, doublereal* tau,
	doublereal* work, integer* info);


CLAPACK_API
int dgsvj0_(char* jobv, integer* m, integer* n, doublereal*
	a, integer* lda, doublereal* d__, doublereal* sva, integer* mv,
	doublereal* v, integer* ldv, doublereal* eps, doublereal* sfmin,
	doublereal* tol, integer* nsweep, doublereal* work, integer* lwork,
	integer* info);


CLAPACK_API
int dgsvj1_(char* jobv, integer* m, integer* n, integer* n1,
	doublereal* a, integer* lda, doublereal* d__, doublereal* sva,
	integer* mv, doublereal* v, integer* ldv, doublereal* eps, doublereal
	* sfmin, doublereal* tol, integer* nsweep, doublereal* work, integer*
	lwork, integer* info);


CLAPACK_API
int dgtcon_(char* norm, integer* n, doublereal* dl,
	doublereal* d__, doublereal* du, doublereal* du2, integer* ipiv,
	doublereal* anorm, doublereal* rcond, doublereal* work, integer*
	iwork, integer* info);


CLAPACK_API
int dgtrfs_(char* trans, integer* n, integer* nrhs,
	doublereal* dl, doublereal* d__, doublereal* du, doublereal* dlf,
	doublereal* df, doublereal* duf, doublereal* du2, integer* ipiv,
	doublereal* b, integer* ldb, doublereal* x, integer* ldx, doublereal*
	ferr, doublereal* berr, doublereal* work, integer* iwork, integer*
	info);


CLAPACK_API
int dgtsv_(integer* n, integer* nrhs, doublereal* dl,
	doublereal* d__, doublereal* du, doublereal* b, integer* ldb, integer
	* info);


CLAPACK_API
int dgtsvx_(char* fact, char* trans, integer* n, integer*
	nrhs, doublereal* dl, doublereal* d__, doublereal* du, doublereal*
	dlf, doublereal* df, doublereal* duf, doublereal* du2, integer* ipiv,
	doublereal* b, integer* ldb, doublereal* x, integer* ldx, doublereal*
	rcond, doublereal* ferr, doublereal* berr, doublereal* work, integer*
	iwork, integer* info);


CLAPACK_API
int dgttrf_(integer* n, doublereal* dl, doublereal* d__,
	doublereal* du, doublereal* du2, integer* ipiv, integer* info);


CLAPACK_API
int dgttrs_(char* trans, integer* n, integer* nrhs,
	doublereal* dl, doublereal* d__, doublereal* du, doublereal* du2,
	integer* ipiv, doublereal* b, integer* ldb, integer* info);


CLAPACK_API
int dgtts2_(integer* itrans, integer* n, integer* nrhs,
	doublereal* dl, doublereal* d__, doublereal* du, doublereal* du2,
	integer* ipiv, doublereal* b, integer* ldb);


CLAPACK_API
int dhgeqz_(char* job, char* compq, char* compz, integer* n,
	integer* ilo, integer* ihi, doublereal* h__, integer* ldh, doublereal
	* t, integer* ldt, doublereal* alphar, doublereal* alphai, doublereal*
	beta, doublereal* q, integer* ldq, doublereal* z__, integer* ldz,
	doublereal* work, integer* lwork, integer* info);


CLAPACK_API
int dhsein_(char* side, char* eigsrc, char* initv, logical*
	select, integer* n, doublereal* h__, integer* ldh, doublereal* wr,
	doublereal* wi, doublereal* vl, integer* ldvl, doublereal* vr,
	integer* ldvr, integer* mm, integer* m, doublereal* work, integer*
	ifaill, integer* ifailr, integer* info);


CLAPACK_API
int dhseqr_(char* job, char* compz, integer* n, integer* ilo,
	integer* ihi, doublereal* h__, integer* ldh, doublereal* wr,
	doublereal* wi, doublereal* z__, integer* ldz, doublereal* work,
	integer* lwork, integer* info);


CLAPACK_API
int dlabrd_(integer* m, integer* n, integer* nb, doublereal*
	a, integer* lda, doublereal* d__, doublereal* e, doublereal* tauq,
	doublereal* taup, doublereal* x, integer* ldx, doublereal* y, integer
	* ldy);


CLAPACK_API
int dlacn2_(integer* n, doublereal* v, doublereal* x,
	integer* isgn, doublereal* est, integer* kase, integer* isave);


CLAPACK_API
int dlacon_(integer* n, doublereal* v, doublereal* x,
	integer* isgn, doublereal* est, integer* kase);


CLAPACK_API
int dlaein_(logical* rightv, logical* noinit, integer* n,
	doublereal* h__, integer* ldh, doublereal* wr, doublereal* wi,
	doublereal* vr, doublereal* vi, doublereal* b, integer* ldb,
	doublereal* work, doublereal* eps3, doublereal* smlnum, doublereal*
	bignum, integer* info);


CLAPACK_API
int dlaexc_(logical* wantq, integer* n, doublereal* t,
	integer* ldt, doublereal* q, integer* ldq, integer* j1, integer* n1,
	integer* n2, doublereal* work, integer* info);


CLAPACK_API
int dlag2_(doublereal* a, integer* lda, doublereal* b,
	integer* ldb, doublereal* safmin, doublereal* scale1, doublereal*
	scale2, doublereal* wr1, doublereal* wr2, doublereal* wi);


CLAPACK_API
int dlag2s_(integer* m, integer* n, doublereal* a, integer*
	lda, real* sa, integer* ldsa, integer* info);


CLAPACK_API
int dlags2_(logical* upper, doublereal* a1, doublereal* a2,
	doublereal* a3, doublereal* b1, doublereal* b2, doublereal* b3,
	doublereal* csu, doublereal* snu, doublereal* csv, doublereal* snv,
	doublereal* csq, doublereal* snq);


CLAPACK_API
int dlagtm_(char* trans, integer* n, integer* nrhs,
	doublereal* alpha, doublereal* dl, doublereal* d__, doublereal* du,
	doublereal* x, integer* ldx, doublereal* beta, doublereal* b, integer
	* ldb);


CLAPACK_API
int dlagv2_(doublereal* a, integer* lda, doublereal* b,
	integer* ldb, doublereal* alphar, doublereal* alphai, doublereal*
	beta, doublereal* csl, doublereal* snl, doublereal* csr, doublereal*
	snr);


CLAPACK_API
int dlahqr_(logical* wantt, logical* wantz, integer* n,
	integer* ilo, integer* ihi, doublereal* h__, integer* ldh, doublereal
	* wr, doublereal* wi, integer* iloz, integer* ihiz, doublereal* z__,
	integer* ldz, integer* info);


CLAPACK_API
int dlahr2_(integer* n, integer* k, integer* nb, doublereal*
	a, integer* lda, doublereal* tau, doublereal* t, integer* ldt,
	doublereal* y, integer* ldy);


CLAPACK_API
int dlahrd_(integer* n, integer* k, integer* nb, doublereal*
	a, integer* lda, doublereal* tau, doublereal* t, integer* ldt,
	doublereal* y, integer* ldy);


CLAPACK_API
int dlaic1_(integer* job, integer* j, doublereal* x,
	doublereal* sest, doublereal* w, doublereal* gamma, doublereal*
	sestpr, doublereal* s, doublereal* c__);


CLAPACK_API
int dlaln2_(logical* ltrans, integer* na, integer* nw,
	doublereal* smin, doublereal* ca, doublereal* a, integer* lda,
	doublereal* d1, doublereal* d2, doublereal* b, integer* ldb,
	doublereal* wr, doublereal* wi, doublereal* x, integer* ldx,
	doublereal* scale, doublereal* xnorm, integer* info);


CLAPACK_API
int dlals0_(integer* icompq, integer* nl, integer* nr,
	integer* sqre, integer* nrhs, doublereal* b, integer* ldb, doublereal
	* bx, integer* ldbx, integer* perm, integer* givptr, integer* givcol,
	integer* ldgcol, doublereal* givnum, integer* ldgnum, doublereal*
	poles, doublereal* difl, doublereal* difr, doublereal* z__, integer*
	k, doublereal* c__, doublereal* s, doublereal* work, integer* info);


CLAPACK_API
int dlalsa_(integer* icompq, integer* smlsiz, integer* n,
	integer* nrhs, doublereal* b, integer* ldb, doublereal* bx, integer*
	ldbx, doublereal* u, integer* ldu, doublereal* vt, integer* k,
	doublereal* difl, doublereal* difr, doublereal* z__, doublereal*
	poles, integer* givptr, integer* givcol, integer* ldgcol, integer*
	perm, doublereal* givnum, doublereal* c__, doublereal* s, doublereal*
	work, integer* iwork, integer* info);


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


CLAPACK_API
int dlanv2_(doublereal* a, doublereal* b, doublereal* c__,
	doublereal* d__, doublereal* rt1r, doublereal* rt1i, doublereal* rt2r,
	doublereal* rt2i, doublereal* cs, doublereal* sn);


CLAPACK_API
int dlapll_(integer* n, doublereal* x, integer* incx,
	doublereal* y, integer* incy, doublereal* ssmin);


CLAPACK_API
int dlapmt_(logical* forwrd, integer* m, integer* n,
	doublereal* x, integer* ldx, integer* k);


CLAPACK_API
int dlaqgb_(integer* m, integer* n, integer* kl, integer* ku,
	doublereal* ab, integer* ldab, doublereal* r__, doublereal* c__,
	doublereal* rowcnd, doublereal* colcnd, doublereal* amax, char* equed);


CLAPACK_API
int dlaqge_(integer* m, integer* n, doublereal* a, integer*
	lda, doublereal* r__, doublereal* c__, doublereal* rowcnd, doublereal
	* colcnd, doublereal* amax, char* equed);


CLAPACK_API
int dlaqp2_(integer* m, integer* n, integer* offset,
	doublereal* a, integer* lda, integer* jpvt, doublereal* tau,
	doublereal* vn1, doublereal* vn2, doublereal* work);


CLAPACK_API
int dlaqps_(integer* m, integer* n, integer* offset, integer
	* nb, integer* kb, doublereal* a, integer* lda, integer* jpvt,
	doublereal* tau, doublereal* vn1, doublereal* vn2, doublereal* auxv,
	doublereal* f, integer* ldf);


CLAPACK_API
int dlaqr0_(logical* wantt, logical* wantz, integer* n,
	integer* ilo, integer* ihi, doublereal* h__, integer* ldh, doublereal
	* wr, doublereal* wi, integer* iloz, integer* ihiz, doublereal* z__,
	integer* ldz, doublereal* work, integer* lwork, integer* info);


CLAPACK_API
int dlaqr1_(integer* n, doublereal* h__, integer* ldh,
	doublereal* sr1, doublereal* si1, doublereal* sr2, doublereal* si2,
	doublereal* v);


CLAPACK_API
int dlaqr2_(logical* wantt, logical* wantz, integer* n,
	integer* ktop, integer* kbot, integer* nw, doublereal* h__, integer*
	ldh, integer* iloz, integer* ihiz, doublereal* z__, integer* ldz,
	integer* ns, integer* nd, doublereal* sr, doublereal* si, doublereal*
	v, integer* ldv, integer* nh, doublereal* t, integer* ldt, integer*
	nv, doublereal* wv, integer* ldwv, doublereal* work, integer* lwork);


CLAPACK_API
int dlaqr3_(logical* wantt, logical* wantz, integer* n,
	integer* ktop, integer* kbot, integer* nw, doublereal* h__, integer*
	ldh, integer* iloz, integer* ihiz, doublereal* z__, integer* ldz,
	integer* ns, integer* nd, doublereal* sr, doublereal* si, doublereal*
	v, integer* ldv, integer* nh, doublereal* t, integer* ldt, integer*
	nv, doublereal* wv, integer* ldwv, doublereal* work, integer* lwork);

CLAPACK_API
int dlaqr4_(logical* wantt, logical* wantz, integer* n,
	integer* ilo, integer* ihi, doublereal* h__, integer* ldh, doublereal
	* wr, doublereal* wi, integer* iloz, integer* ihiz, doublereal* z__,
	integer* ldz, doublereal* work, integer* lwork, integer* info);


CLAPACK_API
int dlaqr5_(logical* wantt, logical* wantz, integer* kacc22,
	integer* n, integer* ktop, integer* kbot, integer* nshfts, doublereal
	* sr, doublereal* si, doublereal* h__, integer* ldh, integer* iloz,
	integer* ihiz, doublereal* z__, integer* ldz, doublereal* v, integer*
	ldv, doublereal* u, integer* ldu, integer* nv, doublereal* wv,
	integer* ldwv, integer* nh, doublereal* wh, integer* ldwh);


CLAPACK_API
int dlaqsb_(char* uplo, integer* n, integer* kd, doublereal*
	ab, integer* ldab, doublereal* s, doublereal* scond, doublereal* amax,
	char* equed);


CLAPACK_API
int dlaqsp_(char* uplo, integer* n, doublereal* ap,
	doublereal* s, doublereal* scond, doublereal* amax, char* equed);


CLAPACK_API
int dlaqsy_(char* uplo, integer* n, doublereal* a, integer*
	lda, doublereal* s, doublereal* scond, doublereal* amax, char* equed);


CLAPACK_API
int dlaqtr_(logical* ltran, logical* lreal, integer* n,
	doublereal* t, integer* ldt, doublereal* b, doublereal* w, doublereal
	* scale, doublereal* x, doublereal* work, integer* info);


CLAPACK_API
int dlar1v_(integer* n, integer* b1, integer* bn, doublereal
	* lambda, doublereal* d__, doublereal* l, doublereal* ld, doublereal*
	lld, doublereal* pivmin, doublereal* gaptol, doublereal* z__, logical
	* wantnc, integer* negcnt, doublereal* ztz, doublereal* mingma,
	integer* r__, integer* isuppz, doublereal* nrminv, doublereal* resid,
	doublereal* rqcorr, doublereal* work);


CLAPACK_API
int dlar2v_(integer* n, doublereal* x, doublereal* y,
	doublereal* z__, integer* incx, doublereal* c__, doublereal* s,
	integer* incc);


CLAPACK_API
int dlarf_(char* side, integer* m, integer* n, doublereal* v,
	integer* incv, doublereal* tau, doublereal* c__, integer* ldc,
	doublereal* work);


CLAPACK_API
int dlarfb_(char* side, char* trans, char* direct, char*
	storev, integer* m, integer* n, integer* k, doublereal* v, integer*
	ldv, doublereal* t, integer* ldt, doublereal* c__, integer* ldc,
	doublereal* work, integer* ldwork);


CLAPACK_API
int dlarfg_(integer* n, doublereal* alpha, doublereal* x,
	integer* incx, doublereal* tau);


CLAPACK_API
int dlarfp_(integer* n, doublereal* alpha, doublereal* x,
	integer* incx, doublereal* tau);


CLAPACK_API
int dlarft_(char* direct, char* storev, integer* n, integer*
	k, doublereal* v, integer* ldv, doublereal* tau, doublereal* t,
	integer* ldt);


CLAPACK_API
int dlarfx_(char* side, integer* m, integer* n, doublereal*
	v, doublereal* tau, doublereal* c__, integer* ldc, doublereal* work);


CLAPACK_API
int dlargv_(integer* n, doublereal* x, integer* incx,
	doublereal* y, integer* incy, doublereal* c__, integer* incc);


CLAPACK_API
int dlarrv_(integer* n, doublereal* vl, doublereal* vu,
	doublereal* d__, doublereal* l, doublereal* pivmin, integer* isplit,
	integer* m, integer* dol, integer* dou, doublereal* minrgp,
	doublereal* rtol1, doublereal* rtol2, doublereal* w, doublereal* werr,
	doublereal* wgap, integer* iblock, integer* indexw, doublereal* gers,
	doublereal* z__, integer* ldz, integer* isuppz, doublereal* work,
	integer* iwork, integer* info);


CLAPACK_API
int dlartv_(integer* n, doublereal* x, integer* incx,
	doublereal* y, integer* incy, doublereal* c__, doublereal* s, integer
	* incc);


CLAPACK_API
int dlarz_(char* side, integer* m, integer* n, integer* l,
	doublereal* v, integer* incv, doublereal* tau, doublereal* c__,
	integer* ldc, doublereal* work);


CLAPACK_API
int dlarzb_(char* side, char* trans, char* direct, char*
	storev, integer* m, integer* n, integer* k, integer* l, doublereal* v,
	integer* ldv, doublereal* t, integer* ldt, doublereal* c__, integer*
	ldc, doublereal* work, integer* ldwork);


CLAPACK_API
int dlarzt_(char* direct, char* storev, integer* n, integer*
	k, doublereal* v, integer* ldv, doublereal* tau, doublereal* t,
	integer* ldt);


CLAPACK_API
int dlaswp_(integer* n, doublereal* a, integer* lda, integer
	* k1, integer* k2, integer* ipiv, integer* incx);


CLAPACK_API
int dlasy2_(logical* ltranl, logical* ltranr, integer* isgn,
	integer* n1, integer* n2, doublereal* tl, integer* ldtl, doublereal*
	tr, integer* ldtr, doublereal* b, integer* ldb, doublereal* scale,
	doublereal* x, integer* ldx, doublereal* xnorm, integer* info);


CLAPACK_API
int dlasyf_(char* uplo, integer* n, integer* nb, integer* kb,
	doublereal* a, integer* lda, integer* ipiv, doublereal* w, integer*
	ldw, integer* info);


CLAPACK_API
int dlat2s_(char* uplo, integer* n, doublereal* a, integer*
	lda, real* sa, integer* ldsa, integer* info);


CLAPACK_API
int dlatbs_(char* uplo, char* trans, char* diag, char*
	normin, integer* n, integer* kd, doublereal* ab, integer* ldab,
	doublereal* x, doublereal* scale, doublereal* cnorm, integer* info);


CLAPACK_API
int dlatdf_(integer* ijob, integer* n, doublereal* z__,
	integer* ldz, doublereal* rhs, doublereal* rdsum, doublereal* rdscal,
	integer* ipiv, integer* jpiv);


CLAPACK_API
int dlatps_(char* uplo, char* trans, char* diag, char*
	normin, integer* n, doublereal* ap, doublereal* x, doublereal* scale,
	doublereal* cnorm, integer* info);


CLAPACK_API
int dlatrd_(char* uplo, integer* n, integer* nb, doublereal*
	a, integer* lda, doublereal* e, doublereal* tau, doublereal* w,
	integer* ldw);


CLAPACK_API
int dlatrs_(char* uplo, char* trans, char* diag, char*
	normin, integer* n, doublereal* a, integer* lda, doublereal* x,
	doublereal* scale, doublereal* cnorm, integer* info);


CLAPACK_API
int dlatrz_(integer* m, integer* n, integer* l, doublereal*
	a, integer* lda, doublereal* tau, doublereal* work);


CLAPACK_API
int dlatzm_(char* side, integer* m, integer* n, doublereal*
	v, integer* incv, doublereal* tau, doublereal* c1, doublereal* c2,
	integer* ldc, doublereal* work);


CLAPACK_API
int dlauu2_(char* uplo, integer* n, doublereal* a, integer*
	lda, integer* info);


CLAPACK_API
int dlauum_(char* uplo, integer* n, doublereal* a, integer*
	lda, integer* info);


CLAPACK_API
int dopgtr_(char* uplo, integer* n, doublereal* ap,
	doublereal* tau, doublereal* q, integer* ldq, doublereal* work,
	integer* info);


CLAPACK_API
int dopmtr_(char* side, char* uplo, char* trans, integer* m,
	integer* n, doublereal* ap, doublereal* tau, doublereal* c__, integer
	* ldc, doublereal* work, integer* info);


CLAPACK_API
int dorg2l_(integer* m, integer* n, integer* k, doublereal*
	a, integer* lda, doublereal* tau, doublereal* work, integer* info);


CLAPACK_API
int dorg2r_(integer* m, integer* n, integer* k, doublereal*
	a, integer* lda, doublereal* tau, doublereal* work, integer* info);


CLAPACK_API
int dorgbr_(char* vect, integer* m, integer* n, integer* k,
	doublereal* a, integer* lda, doublereal* tau, doublereal* work,
	integer* lwork, integer* info);


CLAPACK_API
int dorghr_(integer* n, integer* ilo, integer* ihi,
	doublereal* a, integer* lda, doublereal* tau, doublereal* work,
	integer* lwork, integer* info);


CLAPACK_API
int dorgl2_(integer* m, integer* n, integer* k, doublereal*
	a, integer* lda, doublereal* tau, doublereal* work, integer* info);


CLAPACK_API
int dorglq_(integer* m, integer* n, integer* k, doublereal*
	a, integer* lda, doublereal* tau, doublereal* work, integer* lwork,
	integer* info);


CLAPACK_API
int dorgql_(integer* m, integer* n, integer* k, doublereal*
	a, integer* lda, doublereal* tau, doublereal* work, integer* lwork,
	integer* info);


CLAPACK_API
int dorgqr_(integer* m, integer* n, integer* k, doublereal*
	a, integer* lda, doublereal* tau, doublereal* work, integer* lwork,
	integer* info);


CLAPACK_API
int dorgr2_(integer* m, integer* n, integer* k, doublereal*
	a, integer* lda, doublereal* tau, doublereal* work, integer* info);


CLAPACK_API
int dorgrq_(integer* m, integer* n, integer* k, doublereal*
	a, integer* lda, doublereal* tau, doublereal* work, integer* lwork,
	integer* info);


CLAPACK_API
int dorgtr_(char* uplo, integer* n, doublereal* a, integer*
	lda, doublereal* tau, doublereal* work, integer* lwork, integer* info);


CLAPACK_API
int dorm2l_(char* side, char* trans, integer* m, integer* n,
	integer* k, doublereal* a, integer* lda, doublereal* tau, doublereal*
	c__, integer* ldc, doublereal* work, integer* info);


CLAPACK_API
int dorm2r_(char* side, char* trans, integer* m, integer* n,
	integer* k, doublereal* a, integer* lda, doublereal* tau, doublereal*
	c__, integer* ldc, doublereal* work, integer* info);


CLAPACK_API
int dormbr_(char* vect, char* side, char* trans, integer* m,
	integer* n, integer* k, doublereal* a, integer* lda, doublereal* tau,
	doublereal* c__, integer* ldc, doublereal* work, integer* lwork,
	integer* info);


CLAPACK_API
int dormhr_(char* side, char* trans, integer* m, integer* n,
	integer* ilo, integer* ihi, doublereal* a, integer* lda, doublereal*
	tau, doublereal* c__, integer* ldc, doublereal* work, integer* lwork,
	integer* info);


CLAPACK_API
int dorml2_(char* side, char* trans, integer* m, integer* n,
	integer* k, doublereal* a, integer* lda, doublereal* tau, doublereal*
	c__, integer* ldc, doublereal* work, integer* info);


CLAPACK_API
int dormlq_(char* side, char* trans, integer* m, integer* n,
	integer* k, doublereal* a, integer* lda, doublereal* tau, doublereal*
	c__, integer* ldc, doublereal* work, integer* lwork, integer* info);


CLAPACK_API
int dormql_(char* side, char* trans, integer* m, integer* n,
	integer* k, doublereal* a, integer* lda, doublereal* tau, doublereal*
	c__, integer* ldc, doublereal* work, integer* lwork, integer* info);


CLAPACK_API
int dormqr_(char* side, char* trans, integer* m, integer* n,
	integer* k, doublereal* a, integer* lda, doublereal* tau, doublereal*
	c__, integer* ldc, doublereal* work, integer* lwork, integer* info);


CLAPACK_API
int dormr2_(char* side, char* trans, integer* m, integer* n,
	integer* k, doublereal* a, integer* lda, doublereal* tau, doublereal*
	c__, integer* ldc, doublereal* work, integer* info);


CLAPACK_API
int dormr3_(char* side, char* trans, integer* m, integer* n,
	integer* k, integer* l, doublereal* a, integer* lda, doublereal* tau,
	doublereal* c__, integer* ldc, doublereal* work, integer* info);


CLAPACK_API
int dormrq_(char* side, char* trans, integer* m, integer* n,
	integer* k, doublereal* a, integer* lda, doublereal* tau, doublereal*
	c__, integer* ldc, doublereal* work, integer* lwork, integer* info);


CLAPACK_API
int dormrz_(char* side, char* trans, integer* m, integer* n,
	integer* k, integer* l, doublereal* a, integer* lda, doublereal* tau,
	doublereal* c__, integer* ldc, doublereal* work, integer* lwork,
	integer* info);


CLAPACK_API
int dormtr_(char* side, char* uplo, char* trans, integer* m,
	integer* n, doublereal* a, integer* lda, doublereal* tau, doublereal*
	c__, integer* ldc, doublereal* work, integer* lwork, integer* info);


CLAPACK_API
int dpbcon_(char* uplo, integer* n, integer* kd, doublereal*
	ab, integer* ldab, doublereal* anorm, doublereal* rcond, doublereal*
	work, integer* iwork, integer* info);


CLAPACK_API
int dpbequ_(char* uplo, integer* n, integer* kd, doublereal*
	ab, integer* ldab, doublereal* s, doublereal* scond, doublereal* amax,
	integer* info);


CLAPACK_API
int dpbrfs_(char* uplo, integer* n, integer* kd, integer*
	nrhs, doublereal* ab, integer* ldab, doublereal* afb, integer* ldafb,
	doublereal* b, integer* ldb, doublereal* x, integer* ldx, doublereal*
	ferr, doublereal* berr, doublereal* work, integer* iwork, integer*
	info);


CLAPACK_API
int dpbstf_(char* uplo, integer* n, integer* kd, doublereal*
	ab, integer* ldab, integer* info);


CLAPACK_API
int dpbsv_(char* uplo, integer* n, integer* kd, integer*
	nrhs, doublereal* ab, integer* ldab, doublereal* b, integer* ldb,
	integer* info);


CLAPACK_API
int dpbsvx_(char* fact, char* uplo, integer* n, integer* kd,
	integer* nrhs, doublereal* ab, integer* ldab, doublereal* afb,
	integer* ldafb, char* equed, doublereal* s, doublereal* b, integer*
	ldb, doublereal* x, integer* ldx, doublereal* rcond, doublereal* ferr,
	doublereal* berr, doublereal* work, integer* iwork, integer* info);


CLAPACK_API
int dpbtf2_(char* uplo, integer* n, integer* kd, doublereal*
	ab, integer* ldab, integer* info);


CLAPACK_API
int dpbtrf_(char* uplo, integer* n, integer* kd, doublereal*
	ab, integer* ldab, integer* info);


CLAPACK_API
int dpbtrs_(char* uplo, integer* n, integer* kd, integer*
	nrhs, doublereal* ab, integer* ldab, doublereal* b, integer* ldb,
	integer* info);


CLAPACK_API
int dpftrf_(char* transr, char* uplo, integer* n, doublereal
	* a, integer* info);


CLAPACK_API
int dpftri_(char* transr, char* uplo, integer* n, doublereal
	* a, integer* info);


CLAPACK_API
int dpftrs_(char* transr, char* uplo, integer* n, integer*
	nrhs, doublereal* a, doublereal* b, integer* ldb, integer* info);


CLAPACK_API
int dpocon_(char* uplo, integer* n, doublereal* a, integer*
	lda, doublereal* anorm, doublereal* rcond, doublereal* work, integer*
	iwork, integer* info);


CLAPACK_API
int dpoequ_(integer* n, doublereal* a, integer* lda,
	doublereal* s, doublereal* scond, doublereal* amax, integer* info);


CLAPACK_API
int dpoequb_(integer* n, doublereal* a, integer* lda,
	doublereal* s, doublereal* scond, doublereal* amax, integer* info);


CLAPACK_API
int dporfs_(char* uplo, integer* n, integer* nrhs,
	doublereal* a, integer* lda, doublereal* af, integer* ldaf,
	doublereal* b, integer* ldb, doublereal* x, integer* ldx, doublereal*
	ferr, doublereal* berr, doublereal* work, integer* iwork, integer*
	info);


CLAPACK_API
int dposv_(char* uplo, integer* n, integer* nrhs, doublereal
	* a, integer* lda, doublereal* b, integer* ldb, integer* info);


CLAPACK_API
int dposvx_(char* fact, char* uplo, integer* n, integer*
	nrhs, doublereal* a, integer* lda, doublereal* af, integer* ldaf,
	char* equed, doublereal* s, doublereal* b, integer* ldb, doublereal*
	x, integer* ldx, doublereal* rcond, doublereal* ferr, doublereal*
	berr, doublereal* work, integer* iwork, integer* info);


CLAPACK_API
int dpotf2_(char* uplo, integer* n, doublereal* a, integer*
	lda, integer* info);


CLAPACK_API
int dpotrf_(char* uplo, integer* n, doublereal* a, integer*
	lda, integer* info);


CLAPACK_API
int dpotri_(char* uplo, integer* n, doublereal* a, integer*
	lda, integer* info);


CLAPACK_API
int dpotrs_(char* uplo, integer* n, integer* nrhs,
	doublereal* a, integer* lda, doublereal* b, integer* ldb, integer*
	info);


CLAPACK_API
int dppcon_(char* uplo, integer* n, doublereal* ap,
	doublereal* anorm, doublereal* rcond, doublereal* work, integer*
	iwork, integer* info);


CLAPACK_API
int dppequ_(char* uplo, integer* n, doublereal* ap,
	doublereal* s, doublereal* scond, doublereal* amax, integer* info);


CLAPACK_API
int dpprfs_(char* uplo, integer* n, integer* nrhs,
	doublereal* ap, doublereal* afp, doublereal* b, integer* ldb,
	doublereal* x, integer* ldx, doublereal* ferr, doublereal* berr,
	doublereal* work, integer* iwork, integer* info);


CLAPACK_API
int dppsv_(char* uplo, integer* n, integer* nrhs, doublereal
	* ap, doublereal* b, integer* ldb, integer* info);


CLAPACK_API
int dppsvx_(char* fact, char* uplo, integer* n, integer*
	nrhs, doublereal* ap, doublereal* afp, char* equed, doublereal* s,
	doublereal* b, integer* ldb, doublereal* x, integer* ldx, doublereal*
	rcond, doublereal* ferr, doublereal* berr, doublereal* work, integer*
	iwork, integer* info);


CLAPACK_API
int dpptrf_(char* uplo, integer* n, doublereal* ap, integer*
	info);


CLAPACK_API
int dpptri_(char* uplo, integer* n, doublereal* ap, integer*
	info);


CLAPACK_API
int dpptrs_(char* uplo, integer* n, integer* nrhs,
	doublereal* ap, doublereal* b, integer* ldb, integer* info);


CLAPACK_API
int dpstf2_(char* uplo, integer* n, doublereal* a, integer*
	lda, integer* piv, integer* rank, doublereal* tol, doublereal* work,
	integer* info);


CLAPACK_API
int dpstrf_(char* uplo, integer* n, doublereal* a, integer*
	lda, integer* piv, integer* rank, doublereal* tol, doublereal* work,
	integer* info);


CLAPACK_API
int dptcon_(integer* n, doublereal* d__, doublereal* e,
	doublereal* anorm, doublereal* rcond, doublereal* work, integer* info);


CLAPACK_API
int dpteqr_(char* compz, integer* n, doublereal* d__,
	doublereal* e, doublereal* z__, integer* ldz, doublereal* work,
	integer* info);


CLAPACK_API
int dptrfs_(integer* n, integer* nrhs, doublereal* d__,
	doublereal* e, doublereal* df, doublereal* ef, doublereal* b, integer
	* ldb, doublereal* x, integer* ldx, doublereal* ferr, doublereal* berr,
	doublereal* work, integer* info);


CLAPACK_API
int dptsv_(integer* n, integer* nrhs, doublereal* d__,
	doublereal* e, doublereal* b, integer* ldb, integer* info);


CLAPACK_API
int dptsvx_(char* fact, integer* n, integer* nrhs,
	doublereal* d__, doublereal* e, doublereal* df, doublereal* ef,
	doublereal* b, integer* ldb, doublereal* x, integer* ldx, doublereal*
	rcond, doublereal* ferr, doublereal* berr, doublereal* work, integer*
	info);


CLAPACK_API
int dpttrs_(integer* n, integer* nrhs, doublereal* d__,
	doublereal* e, doublereal* b, integer* ldb, integer* info);


CLAPACK_API
int dptts2_(integer* n, integer* nrhs, doublereal* d__,
	doublereal* e, doublereal* b, integer* ldb);


CLAPACK_API
int drscl_(integer* n, doublereal* sa, doublereal* sx,
	integer* incx);


CLAPACK_API
int dsbev_(char* jobz, char* uplo, integer* n, integer* kd,
	doublereal* ab, integer* ldab, doublereal* w, doublereal* z__,
	integer* ldz, doublereal* work, integer* info);


CLAPACK_API
int dsbevd_(char* jobz, char* uplo, integer* n, integer* kd,
	doublereal* ab, integer* ldab, doublereal* w, doublereal* z__,
	integer* ldz, doublereal* work, integer* lwork, integer* iwork,
	integer* liwork, integer* info);


CLAPACK_API
int dsbevx_(char* jobz, char* range, char* uplo, integer* n,
	integer* kd, doublereal* ab, integer* ldab, doublereal* q, integer*
	ldq, doublereal* vl, doublereal* vu, integer* il, integer* iu,
	doublereal* abstol, integer* m, doublereal* w, doublereal* z__,
	integer* ldz, doublereal* work, integer* iwork, integer* ifail,
	integer* info);


CLAPACK_API
int dsbgst_(char* vect, char* uplo, integer* n, integer* ka,
	integer* kb, doublereal* ab, integer* ldab, doublereal* bb, integer*
	ldbb, doublereal* x, integer* ldx, doublereal* work, integer* info);


CLAPACK_API
int dsbgv_(char* jobz, char* uplo, integer* n, integer* ka,
	integer* kb, doublereal* ab, integer* ldab, doublereal* bb, integer*
	ldbb, doublereal* w, doublereal* z__, integer* ldz, doublereal* work,
	integer* info);


CLAPACK_API
int dsbgvd_(char* jobz, char* uplo, integer* n, integer* ka,
	integer* kb, doublereal* ab, integer* ldab, doublereal* bb, integer*
	ldbb, doublereal* w, doublereal* z__, integer* ldz, doublereal* work,
	integer* lwork, integer* iwork, integer* liwork, integer* info);


CLAPACK_API
int dsbgvx_(char* jobz, char* range, char* uplo, integer* n,
	integer* ka, integer* kb, doublereal* ab, integer* ldab, doublereal*
	bb, integer* ldbb, doublereal* q, integer* ldq, doublereal* vl,
	doublereal* vu, integer* il, integer* iu, doublereal* abstol, integer
	* m, doublereal* w, doublereal* z__, integer* ldz, doublereal* work,
	integer* iwork, integer* ifail, integer* info);


CLAPACK_API
int dsbtrd_(char* vect, char* uplo, integer* n, integer* kd,
	doublereal* ab, integer* ldab, doublereal* d__, doublereal* e,
	doublereal* q, integer* ldq, doublereal* work, integer* info);


CLAPACK_API
int dsfrk_(char* transr, char* uplo, char* trans, integer* n,
	integer* k, doublereal* alpha, doublereal* a, integer* lda,
	doublereal* beta, doublereal* c__);


CLAPACK_API int dsgesv_(integer* n, integer* nrhs, doublereal* a,
	integer* lda, integer* ipiv, doublereal* b, integer* ldb, doublereal*
	x, integer* ldx, doublereal* work, real* swork, integer* iter,
	integer* info);


CLAPACK_API
int dspcon_(char* uplo, integer* n, doublereal* ap, integer*
	ipiv, doublereal* anorm, doublereal* rcond, doublereal* work, integer
	* iwork, integer* info);


CLAPACK_API
int dspev_(char* jobz, char* uplo, integer* n, doublereal*
	ap, doublereal* w, doublereal* z__, integer* ldz, doublereal* work,
	integer* info);


CLAPACK_API int dspevd_(char* jobz, char* uplo, integer* n, doublereal*
	ap, doublereal* w, doublereal* z__, integer* ldz, doublereal* work,
	integer* lwork, integer* iwork, integer* liwork, integer* info);


CLAPACK_API
int dspevx_(char* jobz, char* range, char* uplo, integer* n,
	doublereal* ap, doublereal* vl, doublereal* vu, integer* il, integer*
	iu, doublereal* abstol, integer* m, doublereal* w, doublereal* z__,
	integer* ldz, doublereal* work, integer* iwork, integer* ifail,
	integer* info);


CLAPACK_API
int dspgst_(integer* itype, char* uplo, integer* n,
	doublereal* ap, doublereal* bp, integer* info);


CLAPACK_API
int dspgv_(integer* itype, char* jobz, char* uplo, integer*
	n, doublereal* ap, doublereal* bp, doublereal* w, doublereal* z__,
	integer* ldz, doublereal* work, integer* info);


CLAPACK_API
int dspgvd_(integer* itype, char* jobz, char* uplo, integer*
	n, doublereal* ap, doublereal* bp, doublereal* w, doublereal* z__,
	integer* ldz, doublereal* work, integer* lwork, integer* iwork,
	integer* liwork, integer* info);


CLAPACK_API
int dspgvx_(integer* itype, char* jobz, char* range, char*
	uplo, integer* n, doublereal* ap, doublereal* bp, doublereal* vl,
	doublereal* vu, integer* il, integer* iu, doublereal* abstol, integer
	* m, doublereal* w, doublereal* z__, integer* ldz, doublereal* work,
	integer* iwork, integer* ifail, integer* info);


CLAPACK_API
int dsposv_(char* uplo, integer* n, integer* nrhs,
	doublereal* a, integer* lda, doublereal* b, integer* ldb, doublereal*
	x, integer* ldx, doublereal* work, real* swork, integer* iter,
	integer* info);


CLAPACK_API
int dsprfs_(char* uplo, integer* n, integer* nrhs,
	doublereal* ap, doublereal* afp, integer* ipiv, doublereal* b,
	integer* ldb, doublereal* x, integer* ldx, doublereal* ferr,
	doublereal* berr, doublereal* work, integer* iwork, integer* info);


CLAPACK_API
int dspsv_(char* uplo, integer* n, integer* nrhs, doublereal
	* ap, integer* ipiv, doublereal* b, integer* ldb, integer* info);


CLAPACK_API
int dspsvx_(char* fact, char* uplo, integer* n, integer*
	nrhs, doublereal* ap, doublereal* afp, integer* ipiv, doublereal* b,
	integer* ldb, doublereal* x, integer* ldx, doublereal* rcond,
	doublereal* ferr, doublereal* berr, doublereal* work, integer* iwork,
	integer* info);


CLAPACK_API
int dsptrd_(char* uplo, integer* n, doublereal* ap,
	doublereal* d__, doublereal* e, doublereal* tau, integer* info);


CLAPACK_API
int dsptrf_(char* uplo, integer* n, doublereal* ap, integer*
	ipiv, integer* info);


CLAPACK_API
int dsptri_(char* uplo, integer* n, doublereal* ap, integer*
	ipiv, doublereal* work, integer* info);


CLAPACK_API
int dsptrs_(char* uplo, integer* n, integer* nrhs,
	doublereal* ap, integer* ipiv, doublereal* b, integer* ldb, integer*
	info);


CLAPACK_API
int dstegr_(char* jobz, char* range, integer* n, doublereal*
	d__, doublereal* e, doublereal* vl, doublereal* vu, integer* il,
	integer* iu, doublereal* abstol, integer* m, doublereal* w,
	doublereal* z__, integer* ldz, integer* isuppz, doublereal* work,
	integer* lwork, integer* iwork, integer* liwork, integer* info);


CLAPACK_API
int dstein_(integer* n, doublereal* d__, doublereal* e,
	integer* m, doublereal* w, integer* iblock, integer* isplit,
	doublereal* z__, integer* ldz, doublereal* work, integer* iwork,
	integer* ifail, integer* info);


CLAPACK_API
int dstemr_(char* jobz, char* range, integer* n, doublereal*
	d__, doublereal* e, doublereal* vl, doublereal* vu, integer* il,
	integer* iu, integer* m, doublereal* w, doublereal* z__, integer* ldz,
	integer* nzc, integer* isuppz, logical* tryrac, doublereal* work,
	integer* lwork, integer* iwork, integer* liwork, integer* info);


CLAPACK_API
int dstev_(char* jobz, integer* n, doublereal* d__,
	doublereal* e, doublereal* z__, integer* ldz, doublereal* work,
	integer* info);


CLAPACK_API
int dstevd_(char* jobz, integer* n, doublereal* d__,
	doublereal* e, doublereal* z__, integer* ldz, doublereal* work,
	integer* lwork, integer* iwork, integer* liwork, integer* info);


CLAPACK_API
int dstevr_(char* jobz, char* range, integer* n, doublereal*
	d__, doublereal* e, doublereal* vl, doublereal* vu, integer* il,
	integer* iu, doublereal* abstol, integer* m, doublereal* w,
	doublereal* z__, integer* ldz, integer* isuppz, doublereal* work,
	integer* lwork, integer* iwork, integer* liwork, integer* info);


CLAPACK_API
int dstevx_(char* jobz, char* range, integer* n, doublereal*
	d__, doublereal* e, doublereal* vl, doublereal* vu, integer* il,
	integer* iu, doublereal* abstol, integer* m, doublereal* w,
	doublereal* z__, integer* ldz, doublereal* work, integer* iwork,
	integer* ifail, integer* info);


CLAPACK_API
int dsycon_(char* uplo, integer* n, doublereal* a, integer*
	lda, integer* ipiv, doublereal* anorm, doublereal* rcond, doublereal*
	work, integer* iwork, integer* info);


CLAPACK_API
int dsyequb_(char* uplo, integer* n, doublereal* a, integer*
	lda, doublereal* s, doublereal* scond, doublereal* amax, doublereal*
	work, integer* info);


CLAPACK_API
int dsyev_(char* jobz, char* uplo, integer* n, doublereal* a,
	integer* lda, doublereal* w, doublereal* work, integer* lwork,
	integer* info);


CLAPACK_API
int dsyevd_(char* jobz, char* uplo, integer* n, doublereal*
	a, integer* lda, doublereal* w, doublereal* work, integer* lwork,
	integer* iwork, integer* liwork, integer* info);


CLAPACK_API
int dsyevr_(char* jobz, char* range, char* uplo, integer* n,
	doublereal* a, integer* lda, doublereal* vl, doublereal* vu, integer*
	il, integer* iu, doublereal* abstol, integer* m, doublereal* w,
	doublereal* z__, integer* ldz, integer* isuppz, doublereal* work,
	integer* lwork, integer* iwork, integer* liwork, integer* info);


CLAPACK_API
int dsyevx_(char* jobz, char* range, char* uplo, integer* n,
	doublereal* a, integer* lda, doublereal* vl, doublereal* vu, integer*
	il, integer* iu, doublereal* abstol, integer* m, doublereal* w,
	doublereal* z__, integer* ldz, doublereal* work, integer* lwork,
	integer* iwork, integer* ifail, integer* info);


CLAPACK_API
int dsygs2_(integer* itype, char* uplo, integer* n,
	doublereal* a, integer* lda, doublereal* b, integer* ldb, integer*
	info);


CLAPACK_API
int dsygst_(integer* itype, char* uplo, integer* n,
	doublereal* a, integer* lda, doublereal* b, integer* ldb, integer*
	info);


CLAPACK_API
int dsygv_(integer* itype, char* jobz, char* uplo, integer*
	n, doublereal* a, integer* lda, doublereal* b, integer* ldb,
	doublereal* w, doublereal* work, integer* lwork, integer* info);


CLAPACK_API
int dsygvd_(integer* itype, char* jobz, char* uplo, integer*
	n, doublereal* a, integer* lda, doublereal* b, integer* ldb,
	doublereal* w, doublereal* work, integer* lwork, integer* iwork,
	integer* liwork, integer* info);


CLAPACK_API
int dsygvx_(integer* itype, char* jobz, char* range, char*
	uplo, integer* n, doublereal* a, integer* lda, doublereal* b, integer
	* ldb, doublereal* vl, doublereal* vu, integer* il, integer* iu,
	doublereal* abstol, integer* m, doublereal* w, doublereal* z__,
	integer* ldz, doublereal* work, integer* lwork, integer* iwork,
	integer* ifail, integer* info);


CLAPACK_API
int dsyrfs_(char* uplo, integer* n, integer* nrhs,
	doublereal* a, integer* lda, doublereal* af, integer* ldaf, integer*
	ipiv, doublereal* b, integer* ldb, doublereal* x, integer* ldx,
	doublereal* ferr, doublereal* berr, doublereal* work, integer* iwork,
	integer* info);


CLAPACK_API
int dsysv_(char* uplo, integer* n, integer* nrhs, doublereal
	* a, integer* lda, integer* ipiv, doublereal* b, integer* ldb,
	doublereal* work, integer* lwork, integer* info);


CLAPACK_API
int dsysvx_(char* fact, char* uplo, integer* n, integer*
	nrhs, doublereal* a, integer* lda, doublereal* af, integer* ldaf,
	integer* ipiv, doublereal* b, integer* ldb, doublereal* x, integer*
	ldx, doublereal* rcond, doublereal* ferr, doublereal* berr,
	doublereal* work, integer* lwork, integer* iwork, integer* info);


CLAPACK_API
int dsytd2_(char* uplo, integer* n, doublereal* a, integer*
	lda, doublereal* d__, doublereal* e, doublereal* tau, integer* info);


CLAPACK_API
int dsytf2_(char* uplo, integer* n, doublereal* a, integer*
	lda, integer* ipiv, integer* info);


CLAPACK_API
int dsytrd_(char* uplo, integer* n, doublereal* a, integer*
	lda, doublereal* d__, doublereal* e, doublereal* tau, doublereal*
	work, integer* lwork, integer* info);


CLAPACK_API
int dsytrf_(char* uplo, integer* n, doublereal* a, integer*
	lda, integer* ipiv, doublereal* work, integer* lwork, integer* info);


CLAPACK_API
int dsytri_(char* uplo, integer* n, doublereal* a, integer*
	lda, integer* ipiv, doublereal* work, integer* info);


CLAPACK_API
int dsytrs_(char* uplo, integer* n, integer* nrhs,
	doublereal* a, integer* lda, integer* ipiv, doublereal* b, integer*
	ldb, integer* info);


CLAPACK_API
int dtbcon_(char* norm, char* uplo, char* diag, integer* n,
	integer* kd, doublereal* ab, integer* ldab, doublereal* rcond,
	doublereal* work, integer* iwork, integer* info);


CLAPACK_API
int dtbrfs_(char* uplo, char* trans, char* diag, integer* n,
	integer* kd, integer* nrhs, doublereal* ab, integer* ldab, doublereal
	* b, integer* ldb, doublereal* x, integer* ldx, doublereal* ferr,
	doublereal* berr, doublereal* work, integer* iwork, integer* info);


CLAPACK_API
int dtbtrs_(char* uplo, char* trans, char* diag, integer* n,
	integer* kd, integer* nrhs, doublereal* ab, integer* ldab, doublereal
	* b, integer* ldb, integer* info);


CLAPACK_API
int dtfsm_(char* transr, char* side, char* uplo, char* trans,
	char* diag, integer* m, integer* n, doublereal* alpha, doublereal* a,
	doublereal* b, integer* ldb);


CLAPACK_API
int dtftri_(char* transr, char* uplo, char* diag, integer* n,
	doublereal* a, integer* info);


CLAPACK_API int dtfttp_(char* transr, char* uplo, integer* n, doublereal
	* arf, doublereal* ap, integer* info);


CLAPACK_API
int dtfttr_(char* transr, char* uplo, integer* n, doublereal
	* arf, doublereal* a, integer* lda, integer* info);


CLAPACK_API
int dtgevc_(char* side, char* howmny, logical* select,
	integer* n, doublereal* s, integer* lds, doublereal* p, integer* ldp,
	doublereal* vl, integer* ldvl, doublereal* vr, integer* ldvr, integer
	* mm, integer* m, doublereal* work, integer* info);


CLAPACK_API int dtgex2_(logical* wantq, logical* wantz, integer* n,
	doublereal* a, integer* lda, doublereal* b, integer* ldb, doublereal*
	q, integer* ldq, doublereal* z__, integer* ldz, integer* j1, integer*
	n1, integer* n2, doublereal* work, integer* lwork, integer* info);


CLAPACK_API
int dtgexc_(logical* wantq, logical* wantz, integer* n,
	doublereal* a, integer* lda, doublereal* b, integer* ldb, doublereal*
	q, integer* ldq, doublereal* z__, integer* ldz, integer* ifst,
	integer* ilst, doublereal* work, integer* lwork, integer* info);


CLAPACK_API
int dtgsen_(integer* ijob, logical* wantq, logical* wantz,
	logical* select, integer* n, doublereal* a, integer* lda, doublereal*
	b, integer* ldb, doublereal* alphar, doublereal* alphai, doublereal*
	beta, doublereal* q, integer* ldq, doublereal* z__, integer* ldz,
	integer* m, doublereal* pl, doublereal* pr, doublereal* dif,
	doublereal* work, integer* lwork, integer* iwork, integer* liwork,
	integer* info);


CLAPACK_API
int dtgsja_(char* jobu, char* jobv, char* jobq, integer* m,
	integer* p, integer* n, integer* k, integer* l, doublereal* a,
	integer* lda, doublereal* b, integer* ldb, doublereal* tola,
	doublereal* tolb, doublereal* alpha, doublereal* beta, doublereal* u,
	integer* ldu, doublereal* v, integer* ldv, doublereal* q, integer*
	ldq, doublereal* work, integer* ncycle, integer* info);


CLAPACK_API
int dtgsna_(char* job, char* howmny, logical* select,
	integer* n, doublereal* a, integer* lda, doublereal* b, integer* ldb,
	doublereal* vl, integer* ldvl, doublereal* vr, integer* ldvr,
	doublereal* s, doublereal* dif, integer* mm, integer* m, doublereal*
	work, integer* lwork, integer* iwork, integer* info);


CLAPACK_API
int dtgsy2_(char* trans, integer* ijob, integer* m, integer*
	n, doublereal* a, integer* lda, doublereal* b, integer* ldb,
	doublereal* c__, integer* ldc, doublereal* d__, integer* ldd,
	doublereal* e, integer* lde, doublereal* f, integer* ldf, doublereal*
	scale, doublereal* rdsum, doublereal* rdscal, integer* iwork, integer
	* pq, integer* info);


CLAPACK_API
int dtgsyl_(char* trans, integer* ijob, integer* m, integer*
	n, doublereal* a, integer* lda, doublereal* b, integer* ldb,
	doublereal* c__, integer* ldc, doublereal* d__, integer* ldd,
	doublereal* e, integer* lde, doublereal* f, integer* ldf, doublereal*
	scale, doublereal* dif, doublereal* work, integer* lwork, integer*
	iwork, integer* info);


CLAPACK_API
int dtpcon_(char* norm, char* uplo, char* diag, integer* n,
	doublereal* ap, doublereal* rcond, doublereal* work, integer* iwork,
	integer* info);


CLAPACK_API
int dtprfs_(char* uplo, char* trans, char* diag, integer* n,
	integer* nrhs, doublereal* ap, doublereal* b, integer* ldb,
	doublereal* x, integer* ldx, doublereal* ferr, doublereal* berr,
	doublereal* work, integer* iwork, integer* info);


CLAPACK_API
int dtptri_(char* uplo, char* diag, integer* n, doublereal*
	ap, integer* info);


CLAPACK_API int dtptrs_(char* uplo, char* trans, char* diag, integer* n,
	integer* nrhs, doublereal* ap, doublereal* b, integer* ldb, integer*
	info);


CLAPACK_API
int dtpttf_(char* transr, char* uplo, integer* n, doublereal
	* ap, doublereal* arf, integer* info);


CLAPACK_API
int dtpttr_(char* uplo, integer* n, doublereal* ap,
	doublereal* a, integer* lda, integer* info);


CLAPACK_API
int dtrcon_(char* norm, char* uplo, char* diag, integer* n,
	doublereal* a, integer* lda, doublereal* rcond, doublereal* work,
	integer* iwork, integer* info);


CLAPACK_API
int dtrevc_(char* side, char* howmny, logical* select,
	integer* n, doublereal* t, integer* ldt, doublereal* vl, integer*
	ldvl, doublereal* vr, integer* ldvr, integer* mm, integer* m,
	doublereal* work, integer* info);


CLAPACK_API
int dtrexc_(char* compq, integer* n, doublereal* t, integer*
	ldt, doublereal* q, integer* ldq, integer* ifst, integer* ilst,
	doublereal* work, integer* info);


CLAPACK_API
int dtrrfs_(char* uplo, char* trans, char* diag, integer* n,
	integer* nrhs, doublereal* a, integer* lda, doublereal* b, integer*
	ldb, doublereal* x, integer* ldx, doublereal* ferr, doublereal* berr,
	doublereal* work, integer* iwork, integer* info);


CLAPACK_API
int dtrsen_(char* job, char* compq, logical* select, integer
	* n, doublereal* t, integer* ldt, doublereal* q, integer* ldq,
	doublereal* wr, doublereal* wi, integer* m, doublereal* s, doublereal
	* sep, doublereal* work, integer* lwork, integer* iwork, integer*
	liwork, integer* info);


CLAPACK_API
int dtrsna_(char* job, char* howmny, logical* select,
	integer* n, doublereal* t, integer* ldt, doublereal* vl, integer*
	ldvl, doublereal* vr, integer* ldvr, doublereal* s, doublereal* sep,
	integer* mm, integer* m, doublereal* work, integer* ldwork, integer*
	iwork, integer* info);


CLAPACK_API
int dtrsyl_(char* trana, char* tranb, integer* isgn, integer
	* m, integer* n, doublereal* a, integer* lda, doublereal* b, integer*
	ldb, doublereal* c__, integer* ldc, doublereal* scale, integer* info);


CLAPACK_API
int dtrti2_(char* uplo, char* diag, integer* n, doublereal*
	a, integer* lda, integer* info);


CLAPACK_API
int dtrtri_(char* uplo, char* diag, integer* n, doublereal*
	a, integer* lda, integer* info);


CLAPACK_API
int dtrtrs_(char* uplo, char* trans, char* diag, integer* n,
	integer* nrhs, doublereal* a, integer* lda, doublereal* b, integer*
	ldb, integer* info);


CLAPACK_API
int dtrttf_(char* transr, char* uplo, integer* n, doublereal
	* a, integer* lda, doublereal* arf, integer* info);


CLAPACK_API
int dtrttp_(char* uplo, integer* n, doublereal* a, integer*
	lda, doublereal* ap, integer* info);


CLAPACK_API
int dtzrqf_(integer* m, integer* n, doublereal* a, integer*
	lda, doublereal* tau, integer* info);


CLAPACK_API
int dtzrzf_(integer* m, integer* n, doublereal* a, integer*
	lda, doublereal* tau, doublereal* work, integer* lwork, integer* info);

CLAPACK_API
integer iladlc_(integer* m, integer* n, doublereal* a, integer* lda);

CLAPACK_API
integer iladlr_(integer* m, integer* n, doublereal* a, integer* lda);

//}}@@@ Finished

#pragma endregion

#pragma region DXLASRC -- Double precision real LAPACK routines using extra precision (dlax) - ok
//#       
//set(DXLASRC
//	
//  dgerfsx.c dgesvxx.c
//	dla_gerfsx_extended.c dla_geamv.c dla_gercond.c dla_rpvgrw.c dsysvxx.c dsyrfsx.c
//	dla_syrfsx_extended.c dla_syamv.c dla_syrcond.c dla_syrpvgrw.c
//	dposvxx.c dporfsx.c dla_porfsx_extended.c dla_porcond.c
//	dla_porpvgrw.c dgbsvxx.c dgbrfsx.c dla_gbrfsx_extended.c
//	dla_gbamv.c dla_gbrcond.c dla_gbrpvgrw.c dla_lin_berr.c dlarscl2.c
//	dlascl2.c dla_wwaddw.c)
//

/*
int dgbrfsx_(char* trans, char* equed, integer* n, integer*
	kl, integer* ku, integer* nrhs, doublereal* ab, integer* ldab,
	doublereal* afb, integer* ldafb, integer* ipiv, doublereal* r__,
	doublereal* c__, doublereal* b, integer* ldb, doublereal* x, integer*
	ldx, doublereal* rcond, doublereal* berr, integer* n_err_bnds__,
	doublereal* err_bnds_norm__, doublereal* err_bnds_comp__, integer*
	nparams, doublereal* params, doublereal* work, integer* iwork,
	integer* info);

int dgbsvxx_(char* fact, char* trans, integer* n, integer*
	kl, integer* ku, integer* nrhs, doublereal* ab, integer* ldab,
	doublereal* afb, integer* ldafb, integer* ipiv, char* equed,
	doublereal* r__, doublereal* c__, doublereal* b, integer* ldb,
	doublereal* x, integer* ldx, doublereal* rcond, doublereal* rpvgrw,
	doublereal* berr, integer* n_err_bnds__, doublereal* err_bnds_norm__,
	doublereal* err_bnds_comp__, integer* nparams, doublereal* params,
	doublereal* work, integer* iwork, integer* info);

int dgerfsx_(char* trans, char* equed, integer* n, integer*
	nrhs, doublereal* a, integer* lda, doublereal* af, integer* ldaf,
	integer* ipiv, doublereal* r__, doublereal* c__, doublereal* b,
	integer* ldb, doublereal* x, integer* ldx, doublereal* rcond,
	doublereal* berr, integer* n_err_bnds__, doublereal* err_bnds_norm__,
	doublereal* err_bnds_comp__, integer* nparams, doublereal* params,
	doublereal* work, integer* iwork, integer* info);

int dgesvxx_(char* fact, char* trans, integer* n, integer*
	nrhs, doublereal* a, integer* lda, doublereal* af, integer* ldaf,
	integer* ipiv, char* equed, doublereal* r__, doublereal* c__,
	doublereal* b, integer* ldb, doublereal* x, integer* ldx, doublereal*
	rcond, doublereal* rpvgrw, doublereal* berr, integer* n_err_bnds__,
	doublereal* err_bnds_norm__, doublereal* err_bnds_comp__, integer*
	nparams, doublereal* params, doublereal* work, integer* iwork,
	integer* info);

int dla_gbamv__(integer* trans, integer* m, integer* n,
	integer* kl, integer* ku, doublereal* alpha, doublereal* ab, integer*
	ldab, doublereal* x, integer* incx, doublereal* beta, doublereal* y,
	integer* incy);

doublereal dla_gbrcond__(char* trans, integer* n, integer* kl, integer* ku,
	doublereal* ab, integer* ldab, doublereal* afb, integer* ldafb,
	integer* ipiv, integer* cmode, doublereal* c__, integer* info,
	doublereal* work, integer* iwork, ftnlen trans_len);

int dla_gbrfsx_extended__(integer* prec_type__, integer*
	trans_type__, integer* n, integer* kl, integer* ku, integer* nrhs,
	doublereal* ab, integer* ldab, doublereal* afb, integer* ldafb,
	integer* ipiv, logical* colequ, doublereal* c__, doublereal* b,
	integer* ldb, doublereal* y, integer* ldy, doublereal* berr_out__,
	integer* n_norms__, doublereal* errs_n__, doublereal* errs_c__,
	doublereal* res, doublereal* ayb, doublereal* dy, doublereal*
	y_tail__, doublereal* rcond, integer* ithresh, doublereal* rthresh,
	doublereal* dz_ub__, logical* ignore_cwise__, integer* info);

doublereal dla_gbrpvgrw__(integer* n, integer* kl, integer* ku, integer*
	ncols, doublereal* ab, integer* ldab, doublereal* afb, integer* ldafb);


int dla_geamv__(integer* trans, integer* m, integer* n,
	doublereal* alpha, doublereal* a, integer* lda, doublereal* x,
	integer* incx, doublereal* beta, doublereal* y, integer* incy);

CLAPACK_API
doublereal dla_gercond__(char* trans, integer* n, doublereal* a, integer* lda,
	doublereal* af, integer* ldaf, integer* ipiv, integer* cmode,
	doublereal* c__, integer* info, doublereal* work, integer* iwork,
	ftnlen trans_len);


int dla_gerfsx_extended__(integer* prec_type__, integer*
	trans_type__, integer* n, integer* nrhs, doublereal* a, integer* lda,
	doublereal* af, integer* ldaf, integer* ipiv, logical* colequ,
	doublereal* c__, doublereal* b, integer* ldb, doublereal* y, integer*
	ldy, doublereal* berr_out__, integer* n_norms__, doublereal* errs_n__,
	doublereal* errs_c__, doublereal* res, doublereal* ayb, doublereal*
	dy, doublereal* y_tail__, doublereal* rcond, integer* ithresh,
	doublereal* rthresh, doublereal* dz_ub__, logical* ignore_cwise__,
	integer* info);


int dla_lin_berr__(integer* n, integer* nz, integer* nrhs,
	doublereal* res, doublereal* ayb, doublereal* berr);

doublereal dla_porcond__(char* uplo, integer* n, doublereal* a, integer* lda,
	doublereal* af, integer* ldaf, integer* cmode, doublereal* c__,
	integer* info, doublereal* work, integer* iwork, ftnlen uplo_len);

int dla_porfsx_extended__(integer* prec_type__, char* uplo,
	integer* n, integer* nrhs, doublereal* a, integer* lda, doublereal*
	af, integer* ldaf, logical* colequ, doublereal* c__, doublereal* b,
	integer* ldb, doublereal* y, integer* ldy, doublereal* berr_out__,
	integer* n_norms__, doublereal* errs_n__, doublereal* errs_c__,
	doublereal* res, doublereal* ayb, doublereal* dy, doublereal*
	y_tail__, doublereal* rcond, integer* ithresh, doublereal* rthresh,
	doublereal* dz_ub__, logical* ignore_cwise__, integer* info, ftnlen
	uplo_len);

doublereal dla_porpvgrw__(char* uplo, integer* ncols, doublereal* a, integer*
	lda, doublereal* af, integer* ldaf, doublereal* work, ftnlen uplo_len);

doublereal dla_rpvgrw__(integer* n, integer* ncols, doublereal* a, integer*
	lda, doublereal* af, integer* ldaf);

int dla_syamv__(integer* uplo, integer* n, doublereal* alpha,
	doublereal* a, integer* lda, doublereal* x, integer* incx,
	doublereal* beta, doublereal* y, integer* incy);

doublereal dla_syrcond__(char* uplo, integer* n, doublereal* a, integer* lda,
	doublereal* af, integer* ldaf, integer* ipiv, integer* cmode,
	doublereal* c__, integer* info, doublereal* work, integer* iwork,
	ftnlen uplo_len);

int dla_syrfsx_extended__(integer* prec_type__, char* uplo,
	integer* n, integer* nrhs, doublereal* a, integer* lda, doublereal*
	af, integer* ldaf, integer* ipiv, logical* colequ, doublereal* c__,
	doublereal* b, integer* ldb, doublereal* y, integer* ldy, doublereal*
	berr_out__, integer* n_norms__, doublereal* errs_n__, doublereal*
	errs_c__, doublereal* res, doublereal* ayb, doublereal* dy,
	doublereal* y_tail__, doublereal* rcond, integer* ithresh, doublereal
	* rthresh, doublereal* dz_ub__, logical* ignore_cwise__, integer* info,
	ftnlen uplo_len);

doublereal dla_syrpvgrw__(char* uplo, integer* n, integer* info, doublereal*
	a, integer* lda, doublereal* af, integer* ldaf, integer* ipiv,
	doublereal* work, ftnlen uplo_len);

int dla_wwaddw__(integer* n, doublereal* x, doublereal* y,
	doublereal* w);

int dlarscl2_(integer* m, integer* n, doublereal* d__,
	doublereal* x, integer* ldx);

int dlascl2_(integer* m, integer* n, doublereal* d__,
	doublereal* x, integer* ldx);

int dporfsx_(char* uplo, char* equed, integer* n, integer*
	nrhs, doublereal* a, integer* lda, doublereal* af, integer* ldaf,
	doublereal* s, doublereal* b, integer* ldb, doublereal* x, integer*
	ldx, doublereal* rcond, doublereal* berr, integer* n_err_bnds__,
	doublereal* err_bnds_norm__, doublereal* err_bnds_comp__, integer*
	nparams, doublereal* params, doublereal* work, integer* iwork,
	integer* info);

int dposvxx_(char* fact, char* uplo, integer* n, integer*
	nrhs, doublereal* a, integer* lda, doublereal* af, integer* ldaf,
	char* equed, doublereal* s, doublereal* b, integer* ldb, doublereal*
	x, integer* ldx, doublereal* rcond, doublereal* rpvgrw, doublereal*
	berr, integer* n_err_bnds__, doublereal* err_bnds_norm__, doublereal*
	err_bnds_comp__, integer* nparams, doublereal* params, doublereal*
	work, integer* iwork, integer* info);

int dsyrfsx_(char* uplo, char* equed, integer* n, integer*
	nrhs, doublereal* a, integer* lda, doublereal* af, integer* ldaf,
	integer* ipiv, doublereal* s, doublereal* b, integer* ldb, doublereal
	* x, integer* ldx, doublereal* rcond, doublereal* berr, integer*
	n_err_bnds__, doublereal* err_bnds_norm__, doublereal*
	err_bnds_comp__, integer* nparams, doublereal* params, doublereal*
	work, integer* iwork, integer* info);

int dsysvxx_(char* fact, char* uplo, integer* n, integer*
	nrhs, doublereal* a, integer* lda, doublereal* af, integer* ldaf,
	integer* ipiv, char* equed, doublereal* s, doublereal* b, integer*
	ldb, doublereal* x, integer* ldx, doublereal* rcond, doublereal*
	rpvgrw, doublereal* berr, integer* n_err_bnds__, doublereal*
	err_bnds_norm__, doublereal* err_bnds_comp__, integer* nparams,
	doublereal* params, doublereal* work, integer* iwork, integer* info);
*/



#pragma endregion


// S
#pragma region SLASRC -- Single precision real LAPACK routines (sla) - ok

CLAPACK_API
integer ilaslc_(integer* m, integer* n, real* a, integer* lda);

CLAPACK_API
integer ilaslr_(integer* m, integer* n, real* a, integer* lda);

CLAPACK_API
doublereal scsum1_(integer* n, complex* cx, integer* incx);

CLAPACK_API int sgbbrd_(char* vect, integer* m, integer* n, integer* ncc,
	integer* kl, integer* ku, real* ab, integer* ldab, real* d__, real*
	e, real* q, integer* ldq, real* pt, integer* ldpt, real* c__, integer
	* ldc, real* work, integer* info);

CLAPACK_API int sgbcon_(char* norm, integer* n, integer* kl, integer* ku,
	real* ab, integer* ldab, integer* ipiv, real* anorm, real* rcond,
	real* work, integer* iwork, integer* info);

CLAPACK_API int sgbequ_(integer* m, integer* n, integer* kl, integer* ku,
	real* ab, integer* ldab, real* r__, real* c__, real* rowcnd, real*
	colcnd, real* amax, integer* info);

CLAPACK_API int sgbequb_(integer* m, integer* n, integer* kl, integer*
	ku, real* ab, integer* ldab, real* r__, real* c__, real* rowcnd, real
	* colcnd, real* amax, integer* info);

CLAPACK_API int sgbrfs_(char* trans, integer* n, integer* kl, integer*
	ku, integer* nrhs, real* ab, integer* ldab, real* afb, integer* ldafb,
	integer* ipiv, real* b, integer* ldb, real* x, integer* ldx, real*
	ferr, real* berr, real* work, integer* iwork, integer* info);

CLAPACK_API int sgbsv_(integer* n, integer* kl, integer* ku, integer*
	nrhs, real* ab, integer* ldab, integer* ipiv, real* b, integer* ldb,
	integer* info);

CLAPACK_API int sgbsvx_(char* fact, char* trans, integer* n, integer* kl,
	integer* ku, integer* nrhs, real* ab, integer* ldab, real* afb,
	integer* ldafb, integer* ipiv, char* equed, real* r__, real* c__,
	real* b, integer* ldb, real* x, integer* ldx, real* rcond, real* ferr,
	real* berr, real* work, integer* iwork, integer* info);

CLAPACK_API int sgbtf2_(integer* m, integer* n, integer* kl, integer* ku,
	real* ab, integer* ldab, integer* ipiv, integer* info);

CLAPACK_API int sgbtrf_(integer* m, integer* n, integer* kl, integer* ku,
	real* ab, integer* ldab, integer* ipiv, integer* info);

CLAPACK_API int sgbtrs_(char* trans, integer* n, integer* kl, integer*
	ku, integer* nrhs, real* ab, integer* ldab, integer* ipiv, real* b,
	integer* ldb, integer* info);

CLAPACK_API int sgebak_(char* job, char* side, integer* n, integer* ilo,
	integer* ihi, real* scale, integer* m, real* v, integer* ldv, integer
	* info);

CLAPACK_API int sgebal_(char* job, integer* n, real* a, integer* lda,
	integer* ilo, integer* ihi, real* scale, integer* info);

CLAPACK_API int sgebd2_(integer* m, integer* n, real* a, integer* lda,
	real* d__, real* e, real* tauq, real* taup, real* work, integer* info);

CLAPACK_API int sgebrd_(integer* m, integer* n, real* a, integer* lda,
	real* d__, real* e, real* tauq, real* taup, real* work, integer*
	lwork, integer* info);

CLAPACK_API int sgecon_(char* norm, integer* n, real* a, integer* lda,
	real* anorm, real* rcond, real* work, integer* iwork, integer* info);

CLAPACK_API int sgeequ_(integer* m, integer* n, real* a, integer* lda,
	real* r__, real* c__, real* rowcnd, real* colcnd, real* amax, integer
	* info);

CLAPACK_API int sgeequb_(integer* m, integer* n, real* a, integer* lda,
	real* r__, real* c__, real* rowcnd, real* colcnd, real* amax, integer
	* info);

CLAPACK_API int sgees_(char* jobvs, char* sort, L_fp select, integer* n,
	real* a, integer* lda, integer* sdim, real* wr, real* wi, real* vs,
	integer* ldvs, real* work, integer* lwork, logical* bwork, integer*
	info);

CLAPACK_API int sgeesx_(char* jobvs, char* sort, L_fp select, char*
	sense, integer* n, real* a, integer* lda, integer* sdim, real* wr,
	real* wi, real* vs, integer* ldvs, real* rconde, real* rcondv, real*
	work, integer* lwork, integer* iwork, integer* liwork, logical* bwork,
	integer* info);

CLAPACK_API int sgeev_(char* jobvl, char* jobvr, integer* n, real* a,
	integer* lda, real* wr, real* wi, real* vl, integer* ldvl, real* vr,
	integer* ldvr, real* work, integer* lwork, integer* info);

CLAPACK_API int sgeevx_(char* balanc, char* jobvl, char* jobvr, char*
	sense, integer* n, real* a, integer* lda, real* wr, real* wi, real*
	vl, integer* ldvl, real* vr, integer* ldvr, integer* ilo, integer*
	ihi, real* scale, real* abnrm, real* rconde, real* rcondv, real* work,
	integer* lwork, integer* iwork, integer* info);

CLAPACK_API int sgegs_(char* jobvsl, char* jobvsr, integer* n, real* a,
	integer* lda, real* b, integer* ldb, real* alphar, real* alphai, real
	* beta, real* vsl, integer* ldvsl, real* vsr, integer* ldvsr, real*
	work, integer* lwork, integer* info);

CLAPACK_API int sgegv_(char* jobvl, char* jobvr, integer* n, real* a,
	integer* lda, real* b, integer* ldb, real* alphar, real* alphai, real
	* beta, real* vl, integer* ldvl, real* vr, integer* ldvr, real* work,
	integer* lwork, integer* info);

CLAPACK_API int sgehd2_(integer* n, integer* ilo, integer* ihi, real* a,
	integer* lda, real* tau, real* work, integer* info);

CLAPACK_API int sgehrd_(integer* n, integer* ilo, integer* ihi, real* a,
	integer* lda, real* tau, real* work, integer* lwork, integer* info);

CLAPACK_API int sgejsv_(char* joba, char* jobu, char* jobv, char* jobr,
	char* jobt, char* jobp, integer* m, integer* n, real* a, integer* lda,
	real* sva, real* u, integer* ldu, real* v, integer* ldv, real* work,
	integer* lwork, integer* iwork, integer* info);

CLAPACK_API int sgelq2_(integer* m, integer* n, real* a, integer* lda,
	real* tau, real* work, integer* info);

CLAPACK_API int sgelqf_(integer* m, integer* n, real* a, integer* lda,
	real* tau, real* work, integer* lwork, integer* info);

CLAPACK_API int sgels_(char* trans, integer* m, integer* n, integer*
	nrhs, real* a, integer* lda, real* b, integer* ldb, real* work,
	integer* lwork, integer* info);

CLAPACK_API int sgelsd_(integer* m, integer* n, integer* nrhs, real* a,
	integer* lda, real* b, integer* ldb, real* s, real* rcond, integer*
	rank, real* work, integer* lwork, integer* iwork, integer* info);

CLAPACK_API int sgelss_(integer* m, integer* n, integer* nrhs, real* a,
	integer* lda, real* b, integer* ldb, real* s, real* rcond, integer*
	rank, real* work, integer* lwork, integer* info);

CLAPACK_API int sgelsx_(integer* m, integer* n, integer* nrhs, real* a,
	integer* lda, real* b, integer* ldb, integer* jpvt, real* rcond,
	integer* rank, real* work, integer* info);

CLAPACK_API int sgelsy_(integer* m, integer* n, integer* nrhs, real* a,
	integer* lda, real* b, integer* ldb, integer* jpvt, real* rcond,
	integer* rank, real* work, integer* lwork, integer* info);

CLAPACK_API int sgeql2_(integer* m, integer* n, real* a, integer* lda,
	real* tau, real* work, integer* info);

CLAPACK_API int sgeqlf_(integer* m, integer* n, real* a, integer* lda,
	real* tau, real* work, integer* lwork, integer* info);

CLAPACK_API int sgeqp3_(integer* m, integer* n, real* a, integer* lda,
	integer* jpvt, real* tau, real* work, integer* lwork, integer* info);

CLAPACK_API int sgeqpf_(integer* m, integer* n, real* a, integer* lda,
	integer* jpvt, real* tau, real* work, integer* info);

CLAPACK_API int sgeqr2_(integer* m, integer* n, real* a, integer* lda,
	real* tau, real* work, integer* info);

CLAPACK_API int sgeqrf_(integer* m, integer* n, real* a, integer* lda,
	real* tau, real* work, integer* lwork, integer* info);

CLAPACK_API int sgerfs_(char* trans, integer* n, integer* nrhs, real* a,
	integer* lda, real* af, integer* ldaf, integer* ipiv, real* b,
	integer* ldb, real* x, integer* ldx, real* ferr, real* berr, real*
	work, integer* iwork, integer* info);

CLAPACK_API int sgerq2_(integer* m, integer* n, real* a, integer* lda,
	real* tau, real* work, integer* info);

CLAPACK_API int sgerqf_(integer* m, integer* n, real* a, integer* lda,
	real* tau, real* work, integer* lwork, integer* info);

CLAPACK_API int sgesc2_(integer* n, real* a, integer* lda, real* rhs,
	integer* ipiv, integer* jpiv, real* scale);

CLAPACK_API int sgesdd_(char* jobz, integer* m, integer* n, real* a,
	integer* lda, real* s, real* u, integer* ldu, real* vt, integer* ldvt,
	real* work, integer* lwork, integer* iwork, integer* info);

CLAPACK_API int sgesv_(integer* n, integer* nrhs, real* a, integer* lda,
	integer* ipiv, real* b, integer* ldb, integer* info);

CLAPACK_API int sgesvd_(char* jobu, char* jobvt, integer* m, integer* n,
	real* a, integer* lda, real* s, real* u, integer* ldu, real* vt,
	integer* ldvt, real* work, integer* lwork, integer* info);

CLAPACK_API int sgesvj_(char* joba, char* jobu, char* jobv, integer* m,
	integer* n, real* a, integer* lda, real* sva, integer* mv, real* v,
	integer* ldv, real* work, integer* lwork, integer* info);

CLAPACK_API int sgesvx_(char* fact, char* trans, integer* n, integer*
	nrhs, real* a, integer* lda, real* af, integer* ldaf, integer* ipiv,
	char* equed, real* r__, real* c__, real* b, integer* ldb, real* x,
	integer* ldx, real* rcond, real* ferr, real* berr, real* work,
	integer* iwork, integer* info);

CLAPACK_API int sgetc2_(integer* n, real* a, integer* lda, integer* ipiv,
	integer* jpiv, integer* info);

CLAPACK_API int sgetf2_(integer* m, integer* n, real* a, integer* lda,
	integer* ipiv, integer* info);

CLAPACK_API int sgetrf_(integer* m, integer* n, real* a, integer* lda,
	integer* ipiv, integer* info);

CLAPACK_API int sgetri_(integer* n, real* a, integer* lda, integer* ipiv,
	real* work, integer* lwork, integer* info);

CLAPACK_API int sgetrs_(char* trans, integer* n, integer* nrhs, real* a,
	integer* lda, integer* ipiv, real* b, integer* ldb, integer* info);

CLAPACK_API int sggbak_(char* job, char* side, integer* n, integer* ilo,
	integer* ihi, real* lscale, real* rscale, integer* m, real* v,
	integer* ldv, integer* info);

CLAPACK_API int sggbal_(char* job, integer* n, real* a, integer* lda,
	real* b, integer* ldb, integer* ilo, integer* ihi, real* lscale, real
	* rscale, real* work, integer* info);

CLAPACK_API int sgges_(char* jobvsl, char* jobvsr, char* sort, L_fp
	selctg, integer* n, real* a, integer* lda, real* b, integer* ldb,
	integer* sdim, real* alphar, real* alphai, real* beta, real* vsl,
	integer* ldvsl, real* vsr, integer* ldvsr, real* work, integer* lwork,
	logical* bwork, integer* info);

CLAPACK_API int sggesx_(char* jobvsl, char* jobvsr, char* sort, L_fp
	selctg, char* sense, integer* n, real* a, integer* lda, real* b,
	integer* ldb, integer* sdim, real* alphar, real* alphai, real* beta,
	real* vsl, integer* ldvsl, real* vsr, integer* ldvsr, real* rconde,
	real* rcondv, real* work, integer* lwork, integer* iwork, integer*
	liwork, logical* bwork, integer* info);

CLAPACK_API int sggev_(char* jobvl, char* jobvr, integer* n, real* a,
	integer* lda, real* b, integer* ldb, real* alphar, real* alphai, real
	* beta, real* vl, integer* ldvl, real* vr, integer* ldvr, real* work,
	integer* lwork, integer* info);

CLAPACK_API int sggevx_(char* balanc, char* jobvl, char* jobvr, char*
	sense, integer* n, real* a, integer* lda, real* b, integer* ldb, real
	* alphar, real* alphai, real* beta, real* vl, integer* ldvl, real* vr,
	integer* ldvr, integer* ilo, integer* ihi, real* lscale, real* rscale,
	real* abnrm, real* bbnrm, real* rconde, real* rcondv, real* work,
	integer* lwork, integer* iwork, logical* bwork, integer* info);

CLAPACK_API int sggglm_(integer* n, integer* m, integer* p, real* a,
	integer* lda, real* b, integer* ldb, real* d__, real* x, real* y,
	real* work, integer* lwork, integer* info);

CLAPACK_API int sgghrd_(char* compq, char* compz, integer* n, integer*
	ilo, integer* ihi, real* a, integer* lda, real* b, integer* ldb, real
	* q, integer* ldq, real* z__, integer* ldz, integer* info);

CLAPACK_API int sgglse_(integer* m, integer* n, integer* p, real* a,
	integer* lda, real* b, integer* ldb, real* c__, real* d__, real* x,
	real* work, integer* lwork, integer* info);

CLAPACK_API int sggqrf_(integer* n, integer* m, integer* p, real* a,
	integer* lda, real* taua, real* b, integer* ldb, real* taub, real*
	work, integer* lwork, integer* info);

CLAPACK_API int sggrqf_(integer* m, integer* p, integer* n, real* a,
	integer* lda, real* taua, real* b, integer* ldb, real* taub, real*
	work, integer* lwork, integer* info);

CLAPACK_API int sggsvd_(char* jobu, char* jobv, char* jobq, integer* m,
	integer* n, integer* p, integer* k, integer* l, real* a, integer* lda,
	real* b, integer* ldb, real* alpha, real* beta, real* u, integer*
	ldu, real* v, integer* ldv, real* q, integer* ldq, real* work,
	integer* iwork, integer* info);

CLAPACK_API int sggsvp_(char* jobu, char* jobv, char* jobq, integer* m,
	integer* p, integer* n, real* a, integer* lda, real* b, integer* ldb,
	real* tola, real* tolb, integer* k, integer* l, real* u, integer* ldu,
	real* v, integer* ldv, real* q, integer* ldq, integer* iwork, real*
	tau, real* work, integer* info);

CLAPACK_API int sgsvj0_(char* jobv, integer* m, integer* n, real* a,
	integer* lda, real* d__, real* sva, integer* mv, real* v, integer*
	ldv, real* eps, real* sfmin, real* tol, integer* nsweep, real* work,
	integer* lwork, integer* info);

CLAPACK_API int sgsvj1_(char* jobv, integer* m, integer* n, integer* n1,
	real* a, integer* lda, real* d__, real* sva, integer* mv, real* v,
	integer* ldv, real* eps, real* sfmin, real* tol, integer* nsweep,
	real* work, integer* lwork, integer* info);

CLAPACK_API int sgtcon_(char* norm, integer* n, real* dl, real* d__,
	real* du, real* du2, integer* ipiv, real* anorm, real* rcond, real*
	work, integer* iwork, integer* info);

CLAPACK_API int sgtrfs_(char* trans, integer* n, integer* nrhs, real* dl,
	real* d__, real* du, real* dlf, real* df, real* duf, real* du2,
	integer* ipiv, real* b, integer* ldb, real* x, integer* ldx, real*
	ferr, real* berr, real* work, integer* iwork, integer* info);

CLAPACK_API int sgtsv_(integer* n, integer* nrhs, real* dl, real* d__,
	real* du, real* b, integer* ldb, integer* info);

CLAPACK_API int sgtsvx_(char* fact, char* trans, integer* n, integer*
	nrhs, real* dl, real* d__, real* du, real* dlf, real* df, real* duf,
	real* du2, integer* ipiv, real* b, integer* ldb, real* x, integer*
	ldx, real* rcond, real* ferr, real* berr, real* work, integer* iwork,
	integer* info);

CLAPACK_API int sgttrf_(integer* n, real* dl, real* d__, real* du, real*
	du2, integer* ipiv, integer* info);

CLAPACK_API int sgttrs_(char* trans, integer* n, integer* nrhs, real* dl,
	real* d__, real* du, real* du2, integer* ipiv, real* b, integer* ldb,
	integer* info);

CLAPACK_API int sgtts2_(integer* itrans, integer* n, integer* nrhs, real
	* dl, real* d__, real* du, real* du2, integer* ipiv, real* b, integer*
	ldb);

CLAPACK_API int shgeqz_(char* job, char* compq, char* compz, integer* n,
	integer* ilo, integer* ihi, real* h__, integer* ldh, real* t, integer
	* ldt, real* alphar, real* alphai, real* beta, real* q, integer* ldq,
	real* z__, integer* ldz, real* work, integer* lwork, integer* info);

CLAPACK_API int shsein_(char* side, char* eigsrc, char* initv, logical*
	select, integer* n, real* h__, integer* ldh, real* wr, real* wi, real
	* vl, integer* ldvl, real* vr, integer* ldvr, integer* mm, integer* m,
	real* work, integer* ifaill, integer* ifailr, integer* info);

CLAPACK_API int shseqr_(char* job, char* compz, integer* n, integer* ilo,
	integer* ihi, real* h__, integer* ldh, real* wr, real* wi, real* z__,
	integer* ldz, real* work, integer* lwork, integer* info);

CLAPACK_API int slabrd_(integer* m, integer* n, integer* nb, real* a,
	integer* lda, real* d__, real* e, real* tauq, real* taup, real* x,
	integer* ldx, real* y, integer* ldy);

CLAPACK_API int slacn2_(integer* n, real* v, real* x, integer* isgn,
	real* est, integer* kase, integer* isave);

CLAPACK_API int slacon_(integer* n, real* v, real* x, integer* isgn,
	real* est, integer* kase);

CLAPACK_API int slaein_(logical* rightv, logical* noinit, integer* n,
	real* h__, integer* ldh, real* wr, real* wi, real* vr, real* vi, real
	* b, integer* ldb, real* work, real* eps3, real* smlnum, real* bignum,
	integer* info);

CLAPACK_API int slaexc_(logical* wantq, integer* n, real* t, integer*
	ldt, real* q, integer* ldq, integer* j1, integer* n1, integer* n2,
	real* work, integer* info);

CLAPACK_API int slag2_(real* a, integer* lda, real* b, integer* ldb,
	real* safmin, real* scale1, real* scale2, real* wr1, real* wr2, real*
	wi);

CLAPACK_API int slag2d_(integer* m, integer* n, real* sa, integer* ldsa,
	doublereal* a, integer* lda, integer* info);

CLAPACK_API int slags2_(logical* upper, real* a1, real* a2, real* a3,
	real* b1, real* b2, real* b3, real* csu, real* snu, real* csv, real*
	snv, real* csq, real* snq);

CLAPACK_API int slagtm_(char* trans, integer* n, integer* nrhs, real*
	alpha, real* dl, real* d__, real* du, real* x, integer* ldx, real*
	beta, real* b, integer* ldb);

CLAPACK_API int slagv2_(real* a, integer* lda, real* b, integer* ldb,
	real* alphar, real* alphai, real* beta, real* csl, real* snl, real*
	csr, real* snr);

CLAPACK_API int slahqr_(logical* wantt, logical* wantz, integer* n,
	integer* ilo, integer* ihi, real* h__, integer* ldh, real* wr, real*
	wi, integer* iloz, integer* ihiz, real* z__, integer* ldz, integer*
	info);

CLAPACK_API int slahr2_(integer* n, integer* k, integer* nb, real* a,
	integer* lda, real* tau, real* t, integer* ldt, real* y, integer* ldy);

CLAPACK_API int slahrd_(integer* n, integer* k, integer* nb, real* a,
	integer* lda, real* tau, real* t, integer* ldt, real* y, integer* ldy);

CLAPACK_API int slaic1_(integer* job, integer* j, real* x, real* sest,
	real* w, real* gamma, real* sestpr, real* s, real* c__);

CLAPACK_API int slaln2_(logical* ltrans, integer* na, integer* nw, real*
	smin, real* ca, real* a, integer* lda, real* d1, real* d2, real* b,
	integer* ldb, real* wr, real* wi, real* x, integer* ldx, real* scale,
	real* xnorm, integer* info);

CLAPACK_API int slals0_(integer* icompq, integer* nl, integer* nr,
	integer* sqre, integer* nrhs, real* b, integer* ldb, real* bx,
	integer* ldbx, integer* perm, integer* givptr, integer* givcol,
	integer* ldgcol, real* givnum, integer* ldgnum, real* poles, real*
	difl, real* difr, real* z__, integer* k, real* c__, real* s, real*
	work, integer* info);

CLAPACK_API int slalsa_(integer* icompq, integer* smlsiz, integer* n,
	integer* nrhs, real* b, integer* ldb, real* bx, integer* ldbx, real*
	u, integer* ldu, real* vt, integer* k, real* difl, real* difr, real*
	z__, real* poles, integer* givptr, integer* givcol, integer* ldgcol,
	integer* perm, real* givnum, real* c__, real* s, real* work, integer*
	iwork, integer* info);

CLAPACK_API int slalsd_(char* uplo, integer* smlsiz, integer* n, integer
	* nrhs, real* d__, real* e, real* b, integer* ldb, real* rcond,
	integer* rank, real* work, integer* iwork, integer* info);

CLAPACK_API doublereal slangb_(char* norm, integer* n, integer* kl, integer* ku, real* ab,
	integer* ldab, real* work);

CLAPACK_API doublereal slange_(char* norm, integer* m, integer* n, real* a, integer* lda,
	real* work);

CLAPACK_API doublereal slangt_(char* norm, integer* n, real* dl, real* d__, real* du);

CLAPACK_API doublereal slanhs_(char* norm, integer* n, real* a, integer* lda, real* work);

CLAPACK_API doublereal slansb_(char* norm, char* uplo, integer* n, integer* k, real* ab,
	integer* ldab, real* work);

CLAPACK_API doublereal slansf_(char* norm, char* transr, char* uplo, integer* n, real* a,
	real* work);

CLAPACK_API doublereal slansp_(char* norm, char* uplo, integer* n, real* ap, real* work);

CLAPACK_API doublereal slansy_(char* norm, char* uplo, integer* n, real* a, integer* lda,
	real* work);

CLAPACK_API doublereal slantb_(char* norm, char* uplo, char* diag, integer* n, integer* k,
	real* ab, integer* ldab, real* work);

CLAPACK_API doublereal slantp_(char* norm, char* uplo, char* diag, integer* n, real* ap,
	real* work);

CLAPACK_API doublereal slantr_(char* norm, char* uplo, char* diag, integer* m, integer* n,
	real* a, integer* lda, real* work);

CLAPACK_API int slanv2_(real* a, real* b, real* c__, real* d__, real*
	rt1r, real* rt1i, real* rt2r, real* rt2i, real* cs, real* sn);

CLAPACK_API int slapll_(integer* n, real* x, integer* incx, real* y,
	integer* incy, real* ssmin);

CLAPACK_API int slapmt_(logical* forwrd, integer* m, integer* n, real* x,
	integer* ldx, integer* k);

CLAPACK_API int slaqgb_(integer* m, integer* n, integer* kl, integer* ku,
	real* ab, integer* ldab, real* r__, real* c__, real* rowcnd, real*
	colcnd, real* amax, char* equed);

CLAPACK_API int slaqge_(integer* m, integer* n, real* a, integer* lda,
	real* r__, real* c__, real* rowcnd, real* colcnd, real* amax, char*
	equed);

CLAPACK_API int slaqp2_(integer* m, integer* n, integer* offset, real* a,
	integer* lda, integer* jpvt, real* tau, real* vn1, real* vn2, real*
	work);

CLAPACK_API int slaqps_(integer* m, integer* n, integer* offset, integer
	* nb, integer* kb, real* a, integer* lda, integer* jpvt, real* tau,
	real* vn1, real* vn2, real* auxv, real* f, integer* ldf);

CLAPACK_API int slaqr0_(logical* wantt, logical* wantz, integer* n,
	integer* ilo, integer* ihi, real* h__, integer* ldh, real* wr, real*
	wi, integer* iloz, integer* ihiz, real* z__, integer* ldz, real* work,
	integer* lwork, integer* info);

CLAPACK_API int slaqr1_(integer* n, real* h__, integer* ldh, real* sr1,
	real* si1, real* sr2, real* si2, real* v);

CLAPACK_API int slaqr2_(logical* wantt, logical* wantz, integer* n,
	integer* ktop, integer* kbot, integer* nw, real* h__, integer* ldh,
	integer* iloz, integer* ihiz, real* z__, integer* ldz, integer* ns,
	integer* nd, real* sr, real* si, real* v, integer* ldv, integer* nh,
	real* t, integer* ldt, integer* nv, real* wv, integer* ldwv, real*
	work, integer* lwork);

CLAPACK_API int slaqr3_(logical* wantt, logical* wantz, integer* n,
	integer* ktop, integer* kbot, integer* nw, real* h__, integer* ldh,
	integer* iloz, integer* ihiz, real* z__, integer* ldz, integer* ns,
	integer* nd, real* sr, real* si, real* v, integer* ldv, integer* nh,
	real* t, integer* ldt, integer* nv, real* wv, integer* ldwv, real*
	work, integer* lwork);

CLAPACK_API int slaqr4_(logical* wantt, logical* wantz, integer* n,
	integer* ilo, integer* ihi, real* h__, integer* ldh, real* wr, real*
	wi, integer* iloz, integer* ihiz, real* z__, integer* ldz, real* work,
	integer* lwork, integer* info);

CLAPACK_API int slaqr5_(logical* wantt, logical* wantz, integer* kacc22,
	integer* n, integer* ktop, integer* kbot, integer* nshfts, real* sr,
	real* si, real* h__, integer* ldh, integer* iloz, integer* ihiz, real
	* z__, integer* ldz, real* v, integer* ldv, real* u, integer* ldu,
	integer* nv, real* wv, integer* ldwv, integer* nh, real* wh, integer*
	ldwh);

CLAPACK_API int slaqsb_(char* uplo, integer* n, integer* kd, real* ab,
	integer* ldab, real* s, real* scond, real* amax, char* equed);

CLAPACK_API int slaqsp_(char* uplo, integer* n, real* ap, real* s, real*
	scond, real* amax, char* equed);

CLAPACK_API int slaqsy_(char* uplo, integer* n, real* a, integer* lda,
	real* s, real* scond, real* amax, char* equed);

CLAPACK_API int slaqtr_(logical* ltran, logical* lreal, integer* n, real
	* t, integer* ldt, real* b, real* w, real* scale, real* x, real* work,
	integer* info);

CLAPACK_API int slar1v_(integer* n, integer* b1, integer* bn, real*
	lambda, real* d__, real* l, real* ld, real* lld, real* pivmin, real*
	gaptol, real* z__, logical* wantnc, integer* negcnt, real* ztz, real*
	mingma, integer* r__, integer* isuppz, real* nrminv, real* resid,
	real* rqcorr, real* work);

CLAPACK_API int slar2v_(integer* n, real* x, real* y, real* z__, integer
	* incx, real* c__, real* s, integer* incc);

CLAPACK_API int slarf_(char* side, integer* m, integer* n, real* v,
	integer* incv, real* tau, real* c__, integer* ldc, real* work);

CLAPACK_API int slarfb_(char* side, char* trans, char* direct, char*
	storev, integer* m, integer* n, integer* k, real* v, integer* ldv,
	real* t, integer* ldt, real* c__, integer* ldc, real* work, integer*
	ldwork);

CLAPACK_API int slarfg_(integer* n, real* alpha, real* x, integer* incx,
	real* tau);

CLAPACK_API int slarfp_(integer* n, real* alpha, real* x, integer* incx,
	real* tau);

CLAPACK_API int slarft_(char* direct, char* storev, integer* n, integer*
	k, real* v, integer* ldv, real* tau, real* t, integer* ldt);

CLAPACK_API int slarfx_(char* side, integer* m, integer* n, real* v,
	real* tau, real* c__, integer* ldc, real* work);

CLAPACK_API int slargv_(integer* n, real* x, integer* incx, real* y,
	integer* incy, real* c__, integer* incc);

CLAPACK_API int slarrv_(integer* n, real* vl, real* vu, real* d__, real*
	l, real* pivmin, integer* isplit, integer* m, integer* dol, integer*
	dou, real* minrgp, real* rtol1, real* rtol2, real* w, real* werr,
	real* wgap, integer* iblock, integer* indexw, real* gers, real* z__,
	integer* ldz, integer* isuppz, real* work, integer* iwork, integer*
	info);

CLAPACK_API int slarscl2_(integer* m, integer* n, real* d__, real* x,
	integer* ldx);

CLAPACK_API int slartv_(integer* n, real* x, integer* incx, real* y,
	integer* incy, real* c__, real* s, integer* incc);

CLAPACK_API int slarz_(char* side, integer* m, integer* n, integer* l,
	real* v, integer* incv, real* tau, real* c__, integer* ldc, real*
	work);

CLAPACK_API int slarzb_(char* side, char* trans, char* direct, char*
	storev, integer* m, integer* n, integer* k, integer* l, real* v,
	integer* ldv, real* t, integer* ldt, real* c__, integer* ldc, real*
	work, integer* ldwork);

CLAPACK_API int slarzt_(char* direct, char* storev, integer* n, integer*
	k, real* v, integer* ldv, real* tau, real* t, integer* ldt);

CLAPACK_API int slascl2_(integer* m, integer* n, real* d__, real* x,
	integer* ldx);

CLAPACK_API int slaswp_(integer* n, real* a, integer* lda, integer* k1,
	integer* k2, integer* ipiv, integer* incx);

CLAPACK_API int slasy2_(logical* ltranl, logical* ltranr, integer* isgn,
	integer* n1, integer* n2, real* tl, integer* ldtl, real* tr, integer*
	ldtr, real* b, integer* ldb, real* scale, real* x, integer* ldx, real
	* xnorm, integer* info);

CLAPACK_API int slasyf_(char* uplo, integer* n, integer* nb, integer* kb,
	real* a, integer* lda, integer* ipiv, real* w, integer* ldw, integer
	* info);

CLAPACK_API int slatbs_(char* uplo, char* trans, char* diag, char*
	normin, integer* n, integer* kd, real* ab, integer* ldab, real* x,
	real* scale, real* cnorm, integer* info);

CLAPACK_API int slatdf_(integer* ijob, integer* n, real* z__, integer*
	ldz, real* rhs, real* rdsum, real* rdscal, integer* ipiv, integer*
	jpiv);

CLAPACK_API int slatps_(char* uplo, char* trans, char* diag, char*
	normin, integer* n, real* ap, real* x, real* scale, real* cnorm,
	integer* info);

CLAPACK_API int slatrd_(char* uplo, integer* n, integer* nb, real* a,
	integer* lda, real* e, real* tau, real* w, integer* ldw);

CLAPACK_API int slatrs_(char* uplo, char* trans, char* diag, char*
	normin, integer* n, real* a, integer* lda, real* x, real* scale, real
	* cnorm, integer* info);

CLAPACK_API int slatrz_(integer* m, integer* n, integer* l, real* a,
	integer* lda, real* tau, real* work);

CLAPACK_API int slatzm_(char* side, integer* m, integer* n, real* v,
	integer* incv, real* tau, real* c1, real* c2, integer* ldc, real*
	work);

CLAPACK_API int slauu2_(char* uplo, integer* n, real* a, integer* lda,
	integer* info);

CLAPACK_API int slauum_(char* uplo, integer* n, real* a, integer* lda,
	integer* info);

CLAPACK_API int sopgtr_(char* uplo, integer* n, real* ap, real* tau,
	real* q, integer* ldq, real* work, integer* info);

CLAPACK_API int sopmtr_(char* side, char* uplo, char* trans, integer* m,
	integer* n, real* ap, real* tau, real* c__, integer* ldc, real* work,
	integer* info);

CLAPACK_API int sorg2l_(integer* m, integer* n, integer* k, real* a,
	integer* lda, real* tau, real* work, integer* info);

CLAPACK_API int sorg2r_(integer* m, integer* n, integer* k, real* a,
	integer* lda, real* tau, real* work, integer* info);

CLAPACK_API int sorgbr_(char* vect, integer* m, integer* n, integer* k,
	real* a, integer* lda, real* tau, real* work, integer* lwork, integer
	* info);

CLAPACK_API int sorghr_(integer* n, integer* ilo, integer* ihi, real* a,
	integer* lda, real* tau, real* work, integer* lwork, integer* info);

CLAPACK_API int sorgl2_(integer* m, integer* n, integer* k, real* a,
	integer* lda, real* tau, real* work, integer* info);

CLAPACK_API int sorglq_(integer* m, integer* n, integer* k, real* a,
	integer* lda, real* tau, real* work, integer* lwork, integer* info);

CLAPACK_API int sorgql_(integer* m, integer* n, integer* k, real* a,
	integer* lda, real* tau, real* work, integer* lwork, integer* info);

CLAPACK_API int sorgqr_(integer* m, integer* n, integer* k, real* a,
	integer* lda, real* tau, real* work, integer* lwork, integer* info);

CLAPACK_API int sorgr2_(integer* m, integer* n, integer* k, real* a,
	integer* lda, real* tau, real* work, integer* info);

CLAPACK_API int sorgrq_(integer* m, integer* n, integer* k, real* a,
	integer* lda, real* tau, real* work, integer* lwork, integer* info);

CLAPACK_API int sorgtr_(char* uplo, integer* n, real* a, integer* lda,
	real* tau, real* work, integer* lwork, integer* info);

CLAPACK_API int sorm2l_(char* side, char* trans, integer* m, integer* n,
	integer* k, real* a, integer* lda, real* tau, real* c__, integer* ldc,
	real* work, integer* info);

CLAPACK_API int sorm2r_(char* side, char* trans, integer* m, integer* n,
	integer* k, real* a, integer* lda, real* tau, real* c__, integer* ldc,
	real* work, integer* info);

CLAPACK_API int sormbr_(char* vect, char* side, char* trans, integer* m,
	integer* n, integer* k, real* a, integer* lda, real* tau, real* c__,
	integer* ldc, real* work, integer* lwork, integer* info);

CLAPACK_API int sormhr_(char* side, char* trans, integer* m, integer* n,
	integer* ilo, integer* ihi, real* a, integer* lda, real* tau, real*
	c__, integer* ldc, real* work, integer* lwork, integer* info);

CLAPACK_API int sorml2_(char* side, char* trans, integer* m, integer* n,
	integer* k, real* a, integer* lda, real* tau, real* c__, integer* ldc,
	real* work, integer* info);

CLAPACK_API int sormlq_(char* side, char* trans, integer* m, integer* n,
	integer* k, real* a, integer* lda, real* tau, real* c__, integer* ldc,
	real* work, integer* lwork, integer* info);

CLAPACK_API int sormql_(char* side, char* trans, integer* m, integer* n,
	integer* k, real* a, integer* lda, real* tau, real* c__, integer* ldc,
	real* work, integer* lwork, integer* info);

CLAPACK_API int sormqr_(char* side, char* trans, integer* m, integer* n,
	integer* k, real* a, integer* lda, real* tau, real* c__, integer* ldc,
	real* work, integer* lwork, integer* info);

CLAPACK_API int sormr2_(char* side, char* trans, integer* m, integer* n,
	integer* k, real* a, integer* lda, real* tau, real* c__, integer* ldc,
	real* work, integer* info);

CLAPACK_API int sormr3_(char* side, char* trans, integer* m, integer* n,
	integer* k, integer* l, real* a, integer* lda, real* tau, real* c__,
	integer* ldc, real* work, integer* info);

CLAPACK_API int sormrq_(char* side, char* trans, integer* m, integer* n,
	integer* k, real* a, integer* lda, real* tau, real* c__, integer* ldc,
	real* work, integer* lwork, integer* info);

CLAPACK_API int sormrz_(char* side, char* trans, integer* m, integer* n,
	integer* k, integer* l, real* a, integer* lda, real* tau, real* c__,
	integer* ldc, real* work, integer* lwork, integer* info);

CLAPACK_API int sormtr_(char* side, char* uplo, char* trans, integer* m,
	integer* n, real* a, integer* lda, real* tau, real* c__, integer* ldc,
	real* work, integer* lwork, integer* info);

CLAPACK_API int spbcon_(char* uplo, integer* n, integer* kd, real* ab,
	integer* ldab, real* anorm, real* rcond, real* work, integer* iwork,
	integer* info);

CLAPACK_API int spbequ_(char* uplo, integer* n, integer* kd, real* ab,
	integer* ldab, real* s, real* scond, real* amax, integer* info);

CLAPACK_API int spbrfs_(char* uplo, integer* n, integer* kd, integer*
	nrhs, real* ab, integer* ldab, real* afb, integer* ldafb, real* b,
	integer* ldb, real* x, integer* ldx, real* ferr, real* berr, real*
	work, integer* iwork, integer* info);

CLAPACK_API int spbstf_(char* uplo, integer* n, integer* kd, real* ab,
	integer* ldab, integer* info);

CLAPACK_API int spbsv_(char* uplo, integer* n, integer* kd, integer*
	nrhs, real* ab, integer* ldab, real* b, integer* ldb, integer* info);

CLAPACK_API int spbsvx_(char* fact, char* uplo, integer* n, integer* kd,
	integer* nrhs, real* ab, integer* ldab, real* afb, integer* ldafb,
	char* equed, real* s, real* b, integer* ldb, real* x, integer* ldx,
	real* rcond, real* ferr, real* berr, real* work, integer* iwork,
	integer* info);

CLAPACK_API int spbtf2_(char* uplo, integer* n, integer* kd, real* ab,
	integer* ldab, integer* info);

CLAPACK_API int spbtrf_(char* uplo, integer* n, integer* kd, real* ab,
	integer* ldab, integer* info);

CLAPACK_API int spbtrs_(char* uplo, integer* n, integer* kd, integer*
	nrhs, real* ab, integer* ldab, real* b, integer* ldb, integer* info);

CLAPACK_API int spftrf_(char* transr, char* uplo, integer* n, real* a,
	integer* info);

CLAPACK_API int spftri_(char* transr, char* uplo, integer* n, real* a,
	integer* info);

CLAPACK_API int spftrs_(char* transr, char* uplo, integer* n, integer*
	nrhs, real* a, real* b, integer* ldb, integer* info);

CLAPACK_API int spocon_(char* uplo, integer* n, real* a, integer* lda,
	real* anorm, real* rcond, real* work, integer* iwork, integer* info);

CLAPACK_API int spoequ_(integer* n, real* a, integer* lda, real* s, real
	* scond, real* amax, integer* info);

CLAPACK_API int spoequb_(integer* n, real* a, integer* lda, real* s,
	real* scond, real* amax, integer* info);

CLAPACK_API int sporfs_(char* uplo, integer* n, integer* nrhs, real* a,
	integer* lda, real* af, integer* ldaf, real* b, integer* ldb, real* x,
	integer* ldx, real* ferr, real* berr, real* work, integer* iwork,
	integer* info);

CLAPACK_API int sposv_(char* uplo, integer* n, integer* nrhs, real* a,
	integer* lda, real* b, integer* ldb, integer* info);

CLAPACK_API int sposvx_(char* fact, char* uplo, integer* n, integer*
	nrhs, real* a, integer* lda, real* af, integer* ldaf, char* equed,
	real* s, real* b, integer* ldb, real* x, integer* ldx, real* rcond,
	real* ferr, real* berr, real* work, integer* iwork, integer* info);

CLAPACK_API int spotf2_(char* uplo, integer* n, real* a, integer* lda,
	integer* info);

CLAPACK_API int spotrf_(char* uplo, integer* n, real* a, integer* lda,
	integer* info);

CLAPACK_API int spotri_(char* uplo, integer* n, real* a, integer* lda,
	integer* info);

CLAPACK_API int spotrs_(char* uplo, integer* n, integer* nrhs, real* a,
	integer* lda, real* b, integer* ldb, integer* info);

CLAPACK_API int sppcon_(char* uplo, integer* n, real* ap, real* anorm,
	real* rcond, real* work, integer* iwork, integer* info);

CLAPACK_API int sppequ_(char* uplo, integer* n, real* ap, real* s, real*
	scond, real* amax, integer* info);

CLAPACK_API int spprfs_(char* uplo, integer* n, integer* nrhs, real* ap,
	real* afp, real* b, integer* ldb, real* x, integer* ldx, real* ferr,
	real* berr, real* work, integer* iwork, integer* info);

CLAPACK_API int sppsv_(char* uplo, integer* n, integer* nrhs, real* ap,
	real* b, integer* ldb, integer* info);

CLAPACK_API int sppsvx_(char* fact, char* uplo, integer* n, integer*
	nrhs, real* ap, real* afp, char* equed, real* s, real* b, integer*
	ldb, real* x, integer* ldx, real* rcond, real* ferr, real* berr, real
	* work, integer* iwork, integer* info);

CLAPACK_API int spptrf_(char* uplo, integer* n, real* ap, integer* info);

CLAPACK_API int spptri_(char* uplo, integer* n, real* ap, integer* info);

CLAPACK_API int spptrs_(char* uplo, integer* n, integer* nrhs, real* ap,
	real* b, integer* ldb, integer* info);

CLAPACK_API int spstf2_(char* uplo, integer* n, real* a, integer* lda,
	integer* piv, integer* rank, real* tol, real* work, integer* info);

CLAPACK_API int spstrf_(char* uplo, integer* n, real* a, integer* lda,
	integer* piv, integer* rank, real* tol, real* work, integer* info);

CLAPACK_API int sptcon_(integer* n, real* d__, real* e, real* anorm,
	real* rcond, real* work, integer* info);

CLAPACK_API int spteqr_(char* compz, integer* n, real* d__, real* e,
	real* z__, integer* ldz, real* work, integer* info);

CLAPACK_API int sptrfs_(integer* n, integer* nrhs, real* d__, real* e,
	real* df, real* ef, real* b, integer* ldb, real* x, integer* ldx,
	real* ferr, real* berr, real* work, integer* info);

CLAPACK_API int sptsv_(integer* n, integer* nrhs, real* d__, real* e,
	real* b, integer* ldb, integer* info);

CLAPACK_API int sptsvx_(char* fact, integer* n, integer* nrhs, real* d__,
	real* e, real* df, real* ef, real* b, integer* ldb, real* x, integer
	* ldx, real* rcond, real* ferr, real* berr, real* work, integer* info);

CLAPACK_API int spttrs_(integer* n, integer* nrhs, real* d__, real* e,
	real* b, integer* ldb, integer* info);

CLAPACK_API int sptts2_(integer* n, integer* nrhs, real* d__, real* e,
	real* b, integer* ldb);

CLAPACK_API int srscl_(integer* n, real* sa, real* sx, integer* incx);

CLAPACK_API int ssbev_(char* jobz, char* uplo, integer* n, integer* kd,
	real* ab, integer* ldab, real* w, real* z__, integer* ldz, real* work,
	integer* info);

CLAPACK_API int ssbevd_(char* jobz, char* uplo, integer* n, integer* kd,
	real* ab, integer* ldab, real* w, real* z__, integer* ldz, real* work,
	integer* lwork, integer* iwork, integer* liwork, integer* info);

CLAPACK_API int ssbevx_(char* jobz, char* range, char* uplo, integer* n,
	integer* kd, real* ab, integer* ldab, real* q, integer* ldq, real* vl,
	real* vu, integer* il, integer* iu, real* abstol, integer* m, real*
	w, real* z__, integer* ldz, real* work, integer* iwork, integer*
	ifail, integer* info);

CLAPACK_API int ssbgst_(char* vect, char* uplo, integer* n, integer* ka,
	integer* kb, real* ab, integer* ldab, real* bb, integer* ldbb, real*
	x, integer* ldx, real* work, integer* info);

CLAPACK_API int ssbgv_(char* jobz, char* uplo, integer* n, integer* ka,
	integer* kb, real* ab, integer* ldab, real* bb, integer* ldbb, real*
	w, real* z__, integer* ldz, real* work, integer* info);

CLAPACK_API int ssbgvd_(char* jobz, char* uplo, integer* n, integer* ka,
	integer* kb, real* ab, integer* ldab, real* bb, integer* ldbb, real*
	w, real* z__, integer* ldz, real* work, integer* lwork, integer*
	iwork, integer* liwork, integer* info);

CLAPACK_API int ssbgvx_(char* jobz, char* range, char* uplo, integer* n,
	integer* ka, integer* kb, real* ab, integer* ldab, real* bb, integer*
	ldbb, real* q, integer* ldq, real* vl, real* vu, integer* il, integer
	* iu, real* abstol, integer* m, real* w, real* z__, integer* ldz, real
	* work, integer* iwork, integer* ifail, integer* info);

CLAPACK_API int ssbtrd_(char* vect, char* uplo, integer* n, integer* kd,
	real* ab, integer* ldab, real* d__, real* e, real* q, integer* ldq,
	real* work, integer* info);

CLAPACK_API int ssfrk_(char* transr, char* uplo, char* trans, integer* n,
	integer* k, real* alpha, real* a, integer* lda, real* beta, real*
	c__);

CLAPACK_API int sspcon_(char* uplo, integer* n, real* ap, integer* ipiv,
	real* anorm, real* rcond, real* work, integer* iwork, integer* info);

CLAPACK_API int sspev_(char* jobz, char* uplo, integer* n, real* ap,
	real* w, real* z__, integer* ldz, real* work, integer* info);

CLAPACK_API int sspevd_(char* jobz, char* uplo, integer* n, real* ap,
	real* w, real* z__, integer* ldz, real* work, integer* lwork, integer
	* iwork, integer* liwork, integer* info);

CLAPACK_API int sspevx_(char* jobz, char* range, char* uplo, integer* n,
	real* ap, real* vl, real* vu, integer* il, integer* iu, real* abstol,
	integer* m, real* w, real* z__, integer* ldz, real* work, integer*
	iwork, integer* ifail, integer* info);

CLAPACK_API int sspgst_(integer* itype, char* uplo, integer* n, real* ap,
	real* bp, integer* info);

CLAPACK_API int sspgv_(integer* itype, char* jobz, char* uplo, integer*
	n, real* ap, real* bp, real* w, real* z__, integer* ldz, real* work,
	integer* info);

CLAPACK_API int sspgvd_(integer* itype, char* jobz, char* uplo, integer*
	n, real* ap, real* bp, real* w, real* z__, integer* ldz, real* work,
	integer* lwork, integer* iwork, integer* liwork, integer* info);

CLAPACK_API int sspgvx_(integer* itype, char* jobz, char* range, char*
	uplo, integer* n, real* ap, real* bp, real* vl, real* vu, integer* il,
	integer* iu, real* abstol, integer* m, real* w, real* z__, integer*
	ldz, real* work, integer* iwork, integer* ifail, integer* info);

CLAPACK_API int ssprfs_(char* uplo, integer* n, integer* nrhs, real* ap,
	real* afp, integer* ipiv, real* b, integer* ldb, real* x, integer*
	ldx, real* ferr, real* berr, real* work, integer* iwork, integer*
	info);

CLAPACK_API int sspsv_(char* uplo, integer* n, integer* nrhs, real* ap,
	integer* ipiv, real* b, integer* ldb, integer* info);

CLAPACK_API int sspsvx_(char* fact, char* uplo, integer* n, integer*
	nrhs, real* ap, real* afp, integer* ipiv, real* b, integer* ldb, real
	* x, integer* ldx, real* rcond, real* ferr, real* berr, real* work,
	integer* iwork, integer* info);

CLAPACK_API int ssptrd_(char* uplo, integer* n, real* ap, real* d__,
	real* e, real* tau, integer* info);

CLAPACK_API int ssptrf_(char* uplo, integer* n, real* ap, integer* ipiv,
	integer* info);

CLAPACK_API int ssptri_(char* uplo, integer* n, real* ap, integer* ipiv,
	real* work, integer* info);

CLAPACK_API int ssptrs_(char* uplo, integer* n, integer* nrhs, real* ap,
	integer* ipiv, real* b, integer* ldb, integer* info);

CLAPACK_API int sstegr_(char* jobz, char* range, integer* n, real* d__,
	real* e, real* vl, real* vu, integer* il, integer* iu, real* abstol,
	integer* m, real* w, real* z__, integer* ldz, integer* isuppz, real*
	work, integer* lwork, integer* iwork, integer* liwork, integer* info);

CLAPACK_API int sstein_(integer* n, real* d__, real* e, integer* m, real
	* w, integer* iblock, integer* isplit, real* z__, integer* ldz, real*
	work, integer* iwork, integer* ifail, integer* info);

CLAPACK_API int sstemr_(char* jobz, char* range, integer* n, real* d__,
	real* e, real* vl, real* vu, integer* il, integer* iu, integer* m,
	real* w, real* z__, integer* ldz, integer* nzc, integer* isuppz,
	logical* tryrac, real* work, integer* lwork, integer* iwork, integer*
	liwork, integer* info);

CLAPACK_API int sstev_(char* jobz, integer* n, real* d__, real* e, real*
	z__, integer* ldz, real* work, integer* info);

CLAPACK_API int sstevd_(char* jobz, integer* n, real* d__, real* e, real
	* z__, integer* ldz, real* work, integer* lwork, integer* iwork,
	integer* liwork, integer* info);

CLAPACK_API int sstevr_(char* jobz, char* range, integer* n, real* d__,
	real* e, real* vl, real* vu, integer* il, integer* iu, real* abstol,
	integer* m, real* w, real* z__, integer* ldz, integer* isuppz, real*
	work, integer* lwork, integer* iwork, integer* liwork, integer* info);

CLAPACK_API int sstevx_(char* jobz, char* range, integer* n, real* d__,
	real* e, real* vl, real* vu, integer* il, integer* iu, real* abstol,
	integer* m, real* w, real* z__, integer* ldz, real* work, integer*
	iwork, integer* ifail, integer* info);

CLAPACK_API int ssycon_(char* uplo, integer* n, real* a, integer* lda,
	integer* ipiv, real* anorm, real* rcond, real* work, integer* iwork,
	integer* info);

CLAPACK_API int ssyequb_(char* uplo, integer* n, real* a, integer* lda,
	real* s, real* scond, real* amax, real* work, integer* info);

CLAPACK_API int ssyev_(char* jobz, char* uplo, integer* n, real* a,
	integer* lda, real* w, real* work, integer* lwork, integer* info);

CLAPACK_API int ssyevd_(char* jobz, char* uplo, integer* n, real* a,
	integer* lda, real* w, real* work, integer* lwork, integer* iwork,
	integer* liwork, integer* info);

CLAPACK_API int ssyevr_(char* jobz, char* range, char* uplo, integer* n,
	real* a, integer* lda, real* vl, real* vu, integer* il, integer* iu,
	real* abstol, integer* m, real* w, real* z__, integer* ldz, integer*
	isuppz, real* work, integer* lwork, integer* iwork, integer* liwork,
	integer* info);

CLAPACK_API int ssyevx_(char* jobz, char* range, char* uplo, integer* n,
	real* a, integer* lda, real* vl, real* vu, integer* il, integer* iu,
	real* abstol, integer* m, real* w, real* z__, integer* ldz, real*
	work, integer* lwork, integer* iwork, integer* ifail, integer* info);

CLAPACK_API int ssygs2_(integer* itype, char* uplo, integer* n, real* a,
	integer* lda, real* b, integer* ldb, integer* info);

CLAPACK_API int ssygst_(integer* itype, char* uplo, integer* n, real* a,
	integer* lda, real* b, integer* ldb, integer* info);

CLAPACK_API int ssygv_(integer* itype, char* jobz, char* uplo, integer*
	n, real* a, integer* lda, real* b, integer* ldb, real* w, real* work,
	integer* lwork, integer* info);

CLAPACK_API int ssygvd_(integer* itype, char* jobz, char* uplo, integer*
	n, real* a, integer* lda, real* b, integer* ldb, real* w, real* work,
	integer* lwork, integer* iwork, integer* liwork, integer* info);

CLAPACK_API int ssygvx_(integer* itype, char* jobz, char* range, char*
	uplo, integer* n, real* a, integer* lda, real* b, integer* ldb, real*
	vl, real* vu, integer* il, integer* iu, real* abstol, integer* m,
	real* w, real* z__, integer* ldz, real* work, integer* lwork, integer
	* iwork, integer* ifail, integer* info);

CLAPACK_API int ssyrfs_(char* uplo, integer* n, integer* nrhs, real* a,
	integer* lda, real* af, integer* ldaf, integer* ipiv, real* b,
	integer* ldb, real* x, integer* ldx, real* ferr, real* berr, real*
	work, integer* iwork, integer* info);

CLAPACK_API int ssysv_(char* uplo, integer* n, integer* nrhs, real* a,
	integer* lda, integer* ipiv, real* b, integer* ldb, real* work,
	integer* lwork, integer* info);

CLAPACK_API int ssysvx_(char* fact, char* uplo, integer* n, integer*
	nrhs, real* a, integer* lda, real* af, integer* ldaf, integer* ipiv,
	real* b, integer* ldb, real* x, integer* ldx, real* rcond, real* ferr,
	real* berr, real* work, integer* lwork, integer* iwork, integer*
	info);

CLAPACK_API int ssytd2_(char* uplo, integer* n, real* a, integer* lda,
	real* d__, real* e, real* tau, integer* info);

CLAPACK_API int ssytf2_(char* uplo, integer* n, real* a, integer* lda,
	integer* ipiv, integer* info);

CLAPACK_API int ssytrd_(char* uplo, integer* n, real* a, integer* lda,
	real* d__, real* e, real* tau, real* work, integer* lwork, integer*
	info);

CLAPACK_API int ssytrf_(char* uplo, integer* n, real* a, integer* lda,
	integer* ipiv, real* work, integer* lwork, integer* info);

CLAPACK_API int ssytri_(char* uplo, integer* n, real* a, integer* lda,
	integer* ipiv, real* work, integer* info);

CLAPACK_API int ssytrs_(char* uplo, integer* n, integer* nrhs, real* a,
	integer* lda, integer* ipiv, real* b, integer* ldb, integer* info);

CLAPACK_API int stbcon_(char* norm, char* uplo, char* diag, integer* n,
	integer* kd, real* ab, integer* ldab, real* rcond, real* work,
	integer* iwork, integer* info);

CLAPACK_API int stbrfs_(char* uplo, char* trans, char* diag, integer* n,
	integer* kd, integer* nrhs, real* ab, integer* ldab, real* b, integer
	* ldb, real* x, integer* ldx, real* ferr, real* berr, real* work,
	integer* iwork, integer* info);

CLAPACK_API int stbtrs_(char* uplo, char* trans, char* diag, integer* n,
	integer* kd, integer* nrhs, real* ab, integer* ldab, real* b, integer
	* ldb, integer* info);

CLAPACK_API int stfsm_(char* transr, char* side, char* uplo, char* trans,
	char* diag, integer* m, integer* n, real* alpha, real* a, real* b,
	integer* ldb);

CLAPACK_API int stftri_(char* transr, char* uplo, char* diag, integer* n,
	real* a, integer* info);

CLAPACK_API int stfttp_(char* transr, char* uplo, integer* n, real* arf,
	real* ap, integer* info);

CLAPACK_API int stfttr_(char* transr, char* uplo, integer* n, real* arf,
	real* a, integer* lda, integer* info);

CLAPACK_API int stgevc_(char* side, char* howmny, logical* select,
	integer* n, real* s, integer* lds, real* p, integer* ldp, real* vl,
	integer* ldvl, real* vr, integer* ldvr, integer* mm, integer* m, real
	* work, integer* info);

CLAPACK_API int stgex2_(logical* wantq, logical* wantz, integer* n, real
	* a, integer* lda, real* b, integer* ldb, real* q, integer* ldq, real*
	z__, integer* ldz, integer* j1, integer* n1, integer* n2, real* work,
	integer* lwork, integer* info);

CLAPACK_API int stgexc_(logical* wantq, logical* wantz, integer* n, real
	* a, integer* lda, real* b, integer* ldb, real* q, integer* ldq, real*
	z__, integer* ldz, integer* ifst, integer* ilst, real* work, integer*
	lwork, integer* info);

CLAPACK_API int stgsen_(integer* ijob, logical* wantq, logical* wantz,
	logical* select, integer* n, real* a, integer* lda, real* b, integer*
	ldb, real* alphar, real* alphai, real* beta, real* q, integer* ldq,
	real* z__, integer* ldz, integer* m, real* pl, real* pr, real* dif,
	real* work, integer* lwork, integer* iwork, integer* liwork, integer*
	info);

CLAPACK_API int stgsja_(char* jobu, char* jobv, char* jobq, integer* m,
	integer* p, integer* n, integer* k, integer* l, real* a, integer* lda,
	real* b, integer* ldb, real* tola, real* tolb, real* alpha, real*
	beta, real* u, integer* ldu, real* v, integer* ldv, real* q, integer*
	ldq, real* work, integer* ncycle, integer* info);

CLAPACK_API int stgsna_(char* job, char* howmny, logical* select,
	integer* n, real* a, integer* lda, real* b, integer* ldb, real* vl,
	integer* ldvl, real* vr, integer* ldvr, real* s, real* dif, integer*
	mm, integer* m, real* work, integer* lwork, integer* iwork, integer*
	info);

CLAPACK_API int stgsy2_(char* trans, integer* ijob, integer* m, integer*
	n, real* a, integer* lda, real* b, integer* ldb, real* c__, integer*
	ldc, real* d__, integer* ldd, real* e, integer* lde, real* f, integer
	* ldf, real* scale, real* rdsum, real* rdscal, integer* iwork, integer
	* pq, integer* info);

CLAPACK_API int stgsyl_(char* trans, integer* ijob, integer* m, integer*
	n, real* a, integer* lda, real* b, integer* ldb, real* c__, integer*
	ldc, real* d__, integer* ldd, real* e, integer* lde, real* f, integer
	* ldf, real* scale, real* dif, real* work, integer* lwork, integer*
	iwork, integer* info);

CLAPACK_API int stpcon_(char* norm, char* uplo, char* diag, integer* n,
	real* ap, real* rcond, real* work, integer* iwork, integer* info);

CLAPACK_API int stprfs_(char* uplo, char* trans, char* diag, integer* n,
	integer* nrhs, real* ap, real* b, integer* ldb, real* x, integer* ldx,
	real* ferr, real* berr, real* work, integer* iwork, integer* info);

CLAPACK_API int stptri_(char* uplo, char* diag, integer* n, real* ap,
	integer* info);

CLAPACK_API int stptrs_(char* uplo, char* trans, char* diag, integer* n,
	integer* nrhs, real* ap, real* b, integer* ldb, integer* info);

CLAPACK_API int stpttf_(char* transr, char* uplo, integer* n, real* ap,
	real* arf, integer* info);

CLAPACK_API int stpttr_(char* uplo, integer* n, real* ap, real* a,
	integer* lda, integer* info);

CLAPACK_API int strcon_(char* norm, char* uplo, char* diag, integer* n,
	real* a, integer* lda, real* rcond, real* work, integer* iwork,
	integer* info);

CLAPACK_API int strevc_(char* side, char* howmny, logical* select,
	integer* n, real* t, integer* ldt, real* vl, integer* ldvl, real* vr,
	integer* ldvr, integer* mm, integer* m, real* work, integer* info);

CLAPACK_API int strexc_(char* compq, integer* n, real* t, integer* ldt,
	real* q, integer* ldq, integer* ifst, integer* ilst, real* work,
	integer* info);

CLAPACK_API int strrfs_(char* uplo, char* trans, char* diag, integer* n,
	integer* nrhs, real* a, integer* lda, real* b, integer* ldb, real* x,
	integer* ldx, real* ferr, real* berr, real* work, integer* iwork,
	integer* info);

CLAPACK_API int strsen_(char* job, char* compq, logical* select, integer
	* n, real* t, integer* ldt, real* q, integer* ldq, real* wr, real* wi,
	integer* m, real* s, real* sep, real* work, integer* lwork, integer*
	iwork, integer* liwork, integer* info);

CLAPACK_API int strsna_(char* job, char* howmny, logical* select,
	integer* n, real* t, integer* ldt, real* vl, integer* ldvl, real* vr,
	integer* ldvr, real* s, real* sep, integer* mm, integer* m, real*
	work, integer* ldwork, integer* iwork, integer* info);

CLAPACK_API int strsyl_(char* trana, char* tranb, integer* isgn, integer
	* m, integer* n, real* a, integer* lda, real* b, integer* ldb, real*
	c__, integer* ldc, real* scale, integer* info);

CLAPACK_API int strti2_(char* uplo, char* diag, integer* n, real* a,
	integer* lda, integer* info);

CLAPACK_API int strtri_(char* uplo, char* diag, integer* n, real* a,
	integer* lda, integer* info);

CLAPACK_API int strtrs_(char* uplo, char* trans, char* diag, integer* n,
	integer* nrhs, real* a, integer* lda, real* b, integer* ldb, integer*
	info);

CLAPACK_API int strttf_(char* transr, char* uplo, integer* n, real* a,
	integer* lda, real* arf, integer* info);

CLAPACK_API int strttp_(char* uplo, integer* n, real* a, integer* lda,
	real* ap, integer* info);

CLAPACK_API int stzrqf_(integer* m, integer* n, real* a, integer* lda,
	real* tau, integer* info);

CLAPACK_API int stzrzf_(integer* m, integer* n, real* a, integer* lda,
	real* tau, real* work, integer* lwork, integer* info);


//}} finished

#pragma endregion

#pragma region SXLASRC -- Single precision real LAPACK routines using extra precision (slax) - ok
//#       
//set(SXLASRC sgesvxx.c sgerfsx.c
//	sla_gerfsx_extended.c sla_geamv.c
//	sla_gercond.c sla_rpvgrw.c ssysvxx.c ssyrfsx.c
//	sla_syrfsx_extended.c sla_syamv.c sla_syrcond.c sla_syrpvgrw.c
//	sposvxx.c sporfsx.c sla_porfsx_extended.c sla_porcond.c
//	sla_porpvgrw.c sgbsvxx.c sgbrfsx.c sla_gbrfsx_extended.c
//	sla_gbamv.c sla_gbrcond.c sla_gbrpvgrw.c sla_lin_berr.c slarscl2.c
//	slascl2.c sla_wwaddw.c)

/*

int sgbrfsx_(char *trans, char *equed, integer *n, integer *
	kl, integer *ku, integer *nrhs, real *ab, integer *ldab, real *afb,
	integer *ldafb, integer *ipiv, real *r__, real *c__, real *b, integer
	*ldb, real *x, integer *ldx, real *rcond, real *berr, integer *
	n_err_bnds__, real *err_bnds_norm__, real *err_bnds_comp__, integer *
	nparams, real *params, real *work, integer *iwork, integer *info);

int sgbsvxx_(char *fact, char *trans, integer *n, integer *
	kl, integer *ku, integer *nrhs, real *ab, integer *ldab, real *afb,
	integer *ldafb, integer *ipiv, char *equed, real *r__, real *c__,
	real *b, integer *ldb, real *x, integer *ldx, real *rcond, real *
	rpvgrw, real *berr, integer *n_err_bnds__, real *err_bnds_norm__,
	real *err_bnds_comp__, integer *nparams, real *params, real *work,
	integer *iwork, integer *info);

int sgerfsx_(char* trans, char* equed, integer* n, integer*
	nrhs, real* a, integer* lda, real* af, integer* ldaf, integer* ipiv,
	real* r__, real* c__, real* b, integer* ldb, real* x, integer* ldx,
	real* rcond, real* berr, integer* n_err_bnds__, real* err_bnds_norm__,
	real* err_bnds_comp__, integer* nparams, real* params, real* work,
	integer* iwork, integer* info);

int sgesvxx_(char* fact, char* trans, integer* n, integer*
	nrhs, real* a, integer* lda, real* af, integer* ldaf, integer* ipiv,
	char* equed, real* r__, real* c__, real* b, integer* ldb, real* x,
	integer* ldx, real* rcond, real* rpvgrw, real* berr, integer*
	n_err_bnds__, real* err_bnds_norm__, real* err_bnds_comp__, integer*
	nparams, real* params, real* work, integer* iwork, integer* info);

int sla_gbamv__(integer *trans, integer *m, integer *n,
	integer *kl, integer *ku, real *alpha, real *ab, integer *ldab, real *
	x, integer *incx, real *beta, real *y, integer *incy);

doublereal sla_gbrcond__(char *trans, integer *n, integer *kl, integer *ku,
	real *ab, integer *ldab, real *afb, integer *ldafb, integer *ipiv,
	integer *cmode, real *c__, integer *info, real *work, integer *iwork,
	ftnlen trans_len);

int sla_gbrfsx_extended__(integer *prec_type__, integer *
	trans_type__, integer *n, integer *kl, integer *ku, integer *nrhs,
	real *ab, integer *ldab, real *afb, integer *ldafb, integer *ipiv,
	logical *colequ, real *c__, real *b, integer *ldb, real *y, integer *
	ldy, real *berr_out__, integer *n_norms__, real *errs_n__, real *
	errs_c__, real *res, real *ayb, real *dy, real *y_tail__, real *rcond,
	 integer *ithresh, real *rthresh, real *dz_ub__, logical *
	ignore_cwise__, integer *info);

doublereal sla_gbrpvgrw__(integer *n, integer *kl, integer *ku, integer *
	ncols, real *ab, integer *ldab, real *afb, integer *ldafb);

int sla_geamv__(integer *trans, integer *m, integer *n, real
	*alpha, real *a, integer *lda, real *x, integer *incx, real *beta,
	real *y, integer *incy);

doublereal sla_gercond__(char *trans, integer *n, real *a, integer *lda, real
	*af, integer *ldaf, integer *ipiv, integer *cmode, real *c__, integer
	*info, real *work, integer *iwork, ftnlen trans_len);

int sla_gerfsx_extended__(integer *prec_type__, integer *
	trans_type__, integer *n, integer *nrhs, real *a, integer *lda, real *
	af, integer *ldaf, integer *ipiv, logical *colequ, real *c__, real *b,
	 integer *ldb, real *y, integer *ldy, real *berr_out__, integer *
	n_norms__, real *errs_n__, real *errs_c__, real *res, real *ayb, real
	*dy, real *y_tail__, real *rcond, integer *ithresh, real *rthresh,
	real *dz_ub__, logical *ignore_cwise__, integer *info);

int sla_lin_berr__(integer *n, integer *nz, integer *nrhs,
	real *res, real *ayb, real *berr);

doublereal sla_porcond__(char *uplo, integer *n, real *a, integer *lda, real *
	af, integer *ldaf, integer *cmode, real *c__, integer *info, real *
	work, integer *iwork, ftnlen uplo_len);

int sla_porfsx_extended__(integer *prec_type__, char *uplo,
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

int sla_syamv__(integer *uplo, integer *n, real *alpha, real
	*a, integer *lda, real *x, integer *incx, real *beta, real *y,
	integer *incy);

doublereal sla_syrcond__(char *uplo, integer *n, real *a, integer *lda, real *
	af, integer *ldaf, integer *ipiv, integer *cmode, real *c__, integer *
	info, real *work, integer *iwork, ftnlen uplo_len);

int sla_syrfsx_extended__(integer *prec_type__, char *uplo,
	integer *n, integer *nrhs, real *a, integer *lda, real *af, integer *
	ldaf, integer *ipiv, logical *colequ, real *c__, real *b, integer *
	ldb, real *y, integer *ldy, real *berr_out__, integer *n_norms__,
	real *errs_n__, real *errs_c__, real *res, real *ayb, real *dy, real *
	y_tail__, real *rcond, integer *ithresh, real *rthresh, real *dz_ub__,
	 logical *ignore_cwise__, integer *info, ftnlen uplo_len);

doublereal sla_syrpvgrw__(char *uplo, integer *n, integer *info, real *a,
	integer *lda, real *af, integer *ldaf, integer *ipiv, real *work,
	ftnlen uplo_len);

int sla_wwaddw__(integer *n, real *x, real *y, real *w);

int sporfsx_(char *uplo, char *equed, integer *n, integer *
	nrhs, real *a, integer *lda, real *af, integer *ldaf, real *s, real *
	b, integer *ldb, real *x, integer *ldx, real *rcond, real *berr,
	integer *n_err_bnds__, real *err_bnds_norm__, real *err_bnds_comp__,
	integer *nparams, real *params, real *work, integer *iwork, integer *
	info);

int sposvxx_(char *fact, char *uplo, integer *n, integer *
	nrhs, real *a, integer *lda, real *af, integer *ldaf, char *equed,
	real *s, real *b, integer *ldb, real *x, integer *ldx, real *rcond,
	real *rpvgrw, real *berr, integer *n_err_bnds__, real *
	err_bnds_norm__, real *err_bnds_comp__, integer *nparams, real *
	params, real *work, integer *iwork, integer *info);

int ssyrfsx_(char *uplo, char *equed, integer *n, integer *
	nrhs, real *a, integer *lda, real *af, integer *ldaf, integer *ipiv,
	real *s, real *b, integer *ldb, real *x, integer *ldx, real *rcond,
	real *berr, integer *n_err_bnds__, real *err_bnds_norm__, real *
	err_bnds_comp__, integer *nparams, real *params, real *work, integer *
	iwork, integer *info);

int ssysvxx_(char *fact, char *uplo, integer *n, integer *
	nrhs, real *a, integer *lda, real *af, integer *ldaf, integer *ipiv,
	char *equed, real *s, real *b, integer *ldb, real *x, integer *ldx,
	real *rcond, real *rpvgrw, real *berr, integer *n_err_bnds__, real *
	err_bnds_norm__, real *err_bnds_comp__, integer *nparams, real *
	params, real *work, integer *iwork, integer *info);



*/

#pragma endregion


// Z
#pragma region ZLASRC -- Double precision complex LAPACK routines (zla) - ok

CLAPACK_API
doublereal dzsum1_(integer* n, doublecomplex* cx, integer* incx);

CLAPACK_API
integer ilazlc_(integer* m, integer* n, doublecomplex* a, integer* lda);

CLAPACK_API
integer ilazlr_(integer* m, integer* n, doublecomplex* a, integer* lda);

CLAPACK_API
integer izmax1_(integer* n, doublecomplex* cx, integer* incx);

CLAPACK_API int zbdsqr_(char* uplo, integer* n, integer* ncvt, integer*
	nru, integer* ncc, doublereal* d__, doublereal* e, doublecomplex* vt,
	integer* ldvt, doublecomplex* u, integer* ldu, doublecomplex* c__,
	integer* ldc, doublereal* rwork, integer* info);

CLAPACK_API int zcgesv_(integer* n, integer* nrhs, doublecomplex* a,
	integer* lda, integer* ipiv, doublecomplex* b, integer* ldb,
	doublecomplex* x, integer* ldx, doublecomplex* work, complex* swork,
	doublereal* rwork, integer* iter, integer* info);

CLAPACK_API int zcposv_(char* uplo, integer* n, integer* nrhs,
	doublecomplex* a, integer* lda, doublecomplex* b, integer* ldb,
	doublecomplex* x, integer* ldx, doublecomplex* work, complex* swork,
	doublereal* rwork, integer* iter, integer* info);

CLAPACK_API int zdrscl_(integer* n, doublereal* sa, doublecomplex* sx,
	integer* incx);

CLAPACK_API int zgbcon_(char* norm, integer* n, integer* kl, integer* ku,
	doublecomplex* ab, integer* ldab, integer* ipiv, doublereal* anorm,
	doublereal* rcond, doublecomplex* work, doublereal* rwork, integer*
	info);

CLAPACK_API int zgbequ_(integer* m, integer* n, integer* kl, integer* ku,
	doublecomplex* ab, integer* ldab, doublereal* r__, doublereal* c__,
	doublereal* rowcnd, doublereal* colcnd, doublereal* amax, integer*
	info);

CLAPACK_API int zgbequb_(integer* m, integer* n, integer* kl, integer*
	ku, doublecomplex* ab, integer* ldab, doublereal* r__, doublereal*
	c__, doublereal* rowcnd, doublereal* colcnd, doublereal* amax,
	integer* info);

CLAPACK_API int zgbrfs_(char* trans, integer* n, integer* kl, integer*
	ku, integer* nrhs, doublecomplex* ab, integer* ldab, doublecomplex*
	afb, integer* ldafb, integer* ipiv, doublecomplex* b, integer* ldb,
	doublecomplex* x, integer* ldx, doublereal* ferr, doublereal* berr,
	doublecomplex* work, doublereal* rwork, integer* info);

CLAPACK_API int zgbsv_(integer* n, integer* kl, integer* ku, integer*
	nrhs, doublecomplex* ab, integer* ldab, integer* ipiv, doublecomplex*
	b, integer* ldb, integer* info);

CLAPACK_API int zgbsvx_(char* fact, char* trans, integer* n, integer* kl,
	integer* ku, integer* nrhs, doublecomplex* ab, integer* ldab,
	doublecomplex* afb, integer* ldafb, integer* ipiv, char* equed,
	doublereal* r__, doublereal* c__, doublecomplex* b, integer* ldb,
	doublecomplex* x, integer* ldx, doublereal* rcond, doublereal* ferr,
	doublereal* berr, doublecomplex* work, doublereal* rwork, integer*
	info);

CLAPACK_API int zgbtf2_(integer* m, integer* n, integer* kl, integer* ku,
	doublecomplex* ab, integer* ldab, integer* ipiv, integer* info);

CLAPACK_API int zgbtrf_(integer* m, integer* n, integer* kl, integer* ku,
	doublecomplex* ab, integer* ldab, integer* ipiv, integer* info);

CLAPACK_API int zgbtrs_(char* trans, integer* n, integer* kl, integer*
	ku, integer* nrhs, doublecomplex* ab, integer* ldab, integer* ipiv,
	doublecomplex* b, integer* ldb, integer* info);

CLAPACK_API int zgebak_(char* job, char* side, integer* n, integer* ilo,
	integer* ihi, doublereal* scale, integer* m, doublecomplex* v,
	integer* ldv, integer* info);

CLAPACK_API int zgebal_(char* job, integer* n, doublecomplex* a, integer
	* lda, integer* ilo, integer* ihi, doublereal* scale, integer* info);

CLAPACK_API int zgebd2_(integer* m, integer* n, doublecomplex* a,
	integer* lda, doublereal* d__, doublereal* e, doublecomplex* tauq,
	doublecomplex* taup, doublecomplex* work, integer* info);

CLAPACK_API int zgebrd_(integer* m, integer* n, doublecomplex* a,
	integer* lda, doublereal* d__, doublereal* e, doublecomplex* tauq,
	doublecomplex* taup, doublecomplex* work, integer* lwork, integer*
	info);

CLAPACK_API int zgecon_(char* norm, integer* n, doublecomplex* a,
	integer* lda, doublereal* anorm, doublereal* rcond, doublecomplex*
	work, doublereal* rwork, integer* info);

CLAPACK_API int zgeequ_(integer* m, integer* n, doublecomplex* a,
	integer* lda, doublereal* r__, doublereal* c__, doublereal* rowcnd,
	doublereal* colcnd, doublereal* amax, integer* info);

CLAPACK_API int zgeequb_(integer* m, integer* n, doublecomplex* a,
	integer* lda, doublereal* r__, doublereal* c__, doublereal* rowcnd,
	doublereal* colcnd, doublereal* amax, integer* info);

CLAPACK_API int zgees_(char* jobvs, char* sort, L_fp select, integer* n,
	doublecomplex* a, integer* lda, integer* sdim, doublecomplex* w,
	doublecomplex* vs, integer* ldvs, doublecomplex* work, integer* lwork,
	doublereal* rwork, logical* bwork, integer* info);

CLAPACK_API int zgeesx_(char* jobvs, char* sort, L_fp select, char*
	sense, integer* n, doublecomplex* a, integer* lda, integer* sdim,
	doublecomplex* w, doublecomplex* vs, integer* ldvs, doublereal*
	rconde, doublereal* rcondv, doublecomplex* work, integer* lwork,
	doublereal* rwork, logical* bwork, integer* info);

CLAPACK_API int zgeev_(char* jobvl, char* jobvr, integer* n,
	doublecomplex* a, integer* lda, doublecomplex* w, doublecomplex* vl,
	integer* ldvl, doublecomplex* vr, integer* ldvr, doublecomplex* work,
	integer* lwork, doublereal* rwork, integer* info);

CLAPACK_API int zgeevx_(char* balanc, char* jobvl, char* jobvr, char*
	sense, integer* n, doublecomplex* a, integer* lda, doublecomplex* w,
	doublecomplex* vl, integer* ldvl, doublecomplex* vr, integer* ldvr,
	integer* ilo, integer* ihi, doublereal* scale, doublereal* abnrm,
	doublereal* rconde, doublereal* rcondv, doublecomplex* work, integer*
	lwork, doublereal* rwork, integer* info);

CLAPACK_API int zgegs_(char* jobvsl, char* jobvsr, integer* n,
	doublecomplex* a, integer* lda, doublecomplex* b, integer* ldb,
	doublecomplex* alpha, doublecomplex* beta, doublecomplex* vsl,
	integer* ldvsl, doublecomplex* vsr, integer* ldvsr, doublecomplex*
	work, integer* lwork, doublereal* rwork, integer* info);

CLAPACK_API int zgegv_(char* jobvl, char* jobvr, integer* n,
	doublecomplex* a, integer* lda, doublecomplex* b, integer* ldb,
	doublecomplex* alpha, doublecomplex* beta, doublecomplex* vl, integer
	* ldvl, doublecomplex* vr, integer* ldvr, doublecomplex* work, integer
	* lwork, doublereal* rwork, integer* info);

CLAPACK_API int zgehd2_(integer* n, integer* ilo, integer* ihi,
	doublecomplex* a, integer* lda, doublecomplex* tau, doublecomplex*
	work, integer* info);

CLAPACK_API int zgehrd_(integer* n, integer* ilo, integer* ihi,
	doublecomplex* a, integer* lda, doublecomplex* tau, doublecomplex*
	work, integer* lwork, integer* info);

CLAPACK_API int zgelq2_(integer* m, integer* n, doublecomplex* a,
	integer* lda, doublecomplex* tau, doublecomplex* work, integer* info);

CLAPACK_API int zgelqf_(integer* m, integer* n, doublecomplex* a,
	integer* lda, doublecomplex* tau, doublecomplex* work, integer* lwork,
	integer* info);

CLAPACK_API int zgels_(char* trans, integer* m, integer* n, integer*
	nrhs, doublecomplex* a, integer* lda, doublecomplex* b, integer* ldb,
	doublecomplex* work, integer* lwork, integer* info);

CLAPACK_API int zgelsd_(integer* m, integer* n, integer* nrhs,
	doublecomplex* a, integer* lda, doublecomplex* b, integer* ldb,
	doublereal* s, doublereal* rcond, integer* rank, doublecomplex* work,
	integer* lwork, doublereal* rwork, integer* iwork, integer* info);

CLAPACK_API int zgelss_(integer* m, integer* n, integer* nrhs,
	doublecomplex* a, integer* lda, doublecomplex* b, integer* ldb,
	doublereal* s, doublereal* rcond, integer* rank, doublecomplex* work,
	integer* lwork, doublereal* rwork, integer* info);

CLAPACK_API int zgelsx_(integer* m, integer* n, integer* nrhs,
	doublecomplex* a, integer* lda, doublecomplex* b, integer* ldb,
	integer* jpvt, doublereal* rcond, integer* rank, doublecomplex* work,
	doublereal* rwork, integer* info);

CLAPACK_API int zgelsy_(integer* m, integer* n, integer* nrhs,
	doublecomplex* a, integer* lda, doublecomplex* b, integer* ldb,
	integer* jpvt, doublereal* rcond, integer* rank, doublecomplex* work,
	integer* lwork, doublereal* rwork, integer* info);

CLAPACK_API int zgeql2_(integer* m, integer* n, doublecomplex* a,
	integer* lda, doublecomplex* tau, doublecomplex* work, integer* info);

CLAPACK_API int zgeqlf_(integer* m, integer* n, doublecomplex* a,
	integer* lda, doublecomplex* tau, doublecomplex* work, integer* lwork,
	integer* info);

CLAPACK_API int zgeqp3_(integer* m, integer* n, doublecomplex* a,
	integer* lda, integer* jpvt, doublecomplex* tau, doublecomplex* work,
	integer* lwork, doublereal* rwork, integer* info);

CLAPACK_API int zgeqpf_(integer* m, integer* n, doublecomplex* a,
	integer* lda, integer* jpvt, doublecomplex* tau, doublecomplex* work,
	doublereal* rwork, integer* info);

CLAPACK_API int zgeqr2_(integer* m, integer* n, doublecomplex* a,
	integer* lda, doublecomplex* tau, doublecomplex* work, integer* info);

CLAPACK_API int zgeqrf_(integer* m, integer* n, doublecomplex* a,
	integer* lda, doublecomplex* tau, doublecomplex* work, integer* lwork,
	integer* info);

CLAPACK_API int zgerfs_(char* trans, integer* n, integer* nrhs,
	doublecomplex* a, integer* lda, doublecomplex* af, integer* ldaf,
	integer* ipiv, doublecomplex* b, integer* ldb, doublecomplex* x,
	integer* ldx, doublereal* ferr, doublereal* berr, doublecomplex* work,
	doublereal* rwork, integer* info);

CLAPACK_API int zgerq2_(integer* m, integer* n, doublecomplex* a,
	integer* lda, doublecomplex* tau, doublecomplex* work, integer* info);

CLAPACK_API int zgerqf_(integer* m, integer* n, doublecomplex* a,
	integer* lda, doublecomplex* tau, doublecomplex* work, integer* lwork,
	integer* info);

CLAPACK_API int zgesc2_(integer* n, doublecomplex* a, integer* lda,
	doublecomplex* rhs, integer* ipiv, integer* jpiv, doublereal* scale);

CLAPACK_API int zgesdd_(char* jobz, integer* m, integer* n,
	doublecomplex* a, integer* lda, doublereal* s, doublecomplex* u,
	integer* ldu, doublecomplex* vt, integer* ldvt, doublecomplex* work,
	integer* lwork, doublereal* rwork, integer* iwork, integer* info);

CLAPACK_API int zgesv_(integer* n, integer* nrhs, doublecomplex* a,
	integer* lda, integer* ipiv, doublecomplex* b, integer* ldb, integer*
	info);

CLAPACK_API int zgesvd_(char* jobu, char* jobvt, integer* m, integer* n,
	doublecomplex* a, integer* lda, doublereal* s, doublecomplex* u,
	integer* ldu, doublecomplex* vt, integer* ldvt, doublecomplex* work,
	integer* lwork, doublereal* rwork, integer* info);

CLAPACK_API int zgesvx_(char* fact, char* trans, integer* n, integer*
	nrhs, doublecomplex* a, integer* lda, doublecomplex* af, integer*
	ldaf, integer* ipiv, char* equed, doublereal* r__, doublereal* c__,
	doublecomplex* b, integer* ldb, doublecomplex* x, integer* ldx,
	doublereal* rcond, doublereal* ferr, doublereal* berr, doublecomplex*
	work, doublereal* rwork, integer* info);

CLAPACK_API int zgetc2_(integer* n, doublecomplex* a, integer* lda,
	integer* ipiv, integer* jpiv, integer* info);

CLAPACK_API int zgetf2_(integer* m, integer* n, doublecomplex* a,
	integer* lda, integer* ipiv, integer* info);

CLAPACK_API int zgetrf_(integer* m, integer* n, doublecomplex* a,
	integer* lda, integer* ipiv, integer* info);

CLAPACK_API int zgetri_(integer* n, doublecomplex* a, integer* lda,
	integer* ipiv, doublecomplex* work, integer* lwork, integer* info);

CLAPACK_API int zgetrs_(char* trans, integer* n, integer* nrhs,
	doublecomplex* a, integer* lda, integer* ipiv, doublecomplex* b,
	integer* ldb, integer* info);

CLAPACK_API int zggbak_(char* job, char* side, integer* n, integer* ilo,
	integer* ihi, doublereal* lscale, doublereal* rscale, integer* m,
	doublecomplex* v, integer* ldv, integer* info);

CLAPACK_API int zggbal_(char* job, integer* n, doublecomplex* a, integer
	* lda, doublecomplex* b, integer* ldb, integer* ilo, integer* ihi,
	doublereal* lscale, doublereal* rscale, doublereal* work, integer*
	info);

CLAPACK_API int zgges_(char* jobvsl, char* jobvsr, char* sort, L_fp
	selctg, integer* n, doublecomplex* a, integer* lda, doublecomplex* b,
	integer* ldb, integer* sdim, doublecomplex* alpha, doublecomplex*
	beta, doublecomplex* vsl, integer* ldvsl, doublecomplex* vsr, integer
	* ldvsr, doublecomplex* work, integer* lwork, doublereal* rwork,
	logical* bwork, integer* info);

CLAPACK_API int zggesx_(char* jobvsl, char* jobvsr, char* sort, L_fp
	selctg, char* sense, integer* n, doublecomplex* a, integer* lda,
	doublecomplex* b, integer* ldb, integer* sdim, doublecomplex* alpha,
	doublecomplex* beta, doublecomplex* vsl, integer* ldvsl,
	doublecomplex* vsr, integer* ldvsr, doublereal* rconde, doublereal*
	rcondv, doublecomplex* work, integer* lwork, doublereal* rwork,
	integer* iwork, integer* liwork, logical* bwork, integer* info);

CLAPACK_API int zggev_(char* jobvl, char* jobvr, integer* n,
	doublecomplex* a, integer* lda, doublecomplex* b, integer* ldb,
	doublecomplex* alpha, doublecomplex* beta, doublecomplex* vl, integer
	* ldvl, doublecomplex* vr, integer* ldvr, doublecomplex* work, integer
	* lwork, doublereal* rwork, integer* info);

CLAPACK_API int zggevx_(char* balanc, char* jobvl, char* jobvr, char*
	sense, integer* n, doublecomplex* a, integer* lda, doublecomplex* b,
	integer* ldb, doublecomplex* alpha, doublecomplex* beta,
	doublecomplex* vl, integer* ldvl, doublecomplex* vr, integer* ldvr,
	integer* ilo, integer* ihi, doublereal* lscale, doublereal* rscale,
	doublereal* abnrm, doublereal* bbnrm, doublereal* rconde, doublereal*
	rcondv, doublecomplex* work, integer* lwork, doublereal* rwork,
	integer* iwork, logical* bwork, integer* info);

CLAPACK_API int zggglm_(integer* n, integer* m, integer* p,
	doublecomplex* a, integer* lda, doublecomplex* b, integer* ldb,
	doublecomplex* d__, doublecomplex* x, doublecomplex* y, doublecomplex
	* work, integer* lwork, integer* info);

CLAPACK_API int zgghrd_(char* compq, char* compz, integer* n, integer*
	ilo, integer* ihi, doublecomplex* a, integer* lda, doublecomplex* b,
	integer* ldb, doublecomplex* q, integer* ldq, doublecomplex* z__,
	integer* ldz, integer* info);

CLAPACK_API int zgglse_(integer* m, integer* n, integer* p,
	doublecomplex* a, integer* lda, doublecomplex* b, integer* ldb,
	doublecomplex* c__, doublecomplex* d__, doublecomplex* x,
	doublecomplex* work, integer* lwork, integer* info);

CLAPACK_API int zggqrf_(integer* n, integer* m, integer* p,
	doublecomplex* a, integer* lda, doublecomplex* taua, doublecomplex* b,
	integer* ldb, doublecomplex* taub, doublecomplex* work, integer*
	lwork, integer* info);

CLAPACK_API int zggrqf_(integer* m, integer* p, integer* n,
	doublecomplex* a, integer* lda, doublecomplex* taua, doublecomplex* b,
	integer* ldb, doublecomplex* taub, doublecomplex* work, integer*
	lwork, integer* info);

CLAPACK_API int zggsvd_(char* jobu, char* jobv, char* jobq, integer* m,
	integer* n, integer* p, integer* k, integer* l, doublecomplex* a,
	integer* lda, doublecomplex* b, integer* ldb, doublereal* alpha,
	doublereal* beta, doublecomplex* u, integer* ldu, doublecomplex* v,
	integer* ldv, doublecomplex* q, integer* ldq, doublecomplex* work,
	doublereal* rwork, integer* iwork, integer* info);

CLAPACK_API int zggsvp_(char* jobu, char* jobv, char* jobq, integer* m,
	integer* p, integer* n, doublecomplex* a, integer* lda, doublecomplex
	* b, integer* ldb, doublereal* tola, doublereal* tolb, integer* k,
	integer* l, doublecomplex* u, integer* ldu, doublecomplex* v, integer
	* ldv, doublecomplex* q, integer* ldq, integer* iwork, doublereal*
	rwork, doublecomplex* tau, doublecomplex* work, integer* info);

CLAPACK_API int zgtcon_(char* norm, integer* n, doublecomplex* dl,
	doublecomplex* d__, doublecomplex* du, doublecomplex* du2, integer*
	ipiv, doublereal* anorm, doublereal* rcond, doublecomplex* work,
	integer* info);

CLAPACK_API int zgtrfs_(char* trans, integer* n, integer* nrhs,
	doublecomplex* dl, doublecomplex* d__, doublecomplex* du,
	doublecomplex* dlf, doublecomplex* df, doublecomplex* duf,
	doublecomplex* du2, integer* ipiv, doublecomplex* b, integer* ldb,
	doublecomplex* x, integer* ldx, doublereal* ferr, doublereal* berr,
	doublecomplex* work, doublereal* rwork, integer* info);

CLAPACK_API int zgtsv_(integer* n, integer* nrhs, doublecomplex* dl,
	doublecomplex* d__, doublecomplex* du, doublecomplex* b, integer* ldb,
	integer* info);

CLAPACK_API int zgtsvx_(char* fact, char* trans, integer* n, integer*
	nrhs, doublecomplex* dl, doublecomplex* d__, doublecomplex* du,
	doublecomplex* dlf, doublecomplex* df, doublecomplex* duf,
	doublecomplex* du2, integer* ipiv, doublecomplex* b, integer* ldb,
	doublecomplex* x, integer* ldx, doublereal* rcond, doublereal* ferr,
	doublereal* berr, doublecomplex* work, doublereal* rwork, integer*
	info);

CLAPACK_API int zgttrf_(integer* n, doublecomplex* dl, doublecomplex*
	d__, doublecomplex* du, doublecomplex* du2, integer* ipiv, integer*
	info);

CLAPACK_API int zgttrs_(char* trans, integer* n, integer* nrhs,
	doublecomplex* dl, doublecomplex* d__, doublecomplex* du,
	doublecomplex* du2, integer* ipiv, doublecomplex* b, integer* ldb,
	integer* info);

CLAPACK_API int zgtts2_(integer* itrans, integer* n, integer* nrhs,
	doublecomplex* dl, doublecomplex* d__, doublecomplex* du,
	doublecomplex* du2, integer* ipiv, doublecomplex* b, integer* ldb);

CLAPACK_API int zhbev_(char* jobz, char* uplo, integer* n, integer* kd,
	doublecomplex* ab, integer* ldab, doublereal* w, doublecomplex* z__,
	integer* ldz, doublecomplex* work, doublereal* rwork, integer* info);

CLAPACK_API int zhbevd_(char* jobz, char* uplo, integer* n, integer* kd,
	doublecomplex* ab, integer* ldab, doublereal* w, doublecomplex* z__,
	integer* ldz, doublecomplex* work, integer* lwork, doublereal* rwork,
	integer* lrwork, integer* iwork, integer* liwork, integer* info);

CLAPACK_API int zhbevx_(char* jobz, char* range, char* uplo, integer* n,
	integer* kd, doublecomplex* ab, integer* ldab, doublecomplex* q,
	integer* ldq, doublereal* vl, doublereal* vu, integer* il, integer*
	iu, doublereal* abstol, integer* m, doublereal* w, doublecomplex* z__,
	integer* ldz, doublecomplex* work, doublereal* rwork, integer* iwork,
	integer* ifail, integer* info);

CLAPACK_API int zhbgst_(char* vect, char* uplo, integer* n, integer* ka,
	integer* kb, doublecomplex* ab, integer* ldab, doublecomplex* bb,
	integer* ldbb, doublecomplex* x, integer* ldx, doublecomplex* work,
	doublereal* rwork, integer* info);

CLAPACK_API int zhbgv_(char* jobz, char* uplo, integer* n, integer* ka,
	integer* kb, doublecomplex* ab, integer* ldab, doublecomplex* bb,
	integer* ldbb, doublereal* w, doublecomplex* z__, integer* ldz,
	doublecomplex* work, doublereal* rwork, integer* info);

CLAPACK_API int zhbgvd_(char* jobz, char* uplo, integer* n, integer* ka,
	integer* kb, doublecomplex* ab, integer* ldab, doublecomplex* bb,
	integer* ldbb, doublereal* w, doublecomplex* z__, integer* ldz,
	doublecomplex* work, integer* lwork, doublereal* rwork, integer*
	lrwork, integer* iwork, integer* liwork, integer* info);

CLAPACK_API int zhbgvx_(char* jobz, char* range, char* uplo, integer* n,
	integer* ka, integer* kb, doublecomplex* ab, integer* ldab,
	doublecomplex* bb, integer* ldbb, doublecomplex* q, integer* ldq,
	doublereal* vl, doublereal* vu, integer* il, integer* iu, doublereal*
	abstol, integer* m, doublereal* w, doublecomplex* z__, integer* ldz,
	doublecomplex* work, doublereal* rwork, integer* iwork, integer*
	ifail, integer* info);

CLAPACK_API int zhbtrd_(char* vect, char* uplo, integer* n, integer* kd,
	doublecomplex* ab, integer* ldab, doublereal* d__, doublereal* e,
	doublecomplex* q, integer* ldq, doublecomplex* work, integer* info);

CLAPACK_API int zhecon_(char* uplo, integer* n, doublecomplex* a,
	integer* lda, integer* ipiv, doublereal* anorm, doublereal* rcond,
	doublecomplex* work, integer* info);

CLAPACK_API int zheequb_(char* uplo, integer* n, doublecomplex* a,
	integer* lda, doublereal* s, doublereal* scond, doublereal* amax,
	doublecomplex* work, integer* info);

CLAPACK_API int zheev_(char* jobz, char* uplo, integer* n, doublecomplex
	* a, integer* lda, doublereal* w, doublecomplex* work, integer* lwork,
	doublereal* rwork, integer* info);

CLAPACK_API int zheevd_(char* jobz, char* uplo, integer* n,
	doublecomplex* a, integer* lda, doublereal* w, doublecomplex* work,
	integer* lwork, doublereal* rwork, integer* lrwork, integer* iwork,
	integer* liwork, integer* info);

CLAPACK_API int zheevr_(char* jobz, char* range, char* uplo, integer* n,
	doublecomplex* a, integer* lda, doublereal* vl, doublereal* vu,
	integer* il, integer* iu, doublereal* abstol, integer* m, doublereal*
	w, doublecomplex* z__, integer* ldz, integer* isuppz, doublecomplex*
	work, integer* lwork, doublereal* rwork, integer* lrwork, integer*
	iwork, integer* liwork, integer* info);

CLAPACK_API int zheevx_(char* jobz, char* range, char* uplo, integer* n,
	doublecomplex* a, integer* lda, doublereal* vl, doublereal* vu,
	integer* il, integer* iu, doublereal* abstol, integer* m, doublereal*
	w, doublecomplex* z__, integer* ldz, doublecomplex* work, integer*
	lwork, doublereal* rwork, integer* iwork, integer* ifail, integer*
	info);

CLAPACK_API int zhegs2_(integer* itype, char* uplo, integer* n,
	doublecomplex* a, integer* lda, doublecomplex* b, integer* ldb,
	integer* info);

CLAPACK_API int zhegst_(integer* itype, char* uplo, integer* n,
	doublecomplex* a, integer* lda, doublecomplex* b, integer* ldb,
	integer* info);

CLAPACK_API int zhegv_(integer* itype, char* jobz, char* uplo, integer*
	n, doublecomplex* a, integer* lda, doublecomplex* b, integer* ldb,
	doublereal* w, doublecomplex* work, integer* lwork, doublereal* rwork,
	integer* info);

CLAPACK_API int zhegvd_(integer* itype, char* jobz, char* uplo, integer*
	n, doublecomplex* a, integer* lda, doublecomplex* b, integer* ldb,
	doublereal* w, doublecomplex* work, integer* lwork, doublereal* rwork,
	integer* lrwork, integer* iwork, integer* liwork, integer* info);

CLAPACK_API int zhegvx_(integer* itype, char* jobz, char* range, char*
	uplo, integer* n, doublecomplex* a, integer* lda, doublecomplex* b,
	integer* ldb, doublereal* vl, doublereal* vu, integer* il, integer*
	iu, doublereal* abstol, integer* m, doublereal* w, doublecomplex* z__,
	integer* ldz, doublecomplex* work, integer* lwork, doublereal* rwork,
	integer* iwork, integer* ifail, integer* info);

CLAPACK_API int zherfs_(char* uplo, integer* n, integer* nrhs,
	doublecomplex* a, integer* lda, doublecomplex* af, integer* ldaf,
	integer* ipiv, doublecomplex* b, integer* ldb, doublecomplex* x,
	integer* ldx, doublereal* ferr, doublereal* berr, doublecomplex* work,
	doublereal* rwork, integer* info);

CLAPACK_API int zhesv_(char* uplo, integer* n, integer* nrhs,
	doublecomplex* a, integer* lda, integer* ipiv, doublecomplex* b,
	integer* ldb, doublecomplex* work, integer* lwork, integer* info);

CLAPACK_API int zhesvx_(char* fact, char* uplo, integer* n, integer*
	nrhs, doublecomplex* a, integer* lda, doublecomplex* af, integer*
	ldaf, integer* ipiv, doublecomplex* b, integer* ldb, doublecomplex* x,
	integer* ldx, doublereal* rcond, doublereal* ferr, doublereal* berr,
	doublecomplex* work, integer* lwork, doublereal* rwork, integer* info);

CLAPACK_API int zhetd2_(char* uplo, integer* n, doublecomplex* a,
	integer* lda, doublereal* d__, doublereal* e, doublecomplex* tau,
	integer* info);

CLAPACK_API int zhetf2_(char* uplo, integer* n, doublecomplex* a,
	integer* lda, integer* ipiv, integer* info);

CLAPACK_API int zhetrd_(char* uplo, integer* n, doublecomplex* a,
	integer* lda, doublereal* d__, doublereal* e, doublecomplex* tau,
	doublecomplex* work, integer* lwork, integer* info);

CLAPACK_API int zhetrf_(char* uplo, integer* n, doublecomplex* a,
	integer* lda, integer* ipiv, doublecomplex* work, integer* lwork,
	integer* info);

CLAPACK_API int zhetri_(char* uplo, integer* n, doublecomplex* a,
	integer* lda, integer* ipiv, doublecomplex* work, integer* info);

CLAPACK_API int zhetrs_(char* uplo, integer* n, integer* nrhs,
	doublecomplex* a, integer* lda, integer* ipiv, doublecomplex* b,
	integer* ldb, integer* info);

CLAPACK_API int zhfrk_(char* transr, char* uplo, char* trans, integer* n,
	integer* k, doublereal* alpha, doublecomplex* a, integer* lda,
	doublereal* beta, doublecomplex* c__);

CLAPACK_API int zhgeqz_(char* job, char* compq, char* compz, integer* n,
	integer* ilo, integer* ihi, doublecomplex* h__, integer* ldh,
	doublecomplex* t, integer* ldt, doublecomplex* alpha, doublecomplex*
	beta, doublecomplex* q, integer* ldq, doublecomplex* z__, integer*
	ldz, doublecomplex* work, integer* lwork, doublereal* rwork, integer*
	info);

CLAPACK_API int zhpcon_(char* uplo, integer* n, doublecomplex* ap,
	integer* ipiv, doublereal* anorm, doublereal* rcond, doublecomplex*
	work, integer* info);

CLAPACK_API int zhpev_(char* jobz, char* uplo, integer* n, doublecomplex
	* ap, doublereal* w, doublecomplex* z__, integer* ldz, doublecomplex*
	work, doublereal* rwork, integer* info);

CLAPACK_API int zhpevd_(char* jobz, char* uplo, integer* n,
	doublecomplex* ap, doublereal* w, doublecomplex* z__, integer* ldz,
	doublecomplex* work, integer* lwork, doublereal* rwork, integer*
	lrwork, integer* iwork, integer* liwork, integer* info);

CLAPACK_API int zhpevx_(char* jobz, char* range, char* uplo, integer* n,
	doublecomplex* ap, doublereal* vl, doublereal* vu, integer* il,
	integer* iu, doublereal* abstol, integer* m, doublereal* w,
	doublecomplex* z__, integer* ldz, doublecomplex* work, doublereal*
	rwork, integer* iwork, integer* ifail, integer* info);

CLAPACK_API int zhpgst_(integer* itype, char* uplo, integer* n,
	doublecomplex* ap, doublecomplex* bp, integer* info);

CLAPACK_API int zhpgv_(integer* itype, char* jobz, char* uplo, integer*
	n, doublecomplex* ap, doublecomplex* bp, doublereal* w, doublecomplex
	* z__, integer* ldz, doublecomplex* work, doublereal* rwork, integer*
	info);

CLAPACK_API int zhpgvd_(integer* itype, char* jobz, char* uplo, integer*
	n, doublecomplex* ap, doublecomplex* bp, doublereal* w, doublecomplex
	* z__, integer* ldz, doublecomplex* work, integer* lwork, doublereal*
	rwork, integer* lrwork, integer* iwork, integer* liwork, integer*
	info);

CLAPACK_API int zhpgvx_(integer* itype, char* jobz, char* range, char*
	uplo, integer* n, doublecomplex* ap, doublecomplex* bp, doublereal*
	vl, doublereal* vu, integer* il, integer* iu, doublereal* abstol,
	integer* m, doublereal* w, doublecomplex* z__, integer* ldz,
	doublecomplex* work, doublereal* rwork, integer* iwork, integer*
	ifail, integer* info);

CLAPACK_API int zhprfs_(char* uplo, integer* n, integer* nrhs,
	doublecomplex* ap, doublecomplex* afp, integer* ipiv, doublecomplex*
	b, integer* ldb, doublecomplex* x, integer* ldx, doublereal* ferr,
	doublereal* berr, doublecomplex* work, doublereal* rwork, integer*
	info);

CLAPACK_API int zhpsv_(char* uplo, integer* n, integer* nrhs,
	doublecomplex* ap, integer* ipiv, doublecomplex* b, integer* ldb,
	integer* info);

CLAPACK_API int zhpsvx_(char* fact, char* uplo, integer* n, integer*
	nrhs, doublecomplex* ap, doublecomplex* afp, integer* ipiv,
	doublecomplex* b, integer* ldb, doublecomplex* x, integer* ldx,
	doublereal* rcond, doublereal* ferr, doublereal* berr, doublecomplex*
	work, doublereal* rwork, integer* info);

CLAPACK_API int zhptrd_(char* uplo, integer* n, doublecomplex* ap,
	doublereal* d__, doublereal* e, doublecomplex* tau, integer* info);

CLAPACK_API int zhptrf_(char* uplo, integer* n, doublecomplex* ap,
	integer* ipiv, integer* info);

CLAPACK_API int zhptri_(char* uplo, integer* n, doublecomplex* ap,
	integer* ipiv, doublecomplex* work, integer* info);

CLAPACK_API int zhptrs_(char* uplo, integer* n, integer* nrhs,
	doublecomplex* ap, integer* ipiv, doublecomplex* b, integer* ldb,
	integer* info);

CLAPACK_API int zhsein_(char* side, char* eigsrc, char* initv, logical*
	select, integer* n, doublecomplex* h__, integer* ldh, doublecomplex*
	w, doublecomplex* vl, integer* ldvl, doublecomplex* vr, integer* ldvr,
	integer* mm, integer* m, doublecomplex* work, doublereal* rwork,
	integer* ifaill, integer* ifailr, integer* info);

CLAPACK_API int zhseqr_(char* job, char* compz, integer* n, integer* ilo,
	integer* ihi, doublecomplex* h__, integer* ldh, doublecomplex* w,
	doublecomplex* z__, integer* ldz, doublecomplex* work, integer* lwork,
	integer* info);

CLAPACK_API int zlabrd_(integer* m, integer* n, integer* nb,
	doublecomplex* a, integer* lda, doublereal* d__, doublereal* e,
	doublecomplex* tauq, doublecomplex* taup, doublecomplex* x, integer*
	ldx, doublecomplex* y, integer* ldy);

CLAPACK_API int zlacgv_(integer* n, doublecomplex* x, integer* incx);

CLAPACK_API int zlacn2_(integer* n, doublecomplex* v, doublecomplex* x,
	doublereal* est, integer* kase, integer* isave);

CLAPACK_API int zlacon_(integer* n, doublecomplex* v, doublecomplex* x,
	doublereal* est, integer* kase);

CLAPACK_API int zlacp2_(char* uplo, integer* m, integer* n, doublereal*
	a, integer* lda, doublecomplex* b, integer* ldb);

CLAPACK_API int zlacpy_(char* uplo, integer* m, integer* n,
	doublecomplex* a, integer* lda, doublecomplex* b, integer* ldb);

CLAPACK_API int zlacrm_(integer* m, integer* n, doublecomplex* a,
	integer* lda, doublereal* b, integer* ldb, doublecomplex* c__,
	integer* ldc, doublereal* rwork);

CLAPACK_API int zlacrt_(integer* n, doublecomplex* cx, integer* incx,
	doublecomplex* cy, integer* incy, doublecomplex* c__, doublecomplex*
	s);

/* Double Complex */ VOID zladiv_(doublecomplex* ret_val, doublecomplex* x,
	doublecomplex* y);

CLAPACK_API int zlaed0_(integer* qsiz, integer* n, doublereal* d__,
	doublereal* e, doublecomplex* q, integer* ldq, doublecomplex* qstore,
	integer* ldqs, doublereal* rwork, integer* iwork, integer* info);

CLAPACK_API int zlaed7_(integer* n, integer* cutpnt, integer* qsiz,
	integer* tlvls, integer* curlvl, integer* curpbm, doublereal* d__,
	doublecomplex* q, integer* ldq, doublereal* rho, integer* indxq,
	doublereal* qstore, integer* qptr, integer* prmptr, integer* perm,
	integer* givptr, integer* givcol, doublereal* givnum, doublecomplex*
	work, doublereal* rwork, integer* iwork, integer* info);

CLAPACK_API int zlaed8_(integer* k, integer* n, integer* qsiz,
	doublecomplex* q, integer* ldq, doublereal* d__, doublereal* rho,
	integer* cutpnt, doublereal* z__, doublereal* dlamda, doublecomplex*
	q2, integer* ldq2, doublereal* w, integer* indxp, integer* indx,
	integer* indxq, integer* perm, integer* givptr, integer* givcol,
	doublereal* givnum, integer* info);

CLAPACK_API int zlaein_(logical* rightv, logical* noinit, integer* n,
	doublecomplex* h__, integer* ldh, doublecomplex* w, doublecomplex* v,
	doublecomplex* b, integer* ldb, doublereal* rwork, doublereal* eps3,
	doublereal* smlnum, integer* info);

CLAPACK_API int zlaesy_(doublecomplex* a, doublecomplex* b,
	doublecomplex* c__, doublecomplex* rt1, doublecomplex* rt2,
	doublecomplex* evscal, doublecomplex* cs1, doublecomplex* sn1);

CLAPACK_API int zlaev2_(doublecomplex* a, doublecomplex* b,
	doublecomplex* c__, doublereal* rt1, doublereal* rt2, doublereal* cs1,
	doublecomplex* sn1);

CLAPACK_API int zlag2c_(integer* m, integer* n, doublecomplex* a,
	integer* lda, complex* sa, integer* ldsa, integer* info);

CLAPACK_API int zlags2_(logical* upper, doublereal* a1, doublecomplex*
	a2, doublereal* a3, doublereal* b1, doublecomplex* b2, doublereal* b3,
	doublereal* csu, doublecomplex* snu, doublereal* csv, doublecomplex*
	snv, doublereal* csq, doublecomplex* snq);

CLAPACK_API int zlagtm_(char* trans, integer* n, integer* nrhs,
	doublereal* alpha, doublecomplex* dl, doublecomplex* d__,
	doublecomplex* du, doublecomplex* x, integer* ldx, doublereal* beta,
	doublecomplex* b, integer* ldb);

CLAPACK_API int zlahef_(char* uplo, integer* n, integer* nb, integer* kb,
	doublecomplex* a, integer* lda, integer* ipiv, doublecomplex* w,
	integer* ldw, integer* info);

CLAPACK_API int zlahqr_(logical* wantt, logical* wantz, integer* n,
	integer* ilo, integer* ihi, doublecomplex* h__, integer* ldh,
	doublecomplex* w, integer* iloz, integer* ihiz, doublecomplex* z__,
	integer* ldz, integer* info);

CLAPACK_API int zlahr2_(integer* n, integer* k, integer* nb,
	doublecomplex* a, integer* lda, doublecomplex* tau, doublecomplex* t,
	integer* ldt, doublecomplex* y, integer* ldy);

CLAPACK_API int zlahrd_(integer* n, integer* k, integer* nb,
	doublecomplex* a, integer* lda, doublecomplex* tau, doublecomplex* t,
	integer* ldt, doublecomplex* y, integer* ldy);

CLAPACK_API int zlaic1_(integer* job, integer* j, doublecomplex* x,
	doublereal* sest, doublecomplex* w, doublecomplex* gamma, doublereal*
	sestpr, doublecomplex* s, doublecomplex* c__);

CLAPACK_API int zlals0_(integer* icompq, integer* nl, integer* nr,
	integer* sqre, integer* nrhs, doublecomplex* b, integer* ldb,
	doublecomplex* bx, integer* ldbx, integer* perm, integer* givptr,
	integer* givcol, integer* ldgcol, doublereal* givnum, integer* ldgnum,
	doublereal* poles, doublereal* difl, doublereal* difr, doublereal*
	z__, integer* k, doublereal* c__, doublereal* s, doublereal* rwork,
	integer* info);

CLAPACK_API int zlalsa_(integer* icompq, integer* smlsiz, integer* n,
	integer* nrhs, doublecomplex* b, integer* ldb, doublecomplex* bx,
	integer* ldbx, doublereal* u, integer* ldu, doublereal* vt, integer*
	k, doublereal* difl, doublereal* difr, doublereal* z__, doublereal*
	poles, integer* givptr, integer* givcol, integer* ldgcol, integer*
	perm, doublereal* givnum, doublereal* c__, doublereal* s, doublereal*
	rwork, integer* iwork, integer* info);

CLAPACK_API int zlalsd_(char* uplo, integer* smlsiz, integer* n, integer
	* nrhs, doublereal* d__, doublereal* e, doublecomplex* b, integer* ldb,
	doublereal* rcond, integer* rank, doublecomplex* work, doublereal*
	rwork, integer* iwork, integer* info);

CLAPACK_API
doublereal zlangb_(char* norm, integer* n, integer* kl, integer* ku,
	doublecomplex* ab, integer* ldab, doublereal* work);

CLAPACK_API
doublereal zlange_(char* norm, integer* m, integer* n, doublecomplex* a,
	integer* lda, doublereal* work);

CLAPACK_API
doublereal zlangt_(char* norm, integer* n, doublecomplex* dl, doublecomplex*
	d__, doublecomplex* du);

CLAPACK_API
doublereal zlanhb_(char* norm, char* uplo, integer* n, integer* k,
	doublecomplex* ab, integer* ldab, doublereal* work);

CLAPACK_API
doublereal zlanhe_(char* norm, char* uplo, integer* n, doublecomplex* a,
	integer* lda, doublereal* work);

CLAPACK_API
doublereal zlanhf_(char* norm, char* transr, char* uplo, integer* n,
	doublecomplex* a, doublereal* work);

CLAPACK_API
doublereal zlanhp_(char* norm, char* uplo, integer* n, doublecomplex* ap,
	doublereal* work);

CLAPACK_API
doublereal zlanhs_(char* norm, integer* n, doublecomplex* a, integer* lda,
	doublereal* work);

CLAPACK_API
doublereal zlanht_(char* norm, integer* n, doublereal* d__, doublecomplex* e);

CLAPACK_API
doublereal zlansb_(char* norm, char* uplo, integer* n, integer* k,
	doublecomplex* ab, integer* ldab, doublereal* work);

CLAPACK_API
doublereal zlansp_(char* norm, char* uplo, integer* n, doublecomplex* ap,
	doublereal* work);

CLAPACK_API
doublereal zlansy_(char* norm, char* uplo, integer* n, doublecomplex* a,
	integer* lda, doublereal* work);

CLAPACK_API
doublereal zlantb_(char* norm, char* uplo, char* diag, integer* n, integer* k,
	doublecomplex* ab, integer* ldab, doublereal* work);

CLAPACK_API
doublereal zlantp_(char* norm, char* uplo, char* diag, integer* n,
	doublecomplex* ap, doublereal* work);

CLAPACK_API
doublereal zlantr_(char* norm, char* uplo, char* diag, integer* m, integer* n,
	doublecomplex* a, integer* lda, doublereal* work);

CLAPACK_API int zlapll_(integer* n, doublecomplex* x, integer* incx,
	doublecomplex* y, integer* incy, doublereal* ssmin);

CLAPACK_API int zlapmt_(logical* forwrd, integer* m, integer* n,
	doublecomplex* x, integer* ldx, integer* k);

CLAPACK_API int zlaqgb_(integer* m, integer* n, integer* kl, integer* ku,
	doublecomplex* ab, integer* ldab, doublereal* r__, doublereal* c__,
	doublereal* rowcnd, doublereal* colcnd, doublereal* amax, char* equed);

CLAPACK_API int zlaqge_(integer* m, integer* n, doublecomplex* a,
	integer* lda, doublereal* r__, doublereal* c__, doublereal* rowcnd,
	doublereal* colcnd, doublereal* amax, char* equed);

CLAPACK_API int zlaqhb_(char* uplo, integer* n, integer* kd,
	doublecomplex* ab, integer* ldab, doublereal* s, doublereal* scond,
	doublereal* amax, char* equed);

CLAPACK_API int zlaqhe_(char* uplo, integer* n, doublecomplex* a,
	integer* lda, doublereal* s, doublereal* scond, doublereal* amax,
	char* equed);

CLAPACK_API int zlaqhp_(char* uplo, integer* n, doublecomplex* ap,
	doublereal* s, doublereal* scond, doublereal* amax, char* equed);

CLAPACK_API int zlaqp2_(integer* m, integer* n, integer* offset,
	doublecomplex* a, integer* lda, integer* jpvt, doublecomplex* tau,
	doublereal* vn1, doublereal* vn2, doublecomplex* work);

CLAPACK_API int zlaqps_(integer* m, integer* n, integer* offset, integer
	* nb, integer* kb, doublecomplex* a, integer* lda, integer* jpvt,
	doublecomplex* tau, doublereal* vn1, doublereal* vn2, doublecomplex*
	auxv, doublecomplex* f, integer* ldf);

CLAPACK_API int zlaqr0_(logical* wantt, logical* wantz, integer* n,
	integer* ilo, integer* ihi, doublecomplex* h__, integer* ldh,
	doublecomplex* w, integer* iloz, integer* ihiz, doublecomplex* z__,
	integer* ldz, doublecomplex* work, integer* lwork, integer* info);

CLAPACK_API int zlaqr1_(integer* n, doublecomplex* h__, integer* ldh,
	doublecomplex* s1, doublecomplex* s2, doublecomplex* v);

CLAPACK_API int zlaqr2_(logical* wantt, logical* wantz, integer* n,
	integer* ktop, integer* kbot, integer* nw, doublecomplex* h__,
	integer* ldh, integer* iloz, integer* ihiz, doublecomplex* z__,
	integer* ldz, integer* ns, integer* nd, doublecomplex* sh,
	doublecomplex* v, integer* ldv, integer* nh, doublecomplex* t,
	integer* ldt, integer* nv, doublecomplex* wv, integer* ldwv,
	doublecomplex* work, integer* lwork);

CLAPACK_API int zlaqr3_(logical* wantt, logical* wantz, integer* n,
	integer* ktop, integer* kbot, integer* nw, doublecomplex* h__,
	integer* ldh, integer* iloz, integer* ihiz, doublecomplex* z__,
	integer* ldz, integer* ns, integer* nd, doublecomplex* sh,
	doublecomplex* v, integer* ldv, integer* nh, doublecomplex* t,
	integer* ldt, integer* nv, doublecomplex* wv, integer* ldwv,
	doublecomplex* work, integer* lwork);

CLAPACK_API int zlaqr4_(logical* wantt, logical* wantz, integer* n,
	integer* ilo, integer* ihi, doublecomplex* h__, integer* ldh,
	doublecomplex* w, integer* iloz, integer* ihiz, doublecomplex* z__,
	integer* ldz, doublecomplex* work, integer* lwork, integer* info);

CLAPACK_API int zlaqr5_(logical* wantt, logical* wantz, integer* kacc22,
	integer* n, integer* ktop, integer* kbot, integer* nshfts,
	doublecomplex* s, doublecomplex* h__, integer* ldh, integer* iloz,
	integer* ihiz, doublecomplex* z__, integer* ldz, doublecomplex* v,
	integer* ldv, doublecomplex* u, integer* ldu, integer* nv,
	doublecomplex* wv, integer* ldwv, integer* nh, doublecomplex* wh,
	integer* ldwh);

CLAPACK_API int zlaqsb_(char* uplo, integer* n, integer* kd,
	doublecomplex* ab, integer* ldab, doublereal* s, doublereal* scond,
	doublereal* amax, char* equed);

CLAPACK_API int zlaqsp_(char* uplo, integer* n, doublecomplex* ap,
	doublereal* s, doublereal* scond, doublereal* amax, char* equed);

CLAPACK_API int zlaqsy_(char* uplo, integer* n, doublecomplex* a,
	integer* lda, doublereal* s, doublereal* scond, doublereal* amax,
	char* equed);

CLAPACK_API int zlar1v_(integer* n, integer* b1, integer* bn, doublereal
	* lambda, doublereal* d__, doublereal* l, doublereal* ld, doublereal*
	lld, doublereal* pivmin, doublereal* gaptol, doublecomplex* z__,
	logical* wantnc, integer* negcnt, doublereal* ztz, doublereal* mingma,
	integer* r__, integer* isuppz, doublereal* nrminv, doublereal* resid,
	doublereal* rqcorr, doublereal* work);

CLAPACK_API int zlar2v_(integer* n, doublecomplex* x, doublecomplex* y,
	doublecomplex* z__, integer* incx, doublereal* c__, doublecomplex* s,
	integer* incc);

CLAPACK_API int zlarcm_(integer* m, integer* n, doublereal* a, integer*
	lda, doublecomplex* b, integer* ldb, doublecomplex* c__, integer* ldc,
	doublereal* rwork);

CLAPACK_API int zlarf_(char* side, integer* m, integer* n, doublecomplex
	* v, integer* incv, doublecomplex* tau, doublecomplex* c__, integer*
	ldc, doublecomplex* work);

CLAPACK_API int zlarfb_(char* side, char* trans, char* direct, char*
	storev, integer* m, integer* n, integer* k, doublecomplex* v, integer
	* ldv, doublecomplex* t, integer* ldt, doublecomplex* c__, integer*
	ldc, doublecomplex* work, integer* ldwork);

CLAPACK_API int zlarfg_(integer* n, doublecomplex* alpha, doublecomplex*
	x, integer* incx, doublecomplex* tau);

CLAPACK_API int zlarfp_(integer* n, doublecomplex* alpha, doublecomplex*
	x, integer* incx, doublecomplex* tau);

CLAPACK_API int zlarft_(char* direct, char* storev, integer* n, integer*
	k, doublecomplex* v, integer* ldv, doublecomplex* tau, doublecomplex*
	t, integer* ldt);

CLAPACK_API int zlarfx_(char* side, integer* m, integer* n,
	doublecomplex* v, doublecomplex* tau, doublecomplex* c__, integer*
	ldc, doublecomplex* work);

CLAPACK_API int zlargv_(integer* n, doublecomplex* x, integer* incx,
	doublecomplex* y, integer* incy, doublereal* c__, integer* incc);

CLAPACK_API int zlarnv_(integer* idist, integer* iseed, integer* n,
	doublecomplex* x);

CLAPACK_API int zlarrv_(integer* n, doublereal* vl, doublereal* vu,
	doublereal* d__, doublereal* l, doublereal* pivmin, integer* isplit,
	integer* m, integer* dol, integer* dou, doublereal* minrgp,
	doublereal* rtol1, doublereal* rtol2, doublereal* w, doublereal* werr,
	doublereal* wgap, integer* iblock, integer* indexw, doublereal* gers,
	doublecomplex* z__, integer* ldz, integer* isuppz, doublereal* work,
	integer* iwork, integer* info);

CLAPACK_API int zlartg_(doublecomplex* f, doublecomplex* g, doublereal*
	cs, doublecomplex* sn, doublecomplex* r__);

CLAPACK_API int zlartv_(integer* n, doublecomplex* x, integer* incx,
	doublecomplex* y, integer* incy, doublereal* c__, doublecomplex* s,
	integer* incc);

CLAPACK_API int zlarz_(char* side, integer* m, integer* n, integer* l,
	doublecomplex* v, integer* incv, doublecomplex* tau, doublecomplex*
	c__, integer* ldc, doublecomplex* work);

CLAPACK_API int zlarzb_(char* side, char* trans, char* direct, char*
	storev, integer* m, integer* n, integer* k, integer* l, doublecomplex
	* v, integer* ldv, doublecomplex* t, integer* ldt, doublecomplex* c__,
	integer* ldc, doublecomplex* work, integer* ldwork);

CLAPACK_API int zlarzt_(char* direct, char* storev, integer* n, integer*
	k, doublecomplex* v, integer* ldv, doublecomplex* tau, doublecomplex*
	t, integer* ldt);

CLAPACK_API int zlascl_(char* type__, integer* kl, integer* ku,
	doublereal* cfrom, doublereal* cto, integer* m, integer* n,
	doublecomplex* a, integer* lda, integer* info);

CLAPACK_API int zlascl2_(integer* m, integer* n, doublereal* d__,
	doublecomplex* x, integer* ldx);

CLAPACK_API int zlaset_(char* uplo, integer* m, integer* n,
	doublecomplex* alpha, doublecomplex* beta, doublecomplex* a, integer*
	lda);

CLAPACK_API int zlasr_(char* side, char* pivot, char* direct, integer* m,
	integer* n, doublereal* c__, doublereal* s, doublecomplex* a,
	integer* lda);

CLAPACK_API int zlassq_(integer* n, doublecomplex* x, integer* incx,
	doublereal* scale, doublereal* sumsq);

CLAPACK_API int zlaswp_(integer* n, doublecomplex* a, integer* lda,
	integer* k1, integer* k2, integer* ipiv, integer* incx);

CLAPACK_API int zlasyf_(char* uplo, integer* n, integer* nb, integer* kb,
	doublecomplex* a, integer* lda, integer* ipiv, doublecomplex* w,
	integer* ldw, integer* info);

CLAPACK_API int zlat2c_(char* uplo, integer* n, doublecomplex* a,
	integer* lda, complex* sa, integer* ldsa, integer* info);

CLAPACK_API int zlatbs_(char* uplo, char* trans, char* diag, char*
	normin, integer* n, integer* kd, doublecomplex* ab, integer* ldab,
	doublecomplex* x, doublereal* scale, doublereal* cnorm, integer* info);

CLAPACK_API int zlatdf_(integer* ijob, integer* n, doublecomplex* z__,
	integer* ldz, doublecomplex* rhs, doublereal* rdsum, doublereal*
	rdscal, integer* ipiv, integer* jpiv);

CLAPACK_API int zlatps_(char* uplo, char* trans, char* diag, char*
	normin, integer* n, doublecomplex* ap, doublecomplex* x, doublereal*
	scale, doublereal* cnorm, integer* info);

CLAPACK_API int zlatrd_(char* uplo, integer* n, integer* nb,
	doublecomplex* a, integer* lda, doublereal* e, doublecomplex* tau,
	doublecomplex* w, integer* ldw);

CLAPACK_API int zlatrs_(char* uplo, char* trans, char* diag, char*
	normin, integer* n, doublecomplex* a, integer* lda, doublecomplex* x,
	doublereal* scale, doublereal* cnorm, integer* info);

CLAPACK_API int zlatrz_(integer* m, integer* n, integer* l,
	doublecomplex* a, integer* lda, doublecomplex* tau, doublecomplex*
	work);

CLAPACK_API int zlatzm_(char* side, integer* m, integer* n,
	doublecomplex* v, integer* incv, doublecomplex* tau, doublecomplex*
	c1, doublecomplex* c2, integer* ldc, doublecomplex* work);

CLAPACK_API int zlauu2_(char* uplo, integer* n, doublecomplex* a,
	integer* lda, integer* info);

CLAPACK_API int zlauum_(char* uplo, integer* n, doublecomplex* a,
	integer* lda, integer* info);

CLAPACK_API int zpbcon_(char* uplo, integer* n, integer* kd,
	doublecomplex* ab, integer* ldab, doublereal* anorm, doublereal*
	rcond, doublecomplex* work, doublereal* rwork, integer* info);

CLAPACK_API int zpbequ_(char* uplo, integer* n, integer* kd,
	doublecomplex* ab, integer* ldab, doublereal* s, doublereal* scond,
	doublereal* amax, integer* info);

CLAPACK_API int zpbrfs_(char* uplo, integer* n, integer* kd, integer*
	nrhs, doublecomplex* ab, integer* ldab, doublecomplex* afb, integer*
	ldafb, doublecomplex* b, integer* ldb, doublecomplex* x, integer* ldx,
	doublereal* ferr, doublereal* berr, doublecomplex* work, doublereal*
	rwork, integer* info);

CLAPACK_API int zpbstf_(char* uplo, integer* n, integer* kd,
	doublecomplex* ab, integer* ldab, integer* info);

CLAPACK_API int zpbsv_(char* uplo, integer* n, integer* kd, integer*
	nrhs, doublecomplex* ab, integer* ldab, doublecomplex* b, integer*
	ldb, integer* info);

CLAPACK_API int zpbsvx_(char* fact, char* uplo, integer* n, integer* kd,
	integer* nrhs, doublecomplex* ab, integer* ldab, doublecomplex* afb,
	integer* ldafb, char* equed, doublereal* s, doublecomplex* b, integer
	* ldb, doublecomplex* x, integer* ldx, doublereal* rcond, doublereal*
	ferr, doublereal* berr, doublecomplex* work, doublereal* rwork,
	integer* info);

CLAPACK_API int zpbtf2_(char* uplo, integer* n, integer* kd,
	doublecomplex* ab, integer* ldab, integer* info);

CLAPACK_API int zpbtrf_(char* uplo, integer* n, integer* kd,
	doublecomplex* ab, integer* ldab, integer* info);

CLAPACK_API int zpbtrs_(char* uplo, integer* n, integer* kd, integer*
	nrhs, doublecomplex* ab, integer* ldab, doublecomplex* b, integer*
	ldb, integer* info);

CLAPACK_API int zpftrf_(char* transr, char* uplo, integer* n,
	doublecomplex* a, integer* info);

CLAPACK_API int zpftri_(char* transr, char* uplo, integer* n,
	doublecomplex* a, integer* info);

CLAPACK_API int zpftrs_(char* transr, char* uplo, integer* n, integer*
	nrhs, doublecomplex* a, doublecomplex* b, integer* ldb, integer* info);

CLAPACK_API int zpocon_(char* uplo, integer* n, doublecomplex* a,
	integer* lda, doublereal* anorm, doublereal* rcond, doublecomplex*
	work, doublereal* rwork, integer* info);

CLAPACK_API int zpoequ_(integer* n, doublecomplex* a, integer* lda,
	doublereal* s, doublereal* scond, doublereal* amax, integer* info);

CLAPACK_API int zpoequb_(integer* n, doublecomplex* a, integer* lda,
	doublereal* s, doublereal* scond, doublereal* amax, integer* info);

CLAPACK_API int zporfs_(char* uplo, integer* n, integer* nrhs,
	doublecomplex* a, integer* lda, doublecomplex* af, integer* ldaf,
	doublecomplex* b, integer* ldb, doublecomplex* x, integer* ldx,
	doublereal* ferr, doublereal* berr, doublecomplex* work, doublereal*
	rwork, integer* info);

CLAPACK_API int zposv_(char* uplo, integer* n, integer* nrhs,
	doublecomplex* a, integer* lda, doublecomplex* b, integer* ldb,
	integer* info);

CLAPACK_API int zposvx_(char* fact, char* uplo, integer* n, integer*
	nrhs, doublecomplex* a, integer* lda, doublecomplex* af, integer*
	ldaf, char* equed, doublereal* s, doublecomplex* b, integer* ldb,
	doublecomplex* x, integer* ldx, doublereal* rcond, doublereal* ferr,
	doublereal* berr, doublecomplex* work, doublereal* rwork, integer*
	info);

CLAPACK_API int zpotf2_(char* uplo, integer* n, doublecomplex* a,
	integer* lda, integer* info);

CLAPACK_API int zpotrf_(char* uplo, integer* n, doublecomplex* a,
	integer* lda, integer* info);

CLAPACK_API int zpotri_(char* uplo, integer* n, doublecomplex* a,
	integer* lda, integer* info);

CLAPACK_API int zpotrs_(char* uplo, integer* n, integer* nrhs,
	doublecomplex* a, integer* lda, doublecomplex* b, integer* ldb,
	integer* info);

CLAPACK_API int zppcon_(char* uplo, integer* n, doublecomplex* ap,
	doublereal* anorm, doublereal* rcond, doublecomplex* work, doublereal
	* rwork, integer* info);

CLAPACK_API int zppequ_(char* uplo, integer* n, doublecomplex* ap,
	doublereal* s, doublereal* scond, doublereal* amax, integer* info);

CLAPACK_API int zpprfs_(char* uplo, integer* n, integer* nrhs,
	doublecomplex* ap, doublecomplex* afp, doublecomplex* b, integer* ldb,
	doublecomplex* x, integer* ldx, doublereal* ferr, doublereal* berr,
	doublecomplex* work, doublereal* rwork, integer* info);

CLAPACK_API int zppsv_(char* uplo, integer* n, integer* nrhs,
	doublecomplex* ap, doublecomplex* b, integer* ldb, integer* info);

CLAPACK_API int zppsvx_(char* fact, char* uplo, integer* n, integer*
	nrhs, doublecomplex* ap, doublecomplex* afp, char* equed, doublereal*
	s, doublecomplex* b, integer* ldb, doublecomplex* x, integer* ldx,
	doublereal* rcond, doublereal* ferr, doublereal* berr, doublecomplex*
	work, doublereal* rwork, integer* info);

CLAPACK_API int zpptrf_(char* uplo, integer* n, doublecomplex* ap,
	integer* info);

CLAPACK_API int zpptri_(char* uplo, integer* n, doublecomplex* ap,
	integer* info);

CLAPACK_API int zpptrs_(char* uplo, integer* n, integer* nrhs,
	doublecomplex* ap, doublecomplex* b, integer* ldb, integer* info);

CLAPACK_API int zpstf2_(char* uplo, integer* n, doublecomplex* a,
	integer* lda, integer* piv, integer* rank, doublereal* tol,
	doublereal* work, integer* info);

CLAPACK_API int zpstrf_(char* uplo, integer* n, doublecomplex* a,
	integer* lda, integer* piv, integer* rank, doublereal* tol,
	doublereal* work, integer* info);

CLAPACK_API int zptcon_(integer* n, doublereal* d__, doublecomplex* e,
	doublereal* anorm, doublereal* rcond, doublereal* rwork, integer*
	info);

CLAPACK_API int zpteqr_(char* compz, integer* n, doublereal* d__,
	doublereal* e, doublecomplex* z__, integer* ldz, doublereal* work,
	integer* info);

CLAPACK_API int zptrfs_(char* uplo, integer* n, integer* nrhs,
	doublereal* d__, doublecomplex* e, doublereal* df, doublecomplex* ef,
	doublecomplex* b, integer* ldb, doublecomplex* x, integer* ldx,
	doublereal* ferr, doublereal* berr, doublecomplex* work, doublereal*
	rwork, integer* info);

CLAPACK_API int zptsv_(integer* n, integer* nrhs, doublereal* d__,
	doublecomplex* e, doublecomplex* b, integer* ldb, integer* info);

CLAPACK_API int zptsvx_(char* fact, integer* n, integer* nrhs,
	doublereal* d__, doublecomplex* e, doublereal* df, doublecomplex* ef,
	doublecomplex* b, integer* ldb, doublecomplex* x, integer* ldx,
	doublereal* rcond, doublereal* ferr, doublereal* berr, doublecomplex*
	work, doublereal* rwork, integer* info);

CLAPACK_API int zpttrf_(integer* n, doublereal* d__, doublecomplex* e,
	integer* info);

CLAPACK_API int zpttrs_(char* uplo, integer* n, integer* nrhs,
	doublereal* d__, doublecomplex* e, doublecomplex* b, integer* ldb,
	integer* info);

CLAPACK_API int zptts2_(integer* iuplo, integer* n, integer* nrhs,
	doublereal* d__, doublecomplex* e, doublecomplex* b, integer* ldb);

CLAPACK_API int zrot_(integer* n, doublecomplex* cx, integer* incx,
	doublecomplex* cy, integer* incy, doublereal* c__, doublecomplex* s);

CLAPACK_API int zspcon_(char* uplo, integer* n, doublecomplex* ap,
	integer* ipiv, doublereal* anorm, doublereal* rcond, doublecomplex*
	work, integer* info);

CLAPACK_API int zspmv_(char* uplo, integer* n, doublecomplex* alpha,
	doublecomplex* ap, doublecomplex* x, integer* incx, doublecomplex*
	beta, doublecomplex* y, integer* incy);

CLAPACK_API int zspr_(char* uplo, integer* n, doublecomplex* alpha,
	doublecomplex* x, integer* incx, doublecomplex* ap);

CLAPACK_API int zsprfs_(char* uplo, integer* n, integer* nrhs,
	doublecomplex* ap, doublecomplex* afp, integer* ipiv, doublecomplex*
	b, integer* ldb, doublecomplex* x, integer* ldx, doublereal* ferr,
	doublereal* berr, doublecomplex* work, doublereal* rwork, integer*
	info);

CLAPACK_API int zspsv_(char* uplo, integer* n, integer* nrhs,
	doublecomplex* ap, integer* ipiv, doublecomplex* b, integer* ldb,
	integer* info);

CLAPACK_API int zspsvx_(char* fact, char* uplo, integer* n, integer*
	nrhs, doublecomplex* ap, doublecomplex* afp, integer* ipiv,
	doublecomplex* b, integer* ldb, doublecomplex* x, integer* ldx,
	doublereal* rcond, doublereal* ferr, doublereal* berr, doublecomplex*
	work, doublereal* rwork, integer* info);

CLAPACK_API int zsptrf_(char* uplo, integer* n, doublecomplex* ap,
	integer* ipiv, integer* info);

CLAPACK_API int zsptri_(char* uplo, integer* n, doublecomplex* ap,
	integer* ipiv, doublecomplex* work, integer* info);

CLAPACK_API int zsptrs_(char* uplo, integer* n, integer* nrhs,
	doublecomplex* ap, integer* ipiv, doublecomplex* b, integer* ldb,
	integer* info);

CLAPACK_API int zstedc_(char* compz, integer* n, doublereal* d__,
	doublereal* e, doublecomplex* z__, integer* ldz, doublecomplex* work,
	integer* lwork, doublereal* rwork, integer* lrwork, integer* iwork,
	integer* liwork, integer* info);

CLAPACK_API int zstegr_(char* jobz, char* range, integer* n, doublereal*
	d__, doublereal* e, doublereal* vl, doublereal* vu, integer* il,
	integer* iu, doublereal* abstol, integer* m, doublereal* w,
	doublecomplex* z__, integer* ldz, integer* isuppz, doublereal* work,
	integer* lwork, integer* iwork, integer* liwork, integer* info);

CLAPACK_API int zstein_(integer* n, doublereal* d__, doublereal* e,
	integer* m, doublereal* w, integer* iblock, integer* isplit,
	doublecomplex* z__, integer* ldz, doublereal* work, integer* iwork,
	integer* ifail, integer* info);

CLAPACK_API int zstemr_(char* jobz, char* range, integer* n, doublereal*
	d__, doublereal* e, doublereal* vl, doublereal* vu, integer* il,
	integer* iu, integer* m, doublereal* w, doublecomplex* z__, integer*
	ldz, integer* nzc, integer* isuppz, logical* tryrac, doublereal* work,
	integer* lwork, integer* iwork, integer* liwork, integer* info);

CLAPACK_API int zsteqr_(char* compz, integer* n, doublereal* d__,
	doublereal* e, doublecomplex* z__, integer* ldz, doublereal* work,
	integer* info);

CLAPACK_API int zsycon_(char* uplo, integer* n, doublecomplex* a,
	integer* lda, integer* ipiv, doublereal* anorm, doublereal* rcond,
	doublecomplex* work, integer* info);

CLAPACK_API int zsyequb_(char* uplo, integer* n, doublecomplex* a,
	integer* lda, doublereal* s, doublereal* scond, doublereal* amax,
	doublecomplex* work, integer* info);

CLAPACK_API int zsymv_(char* uplo, integer* n, doublecomplex* alpha,
	doublecomplex* a, integer* lda, doublecomplex* x, integer* incx,
	doublecomplex* beta, doublecomplex* y, integer* incy);

CLAPACK_API int zsyr_(char* uplo, integer* n, doublecomplex* alpha,
	doublecomplex* x, integer* incx, doublecomplex* a, integer* lda);

CLAPACK_API int zsyrfs_(char* uplo, integer* n, integer* nrhs,
	doublecomplex* a, integer* lda, doublecomplex* af, integer* ldaf,
	integer* ipiv, doublecomplex* b, integer* ldb, doublecomplex* x,
	integer* ldx, doublereal* ferr, doublereal* berr, doublecomplex* work,
	doublereal* rwork, integer* info);

CLAPACK_API int zsysv_(char* uplo, integer* n, integer* nrhs,
	doublecomplex* a, integer* lda, integer* ipiv, doublecomplex* b,
	integer* ldb, doublecomplex* work, integer* lwork, integer* info);

CLAPACK_API int zsysvx_(char* fact, char* uplo, integer* n, integer*
	nrhs, doublecomplex* a, integer* lda, doublecomplex* af, integer*
	ldaf, integer* ipiv, doublecomplex* b, integer* ldb, doublecomplex* x,
	integer* ldx, doublereal* rcond, doublereal* ferr, doublereal* berr,
	doublecomplex* work, integer* lwork, doublereal* rwork, integer* info);

CLAPACK_API int zsytf2_(char* uplo, integer* n, doublecomplex* a,
	integer* lda, integer* ipiv, integer* info);

CLAPACK_API int zsytrf_(char* uplo, integer* n, doublecomplex* a,
	integer* lda, integer* ipiv, doublecomplex* work, integer* lwork,
	integer* info);

CLAPACK_API int zsytri_(char* uplo, integer* n, doublecomplex* a,
	integer* lda, integer* ipiv, doublecomplex* work, integer* info);

CLAPACK_API int zsytrs_(char* uplo, integer* n, integer* nrhs,
	doublecomplex* a, integer* lda, integer* ipiv, doublecomplex* b,
	integer* ldb, integer* info);

CLAPACK_API int ztbcon_(char* norm, char* uplo, char* diag, integer* n,
	integer* kd, doublecomplex* ab, integer* ldab, doublereal* rcond,
	doublecomplex* work, doublereal* rwork, integer* info);

CLAPACK_API int ztbrfs_(char* uplo, char* trans, char* diag, integer* n,
	integer* kd, integer* nrhs, doublecomplex* ab, integer* ldab,
	doublecomplex* b, integer* ldb, doublecomplex* x, integer* ldx,
	doublereal* ferr, doublereal* berr, doublecomplex* work, doublereal*
	rwork, integer* info);

CLAPACK_API int ztbtrs_(char* uplo, char* trans, char* diag, integer* n,
	integer* kd, integer* nrhs, doublecomplex* ab, integer* ldab,
	doublecomplex* b, integer* ldb, integer* info);

CLAPACK_API int ztfsm_(char* transr, char* side, char* uplo, char* trans,
	char* diag, integer* m, integer* n, doublecomplex* alpha,
	doublecomplex* a, doublecomplex* b, integer* ldb);

CLAPACK_API int ztftri_(char* transr, char* uplo, char* diag, integer* n,
	doublecomplex* a, integer* info);

CLAPACK_API int ztfttp_(char* transr, char* uplo, integer* n,
	doublecomplex* arf, doublecomplex* ap, integer* info);

CLAPACK_API int ztfttr_(char* transr, char* uplo, integer* n,
	doublecomplex* arf, doublecomplex* a, integer* lda, integer* info);

CLAPACK_API int ztgevc_(char* side, char* howmny, logical* select,
	integer* n, doublecomplex* s, integer* lds, doublecomplex* p, integer
	* ldp, doublecomplex* vl, integer* ldvl, doublecomplex* vr, integer*
	ldvr, integer* mm, integer* m, doublecomplex* work, doublereal* rwork,
	integer* info);

CLAPACK_API int ztgex2_(logical* wantq, logical* wantz, integer* n,
	doublecomplex* a, integer* lda, doublecomplex* b, integer* ldb,
	doublecomplex* q, integer* ldq, doublecomplex* z__, integer* ldz,
	integer* j1, integer* info);

CLAPACK_API int ztgexc_(logical* wantq, logical* wantz, integer* n,
	doublecomplex* a, integer* lda, doublecomplex* b, integer* ldb,
	doublecomplex* q, integer* ldq, doublecomplex* z__, integer* ldz,
	integer* ifst, integer* ilst, integer* info);

CLAPACK_API int ztgsen_(integer* ijob, logical* wantq, logical* wantz,
	logical* select, integer* n, doublecomplex* a, integer* lda,
	doublecomplex* b, integer* ldb, doublecomplex* alpha, doublecomplex*
	beta, doublecomplex* q, integer* ldq, doublecomplex* z__, integer*
	ldz, integer* m, doublereal* pl, doublereal* pr, doublereal* dif,
	doublecomplex* work, integer* lwork, integer* iwork, integer* liwork,
	integer* info);

CLAPACK_API int ztgsja_(char* jobu, char* jobv, char* jobq, integer* m,
	integer* p, integer* n, integer* k, integer* l, doublecomplex* a,
	integer* lda, doublecomplex* b, integer* ldb, doublereal* tola,
	doublereal* tolb, doublereal* alpha, doublereal* beta, doublecomplex*
	u, integer* ldu, doublecomplex* v, integer* ldv, doublecomplex* q,
	integer* ldq, doublecomplex* work, integer* ncycle, integer* info);

CLAPACK_API int ztgsna_(char* job, char* howmny, logical* select,
	integer* n, doublecomplex* a, integer* lda, doublecomplex* b, integer
	* ldb, doublecomplex* vl, integer* ldvl, doublecomplex* vr, integer*
	ldvr, doublereal* s, doublereal* dif, integer* mm, integer* m,
	doublecomplex* work, integer* lwork, integer* iwork, integer* info);

CLAPACK_API int ztgsy2_(char* trans, integer* ijob, integer* m, integer*
	n, doublecomplex* a, integer* lda, doublecomplex* b, integer* ldb,
	doublecomplex* c__, integer* ldc, doublecomplex* d__, integer* ldd,
	doublecomplex* e, integer* lde, doublecomplex* f, integer* ldf,
	doublereal* scale, doublereal* rdsum, doublereal* rdscal, integer*
	info);

CLAPACK_API int ztgsyl_(char* trans, integer* ijob, integer* m, integer*
	n, doublecomplex* a, integer* lda, doublecomplex* b, integer* ldb,
	doublecomplex* c__, integer* ldc, doublecomplex* d__, integer* ldd,
	doublecomplex* e, integer* lde, doublecomplex* f, integer* ldf,
	doublereal* scale, doublereal* dif, doublecomplex* work, integer*
	lwork, integer* iwork, integer* info);

CLAPACK_API int ztpcon_(char* norm, char* uplo, char* diag, integer* n,
	doublecomplex* ap, doublereal* rcond, doublecomplex* work, doublereal
	* rwork, integer* info);

CLAPACK_API int ztprfs_(char* uplo, char* trans, char* diag, integer* n,
	integer* nrhs, doublecomplex* ap, doublecomplex* b, integer* ldb,
	doublecomplex* x, integer* ldx, doublereal* ferr, doublereal* berr,
	doublecomplex* work, doublereal* rwork, integer* info);

CLAPACK_API int ztptri_(char* uplo, char* diag, integer* n,
	doublecomplex* ap, integer* info);

CLAPACK_API int ztptrs_(char* uplo, char* trans, char* diag, integer* n,
	integer* nrhs, doublecomplex* ap, doublecomplex* b, integer* ldb,
	integer* info);

CLAPACK_API int ztpttf_(char* transr, char* uplo, integer* n,
	doublecomplex* ap, doublecomplex* arf, integer* info);

CLAPACK_API int ztpttr_(char* uplo, integer* n, doublecomplex* ap,
	doublecomplex* a, integer* lda, integer* info);

CLAPACK_API int ztrcon_(char* norm, char* uplo, char* diag, integer* n,
	doublecomplex* a, integer* lda, doublereal* rcond, doublecomplex*
	work, doublereal* rwork, integer* info);

CLAPACK_API int ztrevc_(char* side, char* howmny, logical* select,
	integer* n, doublecomplex* t, integer* ldt, doublecomplex* vl,
	integer* ldvl, doublecomplex* vr, integer* ldvr, integer* mm, integer
	* m, doublecomplex* work, doublereal* rwork, integer* info);

CLAPACK_API int ztrexc_(char* compq, integer* n, doublecomplex* t,
	integer* ldt, doublecomplex* q, integer* ldq, integer* ifst, integer*
	ilst, integer* info);

CLAPACK_API int ztrrfs_(char* uplo, char* trans, char* diag, integer* n,
	integer* nrhs, doublecomplex* a, integer* lda, doublecomplex* b,
	integer* ldb, doublecomplex* x, integer* ldx, doublereal* ferr,
	doublereal* berr, doublecomplex* work, doublereal* rwork, integer*
	info);

CLAPACK_API int ztrsen_(char* job, char* compq, logical* select, integer
	* n, doublecomplex* t, integer* ldt, doublecomplex* q, integer* ldq,
	doublecomplex* w, integer* m, doublereal* s, doublereal* sep,
	doublecomplex* work, integer* lwork, integer* info);

CLAPACK_API int ztrsna_(char* job, char* howmny, logical* select,
	integer* n, doublecomplex* t, integer* ldt, doublecomplex* vl,
	integer* ldvl, doublecomplex* vr, integer* ldvr, doublereal* s,
	doublereal* sep, integer* mm, integer* m, doublecomplex* work,
	integer* ldwork, doublereal* rwork, integer* info);

CLAPACK_API int ztrsyl_(char* trana, char* tranb, integer* isgn, integer
	* m, integer* n, doublecomplex* a, integer* lda, doublecomplex* b,
	integer* ldb, doublecomplex* c__, integer* ldc, doublereal* scale,
	integer* info);

CLAPACK_API int ztrti2_(char* uplo, char* diag, integer* n,
	doublecomplex* a, integer* lda, integer* info);

CLAPACK_API int ztrtri_(char* uplo, char* diag, integer* n,
	doublecomplex* a, integer* lda, integer* info);

CLAPACK_API int ztrtrs_(char* uplo, char* trans, char* diag, integer* n,
	integer* nrhs, doublecomplex* a, integer* lda, doublecomplex* b,
	integer* ldb, integer* info);

CLAPACK_API int ztrttf_(char* transr, char* uplo, integer* n,
	doublecomplex* a, integer* lda, doublecomplex* arf, integer* info);

CLAPACK_API int ztrttp_(char* uplo, integer* n, doublecomplex* a,
	integer* lda, doublecomplex* ap, integer* info);

CLAPACK_API int ztzrqf_(integer* m, integer* n, doublecomplex* a,
	integer* lda, doublecomplex* tau, integer* info);

CLAPACK_API int ztzrzf_(integer* m, integer* n, doublecomplex* a,
	integer* lda, doublecomplex* tau, doublecomplex* work, integer* lwork,
	integer* info);

CLAPACK_API int zung2l_(integer* m, integer* n, integer* k,
	doublecomplex* a, integer* lda, doublecomplex* tau, doublecomplex*
	work, integer* info);

CLAPACK_API int zung2r_(integer* m, integer* n, integer* k,
	doublecomplex* a, integer* lda, doublecomplex* tau, doublecomplex*
	work, integer* info);

CLAPACK_API int zungbr_(char* vect, integer* m, integer* n, integer* k,
	doublecomplex* a, integer* lda, doublecomplex* tau, doublecomplex*
	work, integer* lwork, integer* info);

CLAPACK_API int zunghr_(integer* n, integer* ilo, integer* ihi,
	doublecomplex* a, integer* lda, doublecomplex* tau, doublecomplex*
	work, integer* lwork, integer* info);

CLAPACK_API int zungl2_(integer* m, integer* n, integer* k,
	doublecomplex* a, integer* lda, doublecomplex* tau, doublecomplex*
	work, integer* info);

CLAPACK_API int zunglq_(integer* m, integer* n, integer* k,
	doublecomplex* a, integer* lda, doublecomplex* tau, doublecomplex*
	work, integer* lwork, integer* info);

CLAPACK_API int zungql_(integer* m, integer* n, integer* k,
	doublecomplex* a, integer* lda, doublecomplex* tau, doublecomplex*
	work, integer* lwork, integer* info);

CLAPACK_API int zungqr_(integer* m, integer* n, integer* k,
	doublecomplex* a, integer* lda, doublecomplex* tau, doublecomplex*
	work, integer* lwork, integer* info);

CLAPACK_API int zungr2_(integer* m, integer* n, integer* k,
	doublecomplex* a, integer* lda, doublecomplex* tau, doublecomplex*
	work, integer* info);

CLAPACK_API int zungrq_(integer* m, integer* n, integer* k,
	doublecomplex* a, integer* lda, doublecomplex* tau, doublecomplex*
	work, integer* lwork, integer* info);

CLAPACK_API int zungtr_(char* uplo, integer* n, doublecomplex* a,
	integer* lda, doublecomplex* tau, doublecomplex* work, integer* lwork,
	integer* info);

CLAPACK_API int zunm2l_(char* side, char* trans, integer* m, integer* n,
	integer* k, doublecomplex* a, integer* lda, doublecomplex* tau,
	doublecomplex* c__, integer* ldc, doublecomplex* work, integer* info);

CLAPACK_API int zunm2r_(char* side, char* trans, integer* m, integer* n,
	integer* k, doublecomplex* a, integer* lda, doublecomplex* tau,
	doublecomplex* c__, integer* ldc, doublecomplex* work, integer* info);

CLAPACK_API int zunmbr_(char* vect, char* side, char* trans, integer* m,
	integer* n, integer* k, doublecomplex* a, integer* lda, doublecomplex
	* tau, doublecomplex* c__, integer* ldc, doublecomplex* work, integer*
	lwork, integer* info);

CLAPACK_API int zunmhr_(char* side, char* trans, integer* m, integer* n,
	integer* ilo, integer* ihi, doublecomplex* a, integer* lda,
	doublecomplex* tau, doublecomplex* c__, integer* ldc, doublecomplex*
	work, integer* lwork, integer* info);

CLAPACK_API int zunml2_(char* side, char* trans, integer* m, integer* n,
	integer* k, doublecomplex* a, integer* lda, doublecomplex* tau,
	doublecomplex* c__, integer* ldc, doublecomplex* work, integer* info);

CLAPACK_API int zunmlq_(char* side, char* trans, integer* m, integer* n,
	integer* k, doublecomplex* a, integer* lda, doublecomplex* tau,
	doublecomplex* c__, integer* ldc, doublecomplex* work, integer* lwork,
	integer* info);

CLAPACK_API int zunmql_(char* side, char* trans, integer* m, integer* n,
	integer* k, doublecomplex* a, integer* lda, doublecomplex* tau,
	doublecomplex* c__, integer* ldc, doublecomplex* work, integer* lwork,
	integer* info);

CLAPACK_API int zunmqr_(char* side, char* trans, integer* m, integer* n,
	integer* k, doublecomplex* a, integer* lda, doublecomplex* tau,
	doublecomplex* c__, integer* ldc, doublecomplex* work, integer* lwork,
	integer* info);

CLAPACK_API int zunmr2_(char* side, char* trans, integer* m, integer* n,
	integer* k, doublecomplex* a, integer* lda, doublecomplex* tau,
	doublecomplex* c__, integer* ldc, doublecomplex* work, integer* info);

CLAPACK_API int zunmr3_(char* side, char* trans, integer* m, integer* n,
	integer* k, integer* l, doublecomplex* a, integer* lda, doublecomplex
	* tau, doublecomplex* c__, integer* ldc, doublecomplex* work, integer*
	info);

CLAPACK_API int zunmrq_(char* side, char* trans, integer* m, integer* n,
	integer* k, doublecomplex* a, integer* lda, doublecomplex* tau,
	doublecomplex* c__, integer* ldc, doublecomplex* work, integer* lwork,
	integer* info);

CLAPACK_API int zunmrz_(char* side, char* trans, integer* m, integer* n,
	integer* k, integer* l, doublecomplex* a, integer* lda, doublecomplex
	* tau, doublecomplex* c__, integer* ldc, doublecomplex* work, integer*
	lwork, integer* info);

CLAPACK_API int zunmtr_(char* side, char* uplo, char* trans, integer* m,
	integer* n, doublecomplex* a, integer* lda, doublecomplex* tau,
	doublecomplex* c__, integer* ldc, doublecomplex* work, integer* lwork,
	integer* info);

CLAPACK_API int zupgtr_(char* uplo, integer* n, doublecomplex* ap,
	doublecomplex* tau, doublecomplex* q, integer* ldq, doublecomplex*
	work, integer* info);

CLAPACK_API int zupmtr_(char* side, char* uplo, char* trans, integer* m,
	integer* n, doublecomplex* ap, doublecomplex* tau, doublecomplex* c__,
	integer* ldc, doublecomplex* work, integer* info);


//}} finished

#pragma endregion

#pragma region ZXLASRC -- Double precision complex LAPACK routines using extra precision (zlax)
//#       
//set(ZXLASRC
//	zherfsx.c zhesvxx.c 
//	zgbrfsx.c zgbsvxx.c 
//	zgerfsx.c zgesvxx.c
//	zporfsx.c zposvxx.c 
//	zsyrfsx.c zsysvxx.c 
//  zlarscl2.c zlascl2.c zla_wwaddw.c
//	zla_gerfsx_extended.c zla_geamv.c
//	zla_gercond_c.c zla_gercond_x.c zla_rpvgrw.c
//	zla_syrfsx_extended.c zla_syamv.c zla_syrcond_c.c zla_syrcond_x.c
//	zla_syrpvgrw.c zla_porfsx_extended.c
//	zla_porcond_c.c zla_porcond_x.c zla_porpvgrw.c
//	zla_gbrfsx_extended.c zla_gbamv.c zla_gbrcond_c.c zla_gbrcond_x.c
//	zla_gbrpvgrw.c zla_herfsx_extended.c
//	zla_heamv.c zla_hercond_c.c zla_hercond_x.c zla_herpvgrw.c
//	zla_lin_berr.c )


CLAPACK_API
int zlarscl2_(integer* m, integer* n, doublereal* d__,
	doublecomplex* x, integer* ldx);


/*

int zgbbrd_(char *vect, integer *m, integer *n, integer *ncc,
	 integer *kl, integer *ku, doublecomplex *ab, integer *ldab,
	doublereal *d__, doublereal *e, doublecomplex *q, integer *ldq,
	doublecomplex *pt, integer *ldpt, doublecomplex *c__, integer *ldc,
	doublecomplex *work, doublereal *rwork, integer *info);

int zgbrfsx_(char *trans, char *equed, integer *n, integer *
	kl, integer *ku, integer *nrhs, doublecomplex *ab, integer *ldab,
	doublecomplex *afb, integer *ldafb, integer *ipiv, doublereal *r__,
	doublereal *c__, doublecomplex *b, integer *ldb, doublecomplex *x,
	integer *ldx, doublereal *rcond, doublereal *berr, integer *
	n_err_bnds__, doublereal *err_bnds_norm__, doublereal *
	err_bnds_comp__, integer *nparams, doublereal *params, doublecomplex *
	work, doublereal *rwork, integer *info);

int zgbsvxx_(char *fact, char *trans, integer *n, integer *
	kl, integer *ku, integer *nrhs, doublecomplex *ab, integer *ldab,
	doublecomplex *afb, integer *ldafb, integer *ipiv, char *equed,
	doublereal *r__, doublereal *c__, doublecomplex *b, integer *ldb,
	doublecomplex *x, integer *ldx, doublereal *rcond, doublereal *rpvgrw,
	 doublereal *berr, integer *n_err_bnds__, doublereal *err_bnds_norm__,
	 doublereal *err_bnds_comp__, integer *nparams, doublereal *params,
	doublecomplex *work, doublereal *rwork, integer *info);

int zgerfsx_(char *trans, char *equed, integer *n, integer *
	nrhs, doublecomplex *a, integer *lda, doublecomplex *af, integer *
	ldaf, integer *ipiv, doublereal *r__, doublereal *c__, doublecomplex *
	b, integer *ldb, doublecomplex *x, integer *ldx, doublereal *rcond,
	doublereal *berr, integer *n_err_bnds__, doublereal *err_bnds_norm__,
	doublereal *err_bnds_comp__, integer *nparams, doublereal *params,
	doublecomplex *work, doublereal *rwork, integer *info);

int zgesvxx_(char *fact, char *trans, integer *n, integer *
	nrhs, doublecomplex *a, integer *lda, doublecomplex *af, integer *
	ldaf, integer *ipiv, char *equed, doublereal *r__, doublereal *c__,
	doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx,
	doublereal *rcond, doublereal *rpvgrw, doublereal *berr, integer *
	n_err_bnds__, doublereal *err_bnds_norm__, doublereal *
	err_bnds_comp__, integer *nparams, doublereal *params, doublecomplex *
	work, doublereal *rwork, integer *info);

int zherfsx_(char *uplo, char *equed, integer *n, integer *
	nrhs, doublecomplex *a, integer *lda, doublecomplex *af, integer *
	ldaf, integer *ipiv, doublereal *s, doublecomplex *b, integer *ldb,
	doublecomplex *x, integer *ldx, doublereal *rcond, doublereal *berr,
	integer *n_err_bnds__, doublereal *err_bnds_norm__, doublereal *
	err_bnds_comp__, integer *nparams, doublereal *params, doublecomplex *
	work, doublereal *rwork, integer *info);

int zhesvxx_(char *fact, char *uplo, integer *n, integer *
	nrhs, doublecomplex *a, integer *lda, doublecomplex *af, integer *
	ldaf, integer *ipiv, char *equed, doublereal *s, doublecomplex *b,
	integer *ldb, doublecomplex *x, integer *ldx, doublereal *rcond,
	doublereal *rpvgrw, doublereal *berr, integer *n_err_bnds__,
	doublereal *err_bnds_norm__, doublereal *err_bnds_comp__, integer *
	nparams, doublereal *params, doublecomplex *work, doublereal *rwork,
	integer *info);

int zla_gbamv__(integer *trans, integer *m, integer *n,
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

int zla_gbrfsx_extended__(integer *prec_type__, integer *
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

int zla_geamv__(integer *trans, integer *m, integer *n,
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

int zla_gerfsx_extended__(integer *prec_type__, integer *
	trans_type__, integer *n, integer *nrhs, doublecomplex *a, integer *
	lda, doublecomplex *af, integer *ldaf, integer *ipiv, logical *colequ,
	 doublereal *c__, doublecomplex *b, integer *ldb, doublecomplex *y,
	integer *ldy, doublereal *berr_out__, integer *n_norms__, doublereal *
	errs_n__, doublereal *errs_c__, doublecomplex *res, doublereal *ayb,
	doublecomplex *dy, doublecomplex *y_tail__, doublereal *rcond,
	integer *ithresh, doublereal *rthresh, doublereal *dz_ub__, logical *
	ignore_cwise__, integer *info);

int zla_heamv__(integer *uplo, integer *n, doublereal *alpha,
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

int zla_herfsx_extended__(integer *prec_type__, char *uplo,
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

int zla_lin_berr__(integer *n, integer *nz, integer *nrhs,
	doublecomplex *res, doublereal *ayb, doublereal *berr);

doublereal zla_porcond_c__(char *uplo, integer *n, doublecomplex *a, integer *
	lda, doublecomplex *af, integer *ldaf, doublereal *c__, logical *
	capply, integer *info, doublecomplex *work, doublereal *rwork, ftnlen
	uplo_len);

doublereal zla_porcond_x__(char *uplo, integer *n, doublecomplex *a, integer *
	lda, doublecomplex *af, integer *ldaf, doublecomplex *x, integer *
	info, doublecomplex *work, doublereal *rwork, ftnlen uplo_len);

int zla_porfsx_extended__(integer *prec_type__, char *uplo,
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

int zla_syamv__(integer *uplo, integer *n, doublereal *alpha,
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

int zla_syrfsx_extended__(integer *prec_type__, char *uplo,
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

int zla_wwaddw__(integer *n, doublecomplex *x, doublecomplex
	*y, doublecomplex *w);

int zporfsx_(char *uplo, char *equed, integer *n, integer *
	nrhs, doublecomplex *a, integer *lda, doublecomplex *af, integer *
	ldaf, doublereal *s, doublecomplex *b, integer *ldb, doublecomplex *x,
	 integer *ldx, doublereal *rcond, doublereal *berr, integer *
	n_err_bnds__, doublereal *err_bnds_norm__, doublereal *
	err_bnds_comp__, integer *nparams, doublereal *params, doublecomplex *
	work, doublereal *rwork, integer *info);

int zposvxx_(char *fact, char *uplo, integer *n, integer *
	nrhs, doublecomplex *a, integer *lda, doublecomplex *af, integer *
	ldaf, char *equed, doublereal *s, doublecomplex *b, integer *ldb,
	doublecomplex *x, integer *ldx, doublereal *rcond, doublereal *rpvgrw,
	 doublereal *berr, integer *n_err_bnds__, doublereal *err_bnds_norm__,
	 doublereal *err_bnds_comp__, integer *nparams, doublereal *params,
	doublecomplex *work, doublereal *rwork, integer *info);

int zsyrfsx_(char *uplo, char *equed, integer *n, integer *
	nrhs, doublecomplex *a, integer *lda, doublecomplex *af, integer *
	ldaf, integer *ipiv, doublereal *s, doublecomplex *b, integer *ldb,
	doublecomplex *x, integer *ldx, doublereal *rcond, doublereal *berr,
	integer *n_err_bnds__, doublereal *err_bnds_norm__, doublereal *
	err_bnds_comp__, integer *nparams, doublereal *params, doublecomplex *
	work, doublereal *rwork, integer *info);

int zsysvxx_(char *fact, char *uplo, integer *n, integer *
	nrhs, doublecomplex *a, integer *lda, doublecomplex *af, integer *
	ldaf, integer *ipiv, char *equed, doublereal *s, doublecomplex *b,
	integer *ldb, doublecomplex *x, integer *ldx, doublereal *rcond,
	doublereal *rpvgrw, doublereal *berr, integer *n_err_bnds__,
	doublereal *err_bnds_norm__, doublereal *err_bnds_comp__, integer *
	nparams, doublereal *params, doublecomplex *work, doublereal *rwork,
	integer *info);




*/

#pragma endregion



