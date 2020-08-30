#pragma once

#include "_f2c.h"

CBLAS_API
doublereal dasum_(integer* n,
	doublereal* dx, integer* incx);

/* Subroutine */
CBLAS_API
int daxpy_(integer* n,
	doublereal* da,
	doublereal* dx,	integer* incx,
	doublereal* dy, integer* incy);


CBLAS_API
doublereal dcabs1_(doublecomplex* z__);


/* Subroutine */
CBLAS_API
int dcopy_(integer* n,
	doublereal* dx, integer* incx,
	doublereal* dy, integer* incy);

CBLAS_API
doublereal ddot_(integer* n,
	doublereal* dx, integer* incx,
	doublereal* dy,	integer* incy);






