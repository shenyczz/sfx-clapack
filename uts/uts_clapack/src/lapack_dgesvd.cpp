#include "pch.h"
#include "CppUnitTest.h"
using namespace Microsoft::VisualStudio::CppUnitTestFramework;
#include "clapack.h"


/*
DGESVD computes the singular value decomposition (SVD) of a real
M-by-N matrix A, optionally computing the left and/or right singular
vectors. The SVD is written
       A = U * SIGMA * transpose(V)
where SIGMA is an M-by-N matrix which is zero except for its min(m,n) diagonal elements,
      U is an M-by-M orthogonal matrix(Õý½»¾ØÕó),
	  V is an N-by-N orthogonal matrix.

The diagonal elements of SIGMA are the singular values of A; they are real and non-negative,
and are returned in descending order.
The first min(m,n) columns of U and V are the left and right singular vectors of A.

Note that the routine returns V**T, not V.

Arguments
=========
	JOBU	(input) CHARACTER*1
			Specifies options for computing all or part of the matrix U:
			= 'A':  all M columns of U are returned in array U:
			= 'S':  the first min(m,n) columns of U (the left singular
			        vectors) are returned in the array U
			= 'O':  the first min(m,n) columns of U (the left singular
			        vectors) are overwritten on the array A;
			= 'N':  no columns of U (no left singular vectors) are
                    computed.




*/
void lapack_dgesvd()
{
	const int M = 5, N = 3;
	double A[] = { 18.91, 14.91, -6.15, -18.15, 27.5, -1.59, -1.59,
		-2.25,  -1.59, -2.25, -1.59, 1.59, 0.0, 1.59, 0.0 };

	double u[M * M];
	double vt[N * N];

	double s[M * N];

	const int nwork = max(3 * min(M, N) + max(M, N), 5 * min(M, N));
	double work[nwork];

	//int dgesvd_(char* jobu, char* jobvt, integer * m, integer * n,
	//	doublereal* a, integer * lda, doublereal * s, doublereal * u, integer *
	//	ldu, doublereal * vt, integer * ldvt, doublereal * work, integer * lwork,
	//	integer * info);

	char jobu = 'A', jobvt = 'A';
	//int lda,
	//int ret = dgesvd_(&jobu, &jobvt, (integer*)&M, (integer*)&N, a, 0, s, u, 0, vt, 0, work, 0, 0);

	//dgesvd_
	//org.netlib.util.intW info = new org.netlib.util.intW(2);
	//Dgesvd.dgesvd("A", "A", M, N, m, 0, M, s, 0, u, 0, M, vt, 0, N, work, 0, work.length, info);
	//System.out.println("info = " + info.val);


	return;
}
