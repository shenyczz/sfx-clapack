#include "pch.h"
#include "CppUnitTest.h"

#include "cblas.h"
#pragma comment(lib,"cblas.lib")

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace utscblas
{
	TEST_CLASS(utscblas)
	{
	public:
		
		TEST_METHOD(TestMethod1)
		{
			blas_dasum();
			blas_ddot();
		}



		void blas_dasum()
		{
			doublereal dx[] = { 1.1, 2.2, 3.3, 4.4 };
			doublereal dy[] = { 1.1, 2.2, 3.3, 4.4 };

			integer n = 4;
			integer incx = 1;
			integer incy = 1;

			auto sum = dasum_(&n, dx, &incx);

			Assert::AreEqual(sum, 11.0);
		}

		void blas_ddot()
		{
			doublereal dx[] = { 1.1, 2.2, 3.3, 4.4 };
			doublereal dy[] = { 1.1, 2.2, 3.3, 4.4 };

			integer n = 4;
			integer incx = 1;
			integer incy = 1;

			auto dot = ddot_(&n, dx, &incx, dy, &incy);
			//auto sum = dasum_(&n, dx, &incx);

			Assert::AreEqual(dot, 36.3);
		}


		//}}@@@
	};




}
