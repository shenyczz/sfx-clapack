#include "pch.h"
#include "CppUnitTest.h"

#include "clapack.h"

#pragma comment(lib,"cblas.lib")
#pragma comment(lib,"clapack.lib")

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace utsclapack
{
	TEST_CLASS(utsclapack)
	{
	public:
		
		TEST_METHOD(TestMethod1)
		{
			blas_dasum();
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

	};
}
