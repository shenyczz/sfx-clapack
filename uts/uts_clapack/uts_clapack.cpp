#include "pch.h"
#include "CppUnitTest.h"
using namespace Microsoft::VisualStudio::CppUnitTestFramework;

#include "clapack.h"

//#pragma comment(lib,"cblas.lib")
#pragma comment(lib,"clapack.lib")


namespace utsclapack
{
	TEST_CLASS(utsclapack)
	{
	public:
		
		TEST_METHOD(TCLAPACK)
		{
			void lapack_dgesvd();
			lapack_dgesvd();
		}


		//}}@@@
	};
}
