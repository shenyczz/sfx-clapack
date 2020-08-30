#include <iostream>
#include "cblas.h"

void ddot_test()
{
    doublereal dx[] = { 1.1, 2.2, 3.3, 4.4 };
    doublereal dy[] = { 1.1, 2.2, 3.3, 4.4 };

    integer n = 4;
    integer incx = 1;
    integer incy = 1;

    auto dot = ddot_(&n, dx, &incx, dy, &incy);
    //auto sum = dasum_(&n, dx, &incx);

    std::cout << "ddot_ = " << dot << std::endl;
}

