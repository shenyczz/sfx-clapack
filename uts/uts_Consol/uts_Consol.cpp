// uts_Consol.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>

#pragma comment(lib,"cblas.lib")
#pragma comment(lib,"clapack.lib")


int main()
{
    extern void ddot_test();
    ddot_test();

    std::cout << "Hello World!\n";
}

