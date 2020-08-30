using System.Runtime.InteropServices;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace uts_clapack_all
{
    [TestClass]
    public class UnitTestClapack
    {
        [TestMethod]
        public void TestMethod1()
        {
            var v1 = fn_clapack();
            Assert.AreEqual(v1, 100);
        }

        [DllImport("x86_clapack.dll")]
        static extern int fn_clapack();


        //}}@@@
    }
}
