
#ifndef LINALGTEST_H_
#define LINALGTEST_H_

#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>

class LinAlgTest : public CppUnit::TestCase {
 private:
  CPPUNIT_TEST_SUITE(LinAlgTest);
  CPPUNIT_TEST(testMoorePenroseInverse);
  CPPUNIT_TEST(testYlm);
  CPPUNIT_TEST(testFactorial);
  CPPUNIT_TEST(testTensor2Euler);
  CPPUNIT_TEST(testdd2Mndb);
  CPPUNIT_TEST(testPolar2Eigen);
  CPPUNIT_TEST(testEigen2Polar);
  CPPUNIT_TEST(testSphere2B);
  CPPUNIT_TEST(testSv2St);
  CPPUNIT_TEST(testlinear_Regression);
  CPPUNIT_TEST(testaxb3d);
  CPPUNIT_TEST(testS3_permutation_matrix);
  CPPUNIT_TEST(testS3_permutation_index);
  CPPUNIT_TEST(testS3_permutation_order);
  CPPUNIT_TEST(testS3_permutation_parity);
  CPPUNIT_TEST(testderived_Leibnitz_determinante);
  CPPUNIT_TEST_SUITE_END();

 public:
  LinAlgTest();
  void testMoorePenroseInverse();
  void testYlm();
  void testFactorial();
  void testTensor2Euler();
  void testdd2Mndb();
  void testSphere2B();
  void testPolar2Eigen();
  void testEigen2Polar();
  void testSv2St();
  void testlinear_Regression();
  void testaxb3d();
  void testS3_permutation_matrix();
  void testS3_permutation_index();
  void testS3_permutation_order();
  void testS3_permutation_parity();
  void testderived_Leibnitz_determinante();
};

#endif /* LINALGTEST_H_ */
