
#ifndef SIMPLEXTEST_HPP_
#define SIMPLEXTEST_HPP_

#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>

class SimplexTest : public CppUnit::TestCase {
 private:
  CPPUNIT_TEST_SUITE(SimplexTest);
  // TODO(vsz) move the first two to some math util module?
  CPPUNIT_TEST(testNLoptExample);
  CPPUNIT_TEST(testNLoptRosenbrock2D_NM);
  CPPUNIT_TEST(testSimplexNLoptRosenbrock2D_NM);
  CPPUNIT_TEST_SUITE_END();

 public:
  SimplexTest();
  ~SimplexTest();
  void testNLoptExample();
  void testNLoptRosenbrock2D_NM();
  void testSimplexNLoptRosenbrock2D_NM();
};


#endif /* SIMPLEXTEST_HPP_ */
