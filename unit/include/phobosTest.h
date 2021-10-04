
#ifndef UNIT_INCLUDE_PHOBOSTEST_H_
#define UNIT_INCLUDE_PHOBOSTEST_H_

#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>

// use lower case here to be consistent with file name
class phobosTest : public CppUnit::TestCase {
 private:
  CPPUNIT_TEST_SUITE(phobosTest);
  CPPUNIT_TEST(testPhobosRosenbrock);
  CPPUNIT_TEST(testPhobosBeale);
  CPPUNIT_TEST(testPhobosRastriginNumerical);
  CPPUNIT_TEST_SUITE_END();

 public:
  phobosTest();
  void testPhobosRosenbrock();
  void testPhobosRastriginNumerical();
  void testPhobosBeale();
};


#endif /* UNIT_INCLUDE_PHOBOSTEST_H_ */
