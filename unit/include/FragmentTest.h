
#ifndef FRAGMENTTEST_H_
#define FRAGMENTTEST_H_

#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>

#ifndef INCLUDE_FRAGMENT_HPP_
class Fragment;
#endif

class FragmentTest : public CppUnit::TestCase, public Fragment {
 private:
  CPPUNIT_TEST_SUITE(FragmentTest);
  CPPUNIT_TEST(testcalculate_distance);
  CPPUNIT_TEST(testcalculate_distance_derivative);
  CPPUNIT_TEST(testcalculate_angle);
  CPPUNIT_TEST(testcalculate_angle_derivative);
  CPPUNIT_TEST(testcalculate_rdc_angle);
  CPPUNIT_TEST(testcalculate_rdc_angle_derivative);
  CPPUNIT_TEST(mpi_ff_test);
  CPPUNIT_TEST(testperform_optimization_step);
  CPPUNIT_TEST_SUITE_END();

 public:
  FragmentTest();
  void testcalculate_distance();
  void testcalculate_distance_derivative();
  void testcalculate_angle();
  void testcalculate_angle_derivative();
  void testcalculate_rdc_angle();
  void testcalculate_rdc_angle_derivative();
  void mpi_ff_test();
  void testperform_optimization_step();
};

#endif /* FRAGMENTTEST_H_ */
