
#ifndef ATOMTEST_H_
#define ATOMTEST_H_

#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>

#ifndef ATOM_HPP_
class Atom;
#endif


class AtomTest : public CppUnit::TestCase {
 private:
  CPPUNIT_TEST_SUITE(AtomTest);
  CPPUNIT_TEST(testAtomCoordinates);
  CPPUNIT_TEST(testAtomLabels);
  CPPUNIT_TEST(testAtomIsotopes);
  CPPUNIT_TEST(testSetGamma);
  CPPUNIT_TEST(testGetDistance);
  CPPUNIT_TEST_SUITE_END();
  Atom *A1, *A2, *A3, *A4;

 public:
  AtomTest();
  void testAtomCoordinates();
  void testAtomLabels();
  void testAtomIsotopes();
  void testSetGamma();
  void testGetDistance();
};

#endif /* ATOMTEST_H_ */
