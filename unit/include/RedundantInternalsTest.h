
#ifndef REDUNDANTINTERNALSTEST_H_
#define REDUNDANTINTERNALSTEST_H_

#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>

#ifndef MOLECULE_HPP_
class Molecule;
#endif

#ifndef STRUCTURE_SIMULATOR_HPP_
class StructureSimulator;
#endif

class RedundantInternalsTest : public CppUnit::TestCase {
 private:
  CPPUNIT_TEST_SUITE(RedundantInternalsTest);
  CPPUNIT_TEST(testRedundantSetup);
  CPPUNIT_TEST_SUITE_END();
  Molecule *CurrMol;
  StructureSimulator *Sim;

 public:
  RedundantInternalsTest();
  ~RedundantInternalsTest();
  void testRedundantSetup();
};

#endif /* REDUNDANTINTERNALSTEST_H_ */
