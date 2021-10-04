#ifndef STRUCTURESIMULATORTEST_HPP_
#define STRUCTURESIMULATORTEST_HPP_

#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>

#ifndef MOLECULE_HPP_
class Molecule;
#endif

#ifndef STRUCTURE_SIMULATOR_HPP_
class StructureSimulator;
#endif

class StructureSimulatorTest : public CppUnit::TestCase {
 private:
  CPPUNIT_TEST_SUITE(StructureSimulatorTest);
  CPPUNIT_TEST(testSimplex);
  CPPUNIT_TEST_SUITE_END();
  Molecule *CurrMol;
  StructureSimulator *Simulator;

 public:
  StructureSimulatorTest();
  ~StructureSimulatorTest();
  void testSimplex();
};


#endif /* STRUCTURESIMULATORTEST_HPP_ */
