#ifndef REMOVEPLANARITYTEST_H_
#define REMOVEPLANARITYTEST_H_

#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>

#ifndef MOLECULE_HPP_
class Molecule;
#endif

#ifndef STRUCTURE_SIMULATOR_HPP_
class StructureSimulator;
#endif

class RemovePlanarityTest : public CppUnit::TestCase {
 private:
  CPPUNIT_TEST_SUITE(RemovePlanarityTest);
  CPPUNIT_TEST(testRemovePlanarity);
  CPPUNIT_TEST_SUITE_END();
  Molecule *CurrMol;
  StructureSimulator *Sim;

 public:
  RemovePlanarityTest();
  ~RemovePlanarityTest();
  void testRemovePlanarity();
};

#endif /* REMOVEPLANARITYTEST_H_ */
