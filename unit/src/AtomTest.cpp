
#include <Atom.hpp>
#include <AtomTest.h>
#include <Declarations.hpp>
#include <Molecule.hpp>
#include <StanIPC.h>


#ifndef ALLOWED_ERROR
#define ALLOWED_ERROR 1e-3
#endif

#define GAMMA_1H  26.7513e7
#define GAMMA_2H  4.1065e7
#define GAMMA_13C 6.7262e7
#define GAMMA_15N -2.7116e7

#define x_A1 1.38
#define y_A1 1.69
#define z_A1 -4.87

#define x_A2 -1.11
#define y_A2 1.51
#define z_A2 -3.91

#define x_A3 1.11
#define y_A3 -1.51
#define z_A3 3.75

#define r12 2.6747
#define r13 9.1988
#define r23 8.5279

CPPUNIT_TEST_SUITE_REGISTRATION(AtomTest);

AtomTest::AtomTest()
{
  // Setup some atoms using the element + label convention

  A1              = new Atom();
  std::string a1l = "C1";
  atomLabelParser(a1l, A1);
  A1->setIndex(0);
  A1->setrdcIndex(2);
  A1->setCoordinates(x_A1, y_A1, z_A1, StructureOptions::Optimized);
  // Prepare one atom for the isotope test

  A2              = new Atom();
  std::string a2l = "2H1a";
  atomLabelParser(a2l, A2);
  A2->setIndex(1);
  A2->setrdcIndex(1);
  A2->setCoordinates(x_A2, y_A2, z_A2, StructureOptions::Optimized);

  A3              = new Atom();
  std::string a3l = "15N2";
  atomLabelParser(a3l, A3);
  A3->setIndex(3);
  A3->setrdcIndex(3);
  A3->setCoordinates(x_A3, y_A3, z_A3, StructureOptions::Optimized);

  A4              = new Atom();
  std::string a4l = "H2";
  atomLabelParser(a4l, A4);
  A4->setIndex(4);
  A4->setrdcIndex(4);
}


void
AtomTest::testAtomCoordinates()
{
  // Test the atom coordinate getter
  CPPUNIT_ASSERT_EQUAL(x_A1,
                       A1->getCoordinates(StructureOptions::Optimized)->x);
  CPPUNIT_ASSERT_EQUAL(y_A1,
                       A1->getCoordinates(StructureOptions::Optimized)->y);
  CPPUNIT_ASSERT_EQUAL(z_A1,
                       A1->getCoordinates(StructureOptions::Optimized)->z);
  CPPUNIT_ASSERT_EQUAL(0.0, A1->getCoordinates(StructureOptions::Initial)->x);
  CPPUNIT_ASSERT_EQUAL(0.0, A1->getCoordinates(StructureOptions::Initial)->y);
  CPPUNIT_ASSERT_EQUAL(0.0, A1->getCoordinates(StructureOptions::Initial)->z);

  A1->setCoordinates(0.1, 0.2, 0.3, StructureOptions::Initial);
  A1->retainCoordinates();
  CPPUNIT_ASSERT_EQUAL(0.1, A1->getCoordinates(StructureOptions::Optimized)->x);
  CPPUNIT_ASSERT_EQUAL(0.2, A1->getCoordinates(StructureOptions::Optimized)->y);
  CPPUNIT_ASSERT_EQUAL(0.3, A1->getCoordinates(StructureOptions::Optimized)->z);
  CPPUNIT_ASSERT_EQUAL(0.1, A1->getCoordinates(StructureOptions::Initial)->x);
  CPPUNIT_ASSERT_EQUAL(0.2, A1->getCoordinates(StructureOptions::Initial)->y);
  CPPUNIT_ASSERT_EQUAL(0.3, A1->getCoordinates(StructureOptions::Initial)->z);

  // Reset the coordinates to the initial value;
  A1->setCoordinates(.0, .0, .0, StructureOptions::Initial);
  A1->setCoordinates(x_A1, y_A1, z_A1, StructureOptions::Optimized);
}

void
AtomTest::testAtomLabels()
{
  CPPUNIT_ASSERT(A1->getIdentifier() == "C1");
  CPPUNIT_ASSERT(A1->getLabel() == "1");
  CPPUNIT_ASSERT(A2->getIdentifier() == "H1a");
  CPPUNIT_ASSERT(A2->getLabel() == "1a");
}

void
AtomTest::testAtomIsotopes()
{
  CPPUNIT_ASSERT_EQUAL(6, A1->getZ());
  CPPUNIT_ASSERT_EQUAL(13, A1->getA());
  CPPUNIT_ASSERT(A1->getElement() == "C");
  CPPUNIT_ASSERT_EQUAL(1, A2->getZ());
  CPPUNIT_ASSERT_EQUAL(2, A2->getA());
  CPPUNIT_ASSERT(A2->getElement() == "H");
}

void
AtomTest::testSetGamma()
{
  // Difference should be normalized on gamma since the values are around 1e7
  CPPUNIT_ASSERT(fabs((A1->getGamma() - GAMMA_13C) / GAMMA_13C) <
                 ALLOWED_ERROR);
  CPPUNIT_ASSERT(fabs((A2->getGamma() - GAMMA_2H) / GAMMA_2H) < ALLOWED_ERROR);
  CPPUNIT_ASSERT(fabs((A3->getGamma() - GAMMA_15N) / GAMMA_15N) <
                 ALLOWED_ERROR);
  CPPUNIT_ASSERT(fabs((A4->getGamma() - GAMMA_1H) / GAMMA_1H) < ALLOWED_ERROR);
}

void
AtomTest::testGetDistance()
{
  CPPUNIT_ASSERT(fabs(r12 - A1->getDistance(A2, StructureOptions::Optimized)) <
                 ALLOWED_ERROR);
  CPPUNIT_ASSERT(fabs(r13 - A1->getDistance(A3, StructureOptions::Optimized)) <
                 ALLOWED_ERROR);
  CPPUNIT_ASSERT(fabs(r23 - A2->getDistance(A3, StructureOptions::Optimized)) <
                 ALLOWED_ERROR);
}
