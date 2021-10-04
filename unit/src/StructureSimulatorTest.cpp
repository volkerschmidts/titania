#include <Atom.hpp>
#include <LinAlg.hpp>
#include <Molecule.hpp>
#include <Output.hpp>
#include <Simplex.hpp>
#include <StanIPC.h>
#include <Structure.hpp>
#include <StructureSimulator.hpp>
#include <StructureSimulatorTest.hpp>
// there is still some weird include ordering dependence hidden in
// Declarations.hpp
// clang-format off
#include <Declarations.hpp>
// clang-format on
#include <cmath>
#include <main.h>
#include <nlopt.hpp>
#include <vector>

CPPUNIT_TEST_SUITE_REGISTRATION(StructureSimulatorTest);

StructureSimulatorTest::StructureSimulatorTest()
{
  CurrMol = ipcStan();

  Flags flags;
  Simulator = new StructureSimulator(CurrMol->getHeadStruc(), flags);
}

StructureSimulatorTest::~StructureSimulatorTest()
{
  delete CurrMol;
  delete Simulator;
}

void
StructureSimulatorTest::testSimplex()
{
  Structure *CurrStruc = CurrMol->getHeadStruc();

  // setup molecule with distorted geometry
  unsigned int displace_atom_index = 4; // C4
  Eigen::Vector3d displacement(2.0, 2.0, 2.0);
  for (Atom *CurrAtom = CurrMol->getHeadStruc()->getHeadAtom(); CurrAtom;
       CurrAtom       = CurrAtom->getNext())
  {
    if (CurrAtom->getIndex() == displace_atom_index)
    {
      CurrAtom->shiftCoordinates(displacement, StructureOptions::Initial);
    }
  }

  std::string titania_unit_path = TITANIA_UNIT_TESTS_DIR__;
  BasicInformation baseInformation;
  baseInformation.xyzFile.open(
      (titania_unit_path +
       "/systems/cppunit_input.tna.out.StructureSimulatorSimplex.xyz"),
      std::ios::binary | std::ios::out | std::ios::trunc);
  outputXYZ(CurrStruc->getHeadAtom(), baseInformation,
            StructureOptions::Initial, "distorted structure");

  bool reduced        = false;
  unsigned int NMAX   = 1e6;
  double initial_step = 0.5;
  double Etol         = 1e-6;
  unsigned int nfunc =
      Simulator->startSimplexMinimizer(NMAX, initial_step, Etol, reduced);

  outputXYZ(CurrStruc->getHeadAtom(), baseInformation,
            StructureOptions::Initial, "re-optimized structure");

  double E_final = Simulator->totalEnergy();
  Tensor *ITensor_final =
      CurrStruc->CalculateInertiaTensor(StructureOptions::Initial);
  Eigen::Vector3d MomentsOfInertia_final = ITensor_final->getEvals();

  // this are a very relaxed tolerances -- but any geometry with a total energy
  // and Moments of Inertia in these ranges should be a good enough starting
  // point for TITANIA
  double E_expected = 95.0;
  Eigen::Vector3d MomentsOfInertia_expected(306.0, 523.0, 605.0);

  CPPUNIT_ASSERT_DOUBLES_EQUAL(E_expected, E_final, 5);
  for (unsigned int i = 0; i < 3; ++i)
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(MomentsOfInertia_expected(i),
                                 MomentsOfInertia_final(i), 5);
  }
}
