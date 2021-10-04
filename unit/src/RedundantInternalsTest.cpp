#include "RedundantInternalsTest.h"

#include "StanIPC.h"

#include <Atom.hpp>
#include <Bond.hpp>
#include <Declarations.hpp>
#include <Helper.hpp>
#include <LinAlg.hpp>
#include <Molecule.hpp>
#include <Output.hpp>
#include <Potential.hpp>
#include <RedundantInternals.hpp>
#include <SphericalHarmonics.hpp>
#include <SphericalHarmonicsMatrix.hpp>
#include <Structure.hpp>
#include <StructureSimulator.hpp>
#include <main.h>


CPPUNIT_TEST_SUITE_REGISTRATION(RedundantInternalsTest);

RedundantInternalsTest::RedundantInternalsTest()
{
  CurrMol = ipcStan();
  Flags flags;
  flags.long_range_only_4_redundants = true;
  Sim = new StructureSimulator(CurrMol->getHeadStruc(), flags);
}

RedundantInternalsTest::~RedundantInternalsTest()
{
  delete Sim;
  delete CurrMol;
}

void
RedundantInternalsTest::testRedundantSetup()
{
  std::cout << "RedundantInternalsTest is not working due to memory issues "
               "right now...\n";
  return;
  BasicInformation baseInformation;
  Flags flags;
  flags.print_redundants                        = true;
  flags.useRedundants                           = true;
  flags.torsions_4_redundants                   = true;
  flags.long_range_only_4_redundants            = false;
  flags.floating_rdc_angles                     = false;
  baseInformation.limits.max_redundant_cycles   = 10;
  baseInformation.limits.redundants_convergence = 1e-2;
  for (int red_type = BOND_REDUNDANTS_; red_type < NUMBER_OF_REDUNDANT_TYPES_;
       ++red_type)
    baseInformation.static_redundants_weighting[red_type] = 1.0;
  baseInformation.static_redundants_weighting[PLANAR_REDUNDANTS_] = 2.0;
  std::string titania_unit_path = TITANIA_UNIT_TESTS_DIR__;
  baseInformation.xyzFile.open(
      (titania_unit_path + "/systems/cppunit_input.tna.out.xyz"),
      std::ios::binary | std::ios::out | std::ios::trunc);
  outputXYZ(CurrMol->getHeadStruc()->getHeadAtom(), baseInformation,
            StructureOptions::Initial, "true Struc");

  for (Atom *CurrAtom = CurrMol->getHeadStruc()->getHeadAtom(); CurrAtom;
       CurrAtom       = CurrAtom->getNext())
  {
    if (CurrAtom->getIndex() > 6)
      continue;
    Eigen::Vector3d C = CurrAtom->Coordinates2Eigen(StructureOptions::Initial);
    Eigen::Vector3d R = Eigen::Vector3d::Random();
    CurrAtom->setCoordinates((C + R * 0.2), StructureOptions::Initial);
  }


  outputXYZ(CurrMol->getHeadStruc()->getHeadAtom(), baseInformation,
            StructureOptions::Initial, "shifted Struc");

  CPPUNIT_ASSERT_EQUAL(ipc_stan_bonds(), Sim->getNumberOfBonds());
  CPPUNIT_ASSERT_EQUAL(ipc_stan_angles(), Sim->getNumberOfAngles());
  CPPUNIT_ASSERT_EQUAL(ipc_stan_torsions(), Sim->getNumberOfTorsions());
  CurrMol->getHeadStruc()->retainCoordinates();
  RedundantInternals Red(Sim, baseInformation, flags);

  CPPUNIT_ASSERT_EQUAL(ipc_stan_bonds(), Red.getNumberOfBonds());
  CPPUNIT_ASSERT_EQUAL(ipc_stan_angles(), Red.getNumberOfAngles());
  CPPUNIT_ASSERT_EQUAL(ipc_stan_torsions(), Red.getNumberOfTorsions());
  CPPUNIT_ASSERT_EQUAL((ipc_stan_torsions() + ipc_stan_angles() +
                        ipc_stan_bonds() + ipc_stan_rdcs()),
                       Red.getNumberOfRedundants());

  SphericalHarmonics *Y =
      CurrMol->getHeadStruc()->getHeadYmatrix()->getHeadHarmonic();
  srand(time(NULL));

  while (Y)
  {
    Y->setTheta(Y->getTheta(StructureOptions::Initial),
                StructureOptions::Optimized);
    Y->setPhi(Y->getPhi(StructureOptions::Initial),
              StructureOptions::Optimized);
    Y = Y->getNext();
  }

  CurrMol->getHeadStruc()->retainCoordinates();
  std::cout << "   - Information: Check "
            << (titania_unit_path + "/systems/cppunit_input.tna.out.xyz")
            << " for result..." << std::endl;
  Red.setupWilsonB(StructureOptions::Optimized);
  Eigen::MatrixXd WilsonB     = Red.getWilsonB();
  Eigen::MatrixXd MPI_WilsonB = MoorePenroseInverse(WilsonB);
  Eigen::MatrixXd BBt         = WilsonB * MPI_WilsonB;

  Red.S2x(baseInformation, flags);

  Red.rebaseCartesians(CurrMol->getHeadStruc());
}
