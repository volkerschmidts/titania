#include "RemovePlanarityTest.hpp"

#include "StanIPC.h"
#include "main.h"

#include <Atom.hpp>
#include <Declarations.hpp>
#include <LinAlg.hpp>
#include <Molecule.hpp>
#include <Output.hpp>
#include <Structure.hpp>
#include <StructureSimulator.hpp>

CPPUNIT_TEST_SUITE_REGISTRATION(RemovePlanarityTest);

RemovePlanarityTest::RemovePlanarityTest()
{
  CurrMol = ipcStan();
  for (Atom *CurrAtom = CurrMol->getHeadStruc()->getHeadAtom(); CurrAtom;
       CurrAtom       = CurrAtom->getNext())
  {
    Eigen::Vector3d coord =
        CurrAtom->Coordinates2Eigen(StructureOptions::Initial);
    // flatten initial coordinates along z-axis
    coord(2) = 0.0;
    CurrAtom->setCoordinates(coord, StructureOptions::Initial);
  }

  Flags flags;
  Sim = new StructureSimulator(CurrMol->getHeadStruc(), flags);
}

RemovePlanarityTest::~RemovePlanarityTest()
{
  delete Sim;
  delete CurrMol;
}

void
RemovePlanarityTest::testRemovePlanarity()
{
  // fix random seed to something we know works, otherwise the
  // constraints might be too strict for a randomly large displacement
  srand(1632403515);
  BasicInformation baseInformation;
  Structure *CurrStruc = CurrMol->getHeadStruc();

  std::string titania_unit_path = TITANIA_UNIT_TESTS_DIR__;
  baseInformation.xyzFile.open(
      (titania_unit_path +
       "/systems/cppunit_input.tna.out.removePlanarity.xyz"),
      std::ios::binary | std::ios::out | std::ios::trunc);
  outputXYZ(CurrStruc->getHeadAtom(), baseInformation,
            StructureOptions::Initial, "flattened structure");

  removePlanarity(CurrStruc, *Sim, baseInformation);

  // For the IPC geometry with all z-coordinates set to 0.0 (see constructor
  // above), the initial eigensystem of the inertia tensor is:
  //
  // Moments of Inertia
  // 174.659920 415.561870 590.221790
  //
  // After the coordinates are mangled by the Simplex Minimizer called by
  // removePlanarity() the moments of inertia will no longer satisfy the
  // planarity condition (I_z = Ix + Iy).
  // The actual values of the moments of inertia will be different for each
  // diastereoisomer, and as the function uses rand() to mutate the initial
  // planar geometry, we don't know at which diastereoisomer the minimizer
  // will arrive.
  //
  // The threshold in this test is set to a relaxed / optimized
  // value of 200. For a proper, relaxed IPC geometry, a value around
  // I_z - (Ix + Iy) = -232 is expected. With the limited constraints imposed
  // in removePlanarity(), this may not be reached. Manual calculations with
  // the current implementation and the above seed yield a value of -221.068025.
  //
  double MINIMUM_PLANARITY_VIOLATION = 200;

  Tensor *ITensor =
      CurrStruc->CalculateInertiaTensor(StructureOptions::Initial);
  Eigen::Vector3d MomentsOfInertia = ITensor->getEvals();

  outputXYZ(CurrStruc->getHeadAtom(), baseInformation,
            StructureOptions::Initial, "de-flattened structure");

  CPPUNIT_ASSERT(
      fabs(MomentsOfInertia(2) - (MomentsOfInertia(1) + MomentsOfInertia(0))) >
      MINIMUM_PLANARITY_VIOLATION);
}
