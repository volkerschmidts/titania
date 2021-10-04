
#include <Atom.hpp>
#include <Bond.hpp>
#include <Declarations.hpp>
#include <IPC>
#include <Molecule.hpp>
#include <Parser/Parser.hpp>
#include <StanIPC.h>
#include <Structure.hpp>
#include <iostream>
#include <string>

Molecule *
ipcStan()
{
  std::string titania_unit_path = TITANIA_UNIT_TESTS_DIR__;
  BasicInformation baseInformation;
  Flags flags;
  std::fstream input;
  baseInformation.limits.phobos_opts_eckart =
      (double *) malloc(4 * sizeof(double));
  baseInformation.limits.phobos_opts_sphericals =
      (double *) malloc(4 * sizeof(double));
  baseInformation.inputFileName =
      (titania_unit_path + "/systems/cppunit_input.tna");
  input.open(baseInformation.inputFileName, std::ios::binary | std::ios::in);
  InputFile reader(input);
  reader.checkFile();
  Molecule *M = new Molecule("");
  M->loadInput(reader, baseInformation, flags);
  input.close(), input.clear();

  setupBonds(*M, flags);
  M->determineMolecularMass();
  M->getHeadStruc()->initializeRDCindex(*M);
  initializeSphericalHarmonics(*M, M->getHeadStruc());
  free(baseInformation.limits.phobos_opts_eckart);
  free(baseInformation.limits.phobos_opts_sphericals);
  return M;
}

unsigned int
ipc_stan_bonds()
{
  return NUMBOND_TITANIA_TEST_IPC__;
}

unsigned int
ipc_stan_angles()
{
  return NUMANGL_TITANIA_TEST_IPC__;
}

unsigned int
ipc_stan_torsions()
{
  return NUMTORS_TITANIA_TEST_IPC__;
}

unsigned int
ipc_stan_rdcs()
{
  return NUMRDCS_TITANIA_TEST_IPC__;
}
void
outputStanMol(Molecule *M)
{
  std::cout << "Label: " << M->getLabel() << std::endl;
  std::cout << "Atoms: " << M->getNOA() << std::endl;

  Structure *S = M->getHeadStruc();
  std::cout << "Struc: " << S->getLabel() << std::endl;

  Atom *A = S->getHeadAtom();
  while (A)
  {
    std::cout << A->getIndex() << "\t" << A->getElement() << "\t"
              << A->getLabel() << "\t" << A->getA() << "\t" << A->getZ() << "\t"
              << A->getCoordinates(StructureOptions::Initial)->x << "\t"
              << A->getCoordinates(StructureOptions::Initial)->y << "\t"
              << A->getCoordinates(StructureOptions::Initial)->z << std::endl;
    A = A->getNext();
  }
  A = NULL;
  S = NULL;
  M = NULL;
}
