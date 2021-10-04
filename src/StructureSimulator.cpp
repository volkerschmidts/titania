
#include <Atom.hpp>
#include <Bond.hpp>
#include <Declarations.hpp>
#include <Molecule.hpp>
#include <Potential.hpp>
#include <Simplex.hpp>
#include <SphericalHarmonics.hpp>
#include <SphericalHarmonicsMatrix.hpp>
#include <Structure.hpp>
#include <StructureSimulator.hpp>
#include <iostream>
#include <mmff94.ff>
#include <omp.h>

unsigned int NumberOfBond_ff    = (sizeof(b_ff) / sizeof(bond_ff));
unsigned int NumberOfAngle_ff   = (sizeof(a_ff) / sizeof(angle_ff));
unsigned int NumberOfSB_ff      = (sizeof(sb_ff) / sizeof(stretch_bend_ff));
unsigned int NumberOfTorsion_ff = (sizeof(t_ff) / sizeof(torsion_ff));
unsigned int NumberOfBCI_ff =
    (sizeof(bci_ff) / sizeof(bond_charge_increment_ff));

StructureSimulator::StructureSimulator()
{
  HeadPotential = NULL;
  CurrStruc     = NULL;
  ListToOpt     = NULL;
  NumberToOpt = NumberOfAtoms = NumberOfBonds = NumberOfAngles =
      NumberOfTorsions = NumberOfPlanarCenters = 0;
  defOpt                                       = StructureOptions::Initial;
  monitoring                                   = false;
}

StructureSimulator::StructureSimulator(
    Structure *s,
    Flags &flags,
    StructureOptions opt) //= StructureOptions::Initial )
{
#pragma omp critical
  *std::cin.tie() << "\t\tLoading MMFF94 force field...\n";
  int bI, bI2, aI, sbI, dI;
  unsigned int ai;
  CurrStruc = s;
  bI = bI2 = aI = sbI = dI = -1;
  this->setupPartialCharge(flags);

  monitoring = false;

  NumberToOpt           = 0;
  NumberOfAtoms         = CurrStruc->getNOA();
  NumberOfBonds         = s->getNumberOfBonds();
  NumberOfAngles        = 0;
  NumberOfTorsions      = 0;
  NumberOfPlanarCenters = 0;
  defOpt                = opt;

  Atom *A1, *A2, *A3, *A4;
  A2 = A3 = A4 = NULL;

  A1 = s->getHeadAtom();
  while (A1)
  {
    if (!A1->getNumberOfHarmonics())
      ++NumberToOpt;
    A1 = A1->getNext();
  }

  ListToOpt = (Atom **) malloc(NumberToOpt * sizeof(Atom *));

  ai = 0;
  A1 = s->getHeadAtom();
  while (A1)
  {
    if (!A1->getNumberOfHarmonics())
    {
      ListToOpt[ai] = A1;
      ++ai;
    }
    A1 = A1->getNext();
  }

  Bond **ListOfBonds = CurrStruc->getListOfBonds();
  HeadPotential      = NULL;
  Potential *CurrP   = HeadPotential;

  for (unsigned int b = 0; b < NumberOfBonds; ++b)
  {
    if (ListOfBonds[b]->getAtom1()->getAtomType() <
        ListOfBonds[b]->getAtom2()->getAtomType())
    {
      A1 = ListOfBonds[b]->getAtom1();
      A2 = ListOfBonds[b]->getAtom2();
    }
    else
    {
      A2 = ListOfBonds[b]->getAtom1();
      A1 = ListOfBonds[b]->getAtom2();
    }
    bI = getStretchIndex(A1, A2);
    if (bI < 0)
    {
      bI = 0;
      if (flags.printWarnings)
      {
#pragma omp critical
        *std::cin.tie()
            << "WARNING:\tFailed to find a proper stretch force field: "
            << A1->getIdentifier() << " (" << A1->getAtomType() << "/"
            << A1->getHybridisation() << ") " << A2->getIdentifier() << " ("
            << A2->getAtomType() << "/" << A2->getHybridisation() << ")"
            << std::endl;
      }
    }
    if (HeadPotential == NULL)
    {
      HeadPotential = new Potential(NULL, PotentialType::Stretch, A1, A2, NULL,
                                    NULL, bI, &b_ff[bI]);
      CurrP         = HeadPotential;
    }
    else
    {
      CurrP = CurrP->addNext(PotentialType::Stretch, A1, A2, NULL, NULL, bI,
                             &b_ff[bI]);
    }
  }

  Potential *P1 = HeadPotential;
  Potential *P2; // = B1->getNext();

  while (P1)
  {
    P2 = P1->getNext();
    while (P2 && P2->getType() == PotentialType::Stretch)
    {
      aI = check4Bend(P1, P2, &A1, &A2, &A3);
      if (aI == -2)
      {
        P2 = P2->getNext();
        continue;
      }
      if (P1->getAtom(0) != A1 || P1->getAtom(1) != A1)
      {
        bI  = P1->getIndex();
        bI2 = P2->getIndex();
      }
      else
      {
        bI = P2->getIndex(), bI2 = P1->getIndex();
      } // Has to be done due to Fijk vs Fkji in sb_ff!

      if (aI == -1)
      {
        aI = 0;
        if (flags.printWarnings)
        {
#pragma omp critical
          *std::cin.tie()
              << "WARNING:\tFailed to find a proper bend force field: "
              << A1->getIdentifier() << " (" << A1->getAtomType() << "/"
              << A1->getHybridisation() << ") " << A2->getIdentifier() << " ("
              << A2->getAtomType() << "/" << A2->getHybridisation() << ") "
              << A3->getIdentifier() << " (" << A3->getAtomType() << "/"
              << A3->getHybridisation() << ") " << std::endl;
        }
      }

      CurrP =
          CurrP->addNext(PotentialType::Bend, A1, A2, A3, NULL, aI, &a_ff[aI]);
      ++NumberOfAngles;

      sbI = getStretchBendIndex(A1, A2, A3);
      if (sbI < 0)
      {
        sbI = 0;
        if (flags.printWarnings)
        {
#pragma omp critical
          *std::cin.tie()
              << "WARNING:\tFailed to find a proper stretch-bend force field: "
              << A1->getIdentifier() << " (" << A1->getAtomType() << "/"
              << A1->getHybridisation() << ") " << A2->getIdentifier() << " ("
              << A2->getAtomType() << "/" << A2->getHybridisation() << ") "
              << A3->getIdentifier() << " (" << A3->getAtomType() << "/"
              << A3->getHybridisation() << ") " << std::endl;
        }
      }
      CurrP = CurrP->addNext(PotentialType::StretchBend, A1, A2, A3, NULL, sbI,
                             &sb_ff[sbI], &b_ff[bI], &b_ff[bI2], &a_ff[aI]);
      P2    = P2->getNext();
    }
    P1 = P1->getNext();
  }

  P1 = HeadPotential;
  P1 = P1->getNext(PotentialType::Bend);

  while (P1)
  {
    P2 = P1->getNext(PotentialType::Bend);
    while (P2)
    {
      dI = check4Dihedral(P1, P2, &A1, &A2, &A3, &A4);
      if (dI == (-2))
      {
        P2 = P2->getNext(PotentialType::Bend);
        continue;
      }
      else if (dI == (-1))
      {
        dI = 0;
        if (flags.printWarnings)
        {
#pragma omp critical
          *std::cin.tie()
              << "WARNING:\tFailed to find a proper dihedral force field: "
              << A1->getIdentifier() << " (" << A1->getAtomType() << "/"
              << A1->getHybridisation() << ") " << A2->getIdentifier() << " ("
              << A2->getAtomType() << "/" << A2->getHybridisation() << ") "
              << A3->getIdentifier() << " (" << A3->getAtomType() << "/"
              << A3->getHybridisation() << ") " << A4->getIdentifier() << " ("
              << A4->getAtomType() << "/" << A4->getHybridisation() << ") "
              << std::endl;
        }
      }
      CurrP = CurrP->addNext(PotentialType::Dihedral, A1, A2, A3, A4, dI,
                             &t_ff[dI]);
      ++NumberOfTorsions;
      P2 = P2->getNext(PotentialType::Bend);
    }
    P1 = P1->getNext(PotentialType::Bend);
  }
  load_lr_Energy(s, &CurrP);

  for (A1 = CurrStruc->getHeadAtom(); A1; A1 = A1->getNext())
  {
    if (check4oop(&A2, &A3, &A4, &A1))
    {
      //         dI = getoopIndex( A2, A3, A4, A1);
      // getoopIndex is not implemented yet. So we can just set dI to 0
      dI = 0;
      CurrP =
          CurrP->addNext(PotentialType::oop, A2, A3, A4, A1, dI, &oop_ff[dI]);
      ++NumberOfPlanarCenters;
    }
  }

  if (flags.use_initial_holonomics)
  {
    for (CurrP = HeadPotential; CurrP; CurrP = CurrP->getNext())
    {
      CurrP->set_initial_holonomics();
    }
  }
  CurrP       = NULL;
  ListOfBonds = NULL;
  A1 = A2 = A3 = A4 = NULL;
}

StructureSimulator::~StructureSimulator()
{
  if (monitoringFile.is_open())
  {
    monitoringFile.close();
    monitoringFile.clear();
  }
  unsigned int a;
  for (a = 0; a < NumberToOpt; ++a)
    ListToOpt[a] = NULL;

  while (HeadPotential->getNext())
  {
    HeadPotential = HeadPotential->getNext();
    delete HeadPotential->getPrev();
  }
  delete HeadPotential;
  CurrStruc = NULL;
  free(ListToOpt);
}

void
StructureSimulator::newStructure(Structure *s, StructureOptions opt)
{
  defOpt = opt;
  Atom *CurrAtom;
  CurrAtom  = s->getHeadAtom();
  CurrStruc = s;
  for (unsigned a = 0; a < NumberToOpt; ++a)
  {
    if (ListToOpt[a])
    {
      atomByrdcIndex(&CurrAtom, ListToOpt[a]->getrdcIndex());
      ListToOpt[a] = CurrAtom;
    }
  }

  Potential *CurrP = HeadPotential;
  while (CurrP)
  {
    CurrP->resetStructure(s);
    CurrP = CurrP->getNext();
  }

  CurrAtom = NULL;
  delete CurrAtom;
}

// making this a bound method of the StructureSimulator class would simplify all
// these lookups but messes up the pointer assignment to the Simplex::vfunc
// below ...
double
totalEfunc(const std::vector<double> &x,
           std::vector<double> &grad,
           void *func_data)
{
  Efunc_par *pars = reinterpret_cast<Efunc_par *>(func_data);

  if (!grad.empty())
  {
    // pass, Nelder-Mead doesn't use gradients
  }

  Eigen::VectorXd point =
      Eigen::Map<const Eigen::VectorXd, Eigen::Unaligned>(x.data(), x.size());
  pars->Simulator->setCoordinates(point, pars->reduced);

  if (pars->Simulator->getmonitoring() && !(pars->count % 100))
  {
    pars->Simulator->printStructure(pars->count);
  }

  return pars->Simulator->totalEnergy();
}

unsigned int
StructureSimulator::startSimplexMinimizer(bool reduced)
{
  unsigned int NMAX = 1000;
  return this->startSimplexMinimizer(NMAX, reduced);
}

unsigned int
StructureSimulator::startSimplexMinimizer(unsigned int steps, bool reduced)
{
  double del  = 0.01;
  double ftol = 5e-5;
  return this->startSimplexMinimizer(steps, del, ftol, reduced);
}

unsigned int
StructureSimulator::startSimplexMinimizer(unsigned int steps,
                                          double del,
                                          double ftol,
                                          bool reduced)
{
  if (monitoring && !monitoringFile.is_open())
    monitoring = false;
  unsigned int ai;
  unsigned int ndim;
  if (!reduced)
    ndim = 3 * NumberOfAtoms;
  else
    ndim = 3 * NumberToOpt;

  Eigen::VectorXd C = Eigen::VectorXd::Zero(ndim);

  Atom *CurrAtom = CurrStruc->getHeadAtom();
  Coordinates *Coord;

  ai = 0;
  while (CurrAtom)
  {
    if (reduced && CurrAtom->getNumberOfHarmonics())
    {
      CurrAtom = CurrAtom->getNext();
      continue;
    }
    Coord         = CurrAtom->getCoordinates(defOpt);
    C(3 * ai)     = Coord->x;
    C(3 * ai + 1) = Coord->y;
    C(3 * ai + 2) = Coord->z;
    ++ai; // TODO(vsz) shouldn't this check if (ai > NumberToOpt)?
    CurrAtom = CurrAtom->getNext();
  }

  Simplex simplex(C);
  simplex.vfunc = totalEfunc;
  Efunc_par func_par;
  func_par.count     = 0;
  func_par.reduced   = reduced;
  func_par.Simulator = this;
  simplex.func_par   = &func_par;

  simplex.setinitialstep(del);
  simplex.setftol(ftol);
  simplex.setNMAX(steps);
  Eigen::VectorXd result; // TODO(vsz) result unused?

  result = simplex.minimize();
  return simplex.getnfunc();
}

void
StructureSimulator::setMonitorFile(std::string name)
{
  if (monitoring == false)
  {
    monitoring = true;
  }
  if (monitoringFile.is_open()) // There is currently a file, close it previous
                                // to opening a new one
  {
    monitoringFile.close();
    monitoringFile.clear();
    *std::cin.tie() << "Information:\tThe structure simulator function closed "
                       "the priviously opened file...\n";
  }
  monitoringFile.open(name, std::ios::binary | std::ios::out | std::ios::trunc);
}

void
StructureSimulator::monitorMinimizer(bool m)
{
  monitoring = m;
  if (monitoring)
    return;

  if (monitoringFile.is_open())
  {
    monitoringFile.close();
    monitoringFile.clear();
  }
}

void
StructureSimulator::setCoordinates(const Eigen::VectorXd &C, const bool reduced)
{
  unsigned int ai;

  Atom *CurrAtom = CurrStruc->getHeadAtom();
  Coordinates *Coord;
  ai = 0;
  while (CurrAtom)
  {
    if (reduced && CurrAtom->getNumberOfHarmonics())
    {
      CurrAtom = CurrAtom->getNext();
      continue;
    }
    Coord    = CurrAtom->getCoordinates(defOpt);
    Coord->x = C(3 * ai);
    Coord->y = C(3 * ai + 1);
    Coord->z = C(3 * ai + 2);
    ++ai;
    CurrAtom = CurrAtom->getNext();
  }
  CurrAtom = NULL;
  Coord    = NULL;
  delete CurrAtom;
  delete Coord;
}

void
StructureSimulator::setupPartialCharge(Flags &flags)
{
  unsigned int b, ci;
  int i1, i2;
  double cpci, pci = .0;

  Atom *CurrAtom, *BondAtom;

  for (CurrAtom = CurrStruc->getHeadAtom(); CurrAtom;
       CurrAtom = CurrAtom->getNext())
  {
    pci = CurrAtom->getPartialCharge();
    for (b = 0; b < CurrAtom->getNumberOfBonds(); ++b)
    {
      BondAtom = CurrAtom->getBond(b)->getBondpartner(CurrAtom);
      if (CurrAtom->getAtomType() <= BondAtom->getAtomType())
      {
        cpci = 1.0;
        i1   = CurrAtom->getAtomType();
        i2   = BondAtom->getAtomType();
      }
      else
      {
        cpci = -1.0;
        i2   = CurrAtom->getAtomType();
        i1   = BondAtom->getAtomType();
      }

      for (ci = 0; ci < NumberOfBCI_ff; ++ci)
      {
        if (i1 == bci_ff[ci].Atom_type_1 && i2 == bci_ff[ci].Atom_type_2)
        {
          cpci *= bci_ff[ci].bci;
          break;
        }
      }
      if (ci == NumberOfBCI_ff)
      {
        cpci = .0;
        if (flags.printWarnings)
          *std::cin.tie() << "WARNING:\tBCI for " << CurrAtom->getIdentifier()
                          << "-" << BondAtom->getIdentifier()
                          << " bond not found\n";
        continue;
      }
      pci += cpci;
    }
    CurrAtom->setPartialCharge(pci);
  }

  BondAtom = NULL;
}

Potential *
StructureSimulator::getPotential(PotentialType t)
{
  if (t == PotentialType::Stretch)
    return HeadPotential;
  else
    return HeadPotential->getNext(t);
}

double
StructureSimulator::totalEnergy()
{
  double Bend, Stretch, StretchBend, Dihedral, VDW, TE;
  Bend = Stretch = StretchBend = Dihedral = VDW = .0;


  //   #pragma omp parallel shared(Stretch,Bend,StretchBend,Dihedral,VDW)
  {
    //      #pragma omp single
    {
#pragma omp task shared(Stretch)
      {
        Stretch = bondEnergy();
      }
#pragma omp task shared(Bend)
      {
        Bend = angleEnergy();
      }
#pragma omp task shared(StretchBend)
      {
        StretchBend = stretch_bend_Energy();
      }
#pragma omp task shared(Dihedral)
      {
        Dihedral = torsionEnergy();
      }
#pragma omp task shared(VDW)
      {
        VDW = lrEnergy();
      }
    }
  }
#pragma omp taskwait
  TE = Stretch + StretchBend + Bend + Dihedral + VDW;

  return TE;
}

double
StructureSimulator::AtomEnergy(Atom *CurrAtom, StructureOptions opt)
{
  double tmp, TE = .0;
  int N            = 0;
  Potential *CurrP = HeadPotential->getNext(PotentialType::Bend);
  while (CurrP)
  {
    if (CurrP->getAtom(0) == CurrAtom || CurrP->getAtom(2) == CurrAtom)
    {
      tmp = CurrP->getEnergy(opt);
      if (std::isnan(tmp))
        tmp = .0;
      TE += tmp;
      ++N;
    }
    CurrP = CurrP->getNext(PotentialType::Bend);
  }
  CurrP = NULL;
  TE    = TE / ((double) N);
  CurrAtom->setEnergy(TE);
  return TE;
}

int
StructureSimulator::optimizeAngles(Atom *CurrAtom)
{
  int optAngle = 0;
  double tmp, TEi, TEo, x, y, z;

  Atom *Partner;
  Coordinates *C, *P;

  StructureOptions opt = StructureOptions::Optimized;

  TEi = TEo = .0;

  Partner = NULL;
  C       = CurrAtom->getCoordinates(opt);

  // Recalculate full Energy and find inversion partner
  Potential *CurrP = HeadPotential->getNext(PotentialType::Bend);
  while (CurrP)
  {
    if (CurrP->getAtom(0) == CurrAtom || CurrP->getAtom(2) == CurrAtom)
    {
      tmp = CurrP->getEnergy(opt);
      if (std::isnan(tmp))
        tmp = .0;
      TEi += tmp;

      if (Partner)
      {
        if (CurrP->getAtom(0) == CurrAtom &&
            CurrP->getAtom(1)->getrdcIndex() < Partner->getrdcIndex())
          Partner = CurrP->getAtom(1);
        else if (CurrP->getAtom(1) == CurrAtom &&
                 CurrP->getAtom(0)->getrdcIndex() < Partner->getrdcIndex())
          Partner = CurrP->getAtom(1);
        else
        {
          CurrP = CurrP->getNext(PotentialType::Bend);
          continue;
        }
      }
      else
      {
        if (CurrP->getAtom(0) == CurrAtom)
          Partner = CurrP->getAtom(1);
        else
          Partner = CurrP->getAtom(0);
      }
    }
    CurrP = CurrP->getNext(PotentialType::Bend);
  }

  // Setup and perform inversion

  P = Partner->getCoordinates(opt);

  // Defintion: C = P + r   <-> r = C - P

  x = C->x - P->x;
  y = C->y - P->y;
  z = C->z - P->z;

  // Invert so C = P - r

  C->x = P->x - x;
  C->y = P->y - y;
  C->z = P->z - z;

  CurrP = HeadPotential->getNext(PotentialType::Bend);
  while (CurrP)
  {
    if (CurrP->getAtom(0) == CurrAtom || CurrP->getAtom(2) == CurrAtom)
    {
      tmp = CurrP->getEnergy(opt);
      if (std::isnan(tmp))
        tmp = .0;
      TEo += tmp;
    }
    CurrP = CurrP->getNext(PotentialType::Bend);
  }

  if (TEo > TEi)
  {
    C->x     = P->x + x;
    C->y     = P->y + y;
    C->z     = P->z + z;
    optAngle = 0;
  }
  P = C   = NULL;
  Partner = NULL;
  delete C;
  delete P;
  delete Partner;
  return optAngle;
}

int
getStretchIndex(Atom *A1, Atom *A2)
{
  int bI = 0;
  Atom *a1, *a2;
  if (A1->getAtomType() < A2->getAtomType())
  {
    a1 = A1;
    a2 = A2;
  }
  else
  {
    a1 = A2;
    a2 = A1;
  }
  for (unsigned int i = 0; i < NumberOfBond_ff; ++i, ++bI)
  {
    if (b_ff[bI].Atom_type_1 != a1->getAtomType())
      continue;
    if (b_ff[bI].Atom_type_2 != a2->getAtomType())
      continue;
    return bI;
  }
  return -1;
}

int
check4Bend(Potential *P1, Potential *P2, Atom **A1, Atom **A2, Atom **A3)
{
  Atom *A11 = P1->getAtom(0);
  Atom *A12 = P1->getAtom(1);
  Atom *A21 = P2->getAtom(0);
  Atom *A22 = P2->getAtom(1);
  int res   = -2;
  if (A11 == A22)
  {
    res   = getBendIndex(A12, A11, A21);
    A1[0] = A12;
    A2[0] = A11;
    A3[0] = A21;
  }
  if (A11 == A21)
  {
    res   = getBendIndex(A12, A11, A22);
    A1[0] = A12;
    A2[0] = A11;
    A3[0] = A22;
  }
  if (A12 == A22)
  {
    res   = getBendIndex(A11, A12, A21);
    A1[0] = A11;
    A2[0] = A12;
    A3[0] = A21;
  }
  if (A12 == A21)
  {
    res   = getBendIndex(A11, A12, A22);
    A1[0] = A11;
    A2[0] = A12;
    A3[0] = A22;
  }
  if (A1[0]->getAtomType() > A3[0]->getAtomType())
  {
    A11   = A1[0];
    A1[0] = A3[0];
    A3[0] = A11;
  }
  A11 = A12 = A21 = A22 = NULL;
  return res;
}

int
getBendIndex(Atom *A1, Atom *A2, Atom *A3)
{
  int bI = 0;
  Atom *a1, *a2, *a3;
  if (A1->getAtomType() < A3->getAtomType())
  {
    a1 = A1;
    a2 = A2;
    a3 = A3;
  }
  else
  {
    a1 = A3;
    a2 = A2;
    a3 = A1;
  }
  for (unsigned int i = 0; i < NumberOfAngle_ff; ++i, ++bI)
  {
    if (a_ff[bI].Atom_type_1 != a1->getAtomType())
      continue;
    if (a_ff[bI].Atom_type_2 != a2->getAtomType())
      continue;
    if (a_ff[bI].Atom_type_3 != a3->getAtomType())
      continue;
    return bI;
  }
  return -1;
}

int
getStretchBendIndex(Atom *A1, Atom *A2, Atom *A3)
{
  int bI = 0;
  Atom *a1, *a2, *a3;
  if (A1->getAtomType() < A3->getAtomType())
  {
    a1 = A1;
    a2 = A2;
    a3 = A3;
  }
  else
  {
    a1 = A3;
    a2 = A2;
    a3 = A1;
  }
  for (unsigned int i = 0; i < NumberOfSB_ff; ++i, ++bI)
  {
    if (sb_ff[bI].Atom_type_1 != a1->getAtomType())
      continue;
    if (sb_ff[bI].Atom_type_2 != a2->getAtomType())
      continue;
    if (sb_ff[bI].Atom_type_3 != a3->getAtomType())
      continue;
    return bI;
  }
  return -1;
}

int
check4Dihedral(Potential *P1,
               Potential *P2,
               Atom **A1,
               Atom **A2,
               Atom **A3,
               Atom **A4)
{
  Atom *A11 = P1->getAtom(0);
  Atom *A12 = P1->getAtom(1);
  Atom *A13 = P1->getAtom(2);
  Atom *A21 = P2->getAtom(0);
  Atom *A22 = P2->getAtom(1);
  Atom *A23 = P2->getAtom(2);
  int res   = -2;
  if ((A11 == A22) && (A12 == A21))
  {
    if (A12->getAtomType() < A11->getAtomType())
    {
      res   = getDihedralIndex(A13, A12, A11, A23);
      A1[0] = A13;
      A2[0] = A12;
      A3[0] = A11;
      A4[0] = A23;
    }
    else
    {
      res   = getDihedralIndex(A23, A11, A12, A13);
      A1[0] = A23;
      A2[0] = A11;
      A3[0] = A12;
      A4[0] = A13;
    }
  }
  else if ((A11 == A22) && (A12 == A23))
  {
    if (A12->getAtomType() < A11->getAtomType())
    {
      res   = getDihedralIndex(A13, A12, A11, A21);
      A1[0] = A13;
      A2[0] = A12;
      A3[0] = A11;
      A4[0] = A21;
    }
    else
    {
      res   = getDihedralIndex(A21, A11, A12, A13);
      A1[0] = A21;
      A2[0] = A11;
      A3[0] = A12;
      A4[0] = A13;
    }
  }
  else if ((A12 == A21) && (A13 == A22))
  {
    if (A12->getAtomType() < A23->getAtomType())
    {
      res   = getDihedralIndex(A11, A12, A13, A23);
      A1[0] = A11;
      A2[0] = A12;
      A3[0] = A13;
      A4[0] = A23;
    }
    else
    {
      res   = getDihedralIndex(A23, A13, A12, A11);
      A1[0] = A23;
      A2[0] = A13;
      A3[0] = A12;
      A4[0] = A11;
    }
  }
  else if ((A12 == A23) && (A13 == A22))
  {
    if (A12->getAtomType() < A13->getAtomType())
    {
      res   = getDihedralIndex(A11, A12, A13, A21);
      A1[0] = A11;
      A2[0] = A12;
      A3[0] = A13;
      A4[0] = A21;
    }
    else
    {
      res   = getDihedralIndex(A21, A13, A12, A11);
      A1[0] = A21;
      A2[0] = A13;
      A3[0] = A12;
      A4[0] = A11;
    }
  }
  A11 = A12 = A13 = A21 = A22 = A23 = NULL;
  return res;
}

int
getDihedralIndex(Atom *a1, Atom *a2, Atom *a3, Atom *a4)
{
  unsigned int t1, t2, t3, t4, bI;
  t1 = t2 = t3 = t4 = 0;
  while (t_ff[t2].Atom_type_2 != a2->getAtomType())
  {
    if (++t2 == NumberOfTorsion_ff)
      return -2;
    if (t_ff[t2].Atom_type_2 > a2->getAtomType())
      return -2;
  }
  while (t_ff[(t3 + t2)].Atom_type_3 != a3->getAtomType())
  {
    if ((++t3 + t2) == NumberOfTorsion_ff)
      return -2;
    if (t_ff[t3 + t2].Atom_type_3 > a3->getAtomType())
      return -2;
  }
  bI = t2 + t3;
  while (t_ff[(t1 + bI)].Atom_type_1 != a1->getAtomType())
  {
    if ((++t1 + bI) > NumberOfTorsion_ff)
      return bI;
    if (t_ff[(t1 + bI)].Atom_type_1 > a1->getAtomType())
    {
      t1 = 0;
      break;
    }
  }
  bI += t1;
  while (t_ff[(t4 + bI)].Atom_type_4 != a4->getAtomType())
  {
    if ((++t4 + bI) > NumberOfTorsion_ff)
      return bI;
    if (t_ff[(t4 + bI)].Atom_type_1 > a1->getAtomType())
    {
      t4 = 0;
      break;
    }
  }
  return (t4 + bI);
}

int
check4oop(Atom **a1, Atom **a2, Atom **a3, Atom **c)
{
  if ((c[0]->hasChiralVolume() &&
       (fabs(c[0]->getChiralVolume(StructureOptions::Initial, true)) > 1e-1)) ||
      c[0]->getNumberOfBonds() != 3)
    return 0;
  a1[0] = c[0]->getBondpartner(0);
  a2[0] = c[0]->getBondpartner(1);
  a3[0] = c[0]->getBondpartner(2);

  return 1;
}

int
getoopIndex(Atom *a1, Atom *a2, Atom *a3, Atom *c)
{
  if (a1 && a2 && a3 && c)
    std::cout << "Information: getoopIndex is not fully implemented\n";
  return 0;
}

double
StructureSimulator::bondEnergy()
{
  Potential *CurrP;
  double stretch = .0;
  for (CurrP = HeadPotential; CurrP;
       CurrP = CurrP->getNext(PotentialType::Stretch))
  {
    stretch += CurrP->getEnergy(defOpt);
  }
  return stretch;
}

double
StructureSimulator::angleEnergy()
{
  Potential *CurrP;
  double bend = .0;
  for (CurrP = HeadPotential->getNext(PotentialType::Bend); CurrP;
       CurrP = CurrP->getNext(PotentialType::Bend))
  {
    {
      bend += CurrP->getEnergy(defOpt);
    }
  }
  return bend;
}

double
StructureSimulator::stretch_bend_Energy()
{
  Potential *CurrP;
  double stretchbend = .0;
  for (CurrP = HeadPotential->getNext(PotentialType::StretchBend); CurrP;
       CurrP = CurrP->getNext(PotentialType::StretchBend))
  {
    {
      stretchbend += CurrP->getEnergy(defOpt);
    }
  }
  return stretchbend;
}

double
StructureSimulator::torsionEnergy()
{
  Potential *CurrP;
  double dihedral = .0;
  for (CurrP = HeadPotential->getNext(PotentialType::Dihedral); CurrP;
       CurrP = CurrP->getNext(PotentialType::Dihedral))
  {
    {
      dihedral += CurrP->getEnergy(defOpt);
    }
  }
  return dihedral;
}

double
StructureSimulator::lrEnergy()
{
  Potential *CurrP;
  double vdw, vdwSum = .0;
  for (CurrP = HeadPotential->getNext(PotentialType::vanderWaals); CurrP;
       CurrP = CurrP->getNext(PotentialType::vanderWaals))
  {
    vdw = CurrP->getEnergy(defOpt);
    vdwSum += vdw;
  }
  return vdwSum;
}

void
load_lr_Energy(Structure *s, Potential **CurrP)
{
  unsigned int s1, s2, s3;
  double int14;
  bool defined;

  Atom *A = s->getHeadAtom();
  Atom *CurrAtom;
  Atom *A2, *A3, *A4;
  while (A)
  {
    CurrAtom = A->getNext();
    while (CurrAtom)
    {
      int14   = 1.0;
      defined = false;
      for (s1 = 0; s1 < A->getNumberOfBonds(); ++s1)
      {
        A2 = A->getBond(s1)->getBondpartner(A);
        if (A2 == CurrAtom)
        {
          defined = true;
          break;
        }
        for (s2 = 0; s2 < A2->getNumberOfBonds(); ++s2)
        {
          A3 = A2->getBond(s2)->getBondpartner(A2);
          if (A3 == CurrAtom)
          {
            defined = true;
            break;
          }
          else if (A3 == A)
            continue;
          for (s3 = 0; s3 < A3->getNumberOfBonds(); ++s3)
          {
            A4 = A3->getBond(s3)->getBondpartner(A3);
            if (A4 == CurrAtom)
            {
              int14 = 0.75;
              break;
            }
          }
        }
      }
      if (defined)
      {
        CurrAtom = CurrAtom->getNext();
        continue;
      }
      else
      {
        CurrP[0] = CurrP[0]->addNext(
            PotentialType::vanderWaals, A, CurrAtom, NULL, NULL, 0,
            &vdW_ff[A->getAtomType() - 1], &vdW_ff[CurrAtom->getAtomType() - 1],
            &int14, NULL);
      }
      CurrAtom = CurrAtom->getNext();
    }
    A = A->getNext();
  }
  CurrAtom = NULL;
  A2       = NULL;
  A3       = NULL;

  delete CurrAtom;
  delete A2;
  delete A3;
}

void
StructureSimulator::printStructure(unsigned int step)
{
  Atom *CurrAtom = CurrStruc->getHeadAtom();
  Coordinates *C;
  monitoringFile << NumberOfAtoms << std::endl
                 << CurrStruc->getLabel() << "_" << step << std::endl;
  while (CurrAtom)
  {
    C = CurrAtom->getCoordinates(defOpt);
    if (CurrAtom->getLabel().back() == 'a')
      monitoringFile << "T  " << C->x << " " << C->y << " " << C->z
                     << std::endl;
    else
      monitoringFile << CurrAtom->getElement() << "  " << C->x << " " << C->y
                     << " " << C->z << std::endl;

    CurrAtom = CurrAtom->getNext();
  }
  CurrAtom = NULL;
  C        = NULL;
  delete CurrAtom;
  delete C;
}

void
removePlanarity(Structure *CurrStruc,
                StructureSimulator &MainSimulation,
                BasicInformation &baseInformation)
{
  if (!((baseInformation.secondaryInput & StructureInputType::planarX) |
        (baseInformation.secondaryInput & StructureInputType::planarY) |
        (baseInformation.secondaryInput & StructureInputType::planarZ)))
    return;

  *std::cin.tie() << "Rearraning Atoms due to planar input\n";
  Atom *CurrAtom;
  unsigned int NOA = CurrStruc->getNOA();
  CurrAtom         = CurrStruc->getHeadAtom();
  while (CurrAtom)
  {
    if (baseInformation.secondaryInput & StructureInputType::planarX)
      CurrAtom->getCoordinates(StructureOptions::Initial)->x =
          (double) (100 - rand() % 200) / 10.0;
    else if (baseInformation.secondaryInput & StructureInputType::planarY)
      CurrAtom->getCoordinates(StructureOptions::Initial)->y =
          (double) (100 - rand() % 200) / 10.0;
    else if (baseInformation.secondaryInput & StructureInputType::planarZ)
      CurrAtom->getCoordinates(StructureOptions::Initial)->z =
          (double) (100 - rand() % 200) / 10.0;

    CurrAtom = CurrAtom->getNext();
  }

  MainSimulation.setStanStructureOptions(StructureOptions::Initial);
  MainSimulation.monitorMinimizer(true);
  MainSimulation.setMonitorFile((baseInformation.outputFileName + ".p23"));
  MainSimulation.startSimplexMinimizer(10 * NOA, 0.5, 1e-1, false);
  MainSimulation.startSimplexMinimizer(100 * NOA, 0.5, 1e-2, false);
  MainSimulation.startSimplexMinimizer(1000 * NOA, 0.5, 1e-3, false);
  MainSimulation.startSimplexMinimizer(10000 * NOA, 0.5, 1e-4, false);
  *std::cin.tie() << std::endl;
  MainSimulation.monitorMinimizer(false);
  for (SphericalHarmonics *CurrHarmonic =
           CurrStruc->getHeadYmatrix()->getHeadHarmonic();
       CurrHarmonic; CurrHarmonic = CurrHarmonic->getNext())
  {
    CurrHarmonic->recalculateAngles(StructureOptions::Initial);
  }
}

int
maxEnergyAtom(Atom *CurrAtom, StructureSimulator &MainSimulation)
{
  int optAngle = 1;
  double maxE, currE;
  maxE = currE = .0;
  Atom *maxA   = NULL;
  while (CurrAtom)
  {
    if (CurrAtom->getNumberOfHarmonics())
      currE = CurrAtom->getEnergy();
    else
      currE = .0;
    if (currE > maxE)
    {
      maxA = CurrAtom;
      maxE = currE;
    }
    CurrAtom = CurrAtom->getNext();
  }
  if (maxA)
  {
    optAngle = MainSimulation.optimizeAngles(maxA);
  }
  maxA = NULL;
  return optAngle;
}
