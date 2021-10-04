
#include <Atom.hpp>
#include <Bond.hpp>
#include <Declarations.hpp>
#include <LinAlg.hpp>
#include <Molecule.hpp>
#include <Potential.hpp>
#include <Properties.hpp>
#include <RDCdata.hpp>
#include <SphericalHarmonics.hpp>
#include <Structure.hpp>
#include <StructureSimulator.hpp>
#include <eigen3/Eigen/Geometry>
#include <iostream>

Atom::Atom()
{
  element           = "";
  label             = "";
  chiralVolumeAtoms = "";

  z                   = 0;
  A                   = 0;
  AtomType            = 0;
  distance_violations = 0;

  inputIndex        = 0;
  rdcIndex          = 0;
  NumberOfBonds     = 0;
  NumberOfBondRDCs  = 0;
  NumberOfHarmonics = 0;

  gamma               = .0;
  mass                = .0;
  partialCharge       = .0;
  Energy              = .0;
  initialChiralVolume = optimizedChiralVolume = .0;
  vdW_Potential                               = .0;

  fixed = invertable = increased_violations = false;

  nucleus        = new struct NucleusList;
  coordinates    = new struct Coordinates;
  optimized      = new struct Coordinates;
  coordinates->x = coordinates->y = coordinates->z = .0;
  optimized->x = optimized->y = optimized->z = .0;

  hybr            = Hybridisation::undefined;
  chiral_index[0] = chiral_index[1] = chiral_index[2] = 0;

  parent    = NULL;
  prevAtom  = NULL;
  nextAtom  = NULL;
  harmonics = NULL;
  bonds     = NULL;
}

Atom::Atom(Atom *a)
{
  element           = "";
  label             = "";
  chiralVolumeAtoms = "";

  z                   = 0;
  A                   = 0;
  AtomType            = 0;
  distance_violations = 0;

  inputIndex        = 0;
  rdcIndex          = 0;
  NumberOfBonds     = 0;
  NumberOfBondRDCs  = 0;
  NumberOfHarmonics = 0;

  mass                = .0;
  gamma               = .0;
  partialCharge       = .0;
  Energy              = .0;
  initialChiralVolume = optimizedChiralVolume = .0;
  vdW_Potential                               = .0;

  fixed = invertable = increased_violations = false;

  nucleus        = new struct NucleusList;
  coordinates    = new struct Coordinates;
  optimized      = new struct Coordinates;
  coordinates->x = coordinates->y = coordinates->z = .0;
  optimized->x = optimized->y = optimized->z = .0;
  hybr                                       = Hybridisation::undefined;
  chiral_index[0] = chiral_index[1] = chiral_index[2] = 0;

  parent    = a->parent;
  prevAtom  = a;
  nextAtom  = NULL;
  harmonics = NULL;
  bonds     = NULL;
}

Atom::~Atom()
{
  parent   = NULL;
  prevAtom = NULL;
  nextAtom = NULL;
  nucleus  = NULL;

  for (unsigned int i = 0; (bonds && i < NumberOfBonds); ++i)
    bonds[i] = NULL;
  for (unsigned int i = 0; (harmonics && i < NumberOfHarmonics); ++i)
    harmonics[i] = NULL;

  free(bonds);
  free(harmonics);
  delete coordinates;
  delete optimized;
}

void
Atom::loadBondorder()
{
  unsigned int b;
  int db = 0, at = 0;
  bool conjSet   = false;
  bool chinonSys = false;
  Atom *BondAtom;
  Bond *CurrBond;
  Molecule *CurrMol = this->getParent()->getParent();

  if (hybr == Hybridisation::undefined)
    at = 0;

  if (hybr == Hybridisation::sp3 || hybr == Hybridisation::s)
  {
    switch (z)
    {
      case 1:
        at = 5;
        break; // H
      case 6:
        at = 1;
        break; // C
      case 7:
        at = 8;
        break; // N
      case 8:
        at = 6;
        break; // O
      case 16:
        at = 15;
        break; // S
      default:
        std::cerr << "ERROR:\tNo force field defined for "
                  << this->getIdentifier()
                  << "...\n\t\tContact the developer or add an own one...\n";
        break;
    }
  }
  for (b = 0; b < NumberOfBonds; ++b)
  {
    CurrBond = bonds[b];
    if (!CurrBond)
    {
      std::cerr << "ERROR:\tMisslinking in bonds of " << this->getIdentifier()
                << "...\n";
      return;
    }
    BondAtom = CurrBond->getBondpartner(this);
    if (CurrBond->getOrder())
    {
      if (hybr == Hybridisation::sp2 && BondAtom->hybr == Hybridisation::sp2)
        ++db;
      continue;
    }

    if ((hybr == Hybridisation::s ||
         hybr ==
             Hybridisation::sp3)) // && ( BondAtom->hybr == Hybridisation::s ||
                                  // BondAtom->hybr == Hybridisation::sp3 ) )
    {
      CurrBond->setOrder(1);
      continue;
    }
    else if (hybr == Hybridisation::sp2 && NumberOfBonds == 1)
    {
      if (z == 8)
      {
        AtomType           = 7;
        BondAtom->AtomType = 3;
      }
      CurrBond->setOrder(2);
      CurrMol->raiseNumberOfDB();
      continue;
    }
    else if (hybr == Hybridisation::sp1)
    {
      AtomType = 4;
      if (BondAtom->hybr == Hybridisation::sp1)
      {
        if (z == 6)
        {
          BondAtom->AtomType = 4;
        }
        CurrBond->setOrder(3);
        CurrMol->raiseNumberOfDB();
        CurrMol->raiseNumberOfDB();
      }
      else if (BondAtom->hybr == Hybridisation::sp2)
      {
        if (z == 6 && BondAtom->z == 6)
        {
          BondAtom->AtomType = 2;
        }
        else if (z == 6 && BondAtom->z == 8)
        {
          BondAtom->AtomType = 7;
        }
        CurrBond->setOrder(2);
        CurrMol->raiseNumberOfDB();
      }
    }
    else if (hybr == Hybridisation::sp2 && BondAtom->z == 8)
      AtomType = 3;
    else if (hybr == Hybridisation::sp2 && BondAtom->hybr == Hybridisation::sp2)
      ++db;
  }
  for (b = 0; db >= 2 && b < NumberOfBonds; ++b)
  {
    CurrBond = bonds[b];
    if (CurrBond->getBondpartner(this)->z == 8 ||
        CurrBond->getBondpartner(this)->AtomType == 3)
    {
      chinonSys = true;
      break;
    }
  }
  for (b = 0; db >= 2 && b < NumberOfBonds; ++b)
  {
    CurrBond = bonds[b];
    if (CurrBond->getBondpartner(this)->z == 8)
      CurrBond->getBondpartner(this)->AtomType = 7;
    else if (chinonSys && CurrBond->getBondpartner(this)->z == 6 &&
             (CurrBond->getBondpartner(this)->AtomType == 37 ||
              CurrBond->getBondpartner(this)->AtomType == 0))
      CurrBond->getBondpartner(this)->AtomType = 2;
  }


  for (b = 0; b < NumberOfBonds; ++b)
  {
    CurrBond = bonds[b];
    BondAtom = CurrBond->getBondpartner(this);
    if (CurrBond->getOrder())
      continue;
    else if (db == 2 && !chinonSys)
    {
      if (CurrBond->getOrder() == 2)
      {
        AtomType = BondAtom->AtomType = 37;
        conjSet                       = true;
        continue;
      }
    }

    if (db == 1 && BondAtom->hybr == Hybridisation::sp2 &&
        !CurrBond->getOrder())
    {
      if (z == 6 && BondAtom->z == 6)
      {
        AtomType           = 2;
        BondAtom->AtomType = 2;
      }
      else if (z == 6 && BondAtom->z == 8)
      {
        AtomType           = 3;
        BondAtom->AtomType = 7;
      }
      else if (BondAtom->z == 6 && z == 8)
      {
        AtomType           = 7;
        BondAtom->AtomType = 3;
      }
      CurrBond->setOrder(2);
      CurrMol->raiseNumberOfDB();
      continue;
    }
    else if (db == 1)
    {
      CurrBond->setOrder(1);
      if (BondAtom->z == 1)
        BondAtom->AtomType = 5;
      else if (BondAtom->z == 6)
        BondAtom->AtomType = 1;
      continue;
    }
    else if (db == 2 && BondAtom->hybr == Hybridisation::sp2)
    {
      if (!conjSet && !chinonSys)
      {
        AtomType = BondAtom->AtomType = 37, CurrBond->setOrder(2);
        conjSet                       = true;
        CurrMol->raiseNumberOfDB();
      }
      else
      {
        CurrBond->setOrder(1);
      }
    }
  }

  if (!AtomType)
    AtomType = at;
  CurrMol  = NULL;
  BondAtom = NULL;
  CurrBond = NULL;
  delete CurrMol;
  delete BondAtom;
  delete CurrBond;
}

void
Atom::loadAtomType(Flags &flags)
{
  unsigned int b;
  Atom *BondAtom;

  if (z == 1)
  {
    switch (bonds[0]->getBondpartner(this)->AtomType)
    {
      case 1:
        AtomType = 5;
        return;
      case 2:
        AtomType = 5;
        return;
      case 3:
        AtomType = 5;
        return;
      case 4:
        AtomType = 28;
        return;
      case 6: {
        BondAtom = bonds[0]->getBondpartner(this);
        for (b = 0; b < BondAtom->getNumberOfBonds(); ++b)
        {
          if (BondAtom->getBondpartner(b) != this)
          {
            if (BondAtom->getBondpartner(b)->getAtomType() == 37)
              AtomType = 29;
            else
              AtomType = 21;
            BondAtom = NULL;
            return;
          }
        }
      }
      break;
      case 8:
        AtomType = 23;
        return;
      case 10:
        AtomType = 28;
        return;
      case 20:
        AtomType = 5;
        return;
      case 22:
        AtomType = 5;
        return;
      case 30:
        AtomType = 5;
        return;
      case 34:
        AtomType = 36;
        return;
      case 37:
        AtomType = 5;
        return;
      case 43:
        AtomType = 28;
        return;
      case 56:
        AtomType = 36;
        return;
      default:
        std::cerr << "ERROR:\tNo atom type found for " << this->getIdentifier()
                  << "...\n\t\tAdd own force field and atom type or contact "
                     "developer to add one if geometry wont converge...\n";
        return;
    }
  }
  else if (z == 6)
  {
    if (AtomType == 37) // aromatic
    {
      for (b = 0; b < NumberOfBonds; ++b)
      {
        BondAtom = bonds[b]->getBondpartner(this);
        if (BondAtom->z == 1)
        {
          BondAtom->AtomType = 5;
          continue;
        }
        else if (BondAtom->z == 6 || BondAtom->z == 7 || BondAtom->z == 8)
        {
          if (BondAtom->hybr == Hybridisation::sp2)
          {
            if (BondAtom->AtomType == 2)
              BondAtom->AtomType = 37; // DB
            else if (BondAtom->AtomType == 3)
            {
              std::cout << "WARNING:\tAtom type of "
                        << BondAtom->getIdentifier() << " is changed from "
                        << BondAtom->AtomType << " to " << 37 << std::endl;
              BondAtom->AtomType = 37;
            } // Carbonyl
            else if (BondAtom->AtomType != 37 && BondAtom->z == 6)
            {
              if (flags.printWarnings)
                std::cout << "Warning:\tUnknown type of "
                          << BondAtom->getIdentifier() << " current type ("
                          << BondAtom->AtomType << ")" << std::endl;
              BondAtom->AtomType = 2;
            }
            continue;
          }
          else if (BondAtom->hybr == Hybridisation::sp3)
          {
            if (BondAtom->AtomType == 1)
              BondAtom->AtomType = 1;
            else if (BondAtom->AtomType == 8)
            {
              BondAtom->AtomType = 10;
              BondAtom->hybr     = Hybridisation::sp2;
            }
          }
        }
      }
    }
    else if (AtomType == 3)
    {
      for (b = 0; b < NumberOfBonds; ++b)
      {
        BondAtom = bonds[b]->getBondpartner(this);
        if (BondAtom->z == 6 && BondAtom->NumberOfBonds == 3)
        {
          AtomType = 37;
          if (BondAtom->AtomType == 2)
            BondAtom->AtomType = 37;
        }
        else if (BondAtom->z == 6 || BondAtom->z == 7 || BondAtom->z == 8)
        {
          if (BondAtom->AtomType == 8)
          {
            BondAtom->AtomType = 10;
            BondAtom->hybr     = Hybridisation::sp2;
          }
        }
      }
    }
  }
  else if (z == 7)
  {
    for (b = 0; b < NumberOfBonds; ++b)
    {
      if (hybr == Hybridisation::sp3 &&
          bonds[b]->getBondpartner(this)->getHybridisation() ==
              Hybridisation::sp2)
      {
        hybr = Hybridisation::sp2;
        break;
      }
    }
    if (hybr == Hybridisation::sp3)
    {
      if (NumberOfBonds)
        AtomType = 34;
    }
    else if (hybr == Hybridisation::sp2)
    {
      for (b = 0; b < NumberOfBonds; ++b)
      {
        BondAtom = bonds[b]->getBondpartner(this);
        if (BondAtom->getZ() == 6 && BondAtom->AtomType == 0)
        {
          bonds[b]->getBondpartner(this)->AtomType = 57;
          AtomType                                 = 56;
          break;
        }
        else if (BondAtom->getZ() == 1)
          BondAtom->AtomType = 28;
        else if (BondAtom->AtomType == 57)
        {
          AtomType = 56;
          break;
        }
      }
      for (b = 0; b < NumberOfBonds; ++b)
      {
        BondAtom = bonds[b]->getBondpartner(this);
        if (BondAtom->getZ() == 1 && AtomType == 56)
        {
          bonds[b]->getBondpartner(this)->AtomType = 36;
        }
      }
    }
  }
  else if (z == 8)
  {
    if (hybr == Hybridisation::sp2)
    {
      if (AtomType == 32)
        ;
      else if (NumberOfBonds == 2)
        AtomType = 6;
      else
        AtomType = 7;
    }
  }
  else if (z == 9)
  {
    AtomType = 11;
  }
  else if (z == 16)
  {
    int NOO = 0;
    for (b = 0; b < NumberOfBonds; ++b)
    {
      BondAtom = bonds[b]->getBondpartner(this);
      if (BondAtom->z == 6 || BondAtom->z == 7)
        ++NOO;
    }
    if (NOO >= 2)
      AtomType = 18;
    else if (NOO == 1)
      AtomType = 17;
    else
      AtomType = 15;

    for (b = 0; b < NumberOfBonds; ++b)
    {
      if (AtomType < 17)
        break;
      BondAtom = bonds[b]->getBondpartner(this);
      if (BondAtom->z == 8)
      {
        BondAtom->AtomType = 32;
      }
      else if (BondAtom->z == 7)
      {
        BondAtom->AtomType = 43;
      }
    }
  }
  else if (z == 17)
  {
    AtomType = 12;
  }
  else if (z == 35)
  {
    AtomType = 13;
  }
  else if (z == 53)
  {
    AtomType = 14;
  }
  BondAtom = NULL;
  delete BondAtom;
}

/*****************************************************
 *                      setOC                        *
 *                                                   *
 *     Sets the optimized coordinates of a atom.     *
 *                                                   *
 *****************************************************/

void
Atom::setOC(double X, double Y, double Z)
{
  optimized->x = X;
  optimized->y = Y;
  optimized->z = Z;
}

/*****************************************************
 *                      setIC                        *
 *                                                   *
 *      Sets the initial coordinates of a atom.      *
 *                                                   *
 *****************************************************/

void
Atom::setIC(double X, double Y, double Z)
{
  coordinates->x = X;
  coordinates->y = Y;
  coordinates->z = Z;
}

/*****************************************************
 *                  setCoordinates                   *
 *                                                   *
 *      Public (overloaded) functions to set the     *
 *   Coordinates of the atoms. These functions call  *
 *       each other in some kind of hirachy:         *
 *                                                   *
 * 1. Eigen type calls                               *
 * 2. doubel type                                    *
 * 3. private specialized 'setIC' or 'setOC'         *
 *                                                   *
 *****************************************************/

void
Atom::setCoordinates(double nX,
                     double nY,
                     double nZ,
                     enum StructureOptions options)
{
  if (!((options & StructureOptions::Initial) ||
        (options & StructureOptions::Optimized)))
  {
    std::cerr << "ERROR:\nUnknown type of coordinates were "
                 "requested...\n\t\tReturning initial coordinates...\n";
    options = StructureOptions::Initial;
  }

  if (options & StructureOptions::Initial)
    this->setIC(nX, nY, nZ);
  else if (options & StructureOptions::Optimized)
    this->setOC(nX, nY, nZ);
}

void
Atom::setCoordinates(Eigen::Vector3d r, enum StructureOptions options)
{
  this->setCoordinates(r(STANDARD_AXIS_X_), r(STANDARD_AXIS_Y_),
                       r(STANDARD_AXIS_Z_), options);
}

void
Atom::setCoordinates(double *r, StructureOptions options)
{
  this->setCoordinates(r[0], r[1], r[2], options);
}

/*****************************************************
 *                     setGamma                      *
 *                                                   *
 *   In the list with nuclear information only the   *
 *   magnetic moment and nuclear spin is saved. By   *
 *   this one has to calculate gamma from it. This   *
 *   function normally should just be called once.   *
 *                                                   *
 *****************************************************/

void
Atom::setGamma()
{
  if (nucleus->magneticmoment == 0)
    gamma = 0.0;
  else
    gamma = nucleus->magneticmoment * MU_N_ / H_BAR_ / nucleus->nuclearspin;
}

/*****************************************************
 *                     setNext                       *
 *                                                   *
 * Initializes and links a standard atom template to *
 *         the tail atom. Be carefull, setNext       *
 *    automatically links the current Atom to the    *
 * TailAtom of the parent Structure to prevent false *
 *                    linking!                       *
 *                                                   *
 *****************************************************/

Atom *
Atom::setNext()
{
  Atom *tmp                       = new Atom(this);
  parent->getTailAtom()->nextAtom = tmp;
  tmp->prevAtom                   = parent->getTailAtom();
  parent->setTailAtom(tmp);
  return tmp;
}

/*****************************************************
 *                 initializeBonds                   *
 *                                                   *
 * Initializes and clears the bonds pointer of this. *
 *                                                   *
 *****************************************************/

void
Atom::initializeBonds()
{
  if (NumberOfBonds && !bonds)
  {
    bonds = (Bond **) malloc(NumberOfBonds * sizeof(Bond *));
  }
  else
    return;
  for (unsigned int i = 0; i < NumberOfBonds; ++i)
  {
    bonds[i] = NULL;
  }
}

/*****************************************************
 *                      getBond                      *
 *                                                   *
 *   Returns the bond between this and a requested   *
 *              (possible) partner atom.             *
 *                                                   *
 *****************************************************/

Bond *
Atom::getBond(Atom *a) const
{
  unsigned int i;
  for (i = 0; i < NumberOfBonds; ++i)
  {
    if (bonds[i]->getBondpartner(this) == a)
      return bonds[i];
  }
  return NULL;
}

/*****************************************************
 *                  getBondpartner                   *
 *                                                   *
 *  Returns the bond partner for this at bond with   *
 *           the respective bond index.              *
 *                                                   *
 *****************************************************/


Atom *
Atom::getBondpartner(unsigned int i) const
{
  if (i < NumberOfBonds)
  {
    return bonds[i]->getBondpartner(this);
  }
  else
    return NULL;
}


/*****************************************************
 *                  addBondpartner                   *
 *                                                   *
 *  Initializes a new bond, links this and atom a to *
 *   it and does the same thing directly to atom a!  *
 *        By this the linking is accelerated.        *
 *                                                   *
 *****************************************************/

void
Atom::addBondpartner(Atom *a)
{
  unsigned int i;
  for (i = 0; i < NumberOfBonds; ++i)
  {
    if (bonds[i] == NULL)
    {
      bonds[i] = new Bond(this, a);
      parent->NBPP();
      a->addBond(bonds[i]);
      return;
    }
    else if (bonds[i]->getBondpartner(this) == a)
    {
      return;
    }
  }
  std::cerr << "ERROR:\tTried to add a fifth or higher bond partner ("
            << a->element + a->label << ") to atom " << element + label
            << " which should have " << NumberOfBonds << " partners...\n";
}

/*****************************************************
 *                     addBond                       *
 *                                                   *
 * Initializes a new bond. This is something like an *
 *      acceleration function for addBondpartner.    *
 *                                                   *
 *****************************************************/

void
Atom::addBond(Bond *b)
{
  unsigned int i;
  for (i = 0; i < NumberOfBonds; ++i)
  {
    if (bonds[i] == b)
    {
      return;
    }
    else if (bonds[i] == NULL)
    {
      bonds[i] = b;
      return;
    }
  }
  std::cerr << "ERROR:\tTried to add a fifth or higher bond partner ("
            << b->getBondpartner(this)->getIdentifier() << ") to atom "
            << this->getIdentifier() << " which should have " << NumberOfBonds
            << " partners...\n";
}

void
Atom::set_chiral_volume_indizes()
{
  unsigned int i, j;
  int k, NOHs;

  for (i = j = 0, NOHs = k = 0; i < NumberOfBonds; ++i, ++k)
  {
    if (bonds[i]->getBondpartner(this)->getZ() == 1)
      ++NOHs;
    if (bonds[i]->getHarmonic())
    {
      chiral_index[j++] = k;
      if (j == CHIRAL_CENTER_ATOMS_)
        return;
    }
  }
  if (NOHs > 2)
  {
    invertable = true;
    return;
  }
  for (i = 0, k = 0; i < NumberOfBonds; ++i, ++k)
  {
    if (bonds[i]->getHarmonic())
    {
      continue;
    }
    else
    {
      chiral_index[j++] = k;
      if (j == CHIRAL_CENTER_ATOMS_)
        return;
    }
  }
}

double
Atom::getChiralVolume(StructureOptions opt, bool norm)
{
  if (NumberOfBonds < CHIRAL_CENTER_ATOMS_)
    return .0;
  else if (invertable)
    return .0;
  else if (chiral_index[0] == chiral_index[1] ||
           chiral_index[1] == chiral_index[2] ||
           chiral_index[0] == chiral_index[2])
    set_chiral_volume_indizes();
  else if (opt == StructureOptions::Initial && initialChiralVolume)
    return initialChiralVolume;


  double Vchiral = .0;

  Atom *A, *B, *C;
  A = bonds[chiral_index[0]]->getBondpartner(this);
  B = bonds[chiral_index[1]]->getBondpartner(this);
  C = bonds[chiral_index[2]]->getBondpartner(this);

  Eigen::Vector3d a = A->Coordinates2Eigen(opt) - this->Coordinates2Eigen(opt);
  Eigen::Vector3d b = B->Coordinates2Eigen(opt) - this->Coordinates2Eigen(opt);
  Eigen::Vector3d c = C->Coordinates2Eigen(opt) - this->Coordinates2Eigen(opt);
  if (norm)
  {
    a = a / a.norm();
    b = b / b.norm();
    c = c / c.norm();
  }
  Eigen::Vector3d bxc = b.cross(c);

  Vchiral = (a.transpose() * bxc);

  A = B = C = NULL;

  if (opt == StructureOptions::Initial && initialChiralVolume)
    initialChiralVolume = Vchiral;
  else if (opt == StructureOptions::Optimized && optimizedChiralVolume)
    optimizedChiralVolume = Vchiral;

  return Vchiral;
}

std::string
Atom::printChiralVolumeAtoms()
{
  if (NumberOfBonds < CHIRAL_CENTER_ATOMS_)
    return "";
  else if (invertable)
    return "";
  else if (chiral_index[0] == chiral_index[1] ||
           chiral_index[1] == chiral_index[2] ||
           chiral_index[0] == chiral_index[2])
    set_chiral_volume_indizes();
  if (chiralVolumeAtoms == "")
  {
    for (int i = 0; i < CHIRAL_CENTER_ATOMS_; ++i)
    {
      chiralVolumeAtoms +=
          ((i == 0 ? "" : "-") +
           bonds[chiral_index[i]]->getBondpartner(this)->getIdentifier());
    }
  }
  return chiralVolumeAtoms;
}

const char *
Atom::printChiralVolumeAtoms(int length)
{
  if (NumberOfBonds < CHIRAL_CENTER_ATOMS_)
    return "";
  else if (invertable)
    return "";
  else if (chiral_index[0] == chiral_index[1] ||
           chiral_index[1] == chiral_index[2] ||
           chiral_index[0] == chiral_index[2])
    set_chiral_volume_indizes();
  char *buf    = (char *) malloc(STANDARD_BUFFER_SIZE * sizeof(char));
  char *retBuf = (char *) malloc(STANDARD_BUFFER_SIZE * sizeof(char));
  for (int i = 0; i < STANDARD_BUFFER_SIZE; ++i)
    buf[i] = 0;
  for (int i = 0; i < CHIRAL_CENTER_ATOMS_; ++i)
  {
    snprintf(
        retBuf, STANDARD_BUFFER_SIZE * sizeof(char), "%s %*s", buf, length,
        bonds[chiral_index[i]]->getBondpartner(this)->getIdentifier().c_str());
    strcpy(buf, retBuf);
  }
  free(buf);
  return retBuf;
}


double
Atom::check_Chiral_Volume_validity(StructureOptions opt)
{
  if (NumberOfBonds <= 3)
    return .0;
  Atom **chiral_list = (Atom **) malloc(3 * sizeof(Atom *));
  unsigned int kick, atom, index, permutation;
  Eigen::Vector3d a, b, c, bxc;
  double Vchiral, mean_V = .0, stddev_V = .0;
  for (kick = 0; kick < 4; ++kick)
  {
    for (atom = index = 0; atom < NumberOfBonds; ++atom)
    {
      if (atom == kick)
        continue;
      chiral_list[index] = getBondpartner(atom);
      ++index;
    }
    for (permutation = 1; permutation <= 6; ++permutation)
    {
      a = chiral_list[get_S3_permuted_index(permutation, 0)]->Coordinates2Eigen(
              opt) -
          this->Coordinates2Eigen(opt);
      b = chiral_list[get_S3_permuted_index(permutation, 1)]->Coordinates2Eigen(
              opt) -
          this->Coordinates2Eigen(opt);
      c = chiral_list[get_S3_permuted_index(permutation, 2)]->Coordinates2Eigen(
              opt) -
          this->Coordinates2Eigen(opt);
      bxc = b.cross(c);

      Vchiral = a.transpose() * bxc;
      Vchiral /= (a.norm() * bxc.norm());
      mean_V += fabs(Vchiral);
    }
  }
  mean_V /= 24.0;
  for (kick = 0; kick < 4; ++kick)
  {
    for (atom = index = 0; atom < NumberOfBonds; ++atom)
    {
      if (atom == kick)
        continue;
      chiral_list[index] = getBondpartner(atom);
      ++index;
    }
    for (permutation = 1; permutation <= 6; ++permutation)
    {
      a = chiral_list[get_S3_permuted_index(permutation, 0)]->Coordinates2Eigen(
              opt) -
          this->Coordinates2Eigen(opt);
      b = chiral_list[get_S3_permuted_index(permutation, 1)]->Coordinates2Eigen(
              opt) -
          this->Coordinates2Eigen(opt);
      c = chiral_list[get_S3_permuted_index(permutation, 2)]->Coordinates2Eigen(
              opt) -
          this->Coordinates2Eigen(opt);
      bxc = b.cross(c);

      Vchiral = a.transpose() * bxc;
      Vchiral /= (a.norm() * bxc.norm());
      stddev_V += pow(mean_V - fabs(Vchiral), 2);
    }
  }
  stddev_V = sqrt(stddev_V / 24.0);
  for (index = 0; index < 3; ++index)
    chiral_list[index] = NULL;
  return stddev_V;
}


/*****************************************************
 *               initializeHarmonics                 *
 *                                                   *
 * Allocate the full harmonics pointer and reset all *
 *               subpointers to NULL.                *
 *                                                   *
 *  Don't ask for the clearHarmonics function! Even  *
 *     though the pointers were just cleared, the    *
 *    program ran into memory error when accessing   *
 *                  the pointers!!!                  *
 *                                                   *
 *****************************************************/

void
Atom::initializeHarmonics()
{
  if (NumberOfHarmonics)
  {
    harmonics = (SphericalHarmonics **) malloc(NumberOfHarmonics *
                                               sizeof(SphericalHarmonics *));
    for (unsigned int i = 0; i < NumberOfHarmonics; ++i)
      harmonics[i] = NULL;
  }
  this->clearHarmonics();
}

void
Atom::clearHarmonics()
{
  unsigned int i;
  if (NumberOfBonds)
  {
    for (i = 0; i < NumberOfBonds; ++i)
    {
      bonds[i]->setHarmonic(NULL);
    }
  }
  if (NumberOfHarmonics)
  {
    for (i = 0; i < NumberOfHarmonics; ++i)
    {
      harmonics[i] = NULL;
    }
  }
}

/*****************************************************
 *                  addBondharmonic                  *
 *                                                   *
 *        Link a harmonic to this and to bond!       *
 *                                                   *
 *****************************************************/

void
Atom::addBondharmonic(SphericalHarmonics *h)
{
  unsigned int i;
  if (h->getRDC()->getRange() == 1)
  {
    for (i = 0; i < NumberOfBonds; ++i)
    {
      if (bonds[i]->getHarmonic() == h)
      {
        break;
      }
      else if (bonds[i]->getBondpartner(this) == h->getAtom1() ||
               bonds[i]->getBondpartner(this) == h->getAtom2())
      {
        bonds[i]->setHarmonic(h);
      }
    }
  }
  for (i = 0; i < NumberOfHarmonics; ++i)
  {
    if (harmonics[i] == h)
    {
      return;
    }
    else if (harmonics[i] == NULL)
    {
      harmonics[i] = h;
      return;
    }
  }
}

/*****************************************************
 *                  getBondharmonic                  *
 *                                                   *
 *  Returns the harmonic associated with bond i, and *
 *   only if there is an harmonic (and hence rdc)!   *
 *                                                   *
 *****************************************************/

SphericalHarmonics *
Atom::getBondharmonic(unsigned int i) const
{
  if (i < NumberOfBonds)
  {
    return bonds[i]->getHarmonic();
  }
  else
    return NULL;
}

/*****************************************************
 *                    getHarmonic                    *
 *                                                   *
 *   Returns the harmonic i, independent of bonds!   *
 *                                                   *
 *****************************************************/

SphericalHarmonics *
Atom::getHarmonic(unsigned int i) const
{
  if (i < NumberOfHarmonics)
  {
    return harmonics[i];
  }
  else
    return NULL;
}

SphericalHarmonics *
Atom::getHarmonic(Atom *partner) const
{
  for (unsigned int i = 0; i < NumberOfHarmonics; ++i)
  {
    if ((harmonics[i]->getAtom1() == partner &&
         harmonics[i]->getAtom2() == this) ||
        (harmonics[i]->getAtom2() == partner &&
         harmonics[i]->getAtom1() == this))
      return harmonics[i];
  }
  return NULL;
}


Atom *
Atom::initializeOptimizedAtoms(Structure *preStruc, Structure *p)
{
  unsigned int i;
  Atom *CurrAtom = preStruc->getHeadAtom();
  Atom *preAtom;
  Atom *tmpAtom;

  parent     = p;
  inputIndex = CurrAtom->inputIndex;
  rdcIndex   = CurrAtom->rdcIndex;
  AtomType   = CurrAtom->AtomType;
  if (CurrAtom->increased_violations)
    distance_violations = CurrAtom->distance_violations;
  element           = CurrAtom->element;
  label             = CurrAtom->label;
  chiralVolumeAtoms = CurrAtom->chiralVolumeAtoms;

  findNucleusByName(*this);

  hybr            = CurrAtom->hybr;
  chiral_index[0] = CurrAtom->chiral_index[0];
  chiral_index[1] = CurrAtom->chiral_index[1];
  chiral_index[2] = CurrAtom->chiral_index[2];

  invertable = CurrAtom->invertable;

  coordinates->x = CurrAtom->optimized->x;
  coordinates->y = CurrAtom->optimized->y;
  coordinates->z = CurrAtom->optimized->z;

  gamma                 = CurrAtom->gamma;
  mass                  = CurrAtom->mass;
  partialCharge         = CurrAtom->partialCharge;
  Energy                = CurrAtom->Energy;
  initialChiralVolume   = CurrAtom->optimizedChiralVolume;
  optimizedChiralVolume = .0;

  NumberOfHarmonics = CurrAtom->NumberOfHarmonics;
  NumberOfBonds     = CurrAtom->NumberOfBonds;

  bonds = (Bond **) malloc(NumberOfBonds * sizeof(Bond *));
  for (i = 0; i < NumberOfBonds; ++i)
    bonds[i] = NULL;

  harmonics = (SphericalHarmonics **) malloc(NumberOfHarmonics *
                                             sizeof(SphericalHarmonics *));
  for (i = 0; i < NumberOfHarmonics; ++i)
    harmonics[i] = NULL;

  CurrAtom            = CurrAtom->getNext();
  Atom *newAtom       = new Atom(this);
  nextAtom            = newAtom;
  newAtom->inputIndex = CurrAtom->inputIndex;
  newAtom->rdcIndex   = CurrAtom->rdcIndex;
  newAtom->AtomType   = CurrAtom->AtomType;
  if (CurrAtom->increased_violations)
    newAtom->distance_violations = CurrAtom->distance_violations;
  newAtom->element           = CurrAtom->element;
  newAtom->label             = CurrAtom->label;
  newAtom->chiralVolumeAtoms = CurrAtom->chiralVolumeAtoms;


  findNucleusByName(*newAtom);

  newAtom->hybr            = CurrAtom->hybr;
  newAtom->chiral_index[0] = CurrAtom->chiral_index[0];
  newAtom->chiral_index[1] = CurrAtom->chiral_index[1];
  newAtom->chiral_index[2] = CurrAtom->chiral_index[2];

  newAtom->invertable = CurrAtom->invertable;

  newAtom->coordinates->x = CurrAtom->optimized->x;
  newAtom->coordinates->y = CurrAtom->optimized->y;
  newAtom->coordinates->z = CurrAtom->optimized->z;

  newAtom->gamma               = CurrAtom->gamma;
  newAtom->mass                = CurrAtom->mass;
  newAtom->partialCharge       = CurrAtom->partialCharge;
  newAtom->Energy              = CurrAtom->Energy;
  newAtom->initialChiralVolume = CurrAtom->optimizedChiralVolume;

  newAtom->NumberOfHarmonics = CurrAtom->NumberOfHarmonics;
  newAtom->NumberOfBonds     = CurrAtom->NumberOfBonds;

  newAtom->bonds = (Bond **) malloc(newAtom->NumberOfBonds * sizeof(Bond *));
  for (i = 0; i < newAtom->NumberOfBonds; ++i)
    newAtom->bonds[i] = NULL;

  newAtom->harmonics = (SphericalHarmonics **) malloc(
      newAtom->NumberOfHarmonics * sizeof(SphericalHarmonics *));
  for (i = 0; i < newAtom->NumberOfHarmonics; ++i)
    newAtom->harmonics[i] = NULL;

  while (CurrAtom->getNext())
  {
    CurrAtom                    = CurrAtom->getNext();
    newAtom->nextAtom           = new Atom(newAtom);
    newAtom->nextAtom->prevAtom = newAtom;
    newAtom                     = newAtom->nextAtom;

    newAtom->inputIndex = CurrAtom->inputIndex;
    newAtom->rdcIndex   = CurrAtom->rdcIndex;
    newAtom->AtomType   = CurrAtom->AtomType;
    if (CurrAtom->increased_violations)
      newAtom->distance_violations = CurrAtom->distance_violations;
    newAtom->element           = CurrAtom->element;
    newAtom->label             = CurrAtom->label;
    newAtom->chiralVolumeAtoms = CurrAtom->chiralVolumeAtoms;

    findNucleusByName(*newAtom);

    newAtom->hybr            = CurrAtom->hybr;
    newAtom->chiral_index[0] = CurrAtom->chiral_index[0];
    newAtom->chiral_index[1] = CurrAtom->chiral_index[1];
    newAtom->chiral_index[2] = CurrAtom->chiral_index[2];

    newAtom->invertable = CurrAtom->invertable;

    newAtom->coordinates->x = CurrAtom->optimized->x;
    newAtom->coordinates->y = CurrAtom->optimized->y;
    newAtom->coordinates->z = CurrAtom->optimized->z;

    newAtom->gamma               = CurrAtom->gamma;
    newAtom->mass                = CurrAtom->mass;
    newAtom->partialCharge       = CurrAtom->partialCharge;
    newAtom->Energy              = CurrAtom->Energy;
    newAtom->initialChiralVolume = CurrAtom->optimizedChiralVolume;

    newAtom->NumberOfHarmonics = CurrAtom->NumberOfHarmonics;
    newAtom->NumberOfBonds     = CurrAtom->NumberOfBonds;

    newAtom->bonds = (Bond **) malloc(newAtom->NumberOfBonds * sizeof(Bond *));
    for (i = 0; i < newAtom->NumberOfBonds; ++i)
      newAtom->bonds[i] = NULL;

    newAtom->harmonics = (SphericalHarmonics **) malloc(
        newAtom->NumberOfHarmonics * sizeof(SphericalHarmonics *));
    for (i = 0; i < newAtom->NumberOfHarmonics; ++i)
      newAtom->harmonics[i] = NULL;
  }

  /* Add Bondpartner copy */
  CurrAtom = p->getHeadAtom();
  preAtom  = preStruc->getHeadAtom();

  while (CurrAtom)
  {
    for (i = 0; i < preAtom->NumberOfBonds; ++i)
    {
      if (preAtom->getBondpartner(i))
      {
        tmpAtom = p->getHeadAtom();
        atomByIdentifier(preAtom->getBondpartner(i)->getIdentifier(), &tmpAtom);
        CurrAtom->addBondpartner(tmpAtom);
        tmpAtom->addBondpartner(CurrAtom);
        CurrAtom->getBond(tmpAtom)->setLength(preAtom->getBond(i)->getLength());
      }
      else
        break;
    }
    CurrAtom = CurrAtom->getNext();
    preAtom  = preAtom->getNext();
  }

  tmpAtom  = NULL;
  CurrAtom = NULL;
  preAtom  = NULL;
  return newAtom;
}

/*****************************************************
 *                  getCoordinates                   *
 *                                                   *
 *           Returns the Coordinates of this.        *
 *                                                   *
 *****************************************************/

struct Coordinates *
Atom::getCoordinates(enum StructureOptions options) const
{
  if (options & StructureOptions::Initial)
    return coordinates;
  else if (options & StructureOptions::Optimized)
    return optimized;
  else
  {
    std::cerr << "ERROR:\nUnknown type of coordinates were "
                 "requested...\n\t\tReturning initial coordinates...\n";
    return coordinates;
  }
}

/*****************************************************
 *                       getIC                       *
 *                                                   *
 *  Translates the initial Coordinates to a vector   *
 *   representation. This is a private sub routine   *
 *  and should acutally never be called besides the  *
 *     faster coordinates access in more general     *
 *                     functions.                    *
 *                                                   *
 *****************************************************/

Eigen::Vector3d
Atom::getIC()
{
  Eigen::Vector3d initialCoordinatesVector;

  initialCoordinatesVector(STANDARD_AXIS_X_) = coordinates->x;
  initialCoordinatesVector(STANDARD_AXIS_Y_) = coordinates->y;
  initialCoordinatesVector(STANDARD_AXIS_Z_) = coordinates->z;

  return initialCoordinatesVector;
}

/*****************************************************
 *                       getOC                       *
 *                                                   *
 * Translates the optimized Coordinates to a vector  *
 *   representation. This is a private sub routine   *
 *  and should acutally never be called besides the  *
 *     faster coordinates access in more general     *
 *                     functions.                    *
 *                                                   *
 *****************************************************/

Eigen::Vector3d
Atom::getOC()
{
  Eigen::Vector3d optimizedCoordinatesVector;

  optimizedCoordinatesVector(STANDARD_AXIS_X_) = optimized->x;
  optimizedCoordinatesVector(STANDARD_AXIS_Y_) = optimized->y;
  optimizedCoordinatesVector(STANDARD_AXIS_Z_) = optimized->z;

  return optimizedCoordinatesVector;
}

/*****************************************************
 *                  Coordinates2Eigen                *
 *                                                   *
 * Translates the optimized Coordinates to a vector  *
 *                   representation.                 *
 *                                                   *
 *****************************************************/

Eigen::Vector3d
Atom::Coordinates2Eigen(enum StructureOptions options)
{
  Eigen::Vector3d tmp = Eigen::Vector3d::Zero();
  if (!((options & StructureOptions::Initial) ||
        (options & StructureOptions::Optimized)))
  {
    std::cerr << "ERROR:\tUnknown type of coordinates were "
                 "requested...\n\t\tReturning initial Coordinates...\n";
    options = StructureOptions::Initial;
  }

  if (options & StructureOptions::Initial)
    return this->getIC();
  else if (options & StructureOptions::Optimized)
    return this->getOC();
  return tmp;
}

/*****************************************************
 *                   shiftCoordinates                *
 *                                                   *
 *   Shifts the coordinates of an atom by a defined  *
 *  vector V. This is done by summation! Be careful  *
 *                    about that!                    *
 *                                                   *
 *****************************************************/

void
Atom::shiftCoordinates(Eigen::Vector3d &shiftVector, StructureOptions opt)
{
  Eigen::Vector3d shiftedCoordinates = Coordinates2Eigen(opt) + shiftVector;
  setCoordinates(shiftedCoordinates, opt);
}

/*****************************************************
 *                  rotateCoordinates                *
 *                                                   *
 *  Rotates the coordinates of an atom by a defined  *
 * transformation matrix R. This will also allow non *
 *   linear transformations! Be careful about that!  *
 *                                                   *
 *****************************************************/

void
Atom::rotateCoordinates(Eigen::Matrix3d &rotationMatrix, StructureOptions opt)
{
  Eigen::Vector3d rotatedCoordinates = rotationMatrix * Coordinates2Eigen(opt);
  setCoordinates(rotatedCoordinates, opt);
}

/*****************************************************
 *                  retainCoordinates                *
 *                                                   *
 *   Copies the initial coodinates to the optimized  *
 *    ones. This is actually usefull when making a   *
 *      function more useable by StructureOptions    *
 *                      arguments.                   *
 *                                                   *
 *****************************************************/

void
Atom::retainCoordinates()
{
  optimized->x = coordinates->x;
  optimized->y = coordinates->y;
  optimized->z = coordinates->z;
}

/*****************************************************
 *                     getDistance                   *
 *                                                   *
 *     Return the distance between this and a in     *
 *  Angstrom. This function uses the Eigen function  *
 *   norm, which returns a double representing the   *
 *                 length of a vector.               *
 *                                                   *
 *****************************************************/

double
Atom::getDistance(Atom *referenceAtom, StructureOptions opt)
{
  return ((referenceAtom->Coordinates2Eigen(opt) - this->Coordinates2Eigen(opt))
              .norm());
}

double
Atom::getAngle(Atom *A1, Atom *A2, StructureOptions opt)
{
  Eigen::Vector3d V1 =
      A1->Coordinates2Eigen(opt) - this->Coordinates2Eigen(opt);
  Eigen::Vector3d V2 = A1->Coordinates2Eigen(opt) - A2->Coordinates2Eigen(opt);
  double V1_n        = V1.norm();
  double V2_n        = V2.norm();
  if (V1_n > .0)
    V1 /= (V1.norm());
  if (V2_n > .0)
    V2 /= (V2.norm());

  return (acos(V1.transpose() * V2));
}


double
Atom::calculateVdW_Potential(StructureSimulator &MainSim, StructureOptions opt)
{
  Potential *CurrP = MainSim.getPotential(PotentialType::Stretch);
  CurrP            = CurrP->getNext(this, PotentialType::vanderWaals);
  double P         = .0;
  while (CurrP)
  {
    P += CurrP->getEnergy(opt);
    CurrP = CurrP->getNext(this, PotentialType::vanderWaals);
  }

  if (opt & StructureOptions::Optimized)
  {
    vdW_Potential = P;
  }
  return P;
}

double
Atom::getVdW_Potential()
{
  return vdW_Potential;
}

/*****************************************************
 *                  atomByIdentifier                 *
 *                                                   *
 *    Shifts the adress in A to the Atom with the    *
 *  Identifier (Element + Label) requested of to a   *
 *                    NULL pointer!                  *
 *                                                   *
 *****************************************************/

void
atomByIdentifier(std::string ident, Atom **atom)
{
  atom[0] = atom[0]->getParent()->getHeadAtom();
  while (atom[0] && ident.compare(atom[0]->getElement() + atom[0]->getLabel()))
  {
    atom[0] = atom[0]->getNext();
  }
}

/*****************************************************
 *                    atomByIndex                    *
 *                                                   *
 *    Shifts the adress in A to the Atom with the    *
 *     rdc index requested of to a NULL pointer!     *
 *                                                   *
 *****************************************************/

void
atomByrdcIndex(Atom **atom, unsigned int index)
{
  atom[0] = atom[0]->getParent()->getHeadAtom();
  while (atom[0] && atom[0]->getrdcIndex() != index)
  {
    atom[0] = atom[0]->getNext();
  }
}
