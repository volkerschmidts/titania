
#include <Atom.hpp>
#include <Bond.hpp>
#include <Declarations.hpp>
#include <Molecule.hpp>
#include <Properties.hpp>
#include <RDCdata.hpp>
#include <RDCset.hpp>
#include <Structure.hpp>
#include <iostream>


Bond::Bond()
{
  order    = 0;
  length   = 0;
  atom1    = NULL;
  atom2    = NULL;
  harmonic = NULL;
}

Bond::Bond(Atom *a1, Atom *a2, SphericalHarmonics *Y)
{
  order    = 0;
  length   = 0;
  atom1    = a1;
  atom2    = a2;
  harmonic = Y;

  Eigen::Vector3d r = a1->Coordinates2Eigen(StructureOptions::Initial) -
                      a2->Coordinates2Eigen(StructureOptions::Initial);
  length = r.norm();
}

Bond::~Bond()
{
  atom1    = NULL;
  atom2    = NULL;
  harmonic = NULL;
}

void
Bond::copyOrder()
{
  Atom *iniA1 = atom1->getParent()->getParent()->getHeadStruc()->getHeadAtom();
  Atom *iniA2 = atom2->getParent()->getParent()->getHeadStruc()->getHeadAtom();
  atomByrdcIndex(&iniA1, atom1->getrdcIndex());
  atomByrdcIndex(&iniA2, atom2->getrdcIndex());

  if (iniA1->getBond(iniA2)->getOrder() == iniA2->getBond(iniA1)->getOrder())
    order = iniA1->getBond(iniA2)->getOrder();
  else
    std::cerr << "Something went wrong in copyOrder of "
              << atom1->getIdentifier() << " " << atom2->getIdentifier()
              << std::endl;
  iniA1 = iniA2 = NULL;
  delete iniA1;
  delete iniA2;
}

/* End of Bond implementation */

void
setupBonds(Molecule &CurrMol, Flags &flags)
{
  unsigned int NB    = CurrMol.getHeadStruc()->getNumberOfBonds();
  Bond **ListOfBonds = CurrMol.getHeadStruc()->initializeListOfBonds();
  unsigned int b, bondI;
  Atom *CurrAtom;
  Structure *CurrStruc;
  RDCset *CurrSet;
  RDCdata *CurrData;
  CurrStruc = CurrMol.getHeadStruc();
  CurrAtom  = CurrStruc->getHeadAtom();
  while (CurrAtom)
  {
    for (bondI = 0; bondI < CurrAtom->getNumberOfBonds(); ++bondI)
    {
      for (b = 0; b < NB; ++b)
      {
        if (CurrAtom->getBond(bondI) == ListOfBonds[b])
          break;
        else if (ListOfBonds[b] == NULL)
        {
          ListOfBonds[b] = CurrAtom->getBond(bondI);
          break;
        }
      }
    }
    CurrAtom = CurrAtom->getNext();
  }

  CurrSet = CurrMol.getHeadSet();

  while (CurrSet)
  {
    CurrData = CurrSet->getHeadData();
    while (CurrData)
    {
      CurrData->determineRange(CurrMol);
      CurrData->determineEffectiveWeight(flags);
      CurrData = CurrData->getNext();
    }
    CurrSet = CurrSet->getNext();
  }
  CurrMol.initializeWmatrix();

  CurrAtom = CurrStruc->getHeadAtom();
  while (CurrAtom)
  {
    loadHybridisation(CurrAtom);
    CurrAtom = CurrAtom->getNext();
  }
  CurrAtom = CurrStruc->getHeadAtom();
  while (CurrAtom)
  {
    CurrAtom->loadBondorder();
    CurrAtom = CurrAtom->getNext();
  }
  CurrAtom = CurrStruc->getHeadAtom();
  while (CurrAtom)
  {
    CurrAtom->loadAtomType(flags);
    CurrAtom = CurrAtom->getNext();
  }

  CurrStruc = NULL;
  CurrAtom  = NULL;
  CurrSet   = NULL;
  CurrData  = NULL;
}
