
#include <Atom.hpp>
#include <Declarations.hpp>
#include <Molecule.hpp>
#include <RDCdata.hpp>
#include <RDCset.hpp>
#include <Structure.hpp>

RDCset::RDCset(RDCset *s)
{
  label                       = "";
  inputIndex                  = 0;
  HeadData                    = s->HeadData;
  TailData                    = NULL;
  prevSet                     = s;
  nextSet                     = NULL;
  parent                      = s->parent;
  NumberOfRDCs                = 0;
  SECONDA_gap_1_5_sensitivity = SECONDA_gap_5_6_sensitivity = .0;
  SECONDA_sensitivity_rank                                  = 0;
}

RDCset::RDCset()
{
  label                       = "";
  inputIndex                  = 0;
  HeadData                    = NULL;
  TailData                    = NULL;
  prevSet                     = NULL;
  nextSet                     = NULL;
  parent                      = NULL;
  NumberOfRDCs                = 0;
  SECONDA_gap_1_5_sensitivity = SECONDA_gap_5_6_sensitivity = .0;
  SECONDA_sensitivity_rank                                  = 0;
}

RDCset::RDCset(Molecule &M, Molecule &P, Flags &flags)
{
  HeadData = TailData = NULL;
  nextSet = prevSet = NULL;
  RDCset *pre       = P.getHeadSet();
  RDCset *tmp       = this;

  RDCdata *pRDC, *tRDC;

  Atom *A1, *A2;

  M.setHeadSet(this);
  while (pre)
  {
    tRDC            = tmp->addHeadData();
    pRDC            = pre->getHeadData();
    tmp->label      = pre->label;
    tmp->inputIndex = pre->inputIndex;
    parent          = &M;
    NumberOfRDCs    = pre->NumberOfRDCs;

    while (pRDC)
    {
      A2 = A1 = M.getTailStruc()->getHeadAtom();
      atomByrdcIndex(&A1, pRDC->getAtom1()->getrdcIndex());
      atomByrdcIndex(&A2, pRDC->getAtom2()->getrdcIndex());
      tRDC->setAtom1(A1);
      tRDC->setAtom2(A2);
      tRDC->setD(pRDC->getD());
      tRDC->setD_scaled(pRDC->getD_scaled());
      tRDC->setD_norm(pRDC->getD_norm());
      tRDC->setDeltaD(pRDC->getDeltaD());
      tRDC->setWeight(pRDC->getWeight());
      tRDC->setRange(pRDC->getRange());
      tRDC->setIndex(pRDC->getInputIndex());
      tRDC->setKappa();
      tRDC->determineEffectiveWeight(flags);
      if (pRDC->getNext())
        tRDC = tRDC->appendData();
      pRDC = pRDC->getNext();
    }

    if (pre->getNext())
      tmp = tmp->appendSet();
    pre = pre->getNext();
  }
  SECONDA_gap_1_5_sensitivity = SECONDA_gap_5_6_sensitivity = .0;
  SECONDA_sensitivity_rank                                  = 0;
}

RDCset::~RDCset()
{
  while (HeadData->getNext())
  {
    HeadData = HeadData->getNext();
    delete HeadData->getPrev();
  }
  TailData = NULL;
  delete HeadData;

  prevSet = NULL;
  nextSet = NULL;
  parent  = NULL;
}

std::string
RDCset::getLabel(int sub, bool extend)
{
  if (sub > 0 && sub < static_cast<int>(label.size()))
  {
    std::string result = label.substr(0, sub);
    if (extend && sub < static_cast<int>(label.size() - 3))
      result = label.substr(0, sub - 3) + "...";
    return result;
  }
  else
    return label;
}


RDCdata *
RDCset::addHeadData()
{
  RDCdata *tmp = new RDCdata();
  HeadData = TailData = tmp;
  tmp->setParent(this);
  return tmp;
}

RDCdata *
RDCdata::copyRDC(RDCdata *prev)
{
  RDCdata *copy    = new RDCdata();
  copy->inputIndex = inputIndex;
  copy->rdcAtom1   = rdcAtom1;
  copy->rdcAtom2   = rdcAtom2;
  copy->kappa      = kappa;
  copy->parent     = parent;
  if (prev)
    prev->nextData = copy;
  copy->prevData = prev;
  return copy;
}

RDCset *
RDCset::appendSet()
{
  RDCset *tmp = new RDCset(this);
  parent->setTailSet(tmp);
  nextSet = tmp;
  return tmp;
}

void
RDCset::copyRDCs(RDCset *s)
{
  RDCdata *CurrRDC, *sRDC, *nRDC;
  if (TailData)
  {
    CurrRDC = TailData;
    while (CurrRDC)
    {
      if (CurrRDC->getPrev())
      {
        CurrRDC = CurrRDC->getPrev();
        delete CurrRDC->getNext();
      }
      else
        delete CurrRDC;
    }
    HeadData = TailData = NULL;
  }
  sRDC     = s->getHeadData();
  CurrRDC  = sRDC->copyRDC(NULL);
  HeadData = CurrRDC;
  sRDC     = sRDC->getNext();
  while (sRDC)
  {
    nRDC = sRDC->copyRDC(CurrRDC);
    CurrRDC->setNext(nRDC);
    CurrRDC = nRDC;
    sRDC    = sRDC->getNext();
  }
  TailData = CurrRDC;
  CurrRDC  = NULL;
  sRDC     = NULL;
  nRDC     = NULL;
  delete CurrRDC;
  delete sRDC;
  delete nRDC;
}

double
RDCset::get_sigma_square()
{
  double sigmasq = .0;

  for (RDCdata *RDC = HeadData; RDC; RDC = RDC->getNext())
  {
    sigmasq += (RDC->getDeltaD() * RDC->getDeltaD());
  }
  return sigmasq;
}

Eigen::MatrixXd
RDCset::getWeighting()
{
  unsigned int i;
  RDCdata *CurrData;

  weighting = Eigen::MatrixXd::Identity(NumberOfRDCs, NumberOfRDCs);

  for (i = 0, CurrData = HeadData; ((i < NumberOfRDCs) && CurrData);
       ++i, CurrData   = CurrData->getNext())
  {
    weighting(i, i) = CurrData->getEffectiveWeight();
  }
  return weighting;
}
