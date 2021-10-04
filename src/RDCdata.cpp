
#include <Atom.hpp>
#include <Declarations.hpp>
#include <Helper.hpp>
#include <LinAlg.hpp>
#include <Molecule.hpp>
#include <RDCdata.hpp>
#include <RDCset.hpp>
#include <Structure.hpp>
#include <iomanip>
#include <iostream>


RDCdata::RDCdata(RDCdata *d)
{
  inputIndex = (d->inputIndex + 1);
  rdcAtom1   = NULL;
  rdcAtom2   = NULL;
  D          = .0;
  D_norm     = .0;
  D_scaled   = .0;
  deltaD     = .0;
  weight     = .0;
  effWeight  = .0;
  kappa      = .0;
  prevData   = d;
  nextData   = NULL;
  prevSet = nextSet = NULL;
  parent            = d->parent;
  defined           = true;
  range             = 1;
}

RDCdata::RDCdata()
{
  inputIndex = 0;
  rdcAtom1   = NULL;
  rdcAtom2   = NULL;
  D          = .0;
  D_norm     = .0;
  D_scaled   = .0;
  deltaD     = .0;
  weight     = .0;
  effWeight  = .0;
  kappa      = .0;
  prevData = nextData = NULL;
  prevSet = nextSet = NULL;
  parent            = NULL;
  defined           = true;
  range             = 1;
}

RDCdata::~RDCdata()
{
  rdcAtom1 = NULL;
  rdcAtom2 = NULL;
  prevData = NULL;
  nextData = NULL;
  prevSet  = NULL;
  nextSet  = NULL;
  parent   = NULL;
}

void
RDCdata::determineEffectiveWeight(Flags &flags)
{
  if (!defined)
    weight = effWeight = .0;
  else
  {
    effWeight = weight;
    if (flags.errorWeightInSVD && defined)
      effWeight /= deltaD;
  }
}

Eigen::MatrixXd
RDCdata::getVectorWeights()
{
  unsigned int i, NOS = parent->getParent()->getNORsets();
  RDCdata *CurrData;
  weighting = Eigen::MatrixXd::Identity(NOS, NOS);
  for (i = 0, CurrData = this; ((i < NOS) && CurrData);
       ++i, CurrData   = CurrData->getNextSetData())
  {
    weighting(i, i) = CurrData->getEffectiveWeight();
  }
  return weighting;
}

void
RDCdata::setKappa()
{
  kappa = -(rdcAtom1->getGamma() * rdcAtom2->getGamma() * MU_0_ * H_BAR_) /
          (8.0 * PI_ * PI_);
}

RDCdata *
RDCdata::appendData()
{
  RDCdata *tmp = new RDCdata(this);
  parent->setTailData(tmp);
  nextData = tmp;
  return tmp;
}

double
RDCdata::getDistance(enum StructureOptions options) const
{
  double r = .0;
  Atom *a1, *a2;
  if (range == 2)
  {
    a1 = rdcAtom1->getParent()->getParent()->getHeadStruc()->getAtomByIndex(
        rdcAtom1->getIndex());
    a2 = rdcAtom2->getParent()->getParent()->getHeadStruc()->getAtomByIndex(
        rdcAtom2->getIndex());
    r = a1->getDistance(a2, StructureOptions::Initial);
  }
  else
  {
    a1 = rdcAtom1->getParent()->getParent()->getTailStruc()->getAtomByIndex(
        rdcAtom1->getIndex());
    a2 = rdcAtom2->getParent()->getParent()->getTailStruc()->getAtomByIndex(
        rdcAtom2->getIndex());
    r = a1->getDistance(a2, options);
  }
  a1 = a2 = NULL;
  delete a1;
  delete a2;
  return (r * 1e-10);
}

void
RDCdata::determineRange(Molecule &CurrMol)
{
  Atom *List[CurrMol.getNOA()];
  unsigned int i, j, k;
  unsigned int Range[CurrMol.getNOA()];
  range = 0;
  for (i = 0; i < CurrMol.getNOA(); ++i)
  {
    List[i]  = NULL;
    Range[i] = CurrMol.getNOA(); // > max possible range for this molecule
  }

  List[0] = rdcAtom1;
  for (j = 0, i = 1; j < rdcAtom1->getNumberOfBonds(); ++j, ++i)
  {
    List[i]  = rdcAtom1->getBondpartner(j);
    Range[i] = 1;
    if (rdcAtom1->getBondpartner(j) == rdcAtom2)
    {
      range = 1;
    }
  }

  i = 1;
  while (List[i] && !range && i < CurrMol.getNOA())
  {
    if (List[i] == rdcAtom2)
    {
      range = Range[i];
      break;
    }
    else
    {
      for (j = 0; j < List[i]->getNumberOfBonds(); ++j)
      {
        k = 1;
        if (List[i]->getBondpartner(j) == rdcAtom1)
        {
          continue;
        }
        while (List[k])
        {
          if (List[k] == List[i]->getBondpartner(j))
          {
            if (Range[k] <= (Range[i] + 1))
              break;
            else if (Range[k] > (Range[i] + 1))
            {
              Range[k] = Range[i] + 1;
              break;
            }
          }
          ++k;
        }
        List[k]  = List[i]->getBondpartner(j);
        Range[k] = Range[i] + 1;
      }
    }
    ++i;
  }
  for (i = 0; i < CurrMol.getNOA(); ++i)
    List[i] = NULL;
}


void
initializeRDCmatrix(Molecule &CurrMol,
                    BasicInformation &baseInformation,
                    Flags &flags)
{
  Eigen::MatrixXd rdcMatrix = Eigen::MatrixXd::Zero(
      baseInformation.NumberOfRDCs, baseInformation.NumberOfSets);
  Eigen::MatrixXd rdcMatrixScaled = rdcMatrix;
  Eigen::MatrixXd rdcMatrixNorm   = rdcMatrix;
  RDCset *CurrSet                 = CurrMol.getHeadSet();
  RDCdata *CurrData;
  unsigned int i, j;
  for (i = 0; i < CurrMol.getNORsets() && CurrSet; ++i)
  {
    CurrData = CurrSet->getHeadData();
    for (j = 0; j < CurrMol.getNOR() && CurrData; ++j)
    {
      rdcMatrix(j, i) = CurrData->getD();
      CurrData        = CurrData->getNext();
    }
    CurrSet = CurrSet->getNext();
  }

  /* Fill the rdcMatrix with the respective elements found in the input. */

  CurrMol.setRDCMatrix(rdcMatrix);
  ScaleRDCMatrix(rdcMatrix, rdcMatrixScaled, rdcMatrixNorm, CurrMol,
                 baseInformation, StructureOptions::Initial);

  CurrMol.getHeadStruc()->setRDCmatrix(rdcMatrixOptions::Unscaled);
  CurrMol.getHeadStruc()->updateRDCmatrix(baseInformation, flags,
                                          StructureOptions::Initial);

  CurrSet  = NULL;
  CurrData = NULL;
  delete CurrSet;
  delete CurrData;
}


int
predictRDCs(Molecule &CurrMol, BasicInformation &baseInformation)
{
  if (baseInformation.numberOfFills < 0)
    baseInformation.numberOfFills = 0;
  *std::cin.tie() << "\nInformation:\tStarting predict rdcs...\n\t\t"
                  << baseInformation.numberOfFills
                  << " datasets will be created by random guesses...\n";

  int Aindex, i;

  Atom *CurrAtom, *PartnerAtom;
  std::fstream output;
  Eigen::MatrixXcd F = Eigen::MatrixXcd::Zero(5, 1);

  for (Aindex = 0; Aindex < (baseInformation.numberOfAlignments +
                             baseInformation.numberOfFills);
       ++Aindex)
  {
    // Check if random alignment condition is requested.
    if (Aindex >= baseInformation.numberOfAlignments)
    {
      // Allocate space for information on randomly generated alignment
      // conditions.
      baseInformation.alignments[Aindex] =
          (double *) malloc(5 * sizeof(double));
      // Add additional randomization by running rand() for x times.
      for (i = (rand() % 100); i < (rand() % 1000); ++i)
        ;
      // Generate Azz < 1E-2
      baseInformation.alignments[Aindex][0] = (rand() % 100000) / 1e8;
      // Generate R [0, 2/3].
      baseInformation.alignments[Aindex][1] =
          2.0 / 3.0 - ((rand() % 10000) / (1e4 / (2.0 / 3.0)));
      // Generate alpha [0°, 180°].
      baseInformation.alignments[Aindex][2] = (rand() % 18000) / 100.0;
      // Chose way of calculating beta. (rand[0°, 360°] - 180° or 180° -
      // rand[0°, 360°])
      if (rand() % 2)
      {
        // Generate beta [ -180°, 180° ]
        baseInformation.alignments[Aindex][3] =
            ((rand() % 36000) / 100.0) - 180.0;
      }
      else
      {
        // Generate beta [ -180°, 180° ]
        baseInformation.alignments[Aindex][3] =
            180.0 - (rand() % 36000) / 100.0;
      }
      // Chose way of calculating gamma. (rand[0°, 360°] - 180° or 180° -
      // rand[0°, 360°])
      if (rand() % 2)
      {
        // Generate gamma [ -180°, 180° ]
        baseInformation.alignments[Aindex][4] =
            ((rand() % 36000) / 100.0) - 180.0;
      }
      else
      {
        // Generate gamma [ -180°, 180° ]
        baseInformation.alignments[Aindex][4] =
            180.0 - (rand() % 36000) / 100.0;
      }
      *std::cin.tie() << "\t\tStarting predict set #" << Aindex
                      << " based on random orientation\n";
    } // end of if ( Aindex >= baseInformation.numberOfAlignments )
    else
    {
      *std::cin.tie() << "\t\tStarting predefined set #" << Aindex << std::endl;
    }

    // Build full wigner rotation matrix.
    for (int i = 0; i < 5; ++i)
    {
      F(i, 0) =
          baseInformation.alignments[Aindex][0] * sqrt(4.0 * PI_ / 5.0) *
          (D2Mm(i - 2, 0, -deg2rad(baseInformation.alignments[Aindex][4]),
                -deg2rad(baseInformation.alignments[Aindex][3]),
                -deg2rad(baseInformation.alignments[Aindex][2])) +
           sqrt(3.0 / 8.0) * baseInformation.alignments[Aindex][1] *
               (D2Mm(i - 2, 2, -deg2rad(baseInformation.alignments[Aindex][4]),
                     -deg2rad(baseInformation.alignments[Aindex][3]),
                     -deg2rad(baseInformation.alignments[Aindex][2])) +
                D2Mm(i - 2, -2, -deg2rad(baseInformation.alignments[Aindex][4]),
                     -deg2rad(baseInformation.alignments[Aindex][3]),
                     -deg2rad(baseInformation.alignments[Aindex][2]))));
    }
    // Open output file for current alignment condition.
    output.open(baseInformation.workingDirectory + CurrMol.getLabel() +
                    std::to_string(Aindex) + ".rdc",
                std::ios::binary | std::ios::out | std::ios::trunc);
    CurrAtom = CurrMol.getHeadStruc()->getHeadAtom();
    // Save the uesed alignment conditions.
    output << "Azz = " << baseInformation.alignments[Aindex][0] << ";"
           << " R = " << baseInformation.alignments[Aindex][1] << ";"
           << " alpha = " << baseInformation.alignments[Aindex][2] << ";"
           << " beta = " << baseInformation.alignments[Aindex][3] << ";"
           << " gamma = " << baseInformation.alignments[Aindex][4] << std::endl;

    CurrAtom = CurrMol.getHeadStruc()->getHeadAtom();
    // Loop all atoms...
    while (CurrAtom)
    {
      PartnerAtom = CurrAtom->getNext();
      // ... and generate all combinations...
      while (PartnerAtom)
      {
        // ... calculate and output all combinations with the respective RDCs.
        output << CurrAtom->getIdentifier() << ";"
               << PartnerAtom->getIdentifier() << ";" << std::setprecision(5)
               << Dcalc(CurrAtom, PartnerAtom, F) << ";0.5;1.0\n";
        PartnerAtom = PartnerAtom->getNext();
      }
      CurrAtom = CurrAtom->getNext();
    }
    // Finish the output of current set.
    output.close();
    output.clear();
  } // for ( Aindex = 0; Aindex < ( baseInformation.numberOfAlignments +
    // baseInformation.numberOfFills ) ; ++Aindex )
  *std::cin.tie() << "\nInformation:\tPrediction of RDCs is finished...\n\n";
  return 0;
}

double
Dcalc(Atom *a, Atom *b, const Eigen::MatrixXcd &F)
{
  Coordinates *Ca, *Cb;
  Ca                 = a->getCoordinates(StructureOptions::Initial);
  Cb                 = b->getCoordinates(StructureOptions::Initial);
  double r           = a->getDistance(b, StructureOptions::Initial);
  double theta       = acos((Ca->z - Cb->z) / r);
  double phi         = atan2(((Ca->y - Cb->y) / r), ((Ca->x - Cb->x) / r));
  Eigen::MatrixXcd Y = Eigen::MatrixXcd::Zero(1, 5);
  for (int i = 0; i < 5; ++i)
    Y(0, i) = Ylm(2, i - 2, theta, phi);
  double kappa =
      -(a->getGamma() * b->getGamma() * MU_0_ * H_BAR_) / (8.0 * PI_ * PI_);
  Eigen::MatrixXd tmp = (Y * F).real();
  double D = ((kappa / pow((r * pow(10.0, -10.0)), 3.0)) * tmp(0, 0));
  Ca       = NULL;
  Cb       = NULL;
  delete Ca;
  delete Cb;
  return D;
}

double
Dcalc(const Eigen::MatrixXcd &Y,
      const Eigen::MatrixXcd &F,
      const double kappa,
      const double r)
{
  Eigen::MatrixXd tmp = (Y * F).real();
  double D            = ((kappa / pow(r, 3.0)) * tmp(0, 0));
  return D;
}

Eigen::MatrixXd
Dcalc(Eigen::MatrixXd &b, Eigen::MatrixXd &C, Eigen::MatrixXd &A)
{
  Eigen::MatrixXd dcalc;
  dcalc = b * C * A;
  return dcalc;
}

int
linkCascade(RDCset *S1, RDCdata *D1, RDCset *S2, RDCdata *D2)
{
  int state1, state2;
  if (D1 && D2)
  {
    D1->setNextSetData(D2);
    state1 = state2 = GOOD_STATE;
    if (D1 == S1->getHeadData() && D2 == S2->getHeadData() && S2->getNext())
      state1 = linkCascade(S2, D2, S2->getNext(), NULL);
    if (D1 != S1->getTailData() && D2 != S2->getTailData())
      state2 = linkCascade(S1, D1->getNext(), S2, D2->getNext());
    return (state1 | state2);
  }
  if (S1 && S2)
  {
    return linkCascade(S1, S1->getHeadData(), S2, S2->getHeadData());
  }
  else
  {
    std::cerr << "ERROR:\tMisslinking on RDCs: ";
    if (S1)
    {
      std::cerr << S1->getLabel() << " (";
      if (D1)
      {
        std::cerr << D1->getAtom1()->getIdentifier() << " - "
                  << D1->getAtom2()->getIdentifier() << ") ";
      }
      else
        std::cerr << "unknown RDC pair) ";
    }
    else
      std::cerr << "unkown set 1 ";
    std::cerr << " and ";
    if (S2)
    {
      std::cerr << S2->getLabel() << " (";
      if (D2)
      {
        std::cerr << D2->getAtom1()->getIdentifier() << " - "
                  << D2->getAtom2()->getIdentifier() << ") ";
      }
      else
        std::cerr << "unknown RDC pair) ";
    }
    else
      std::cerr << "unkown set 2";
    std::cerr << "...\n";
    return TITANIA_RDC_LINKING_ERROR;
  }
}
