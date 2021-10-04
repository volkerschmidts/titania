
#include <Atom.hpp>
#include <Declarations.hpp>
#include <Helper.hpp>
#include <Molecule.hpp>
#include <Output.hpp>
#include <Parser/Parser.hpp>
#include <Properties.hpp>
#include <RDCdata.hpp>
#include <RDCset.hpp>
#include <Structure.hpp>
#include <iostream>
#include <omp.h>

Molecule::Molecule(std::string newLabel)
{
  label                     = newLabel;
  NumberOfAtoms             = 0;
  NumberOfRDCatoms          = 0;
  NumberOfRDCs              = 0;
  NumberOfRDCsets           = 0;
  NumberOfStructs           = 0;
  NumberOfDoubleBonds       = 0;
  NumberOfLongRangeRDCs     = 0;
  MolecularMass             = .0;
  DBE                       = 0;
  inputB                    = Eigen::MatrixXd::Zero(1, 1);
  kappa                     = Eigen::MatrixXd::Zero(1, 1);
  HeadSet                   = NULL;
  TailSet                   = NULL;
  HeadStruc                 = NULL;
  TailStruc                 = NULL;
  weights                   = Eigen::MatrixXd::Zero(1, 1);
  input_Transformation      = Eigen::MatrixXd::Identity(4, 4);
  input_InertiaTensor       = Eigen::Matrix3d::Zero();
  input_InertiaTensorEigVal = Eigen::Vector3d::Zero();
  input_InertiaTensorEigVec = Eigen::Matrix3d::Zero();
}

Molecule::~Molecule()
{
  *std::cin.tie() << "\nInformation:\tStarting to clean up " << label
                  << "...\n";

  while (HeadSet && HeadSet->getNext())
  {
    HeadSet = HeadSet->getNext();
    delete HeadSet->getPrev();
  }
  delete HeadSet;
  TailSet = NULL;

  *std::cin.tie() << "\t\tRDC data were deleted...\n";

  while (HeadStruc->getNext())
  {
    HeadStruc = HeadStruc->getNext();
    delete HeadStruc->getPrev();
  }
  TailStruc = NULL;
  delete HeadStruc;

  *std::cin.tie() << "\t\tStructures (containing atoms & spherical harmonics) "
                     "were deleted...\n";
}

int
Molecule::loadStructure(InputFile &inf,
                        StructureInput *si,
                        BasicInformation &bi)
{
  unsigned int i, j, NOB;
  double *coordinates;

  std::string atomlabel, bondlabel;
  std::vector<std::string> Keys = si->getKeys();

  Atom *CurrAtom, *BondAtom;

  label         = si->getMolecule();
  NumberOfAtoms = si->getNOA();

  /******************/
  /* ERROR handling */
  /******************/

  if (!(Keys.size()))
  {
    std::cerr << "\nError:\tNo structure was defined in the input file.\n"
              << "\tTITANIA will terminate. Check your input file and the "
                 "manual for more information.\n";
    si          = NULL;
    bi.errorKey = "No structure was properly defined in the input. Check the "
                  "manual for more information.\n";
    bi.state = ERROR_CORRUPTED_STRUCTURE_INPUT;
    return ERROR_CORRUPTED_STRUCTURE_INPUT;
  }

  /******************************/
  /* Handle the structure input */
  /******************************/

  HeadStruc = addHeadStruc();
  HeadStruc->setLabel(si->getLabel());
  CurrAtom = HeadStruc->setHeadAtom();


  /*************************/
  /* Parse the coordinates */
  /*************************/

  for (i = 0; i < Keys.size(); ++i)
  {
    if (i)
      CurrAtom = CurrAtom->setNext();
    if (!CurrAtom)
    {
      std::cerr << "\nERROR:\tTITANIA could not defined atom " << Keys.at(i)
                << " properly.\n"
                << "\tTerminating process.\n";
      bi.state    = ERROR_CORRUPTED_STRUCTURE_INPUT;
      bi.errorKey = "Internal error on setting next atom. Check function "
                    "\"Molecule::loadInput\".\n";
      CurrAtom = NULL;
      si       = NULL;
      return ERROR_CORRUPTED_STRUCTURE_INPUT;
    }
    CurrAtom->setIndex(i + 1);
    atomlabel = Keys.at(i);

    switch (atomLabelParser(atomlabel, CurrAtom))
    {
      case -1: {
        unsigned int line = 0;
        inf.getKeyLine(atomlabel, bi.errorKey, line);
        CurrAtom = NULL;
        si       = NULL;
        return ERROR_CORRUPTED_STRUCTURE_INPUT;
        break;
      }
      case ERROR_IDENTIFIER_MISSING: {
        bi.errorKey = "ERROR:\tThe structure was not defined "
                      "properly...\n\tThe atom with input index ";
        bi.errorKey += std::to_string(i + 1);
        bi.errorKey +=
            " was just defined via element number...\n\tTITANIA needs the "
            "element symbol in combination with an identifier...\n";
        bi.state = ERROR_IDENTIFIER_MISSING;
        CurrAtom = NULL;
        si       = NULL;
        return ERROR_IDENTIFIER_MISSING;
        break;
      }
      default:
        break;
    }
    coordinates = NULL;
    NOB         = 0;

    if (!si->setAtomInformation(atomlabel, &coordinates, NOB))
    {
      unsigned int line = 0;
      std::string tmp   = "";
      si->checkDiff(atomlabel, tmp);
      bi.errorKey = ("ERROR:\tAtom " + atomlabel +
                     " was not defined properly. It does not match the "
                     "following lines:\n\n");
      inf.getKeyLine(tmp, bi.errorKey, line);
      bi.errorKey +=
          ("\n\t" + atomlabel + " can be found in the following lines:\n\n");
      line = 0;
      inf.getKeyLine(atomlabel, bi.errorKey, line);
      CurrAtom = NULL;
      si       = NULL;
      return ERROR_CORRUPTED_STRUCTURE_INPUT;
    }
    if (coordinates)
    {
      CurrAtom->setCoordinates(coordinates[0], coordinates[1], coordinates[2],
                               StructureOptions::Initial);
      bi.structureInput = StructureInputType::xyzCoordinates;
    }
    CurrAtom->setNumberOfBonds(NOB);
    CurrAtom->initializeBonds();
    if (!NOB)
    {
      std::cerr << "\nERROR:\tTITANIA the connectivity of the molecule in "
                   "order to optimize the structuer.\n"
                << "\tFor more information check the manual.\n";
      CurrAtom = NULL;
      si       = NULL;
      return ERROR_CORRUPTED_STRUCTURE_INPUT;
    }
  }

  BondAtom = NULL;
  for (i = 0; i < Keys.size(); ++i)
  {
    CurrAtom  = HeadStruc->getHeadAtom();
    atomlabel = Keys.at(i);
    atomByIdentifier(atomlabel, &CurrAtom);
    for (j = 0; j < CurrAtom->getNumberOfBonds(); ++j)
    {
      BondAtom  = HeadStruc->getHeadAtom();
      bondlabel = si->getBondPartner(atomlabel, j);
      if (bondlabel == "")
      {
        std::cerr << "\nERROR:\tThe connectivity of atom " << atomlabel
                  << " was not defined properly.\n"
                  << "\tTITANIA will be terminated. For more information check "
                     "the manual.\n";
        CurrAtom = BondAtom = NULL;
        si                  = NULL;
        return ERROR_CORRUPTED_STRUCTURE_INPUT;
      }
      atomByIdentifier(bondlabel, &BondAtom);
      if (BondAtom)
        CurrAtom->addBondpartner(BondAtom);
      else
      {
        unsigned int line = 0;
        std::string dump  = "";
        std::string con   = "onnectivity";
        inf.getKeyLine(con, dump, line);
        bi.errorKey +=
            ("ERROR:\tConnectivity was not defined properly. Atom \"" +
             bondlabel + "\" was not found.\n");
        inf.getKeyLine(bondlabel, bi.errorKey, line);
        return ERROR_CORRUPTED_STRUCTURE_INPUT;
      }
    }
  }

  coordinates = NULL;
  CurrAtom = BondAtom = NULL;
  return 0;
}

int
Molecule::loadRDCs(InputFile &inf, RDCinput *ri, BasicInformation &bi)
{
  unsigned int i, j;
  double *rdctable;

  std::string a1, a2, id;
  std::vector<std::string> Keys;
  std::vector<std::string> allRDCs;

  RDCset *CurrSet  = HeadSet;
  RDCdata *CurrRDC = NULL;
  Atom *CurrAtom, *BondAtom;
  /* Generate full rdc list */

  while (ri)
  {
    if (!HeadSet)
      addHeadSet();
    else
      TailSet->appendSet();

    TailSet->setLabel(ri->getLabel());
    TailSet->setIndex(raiseNumberOfRDCsets());

    Keys = ri->getKeys();
    i    = 0;
    if (allRDCs.size() == 0)
      allRDCs.push_back(Keys.at(i++));
    for (; i < Keys.size(); ++i)
    {
      for (j = 0; j < allRDCs.size(); ++j)
      {
        if (!allRDCs.at(j).compare(Keys.at(i)))
          break;
      }
      if (j == allRDCs.size())
        allRDCs.push_back(Keys.at(i));
    }
    ri = ri->getNext();
  }

  ri = inf.getRI();


  /* Initialize the HeadSet */
  CurrSet = HeadSet;
  for (i = 0; i < allRDCs.size(); ++i)
  {
    if (i == 0)
      CurrRDC = CurrSet->addHeadData();
    else
      CurrRDC = CurrRDC->appendData();

    a1 = a2 = id = "";
    id           = allRDCs.at(i);

    j = 0;
    while (id.at(j) != ';')
      ++j;

    a1 = id.substr(0, (j));
    a2 = id.substr((j + 1), (id.size() - j));

    BondAtom = CurrAtom = HeadStruc->getHeadAtom();
    atomByIdentifier(a1, &CurrAtom);
    atomByIdentifier(a2, &BondAtom);
    if (!(CurrAtom && BondAtom))
    {
      std::string set   = ri->getLabel();
      unsigned int line = 0;
      bi.errorKey += "ERROR:\tThe RDC of " + a1 + " and " + a2 +
                     " is not defined properly. Check the respective set\n";
      inf.getKeyLine(set, bi.errorKey, line);
      return ERROR_CORRUPTED_RDC_INPUT;
    }

    CurrRDC->setAtom1(CurrAtom);
    BondAtom->increaseNumberOfHarmonics();
    CurrRDC->setAtom2(BondAtom);
    CurrAtom->increaseNumberOfHarmonics();
    raiseNOR();
  }

  while (ri)
  {
    for (i = 0; i < allRDCs.size(); ++i)
    {
      a1 = a2 = id = "";
      id           = allRDCs.at(i);

      j = 0;
      while (id.at(j) != ';')
        ++j;

      a1 = id.substr(0, (j));
      a2 = id.substr((j + 1), (id.size() - j));

      BondAtom = CurrAtom = HeadStruc->getHeadAtom();
      atomByIdentifier(a1, &CurrAtom);
      atomByIdentifier(a2, &BondAtom);
      if (!(CurrAtom && BondAtom))
      {
        std::string set   = ri->getLabel();
        unsigned int line = 0;
        bi.errorKey += "ERROR:\tThe RDC of " + a1 + " and " + a2 +
                       " is not defined properly. Check the respective set\n";
        inf.getKeyLine(set, bi.errorKey, line);
        return ERROR_CORRUPTED_RDC_INPUT;
      }
      CurrRDC = CurrSet->getHeadData();

      while (CurrRDC->getAtom1() != CurrAtom || CurrRDC->getAtom2() != BondAtom)
        CurrRDC = CurrRDC->getNext();

      rdctable = ri->getValuesD(id);
      if (!rdctable)
      {
        std::string set = ri->getLabel();
        CurrRDC->setUndefined();
        CurrRDC->getAtom1()->getParent()->setUndefined();
        CurrRDC->setWeight(.0);
        *std::cin.tie() << "Information:\tThe RDC of " + a1 + " and " + a2 +
                               " in set "
                        << set << " is not defined properly.\n";
      }
      else
      {
        CurrRDC->setD(rdctable[0]);
        CurrRDC->setDeltaD(rdctable[1]);
        CurrRDC->setWeight(rdctable[2]);
      }
      CurrRDC->setParent(CurrSet);
      CurrRDC = CurrRDC->getNext();
    }
    ri = ri->getNext();
    CurrSet->setNOR(getNOR());
    CurrSet = CurrSet->getNext();
    if (CurrSet)
      CurrSet->copyRDCs(HeadSet);
  }


  CurrSet  = NULL;
  CurrRDC  = NULL;
  CurrAtom = BondAtom = NULL;
  rdctable            = NULL;

  return linkSets();
}

int
Molecule::loadKeywords(KeywordInput *ki, BasicInformation &bi, Flags &flags)
{
  int a               = 0;
  std::string keyline = "";
  unsigned int keyword, i, j;
  keyword = 0;

  std::vector<std::string> Keys = ki->getKeys();

  if (ki->contains("alignCount"))
  {
    a = ki->getValueI("predictrdcs", 1);
    if (a == 0)
      a = ki->getValueI("alignCount", 0);
    bi.numberOfAlignments = a;
    bi.alignments         = (double **) malloc(a * sizeof(double *));
    for (i = 0; i < (unsigned int) a; ++i)
      bi.alignments[i] = NULL;
    a = 0;
  }

  for (i = 0; i < Keys.size(); ++i)
  {
    keyline = Keys.at(i);
    for (j = 0; j < keyline.size(); ++j)
      keyline.at(j) = std::tolower(keyline.at(j));
    keyword = identifyKey(keyline);
    switch (keyword)
    {
      case 0:
        break;
      case 1:
        flags.echo = (ki->getValueI(keyline, 0) > 0 ? true : false);
        break;
      case 8:
        bi.limits.max_lm_iterations = ki->getValueI(Keys.at(i), 0);
        break;
      case 9:
        bi.limits.max_titania_iterations = ki->getValueI(Keys.at(i), 0);
        break;
      case 10:
        break; // just take the standard value here.
      case 11:
        flags.skipSCRM = ki->getValueI(Keys.at(i), 0);
        break;
      case 12: {
        bi.predictRDCs    = ki->getValueI(Keys.at(i), 0);
        int numberOfFills = 0;
        if (bi.predictRDCs == 0)
          break;
        numberOfFills = ki->getValueI(Keys.at(i), 1);
        if (numberOfFills)
          numberOfFills -= (ki->getValueI("alignCount", 0));
        if (bi.numberOfAlignments)
          bi.numberOfAlignments -= numberOfFills;
        bi.numberOfFills = numberOfFills;
        flags.skipSCRM   = true;
        if (bi.alignments == NULL)
        {
          bi.alignments = (double **) malloc(numberOfFills * sizeof(double *));
          for (i = 0; i < (unsigned int) numberOfFills; ++i)
            bi.alignments[i] = NULL;
        }
      }
      break;
      case 13:
        while (bi.alignments[a])
          ++a;
        bi.alignments[a] = ki->getValuesD(Keys.at(i));
        break;
      case 14:
        break; // Deactivated
      case 17:
        break; // silent is handled previously
      case 18:
        bi.limits.Q_factor_convergence = ki->getValueD(keyline, 0);
        break;
      case 19:
        flags.titania2hotfcht = (ki->getValueI(keyline, 0) > 0 ? true : false);
        break;
      case 20:
        flags.outputAli = (ki->getValueI(keyline, 0) > 0 ? true : false);
        break;
      case 21:
        if (bi.numOfThreads)
          break;
        else
        {
          bi.numOfThreads = ki->getValueI(keyline, 0);
          omp_set_num_threads(bi.numOfThreads);
          break;
        }
      case 22:
        flags.monteCarloBootstrapping =
            (ki->getValueI(keyline, 0) > 0 ? true : false);
        break;
      case 23:
        flags.plotKappaQ = (ki->getValueI(keyline, 0) > 0 ? true : false);
        break;
      case 24:
        bi.gnuOutputFormat = ki->getValueS(keyline, 0);
        break;
      case 25:
        flags.numericalGradients =
            (ki->getValueI(keyline, 0) > 0 ? true : false);
        break;
      case 26:
        flags.outputLM = (ki->getValueI(keyline, 0) > 0 ? true : false);
        break;
      case 27:
        bi.limits.redundants_convergence = ki->getValueD(keyline, 0);
        break;
      case 28:
        bi.limits.alignment_mean_convergence = ki->getValueD(keyline, 0);
        break;
      case 29:
        bi.limits.alignment_sigm_convergence = ki->getValueD(keyline, 0);
        break;
      case 30:
        bi.limits.sphericals_mean_convergence = ki->getValueD(keyline, 0);
        break;
      case 31:
        bi.limits.sphericals_sigm_convergence = ki->getValueD(keyline, 0);
        break;
      case 32:
        bi.limits.sphericals_spread_convergence = ki->getValueD(keyline, 0);
        break;
      case 33:
        flags.plotTrajectory = (ki->getValueI(keyline, 0) > 0 ? true : false);
        break;
      case 34:
        flags.recalculateRDCs = (ki->getValueI(keyline, 0) > 0 ? true : false);
        break;
      case 35:
        flags.printWarnings = (ki->getValueI(keyline, 0) > 0 ? true : false);
        break;
      case 36: {
        for (int lmopt = 0; lmopt < 4; ++lmopt)
          bi.limits.phobos_opts_sphericals[lmopt] =
              bi.limits.phobos_opts_eckart[lmopt] =
                  ki->getValueD(keyline, lmopt);
        break;
      }
      case 37: {
        for (int lmopt = 0; lmopt < 4; ++lmopt)
          bi.limits.phobos_opts_eckart[lmopt] = ki->getValueD(keyline, lmopt);
        break;
      }
      case 38: {
        for (int lmopt = 0; lmopt < 4; ++lmopt)
          bi.limits.phobos_opts_sphericals[lmopt] =
              ki->getValueD(keyline, lmopt);
        break;
      }
      case 39:
        flags.skipEckart = (ki->getValueI(keyline, 0) > 0 ? true : false);
        break;
      case 40:
        flags.print_redundants = (ki->getValueI(keyline, 0) > 0 ? true : false);
        break;
      case 41:
        flags.titania2titania = (ki->getValueI(keyline, 0) > 0 ? true : false);
        break;
      case 42:
        flags.normChiralVolume = (ki->getValueI(keyline, 0) > 0 ? true : false);
        break;
      case 43:
        bi.overOptimization = ki->getValueI(keyline, 0);
        break;
      case 44:
        flags.bigData = (ki->getValueI(keyline, 0) > 0 ? true : false);
        break;
      case 45:
        flags.plotRDCdynamics = (ki->getValueI(keyline, 0) > 0 ? true : false);
        break;
      case 46:
        flags.scaleWithSoverall =
            (ki->getValueI(keyline, 0) > 0 ? true : false);
        break;
      case 47:
        flags.monteCarloOutput = (ki->getValueI(keyline, 0) > 0 ? true : false);
        break;
      case 48:
        flags.plotMonteCarlo = (ki->getValueI(keyline, 0) > 0 ? true : false);
        break;
      case 49:
        flags.calculateFullMatrix =
            (ki->getValueI(keyline, 0) > 0 ? true : false);
        break;
      case 50:
        flags.errorWeightInSVD = (ki->getValueI(keyline, 0) > 0 ? true : false);
        break;
      case 51:
        bi.useRedundantsOnlyAfter = ki->getValueI(keyline, 0);
        break;
      case 52:
        bi.limits.max_redundant_cycles = ki->getValueI(keyline, 0);
        break;
      case 53:
        bi.limits.redundants_validity = ki->getValueD(keyline, 0);
        break;
      case 54:
        bi.redundants_damping    = ki->getValueD(keyline, 0);
        flags.redundants_damping = true;
        break;
      case 55:
        flags.use_gpu = (ki->getValueI(keyline, 0) == 1 ? true : false);
        break;
      case 56:
        bi.redundants_distance_optimization = ki->getValueI(keyline, 0);
        break;
      case 57:
        bi.lowerInversionAfter =
            (ki->getValueI(keyline, 0) <= 0 ? -1 : ki->getValueI(keyline, 0));
        break;
      case 58:
        bi.memoryPurge = ((int) (ki->getValueD(keyline, 0) * 10));
        if (bi.memoryPurge >= 10)
          bi.memoryPurge = 9;
        else if (bi.memoryPurge < 0)
          bi.memoryPurge = 0;
        break;
      case 59:
        bi.static_redundants_weighting[BOND_REDUNDANTS_] =
            fabs(ki->getValueD(keyline, 0));
        break;
      case 60:
        bi.static_redundants_weighting[ANGLE_REDUNDANTS_] =
            fabs(ki->getValueD(keyline, 0));
        break;
      case 61:
        bi.static_redundants_weighting[TORSION_REDUNDANTS_] =
            fabs(ki->getValueD(keyline, 0));
        break;
      case 62:
        bi.static_redundants_weighting[RDC_REDUNDANTS_] =
            fabs(ki->getValueD(keyline, 0));
        break;
      case 63:
        bi.static_redundants_weighting[PLANAR_REDUNDANTS_] =
            fabs(ki->getValueD(keyline, 0));
        break;
      case 64:
        bi.static_redundants_weighting[DISTANCE_REDUNDANTS_] =
            fabs(ki->getValueD(keyline, 0));
        break;
      case 65:
        flags.floating_rdc_angles =
            (ki->getValueI(keyline, 0) == 0 ? false : true);
        break;
      case 66:
        flags.use_initial_holonomics =
            (ki->getValueI(keyline, 0) == 0 ? false : true);
        break;
      case 67:
        flags.SECONDAreducedCovariance =
            (ki->getValueI(keyline, 0) == 0 ? false : true);
        break;
      case 666:
        break;
      default:
        break;
    }
  }

  if (bi.useRedundantsOnlyAfter == (-1))
    bi.useRedundantsOnlyAfter =
        (bi.overOptimization + bi.limits.max_titania_iterations);

  return 0;
}

int
Molecule::loadInput(InputFile &inf, BasicInformation &bi, Flags &flags)
{
  unsigned int i, j, errorKey;

  std::string kl;
  unsigned int kw;
  KeywordInput *ki              = inf.getKI();
  std::vector<std::string> Keys = ki->getKeys();

  // First check for communication
  for (i = 0; i < Keys.size(); ++i)
  {
    kl = Keys.at(i);
    for (j = 0; j < kl.size(); ++j)
      kl.at(j) = std::tolower(kl.at(j));
    kw = identifyKey(kl);
    switch (kw)
    {
      case 17:
        flags.silent   = true;
        bi.comFileName = ki->getValueS(kl, 0);
        if (bi.comFileName == "1")
        {
          bi.comFileName = "/dev/null";
        }
        else if (bi.comFileName == "0")
        {
          bi.comFileName = "";
          flags.silent   = false;
        }
        // Untie the output stream if flag "silent" was used
        if (flags.silent)
          shutup(bi);
        break;
      default:
        break;
    }
  }

  errorKey = loadStructure(inf, inf.getSI(), bi);
  if (errorKey)
    return errorKey;
  errorKey = loadRDCs(inf, inf.getRI(), bi);
  if (errorKey)
    return errorKey;
  errorKey = loadKeywords(ki, bi, flags);
  if (errorKey)
    return errorKey;

  ki = NULL;

  bi.NumberOfSets  = getNORsets();
  bi.NumberOfRDCs  = getNOR();
  bi.NumberOfAtoms = getNOA();

  return 0;
}

Structure *
Molecule::addHeadStruc()
{
  TailStruc = HeadStruc = new Structure();
  HeadStruc->setParent(this);
  raiseNumberOfStrucs();
  return HeadStruc;
}
void
Molecule::removeStructures(int amount)
{
  if (amount > 0)
  {
    Structure *anchor = HeadStruc->getNext();
    /*      Atom *oa1, *oa2, *na1, *na2;
          for ( int i = 0; i < amount; ++i ) newHead = newHead->getNext();
          for ( RDCdata *data = HeadSet->getHeadData(); data; data =
       data->getNext() )
          {
             oa1 = data->getAtom1(); na1 = newHead->getAtomByIndex (
       oa1->getIndex() ); oa2 = data->getAtom2(); na2 = newHead->getAtomByIndex
       ( oa2->getIndex() ); for ( RDCdata *set = data; set; set =
       set->getNextSetData() )
             {
                set->setAtom1(na1);
                set->setAtom2(na2);
             }
          }
          oa1 = oa2 = na1 = na2 = NULL;*/
    for (int i = 0; i < amount; ++i)
    {
      anchor = anchor->getNext();
      delete anchor->getPrev();
    }
    HeadStruc->setNext(anchor);
    anchor->setPrev(HeadStruc);
    anchor = NULL;
  }
}


void
Molecule::checkPlanarity(BasicInformation &baseInformation)
{
  double x, y, z;
  x = y = z = .0;

  Atom *CurrAtom = HeadStruc->getHeadAtom();
  Coordinates *C;

  while (CurrAtom)
  {
    C = CurrAtom->getCoordinates(StructureOptions::Initial);
    x += (pow(C->x, 2.0));
    y += (pow(C->y, 2.0));
    z += (pow(C->z, 2.0));
    CurrAtom = CurrAtom->getNext();
  }

  if (x < 1e-3)
    baseInformation.secondaryInput = StructureInputType::planarX;
  else if (y < 1e-3)
    baseInformation.secondaryInput = StructureInputType::planarY;
  else if (z < 1e-3)
    baseInformation.secondaryInput = StructureInputType::planarZ;
  C = NULL;
}

Eigen::MatrixXd
Molecule::getKappaMatrix()
{
  /*
   * Kappa is independent of r(A-B)! So this
   * matrix should be generated once.
   */

  if (kappa.rows() == NumberOfRDCs)
    return kappa;

  kappa = Eigen::MatrixXd::Zero(NumberOfRDCs, NumberOfRDCs);

  RDCdata *r;
  unsigned int i;

  for (i = 0, r = HeadSet->getHeadData(); (r && i < NumberOfRDCs);
       r = r->getNext(), ++i)
    kappa(i, i) = r->getKappa();

  return kappa;
}

void
Molecule::initializeWmatrix()
{
  RDCset *CurrSet;
  RDCdata *CurrData;
  unsigned int r, s;
  weights = Eigen::MatrixXd::Zero(NumberOfRDCs, NumberOfRDCsets);
  for (CurrSet = HeadSet, s = 0; (CurrSet && (s < NumberOfRDCsets));
       CurrSet = CurrSet->getNext(), ++s)
  {
    for (CurrData = CurrSet->getHeadData(), r = 0;
         (CurrSet && (r < NumberOfRDCs)); CurrData = CurrData->getNext(), ++r)
    {
      if (CurrData->isUndefined())
        weights(r, s) = .0;
      else
        weights(r, s) = CurrData->getEffectiveWeight();
    }
  }
}

Eigen::MatrixXd
Molecule::getWmatrix()
{
  return weights;
}

void
Molecule::setInputInertiaTensor(Eigen::Matrix3d IT)
{
  static bool set = false;
  if (set)
    return;
  for (int i = 0; i < NUMBER_OF_AXIS_; ++i)
  {
    for (int j = 0; j < NUMBER_OF_AXIS_; ++j)
      input_InertiaTensor(i, j) = IT(i, j);
  }
  set = true;
}

void
Molecule::setInputInertiaTensorEigVal(Eigen::Vector3d EigVal)
{
  static bool set = false;
  if (set)
    return;
  for (int i = 0; i < NUMBER_OF_AXIS_; ++i)
  {
    input_InertiaTensorEigVal(i) = EigVal(i);
  }
  set = true;
}

void
Molecule::setInputInertiaTensorEigVec(Eigen::Matrix3d EigVec)
{
  static bool set = false;
  if (set)
    return;
  for (int i = 0; i < NUMBER_OF_AXIS_; ++i)
  {
    for (int j = 0; j < NUMBER_OF_AXIS_; ++j)
      input_InertiaTensorEigVec(i, j) = EigVec(i, j);
  }
  set = true;
}

void
Molecule::setInputTransformation(Eigen::Matrix3d Rotation)
{
  static bool set = false;
  if (set)
    return;
  for (int i = 0; i < NUMBER_OF_AXIS_; ++i)
  {
    for (int j = 0; j < NUMBER_OF_AXIS_; ++j)
      input_Transformation(i, j) = Rotation(i, j);
  }
  set = true;
}

void
Molecule::setInputTransformation(Eigen::Vector3d Translation)
{
  static bool set = false;
  if (set)
    return;
  for (int i = 0; i < NUMBER_OF_AXIS_; ++i)
    input_Transformation(i, NUMBER_OF_AXIS_) = Translation(i);
  set = true;
}

RDCset *
Molecule::addHeadSet()
{
  TailSet = HeadSet = new RDCset();
  HeadSet->setParent(this);
  return HeadSet;
}

int
Molecule::linkSets()
{
  if (HeadSet)
    return linkCascade(HeadSet, NULL, HeadSet->getNext(), NULL);
  else
    return GOOD_STATE;
}

void
Molecule::setInputBMatrix(Structure *CurrStruc)
{
  CurrStruc->determineCosineMatrix(StructureOptions::Initial);
  inputB = CurrStruc->getCosineMatrix(StructureOptions::Initial);
}

double
Molecule::determineMolecularMass()
{
  MolecularMass = .0;
  Atom *CurrAtom;
  for (CurrAtom = HeadStruc->getHeadAtom(); CurrAtom;
       CurrAtom = CurrAtom->getNext())
    MolecularMass += CurrAtom->getMass();

  return MolecularMass;
}

void
Molecule::countLongRangeRDCs()
{
  if (!HeadSet)
    return;
  RDCdata *R;
  for (R = HeadSet->getHeadData(); R; R = R->getNext())
  {
    if (R->getRange() > 1)
      ++NumberOfLongRangeRDCs;
  }
}

/**************************************************************
 *                                                            *
 *                      atomLabelParser                       *
 *                                                            *
 * Parses atom labels containing the definition of an special *
 *   isotope, element symbols, labels or just of the atomic   *
 *     number and stores the information to struct Atom.      *
 *  The return value represents the index of the first C_str  *
 *  element after the element. From this point one can parse  *
 *            for i.e. coordinates, bonded atoms.             *
 *                                                            *
 * The return value of -1 means, that the syntax of the label *
 *                   is not interpretable.                    *
 *   ERROR_IDENTIFIER_MISSING (8) means that just an element  *
 *       number was defined. This is not implemented.         *
 *  0 means, that the label was commented using # and should  *
 *                        be skipped.                         *
 *                                                            *
 **************************************************************/

int
atomLabelParser(std::string argument, Atom *CurrAtom)
{
  std::string isotope;
  bool inIso, inAt;
  unsigned int i, j;
  std::string key    = "";
  std::string tmpKey = "";
  std::string atom   = "";
  int mass           = 0;

  inIso = true;
  inAt  = false;

  j = 0;

  for (i = 0; i < argument.size();
       ++i) /* Check for isotope definition on argument */
  {
    // Assume, that there is an isotope definition.
    if (inIso)
    {
      // If identifier starts with a number, there actually is an definition.
      if (isdigit(argument.at(i)))
        continue;
      // Else there isnt (Cpt. Obvious)
      else
        inIso = false;
    }
    // If we did not find any additional number prior to characters we are in
    // the actual identifier.
    else if (!inAt)
    {
      if (i == j) // This should actually never happen since we assume first
                  // entry to be a number.
      {
        break;
      }
      else
      {
        // Save the substring and convert it to an int.
        isotope = argument.substr(j, i);
        // If substring is no int -> mass = 0
        mass = atoi(isotope.c_str());
        break;
      }
    }
  }

  // We now can be sure to be inside the identifier
  inAt = true;

  // But since we are !inside! the identifier we have to start at i-1
  for (j = --i; i < argument.size(); ++i) /* Parse the element symbol */
  {
    if (i == j)
    {
      continue; /* The AtomIdentifier has to be at least one char long! */
    }
    if (inAt) // This allways should be true
    {
      if (isdigit(argument.at(
              i))) // First digit means we have definetly reached the label
      {
        atom = argument.substr(j, (i - j));
        break;
      }
    }
  }
  if (!atom.size())
  {
    return ERROR_IDENTIFIER_MISSING;
  }
  while (!isElement(atom)) // As long as the element symbol is not truncated
  {
    --i;
    atom = argument.substr(j, (i - j)); // shorten the label.
    if ((i - j) == 0)
      break;
  }


  if (mass)
    CurrAtom->setA(mass);
  if (atom.size())
    CurrAtom->setElement(atom);
  if (findNucleusByName(*(CurrAtom)))
    return -1;

  // Everything after the element symbol is the label.
  key = argument.substr(i, (argument.size() - i));

  if (key.size())
    CurrAtom->setLabel(key);

  return 0;
}

int
initializeChiralVolumes(Molecule &CurrMol, BasicInformation &baseInformation)
{
  for (Atom *CurrAtom = CurrMol.getHeadStruc()->getHeadAtom(); CurrAtom;
       CurrAtom       = CurrAtom->getNext())
  {
    if (CurrAtom->hasChiralVolume())
      ++baseInformation.NumberOfChiralCenters;
  }
  baseInformation.chiral_volumes =
      (double *) malloc(baseInformation.NumberOfChiralCenters *
                        (baseInformation.limits.max_titania_iterations +
                         baseInformation.overOptimization + 1) *
                        sizeof(double));
  for (unsigned int i = 0;
       i < baseInformation.NumberOfChiralCenters *
               (baseInformation.limits.max_titania_iterations +
                baseInformation.overOptimization);
       ++i)
    baseInformation.chiral_volumes[i] = .0;

  return 0;
}
