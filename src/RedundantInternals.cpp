
#include <Atom.hpp>
#include <Bond.hpp>
#include <Declarations.hpp>
#include <Helper.hpp>
#include <LinAlg.hpp>
#include <Output.hpp>
#include <Potential.hpp>
#include <RedundantInternals.hpp>
#include <SphericalHarmonics.hpp>
#include <SphericalHarmonicsMatrix.hpp>
#include <Structure.hpp>
#include <StructureSimulator.hpp>
#include <eigen3/Eigen/Geometry>
#include <iostream>

/****************************************************
 *                RedundantInterals                 *
 *                                                  *
 *     Starting the definition of the redundant     *
 *      internal coordinates. This class will       *
 *     define everything needed to interconvert     *
 *  between cartesian coordinates and (redundant)   *
 *   internal coordinates. Additionally the first   *
 *   derivative of the internal coordinates with    *
 *  respect to the cartesian coordinates (dqi/dxi)  *
 *  will automatically computed. This will be used  *
 * for the interconvertion of the redudant internal *
 *    coordinates and the cartesian coordinates.    *
 *                                                  *
 *                                                  *
 ****************************************************/
/*
 * The implementation is based on the publications:
 * C. Peng, P. Y. Ayala, H. B. Schlegel, M. J. Frisch, Journal of Computational
 * Chemistry 1996, 17, 49–56,
 * DOI 10.1002/(SICI)1096-987X(19960115)17:1<49::AID-JCC5>3.0.CO;2-0. V. Bakken,
 * T. Helgaker, Journal of Chemical Physics 2002, 117, 9160–9174,
 * DOI 10.1063/1.1515483.
 */

RedundantInternals::RedundantInternals(BasicInformation &baseInformation,
                                       Flags &flags)
{
  InternalsList = NULL;
  ListOfAtoms   = NULL;
  CurrStruc     = NULL;
  q0            = Eigen::VectorXd::Zero(1);
  qk            = Eigen::VectorXd::Zero(1);
  Sq            = Eigen::VectorXd::Zero(1);
  sq            = Eigen::VectorXd::Zero(1);
  x0            = Eigen::VectorXd::Zero(1);
  xk            = Eigen::VectorXd::Zero(1);
  x             = Eigen::VectorXd::Zero(1);
  dx_redundants = Eigen::VectorXd::Zero(1);
  dx_distances  = Eigen::VectorXd::Zero(1);
  WilsonB       = Eigen::MatrixXd::Zero(1, 1);

  maxCycles               = baseInformation.limits.max_redundant_cycles;
  rdc_inversion_threshold = STAN_RDC_INVERSION_THRESHOLD;
  convergence_limit       = baseInformation.limits.redundants_convergence;
  restarted               = false;
  floating_rdc_angles     = flags.floating_rdc_angles;
  use_Torsions            = flags.torsions_4_redundants;
  use_distances =
      ((baseInformation.redundants_distance_optimization == FULL_DISTANCES_) ||
       (baseInformation.redundants_distance_optimization ==
        REDUNDANT_DISTANCES_));
  use_long_range_only = flags.long_range_only_4_redundants;

  damping = 1.0;
  static_bond_weighting =
      baseInformation.static_redundants_weighting[BOND_REDUNDANTS_];
  static_angle_weighting =
      baseInformation.static_redundants_weighting[ANGLE_REDUNDANTS_];
  static_torsion_weighting =
      baseInformation.static_redundants_weighting[TORSION_REDUNDANTS_];
  static_rdc_weighting =
      baseInformation.static_redundants_weighting[RDC_REDUNDANTS_];
  static_planar_weighting =
      baseInformation.static_redundants_weighting[PLANAR_REDUNDANTS_];
  static_distance_weighting =
      baseInformation.static_redundants_weighting[DISTANCE_REDUNDANTS_];
  rmsd_q = rmsd_x = .0;
  stop_crit       = .0;

  cycles             = 0;
  NumberOfRedundants = 0;
  NumberOfBonds      = 0;
  NumberOfAngles = cNumberOfAngles = 0;
  NumberOfTorsions = cNumberOfTorsions = 0;
  NumberOfAtoms                        = 0;
  NumberOfCoordinates                  = 0;
  NumberOfRDCAngles = cNumberOfRDCAngles = 0;
  NumberOfChiralVolumes                  = 0;
  NumberOfDistances                      = 0;
}

RedundantInternals::RedundantInternals(StructureSimulator *StrucSim,
                                       BasicInformation &baseInformation,
                                       Flags &flags)
{
  // First do the standard definitions
  maxCycles               = baseInformation.limits.max_redundant_cycles;
  rdc_inversion_threshold = STAN_RDC_INVERSION_THRESHOLD;
  convergence_limit       = baseInformation.limits.redundants_convergence;
  restarted               = false;
  floating_rdc_angles     = flags.floating_rdc_angles;
  use_Torsions            = flags.torsions_4_redundants;
  use_distances =
      ((baseInformation.redundants_distance_optimization == FULL_DISTANCES_) ||
       (baseInformation.redundants_distance_optimization ==
        REDUNDANT_DISTANCES_));
  use_long_range_only = flags.long_range_only_4_redundants;
  cycles              = 0;
  CurrStruc           = StrucSim->getCurrStruc();

  damping = 1.0;
  static_bond_weighting =
      baseInformation.static_redundants_weighting[BOND_REDUNDANTS_];
  static_angle_weighting =
      baseInformation.static_redundants_weighting[ANGLE_REDUNDANTS_];
  static_torsion_weighting =
      baseInformation.static_redundants_weighting[TORSION_REDUNDANTS_];
  static_rdc_weighting =
      baseInformation.static_redundants_weighting[RDC_REDUNDANTS_];
  static_planar_weighting =
      baseInformation.static_redundants_weighting[PLANAR_REDUNDANTS_];
  static_distance_weighting =
      baseInformation.static_redundants_weighting[DISTANCE_REDUNDANTS_];
  rmsd_q = rmsd_x = .0;
  stop_crit       = .0;

  // Now count the number of redundant internal coordinates
  countRedundants(StrucSim);

  InternalsList = (InternalCoordinate *) malloc(NumberOfRedundants *
                                                sizeof(InternalCoordinate));
  for (unsigned int i = 0; i < NumberOfRedundants; ++i)
    InternalsList[i].participants = NULL;
  ListOfAtoms   = NULL;
  q0            = Eigen::VectorXd::Zero(NumberOfRedundants);
  qk            = Eigen::VectorXd::Zero(NumberOfRedundants);
  Sq            = Eigen::VectorXd::Zero(NumberOfRedundants);
  sq            = Eigen::VectorXd::Zero(NumberOfRedundants);
  x0            = Eigen::VectorXd::Zero(NumberOfCoordinates);
  xk            = Eigen::VectorXd::Zero(NumberOfCoordinates);
  x             = Eigen::VectorXd::Zero(NumberOfCoordinates);
  dx_redundants = Eigen::VectorXd::Zero(NumberOfCoordinates);
  dx_distances  = Eigen::VectorXd::Zero(NumberOfCoordinates);
  WilsonB = Eigen::MatrixXd::Zero(NumberOfRedundants, NumberOfCoordinates);

#pragma omp task shared(StrucSim)
  this->setupRedundants(StrucSim);

#pragma omp task shared(CurrStruc)
  this->setupCartesians(CurrStruc);

#pragma omp taskwait
}

RedundantInternals::~RedundantInternals()
{
  if (InternalsList)
    free(InternalsList);
  for (unsigned int i = 0; (ListOfAtoms && i < NumberOfAtoms); ++i)
    ListOfAtoms[i] = NULL;
  free(ListOfAtoms);
  CurrStruc = NULL;
}

void
RedundantInternals::countRedundants(StructureSimulator *StrucSim)
{
  NumberOfAtoms       = StrucSim->getNumberOfAtoms();
  NumberOfDistances   = getSizeOfUpperTriangle(NumberOfAtoms);
  NumberOfCoordinates = (3 * NumberOfAtoms);
  NumberOfBonds       = StrucSim->getNumberOfBonds();
  NumberOfAngles      = StrucSim->getNumberOfAngles();
  if (use_Torsions)
    NumberOfTorsions = StrucSim->getNumberOfTorsions();
  else
    NumberOfTorsions = 0;

  if (use_long_range_only)
    NumberOfRDCAngles = CurrStruc->getNumberOfLongRangeRDCs();
  else
    NumberOfRDCAngles = CurrStruc->getNOR();
  NumberOfChiralVolumes = StrucSim->getNumberOfPlanarCenters();

  NumberOfRedundants = NumberOfBonds + NumberOfAngles + NumberOfTorsions +
                       NumberOfRDCAngles + NumberOfChiralVolumes;
  cNumberOfAngles    = NumberOfAngles + NumberOfBonds;
  cNumberOfTorsions  = cNumberOfAngles + NumberOfTorsions;
  cNumberOfRDCAngles = cNumberOfTorsions + NumberOfRDCAngles;
}

void
RedundantInternals::setupRedundants(StructureSimulator *StrucSim)
{
  unsigned int index = 0;
  this->setupDistances(StrucSim, index);
  this->setupAngles(StrucSim, index);
  this->setupTorsions(StrucSim, index);
  this->setupRDCs(index);
  this->setupChiralVolumes(StrucSim, index);
}

void
RedundantInternals::setupSq(RedundantCalculation rc)
{
  unsigned int index = 0;
  Sq                 = Eigen::VectorXd::Zero(NumberOfRedundants);
  for (; index < NumberOfBonds; ++index)
  {
    Sq(index) = InternalsList[index].potential->get_equilibrium();
    InternalsList[index].index = index;
  }
  for (; index < cNumberOfAngles; ++index)
  {
    Sq(index)                  = this->calculateAngleQ(InternalsList[index],
                                      StructureOptions::Optimized, rc);
    InternalsList[index].index = index;
  }
  for (; index < cNumberOfTorsions; ++index)
  {
    Sq(index)                  = this->calculateTorsionQ(InternalsList[index],
                                        StructureOptions::Optimized, rc);
    InternalsList[index].index = index;
  }
  for (; index < cNumberOfRDCAngles; ++index)
  {
    Sq(index)                  = .0;
    InternalsList[index].index = index;
  }
  for (; index < NumberOfRedundants; ++index)
  {
    Sq(index)                  = .0;
    InternalsList[index].index = index;
  }
}

void
RedundantInternals::setupSq_distances()
{
  sq_distances = Eigen::VectorXd::Zero(NumberOfDistances);
  Eigen::MatrixXd distance_matrix =
      CurrStruc->getDistances(StructureOptions::Optimized);
  unsigned int a, b, d;
  double opt_dist;
  for (a = d = 0; a < NumberOfAtoms && d < NumberOfDistances; ++a)
  {
    for (b = (a + 1); b < NumberOfAtoms && d < NumberOfDistances; ++b, ++d)
    {
      if (ListOfAtoms[a]->getBond(ListOfAtoms[b]))
        sq_distances(d) = .0;
      else
      {
        opt_dist        = get_optimal_distance(ListOfAtoms[a]->getA(),
                                        ListOfAtoms[b]->getA());
        sq_distances(d) = (distance_matrix(a, b) > opt_dist ?
                               .0 :
                               (opt_dist - distance_matrix(a, b)));
        if (ListOfAtoms[a]->getHarmonic(ListOfAtoms[b]))
          sq_distances(d) *= 0.1;
        if (sq_distances(d) != .0)
        {
          ListOfAtoms[a]->increase_distance_violations();
          ListOfAtoms[b]->increase_distance_violations();
        }
        sq_distances(d) *= static_distance_weighting;
      }
    }
  }
}

void
RedundantInternals::setupDistances(StructureSimulator *StrucSim,
                                   unsigned int &index)
{
  Potential *CurrP = StrucSim->getPotential(PotentialType::Stretch);
  unsigned int i, j;
  for (i = 0; i < NumberOfBonds && index < NumberOfRedundants; ++i, ++index)
  {
    InternalsList[index].type = RedundantType::Distance;
    if (InternalsList[index].participants == NULL)
      InternalsList[index].participants = (Atom **) malloc(2 * sizeof(Atom *));
    for (j = 0; j < 2; ++j)
      InternalsList[index].participants[j] = CurrP->getAtom(j);
    InternalsList[index].value = q0(index) = this->calculateDistanceQ(
        InternalsList[index], StructureOptions::Optimized);
    InternalsList[index].potential = CurrP;
    CurrP                          = CurrP->getNext(PotentialType::Stretch);
  }
  CurrP = NULL;
}

void
RedundantInternals::setupAngles(StructureSimulator *StrucSim,
                                unsigned int &index)
{
  Potential *CurrP = StrucSim->getPotential(PotentialType::Bend);
  unsigned i, j;
  for (i = 0; i < NumberOfAngles && index < NumberOfRedundants; ++i, ++index)
  {
    InternalsList[index].type = RedundantType::Angle;
    if (InternalsList[index].participants == NULL)
      InternalsList[index].participants = (Atom **) malloc(3 * sizeof(Atom *));
    for (j = 0; j < 3; ++j)
      InternalsList[index].participants[j] = CurrP->getAtom(j);
    InternalsList[index].value = q0(index) = this->calculateAngleQ(
        InternalsList[index], StructureOptions::Optimized);
    InternalsList[index].potential = CurrP;

    CurrP = CurrP->getNext(PotentialType::Bend);
  }
  CurrP = NULL;
}

void
RedundantInternals::setupRDCs(unsigned int &index)
{
  SphericalHarmonics *CurrY = CurrStruc->getHeadYmatrix()->getHeadHarmonic();
  for (; CurrY && index < NumberOfRedundants; CurrY = CurrY->getNext())
  {
    if (use_long_range_only && CurrY->getRange() < 2)
      continue;
    InternalsList[index].type     = RedundantType::RDC;
    InternalsList[index].harmonic = CurrY;
    if (InternalsList[index].participants == NULL)
      InternalsList[index].participants = (Atom **) malloc(2 * sizeof(Atom *));
    InternalsList[index].participants[0] = CurrY->getAtom1();
    InternalsList[index].participants[1] = CurrY->getAtom2();

    InternalsList[index].value = q0(index) = this->calculateRDCangleQ(
        InternalsList[index], StructureOptions::Optimized);
    ++index;
  }
}

void
RedundantInternals::setupChiralVolumes(StructureSimulator *StrucSim,
                                       unsigned int &index)
{
  Potential *CurrP = StrucSim->getPotential(PotentialType::oop);
  for (; CurrP && index < NumberOfRedundants;
       CurrP = CurrP->getNext(PotentialType::oop))
  {
    InternalsList[index].type = RedundantType::oop;
    if (InternalsList[index].participants == NULL)
      InternalsList[index].participants = (Atom **) malloc(4 * sizeof(Atom *));

    for (unsigned int i = 0; i < 4; ++i)
      InternalsList[index].participants[i] = CurrP->getAtom(i);

    InternalsList[index].value = q0(index) = this->calculateChiralVQ(
        InternalsList[index], StructureOptions::Optimized);
    ++index;
  }
}

void
RedundantInternals::setupTorsions(StructureSimulator *StrucSim,
                                  unsigned int &index)
{
  Potential *CurrP = StrucSim->getPotential(PotentialType::Dihedral);
  unsigned int i, j;
  for (i = 0; i < NumberOfTorsions && index < NumberOfRedundants; ++i, ++index)
  {
    InternalsList[index].type = RedundantType::Torsion;
    if (InternalsList[index].participants == NULL)
      InternalsList[index].participants = (Atom **) malloc(4 * sizeof(Atom *));
    for (j = 0; j < 4; ++j)
      InternalsList[index].participants[j] = CurrP->getAtom(j);
    InternalsList[index].value = q0(index) = this->calculateTorsionQ(
        InternalsList[index], StructureOptions::Optimized);
    InternalsList[index].potential = CurrP;
    CurrP                          = CurrP->getNext();
  }
  CurrP = NULL;
}

void
RedundantInternals::setupCartesians(Structure *Struc)
{
  Atom *A = Struc->getHeadAtom();

  unsigned int i = 0, j = 0, k;
  if (ListOfAtoms == NULL)
    ListOfAtoms = (Atom **) malloc(NumberOfAtoms * sizeof(Atom *));
  Eigen::Vector3d coord;
  while (A)
  {
    ListOfAtoms[j++] = A;
    coord            = A->Coordinates2Eigen(StructureOptions::Optimized);
    for (k = 0; k < NUMBER_OF_AXIS_; ++k, ++i)
      x0(i) = coord(k);

    A = A->getNext();
    if (i == NumberOfCoordinates || j == NumberOfAtoms)
      break;
  }
  A = NULL;
}

void
RedundantInternals::rebaseCartesians(Structure *Struc)
{
  Atom *A        = Struc->getHeadAtom();
  unsigned int i = 0;
  Eigen::Vector3d coord;
  while (A)
  {
    coord(STANDARD_AXIS_X_) = x(i++);
    coord(STANDARD_AXIS_Y_) = x(i++);
    coord(STANDARD_AXIS_Z_) = x(i++);
#pragma omp task firstprivate(coord, A)
    A->setCoordinates(coord, StructureOptions::Optimized);

    A = A->getNext();
    if (i == NumberOfCoordinates)
      break;
  }
  A = NULL;
#pragma omp taskwait
}

void
RedundantInternals::setupWilsonB(StructureOptions opt, RedundantCalculation rc)
{
  unsigned int i, b;
  i = b = 0;
  //   #pragma omp parallel for
  //   shared(NumberOfCoordinates,NumberOfRedundants,WilsonB)
  for (i = 0; i < NumberOfCoordinates; ++i)
  {
    //      #pragma omp parallel for
    //      shared(NumberOfCoordinates,NumberOfRedundants,WilsonB,i)
    for (b = 0; b < NumberOfRedundants; ++b)
    {
      WilsonB(b, i) = this->getDerivative(i, b, opt, rc);
    }
  }
}

void
RedundantInternals::setupWilsonB_distances(StructureOptions opt,
                                           RedundantCalculation rc)
{
  unsigned int a, b, d, d_atom, d_axis, distance_index;
  distance_index = 0;
  for (a = 0; a < NumberOfAtoms; ++a)
  {
    for (b = (a + 1); b < NumberOfAtoms; ++b, ++distance_index)
    {
      for (d = 0; d < NumberOfCoordinates; ++d)
      {
        d_atom = d / NUMBER_OF_AXIS_;
        d_axis = d % NUMBER_OF_AXIS_;
        WilsonB_distances(distance_index, d) =
            drdx(a, b, d_atom, d_axis, opt, rc);
      }
    }
  }
}

double
RedundantInternals::getDerivative(const unsigned int i,
                                  const unsigned int b,
                                  StructureOptions opt,
                                  RedundantCalculation rc)
{
  switch (InternalsList[b].type)
  {
    case RedundantType::Distance:
      return this->drdx(i, b, opt, rc);
    case RedundantType::Angle:
      return this->dadx(i, b, opt, rc);
    case RedundantType::Torsion:
      return this->dtdx(i, b, opt, rc);
    case RedundantType::RDC:
      return this->drdcdx(i, b, opt, rc);
    case RedundantType::oop:
      return this->dVcdx(i, b, opt);
    default:
      return .0;
  }
}

double
RedundantInternals::drdx(const unsigned int i,
                         const unsigned int b,
                         StructureOptions opt,
                         RedundantCalculation rc)
{
  unsigned int index = i / 3;
  unsigned int axis  = i % 3;
  unsigned int A     = (InternalsList[b].participants[0]->getIndex() - 1);
  unsigned int B     = (InternalsList[b].participants[1]->getIndex() - 1);
  return drdx(A, B, index, axis, opt, rc);
}

double
RedundantInternals::drdx(const unsigned int A,
                         const unsigned int B,
                         const unsigned int D,
                         const unsigned int axis,
                         StructureOptions opt,
                         RedundantCalculation rc)
{
  Eigen::VectorXd u = getVec(ListOfAtoms[B], ListOfAtoms[A], opt, rc);
  u                 = u / u.norm();
  double der =
      (zetaFunction(ListOfAtoms[D], ListOfAtoms[B], ListOfAtoms[A]) * u(axis));
  if (der != der)
    return .0;
  else
    return der;
}

double
RedundantInternals::dadx(const unsigned int i,
                         const unsigned int b,
                         StructureOptions opt,
                         RedundantCalculation rc)
{
  unsigned int index = i / 3;
  unsigned int j     = i % 3;

  double zamo =
      zetaFunction(ListOfAtoms[(index)], InternalsList[b].participants[2],
                   InternalsList[b].participants[1]);
  double zano =
      zetaFunction(ListOfAtoms[(index)], InternalsList[b].participants[0],
                   InternalsList[b].participants[1]);
  if (!(zamo + zano))
    return .0;
  Eigen::Vector3d ub =
      RedundantInternals::getVec(InternalsList[b].participants[2],
                                 InternalsList[b].participants[1], opt, rc);
  Eigen::Vector3d vb =
      RedundantInternals::getVec(InternalsList[b].participants[0],
                                 InternalsList[b].participants[1], opt, rc);
  Eigen::Vector3d wb =
      RedundantInternals::setupWvector(InternalsList[b], opt, rc);

  double lu = ub.norm();
  double lv = vb.norm();
  double lw = wb.norm();

  Eigen::Vector3d u = ub / lu;
  Eigen::Vector3d v = vb / lv;
  Eigen::Vector3d w = wb / lw;

  Eigen::Vector3d uxw = u.cross(w);
  Eigen::Vector3d wxv = w.cross(v);
  double der          = zamo * (uxw(j) / lu) + zano * (wxv(j) / lv);
  if (der != der)
    return .0;
  else
    return der;
}

double
RedundantInternals::drdcdx(const unsigned int i,
                           const unsigned int b,
                           StructureOptions opt,
                           RedundantCalculation rc)
{
  unsigned int index = i / 3;
  unsigned int j     = i % 3;

  double zamo =
      zetaFunction(ListOfAtoms[index], InternalsList[b].participants[0],
                   InternalsList[b].participants[1]);
  //   double zano = 0;

  Eigen::Vector3d ub =
      RedundantInternals::getVec(InternalsList[b].participants[0],
                                 InternalsList[b].participants[1], opt, rc);
  Eigen::Vector3d vb = Polar2Eigen(InternalsList[b].harmonic->getTheta(opt),
                                   InternalsList[b].harmonic->getPhi(opt));

  double lu = ub.norm();
  double lv = vb.norm();

  Eigen::Vector3d u  = ub / lu;
  Eigen::Vector3d v  = vb / lv;
  Eigen::Vector3d wb = u.cross(v);
  double lw          = wb.norm();
  Eigen::Vector3d w  = wb / lw;

  Eigen::Vector3d uxw = u.cross(w);
  double der          = zamo * (uxw(j) / lu);
  if (der != der || lw < 1e-10)
    return .0;
  else
    return der;
}

double
RedundantInternals::dVcdx(const unsigned int i,
                          const unsigned int b,
                          StructureOptions opt) //, RedundantCalculation rc )
{
  unsigned int index                = i / 3;
  unsigned int axis_to_derive       = i % 3;
  unsigned int coordinate_to_derive = 4;

  for (unsigned int k = 0; k < 4; ++k)
  {
    if (ListOfAtoms[index] == InternalsList[b].participants[k])
    {
      coordinate_to_derive = k;
      break;
    }
  }
  if (coordinate_to_derive == 4)
    return .0;
  //   else std::cout << "Derivative: R(" << axis_to_derive << "," <<
  //   coordinate_to_derive << ") of " << ListOfAtoms[index]->getIdentifier()
  //                  << " at " <<
  //                  InternalsList[b].participants[3]->getIdentifier() << "\n";

  Eigen::Vector3d ub = RedundantInternals::getVec(
      InternalsList[b].participants[0], InternalsList[b].participants[3], opt,
      RedundantCalculation::QInitial); // Ligand, Central, allways Q
  Eigen::Vector3d vb = RedundantInternals::getVec(
      InternalsList[b].participants[1], InternalsList[b].participants[3], opt,
      RedundantCalculation::QInitial); // Ligand, Central, allways Q
  Eigen::Vector3d wb = RedundantInternals::getVec(
      InternalsList[b].participants[2], InternalsList[b].participants[3], opt,
      RedundantCalculation::QInitial); // Ligand, Central, allways Q
  Eigen::Matrix3d Det;

  Det << ub, vb, wb;


  double der = get_derived_Lebinitz_determinante(Det, axis_to_derive,
                                                 coordinate_to_derive);
  if (der != der)
    return .0;
  else
    return der;
}

double
RedundantInternals::dtdx(const unsigned int i,
                         const unsigned int b,
                         StructureOptions opt,
                         RedundantCalculation rc)
{
  unsigned int index = i / 3;
  unsigned int j     = i % 3;

  Eigen::Vector3d ub =
      RedundantInternals::getVec(InternalsList[b].participants[3],
                                 InternalsList[b].participants[2], opt, rc);
  Eigen::Vector3d wb =
      RedundantInternals::getVec(InternalsList[b].participants[1],
                                 InternalsList[b].participants[2], opt, rc);
  Eigen::Vector3d vb =
      RedundantInternals::getVec(InternalsList[b].participants[0],
                                 InternalsList[b].participants[1], opt, rc);

  double lu = ub.norm();
  double lw = wb.norm();
  double lv = vb.norm();

  Eigen::Vector3d u = ub / lu;
  Eigen::Vector3d w = wb / lw;
  Eigen::Vector3d v = vb / lv;

  double zamo =
      zetaFunction(ListOfAtoms[index], InternalsList[b].participants[3],
                   InternalsList[b].participants[2]);
  double zapn =
      zetaFunction(ListOfAtoms[index], InternalsList[b].participants[1],
                   InternalsList[b].participants[0]);
  double zaop =
      zetaFunction(ListOfAtoms[index], InternalsList[b].participants[2],
                   InternalsList[b].participants[1]);

  Eigen::Vector3d uxw = u.cross(w);
  Eigen::Vector3d vxw = v.cross(w);

  double cospu = u.transpose() * w;
  double cospv = -1.0 * v.transpose() * w;

  double der = zamo * (uxw(j) / (lu * pow(sin(acos(cospu)), 2.0))) +
               zapn * (vxw(j) / (lv * pow(sin(acos(cospv)), 2.0))) +
               zaop * ((uxw(j) * cospu / (lw * pow(sin(acos(cospu)), 2))) +
                       (vxw(j) * cospv / (lw * pow(sin(acos(cospv)), 2))));
  if (der != der)
    return .0;
  else
    return der;
}

int
RedundantInternals::newStructure(Structure *newS,
                                 StructureSimulator *StrucSim,
                                 StructureOptions options)
{
  // First do the standard definitions
  restarted = false;
  cycles    = 0;
  CurrStruc = newS;

  // Clear vectors / matrices
  q0            = Eigen::VectorXd::Zero(NumberOfRedundants);
  qk            = Eigen::VectorXd::Zero(NumberOfRedundants);
  Sq            = Eigen::VectorXd::Zero(NumberOfRedundants);
  sq            = Eigen::VectorXd::Zero(NumberOfRedundants);
  x0            = Eigen::VectorXd::Zero(NumberOfCoordinates);
  xk            = Eigen::VectorXd::Zero(NumberOfCoordinates);
  x             = Eigen::VectorXd::Zero(NumberOfCoordinates);
  dx_redundants = Eigen::VectorXd::Zero(NumberOfCoordinates);
  dx_distances  = Eigen::VectorXd::Zero(NumberOfCoordinates);
  WilsonB = Eigen::MatrixXd::Zero(NumberOfRedundants, NumberOfCoordinates);

  //   #pragma omp task shared(StrucSim)
  //      this->setupRedundants ( StrucSim );

#pragma omp task shared(StrucSim)
  this->setupRedundants(StrucSim);

#pragma omp task shared(CurrStruc)
  this->setupCartesians(CurrStruc);

#pragma omp taskwait
  return 0;
}

int
RedundantInternals::S2x(BasicInformation &bi, Flags &flags)
{
  x = x0;
  Eigen::MatrixXd MPI_WilsonB, MPI_WilsonB_distances;
  if (flags.redundants_damping)
    damping = calculate_damping(bi);

#pragma omp task
  setupSq_distances();

  if (use_distances)
  {
#pragma omp task shared(MPI_WilsonB_distances)
    {
      WilsonB_distances =
          Eigen::MatrixXd::Zero(NumberOfDistances, NumberOfCoordinates);
      setupWilsonB_distances(StructureOptions::Optimized);
      MPI_WilsonB_distances =
          MoorePenroseInverse<double>(WilsonB_distances, 0.5, flags.use_gpu);
    }
  }
  else
  {
    dx_distances = Eigen::VectorXd::Zero(NumberOfCoordinates);
  }

#pragma omp task shared(MPI_WilsonB)
  {
    this->setupWilsonB(StructureOptions::Optimized);
    MPI_WilsonB = MoorePenroseInverse<double>(WilsonB, 0.5, flags.use_gpu);
  }

#pragma omp task shared(sq, Sq, q0, x, x0)
  {
    this->setupSq();
    sq = Sq - q0;
  }
#pragma omp taskwait

  if (use_distances)
  {
#pragma omp task shared(sq, WilsonB, MPI_WilsonB)
    sq = WilsonB * MPI_WilsonB * sq;
#pragma omp task shared(sq_distances, WilsonB_distances, MPI_WilsonB_distances)
    sq_distances = WilsonB_distances * MPI_WilsonB_distances * sq_distances;
#pragma omp taskwait
  }
  else
    sq = WilsonB * MPI_WilsonB * sq;

  Eigen::VectorXd SQ = sq;

  double pre    = .0;
  double sqNorm = sq.norm();

  for (cycles = 0; cycles < maxCycles; ++cycles)
  {
#pragma omp task shared(dx_redundants, MPI_WilsonB, sq)
    dx_redundants = MPI_WilsonB * sq;

    if (use_distances)
    {
#pragma omp task shared(dx_distances, MPI_WilsonB_distances, sq_distances)
      dx_distances = MPI_WilsonB_distances * sq_distances;
    }

#pragma omp taskwait
    x = x + dx_redundants;
    if (use_distances)
      x = x + dx_distances;

    this->rebaseCartesians(CurrStruc);

    if (flags.print_redundants)
    {
      std::string cycle = std::to_string(cycles);
      outputXYZ(CurrStruc->getHeadAtom(), bi, StructureOptions::Optimized,
                ("Redundant_" + cycle));
    }

    if (use_distances)
    {
#pragma omp task shared(MPI_WilsonB_distances)
      {
        setupWilsonB_distances(StructureOptions::Optimized);
        MPI_WilsonB_distances =
            MoorePenroseInverse<double>(WilsonB_distances, 0.5, flags.use_gpu);
      }
#pragma omp task
      setupSq_distances();
    }

#pragma omp task shared(WilsonB, MPI_WilsonB)
    {
      this->setupWilsonB(StructureOptions::Optimized);
      MPI_WilsonB = MoorePenroseInverse<double>(WilsonB, 0.5, flags.use_gpu);
    }

// Recalculate q
#pragma omp task shared(sq, SQ, qk, q0, damping)
    {
      this->recalcQk(StructureOptions::Optimized);
      sq = this->calculate_sq(SQ, damping);
    }

#pragma omp taskwait

    if (use_distances)
    {
#pragma omp task shared(sq, WilsonB, MPI_WilsonB)
      sq = WilsonB * MPI_WilsonB * sq;
#pragma omp task shared(sq_distances, WilsonB_distances, MPI_WilsonB_distances)
      sq_distances = WilsonB_distances * MPI_WilsonB_distances * sq_distances;
#pragma omp taskwait
    }
    else
      sq = WilsonB * MPI_WilsonB * sq;

    sqNorm    = sq.norm();
    stop_crit = fabs(pre - sqNorm);
    if (stop_crit < convergence_limit)
      break;
    else
      pre = sqNorm;
  }

#pragma omp task firstprivate(CurrStruc) shared(bi)
  {
    double rmsd = sq.segment(cNumberOfTorsions, NumberOfRDCAngles).norm();
    CurrStruc->setRDC_rmsd(rmsd);
  }

  int state = check_validity(bi, flags);
  if (state)
    return state;
#pragma omp task shared(dx_redundants, dx_distances, sq, rmsd_x, rmsd_q, \
                        NumberOfCoordinates, NumberOfRedundants)
  {
    rmsd_x = (dx_distances + dx_redundants).norm() / sqrt(NumberOfCoordinates);
    rmsd_q = sq.norm() / sqrt(NumberOfRedundants);
  }
  this->rebaseCartesians(CurrStruc);

  if (!flags.print_redundants && flags.skipEckart)
    outputXYZ(CurrStruc->getHeadAtom(), bi, StructureOptions::Optimized,
              "Redundant-opt");
  CurrStruc->setPerformedRedundantSteps(cycles);
  return GOOD_STATE;
}

double
RedundantInternals::calculateDistanceQ(const InternalCoordinate &q,
                                       StructureOptions opt,
                                       RedundantCalculation rc)
{
  Eigen::Vector3d u =
      this->getVec(q.participants[1], q.participants[0], opt, rc);
  double d = u.norm();
  return d;
}

double
RedundantInternals::calculateAngleQ(const InternalCoordinate &q,
                                    StructureOptions opt,
                                    RedundantCalculation rc)
{
  if (rc == RedundantCalculation::SqRDCs)
  {
    //      if ( !( q.participants[0]->getBond(q.participants[1])->getHarmonic()
    //      && q.participants[2]->getBond(q.participants[1])->getHarmonic() ) )
    return deg2rad(q.potential->get_equilibrium());
  }
  else if (rc == RedundantCalculation::Standard)
    return deg2rad(q.potential->get_equilibrium());
  double a = .0;
  Eigen::Vector3d ub =
      RedundantInternals::getVec(q.participants[2], q.participants[1], opt, rc);
  Eigen::Vector3d vb =
      RedundantInternals::getVec(q.participants[0], q.participants[1], opt, rc);
  Eigen::Vector3d u = ub / ub.norm();
  Eigen::Vector3d v = vb / vb.norm();
  a                 = acos(u.transpose() * v);

  return a;
}

double
RedundantInternals::calculateRDCangleQ(const InternalCoordinate &q,
                                       StructureOptions opt,
                                       RedundantCalculation rc)
{
  double a = .0, t, p;
  Eigen::Vector3d ub =
      RedundantInternals::getVec(q.participants[0], q.participants[1], opt, rc);
  Eigen::Vector3d vb;
  t                 = q.harmonic->getTheta(opt);
  p                 = q.harmonic->getPhi(opt);
  vb(0)             = sin(t) * cos(p);
  vb(1)             = sin(t) * sin(p);
  vb(2)             = cos(t);
  Eigen::Vector3d u = ub / ub.norm();

  a = acos(u.transpose() * vb);
  if (a != a)
    return .0;
  while (a >= PI_)
    a -= PI_;

  if (q.harmonic->getRange() == 1 && ub.norm() < rdc_inversion_threshold)
  {
    return -(a + PI_);
  }
  return (a < 1e-5 ? .0 : a);
}

double
RedundantInternals::calculateChiralVQ(
    const InternalCoordinate &q,
    StructureOptions opt //,
    //                                               RedundantCalculation rc
)
{
  double V;
  Eigen::Vector3d ub = RedundantInternals::getVec(
      q.participants[0], q.participants[3], opt,
      RedundantCalculation::QInitial); // Ligand, Central, allways Q
  Eigen::Vector3d vb = RedundantInternals::getVec(
      q.participants[1], q.participants[3], opt,
      RedundantCalculation::QInitial); // Ligand, Central, allways Q
  Eigen::Vector3d wb = RedundantInternals::getVec(
      q.participants[2], q.participants[3], opt,
      RedundantCalculation::QInitial); // Ligand, Central, allways Q
  Eigen::Matrix3d Det;

  Det << ub, vb, wb;
  V = Det.determinant();
  return V;
}

Eigen::Vector3d
RedundantInternals::setupWvector(const InternalCoordinate &q,
                                 StructureOptions opt,
                                 RedundantCalculation rc)
{
  Eigen::Vector3d w = Eigen::Vector3d::Zero();
  double a          = .0;
  Eigen::Vector3d ub =
      RedundantInternals::getVec(q.participants[2], q.participants[1], opt, rc);
  Eigen::Vector3d vb =
      RedundantInternals::getVec(q.participants[0], q.participants[1], opt, rc);
  Eigen::Vector3d u = ub / ub.norm();
  Eigen::Vector3d v = vb / vb.norm();
  a                 = acos(u.transpose() * v);

  if (a == .0 || a == PI_)
  {
    vb(0) = 1.0;
    vb(1) = -1.0;
    vb(2) = 1.0;
    v     = vb / vb.norm();
    a     = acos(u.transpose() * v);
    if (a == .0 || a == PI_)
    {
      vb(0) = -1.0, vb(1) = 1.0;
      v = vb / vb.norm();
    }
  }
  w = u.cross(v);
  return w;
}

Eigen::Vector3d
RedundantInternals::getVec(Atom *n,
                           Atom *m,
                           StructureOptions opt,
                           RedundantCalculation rc)
{
  if (rc == RedundantCalculation::QInitial ||
      n->getBond(m)->getHarmonic() == NULL)
    return (n->Coordinates2Eigen(opt) - m->Coordinates2Eigen(opt));

  Eigen::Vector3d V = Polar2Eigen(n->getBond(m)->getHarmonic()->getTheta(opt),
                                  n->getBond(m)->getHarmonic()->getPhi(opt));

  Eigen::Vector3d helpb = n->Coordinates2Eigen(opt) - m->Coordinates2Eigen(opt);
  Eigen::Vector3d help  = helpb / helpb.norm();
  if ((V.transpose() * help) < .0)
    V = (-1.0) * V;
  V = V * helpb.norm();
  return V;
}

Eigen::VectorXd
RedundantInternals::calculate_sq(Eigen::VectorXd &SQ, double &d)
{
  sq = SQ - (qk - q0);
  sq = apply_damping(d);
  return sq;
}

Eigen::VectorXd
RedundantInternals::apply_damping(double d)
{
  unsigned int i = 0;

  sq.head(NumberOfBonds) *= (d * static_bond_weighting);
  for (i = NumberOfBonds; i < cNumberOfAngles; ++i)
  {
    if (InternalsList[i]
                .participants[0]
                ->getBond(InternalsList[i].participants[1])
                ->getHarmonic() == NULL &&
        InternalsList[i]
                .participants[2]
                ->getBond(InternalsList[i].participants[1])
                ->getHarmonic() == NULL)
      sq(i) *= (d * static_angle_weighting);
    else if (floating_rdc_angles &&
             (InternalsList[i]
                      .participants[0]
                      ->getBond(InternalsList[i].participants[1])
                      ->getHarmonic() == NULL ||
              InternalsList[i]
                      .participants[2]
                      ->getBond(InternalsList[i].participants[1])
                      ->getHarmonic() == NULL))
      sq(i) = .0;
    else if (InternalsList[i].participants[0]->getDistance(
                 InternalsList[i].participants[2],
                 StructureOptions::Optimized) < 0.5)
      sq(i) = .0;
  }

  if (use_Torsions)
    sq.segment(cNumberOfAngles, NumberOfTorsions) *=
        (d * static_torsion_weighting);

  sq.segment(cNumberOfTorsions, NumberOfRDCAngles) *=
      (d * static_rdc_weighting);
  sq.tail(NumberOfChiralVolumes) *= (d * static_planar_weighting);

  return sq;
}

void
RedundantInternals::recalcQk(StructureOptions opt)
{
  qk                 = q0;
  unsigned int index = 0;
  unsigned int i;

  for (i = 0; i < NumberOfBonds && index < NumberOfRedundants; ++i, ++index)
  {
    qk(index) = this->calculateDistanceQ(InternalsList[index], opt);
  }

  for (i = 0; i < NumberOfAngles && index < NumberOfRedundants; ++i, ++index)
  {
    qk(index) = this->calculateAngleQ(InternalsList[index], opt);
  }

  for (i = 0; i < NumberOfTorsions && index < NumberOfRedundants; ++i, ++index)
  {
    qk(index) = this->calculateTorsionQ(InternalsList[index], opt);
  }

  for (i = 0; i < NumberOfRDCAngles && index < NumberOfRedundants; ++i, ++index)
  {
    qk(index) = this->calculateRDCangleQ(InternalsList[index], opt);
  }

  for (i = 0; i < NumberOfChiralVolumes && index < NumberOfRedundants;
       ++i, ++index)
  {
    qk(index) = this->calculateChiralVQ(InternalsList[index], opt);
  }
}

double
RedundantInternals::calculate_damping(BasicInformation &baseInformation)
{
  double step     = (double) baseInformation.numOfOptSteps;
  double maxSteps = (double) baseInformation.limits.max_titania_iterations;
  double fac, damping, exp_1, damper, x;
  x       = step / maxSteps;
  damper  = 3.5;
  fac     = baseInformation.redundants_damping;
  exp_1   = exp(damper);
  damping = exp(damper * x) * exp(damper * x) /
            (fac + exp(damper * x) * exp(damper * x));

  const double norm = 1.0 / (exp_1 * exp_1 / (fac + exp_1 * exp_1));
  return (damping * norm);
}

int
RedundantInternals::check_validity(BasicInformation &baseInformation,
                                   Flags &flags)
{
  double rg_ratio, rg_i, rg_o;
  rg_i     = CurrStruc->get_radius_of_Gyration(StructureOptions::Initial);
  rg_o     = CurrStruc->get_radius_of_Gyration(StructureOptions::Optimized);
  rg_ratio = rg_o / rg_i;

  if (rg_ratio > baseInformation.limits.redundants_validity)
  {
    CurrStruc->retainCoordinates();

    if (restarted)
    {
      std::string runtimeChange = "TITANIA is going to turn off redundant "
                                  "internal coordinates for this step...\n";
      if (flags.printWarnings)
        *std::cin.tie() << runtimeChange;
      baseInformation.runtime_changes += runtimeChange;
      return ERROR_INVALID_REDUNDANTS_STEP;
    }


    maxCycles         = 5;
    convergence_limit = 0.1;
    restarted         = true;

    if (baseInformation.runtime_changes.find("extraordinary") !=
        std::string::npos)
    {
      std::string runtimeChange = "Repeated extraordinary increase (";
      runtimeChange += std::to_string(rg_ratio);
      runtimeChange += "-fold) in iteration ";
      runtimeChange += std::to_string(baseInformation.numOfOptSteps);
      runtimeChange += "...\n";
      if (flags.printWarnings)
        *std::cin.tie() << runtimeChange;
      baseInformation.runtime_changes += runtimeChange;
    }
    else
    {
      baseInformation.runtime_changes +=
          "TITANIA detected an extraordinary increase (";
      baseInformation.runtime_changes += std::to_string(rg_ratio);
      baseInformation.runtime_changes +=
          "-fold) in the radius of gyration in iteration ";
      baseInformation.runtime_changes +=
          std::to_string(baseInformation.numOfOptSteps);
      baseInformation.runtime_changes += "...\n";
      baseInformation.runtime_changes +=
          "This normally is due to too strict optimization values for "
          "redundant internal coordinates...\n";
      baseInformation.runtime_changes +=
          "The new values are:\n   * convergence limit: ";
      baseInformation.runtime_changes += std::to_string(convergence_limit);
      baseInformation.runtime_changes += "\n   * max interations: ";
      baseInformation.runtime_changes += std::to_string(maxCycles);
      baseInformation.runtime_changes += "\n";
      if (flags.printWarnings)
        *std::cin.tie() << baseInformation.runtime_changes;
    }
    recalcQk(StructureOptions::Optimized);
    this->setupCartesians(CurrStruc);
    return S2x(baseInformation, flags);
  }
  return GOOD_STATE;
}

double
RedundantInternals::calculateTorsionQ(const InternalCoordinate &q,
                                      StructureOptions opt,
                                      RedundantCalculation rc)
{
  double t = .0;

  Eigen::Vector3d ub =
      RedundantInternals::getVec(q.participants[3], q.participants[2], opt, rc);
  Eigen::Vector3d wb =
      RedundantInternals::getVec(q.participants[1], q.participants[2], opt, rc);
  Eigen::Vector3d vb =
      RedundantInternals::getVec(q.participants[0], q.participants[1], opt, rc);

  Eigen::Vector3d u = ub / ub.norm();
  Eigen::Vector3d w = wb / wb.norm();
  Eigen::Vector3d v = vb / vb.norm();

  t = ((u.cross(w)).transpose()) * v.cross(w);
  t /= (sqrt(1.0 - pow(u.transpose() * w, 2)) *
        sqrt(1.0 - pow(v.transpose() * w, 2)));

  return acos(t);
}

void
RedundantInternals::refactorInternals(Eigen::VectorXd &q)
{
  unsigned int i;
  for (i = NumberOfBonds; i < cNumberOfAngles; ++i)
  {
    if (q(i) >= (2.0 * PI_))
    {
      std::cout << i << ": " << q(i) << std::endl;
      q(i) = q(i) - 2.0 * PI_;
      --i;
    }
  }
}

double
get_optimal_distance(int A, int B)
{
  int a, b;
  if (A > B)
  {
    a = A;
    b = B;
  }
  else
  {
    a = B;
    b = A;
  }

  if (a == 1 && b == 1)
    return 1.3; // 5;
  else if (a == 6 && b == 1)
    return 1.6; // 8;
  else if (a == 7 && b == 1)
    return 1.4; // 8;
  else if (a == 8 && b == 1)
    return 1.4; // 8;
  else if (a == 6 && b == 6)
    return 1.8; // 2.2;
  else if (a == 7 && b == 6)
    return 1.8; // 2.2;
  else if (a == 8 && b == 6)
    return 1.9; // 2.3;
  else
    return 1.9; // 2.3;
}
