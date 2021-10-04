
#ifndef REDUNDANTINTERNALS_HPP_
#define REDUNDANTINTERNALS_HPP_

#include <eigen3/Eigen/Core>

#ifndef STAN_MAX_REDUNDANT_CYCLES
#define STAN_MAX_REDUNDANT_CYCLES 5
#endif

#ifndef STAN_RDC_INVERSION_THRESHOLD
#define STAN_RDC_INVERSION_THRESHOLD 0.8
#endif

#ifndef ATOM_HPP_
class Atom;
#endif

#ifndef STRUCTURE_HPP_
class Structure;
#endif

#ifndef STRUCTURE_SIMULATOR_HPP
class StructureSimulator;
#endif

#ifndef REDUNANDTTYPE_TITANIA_
enum struct RedundantType : unsigned int;
#endif

#ifndef INTERNALCOORDINATE_TITANIA_
struct InternalCoordinate;
#endif

#ifndef REDUNDANT_CALCULATION_TITANIA_
#include <RedundantCalculation.hpp>
#endif

#ifndef STRUCTUREOPTIONS_TITANIA_
#include <StructureOptions.hpp>
#endif

#ifndef BASIC_INFORMATION_TITANIA_
struct BasicInformation;
#endif

#ifndef FLAGS_TITANIA_
struct Flags;
#endif

class RedundantInternals {
 private:
  // General Variables

  Structure *CurrStruc;

  // Definition of internal

  InternalCoordinate *InternalsList; // Types of the redundant coordinates
  Eigen::VectorXd q0;                // Initial redundant internals
  Eigen::VectorXd qk;                // Redundant internals of step k
  unsigned int NumberOfRedundants;
  unsigned int NumberOfBonds;
  unsigned int NumberOfAngles;
  unsigned int cNumberOfAngles;
  unsigned int NumberOfTorsions;
  unsigned int cNumberOfTorsions;
  unsigned int NumberOfRDCAngles;
  unsigned int cNumberOfRDCAngles;
  unsigned int NumberOfChiralVolumes;
  unsigned int NumberOfDistances;

  // Definition of cartesian coordinates

  Atom **ListOfAtoms;
  Eigen::VectorXd x0; // Initial cartesian coordinates
  Eigen::VectorXd xk; // Cartesian coordinates of step k
  Eigen::VectorXd x;  // Current cartesian coordinates
  Eigen::VectorXd dx_redundants;
  Eigen::VectorXd dx_distances;
  unsigned int NumberOfAtoms;
  unsigned int NumberOfCoordinates; // 3 * NumberOfAtoms

  // Definition of interconversion variables

  Eigen::MatrixXd WilsonB; // Transformation matrix according
  Eigen::MatrixXd WilsonB_distances;
  Eigen::VectorXd Sq; // Optimized redundant internals
  Eigen::VectorXd
      sq; // Difference between current and optimized redundant internals
  Eigen::VectorXd sq_distances;

  // Monitoring of the interconversion

  unsigned int
      maxCycles; // Defines the number of maximum cycles for the interconversion
  unsigned int
      cycles; // Number of cycles to obtain optimized cartesian coordinates

  double rdc_inversion_threshold;
  double convergence_limit;
  bool use_Torsions;
  bool use_distances;
  bool use_long_range_only;
  bool restarted;
  bool floating_rdc_angles;

  double damping;
  double static_bond_weighting;
  double static_angle_weighting;
  double static_torsion_weighting;
  double static_rdc_weighting;
  double static_planar_weighting;
  double static_distance_weighting;

  double rmsd_x;
  double rmsd_q;

  double stop_crit;

  // Internal functions
  void countRedundants(StructureSimulator *);
  void setupRedundants(StructureSimulator *);
  void setupCartesians(Structure *);
  void setupSq(RedundantCalculation rc = RedundantCalculation::SqRDCs);
  void setupSq_distances(); // RedundantCalculation rc =
                            // RedundantCalculation::SqRDCs );
  // void setupAtoms (Structure* ); // Automatically run while setup Cartesians

  void setupDistances(StructureSimulator *, unsigned int &);
  void setupAngles(StructureSimulator *, unsigned int &);
  void setupTorsions(StructureSimulator *, unsigned int &);
  void setupRDCs(unsigned int &);
  void setupChiralVolumes(StructureSimulator *, unsigned int &);

  double
  calculateDistanceQ(const InternalCoordinate &,
                     StructureOptions,
                     RedundantCalculation rc = RedundantCalculation::QInitial);
  double
  calculateAngleQ(const InternalCoordinate &,
                  StructureOptions,
                  RedundantCalculation rc = RedundantCalculation::QInitial);
  double
  calculateRDCangleQ(const InternalCoordinate &,
                     StructureOptions,
                     RedundantCalculation rc = RedundantCalculation::QInitial);
  double
  calculateTorsionQ(const InternalCoordinate &,
                    StructureOptions,
                    RedundantCalculation rc = RedundantCalculation::QInitial);
  double
  calculateChiralVQ(const InternalCoordinate &,
                    StructureOptions); //, RedundantCalculation rc =
                                       // RedundantCalculation::QInitial );
  void refactorInternals(Eigen::VectorXd &);

  static Eigen::Vector3d
  setupWvector(const InternalCoordinate &,
               StructureOptions,
               RedundantCalculation rc = RedundantCalculation::QInitial);
  static Eigen::Vector3d
  getVec(Atom *,
         Atom *,
         StructureOptions,
         RedundantCalculation rc = RedundantCalculation::QInitial);
  Eigen::VectorXd calculate_sq(Eigen::VectorXd &, double &);
  Eigen::VectorXd apply_damping(double);
  void recalcQk(StructureOptions);
  double calculate_damping(BasicInformation &);

  int check_validity(BasicInformation &, Flags &);

 public:
  RedundantInternals(BasicInformation &, Flags &);
  RedundantInternals(StructureSimulator *, BasicInformation &, Flags &);
  ~RedundantInternals();

  unsigned int getNumberOfAtoms()
  {
    return NumberOfAtoms;
  }
  unsigned int getNumberOfCoordinates()
  {
    return NumberOfCoordinates;
  }

  unsigned int getNumberOfBonds()
  {
    return NumberOfBonds;
  }
  unsigned int getNumberOfAngles()
  {
    return NumberOfAngles;
  }
  unsigned int getNumberOfTorsions()
  {
    return NumberOfTorsions;
  }
  unsigned int getNumberOfChiralVolumes()
  {
    return NumberOfChiralVolumes;
  }
  unsigned int getNumberOfRedundants()
  {
    return NumberOfRedundants;
  }
  unsigned int getNumberOfDistance()
  {
    return NumberOfDistances;
  }

  unsigned int getMaxCycles() const
  {
    return maxCycles;
  }
  unsigned int getCycles() const
  {
    return cycles;
  }

  double getDamping() const
  {
    return damping;
  }

  double getRMSDx() const
  {
    return rmsd_x;
  }
  double getRMSDq() const
  {
    return rmsd_q;
  }

  double getStopCrit() const
  {
    return stop_crit;
  }

  void setupWilsonB(StructureOptions,
                    RedundantCalculation rc = RedundantCalculation::QInitial);
  void setupWilsonB_distances(
      StructureOptions opt,
      RedundantCalculation rc = RedundantCalculation::QInitial);
  double
  getDerivative(unsigned int,
                unsigned int,
                StructureOptions,
                RedundantCalculation rc = RedundantCalculation::QInitial);
  double drdx(const unsigned int,
              const unsigned int,
              StructureOptions,
              RedundantCalculation rc = RedundantCalculation::QInitial);
  double drdx(const unsigned int,
              const unsigned int,
              const unsigned int,
              const unsigned int,
              StructureOptions,
              RedundantCalculation rc = RedundantCalculation::QInitial);
  double dadx(const unsigned int,
              const unsigned int,
              StructureOptions,
              RedundantCalculation rc = RedundantCalculation::QInitial);
  double dtdx(const unsigned int,
              const unsigned int,
              StructureOptions,
              RedundantCalculation rc = RedundantCalculation::QInitial);
  double drdcdx(const unsigned int,
                const unsigned int,
                StructureOptions,
                RedundantCalculation rc = RedundantCalculation::QInitial);
  double dVcdx(const unsigned int,
               const unsigned int,
               StructureOptions); //, RedundantCalculation rc =
                                  // RedundantCalculation::QInitial );

  Eigen::MatrixXd getWilsonB()
  {
    return WilsonB;
  }

  int newStructure(Structure *, StructureSimulator *, StructureOptions);
  int S2x(BasicInformation &, Flags &);

  InternalCoordinate *getInternal(unsigned int, RedundantType);
  void rebaseCartesians(Structure *Struc);
};

double get_optimal_distance(int, int);
#endif
