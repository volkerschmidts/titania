
#ifndef DECLARATIONS_HPP_
#define DECLARATIONS_HPP_

//#define WEIGHTED_COVARIANCE_MATRIX      1

#include <fstream>

// Error keys
#define GOOD_STATE                      0
#define ERROR_CORRUPTED_STRUCTURE_INPUT 1
#define ERROR_UNKNOWN_KEYWORDS          2
#define OUTPUT_HELP                     3
#define OUTPUT_VERSION                  4
#define ERROR_ECKART_TRANSFORMATION     5
#define ERROR_CORRUPTED_RDC_INPUT       6
#define OUTPUT_BASE_SETTINGS            7
#define ERROR_IDENTIFIER_MISSING        8
#define ERROR_FILE_DOES_NOT_EXIST       9
#define ERROR_STRUCTURE_EXPLOSION       10
#define TITANIA_MULTI_JOB               11
#define TITANIA_RDC_LINKING_ERROR       12
#define ERROR_INPUT_EQ_OUTPUT           13
#define ERROR_NO_INPUT_FILE             14
#define ERROR_INVALID_REDUNDANTS_STEP   15
#define ERROR_OUTPUT_OPEN               16
#define ERROR_TRJ_OPEN                  17

// LM information for F/Y optimization
#define NUMBER_OF_SPHERICAL_COORDINATES 2
#define NUMBER_OF_SPHERICAL_FUNCTIONS   2
#define PHOBOS_SPHERICAL_COVAR_SIZE     4
#define PHOBOS_THETA                    1
#define PHOBOS_PHI                      0

// LM information for Eckart transformation
#define ECKART_THETA 0
#define ECKART_PHI   1
#define ECKART_CHI   2

#define NUMBER_OF_ECKART_ANGLES   3
#define NUMBER_OF_ECKART_FUNTIONS 3
#define ECKART_COVAR_SIZE         9

// Redundant internal coordinates
#define BOND_REDUNDANTS_           0
#define ANGLE_REDUNDANTS_          1
#define TORSION_REDUNDANTS_        2
#define RDC_REDUNDANTS_            3
#define PLANAR_REDUNDANTS_         4
#define DISTANCE_REDUNDANTS_       5
#define NUMBER_OF_REDUNDANT_TYPES_ 6

#define SKIP_DISTANCES_      0
#define FULL_DISTANCES_      1
#define REDUNDANT_DISTANCES_ 2
#define INVERSION_DISTANCES_ 3

// Some standard io information
#define STANDARD_DIGITS      9
#define STANDARD_BUFFER_SIZE 256
#define CLEAN_LINE           50

// Axis definitions
#ifndef AXESDEFINITIONS_TITANIA_
#include <AxesDefinitions.hpp>
#endif

// Number of X elements
#define SAUPE_ELEMENTS_             5
#define SAUPE_TENSOR_DIM_           3
#define SAUPE_ZZ_                   0
#define SAUPE_XX_YY_                1
#define SAUPE_XY_                   2
#define SAUPE_XZ_                   3
#define SAUPE_YZ_                   4
#define FULL_SAUPE_ELEMENTS_        6
#define WIGNER_ELEMENTS_            5
#define COSINE_ELEMENTS_            5
#define HARMONIC_ELEMENTS_          5
#define NUM_EULER_ANGLES_           3
#define EULER_ALPHA_                0
#define EULER_BETA_                 1
#define EULER_GAMMA_                2
#define NUM_EULER_PERMUTATIONS_     4
#define NUM_EIGENVALUES_            3
#define NUM_EIGENVECTORS_           9
#define EIGENVECTOR_ELEMENTS_       3
#define MIN_SECONDA_SETS_           6
#define NUMBER_OF_TRIG_COORDINATES_ 4
#define LARGE_MATRIX_LIMIT_         2000

// Leibnitz definition (Leibnitz form of determinants)
#define NUMBER_OF_S3_PERMUTATIONS_ 6
#define LIGAND_1_                  0
#define LIGAND_2_                  1
#define LIGAND_3_                  2
#define CENTRAL_ATOM_              3

// Define string lengths for output
#define SHORT_STRING_       4
#define MEDIUM_STRING_      6
#define LONG_STRING_        12
#define EXTRA_STRING_       16
#define EXTRA_EXTRA_STRING_ 28

// Define number formats for output
#define SHORT_NUMBER_  3
#define MEDIUM_NUMBER_ 8
#define LONG_NUMBER_   12
#define EXTRA_NUMBER_  16

#define LOW_PRECISION_    1
#define MEDIUM_PRECISION_ 3
#define HIGH_PRECISION_   5

#ifndef TITANIA_OUTPUT_FORMAT_
#define TITANIA_OUTPUT_FORMAT_
struct Output_Format {
  int string_s;
  int string_m;
  int string_l;
  int string_xl;
  int string_xxl;

  int number_s;
  int number_m;
  int number_l;
  int number_xl;

  int prec_s;
  int prec_m;
  int prec_l;

  Output_Format()
  {
    string_s   = SHORT_STRING_;
    string_m   = MEDIUM_STRING_;
    string_l   = LONG_STRING_;
    string_xl  = EXTRA_STRING_;
    string_xxl = EXTRA_EXTRA_STRING_;

    number_s  = SHORT_NUMBER_;
    number_m  = MEDIUM_NUMBER_;
    number_l  = LONG_NUMBER_;
    number_xl = EXTRA_NUMBER_;

    prec_s = LOW_PRECISION_;
    prec_m = MEDIUM_PRECISION_;
    prec_l = HIGH_PRECISION_;
  }
};
#endif

#define DO_PRAGMA(x) _Pragma(#x)
#define TODO(x)      DO_PRAGMA(message("TODO - " #x))
#define getName(var) #var
#define caller       std::cout << __PRETTY_FUNCTION__ << std::endl;

#include <eigen3/Eigen/Core>

#ifndef PI_
#define PI_ 3.14159265358979323846264338327950288
constexpr double PI_HALF_      = (PI_ / 2.0);
constexpr double PI_TIMES_TWO_ = (PI_ * 2.0);
#endif

#ifndef MU_N_
#define MU_N_ 5.05078369821105E-27
#endif

#ifndef H_
#define H_ 6.626070040E-34
#endif


#ifndef H_BAR_
constexpr double H_BAR_ = (H_ / PI_TIMES_TWO_);
#endif

#ifndef MU_0_
constexpr double MU_0_ = (4.0 * PI_ * 1E-7);
#endif

#ifndef ZERO_CUTOFF_
#define ZERO_CUTOFF_ 1e-9
#endif

#ifndef RAD_2_DEG_TITANIA_
#define RAD_2_DEG_TITANIA_
inline double
rad2deg(double a)
{
  return (a * 180.0 / PI_);
}
inline double
deg2rad(double a)
{
  return (a * PI_ / 180.0);
}
inline Eigen::MatrixXd
rad2deg(Eigen::MatrixXd a)
{
  return (a * 180.0 / PI_);
}
inline Eigen::MatrixXd
deg2rad(Eigen::MatrixXd a)
{
  return (a * PI_ / 180.0);
}
inline Eigen::VectorXd
rad2deg(Eigen::VectorXd a)
{
  return (a * 180.0 / PI_);
}
inline Eigen::VectorXd
deg2rad(Eigen::VectorXd a)
{
  return (a * PI_ / 180.0);
}
inline Eigen::Vector3d
rad2deg(Eigen::Vector3d a)
{
  return (a * 180.0 / PI_);
}
inline Eigen::Vector3d
deg2rad(Eigen::Vector3d a)
{
  return (a * PI_ / 180.0);
}
#endif


#ifndef ATOM_HPP_
class Atom;
#endif

#ifndef BOND_HPP_
class Bond;
#endif

#ifndef POTENTIAL_HPP_
class Potential;
#endif

#ifndef SPHERICALHARMONICS_HPP_
class SphericalHarmonics;
#endif

#ifndef BASIC_INFORMATION_TITANIA_
struct BasicInformation;
#endif

#ifndef AXIS_TITANIA_
#define AXIS_TITANIA_
enum struct Axis : unsigned int
{
  undefined = 0,
  x,
  y,
  z,
  Xi,
  Eta,
  Zeta,
};
#endif

#ifndef ECKART_OPT_PARAMETERS_TITANIA_
#define ECKART_OPT_PARAMETERS_TITANIA_
struct Eckart_Opt_Parameters {
  Eigen::Matrix3d AlphaTau;
  BasicInformation *baseInformation;
};
#endif

#ifndef OUTPUT_OPTIONS_TITANIA_
#include <OutputOptions.hpp>
#endif

#ifndef ALIGNMENT_OUTPUT_TITANIA_
#define ALIGNMENT_OUTPUT_TITANIA_
struct AlignmentOutput {
  Eigen::MatrixXd rdc_calc;
  //   Eigen::MatrixXd rdc_calc_uw;
  Eigen::MatrixXd SaupeTensor;
  //   Eigen::MatrixXd SaupeTensor_uw;
  Eigen::MatrixXd SaupeEigenValues;
  //   Eigen::MatrixXd SaupeEigenValues_uw;
  Eigen::MatrixXd SaupeEigenVectors;
  //   Eigen::MatrixXd SaupeEigenVectors_uw;
  Eigen::MatrixXd EulerAngles;
  //   Eigen::MatrixXd EulerAngles_uw;
  Eigen::MatrixXd w;
  //   Eigen::MatrixXd uw;
  Eigen::VectorXd Q_factor;
  //   Eigen::VectorXd Q_factor_uw;

  AlignmentOutput(){};
  AlignmentOutput(unsigned int NOR, unsigned int NOS)
  {
    /*      rdc_calc_uw = rdc_calc_w = Eigen::MatrixXd::Zero(NOR,NOS);
          SaupeTensor_uw = SaupeTensor_w = Eigen::MatrixXd::Zero(9,NOS);
          SaupeEigenValues_uw = SaupeEigenValues_w =
       Eigen::MatrixXd::Zero(3,NOS); SaupeEigenVectors_uw = SaupeEigenVectors_w
       = Eigen::MatrixXd::Zero(9,NOS); EulerAngles_uw = EulerAngles_w =
       Eigen::MatrixXd::Zero(3,NOS); Q_factor_uw = Q_factor_w =
       Eigen::VectorXd::Zero(NOS); uw = Eigen::MatrixXd::Identity(NOR, NOR);*/
    rdc_calc          = Eigen::MatrixXd::Zero(NOR, NOS);
    SaupeTensor       = Eigen::MatrixXd::Zero(9, NOS);
    SaupeEigenValues  = Eigen::MatrixXd::Zero(3, NOS);
    SaupeEigenVectors = Eigen::MatrixXd::Zero(9, NOS);
    EulerAngles       = Eigen::MatrixXd::Zero(3, NOS);
    Q_factor          = Eigen::VectorXd::Zero(NOS);
    w                 = Eigen::MatrixXd::Identity(NOR, NOR);
  }
};
#endif

#ifndef HYBRIDISATION_TITANIA_
#define HYBRIDISATION_TITANIA_
enum struct Hybridisation : unsigned int
{
  undefined = 0,
  s         = 1,
  sp3       = 1 << 1,
  sp2       = 1 << 2,
  sp1       = 1 << 3,
  coo       = 1 << 4,
};

inline std::ostream &
operator<<(std::ostream &os, Hybridisation &inp)
{
  switch (((unsigned int) inp))
  {
    case (0):
      os << "undefined";
      break;
    case (1):
      os << "s";
      break;
    case (1 << 1):
      os << "sp3";
      break;
    case (1 << 2):
      os << "sp2";
      break;
    case (1 << 3):
      os << "sp";
      break;
    case (1 << 4):
      os << "coordinated";
      break;
    default:
      os << "undefined hybridisation";
  }
  return os;
}
#endif

/* Structure to define the requested rdc matrix. Human reading */
#ifndef RDCMATRIXOPTIONS_TITANIA_
#define RDCMATRIXOPTIONS_TITANIA_
enum struct rdcMatrixOptions : unsigned int
{
  Unscaled = 1,
  Scaled   = 1 << 1,
  Norm     = 1 << 2,
};
#endif

#ifndef STRUCTUREOPTIONS_TITANIA_
#include <StructureOptions.hpp>
#endif /* STRUCTUREOPTIONS_TITANIA_ */

#ifndef FRAGMENT_FLAG_TITANIA_
#include <FragmentFlag.hpp>
#endif

#ifndef REDUNDANTTYPE_TITANIA_
#define REDUNANDTTYPE_TITANIA_
enum struct RedundantType : unsigned int
{
  Undefined = 0,
  Distance  = 1,
  Angle     = 1 << 1,
  Torsion   = 1 << 2,
  RDC       = 1 << 3,
  oop       = 1 << 4,
};
#endif

#ifndef REDUNDANT_CALCULATION_TITANIA_
#include <RedundantCalculation.hpp>
#endif

#ifndef STRUCTURE_TYPE_TITANIA_
#define STRUCTURE_TYPE_TITANIA_
enum struct StructureInputType : unsigned int
{
  undefined      = 0,
  xyzCoordinates = 1,
  connectivity   = 1 << 2, // 2
  zMatrix        = 1 << 3,
  cosineMatrix   = 1 << 4,
  planarX        = 1 << 5,
  planarY        = 1 << 6,
  planarZ        = 1 << 7,
};

inline StructureInputType
operator|(StructureInputType l, StructureInputType r)
{
  return static_cast<StructureInputType>(static_cast<unsigned int>(l) |
                                         static_cast<unsigned int>(r));
}

inline bool
operator&(StructureInputType l, StructureInputType r)
{
  return static_cast<bool>(static_cast<unsigned int>(l) &
                           static_cast<unsigned int>(r));
}

#endif

#ifndef COORDINATES_TITANIA_
#define COORDINATES_TITANIA_
struct Coordinates {
  double x = .0;
  double y = .0;
  double z = .0;
};
#endif

#ifndef INTERNALCOORDINATE_TITANIA_
#define INTERNALCOORDINATE_TITANIA_
struct InternalCoordinate {
  RedundantType type;
  Atom **participants;
  Potential *potential;
  SphericalHarmonics *harmonic;
  double value;
  unsigned int index;
  InternalCoordinate() :
      type(RedundantType::Undefined), participants(NULL), potential(NULL),
      harmonic(NULL), value(.0), index(0)
  {
    ;
  }
  ~InternalCoordinate()
  {
    if (type == RedundantType::Undefined)
      return;
    else if (type == RedundantType::Distance)
      participants[0] = participants[1] = NULL;
    else if (type == RedundantType::Angle)
      participants[0] = participants[1] = participants[2] = NULL;
    else if (type == RedundantType::Torsion)
      participants[0] = participants[1] = participants[2] = participants[3] =
          NULL;
    else if (type == RedundantType::oop)
      participants[0] = participants[1] = participants[2] = participants[3] =
          NULL;
    potential = NULL;
    harmonic  = NULL;
    free(participants);
  }
};
#endif

#ifndef BOND_ANGLE_TITANIA_
#define BOND_ANGLE_TITANIA_
struct BondAngle {
  Bond *bond1;
  Bond *bond2;
  double angle;
};
#endif

#ifndef Y_PARAMETERS_TITANIA_
#define Y_PARAMETERS_TITANIA_
struct Y_parameters {
  Eigen::MatrixXcd Yref;
  unsigned int rdc;
};
#endif

#ifndef LM_STRUC_MIN_TITANIA_
#define LM_STRUC_MIN_TITANIA_
struct lm_strucmin {
  Eigen::MatrixXd Vecs;
  Eigen::MatrixXd Calc;
  Eigen::MatrixXd Angles;
};
#endif

#ifndef POTENTIAL_TYPE_TITANIA_
#define POTENTIAL_TYPE_TITANIA_
enum struct PotentialType : unsigned int
{
  Undefined   = 0,
  Stretch     = 1,
  Bend        = 1 << 1,
  StretchBend = 1 << 2,
  Dihedral    = 1 << 3,
  vanderWaals = 1 << 4,
  oop         = 1 << 5,
};

inline PotentialType
operator|(PotentialType l, PotentialType r)
{
  return static_cast<PotentialType>(static_cast<unsigned int>(l) |
                                    static_cast<unsigned int>(r));
}

inline PotentialType
operator&(PotentialType l, PotentialType r)
{
  return static_cast<PotentialType>(static_cast<unsigned int>(l) &
                                    static_cast<unsigned int>(r));
}

inline std::ostream &
operator<<(std::ostream &os, PotentialType inp)
{
  switch (((unsigned int) inp))
  {
    case (0):
      os << "undefined potential";
      break;
    case (1):
      os << "stretch potential";
      break;
    case (1 << 1):
      os << "bend potential";
      break;
    case (1 << 2):
      os << "stretch-bend potential";
      break;
    case (1 << 3):
      os << "dihderal potential";
      break;
    case (1 << 4):
      os << "long range (van der Waals) potential";
      break;
    case (1 << 5):
      os << "out of plane potential";
      break;
    default:
      os << "undefined potential";
  }
  return os;
}

#endif

#ifndef BOND_FF_TITANIA_
#define BOND_FF_TITANIA_
struct bond_ff {
  int Bond_type;
  int Atom_type_1;
  int Atom_type_2;
  double k;  /* force constant */
  double r0; /* equilibrium bond length */
};
#endif

#ifndef ANGLE_FF_TITANIA_
#define ANGLE_FF_TITANIA_
struct angle_ff {
  int Bond_type;
  int Atom_type_1;
  int Atom_type_2;
  int Atom_type_3;
  double k;  /* force constant */
  double a0; /* equilibrium bond angle */
};
#endif

#ifndef STRETCH_BEND_FF_TITANIA_
#define STRETCH_BEND_FF_TITANIA_
struct stretch_bend_ff {
  int Bond_type;
  int Atom_type_1; // i
  int Atom_type_2; // j
  int Atom_type_3; // k
  double Fijk;
  double Fkij;
};
#endif

#ifndef TORSION_FF_TITANIA_
#define TORSION_FF_TITANIA_
struct torsion_ff {
  int Bond_type;
  int Atom_type_1;
  int Atom_type_2;
  int Atom_type_3;
  int Atom_type_4;
  double V1; /* dihedral phase */
  double V2; /* force constant */
  double V3; /* dihedral multiplicity */
};
#endif

#ifndef VAN_DER_VAALS_FF_TITANIA_
#define VAN_DER_VAALS_FF_TITANIA_
struct vanderWaals_ff {
  int Atom_type;
  double alpha;
  double N;
  double A;
  double G;
};
#endif

#ifndef BOND_CHARGE_INC_FF_TITANIA_
#define BOND_CHARGE_INC_FF_TITANIA_
struct bond_charge_increment_ff {
  int Bond_type;
  int Atom_type_1;
  int Atom_type_2;
  double bci;
};
#endif

#ifndef OUT_OF_PLANE_FF_TITANIA_
#define OUT_OF_PLANE_FF_TITANIA_
struct out_of_plane_ff {
  int Atom_type_1; // Ligand 1
  int Atom_type_2; // Central atom
  int Atom_type_3; // Ligand 2
  int Atom_type_4; // Ligand 3
  double koop;     // Force constant for oop's:
                   // a( r(2-1) - A(2-3-4) )
                   // a( r(2-3) - A(2-1-4) )
                   // a( r(2-4) - A(2-1-3) )
};
#endif

#ifndef BOND_E_TITANIA_
#define BOND_E_TITANIA_
struct bond_E {
  Atom *Atom_1;
  Atom *Atom_2;
  int index;
  double E_new;
  double E_old;
};
#endif

#ifndef ANGLE_E_TITANIA_
#define ANGLE_E_TITANIA_
struct angle_E {
  Atom *Atom_1;
  Atom *Atom_2;
  Atom *Atom_3;
  int index;
  double E_new;
  double E_old;
};
#endif

#ifndef STRETCH_BEND_E_TITANIA_
#define STRETCH_BEND_E_TITANIA_
struct stretch_bend_E {
  Atom *Atom_1;
  Atom *Atom_2;
  Atom *Atom_3;
  int index_sb;
  int index_b_1;
  int index_b_2;
  int index_a;
  double E_new;
  double E_old;
};
#endif

#ifndef TORSION_E_TITANIA_
#define TORSION_E_TITANIA_
struct torsion_E {
  Atom *Atom_1;
  Atom *Atom_2;
  Atom *Atom_3;
  Atom *Atom_4;
  int index;
  double E_new;
  double E_old;
};
#endif

#ifndef MC_RUN_INFORMATION_TITANIA_
#define MC_RUN_INFORMATION_TITANIA_
struct MC_Run_Information {
  unsigned int max_steps;
  unsigned int steps;
  double convergence;
  bool console_output;
};
#endif

#ifndef MC_STOP_TITANIA_
#define MC_STOP_TITANIA_
enum struct MC_stop : unsigned int
{
  undefined = 0,

  pMean = 1,
  pSigm = 1 << 1,
  RMean = 1 << 2,

  AMean = 1 << 3,
  ASigm = 1 << 4,

  QFac = 1 << 5,

  MaxIter = 1 << 6,
};

inline MC_stop
operator|(MC_stop l, MC_stop r)
{
  return static_cast<MC_stop>(static_cast<unsigned int>(l) |
                              static_cast<unsigned int>(r));
}

inline bool
operator&(MC_stop l, MC_stop r)
{
  return static_cast<bool>(static_cast<unsigned int>(l) &
                           static_cast<unsigned int>(r));
}

inline std::ostream &
operator<<(std::ostream &os, MC_stop &inp)
{
  switch (((unsigned int) inp))
  {
    case (0):
      os << "undefined";
      break;
    case (1):
      os << "change of mean spherical coordinates";
      break;
    case (1 << 1):
      os << "chance of spherical coordinates standard deviation";
      break;
    case (1 << 2):
      os << "change circular distribution R^2";
      break;
    case (1 << 3):
      os << "change of mean alignment tensor";
      break;
    case (1 << 4):
      os << "change of alignment tensor standard deviation";
      break;
    case (1 << 5):
      os << " change of Q-factor";
      break;
    case (1 << 6):
      os << "max TITANIA iterations";
      break;
    default:
      os << "undefined option";
  }
  return os;
}

#endif

#ifndef MC_OUTPUT_TITANIA_
#include <MC_Output.hpp>
#endif

#ifndef SECONDA_OUTPUT_TITANIA_
#define SECONDA_OUTPUT_TITANIA_
struct SECONDA_Output {
  Eigen::VectorXd kappa_q;
  Eigen::VectorXd Covariance_Evals;
  Eigen::MatrixXd Covariance_Evecs;
  Eigen::VectorXd cumulative_variance;
  Eigen::VectorXd Tolman_singular_values;
  Eigen::MatrixXd RDC_R_2;
  Eigen::MatrixXd RDC_m;
  Eigen::MatrixXd RDC_b;
  Eigen::MatrixXd twoSetSingularValues;
  Eigen::VectorXd heterogeneity;
  double rho_5_6;
  double rho_mean_6;
  double covar_condition_number;
  double Tolman_condition_number;
  unsigned int covar_rank;
  unsigned int Tolman_rank;
  unsigned int Tolman_CN_rank;
  bool collected;
  SECONDA_Output()
  {
    rho_5_6 = rho_mean_6   = .0;
    covar_condition_number = Tolman_condition_number = .0;
    covar_rank = Tolman_rank = Tolman_CN_rank = 0;
    collected                                 = false;
  }
};
#endif

/* Define Structure inputs for better readability for humans... */

#ifndef LIMITS_TITANIA_
#define LIMITS_TITANIA_
struct Limits {
  int max_lm_iterations;      // User defined number of max iterations per lm
                              // optimization cylce.
  int max_titania_iterations; // User defined number of max titania iterations
  int max_redundant_cycles;
  double Q_factor_convergence; // Difference in Q between to interations to stop
                               // the optimization
  double redundants_convergence;
  double redundants_validity;
  double alignment_mean_convergence;
  double alignment_sigm_convergence;
  double sphericals_mean_convergence;
  double sphericals_sigm_convergence;
  double sphericals_spread_convergence;
  double zero_cutoff;
  double *phobos_opts_sphericals;
  double *phobos_opts_eckart;
};

#endif


#ifndef STOP_CRITERIA_TITANIA_
#define STOP_CRITERIA_TITANIA_
struct StopCriteria {
  double Q_factor_convergence;
  double redundants_convergence;
  double alignment_mean_convergence;
  double alignment_sigm_convergence;
  double sphericals_mean_convergence;
  double sphericals_sigm_convergence;
  double sphericals_spread_convergence;
};

#endif


union tmp_data {
  int int_;
  unsigned int unsigned_int_;
  float float_;
  double double_;
  char char_;
};

#ifndef BASIC_INFORMATION_TITANIA_
#define BASIC_INFORMATION_TITANIA_
struct BasicInformation {
  tmp_data tmp_information;
  int numberOfInputs;
  int numberOfStructures; /* Right now default = 1. Optimized Structures are
                             neglected... TODO */
  //   int LMDeterminationSteps;                 /* User defined criterion the
  //   maximum amount of LM optimizations */
  int numOfThreads; /* Requested amount of cores by user to run programm on */
  int numOfThreadsUsed;   /* Used amount of threads for TITANIA */
  int maxThreads;         /* Max amount of cores on the used system */
  int predictRDCs;        /* How to predict rdcs */
  int numberOfAlignments; /* Number of sets to be predicted */
  int numberOfFills;      /* Number of sets from random alignments */
  int state;
  int nestingDepth;
  int phobos_worksize_Sphericals;
  int phobos_worksize_Eckart;
  int numOfOptSteps;
  int full_mc_steps;
  int overOptimization;
  int useRedundantsOnlyAfter;
  int memoryPurge;
  int lowerInversionAfter;
  int redundants_distance_optimization;
  int over_titania_iterations;
  unsigned int digits;
  unsigned int NumberOfSets;
  unsigned int NumberOfRDCs;
  unsigned int NumberOfAtoms;
  unsigned int NumberOfChiralCenters;
  double normFactRDC; /* Normalizatoni factor based on user settings */
  double MCvariation; /* Difference in Q between to interations to stop the
                         optimization */
  double numericalDeltaX;
  double redundants_damping;
  double timeStamp;
  double fulltime; /* Saves the times used for full program execution */
  double SCRMtime; /* Saves the times used only for SCRM cycles */
  double structureTime;
  double MCtime;     /* Saves the times used only for Monte Carlo SCRM cycles */
  double MMFF94time; /* Saves the times used only for MMFF94 structure
                        optimization */
  double Systemtime; /* Everything else */
  double inputParseTime; /* Saves the times used only for input parsing */
  double static_redundants_weighting[NUMBER_OF_REDUNDANT_TYPES_];
  double **phobos_workspace_Sphericals;
  double *phobos_workspace_Eckart;
  double **alignments; /* Alignment information for prediction */
  double *q_factors;
  double *chiral_volumes;
  time_t starting_time;
  char hostname[STANDARD_BUFFER_SIZE];
  char username[STANDARD_BUFFER_SIZE];
  char rundir[STANDARD_BUFFER_SIZE];
  FILE *output; /* Every methode should be able to right in output file -> saved
                   in baseInformation */
  FILE *trjFile;
  std::fstream debugFile; /* Every methode should be able to right in debug file
                             -> saved in baseInformation */
  std::fstream xyzFile; /* Every methode should be able to right in xyz file ->
                           saved in baseInformation */
  std::fstream comFile; /* If TITANIA is set to silent mode the communication is
                           printed to a file */
  std::ostream *console;        /* Used for the output of communication */
  std::string workingDirectory; /* Separate the directory of the input file from
                                   the filename */
  std::string inputFileName;    /* Defines the name of input file */
  std::string outputFileName;   /* Defines the name of output file */
  std::string hotFCHTbase;
  std::string xyzFileName;   /* Defines the name of xyz structure file */
  std::string debugFileName; /* File name of the debug file */
  std::string comFileName;   /* File name of the communication file */
  std::string trjFileName;
  std::string errorKey; /* String to save the error message */
  std::string gnuOutputFormat;
  std::string structureLabelExtension;
  std::string runtime_changes;
  std::string gpu_device_name;
  StructureInputType
      structureInput; /* Saves the primary (first) type of structure input */
  StructureInputType
      secondaryInput; /* Saves the second type of structure input (e.g. one can
                         define xyz and connectivity) */
  Output_Format output_format;
  MC_Output
      MonteCarloOutput; /* Saves all results of the Monte Carlo Bootstrapping */
  SECONDA_Output SECONDA_output_initial;
  SECONDA_Output SECONDA_output_final;
  Limits limits;
  StopCriteria stop_crit;
  MC_stop MC_stop_reason;
};
#endif

/* Additional Flags to basicInformation */
#ifndef FLAGS_TITANIA_
#define FLAGS_TITANIA_
struct Flags {
  bool echo;

  bool skipDIDC;
  bool skipSCRM;
  bool skipEckart;

  bool scaleDmatrix;
  bool normDmatrix;
  bool scaleWithSoverall;

  bool noInput;  /* Should be true right after start. If not programm terminates
                  */
  bool noOutput; /* If false a default output file is generated based on
                    inputfilename.out */
  bool appendOutput; /* If an output is defined one can append the output to it
                        or truncate it */
  bool debug;        /* Possibility to generate a debug file containing way more
                        information TODO */
  bool silent;
  bool print_redundants;
  bool planar_input;
  bool printWarnings;
  bool titania2titania;
  bool bigData;
  bool monteCarloOutput;

  bool titania2hotfcht;
  bool outputAli;
  bool outputLM;
  bool minimalOutput;

  bool useRedundants;
  bool redundants_damping;
  bool monteCarloBootstrapping;
  bool calculateFullMatrix;
  bool errorWeightInSVD;
  bool torsions_4_redundants;
  bool long_range_only_4_redundants;
  bool floating_rdc_angles;
  bool use_initial_holonomics;

  bool weightMonteCarlo;
  bool numericalGradients;
  bool recalculateRDCs;
  bool normChiralVolume;
  bool SECONDAreducedCovariance;

  bool plotKappaQ;
  bool plotTrajectory;
  bool plotRDCdynamics;
  bool plotMonteCarlo;
  bool plotRDCrmsd;

  bool converged;

  bool use_gpu;
  bool has_gpu;
};
#endif

#ifndef TITANIA_CONFIGURATION_
#define TITANIA_CONFIGURATION_
struct Mars_Configuration {
  std::string titania_directory;
};
#endif

#endif
