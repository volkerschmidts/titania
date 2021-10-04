#ifndef INCLUDE_FRAGMENT_HPP_
#define INCLUDE_FRAGMENT_HPP_

#include <vector>

#ifndef INTERNALCOORDINATE_TITANIA_
struct InternalCoordinate;
#endif

#ifndef POTENTIAL_TYPE_TITANIA_
enum struct PotentialType : unsigned int;
#endif

#ifndef REDUNANDTTYPE_TITANIA_
enum struct RedundantType : unsigned int;
#endif

#ifndef ATOM_HPP_
class Atom;
#endif

#ifndef FRAGMENT_FLAG_TITANIA_
#include <FragmentFlag.hpp>
#endif

#ifndef STRUCTUREOPTIONS_TITANIA_
enum struct StructureOptions : unsigned int;
#endif

#ifndef POTENTIAL_HPP_
class Potential;
#endif // POTENTIAL_HPP_

#ifndef SPHERICALHARMONICS_HPP_
class SphericalHarmonics;
#endif

#ifndef BASIC_INFORMATION_TITANIA_
struct BasicInformation;
#endif

class Fragment {
 private:
  unsigned int NumberOfAtoms;
  unsigned int NumberOfCoordinates;

  unsigned int NumberOfBonds;
  unsigned int NumberOfAngles;
  unsigned int NumberOfTorsions;
  unsigned int NumberOfRDCs;
  unsigned int NumberOfDistances;
  unsigned int NumberOfRedundants;

  double *initial_Internal_Coordinates;
  double *current_Internal_Coordinates;
  double *optimal_Internal_Coordinates;
  double *difference_Internal_Coordinates;
  double *current_Cartesian_Coordinates;
  double *difference_Cartesian_Coordinates;
  double *Wilson_Matrix;
  double *inverse_Wilson_Matrix;

  bool use_torsion_angles;
  bool was_combined;

  Fragment *nextFragment;
  Fragment *prevFragment;
  Atom **list_of_Atoms;
  InternalCoordinate *InternalsList;

  std::vector<Atom *> current_Fragment;
  std::vector<Atom *> to_check_Atoms;
  std::vector<Atom *> optimization_Atoms;

  FragmentFlag fragment_type;

 protected: // define as protected to make them available for cpp unit tests
  static FragmentFlag check_Fragment_type(Atom *);
  static unsigned int build_Fragment(Fragment *,
                                     std::vector<Atom *> &,
                                     std::vector<Atom *> &,
                                     FragmentFlag);
  static Fragment *is_in_Fragment(Fragment *, Atom *);
  static bool is_in_Fragment(Fragment *, Fragment *, Atom *);
  static void add_small_subfragment(Fragment *, std::vector<Atom *> &);
  static void relink_Fragments(Fragment *, Fragment *);

  bool check_Potential(Potential *, bool, bool);
  static bool check_Stretch(Fragment *, Potential *, bool, bool);
  static bool check_Bend(Fragment *, Potential *, bool, bool);
  static bool check_Torsion(Fragment *, Potential *, bool, bool);
  bool check_Harmonic(SphericalHarmonics *);

  static double calculate_value(InternalCoordinate *,
                                std::vector<Atom *>,
                                double *,
                                StructureOptions);
  double calculate_value(unsigned int,
                         std::vector<Atom *>,
                         double *,
                         StructureOptions);
  static double calculate_distance(double *, double *);
  static double calculate_angle(double *, double *, double *);
  static double calculate_torsion(double *, double *, double *, double *);
  static double calculate_rdc_angle(double *, double *, double &, double &);

  static double *get_w_vector(double *, double *);
  static double calculate_derivative(InternalCoordinate *,
                                     std::vector<Atom *>,
                                     Atom *,
                                     unsigned int,
                                     double *);
  static double
  calculate_distance_derivative(double *, double *, double *, unsigned int);
  static double calculate_angle_derivative(double *,
                                           double *,
                                           double *,
                                           double *,
                                           unsigned int);
  static double calculate_torsion_derivative(double *,
                                             double *,
                                             double *,
                                             double *,
                                             double *,
                                             unsigned int);
  static double calculate_rdc_angle_derivative(double *,
                                               double *,
                                               double *,
                                               double,
                                               double,
                                               unsigned int);

  static double *get_vector(double *, double *);
  static double *get_unit_vector(double *, double *, double &);

  void rebase_Coordinates();

  static RedundantType Potential_2_Redundant_Type(PotentialType);
  static int NumberOfParticipants(RedundantType);
  void reserve_Memory();

 public:
  Fragment() :
      NumberOfAtoms(0), NumberOfCoordinates(0), NumberOfBonds(0),
      NumberOfAngles(0), NumberOfTorsions(0), NumberOfRDCs(0),
      NumberOfDistances(0), NumberOfRedundants(0),
      initial_Internal_Coordinates(NULL), current_Internal_Coordinates(NULL),
      optimal_Internal_Coordinates(NULL), difference_Internal_Coordinates(NULL),
      current_Cartesian_Coordinates(NULL),
      difference_Cartesian_Coordinates(NULL), Wilson_Matrix(NULL),
      inverse_Wilson_Matrix(NULL), use_torsion_angles(false),
      was_combined(false), nextFragment(NULL), prevFragment(NULL),
      list_of_Atoms(NULL), InternalsList(NULL), current_Fragment(0, NULL),
      to_check_Atoms(0, NULL), optimization_Atoms(0, NULL),
      fragment_type(FragmentFlag::Undefined)
  {
    ;
  };
  Fragment(Fragment *, Atom *);
  ~Fragment();
  Fragment *getNext()
  {
    return nextFragment;
  }
  Fragment *getPrev()
  {
    return prevFragment;
  }

  void perform_optimization_step();

  Atom *getHeadAtom()
  {
    if (current_Fragment.size())
      return current_Fragment.at(0);
    else
      return NULL;
  }

  Fragment *check_for_fusing();
  void load_Potentials(Potential *, SphericalHarmonics *, StructureOptions);
  unsigned int
  save_Potential(Potential *, SphericalHarmonics *, unsigned int &);

  unsigned int getNumberOfAtoms() const
  {
    return NumberOfAtoms;
  }
  FragmentFlag getFragmentFlag() const
  {
    return fragment_type;
  }
  void combine_Fragments(Fragment *);

  void load_Wilson_B();
  //      void printCoordinates (BasicInformation&);
  //      void printJmol();
};

void svdMoorePenroseInverse(double *, int, int, double *, double thresh = 1e-1);


#endif /* INCLUDE_FRAGMENT_HPP_ */
