
#ifndef STRUCTURE_SIMULATOR_HPP_
#define STRUCTURE_SIMULATOR_HPP_

#include <eigen3/Eigen/Core>
#include <fstream>
#include <string>
#include <vector>

#ifndef POTENTIAL_TYPE_TITANIA_
enum struct PotentialType : unsigned int;
#endif

#ifndef BASIC_INFORMATION_TITANIA_
struct BasicInformation;
#endif

#ifndef ATOM_HPP_
class Atom;
#endif

#ifndef POTENTIAL_HPP_
class Potential;
#endif

#ifndef STRUCTURE_HPP_
class Structure;
#endif

#ifndef STRUCTUREOPTIONS_TITANIA_
#include <StructureOptions.hpp>
#endif

#ifndef FLAGS_TITANIA_
struct Flags;
#endif

class StructureSimulator {
 private:
  /* Common structure stuff */
  unsigned int NumberToOpt;
  unsigned int NumberOfAtoms;
  unsigned int NumberOfBonds;
  unsigned int NumberOfAngles;
  unsigned int NumberOfTorsions;
  unsigned int NumberOfPlanarCenters;

  Structure *CurrStruc;
  StructureOptions defOpt;

  Potential *HeadPotential;

  /* List of Atoms (to optimize), all bonds / angles / torsions */
  Atom **ListToOpt;

  bool monitoring;
  std::fstream monitoringFile;

 public:
  StructureSimulator();
  StructureSimulator(Structure *,
                     Flags &flags,
                     StructureOptions opt = StructureOptions::Initial);
  ~StructureSimulator();

  void newStructure(Structure *, StructureOptions);

  void setupPartialCharge(Flags &flags);

  /* Accessing */

  unsigned int getNumberOfAtoms()
  {
    return NumberOfAtoms;
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
  unsigned int getNumberOfPlanarCenters()
  {
    return NumberOfPlanarCenters;
  }

  Potential *getPotential(PotentialType);

  Structure *getCurrStruc()
  {
    return CurrStruc;
  }
  /******Energy Calculations******/

  double totalEnergy();

  double bondEnergy();
  double angleEnergy();
  double stretch_bend_Energy();
  double torsionEnergy();
  double lrEnergy();

  double AtomEnergy(Atom *, StructureOptions);

  /******End of Energy Calculations*******/

  /******Angle Optimization******/

  int optimizeAngles(Atom *);

  /******End of Angle Optimization******/

  /******Simplex Stuff******/

  unsigned int startSimplexMinimizer(unsigned int, double, double, bool = true);
  unsigned int startSimplexMinimizer(unsigned int, bool = true);
  unsigned int startSimplexMinimizer(bool = true);

  void setCoordinates(const Eigen::VectorXd &, const bool);

  bool getmonitoring()
  {
    return monitoring;
  }
  void setMonitorFile(
      std::string n); // { monitoringFile.open(n,
                      // std::ios::binary|std::ios::out|std::ios::trunc); }
  void monitorMinimizer(bool m); // { monitoring = m; }

  unsigned int getNumberToOpt() const
  {
    return NumberToOpt;
  }

  void setStanStructureOptions(StructureOptions opt)
  {
    defOpt = opt;
  }

  /******End Simplex Stuff******/

  void printStructure(unsigned int);
};

int getStretchIndex(Atom *, Atom *);
int check4Bend(Potential *, Potential *, Atom **, Atom **, Atom **);
int getBendIndex(Atom *, Atom *, Atom *);
int getStretchBendIndex(Atom *, Atom *, Atom *);
int
check4Dihedral(Potential *, Potential *, Atom **, Atom **, Atom **, Atom **);
int getDihedralIndex(Atom *, Atom *, Atom *, Atom *);
int check4oop(Atom **, Atom **, Atom **, Atom **);
int getoopIndex(Atom *, Atom *, Atom *, Atom *);

void load_lr_Energy(Structure *, Potential **);

void removePlanarity(Structure *, StructureSimulator &, BasicInformation &);

int maxEnergyAtom(Atom *, StructureSimulator &);

typedef struct {
  bool reduced;
  StructureSimulator *Simulator;
  unsigned int count;
} Efunc_par;
double totalEfunc(const std::vector<double> &, std::vector<double> &, void *);

#endif
