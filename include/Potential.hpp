
#ifndef POTENTIAL_HPP_
#define POTENTIAL_HPP_


#ifndef STRUCTUREOPTIONS_TITANIA_
enum struct StructureOptions : unsigned int;
#endif

#ifndef ATOM_HPP_
class Atom;
#endif

#ifndef STRUCTURE_HPP_
class Structure;
#endif

#ifndef POTENTIAL_TYPE_TITANIA_
enum struct PotentialType : unsigned int;
#endif

#ifndef BOND_FF_TITANIA_
struct bond_ff;
#endif

#ifndef ANGLE_FF_TITANIA_
struct angle_ff;
#endif

#ifndef STRETCH_BEND_FF_TITANIA_
struct stretch_bend_ff;
#endif

#ifndef TORSION_FF_TITANIA_
struct torsion_ff;
#endif

#ifndef VAN_DER_VAALS_FF_TITANIA_
struct vanderWaals_ff;
#endif

#ifndef OUT_OF_PLANE_FF_TITANIA_
struct out_of_plane_ff;
#endif

class Potential {
 private:
  PotentialType type;
  unsigned int index; // Save the index of the respective
                      // MMFF94 entry according to the file
                      // mmff94.ff
  unsigned int NumberOfAtoms;

  // Try to cover full force field
  // in universal variables
  double eq;
  double eq2;
  double eq3;
  double V1;
  double V2;
  double V3;

  double currentState;  // Save current r/theta/phi
  double currentEnergy; // Save current potential

  Atom **atoms; // Participating atoms

  Potential *next;
  Potential *prev;

 public:
  Potential();
  Potential(Potential *,
            PotentialType,
            Atom *,
            Atom *,
            Atom *,
            Atom *,
            unsigned int,
            void *,
            void *b_ff_1 = NULL,
            void *b_ff_2 = NULL,
            void *a_ff   = NULL);
  ~Potential();

  Potential *addNext(PotentialType,
                     Atom *,
                     Atom *,
                     Atom *,
                     Atom *,
                     unsigned int,
                     void *,
                     void *b_ff_1 = NULL,
                     void *b_ff_2 = NULL,
                     void *a_ff   = NULL);
  Potential *getNext();
  Potential *getPrev()
  {
    return prev;
  }
  Potential *getNext(PotentialType);
  Potential *getNext(Atom *);
  Potential *getNext(Atom *, PotentialType);

  void load_stretch_E(bond_ff *, Atom *, Atom *);
  void load_bend_E(angle_ff *, Atom *, Atom *, Atom *);
  void load_stretch_bend_E(stretch_bend_ff *,
                           Atom *,
                           Atom *,
                           Atom *,
                           bond_ff *,
                           bond_ff *,
                           angle_ff *);
  void load_dihedral_E(torsion_ff *, Atom *, Atom *, Atom *, Atom *);
  void load_van_der_waals_E(vanderWaals_ff *,
                            vanderWaals_ff *,
                            Atom *,
                            Atom *,
                            double *);
  void load_out_of_plane_E(out_of_plane_ff *, Atom *, Atom *, Atom *, Atom *);

  double collectEnergys();
  double getEnergy(StructureOptions);
  double get_stretch_E(StructureOptions);
  double get_bend_E(StructureOptions);
  double get_stretch_bend_E(StructureOptions);
  double get_dihedral_E(StructureOptions);
  double get_van_der_waals_E(StructureOptions);
  double get_equilibrium()
  {
    return eq;
  }
  double get_value(StructureOptions);
  double get_distance(StructureOptions);
  double get_angle(StructureOptions);
  double get_dihedral(StructureOptions);
  void set_initial_holonomics();
  void resetStructure(Structure *);

  unsigned int getIndex()
  {
    return index;
  }
  Atom *getAtom(unsigned int i)
  {
    if (atoms && i < NumberOfAtoms)
      return atoms[i];
    else
      return NULL;
  }
  PotentialType getType()
  {
    return type;
  }
  void print();
};

#endif
