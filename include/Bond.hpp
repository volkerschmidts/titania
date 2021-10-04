
#ifndef BOND_HPP_
#define BOND_HPP_

#ifndef FLAGS_TITANIA_
struct Flags;
#endif

#ifndef MOLECULE_HPP_
class Molecule;
#endif

#ifndef ATOM_HPP_
class Atom;
#endif

#ifndef SPHERICALHARMONICS_HPP_
class SphericalHarmonics;
#endif

class Bond {
 private:
  unsigned int order;           // Bond order
  double length;                // Bondlength
  Atom *atom1;                  // Atom 1 participating at this bond
  Atom *atom2;                  // Atom 2 participating at this bond
  SphericalHarmonics *harmonic; // Spherical harmonic which describes this bond
                                // (only if this is known by rdcs)

 public:
  Bond();
  Bond(Atom *, Atom *, SphericalHarmonics * = NULL);
  ~Bond();

  Atom *getAtom1() const
  {
    return atom1;
  }
  Atom *getAtom2() const
  {
    return atom2;
  }
  Atom *getBondpartner(const Atom *a) const
  {
    return (a == atom1 ? atom2 : atom1);
  }
  Atom *getBondpartner(Atom *a) const
  {
    return (a == atom1 ? atom2 : atom1);
  }

  void setHarmonic(SphericalHarmonics *s)
  {
    harmonic = s;
  }
  SphericalHarmonics *getHarmonic() const
  {
    return harmonic;
  };

  int setLength(double l)
  {
    if (l > .0)
    {
      length = l;
      return 0;
    }
    else
      return 1;
  }
  double getLength() const
  {
    return length;
  }

  int setOrder(unsigned int o)
  {
    order = o;
    return 0;
  }
  int setOrder(int o)
  {
    if (o > 0)
    {
      order = (unsigned int) o;
      return 0;
    }
    else
      return 1;
  }
  unsigned int getOrder() const
  {
    return order;
  }
  void copyOrder();
};

void setupBonds(Molecule &, Flags &flags);

#endif
