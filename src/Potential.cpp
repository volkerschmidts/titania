
#include <Atom.hpp>
#include <Declarations.hpp>
#include <LinAlg.hpp>
#include <Potential.hpp>
#include <Structure.hpp>
#include <iostream>

// Constants for stretching potentials in mmff94
#define STRETCH_BASE_POTENTIAL__      71.9663
#define CUBIC_BASE_STRETCH_CONSTANT__ -2.0
#define CUBIC_STRETCH_CONSTANT__      (7.0 / 3.0)

// Constants for bending potentials in mmff94
#define BEND_BASE_POTENTIAL__      0.021922
#define CUBIC_BASE_BEND_CONSTANT__ -0.007

// Constants for stretch-bending potentials in mmff94
#define STRETCH_BEND_BASE_POTENTIAL__ 2.51210

// Constants for van-der-Waals potentials in mmff94
#define CUTOFF                            7.0
#define ELECTROSTATIC_BUFFERING_CONSTANT_ 0.05
#define DIELECTRIC_CONSTANT_              1.0
#define DISTANCE_DEPENDENT_EXPONENT_      1.0
#define ELECTROSTATIC_BASE_CONSTANT_      322.0716


Potential::Potential()
{
  index         = 0;
  NumberOfAtoms = 0;

  eq  = .0;
  eq2 = .0;
  eq3 = .0;
  V1  = .0;
  V2  = .0;
  V3  = .0;

  currentState  = .0;
  currentEnergy = .0;

  atoms = NULL;
  next  = NULL;
  prev  = NULL;
  type  = PotentialType::Undefined;
}

Potential::Potential(Potential *p,
                     PotentialType t,
                     Atom *a1,
                     Atom *a2,
                     Atom *a3,
                     Atom *a4,
                     unsigned int i,
                     void *ff,
                     void *b_ff_1,
                     void *b_ff_2,
                     void *a_ff)
{
  // set standard values
  eq = eq2 = eq3 = V1 = V2 = V3 = currentState = currentEnergy = .0;
  atoms                                                        = NULL;
  index                                                        = i;
  type                                                         = t;

  switch (type)
  {
    case PotentialType::Undefined: {
      NumberOfAtoms = 0;
      break;
    }
    case PotentialType::Stretch: {
      NumberOfAtoms        = 2;
      struct bond_ff *dptr = (struct bond_ff *) ff;
      load_stretch_E(dptr, a1, a2);
      break;
    }
    case PotentialType::Bend: {
      NumberOfAtoms         = 3;
      struct angle_ff *dptr = (struct angle_ff *) ff;
      load_bend_E(dptr, a1, a2, a3);
      break;
    }
    case PotentialType::StretchBend: {
      NumberOfAtoms                = 3;
      struct stretch_bend_ff *dptr = (struct stretch_bend_ff *) ff;
      load_stretch_bend_E(dptr, a1, a2, a3, (struct bond_ff *) b_ff_1,
                          (struct bond_ff *) b_ff_2, (struct angle_ff *) a_ff);
      break;
    }
    case PotentialType::Dihedral: {
      NumberOfAtoms           = 4;
      struct torsion_ff *dptr = (struct torsion_ff *) ff;
      load_dihedral_E(dptr, a1, a2, a3, a4);
      break;
    }
    case PotentialType::vanderWaals: {
      NumberOfAtoms = 2;
      load_van_der_waals_E((struct vanderWaals_ff *) ff,
                           (struct vanderWaals_ff *) b_ff_1, a1, a2,
                           (double *) b_ff_2);
      break;
    }
    case PotentialType::oop: {
      NumberOfAtoms = 4;
      load_out_of_plane_E((struct out_of_plane_ff *) ff, a1, a2, a3, a4);
      break;
    }
    default:
      break;
  }

  next = NULL;
  prev = p;
}

Potential::~Potential()
{
  for (unsigned int a = 0; a < NumberOfAtoms; ++a)
    atoms[a] = NULL;
  next = NULL;
  prev = NULL;
  free(atoms);
}

Potential *
Potential::addNext(PotentialType t,
                   Atom *a1,
                   Atom *a2,
                   Atom *a3,
                   Atom *a4,
                   unsigned int i,
                   void *ff,
                   void *b_ff_1,
                   void *b_ff_2,
                   void *a_ff)
{
  next = new Potential(this, t, a1, a2, a3, a4, i, ff, b_ff_1, b_ff_2, a_ff);
  return next;
}

Potential *
Potential::getNext()
{
  return next;
}

Potential *
Potential::getNext(PotentialType t)
{
  Potential *tmp = next;
  while (tmp && tmp->type != t)
  {
    tmp = tmp->next;
  }
  return tmp;
}

Potential *
Potential::getNext(Atom *A)
{
  Potential *tmp = next;
  unsigned int a;
  while (tmp)
  {
    for (a = 0; a < NumberOfAtoms; ++a)
    {
      if (atoms[a] == A)
        return tmp;
    }
    tmp = tmp->getNext();
  }
  return NULL;
}

Potential *
Potential::getNext(Atom *A, PotentialType pt)
{
  Potential *tmp = next;
  unsigned int a;
  while (tmp)
  {
    for (a = 0; (tmp->type == pt) && (a < tmp->NumberOfAtoms); ++a)
    {
      if (tmp->atoms[a] == A)
        return tmp;
    }
    tmp = tmp->getNext(pt);
  }
  return NULL;
}

void
Potential::load_stretch_E(bond_ff *ff, Atom *A1, Atom *A2)
{
  if (ff->r0)
    eq = ff->r0;
  else
  {
    eq = A1->getDistance(A2, StructureOptions::Initial);
  }
  V1       = ff->k;
  atoms    = (Atom **) malloc(NumberOfAtoms * sizeof(Atom *));
  atoms[0] = A1;
  atoms[1] = A2;
}

void
Potential::load_bend_E(angle_ff *ff, Atom *A1, Atom *A2, Atom *A3)
{
  if (ff->a0)
    eq = ff->a0;
  else
  {
    eq = rad2deg(A1->getAngle(A2, A3, StructureOptions::Initial));
  }
  V1       = ff->k;
  atoms    = (Atom **) malloc(NumberOfAtoms * sizeof(Atom *));
  atoms[0] = A1;
  atoms[1] = A2;
  atoms[2] = A3;
}

void
Potential::load_stretch_bend_E(stretch_bend_ff *ff,
                               Atom *A1,
                               Atom *A2,
                               Atom *A3,
                               bond_ff *b_ff_1,
                               bond_ff *b_ff_2,
                               angle_ff *a_ff)
{
  V1 = ff->Fijk; // F123
  V2 = ff->Fkij; // F321

  if (b_ff_1->r0)
    eq2 = b_ff_1->r0;
  else
  {
    eq2 = A1->getDistance(A2, StructureOptions::Initial);
  }

  if (b_ff_2->r0)
    eq3 = b_ff_2->r0;
  else
  {
    eq3 = A3->getDistance(A2, StructureOptions::Initial);
  }

  if (a_ff->a0)
    eq = a_ff->a0;
  else
  {
    eq = rad2deg(A1->getAngle(A2, A3, StructureOptions::Initial));
  }

  atoms    = (Atom **) malloc(NumberOfAtoms * sizeof(Atom *));
  atoms[0] = A1;
  atoms[1] = A2;
  atoms[2] = A3;
}

void
Potential::load_dihedral_E(torsion_ff *ff,
                           Atom *A1,
                           Atom *A2,
                           Atom *A3,
                           Atom *A4)
{
  V1       = ff->V1;
  V2       = ff->V2;
  V3       = ff->V3;
  atoms    = (Atom **) malloc(NumberOfAtoms * sizeof(Atom *));
  atoms[0] = A1;
  atoms[1] = A2;
  atoms[2] = A3;
  atoms[3] = A4;
}

void
Potential::load_van_der_waals_E(vanderWaals_ff *ff1,
                                vanderWaals_ff *ff2,
                                Atom *A1,
                                Atom *A2,
                                double *int14)
{
  double Aa, AN, AA, AG, Aa2, AN2, AA2, AG2, R11, R22, gamma12;

  // Load all values
  Aa  = ff1->alpha;
  Aa2 = ff2->alpha;
  AN  = ff1->N;
  AN2 = ff2->N;
  AA  = ff1->A;
  AA2 = ff2->A;
  AG  = ff1->G;
  AG2 = ff2->G;

  eq2 = A1->getPartialCharge();
  eq3 = A2->getPartialCharge();

  R11 = AA * pow(Aa, 0.25);
  R22 = AA2 * pow(Aa2, 0.25);

  gamma12 = (R11 - R22) / (R11 + R22);

  // Save R12
  eq = 0.5 * (R11 + R22) * (1.0 + 0.2 * (1.0 - exp(-12.0 * gamma12 * gamma12)));
  V3 = pow(eq, 7.0);
  // Save epsilon12
  V1 = (181.16 * AG * AG2 * Aa * Aa2) /
       ((sqrt(Aa / AN) + sqrt(Aa2 / AN2)) * pow(eq, 6.0));

  V2 = int14[0];

  atoms    = (Atom **) malloc(NumberOfAtoms * sizeof(Atom *));
  atoms[0] = A1;
  atoms[1] = A2;
}

void
Potential::load_out_of_plane_E(out_of_plane_ff *ff,
                               Atom *A1,
                               Atom *A2,
                               Atom *A3,
                               Atom *A4)
{
  atoms    = (Atom **) malloc(NumberOfAtoms * sizeof(Atom *));
  atoms[0] = A1; //
  atoms[1] = A2; //  LIGANDS
  atoms[2] = A3; //
  atoms[3] = A4; // Central atom!
  V1       = ff->koop;
}

double
Potential::collectEnergys()
{
  Potential *CurrP = this;
  double totalE    = .0;
  while (CurrP->prev)
    CurrP = CurrP->prev;
  while (CurrP)
  {
    totalE += CurrP->currentEnergy;
    CurrP = CurrP->next;
  }
  return totalE;
}

double
Potential::getEnergy(StructureOptions opt)
{
  switch (type)
  {
    case PotentialType::Undefined:
      return .0;
    case PotentialType::Stretch:
      return get_stretch_E(opt);
    case PotentialType::Bend:
      return get_bend_E(opt);
    case PotentialType::StretchBend:
      return get_stretch_bend_E(opt);
    case PotentialType::Dihedral:
      return get_dihedral_E(opt);
    case PotentialType::vanderWaals:
      return get_van_der_waals_E(opt);
    default:
      return .0;
  }
}

double
Potential::get_stretch_E(StructureOptions opt)
{
  double r = get_distance(opt);
  //   if ( r == currentState ) return currentEnergy;
  currentState = r;
  double dr    = r - eq;
  double dr2   = dr * dr;
  double Vc    = STRETCH_BASE_POTENTIAL__ * V1 * dr2;
  double Vl =
      1.0 + CUBIC_BASE_STRETCH_CONSTANT__ * dr + CUBIC_STRETCH_CONSTANT__ * dr2;
  currentEnergy = Vc * Vl;
  return currentEnergy;
}

double
Potential::get_bend_E(StructureOptions opt)
{
  double alpha = get_angle(opt);
  //   if ( alpha == currentState ) return currentEnergy;
  currentState  = alpha;
  double da     = rad2deg(alpha) - eq;
  double Vc     = BEND_BASE_POTENTIAL__ * V1 * da * da;
  double Vl     = 1 + CUBIC_BASE_BEND_CONSTANT__ * da;
  currentEnergy = Vc * Vl;
  return currentEnergy;
}

double
Potential::get_stretch_bend_E(StructureOptions opt)
{
  double R1    = atoms[0]->getDistance(atoms[1], opt);
  double R2    = atoms[2]->getDistance(atoms[1], opt);
  double alpha = rad2deg(get_angle(opt));

  double dr1 = R1 - eq2;
  double dr2 = R2 - eq3;
  double da  = alpha - eq;

  currentEnergy = STRETCH_BEND_BASE_POTENTIAL__ * (V1 * dr1 + V2 * dr2) * da;

  return .0;
  return currentEnergy;
}

double
Potential::get_dihedral_E(StructureOptions opt)
{
  double phi = get_dihedral(opt);

  //   if ( phi == currentState ) return currentEnergy;
  currentState  = phi;
  double c1     = 1.0 + cos(phi);
  double c2     = 1.0 - cos(2.0 * phi);
  double c3     = 1.0 + cos(3.0 * phi);
  currentEnergy = 0.5 * (V1 * c1 + V2 * c2 + V3 * c3);
  return currentEnergy;
}

double
Potential::get_van_der_waals_E(StructureOptions opt)
{
  currentState = atoms[0]->getDistance(atoms[1], opt);
  if (currentState > CUTOFF)
  {
    currentEnergy = .0;
    return .0;
  }

  // !!!!!!!!!!!!!!!!!!!! //
  // van der Waals Term : //
  // !!!!!!!!!!!!!!!!!!!! //

  currentEnergy = (V1 * pow((1.07 * eq / (currentState + 0.07 * eq)), 7.0) *
                   (1.12 * V3 / (pow(currentState, 7.0) + 0.12 * V3) - 2.0));

  // !!!!!!!!!!!!!!!!!!!!! //
  // electro static Term : //
  // !!!!!!!!!!!!!!!!!!!!! //

  currentEnergy += (V2 * ELECTROSTATIC_BASE_CONSTANT_ * eq2 * eq3 /
                    (DIELECTRIC_CONSTANT_ *
                     pow((currentState + ELECTROSTATIC_BUFFERING_CONSTANT_),
                         DISTANCE_DEPENDENT_EXPONENT_)));
  return currentEnergy;
}

double
Potential::get_value(StructureOptions opt)
{
  switch (type)
  {
    case PotentialType::Undefined:
      return .0;
    case PotentialType::Stretch:
      return get_distance(opt);
    case PotentialType::Bend:
      return get_angle(opt);
    case PotentialType::Dihedral:
      return get_dihedral(opt);
    case PotentialType::vanderWaals:
      return get_distance(opt);
    default:
      return .0;
  }
}

double
Potential::get_distance(StructureOptions opt)
{
  return atoms[0]->getDistance(atoms[1], opt);
}

double
Potential::get_angle(StructureOptions opt)
{
  return (atoms[0]->getAngle(atoms[1], atoms[2], opt));
}

double
Potential::get_dihedral(StructureOptions opt)
{
  return calculateDihedralAngle(atoms[0], atoms[1], atoms[2], atoms[3], opt);
}

void
Potential::set_initial_holonomics()
{
  switch (type)
  {
    case PotentialType::Undefined:
      break;
    case PotentialType::Stretch:
      eq = get_distance(StructureOptions::Initial);
      break;
    case PotentialType::Bend:
      eq = rad2deg(get_angle(StructureOptions::Initial));
      break;
    case PotentialType::StretchBend:
      eq  = rad2deg(get_angle(StructureOptions::Initial));
      eq2 = atoms[0]->getDistance(atoms[1], StructureOptions::Initial);
      eq3 = atoms[2]->getDistance(atoms[1], StructureOptions::Initial);
      break;
    default:
      break;
  }
}

void
Potential::resetStructure(Structure *s)
{
  Atom *CurrAtom = s->getHeadAtom();
  for (unsigned int i = 0; i < NumberOfAtoms; ++i)
  {
    CurrAtom = s->getHeadAtom();
    atomByrdcIndex(&CurrAtom, atoms[i]->getrdcIndex());
    atoms[i] = CurrAtom;
  }
  CurrAtom = NULL;
}

void
Potential::print()
{
  Potential *tmp = this;
  unsigned int a = 0;
  while (tmp)
  {
    std::cout << a++ << "  ";
    for (unsigned int i = 0; i < tmp->NumberOfAtoms; ++i)
      std::cout << tmp->atoms[i]->getIdentifier() << "  ";
    switch (tmp->type)
    {
      case PotentialType::Undefined:
        std::cout << "undefined\n";
        break;
      case PotentialType::Stretch:
        std::cout << "stretch\n";
        break;
      case PotentialType::Bend:
        std::cout << "bend\n";
        break;
      case PotentialType::StretchBend:
        std::cout << "stretch-bend\n";
        break;
      case PotentialType::Dihedral:
        std::cout << "dihedral\n";
        break;
      case PotentialType::vanderWaals:
        std::cout << "vdw\n";
        break;
      default:
        break;
    }
    tmp = tmp->getNext();
  }
}
