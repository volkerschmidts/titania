
#include <Atom.hpp>
#include <Declarations.hpp>
#include <Fragment.hpp>
#include <Helper.hpp>
#include <LinAlg.hpp>
#include <Potential.hpp>
#include <SphericalHarmonics.hpp>
#include <algorithm> //contains
#include <cblas.h>
#include <lapacke.h>

Fragment::Fragment(Fragment *prev, Atom *A)
{
  // Set all types of redundant coordinates to 0
  NumberOfCoordinates = NumberOfBonds = NumberOfAngles = NumberOfTorsions =
      NumberOfRDCs = NumberOfDistances = NumberOfRedundants = 0;

  // Set all arrays to NULL
  initial_Internal_Coordinates     = current_Internal_Coordinates =
      optimal_Internal_Coordinates = difference_Internal_Coordinates = NULL;
  current_Cartesian_Coordinates = difference_Cartesian_Coordinates = NULL;
  Wilson_Matrix = inverse_Wilson_Matrix = NULL;
  use_torsion_angles                    = true;
  was_combined                          = false;

  // For later checks set all links to NULL;
  prevFragment  = prev;
  nextFragment  = NULL;
  list_of_Atoms = NULL;
  InternalsList = NULL;

  // Add the first atom to the current fragment list.
  current_Fragment.push_back(A);

  // define the current fragment type
  fragment_type = Fragment::check_Fragment_type(A);

  // Build the current Fragment
  NumberOfAtoms = Fragment::build_Fragment(this, current_Fragment,
                                           to_check_Atoms, fragment_type);


  Fragment *tmp = this;
  for (unsigned int i = 0; i < to_check_Atoms.size(); ++i)
  {
    if (Fragment::is_in_Fragment(tmp, to_check_Atoms.at(i)))
      continue;
    else
      tmp->nextFragment = new Fragment(tmp, to_check_Atoms.at(i));

    while (tmp->nextFragment)
      tmp = tmp->nextFragment;
  }
  tmp = NULL;
}

Fragment::~Fragment()
{
  free(initial_Internal_Coordinates);
  free(current_Internal_Coordinates);
  free(optimal_Internal_Coordinates);
  free(difference_Internal_Coordinates);
  free(current_Cartesian_Coordinates);
  free(Wilson_Matrix);
  free(inverse_Wilson_Matrix);
  free(InternalsList);
  for (unsigned int i = 0; (list_of_Atoms && i < NumberOfAtoms); ++i)
    list_of_Atoms[i] = NULL;
  free(list_of_Atoms);
  current_Fragment.clear();
  to_check_Atoms.clear();
  prevFragment = NULL;
}

FragmentFlag
Fragment::check_Fragment_type(Atom *A)
{
  return (A->getNumberOfHarmonics() ? FragmentFlag::RDC : FragmentFlag::NonRDC);
}

unsigned int
Fragment::build_Fragment(Fragment *CurrFrag,
                         std::vector<Atom *> &current_Fragment_List,
                         std::vector<Atom *> &to_check_Atoms_List,
                         FragmentFlag fragType)
{
  Atom *CurrAtom, *BondAtom;
  unsigned int fragment_index = 0;
  unsigned int fragment_size  = current_Fragment_List.size();
  unsigned int to_check_index;
  unsigned int to_check_size = 0;
  unsigned int bond_index;

  for (; fragment_index < fragment_size; ++fragment_index)
  {
    CurrAtom = current_Fragment_List.at(fragment_index);
    for (bond_index = 0; bond_index < CurrAtom->getNumberOfBonds();
         ++bond_index)
    {
      BondAtom = CurrAtom->getBondpartner(bond_index);

      // Check if Atom is allready in another Fragment
      //         if (
      //         std::find(current_Fragment_List.begin(),current_Fragment_List.end(),BondAtom)
      //         != current_Fragment_List.end() ) continue;
      if (Fragment::is_in_Fragment(CurrFrag, BondAtom))
        continue;

      // Check if Atom is of the right type
      if (Fragment::check_Fragment_type(BondAtom) == fragType)
      {
        //            if ( Fragment::is_in_Fragment( CurrFrag, BondAtom ) )
        //            continue; if (
        //            (std::find(current_Fragment_List.begin(),current_Fragment_List.end(),BondAtom)
        //            == current_Fragment_List.end() ))
        //            {
        current_Fragment_List.push_back(BondAtom);
        fragment_size = current_Fragment_List.size();
        //            }
      }
      else if (fragType == FragmentFlag::RDC && BondAtom->getZ() == 1)
      {
        current_Fragment_List.push_back(BondAtom);
        fragment_size = current_Fragment_List.size();
      }
      else
      {
        if ((std::find(to_check_Atoms_List.begin(), to_check_Atoms_List.end(),
                       BondAtom) == to_check_Atoms_List.end()))
        {
          to_check_Atoms_List.push_back(BondAtom);
          to_check_size = to_check_Atoms_List.size();
        }
      }
    }

    if (fragment_index == (fragment_size - 1))
    {
      for (to_check_index = 0; to_check_index < to_check_size; ++to_check_index)
      {
        CurrAtom = to_check_Atoms_List.at(to_check_index);
        for (bond_index = 0; bond_index < CurrAtom->getNumberOfBonds();
             ++bond_index)
        {
          BondAtom = CurrAtom->getBondpartner(bond_index);
          if (fragType == FragmentFlag::RDC &&
              Fragment::check_Fragment_type(BondAtom) == fragType)
          {
            if ((std::find(current_Fragment_List.begin(),
                           current_Fragment_List.end(), BondAtom) !=
                 current_Fragment_List.end())) // does not contain
            {
              current_Fragment_List.push_back(CurrAtom);
              fragment_size = current_Fragment_List.size();
              to_check_Atoms_List.erase(to_check_Atoms_List.begin() +
                                        to_check_index);
              to_check_size  = to_check_Atoms_List.size();
              to_check_index = to_check_size;
              break;
            }
          }
          else
          {
            if ((std::find(to_check_Atoms_List.begin(),
                           to_check_Atoms_List.end(), BondAtom) !=
                 to_check_Atoms_List.end())) // does contain
            {
              current_Fragment_List.push_back(CurrAtom);
              current_Fragment_List.push_back(BondAtom);
              fragment_size = current_Fragment_List.size();
              to_check_Atoms_List.erase(to_check_Atoms_List.begin() +
                                        to_check_index);
              to_check_size = to_check_Atoms_List.size();
              for (; to_check_index < to_check_size; ++to_check_index)
              {
                if (to_check_Atoms_List.at(to_check_index) == BondAtom)
                {
                  to_check_Atoms_List.erase(to_check_Atoms_List.begin() +
                                            to_check_index);
                  to_check_size  = to_check_Atoms_List.size();
                  to_check_index = to_check_size;
                }
              }
              break;
            }
          }
        }
      }
    }
  }
  BondAtom = CurrAtom = NULL;
  return fragment_size;
}

Fragment *
Fragment::is_in_Fragment(Fragment *prev, Atom *CurrAtom)
{
  while (prev)
  {
    //      if (
    //      std::find(prev->current_Fragment.begin(),prev->current_Fragment.end(),CurrAtom)
    //      != prev->current_Fragment.end() ) return prev;
    if (vector_contains(prev->current_Fragment, CurrAtom))
      return prev;
    prev = prev->prevFragment;
  }
  return NULL;
}

bool
Fragment::is_in_Fragment(Fragment *Frag_1, Fragment *Frag_2, Atom *border)
{
  if (vector_contains(Frag_1->optimization_Atoms, border) &&
      vector_contains(Frag_2->optimization_Atoms, border))
    return true;
  return false;
}

void
Fragment::add_small_subfragment(Fragment *bigFrag,
                                std::vector<Atom *> &to_be_added)
{
  for (unsigned int i = 0; i < to_be_added.size(); ++i)
  {
    bigFrag->current_Fragment.push_back(to_be_added.at(i));
  }
}

void
Fragment::relink_Fragments(Fragment *Head, Fragment *Tail)
{
  if (Head)
    Head->nextFragment = Tail;
  if (Tail)
    Tail->prevFragment = Head;
}

bool
Fragment::check_Potential(Potential *CurrP, bool fullPotential, bool count)
{
  switch (CurrP->getType())
  {
    case (PotentialType::Stretch):
      if (Fragment::check_Stretch(this, CurrP, fullPotential, count))
      {
        if (count)
          ++NumberOfBonds;
        return true;
      }
      return false;
    case (PotentialType::Bend):
      if (Fragment::check_Bend(this, CurrP, fullPotential, count))
      {
        if (count)
          ++NumberOfAngles;
        return true;
      }
      return false;
    case (PotentialType::Dihedral):
      if (use_torsion_angles &&
          Fragment::check_Torsion(this, CurrP, fullPotential, count))
      {
        if (count)
          ++NumberOfTorsions;
        return true;
      }
      return false;
    default:
      return false;
  }
  return false;
}

bool
Fragment::check_Stretch(Fragment *CurrF,
                        Potential *CurrP,
                        bool fullPotential,
                        bool count)
{
  if (!fullPotential &&
      vector_contains(CurrF->current_Fragment, CurrP->getAtom(0)))
  {
    if (count &&
        !(vector_contains(CurrF->current_Fragment, CurrP->getAtom(1)) ||
          vector_contains(CurrF->to_check_Atoms, CurrP->getAtom(1))))
      CurrF->to_check_Atoms.push_back(CurrP->getAtom(1));
    return true;
  }
  else if (!fullPotential &&
           vector_contains(CurrF->current_Fragment, CurrP->getAtom(1)))
  {
    if (count &&
        !(vector_contains(CurrF->current_Fragment, CurrP->getAtom(0)) ||
          vector_contains(CurrF->to_check_Atoms, CurrP->getAtom(0))))
      CurrF->to_check_Atoms.push_back(CurrP->getAtom(0));
    return true;
  }
  else if (fullPotential &&
           vector_contains(CurrF->optimization_Atoms, CurrP->getAtom(0)) &&
           vector_contains(CurrF->optimization_Atoms, CurrP->getAtom(1)))
  {
    return true;
  }

  return false;
}

bool
Fragment::check_Bend(Fragment *CurrF,
                     Potential *CurrP,
                     bool fullPotential,
                     bool count)
{
  if (!fullPotential &&
      vector_contains(CurrF->current_Fragment, CurrP->getAtom(0)))
  {
    if (count &&
        !(vector_contains(CurrF->current_Fragment, CurrP->getAtom(1)) ||
          vector_contains(CurrF->to_check_Atoms, CurrP->getAtom(1))))
      CurrF->to_check_Atoms.push_back(CurrP->getAtom(1));
    if (count &&
        !(vector_contains(CurrF->current_Fragment, CurrP->getAtom(2)) ||
          vector_contains(CurrF->to_check_Atoms, CurrP->getAtom(2))))
      CurrF->to_check_Atoms.push_back(CurrP->getAtom(2));
    return true;
  }
  else if (!fullPotential &&
           vector_contains(CurrF->current_Fragment, CurrP->getAtom(1)))
  {
    if (count &&
        !(vector_contains(CurrF->current_Fragment, CurrP->getAtom(0)) ||
          vector_contains(CurrF->to_check_Atoms, CurrP->getAtom(0))))
      CurrF->to_check_Atoms.push_back(CurrP->getAtom(0));
    if (count &&
        !(vector_contains(CurrF->current_Fragment, CurrP->getAtom(2)) ||
          vector_contains(CurrF->to_check_Atoms, CurrP->getAtom(2))))
      CurrF->to_check_Atoms.push_back(CurrP->getAtom(2));
    return true;
  }
  else if (!fullPotential &&
           vector_contains(CurrF->current_Fragment, CurrP->getAtom(2)))
  {
    if (count &&
        !(vector_contains(CurrF->current_Fragment, CurrP->getAtom(0)) ||
          vector_contains(CurrF->to_check_Atoms, CurrP->getAtom(0))))
      CurrF->to_check_Atoms.push_back(CurrP->getAtom(0));
    if (count &&
        !(vector_contains(CurrF->current_Fragment, CurrP->getAtom(1)) ||
          vector_contains(CurrF->to_check_Atoms, CurrP->getAtom(1))))
      CurrF->to_check_Atoms.push_back(CurrP->getAtom(1));
    return true;
  }
  else if (fullPotential &&
           vector_contains(CurrF->optimization_Atoms, CurrP->getAtom(0)) &&
           vector_contains(CurrF->optimization_Atoms, CurrP->getAtom(1)) &&
           vector_contains(CurrF->optimization_Atoms, CurrP->getAtom(2)))
  {
    return true;
  }
  return false;
}

bool
Fragment::check_Torsion(Fragment *CurrF,
                        Potential *CurrP,
                        bool fullPotential,
                        bool count)
{
  if (!fullPotential &&
      vector_contains(CurrF->current_Fragment, CurrP->getAtom(0)))
  {
    if (count &&
        !(vector_contains(CurrF->current_Fragment, CurrP->getAtom(1)) ||
          vector_contains(CurrF->to_check_Atoms, CurrP->getAtom(1))))
      CurrF->to_check_Atoms.push_back(CurrP->getAtom(1));
    if (count &&
        !(vector_contains(CurrF->current_Fragment, CurrP->getAtom(2)) ||
          vector_contains(CurrF->to_check_Atoms, CurrP->getAtom(2))))
      CurrF->to_check_Atoms.push_back(CurrP->getAtom(2));
    if (count &&
        !(vector_contains(CurrF->current_Fragment, CurrP->getAtom(3)) ||
          vector_contains(CurrF->to_check_Atoms, CurrP->getAtom(3))))
      CurrF->to_check_Atoms.push_back(CurrP->getAtom(3));
    return true;
  }
  else if (!fullPotential &&
           vector_contains(CurrF->current_Fragment, CurrP->getAtom(1)))
  {
    if (count &&
        !(vector_contains(CurrF->current_Fragment, CurrP->getAtom(0)) ||
          vector_contains(CurrF->to_check_Atoms, CurrP->getAtom(0))))
      CurrF->to_check_Atoms.push_back(CurrP->getAtom(0));
    if (count &&
        !(vector_contains(CurrF->current_Fragment, CurrP->getAtom(2)) ||
          vector_contains(CurrF->to_check_Atoms, CurrP->getAtom(2))))
      CurrF->to_check_Atoms.push_back(CurrP->getAtom(2));
    if (count &&
        !(vector_contains(CurrF->current_Fragment, CurrP->getAtom(3)) ||
          vector_contains(CurrF->to_check_Atoms, CurrP->getAtom(3))))
      CurrF->to_check_Atoms.push_back(CurrP->getAtom(3));
    return true;
  }
  else if (!fullPotential &&
           vector_contains(CurrF->current_Fragment, CurrP->getAtom(2)))
  {
    if (count &&
        !(vector_contains(CurrF->current_Fragment, CurrP->getAtom(0)) ||
          vector_contains(CurrF->to_check_Atoms, CurrP->getAtom(0))))
      CurrF->to_check_Atoms.push_back(CurrP->getAtom(0));
    if (count &&
        !(vector_contains(CurrF->current_Fragment, CurrP->getAtom(1)) ||
          vector_contains(CurrF->to_check_Atoms, CurrP->getAtom(1))))
      CurrF->to_check_Atoms.push_back(CurrP->getAtom(1));
    if (count &&
        !(vector_contains(CurrF->current_Fragment, CurrP->getAtom(3)) ||
          vector_contains(CurrF->to_check_Atoms, CurrP->getAtom(3))))
      CurrF->to_check_Atoms.push_back(CurrP->getAtom(3));
    return true;
  }
  else if (!fullPotential &&
           vector_contains(CurrF->current_Fragment, CurrP->getAtom(3)))
  {
    if (count &&
        !(vector_contains(CurrF->current_Fragment, CurrP->getAtom(0)) ||
          vector_contains(CurrF->to_check_Atoms, CurrP->getAtom(0))))
      CurrF->to_check_Atoms.push_back(CurrP->getAtom(0));
    if (count &&
        !(vector_contains(CurrF->current_Fragment, CurrP->getAtom(1)) ||
          vector_contains(CurrF->to_check_Atoms, CurrP->getAtom(1))))
      CurrF->to_check_Atoms.push_back(CurrP->getAtom(1));
    if (count &&
        !(vector_contains(CurrF->current_Fragment, CurrP->getAtom(2)) ||
          vector_contains(CurrF->to_check_Atoms, CurrP->getAtom(2))))
      return true;
    CurrF->to_check_Atoms.push_back(CurrP->getAtom(2));
  }
  else if (fullPotential &&
           vector_contains(CurrF->optimization_Atoms, CurrP->getAtom(0)) &&
           vector_contains(CurrF->optimization_Atoms, CurrP->getAtom(1)) &&
           vector_contains(CurrF->optimization_Atoms, CurrP->getAtom(2)) &&
           vector_contains(CurrF->optimization_Atoms, CurrP->getAtom(3)))
  {
    return true;
  }
  return false;
}

bool
Fragment::check_Harmonic(SphericalHarmonics *CurrY)
{
  if ((vector_contains(current_Fragment, CurrY->getAtom1()) &&
       vector_contains(current_Fragment, CurrY->getAtom2())))
    return true;
  return false;
}

double
Fragment::calculate_value(InternalCoordinate *CurrCoordinate,
                          std::vector<Atom *> atoms,
                          double *CartCoordinates,
                          StructureOptions opt)
{
  switch (CurrCoordinate->type)
  {
    case RedundantType::Undefined:
      return .0;
    case RedundantType::Distance: {
      unsigned int index_a1, index_a2;
      index_a1 = 3 * get_vector_index(atoms, CurrCoordinate->participants[0]);
      index_a2 = 3 * get_vector_index(atoms, CurrCoordinate->participants[1]);
      return Fragment::calculate_distance((CartCoordinates + index_a1),
                                          (CartCoordinates + index_a2));
    }
    case RedundantType::Angle: {
      unsigned int index_a1, index_a2, index_a3;
      index_a1 = 3 * get_vector_index(atoms, CurrCoordinate->participants[0]);
      index_a2 = 3 * get_vector_index(atoms, CurrCoordinate->participants[1]);
      index_a3 = 3 * get_vector_index(atoms, CurrCoordinate->participants[2]);
      return Fragment::calculate_angle((CartCoordinates + index_a1),
                                       (CartCoordinates + index_a2),
                                       (CartCoordinates + index_a3));
    }
    case RedundantType::Torsion: {
      unsigned int index_a1, index_a2, index_a3, index_a4;
      index_a1 = 3 * get_vector_index(atoms, CurrCoordinate->participants[0]);
      index_a2 = 3 * get_vector_index(atoms, CurrCoordinate->participants[1]);
      index_a3 = 3 * get_vector_index(atoms, CurrCoordinate->participants[2]);
      index_a4 = 3 * get_vector_index(atoms, CurrCoordinate->participants[3]);
      return Fragment::calculate_torsion(
          (CartCoordinates + index_a1), (CartCoordinates + index_a2),
          (CartCoordinates + index_a3), (CartCoordinates + index_a4));
    }
    case RedundantType::RDC: {
      unsigned int index_a1, index_a2;
      double theta = CurrCoordinate->harmonic->getTheta(opt);
      double phi   = CurrCoordinate->harmonic->getPhi(opt);
      index_a1 = 3 * get_vector_index(atoms, CurrCoordinate->participants[0]);
      index_a2 = 3 * get_vector_index(atoms, CurrCoordinate->participants[1]);

      return Fragment::calculate_rdc_angle((CartCoordinates + index_a1),
                                           (CartCoordinates + index_a2), theta,
                                           phi);
    }
    default:
      return .0;
  }
  return .0;
}

double
Fragment::calculate_value(unsigned int IntCoord,
                          std::vector<Atom *> atoms,
                          double *CartCoordinates,
                          StructureOptions opt)
{
  return (Fragment::calculate_value((InternalsList + IntCoord), atoms,
                                    CartCoordinates, opt));
}

double
Fragment::calculate_distance(double *n, double *m)
{
  double R;
  double *u = Fragment::get_unit_vector(n, m, R);
  free(u);
  if (R == R)
    return R;
  else
    return .0;
}

double
Fragment::calculate_angle(double *n, double *o, double *m)
{
  double Ru, Rv, a;

  double *u = Fragment::get_unit_vector(o, m, Ru);
  double *v = Fragment::get_unit_vector(o, n, Rv);

  // compute u' * v
  cblas_dgemv(CblasRowMajor, CblasTrans, NUMBER_OF_AXIS_, 1, 1.0, u, 1, v, 1, 0,
              &a, 1);

  a = acos(a);
  free(u);
  free(v);
  if (a == a)
    return a;
  else
    return .0;
}

double
Fragment::calculate_torsion(double *m, double *o, double *p, double *n)
{
  return .0;
  double *u = (double *) malloc(double_vector_size);
  double *v = (double *) malloc(double_vector_size);
  double *w = (double *) malloc(double_vector_size);
  double *uxw, *vxw;

  double Ru, Rv, Rw, udw, vdw, sinUsinV, a;
  Ru = Rv = Rw = udw = vdw = sinUsinV = a = .0;
  // Compute vector u',v',w'
  for (int axis = 0; axis < NUMBER_OF_AXIS_; ++axis)
  {
    u[axis] = (m[axis] - o[axis]);
    Ru += (u[axis] * u[axis]);
    v[axis] = (n[axis] - p[axis]);
    Rv += (v[axis] * v[axis]);
    w[axis] = (p[axis] - o[axis]);
    Rw += (w[axis] * w[axis]);
  }
  Ru = sqrt(Ru);
  Rv = sqrt(Rv);
  Rw = sqrt(Rw);
  // Normalize vector u',v',w' to obtain u,v,w
  // Additionally compute dot products
  for (int axis = 0; axis < NUMBER_OF_AXIS_; ++axis)
  {
    u[axis] /= Ru;
    v[axis] /= Rv;
    w[axis] /= Rw;
    udw += (u[axis] * w[axis]);
    vdw += (v[axis] * w[axis]);
  }

  // Compute the cross products
  uxw = axb3d(u, w);
  vxw = axb3d(v, w);

  // Compute the product of the sin(Bondangles)
  sinUsinV = sqrt(1.0 - (vdw * vdw)) * sqrt(1.0 - (udw * udw));

  // Compute the cos(dihedralangle)
  for (int axis = 0; axis < NUMBER_OF_AXIS_; ++axis)
  {
    a += (uxw[axis] * vxw[axis]);
  }

  a = acos(a / (sinUsinV));
  free(u);
  free(v);
  free(w);
  free(uxw);
  free(vxw);
  if (a == a)
    return a;
  else
    return .0;
}

double
Fragment::calculate_rdc_angle(double *n, double *m, double &theta, double &phi)
{
  double Ru, a;
  double *u = Fragment::get_unit_vector(n, m, Ru);
  double *v = (double *) malloc(double_vector_size);

  v[STANDARD_AXIS_X_] = sin(theta) * cos(phi);
  v[STANDARD_AXIS_Y_] = sin(theta) * sin(phi);
  v[STANDARD_AXIS_Z_] = cos(theta);

  // compute u' * v
  cblas_dgemv(CblasRowMajor, CblasTrans, NUMBER_OF_AXIS_, 1, 1.0, u, 1, v, 1, 0,
              &a, 1);

  a = acos(a);
  free(u);
  free(v);

  if (a == a)
  {
    if (a > PI_HALF_)
    {
      a     = PI_ - a;
      theta = (PI_ - theta);
      phi   = (phi - PI_);
    }
    return (a < 1e-5 ? .0 : a);
  }
  else
    return .0;
}

double *
Fragment::get_w_vector(double *u, double *v)
{
  static double null[3] = {0.0, 0.0, 0.0};
  static double vec1[3] = {1.0, -1.0, 1.0};
  static double vec2[3] = {-1.0, 1.0, 1.0};
  double *w;
  if (Fragment::calculate_angle(u, &(null[0]), v) == .0)
  {
    if (Fragment::calculate_angle(u, &(null[0]), &(vec1[0])) == .0)
    {
      w = axb3d(u, &(vec2[0]));
    }
    else
      w = axb3d(u, &(vec1[0]));
  }
  else
    w = axb3d(u, v);
  return w;
}

double
Fragment::calculate_derivative(InternalCoordinate *CurrCoordinate,
                               std::vector<Atom *> atoms,
                               Atom *atom,
                               unsigned int axis,
                               double *CartCoordinates)
{
  switch (CurrCoordinate->type)
  {
    case RedundantType::Distance: {
      unsigned int index_a1, index_a2, derivative_index;
      index_a1 = 3 * get_vector_index(atoms, CurrCoordinate->participants[0]);
      index_a2 = 3 * get_vector_index(atoms, CurrCoordinate->participants[1]);
      derivative_index = 3 * get_vector_index(atoms, atom);
      return Fragment::calculate_distance_derivative(
          (CartCoordinates + index_a1), (CartCoordinates + index_a2),
          (CartCoordinates + derivative_index), axis);
    }
    case RedundantType::Angle: {
      unsigned int index_a1, index_a2, index_a3, derivative_index;
      index_a1 = 3 * get_vector_index(atoms, CurrCoordinate->participants[0]);
      index_a2 = 3 * get_vector_index(atoms, CurrCoordinate->participants[1]);
      index_a3 = 3 * get_vector_index(atoms, CurrCoordinate->participants[2]);
      derivative_index = 3 * get_vector_index(atoms, atom);
      return Fragment::calculate_angle_derivative(
          (CartCoordinates + index_a1), (CartCoordinates + index_a2),
          (CartCoordinates + index_a3), (CartCoordinates + derivative_index),
          axis);
    }
    case RedundantType::Torsion: {
      unsigned int index_a1, index_a2, index_a3, index_a4, derivative_index;
      index_a1 = 3 * get_vector_index(atoms, CurrCoordinate->participants[0]);
      index_a2 = 3 * get_vector_index(atoms, CurrCoordinate->participants[1]);
      index_a3 = 3 * get_vector_index(atoms, CurrCoordinate->participants[2]);
      index_a4 = 3 * get_vector_index(atoms, CurrCoordinate->participants[3]);
      derivative_index = 3 * get_vector_index(atoms, atom);
      return Fragment::calculate_torsion_derivative(
          (CartCoordinates + index_a1), (CartCoordinates + index_a2),
          (CartCoordinates + index_a3), (CartCoordinates + index_a4),
          (CartCoordinates + derivative_index), axis);
    }
    case RedundantType::RDC: {
      unsigned int index_a1, index_a2, derivative_index;
      double theta =
          CurrCoordinate->harmonic->getTheta(StructureOptions::Optimized);
      double phi =
          CurrCoordinate->harmonic->getPhi(StructureOptions::Optimized);
      index_a1 = 3 * get_vector_index(atoms, CurrCoordinate->participants[0]);
      index_a2 = 3 * get_vector_index(atoms, CurrCoordinate->participants[1]);
      derivative_index = 3 * get_vector_index(atoms, atom);
      return Fragment::calculate_rdc_angle_derivative(
          (CartCoordinates + index_a1), (CartCoordinates + index_a2),
          (CartCoordinates + derivative_index), theta, phi, axis);
    }
    default:
      return .0;
  }
  return .0;
}

double
Fragment::calculate_distance_derivative(double *n,
                                        double *m,
                                        double *a,
                                        unsigned int axis)
{
  double zeta_amn = zetaFunction(a, m, n);
  if (zeta_amn == .0)
    return zeta_amn;
  double R;
  double *u = Fragment::get_unit_vector(n, m, R);
  zeta_amn *= u[axis]; // calcualte derivative
  free(u);

  return zeta_amn; // R(axis)*zeta_amn;
}

double
Fragment::calculate_angle_derivative(double *n,
                                     double *o,
                                     double *m,
                                     double *a,
                                     unsigned int axis)
{
  double zeta_amo = zetaFunction(a, m, o);
  double zeta_ano = zetaFunction(a, n, o);

  if ((zeta_amo == .0) && (zeta_ano == .0))
    return .0;

  double Ru, Rv, Rw, derivative;
  double *w, *uxw, *wxv;
  double *u = Fragment::get_unit_vector(o, m, Ru);
  double *v = Fragment::get_unit_vector(o, n, Rv);

  w  = Fragment::get_w_vector(u, v);
  Rw = cblas_dnrm2(NUMBER_OF_AXIS_, w, 1);        // calcualte ||w||
  cblas_dscal(NUMBER_OF_AXIS_, (1.0 / Rw), w, 1); // normalize w

  uxw = axb3d(u, w);
  wxv = axb3d(w, v);

  derivative = zeta_amo * (uxw[axis] / Ru) + zeta_ano * (wxv[axis] / Rv);
  free(u);
  free(v);
  free(w);
  free(uxw);
  free(wxv);

  if (derivative == derivative)
    return derivative;
  else
    return .0;
}

double
Fragment::calculate_torsion_derivative(double *n,
                                       double *p,
                                       double *o,
                                       double *m,
                                       double *a,
                                       unsigned int axis)
{
  double zeta_amo = zetaFunction(a, m, o);
  double zeta_apn = zetaFunction(a, p, n);
  double zeta_aop = zetaFunction(a, o, p);
  if ((zeta_aop == .0) && (zeta_apn == .0) && (zeta_amo == .0))
    return .0;

  double Ru, Rv, Rw, derivative;
  double *u = Fragment::get_unit_vector(o, m, Ru);
  double *v = Fragment::get_unit_vector(p, n, Rv);
  double *w = Fragment::get_unit_vector(o, p, Rw);

  double *uxw    = axb3d(u, w);
  double *vxw    = axb3d(v, w);
  double phi_u   = Fragment::calculate_angle(m, o, p);
  double phi_v   = Fragment::calculate_angle(o, p, n);
  double sin_u_2 = sin(phi_u) * sin(phi_u);
  double sin_v_2 = sin(phi_v) * sin(phi_v);
  double term_u, term_v, cross_term;

  term_u = zeta_amo * (uxw[axis] / (Ru * sin_u_2));
  term_v = zeta_apn * (vxw[axis] / (Rv * sin_v_2));

  cross_term = zeta_aop * ((uxw[axis] * cos(phi_u) / (Rw * sin_u_2)) -
                           (vxw[axis] * cos(phi_v) / (Rw * sin_v_2)));

  derivative = term_u + term_v + cross_term;

  free(uxw);
  free(vxw);
  free(u);
  free(v);
  free(w);

  if (derivative == derivative)
    return cross_term;
  else
    return .0;
}

double
Fragment::calculate_rdc_angle_derivative(double *n,
                                         double *m,
                                         double *a,
                                         double theta,
                                         double phi,
                                         unsigned int axis)
{
  double zeta_amn = zetaFunction(a, m, n);

  if (zeta_amn == .0)
    return .0;

  Fragment::calculate_rdc_angle(n, m, theta, phi);

  double Ru, Rw, derivative;
  double *w, *uxw;
  double *u = Fragment::get_unit_vector(n, m, Ru);
  double *v = (double *) malloc(double_vector_size);

  // calculate vector v
  v[STANDARD_AXIS_X_] = sin(theta) * cos(phi);
  v[STANDARD_AXIS_Y_] = sin(theta) * sin(phi);
  v[STANDARD_AXIS_Z_] = cos(theta);

  w  = Fragment::get_w_vector(u, v);
  Rw = cblas_dnrm2(NUMBER_OF_AXIS_, w, 1);        // calcualte ||w||
  cblas_dscal(NUMBER_OF_AXIS_, (1.0 / Rw), w, 1); // normalize w

  uxw = axb3d(u, w);

  derivative = zeta_amn * (uxw[axis] / Ru);

  free(u);
  free(v);
  free(w);
  free(uxw);

  if (derivative == derivative)
    return derivative;
  else
    return .0;
}

double *
Fragment::get_vector(double *a, double *b)
{
  double *c = (double *) malloc(double_vector_size);
  cblas_dcopy(NUMBER_OF_AXIS_, b, 1, c, 1);         // copy b in c
  cblas_daxpy(NUMBER_OF_AXIS_, (-1.0), a, 1, c, 1); // calculate u = m - o
  return c;
}

double *
Fragment::get_unit_vector(double *a, double *b, double &R)
{
  double *c = Fragment::get_vector(a, b);
  R         = cblas_dnrm2(NUMBER_OF_AXIS_, c, 1); // calcualte ||u||
  cblas_dscal(NUMBER_OF_AXIS_, (1.0 / R), c, 1);  // normalize u
  return c;
}

void
Fragment::rebase_Coordinates()
{
  unsigned int atom_index, coordinate_index;
  for (atom_index = coordinate_index = 0;
       ((atom_index < current_Fragment.size()) &&
        (coordinate_index < NumberOfCoordinates));
       ++atom_index, (coordinate_index += 3))
  {
    current_Fragment.at(atom_index)
        ->setCoordinates((current_Cartesian_Coordinates + coordinate_index),
                         StructureOptions::Initial);
  }
}

RedundantType
Fragment::Potential_2_Redundant_Type(PotentialType pt)
{
  switch (pt)
  {
    case PotentialType::Undefined:
      return RedundantType::Undefined;
    case PotentialType::Stretch:
      return RedundantType::Distance;
    case PotentialType::Bend:
      return RedundantType::Angle;
    case PotentialType::Dihedral:
      return RedundantType::Torsion;
    case PotentialType::StretchBend:
      return RedundantType::Undefined;
    case PotentialType::vanderWaals:
      return RedundantType::Undefined;
    default:
      return RedundantType::Undefined;
  }
}

int
Fragment::NumberOfParticipants(RedundantType rt)
{
  switch (rt)
  {
    case RedundantType::Undefined:
      return 0;
    case RedundantType::Distance:
      return 2;
    case RedundantType::Angle:
      return 3;
    case RedundantType::Torsion:
      return 4;
    case RedundantType::RDC:
      return 2;
    default:
      return 0;
  }
}

void
Fragment::reserve_Memory()
{
  initial_Internal_Coordinates =
      (double *) malloc(NumberOfRedundants * sizeof(double));
  current_Internal_Coordinates =
      (double *) malloc(NumberOfRedundants * sizeof(double));
  optimal_Internal_Coordinates =
      (double *) malloc(NumberOfRedundants * sizeof(double));
  difference_Internal_Coordinates =
      (double *) malloc(NumberOfRedundants * sizeof(double));
  current_Cartesian_Coordinates =
      (double *) malloc(NumberOfCoordinates * sizeof(double));
  difference_Cartesian_Coordinates =
      (double *) malloc(NumberOfCoordinates * sizeof(double));
  Wilson_Matrix = (double *) malloc(NumberOfRedundants * NumberOfCoordinates *
                                    sizeof(double));
  inverse_Wilson_Matrix = (double *) malloc(
      NumberOfRedundants * NumberOfCoordinates * sizeof(double));
  list_of_Atoms = (Atom **) malloc(NumberOfAtoms * sizeof(Atom *));
  InternalsList = (InternalCoordinate *) malloc(NumberOfRedundants *
                                                sizeof(InternalCoordinate));
}

void
Fragment::perform_optimization_step()
{
  for (unsigned int i = 0; i < NumberOfRedundants; ++i)
  {
    current_Internal_Coordinates[i] = Fragment::calculate_value(
        (InternalsList + i), optimization_Atoms, current_Cartesian_Coordinates,
        StructureOptions::Initial);
  }
  cblas_dcopy(NumberOfRedundants, optimal_Internal_Coordinates, 1,
              difference_Internal_Coordinates, 1); // copy m in u
  cblas_daxpy(NumberOfRedundants, (-1.0), current_Internal_Coordinates, 1,
              difference_Internal_Coordinates, 1); // calculate u = m - n

  for (unsigned int i = 0; i < NumberOfCoordinates; ++i)
    difference_Cartesian_Coordinates[i] = .0;

  this->load_Wilson_B();

  double nullspacer[NumberOfRedundants * NumberOfRedundants];

  for (unsigned int i = 0; i < NumberOfRedundants * NumberOfRedundants; ++i)
    nullspacer[i] = .0;

  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, NumberOfCoordinates,
              NumberOfRedundants, NumberOfCoordinates, 1.0, Wilson_Matrix,
              NumberOfCoordinates, inverse_Wilson_Matrix, NumberOfRedundants,
              0.0, nullspacer, NumberOfRedundants);
  bool run = true;

  if (run)
  {
    cblas_dgemv(CblasRowMajor, CblasNoTrans, NumberOfRedundants,
                NumberOfRedundants, 1.0, nullspacer, NumberOfRedundants,
                difference_Internal_Coordinates, 1, 0.0,
                current_Internal_Coordinates, 1);
  }

  for (unsigned int i = 0; i < NumberOfRedundants && run; ++i)
  {
    //      std::cout << i << ":   " << current_Internal_Coordinates[i] << "  "
    //      << difference_Internal_Coordinates[i] << std::endl;
    difference_Internal_Coordinates[i] = current_Internal_Coordinates[i];
  }

  cblas_dgemv(CblasRowMajor, CblasNoTrans, NumberOfCoordinates,
              NumberOfRedundants, 1.0, inverse_Wilson_Matrix,
              NumberOfRedundants, difference_Internal_Coordinates, 1, 0.0,
              difference_Cartesian_Coordinates, 1);
  cblas_daxpy(NumberOfCoordinates, (1.0), difference_Cartesian_Coordinates, 1,
              current_Cartesian_Coordinates, 1); // calculate u = m - n
}

Fragment *
Fragment::check_for_fusing()
{
  Fragment *CurrentFragment = this;
  Fragment *output_fragment = this;
  Fragment *lastFragment    = this;
  Fragment *linker;
  while (lastFragment->nextFragment)
    lastFragment = lastFragment->nextFragment;

  while (CurrentFragment)
  {
    if (CurrentFragment->NumberOfAtoms <= 8)
    {
      for (unsigned int i = 0; i < CurrentFragment->current_Fragment.size();
           ++i)
      {
        for (unsigned int j = 0;
             j < CurrentFragment->current_Fragment.at(i)->getNumberOfBonds();
             ++j)
        {
          linker = Fragment::is_in_Fragment(
              lastFragment,
              CurrentFragment->current_Fragment.at(i)->getBondpartner(j));
          // If this is the first fragment it cant be added to another one.
          if (linker && linker != CurrentFragment)
          {
            if (CurrentFragment != this)
            {
              Fragment::add_small_subfragment(
                  linker, CurrentFragment->current_Fragment);
              j = CurrentFragment->current_Fragment.at(i)->getNumberOfBonds();
              i = CurrentFragment->current_Fragment.size();
              CurrentFragment->current_Fragment.clear();
              CurrentFragment->to_check_Atoms.clear();
              Fragment::relink_Fragments(CurrentFragment->prevFragment,
                                         CurrentFragment->nextFragment);
              linker                        = CurrentFragment->prevFragment;
              CurrentFragment->prevFragment = CurrentFragment->nextFragment =
                  NULL;
              delete CurrentFragment;
              CurrentFragment = linker;
              break;
            }
            else
            {
              Fragment::add_small_subfragment(CurrentFragment,
                                              linker->current_Fragment);
              j = CurrentFragment->current_Fragment.at(i)->getNumberOfBonds();
              i = CurrentFragment->current_Fragment.size();
              linker->current_Fragment.clear();
              linker->to_check_Atoms.clear();
              Fragment::relink_Fragments(linker->prevFragment,
                                         linker->nextFragment);
              linker->nextFragment = linker->prevFragment = NULL;
              delete linker;
              break;
            }
          }
        }
      }
    }
    CurrentFragment = CurrentFragment->nextFragment;
  }
  CurrentFragment = linker = lastFragment = NULL;
  return output_fragment;
}

void
Fragment::load_Potentials(Potential *CurrP,
                          SphericalHarmonics *CurrY,
                          //                                Flags &flags,
                          StructureOptions opt)
{
  unsigned int atom_index, list_index;

  //   if ( fragment_type == FragmentFlag::RDC ) use_torsion_angles =
  //   flags.torsions_4_redundants; else
  use_torsion_angles = true;

  SphericalHarmonics *HeadY = CurrY;
  Potential *HeadP          = CurrP;

  to_check_Atoms.clear();

  for (; CurrP; CurrP = CurrP->getNext())
    check_Potential(CurrP, false, true);

  for (; CurrY; CurrY = CurrY->getNext())
  {
    check_Harmonic(CurrY);
  }
  NumberOfAtoms = to_check_Atoms.size() + current_Fragment.size();
  optimization_Atoms.reserve(NumberOfAtoms);
  optimization_Atoms.insert(optimization_Atoms.end(), current_Fragment.begin(),
                            current_Fragment.end());
  optimization_Atoms.insert(optimization_Atoms.end(), to_check_Atoms.begin(),
                            to_check_Atoms.end());
  NumberOfBonds = NumberOfAngles = NumberOfTorsions = 0;

  for (CurrP = HeadP; CurrP; CurrP = CurrP->getNext())
    check_Potential(CurrP, true, true);

  for (CurrY = HeadY; CurrY; CurrY = CurrY->getNext())
  {
    if (check_Harmonic(CurrY))
      ++NumberOfRDCs;
  }

  NumberOfRedundants = NumberOfBonds + NumberOfAngles + NumberOfTorsions +
                       NumberOfRDCs + NumberOfDistances;
  NumberOfCoordinates = 3 * NumberOfAtoms;

  reserve_Memory();

  for (atom_index = 0; atom_index < NumberOfAtoms; ++atom_index)
  {
    list_of_Atoms[atom_index] = optimization_Atoms.at(atom_index);
    current_Cartesian_Coordinates[3 * atom_index + STANDARD_AXIS_X_] =
        list_of_Atoms[atom_index]->Coordinates2Eigen(opt)(STANDARD_AXIS_X_);
    current_Cartesian_Coordinates[3 * atom_index + STANDARD_AXIS_Y_] =
        list_of_Atoms[atom_index]->Coordinates2Eigen(opt)(STANDARD_AXIS_Y_);
    current_Cartesian_Coordinates[3 * atom_index + STANDARD_AXIS_Z_] =
        list_of_Atoms[atom_index]->Coordinates2Eigen(opt)(STANDARD_AXIS_Z_);
    difference_Cartesian_Coordinates[3 * atom_index + STANDARD_AXIS_X_] = .0;
    difference_Cartesian_Coordinates[3 * atom_index + STANDARD_AXIS_Y_] = .0;
    difference_Cartesian_Coordinates[3 * atom_index + STANDARD_AXIS_Z_] = .0;
  }

  for (CurrP = HeadP, list_index = 0; CurrP; CurrP = CurrP->getNext())
  {
    if (check_Potential(CurrP, true, false))
    {
      if (CurrP->getType() == PotentialType::Bend)
        optimal_Internal_Coordinates[list_index] =
            deg2rad(CurrP->get_equilibrium());
      else
        optimal_Internal_Coordinates[list_index] = CurrP->get_equilibrium();
      list_index = save_Potential(CurrP, NULL, list_index);
    }
  }
  for (CurrY = HeadY; CurrY; CurrY = CurrY->getNext())
  {
    if (check_Harmonic(CurrY))
    {
      optimal_Internal_Coordinates[list_index] = .0;
      list_index = save_Potential(NULL, CurrY, list_index);
    }
  }

  HeadP = CurrP = NULL;
  HeadY = CurrY = NULL;
}

unsigned int
Fragment::save_Potential(Potential *CurrP,
                         SphericalHarmonics *CurrY,
                         unsigned int &index)
{
  // First check if the current potential is an structure
  // or rdc based potential.
  if (CurrP)
    InternalsList[index].type =
        Fragment::Potential_2_Redundant_Type(CurrP->getType());
  else if (CurrY)
    InternalsList[index].type = RedundantType::RDC;
  else
    InternalsList[index].type = RedundantType::Undefined;

  // Reserve the respectiv memory space for atom pointers.
  InternalsList[index].participants = (Atom **) malloc(
      Fragment::NumberOfParticipants(InternalsList[index].type) *
      sizeof(Atom *));

  // Save the participating atoms
  if (CurrP)
  {
    for (int j = 0;
         j < Fragment::NumberOfParticipants(InternalsList[index].type); ++j)
      InternalsList[index].participants[j] = CurrP->getAtom(j);
  }
  else if (CurrY)
  {
    InternalsList[index].participants[0] = CurrY->getAtom1();
    InternalsList[index].participants[1] = CurrY->getAtom2();
  }

  // Save the index and the respective object (Y/P)
  InternalsList[index].index     = index;
  InternalsList[index].potential = CurrP;
  InternalsList[index].harmonic  = CurrY;

  // Save the current value as initial and current value
  InternalsList[index].value = initial_Internal_Coordinates[index] =
      current_Internal_Coordinates[index] =
          Fragment::calculate_value((InternalsList + index), optimization_Atoms,
                                    current_Cartesian_Coordinates,
                                    StructureOptions::Initial);

  // Calculate the current difference
  if (InternalsList[index].type == RedundantType::Torsion)
  {
    optimal_Internal_Coordinates[index] = current_Internal_Coordinates[index];
  }

  // Increment the current index by one.
  // By this implementation it is possible
  // to increment by using the return value
  // or to retain it by just letting the
  // return value run into nirvana.
  return (index + 1);
}

void
Fragment::combine_Fragments(Fragment *bonded_Fragment)
{
  Atom *CurrAtom;
  std::vector<Atom *> overlaped_atoms;
  unsigned int overlap = 0;
  unsigned int atom_index, axis_index, current_coordinate;
  double sub_fragment_mass;

  // Check overlapping atoms
  for (atom_index = 0; atom_index < NumberOfAtoms; ++atom_index)
  {
    CurrAtom = optimization_Atoms.at(atom_index);
    if (Fragment::is_in_Fragment(this, bonded_Fragment, CurrAtom))
    {
      overlaped_atoms.push_back(CurrAtom);
      ++overlap;
    }
  }

  // allocate memory for center of mass
  double *COM_1 = (double *) malloc(double_vector_size);
  double *COM_2 = (double *) malloc(double_vector_size);

  // make sure memory is clean
  sub_fragment_mass = .0;
  for (axis_index = 0; axis_index < NUMBER_OF_AXIS_; ++axis_index)
    COM_1[axis_index] = COM_2[axis_index] = .0;


  for (atom_index = 0; atom_index < overlap; ++atom_index)
  {
    CurrAtom = overlaped_atoms.at(atom_index);
    sub_fragment_mass += CurrAtom->getMass();

    // Add mass weighted coordinates to center of mass for MainFragment
    current_coordinate = 3 * get_vector_index(optimization_Atoms, CurrAtom);
    cblas_daxpy(NUMBER_OF_AXIS_, CurrAtom->getMass(),
                (current_Cartesian_Coordinates + current_coordinate), 1, COM_1,
                1);

    // Add mass weighted coordinates to center of mass for BondedFragment
    current_coordinate =
        3 * get_vector_index(bonded_Fragment->optimization_Atoms, CurrAtom);
    cblas_daxpy(
        NUMBER_OF_AXIS_, CurrAtom->getMass(),
        (bonded_Fragment->current_Cartesian_Coordinates + current_coordinate),
        1, COM_2, 1);
  }

  // Divide by fragment mass
  cblas_dscal(NUMBER_OF_AXIS_, (1.0 / sub_fragment_mass), COM_1, 1);
  cblas_dscal(NUMBER_OF_AXIS_, (1.0 / sub_fragment_mass), COM_2, 1);

  // Calculate difference
  cblas_daxpy(NUMBER_OF_AXIS_, (-1.0), COM_2, 1, COM_1, 1);

  // And shift BONDEDFragment by the difference
  for (atom_index = 0; atom_index < bonded_Fragment->NumberOfAtoms;
       ++atom_index)
  {
    cblas_daxpy(
        NUMBER_OF_AXIS_, 1.0, COM_1, 1,
        (bonded_Fragment->current_Cartesian_Coordinates + 3 * atom_index), 1);
  }

  free(COM_1);
  free(COM_2);
  bonded_Fragment->rebase_Coordinates();
  this->rebase_Coordinates();
}

void
Fragment::load_Wilson_B()
{
  // WILSON_B(NUMBER_OF_REDUNDANTS x NUMBER_OF_CARTESIAN_COORDINATES)
  unsigned int Wilson_B_Index, axis, Atom_Index, Redundant_Index;
  Atom *Current_Atom;
  InternalCoordinate *Current_Internal_Coordinate;
  Wilson_B_Index = 0;
  for (Redundant_Index = 0; Redundant_Index < NumberOfRedundants;
       ++Redundant_Index)
  {
    Current_Internal_Coordinate = &(InternalsList[Redundant_Index]);
    for (Atom_Index = 0; Atom_Index < NumberOfAtoms; ++Atom_Index)
    {
      Current_Atom = optimization_Atoms.at(Atom_Index);
      for (axis = 0; axis < NUMBER_OF_AXIS_; ++axis, ++Wilson_B_Index)
      {
        Wilson_Matrix[Wilson_B_Index] = Fragment::calculate_derivative(
            Current_Internal_Coordinate, optimization_Atoms, Current_Atom, axis,
            current_Cartesian_Coordinates);
      }
    }
  }
  if (fragment_type == FragmentFlag::RDC)
    svdMoorePenroseInverse(Wilson_Matrix, (int) NumberOfRedundants,
                           (int) NumberOfCoordinates, inverse_Wilson_Matrix);
  else
    svdMoorePenroseInverse(Wilson_Matrix, (int) NumberOfRedundants,
                           (int) NumberOfCoordinates, inverse_Wilson_Matrix);
}

// void Fragment::printCoordinates (BasicInformation &baseInformation)
//{
//   int j = 0;
//
//   fprintf(baseInformation.trjFile, "%d\n%s\n", NumberOfAtoms,
//   optimization_Atoms.at(0)->getIdentifier().c_str()); for ( unsigned int i =
//   0; i < NumberOfAtoms; ++i )
//   {
//      fprintf(baseInformation.trjFile, "%3.s    %3.5f   %3.5f   %3.5f\n",
//                                        optimization_Atoms.at(i)->getElement().c_str(),
//                                        current_Cartesian_Coordinates[j++],
//                                        current_Cartesian_Coordinates[j++],
//                                        current_Cartesian_Coordinates[j++]);
//   }
//}

// void Fragment::printJmol()
//{
//   std::cout << "NumberOfBonds: " << NumberOfBonds << std::endl
//             << "NumberOfAngles: " << NumberOfAngles << std::endl
//             << "NumberOfTorsions: " << NumberOfTorsions << std::endl
//             << "NumberOfRDCs: " << NumberOfRDCs << std::endl
//             << "NumberOfRedundants: " << NumberOfRedundants << std::endl;
//   std::cout << "NumberOfAtoms: " << NumberOfAtoms << std::endl
//             << "NumberOfCoordinates: " << NumberOfCoordinates << std::endl;
//
//   unsigned int start = NumberOfBonds+NumberOfAngles;
//   unsigned int end = NumberOfBonds+NumberOfAngles+NumberOfTorsions;
//   for ( unsigned int list_index = start; list_index < end; ++list_index )
//   {
//      std::cout << InternalsList[list_index].participants[0]->getIdentifier();
//      std::cout << "-" <<
//      InternalsList[list_index].participants[1]->getIdentifier(); std::cout <<
//      "-" << InternalsList[list_index].participants[2]->getIdentifier();
//      std::cout << "-" <<
//      InternalsList[list_index].participants[3]->getIdentifier(); std::cout <<
//      ":  "; std::cout << rad2deg(optimal_Internal_Coordinates[list_index]) <<
//      "  " << rad2deg(current_Internal_Coordinates[list_index]) << "  " <<
//      difference_Internal_Coordinates[list_index] << std::endl;
//   }
//
//   if ( current_Fragment.size() ) std::cout << "select " <<
//   current_Fragment.at(0)->getIdentifier(); for ( unsigned int i = 1; i <
//   current_Fragment.size(); ++i )
//   {
//      std::cout << ", " << current_Fragment.at(i)->getIdentifier();
//   }
//   if ( to_check_Atoms.size() ) std::cout << std::endl << "select add " <<
//   to_check_Atoms.at(0)->getIdentifier(); for ( unsigned int i = 1; i <
//   to_check_Atoms.size(); ++i )
//   {
//      std::cout << ", " << to_check_Atoms.at(i)->getIdentifier();
//   }
//}

/*********************************************************************
 *                                                                   *
 * svdMoorePenrose:                                                  *
 *                                                                   *
 * Routine to calculate the generalized (pseudo-)inverse of an       *
 * m-by-n matrix A                                                   *
 *                                                                   *
 * The calculation proceeds by singular value decomposition of A and *
 * uses the DGESVD routine from LAPACK. The output matrices of that  *
 * routine are then transposed and inverted as required and finally  *
 * matrix multiplied to calculate the inverse.                       *
 *                                                                   *
 * Input is the m-by-n matrix a. Upon return the n-by-m (!) matrix b *
 * contains the result.                                              *
 *                                                                   *
 *********************************************************************/

void
svdMoorePenroseInverse(double *a, int m, int n, double *ainv, double thresh)
{
  // define leading dimensions of the matrices, so we can use the same syntax
  // as lapack (see e.g. the example on software.intel.com)
  int info = 0, ldu = m, ldvt = n, lda = n;
  int number_of_elements_a = n * m;
  double u[ldu * m], w[n], vt[ldvt * n], superb[(m > n ? n : m) - 1];
  double aa[m * n];

  // as the LAPACK code scrambles a upon calculation, we transfer the initial
  // a-values to aa

  cblas_dcopy(number_of_elements_a, a, 1, &(aa[0]), 1); // copy m in u


  // LAPACK version of Double GEnereal matrix Singular Value Decomposition
  // we have C-style arrays, so set LAPACK_ROW_MAJOR for proper matrix memory
  // management, upon exit store the left singular vectors in u, the singular
  // values in w, the transpose of the right singular vectors in v and
  // ignore superb.
  info = LAPACKE_dgesvd(CblasRowMajor, 'A', 'A', m, n, aa, lda, w, u, ldu, vt,
                        ldvt, superb);

  if (info)
    printf("RDC warning: SVD -- the computed matrix may be singular (%d)!\n",
           info);

  // If singular values w[i] are really small the back-calculation
  // (using 1/w[i]) gives wrong results. So if a w[i] is less than some
  // threshold, the value in the resulting array is set to zero.
  double winv[n * m], tmp[n * m];
  //   thresh = 0.5;//max(LAPACKE_dlamch('P')*w[0],LAPACKE_dlamch('S'));

  // setup an n-by-m matrix of 'diagonal' values 1/w[i] for back-calculation
  //   double t1, t2;
  //   int ti = 0;
  //   t2 = .0;
  //   t1 = w[0];

  for (int i = 0; i < n; i++)
  {
    // zero values
    for (int j = 0; j < m; j++)
      winv[i * m + j] = 0.0;

    // properly set value of the diagonal elements
    if (w[i] > thresh)
      winv[i * m + i] = 1.0 / w[i];
    //      else if ( t2 == .0 ){ t2 = w[i-1]; ti = i-1; }
  }
  //   std::cout << "w(0)/w(" << ti << ") = " << t1 << "/" << t2  << " = " <<
  //   t1/t2 << std::endl;

  // matrix multiply to calculate pinv(A) = v * pinv(w) * u_t
  // (if a is an m-by-n matrix, ainv has dimensions n-by-m!)
  //
  // tmp = winv * u_t
  // as dgesvd returns u, we transpose that here
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n, m, n, 1.0, winv, m, u,
              m, 0.0, tmp, m);

  // ainv = v * tmp
  // as dgesvd returns v_t, we have to transpose that
  cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, n, m, n, 1.0, vt, n, tmp,
              m, 0.0, ainv, m);
}
