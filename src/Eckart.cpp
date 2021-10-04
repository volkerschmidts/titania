
#include <Atom.hpp>
#include <Declarations.hpp>
#include <Eckart.hpp>
#include <Structure.hpp>
#include <iostream>
#include <phobos.h>

/*
 * Eckart frames
 *
 * These implementations follow P.R. Bunker, P. Jensen, Molecular Symmetry and
 * Spectroscopy, NRC Research Press, 2006, chaper 10.2.
 */

void
build_Eckart_AlphaTau(Structure *RefStruc,
                      Structure *CurrStruc,
                      Eigen::Matrix3d &AlphaTau,
                      StructureOptions options)
{
  AlphaTau(0, 0) =
      Eckart_AlphaTau(Axis::x, Axis::Xi, RefStruc, CurrStruc, options);
  AlphaTau(0, 1) =
      Eckart_AlphaTau(Axis::x, Axis::Eta, RefStruc, CurrStruc, options);
  AlphaTau(0, 2) =
      Eckart_AlphaTau(Axis::x, Axis::Zeta, RefStruc, CurrStruc, options);

  AlphaTau(1, 0) =
      Eckart_AlphaTau(Axis::y, Axis::Xi, RefStruc, CurrStruc, options);
  AlphaTau(1, 1) =
      Eckart_AlphaTau(Axis::y, Axis::Eta, RefStruc, CurrStruc, options);
  AlphaTau(1, 2) =
      Eckart_AlphaTau(Axis::y, Axis::Zeta, RefStruc, CurrStruc, options);

  AlphaTau(2, 0) =
      Eckart_AlphaTau(Axis::z, Axis::Xi, RefStruc, CurrStruc, options);
  AlphaTau(2, 1) =
      Eckart_AlphaTau(Axis::z, Axis::Eta, RefStruc, CurrStruc, options);
  AlphaTau(2, 2) =
      Eckart_AlphaTau(Axis::z, Axis::Zeta, RefStruc, CurrStruc, options);
}

void
build_Eckart_Rotation(double *bestFit, Eigen::Matrix3d &Rotation)
{
  build_Eckart_Rotation(bestFit[ECKART_THETA], bestFit[ECKART_PHI],
                        bestFit[ECKART_CHI], Rotation);
}

void
build_Eckart_Rotation(double theta,
                      double phi,
                      double chi,
                      Eigen::Matrix3d &Rotation)
{
  Rotation(0, 0) = Eckart_Lambda(Axis::x, Axis::Xi, theta, phi, chi);
  Rotation(0, 1) = Eckart_Lambda(Axis::x, Axis::Eta, theta, phi, chi);
  Rotation(0, 2) = Eckart_Lambda(Axis::x, Axis::Zeta, theta, phi, chi);

  Rotation(1, 0) = Eckart_Lambda(Axis::y, Axis::Xi, theta, phi, chi);
  Rotation(1, 1) = Eckart_Lambda(Axis::y, Axis::Eta, theta, phi, chi);
  Rotation(1, 2) = Eckart_Lambda(Axis::y, Axis::Zeta, theta, phi, chi);

  Rotation(2, 0) = Eckart_Lambda(Axis::z, Axis::Xi, theta, phi, chi);
  Rotation(2, 1) = Eckart_Lambda(Axis::z, Axis::Eta, theta, phi, chi);
  Rotation(2, 2) = Eckart_Lambda(Axis::z, Axis::Zeta, theta, phi, chi);
}

void
get_Eckart_Angles(Structure *RefStruc,
                  Structure *CurrStruc,
                  BasicInformation &baseInformation,
                  StructureOptions options)
{
  int i;

  Eckart_Opt_Parameters data;
  build_Eckart_AlphaTau(RefStruc, CurrStruc, data.AlphaTau, options);
  data.baseInformation = &baseInformation;

  for (i = 0; i < baseInformation.phobos_worksize_Eckart; ++i)
    baseInformation.phobos_workspace_Eckart[i] = .0;

  double *bestFit   = baseInformation.phobos_workspace_Eckart;
  double *measure   = bestFit + NUMBER_OF_ECKART_ANGLES;
  double *workspace = measure + NUMBER_OF_ECKART_FUNTIONS;
  double *covar =
      bestFit + baseInformation.phobos_worksize_Eckart - ECKART_COVAR_SIZE;

  phobos(Eckart_Angles, NULL, bestFit, measure, NUMBER_OF_ECKART_ANGLES,
         NUMBER_OF_ECKART_FUNTIONS, baseInformation.limits.max_lm_iterations,
         baseInformation.limits.phobos_opts_eckart, NULL, workspace, covar,
         (void *) &data);
  covar = workspace = measure = bestFit = NULL;
  data.baseInformation                  = NULL;
}

double
Eckart_Lambda(Axis axis_ref,
              Axis axis_curr,
              double theta,
              double phi,
              double chi)
{
  switch (axis_ref)
  {
    case (Axis::x): {
      switch (axis_curr)
      {
        case (Axis::Xi):
          return (cos(theta) * cos(phi) * cos(chi) - sin(phi) * sin(chi));
        case (Axis::Eta):
          return (cos(theta) * sin(phi) * cos(chi) + cos(phi) * sin(chi));
        case (Axis::Zeta):
          return (-sin(theta) * cos(chi));
        default:
          return .0;
      }
      break;
    }
    case (Axis::y): {
      switch (axis_curr)
      {
        case (Axis::Xi):
          return (-cos(theta) * cos(phi) * sin(chi) - sin(phi) * cos(chi));
        case (Axis::Eta):
          return (-cos(theta) * sin(phi) * sin(chi) + cos(phi) * cos(chi));
        case (Axis::Zeta):
          return (sin(theta) * sin(chi));
        default:
          return .0;
      }
      break;
    }
    case (Axis::z): {
      switch (axis_curr)
      {
        case (Axis::Xi):
          return (sin(theta) * cos(phi));
        case (Axis::Eta):
          return (sin(theta) * sin(phi));
        case (Axis::Zeta):
          return (cos(theta));
        default:
          return .0;
      }
      break;
    }
    default:
      return .0;
  }
  return .0;
}

double
Eckart_AlphaTau(Axis axis_ref,
                Axis axis_curr,
                Structure *RefStruc,
                Structure *CurrStruc,
                StructureOptions options)
{
  double m, AlphaTau = .0;
  int index_ref = 0, index_curr = 0;
  Eigen::Vector3d RefCoord, CurrCoord;
  Atom *CurrAtom, *RefAtom;
  switch (axis_ref)
  {
    case (Axis::x):
      index_ref = 0;
      break;
    case (Axis::y):
      index_ref = 1;
      break;
    case (Axis::z):
      index_ref = 2;
      break;
    default:
      return .0;
  }

  switch (axis_curr)
  {
    case (Axis::Xi):
      index_curr = 0;
      break;
    case (Axis::Eta):
      index_curr = 1;
      break;
    case (Axis::Zeta):
      index_curr = 2;
      break;
    default:
      return .0;
  }
  for (CurrAtom = CurrStruc->getHeadAtom(), RefAtom = RefStruc->getHeadAtom();
       CurrAtom && RefAtom;
       CurrAtom = CurrAtom->getNext(), RefAtom = RefAtom->getNext())
  {
    m         = CurrAtom->getMass();
    CurrCoord = CurrAtom->Coordinates2Eigen(options);
    RefCoord  = RefAtom->Coordinates2Eigen(options);
    AlphaTau += (m * CurrCoord(index_curr) * RefCoord(index_ref));
  }
  return AlphaTau;
}

void
Eckart_Angles(double *p, double *x, int m, int n, void *data)
{
  Eckart_Opt_Parameters *dptr;
  dptr = (Eckart_Opt_Parameters *) data;
  if (m != n)
  {
    dptr->baseInformation->errorKey +=
        ("ERROR:\tSomething strange happened in " +
         std::string(__PRETTY_FUNCTION__) +
         "...\n\tTITANIA will ignore this error and stop running " +
         "Eckart transformations...\n");
    dptr->baseInformation->state = ERROR_ECKART_TRANSFORMATION;
  }
  Eigen::Matrix3d AT = dptr->AlphaTau;

  x[STANDARD_AXIS_X_] =
      fabs(AT(STANDARD_AXIS_X_, STANDARD_AXIS_X_) *
               Eckart_Lambda(Axis::y, Axis::Xi, p[ECKART_THETA], p[ECKART_PHI],
                             p[ECKART_CHI]) +
           AT(STANDARD_AXIS_X_, STANDARD_AXIS_Y_) *
               Eckart_Lambda(Axis::y, Axis::Eta, p[ECKART_THETA], p[ECKART_PHI],
                             p[ECKART_CHI]) +
           AT(STANDARD_AXIS_X_, STANDARD_AXIS_Z_) *
               Eckart_Lambda(Axis::y, Axis::Zeta, p[ECKART_THETA],
                             p[ECKART_PHI], p[ECKART_CHI]) -
           AT(STANDARD_AXIS_Y_, STANDARD_AXIS_X_) *
               Eckart_Lambda(Axis::x, Axis::Xi, p[ECKART_THETA], p[ECKART_PHI],
                             p[ECKART_CHI]) -
           AT(STANDARD_AXIS_Y_, STANDARD_AXIS_Y_) *
               Eckart_Lambda(Axis::x, Axis::Eta, p[ECKART_THETA], p[ECKART_PHI],
                             p[ECKART_CHI]) -
           AT(STANDARD_AXIS_Y_, STANDARD_AXIS_Z_) *
               Eckart_Lambda(Axis::x, Axis::Zeta, p[ECKART_THETA],
                             p[ECKART_PHI], p[ECKART_CHI]));

  x[STANDARD_AXIS_Y_] =
      fabs(AT(STANDARD_AXIS_Y_, STANDARD_AXIS_X_) *
               Eckart_Lambda(Axis::z, Axis::Xi, p[ECKART_THETA], p[ECKART_PHI],
                             p[ECKART_CHI]) +
           AT(STANDARD_AXIS_Y_, STANDARD_AXIS_Y_) *
               Eckart_Lambda(Axis::z, Axis::Eta, p[ECKART_THETA], p[ECKART_PHI],
                             p[ECKART_CHI]) +
           AT(STANDARD_AXIS_Y_, STANDARD_AXIS_Z_) *
               Eckart_Lambda(Axis::z, Axis::Zeta, p[ECKART_THETA],
                             p[ECKART_PHI], p[ECKART_CHI]) -
           AT(STANDARD_AXIS_Z_, STANDARD_AXIS_X_) *
               Eckart_Lambda(Axis::y, Axis::Xi, p[ECKART_THETA], p[ECKART_PHI],
                             p[ECKART_CHI]) -
           AT(STANDARD_AXIS_Z_, STANDARD_AXIS_Y_) *
               Eckart_Lambda(Axis::y, Axis::Eta, p[ECKART_THETA], p[ECKART_PHI],
                             p[ECKART_CHI]) -
           AT(STANDARD_AXIS_Z_, STANDARD_AXIS_Z_) *
               Eckart_Lambda(Axis::y, Axis::Zeta, p[ECKART_THETA],
                             p[ECKART_PHI], p[ECKART_CHI]));

  x[STANDARD_AXIS_Z_] =
      fabs(AT(STANDARD_AXIS_Z_, STANDARD_AXIS_X_) *
               Eckart_Lambda(Axis::x, Axis::Xi, p[ECKART_THETA], p[ECKART_PHI],
                             p[ECKART_CHI]) +
           AT(STANDARD_AXIS_Z_, STANDARD_AXIS_Y_) *
               Eckart_Lambda(Axis::x, Axis::Eta, p[ECKART_THETA], p[ECKART_PHI],
                             p[ECKART_CHI]) +
           AT(STANDARD_AXIS_Z_, STANDARD_AXIS_Z_) *
               Eckart_Lambda(Axis::x, Axis::Zeta, p[ECKART_THETA],
                             p[ECKART_PHI], p[ECKART_CHI]) -
           AT(STANDARD_AXIS_X_, STANDARD_AXIS_X_) *
               Eckart_Lambda(Axis::z, Axis::Xi, p[ECKART_THETA], p[ECKART_PHI],
                             p[ECKART_CHI]) -
           AT(STANDARD_AXIS_X_, STANDARD_AXIS_Y_) *
               Eckart_Lambda(Axis::z, Axis::Eta, p[ECKART_THETA], p[ECKART_PHI],
                             p[ECKART_CHI]) -
           AT(STANDARD_AXIS_X_, STANDARD_AXIS_Z_) *
               Eckart_Lambda(Axis::z, Axis::Zeta, p[ECKART_THETA],
                             p[ECKART_PHI], p[ECKART_CHI]));
  dptr = NULL;
}

double *
checkEckart(Structure *RefStruc,
            Structure *CheckStruc,
            enum StructureOptions options)
{
  double *eckart_dev = (double *) malloc(3 * sizeof(double));
  double m;
  Atom *CurrAtom, *RefAtom;
  Eigen::Vector3d CurrCoord, RefCoord;
  eckart_dev[STANDARD_AXIS_X_]     = eckart_dev[STANDARD_AXIS_Y_] =
      eckart_dev[STANDARD_AXIS_Z_] = .0;
  for (CurrAtom = CheckStruc->getHeadAtom(), RefAtom = RefStruc->getHeadAtom();
       CurrAtom && RefAtom;
       CurrAtom = CurrAtom->getNext(), RefAtom = RefAtom->getNext())
  {
    CurrCoord = CurrAtom->Coordinates2Eigen(options);
    RefCoord  = RefAtom->Coordinates2Eigen(options);
    m         = CurrAtom->getMass();
    eckart_dev[STANDARD_AXIS_X_] +=
        (m * (RefCoord(STANDARD_AXIS_X_) * CurrCoord(STANDARD_AXIS_Y_) -
              RefCoord(STANDARD_AXIS_Y_) * CurrCoord(STANDARD_AXIS_X_)));
    eckart_dev[STANDARD_AXIS_Y_] +=
        (m * (RefCoord(STANDARD_AXIS_Y_) * CurrCoord(STANDARD_AXIS_Z_) -
              RefCoord(STANDARD_AXIS_Z_) * CurrCoord(STANDARD_AXIS_Y_)));
    eckart_dev[STANDARD_AXIS_Z_] +=
        (m * (RefCoord(STANDARD_AXIS_Z_) * CurrCoord(STANDARD_AXIS_X_) -
              RefCoord(STANDARD_AXIS_X_) * CurrCoord(STANDARD_AXIS_Z_)));
  }
  return eckart_dev;
}

void
Rotate_2_Eckart_Frame(Structure *ReferenceStructure,
                      Structure *Structure2Rotate,
                      BasicInformation &baseInformation,
                      StructureOptions options)
{
  Structure2Rotate->Shift2CenterOfMass(options);
  get_Eckart_Angles(ReferenceStructure, Structure2Rotate, baseInformation,
                    options);

  Eigen::Matrix3d Rotation;
  build_Eckart_Rotation(baseInformation.phobos_workspace_Eckart, Rotation);
  for (Atom *CurrAtom = Structure2Rotate->getHeadAtom(); CurrAtom;
       CurrAtom       = CurrAtom->getNext())
  {
    CurrAtom->rotateCoordinates(Rotation, options);
  }
}
