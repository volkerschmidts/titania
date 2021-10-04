
#include <Atom.hpp>
#include <Declarations.hpp>
#include <IndependencyAnalysis.hpp>
#include <LinAlg.hpp>
#include <Molecule.hpp>
#include <RDCdata.hpp>
#include <RDCset.hpp>
#include <SphericalHarmonics.hpp>
#include <SphericalHarmonicsMatrix.hpp>
#include <Structure.hpp>
#include <eigen3/Eigen/Eigenvalues>
#include <iostream>
#include <vector>

void
sort_evals(Eigen::VectorXd &Evals, Eigen::MatrixXd &Evecs)
{
  int l = Evals.size();
  std::vector<std::tuple<double, Eigen::VectorXd>> eigen_vectors_and_values;
  for (int i = 0; i < l; ++i)
  {
    std::tuple<double, Eigen::VectorXd> vec_and_val(Evals(i), Evecs.col(i));
    eigen_vectors_and_values.push_back(vec_and_val);
  }
  std::sort(eigen_vectors_and_values.begin(), eigen_vectors_and_values.end(),
            [&](const std::tuple<double, Eigen::VectorXd> &a,
                const std::tuple<double, Eigen::VectorXd> &b) -> bool {
              return std::get<0>(a) > std::get<0>(b);
            });
  for (int i = 0; i < l; ++i)
  {
    Evals(i)     = std::get<0>(eigen_vectors_and_values.at(i));
    Evecs.col(i) = std::get<1>(eigen_vectors_and_values.at(i));
  }
}

SECONDA_Output
SECONDA_analysis(Molecule &CurrMol,
                 Structure *CurrStruc,
                 BasicInformation &baseInformation,
                 Flags &flags)
{
#pragma omp critical
  {
    *std::cin.tie() << "\t\tStarting SECONDA analysis for "
                    << CurrMol.getLabel() << std::endl;
  }

  // some general variables
  Eigen::MatrixXd RDCnorm = CurrStruc->getRDCmatrix(rdcMatrixOptions::Norm);
  Eigen::MatrixXd RDC     = CurrStruc->getRDCmatrix(rdcMatrixOptions::Unscaled);
  Eigen::MatrixXd R_2, m, b;
  full_Regression(RDC, R_2, m, b);

  // Calculate the covariance matrix
  Eigen::MatrixXd covarianceMatrix =
      SECONDA_CovarianceMatrix(RDCnorm, CurrMol, flags);

  // Calculate eigen system

  // Compute the eigen system of the covariance matrix
  Eigen::EigenSolver<Eigen::MatrixXd> es(covarianceMatrix);

  // Read eigen vectos (column)
  Eigen::MatrixXd Evecs = es.eigenvectors().real();

  // Read eigen values
  Eigen::VectorXd Evals = es.eigenvalues().real();
  sort_evals(Evals, Evecs);
  for (int i = 0; i < Evals.rows(); ++i)
    Evals(i) = (Evals(i) < baseInformation.limits.zero_cutoff ? 0.0 : Evals(i));

  // Calculate Kappa
  Eigen::VectorXd KappaQ = SECONDA_Kappa(Evecs);

  Eigen::VectorXd heterogeneity = SECONDA_a_2(Evecs, Evals);

  return generate_Seconda_Output(RDCnorm, Evals, Evecs, KappaQ, R_2, m, b,
                                 heterogeneity, baseInformation);
}

Eigen::MatrixXd
SECONDA_CovarianceMatrix(Eigen::MatrixXd &RDCnorm,
                         Molecule &CurrMol,
                         Flags &flags)
{ // TODO: Implement if statement
  if (flags.SECONDAreducedCovariance)
  {
    return SECONDA_CovarianceMatrix_reduced(RDCnorm, CurrMol);
  }
  else
  {
    return SECONDA_CovarianceMatrix_Hus(RDCnorm);
  }
}

//#ifdef WEIGHTED_COVARIANCE_MATRIX
Eigen::MatrixXd
SECONDA_CovarianceMatrix_reduced(Eigen::MatrixXd &RDCnorm, Molecule &CurrMol)
{
  // Get initial number of RDCs per set and number of sets.
  const unsigned int NOR = CurrMol.getNOR();
  const unsigned int NOS = CurrMol.getNORsets();
  unsigned int idx_s, idx_r;
  unsigned int i, j;

  // Initialize the vectors to store the mean RDCs.
  Eigen::VectorXd D_K_reduced = Eigen::VectorXd::Zero(NOR);
  Eigen::VectorXd D_M_reduced = Eigen::VectorXd::Zero(NOS);

  // Initialize the vectors to store the reduced number of RDCs per set and
  // number of sets.
  Eigen::Matrix<unsigned int, Eigen::Dynamic, 1> NOS_reduced =
      Eigen::Matrix<unsigned int, Eigen::Dynamic, 1>::Zero(NOR);
  Eigen::Matrix<unsigned int, Eigen::Dynamic, 1> NOR_reduced =
      Eigen::Matrix<unsigned int, Eigen::Dynamic, 1>::Zero(NOS);

  // Initialize matrix that counts the number of elements used for the
  // individual covariance matrix elements.
  Eigen::Matrix<unsigned int, Eigen::Dynamic, Eigen::Dynamic> ELEM_reduced =
      Eigen::Matrix<unsigned int, Eigen::Dynamic, Eigen::Dynamic>::Zero(NOR,
                                                                        NOR);

  // Initialize the vector of weights and covariance matrix.
  Eigen::VectorXd weights = Eigen::VectorXd::Zero(NOS);
  Eigen::MatrixXd covar   = Eigen::MatrixXd::Zero(NOR, NOR);

  RDCset *CurrSet    = CurrMol.getHeadSet();
  RDCdata *CurrRDC_1 = NULL, *CurrRDC_2 = NULL;

  // Calculate mean RDCs (row and column mean)
  // Loop all RDC sets
  for (idx_s = 0; idx_s < NOS && CurrSet; ++idx_s, CurrSet = CurrSet->getNext())
  {
    CurrRDC_1 = CurrSet->getHeadData();
    // Loop all RDCs
    for (idx_r = 0; idx_r < NOR && CurrRDC_1;
         ++idx_r, CurrRDC_1 = CurrRDC_1->getNext())
    {
      if (CurrRDC_1->isUndefined())
        continue;
      ++NOR_reduced(idx_s);
      ++NOS_reduced(idx_r);
      D_K_reduced(idx_r) += RDCnorm(idx_r, idx_s);
      D_M_reduced(idx_s) += RDCnorm(idx_r, idx_s);
    }
  }

  for (idx_s = 0; idx_s < NOS; ++idx_s)
    D_M_reduced(idx_s) /= NOR_reduced(idx_s);
  for (idx_r = 0; idx_r < NOR; ++idx_r)
    D_K_reduced(idx_r) /= NOS_reduced(idx_r);

  // Calculate weights  //TODO Weighting auch in weighting einsetzen // I dont
  // know what this "todo" should tell me... I think it was fixed.
  for (idx_s = 0, CurrSet = CurrMol.getHeadSet(); idx_s < NOS && CurrSet;
       ++idx_s, CurrSet   = CurrSet->getNext())
  {
    for (idx_r = 0, CurrRDC_1                         = CurrSet->getHeadData();
         idx_r < NOR && CurrRDC_1; ++idx_r, CurrRDC_1 = CurrRDC_1->getNext())
    {
      if (CurrRDC_1->isUndefined())
        continue;
      weights(idx_s) += pow((RDCnorm(idx_r, idx_s) - D_M_reduced(idx_s)), 2.0);
    }
    weights(idx_s) =
        1.0 / ((1.0 / (NOR_reduced(idx_s) - 1.0)) * weights(idx_s));
  }

  CurrSet = CurrMol.getHeadSet();
  for (idx_s = 0; idx_s < NOS && CurrSet; ++idx_s, CurrSet = CurrSet->getNext())
  {
    CurrRDC_1 = CurrSet->getHeadData();
    for (i = 0; i < NOR && CurrRDC_1; ++i, CurrRDC_1 = CurrRDC_1->getNext())
    {
      if (CurrRDC_1->isUndefined())
        continue;
      CurrRDC_2 = CurrRDC_1;
      for (j = i; j < NOR && CurrRDC_2; ++j, CurrRDC_2 = CurrRDC_2->getNext())
      {
        if (CurrRDC_2->isUndefined())
          continue;
        ++ELEM_reduced(i, j);
        covar(i, j) += (weights(idx_s) * (RDCnorm(i, idx_s) - D_K_reduced(i)) *
                        (RDCnorm(j, idx_s) - D_K_reduced(j)));
      }
    }
  }

  // Normalize the covar elements and fill the matrix.
  for (i = 0; i < NOR; ++i)
  {
    for (j = i; j < NOR; ++j)
    {
      covar(i, j) /= (static_cast<double>(ELEM_reduced(i, j) - 1));
      covar(j, i) = covar(i, j);
    }
  }
  return covar;
}

Eigen::MatrixXd
SECONDA_CovarianceMatrix_Hus(Eigen::MatrixXd &RDCnorm)
{
  unsigned int i, j, m, v;

  unsigned int NOR = RDCnorm.rows();
  unsigned int NOS = RDCnorm.cols();

  double one_over_NOS_minus_1 = 1.0 / (((double) NOS) - 1.0);
  double NOR_minus_1          = ((double) NOR - 1.0);

  Eigen::MatrixXd covarianceMatrix = Eigen::MatrixXd::Zero(NOR, NOR);
  Eigen::VectorXd meanD            = Eigen::VectorXd::Zero(NOS);
  Eigen::VectorXd wD               = Eigen::VectorXd::Zero(NOS);

  for (m = 0; m < NOS; ++m)
    meanD(m) = RDCnorm.col(m).mean();

  for (m = 0; m < NOS; ++m)
  {
    for (v = 0; v < NOR; ++v)
    {
      wD(m) += pow((RDCnorm(v, m) - meanD(m)), 2.0);
    }
    wD(m) = (1.0 / wD(m));
  }

  wD *= NOR_minus_1;

  meanD = Eigen::VectorXd::Zero(NOR);

  for (m = 0; m < NOR; ++m)
    meanD(m) = RDCnorm.row(m).mean();

  for (i = 0; i < NOR; ++i)
  {
    for (j = i; j < NOR; ++j)
    {
      for (m = 0; m < NOS; ++m)
      {
        covarianceMatrix(i, j) +=
            (wD(m) * (RDCnorm(i, m) - meanD(i)) * (RDCnorm(j, m) - meanD(j)));
      }
      covarianceMatrix(j, i) = covarianceMatrix(i, j);
    }
  }

  covarianceMatrix *= one_over_NOS_minus_1;

  return covarianceMatrix;
}

Eigen::VectorXd
SECONDA_a_2(Eigen::MatrixXd &Evecs, Eigen::VectorXd &Evals)
{
  Eigen::VectorXd a_2 = Eigen::VectorXd::Zero(Evecs.cols());
  for (unsigned int j = 0; j < Evecs.cols(); ++j)
  {
    for (unsigned int q = 5; q < Evecs.cols(); ++q)
    {
      a_2(j) += (Evecs(j, q) * Evecs(j, q) * Evals(q));
    }
  }
  return a_2;
}

Eigen::VectorXd
SECONDA_Kappa(Eigen::MatrixXd &Evecs)
{
  unsigned int NOR       = Evecs.cols();
  double centum_over_NOR = (100.0 / ((double) NOR));
  unsigned int j, q;
  double efac;

  Eigen::VectorXd KappaQ = Eigen::VectorXd::Zero(NOR);

  for (q = 0; q < NOR; ++q)
  {
    efac = .0;
    for (j = 0; j < NOR; ++j)
    {
      efac += (pow(Evecs(j, q), 2.0) * log(pow(Evecs(j, q), 2.0)));
    }
    KappaQ(q) = exp(-efac) * centum_over_NOR;
  }

  return KappaQ;
}

double
SECONDA_rho(Eigen::MatrixXd &Evals)
{
  double rho = Evals(4, 0) / Evals(5, 0);

  return rho;
}

double
getConditionNumber(Eigen::MatrixXd &Mat,
                   unsigned int &rank,
                   double cutoff,
                   unsigned int base)
{
  Eigen::VectorXd C = getSingularValues(Mat);
  if (rank == 0 || rank > C.rows())
    rank = C.rows();
  while (rank > 0 && C(rank - 1) < cutoff)
    --rank;
  if (rank == 0)
    return .0;

  return (C(base) / C(rank - 1));
}

Eigen::VectorXd
getSingularValues(Eigen::MatrixXd &Mat)
{
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(Mat, Eigen::ComputeFullU |
                                                 Eigen::ComputeFullV);
  return svd.singularValues();
}

Eigen::MatrixXd
getTwoSetSingularValues(Eigen::MatrixXd &D)
{
  Eigen::MatrixXd SV  = Eigen::MatrixXd::Zero(D.cols(), D.cols());
  Eigen::MatrixXd D2C = Eigen::MatrixXd::Zero(D.rows(), 2);
  Eigen::VectorXd sigma;
  int i, j;
  for (i = 0; i < D.cols(); ++i)
  {
    for (j = (i + 1); j < D.cols(); ++j)
    {
      D2C.col(0) = D.col(i);
      D2C.col(1) = D.col(j);
      sigma      = getSingularValues(D2C);
      SV(i, j)   = sigma(0);
      SV(j, i)   = sigma(1);
    }
  }

  return SV;
}

SECONDA_Output
generate_Seconda_Output(Eigen::MatrixXd &RDCnorm,
                        Eigen::VectorXd &Evals,
                        Eigen::MatrixXd &Evecs,
                        Eigen::VectorXd &KappaQ,
                        Eigen::MatrixXd &R_2,
                        Eigen::MatrixXd &m,
                        Eigen::MatrixXd &b,
                        Eigen::VectorXd &a_2,
                        BasicInformation &baseInformation)
{
  SECONDA_Output out;
  out.Tolman_singular_values = getSingularValues(RDCnorm);
  out.kappa_q                = KappaQ;
  out.Covariance_Evals       = Evals;
  out.Covariance_Evecs       = Evecs;
  out.cumulative_variance    = Eigen::VectorXd::Zero(Evals.rows());

  double esum = Evals.sum();
  double csum = .0;

  // Get rank of covariance matrix
  out.covar_rank = Evals.size();
  while (out.covar_rank &&
         Evals(out.covar_rank - 1) < baseInformation.limits.zero_cutoff)
    --out.covar_rank;

  // Calculate cumulative sums
  for (int i = 0; i < out.cumulative_variance.rows(); ++i)
  {
    csum += Evals(i);
    out.cumulative_variance(i) = 100 * csum / esum;
  }

  // Calculate rho gaps
  out.rho_5_6    = (Evals(5) == .0 ? .0 : Evals(4) / Evals(5));
  out.rho_mean_6 = .0;

  for (int i = 0; i < 5; ++i)
    out.rho_mean_6 += Evals(i);
  if (Evals(5) > .0)
    out.rho_mean_6 *= (0.2 / Evals(5));
  else
    out.rho_mean_6 = .0;

  // Calculate "Condition number" (1-5 gap)
  out.covar_condition_number =
      Evals(0) / Evals((out.covar_rank > 4 ? 4 : (out.covar_rank - 1)));

  // Get Tolman Condition number and rank
  out.Tolman_CN_rank = 5;
  out.Tolman_rank    = 0;
  // Determine full rank
  getConditionNumber(RDCnorm, out.Tolman_rank,
                     baseInformation.limits.zero_cutoff, 0);
  // Determine 1-5 gap
  out.Tolman_condition_number = getConditionNumber(
      RDCnorm, out.Tolman_CN_rank, baseInformation.limits.zero_cutoff, 0);

  out.RDC_R_2              = R_2;
  out.RDC_m                = m;
  out.RDC_b                = b;
  out.twoSetSingularValues = getTwoSetSingularValues(RDCnorm);
  out.heterogeneity        = a_2;
  out.collected            = true;
  return out;
}

void
full_Regression(Eigen::MatrixXd &RDCnorm,
                Eigen::MatrixXd &R,
                Eigen::MatrixXd &m,
                Eigen::MatrixXd &b)
{
  unsigned int NOS = RDCnorm.cols();
  unsigned int i, j, I, J;

  R = Eigen::MatrixXd::Zero((NOS - 1), (NOS - 1));
  m = Eigen::MatrixXd::Zero((NOS - 1), (NOS - 1));
  b = Eigen::MatrixXd::Zero((NOS - 1), (NOS - 1));
  Eigen::MatrixXd x, y;

  for (I = i = 0; i < NOS; ++i, ++I)
  {
    x = RDCnorm.col(i);
    for (J = I, j = (i + 1); j < NOS; ++j, ++J)
    {
      y       = RDCnorm.col(j);
      R(I, J) = linear_Regression(x, y, m(I, J), b(I, J));
    }
  }
}

void
rdc_vector_analysis(Structure *CurrStruc, StructureOptions opt)
{
  unsigned int NOR = CurrStruc->getParent()->getNOR(), i = 0;
  int min, max = min = 0;

  Eigen::MatrixXd r_mat = Eigen::MatrixXd::Zero(NOR, 3), r_cur, C;
  Eigen::Vector3d eval;

  for (SphericalHarmonics *Y = CurrStruc->getHeadYmatrix()->getHeadHarmonic();
       Y; Y                  = Y->getNext(), ++i)
  {
    r_cur = Y->getAtom1()->Coordinates2Eigen(opt) -
            Y->getAtom2()->Coordinates2Eigen(opt);
    r_mat.row(i) = (r_cur.transpose() / r_cur.norm());
  }

  C = r_mat.transpose() * r_mat;
  Eigen::EigenSolver<Eigen::MatrixXd> es(C, false);
  eval = es.eigenvalues().real();
  eval.maxCoeff(&max);
  eval.minCoeff(&min);

  CurrStruc->set_rdc_vector_sampling(
      sqrt(eval(max) / eval(min)),
      ((eval((3 - min - max)) - eval(min)) / eval(max)), opt);
}


void
SECONDA_sensitivity(Molecule &CurrMol,
                    Structure *CurrStruc,
                    BasicInformation &baseInformation,
                    Flags &flags)
{
  unsigned int NOS = CurrMol.getNORsets(), NOR = CurrMol.getNOR();
  unsigned int currSet, set, set_full;
  // some general variables
  Eigen::MatrixXd RDCnorm_full =
      CurrStruc->getRDCmatrix(rdcMatrixOptions::Norm);
  Eigen::MatrixXd RDCnorm = Eigen::MatrixXd::Zero(NOR, (NOS - 1));
  RDCset *CurrSet         = CurrMol.getHeadSet();
  for (currSet = 0; currSet < NOS; ++currSet)
  {
    for (set_full = set = 0; set_full < NOS; ++set_full)
    {
      if (set_full == currSet)
        continue;
      RDCnorm.col(set) = RDCnorm_full.col(set_full);
      ++set;
    }

    // Calculate the covariance matrix
    Eigen::MatrixXd covarianceMatrix =
        SECONDA_CovarianceMatrix(RDCnorm, CurrMol, flags);

    // Compute the eigen system of the covariance matrix
    Eigen::EigenSolver<Eigen::MatrixXd> es(covarianceMatrix);

    // Read eigen values
    Eigen::VectorXd Evals = es.eigenvalues().real();
    Eigen::MatrixXd Evecs = es.eigenvectors().real();
    sort_evals(Evals, Evecs);
    for (int i = 0; i < Evals.size(); ++i)
      Evals(i) =
          (Evals(i) < baseInformation.limits.zero_cutoff ? .0 : Evals(i));
    int lead_index = 4;
    while (Evals(lead_index) == .0 && lead_index)
      --lead_index;
    CurrSet->set_SECONDA_gap_5_6_sensitivity(
        (Evals(5) ? Evals(4) / Evals(5) : .0));
    CurrSet->set_SECONDA_gap_1_5_sensitivity(Evals(0) / Evals(lead_index));
    CurrSet->set_SECONDA_sensitivity_rank(lead_index + 1);
    CurrSet = CurrSet->getNext();
  }
}
