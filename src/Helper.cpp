
#include <Atom.hpp>
#include <Declarations.hpp>
#include <Helper.hpp>
#include <LinAlg.hpp>
#include <Molecule.hpp>
#include <Parser/Parser.hpp>
#include <RDCdata.hpp>
#include <RDCset.hpp>
#include <SphericalHarmonics.hpp>
#include <Structure.hpp>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <omp.h>
#include <unistd.h>

#define _ANS_SIZE_ 50

unsigned int
identifyKey(std::string s)
{
  unsigned int i;
  /* To be skipped */
  if (!s.compare("aligncount"))
    return 666;

  for (i = 0; i < NumberOfKeywords; ++i)
  {
    if (s.compare(0, ListOfKeywords[i].keyword.size(),
                  ListOfKeywords[i].keyword))
      continue;
    return i;
  }
  return 0;
}

template<typename T>
double
zetaFunction(T a, T m, T n)
{
  return (deltaFunction(a, m) - deltaFunction(a, n));
}

template double zetaFunction(Atom *, Atom *, Atom *);
template double zetaFunction(double *, double *, double *);

template<typename T>
double
deltaFunction(T a, T b)
{
  return (a == b ? 1.0 : 0.0);
}

template double deltaFunction(Atom *, Atom *);
template double deltaFunction(double *, double *);

template<typename T>
T
max(T a, T b)
{
  return (a > b ? a : b);
}

template double max(double, double);

template<typename T>
T
min(T a, T b)
{
  return (a < b ? a : b);
}

template double min(double, double);

double
getPolar2Beta(double theta1, double phi1, double theta2, double phi2)
{
  double beta        = .0;
  Eigen::Vector3d V1 = Polar2Eigen(theta1, phi1);
  Eigen::Vector3d V2 = Polar2Eigen(theta2, phi2);
  beta               = acos(V1.transpose() * V2);
  return beta;
}

double
getDpolar2Dbeta(double theta1,
                double dT1,
                double phi1,
                double dP1,
                double theta2,
                double dT2,
                double phi2,
                double dP2)
{
  double dB, dcB = .0;
  double dBdT1, dBdT2, dBdP1, dBdP2;
  double cosT1 = cos(theta1);
  double cosT2 = cos(theta2);
  double cosP1 = cos(phi1);
  double cosP2 = cos(phi2);
  double sinT1 = sin(theta1);
  double sinT2 = sin(theta2);
  double sinP1 = sin(phi1);
  double sinP2 = sin(phi2);

  // Stepwise computing of the partial derivatives of beta in respect to the
  // polar angles

  // D(cos(b)) / D(t1) =
  dBdT1 = cosT1 * cosP1 * sinT2 * cosP2 + cosT1 * sinP1 * sinT2 * sinP2 -
          sinT1 * cosT2;
  // D(cos(b)) / D(t2) =
  dBdT2 = sinT1 * cosP1 * cosT2 * cosP2 + sinT1 * sinP1 * cosT2 * sinP2 -
          cosT1 * sinT2;
  // D(cos(b)) / D(p1) =
  dBdP1 = -sinT1 * sinP1 * sinT2 * cosP2 + sinT1 * cosP1 * sinT2 * sinP2;
  // D(cos(b)) / D(p2) =
  dBdP2 = -sinT1 * cosP1 * sinT2 * sinP2 + sinT1 * sinP1 * sinT2 * cosP2;
  dBdT1 *= dT1;
  dBdP1 *= dP1;
  dBdT2 *= dT2;
  dBdP2 *= dP2;

  dcB = sqrt(pow(dBdT1, 2.0) + pow(dBdT2, 2.0) + pow(dBdP1, 2.0) +
             pow(dBdP2, 2.0));

  // D( arcos(dcB) ) / D( dcB ) = 1 / sqrt ( 1.0 - pow ( dcB, 2.0 ) );
  dB = dcB /
       (sqrt(1.0 - pow(cos(getPolar2Beta(theta1, phi1, theta2, phi2)), 2.0)));
  return dB;
}


void
addTimeSlot(BasicInformation &baseInformation, double &slot)
{
  slot += (omp_get_wtime() - baseInformation.timeStamp);
  baseInformation.timeStamp = omp_get_wtime();
}

void
checkWorkingDir(BasicInformation &baseInformation, Flags &flags)
{
  if (baseInformation.inputFileName.find("/") != std::string::npos)
  {
    unsigned int i = baseInformation.inputFileName.length() - 1;
    while (baseInformation.inputFileName.at(i) != '/')
      --i;
    baseInformation.workingDirectory =
        baseInformation.inputFileName.substr(0, ++i);
  }
  flags.noInput = false;
  baseInformation.numberOfInputs++;
}

std::string
getCurrentDir(BasicInformation &baseInformation)
{
  return (getcwd(baseInformation.rundir, STANDARD_BUFFER_SIZE) ?
              std::string(baseInformation.rundir) :
              std::string(""));
}

void
printKeywords(BasicInformation &baseInformation, Flags &flags)
{
  std::cout << " *  echo: " << flags.echo << std::endl;
  std::cout << " *  lmmaxiterations: "
            << baseInformation.limits.max_lm_iterations << std::endl;
  std::cout << " *  maxtitaniaiterations: "
            << baseInformation.limits.max_titania_iterations << std::endl;
  std::cout << " *  skipscrm: " << flags.skipSCRM << std::endl;
  std::cout << " *  silent: " << flags.silent << std::endl;
  std::cout << " *  qfactorconvergence: "
            << baseInformation.limits.Q_factor_convergence << std::endl;
  std::cout << " *  titania2hotfcht: " << flags.titania2hotfcht << std::endl;
  std::cout << " *  outputali: " << flags.outputAli << std::endl;
  std::cout << " *  threads: " << baseInformation.numOfThreads << std::endl;
  std::cout << " *  montecarlobootstrapping: " << flags.monteCarloBootstrapping
            << std::endl;
  std::cout << " *  plotkappaq: " << flags.plotKappaQ << std::endl;
  std::cout << " *  plotTrajectory: " << flags.plotTrajectory << std::endl;
  std::cout << " *  gnuoutputformat: " << baseInformation.gnuOutputFormat
            << std::endl;
  std::cout << " *  numericalgradients: " << flags.numericalGradients
            << std::endl;
  std::cout << " *  outputlm: " << flags.outputLM << std::endl;
  std::cout << " *  redundantsconvergence: "
            << baseInformation.limits.redundants_convergence << std::endl;
  std::cout << " *  meanalignmentconvergence: "
            << baseInformation.limits.alignment_mean_convergence << std::endl;
  std::cout << " *  sigmaalignmentconvergence: "
            << baseInformation.limits.alignment_sigm_convergence << std::endl;
  std::cout << " *  meanangleconvergence: "
            << baseInformation.limits.sphericals_mean_convergence << std::endl;
  std::cout << " *  sigmaangleconvergence: "
            << baseInformation.limits.sphericals_sigm_convergence << std::endl;
  std::cout << " *  spreadangleconvergence: "
            << baseInformation.limits.sphericals_spread_convergence
            << std::endl;
}

int
loadFCHTname(BasicInformation &bI)
{
  std::string out = std::string(".out");
  int outlength   = out.size();
  if (bI.outputFileName.substr(bI.outputFileName.size() - outlength,
                               outlength) == out)
    bI.hotFCHTbase =
        bI.outputFileName.substr(0, (bI.outputFileName.size() - outlength));
  else
    bI.hotFCHTbase = bI.outputFileName;
  return 0;
}

/*int generate_rdc_rmsd_plot ( BasicInformation &bI )
{

   return 0;
}*/

template<typename T>
bool
vector_contains(std::vector<T> vec, T el)
{
  return (std::find(vec.begin(), vec.end(), el) != vec.end());
}

template bool vector_contains(std::vector<Atom *>, Atom *);

template<typename T>
int
get_vector_index(std::vector<T> vec, T el)
{
  auto it = std::find(vec.begin(), vec.end(), el);
  if (it == vec.end())
    return -1;
  else
    return (std::distance(vec.begin(), it));
}

template int get_vector_index(std::vector<Atom *>, Atom *);

unsigned int
getSizeOfUpperTriangle(const unsigned int n)
{
  return n * (n - 1) / 2;
}

unsigned long long
getTotalSystemMemory()
{
  long pages     = sysconf(_SC_PHYS_PAGES);
  long page_size = sysconf(_SC_PAGE_SIZE);
  return pages * page_size;
}
