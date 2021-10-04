
#include <Atom.hpp>
#include <Declarations.hpp>
#include <Helper.hpp>
#include <IndependencyAnalysis.hpp>
#include <LinAlg.hpp>
#include <Molecule.hpp>
#include <MonteCarloStatistic.hpp>
#include <Output.hpp>
#include <Parser/Parser.hpp>
#include <RDCdata.hpp>
#include <RDCset.hpp>
#include <RedundantInternals.hpp>
#include <ShoemakeEulerAngles.hpp>
#include <SphericalHarmonics.hpp>
#include <SphericalHarmonicsMatrix.hpp>
#include <Structure.hpp>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <sys/resource.h>
#include <unistd.h>
#include <version.hpp>

#ifdef USE_CUDA
#include <cuda_helper.hpp>
#endif

constexpr auto INFO_BUF = OUTPUT_INFOBOX_STRING_SIZE * sizeof(char);
constexpr auto STAN_BUF = STANDARD_BUFFER_SIZE * sizeof(char);

const int header_width                 = 60;
const int sub_block_size               = 6;
const int seconda_sub_block_size       = 14;
const int number_of_authors            = 3;
const char *authors[number_of_authors] = {
    "Felix A. Roth",
    "Volker Schmidts",
    "Christina M. Thiele",
};

const char *references[] = {
    "J. Meiler, J.J. Prompers, W. Peti, C. Griesinger, R. Brüschweiler, JACS "
    "2001, 123, 6098 - 6107 (https://doi.org/10.1021/ja010002z).",
    "W. Peti, J. Meiler, R. Brüschweiler, C. Griesinger, JACS 2002, 124, 5822 "
    "- 5833 (https://doi.org/10.1021/ja011883c).",
    "N.A. Lakomek, R. Brüschweiler, J. Meiler, C. Griesinger, J. Biomol. NMR "
    "2008, 41, 139 - 155 (https://doi.org/10.1007/s10858-008-9244-4).",
    "J.C. Hus, R. Brüschweiler, J. Biomol. NMR 2002, 24, 123 - 132 "
    "(https://doi.org/10.1023/A:1020927930910).",
    "R. Brüschweiler, J. Chem. Phys. 1995, 102, 3396-3403 "
    "(https://doi.org/10.1063/1.469213).",
    "J.C. Hus, W. Peti, C. Griesinger, R. Brüschweiler, JACS 2003, 125, "
    "5596-5597 (https://doi.org/10.1021/ja029719s).",
    "J.R. Tolman, JACS 2002, 124, 12020-12030 "
    "(https://doi.org/10.1021/ja0261123).",
    "J.A. Losonczi, M. Andrec, M.W.F. Fischer, J.H. Prestegard, J. of Mag.Res. "
    "1999, 138, 334-342 (https://doi.org/10.1006/jmre.1999.1754).",
    "K.V. Mardia, P.E. Jupp, \"Directional Statistics\" 1999, Wiley Series in "
    "Probability and Statistics (https://doi.org/10.1002/9780470316979).",
    "V. Bakken, T. Helgaker, J. Chem. Phys. 2002, 117, 9160-9174 "
    "(https://doi.org/10.1063/1.1515483).",
    "T.A. Halgren, J. Comp. Chem. 1996, 17, 490 - 519 "
    "(https://doi.org/10.1002/(SICI)1096-987X(199604)17:5/"
    "6<490::AID-JCC1>3.0.CO;2-P).",
    "C. Cornilescu, J.L. Marquardt, M. Ottiger, A. Bax, JACS 1998, 120, "
    "6836-6837 (https://doi.org/10.1021/ja9812610).",
};


namespace Paper {
  enum paper
  {
    Meiler = 0,
    Peti,          // = 1,
    Lakomek,       // = 2,
    Hus1,          // = 3,
    Brueschweiler, // = 4,
    Hus2,          // = 5,
    Tolman,        // = 6,
    Losonczi,      // = 7,
    Mardia,        // = 8,
    Bakken,        // = 9,
    Halgren,       // = 10,
    Cornilescu,    // = 11,
  };

  inline const char *print(paper p)
  {
    return references[p];
  }
  inline const char *cite(paper p)
  {
    char *buf = (char *) malloc(STAN_BUF);
    snprintf(buf, STAN_BUF, "[%d]", static_cast<int>(p + 1));
    return buf;
  }
  inline const char *cite_multiple_start(paper p)
  {
    char *buf = (char *) malloc(STAN_BUF);
    snprintf(buf, STAN_BUF, "[%d", static_cast<int>(p + 1));
    return buf;
  }
  inline const char *cite_multiple_middle(paper p)
  {
    char *buf = (char *) malloc(STAN_BUF);
    snprintf(buf, STAN_BUF, ",%d", static_cast<int>(p + 1));
    return buf;
  }
  inline const char *cite_multiple_end(paper p)
  {
    char *buf = (char *) malloc(STAN_BUF);
    snprintf(buf, STAN_BUF, ",%d]", static_cast<int>(p + 1));
    return buf;
  }
} // namespace Paper


void frame_string(FILE *output, const char name[], int pre_space = 3);

const char *
to_cstring(double input, int width, int precision, bool exponential)
{
  char *buf = (char *) malloc(STAN_BUF);
  if (exponential)
    snprintf(buf, STAN_BUF, "%*.*e", width, precision, input);
  else
    snprintf(buf, STAN_BUF, "%*.*f", width, precision, input);
  return buf;
}

const char *
to_cstring(int input, int width)
{
  char *buf = (char *) malloc(STAN_BUF);
  snprintf(buf, STAN_BUF, "%*d", width, input);
  return buf;
}

const char *
to_cstring(unsigned input, int width)
{
  char *buf = (char *) malloc(STAN_BUF);
  snprintf(buf, STAN_BUF, "%*u", width, input);
  return buf;
}

const char *
to_cstring(long int input, int width)
{
  char *buf = (char *) malloc(STAN_BUF);
  snprintf(buf, STAN_BUF, "%*ld", width, input);
  return buf;
}

const char *
to_cstring(long unsigned input, int width)
{
  char *buf = (char *) malloc(STAN_BUF);
  snprintf(buf, STAN_BUF, "%*lu", width, input);
  return buf;
}

const char *
LM_stop(int reason)
{
  switch (reason)
  {
    case 1:
      return "delta(p)";
    case 2:
      return "grad(p)";
    case 3:
      return "epsilon^2(p)";
    case 4:
      return "max iter";
    case 5:
      return "ini. grad(p)";
    default:
      return "unkown";
  }
}

const Eigen::IOFormat big_data(Eigen::StreamPrecision,
                               Eigen::DontAlignCols,
                               ";",
                               "\n",
                               "",
                               "",
                               "",
                               "");

/*******************************************************
 *                                                     *
 * void shutup                                         *
 *                                                     *
 * TITANIA generates a terminal output by default. If the *
 *  user would like to redirect this output to a file  *
 *    or to suppress the output shutup is called and   *
 *             handles the respective request.         *
 *                                                     *
 *******************************************************/

void
shutup(BasicInformation &baseInformation)
{
  if (baseInformation.comFileName == "")
    baseInformation.comFileName = ".null";
  baseInformation.comFile.open(baseInformation.comFileName,
                               std::ios::binary | std::ios::out |
                                   std::ios::trunc);
  std::cin.tie(&baseInformation.comFile);
}

void
frame_string(FILE *output, const char name[], int pre_space)
{
  const int name_length = static_cast<int>(strlen(name));
  const int suf_space   = header_width - name_length - pre_space - 2;
  fprintf(output, " =%*s%s%*s= \n", pre_space, "", name, suf_space, "");
}

void
information_header(FILE *output, const char name[])
{
  fprintf(output, "\n");
  fprintf(output, " %s \n", std::string(header_width, '=').c_str());
  frame_string(output, name);
  fprintf(output, " %s \n", std::string(header_width, '=').c_str());
  fprintf(output, "\n");
  fflush(output);
}

void
add_contributors(FILE *output)
{
  fprintf(output, " %s \n", " Main contributors:");
  for (int i = 0; i < number_of_authors; ++i)
    fprintf(output, "   * %s\n", authors[i]);
  fflush(output);
}

void
add_externals(FILE *output, BasicInformation &baseInformation, Flags &flags)
{
  fprintf(output, " %s\n", "This TITANIA version uses:");
  fprintf(output, "   * %-*s: %s\n", baseInformation.output_format.string_xl,
          "EIGEN3", "Standard vector & matrix operations.");
  fprintf(output, "   * %-*s: %s\n", baseInformation.output_format.string_xl,
          "CBLAS", "Plugin vector & matrix operations.");
  fprintf(output, "   * %-*s: %s\n", baseInformation.output_format.string_xl,
          "LAPACKE", "Plugin linear algebra routines.");
#ifdef USE_CUDA
  if (flags.use_gpu)
  {
    fprintf(output, "   * %-*s: %s\n", baseInformation.output_format.string_xl,
            "CUDART", "CUDA run time implementation.");
    fprintf(output, "   * %-*s: %s\n", baseInformation.output_format.string_xl,
            "CUBLAS", "CUDA vector & matrix operations.");
    fprintf(output, "   * %-*s: %s\n", baseInformation.output_format.string_xl,
            "CUSOLVER", "CUDA linear algebra routines.");
  }
#endif
  if (baseInformation.numOfThreads)
  {
    fprintf(output, "   * %-*s: %s\n", baseInformation.output_format.string_xl,
            "OMP", "Parallelization standard.");
  }
  fflush(output);
}

void
add_compilation_information(FILE *output, BasicInformation &baseInformation)
{
  fprintf(output, "   * %-*s: %s\n", baseInformation.output_format.string_xl,
          "Version", TITANIA_VERSION__);
  fprintf(output, "   * %-*s: %s\n", baseInformation.output_format.string_xl,
          "Compile date", __DATE__);
  fprintf(output, "   * %-*s: %s\n", baseInformation.output_format.string_xl,
          "Compile host", COMP_HOST__);
  fprintf(output, "   * %-*s: %s\n", baseInformation.output_format.string_xl,
          "Compiler version", COMP_CXX_VERSION__);
  if (output == baseInformation.output)
    fprintf(output, "\n");
  fflush(output);
}

void
add_host_information(FILE *output,
                     BasicInformation &baseInformation,
                     Flags &flags)
{
  fprintf(output, "   * %-*s: %s\n", baseInformation.output_format.string_xl,
          "Hostname", baseInformation.hostname);
  fprintf(output, "   * %-*s: %s\n", baseInformation.output_format.string_xl,
          "Username", baseInformation.username);
  fprintf(output, "   * %-*s: %.1f\n", baseInformation.output_format.string_xl,
          "Sys. memory (GB)", getTotalSystemMemory() * pow(1024.0, -3));

#ifdef USE_CUDA
  if (flags.use_gpu)
  {
    if (flags.has_gpu)
    {
      fprintf(output, "   * %-*s: %s\n",
              baseInformation.output_format.string_xl, "GPU device",
              baseInformation.gpu_device_name.c_str());
      fprintf(output, "   * %-*s: %.1f\n",
              baseInformation.output_format.string_xl, "GPU memory (MB)",
              get_gpu_memory());
    }
    else
      fprintf(output, "   * No GPU device found... Remove TITANIA_CUDA=1 from "
                      "titania_configuration.cfg file...\n");
  }
#endif
  fflush(output);
}

void
add_input_file_information(FILE *output, BasicInformation &baseInformation)
{
  fprintf(output, " %-*s %s\n", baseInformation.output_format.string_xxl,
          "Input file name:", baseInformation.inputFileName.c_str());
  fprintf(output, " %-*s %s\n", baseInformation.output_format.string_xxl,
          "Working directory:", baseInformation.rundir);
  if (output == baseInformation.output)
    fprintf(output, "\n");
  fflush(output);
}

void
echoInput(FILE *output, InputFile &reader, BasicInformation &baseInformation)
{
  for (int line = 0; line < reader.get_number_of_lines(); ++line)
    fprintf(output, " |%*d> %s\n", baseInformation.output_format.number_s,
            line + 1, reader.get_line(line).c_str());
  fflush(output);
}

void
add_run_information(FILE *output, BasicInformation &baseInformation)
{
  fprintf(output, " %-*s: %s", baseInformation.output_format.string_xl,
          "Calculation started", ctime(&baseInformation.starting_time));
  fflush(output);
}

void
add_flag_information(FILE *output,
                     BasicInformation &baseInformation,
                     Flags &flags)
{
  fprintf(output, " User defined keywords:\n");

  if (flags.numericalGradients)
  {
    fprintf(output, "   * %s\n",
            "Levenberg-Marquardt jacobi matrix was estimated using numerical");
    fprintf(output, "   * %s %.*e\n",
            "gradients. Delta(x) = ", baseInformation.output_format.prec_m,
            baseInformation.numericalDeltaX);
  }

  if (flags.useRedundants)
  {
    if (baseInformation.limits.max_redundant_cycles !=
        STAN_MAX_REDUNDANT_CYCLES)
      fprintf(output, "   * %s: %d\n",
              "Maximum iterations for redundant internal coordinates",
              baseInformation.limits.max_redundant_cycles);

    if (baseInformation.useRedundantsOnlyAfter)
    {
      fprintf(output, "   * %s: ",
              "Redundant internal coordinates were exclusivly used after");
      fprintf(output, "%d %s\n", baseInformation.useRedundantsOnlyAfter,
              "steps");
    }

    if (flags.redundants_damping)
      fprintf(output, "   * %-*s %.*e\n",
              2 * baseInformation.output_format.string_xxl,
              "Used redundant coordinates damping constant",
              baseInformation.output_format.prec_m,
              baseInformation.redundants_damping);

    if (baseInformation.static_redundants_weighting[BOND_REDUNDANTS_] != 1.0)
      fprintf(output, "   * %-*s %.*e\n",
              2 * baseInformation.output_format.string_xxl,
              "Used static bond length damping constant",
              baseInformation.output_format.prec_m,
              baseInformation.static_redundants_weighting[BOND_REDUNDANTS_]);
    if (baseInformation.static_redundants_weighting[ANGLE_REDUNDANTS_] != 1.0)
      fprintf(output, "   * %-*s %.*e\n",
              2 * baseInformation.output_format.string_xxl,
              "Used static bond angle damping constant",
              baseInformation.output_format.prec_m,
              baseInformation.static_redundants_weighting[ANGLE_REDUNDANTS_]);
    if (baseInformation.static_redundants_weighting[TORSION_REDUNDANTS_] != 1.0)
      fprintf(output, "   * %-*s %.*e\n",
              2 * baseInformation.output_format.string_xxl,
              "Used static torsion angle damping constantt",
              baseInformation.output_format.prec_m,
              baseInformation.static_redundants_weighting[TORSION_REDUNDANTS_]);
    if (baseInformation.static_redundants_weighting[RDC_REDUNDANTS_] != 1.0)
      fprintf(output, "   * %-*s %.*e\n",
              2 * baseInformation.output_format.string_xxl,
              "Used static rdc angle damping constant",
              baseInformation.output_format.prec_m,
              baseInformation.static_redundants_weighting[RDC_REDUNDANTS_]);
    if (baseInformation.static_redundants_weighting[PLANAR_REDUNDANTS_] != 1.0)
      fprintf(output, "   * %-*s %.*e\n",
              2 * baseInformation.output_format.string_xxl,
              "Used static planar center damping constant",
              baseInformation.output_format.prec_m,
              baseInformation.static_redundants_weighting[PLANAR_REDUNDANTS_]);
    if (baseInformation.static_redundants_weighting[DISTANCE_REDUNDANTS_] !=
        1.0)
      fprintf(
          output, "   * %-*s %.*e\n",
          2 * baseInformation.output_format.string_xxl,
          "Used static distance damping constant",
          baseInformation.output_format.prec_m,
          baseInformation.static_redundants_weighting[DISTANCE_REDUNDANTS_]);

    if (flags.torsions_4_redundants)
      fprintf(output, "   * %s\n",
              "Redundant internal coordinates used torsion angles");
    if (flags.floating_rdc_angles)
      fprintf(output, "   * %s\n",
              "RDC-vector vs. static-vector angles were floating");
  }
  else
    fprintf(output, "   * %s\n",
            "Structures were generated using vector addtion");

  switch (baseInformation.redundants_distance_optimization)
  {
    case 0:
      fprintf(output, "   * %s\n",
              "No distance information were used during structure refinement");
      break;
    case 1:
      fprintf(
          output, "   * %s\n",
          "Distance information were fully used during structure refinement");
      fprintf(output, "     %s %d %s\n",
              "Exclusion volume induced inversions were reduced after",
              baseInformation.lowerInversionAfter, "iterations");
      break;
    case 2:
      fprintf(output, "   * %s\n",
              "Distance information were used during redundant coordinates "
              "structure refinement");
      break;
    case 3:
      fprintf(output, "   * %s\n",
              "Distance information were only used during vector addition.");
      fprintf(output, "     %s %d %s\n",
              "Exclusion volume induced inversions were reduced after",
              baseInformation.lowerInversionAfter, "iterations");
      break;
    default:
      break;
  }

  if (flags.skipEckart)
    fprintf(output, "   * %s\n", "Eckart transformation was skiped");

  if (!flags.calculateFullMatrix)
  {
    fprintf(output, "   * %s\n",
            "Calculations were performed in vector formalism");
    if (flags.errorWeightInSVD)
      fprintf(output, "   * %s\n",
              "Vector formalism used RDC error weighting in SVD");
  }

  if (flags.recalculateRDCs)
    fprintf(output, "   * %s\n",
            "Undefined RDCs were estimated for every optimization step");
  if (!flags.scaleWithSoverall)
    fprintf(output, "   * %s\n",
            "S(overall) scaling  of spherical harmonics Y was skiped");
  if (baseInformation.overOptimization)
    fprintf(output, "   * %s %d %s\n", "TITANIA used",
            baseInformation.overOptimization, "addtional iterations");
  fflush(output);
}

void
add_metrics(FILE *output, BasicInformation &baseInformation)
{
  fprintf(output, " %-*s %d\n", baseInformation.output_format.string_xxl,
          "Number of RDCs:", baseInformation.NumberOfRDCs);
  fprintf(output, " %-*s %d\n", baseInformation.output_format.string_xxl,
          "Number of RDC sets:", baseInformation.NumberOfSets);
  fprintf(output, " %-*s %d\n", baseInformation.output_format.string_xxl,
          "Number of atoms:", baseInformation.NumberOfAtoms);
  fflush(output);
}

void
add_redundant_metrics(FILE *output,
                      BasicInformation &baseInformation,
                      RedundantInternals &Red)
{
  fprintf(output, " %-*s %u\n", baseInformation.output_format.string_xxl,
          "Number of bonds:", Red.getNumberOfBonds());
  fprintf(output, " %-*s %u\n", baseInformation.output_format.string_xxl,
          "Number of bond angles:", Red.getNumberOfAngles());
  fprintf(output, " %-*s %u\n", baseInformation.output_format.string_xxl,
          "Number of dihedral angles:", Red.getNumberOfTorsions());
  fprintf(output, " %-*s %u\n", baseInformation.output_format.string_xxl,
          "Number of chiral volumes:", Red.getNumberOfChiralVolumes());
  fprintf(output, " %-*s %u\n", baseInformation.output_format.string_xxl,
          "Number of redundant coord.:", Red.getNumberOfRedundants());
  fprintf(output, " %-*s %u\n", baseInformation.output_format.string_xxl,
          "Number of distances:", Red.getNumberOfDistance());
  fflush(output);
}

void
add_rdc_sets(FILE *output,
             Molecule &CurrMol,
             BasicInformation &baseInformation,
             Flags &flags)
{
  fprintf(output, " List of RDC sets\n");
  fprintf(output, " %*s | %s\n", baseInformation.output_format.string_m,
          "Index", "Label");
  for (RDCset *CurrSet = CurrMol.getHeadSet(); CurrSet;
       CurrSet         = CurrSet->getNext())
  {
    fprintf(output, " %*d | %s\n", baseInformation.output_format.string_m,
            CurrSet->getIndex(), CurrSet->getLabel().c_str());
  }

  if (flags.minimalOutput)
    return;
  fprintf(output, "\n List of RDCs\n");
  fprintf(output, "    %-*s %-*s | %s\n",
          baseInformation.output_format.string_m, "Nuc 1",
          baseInformation.output_format.string_m, "Nuc 2",
          "Range (bonds between RDC nuclei)");
  for (RDCdata *CurrData = CurrMol.getHeadSet()->getHeadData(); CurrData;
       CurrData          = CurrData->getNext())
  {
    fprintf(output, "    %-*s %-*s | %*u \n",
            baseInformation.output_format.string_m,
            CurrData->getAtom1()->getIdentifier().c_str(),
            baseInformation.output_format.string_m,
            CurrData->getAtom2()->getIdentifier().c_str(),
            baseInformation.output_format.number_s, CurrData->getRange());
  }
  fflush(output);
}

void
add_structure(FILE *output,
              Structure *CurrStruc,
              BasicInformation &baseInformation,
              Flags &flags,
              StructureOptions options)
{
  if (options == StructureOptions::Initial)
    fprintf(output,
            " Initial coordinates (in PAS of inertia tensor) / 1e-10 m:\n");
  else
    fprintf(output, " Final coordinates %s/ 1e-10 m:\n",
            (flags.skipEckart ? "" : "(in Eckart frame) "));
  Coordinates *X;
  for (Atom *CurrAtom = CurrStruc->getHeadAtom(); CurrAtom;
       CurrAtom       = CurrAtom->getNext())
  {
    X = CurrAtom->getCoordinates(options);
    fprintf(output, "    %-*s %*.*f %*.*f %*.*f\n",
            baseInformation.output_format.string_l,
            CurrAtom->getIdentifier().c_str(),
            baseInformation.output_format.number_l,
            baseInformation.output_format.prec_l, X->x,
            baseInformation.output_format.number_l,
            baseInformation.output_format.prec_l, X->y,
            baseInformation.output_format.number_l,
            baseInformation.output_format.prec_l, X->z);
  }
  X = NULL;
  fprintf(output, "\n");
  fflush(output);
}

void
add_cosine_matrix(FILE *output,
                  Structure *CurrStruc,
                  BasicInformation &baseInformation,
                  Flags &flags,
                  StructureOptions options)
{
  if (options == StructureOptions::Initial)
  {
    fprintf(output, " Normalized cosine matrix of the initial structure (in "
                    "PAS of inertia tensor):\n");
    if (!flags.minimalOutput)
    {
      fprintf(output, " %s \n", std::string(header_width, '*').c_str());
      fprintf(output, " *  %-56s* \n", "Condition number = sig[1] / sig[n]");
      fprintf(output, " *  %-56s* \n",
              "                   with sig[n] largest sig[n] != 0.");
      fprintf(output, " %s \n", std::string(header_width, '*').c_str());
    }
  }
  else
    fprintf(output, " Normalized cosine matrix of the final structure%s:\n\n",
            (flags.skipEckart ? "" : " (in Eckart frame)"));

  fprintf(output, "    %-*s %-*s |   %s%s%s%s%s\n",
          baseInformation.output_format.string_m, "Nuc 1",
          baseInformation.output_format.string_m, "Nuc 2",
          centerCstring("zz", baseInformation.output_format.number_xl),
          centerCstring("xx-yy", baseInformation.output_format.number_xl),
          centerCstring("xy", baseInformation.output_format.number_xl),
          centerCstring("xz", baseInformation.output_format.number_xl),
          centerCstring("yz", baseInformation.output_format.number_xl));

  Eigen::MatrixXd C = CurrStruc->getHeadYmatrix()->determineBmatrix(
      *(CurrStruc->getParent()), CurrStruc, options);
  Eigen::MatrixXd sing = getSingularValues(C);
  unsigned int range, base;
  range = base    = 0;
  const double cn = getConditionNumber(C, range, base);
  for (SphericalHarmonics *Y = CurrStruc->getHeadYmatrix()->getHeadHarmonic();
       Y; Y                  = Y->getNext())
  {
    fprintf(output, "    %-*s %-*s | ", baseInformation.output_format.string_m,
            Y->getAtom1()->getIdentifier().c_str(),
            baseInformation.output_format.string_m,
            Y->getAtom2()->getIdentifier().c_str());
    for (int i = 0; i < COSINE_ELEMENTS_; ++i)
      fprintf(output, "%*.*e", baseInformation.output_format.number_xl,
              baseInformation.output_format.prec_l,
              C(Y->getInputIndex() - 1, i));
    fprintf(output, "\n");
  }
  fprintf(output, "\n");
  fprintf(output, " %-*s   ", baseInformation.output_format.string_xl,
          "Singular values:");
  for (int i = 0; i < COSINE_ELEMENTS_; ++i)
    fprintf(output, "%*.*e", baseInformation.output_format.number_xl,
            baseInformation.output_format.prec_l, sing(i));
  fprintf(output, "\n");
  fprintf(output, " %-*s    %*.*f\n", baseInformation.output_format.string_xl,
          "Condition number:", baseInformation.output_format.number_m,
          baseInformation.output_format.prec_m, cn);
  fprintf(output, " %-*s %*d\n", baseInformation.output_format.string_xl,
          "Rank:", baseInformation.output_format.number_m, range);
  fflush(output);
}

void
add_vector_sampling(FILE *output,
                    Structure *CurrStruc,
                    BasicInformation &baseInformation,
                    Flags &flags,
                    StructureOptions options)
{
  double r, a;
  CurrStruc->get_rdc_vector_sampling(r, a, options);
  fprintf(output, "\n Completeness of vector sampling: %s\n",
          Paper::cite(Paper::Lakomek));
  if (!flags.minimalOutput)
  {
    fprintf(output, " %s \n", std::string(header_width, '*').c_str());
    fprintf(output, " *  %-56s* \n",
            "M(Vec) = R^t * R       with R: normalized RDC vectors.");
    fprintf(output, " *  %-56s* \n",
            "r =  sqrt( m[1]/m[3] ) with m: eigenvalues of M.");
    fprintf(output, " *  %-56s* \n",
            "a = (m[2]-m[3])/m[1]    and m[1] >= m[2] >= m[3]");
    fprintf(output, " %s \n", std::string(header_width, '*').c_str());
  }
  fprintf(output, " r = %*.*f \n", baseInformation.output_format.number_l,
          baseInformation.output_format.prec_l, r);
  fprintf(output, " a = %*.*f \n", baseInformation.output_format.number_l,
          baseInformation.output_format.prec_l, a);
  fflush(output);
}

void
add_SECONDA(FILE *output,
            Molecule &CurrMol,
            BasicInformation &baseInformation,
            Flags &flags,
            StructureOptions options)
{
  int srank = 0;
  fprintf(output, " SECONDA analysis: %s\n", Paper::cite(Paper::Hus1));
  if (!flags.minimalOutput)
  {
    fprintf(output, " %s \n", std::string(header_width, '*').c_str());
    fprintf(output, " *  %-14s%-.3e%-15s* \n",
            "RDCs were normalized to Dmax = ", baseInformation.normFactRDC,
            " Hz");
    fprintf(output, " *  %-56s* \n", "Covariance weighting:");
    fprintf(output, " *  %-56s* \n", "   w(M) = 1 / SUM_i ( D[i] - <D>(M) )^2");
    fprintf(output, " *  %-14s%.1e%-35s* \n", "Eigenvalues < ",
            baseInformation.limits.zero_cutoff, " were removed.");
    fprintf(output, " *  %-34s%s.%*s* \n", "Kappa: Collectivity as defined in ",
            Paper::cite(Paper::Brueschweiler),
            static_cast<int>(21 - strlen(Paper::cite(Paper::Brueschweiler))),
            " ");
    fprintf(output, " *  %-56s* \n", "Cumulative variance of eigenvalue:");
    fprintf(output, " *  %-56s* \n",
            "   CSUM = SUM(EigVal[j])   with n | j <= n.");
    fprintf(output, " *  %-56s* \n", "rho     = EigVal[5] / EigVal[6]");
    fprintf(output, " *  %-56s* \n",
            "rho(av) = 0.2 * SUM (EigVal[i]) / EigVal[6]");
    fprintf(output, " *  %-56s* \n", "                    with i=[1,5].");
    fprintf(output, " *  %-56s* \n",
            "Estimated linear dependency of cov. matrix:");
    fprintf(output, " *  %-56s* \n", "   CDep = EigVal[1] / EigVal[5]");
    fprintf(output, " *  %-56s* \n", "Cumulative sum of heterogeneous modes:");
    fprintf(output, " *  %-56s* \n",
            "   a^2 = SUM(EigVal[q]*|Q[i]>) with q=[6,n]");
    fprintf(output, " *  %-56s* \n", "   (Q: eigenvector with index q)");
    fprintf(output, " %s \n", std::string(header_width, '*').c_str());
  }
  SECONDA_Output out = (options == StructureOptions::Optimized ?
                            baseInformation.SECONDA_output_final :
                            baseInformation.SECONDA_output_initial);

  fprintf(output, " %*s | %*s %*s %*s\n",
          baseInformation.output_format.string_l, "Index",
          baseInformation.output_format.string_l, "Eigenvalues",
          baseInformation.output_format.string_l, "Kappa",
          baseInformation.output_format.string_l, "CSUM / %");
  for (unsigned int i = 0;
       i <= out.covar_rank && i < out.Covariance_Evals.size(); ++i)
  {
    if (out.Covariance_Evals(i))
    {
      fprintf(output, " %*u | %*.*e %*.*f %*.*f \n",
              baseInformation.output_format.number_l, (i + 1),
              baseInformation.output_format.number_l,
              baseInformation.output_format.prec_m, out.Covariance_Evals(i),
              baseInformation.output_format.number_l,
              baseInformation.output_format.prec_m, out.kappa_q(i),
              baseInformation.output_format.number_l,
              baseInformation.output_format.prec_m, out.cumulative_variance(i));
    }
    else
    {
      fprintf(
          output, "%*u -%*ld | <%*.*e %s%*s %*.*f \n",
          baseInformation.output_format.number_m, (i + 1),
          baseInformation.output_format.number_s, out.Covariance_Evals.size(),
          baseInformation.output_format.number_l - 1,
          baseInformation.output_format.prec_m,
          baseInformation.limits.zero_cutoff,
          std::string(baseInformation.output_format.number_l / 2, ' ').c_str(),
          baseInformation.output_format.number_l / 2,
          centerCstring("-", baseInformation.output_format.number_l / 2),
          baseInformation.output_format.number_l,
          baseInformation.output_format.prec_m, out.cumulative_variance(i));
    }
  }

  fprintf(output, "\n %s\n", "Additional Parameter:");
  fprintf(output, " %*s | %*s%s \n", baseInformation.output_format.string_l,
          "rho", baseInformation.output_format.string_l,
          (out.Covariance_Evals(5) == .0 ?
               "undefined" :
               to_cstring(out.rho_5_6, baseInformation.output_format.number_l,
                          baseInformation.output_format.prec_m, true)),
          (out.Covariance_Evals(5) == .0 ?
               "  | (rho is not defined since EigVal[6]" :
               ""));
  fprintf(
      output, " %*s | %*s%s \n", baseInformation.output_format.string_l,
      "rho(av)", baseInformation.output_format.string_l,
      (out.Covariance_Evals(5) == .0 ?
           "undefined" :
           to_cstring(out.rho_mean_6, baseInformation.output_format.number_l,
                      baseInformation.output_format.prec_m, true)),
      (out.Covariance_Evals(5) == .0 ? "  |  is lower than the zero cutoff)" :
                                       ""));
  fprintf(output, " %*s | %*.*e\n", baseInformation.output_format.string_l,
          "CDep", baseInformation.output_format.number_l,
          baseInformation.output_format.prec_m, out.covar_condition_number);
  fprintf(output, " %*s |  %*d\n", baseInformation.output_format.number_l,
          "rank", baseInformation.output_format.number_s, out.covar_rank);
  fprintf(output, "\n");
  fprintf(output, "SECONDA cumulative sum of heterogeneous modes: %s\n",
          Paper::cite(Paper::Hus2));
  fprintf(output, "    %-*s %-*s | %s\n",
          baseInformation.output_format.string_m, "Nuc 1",
          baseInformation.output_format.string_m, "Nuc 2",
          centerCstring("a^2", baseInformation.output_format.string_l));
  RDCdata *CurrRDC = CurrMol.getHeadSet()->getHeadData();
  for (int i = 0; i < out.heterogeneity.rows(); ++i)
  {
    fprintf(output, "    %-*s %-*s | %*.*e \n",
            baseInformation.output_format.string_m,
            CurrRDC->getAtom1()->getIdentifier().c_str(),
            baseInformation.output_format.string_m,
            CurrRDC->getAtom2()->getIdentifier().c_str(),
            baseInformation.output_format.number_l,
            baseInformation.output_format.prec_l, out.heterogeneity(i));
    CurrRDC = CurrRDC->getNext();
  }

  fprintf(output, "\n Jackknife resampling of SECONDA rho: \n");
  if (flags.minimalOutput)
  {
    fprintf(output, " %s \n", std::string(header_width, '*').c_str());
    fprintf(output, " *  %-56s* \n",
            "SECONDA rho when individual RDC sets were removed");
    fprintf(output, " %s \n", std::string(header_width, '*').c_str());
  }
  fprintf(output, " %-*s | %*s %*s %*s %*s\n",
          baseInformation.output_format.string_l, "removed set",
          baseInformation.output_format.string_l, "Cdep[1/n]",
          baseInformation.output_format.string_l, "rho[5/6]",
          baseInformation.output_format.string_m, "n",
          baseInformation.output_format.string_m, "rank");

  for (RDCset *CurrSet = CurrMol.getHeadSet(); CurrSet;
       CurrSet         = CurrSet->getNext())
  {
    srank = CurrSet->get_SECONDA_sensitivity_rank();
    fprintf(
        output, " %-*s | %*.*e %*s %*d %*d",
        baseInformation.output_format.string_l,
        CurrSet->getLabel(baseInformation.output_format.string_l, true).c_str(),
        baseInformation.output_format.number_l,
        baseInformation.output_format.prec_l,
        CurrSet->get_SECONDA_gap_1_5_sensitivity(),
        baseInformation.output_format.number_l,
        (srank < 6 ? "undefined" :
                     to_cstring(CurrSet->get_SECONDA_gap_5_6_sensitivity(),
                                baseInformation.output_format.number_l,
                                baseInformation.output_format.prec_l, true)),
        baseInformation.output_format.string_m, (srank == 5 ? 5 : srank),
        baseInformation.output_format.string_m, srank);
    if (CurrSet->get_SECONDA_sensitivity_rank() !=
        static_cast<int>(out.covar_rank))
      fprintf(output, " (Warning: initial rank [%d] changed during resampling)",
              out.covar_rank);
    fprintf(output, " \n");
  }
  fflush(output);
}

void
add_Tolman_analysis(FILE *output,
                    BasicInformation &baseInformation,
                    Flags &flags,
                    StructureOptions options)
{
  SECONDA_Output out = (options == StructureOptions::Optimized ?
                            baseInformation.SECONDA_output_final :
                            baseInformation.SECONDA_output_initial);

  fprintf(output, "\n SVD of normalized RDC matrix: %s \n",
          Paper::cite(Paper::Tolman));
  if (!flags.minimalOutput)
  {
    fprintf(output, " %s \n", std::string(header_width, '*').c_str());
    fprintf(output, " *  %-56s* \n", "Normalization see SECONDA analysis.");
    fprintf(output, " %s \n", std::string(header_width, '*').c_str());
  }
  fprintf(output, " Singular values: \n");
  for (unsigned int i = 0;
       i <= out.Tolman_rank && i < out.Tolman_singular_values.size(); ++i)
  {
    if (out.Tolman_singular_values(i) > baseInformation.limits.zero_cutoff)
    {
      fprintf(output, " %*u | %*.*e \n", baseInformation.output_format.string_l,
              (i + 1), baseInformation.output_format.number_l,
              baseInformation.output_format.prec_l,
              out.Tolman_singular_values(i));
    }
    else
    {
      fprintf(output, "%*u -%*ld | <%*.*e \n",
              baseInformation.output_format.number_m, (i + 1),
              baseInformation.output_format.number_s,
              out.Tolman_singular_values.size(),
              baseInformation.output_format.number_l - 1,
              baseInformation.output_format.prec_l,
              baseInformation.limits.zero_cutoff);
    }
  }
  fprintf(output, "\n %s\n", "Addional Parameter:");
  fprintf(output, " %*s | %*.*e (used %-2u/%2u gap) \n",
          baseInformation.output_format.string_l, "Con. Num.",
          baseInformation.output_format.number_l,
          baseInformation.output_format.prec_l, out.Tolman_condition_number, 1,
          out.Tolman_CN_rank);
  fprintf(output, " %*s | %*d\n", baseInformation.output_format.string_l,
          "rank", baseInformation.output_format.number_s, out.Tolman_rank);
  fflush(output);
}

void
add_pairwise_analysis(FILE *output,
                      Molecule &CurrMol,
                      BasicInformation &baseInformation,
                      Flags &flags)
{
  int id_1, id_2;
  double sig_1, sig_2;
  RDCset *set_1, *set_2;

  if (!flags.minimalOutput)
  {
    fprintf(output, " %s \n", std::string(header_width, '*').c_str());
    fprintf(output, " *  %-56s* \n",
            "All data were generated by fitting a matrix containing");
    fprintf(output, " *  %-56s* \n", "the two respective RDC sets.");
    fprintf(output, " *  %-56s* \n",
            "R^2: Pearson correlation coefficient of linear fit.");
    fprintf(output, " %s \n", std::string(header_width, '*').c_str());
  }
  fprintf(output, " %*s %*s | %*s %*s %*s %*s   %s \n",
          baseInformation.output_format.string_l + 4, "set 1 [index]",
          baseInformation.output_format.string_l + 4, "set 2 [index]", // |
          baseInformation.output_format.number_m, "slope",
          baseInformation.output_format.number_m, "interc.",
          baseInformation.output_format.number_m, "R^2",
          2 * (baseInformation.output_format.number_m) + 2, "singular values",
          "condition number");
  for (set_1 = CurrMol.getHeadSet(); set_1->getNext(); set_1 = set_1->getNext())
  {
    id_1 = set_1->getIndex() - 1;
    for (set_2 = set_1->getNext(); set_2; set_2 = set_2->getNext())
    {
      id_2  = set_2->getIndex() - 1;
      sig_1 = baseInformation.SECONDA_output_initial.twoSetSingularValues(id_1,
                                                                          id_2);
      sig_2 = baseInformation.SECONDA_output_initial.twoSetSingularValues(id_2,
                                                                          id_1);
      fprintf(output, " %-*s[%2u] %-*s[%2u] |",
              baseInformation.output_format.string_l,
              set_1->getLabel(12, true).c_str(), set_1->getIndex(),
              baseInformation.output_format.string_l,
              set_2->getLabel(12, true).c_str(), set_2->getIndex());
      fprintf(output, " %*.*f %*.*f %*.*f ",
              baseInformation.output_format.number_m,
              baseInformation.output_format.prec_m,
              baseInformation.SECONDA_output_initial.RDC_m(id_1, id_2 - 1),
              baseInformation.output_format.number_m,
              baseInformation.output_format.prec_m,
              baseInformation.SECONDA_output_initial.RDC_b(id_1, id_2 - 1),
              baseInformation.output_format.number_m,
              baseInformation.output_format.prec_m,
              baseInformation.SECONDA_output_initial.RDC_R_2(id_1, id_2 - 1));
      fprintf(output, " %*.*f /%*.*f %*s \n",
              baseInformation.output_format.number_m,
              baseInformation.output_format.prec_m, sig_1,
              baseInformation.output_format.number_m,
              baseInformation.output_format.prec_m, sig_2,
              baseInformation.output_format.number_m + 3,
              ((sig_2) < baseInformation.limits.zero_cutoff ?
                   "undefined" :
                   to_cstring((sig_1 / sig_2),
                              baseInformation.output_format.number_m + 3,
                              baseInformation.output_format.prec_m)));
    }
  }
  set_1 = set_2 = NULL;
  fflush(output);
}

void
add_Qfactors(FILE *output, Molecule &CurrMol, BasicInformation &baseInformation)
{
  unsigned int NOS    = CurrMol.getNORsets();
  RDCset *CurrSet     = CurrMol.getHeadSet();
  unsigned int sblock = 0;
  unsigned int set;
  bool split = (NOS > sub_block_size);
  fprintf(output,
          " Q-factors of individual RDC sets during optimization: %s\n\n",
          Paper::cite(Paper::Cornilescu));
  do
  {
    if (split)
      fprintf(
          output, "%s Block of RDC sets %*d - %*d \n", (sblock ? "\n" : ""),
          baseInformation.output_format.number_s, (sblock + 1),
          baseInformation.output_format.number_s,
          ((sblock + sub_block_size) > NOS ? NOS : (sblock + sub_block_size)));
    (sblock += sub_block_size);
    fprintf(output, " %*s | ", baseInformation.output_format.string_m, "Iter.");
    for (; CurrSet && CurrSet->getIndex() <= sblock;
         CurrSet = CurrSet->getNext())
      fprintf(output, "  %*s[%2u]", baseInformation.output_format.string_l,
              CurrSet->getLabel(baseInformation.output_format.string_l, true)
                  .c_str(),
              CurrSet->getIndex());
    fprintf(output, " \n");

    for (int iteration = 0; iteration <= baseInformation.numOfOptSteps;
         ++iteration)
    {
      fprintf(output, " %*d | ", baseInformation.output_format.string_m,
              iteration);
      set = sblock - sub_block_size;
      for (; set < NOS && set < sblock; ++set)
      {
        fprintf(output, "  %*.*f", baseInformation.output_format.number_xl,
                baseInformation.output_format.prec_l,
                baseInformation.q_factors[set + iteration * NOS]);
      }
      fprintf(output, "\n");
    }

  } while (sblock < NOS);
  fprintf(output, "\n");
  fflush(output);
}

void
add_chiral_volumes(FILE *output,
                   Molecule &CurrMol,
                   BasicInformation &baseInformation,
                   Flags &flags)
{
  int iteration;
  unsigned int NOC          = baseInformation.NumberOfChiralCenters;
  unsigned int center       = 0;
  unsigned int label        = 0;
  unsigned int sblock       = 0;
  unsigned int sblock_start = 0;
  double monitor[NOC][5];
  double V;

  bool split = (NOC > sub_block_size);
  char *buf  = (char *) malloc(INFO_BUF);

  const int max_deviation = 0;
  const int pre_deviation = 1;
  const int full_diff     = 2;
  const int ini_volume    = 3;
  const int final_volume  = 4;
  Atom *CurrAtom, *rmsdAtom;

  fprintf(output, " Chiral volumes of the (pro-)chiral centers during "
                  "optimization / (1e-10 m)^3\n");
  if (!flags.minimalOutput)
  {
    fprintf(output, " %s \n", std::string(header_width, '*').c_str());
    fprintf(output, " *  %-56s* \n",
            "Vc = a' * ( b x c ) with v = r_i - r_center.");
    fprintf(output, " *  %-56s* \n", "Chiral centers: ");
    // Insert chiral volume list
    snprintf(buf, INFO_BUF, " %-6s |   %6s %6s %6s", "center", "a", "b", "c");
    fprintf(output, " *  %-56s* \n", buf);
    for (CurrAtom = CurrMol.getHeadStruc()->getHeadAtom(); CurrAtom;
         CurrAtom = CurrAtom->getNext())
    {
      if (!CurrAtom->hasChiralVolume())
        continue;
      snprintf(buf, INFO_BUF, " %-6s |  %s ", CurrAtom->getIdentifier().c_str(),
               CurrAtom->printChiralVolumeAtoms(6));
      fprintf(output, " *  %-56s* \n", buf);
    }
    fprintf(output, " %s \n", std::string(header_width, '*').c_str());
  }
  CurrAtom = CurrMol.getHeadStruc()->getHeadAtom();
  do
  {
    sblock_start = sblock;
    sblock += sub_block_size;
    if (split)
      fprintf(output, "%s Block of chiral centers %*d - %*d \n",
              (sblock > sub_block_size ? "\n" : ""),
              baseInformation.output_format.number_s, (sblock_start + 1),
              baseInformation.output_format.number_s,
              (sblock > NOC ? NOC : sblock));

    fprintf(output, " %*s | ", baseInformation.output_format.string_m, "Iter.");
    for (; CurrAtom && label < sblock; CurrAtom = CurrAtom->getNext())
    {
      if (CurrAtom->hasChiralVolume())
      {
        fprintf(output, " %*s[%*u]", baseInformation.output_format.string_l,
                CurrAtom->getIdentifier().c_str(),
                baseInformation.output_format.number_s, CurrAtom->getIndex());
        ++label;
      }
    }
    fprintf(output, " \n");

    // Initialize monitor values
    for (center = sblock_start; center < NOC && center < sblock; ++center)
    {
      monitor[center][max_deviation] = .0;
      monitor[center][pre_deviation] = .0;
      monitor[center][ini_volume]    = baseInformation.chiral_volumes[center];
      monitor[center][final_volume] =
          baseInformation
              .chiral_volumes[center + NOC * baseInformation.numOfOptSteps];
      monitor[center][full_diff] =
          monitor[center][final_volume] - monitor[center][ini_volume];
    }

    for (iteration = 0; iteration <= baseInformation.numOfOptSteps; ++iteration)
    {
      fprintf(output, " %*d | ", baseInformation.output_format.string_m,
              iteration);
      for (center = sblock_start; center < NOC && center < sblock; ++center)
      {
        V = baseInformation.chiral_volumes[center + iteration * NOC];
        fprintf(output, "  %*.*f", baseInformation.output_format.number_xl,
                baseInformation.output_format.prec_l, V);

        if (fabs(monitor[center][pre_deviation] - V) >
            fabs(monitor[center][max_deviation]))
        {
          monitor[center][max_deviation] = V - monitor[center][pre_deviation];
        }
        monitor[center][pre_deviation] = V;
      }
      fprintf(output, "\n");
    }
  } while (sblock < NOC);
  fprintf(output, "\n Chiral volume analysis: \n");

  if (!flags.minimalOutput)
  {
    fprintf(output, " %s \n", std::string(header_width, '*').c_str());
    fprintf(output, " *  %-56s* \n", "Chiral volume Vc of:");
    fprintf(output, " *  %-56s* \n", "  0) initial geometry");
    fprintf(output, " *  %-56s* \n", "  n) final geometry");
    fprintf(output, " *  %-56s* \n", "  i) geometry of individual iteration");
    fprintf(output, " *  %-56s* \n", "Total change in chiral volume:");
    fprintf(output, " *  %-56s* \n", "   Delta = Vc[n] - Vc[0]");
    fprintf(output, " *%s* \n", std::string(header_width - 2, ' ').c_str());
    fprintf(output, " *  %-56s* \n", "Estimated probability of an inversion:");
    fprintf(output, " *  %-56s* \n",
            "                            min(|Vc[n]|,|Vc[0]|)  ");
    fprintf(output, " *  %-56s* \n",
            "   Prob. = 50 * |Delta| * ------------------------ ");
    fprintf(output, " *  %-56s* \n",
            "                           max(|Vc[n]|,|Vc[0]|)^2 ");
    fprintf(output, " *%s* \n", std::string(header_width - 2, ' ').c_str());
    fprintf(output, " *  %-56s* \n", "Estimated confidence:");
    fprintf(output, " *  %-56s* \n", "                  min(|Vc[n]|,|Vc[0]|) ");
    fprintf(output, " *  %-56s* \n",
            "   Conf. = 100 * ---------------------- ");
    fprintf(output, " *  %-56s* \n", "                  max(|Vc[n]|,|Vc[0]|) ");
    fprintf(output, " *%s* \n", std::string(header_width - 2, ' ').c_str());
    fprintf(output, " *  %-56s* \n",
            "rmsd = sqrt( SUM(<|Vc|> - |Vc[i]|)^2/N )");
    fprintf(output, " %s \n", std::string(header_width, '*').c_str());
  }
  sblock = center = label = 0;
  CurrAtom                = CurrMol.getHeadStruc()->getHeadAtom();
  rmsdAtom                = CurrMol.getHeadStruc()->getHeadAtom();
  do
  {
    sblock_start = sblock;
    sblock += sub_block_size;
    if (split)
      fprintf(output, "%s Block of chiral centers %*d - %*d \n",
              (sblock > sub_block_size ? "\n" : ""),
              baseInformation.output_format.number_s, (sblock_start + 1),
              baseInformation.output_format.number_s,
              (sblock > NOC ? NOC : sblock));

    fprintf(output, " %9s%-14s | ", "Parameter", " / unit");
    for (label    = sblock_start; CurrAtom && label < sblock;
         CurrAtom = CurrAtom->getNext())
    {
      if (CurrAtom->hasChiralVolume())
      {
        fprintf(output, " %*s[%*u]", baseInformation.output_format.string_l,
                CurrAtom->getIdentifier().c_str(),
                baseInformation.output_format.number_s, CurrAtom->getIndex());
        ++label;
      }
    }
    fprintf(output, " \n");
    fprintf(output, " %9s%-14s | ", "Delta ", " ");
    for (center = sblock_start; center < NOC && center < sblock; ++center)
      fprintf(output, "  %*.*f", baseInformation.output_format.number_xl,
              baseInformation.output_format.prec_l, monitor[center][full_diff]);


    fprintf(output, "\n %9s%-14s |  ", "Prob.", " / %");
    for (center = sblock_start; center < NOC && center < sblock; ++center)
    {
      double Prob = min(fabs(monitor[center][final_volume]),
                        fabs(monitor[center][ini_volume])) *
                    fabs(monitor[center][full_diff]) /
                    pow(max(fabs(monitor[center][final_volume]),
                            fabs(monitor[center][ini_volume])),
                        2) *
                    50.0;

      snprintf(buf, INFO_BUF, " within %u steps.",
               baseInformation.MonteCarloOutput.steps);

      if (Prob > 20.0)
        snprintf(buf, INFO_BUF, "[%.*f]", baseInformation.output_format.prec_l,
                 Prob);
      else
        snprintf(buf, INFO_BUF, "%.*f ", baseInformation.output_format.prec_l,
                 Prob);
      fprintf(output, "  %*s", baseInformation.output_format.number_xl, buf);
    }
    fprintf(output, "\n %9s%-14s | ", "Conf.", " / %");
    for (center = sblock_start; center < NOC && center < sblock; ++center)
    {
      double Conf = 100 *
                    min(fabs(monitor[center][final_volume]),
                        fabs(monitor[center][ini_volume])) /
                    max(fabs(monitor[center][final_volume]),
                        fabs(monitor[center][ini_volume]));
      fprintf(output, "  %*.*f", baseInformation.output_format.number_xl,
              baseInformation.output_format.prec_l, Conf);
    }

    fprintf(output, "\n %9s%-14s | ", "rmsd ", " / (1e-10 m)^3");
    for (label    = sblock_start; rmsdAtom && label < sblock;
         rmsdAtom = rmsdAtom->getNext())
    {
      if (!rmsdAtom->hasChiralVolume())
        continue;
      fprintf(
          output, "  %*.*f", baseInformation.output_format.number_xl,
          baseInformation.output_format.prec_l,
          rmsdAtom->check_Chiral_Volume_validity(StructureOptions::Optimized));
      ++label;
    }
    fprintf(output, "\n");
  } while (sblock < NOC);
  fflush(output);
  free(buf);
}

void
add_TITANIA_angles(FILE *output,
                   Molecule &CurrMol,
                   Structure *CurrStruc,
                   BasicInformation &baseInformation)
{
  SphericalHarmonics *CurrY = CurrStruc->getHeadYmatrix()->getHeadHarmonic();
  SphericalHarmonics *IniY =
      CurrMol.getHeadStruc()->getHeadYmatrix()->getHeadHarmonic();
  double thetaI, phiI, thetaO, phiO;
  thetaI = phiI = thetaO = phiO = .0;


  fprintf(output, " Polar angles of RDC vectors / deg: %s%s%s\n",
          Paper::cite_multiple_start(Paper::Meiler),
          Paper::cite_multiple_middle(Paper::Peti),
          Paper::cite_multiple_end(Paper::Lakomek));

  fprintf(output, "\n   %*s  | ",
          2 * baseInformation.output_format.string_m + 1, "Nuclei list");
  fprintf(output, " %*s  | ", 2 * baseInformation.output_format.string_l + 1,
          "Initial values  ");
  fprintf(output, " %*s  | ", 2 * baseInformation.output_format.string_l + 1,
          "Optimize values ");
  fprintf(output, " %*s \n", baseInformation.output_format.string_m, "Change");

  fprintf(output, "   %-*s %-*s  | ", baseInformation.output_format.string_m,
          "Nuc 1", baseInformation.output_format.string_m, "Nuc 2");
  fprintf(output, "  %*s %*s | ", baseInformation.output_format.string_l,
          "theta", baseInformation.output_format.string_l, "phi");
  fprintf(output, "  %*s %*s | ", baseInformation.output_format.string_l,
          "theta", baseInformation.output_format.string_l, "phi");
  fprintf(output, " %*s \n", baseInformation.output_format.number_m, "alpha");

  for (; IniY && CurrY; CurrY = CurrY->getNext(), IniY = IniY->getNext())
  {
    thetaO = CurrY->getTheta(StructureOptions::Optimized);
    phiO   = CurrY->getPhi(StructureOptions::Optimized);
    thetaI = IniY->getTheta(StructureOptions::Initial);
    phiI   = IniY->getPhi(StructureOptions::Initial);
    A2goodA(thetaO, phiO);
    A2goodA(thetaI, phiI);

    fprintf(output, "    %-*s %-*s | ", baseInformation.output_format.string_m,
            CurrY->getAtom1()->getIdentifier().c_str(),
            baseInformation.output_format.string_m,
            CurrY->getAtom2()->getIdentifier().c_str());
    fprintf(output, "  %*.*f %*.*f | ", baseInformation.output_format.number_l,
            baseInformation.output_format.prec_m, rad2deg(thetaI),
            baseInformation.output_format.number_l,
            baseInformation.output_format.prec_m, rad2deg(phiI));
    fprintf(output, "  %*.*f %*.*f | ", baseInformation.output_format.number_l,
            baseInformation.output_format.prec_m, rad2deg(thetaO),
            baseInformation.output_format.number_l,
            baseInformation.output_format.prec_m, rad2deg(phiO));
    fprintf(output, " %*.*f \n", baseInformation.output_format.number_m,
            baseInformation.output_format.prec_m,
            rad2deg(getPolar2Beta(thetaI, phiI, thetaO, phiO)));
  }

  fprintf(output, "\n");
  fflush(output);
}

void
add_MonteCarlo_angles(FILE *output,
                      Structure *CurrStruc,
                      BasicInformation &baseInformation,
                      Flags &flags,
                      MC_Output mcout)
{
  SphericalHarmonics *CurrY = CurrStruc->getHeadYmatrix()->getHeadHarmonic();

  char *buf = (char *) malloc(INFO_BUF);
  double thetaO, phiO, MCtheta, MCphi, dMCtheta, dMCphi, CyVar, CySigma, R;
  thetaO = phiO = MCtheta = MCphi = dMCtheta = dMCphi = CyVar = CySigma = R =
      .0;

  fprintf(output, " Monte-Carlo Bootstrap of final iteration: %s%s \n",
          Paper::cite_multiple_start(Paper::Losonczi),
          Paper::cite_multiple_end(Paper::Mardia));
  if (!flags.minimalOutput && output != baseInformation.trjFile)
  {
    fprintf(output, " %s \n", std::string(header_width, '*').c_str());
    snprintf(buf, INFO_BUF, " within %u steps.",
             baseInformation.MonteCarloOutput.steps);
    fprintf(output, " *  %-56s* \n",
            "Monte-Carlo Bootstrap results were collected");
    fprintf(output, " *  %-56s* \n", buf);
    fprintf(output, " *  %-56s* \n", "Mean length: ");
    fprintf(output, " *  %-56s* \n", "    R = sqrt( X^2 + Y^2 + Z^2 )");
    fprintf(output, " *  %-56s* \n",
            "    (X,Y,Z) as components of mean directional vectors");
    fprintf(output, " *  %-56s* \n", "circular variance:");
    fprintf(output, " *  %-56s* \n", "   cvar = 1 - R");
    fprintf(output, " *  %-56s* \n", "circular standard deviation:");
    fprintf(output, " *  %-56s* \n", "   csig = sqrt( -2 log (R) )");
    fprintf(output, " %s \n", std::string(header_width, '*').c_str());
  }

  // Header line
  fprintf(output, "   %*s  | ", 2 * baseInformation.output_format.string_m + 1,
          "Nuclei list");
  fprintf(output, " %*s | ", 4 * baseInformation.output_format.number_m + 10,
          centerCstring("Monte-Carlo angles / deg",
                        4 * baseInformation.output_format.number_m + 10));
  fprintf(output, " %*s | ",
          2 * baseInformation.output_format.number_m +
              baseInformation.output_format.string_l + 3,
          centerCstring("Directional statistics of MC",
                        2 * baseInformation.output_format.number_m +
                            baseInformation.output_format.string_l + 3));
  fprintf(output, " %*s \n", 2 * baseInformation.output_format.number_m + 4,
          centerCstring("Iterative vs. MC",
                        2 * baseInformation.output_format.number_m + 4));

  // Description
  fprintf(output, "   %-*s %-*s  | ", baseInformation.output_format.string_m,
          "Nuc 1", baseInformation.output_format.string_m, "Nuc 2");
  fprintf(output, " %*s +- %*s  %*s +- %*s | ",
          baseInformation.output_format.number_m, "theta",
          baseInformation.output_format.number_m, "err",
          baseInformation.output_format.number_m, "phi",
          baseInformation.output_format.number_m, "err");
  fprintf(output, " %*s  %*s %*s | ", baseInformation.output_format.number_m,
          centerCstring("R", baseInformation.output_format.string_m),
          baseInformation.output_format.number_m, "cvar",
          baseInformation.output_format.string_l, "csig/deg");
  fprintf(output, " %*s +- %*s \n", baseInformation.output_format.number_m,
          "alpha", baseInformation.output_format.number_m, "err");

  for (; CurrY; CurrY = CurrY->getNext())
  {
    MCtheta  = mcout.p_mean(0, CurrY->getInputIndex() - 1);
    MCphi    = mcout.p_mean(1, CurrY->getInputIndex() - 1);
    dMCtheta = mcout.p_sigm(0, CurrY->getInputIndex() - 1);
    dMCphi   = mcout.p_sigm(1, CurrY->getInputIndex() - 1);
    thetaO   = CurrY->getTheta(StructureOptions::Optimized);
    phiO     = CurrY->getPhi(StructureOptions::Optimized);
    A2goodA(MCtheta, MCphi);
    R       = mcout.R_mean(CurrY->getInputIndex() - 1);
    CyVar   = 1.0 - R;
    CySigma = sqrt(-2.0 * log(R));
    fprintf(output, "    %-*s %-*s | ", baseInformation.output_format.string_m,
            CurrY->getAtom1()->getIdentifier().c_str(),
            baseInformation.output_format.string_m,
            CurrY->getAtom2()->getIdentifier().c_str());
    fprintf(output, " %*.*f +- %*.*f  %*.*f +- %*.*f | ",
            baseInformation.output_format.number_m,
            baseInformation.output_format.prec_m, rad2deg(MCtheta),
            baseInformation.output_format.number_m,
            baseInformation.output_format.prec_m, rad2deg(dMCtheta),
            baseInformation.output_format.number_m,
            baseInformation.output_format.prec_m, rad2deg(MCphi),
            baseInformation.output_format.number_m,
            baseInformation.output_format.prec_m, rad2deg(dMCphi));
    fprintf(output, " %*.*f  %*.*f %*.*f | ",
            baseInformation.output_format.number_m,
            baseInformation.output_format.prec_l, R,
            baseInformation.output_format.number_m,
            baseInformation.output_format.prec_l, CyVar,
            baseInformation.output_format.number_l,
            baseInformation.output_format.prec_l, rad2deg(CySigma));
    fprintf(output, " %*.*f +- %*.*f \n",
            baseInformation.output_format.number_m,
            baseInformation.output_format.prec_m,
            rad2deg(getPolar2Beta(thetaO, phiO, MCtheta, MCphi)),
            baseInformation.output_format.number_m,
            baseInformation.output_format.prec_m,
            rad2deg(getDpolar2Dbeta(thetaO, .0, phiO, .0, MCtheta, dMCtheta,
                                    MCphi, dMCphi)));
  }

  fprintf(output, "\n");
  fflush(output);
  free(buf);
}

void
add_Srdc(FILE *output,
         Structure *CurrStruc,
         BasicInformation &baseInformation,
         Flags &flags,
         Eigen::MatrixXd &D,
         Eigen::MatrixXd &Dcalc)
{
  SphericalHarmonics *CurrY;
  int i = 0;
  Eigen::MatrixXd D_i, Dcalc_i, sig;
  fprintf(output, " %s: \n", "Motional analysis of RDC vectors");
  if (!flags.minimalOutput && output != baseInformation.trjFile)
  {
    fprintf(output, " %s \n", std::string(header_width, '*').c_str());
    fprintf(output, " *  %-56s* \n",
            "The used spherical harmonics Y are expressed in the");
    fprintf(output, " *  %-56s* \n", " respective vector frame.");
    if (flags.scaleWithSoverall)
      fprintf(output, " *  %-56s* \n",
              "The spherical harmonics Y were normalized by S(overall).");
    fprintf(output, " *  %-56s* \n", "S^2(RDC) = 4pi/5 * SUM_m <Y2,m><Y2,m*>");
    fprintf(output, " *  %-56s* \n", "S^2(ax)  = 4pi/5         <Y2,0><Y2,0*>");
    fprintf(output, " *  %-56s* \n", " ");
    fprintf(output, " *  %-56s* \n", "            SUM_i <Yl,i><Yl,i*>");
    fprintf(output, " *  %-56s* \n", "eta(RDC) = ---------------------");
    fprintf(output, " *  %-56s* \n", "            SUM_j <Yl,j><Yl,j*>");
    fprintf(output, " *  %-56s* \n", " ");
    fprintf(output, " *  %-56s* \n", "           1       Im( <Y2,2> )");
    fprintf(output, " *  %-56s* \n", "phi(RDC) = - atan --------------");
    fprintf(output, " *  %-56s* \n", "           2       Re( <Y2,2> )");
    fprintf(output, " *  %-56s* \n", " ");
    fprintf(output, " *  %-56s* \n", "             SUM_i ( <x> - x[i] )^2");
    fprintf(output, " *  %-56s* \n", "Chi^2    =  ------------------------");
    fprintf(output, " *  %-56s* \n", "             SUM_i delta( x[i] )^2");
    fprintf(output, " *  %-56s* \n",
            " where delta(x[i]) is the RDC error estimated");
    fprintf(output, " *  %-56s* \n", " by the user in the input.");
    fprintf(output, " *  %-56s* \n", "S(overall) = sqrt( 1 / S^2(RDC,max) ) ");
    fprintf(output, " %s \n", std::string(header_width, '*').c_str());
  }
  fprintf(output, "\n %12s = %*.*f %s%s \n\n", "S(overall)",
          baseInformation.output_format.number_l,
          baseInformation.output_format.prec_l,
          CurrStruc->getHeadYmatrix()->getSoverall(),
          Paper::cite_multiple_start(Paper::Meiler),
          Paper::cite_multiple_end(Paper::Peti));
  fprintf(output, " %-*s %-*s %*s %*s %*s %*s %*s %*s\n",
          baseInformation.output_format.string_m, "Nuc 1",
          baseInformation.output_format.string_m, "Nuc 2",
          baseInformation.output_format.string_l, "S^2(RDC)",
          baseInformation.output_format.string_l, "S^2(ax)",
          baseInformation.output_format.string_l, "eta(RDC)",
          baseInformation.output_format.string_l, "phi(aniso)",
          baseInformation.output_format.string_l, "Chi^2",
          baseInformation.output_format.string_l, "Dmax / Hz");
  sig = Eigen::MatrixXd::Ones(1, D.cols());
  for (CurrY = CurrStruc->getHeadYmatrix()->getHeadHarmonic(); CurrY;
       CurrY = CurrY->getNext(), ++i)
  {
    Dcalc_i = Dcalc.row(i);
    D_i     = D.row(i);
    fprintf(
        output, " %-*s %-*s %*.*f %*.*f %*.*f %*.*f %*.*f %*.*f \n",
        baseInformation.output_format.string_m,
        CurrY->getAtom1()->getIdentifier().c_str(),
        baseInformation.output_format.string_m,
        CurrY->getAtom2()->getIdentifier().c_str(),
        baseInformation.output_format.number_l,
        baseInformation.output_format.prec_l, CurrY->getS2rdc(),
        baseInformation.output_format.number_l,
        baseInformation.output_format.prec_l, CurrY->getSaxial(),
        baseInformation.output_format.number_l,
        baseInformation.output_format.prec_l, CurrY->getetardc(),
        baseInformation.output_format.number_l,
        baseInformation.output_format.prec_l, rad2deg(CurrY->getAnisoTheta()),
        baseInformation.output_format.number_l,
        baseInformation.output_format.prec_l,
        Chi_square(D_i, Dcalc_i, CurrY->get_sigma_square()),
        baseInformation.output_format.number_l,
        baseInformation.output_format.prec_s,
        (CurrY->getRDC()->getKappa() /
         pow(CurrY->getRDC()->getDistance(StructureOptions::Optimized), 3)));
  }
  fprintf(output, "\n");
  fflush(output);
}

void
add_stop_criteria(FILE *output, BasicInformation &baseInformation, Flags &flags)
{
  fprintf(output, " Stop criteria:\n");
  char *buf = (char *) malloc(INFO_BUF);

  if (!flags.minimalOutput && output != baseInformation.trjFile)
  {
    fprintf(output, " %s \n", std::string(header_width, '*').c_str());
    fprintf(output, " *  %-56s* \n", "Calculation of stop criteria:");
    fprintf(output, " *  %-56s* \n", "                 X[i,j] - X[i-1,j]");
    fprintf(output, " *  %-56s* \n", " rmsd(x) = sqrt -------------------");
    fprintf(output, " *  %-56s* \n", "                        N");
    fprintf(output, " *  %-56s* \n", " where:");
    if (flags.monteCarloBootstrapping)
    {
      fprintf(output, " *  %-56s* \n",
              "  - X are individual parameter (or their standard");
      fprintf(output, " *  %-56s* \n", "     deviation sigm).");
    }
    else
      fprintf(output, " *  %-56s* \n", "  - X are individual parameter.");
    fprintf(output, " *  %-56s* \n", "  - i are the respective iteration.");
    fprintf(output, " *  %-56s* \n",
            "  - j is the index of the respective element.");
    if (flags.monteCarloBootstrapping)
    {
      fprintf(output, " *  %-56s* \n",
              "The order tensor S and the spherical coordinates and");
      fprintf(output, " *  %-56s* \n",
              " their standard deviations were obtained by the MC");
      fprintf(output, " *  %-56s* \n", " bootstrap.");
      fprintf(output, " *  %-56s* \n",
              "The mean orientation vector length R was obtained by");
      fprintf(output, " *  %-56s* \n", " the MC bootstrap.");
    }
    fprintf(output, " *  %-56s* \n",
            "Q was obtained by back calculating the RDC with the");
    fprintf(output, " *  %-56s* \n", " optimized coordinates.");
  }
  fprintf(output, " %s \n", std::string(header_width, '*').c_str());
  snprintf(buf, INFO_BUF, "%-16s[%*s]: %*s %s", "parameter",
           baseInformation.output_format.string_l, "limit",
           baseInformation.output_format.string_l, "value", "converged");
  fprintf(output, " *  %-56s* \n", buf);
  if (flags.monteCarloBootstrapping)
  {
    fprintf(output, " *%58s* \n", std::string(header_width - 2, '-').c_str());

    snprintf(buf, INFO_BUF, "%-16s[%*.*e]: %*.*e  %6s", "rmsd(S)",
             baseInformation.output_format.number_l,
             baseInformation.output_format.prec_m,
             baseInformation.limits.alignment_mean_convergence,
             baseInformation.output_format.number_l,
             baseInformation.output_format.prec_m,
             baseInformation.stop_crit.alignment_mean_convergence,
             (baseInformation.MC_stop_reason & MC_stop::AMean ? "YES" : "NO"));
    fprintf(output, " *  %-56s* \n", buf);

    snprintf(buf, INFO_BUF, "%-16s[%*.*e]: %*.*e  %6s", "rmsd(sigm[S])",
             baseInformation.output_format.number_l,
             baseInformation.output_format.prec_m,
             baseInformation.limits.alignment_sigm_convergence,
             baseInformation.output_format.number_l,
             baseInformation.output_format.prec_m,
             baseInformation.stop_crit.alignment_sigm_convergence,
             (baseInformation.MC_stop_reason & MC_stop::ASigm ? "YES" : "NO"));
    fprintf(output, " *  %-56s* \n", buf);

    snprintf(buf, INFO_BUF, "%-16s[%*.*e]: %*.*e  %6s", "rmsd(p)",
             baseInformation.output_format.number_l,
             baseInformation.output_format.prec_m,
             baseInformation.limits.sphericals_mean_convergence,
             baseInformation.output_format.number_l,
             baseInformation.output_format.prec_m,
             baseInformation.stop_crit.sphericals_mean_convergence,
             (baseInformation.MC_stop_reason & MC_stop::pMean ? "YES" : "NO"));
    fprintf(output, " *  %-56s* \n", buf);

    snprintf(buf, INFO_BUF, "%-16s[%*.*e]: %*.*e  %6s", "rmsd(sigm[p])",
             baseInformation.output_format.number_l,
             baseInformation.output_format.prec_m,
             baseInformation.limits.sphericals_sigm_convergence,
             baseInformation.output_format.number_l,
             baseInformation.output_format.prec_m,
             baseInformation.stop_crit.sphericals_sigm_convergence,
             (baseInformation.MC_stop_reason & MC_stop::pSigm ? "YES" : "NO"));
    fprintf(output, " *  %-56s* \n", buf);

    snprintf(buf, INFO_BUF, "%-16s[%*.*e]: %*.*e  %6s", "rmsd(R)",
             baseInformation.output_format.number_l,
             baseInformation.output_format.prec_m,
             baseInformation.limits.sphericals_spread_convergence,
             baseInformation.output_format.number_l,
             baseInformation.output_format.prec_m,
             baseInformation.stop_crit.sphericals_spread_convergence,
             (baseInformation.MC_stop_reason & MC_stop::RMean ? "YES" : "NO"));
    fprintf(output, " *  %-56s* \n", buf);
  }
  snprintf(buf, INFO_BUF, "%-16s[%*.*e]: %*.*e  %6s", "delta(Q)",
           baseInformation.output_format.number_l,
           baseInformation.output_format.prec_m,
           baseInformation.limits.Q_factor_convergence,
           baseInformation.output_format.number_l,
           baseInformation.output_format.prec_m,
           baseInformation.stop_crit.Q_factor_convergence,
           (baseInformation.MC_stop_reason & MC_stop::QFac ? "YES" : "NO"));
  fprintf(output, " *  %-56s* \n", buf);

  snprintf(
      buf, INFO_BUF, "%-16s[%*d]: %*s  %6s", "max iter.",
      baseInformation.output_format.number_l,
      baseInformation.limits.max_titania_iterations,
      baseInformation.output_format.number_l, " ",
      (baseInformation.MC_stop_reason & MC_stop::MaxIter ? "REACHED" : "NO"));
  fprintf(output, " *  %-56s* \n", buf);
  fprintf(output, " %s \n", std::string(header_width, '*').c_str());
  fprintf(output, "\n");
  fflush(output);
  free(buf);
}

void
add_generated_files(FILE *output,
                    BasicInformation &baseInformation,
                    Flags &flags,
                    Molecule &CurrMol)
{
  char *buf = (char *) malloc(INFO_BUF);

  fprintf(output, " Generated files:\n");
  if (!flags.minimalOutput)
  {
    fprintf(output, " %s \n", std::string(header_width, '*').c_str());
    fprintf(output, " *  %-56s* \n", "File contents:");
    fprintf(output, " *  %-56s* \n",
            " -.out: General file for input and ouput information.");
    fprintf(output, " *  %-56s* \n",
            " -.trj: Detailed file on individual iteration steps.");
    fprintf(output, " *  %-56s* \n",
            "        Most information can be found here.");
    fprintf(output, " *  %-56s* \n",
            " -.xyz: xyz conform file containing all structures.");
    fprintf(output, " *  %-56s* \n",
            "        Labels were erased to be inline with structure");
    fprintf(output, " *  %-56s* \n", "        displaying software.");
    fprintf(output, " *  %-56s* \n",
            " -.ali: Small and specfic file containing alignment");
    fprintf(output, " *  %-56s* \n",
            "        information on individual RDC sets.");
    fprintf(output, " *  %-56s* \n", " PDF files: ");
    fprintf(output, " *  %-56s* \n",
            " -.pca: SECONDA plots of initial and final structure.");
    fprintf(output, " *  %-56s* \n",
            " -.cvp: Trajectory of the chiral volume and stop ");
    fprintf(output, " *  %-56s* \n", "        criteria.");
    fprintf(output, " *  %-56s* \n", " -.dyn: Plot of S^2(rdc) and S^2(as).");
    fprintf(output, " %s ", std::string(header_width, '*').c_str());
  }
  fprintf(output, "\n");

  // TODO convert all outputs to FILE* and update the respective functions
  fprintf(output, "    - %-*s \n", baseInformation.output_format.string_xxl,
          baseInformation.outputFileName.c_str());
  if (baseInformation.trjFileName != "") //&& baseInformation.trj.good() )
    fprintf(output, "    - %-*s \n", baseInformation.output_format.string_xxl,
            baseInformation.trjFileName.c_str());
  if (baseInformation.xyzFileName != "") //&& baseInformation.xyzFile.good() )
    fprintf(output, "    - %-*s \n", baseInformation.output_format.string_xxl,
            baseInformation.xyzFileName.c_str());

  if (!flags.skipSCRM && flags.bigData)
  {
    snprintf(buf, INFO_BUF, "%s%s", baseInformation.outputFileName.c_str(),
             ".dat");
    fprintf(output, "    - %-*s \n", baseInformation.output_format.string_xxl,
            buf);
  }

  if (baseInformation.debugFileName !=
      "") // && baseInformation.debugFile.good() )
    fprintf(output, "    - %-*s \n", baseInformation.output_format.string_xxl,
            baseInformation.debugFileName.c_str());

  if (!flags.skipSCRM && flags.plotKappaQ)
  {
    snprintf(buf, INFO_BUF, "%s%s", baseInformation.outputFileName.c_str(),
             ".pca.pdf");
    fprintf(output, "    - %-*s \n", baseInformation.output_format.string_xxl,
            buf);
  }
  if (!flags.skipSCRM && baseInformation.numOfOptSteps > 1 &&
      flags.plotTrajectory)
  {
    snprintf(buf, INFO_BUF, "%s%s", baseInformation.outputFileName.c_str(),
             ".cvp.pdf");
    fprintf(output, "    - %-*s \n", baseInformation.output_format.string_xxl,
            buf);
  }
  if (!flags.skipSCRM && flags.plotRDCdynamics)
  {
    snprintf(buf, INFO_BUF, "%s%s", baseInformation.outputFileName.c_str(),
             ".dyn.pdf");
    fprintf(output, "    - %-*s \n", baseInformation.output_format.string_xxl,
            buf);
  }
  if (!flags.skipSCRM && flags.plotMonteCarlo)
  {
    snprintf(buf, INFO_BUF, "%s%s", baseInformation.outputFileName.c_str(),
             ".mcp.pdf");
    fprintf(output, "    - %-*s \n", baseInformation.output_format.string_xxl,
            buf);
  }

  RDCset *CurrSet = CurrMol.getHeadSet();
  while (!flags.skipSCRM && flags.outputAli && CurrSet)
  {
    snprintf(buf, INFO_BUF, "%s%s%s", baseInformation.outputFileName.c_str(),
             CurrSet->getLabel().c_str(), ".ali");
    fprintf(output, "    - %-*s \n", baseInformation.output_format.string_xxl,
            buf);
    CurrSet = CurrSet->getNext();
  }
  CurrSet = CurrMol.getHeadSet();
  while (flags.titania2hotfcht && CurrSet)
  {
    snprintf(buf, INFO_BUF, "%s%s%s", baseInformation.hotFCHTbase.c_str(),
             CurrSet->getLabel().c_str(), ".fcht");
    fprintf(output, "    - %-*s \n", baseInformation.output_format.string_xxl,
            buf);
    CurrSet = CurrSet->getNext();
  }
  fprintf(output, "\n");
  fflush(output);
  free(buf);
}

void
add_runtime_metrics(FILE *output,
                    BasicInformation &baseInformation,
                    Flags &flags)
{
  fprintf(output, " Runtime statistics: \n");
  fprintf(output, "\n");
  fprintf(output, "   TITANIA took   * %*.*f s for full execution...\n",
          baseInformation.output_format.number_l,
          baseInformation.output_format.prec_l, baseInformation.fulltime);
  fprintf(output,
          "                  * %*.*f s (%*.*f%%) for solving orientation "
          "equations...\n",
          baseInformation.output_format.number_l,
          baseInformation.output_format.prec_l, baseInformation.SCRMtime,
          baseInformation.output_format.number_m,
          baseInformation.output_format.prec_m,
          (baseInformation.SCRMtime / baseInformation.fulltime) * 100);

  fprintf(output,
          "                  * %*.*f s (%*.*f%%) for generating the "
          "structures...\n",
          baseInformation.output_format.number_l,
          baseInformation.output_format.prec_l, baseInformation.structureTime,
          baseInformation.output_format.number_m,
          baseInformation.output_format.prec_m,
          (baseInformation.structureTime / baseInformation.fulltime) * 100);
  fprintf(output,
          "                  * %*.*f s (%*.*f%%) for handling holonomic terms "
          "(and MMFF94)...\n",
          baseInformation.output_format.number_l,
          baseInformation.output_format.prec_l, baseInformation.MMFF94time,
          baseInformation.output_format.number_m,
          baseInformation.output_format.prec_m,
          (baseInformation.MMFF94time / baseInformation.fulltime) * 100);
  fprintf(output,
          "                  * %*.*f s (%*.*f%%) for performung Monte-Carlo "
          "bootstraps...\n",
          baseInformation.output_format.number_l,
          baseInformation.output_format.prec_l, baseInformation.MCtime,
          baseInformation.output_format.number_m,
          baseInformation.output_format.prec_m,
          (baseInformation.MCtime / baseInformation.fulltime) * 100);
  fprintf(output,
          "                  * %*.*f s (%*.*f%%) for TITANIA background "
          "operations...\n",
          baseInformation.output_format.number_l,
          baseInformation.output_format.prec_l, baseInformation.Systemtime,
          baseInformation.output_format.number_m,
          baseInformation.output_format.prec_m,
          (baseInformation.Systemtime / baseInformation.fulltime) * 100);

  struct rusage rusage;
  getrusage(RUSAGE_SELF, &rusage);
  fprintf(output, "                  * %*.*Lf MB peak memory (rss)...\n",
          baseInformation.output_format.number_l,
          baseInformation.output_format.prec_s,
          (size_t) (rusage.ru_maxrss) / 1024.0L);

  fprintf(output, "\n");
  fflush(output);
}

void
makeHeader(FILE *output,
           BasicInformation &baseInformation,
           Flags &flags,
           InputFile &reader)
{
  // START of output
  fprintf(output, " %s \n", std::string(header_width, '=').c_str());
  fprintf(output, " =%s= \n", std::string(header_width - 2, ' ').c_str());
  frame_string(output, "TITANIA", 25);
  frame_string(output, "General output (.out)", 18);
  fprintf(output, " =%s= ", std::string(header_width - 2, ' ').c_str());

  // General information part
  information_header(output, "General information");
  add_contributors(output);

  // TODO add_TITANIA_citation (output);

  // Technical information part
  information_header(output, "Technical information");
  add_compilation_information(output, baseInformation);
  add_externals(output, baseInformation, flags);

  if (!flags.minimalOutput)
  {
    // Host information part
    information_header(output, "Host information");
    add_host_information(output, baseInformation, flags);
  }

  // Input information part
  information_header(output, "Input information");
  add_input_file_information(output, baseInformation);
  echoInput(output, reader, baseInformation);

  // Run information part
  information_header(output, "Run information");
  add_run_information(output, baseInformation);
  add_flag_information(output, baseInformation, flags);
  fflush(output);
}

void
makeInputInformation(FILE *output,
                     BasicInformation &baseInformation,
                     Flags &flags,
                     Molecule &CurrMol)
{
  information_header(output, "Metrics");
  add_metrics(output, baseInformation);

  information_header(output, "RDC Lists");
  add_rdc_sets(output, CurrMol, baseInformation, flags);

  information_header(output, "Structure information");
  add_structure(output, CurrMol.getHeadStruc(), baseInformation, flags,
                StructureOptions::Initial);
  add_cosine_matrix(output, CurrMol.getHeadStruc(), baseInformation, flags,
                    StructureOptions::Initial);
  add_vector_sampling(output, CurrMol.getHeadStruc(), baseInformation, flags,
                      StructureOptions::Initial);

  information_header(output, "RDC matrix information");
  add_SECONDA(output, CurrMol, baseInformation, flags,
              StructureOptions::Initial);
  add_Tolman_analysis(output, baseInformation, flags,
                      StructureOptions::Initial);

  information_header(output, "Pairwise correlation of RDC sets");
  add_pairwise_analysis(output, CurrMol, baseInformation, flags);
  fflush(output);
}

void
makeReferences(FILE *output, BasicInformation &baseInformation)
{
  information_header(output, "Full references");
  fprintf(output, "  ** Mathematical models **\n");
  fprintf(output, " %s   %s\n", Paper::cite(Paper::Meiler),
          Paper::print(Paper::Meiler));
  fprintf(output, " %s   %s\n", Paper::cite(Paper::Peti),
          Paper::print(Paper::Peti));
  fprintf(output, "  ** Iterative desciption **\n");
  fprintf(output, " %s   %s\n", Paper::cite(Paper::Lakomek),
          Paper::print(Paper::Lakomek));
  fprintf(output, "  ** Interpretation of RDC data **\n");
  fprintf(output, " %s   %s\n", Paper::cite(Paper::Hus1),
          Paper::print(Paper::Hus1));
  fprintf(output, " %s   %s\n", Paper::cite(Paper::Brueschweiler),
          Paper::print(Paper::Brueschweiler));
  fprintf(output, " %s   %s\n", Paper::cite(Paper::Hus2),
          Paper::print(Paper::Hus2));
  fprintf(output, " %s   %s\n", Paper::cite(Paper::Tolman),
          Paper::print(Paper::Tolman));
  fprintf(output, "  ** Statistical Analysis **\n");
  fprintf(output, " %s   %s\n", Paper::cite(Paper::Losonczi),
          Paper::print(Paper::Losonczi));
  fprintf(output, " %s   %s\n", Paper::cite(Paper::Mardia),
          Paper::print(Paper::Mardia));
  fprintf(output, "  ** Redundant internal coordinates **\n");
  fprintf(output, " %s   %s\n", Paper::cite(Paper::Bakken),
          Paper::print(Paper::Bakken));
  fprintf(output, "  ** Holonomic terms (MMFF94) **\n");
  fprintf(output, " %s   %s\n", Paper::cite(Paper::Halgren),
          Paper::print(Paper::Halgren));
  fprintf(output, "  ** Additional Parameter **\n");
  fprintf(output, " %s   %s\n", Paper::cite(Paper::Cornilescu),
          Paper::print(Paper::Cornilescu));
  fflush(output);
}

void
makeTITANIAresults(FILE *output,
                   BasicInformation &baseInformation,
                   Flags &flags,
                   Molecule &CurrMol,
                   Structure *CurrStruc)
{
  AlignmentOutput *aliOut = NULL;
  Eigen::MatrixXd C = CurrStruc->getCosineMatrix(StructureOptions::Optimized);
  Eigen::MatrixXd rdcs = CurrMol.getRDCmatrix();
  Eigen::MatrixXd rdcs_scaled =
      CurrStruc->getRDCmatrix(rdcMatrixOptions::Scaled);
  Eigen::MatrixXd Dmax = CurrStruc->getDmaxMatrix(StructureOptions::Optimized);
  Eigen::MatrixXd D    = CurrStruc->getRDCmatrix(rdcMatrixOptions::Unscaled);
  // TODO Q-Faktor, Vc output auf N spalten limitieren (siehe Frequenzanalyse
  // von ORCA etc.)

  if (!flags.skipSCRM)
  {
    aliOut = new AlignmentOutput(CurrMol.getNOR(), CurrMol.getNORsets());
    if (flags.calculateFullMatrix)
      aliOut->w = Eigen::MatrixXd::Identity(CurrMol.getNOR(), CurrMol.getNOR());
    else
      aliOut->w = CurrMol.getWmatrix();


    SaupeEigenSystems(C, aliOut->w, rdcs_scaled, aliOut->SaupeTensor,
                      aliOut->SaupeEigenValues, aliOut->EulerAngles,
                      aliOut->SaupeEigenVectors, flags.calculateFullMatrix);
    aliOut->rdc_calc = Dcalc(Dmax, C, aliOut->SaupeTensor);
    if (!flags.calculateFullMatrix)
    {
      for (unsigned int r = 0; r < CurrMol.getNOR(); ++r)
      {
        for (unsigned int s = 0; s < CurrMol.getNORsets(); ++s)
        {
          if (aliOut->w(r, s) == .0)
            aliOut->rdc_calc(r, s) = .0;
        }
      }
    }
    aliOut->Q_factor = DDcalc2Q(rdcs, aliOut->rdc_calc);

#pragma omp task shared(aliOut, CurrMol, baseInformation, flags)
    {
      for (RDCset *CurrSet = CurrMol.getHeadSet(); CurrSet;
           CurrSet         = CurrSet->getNext())
      {
#pragma omp task shared(aliOut, CurrMol, baseInformation, flags) \
    firstprivate(CurrSet)
        generateAlignmentOutput(CurrMol, CurrMol.getTailStruc(), CurrSet,
                                aliOut, baseInformation, flags);
      }
    }
  }

  if (!flags.skipSCRM)
  {
    information_header(output, "Results of TITANIA's iterative optimization");
    add_Qfactors(output, CurrMol, baseInformation);
    add_chiral_volumes(output, CurrMol, baseInformation, flags);

    fprintf(output, "\n %s \n", std::string(header_width, '=').c_str());
    frame_string(output, "Results of TITANIA's iterative optimization");
    frame_string(output, "Structure Parameter");
    fprintf(output, " %s \n\n", std::string(header_width, '=').c_str());

    add_TITANIA_angles(output, CurrMol, CurrStruc, baseInformation);
    if (flags.monteCarloBootstrapping)
      add_MonteCarlo_angles(output, CurrStruc, baseInformation, flags,
                            baseInformation.MonteCarloOutput);

    add_structure(output, CurrStruc, baseInformation, flags,
                  StructureOptions::Optimized);
    add_cosine_matrix(output, CurrStruc, baseInformation, flags,
                      StructureOptions::Optimized);

    fprintf(output, "\n %s \n", std::string(header_width, '=').c_str());
    frame_string(output, "Results of TITANIA's iterative optimization");
    frame_string(output, "Motional Analysis");
    fprintf(output, " %s \n\n", std::string(header_width, '=').c_str());

    add_Srdc(output, CurrStruc, baseInformation, flags, D, aliOut->rdc_calc);
  }
  fflush(output);
}

void
makeRuntimeInformation(FILE *output,
                       BasicInformation &baseInformation,
                       Flags &flags,
                       Molecule &CurrMol)
{
  information_header(output, "TITANIA runtime information");
  add_stop_criteria(output, baseInformation, flags);
  add_generated_files(output, baseInformation, flags, CurrMol);
  add_runtime_metrics(output, baseInformation, flags);
  fflush(output);
}

void
aliout_RDCs(FILE *aliOut,
            BasicInformation &baseInformation,
            Flags &flags,
            RDCdata *CurrRDC,
            AlignmentOutput *alignmentOutput,
            Eigen::MatrixXd &D_exp,
            const unsigned int index,
            const unsigned int NOR)
{
  int temp_length;
  fprintf(aliOut, " RDCs used for the optimization: \n");
  fprintf(aliOut, "\n");
  temp_length = 2 * baseInformation.output_format.string_m +
                3 * baseInformation.output_format.number_m +
                baseInformation.output_format.number_l + 5;
  fprintf(aliOut, "    %*s | ", temp_length,
          centerCstring("Input information", temp_length));
  temp_length = 2 * baseInformation.output_format.number_m + 1;
  fprintf(aliOut, "%*s ", temp_length, centerCstring("Optimized", temp_length));

  if (flags.monteCarloBootstrapping)
  {
    temp_length = 4 * baseInformation.output_format.number_m + 7;
    fprintf(aliOut, "| %*s", temp_length,
            centerCstring("Monte Carlo Bootstrap values", temp_length));
    temp_length = 2 * baseInformation.output_format.number_m + 1;
    fprintf(aliOut, " | %*s ", temp_length,
            centerCstring("MC residuals", temp_length));
  }

  fprintf(aliOut, "\n");

  // RDC table subheader
  fprintf(aliOut, "    %-*s %-*s", baseInformation.output_format.string_m,
          "Nuc 1", baseInformation.output_format.string_m, "Nuc 2");
  fprintf(aliOut, "%*s ", baseInformation.output_format.number_l, "Dmax");
  fprintf(aliOut, "%*s", 3 + 2 * baseInformation.output_format.number_m,
          " ( D +- err ) / Hz");
  fprintf(aliOut, "%*s | ", baseInformation.output_format.number_m, "w(eff)");
  fprintf(aliOut, "%*s ", baseInformation.output_format.number_m, "D(calc)");
  fprintf(aliOut, "%*s ", baseInformation.output_format.number_m, "res.");

  if (flags.monteCarloBootstrapping)
  {
    fprintf(aliOut, "| %*s", 3 + 2 * baseInformation.output_format.number_m,
            "D(inp) / Hz");
    fprintf(aliOut, " %*s", 3 + 2 * baseInformation.output_format.number_m,
            "D(calc) / Hz");
    fprintf(aliOut, " | %*s %*s", baseInformation.output_format.number_m,
            "res(inp)", baseInformation.output_format.number_m, "res(calc)");
  }

  fprintf(aliOut, "\n");

  // RDC table
  for (unsigned int r = 0; r < NOR; ++r)
  {
    fprintf(aliOut, "    %-*s %-*s", baseInformation.output_format.string_m,
            CurrRDC->getAtom1()->getIdentifier().c_str(),
            baseInformation.output_format.string_m,
            CurrRDC->getAtom2()->getIdentifier().c_str());
    fprintf(aliOut, "%*.*f ", baseInformation.output_format.number_l,
            baseInformation.output_format.prec_s,
            CurrRDC->getDmax(StructureOptions::Optimized));
    fprintf(aliOut, "%*.*f +-%*.*f", baseInformation.output_format.number_m,
            baseInformation.output_format.prec_m, CurrRDC->getD(),
            baseInformation.output_format.number_m,
            baseInformation.output_format.prec_m, CurrRDC->getDeltaD());
    fprintf(aliOut, "%*.*f | ", baseInformation.output_format.number_m,
            baseInformation.output_format.prec_m,
            CurrRDC->getEffectiveWeight());
    D_exp(r, 0) = CurrRDC->getD();

    fprintf(aliOut, "%*.*f ", baseInformation.output_format.number_m,
            baseInformation.output_format.prec_m,
            alignmentOutput->rdc_calc(r, index));
    fprintf(aliOut, "%*.*f ", baseInformation.output_format.number_m,
            baseInformation.output_format.prec_m,
            (alignmentOutput->rdc_calc(r, index) - CurrRDC->getD()));

    if (flags.monteCarloBootstrapping)
    {
      fprintf(aliOut, "| %*.*f +-%*.*f", baseInformation.output_format.number_m,
              baseInformation.output_format.prec_m,
              baseInformation.MonteCarloOutput.D_mean(r, index),
              baseInformation.output_format.number_m,
              baseInformation.output_format.prec_m,
              baseInformation.MonteCarloOutput.D_sigm(r, index));
      fprintf(aliOut, " %*.*f +-%*.*f", baseInformation.output_format.number_m,
              baseInformation.output_format.prec_m,
              baseInformation.MonteCarloOutput.D_calc_mean(r, index),
              baseInformation.output_format.number_m,
              baseInformation.output_format.prec_m,
              baseInformation.MonteCarloOutput.D_calc_sigm(r, index));
      fprintf(
          aliOut, " | %*.*f %*.*f", baseInformation.output_format.number_m,
          baseInformation.output_format.prec_m,
          (baseInformation.MonteCarloOutput.D_mean(r, index) - CurrRDC->getD()),
          baseInformation.output_format.number_m,
          baseInformation.output_format.prec_m,
          (baseInformation.MonteCarloOutput.D_calc_mean(r, index) -
           CurrRDC->getD()));
    }

    fprintf(aliOut, "\n");
    CurrRDC = CurrRDC->getNext();
  }
  fprintf(aliOut, "\n");
  fflush(aliOut);
}

void
aliout_correlations(FILE *aliOut,
                    BasicInformation &baseInformation,
                    double Q,
                    double X_2,
                    double m,
                    double b,
                    double R_2)
{
  fprintf(aliOut, " Correlation of experimental and calculated RDCs: \n");
  fprintf(aliOut, "%*s: %*.*f \n", baseInformation.output_format.string_xl,
          "Q-factor", baseInformation.output_format.number_l,
          baseInformation.output_format.prec_m, Q);
  fprintf(aliOut, "%*s: %*.*f \n", baseInformation.output_format.string_xl,
          "Chi^2", baseInformation.output_format.number_l,
          baseInformation.output_format.prec_m, X_2);
  fprintf(aliOut, "\n");

  fprintf(aliOut, " Linear Fit of experimental and calculated RDCs: \n");
  fprintf(aliOut, "%*s: %*.*f \n", baseInformation.output_format.string_xl,
          "Slope", baseInformation.output_format.number_l,
          baseInformation.output_format.prec_m, m);
  fprintf(aliOut, "%*s: %*.*f \n", baseInformation.output_format.string_xl,
          "Offset", baseInformation.output_format.number_l,
          baseInformation.output_format.prec_m, b);
  fprintf(aliOut, "%*s: %*.*f \n", baseInformation.output_format.string_xl,
          "Pearson R^2", baseInformation.output_format.number_l,
          baseInformation.output_format.prec_m, R_2);
  fprintf(aliOut, "\n");
  fflush(aliOut);
}

void
aliout_alignment(FILE *aliOut,
                 BasicInformation &baseInformation,
                 Flags &flags,
                 AlignmentOutput *alignmentOutput,
                 const Eigen::MatrixXd Saupe,
                 const unsigned int index)
{
  double f = 1.0;
  int i, j, temp_length;

  Eigen::MatrixXd S_MC_Gauss;
  Eigen::Vector3d Euler_MC, Euler_MC_sigm;
  Eigen::MatrixXd S_MC      = Eigen::MatrixXd::Zero(3, 3);
  Eigen::MatrixXd S_MC_sigm = Eigen::MatrixXd::Zero(3, 3);
  Eigen::MatrixXd EVal_MC, EVal_MC_sigm;
  Eigen::MatrixXd EVec_MC, EVec_MC_sigm;
  Eigen::Matrix3d EigenVectors;

  // Gather Monte-Carlo Errors
  if (flags.monteCarloBootstrapping)
  {
    S_MC(0, 0) = baseInformation.MonteCarloOutput.Saupe_tensor(0, index);
    S_MC(1, 1) = baseInformation.MonteCarloOutput.Saupe_tensor(1, index);
    S_MC(2, 2) = baseInformation.MonteCarloOutput.Saupe_tensor(2, index);
    S_MC(1, 0) = S_MC(0, 1) =
        baseInformation.MonteCarloOutput.Saupe_tensor(3, index);
    S_MC(2, 0) = S_MC(0, 2) =
        baseInformation.MonteCarloOutput.Saupe_tensor(4, index);
    S_MC(2, 1) = S_MC(1, 2) =
        baseInformation.MonteCarloOutput.Saupe_tensor(5, index);

    S_MC_sigm(0, 0) =
        baseInformation.MonteCarloOutput.Saupe_tensor_sigm(0, index);
    S_MC_sigm(1, 1) =
        baseInformation.MonteCarloOutput.Saupe_tensor_sigm(1, index);
    S_MC_sigm(2, 2) =
        baseInformation.MonteCarloOutput.Saupe_tensor_sigm(2, index);
    S_MC_sigm(1, 0) = S_MC_sigm(0, 1) =
        baseInformation.MonteCarloOutput.Saupe_tensor_sigm(3, index);
    S_MC_sigm(2, 0) = S_MC_sigm(0, 2) =
        baseInformation.MonteCarloOutput.Saupe_tensor_sigm(4, index);
    S_MC_sigm(2, 1) = S_MC_sigm(1, 2) =
        baseInformation.MonteCarloOutput.Saupe_tensor_sigm(5, index);

    S_MC_Gauss = estimateSaupeErrors(
        baseInformation.MonteCarloOutput.Aligns.col(index), S_MC, EVec_MC,
        EVal_MC, Euler_MC,
        baseInformation.MonteCarloOutput.Aligns_sigm.col(index), EVec_MC_sigm,
        EVal_MC_sigm, Euler_MC_sigm);
  }

  fprintf(aliOut, " Information on orientation:\n");
  fprintf(aliOut, "\n");

  // BLOCK: Saupe Vector
  fprintf(aliOut, " Saupe vector:  %*s %*s %*s %*s %*s \n",
          baseInformation.output_format.string_l,
          centerCstring("zz", baseInformation.output_format.string_l),
          baseInformation.output_format.string_l,
          centerCstring("xx-yy", baseInformation.output_format.string_l),
          baseInformation.output_format.string_l,
          centerCstring("xy", baseInformation.output_format.string_l),
          baseInformation.output_format.string_l,
          centerCstring("xz", baseInformation.output_format.string_l),
          baseInformation.output_format.string_l,
          centerCstring("yz", baseInformation.output_format.string_l));
  fprintf(aliOut, "                ");
  for (j = 0; j < SAUPE_ELEMENTS_; ++j)
  {
    if (j == SAUPE_ZZ_)
      f = 1.0;
    else if (j == SAUPE_XX_YY_)
      f = sqrt(3.0);
    else if (j == SAUPE_XY_)
      f /= 2.0;
    fprintf(aliOut, "%*.*e ", baseInformation.output_format.number_l,
            baseInformation.output_format.prec_l,
            f * alignmentOutput->SaupeTensor(j, index));
  }
  fprintf(aliOut, "\n");
  fprintf(aliOut, "\n");
  if (flags.monteCarloBootstrapping)
  {
    fprintf(aliOut, " Monte-Carlo error estimation (%d steps): %s \n",
            baseInformation.MonteCarloOutput.steps,
            Paper::cite(Paper::Losonczi));
    fprintf(aliOut, " MC-mean:       ");
    for (j = 0; j < SAUPE_ELEMENTS_; ++j)
    {
      if (j == SAUPE_ZZ_)
        f = 1.0;
      else if (j == SAUPE_XX_YY_)
        f = sqrt(3.0);
      else if (j == SAUPE_XY_)
        f /= 2.0;
      fprintf(aliOut, "%*.*e ", baseInformation.output_format.number_l,
              baseInformation.output_format.prec_l,
              f * baseInformation.MonteCarloOutput.Aligns(j, index));
    }
    fprintf(aliOut, "\n");
    fprintf(aliOut, " abs. error:    ");
    for (j = 0; j < SAUPE_ELEMENTS_; ++j)
    {
      if (j == SAUPE_ZZ_)
        f = 1.0;
      else if (j == SAUPE_XX_YY_)
        f = sqrt(3.0);
      else if (j == SAUPE_XY_)
        f /= 2.0;
      fprintf(aliOut, "%*.*e ", baseInformation.output_format.number_l,
              baseInformation.output_format.prec_l,
              f * baseInformation.MonteCarloOutput.Aligns_sigm(j, index));
    }
    fprintf(aliOut, "\n");
    fprintf(aliOut, " %%-error:   ");
    for (j = 0; j < SAUPE_ELEMENTS_; ++j)
    {
      if (j == SAUPE_ZZ_)
        f = 1.0;
      else if (j == SAUPE_XX_YY_)
        f = sqrt(3.0);
      else if (j == SAUPE_XY_)
        f /= 2.0;
      double rel_Error =
          fabs(baseInformation.MonteCarloOutput.Aligns_sigm(j, index) /
               baseInformation.MonteCarloOutput.Aligns(j, index)) *
          100.0;
      fprintf(aliOut, "%*.*f ", baseInformation.output_format.number_l,
              baseInformation.output_format.prec_m, rel_Error);
    }
    fprintf(aliOut, "\n");
    fprintf(aliOut, "\n");
  }

  // BLOCK: Saupe matrices
  fprintf(aliOut, " Saupe order matrices: \n");
  fprintf(aliOut, "\n");
  temp_length = 3 * baseInformation.output_format.number_l + 3;
  fprintf(aliOut, "%*s", temp_length,
          centerCstring("Final SVD tensor", temp_length));
  if (flags.monteCarloBootstrapping)
    fprintf(aliOut, " |  %*s", temp_length,
            centerCstring("Final MC tensor", temp_length));
  fprintf(aliOut, "\n");

  for (i = 0; i < NUMBER_OF_AXIS_; ++i)
  {
    for (j = 0; j < NUMBER_OF_AXIS_; ++j)
    {
      fprintf(aliOut, "%*.*e ", baseInformation.output_format.number_l,
              baseInformation.output_format.prec_l, Saupe(i, j));
    }
    if (flags.monteCarloBootstrapping)
    {
      fprintf(aliOut, " |  ");
      for (j = 0; j < NUMBER_OF_AXIS_; ++j)
      {
        fprintf(aliOut, "%*.*e ", baseInformation.output_format.number_l,
                baseInformation.output_format.prec_l, S_MC(i, j));
      }
    }
    fprintf(aliOut, "\n");
  }
  fprintf(aliOut, "\n");

  // BLOCK: Saupe matrices error
  if (flags.monteCarloBootstrapping)
  {
    fprintf(aliOut, "%*s |  %*s \n", temp_length, " ", temp_length,
            centerCstring("Monte-Carlo error", temp_length));
    for (i = 0; i < NUMBER_OF_AXIS_; ++i)
    {
      fprintf(aliOut, "%*s |  ", temp_length, " ");
      for (j = 0; j < NUMBER_OF_AXIS_; ++j)
      {
        fprintf(aliOut, "%*.*e ", baseInformation.output_format.number_l,
                baseInformation.output_format.prec_l, S_MC_sigm(i, j));
      }
      fprintf(aliOut, "\n");
    }
    fprintf(aliOut, "\n");
    fprintf(aliOut, "%*s |  %*s \n", temp_length, " ", temp_length,
            centerCstring("Gauss propagation of errors", temp_length));
    for (i = 0; i < NUMBER_OF_AXIS_; ++i)
    {
      fprintf(aliOut, "%*s |  ", temp_length, " ");
      for (j = 0; j < NUMBER_OF_AXIS_; ++j)
      {
        fprintf(aliOut, "%*.*e ", baseInformation.output_format.number_l,
                baseInformation.output_format.prec_l, S_MC_Gauss(i, j));
      }
      fprintf(aliOut, "\n");
    }
    fprintf(aliOut, "\n");
  }

  // BLOCK: Saupe eigensystem

  fprintf(aliOut, " Saupe eigenvalues (xx, yy, zz): \n");
  fprintf(aliOut, "\n");
  fprintf(aliOut, "%*s", temp_length,
          centerCstring("Final SVD eigenvalues", temp_length));
  if (flags.monteCarloBootstrapping)
    fprintf(aliOut, " |  %*s", temp_length,
            centerCstring("Final MC eigenvalues", temp_length));
  fprintf(aliOut, "\n");
  for (i = 0; i < NUM_EIGENVALUES_; ++i)
  {
    fprintf(aliOut, "%*.*e ", baseInformation.output_format.number_l,
            baseInformation.output_format.prec_l,
            alignmentOutput->SaupeEigenValues(i, index));
  }
  if (flags.monteCarloBootstrapping)
  {
    fprintf(aliOut, " |  ");
    for (i = 0; i < NUM_EIGENVALUES_; ++i)
    {
      fprintf(aliOut, "%*.*e ", baseInformation.output_format.number_l,
              baseInformation.output_format.prec_l, EVal_MC(i));
    }
    fprintf(aliOut, "\n");
    fprintf(aliOut, "\n");
    fprintf(aliOut, "%*s |  %*s \n", temp_length, " ", temp_length,
            centerCstring("Gauss propagation of errors", temp_length));
    fprintf(aliOut, "%*s |  ", temp_length, " ");
    for (i = 0; i < 3; ++i)
    {
      fprintf(aliOut, "%*.*e ", baseInformation.output_format.number_l,
              baseInformation.output_format.prec_l, EVal_MC_sigm(i));
    }
  }
  fprintf(aliOut, "\n");
  fprintf(aliOut, "\n");
  fprintf(aliOut,
          " Saupe eigenvectors (1. col: xx, 2. col: yy, 3. col: zz): \n");
  fprintf(aliOut, "%*s \n", temp_length,
          centerCstring("Final SVD eigenvectors", temp_length));
  for (j = 0; j < NUMBER_OF_AXIS_; ++j)
  {
    for (i = 0; i < NUMBER_OF_AXIS_; ++i)
    {
      fprintf(aliOut, "%*.*e ", baseInformation.output_format.number_l,
              baseInformation.output_format.prec_l,
              alignmentOutput->SaupeEigenVectors(j + i * 3, index));
      EigenVectors(j, i) = alignmentOutput->SaupeEigenVectors(j + i * 3, index);
    }
    fprintf(aliOut, "\n");
  }
  fprintf(aliOut, "\n");
  fprintf(aliOut, "\n");

  // BLOCK: Euler angles
  fprintf(aliOut, " Euler Angles (passive ZYZ rotation) / deg: \n");
  fprintf(aliOut, "\n");
  fprintf(aliOut, "%*s", temp_length,
          centerCstring("Final SVD Euler angles", temp_length));
  if (flags.monteCarloBootstrapping)
    fprintf(aliOut, " |  %*s", temp_length,
            centerCstring("Final MC Euler angles", temp_length));
  fprintf(aliOut, "\n");

  fprintf(aliOut, "%*s %*s %*s ", baseInformation.output_format.number_l,
          "alpha", baseInformation.output_format.number_l, "beta",
          baseInformation.output_format.number_l, "gamma");
  if (flags.monteCarloBootstrapping)
    fprintf(aliOut, " |  %*s %*s %*s ", baseInformation.output_format.number_l,
            "alpha", baseInformation.output_format.number_l, "beta",
            baseInformation.output_format.number_l, "gamma");
  fprintf(aliOut, "\n");

  for (i = 0; i < NUM_EULER_ANGLES_; ++i)
  {
    fprintf(aliOut, "%*.*f ", baseInformation.output_format.number_l,
            baseInformation.output_format.prec_l,
            -rad2deg(alignmentOutput->EulerAngles(i, index)));
  }

  if (flags.monteCarloBootstrapping)
  {
    fprintf(aliOut, " |  ");
    for (i = 0; i < NUM_EULER_ANGLES_; ++i)
    {
      fprintf(aliOut, "%*.*f ", baseInformation.output_format.number_l,
              baseInformation.output_format.prec_l, rad2deg(Euler_MC(i)));
    }
    fprintf(aliOut, "\n");
    fprintf(aliOut, "\n");
    fprintf(aliOut, "%*s |  %*s \n", temp_length, " ", temp_length,
            centerCstring("Gauss propagation of errors", temp_length));
    fprintf(aliOut, "%*s |  ", temp_length, " ");
    for (i = 0; i < 3; ++i)
    {
      fprintf(aliOut, "%*.*f ", baseInformation.output_format.number_l,
              baseInformation.output_format.prec_l, rad2deg(Euler_MC_sigm(i)));
    }
  }
  fprintf(aliOut, "\n");
  fprintf(aliOut, "\n");

  double Da, Dr, Aa, Ar, Da_NH, R, eta, GDO;
  double alpha_zyz, beta_zyz, gamma_zyz;
  double alpha_xyz, beta_xyz, gamma_xyz;
  double quat_x, quat_y, quat_z, quat_w;

  fprintf(aliOut, " Additional tensor parameters: \n");
  fprintf(aliOut, "\n");

  EigVecs2Quat(EigenVectors, quat_x, quat_y, quat_z, quat_w);
  Eul_FromHMatrix(EigenVectors, alpha_zyz, beta_zyz, gamma_zyz, EulOrdZYZr,
                  true);
  Eul_FromHMatrix(EigenVectors, alpha_xyz, beta_xyz, gamma_xyz, EulOrdXYZr,
                  true);
  Aa    = alignmentOutput->SaupeEigenValues(2, index);
  Da    = Aa * 0.5;
  Da_NH = Da * 21585.19;
  Dr    = (alignmentOutput->SaupeEigenValues(0, index) -
        alignmentOutput->SaupeEigenValues(1, index)) /
       3.0;
  Ar  = 2.0 * Dr;
  R   = Ar / Aa;
  eta = 1.5 * R;
  GDO = sqrt(2.0 / 3.0) * alignmentOutput->SaupeEigenValues.col(index).norm();

  fprintf(aliOut, " %*s alpha = %*.*f beta = %*.*f gamma = %*.*f %s \n\n",
          baseInformation.output_format.string_l,
          "Euler(ZYZ,active):", baseInformation.output_format.number_m,
          baseInformation.output_format.prec_m, rad2deg(alpha_zyz),
          baseInformation.output_format.number_m,
          baseInformation.output_format.prec_m, rad2deg(beta_zyz),
          baseInformation.output_format.number_m,
          baseInformation.output_format.prec_m, rad2deg(gamma_zyz),
          "(relative orientation, counterclockwise active "
          "z[gamma]-y'[beta]-z''[alpha])");

  fprintf(aliOut, " %*s alpha = %*.*f beta = %*.*f gamma = %*.*f %s \n\n",
          baseInformation.output_format.string_l,
          "Euler(XYZ,active):", baseInformation.output_format.number_m,
          baseInformation.output_format.prec_m, rad2deg(alpha_xyz),
          baseInformation.output_format.number_m,
          baseInformation.output_format.prec_m, rad2deg(beta_xyz),
          baseInformation.output_format.number_m,
          baseInformation.output_format.prec_m, rad2deg(gamma_xyz),
          "(Cardan rotation with alpha = theta, beta = phi, gamma = psi)");

  fprintf(aliOut, " %-*s %*.*e %s \n\n", baseInformation.output_format.string_l,
          "Da:", baseInformation.output_format.number_l,
          baseInformation.output_format.prec_m, Da,
          "(Axial component of the Saupe tensor: 0.5*SaupeEigenValue[zz])");

  fprintf(aliOut, " %-*s %*.*e %s \n\n", baseInformation.output_format.string_l,
          "Dr:", baseInformation.output_format.number_l,
          baseInformation.output_format.prec_m, Dr,
          "(Rhombic component of the Saupe tensor: "
          "1/3*(SaupeEigenValue[xx]-SaupeEigenValue[yy])");

  fprintf(aliOut, " %-*s %*.*e %s \n\n", baseInformation.output_format.string_l,
          "Aa:", baseInformation.output_format.number_l,
          baseInformation.output_format.prec_m, Aa,
          "(Axial component of the Alignment tensor: SaupeEigenValue[zz])");

  fprintf(aliOut, " %-*s %*.*e %s \n\n", baseInformation.output_format.string_l,
          "Da_NH:", baseInformation.output_format.number_l,
          baseInformation.output_format.prec_m, Da_NH,
          "(Da normalized to one-bond NH dipolar interaction: 21585.19*Da)");

  fprintf(aliOut, " %-*s %*.*e %s \n\n", baseInformation.output_format.string_l,
          "R:", baseInformation.output_format.number_l,
          baseInformation.output_format.prec_m, R,
          "(Rhombicity of the Alignment tensor: Ar/Aa (range:[0,2/3]))");


  fprintf(aliOut, " %-*s %*.*e ", baseInformation.output_format.string_l,
          "eta:", baseInformation.output_format.number_l,
          baseInformation.output_format.prec_m, eta);

  temp_length = baseInformation.output_format.string_l +
                baseInformation.output_format.number_l + 1;
  fprintf(aliOut, "%s \n %*s%s\n\n",
          "(Asymmetry parameter of the Alignment tensor: R/Rmax =", temp_length,
          " ",
          "  (SaupeEigenValue[xx]-SaupeEigenValue[yy])/SaupeEigenValue[zz] = "
          "1.5*R)");

  fprintf(aliOut, " %-*s %*.*e %s \n\n", baseInformation.output_format.string_l,
          "GDO:", baseInformation.output_format.number_l,
          baseInformation.output_format.prec_m, GDO,
          "(Generalized degree of order: sqrt(2/3)*|SaupeEigenValues|)");
}

void
generateAlignmentOutput(Molecule &CurrMol,
                        Structure *CurrStruc,
                        RDCset *CurrSet,
                        AlignmentOutput *alignmentOutput,
                        BasicInformation &baseInformation,
                        Flags &flags)
{
  RDCdata *CurrRDC = CurrSet->getHeadData();
  unsigned int NOR = CurrMol.getNOR();
  double m, b, R_2, Chi_2;
  Eigen::MatrixXd Saupe;
  Eigen::MatrixXd D_exp = Eigen::MatrixXd::Zero(CurrMol.getNOR(), 1);
  Eigen::MatrixXd D_calc;
  char *buf = (char *) malloc(STAN_BUF);
  snprintf(buf, STAN_BUF, "%s.%s.ali", baseInformation.outputFileName.c_str(),
           CurrSet->getLabel().c_str());

  unsigned int index = (CurrSet->getIndex() - 1);

  FILE *aliOut;

  aliOut = fopen(buf, "w");

  // HEADER
  fprintf(aliOut, " %s \n", std::string(header_width, '=').c_str());
  frame_string(aliOut, "This is the alignment information output file for:", 3);
  snprintf(buf, STAN_BUF, "Input file:  %s",
           baseInformation.inputFileName.c_str());
  frame_string(aliOut, buf, 3);
  snprintf(buf, STAN_BUF, "RDC set:     %s", CurrSet->getLabel().c_str());
  frame_string(aliOut, buf, 3);
  fprintf(aliOut, " %s \n", std::string(header_width, '=').c_str());
  fprintf(aliOut, "\n");

  CurrRDC = CurrSet->getHeadData();

  // RDC table
  aliout_RDCs(aliOut, baseInformation, flags, CurrRDC, alignmentOutput, D_exp,
              index, NOR);

  // Calculate alignment stuff
  Saupe  = Sv2St(alignmentOutput->SaupeTensor.col(index));
  D_calc = alignmentOutput->rdc_calc.col(index);
  R_2    = linear_Regression(D_exp, D_calc, m, b);
  Chi_2  = Chi_square(D_exp, D_calc, CurrSet->get_sigma_square());

  // Correlation parameter
  aliout_correlations(aliOut, baseInformation, alignmentOutput->Q_factor(index),
                      Chi_2, m, b, R_2);
  aliout_alignment(aliOut, baseInformation, flags, alignmentOutput, Saupe,
                   index);

  fclose(aliOut);
}

void
output_big_data(Molecule &CurrMol,
                BasicInformation &baseInformation,
                Flags &flags)
{
  Eigen::MatrixXd out;
  int i, j;
  double x, y;
  std::complex<double> z;
  RDCdata *RD;
  RDCset *RS;
  SphericalHarmonics *Y;
  Structure *S;
  std::fstream bd;

  bd.open((baseInformation.outputFileName + ".dat"),
          std::ios::binary | std::ios::out | std::ios::trunc);
  bd << "TITANIA data file;\n;\n";


  bd << "Initial RDC matrix;\n;";
  RS  = CurrMol.getHeadSet();
  RD  = RS->getHeadData();
  out = CurrMol.getRDCmatrix();
  while (RS)
  {
    bd << RS->getLabel() << ";";
    RS = RS->getNext();
  }
  bd << std::endl;
  while (RD)
  {
    bd << RD->getAtom1()->getIdentifier() << RD->getAtom2()->getIdentifier()
       << ";";
    bd << out.row(RD->getInputIndex()).format(big_data) << std::endl;
    RD = RD->getNext();
  }


  bd << ";\nScaled RDC matrix of last structure;\n;";
  RS  = CurrMol.getHeadSet();
  RD  = RS->getHeadData();
  out = CurrMol.getTailStruc()->getRDCmatrix(rdcMatrixOptions::Scaled);
  while (RS)
  {
    bd << RS->getLabel() << ";";
    RS = RS->getNext();
  }
  bd << std::endl;
  while (RD)
  {
    bd << RD->getAtom1()->getIdentifier() << RD->getAtom2()->getIdentifier()
       << ";";
    bd << out.row(RD->getInputIndex()).format(big_data) << std::endl;
    RD = RD->getNext();
  }


  bd << ";\nPolar angles (deg);\n;";
  S = CurrMol.getHeadStruc();
  for (Y = S->getHeadYmatrix()->getHeadHarmonic(); Y; Y = Y->getNext())
  {
    bd << "theta(" << Y->getAtom1()->getIdentifier()
       << Y->getAtom2()->getIdentifier() << ");phi("
       << Y->getAtom1()->getIdentifier() << Y->getAtom2()->getIdentifier()
       << ");";
  }
  bd << std::endl;
  i = 0;
  bd << i++ << ";";
  for (Y = S->getHeadYmatrix()->getHeadHarmonic(); Y; Y = Y->getNext())
  {
    x = Y->getTheta(StructureOptions::Initial);
    y = Y->getPhi(StructureOptions::Initial);
    bd << rad2deg(x) << ";" << rad2deg(y) << ";";
  }
  bd << std::endl;
  while (S)
  {
    bd << i++ << ";";
    for (Y = S->getHeadYmatrix()->getHeadHarmonic(); Y; Y = Y->getNext())
    {
      x = Y->getTheta(StructureOptions::Optimized);
      y = Y->getPhi(StructureOptions::Optimized);
      bd << rad2deg(x) << ";" << rad2deg(y) << ";";
    }
    bd << std::endl;
    S = S->getNext();
  }


  bd << ";\nSaupe vectors;\n";
  RS               = CurrMol.getHeadSet();
  unsigned int NOR = CurrMol.getNOR(), NOS = CurrMol.getNORsets();
  Eigen::MatrixXd ST = Eigen::MatrixXd::Zero(SAUPE_ELEMENTS_, NOS);
  Eigen::MatrixXd w  = Eigen::MatrixXd::Zero(NOR, NOR);
  Eigen::MatrixXd weights;
  if (flags.calculateFullMatrix)
    weights = Eigen::MatrixXd::Identity(NOR, NOR);
  else
    weights = CurrMol.getWmatrix();
  Eigen::MatrixXd D = Eigen::MatrixXd::Zero(NOR, 1);
  Eigen::MatrixXd RDCs;
  while (RS)
  {
    bd << ";\n" << RS->getLabel() << "\n;Szz;Sxx-yy;Sxy;Sxz;Syz\n";
    S = CurrMol.getHeadStruc();
    i = 0;
    bd << i++ << ";";
    if (flags.calculateFullMatrix)
    {
      ST = BD2A_unref(S->getCosineMatrix(StructureOptions::Initial), weights,
                      S->getRDCmatrix(rdcMatrixOptions::Scaled));
    }
    else
    {
      RDCs = S->getRDCmatrix(rdcMatrixOptions::Scaled);
      for (unsigned int s = 0; s < NOS; ++s)
      {
        D = RDCs.col(s);
        for (unsigned int r = 0; r < NOR; ++r)
          w(r, r) = weights(r, s);
        ST.col(s) =
            BD2A_unref(S->getCosineMatrix(StructureOptions::Initial), w, D)
                .col(0);
      }
    }
    for (j = 0; j < 5; ++j)
    {
      double f = 1.0;
      if (j == 1)
        f = sqrt(3.0);
      if (j > 1)
        f = sqrt(3.0) / 2.0;
      bd << std::scientific << f * ST(j, RS->getIndex() - 1) << ";";
    }
    bd << std::endl;
    while (S)
    {
      bd << i++ << ";";

      if (flags.calculateFullMatrix)
      {
        ST = BD2A_unref(S->getCosineMatrix(StructureOptions::Optimized),
                        weights, S->getRDCmatrix(rdcMatrixOptions::Scaled));
      }
      else
      {
        RDCs = S->getRDCmatrix(rdcMatrixOptions::Scaled);
        for (unsigned int s = 0; s < NOS; ++s)
        {
          D = RDCs.col(s);
          for (unsigned int r = 0; r < NOR; ++r)
            w(r, r) = weights(r, s);
          ST.col(s) =
              BD2A_unref(S->getCosineMatrix(StructureOptions::Initial), w, D)
                  .col(0);
        }
      }

      for (j = 0; j < 5; ++j)
      {
        double f = 1.0;
        if (j == 1)
          f = sqrt(3.0);
        if (j > 1)
          f = sqrt(3.0) / 2.0;
        bd << std::scientific << f * ST(j, RS->getIndex() - 1) << ";";
      }
      bd << std::endl;
      S = S->getNext();
    }
    RS = RS->getNext();
  }


  bd << ";\nSaupe eigenvalues;\n";
  RS = CurrMol.getHeadSet();
  while (RS)
  {
    bd << ";\n" << RS->getLabel() << "\n;Sxx;Syy;Szz\n";
    S = CurrMol.getHeadStruc();
    S->getHeadYmatrix()->determineEuler(CurrMol, S, baseInformation, flags,
                                        StructureOptions::Initial);
    i = 0;
    bd << i++ << ";" << std::scientific
       << S->getSaupeEigenValues()
              .col(RS->getIndex() - 1)
              .transpose()
              .format(big_data)
       << std::endl;
    while (S)
    {
      S->getHeadYmatrix()->determineEuler(CurrMol, S, baseInformation, flags,
                                          StructureOptions::Optimized);
      bd << i++ << ";" << std::scientific
         << S->getSaupeEigenValues()
                .col(RS->getIndex() - 1)
                .transpose()
                .format(big_data)
         << std::endl;
      S = S->getNext();
    }
    RS = RS->getNext();
  }


  bd << ";\nEuler angles (ZYZ);\n";
  RS = CurrMol.getHeadSet();
  while (RS)
  {
    bd << ";\n" << RS->getLabel() << "\n;alpha(deg);beta(deg);gamma(deg)\n";
    S = CurrMol.getHeadStruc();
    S->getHeadYmatrix()->determineEuler(CurrMol, S, baseInformation, flags,
                                        StructureOptions::Initial);
    i   = 0;
    out = -S->getSaupeEulerAngles().col(RS->getIndex() - 1).transpose();
    bd << i++ << ";" << std::fixed << rad2deg(out).format(big_data)
       << std::endl;
    while (S)
    {
      S->getHeadYmatrix()->determineEuler(CurrMol, S, baseInformation, flags,
                                          StructureOptions::Optimized);
      out = -S->getSaupeEulerAngles().col(RS->getIndex() - 1).transpose();
      bd << i++ << ";" << std::fixed << rad2deg(out).format(big_data)
         << std::endl;
      S = S->getNext();
    }
    RS = RS->getNext();
  }


  bd << ";\nAxial component S^2(rdc);\n;";
  bd << "S(overall);";
  for (Y = CurrMol.getHeadStruc()->getHeadYmatrix()->getHeadHarmonic(); Y;
       Y = Y->getNext())
  {
    bd << Y->getAtom1()->getIdentifier() << Y->getAtom2()->getIdentifier()
       << ";";
  }
  bd << std::endl;
  for (S = CurrMol.getHeadStruc(); S; S = S->getNext())
  {
    bd << S->getIndex() << ";" << S->getHeadYmatrix()->getSoverall() << ";";
    for (Y = S->getHeadYmatrix()->getHeadHarmonic(); Y; Y = Y->getNext())
    {
      bd << std::fixed << Y->getS2rdc() << ";";
    }
    bd << std::endl;
  }


  bd << ";\nAsymmetry parameter eta(rdc);\n;";
  for (Y = CurrMol.getHeadStruc()->getHeadYmatrix()->getHeadHarmonic(); Y;
       Y = Y->getNext())
  {
    bd << Y->getAtom1()->getIdentifier() << Y->getAtom2()->getIdentifier()
       << ";";
  }
  bd << std::endl;

  for (S = CurrMol.getHeadStruc(); S; S = S->getNext())
  {
    bd << S->getIndex() << ";";
    for (Y = S->getHeadYmatrix()->getHeadHarmonic(); Y; Y = Y->getNext())
    {
      bd << std::fixed << Y->getetardc() << ";";
    }
    bd << std::endl;
  }


  bd << ";\nDirection of anisotropic motion (deg);\n;";
  for (Y = CurrMol.getHeadStruc()->getHeadYmatrix()->getHeadHarmonic(); Y;
       Y = Y->getNext())
  {
    bd << Y->getAtom1()->getIdentifier() << Y->getAtom2()->getIdentifier()
       << ";";
  }
  bd << std::endl;

  for (S = CurrMol.getHeadStruc(); S; S = S->getNext())
  {
    bd << S->getIndex() << ";";
    for (Y = S->getHeadYmatrix()->getHeadHarmonic(); Y; Y = Y->getNext())
    {
      bd << std::fixed << deg2rad(Y->getAnisoTheta()) << ";";
    }
    bd << std::endl;
  }
  bd.close();
  bd.clear();
}

void
outputXYZ(Atom *CurrAtom,
          BasicInformation &bI,
          StructureOptions opt,
          std::string label,
          bool mark)
{
  Coordinates *C;
  bI.xyzFile << CurrAtom->getParent()->getNOA() << std::endl
             << CurrAtom->getParent()->getLabel() << " " << label << std::endl;
  while (CurrAtom)
  {
    C = CurrAtom->getCoordinates(opt);
    if (mark && CurrAtom->getIdentifier().back() == 'a')
      bI.xyzFile << std::setw(6) << std::left << "T " << std::setprecision(5)
                 << std::setw(12) << std::right << std::fixed << C->x
                 << std::setprecision(5) << std::setw(12) << std::right
                 << std::fixed << C->y << std::setprecision(5) << std::setw(12)
                 << std::right << std::fixed << C->z << std::endl;
    else
      bI.xyzFile << std::setw(6) << std::left << CurrAtom->getElement()
                 << std::setprecision(5) << std::setw(12) << std::right
                 << std::fixed << C->x << std::setprecision(5) << std::setw(12)
                 << std::right << std::fixed << C->y << std::setprecision(5)
                 << std::setw(12) << std::right << std::fixed << C->z
                 << std::endl;
    CurrAtom = CurrAtom->getNext();
  }
  C = NULL;
  delete C;
  bI.xyzFile << std::endl;
}

void
titania2hotFCHT(Molecule &CurrMol,
                Structure *CurrStruc,
                RDCset *CurrSet,
                BasicInformation &baseInformation,
                Flags &flags)
{
  std::fstream hotFile;
  hotFile.open(
      (baseInformation.hotFCHTbase + "." + CurrSet->getLabel() + ".fcht"),
      std::ios::binary | std::ios::out);
  hotFile << "# Automatically generated RDC@hotFCHT input file.\n";
  hotFile << "# Used TITANIA analysis of " << CurrMol.getLabel() << " - "
          << CurrMol.getHeadStruc()->getLabel() << std::endl;
  hotFile << "# RDC set: " << CurrSet->getLabel() << std::endl << std::endl;

  // Generate all usefull Keywords
  hotFile << "CalculateOnlyRDCs = 1" << std::endl;
  hotFile << "#RDCMonteCarloOutput = 1" << std::endl;
  hotFile << "#RDCLongNumberFormat = 1" << std::endl;
  hotFile << "RDCScaleWeightsWithDmax = 1" << std::endl;
  if (flags.errorWeightInSVD)
    hotFile << "RDCNoErrorWeightInSVD = 1" << std::endl;
  else
    hotFile << "#RDCNoErrorWeightInSVD = 1" << std::endl;
  hotFile << "#RDCRotateCoordinates = 1" << std::endl;
  hotFile << "#RDCQuadrupolarSplittingFactor = 25.0" << std::endl;
  hotFile << "#RDCRotateCoordinates = 1" << std::endl;
  hotFile << std::endl;
  hotFile << "#Atom 1, Atom 2, D, DeltaD, w" << std::endl;
  hotFile << "ExperimentalRDCsD =" << std::endl;

  RDCdata *RDC = CurrSet->getHeadData();
  while (RDC)
  {
    hotFile.width(6);
    hotFile << std::left << RDC->getAtom1()->getIdentifier();
    hotFile.width(6);
    hotFile << std::left << RDC->getAtom2()->getIdentifier();
    hotFile << std::setprecision(3) << std::setw(12) << std::right << std::fixed
            << RDC->getD();
    hotFile << std::setprecision(3) << std::setw(12) << std::right << std::fixed
            << RDC->getDeltaD();
    hotFile << std::setprecision(3) << std::setw(12) << std::right << std::fixed
            << RDC->getWeight();
    hotFile << std::endl;
    RDC = RDC->getNext();
  }
  hotFile << std::endl;
  hotFile << "xyzcoordinates[" << CurrMol.getLabel() << "] =" << std::endl;
  hotFile << CurrMol.getNOA() << std::endl;
  hotFile << CurrStruc->getLabel() << std::endl;
  Atom *CurrAtom = CurrStruc->getHeadAtom();
  while (CurrAtom)
  {
    hotFile.width(6);
    hotFile << std::left << CurrAtom->getIdentifier();
    hotFile << std::setprecision(8) << std::setw(20) << std::right << std::fixed
            << CurrAtom->getCoordinates(StructureOptions::Optimized)->x;
    hotFile << std::setprecision(8) << std::setw(20) << std::right << std::fixed
            << CurrAtom->getCoordinates(StructureOptions::Optimized)->y;
    hotFile << std::setprecision(8) << std::setw(20) << std::right << std::fixed
            << CurrAtom->getCoordinates(StructureOptions::Optimized)->z;
    hotFile << std::endl;
    CurrAtom = CurrAtom->getNext();
  }

  hotFile.close();
}

void
titania2titania(Structure *CurrStruc, BasicInformation &baseInformation)
{
  std::fstream titaniaFile;
  titaniaFile.open((CurrStruc->getLabel() + ".tna"),
                   std::ios::binary | std::ios::out);

  std::fstream input;
  input.open(baseInformation.inputFileName, std::ios::binary | std::ios::in);

  std::string line, small, baseLabel;

  baseLabel = CurrStruc->getParent()->getHeadStruc()->getLabel();
  int a, noa;
  unsigned int i;
  bool inStruc = false;
  a            = -2;
  noa          = CurrStruc->getParent()->getNOA();

  Atom *CurrAtom = CurrStruc->getHeadAtom();
  Coordinates *C;
  while (std::getline(input, line))
  {
    small = line;
    for (i = 0; i < line.size(); ++i)
      small.at(i) = std::tolower(line.at(i));
    if (small.find("xyzcoordinates") != std::string::npos)
    {
      inStruc = true;
      titaniaFile << line << std::endl;
      continue;
    }
    if (inStruc)
    {
      switch (a)
      {
        case (-2):
          titaniaFile << noa << std::endl;
          break;
        case (-1):
          titaniaFile << CurrStruc->getLabel() << std::endl;
          break;
        default:
          C = CurrAtom->getCoordinates(StructureOptions::Optimized);
          titaniaFile << std::setw(6) << std::left << CurrAtom->getIdentifier()
                      << std::setprecision(5) << std::setw(12) << std::right
                      << std::fixed << C->x << std::setprecision(5)
                      << std::setw(12) << std::right << std::fixed << C->y
                      << std::setprecision(5) << std::setw(12) << std::right
                      << std::fixed << C->z << std::endl;
          CurrAtom = CurrAtom->getNext();
          break;
      }
      if (++a == noa)
        inStruc = false;
    }
    else
    {
      if (small.find(baseLabel) != std::string::npos)
        titaniaFile << CurrStruc->getLabel() << std::endl;
      else
        titaniaFile << line << std::endl;
    }
  }
  titaniaFile.close();
  input.close();
}

void
outputMatrix(Eigen::MatrixXd M,
             unsigned int w,
             unsigned int p,
             std::fstream &out,
             OutputOptions opt)
{
  for (unsigned int i = 0; i < M.rows(); ++i)
  {
    for (unsigned int j = 0; j < M.cols(); ++j)
    {
      switch (opt)
      {
        case (OutputOptions::undefined):
          out << std::setprecision(p) << std::setw(w) << M(i, j);
          break;
        case (OutputOptions::left):
          out << std::setprecision(p) << std::setw(w) << std::left << std::fixed
              << M(i, j);
          break;
        case (OutputOptions::center): {
          std::ostringstream ostr;
          ostr << std::setprecision(p) << std::fixed << M(i, j);
          out << centerString(ostr.str(), w);
          break;
        }
        case (OutputOptions::right):
          out << std::setprecision(p) << std::setw(w) << std::right
              << std::fixed << M(i, j);
          break;
        case (OutputOptions::scientific):
          out << std::setprecision(p) << std::setw(w) << std::right
              << std::scientific << M(i, j);
          break;
        default:
          out << std::setprecision(p) << std::setw(w) << M(i, j);
          break;
      }
    }
    out << std::endl;
  }
}

void
outputVector(Eigen::VectorXd M,
             unsigned int w,
             unsigned int p,
             std::fstream &out,
             OutputOptions opt,
             bool trans)
{
  for (unsigned int i = 0; i < M.rows(); ++i)
  {
    switch (opt)
    {
      case (OutputOptions::undefined):
        out << std::setprecision(p) << std::setw(w) << M(i);
        break;
      case (OutputOptions::left):
        out << std::setprecision(p) << std::setw(w) << std::left << std::fixed
            << M(i);
        break;
      case (OutputOptions::center): {
        std::ostringstream ostr;
        ostr << std::setprecision(p) << std::fixed << M(i);
        out << centerString(ostr.str(), w);
        break;
      }
      case (OutputOptions::right):
        out << std::setprecision(p) << std::setw(w) << std::right << std::fixed
            << M(i);
        break;
      default:
        out << std::setprecision(p) << std::setw(w) << M(i);
        break;
    }
    if (!trans)
      out << std::endl;
  }
}

std::string
centerString(std::string s, int fullSize)
{
  int s_length = s.length();
  if (s_length >= fullSize)
    return s;

  int fill        = fullSize - s_length;
  int pFill       = fill / 2;
  int aFill       = fill - pFill;
  std::string out = std::string(pFill, ' ') + s + std::string(aFill, ' ');
  return out;
}

const char *
centerCstring(const char *to_center, int size)
{
  char *centered = (char *) malloc(STAN_BUF);
  int pre        = (size - strlen(to_center)) / 2;
  snprintf(centered, STAN_BUF, "%*s%s%*s", pre, " ", to_center,
           static_cast<int>(size - strlen(to_center) - pre), " ");
  return centered;
}

std::string
ValueError2String(double V, double E, int preW, int p)
{
  int full = preW + p + 1;
  int fill;
  std::ostringstream ostr;
  std::string out = "";
  ostr << std::setprecision(p) << std::fixed << V;
  fill = full - ostr.str().size();
  if (fill > 0)
    out += std::string(fill, ' ');
  out += (ostr.str() + " +-");
  ostr.str(std::string());
  ostr << std::setprecision(p) << std::fixed << E;
  fill = full - ostr.str().size();
  if (fill > 0)
    out += std::string(fill, ' ');
  out += ostr.str();
  return out;
}

void
outputInitialTransformation(FILE *out,
                            Molecule &CurrMol,
                            BasicInformation &baseInformation)
{
  Eigen::Matrix3d InertiaTensor  = CurrMol.getInputInertiaTensor();
  Eigen::Vector3d InertiaEigVal  = CurrMol.getInputInertiaTensorEval();
  Eigen::Matrix3d InertiaEigVec  = CurrMol.getInputInertiaTensorEvec();
  Eigen::MatrixXd Transformation = CurrMol.getInputTransformation();

  fprintf(
      out,
      " Transformation matrix [input frame -> intertia tensor frame](4D): \n");
  for (int r = 0; r < NUMBER_OF_4D_AXIS_; ++r)
  {
    fprintf(out, "  ");
    for (int c = 0; c < NUMBER_OF_4D_AXIS_; ++c)
      fprintf(out, "%*.*e", baseInformation.output_format.number_xl,
              baseInformation.output_format.prec_l, Transformation(r, c));
    fprintf(out, "\n");
  }

  fprintf(out, "\n Inertia tensor in input frame: g*(1e-10 m)^2 / mol \n");
  for (int r = 0; r < NUMBER_OF_AXIS_; ++r)
  {
    fprintf(out, "  ");
    for (int c = 0; c < NUMBER_OF_AXIS_; ++c)
      fprintf(out, "%*.*e", baseInformation.output_format.number_xl,
              baseInformation.output_format.prec_l, InertiaTensor(r, c));
    fprintf(out, "\n");
  }

  fprintf(out, "\n Eigenvalues of inertia tensor: g*(1e-10 m)^2 / mol \n");
  fprintf(out, "  ");
  for (int c = 0; c < NUMBER_OF_AXIS_; ++c)
    fprintf(out, "%*.*e", baseInformation.output_format.number_xl,
            baseInformation.output_format.prec_l, InertiaEigVal(c));
  fprintf(out, "\n");

  fprintf(out, "\n Eigenvectors of inertia tensor in input frame: \n");
  for (int r = 0; r < NUMBER_OF_AXIS_; ++r)
  {
    fprintf(out, "  ");
    for (int c = 0; c < NUMBER_OF_AXIS_; ++c)
      fprintf(out, "%*.*e", baseInformation.output_format.number_xl,
              baseInformation.output_format.prec_l, InertiaEigVec(r, c));
    fprintf(out, "\n");
  }
  fprintf(out, "\n");
  fflush(out);
}

void
outputInertiaTensor(FILE *out,
                    Structure *CurrStruc,
                    BasicInformation &baseInformation,
                    StructureOptions options)
{
  Eigen::Matrix3d InertiaTensor =
      CurrStruc->CalculateInertiaTensor(options)->getTensor();
  Eigen::Matrix3d InertiaEigVec =
      CurrStruc->getInertiaTensor(options)->getEvecs();
  Eigen::Vector3d InertiaEigVal =
      CurrStruc->getInertiaTensor(options)->getEvals();

  fprintf(out, " Inertia tensor in Eckart frame: g*(1e-10 m)^2 / mol \n");
  for (int r = 0; r < NUMBER_OF_AXIS_; ++r)
  {
    fprintf(out, "  ");
    for (int c = 0; c < NUMBER_OF_AXIS_; ++c)
      fprintf(out, "%*.*e", baseInformation.output_format.number_xl,
              baseInformation.output_format.prec_l, InertiaTensor(r, c));
    fprintf(out, "\n");
  }

  fprintf(out, "\n Eigenvalues of inertia tensor: g*(1e-10 m)^2 / mol \n");
  fprintf(out, "  ");
  for (int c = 0; c < NUMBER_OF_AXIS_; ++c)
    fprintf(out, "%*.*e", baseInformation.output_format.number_xl,
            baseInformation.output_format.prec_l, InertiaEigVal(c));
  fprintf(out, "\n");

  fprintf(out, "\n Eigenvectors of inertia tensor in input frame: \n");
  for (int r = 0; r < NUMBER_OF_AXIS_; ++r)
  {
    fprintf(out, "  ");
    for (int c = 0; c < NUMBER_OF_AXIS_; ++c)
      fprintf(out, "%*.*e", baseInformation.output_format.number_xl,
              baseInformation.output_format.prec_l, InertiaEigVec(r, c));
    fprintf(out, "\n");
  }
  fprintf(out, "\n");
  fflush(out);
}

void
outputTrjSECONDA(FILE *out,
                 Structure *CurrStruc,
                 BasicInformation &baseInformation)
{
  bool initial = (CurrStruc->getIndex() == 0);
  fprintf(out, "\nFull SECONDA information on %s structure:\n",
          (initial ? "initial" : "final"));
  Eigen::VectorXd evals;
  Eigen::MatrixXd evecs;
  unsigned int q, index;
  unsigned int sblock       = 0;
  unsigned int sblock_start = 0;
  unsigned int NOEig =
      baseInformation.SECONDA_output_initial.Covariance_Evals.rows();
  bool split = (NOEig > seconda_sub_block_size);

  if (initial)
  {
    evals = baseInformation.SECONDA_output_initial.Covariance_Evals;
    evecs = baseInformation.SECONDA_output_initial.Covariance_Evecs;
  }
  else
  {
    evals = baseInformation.SECONDA_output_final.Covariance_Evals;
    evecs = baseInformation.SECONDA_output_final.Covariance_Evecs;
  }

  do
  {
    sblock_start = sblock;
    sblock += seconda_sub_block_size;
    fprintf(out, "%s Block of SECONDA eigenvector elements %*d - %*d \n",
            (sblock > seconda_sub_block_size ? "\n" : ""),
            baseInformation.output_format.number_s, (sblock_start + 1),
            baseInformation.output_format.number_s,
            (sblock > NOEig ? NOEig : sblock));

    fprintf(out, " %*s | %*s |", baseInformation.output_format.number_s, "n",
            baseInformation.output_format.number_l, "eigenvalue");
    for (index = sblock_start;
         index < baseInformation.NumberOfRDCs && index < sblock; ++index)
      fprintf(out, " %*s", baseInformation.output_format.number_m,
              centerCstring(to_cstring(index + 1, 0),
                            baseInformation.output_format.number_m));
    fprintf(out, "\n");
    for (q = 0; q < NOEig; ++q)
    {
      fprintf(out, " %*u | ", baseInformation.output_format.number_s, (q + 1));
      fprintf(out, "%*.*e |", baseInformation.output_format.number_l,
              baseInformation.output_format.prec_l, evals(q));
      for (index = sblock_start;
           index < baseInformation.NumberOfRDCs && index < sblock; ++index)
        fprintf(out, " %*.*f", baseInformation.output_format.number_m,
                baseInformation.output_format.prec_l, evecs(q, index));
      fprintf(out, "\n");
    }
  } while (sblock < NOEig);
  fprintf(out, "\n");
  fflush(out);
}

Eigen::MatrixXd
outputAlignments(FILE *out,
                 Structure *CurrStruc,
                 Molecule &CurrMol,
                 int iteration,
                 BasicInformation &baseInformation,
                 Flags &flags,
                 StructureOptions opt)
{
  MC_Output monte_carlo;
  RDCset *CurrSet;
  Eigen::MatrixXd rdcs = CurrStruc->getRDCmatrix(rdcMatrixOptions::Scaled);
  Eigen::MatrixXd w, C, rdc_calc;

  int index        = 0;
  unsigned int NOS = CurrMol.getNORsets();
  unsigned int NOR = CurrMol.getNOR();

  if (flags.calculateFullMatrix)
    w = Eigen::MatrixXd::Identity(CurrMol.getNOR(), CurrMol.getNOR());
  else
    w = CurrMol.getWmatrix();

  C = CurrStruc->getHeadYmatrix()->determineBmatrix(CurrMol, CurrStruc, opt);

  Eigen::MatrixXd A           = Eigen::MatrixXd::Zero(SAUPE_ELEMENTS_, NOS);
  Eigen::MatrixXd EulerAngles = Eigen::MatrixXd::Zero(NUM_EULER_ANGLES_, NOS);
  Eigen::MatrixXd EigenVals   = Eigen::MatrixXd::Zero(NUM_EIGENVALUES_, NOS);
  Eigen::MatrixXd EigenVecs   = Eigen::MatrixXd::Zero(NUM_EIGENVECTORS_, NOS);

  SaupeEigenSystems(C, w, rdcs, A, EigenVals, EulerAngles, EigenVecs,
                    flags.calculateFullMatrix);

  rdc_calc = rdcs = CurrStruc->getRDCmatrix(rdcMatrixOptions::Unscaled);

  if (flags.calculateFullMatrix)
    w = Eigen::MatrixXd::Identity(NOR, NOR);
  else
    w = CurrMol.getWmatrix();

  Eigen::VectorXd Q = CurrStruc->Qfacs(rdcs, rdc_calc, opt, w, flags);
  index             = 0;

  fprintf(out, "\nAlignment conditions:\n");
  fprintf(out, "SVD values: (Euler convention: ZYZ)\n");
  //                  L   Aa +-Aa      R +-R   alpha  beta gamma     Q
  fprintf(out, " %-*s  %*s   %*s   %*s   %*s   %*s   %*s   %*s   %*s\n",
          baseInformation.output_format.string_l,
          centerCstring("Label", baseInformation.output_format.string_l),
          baseInformation.output_format.string_l,
          centerCstring("Aa", baseInformation.output_format.string_l),
          baseInformation.output_format.string_l, " ", // Error Dummy
          baseInformation.output_format.string_l,
          centerCstring("R", baseInformation.output_format.string_l),
          baseInformation.output_format.string_l, " ", // Error Dummy
          2 * baseInformation.output_format.number_m + 3,
          centerCstring("alpha / deg",
                        2 * baseInformation.output_format.number_m + 3),
          2 * baseInformation.output_format.number_m + 3,
          centerCstring("beta / deg",
                        2 * baseInformation.output_format.number_m + 3),
          2 * baseInformation.output_format.number_m + 3,
          centerCstring("gamma / deg",
                        2 * baseInformation.output_format.number_m + 3),
          baseInformation.output_format.string_l, "Q-factor");
  for (CurrSet = CurrMol.getHeadSet(); CurrSet;
       CurrSet = CurrSet->getNext(), ++index)
  {
    baseInformation.q_factors[index + iteration * NOS] = Q(index);
    fprintf(
        out,
        " %-*s: %*.*e   %*s   %*.*e   %*s   %*.*f   %*s   %*.*f   %*s   %*.*f  "
        " %*s   %*.*f\n",
        baseInformation.output_format.string_l,
        CurrSet->getLabel(baseInformation.output_format.string_l, true).c_str(),
        baseInformation.output_format.number_l,
        baseInformation.output_format.prec_m,
        EigenVals(STANDARD_AXIS_Z_, index),
        baseInformation.output_format.number_l, " ",
        baseInformation.output_format.number_l,
        baseInformation.output_format.prec_m,
        S2R(EigenVals(STANDARD_AXIS_X_, index),
            EigenVals(STANDARD_AXIS_Y_, index),
            EigenVals(STANDARD_AXIS_Z_, index)),
        baseInformation.output_format.number_l, " ",
        baseInformation.output_format.number_m,
        baseInformation.output_format.prec_m, -rad2deg(EulerAngles(0, index)),
        baseInformation.output_format.number_m, " ",
        baseInformation.output_format.number_m,
        baseInformation.output_format.prec_m, -rad2deg(EulerAngles(1, index)),
        baseInformation.output_format.number_m, " ",
        baseInformation.output_format.number_m,
        baseInformation.output_format.prec_m, -rad2deg(EulerAngles(2, index)),
        baseInformation.output_format.number_m, " ",
        baseInformation.output_format.number_l,
        baseInformation.output_format.prec_m,
        baseInformation.q_factors[index + iteration * NOS]);
  }

  if (!(flags.monteCarloBootstrapping && opt == StructureOptions::Optimized))
    return rdc_calc;

  fprintf(out, "\nMonte-Carlo values:\n");
  monte_carlo = CurrStruc->get_monte_carlo_output();

  Eigen::MatrixXd EVal_MC, EVal_MC_sigm;
  Eigen::MatrixXd EVec_MC, EVec_MC_sigm;
  Eigen::Vector3d Euler_MC, Euler_MC_sigm;
  Eigen::MatrixXd S_MC      = Eigen::MatrixXd::Zero(3, 3);
  Eigen::MatrixXd S_MC_sigm = Eigen::MatrixXd::Zero(3, 3);

  for (CurrSet = CurrMol.getHeadSet(), index = 0; CurrSet;
       CurrSet = CurrSet->getNext(), ++index)
  {
    S_MC(0, 0) = monte_carlo.Saupe_tensor(0, index);
    S_MC(1, 1) = monte_carlo.Saupe_tensor(1, index);
    S_MC(2, 2) = monte_carlo.Saupe_tensor(2, index);
    S_MC(1, 0) = S_MC(0, 1) = monte_carlo.Saupe_tensor(3, index);
    S_MC(2, 0) = S_MC(0, 2) = monte_carlo.Saupe_tensor(4, index);
    S_MC(2, 1) = S_MC(1, 2) = monte_carlo.Saupe_tensor(5, index);

    S_MC_sigm(0, 0) = monte_carlo.Saupe_tensor_sigm(0, index);
    S_MC_sigm(1, 1) = monte_carlo.Saupe_tensor_sigm(1, index);
    S_MC_sigm(2, 2) = monte_carlo.Saupe_tensor_sigm(2, index);
    S_MC_sigm(1, 0) = S_MC_sigm(0, 1) = monte_carlo.Saupe_tensor_sigm(3, index);
    S_MC_sigm(2, 0) = S_MC_sigm(0, 2) = monte_carlo.Saupe_tensor_sigm(4, index);
    S_MC_sigm(2, 1) = S_MC_sigm(1, 2) = monte_carlo.Saupe_tensor_sigm(5, index);

    estimateSaupeErrors(monte_carlo.Aligns.col(index), S_MC, EVec_MC, EVal_MC,
                        Euler_MC, monte_carlo.Aligns_sigm.col(index),
                        EVec_MC_sigm, EVal_MC_sigm, Euler_MC_sigm);

    fprintf(
        out,
        " %-*s: %*.*e +-%*.*e   %*.*e +-%*.*e   %*.*f +-%*.*f   %*.*f +-%*.*f  "
        " %*.*f +-%*.*f\n",
        baseInformation.output_format.string_l,
        CurrSet->getLabel(baseInformation.output_format.string_l, true).c_str(),
        baseInformation.output_format.number_l,
        baseInformation.output_format.prec_m, EVal_MC(STANDARD_AXIS_Z_, 0),
        baseInformation.output_format.number_l,
        baseInformation.output_format.prec_m, EVal_MC_sigm(STANDARD_AXIS_Z_, 0),
        baseInformation.output_format.number_l,
        baseInformation.output_format.prec_m,
        S2R(EVal_MC(STANDARD_AXIS_X_, 0), EVal_MC(STANDARD_AXIS_Y_, 0),
            EVal_MC(STANDARD_AXIS_Z_, 0)),
        baseInformation.output_format.number_l,
        baseInformation.output_format.prec_m,
        dS2dR(EVal_MC(STANDARD_AXIS_X_, 0), EVal_MC(STANDARD_AXIS_Y_, 0),
              EVal_MC(STANDARD_AXIS_Z_, 0), EVal_MC_sigm(STANDARD_AXIS_X_, 0),
              EVal_MC_sigm(STANDARD_AXIS_Y_, 0),
              EVal_MC_sigm(STANDARD_AXIS_Z_, 0)),
        baseInformation.output_format.number_m,
        baseInformation.output_format.prec_m, rad2deg(Euler_MC(0, 0)),
        baseInformation.output_format.number_m,
        baseInformation.output_format.prec_m, rad2deg(Euler_MC_sigm(0, 0)),
        baseInformation.output_format.number_m,
        baseInformation.output_format.prec_m, rad2deg(Euler_MC(1, 0)),
        baseInformation.output_format.number_m,
        baseInformation.output_format.prec_m, rad2deg(Euler_MC_sigm(0, 0)),
        baseInformation.output_format.number_m,
        baseInformation.output_format.prec_m, rad2deg(Euler_MC(2, 0)),
        baseInformation.output_format.number_m,
        baseInformation.output_format.prec_m, rad2deg(Euler_MC_sigm(0, 0)));
  }

  fprintf(out, "\n");
  fprintf(
      out,
      "%-*s  %*s   %*s   %*s   %*s   %*s \n", //   %*s   %*s   %*s   %*s  \n",
      baseInformation.output_format.string_l, "Saupe vector",
      2 * baseInformation.output_format.string_l + 3,
      centerCstring("Szz", 2 * baseInformation.output_format.string_l + 3),
      2 * baseInformation.output_format.string_l + 3,
      centerCstring("Sxx-Syy", 2 * baseInformation.output_format.string_l + 3),
      2 * baseInformation.output_format.string_l + 3,
      centerCstring("Sxy", 2 * baseInformation.output_format.string_l + 3),
      2 * baseInformation.output_format.string_l + 3,
      centerCstring("Sxz", 2 * baseInformation.output_format.string_l + 3),
      2 * baseInformation.output_format.string_l + 3,
      centerCstring("Syz", 2 * baseInformation.output_format.string_l + 3));

  S_MC      = monte_carlo.Aligns;
  S_MC_sigm = monte_carlo.Aligns_sigm;

  for (CurrSet = CurrMol.getHeadSet(), index = 0; CurrSet;
       CurrSet = CurrSet->getNext(), ++index)
  {
    fprintf(out, " %-*s:", baseInformation.output_format.string_l,
            CurrSet->getLabel(baseInformation.output_format.string_l, true)
                .c_str());
    for (int i = 0; i < SAUPE_ELEMENTS_; ++i)
      fprintf(out, " %*.*e +-%*.*e  ", baseInformation.output_format.number_l,
              baseInformation.output_format.prec_m, S_MC(i, index),
              baseInformation.output_format.number_l,
              baseInformation.output_format.prec_m, S_MC_sigm(i, index));
    fprintf(out, "\n");
  }
  fflush(out);
  return rdc_calc;
}

void
outputPolarAngles(FILE *out,
                  Structure *CurrStruc,
                  BasicInformation &baseInformation,
                  StructureOptions options)
{
  SphericalHarmonics *CurrY;
  bool initial = (options == StructureOptions::Initial);
  double theta, phi;
  theta = phi = .0;
  fprintf(out, " Polar angles of RDC vectors / deg: %s%s%s\n",
          Paper::cite_multiple_start(Paper::Meiler),
          Paper::cite_multiple_middle(Paper::Peti),
          Paper::cite_multiple_end(Paper::Lakomek));

  fprintf(out, "\n   %*s  | ", 2 * baseInformation.output_format.string_m + 1,
          "Nuclei list");
  fprintf(out, " %*s  ", 2 * baseInformation.output_format.string_l + 1,
          (initial ? "Initial values  " : "Optimized values "));
  fprintf(out, "\n");

  fprintf(out, "   %-*s %-*s  | ", baseInformation.output_format.string_m,
          "Nuc 1", baseInformation.output_format.string_m, "Nuc 2");
  fprintf(out, "  %*s %*s ", baseInformation.output_format.string_l, "theta",
          baseInformation.output_format.string_l, "phi");
  fprintf(out, "\n");

  for (CurrY = CurrStruc->getHeadYmatrix()->getHeadHarmonic(); CurrY;
       CurrY = CurrY->getNext())
  {
    theta = CurrY->getTheta(options);
    phi   = CurrY->getPhi(options);
    A2goodA(theta, phi);
    fprintf(out, "    %-*s %-*s | ", baseInformation.output_format.string_m,
            CurrY->getAtom1()->getIdentifier().c_str(),
            baseInformation.output_format.string_m,
            CurrY->getAtom2()->getIdentifier().c_str());
    fprintf(out, "  %*.*f %*.*f ", baseInformation.output_format.number_l,
            baseInformation.output_format.prec_m, rad2deg(theta),
            baseInformation.output_format.number_l,
            baseInformation.output_format.prec_m, rad2deg(phi));
    fprintf(out, "\n");
  }

  fprintf(out, "\n");
  fflush(out);
}

void
outputLevenbergMarquardt(FILE *out,
                         Structure *CurrStruc,
                         BasicInformation &baseInformation,
                         Flags &flags)
{
  SphericalHarmonics *CurrY;
  // TODO add commentaries on LM (only first step)
  if (baseInformation.numOfOptSteps == 1 && !flags.minimalOutput)
  {
    fprintf(out, " %s \n", std::string(header_width, '*').c_str());
    fprintf(out, " *  %-56s* \n", "Gradient: grad(0) -> initial gradient");
    fprintf(out, " *  %-56s* \n", "          grad(f) -> final gradient");
    fprintf(out, " *  %-56s* \n", "Hesse matrix estimation: Jac(trans)*Jac");
    fprintf(out, " *  %-56s* \n", "Change in parameters p: delta(p)");
    fprintf(out, " *  %-56s* \n", "epsilon^2 = [f(p_i)-f(p_i+1)]^2");
    fprintf(out, " *  %-56s* \n", "   with f(p): Y[MF] -> Y[VF]");
    fprintf(out, " %s \n", std::string(header_width, '*').c_str());
  }
  double *info;
  int stop_value_index;
  /******************************************************
   *                     Info map                       *
   *                                                    *
   * info[0]: number of iterations.                     *
   * info[1]: ||epsilon||(2) of p(initial).             *
   * info[2]: ||g||(inf) of p(initial).                 *
   * info[3]: ||epsilon||(2) of p(final).               *
   * info[4]: ||g||(inf) of p(final).                   *
   * info[5]: ||delta(p)||(2) at p(final).              *
   * info[6]: max(diag[hessian]) of p(final).           *
   * info[7]: reason for termination.                   *
   *          - 1: low delta(p)                         *
   *          - 2: low g(p)                             *
   *          - 3: low epsilon^2(p)                     *
   *          - 4: max ierations                        *
   *          - 5: low initial g(p)                     *
   * info[8]: respective criterion for info[7].         *
   *                                                    *
   *                                                    *
   *                                                    *
   ******************************************************/

  fprintf(out,
          " Levenberg Marquardt information on Rotation[Y(MF) -> Y(VF)]: \n");
  fprintf(out, "\n");
  fprintf(out, "   %*s  | ", 2 * baseInformation.output_format.string_m + 1,
          "Nuclei list");
  fprintf(out, "%*s %*s | %*s %*s | %*s \n",
          baseInformation.output_format.string_m, " ",
          2 * baseInformation.output_format.number_l + 1,
          centerCstring("gradients",
                        2 * baseInformation.output_format.number_l + 1),
          baseInformation.output_format.number_l,
          centerCstring("hessian", baseInformation.output_format.number_l),
          baseInformation.output_format.number_l,
          centerCstring("norm", baseInformation.output_format.number_l),
          3 * baseInformation.output_format.string_l + 2,
          centerCstring("stop reason",
                        3 * baseInformation.output_format.string_l + 2));

  fprintf(out, "   %-*s %-*s  | ", baseInformation.output_format.string_m,
          "Nuc 1", baseInformation.output_format.string_m, "Nuc 2");
  fprintf(out, "%*s %*s %*s | ", baseInformation.output_format.string_m,
          "iter.", baseInformation.output_format.number_l, "grad(0)",
          baseInformation.output_format.number_l, "grad(f)");

  fprintf(out, "%*s %*s | ", baseInformation.output_format.number_l,
          centerCstring("max(diag)", baseInformation.output_format.number_l),
          baseInformation.output_format.number_l,
          centerCstring("delta(p)", baseInformation.output_format.number_l));

  fprintf(out, "%*s %*s %*s", baseInformation.output_format.string_l,
          "criterion", baseInformation.output_format.string_l, "value",
          baseInformation.output_format.number_l, "threshhold");

  fprintf(out, "\n");

  for (CurrY = CurrStruc->getHeadYmatrix()->getHeadHarmonic(); CurrY;
       CurrY = CurrY->getNext())
  {
    fprintf(out, "    %-*s %-*s | ", baseInformation.output_format.string_m,
            CurrY->getAtom1()->getIdentifier().c_str(),
            baseInformation.output_format.string_m,
            CurrY->getAtom2()->getIdentifier().c_str());
    info = CurrY->getLMinfo();
    switch (static_cast<int>(info[7]))
    {
      case 1:
        stop_value_index = 5;
        break;
      case 2:
        stop_value_index = 4;
        break;
      case 3:
        stop_value_index = 3;
        break;
      case 4:
        stop_value_index = 0;
        break;
      case 5:
        stop_value_index = 2;
        break;
      default:
        stop_value_index = -1;
        break;
    }
    fprintf(out, "%*d %*.*e %*.*e | %*.*e %*.*e | %*s ",
            baseInformation.output_format.string_m, static_cast<int>(info[0]),
            baseInformation.output_format.number_l,
            baseInformation.output_format.prec_l, info[2],
            baseInformation.output_format.number_l,
            baseInformation.output_format.prec_l, info[4],
            baseInformation.output_format.number_l,
            baseInformation.output_format.prec_l, info[6],
            baseInformation.output_format.number_l,
            baseInformation.output_format.prec_l, info[5],
            baseInformation.output_format.string_l,
            LM_stop(static_cast<int>(info[7])));
    if (stop_value_index)
    {
      double reason = (stop_value_index == 3 ?
                           info[stop_value_index] * info[stop_value_index] :
                           info[stop_value_index]);
      fprintf(out, "%*.*e %*.*e", baseInformation.output_format.number_l,
              baseInformation.output_format.prec_l, reason,
              baseInformation.output_format.number_l,
              baseInformation.output_format.prec_l, info[8]);
    }
    else
      fprintf(out, "%*.d %*.d", baseInformation.output_format.number_l,
              static_cast<int>(info[stop_value_index]),
              baseInformation.output_format.number_l,
              static_cast<int>(info[8]));

    fprintf(out, "\n");
  }
  fprintf(out, "\n");
  info = NULL;
  fflush(out);
}

void
outputRedundantInternalCoordinates(FILE *out,
                                   RedundantInternals &RedOpt,
                                   BasicInformation &baseInformation,
                                   Flags &flags)
{
  double damping;
  damping = RedOpt.getDamping();
  fprintf(out, " Information on redundant internal coordinate optimizer: \n");
  fprintf(out, "\n");
  fprintf(out, " %-*s: %*u (%*u) \n", baseInformation.output_format.string_xxl,
          "Iteration (max. iteration)", baseInformation.output_format.number_s,
          RedOpt.getCycles(), baseInformation.output_format.number_s,
          RedOpt.getMaxCycles());
  fprintf(out, " %-*s: %*.*e \n", baseInformation.output_format.string_xxl,
          "Overall damping", baseInformation.output_format.number_l,
          baseInformation.output_format.prec_l, RedOpt.getDamping());
  fprintf(
      out, " %-*s: %*.*e (%*.*e) \n", baseInformation.output_format.string_xxl,
      "Bond lengths (effective)", baseInformation.output_format.number_l,
      baseInformation.output_format.prec_l,
      baseInformation.static_redundants_weighting[BOND_REDUNDANTS_],
      baseInformation.output_format.number_l,
      baseInformation.output_format.prec_l,
      baseInformation.static_redundants_weighting[BOND_REDUNDANTS_] * damping);
  fprintf(
      out, " %-*s: %*.*e (%*.*e) \n", baseInformation.output_format.string_xxl,
      "Bond angles (effective)", baseInformation.output_format.number_l,
      baseInformation.output_format.prec_l,
      baseInformation.static_redundants_weighting[ANGLE_REDUNDANTS_],
      baseInformation.output_format.number_l,
      baseInformation.output_format.prec_l,
      baseInformation.static_redundants_weighting[ANGLE_REDUNDANTS_] * damping);
  if (flags.torsions_4_redundants)
  {
    fprintf(out, " %-*s: %*.*e (%*.*e) \n",
            baseInformation.output_format.string_xxl,
            "Torsion angles (effective)",
            baseInformation.output_format.number_l,
            baseInformation.output_format.prec_l,
            baseInformation.static_redundants_weighting[TORSION_REDUNDANTS_],
            baseInformation.output_format.number_l,
            baseInformation.output_format.prec_l,
            baseInformation.static_redundants_weighting[TORSION_REDUNDANTS_] *
                damping);
  }
  fprintf(
      out, " %-*s: %*.*e (%*.*e) \n", baseInformation.output_format.string_xxl,
      "RDC angles (effective)", baseInformation.output_format.number_l,
      baseInformation.output_format.prec_l,
      baseInformation.static_redundants_weighting[RDC_REDUNDANTS_],
      baseInformation.output_format.number_l,
      baseInformation.output_format.prec_l,
      baseInformation.static_redundants_weighting[RDC_REDUNDANTS_] * damping);
  fprintf(out, " %-*s: %*.*e (%*.*e) \n",
          baseInformation.output_format.string_xxl,
          "Chiral volumes (effective)", baseInformation.output_format.number_l,
          baseInformation.output_format.prec_l,
          baseInformation.static_redundants_weighting[PLANAR_REDUNDANTS_],
          baseInformation.output_format.number_l,
          baseInformation.output_format.prec_l,
          baseInformation.static_redundants_weighting[PLANAR_REDUNDANTS_] *
              damping);
  fprintf(out, " %-*s: %*.*e \n", baseInformation.output_format.string_xxl,
          "Distances", baseInformation.output_format.number_l,
          baseInformation.output_format.prec_l,
          baseInformation.static_redundants_weighting[DISTANCE_REDUNDANTS_]);
  fprintf(out, " %-*s: %*.*e \n", baseInformation.output_format.string_xxl,
          "Delta x (final)", baseInformation.output_format.number_l,
          baseInformation.output_format.prec_l, RedOpt.getRMSDx());
  fprintf(out, " %-*s: %*.*e \n", baseInformation.output_format.string_xxl,
          "Delta q (final)", baseInformation.output_format.number_l,
          baseInformation.output_format.prec_l, RedOpt.getRMSDq());
  fprintf(out, " %-*s: %*.*e (%*.*e) \n",
          baseInformation.output_format.string_xxl, "Stop crit. (threshold)",
          baseInformation.output_format.number_l,
          baseInformation.output_format.prec_l, RedOpt.getStopCrit(),
          baseInformation.output_format.number_l,
          baseInformation.output_format.prec_l,
          baseInformation.limits.redundants_convergence);
  fprintf(out, "\n");
  fflush(out);
}

void
outputChiralVolume(FILE *out,
                   Structure *CurrStruc,
                   BasicInformation &baseInformation,
                   Flags &flags,
                   StructureOptions options)
{
  Atom *CurrAtom;
  fprintf(out, " Chiral volumes of the (pro-)chiral centers during "
               "optimization / (1e-10 m)^3\n");

  for (CurrAtom = CurrStruc->getHeadAtom(); CurrAtom;
       CurrAtom = CurrAtom->getNext())
  {
    if (!CurrAtom->hasChiralVolume())
      continue;
    fprintf(out, "%10s", CurrAtom->getIdentifier().c_str());
  }
  fprintf(out, "\n");
  int center = 0;
  for (CurrAtom = CurrStruc->getHeadAtom(); CurrAtom;
       CurrAtom = CurrAtom->getNext())
  {
    if (!CurrAtom->hasChiralVolume())
      continue;
    baseInformation.chiral_volumes[baseInformation.numOfOptSteps *
                                       baseInformation.NumberOfChiralCenters +
                                   center] =
        CurrAtom->getChiralVolume(options, flags.normChiralVolume);
    fprintf(out, "%10.4f",
            baseInformation
                .chiral_volumes[baseInformation.numOfOptSteps *
                                    baseInformation.NumberOfChiralCenters +
                                center++]);
  }

  fprintf(out, "\n");
  fflush(out);
}

void
outputStructure2trj(BasicInformation &baseInformation,
                    Flags &flags,
                    Molecule &CurrMol,
                    Structure *CurrStruc,
                    RedundantInternals &RedundantOptimizer,
                    int iteration)
{
  // TODO zusammenfassen von logischen Blöcken (auch .ali)
  // TODO Zitate setzen (auch .ali)
  // .ali todos:
  // TODO Chi^2 zu Q-faktor dazu packen
  // TODO MC Fehler -> alle rechts ausrichten.
  // TODO eta klammer verschieben
  // TODO Euler Winkel beschreiben und ausrichten

  StructureOptions options =
      (iteration ? StructureOptions::Optimized : StructureOptions::Initial);

  /*
   * BLOCK #1:
   *  Header
   */
  if (iteration)
  {
    fprintf(baseInformation.trjFile, "%s\n* Iteration # ",
            std::string(60, '*').c_str());
    fprintf(baseInformation.trjFile, "%-45d", iteration);
    fprintf(baseInformation.trjFile, "*\n* Structure: ");
  }
  else
  {
    fprintf(baseInformation.trjFile, "%s\n* Initial state",
            std::string(60, '*').c_str());
    fprintf(baseInformation.trjFile,
            "%s*\n* Structure: ", std::string(44, ' ').c_str());
  }
  fprintf(baseInformation.trjFile, "%-46s", CurrStruc->getLabel().c_str());
  fprintf(baseInformation.trjFile, "*\n%s\n\n", std::string(60, '*').c_str());

  /*
   * BLOCK #2:
   *  SECONDA (only Iteration 0 and n)
   */
  if (iteration == 0 || baseInformation.SECONDA_output_final.collected)
  {
    outputTrjSECONDA(baseInformation.trjFile, CurrStruc, baseInformation);
  }

  /*
   * BLOCK #3:
   *  Polar angles, Levenberg Marquardt and Monte Carlo
   */
  outputPolarAngles(baseInformation.trjFile, CurrStruc, baseInformation,
                    options);
  if (iteration)
  {
    outputLevenbergMarquardt(baseInformation.trjFile, CurrStruc,
                             baseInformation, flags);
    if (flags.monteCarloBootstrapping)
      add_MonteCarlo_angles(baseInformation.trjFile, CurrStruc, baseInformation,
                            flags, CurrStruc->get_monte_carlo_output());
  }
  if (iteration && flags.useRedundants)
    outputRedundantInternalCoordinates(
        baseInformation.trjFile, RedundantOptimizer, baseInformation, flags);

  /*
   * BLOCK #4:
   *  Coordinates with Transformation matrices and inertia tensor
   */
  add_structure(baseInformation.trjFile, CurrStruc, baseInformation, flags,
                options);
  if (iteration == 0)
    outputInitialTransformation(baseInformation.trjFile, CurrMol,
                                baseInformation);
  else
    outputInertiaTensor(baseInformation.trjFile, CurrStruc, baseInformation,
                        options);

  add_cosine_matrix(baseInformation.trjFile, CurrStruc, baseInformation, flags,
                    options);

  /*
   * Block #6:
   *  Chiral volumes
   */
  outputChiralVolume(baseInformation.trjFile, CurrStruc, baseInformation, flags,
                     options);

  /*
   * Block #7:
   *  Alignment conditions
   */
  Eigen::MatrixXd rdc_calc =
      outputAlignments(baseInformation.trjFile, CurrStruc, CurrMol, iteration,
                       baseInformation, flags, options);


  if (iteration)
  {
    Eigen::MatrixXd D = CurrStruc->getRDCmatrix(rdcMatrixOptions::Unscaled);
    add_Srdc(baseInformation.trjFile, CurrStruc, baseInformation, flags, D,
             rdc_calc);
    if (flags.monteCarloBootstrapping)
      add_stop_criteria(baseInformation.trjFile, baseInformation, flags);

    fprintf(baseInformation.trjFile, "\n%s\n* Finished iteration # ",
            std::string(60, '*').c_str());
    fprintf(baseInformation.trjFile, "%-36d", iteration);
    fprintf(baseInformation.trjFile, "*\n%s", std::string(60, '*').c_str());
  }
  else
  {
    fprintf(baseInformation.trjFile, "\n\n%s\n*%-58.58s*\n%s",
            std::string(60, '*').c_str(), "End of initial state",
            std::string(60, '*').c_str());
    // fprintf( baseInformation.trjFile, "",  );
  }
  fprintf(baseInformation.trjFile, "\n\n");
  fflush(baseInformation.trjFile);
}

void
prepare_trj_file(Molecule &CurrMol,
                 RedundantInternals &RedundantsOptimizer,
                 BasicInformation &baseInformation,
                 Flags &flags)
{
  // START of output
  FILE *output;
  if (baseInformation.trjFile)
    output = baseInformation.trjFile;
  else
  {
    baseInformation.trjFile = fopen(baseInformation.trjFileName.c_str(), "w");
    output                  = baseInformation.trjFile;
  }
  fprintf(output, " %s \n", std::string(header_width, '=').c_str());
  fprintf(output, " =%s= \n", std::string(header_width - 2, ' ').c_str());
  frame_string(output, "TITANIA", 25);
  frame_string(output, "Trajectory output (.trj)", 18);
  fprintf(output, " =%s= ", std::string(header_width - 2, ' ').c_str());

  // General information part
  information_header(output, "General information");
  add_contributors(output);

  // TODO add_TITANIA_citation (output);

  // Technical information part
  information_header(output, "Technical information");
  add_compilation_information(output, baseInformation);

  // Input information part
  information_header(output, "Input information");
  add_input_file_information(output, baseInformation);

  // Run information part
  information_header(output, "Run information");
  add_run_information(output, baseInformation);

  information_header(output, "Metrics");
  add_metrics(output, baseInformation);
  add_redundant_metrics(output, baseInformation, RedundantsOptimizer);

  information_header(output, "RDC Lists");
  add_rdc_sets(output, CurrMol, baseInformation, flags);
  fprintf(output, "\n");
  fflush(output);
  output = NULL;
}
