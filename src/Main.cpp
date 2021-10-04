
#include <Atom.hpp>
#include <Bond.hpp>
#include <DIDC.hpp>
#include <Declarations.hpp>
#include <Eckart.hpp>
#include <Helper.hpp>
#include <IndependencyAnalysis.hpp>
#include <LinAlg.hpp>
#include <Main.hpp>
#include <Molecule.hpp>
#include <MonteCarloStatistic.hpp>
#include <Output.hpp>
#include <Parser/Parser.hpp>
#include <Properties.hpp>
#include <RDCdata.hpp>
#include <RDCset.hpp>
#include <RedundantInternals.hpp>
#include <SCRM.hpp>
#include <SphericalHarmonics.hpp>
#include <SphericalHarmonicsMatrix.hpp>
#include <Structure.hpp>
#include <StructureSimulator.hpp>
#include <iostream>
#include <omp.h>
#include <phobos.h>
#include <unistd.h>
#include <version.hpp>

#ifdef USE_FRAGMENTS_
#include <Fragment.hpp>
#endif

#ifdef USE_CUDA
#include <cuda_helper.hpp>
#endif

#define MC_FINAL_RUN_CONVERGENCE 0.01
#define MC_FINAL_RUN_MAX_STEPS   20000

/*******************************************************
 *    loading is a string array, which displays the    *
 * respective string on the terminal. The intention is *
 *    to show the user, that TITANIA is still alive.   *
 ******************************************************/

std::string loading[4] = {"[\\]", "[|]", "[/]", "[-]"};

/*******************************************************
 *                                                     *
 * void setupMolecule                                  *
 *                                                     *
 *     TITANIA uses some general information on the    *
 *   molecule. These information are generated here.   *
 * Additionally the molecule is prepared by coordinate *
 *                  transformations.                   *
 *                                                     *
 *******************************************************/

void
setupMolecule(Molecule &CurrMol, Flags &flags)
{
  CurrMol.determineMolecularMass();
  // Remove the following two lines to not tranform the
  // structure to the PAS of the intertia tensor.
  // This is needed to reproduce the figure showing the
  // distribution of r(C3-H3)_{MC} in three different
  // frames of reference!
  CurrMol.getHeadStruc()->Shift2CenterOfMass(StructureOptions::Initial);
  CurrMol.getHeadStruc()->Rotate2InertiaPAS(StructureOptions::Initial);
  setupBonds(CurrMol, flags);
  CurrMol.countLongRangeRDCs();
}

/*******************************************************
 *                                                     *
 * int ReadInput                                       *
 *                                                     *
 *     TITANIA uses class InputFile to fully parse,    *
 * interpret and load the user defined input file. All *
 * information are directly stored in CurrMol. If any  *
 *  error occures TITANIA will automatically terminate *
 *                        safely.                      *
 *                                                     *
 * Wiki page exists.                                   *
 *                                                     *
 *******************************************************/

int
ReadInput(Molecule &CurrMol,
          InputFile &reader,
          BasicInformation &baseInformation,
          Flags &flags,
          std::fstream &input)
{
  // Collect all data from the pased data
  baseInformation.state = reader.checkFile();

  // Check and load all information from the data
  if (baseInformation.state == GOOD_STATE)
    baseInformation.state = CurrMol.loadInput(reader, baseInformation, flags);

  input.close(), input.clear();

  // Check if the parsing was successfull
  switch (baseInformation.state)
  {
    case GOOD_STATE:
      break;
    case ERROR_CORRUPTED_STRUCTURE_INPUT:
      killTITANIA(baseInformation, flags);
      return baseInformation.state;
    case ERROR_UNKNOWN_KEYWORDS:
      killTITANIA(baseInformation, flags);
      return baseInformation.state;
    case ERROR_CORRUPTED_RDC_INPUT:
      killTITANIA(baseInformation, flags);
      return baseInformation.state;
    case ERROR_IDENTIFIER_MISSING:
      killTITANIA(baseInformation, flags);
      return baseInformation.state;
    case ERROR_FILE_DOES_NOT_EXIST:
      killTITANIA(baseInformation, flags);
      return baseInformation.state;
    case TITANIA_RDC_LINKING_ERROR:
      killTITANIA(baseInformation, flags);
      return baseInformation.state;
    default:
      break;
  }
  if (flags.plotMonteCarlo)
  {
    flags.plotMonteCarlo = flags.monteCarloOutput =
        flags.monteCarloBootstrapping;
  }
  if (flags.plotRDCdynamics)
    flags.plotRDCdynamics = (flags.skipSCRM ? false : true);
  return baseInformation.state;
}

/*******************************************************
 *                                                     *
 * int initializeProgram                               *
 *                                                     *
 *    TITANIA first loads all standard settings (eg    *
 *    limits). Afterwords terminal arguments will be   *
 *  parsed. grumble will check if the expected values  *
 *                 are set correctly.                  *
 *                                                     *
 *******************************************************/

int
initializeProgram(BasicInformation &baseInformation,
                  Flags &flags,
                  int argc,
                  char *argv[])
{
  /* Load standard values for baseInformation and flags */
  baseSettings(baseInformation, flags);

  /* Check the arguments with which the program was started. The return
     value one denotes, that -h (help) was called, and terminates the
     program right after a short documentation was printed on the terminal. */

  switch (initializeArguments(argc, argv, baseInformation, flags))
  {
    case OUTPUT_HELP:
      return OUTPUT_HELP;
    case OUTPUT_VERSION:
      return OUTPUT_VERSION;
    case OUTPUT_BASE_SETTINGS: {
      // First reset everything to be sure the user gets the base settings.
      baseSettings(baseInformation, flags);
      printKeywords(baseInformation, flags);
      return OUTPUT_BASE_SETTINGS;
    }
    case TITANIA_MULTI_JOB:
      return TITANIA_MULTI_JOB;
    default:
      break;
  }
  return grumble(baseInformation, flags);
}

/*******************************************************
 *                                                     *
 * int checkArguments                                  *
 *                                                     *
 *       TITANIA allows some terminal arguments.       *
 *      checkArguments parses and interprets them.     *
 *                                                     *
 *******************************************************/

int
checkArguments(const char *argv)
{
  if (argv[0] == '-')
  {
    switch (argv[1])
    {
      case 'i':
        return 1; // input
      case 'o':
        return 2; // output
      case 'h':
        return 3; // help
      case 'a':
        return 4; // append to existing
      case 'c':
        return 5; // change number of used threads
      case 'd':
        return 6; // activate debug (not implemented right now)
      case 'r':
        return 7; // redundants
      case 'v':
        return 8; // version
      case 'b':
        return 9; // baseSettings
      case 's':
        return 10;
      case 't':
        return 11; // temporary flag! Only turn on for debugging reasons
    }
  }
  return 0;
}

/*******************************************************
 *                                                     *
 * void outputHelp                                     *
 *                                                     *
 *  This will call a minimal documentation on how to   *
 *      use TITANIA. Keep this simple! The full        *
 *  documentation should be located on TITANIA Wiki.   *
 *                                                     *
 *******************************************************/

void
outputHelp()
{
  *std::cin.tie()
      << "The help function of TITANIA was called.\n\n"
         "The following options are available:\n\n"
         "-i\tDefine an input file for the analysis\n"
         "-o\tDefine an output file for the analysis (default: "
         "<inputfilename.out>)\n"
         "-a\tAppend the output to the defined output file (syntax like -o)\n"
         "-c\tDefine the number of cpus used for the calculation (default: "
         "all)\n"
         "-r\tSuppress the use of redundant internal coordinates (default: "
         "1)\n\n"
         "You will find a wiki page on gitlab which might help you with some\n"
         "additional information.\n"
         "If necessary contact froth@thielelab.de\n";
}

/*******************************************************
 *                                                     *
 * int initialzeArguments                              *
 *                                                     *
 *  All arguments will be parsed, interpreted and the  *
 *        Flags/BasicInformation will be saved.        *
 *                                                     *
 *******************************************************/

int
initializeArguments(int argc,
                    char *argv[],
                    BasicInformation &baseInformation,
                    Flags &flags)
{
  int argt, i, j, k, num_of_arg;
  bool multi_input = false, only_input = true;
  // If there are only 2 arguments (TITANIA <argv[1]>) we
  // can assume argv[1] to be the input file (as long as
  // it does not start with "-"
  if (argc == 2 && argv[1][0] != '-')
  {
    baseInformation.inputFileName = argv[1];
    checkWorkingDir(baseInformation, flags);
    return 0;
  }
  else
    only_input = false;

  // Check the number of "-" starting args.

  for (i = 1, num_of_arg = 0; i < argc; ++i)
  {
    if (checkArguments(argv[i]))
      ++num_of_arg;
  }

  // Save the position[0] and type[1] for each arg[num_of_arg]
  int arg[num_of_arg][2];

  for (i = 1, j = 0; i < argc && j < num_of_arg; ++i)
  {
    argt = checkArguments(argv[i]);
    if (!argt)
      continue;
    arg[j][0] = i;
    arg[j][1] = argt;

    switch (argt)
    {
      case 3:
        outputHelp();
        return OUTPUT_HELP;
      case 4:
        flags.appendOutput = true;
        break;
      case 6:
        flags.debug = true;
        break;
      case 8:
        return OUTPUT_VERSION;
      case 9:
        return OUTPUT_BASE_SETTINGS;
      default:
        break;
    }
    ++j;
  }

  for (i = 0; i < num_of_arg && (arg[i][0] != 1); ++i)
    ;

  multi_input =
      ((((i + 1) < num_of_arg) && ((arg[i + 1][0] - arg[i][0]) > 2)) ||
       (only_input && argc > 2 && num_of_arg == 1));

  if (multi_input)
  {
    std::string baseJob = argv[0];
    for (j = 0; j < num_of_arg; ++j)
    {
      if (arg[j][0] == 1)
        continue;
      baseJob += " ";
      baseJob += argv[arg[j][0]];
      baseJob += " ";
      baseJob += argv[arg[j][0] + 1];
    }
    int next_start = ((i < (num_of_arg - 1)) ? arg[i + 1][0] : argc);
    for (j = 1; j < (next_start - arg[i][0]); ++j)
    {
      std::string individualJob = baseJob + " -i " + argv[arg[i][0] + j] +
                                  " -s " + argv[arg[i][0] + j] + ".multiJob";
      *std::cin.tie() << "STARTING: " << individualJob << std::endl;
      k = system(individualJob.c_str());
      if (k)
        *std::cin.tie() << k << std::endl;
    }
    return TITANIA_MULTI_JOB;
  }

  // Loop all argument types.

  for (j = 0; j < num_of_arg; ++j)
  {
    // Jump to argument type arg[j][1] at position arg[j][0]
    // and interpret the following argument at i = arg[j][0] + 1
    for (i = arg[j][0] + 1; i < argc; ++i)
    {
      // Check that we are still in argv range
      // and did not enter next argument
      if ((j + 1) < num_of_arg && i == arg[j + 1][0])
        break;
      // Check the current type
      switch (arg[j][1])
      {
        case 1:
          if (baseInformation.numberOfInputs && flags.printWarnings)
          {
            std::cerr << "WARNING:\tMore than one input file is provided. "
                         "TITANIA is currently supporting just one input.\n"
                      << "\t\tThe program will run with the first file ("
                      << baseInformation.inputFileName << ") and\n"
                      << "\t\tignore the additional one (" << argv[i] << ")."
                      << std::endl;
            break;
          }
          baseInformation.inputFileName = argv[i];
          checkWorkingDir(baseInformation, flags);
          break;
        case 2:
          if (!flags.noOutput && flags.printWarnings)
          {
            std::cerr << "WARNING:\tMore than one output file was requested. "
                         "TITANIA can just handle one output.\n"
                      << "\t\tThe program will run with the first file ("
                      << baseInformation.outputFileName << ") and\n"
                      << "\t\tignore the additional one (" << argv[i] << ")."
                      << std::endl;
            break;
          }
          baseInformation.outputFileName = argv[i];
          flags.noOutput                 = false;
          break;
        case 4:
          if (!flags.noOutput)
          {
            if (argv[arg[j][0] + 1] && flags.printWarnings)
            {
              std::cerr
                  << "WARNING:\tMore than one output file was requested in "
                     "addition the argument \"-a\".\n"
                  << "\t\tTITANIA can just handle one output. "
                  << "The program will run with the first file ("
                  << baseInformation.outputFileName << ") and\n"
                  << "\t\tignore the additional one (" << argv[i]
                  << "). The new output will be appended to prevent data loss!"
                  << std::endl;
            }
            flags.appendOutput = true;
            break;
          }
          baseInformation.outputFileName = argv[i];
          flags.noOutput                 = false;
          flags.appendOutput             = true;
          break;
        case 5:
          if ((arg[j][0] + 1) == argc || argv[i][0] == 'm')
            baseInformation.numOfThreads = baseInformation.maxThreads;
          else
          {
            int nCores = atoi(argv[i]);
            if (nCores > baseInformation.maxThreads && flags.printWarnings)
            {
              std::cerr << "WARNING:\t" << nCores
                        << " threads were requested but only "
                        << baseInformation.maxThreads << " are available.";
              std::cerr << "\n\t\tUsing " << (baseInformation.maxThreads / 2)
                        << " (number of physical cores) for calculation.\n";
              baseInformation.numOfThreads = (baseInformation.maxThreads / 2);
            }
            baseInformation.numOfThreads = nCores;
          }
          break;
        case 7:
          flags.useRedundants = (atoi(argv[i]) > 0 ? true : false);
          break;
        case 10:
          flags.silent                = true;
          baseInformation.comFileName = (argv[i]);
          shutup(baseInformation);
          break;
        case 11:
          std::cout
              << "Temporary information argument does not access any variable."
              << std::endl;
          baseInformation.tmp_information.double_ = atof(argv[i]);
          break;
      }
    }
  }
  return 0;
}

/*******************************************************
 *                                                     *
 * int grumble                                         *
 *                                                     *
 *  Check if TITANIA recieved all information and was  *
 * able to interpret them correctly. If any point was  *
 *   set wrong TITANIA will automatically terminate    *
 *                       safely.                       *
 *                                                     *
 *******************************************************/

int
grumble(BasicInformation &bI, Flags &f)
{
  if (f.noInput)
  {
    bI.errorKey = "ERROR:\t\tNo input file was specified.\n";
    bI.errorKey += "\t\tUsage: TITANIA -i <filename> or TITANIA <filename>\n";
    bI.errorKey += "\t\tFor more information run TITANIA -h.\n";
    return ERROR_NO_INPUT_FILE;
  }
  if (f.noOutput)
  {
    bI.outputFileName = bI.inputFileName + ".out";
  }
  else if (!(bI.inputFileName.compare(bI.outputFileName)))
  {
    bI.errorKey = "ERROR:\t\tName of input and output file is the same.\n";
    bI.errorKey +=
        "\t\tDefine the output file properly or use default setting.\n";
    return ERROR_INPUT_EQ_OUTPUT;
  }

  if (bI.numOfThreads)
    omp_set_num_threads(bI.numOfThreads);
  else
    omp_set_num_threads(omp_get_num_procs() / 2);

  return GOOD_STATE;
}

/*******************************************************
 *                                                     *
 * void handleMemory                                   *
 *                                                     *
 *       TITANIA uses an own Levenberg-Marquardt       *
 * implementation which is run on several threads. To  *
 *  supress massive allocation and freeing of memory   *
 * space a large chunk is allocated at the beginning.  *
 * In every scenario this memory space is will be set  *
 *                      free again.                    *
 *                                                     *
 *******************************************************/

void
handleMemory(BasicInformation &bI)
{
  int i;
  *std::cin.tie()
      << "\t\tAllocating memory for Levenberg-Marquardt optimizations...\n";
  // PHOBOS_WORKSPACE(n_parameter, n_measurements).

  bI.phobos_worksize_Sphericals = PHOBOS_WORKSIZE(
      NUMBER_OF_SPHERICAL_COORDINATES, NUMBER_OF_SPHERICAL_FUNCTIONS);
  bI.phobos_worksize_Eckart =
      PHOBOS_WORKSIZE(NUMBER_OF_ECKART_ANGLES, NUMBER_OF_ECKART_FUNTIONS);

  bI.phobos_workspace_Sphericals =
      (double **) malloc(bI.numOfThreadsUsed * sizeof(double *));
  bI.phobos_workspace_Eckart =
      (double *) malloc(bI.phobos_worksize_Eckart * sizeof(double *));
  bI.q_factors = (double *) malloc(
      bI.NumberOfSets *
      (bI.limits.max_titania_iterations + bI.overOptimization + 1) *
      sizeof(double));

  for (i = 0; i < bI.numOfThreadsUsed; ++i)
  {
    bI.phobos_workspace_Sphericals[i] = NULL;
  }
}

/*******************************************************
 *                                                     *
 * void freeMemory                                     *
 *                                                     *
 *  Counter part to handle memory. Since this function *
 *  might be called withing any parallel part, one has *
 *               to think about tasks...               *
 *                                                     *
 *******************************************************/

void
freeMemory(BasicInformation &bI)
{
// Make sure that all workspace is unused
// and only freed once.
#pragma omp taskwait
#pragma omp critical
  {
    int i;
    *std::cin.tie() << "\t\tCleaning reserved memory space...\n";
    if (bI.alignments)
    {
      for (i = 0; i < (bI.numberOfAlignments + bI.numberOfFills); ++i)
        free(bI.alignments[i]);
    }

    for (i = 0; i < bI.numOfThreadsUsed && bI.phobos_workspace_Sphericals; ++i)
    {
      free(bI.phobos_workspace_Sphericals[i]);
    }

    if (bI.alignments)
      free(bI.alignments);
    if (bI.phobos_workspace_Sphericals)
      free(bI.phobos_workspace_Sphericals);
    if (bI.phobos_workspace_Eckart)
      free(bI.phobos_workspace_Eckart);
    if (bI.q_factors)
      free(bI.q_factors);
    if (bI.chiral_volumes)
      free(bI.chiral_volumes);
  }
}

/*******************************************************
 *                                                     *
 * void handleOutput                                   *
 *                                                     *
 *     Opens the output, coordinate and debug file.    *
 *  The function will take care of how to truncate or  *
 *              append content of the output.          *
 *                                                     *
 *******************************************************/

int
handleOutput(BasicInformation &bI, Flags &f)
{
  if (f.appendOutput)
  {
    *std::cin.tie() << "\nInformation:\tAppending output to the file "
                    << bI.outputFileName << std::endl;
    if (bI.output)
    {
      bI.state = ERROR_OUTPUT_OPEN;
      return ERROR_OUTPUT_OPEN;
    }
    bI.output = fopen(bI.outputFileName.c_str(), "a");
  }
  else
  {
    if (bI.output)
    {
      bI.state = ERROR_OUTPUT_OPEN;
      return ERROR_OUTPUT_OPEN;
    }
    bI.output = fopen(bI.outputFileName.c_str(), "w");
  }
  if (f.debug) /* File containing lots of runtime information. */
  {
    bI.debugFileName = bI.inputFileName + ".debug";
    bI.debugFile.open(bI.debugFileName,
                      std::ios::binary | std::ios::out | std::ios::trunc);
    //       bI.debugFile << "Debug file was created.\n\n";
  }
  // Prepare the xyz output file
  bI.xyzFileName = bI.outputFileName + ".xyz";
  bI.xyzFile.open(bI.xyzFileName,
                  std::ios::binary | std::ios::out | std::ios::trunc);
  bI.trjFileName = bI.outputFileName + ".trj";
  if (bI.trjFile)
  {
    bI.state = ERROR_TRJ_OPEN;
    return ERROR_TRJ_OPEN;
  }
  bI.trjFile = fopen(bI.trjFileName.c_str(), "w");
  return GOOD_STATE;
}

/*******************************************************
 *                                                     *
 * void killTITANIA                                    *
 *                                                     *
 *       If anything goes wrong on TITANIA runtime     *
 * killTITANIA will be called. This function tries to  *
 *  give the user information on what went wrong and   *
 *         then free the allocated memory.             *
 *                                                     *
 *******************************************************/

void
killTITANIA(BasicInformation &bI, Flags &f)
{
  switch (bI.state)
  {
    case 1:
      break;
    case ERROR_UNKNOWN_KEYWORDS:
      bI.errorKey += "\nERROR:\tTo many keywords were undefined. To ensure a "
                     "clean program execution the current\n\t";
      bI.errorKey += "run was killed. Check your input file and read the "
                     "manual for more information.\n";
      break;
    case OUTPUT_HELP:
      break;
    case OUTPUT_VERSION:
      break;
    case ERROR_STRUCTURE_EXPLOSION:
      bI.errorKey += "\nERROR:\tStructure optimization ran into troubles...\n";
      bI.errorKey += "\tTry running: TITANIA -i ";
      bI.errorKey += bI.inputFileName;
      bI.errorKey += " -r 0\n";
      break;
    case TITANIA_RDC_LINKING_ERROR:
      bI.errorKey += "\nERROR:\tTITANIA ran into problems while trying to link "
                     "RDCs of different media...\n";
      bI.errorKey += "\tThis behaviour might happen due to bad RDC inputs. Try "
                     "to remove single datasets\n";
      bI.errorKey += "\titerativly and try to find the corrupted data set! "
                     "Check this dataset for any\n";
      bI.errorKey += "\tpossible misstake (definition of nuclei, missing RDCs, "
                     "unneeded characters,...).\n";
      bI.errorKey += "\tIf console output is disabled: Run the programm with "
                     "active output! The corrupted\n";
      bI.errorKey += "\tset might be written to the terminal output...\n";
      break;
    case ERROR_OUTPUT_OPEN:
      bI.errorKey += ("\nERROR:\tTITANIA output file " + bI.outputFileName +
                      " could not be opened...\n");
      bI.errorKey += ("\tCheck the file and restart TITANIA...\n");
      break;
    case ERROR_TRJ_OPEN:
      bI.errorKey += ("\nERROR:\tTITANIA trajectory file " + bI.trjFileName +
                      " could not be opened...\n");
      bI.errorKey += ("\tCheck the file and restart TITANIA...\n");
      break;
    default:
      break;
  }
  if (f.debug && bI.debugFile.is_open())
  {
    if (bI.errorKey != "")
      bI.debugFile << bI.errorKey << std::endl;
    bI.debugFile.close();
    bI.debugFile.clear();
  }

  if (bI.output)
  {
    if (bI.errorKey != "")
      fprintf(bI.output, "%s", bI.errorKey.c_str());
    if (fclose(bI.output) == 0)
      bI.output = NULL;
  }
  if (bI.xyzFile.is_open())
  {
    bI.xyzFile.close();
    bI.xyzFile.clear();
  }
  if (bI.trjFile != NULL)
  {
    if (fclose(bI.trjFile) == 0)
      bI.trjFile = NULL;
  }
  if (bI.errorKey != "")
    std::cerr << bI.errorKey << std::endl;
  freeMemory(bI);
}

void
handle_stop(Structure *CurrStruc,
            BasicInformation &baseInformation,
            Flags &flags)
{
  unsigned int i;
  if (CurrStruc->check_Stop(baseInformation, flags))
  {
    if (!flags.silent)
      *std::cin.tie() << "              \r" << std::flush;
    *std::cin.tie() << "\t\tStructure optimization converged in "
                    << baseInformation.numOfOptSteps
                    << " steps due to:" << std::endl;
    for (i = 0; i < 6; ++i)
    {
      if ((unsigned int) baseInformation.MC_stop_reason & 1 << i)
      {
        MC_stop out = (MC_stop) (1 << i);
        *std::cin.tie() << "\t\t* low " << out << std::endl;
      }
    }
    if (baseInformation.MC_stop_reason & MC_stop::MaxIter)
    {
      MC_stop out = MC_stop::MaxIter;
      *std::cin.tie() << "\t\t* " << out << std::endl;
    }
    *std::cin.tie() << "\t\tTITANIA runs " << baseInformation.overOptimization
                    << " extra steps...\n";
  }
}

/*******************************************************
 *                                                     *
 * void baseSettings                                   *
 *                                                     *
 *  This function will handle the standard values for  *
 *               BasicSettings and Flags.              *
 *                                                     *
 * IMPORTAND: Keep the wiki page uptodate!             *
 *                                                     *
 * Wiki page exists.                                   *
 *                                                     *
 *******************************************************/
// TODO remove this function and build a proper constructor for BasicInformation
// and flags
void
baseSettings(BasicInformation &baseInformation, Flags &flags)
{
  baseInformation.numberOfInputs                   = 0;
  baseInformation.numberOfStructures               = 0;
  baseInformation.numOfThreads                     = 0;
  baseInformation.maxThreads                       = omp_get_num_procs();
  baseInformation.predictRDCs                      = 0;
  baseInformation.numberOfAlignments               = 0;
  baseInformation.numberOfFills                    = 0;
  baseInformation.state                            = GOOD_STATE;
  baseInformation.nestingDepth                     = 0;
  baseInformation.phobos_worksize_Sphericals       = 0;
  baseInformation.phobos_worksize_Eckart           = 0;
  baseInformation.numOfOptSteps                    = 0;
  baseInformation.digits                           = STANDARD_DIGITS;
  baseInformation.NumberOfSets                     = 0;
  baseInformation.NumberOfRDCs                     = 0;
  baseInformation.NumberOfAtoms                    = 0;
  baseInformation.NumberOfChiralCenters            = 0;
  baseInformation.full_mc_steps                    = 0;
  baseInformation.overOptimization                 = 2;
  baseInformation.useRedundantsOnlyAfter           = 0;
  baseInformation.memoryPurge                      = 0;
  baseInformation.lowerInversionAfter              = 250;
  baseInformation.redundants_distance_optimization = SKIP_DISTANCES_;

  baseInformation.normFactRDC        = -23000.0;
  baseInformation.MCvariation        = .5;
  baseInformation.numericalDeltaX    = 1e-3;
  baseInformation.redundants_damping = 0.0;
  baseInformation.SCRMtime           = .0;
  baseInformation.structureTime      = .0;
  baseInformation.MMFF94time         = .0;
  baseInformation.MCtime             = .0;
  baseInformation.inputParseTime     = .0;
  baseInformation.Systemtime         = .0;

  for (int red_type = BOND_REDUNDANTS_; red_type < NUMBER_OF_REDUNDANT_TYPES_;
       ++red_type)
    baseInformation.static_redundants_weighting[red_type] = 1.0;
  baseInformation.static_redundants_weighting[PLANAR_REDUNDANTS_] = 2.0;

  baseInformation.alignments                  = NULL;
  baseInformation.phobos_workspace_Sphericals = NULL;
  baseInformation.phobos_workspace_Eckart     = NULL;
  baseInformation.q_factors                   = NULL;
  baseInformation.chiral_volumes              = NULL;

  baseInformation.workingDirectory        = "";
  baseInformation.inputFileName           = "";
  baseInformation.outputFileName          = "";
  baseInformation.debugFileName           = "";
  baseInformation.comFileName             = "";
  baseInformation.errorKey                = "";
  baseInformation.gnuOutputFormat         = "pdf";
  baseInformation.structureLabelExtension = "_opt_";
  baseInformation.runtime_changes         = "";
  baseInformation.gpu_device_name         = "";
  baseInformation.structureInput          = StructureInputType::undefined;
  baseInformation.secondaryInput          = StructureInputType::undefined;
  baseInformation.MC_stop_reason          = MC_stop::undefined;

  baseInformation.output  = NULL;
  baseInformation.trjFile = NULL;

  baseInformation.limits.max_lm_iterations      = 50000;
  baseInformation.limits.max_titania_iterations = 10;
  baseInformation.limits.max_redundant_cycles   = STAN_MAX_REDUNDANT_CYCLES;
  baseInformation.limits.Q_factor_convergence   = 1e-3; // Absolute change!
  baseInformation.limits.redundants_convergence = 1e-3; // Absolute change!
  baseInformation.limits.redundants_validity = 5.0; // Relative (ratio) change!
  baseInformation.limits.alignment_mean_convergence = 1.0; // Percentage change!
  baseInformation.limits.alignment_sigm_convergence = 0.1; // Percentage change!
  baseInformation.limits.sphericals_mean_convergence =
      1.0; // Percentage change!
  baseInformation.limits.sphericals_sigm_convergence =
      0.1; // Percentage change!
  baseInformation.limits.sphericals_spread_convergence =
      0.1; // Percentage change!
  baseInformation.limits.zero_cutoff = ZERO_CUTOFF_;
  baseInformation.limits.phobos_opts_sphericals =
      (double *) malloc(4 * sizeof(double));
  baseInformation.limits.phobos_opts_eckart =
      (double *) malloc(4 * sizeof(double));
  baseInformation.limits.phobos_opts_sphericals[0] =
      STAN_TAU; // #define STAN_TAU  1E-03
  baseInformation.limits.phobos_opts_sphericals[1] =
      STAN_EPSILON; // #define STAN_EPSILON 1E-15
  baseInformation.limits.phobos_opts_sphericals[2] =
      STAN_EPSILON; // #define STAN_EPSILON 1E-15
  baseInformation.limits.phobos_opts_sphericals[3] =
      STAN_EPSILON; // #define STAN_EPSILON 1E-15
  baseInformation.limits.phobos_opts_eckart[0] =
      STAN_TAU; // #define STAN_TAU  1E-03
  baseInformation.limits.phobos_opts_eckart[1] =
      1E-10; // #define STAN_EPSILON 1E-15
  baseInformation.limits.phobos_opts_eckart[2] =
      1E-10; // #define STAN_EPSILON 1E-15
  baseInformation.limits.phobos_opts_eckart[3] =
      1E-10; // #define STAN_EPSILON 1E-15

  baseInformation.stop_crit.Q_factor_convergence          = .0;
  baseInformation.stop_crit.redundants_convergence        = .0;
  baseInformation.stop_crit.alignment_mean_convergence    = .0;
  baseInformation.stop_crit.alignment_sigm_convergence    = .0;
  baseInformation.stop_crit.sphericals_mean_convergence   = .0;
  baseInformation.stop_crit.sphericals_sigm_convergence   = .0;
  baseInformation.stop_crit.sphericals_spread_convergence = .0;


  flags.echo                         = false;
  flags.skipDIDC                     = true;
  flags.scaleDmatrix                 = true;
  flags.normDmatrix                  = true;
  flags.skipSCRM                     = false;
  flags.noInput                      = true;
  flags.noOutput                     = true;
  flags.appendOutput                 = false;
  flags.debug                        = false;
  flags.silent                       = false;
  flags.useRedundants                = true;
  flags.redundants_damping           = false;
  flags.titania2hotfcht              = false;
  flags.outputAli                    = true;
  flags.outputLM                     = false;
  flags.minimalOutput                = false;
  flags.monteCarloBootstrapping      = true;
  flags.plotKappaQ                   = true;
  flags.plotTrajectory               = true;
  flags.numericalGradients           = false;
  flags.weightMonteCarlo             = false;
  flags.print_redundants             = false;
  flags.planar_input                 = false;
  flags.skipEckart                   = false;
  flags.recalculateRDCs              = false;
  flags.printWarnings                = true;
  flags.titania2titania              = false;
  flags.normChiralVolume             = true;
  flags.SECONDAreducedCovariance     = false;
  flags.converged                    = false;
  flags.bigData                      = false;
  flags.plotRDCdynamics              = true;
  flags.scaleWithSoverall            = true;
  flags.monteCarloOutput             = false;
  flags.plotMonteCarlo               = false;
  flags.calculateFullMatrix          = true;
  flags.errorWeightInSVD             = false;
  flags.torsions_4_redundants        = false;
  flags.long_range_only_4_redundants = false;
  flags.floating_rdc_angles          = true;
  flags.use_initial_holonomics       = true;
  flags.plotRDCrmsd                  = false;
#ifdef USE_CUDA
  flags.use_gpu = true;
#else
  flags.use_gpu = false;
#endif
  flags.has_gpu = false;
}

int
main(int argc, char *argv[])
{
  /****************************************************
   *  Start TITANIA with defining standard variables  *
   *                 and setting up srand.            *
   ****************************************************/
  srand(time(NULL));
  struct BasicInformation
      baseInformation; // Saves the non boolean values for the program execution
  struct Flags flags;  // Saves the boolean values for the program execution
  time(&baseInformation.starting_time);
  gethostname(baseInformation.hostname, STANDARD_BUFFER_SIZE);
  getlogin_r(baseInformation.username, STANDARD_BUFFER_SIZE);
  if (getCurrentDir(baseInformation) == "")
  {
    *std::cin.tie() << "Warning:\tUnable to determine current directory...\n";
  }
  baseInformation.timeStamp =
      omp_get_wtime(); // Monitor the time for single blocks
  baseInformation.fulltime =
      baseInformation.timeStamp; // Monitor the full runtime
  std::fstream input;            // Obvious

  // Parse the arguments and load the standard settins of baseInformation and
  // flags
  switch (initializeProgram(baseInformation, flags, argc, argv))
  {
    case OUTPUT_HELP:
      killTITANIA(baseInformation, flags);
      exit(EXIT_SUCCESS);
    case OUTPUT_VERSION:
      *std::cin.tie() << "TITANIA running on version " << TITANIA_VERSION__
                      << std::endl;
      killTITANIA(baseInformation, flags);
      exit(EXIT_SUCCESS);
    case OUTPUT_BASE_SETTINGS:
      killTITANIA(baseInformation, flags);
      exit(EXIT_SUCCESS);
    case TITANIA_MULTI_JOB:
      killTITANIA(baseInformation, flags);
      exit(TITANIA_MULTI_JOB);
    case ERROR_INPUT_EQ_OUTPUT:
      killTITANIA(baseInformation, flags);
      exit(ERROR_INPUT_EQ_OUTPUT);
    case ERROR_NO_INPUT_FILE:
      killTITANIA(baseInformation, flags);
      exit(ERROR_NO_INPUT_FILE);
    default:
      break;
  }

  /****************************************************
   * By now all information from the arguments should *
   *     be loaded. Time to check the input file.     *
   ****************************************************/

  input.open(baseInformation.inputFileName, std::ios::binary | std::ios::in);
  if (!input.good())
  {
    std::cerr << "ERROR:\tThe input file " << baseInformation.inputFileName
              << " does not exist...\n"
                 "\tCheck again defined input file...\n";
    killTITANIA(baseInformation, flags);
    exit(EXIT_FAILURE);
  }


  Molecule CurrMol(""); // Main Molecule for the analysis
  MonteCarloStatistic MC_Structure_Run;


  /**************************************
   * ++++++++++++++++++++++++++++++++++ *
   *           Input parser             *
   * ++++++++++++++++++++++++++++++++++ *
   **************************************/
  addTimeSlot(baseInformation, baseInformation.Systemtime);

  InputFile reader(input, baseInformation.workingDirectory);
  if (ReadInput(CurrMol, reader, baseInformation, flags, input))
    exit(baseInformation.state);

  if (!flags.monteCarloBootstrapping)
    flags.plotTrajectory = false;

  if (CurrMol.getNORsets() < MIN_SECONDA_SETS_)
    flags.plotKappaQ = false;

  addTimeSlot(baseInformation, baseInformation.inputParseTime);
  *std::cin.tie() << "\nInformation:\tInput successfully parsed after after "
                  << ((baseInformation.inputParseTime)) << " s...\n";

  if (flags.titania2hotfcht)
    loadFCHTname(baseInformation);
#ifdef USE_CUDA
  if (flags.use_gpu)
  {
    flags.has_gpu = (get_number_of_gpu_devices() > 0 ? true : false);
    if (flags.has_gpu)
      baseInformation.gpu_device_name = get_gpu_device_name();
    else
      flags.use_gpu = false;
  }
#endif

/**************************************
 * ++++++++++++++++++++++++++++++++++ *
 *        Information collector       *
 * ++++++++++++++++++++++++++++++++++ *
 **************************************/
#pragma omp parallel shared(CurrMol, baseInformation, loading, flags, reader)
  {
    /************************************
     * ++++++++++++++++++++++++++++++++ *
     * Object for the whole program run *
     * ++++++++++++++++++++++++++++++++ *
     ************************************/

    Structure *CurrStruc; // Pointer for the current structure
    Atom *CurrAtom;       // Pointer for a specific Atom
    RDCset *CurrSet;      // Pointer for RDC sets

#pragma omp single
    {
      baseInformation.numOfThreadsUsed = omp_get_num_threads();
      handleMemory(baseInformation);

      if (baseInformation.predictRDCs == 0)
      {
        setupMolecule(CurrMol, flags);
        switch (handleOutput(baseInformation, flags))
        {
          case GOOD_STATE:
            break;
          case ERROR_OUTPUT_OPEN:
            killTITANIA(baseInformation, flags);
            exit(baseInformation.state);
          case ERROR_TRJ_OPEN:
            killTITANIA(baseInformation, flags);
            exit(baseInformation.state);
          default:
            break;
        }
#pragma omp task shared(baseInformation, flags, reader)
        {
          makeHeader(baseInformation.output, baseInformation, flags, reader);
          makeReferences(baseInformation.output, baseInformation);
        }
#pragma omp task shared(CurrMol, baseInformation, flags)
        {
          initializeRDCmatrix(CurrMol, baseInformation, flags);
        }
#pragma omp taskwait
        if (baseInformation.predictRDCs > 0)
          predictRDCs(CurrMol, baseInformation);
      }
      else /* predict rdcs defines only to predict new rdcs and skip rest */
      {
        flags.skipSCRM = flags.skipDIDC = true;
        predictRDCs(CurrMol, baseInformation);
      } /* End of predict rdcs only */

      /* Tolmans DIDC: */

      if (!flags.skipDIDC)
        *std::cin.tie()
            << "ERROR:\tThe DIDC of Tolman (currently) is not supported by "
               "TITANIA... and might be never :( ...\n";

      CurrStruc = CurrMol.getHeadStruc();

      if (CurrMol.getHeadSet())
      {
#pragma omp task shared(CurrMol, baseInformation) firstprivate(CurrStruc)
        baseInformation.SECONDA_output_initial =
            SECONDA_analysis(CurrMol, CurrStruc, baseInformation, flags);
      }
#pragma omp task shared(CurrMol, baseInformation) firstprivate(CurrStruc)
      SECONDA_sensitivity(CurrMol, CurrStruc, baseInformation, flags);

      if (baseInformation.predictRDCs == 0)
      {
        // Setup everything that is needed for a proper analysis.
        // This is used for the .trj output (which is generated even if skipSCRM
        // = true).
        CurrStruc->initializeRDCindex(CurrMol);
        initializeSphericalHarmonics(CurrMol, CurrStruc);
      }

      if (!flags.skipSCRM)
      {
        /*********************************************
         * +++++++++++++++++++++++++++++++++++++++++ *
         * Variables, matrizes and all for SCRM part *
         * +++++++++++++++++++++++++++++++++++++++++ *
         *********************************************/

        // Initialize optimization parameters
        struct Y_parameters Ydata;

#pragma omp critical
        *std::cin.tie()
            << "\nInformation:\tStarting the iterative optimization...\n";

#pragma omp task shared(CurrStruc)
        rdc_vector_analysis(CurrStruc, StructureOptions::Initial);
#pragma omp task firstprivate(CurrStruc) private(CurrAtom)
        {
          for (CurrAtom = CurrStruc->getHeadAtom(); CurrAtom;
               CurrAtom = CurrAtom->getNext())
            CurrAtom->set_chiral_volume_indizes();
          initializeChiralVolumes(CurrMol, baseInformation);
        }
        StructureSimulator MainSimulation(CurrStruc, flags);
        CurrMol.setNORA(
            CurrMol.getNOA() -
            MainSimulation.getNumberToOpt()); // Number of RDC defined atoms
        CurrStruc->retainCoordinates();
        RedundantInternals RedundantsOptimizer(&MainSimulation, baseInformation,
                                               flags);

        // Check if there was a 2 D structure in the input

        if (baseInformation.structureInput ==
            StructureInputType::xyzCoordinates)
          CurrMol.checkPlanarity(baseInformation);

        addTimeSlot(baseInformation, baseInformation.Systemtime);
        removePlanarity(CurrStruc, MainSimulation,
                        baseInformation); // writes coordinates to initial
        addTimeSlot(baseInformation, baseInformation.MMFF94time);

        for (CurrAtom = CurrStruc->getHeadAtom(); CurrAtom;
             CurrAtom = CurrAtom->getNext())
        {
#pragma omp task firstprivate(CurrAtom) shared(MainSimulation)
          {
            MainSimulation.AtomEnergy(CurrAtom, StructureOptions::Initial);
          }
        }

/* Handle the B/C matrix */
#pragma omp task firstprivate(CurrStruc)
        CurrMol.setInputBMatrix(CurrStruc);

#pragma omp task firstprivate(CurrStruc)
        outputXYZ(CurrStruc->getHeadAtom(), baseInformation,
                  StructureOptions::Initial, "Start Structure");

#pragma omp task firstprivate(CurrStruc, CurrAtom)
        {
          int center = 0;
          for (CurrAtom = CurrStruc->getHeadAtom(); CurrAtom;
               CurrAtom = CurrAtom->getNext())
          {
            if (!CurrAtom->hasChiralVolume())
              continue;
            baseInformation.chiral_volumes[center++] =
                CurrAtom->getChiralVolume(StructureOptions::Initial, true);
          }
        }
/*
 * Wait for
 * - Atom energies
 * - Calculation of cosine matrix
 * - Output of initial xyz coordinates.
 */
#pragma omp taskwait

#pragma omp task shared(baseInformation, flags, CurrMol)
        makeInputInformation(baseInformation.output, baseInformation, flags,
                             CurrMol);

        addTimeSlot(baseInformation, baseInformation.Systemtime);
        CurrAtom = CurrStruc->getHeadAtom();
#ifdef USE_FRAGMENTS_
        if (false) // Fragment tests
        {
          Fragment *CurrFragment = new Fragment(NULL, CurrStruc->getHeadAtom());
          CurrFragment           = CurrFragment->check_for_fusing();
          Fragment *HeadFragment = CurrFragment;
          Fragment *MainFragment = HeadFragment;
          unsigned int maxAt     = 0;
          while (CurrFragment)
          {
            CurrFragment->load_Potentials(
                MainSimulation.getPotential(PotentialType::Stretch),
                CurrStruc->getHeadYmatrix()->getHeadHarmonic(), flags,
                StructureOptions::Initial);
            //               CurrFragment->printCoordinates(baseInformation);
            //               CurrFragment->printJmol(); std::cout << std::endl;
            for (int i = 0; i < 5; ++i)
            {
              CurrFragment->perform_optimization_step();
              //                  CurrFragment->printCoordinates(baseInformation);
            }
            CurrFragment = CurrFragment->getNext();
          }
          for (CurrFragment = HeadFragment; CurrFragment;
               CurrFragment = CurrFragment->getNext())
          {
            if (CurrFragment->getNumberOfAtoms() > maxAt)
            {
              maxAt        = CurrFragment->getNumberOfAtoms();
              MainFragment = CurrFragment;
            }
          }
          for (CurrFragment = HeadFragment; CurrFragment;
               CurrFragment = CurrFragment->getNext())
          {
            if (CurrFragment == MainFragment)
              continue;
            CurrFragment->printCoordinates(baseInformation);
            //               MainFragment->printCoordinates(baseInformation);
            MainFragment->combine_Fragments(CurrFragment);
            std::cout << std::endl << std::endl;
            CurrFragment->printCoordinates(baseInformation);
          }
          outputXYZ(CurrAtom, baseInformation, StructureOptions::Initial,
                    "after_Frag");
        } // Fragment tests
#endif    // USE_FRAGMENTS_
#pragma omp task shared(CurrMol, CurrStruc, RedundantsOptimizer, \
                        baseInformation, flags)
        {
          prepare_trj_file(CurrMol, RedundantsOptimizer, baseInformation,
                           flags);
          outputStructure2trj(baseInformation, flags, CurrMol, CurrStruc,
                              RedundantsOptimizer, 0);
        }
        baseInformation.over_titania_iterations =
            baseInformation.limits.max_titania_iterations;
        // Start main optimization part
        for (baseInformation.numOfOptSteps = 1;
             baseInformation.numOfOptSteps <=
             baseInformation.over_titania_iterations;
             ++(baseInformation.numOfOptSteps))
        {
          if (!flags.silent)
            *std::cin.tie() << loading[baseInformation.numOfOptSteps % 4]
                            << "\r" << std::flush;
          /*
           * Compute all refined spherical harmonics, extract polar angles
           * and determine dynamical information.
           */
          computeAngles(CurrMol, CurrStruc, baseInformation, flags, Ydata);
          if (!flags.skipDIDC)
          {
            CurrStruc->toxyz(MainSimulation, RedundantsOptimizer,
                             baseInformation, flags,
                             StructureOptions::Optimized);
            compute_B_Angles(CurrStruc);
          }

          // Transform spherical coordinates to xyz
          // addTimeSlot is handled correctly in toxyz
          CurrStruc->toxyz(MainSimulation, RedundantsOptimizer, baseInformation,
                           flags, StructureOptions::Optimized);
#pragma omp task shared(CurrStruc)
          rdc_vector_analysis(CurrStruc, StructureOptions::Optimized);

          // Check if TITANIA has to stop due to max interations
          if (baseInformation.numOfOptSteps ==
              baseInformation.limits.max_titania_iterations)
            baseInformation.MC_stop_reason =
                baseInformation.MC_stop_reason | MC_stop::MaxIter;
#pragma omp taskwait
          addTimeSlot(baseInformation, baseInformation.Systemtime);

          // Run Monte-Carlo bootstrap on optimized structure
          if (flags.monteCarloBootstrapping)
          {
            baseInformation.full_mc_steps +=
                CurrStruc->run_structure_MC(baseInformation, flags);
            addTimeSlot(baseInformation, baseInformation.MCtime);
            if (baseInformation.state == ERROR_STRUCTURE_EXPLOSION)
            {
              killTITANIA(baseInformation, flags);
              exit(baseInformation.state);
            }
          }

          handle_stop(CurrStruc, baseInformation, flags);

          // If this is the last structure wait for the SECONDA analysis to
          // finish on it...
          if (baseInformation.numOfOptSteps <
              baseInformation.over_titania_iterations)
            outputStructure2trj(baseInformation, flags, CurrMol, CurrStruc,
                                RedundantsOptimizer,
                                baseInformation.numOfOptSteps);

          if (baseInformation.numOfOptSteps !=
              baseInformation.over_titania_iterations)
          {
            CurrStruc =
                CurrStruc->addOptimizedStructure(CurrMol, baseInformation);
            if (baseInformation.memoryPurge &&
                baseInformation.numOfOptSteps > 0 &&
                baseInformation.numOfOptSteps % 10 == 0)
              CurrMol.removeStructures(baseInformation.memoryPurge);
            addTimeSlot(baseInformation, baseInformation.Systemtime);
          }
        } // for ( baseInformation.numOfOptSteps )
        // Reduce numOfOptSteps by one. This is due to the implementation
        // of the for loop and the naming of the optimizes structures.
        // Otherwise One will allways have one step to much in the output
        // file.
        --(baseInformation.numOfOptSteps);
#pragma omp taskwait

        if (!flags.silent)
          *std::cin.tie() << std::string(CLEAN_LINE, ' ') << "\r" << std::flush;

        if (flags.monteCarloBootstrapping)
        {
          MC_Run_Information final_run;
          final_run.convergence    = MC_FINAL_RUN_CONVERGENCE;
          final_run.console_output = true;
          final_run.max_steps      = MC_FINAL_RUN_MAX_STEPS;

          baseInformation.full_mc_steps += MC_Structure_Run.startMonteCarlo(
              CurrMol, CurrStruc, baseInformation, flags, final_run);
          baseInformation.MonteCarloOutput =
              MC_Structure_Run.getStatistics(); // CurrStruc,baseInformation);
          addTimeSlot(baseInformation, baseInformation.MCtime);
        }
        if (CurrMol.getHeadSet())
        {
#pragma omp task shared(CurrMol, CurrStruc, RedundantsOptimizer)
          {
            baseInformation.SECONDA_output_final =
                SECONDA_analysis(CurrMol, CurrStruc, baseInformation, flags);
            outputStructure2trj(baseInformation, flags, CurrMol, CurrStruc,
                                RedundantsOptimizer,
                                baseInformation.numOfOptSteps);
          }
        }
      } /* End of skipSCRM */
      else
      {
        CurrStruc->retainCoordinates();
        if (!baseInformation.predictRDCs)
        {
          rdc_vector_analysis(CurrStruc, StructureOptions::Initial);
#pragma omp taskwait
          makeInputInformation(baseInformation.output, baseInformation, flags,
                               CurrMol);
        }
      }

      if (!flags.skipSCRM)
      {
#pragma omp task shared(CurrStruc)
        rdc_vector_analysis(CurrStruc, StructureOptions::Optimized);
      }

#pragma omp taskwait
      baseInformation.xyzFile.close();
      baseInformation.xyzFile.clear();


      /* Print results to output or shell */
      if (!baseInformation.predictRDCs)
      {
        makeTITANIAresults(baseInformation.output, baseInformation, flags,
                           CurrMol, CurrStruc);
// Most likely this block can be deleted!
// I kept it if any compatiblity problems would arise.
#ifdef ENABLE_OLD_OUTPUT_
        baseInformation.output.open(baseInformation.outputFileName,
                                    std::ios::binary | std::ios::out |
                                        std::ios::app);
        generateOutput(CurrMol, baseInformation, flags);
#endif
      }
      if (flags.monteCarloOutput && !flags.skipSCRM)
      {
        MC_Structure_Run.plot_mc_angles(CurrMol, CurrStruc, baseInformation,
                                        flags);
      }

      for (CurrSet = CurrMol.getHeadSet(); (flags.titania2hotfcht && CurrSet);
           CurrSet = CurrSet->getNext())
      {
#pragma omp task firstprivate(CurrSet) shared(CurrMol)
        titania2hotFCHT(CurrMol, CurrStruc, CurrSet, baseInformation, flags);
      }

      if (flags.titania2titania && !flags.skipSCRM)
      {
        std::cout << "Warning:\tOutput of titania2titania file is currently "
                     "disables...\n";
        if (false) // TODO
        {
#pragma omp task firstprivate(CurrStruc)
          titania2titania(CurrStruc, baseInformation);
        }
      }

// See comment above: This most likely can be deleted
#ifdef ENABLE_OLD_OUTPUT_
      baseInformation.output.close();
      baseInformation.output.clear();
#endif

      *std::cin.tie() << "\nInformation:\tOutput file was closed again...\n";
#ifdef DEACTIVATE
      // Currently there are no external scripts implemented one could wait for.
      if (!flags.silent)
        *std::cin.tie() << "\t\tWaiting for external scripts...\r"
                        << std::flush;
#endif
      if (flags.silent)
      {
        baseInformation.comFile.close();
        if (baseInformation.comFile.good() &&
            baseInformation.comFileName == ".null")
        {
          if (std::remove(baseInformation.comFileName.c_str()))
            std::cerr << "ERROR:\tSomething went wrong with the silent mode. "
                         "Contact the developer or set manually a dump file\n";
        }
        baseInformation.comFile.clear();
      }
#pragma omp taskwait
#ifdef DEACTIVATE
      if (!flags.silent)
        *std::cin.tie() << std::string(CLEAN_LINE, ' ') << "\r" << std::flush;
      *std::cin.tie() << "\t\tAll external scripts finished...\n";
#endif
      addTimeSlot(baseInformation, baseInformation.Systemtime);

      baseInformation.fulltime = omp_get_wtime() - baseInformation.fulltime;
      *std::cin.tie() << "\nInformation:\tComputation time: "
                      << baseInformation.fulltime << " s" << std::endl;

      if (baseInformation.numOfThreads)
        *std::cin.tie() << "\t\tUsed " << baseInformation.numOfThreadsUsed
                        << " threads, requested: "
                        << baseInformation.numOfThreads
                        << ", available: " << baseInformation.maxThreads
                        << std::endl;
      else
        *std::cin.tie() << "\t\tUsed: " << baseInformation.numOfThreadsUsed
                        << " threads, available: " << baseInformation.maxThreads
                        << std::endl;

      if (!baseInformation.predictRDCs)
      {
        makeRuntimeInformation(baseInformation.output, baseInformation, flags,
                               CurrMol);
      }

      // This has to be called in the very end!
      // If this functions are called earlier the MC-results will not
      // be appended to the output files.
      if (flags.monteCarloOutput && !flags.skipSCRM)
        MC_Structure_Run.print_mc_polar_angles(CurrStruc, baseInformation);

      freeMemory(baseInformation);

      CurrSet   = NULL;
      CurrStruc = NULL;
      CurrAtom  = NULL;
    } // #pragma omp single
  }   // #pragma omp parallel shared( CurrMol, baseInformation, loading, flags,
      // reader )
  return EXIT_SUCCESS;
}
