##################################################
#       Standard keywords / flags file for       #
#      Titania performs ITerativ ANalysis of     #
#        Independent Alignments (TITANIA).       #
#                                                #
#   All base settings used by TITANIA are given  #
#   is (first) value in the commentary section.  #
#     If only specific values are valid, the     #
#   respective settings can be find in the same  #
#                    section.                    #
# ################################################


#echo = 0 # [0|1]
# Reprints the whole input in the output.

#LMmaxIterations = 1000  # [50000]
# Sets the maximum iteration for the Levenberg-
# Marquardt minimization.

#LMoptions = 1e-2 1e-10 1e-10 1e-15 # [1e-3 1e-15 1e-15 1e-15]
# Sets the Levenberg-Marquardt minimization 
# parameters. Syntax: 
# LMopts = tau epsilon1 epsilon2 epsilon3

#numericalgradients = 1 # [0|1]
# Forces Levenberg-Marquardt minimizer to
# use numerical gradients. There should be
# no scenario where this keyword is needed.

MaxTITANIAiterations = 1000 #[10]
# Sets the maximum iteration cycles of TITANIA
# (including the determination of the saupe tensor,
# the spherical harmonic and the refinded structure).
# Overoptimization steps (see below) are not included.

overoptimizationsteps = 5 #[2]
# Sets the number of steps that are performed
# after a stop critereon was reached.

qfactorconvergence        = 0.0 # [1e-3]
# Stop criterion: Q-factor change.

MeanAlignmentConvergence  = 1e-1  # [1e0]
# Stop criterion: rmsd(A_i - A_i-1)

SigmaAlignmentConvergence = 1e-2 # [1e-1] 
# Stop criterion: rmsd(sigma(A_i) - sigma(A_i-1))

MeanAngleConvergence      = 1e-1 # [1e0]
# Stop criterion: rmsd(p_i - p_i-1)

SigmaAngleConvergence     = 1e-2 # [1e-1]
# Stop criterion: rmsd(sigma(p_i) - sigma(p_i-1))

SpreadAngleConvergence    = 1e-2 # [1e-1]
# Stop criterion: rmsd(R^2(p_i) - R^2(p_i-1))

MonteCarloBootstrapping = 0 # [1|0]
# En-/diables Monte-Carlo bootstrapping on every
# iteration step. If disabled all alignment
# and angle stop criterions are ignored! 
# This is due the fact, that all information
# are collected by MC-results.

useRedundantsOnlyAfter = 50 # [0]
# Enables rmsd and vector addition based structure 
# refinement previous to redundant internal 
# coordinates for n steps. This can speed up the 
# algorithm significantly, but also can mess up
# the refinement when using poor data or too complex
# molecules.

useDistances = 1 # [0] - 0: skip, 1: full, 2: red, 3: inv
# This enables the use of distance information leading
# to better results (especially for Atoms with few
# RDC information [e.g. methyl groups] or no RDC information
# [e.g. hydroxy groups]). 
# If used with argument 1 or 2 will massivly increase the
# calculation time for large molecules. In this case TITANIA
# might be used with CUDA.

#RedundantsConvergence = 1e-3  # [1e-3]
# Sets the stop criterion for redundant internal
# coordinates. The convergence is an rmsd for
# RDC information.

#redundantCycles = 3 # [5]
# Sets the maximum number of redundant internal
# coordinates optimization cycles. 
# CAUTION: To high values can lead to unpredictible
# behaviour in the structure generation. Low
# values will:
#  a) speed up optimizations for good structures.
#  b) slow down optimizations for poor structures.

redundantsDamping           = 50 # [0|{e|e>1.0}]
# Enables a smoother (and slower) optimization
# using redundant internal coordinates. 
# The damping uses a sigmoid function which 
# is normalized to 1.0 at MaxTITANIAiterations.
# This damping is disabled for rdc and 
# distance information!

StaticBondWeighting         = 1.0 # [1.0]
StaticAngleWeighting        = 1.0 # [1.0]
StaticTorsionWeighting      = 1.0 # [1.0]
StaticRDCWeighting          = 1.0 # [1.0]
StaticChiralVolumeWeighting = 1.5 # [1.0]
StaticDistanceWeighting     = 1.5 # [1.0]
# Static weightings counteract redundants damping
# (or might increase the effect if <1.0). 

#skipEckart = 1 # [0|1]
# This will prevent the Eckart transformation after the
# Structure refinemant via redundant internal coordinates.

InversionStrictness = 150 # [250]
# if useDistances = [1|3] (see below) an inversion is 
# possible after 5 violations. InversionStictness
# defines the iteration after which this violation
# count is increased to 10. This will effect in 
# an slower but more constistent structure refinement.

#printredundants = 1 # [0|1]
# Enables the output of the structures generated
# by every redundant internal coordinate step.
# This might lead to a massive output in the .xyz
# file!

#skipSCRM = 1 # [0|1]
# Prevents TITANIA from running its standard cycle. 
# Can be used to just perform SECONDA.

#TITANIA2hotFCHT = 1 # [0|1]
# Generates RDC@hotFCHT conform input files of the 
# final structure for each RDC set. The RDC@hotFCHT 
# keywords will be set to run the calculation 
# according to TITANIA standard settings.
 
#TITANIA2TITANIA = 1 #[0|1]
# Generates a TITANIA input file with the same 
# settings and datasets like in the current run. 
# The input structure will be replaced with the 
# current structure. 
# CAUTION: In the current version TITANIA can not 
# alter the structure if it is included via an 
# external file! 

outputali = 0 # [1|0]
# Generates an output file for every single alignment 
# medium containing all information on the alignment.

#bigdata = 1 # [0|1]
# This will generate an output file containing a LOT 
# of information. This is not set up properly! 
# Reworks will follow in the future.

PlotKappaQ = 0 # [1|0]
# Generates a SECONDA plot for the RDC sets using
# the initial and final structure (for normalization
# of long range couplings) via gnuplot.

PlotChiralVolume = 0 # [1|0]
# Generates a plot of chiral volumes with respect to 
# the iteration step via gnuplot.

plotrdcdynamics = 0 # [1|0]
# Generates a plot of all S^2 values and S(overall) 
# via gnuplot.

#plotMonteCarlo = 1 # [0|1]
# Generates a plot containing the Monte Carlo results 
# of the polar coordinates (for each RDC vector) in 
# final MC run via gnuplot. 
# Needs:  - MonteCarloBootstrapping = 1 
#         - MonteCarloOutput = 1 (see below)

#MonteCarloOutput = 1 # [0|1]
# Appends the Monte Carlo results of the final MC run
# to the standard output file. 
# Needs:  - MonteCarloBootstrapping = 1 

#gnuOutputFormat = pdf # [pdf|png|svg]
# Defines the standard output format of the automatically
# generated plot files.

#recalculateRDCs = 1 [0|1]
# This can only be used if datasets are incomplete. If a
# RDC set hat missing values:
# recalculateRDCs = 0 -> TITANIA will interpret them as 0
# recalculateRDCs = 1 -> TITANIA will calculate them on the end
#   of every step with the current alignment and optiized 
#   structure.
#
# For additional options see CalculateFullMatrix keyword (below)

#CalculateFullMatrix = 0 # [1|0]
# TITANIA solves on default two systems of equations:
#  w * D = w * B * A ( all in matrix representation )
#  w * D = w * F * Y ( all in matrix representation )
#
# The weighting w cant target a single RDC in a single
# set if this is unknown of poorly determined.
# To do so one can split the respective system of 
# equations and solve them all in vector representation.
# By this one can neglect a single RDC (e.g. a missing one)
# for the alignment calculation and the determination
# of the polar angles.

#errorWeightInSVD = 1 # [0|1]
# Weights every RDC with its respective user defined
# error. This keyword can only be used if CalculateFullMatrix
# is turned off!

#printWarnings = 0 #[1|0]
# This will prevent TITANIA from printing warnings. Not 
# included are errors (which are output via the error
# stream anyway) nor the standard output of TITANIA.

#normChiralVolumes = 0 #[1|0]
# Prevents the normalization of the bond vectors when
# calulating the chiral volumen. This might lead to 
# a higher dispersion of the results but remove the 
# meaning of the values on the first sight.

#ScaleWithSoverall = 1 # [0|1]
# This will lead to a second run on the systems of
# equations (explanation see CalculateFullMatrix) with
# scaling of the RDCs with Soverall of the first run.

useGPU = 1 # [USE_CUDA]
# Enables the use of the GPU. This keyword is only 
# available if TITANIA was compiled with the USE_CUDA flag.
# Otherwise TITANIA is not able to use GPU/CUDA code.

MemoryPurge = 0.8 # [0]
# On large systems (>100 atoms) with a large amount of 
# iterations (>500) TITANIA might use to much RAM. If you
# are using a machine with > 16GB RAM you might not run 
# into troubles. Otherwise you should activate MemoryPurge
# with an float value of 0.5-0.9. This will enable TITANIA
# to free unused memory (respective to the defined value 
# in percent).

#silent = 1 # [0|1|<outputfile>]
# This will stop the output of TITANIA to the console. If
# just set to 1 the output will be redirected to /dev/null.
# If an output is defined the respective file will be used.
# Since TITANIA uses string modifiers in the console output
# this option is prefered over using Linux ">" directive.

#threads = 1
# Defines a standard number of threads for the parallel run.
# This value is ignored if the user defined another value
# via the "-c" argument.

