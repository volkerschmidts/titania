    .___________. __  .___________.    ___      .__   __.  __       ___  
    |           ||  | |           |   /   \     |  \ |  | |  |     /   \  
    `---|  |----`|  | `---|  |----`  /  ^  \    |   \|  | |  |    /  ^  \  
        |  |     |  |     |  |      /  /_\  \   |  . `  | |  |   /  /_\  \  
        |  |     |  |     |  |     /  _____  \  |  |\   | |  |  /  _____  \  
        |__|     |__|     |__|    /__/     \__\ |__| \__| |__| /__/     \__\

TITANIA is a program package to directly extract structure and dynamics parameters
from Residual Dipolar Couplings for weakly aligned organic compounds using an
iterative model-free approach. RDC data from multiple alignment sets are analyzed
to jointly determine the conformation and configuration of organic compounds.

Some applications of the TITANIA program package are (in preparation to be)
published in:
- F. A. Roth, V. Schmidts, C. M. Thiele, "TITANIA: Model Free Interpretation of
    Residual Dipolar Couplings in the context of Organic Compounds",
    ChemRxiv 2021, DOI [10.26434/chemrxiv.14636070.v1](https://doi.org/10.26434/chemrxiv.14636070.v1).
- F. A. Roth, V. Schmidts, J. Rettig, C. M. Thiele, "Incomplete Data Sets in the
    Model Free Analysis of Experimental Residual Dipolar Couplings in Small
    Organic Compounds", ChemRxiv 2021, DOI [10.26434/chemrxiv.14636316.v1](https://doi.org/10.26434/chemrxiv.14636316.v1).

Please cite our work if you use TITANIA in your RDC analysis projects.

## Compiling from Source

TITANIA can be compiled from sources using the CMake build system. The necessary
steps are described in the [CompileNotes](CompileNotes.md). Building has been
tested extensively on recent Debian and Ubuntu systems but should also work
on other modern Linux distributions, provided the dependencies are setup properly.

## Example Data Sets and Scripts

The [standard_inputs](standard_inputs) directory contains some example input
files for common organic compounds which can be run via `TITANIA <inputfile>`.

The [scripts](scripts) directory contains some helper bash, python and gnuplot
scripts that were used to prepare the figures in the publications mentioned above.
