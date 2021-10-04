# Building TITANIA from source

TITANIA uses the `cmake` build system to generate Makefiles. Currently, the
process is only tested on Linux (Debian and Ubuntu) with GNU `make` and the
GNU Compiler Collection toolchain.

## Dependencies

TITANIA's build process is based on CMake and generates standard Makefiles.
TITANIA relies on some external libraries. For runtime, these are [CBLAS/LAPACKE](https://netlib.org/lapack/)
(note that not all BLAS/LAPACK packages may bring the necessary C/C++ interfaces)
and [NLopt](https://github.com/stevengj/nlopt). In cases where these libraries
cannot be installed system-wide, the libraries may be built locally and consumed
in TITANIA's build process.

For building TITANIA from source, an OpenMP-compatible C/C++ compiler is required.
Compiling the CBLAS/LAPACKE dependency additionally requires a Fortran compiler.
Development is extensively tested with the GNU Compiler Collection (GCC 10 and 11).

To compile the TITANIA sources also depends on the [Eigen3](https://eigen.tuxfamily.org/index.php?title=Main_Page)
headers being present. Compiling TITANIA's (unit) test suite requires the
[CPPunit](https://freedesktop.org/wiki/Software/cppunit/) library. Compiling
CPPunit from source requires the libtool and pkg-config packages.

On Debian / Ubuntu based systems, these dependencies can be installed via
```console
$ sudo apt install build-essential cmake libeigen3-dev liblapacke-dev libnlopt-cxx-dev libcppunit-dev
```

## Configure the Build

To start a build, create a target directory (for example `./build`) in which all
generated files are stored. Then configure the build with CMake, potentially
adding user definitions as required. Unless specified otherwise, the default build
type is `RelWithDebugInfo`.
```console
$ cd /path/to/titania/source/directory
$ mkdir build
$ cd build
$ cmake .. -DCMAKE_INSTALL_PREFIX=/my/local/TITANIA/install/path
$ make
$ make install
```

The `make install` command installs the TITANIA binary to `<PREFIX>/bin` with
`PREFIX=/usr/local/` by default. No other binaries, no libraries or headers are
installed otherwise. Once installed, there is no default way to uninstall the
files. CMake recommends to manually check the files listed in the file
`build/install_manifest.txt` and remove them manually.

## Setting Specific Configuration Settings

The default project settings should work out-of-the-box in most cases. To change
the build configuration, use the `cmake-gui` GUI or the `ccmake` curses GUI if
available. Most settings are also accessible as command-line arguments or
environment variables.

Examples:
- to generate a configuration for a 'release' build from the top-level source
  directory you could use:  
  `cmake -S . -B ./build -DCMAKE_BUILD_TYPE=release`
- to build an already configured target in the './build' directory with 6 parallel
  processes you could use:  
  `cmake --build ./build -DCMAKE_BUILD_PARALLEL_LEVEL=5`

On a system where all dependencies are installed, the default configuration should
work out of the box. Some users may wish to fully build the libraries themselves
from source. After building and installing the dependency, simply providing a
`<LIBRARY>_DIR=...` definition to TITANIA's cmake build step should be sufficient
for a successful build. The more general option of appending the target path to
`CMAKE_PREFIX_PATH` should also allow the build system to pick up the dependencies,
but be careful not to override the system's default paths. We also provide an
automated way to download and build the libraries. To enable this automatic process
use the `TITANIA_BUILD_3RDPARTY_<LIBRARY>` cmake flag.

For example, to rely on the system libraries but locally download and build NLopt
from GitHub, one would do:
```console
$ cd /path/to/titania/source/directory
$ mkdir build
$ cmake -S . -B ./build -DTITANIA_BUILD_3RDPARTY_NLOPT=ON
$ cmake --build ./build
```

or pointing cmake to pre-built libraries:
```console
$ cd /path/to/titania/source/directory
$ mkdir build
$ cmake -S . -B ./build -DNLopt_DIR=/path/to/my/workspace/nlopt/build
$ cmake --build ./build
```

## Testing

To build TITANIA's test suite, set the `-DENABLE_TESTING=ON` flag and run `make test`
after compilation. Note that this may take several minutes to finish as some of
the tests internally run expensive structure optimizations. To see individual
test progress run `ctest --verbose` instead. To run single tests run the unit
test binary with the test name as argument (for example: `unit/TITANIA_UNIT_TEST AtomTest`).
```console
$ cd /path/to/titania/source/directory
$ mkdir build
$ cd build
$ cmake .. -DENABLE_TESTING=ON
$ make
$ make test
```
