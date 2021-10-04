
#ifndef SCRM_H_
#define SCRM_H_

#include <eigen3/Eigen/Core>

#ifndef SPHERICALHARMONICS_HPP_
class SphericalHarmonics;
#endif

#ifndef SPHERICALHARMONICSMATRIX_HPP_
class SphericalHarmonicsMatrix;
#endif

#ifndef MOLECULE_HPP_
class Molecule;
#endif

#ifndef STRUCTURE_HPP_
class Structure;
#endif

#ifndef BASIC_INFORMATION_TITANIA_
struct BasicInformation;
#endif

#ifndef FLAGS_TITANIA_
struct Flags;
#endif

#ifndef Y_PARAMETERS_TITANIA_
struct Y_parameters;
#endif

#ifndef STRUCTUREOPTIONS_TITANIA_
enum struct StructureOptions : unsigned int;
#endif

int computeAngles(Molecule &,
                  Structure *,
                  BasicInformation &,
                  Flags &,
                  Y_parameters &);

int SCRM(Structure *, BasicInformation &, Flags &, Y_parameters &, bool);

int computeYref(Molecule &,
                Structure *,
                BasicInformation &,
                Flags &,
                Y_parameters &,
                StructureOptions,
                bool);

void analyzeYtensor(SphericalHarmonics *, Eigen::MatrixXcd);

void SrdcBySoverall(SphericalHarmonicsMatrix *);

#endif
