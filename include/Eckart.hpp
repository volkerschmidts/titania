
#ifndef ECKART_HPP_
#define ECKART_HPP_

#include <eigen3/Eigen/Core>

#ifndef STRUCTURE_HPP_
class Structure;
#endif

#ifndef AXIS_TITANIA_
enum struct Axis : unsigned int;
#endif

#ifndef STRUCTUREOPTIONS_TITANIA_
enum struct StructureOptions : unsigned int;
#endif

#ifndef BASIC_INFORMATION_TITANIA_
struct BasicInformation;
#endif

void build_Eckart_AlphaTau(Structure *,
                           Structure *,
                           Eigen::Matrix3d &,
                           StructureOptions);


void build_Eckart_Rotation(double *, Eigen::Matrix3d &);


void build_Eckart_Rotation(double, double, double, Eigen::Matrix3d &);


void get_Eckart_Angles(Structure *,
                       Structure *,
                       BasicInformation &,
                       StructureOptions);


double Eckart_Lambda(Axis, Axis, double, double, double);


double Eckart_AlphaTau(Axis, Axis, Structure *, Structure *, StructureOptions);


void Eckart_Angles(double *, double *, int, int, void *);


double *checkEckart(Structure *, Structure *, StructureOptions);


void Rotate_2_Eckart_Frame(Structure *,
                           Structure *,
                           BasicInformation &,
                           StructureOptions);

#endif
