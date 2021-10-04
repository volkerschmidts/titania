#ifndef AXESDEFINITIONS_TITANIA_
#define AXESDEFINITIONS_TITANIA_

#define STANDARD_AXIS_X_   0
#define STANDARD_AXIS_Y_   1
#define STANDARD_AXIS_Z_   2
#define STANDARD_AXIS_R_   3
#define NUMBER_OF_AXIS_    3
#define NUMBER_OF_4D_AXIS_ 4
#define STANDARD_THETA_    0
#define STANDARD_PHI_      1
constexpr auto double_vector_size = (NUMBER_OF_AXIS_ * sizeof(double));

#endif