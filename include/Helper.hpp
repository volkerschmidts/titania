
#ifndef HELPER_HPP_
#define HELPER_HPP_

#include <string>
#include <vector>

#ifndef BASIC_INFORMATION_TITANIA_
struct BasicInformation;
#endif

#ifndef FLAGS_TITANIA_
struct Flags;
#endif


unsigned int identifyKey(std::string);

template<typename T>
double zetaFunction(T, T, T);

template<typename T>
double deltaFunction(T, T);

template<typename T>
T max(T, T);
template<typename T>
T min(T, T);

double getPolar2Beta(double, // theta_1
                     double, //   phi_1
                     double, // theta_2
                     double  //   phi_2
);


double getDpolar2Dbeta(double, //        theta_1
                       double, // Delta(theta)_1
                       double, //          phi_1
                       double, // Delta( phi )_1
                       double, //        theta_2
                       double, // Delta(theta)_2
                       double, //          phi_2
                       double  // Delta( phi )_2
);


void addTimeSlot(BasicInformation &, double &);


void checkWorkingDir(BasicInformation &, Flags &);


std::string getCurrentDir(BasicInformation &);


void printKeywords(BasicInformation &, Flags &);


int loadFCHTname(BasicInformation &);

template<typename T>
bool vector_contains(std::vector<T>, T);
template<typename T>
int get_vector_index(std::vector<T>, T);

unsigned int getSizeOfUpperTriangle(const unsigned int);

unsigned long long getTotalSystemMemory();

#endif
