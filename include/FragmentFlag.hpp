#ifndef FRAGMENT_FLAG_TITANIA_
#define FRAGMENT_FLAG_TITANIA_

#include <ostream>

enum struct FragmentFlag : unsigned int
{
  Undefined = 0,
  RDC       = 1,
  NonRDC    = 1 << 1,
};

inline FragmentFlag
operator|(FragmentFlag l, FragmentFlag r)
{
  return static_cast<FragmentFlag>(static_cast<unsigned int>(l) |
                                   static_cast<unsigned int>(r));
}

inline bool
operator&(FragmentFlag l, FragmentFlag r)
{
  return static_cast<bool>(static_cast<unsigned int>(l) &
                           static_cast<unsigned int>(r));
}

inline std::ostream &
operator<<(std::ostream &os, FragmentFlag inp)
{
  switch (((unsigned int) inp))
  {
    case (0):
      os << "undefined fragment";
      break;
    case (1):
      os << "RDC defined fragment";
      break;
    case (1 << 1):
      os << "non-RDC defined fragment";
      break;
    default:
      os << "undefined fragment";
  }
  return os;
}

#endif // FRAGMENT_FLAG_TITANIA_
