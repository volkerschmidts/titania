#ifndef OUTPUT_OPTIONS_TITANIA_
#define OUTPUT_OPTIONS_TITANIA_

enum struct OutputOptions : unsigned int
{
  undefined  = 0,
  left       = 1,
  center     = 1 << 1,
  right      = 1 << 2,
  scientific = 1 << 3,
};
inline OutputOptions
operator|(OutputOptions l, OutputOptions r)
{
  return static_cast<OutputOptions>(static_cast<unsigned int>(l) |
                                    static_cast<unsigned int>(r));
}

inline bool
operator&(OutputOptions l, OutputOptions r)
{
  return static_cast<bool>(static_cast<unsigned int>(l) &
                           static_cast<unsigned int>(r));
}

#endif