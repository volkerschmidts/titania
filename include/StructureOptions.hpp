#ifndef STRUCTUREOPTIONS_TITANIA_
#define STRUCTUREOPTIONS_TITANIA_

enum struct StructureOptions : unsigned int
{
  Undefined = 0,
  Initial   = 1,
  Optimized = 1 << 1,
};

inline StructureOptions
operator|(StructureOptions l, StructureOptions r)
{
  return static_cast<StructureOptions>(static_cast<unsigned int>(l) |
                                       static_cast<unsigned int>(r));
}

inline bool
operator&(StructureOptions l, StructureOptions r)
{
  return static_cast<bool>(static_cast<unsigned int>(l) &
                           static_cast<unsigned int>(r));
}
#endif /* STRUCTUREOPTIONS_TITANIA_ */
