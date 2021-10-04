#ifndef REDUNDANT_CALCULATION_TITANIA_
#define REDUNDANT_CALCULATION_TITANIA_

enum struct RedundantCalculation : unsigned int
{
  Undefined = 0,
  QInitial  = 1,
  SqRDCs    = 1 << 1,
  Standard  = 1 << 2,
};

#endif
