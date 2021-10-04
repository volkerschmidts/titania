
#ifndef RDCSET_HPP_
#define RDCSET_HPP_

#include <eigen3/Eigen/Core>
#include <string>

#ifndef RDCDATA_HPP_
class RDCdata;
#endif

#ifndef MOLECULE_HPP_
class Molecule;
#endif

#ifndef FLAGS_TITANIA_
struct Flags;
#endif

#ifndef STRUCTUREOPTIONS_TITANIA_
enum struct StructureOptions : unsigned int;
#endif

class RDCset {
 private:
  RDCdata *HeadData; /* First rdc in this set */
  RDCdata *TailData; /* Last rdc in this set */
  RDCset *prevSet;   /* Previous rdc set defined by input file */
  RDCset *nextSet;   /* Next rdc set defined by input file */
  Molecule *parent;  /* Molecule which this set belongs to. Since sets do not
                        have to be changed after SCRM optimization this is saved
                        in Molecules! */
  std::string label; /* String representation of this set (e.g. alignment
                        medium, temperature, etc.) */
  unsigned int inputIndex;   /* Input index of this set based on input file */
  unsigned int NumberOfRDCs; /* Number of rdcs in this set (allready in
                                Molecule, needed?) */
  double SECONDA_gap_5_6_sensitivity;
  double SECONDA_gap_1_5_sensitivity;
  int SECONDA_sensitivity_rank;
  Eigen::MatrixXd weighting;

 public:
  /* Standard constructor */
  RDCset();
  /* Copy constructor for faster appending of nodes */
  RDCset(RDCset *s);
  RDCset(Molecule &, Molecule &, Flags &);

  /* Standard destructor */
  ~RDCset();

  /* Accessing label */
  void setLabel(std::string nL)
  {
    label = nL;
  }
  //      std::string getLabel () const { return label; }
  std::string getLabel(int sub = 0, bool extend = false);

  /* Accessing input index of this set */
  void setIndex(unsigned int ii)
  {
    inputIndex = ii;
  }
  unsigned int getIndex() const
  {
    return inputIndex;
  }

  /* Functions to handle nodes of rdc sets*/
  RDCset *appendSet();
  RDCset *getNext() const
  {
    return nextSet;
  };
  RDCset *getPrev() const
  {
    return prevSet;
  };

  /* Accessing rdcs in this set */
  RDCdata *addHeadData();
  void setTailData(RDCdata *t)
  {
    TailData = t;
  };
  RDCdata *getHeadData() const
  {
    return HeadData;
  };
  RDCdata *getTailData() const
  {
    return TailData;
  };
  void copyRDCs(RDCset *);

  /* Accessing Molecule, which defines this sets */
  void setParent(Molecule *p)
  {
    parent = p;
  };
  Molecule *getParent() const
  {
    return parent;
  };


  /* Internal function to obtain number of rdcs */
  unsigned int raiseNOR()
  {
    return ++NumberOfRDCs;
  };
  void setNOR(const unsigned int nor)
  {
    NumberOfRDCs = nor;
  }

  void set_SECONDA_gap_5_6_sensitivity(const double gap)
  {
    SECONDA_gap_5_6_sensitivity = gap;
  }
  void set_SECONDA_gap_1_5_sensitivity(const double gap)
  {
    SECONDA_gap_1_5_sensitivity = gap;
  }
  void set_SECONDA_sensitivity_rank(const int rank)
  {
    SECONDA_sensitivity_rank = rank;
  }
  double get_SECONDA_gap_5_6_sensitivity() const
  {
    return SECONDA_gap_5_6_sensitivity;
  }
  double get_SECONDA_gap_1_5_sensitivity() const
  {
    return SECONDA_gap_1_5_sensitivity;
  }
  int get_SECONDA_sensitivity_rank() const
  {
    return SECONDA_sensitivity_rank;
  }

  double get_sigma_square();

  Eigen::MatrixXd getWeighting();
};

#endif
