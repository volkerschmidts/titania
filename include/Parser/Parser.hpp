#ifndef INPUTPARSER_HPP_
#define INPUTPARSER_HPP_

#include <fstream>
#include <string>
#include <vector>

#define FILESTART_HEX "0x460x690x6C0x650x530x740x610x720x74"
#define FILEEND_HEX   "0x460x690x6C0x650x450x6E0x64"

#define MOLFILE_VERSION "V2"

#ifndef OWNHASH_HPP_
enum class ValueType;
class OwnHash;
#endif

struct Keyword {
  std::string keyword;
  unsigned int Arguments;
  ValueType valueType;
};

#ifndef _LIST_OF_KEYWORDS_
#define _LIST_OF_KEYWORDS_
extern struct Keyword ListOfKeywords[]; //=
extern unsigned int
    NumberOfKeywords; // = sizeof(ListOfKeywords)/sizeof(Keyword);
#endif

#ifndef STRUCTURE_HPP_
class Structure;
#endif

class InputFile;

class StructureInput {
 private:
  std::vector<std::string> coordinates;
  std::vector<std::string> connectivity;
  unsigned int NOA;
  OwnHash *xyzTable;
  OwnHash *conTable;
  std::string label;
  std::string molecule;
  InputFile *parent;

 public:
  StructureInput();
  ~StructureInput();

  std::string getLabel() const
  {
    return label;
  }
  std::string getMolecule() const
  {
    return molecule;
  }
  unsigned int getNOA() const
  {
    return NOA;
  }
  std::vector<std::string> getKeys();
  int setAtomInformation(std::string &, double **, unsigned int &);
  std::string getBondPartner(std::string, unsigned int);
  int checkDiff(std::string &, std::string &);

  void parseXYZ(std::vector<std::string>, unsigned int &);
  void parseXYZ(std::fstream &);
  int parseMOLfile(std::string &, std::string &);
  void parseCon(std::vector<std::string>, unsigned int &);
  void parseCon(std::fstream &);
  void setParent(InputFile *f)
  {
    parent = f;
  }
  void printCoordinates();
  void printConnectivity();
};

class RDCinput {
 private:
  RDCinput *next;
  RDCinput *previous;
  std::vector<std::string> rdcs;
  OwnHash *rdcTable;
  std::string label;

 public:
  RDCinput();
  RDCinput(std::string);
  RDCinput(std::string, std::fstream &);
  RDCinput(std::string, std::vector<std::string>);
  ~RDCinput();

  OwnHash *initializeRDCtable(); /* Check rdcTable == NULL ) */
  OwnHash *getRDCtable()
  {
    return rdcTable;
  }
  std::vector<std::string> getKeys(); // { return keys; }
  std::vector<std::string> getRDCs()
  {
    return rdcs;
  }
  void reorder();
  double *getValuesD(std::string);

  RDCinput *setNext(std::string, std::fstream &); /* Call initializeRDCtable */
  RDCinput *append(std::string, std::fstream &);
  RDCinput *getNext()
  {
    return next;
  }
  RDCinput *getPrev()
  {
    return previous;
  }

  void setLabel(std::string s)
  {
    label = s;
  }
  std::string getLabel()
  {
    return label;
  }
};

class KeywordInput {
 private:
  std::vector<std::string> keys;
  OwnHash *keyTable;
  InputFile *parent;

 public:
  KeywordInput();
  KeywordInput(std::string);
  KeywordInput(std::vector<std::string>);
  KeywordInput(std::fstream &);
  ~KeywordInput();

  std::vector<std::string> getKeys(); // { return keys; }
  int getValueI(std::string, int);
  double *getValuesD(std::string);
  double getValueD(std::string, int);
  std::string getValueS(std::string, int);
  std::vector<std::string> getValuesS(std::string);
  bool contains(std::string);

  void addKeyword(std::string, Keyword);
  void removeEqSgn();
  void setParent(InputFile *f)
  {
    parent = f;
  }
};

class InputFile {
 private:
  StructureInput *structure;
  RDCinput *rdcs;
  KeywordInput *keywords;
  std::string workingDir;
  std::vector<std::string> content;
  std::vector<std::string> commentary;
  std::vector<std::string> unformated;
  std::vector<std::string> mainFile;

 public:
  InputFile(std::fstream &, std::string WD = "");

  int checkFile();

  void removeSpecialChars();

  StructureInput *getSI() const
  {
    return structure;
  }
  RDCinput *getRI() const
  {
    return rdcs;
  }
  KeywordInput *getKI() const
  {
    return keywords;
  }

  void getKeyLine(std::string &, std::string &, unsigned int &);
  void echo(std::fstream &);
  int comms(std::fstream &);

  void printCoordinates()
  {
    structure->printCoordinates();
  }
  void printConnectivity()
  {
    structure->printConnectivity();
  }
  void printKeywords()
  {
    ;
  } // keywords->printKeywords();}
  std::string findLabel(std::string, std::string);
  int includeFile(unsigned int &);

  int get_number_of_lines() const
  {
    return mainFile.size();
  }
  std::string get_line(const int l) const
  {
    return mainFile.at(l);
  }
};

void ta2csvString(std::string &);
int checkLine(std::string);

#endif
