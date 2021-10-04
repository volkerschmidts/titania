
#ifndef OWNHASH_HPP_
#define OWNHASH_HPP_

#include <string>
#include <vector>

class OwnHash;
struct Keyword;

enum struct HashType
{
  undefined    = 0,
  xyz          = 1,
  connectivity = 1 << 1,
  rdc          = 1 << 2,
  keyword      = 1 << 3,
};

enum struct ValueType
{
  undefined  = 0,
  intType    = 1,
  doubleType = 1 << 1,
  stringType = 1 << 2,
};

class String_ptr {
 private:
  std::string content;
  String_ptr *next;
  String_ptr *prev;
  int index;

 public:
  String_ptr()
  {
    content = "";
    prev = next = NULL;
    index       = 0;
  }
  ~String_ptr()
  {
    next = NULL;
    prev = NULL;
  }
  std::string getContent(int);
  String_ptr *setContent(int, std::string);
  String_ptr *getNext()
  {
    return next;
  }
  String_ptr *addNext();
  int getIndex()
  {
    return index;
  }
};

class HashElement {
 private:
  HashElement *next;
  HashElement *previous;
  HashElement *first;
  HashElement *last;
  OwnHash *parent;

  enum HashType type;
  std::string HashIdent;
  double *valueD;
  int *valueI;
  String_ptr *valueS;
  int hashIndex;

 public:
  HashElement();
  HashElement(HashType);
  HashElement(HashType, std::string, Keyword *k = NULL);
  ~HashElement();

  HashElement *getNext();
  HashElement *getPrev();
  HashElement *getFirst();
  HashElement *getLast();

  void setNext(HashElement *); // Just define this -> no reordering of the list!

  void setType(HashType);

  std::string setIdentifier(std::string);
  std::string getIdentifier()
  {
    return HashIdent;
  };

  double *setValueD(int, double);
  double *setValuesD(int, double *);

  int *setValueI(int, int);
  int *setValuesI(int, int *);

  String_ptr *setValueS(int, std::string);

  double *getValueD(int);
  double *getValuesD();

  int *getValueI(int);
  int *getValuesI();

  std::string getValueS(int);
};


class OwnHash {
 private:
  HashType type;
  HashElement *first;
  HashElement *last;
  HashElement *current;

 public:
  OwnHash();
  OwnHash(HashType);
  ~OwnHash();

  HashElement *addElement(std::string, Keyword *k = NULL);
  HashElement *getElement(std::string);
  HashElement *getElement(int);

  std::vector<std::string> getKeys();

  void setValueD(std::string, int, double);
  void setValuesD(std::string, int, double *);
  void setValueI(std::string, int, int);
  void setValuesI(std::string, int, int *);
  void setValueS(std::string, int, std::string);
  void setValueS(std::string, std::string);

  double getValueD(std::string, int);
  double *getValuesD(std::string);
  int getValueI(std::string, int);
  int *getValuesI(std::string);
  std::string getValueS(std::string, int);
  std::vector<std::string> getValuesS(std::string);

  void printrdcTable();
  void printCoordinates();
  void printConn();
  void printKeys();
};

#endif
