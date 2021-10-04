#include <Parser/OwnHash.hpp>
#include <Parser/Parser.hpp>
#include <fstream>
#include <iostream>
#include <vector>


#define MAX_NUM_KEYINTS  1
#define MAX_NUM_KEYDOUBS 5

String_ptr *
String_ptr::addNext()
{
  next = new String_ptr(); //(String_ptr*) malloc ( sizeof ( String_ptr) );
                           ////new String_ptr();
  next->index   = (index + 1);
  next->prev    = this;
  next->content = "";
  return next;
}

std::string
String_ptr::getContent(int i)
{
  if (i == index)
    return content;
  else
  {
    String_ptr *tmp = this;
    while (tmp->index >= i)
      tmp = tmp->prev;
    while (tmp->next)
    {
      tmp = tmp->next;
      if (tmp->index == i)
      {
        return tmp->content;
      }
    }
  }
  return "";
}

String_ptr *
String_ptr::setContent(int i, std::string s)
{
  if (i == index)
  {
    content = s;
    return this;
  }
  else
  {
    String_ptr *tmp = this;
    //      while ( tmp->index >= i ) tmp = tmp->prev;
    while (tmp->next)
    {
      tmp = tmp->next;
      if (tmp->index == i)
      {
        tmp->content = s;
        return tmp;
      }
    }
    tmp->next          = tmp->addNext();
    tmp->next->index   = i;
    tmp->next->content = s;
    return tmp->next;
  }
  return NULL;
}

HashElement::HashElement()
{
  next = previous = NULL;
  first = last = this;
  parent       = NULL;

  type      = HashType::undefined;
  HashIdent = "";
  valueD    = NULL;
  valueI    = NULL;
  valueS    = NULL;

  hashIndex = 0;
}

HashElement::HashElement(enum HashType ht)
{
  next = previous = NULL;
  first = last = this;
  parent       = NULL;
  ;

  type      = ht;
  HashIdent = "";
  valueD    = NULL;
  valueI    = NULL;
  valueS    = NULL;

  if (ht == HashType::xyz || ht == HashType::rdc)
  {
    valueD = (double *) malloc(3 * sizeof(double));
    for (unsigned int i = 0; i < 3; ++i)
      valueD[i] = .0;
  }

  hashIndex = 0;
}

HashElement::HashElement(enum HashType ht, std::string s, Keyword *k)
{
  next = previous = NULL;
  first = last = this;
  parent       = NULL;
  ;

  type      = ht;
  HashIdent = s;
  valueD    = NULL;
  valueI    = NULL;
  valueS    = NULL;

  if (ht == HashType::xyz || ht == HashType::rdc)
  {
    valueD = (double *) malloc(3 * sizeof(double));
    for (unsigned int i = 0; i < 3; ++i)
      valueD[i] = .0;
  }
  else if (ht == HashType::keyword)
  {
    int i, j, d;
    i = j = d = 0;
    if (k)
    {
      if (k->valueType == ValueType::undefined)
      {
        i = MAX_NUM_KEYINTS;
        d = MAX_NUM_KEYDOUBS;
      }
      else if (k->valueType == ValueType::intType)
      {
        i = k->Arguments;
        d = 0;
      }
      else if (k->valueType == ValueType::doubleType)
      {
        d = k->Arguments;
        i = 0;
      }
      else if (k->valueType == ValueType::stringType)
      {
        i = d = 0;
      }
    }
    if (i)
      valueI = (int *) malloc(i * sizeof(int));
    for (j = 0; j < i; ++j)
      valueI[j] = 0;

    if (d)
      valueD = (double *) malloc(d * sizeof(double));
    for (j = 0; j < d; ++j)
      valueD[j] = .0;
  }

  hashIndex = 0;
}

/*HashElement::HashElement ( enum HashType ht, enum ValueType vt, int i )
{

}*/

HashElement::~HashElement()
{
  next = previous = first = last = NULL;
  parent                         = NULL;
  delete valueD;
  delete valueI;
}

HashElement *
HashElement::getNext()
{
  return next;
}

HashElement *
HashElement::getPrev()
{
  return previous;
}

HashElement *
HashElement::getFirst()
{
  return first;
}

HashElement *
HashElement::getLast()
{
  return last;
}

void
HashElement::setNext(HashElement *ne)
{
  next         = ne;
  ne->previous = this;
  ne->first    = first;

  HashElement *tmp = first;

  while (tmp)
  {
    tmp->last = ne;
    tmp       = tmp->next;
  }
  delete tmp;
}

void
HashElement::setType(enum HashType ht)
{
  if (ht != HashType::undefined)
    return;
  if (ht == HashType::xyz || ht == HashType::rdc)
  {
    valueD = (double *) malloc(3 * sizeof(double));
    for (unsigned int i = 0; i < 3; ++i)
      valueD[i] = .0;
  }
  else if (ht == HashType::keyword)
  {
    valueI    = (int *) malloc(sizeof(int));
    valueI[0] = 0;
    valueD    = (double *) malloc(5 * sizeof(double));
    for (unsigned int i = 0; i < 5; ++i)
      valueD[i] = .0;
  }
  else
  {
    valueD = NULL;
    valueI = NULL;
  }
}

std::string
HashElement::setIdentifier(std::string s)
{
  //   if ( HashIdent == "" ) HashIdent = s;
  HashIdent = "";
  HashIdent = s;
  return HashIdent;
}

double *
HashElement::setValueD(int i, double v)
{
  if (!(valueD))
    return NULL;

  if (type == HashType::xyz || type == HashType::rdc)
  {
    if (i >= 3)
      return NULL;
    valueD[i] = v;
    return (valueD + i);
  }
  else if (type == HashType::keyword)
  {
    if (i >= 5)
      return NULL;
    valueD[i] = v;
    return (valueD + i);
  }
  return NULL;
}

double *
HashElement::setValuesD(int i, double *v)
{
  if (!(valueD))
    return NULL;

  if (type == HashType::xyz || type == HashType::rdc)
  {
    if (i != 3)
      return NULL;
    for (int j = 0; j < i; ++j)
      valueD[j] = v[j];
    return valueD;
  }
  else if (type == HashType::keyword)
  {
    if (i != 5)
      return NULL;
    for (int j = 0; j < i; ++j)
      valueD[j] = v[j];
    return valueD;
  }
  return NULL;
}

int *
HashElement::setValueI(int i, int v)
{
  if (!(valueI))
    return NULL;

  if (type == HashType::keyword)
  {
    //      if ( i != 0 ) return NULL;
    valueI[i] = v;
    return valueI;
  }
  return NULL;
}

int *
HashElement::setValuesI(int i, int *v)
{
  if (!(valueI))
    return NULL;

  if (type == HashType::keyword)
  {
    if (i != 0)
      return NULL;
    valueI[i] = v[i];
    return valueI;
  }
  return NULL;
}

String_ptr *
HashElement::setValueS(int i, std::string v)
{
  if (valueS)
    return valueS->setContent(i, v);
  else
  {
    valueS = new String_ptr();
    return valueS->setContent(i, v);
  }
  return NULL;
}

std::string
HashElement::getValueS(int i)
{
  if (valueS)
    return valueS->getContent(i);
  return "";
}

double *
HashElement::getValueD(int i)
{
  if (!(valueD))
    return NULL;

  if (type == HashType::xyz || type == HashType::rdc)
  {
    if (i >= 3)
      return NULL;
    return (valueD + i);
  }
  else if (type == HashType::keyword)
  {
    if (i >= 5)
      return NULL;
    return (valueD + i);
  }
  return NULL;
}

double *
HashElement::getValuesD()
{
  return valueD;
}

int *
HashElement::getValueI(int i)
{
  if (!(valueI))
    return NULL;

  if (type == HashType::keyword)
  {
    //      if ( i != 0 ) return NULL;
    return (valueI + i);
  }
  return NULL;
}

int *
HashElement::getValuesI()
{
  return valueI;
}


/**************************************
 *      Implementation of OwnHash     *
 **************************************/

OwnHash::OwnHash()
{
  type  = HashType::undefined;
  first = last = current = NULL;
}

OwnHash::OwnHash(enum HashType ht)
{
  type  = ht;
  first = last = current = NULL;
}

OwnHash::~OwnHash()
{
  HashElement *tmp = first;
  HashElement *de  = first;
  while (tmp)
  {
    tmp = tmp->getNext();
    delete de;
    de = tmp;
  }
}

HashElement *
OwnHash::addElement(std::string s, Keyword *k)
{
  if (this->getElement(s))
    return current;

  HashElement *tmp = new HashElement(type, s, k);
  if (first)
  {
    last->setNext(tmp);
    last = tmp;
  }
  else
  {
    first = last = tmp;
  }
  return tmp;
}

HashElement *
OwnHash::getElement(std::string s)
{
  if (current && !(current->getIdentifier().compare(s)))
    return current;

  current = first;

  while (current)
  {
    if (current->getIdentifier().compare(s))
      current = current->getNext();
    else
      return current;
  }

  return NULL;
}

HashElement *
OwnHash::getElement(int i)
{
  current = first;
  for (int j = 0; j < i; ++j)
  {
    if (current == NULL)
      return NULL;
    current = current->getNext();
  }
  return current;
}

std::vector<std::string>
OwnHash::getKeys()
{
  std::vector<std::string> k;
  current = first;
  while (current)
  {
    k.push_back(current->getIdentifier());
    current = current->getNext();
  }
  current = first;
  return k;
}

void
OwnHash::setValueD(std::string s, int i, double v)
{
  current = this->getElement(s);
  if (current)
    current->setValueD(i, v);
}

void
OwnHash::setValuesD(std::string s, int i, double *v)
{
  this->getElement(s);
  if (current)
    current->setValuesD(i, v);
}

void
OwnHash::setValueI(std::string s, int i, int v)
{
  this->getElement(s);

  if (current)
    current->setValueI(i, v);
}

void
OwnHash::setValuesI(std::string s, int i, int *v)
{
  this->getElement(s);
  if (current)
    current->setValuesI(i, v);
}

void
OwnHash::setValueS(std::string s, std::string v)
{
  this->getElement(s);
  int i = 0;
  while (getValueS(s, i) != "" && getValueS(s, i) != v)
    ++i;
  if (getValueS(s, i) == v)
    return;
  else
    setValueS(s, i, v);
}

void
OwnHash::setValueS(std::string s, int i, std::string v)
{
  //   String_ptr* tmp;// = NULL;
  this->getElement(s);
  if (current)
    current->setValueS(i, v);
}

double
OwnHash::getValueD(std::string s, int i)
{
  this->getElement(s);
  if (current)
    return *(current->getValueD(i));
  return .0;
}

double *
OwnHash::getValuesD(std::string s)
{
  this->getElement(s);
  if (current)
    return current->getValuesD();
  return NULL;
}

int
OwnHash::getValueI(std::string s, int i)
{
  this->getElement(s);
  if (current)
    return *(current->getValueI(i));
  return .0;
}

int *
OwnHash::getValuesI(std::string s)
{
  this->getElement(s);
  if (current)
    return current->getValuesI();
  return NULL;
}

std::string
OwnHash::getValueS(std::string s, int i)
{
  this->getElement(s);
  if (current)
    return current->getValueS(i);
  return "";
}

std::vector<std::string>
OwnHash::getValuesS(std::string s)
{
  std::vector<std::string> r;
  this->getElement(s);
  if (current)
  {
    unsigned int i = 0;
    while (current->getValueS(i) != "")
    {
      r.push_back(current->getValueS(i++));
    }
  }
  return r;
}

void
OwnHash::printrdcTable()
{
  std::cout << "start to print\n";
  current = first;
  while (current)
  {
    std::cout << current->getIdentifier();
    for (int i = 0; i < 3; ++i)
      std::cout << "  " << *current->getValueD(i);
    std::cout << std::endl;
    current = current->getNext();
  }
}

void
OwnHash::printCoordinates()
{
  std::cout << "start to print\n";
  current = first;
  while (current)
  {
    std::cout << current->getIdentifier();
    for (int i = 0; i < 3; ++i)
      std::cout << "  " << *current->getValueD(i);
    std::cout << std::endl;
    current = current->getNext();
  }
}

void
OwnHash::printConn()
{
  std::cout << "start to print connectivity\n";
  std::string out;
  current = first;
  while (current)
  {
    std::cout << current->getIdentifier();
    int i = 0;
    out   = "";
    while (true)
    {
      out = current->getValueS(i++);
      if (out == "")
        break;
      std::cout << "   " << out;
    }
    std::cout << std::endl;
    current = current->getNext();
  }
}

void
OwnHash::printKeys()
{
  std::cout << "start to print keywors\n";

  double *d;
  int *i;
  int j;
  std::string out;

  current = first;
  while (current)
  {
    j = 0;
    std::cout << current->getIdentifier();
    d = current->getValuesD();
    if (d)
    {
      for (int k = 0; k < 5; ++k)
        std::cout << "  " << d[k];
    }
    i = current->getValuesI();
    if (i)
    {
      std::cout << "  " << i[0];
    }
    while (true)
    {
      out = current->getValueS(j++);
      if (out == "")
        break;
      std::cout << "  " << out;
    }
    std::cout << std::endl;
    current = current->getNext();
  }
  d = NULL;
  i = NULL;
}
