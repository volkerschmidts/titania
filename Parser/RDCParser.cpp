#include <Parser/OwnHash.hpp>
#include <Parser/Parser.hpp>
#include <fstream>
#include <iostream>
#include <vector>

RDCinput::RDCinput()
{
  next = previous = this;
  rdcTable        = NULL; // new OwnHash(HashType::rdc);
  label           = "";
}

RDCinput::RDCinput(std::string s)
{
  next = previous = this;
  rdcTable        = NULL; // new OwnHash(HashType::rdc);
  label           = s;
}

RDCinput::RDCinput(std::string s, std::fstream &f)
{
  next = previous = NULL;
  label           = s;

  std::string il = "";

  while (std::getline(f, il))
  {
    rdcs.push_back(il);
    //      std::cout << il << std::endl;
  }
  rdcTable = this->initializeRDCtable();

  this->reorder();
}

RDCinput::RDCinput(std::string s, std::vector<std::string> v)
{
  next = previous = NULL;
  label           = s;

  rdcs     = v;
  rdcTable = this->initializeRDCtable();
}

RDCinput::~RDCinput()
{
  next = previous = NULL;
  delete rdcTable;
}

std::vector<std::string>
RDCinput::getKeys()
{
  if (rdcTable)
  {
    return rdcTable->getKeys();
  }
  std::vector<std::string> n;
  return n;
}

OwnHash *
RDCinput::initializeRDCtable()
{
  std::string ident    = "";
  OwnHash *table       = new OwnHash(HashType::rdc);
  HashElement *element = NULL;
  for (unsigned int i = 0; i < rdcs.size(); ++i)
  {
    ta2csvString(rdcs.at(i));
    if (rdcs.at(i).size() == 0)
    {
      rdcs.erase(rdcs.begin() + i--);
      continue;
    }
    int j = 0;
    while (rdcs.at(i).at(++j) != 59)
      ;
    while (rdcs.at(i).at(++j) != 59)
      ;
    ident          = rdcs.at(i).substr(0, j);
    element        = table->addElement(ident);
    unsigned int k = ++j;
    for (int e = 0; e < 3; ++e) // Read the 3 values
    {
      j = k;
      while ((k) < rdcs.at(i).size() && rdcs.at(i).at(k) != 59)
        ++k;
      if (rdcs.at(i).substr(j, (k - j)) == "u")
      {
        std::cout << "undefined RDC\n";
        break;
      }
      //         std::cout << atof(rdcs.at(i).substr(j,(k-j)).c_str()) <<
      //         std::endl;
      if (!(element->setValueD(e, atof(rdcs.at(i).substr(j, (k - j)).c_str()))))
      {
        std::cout << "\nTODO:\tImplement the parsing of external files to "
                     "generat error Messages\n";
      }
      ++k;
    }
  }

  return table;
}

RDCinput *
RDCinput::setNext(std::string s, std::fstream &f)
{
  next           = new RDCinput(s, f);
  next->previous = this;
  return next;
}

RDCinput *
RDCinput::append(std::string s, std::fstream &f)
{
  RDCinput *tmp = next;
  if (!next)
  {
    tmp = this;
  }
  else
    while (tmp->next)
      tmp = tmp->next;
  tmp->next = new RDCinput(s, f);
  this->reorder();
  tmp->next->previous = tmp;
  return tmp->next;
}

void
RDCinput::reorder()
{
  unsigned int i, j;
  std::vector<std::string> r = rdcTable->getKeys();
  std::string a1, a2, id;

  for (i = 0; i < r.size(); ++i)
  {
    j  = 0;
    a1 = a2 = "";
    id      = r.at(i);
    while (id.at(j) != ';')
      ++j;
    a1 = id.substr(0, j);
    a2 = id.substr((j + 1), id.size() - j - 1);
    j  = 0;
    while ((j < (a1.size() - 1) && j < (a2.size() - 1)) && a1.at(j) == a2.at(j))
      ++j;
    if (a2.at(j) < a1.at(j))
    {
      rdcTable->getElement(id)->setIdentifier((a2 + ";" + a1));
    }
  }
}

double *
RDCinput::getValuesD(std::string s)
{
  if (rdcTable)
    return rdcTable->getValuesD(s);
  return NULL;
}
