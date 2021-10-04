#include <Parser/OwnHash.hpp>
#include <Parser/Parser.hpp>
#include <fstream>
#include <iostream>
#include <vector>

// Keyword rCount = { "ringCount", 1, ValueType::intType };
Keyword aCount = {"alignCount", 1, ValueType::intType};

KeywordInput::KeywordInput()
{
  keyTable = new OwnHash(HashType::keyword);
  parent   = NULL;
}

KeywordInput::KeywordInput(std::string s)
{
  keyTable = new OwnHash(HashType::keyword);
  ta2csvString(s);
  keys.push_back(s);
  parent = NULL;
}

KeywordInput::KeywordInput(std::vector<std::string> v)
{
  keyTable      = new OwnHash(HashType::keyword);
  std::string s = "";
  for (unsigned int i = 0; i < v.size(); ++i)
  {
    s = v.at(i);
    ta2csvString(s);
    if (s.size() == 0)
      continue;
    keys.push_back(s);
  }
  parent = NULL;
}

KeywordInput::KeywordInput(std::fstream &f)
{
  std::string il = "";
  keyTable       = new OwnHash(HashType::keyword);

  while (std::getline(f, il))
  {
    ta2csvString(il);
    if (il.size() == 0)
      continue;
    keys.push_back(il);
  }
  parent = NULL;
}

KeywordInput::~KeywordInput()
{
  parent = NULL;
}

std::vector<std::string>
KeywordInput::getKeys()
{
  if (keyTable)
  {
    return keyTable->getKeys();
  }
  std::vector<std::string> n;
  return n;
}

int
KeywordInput::getValueI(std::string s, int i)
{
  if (keyTable)
    return keyTable->getValueI(s, i);
  return 0;
}

std::string
KeywordInput::getValueS(std::string s, int i)
{
  if (keyTable)
    return keyTable->getValueS(s, i);
  return "";
}

std::vector<std::string>
KeywordInput::getValuesS(std::string s)
{
  if (keyTable)
    return keyTable->getValuesS(s);
  std::vector<std::string> null;
  return null;
}

double *
KeywordInput::getValuesD(std::string s)
{
  if (keyTable)
    return keyTable->getValuesD(s);
  return NULL;
}

double
KeywordInput::getValueD(std::string s, int i)
{
  if (keyTable)
    return keyTable->getValueD(s, i);
  return .0;
}

bool
KeywordInput::contains(std::string s)
{
  if (keyTable->getElement(s))
    return true;
  else
    return false;
}

void
KeywordInput::addKeyword(std::string s, Keyword k)
{
  static int aligncount = 1;
  unsigned int i, j, l;
  bool check_argument_num = false;
  std::string ident       = "";
  std::string arg         = "";
  ta2csvString(s);
  if (s.size() == 0)
    return;
  for (j = 0; (j < s.size() && s.at(j) != 59); ++j)
    s.at(j) = std::tolower(s.at(j));

  keys.push_back(s);
  // Inititalize the hash element

  i = j = l = 0;
  while (j < s.size() && (s.at(j) != 59 && s.at(j) != 91 && s.at(j) != 93))
    ++j;

  ident = s.substr(l, j);

  if (ident == "alignment")
  {
    if (aligncount == 1)
      keyTable->addElement(aCount.keyword, &aCount);
    keyTable->setValueI(aCount.keyword, 0, aligncount);
    ident += std::to_string(aligncount++);
  }
  else if (ident == "predictrdcs")
    check_argument_num = true;
  keyTable->addElement(ident, &k);

  // Parse the arguments and save them to hash table

  for (i = 0; i < k.Arguments; ++i)
  {
    l = ++j;

    while (j != s.size() && (s.at(j) != 59 && s.at(j) != 91 && s.at(j) != 93))
      ++j;

    arg = s.substr(l, (j - l));

    if (!(atoi(arg.c_str()) || atof(arg.c_str())))
      continue;
    if (k.valueType == ValueType::undefined)
      break;
    else if (k.valueType == ValueType::intType)
      keyTable->setValueI(ident, i, atoi(arg.c_str()));
    else if (k.valueType == ValueType::doubleType)
      keyTable->setValueD(ident, i, atof(arg.c_str()));
    if (check_argument_num)
    {
      check_argument_num = false;
      if (keyTable->getValueI(ident, i) == 0)
      {
        break;
      }
    }
  }
  if (k.valueType == ValueType::stringType)
  {
    i = 0;

    while (true)
    {
      l = ++j;

      while (j < s.size() && (s.at(j) != 59 && s.at(j) != 91 && s.at(j) != 93))
        ++j;

      arg = s.substr(l, (j - l));
      keyTable->setValueS(ident, i++, arg);
      if (j == s.size())
        break;
    }
  }
}

void
KeywordInput::removeEqSgn()
{
  for (unsigned int i = 0; i < keys.size(); ++i)
  {
    for (unsigned int j = 0; j < keys.at(i).size(); ++j)
    {
      if (keys.at(i).at(j) == 61)
        keys.at(i).at(j) = 59;
    }
    ta2csvString(keys.at(i));
  }
}
