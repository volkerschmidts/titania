#include <Parser/OwnHash.hpp>
#include <Parser/Parser.hpp>
#include <fstream>
#include <iostream>
#include <vector>

StructureInput::StructureInput()
{
  NOA      = 0;
  xyzTable = new OwnHash(HashType::xyz);
  conTable = new OwnHash(HashType::connectivity);
  parent   = NULL;
}

StructureInput::~StructureInput()
{
  NOA = 0;
  if (xyzTable)
    delete xyzTable;
  if (conTable)
    delete conTable;
  parent = NULL;
}

std::vector<std::string>
StructureInput::getKeys()
{
  if (xyzTable)
  {
    return xyzTable->getKeys();
  }
  else if (conTable)
  {
    return conTable->getKeys();
  }
  std::vector<std::string> n;
  return n;
}

int
StructureInput::setAtomInformation(std::string &ident,
                                   double **c,
                                   unsigned int &NOB)
{
  if (!xyzTable)
    return 0;

  c[0] = xyzTable->getValuesD(ident);

  NOB = 0;

  if (!conTable->getElement(ident))
    return 0;
  while (conTable->getElement(ident)->getValueS(NOB) != "")
    ++NOB;
  return 1;
}

std::string
StructureInput::getBondPartner(std::string ident, unsigned int b)
{
  return conTable->getValueS(ident, b);
}

void
StructureInput::parseCon(std::vector<std::string> v, unsigned int &i)
{
  unsigned int j, k, l, tmpNOA;
  int e;
  std::string ident = "";
  // Parse molecule name
  j = 0;
  while (!(v.at(i).at(j) == 59 || v.at(i).at(j) == 91 || v.at(i).at(j) == 93))
    ++j;
  k = ++j;
  while ((j + 1) < v.size() &&
         !(v.at(i).at(j) == 59 || v.at(i).at(j) == 91 || v.at(i).at(j) == 93))
    ++j;
  molecule = v.at(i).substr(k, (j - k));

  // parse NOA

  tmpNOA = atoi(v.at(++i).c_str());
  if (!tmpNOA)
    --i; // If no number of atoms are defined go back to last line and continue

  // Parse label (from unformated vector in class InputFile
  label = parent->findLabel(v.at(++i), "connectivity");
  // Start to parse connectivity
  j = 0;
  while (true)
  {
    ++i;
    k = 0;
    e = -1;
    while (true) // one does not know about the amount of bond partners
    {
      l = k;
      while (k < v.at(i).size() && v.at(i).at(k) != 59)
        if (v.at(i).at(k++) == 29)
          break;
      if (k > v.at(i).size())
        break;
      if (e < 0) // Hash identifier
      {
        ident = v.at(i).substr(l, (k - l));
        conTable->addElement(ident);
        ++e;
      }
      else // bonds
      {
        std::string x = v.at(i).substr(l, (k - l));
        if (atoi(x.c_str()))
        {
          x = xyzTable->getElement(atoi(x.c_str()) - 1)->getIdentifier();
        }
        // Add bond partner to current atom
        conTable->setValueS(ident, x);

        // Add Bond partner to list of atoms
        conTable->addElement(x);

        // Add current atom as partner to bond partner
        conTable->setValueS(x, ident);
      }
      ++k;
    }
    if (tmpNOA && ++j == tmpNOA)
      break;
  }
}

void
StructureInput::parseXYZ(std::vector<std::string> v, unsigned int &i)
{
  unsigned int j, k, l;
  std::string ident = "";
  // Parse molecule name
  j = 0;

  while (!(v.at(i).at(j) == 59 || v.at(i).at(j) == 91 || v.at(i).at(j) == 93))
    ++j;
  k = ++j;
  while ((j + 1) < v.size() &&
         !(v.at(i).at(j) == 59 || v.at(i).at(j) == 91 || v.at(i).at(j) == 93))
    ++j;

  molecule = v.at(i).substr(k, (j - k));

  // parse NOA

  NOA = atoi(v.at(++i).c_str());
  if (!NOA)
    --i; // If no number of atoms are defined go back to last line and continue

  // Parse label (from unformated vector in class InputFile
  label = parent->findLabel(v.at(++i), "xyzcoordinates");

  // Start to parse coordinates
  j = 0;
  while (true)
  {
    ++i;
    k = 0;
    for (int e = -1; e < 3; ++e)
    {
      l = k;
      while ((k + 1) < v.at(i).size() && v.at(i).at(++k) != 59)
        ;
      if (e < 0) // Hash identifier
      {
        ident = v.at(i).substr(l, (k - l));
        xyzTable->addElement(ident);
      }
      else // coordinates
      {
        double x = atof(v.at(i).substr(l, (k - l)).c_str());
        xyzTable->setValueD(ident, e, x);
      }
      ++k;
    }
    if (NOA && ++j == NOA)
      break;
  }
}

int
StructureInput::parseMOLfile(std::string &s, std::string &workingDir)
{
  unsigned int i, j, l;
  int k;

  // First parse the filename

  i = s.size() - 1;
  l = 0;
  if (s.at(i) == 59)
    --i;
  while (s.at(i--) != 59)
    ++l;
  i += 2;

  std::string fn = workingDir + s.substr(i, l);
  std::string il;
  std::fstream f(fn, std::ios::in | std::ios::binary);

  if (f.good())
  {
    std::getline(f, molecule); // Parse molecule name
    std::getline(
        f,
        label); // Parse label (obtained from the unique structure identifier)
    std::getline(f, il); // Skip comment line
    std::getline(f, il); // Parse counts line
    ta2csvString(
        il); // Trunkate the counts line to obtain easier to parse string
    //      std::cout << "Molecule: " << molecule << "\nStructure: " << label <<
    //      std::endl; std::cout << il << std::endl;

    i = il.size() - 1;
    l = 1;
    while (il.at(--i) != 59)
      ++l;

    // Check for the right file version (V2000)
    if (il.substr(++i, 2).compare(MOLFILE_VERSION))
    {
      std::cerr << "ERROR:\tThe molfile version (" << il.substr(i, l)
                << ") does not match the implemented one (" << MOLFILE_VERSION
                << "x)\n"
                << "\tTITANIA will be shot down. For more information see the "
                   "manual.\n";
      f.close();
      f.clear();
      return 1;
    }

    // Parse the needed counts (Atoms, Bonds)
    int NOB;
    i = 0;
    l = 0;

    while (il.at(i++) != 59)
      ++l;
    NOA = atoi(il.substr(0, l).c_str());
    l   = 0;
    j   = i;
    while (il.at(i++) != 59)
      ++l;
    NOB = atoi(il.substr(j, l).c_str());
    //      std::cout << "NumberOfAtom: " << NOA << "\nNumberOfBonds: " << NOB
    //      << std::endl;

    unsigned int a;
    int b;
    double *coord = (double *) malloc(3 * sizeof(double));
    std::string ident;

    for (a = 0; a < NOA; ++a)
    {
      if (!std::getline(f, il))
      {
        std::cerr << "ERROR:\tThe molfile does not match standards or is "
                     "corrupted. Check the file!\n"
                  << "\tTITANIA will be shot down. For more information see "
                     "the manual.\n";
        f.close();
        f.clear();
        return 1;
      }
      ta2csvString(il);
      i = 0;
      for (k = 0; k < 3; ++k) // Parse the coodinates
      {
        l = 0;
        j = i;
        while (il.at(i++) != 59)
          ++l;

        coord[k] = atof(il.substr(j, l).c_str());
      }
      l = 0;
      j = i;
      while (il.at(i++) != 59)
        ++l;
      ident = il.substr(j, l) + std::to_string((a + 1));

      xyzTable->addElement(ident);
      xyzTable->setValuesD(ident, 3, coord);
    }
    delete coord;
    // For the most cases the current idents are wrong! So one should first skip
    // the connectivity part and read all proper labels
    std::streampos constart = f.tellg();
    for (b = 0; b < NOB; ++b)
      std::getline(f, il);


    std::string realL;
    HashElement *atom;

    // Now check for labels
    while (std::getline(f, il))
    {
      ta2csvString(il);

      i = il.size() - 1;
      l = 1;
      while (il.at(--i) != 59)
        ++l; // Parse if there is a to update the label
      a = atoi(il.substr(++i, l).c_str()); // Get the atom number

      if (!std::getline(f, realL))
        break;             // Parse the real label
      ta2csvString(realL); // = il;

      atom = xyzTable->getElement(
          (a - 1)); // Hash table is 0 based, atom numbers 1 based

      if (atom)
        realL = atom->setIdentifier(realL); // and update the identifier
    }
    f.clear();
    f.seekg(
        constart); // Since eof might be reached clear and jump to connectivity

    // Worst thing is done!

    int con[4];
    HashElement *bondA;
    for (b = 0; b < NOB; ++b)
    {
      std::getline(f, il);
      ta2csvString(il);
      i = 0;
      //         std::cout << il << std::endl;
      for (k = 0; k < 4; ++k)
      {
        j = i;
        l = 1;
        while (++i < il.size() && il.at(i) != 59)
          ++l; // Parse if there is a to update the label
        con[k] = atoi(il.substr(j, l).c_str()); // Get the atom number
        ++i;
      }
      atom  = xyzTable->getElement(con[0] - 1);
      bondA = xyzTable->getElement(con[1] - 1);
      if (!(atom && bondA))
      {
        std::cerr << "ERROR:\tThe atom indizes in the bond array of the "
                     "provided .mol file are out of bounds.\n"
                  << "\tCheck the file or provide the structure by hand.\n";
        return 1;
      }

      conTable->addElement(atom->getIdentifier());
      conTable->addElement(bondA->getIdentifier());

      for (k = 0; k < con[2]; ++k)
      {
        i = conTable->getValuesS(atom->getIdentifier()).size();
        conTable->setValueS(atom->getIdentifier(), i, bondA->getIdentifier());
        i = conTable->getValuesS(bondA->getIdentifier()).size();
        conTable->setValueS(bondA->getIdentifier(), i, atom->getIdentifier());
      }
    }
    atom  = NULL;
    bondA = NULL;
  }
  return 0;
}

int
StructureInput::checkDiff(std::string &key, std::string &missing)
{
  std::vector<std::string> cKeys;
  std::vector<std::string> xKeys;

  if (conTable)
  {
    cKeys = conTable->getKeys();
  }
  if (xyzTable)
  {
    xKeys = xyzTable->getKeys();
  }
  if (xKeys.size() == 0 || cKeys.size() == 0)
  {
    // TODO write proper error message
    return 1;
  }
  unsigned int i;

  while (true)
  {
    for (i = 0; i < xKeys.size(); ++i)
    {
      if (!xKeys.at(i).compare(cKeys.at(0)))
      {
        xKeys.erase(xKeys.begin() + i);
        cKeys.erase(cKeys.begin());
        break;
      }
    }
    if (i == xKeys.size())
      break;
    if (xKeys.size() == 0 || cKeys.size() == 0)
      break;
  }

  if (cKeys.at(0).compare(key))
    missing = cKeys.at(0);
  return 0;
}

void
StructureInput::printCoordinates()
{
  xyzTable->printCoordinates();
}

void
StructureInput::printConnectivity()
{
  conTable->printConn();
}
