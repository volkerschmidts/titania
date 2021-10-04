#include <Declarations.hpp>
#include <Parser/OwnHash.hpp>
#include <Parser/Parser.hpp>
#include <fstream>
#include <iostream>
#include <vector>


// Never delete Keywords! one can deactivate them in Molecule.cpp
// Never change the numbers/positions of the keywords in the list!
struct Keyword ListOfKeywords[] = {
    /* 000 */ {"undefined", 0, ValueType::undefined},
    /* 001 */ {"echo", 1, ValueType::intType},
    /* 002 */ {"xyzcoordinates", 0, ValueType::undefined},
    /* 003 */ {"connectivity", 0, ValueType::undefined},
    /* 004 */ {"molfile", 0, ValueType::undefined},
    /* 005 */ {"experimentalrdcs", 0, ValueType::undefined},
    /* 006 */ {"csv2rdc", 0, ValueType::undefined},
    /* 007 */ {"include", 0, ValueType::stringType},
    /* 008 */ {"lmmaxiterations", 1, ValueType::intType},
    /* 009 */ {"maxtitaniaiterations", 1, ValueType::intType},
    /* 010 */ {"lmconvergencelimit", 1, ValueType::doubleType}, // deactivated
    /* 011 */ {"skipscrm", 1, ValueType::intType},
    /* 012 */ {"predictrdcs", 2, ValueType::intType},
    /* 013 */ {"alignment", 5, ValueType::doubleType},
    /* 014 */ {"ring", 0, ValueType::stringType}, // deactivated
    /* 015 */ {FILESTART_HEX, 0, ValueType::undefined},
    /* 016 */ {FILEEND_HEX, 0, ValueType::undefined},
    /* 017 */ {"silent", 0, ValueType::stringType},
    /* 018 */ {"qfactorconvergence", 1, ValueType::doubleType},
    /* 019 */ {"titania2hotfcht", 1, ValueType::intType},
    /* 020 */ {"outputali", 1, ValueType::intType},
    /* 021 */ {"threads", 1, ValueType::intType},
    /* 022 */ {"montecarlobootstrapping", 1, ValueType::intType},
    /* 023 */ {"plotkappaq", 1, ValueType::intType},
    /* 024 */ {"gnuoutputformat", 0, ValueType::stringType},
    /* 025 */ {"numericalgradients", 1, ValueType::intType},
    /* 026 */ {"outputlm", 1, ValueType::intType},
    /* 027 */ {"redundantsconvergence", 1, ValueType::doubleType},
    /* 028 */ {"meanalignmentconvergence", 1, ValueType::doubleType},
    /* 029 */ {"sigmaalignmentconvergence", 1, ValueType::doubleType},
    /* 030 */ {"meanangleconvergence", 1, ValueType::doubleType},
    /* 031 */ {"sigmaangleconvergence", 1, ValueType::doubleType},
    /* 032 */ {"spreadangleconvergence", 1, ValueType::doubleType},
    /* 033 */ {"plottrajectory", 1, ValueType::intType},
    /* 034 */ {"recalculaterdcs", 1, ValueType::intType},
    /* 035 */ {"printwarnings", 1, ValueType::intType},
    /* 036 */ {"lmoptions", 4, ValueType::doubleType},
    /* 037 */ {"lmeckartoptions", 4, ValueType::doubleType},
    /* 038 */ {"lmsphericalsoptions", 4, ValueType::doubleType},
    /* 039 */ {"skipeckart", 1, ValueType::intType},
    /* 040 */ {"printredundants", 1, ValueType::intType},
    /* 041 */ {"titania2titania", 1, ValueType::intType},
    /* 042 */ {"normchiralvolumes", 1, ValueType::intType},
    /* 043 */ {"overoptimizationsteps", 1, ValueType::intType},
    /* 044 */ {"bigdata", 1, ValueType::intType},
    /* 045 */ {"plotrdcdynamics", 1, ValueType::intType},
    /* 046 */ {"scalewithsoverall", 1, ValueType::intType},
    /* 047 */ {"montecarlooutput", 1, ValueType::intType},
    /* 048 */ {"plotmontecarlo", 1, ValueType::intType},
    /* 049 */ {"calculatefullmatrix", 1, ValueType::intType},
    /* 050 */ {"errorweightinsvd", 1, ValueType::intType},
    /* 051 */ {"useredundantsonlyafter", 1, ValueType::intType},
    /* 052 */ {"redundantcycles", 1, ValueType::intType},
    /* 053 */ {"redundantsvalidity", 1, ValueType::doubleType},
    /* 054 */ {"redundantsdamping", 1, ValueType::doubleType},
    /* 055 */ {"usegpu", 1, ValueType::intType},
    /* 056 */ {"usedistances", 1, ValueType::intType},
    /* 057 */ {"inversionstrictness", 1, ValueType::intType},
    /* 058 */ {"memorypurge", 1, ValueType::doubleType},
    /* 059 */ {"staticbondweighting", 1, ValueType::doubleType},
    /* 060 */ {"staticangleweighting", 1, ValueType::doubleType},
    /* 061 */ {"statictorsionweighting", 1, ValueType::doubleType},
    /* 062 */ {"staticrdcweighting", 1, ValueType::doubleType},
    /* 063 */ {"staticchiralvolumeweighting", 1, ValueType::doubleType},
    /* 064 */ {"staticdistanceweighting", 1, ValueType::doubleType},
    /* 065 */ {"floatingrdcangles", 1, ValueType::intType},
    /* 066 */ {"useinitialholonomics", 1, ValueType::intType},
    /* 067 */ {"usereducedcovariance", 1, ValueType::intType},
};

unsigned int NumberOfKeywords = sizeof(ListOfKeywords) / sizeof(Keyword);

InputFile::InputFile(std::fstream &f, std::string WD)
{
  workingDir     = WD;
  std::string il = "";
  while (std::getline(f, il))
  {
    unformated.push_back(il);
    mainFile.push_back(il);
    ta2csvString(il);

    if (il.size() == 0)
      continue;
    else if (il.at(0) == 35)
    {
      commentary.push_back(il);
      continue;
    }
    content.push_back(il);
    if (checkLine(il) == 7)
    {
      unsigned int i = content.size() - 1;
      this->includeFile(i);
    }
    if (content.back().size() == 0)
    {
      content.erase(content.begin() + content.size());
    }
  }

  structure = new StructureInput();
  rdcs      = NULL;
  keywords  = new KeywordInput();
  structure->setParent(this);
  keywords->setParent(this);
}

int
InputFile::checkFile()
{
  this->removeSpecialChars();
  int kw;
  kw = 0;
  unsigned int i;
  static int unknowns = 0;

  for (i = 0; i < content.size(); ++i)
  {
    //      std::cout << "in check file: " << content.at(i) << std::endl;
    kw = checkLine(content.at(i));

    switch (kw)
    {
      case 0: // { "undefined" },
      {
        unsigned int uk    = 0;
        std::string errMsg = "";
        std::string tmp    = "";
        while (uk != content.at(i).size() && content.at(i).at(uk) != 59)
          ++uk;
        tmp = content.at(i).substr(0, uk);
        std::cout << content.at(i) << std::endl;
        this->getKeyLine(tmp, errMsg, uk);
        std::cout << "WARNING:\tUnknown keyword \"" << tmp
                  << "\" in the input file. The respective line will be "
                     "skipped...\n\t"
                  << errMsg;
        if (++unknowns >= 3)
          return 2;
        break;
      }
      case 1:
        keywords->addKeyword(content.at(i), ListOfKeywords[kw]);
        break; // { "echo" },
      case 2:  // { "xyzcoordinates" },
      {
        structure->parseXYZ(content, i);
        break;
      }
      case 3: // { "connectivity" },
      {
        structure->parseCon(content, i);
        break;
      }
      case 4:
        if (structure->parseMOLfile(content.at(i), workingDir))
          return 1;
        break; // { "molfile" },
      case 5:
        break; // { "experimentalrdcs" },
      case 6:  // { "csv2rdc" },
      {
        std::string s        = content.at(i);
        std::string label    = "";
        std::string filename = "";
        int j, k;
        j = s.size();
        while (true)
        {
          k = 0;
          while (s.at(--j) != 59)
            ++k;
          filename = s.substr(j + 1, k);
          if (filename.at(0) != 47)
            filename = workingDir + filename;
          if (filename.at(0) != 35)
            break;
        }
        k = 0;
        --j;
        while (--j >= 0 && s.at(j) != 91)
          ++k; // Labels are given in []-brackets
        label = s.substr(j + 1, k);

        std::fstream f(filename, std::ios::in | std::ios::binary);
        if (!f.good())
        {
          std::cerr << "ERROR:\tFile " << filename
                    << " does not exist. Terminating TITANIA.\n";
          unsigned int uk    = 0;
          std::string errMsg = "";
          std::string tmp    = "";
          this->getKeyLine(filename, errMsg, uk);
          std::cerr << errMsg;
          return 1;
        }
        if (rdcs)
          rdcs->append(label, f);
        else
          rdcs = new RDCinput(label, f);
        f.close();
        break;
      }
      case 7: // if ( this->includeFile(i) ) return ERROR_FILE_DOES_NOT_EXIST;
              // // { "includefile" },
        keywords->addKeyword(content.at(i), ListOfKeywords[kw]);
        break;
      case 8:
        keywords->addKeyword(content.at(i), ListOfKeywords[kw]);
        break; // { "lmmaxiterations" },
      case 9:
        keywords->addKeyword(content.at(i), ListOfKeywords[kw]);
        break; // { "maxtitaniaiterations" },
      case 10:
        keywords->addKeyword(content.at(i), ListOfKeywords[kw]);
        break; // { "lmconvergencelimit" },
      case 11:
        keywords->addKeyword(content.at(i), ListOfKeywords[kw]);
        break; // { "skipscrm" },
      case 12:
        keywords->addKeyword(content.at(i), ListOfKeywords[kw]);
        break; // { "predictrdcs" },
      case 13:
        keywords->addKeyword(content.at(i), ListOfKeywords[kw]);
        break; // { "alignment" },
      case 14:
        break; // deactivated keywords->addKeyword ( content.at(i),
               // ListOfKeywords[kw] ); break;  // { "ring" },
      case 15:
        break; // { FILESTART_HEX },
      case 16:
        break; // { FILEEND_HEX },
      case 17:
        keywords->addKeyword(content.at(i), ListOfKeywords[kw]);
        break; // { "silent" },
      case 18:
        keywords->addKeyword(content.at(i), ListOfKeywords[kw]);
        break; // { "qconvergence" },
      case 19:
        keywords->addKeyword(content.at(i), ListOfKeywords[kw]);
        break; // { "titania2hotfcht" },
      case 20:
        keywords->addKeyword(content.at(i), ListOfKeywords[kw]);
        break; // { "outputali" },
      case 21:
        keywords->addKeyword(content.at(i), ListOfKeywords[kw]);
        break; // { "threads" },
      case 22:
        keywords->addKeyword(content.at(i), ListOfKeywords[kw]);
        break; // { "montecarlobootstrap" },
      case 23:
        keywords->addKeyword(content.at(i), ListOfKeywords[kw]);
        break; // { "plotKappaQ" },
      case 24:
        keywords->addKeyword(content.at(i), ListOfKeywords[kw]);
        break; // { "gnuoutputformat" },
      case 25:
        keywords->addKeyword(content.at(i), ListOfKeywords[kw]);
        break; // { "numericalgradients" },
      case 26:
        keywords->addKeyword(content.at(i), ListOfKeywords[kw]);
        break; // { "outputlm" },
      case 27:
        keywords->addKeyword(content.at(i), ListOfKeywords[kw]);
        break; // { "redundantsconvergence" },
      case 28:
        keywords->addKeyword(content.at(i), ListOfKeywords[kw]);
        break; // { "meanalignmentconvergence" },
      case 29:
        keywords->addKeyword(content.at(i), ListOfKeywords[kw]);
        break; // { "sigmaalignmentconvergence" },
      case 30:
        keywords->addKeyword(content.at(i), ListOfKeywords[kw]);
        break; // { "meanangleconvergence" },
      case 31:
        keywords->addKeyword(content.at(i), ListOfKeywords[kw]);
        break; // { "sigmaangleconvergence" },
      case 32:
        keywords->addKeyword(content.at(i), ListOfKeywords[kw]);
        break; // { "spreadangleconvergence" },
      case 33:
        keywords->addKeyword(content.at(i), ListOfKeywords[kw]);
        break; // { "plottrajectory" }
      case 34:
        keywords->addKeyword(content.at(i), ListOfKeywords[kw]);
        break; // { "recalculaterdcs" }
      case 35:
        keywords->addKeyword(content.at(i), ListOfKeywords[kw]);
        break; // { "printwarnings" }
      case 36:
        keywords->addKeyword(content.at(i), ListOfKeywords[kw]);
        break; // { "lmoptions" }
      case 37:
        keywords->addKeyword(content.at(i), ListOfKeywords[kw]);
        break; // { "lmeckartoptions" }
      case 38:
        keywords->addKeyword(content.at(i), ListOfKeywords[kw]);
        break; // { "lmsphericalsoptions" }
      case 39:
        keywords->addKeyword(content.at(i), ListOfKeywords[kw]);
        break; // { "skipeckart" }
      case 40:
        keywords->addKeyword(content.at(i), ListOfKeywords[kw]);
        break; // { "printredundants" },
      case 41:
        keywords->addKeyword(content.at(i), ListOfKeywords[kw]);
        break; // { "titania2titania" },
      case 42:
        keywords->addKeyword(content.at(i), ListOfKeywords[kw]);
        break; // { "normchiralvolumes" },
      case 43:
        keywords->addKeyword(content.at(i), ListOfKeywords[kw]);
        break; // { "overoptimizationsteps" },
      case 44:
        keywords->addKeyword(content.at(i), ListOfKeywords[kw]);
        break; // { "bigdata" },
      case 45:
        keywords->addKeyword(content.at(i), ListOfKeywords[kw]);
        break; // { "plotrdcdynamics" },
      case 46:
        keywords->addKeyword(content.at(i), ListOfKeywords[kw]);
        break; // { "scalewithsovrall" },
      case 47:
        keywords->addKeyword(content.at(i), ListOfKeywords[kw]);
        break; // { "montecarlooutput" },
      case 48:
        keywords->addKeyword(content.at(i), ListOfKeywords[kw]);
        break; // { "plotmontecarlo" },
      case 49:
        keywords->addKeyword(content.at(i), ListOfKeywords[kw]);
        break; // { "calculatefullmatrix" },
      case 50:
        keywords->addKeyword(content.at(i), ListOfKeywords[kw]);
        break; // { "errorweightinsvd" },
      case 51:
        keywords->addKeyword(content.at(i), ListOfKeywords[kw]);
        break; // { "useredundantsonlyafter" },
      case 52:
        keywords->addKeyword(content.at(i), ListOfKeywords[kw]);
        break; // { "redundantcycles" },
      case 53:
        keywords->addKeyword(content.at(i), ListOfKeywords[kw]);
        break; // { "redundantsvalidity" },
      case 54:
        keywords->addKeyword(content.at(i), ListOfKeywords[kw]);
        break; // { "redundantsdamping" },
      case 55:
        keywords->addKeyword(content.at(i), ListOfKeywords[kw]);
        break; // { "usegpu" },
      case 56:
        keywords->addKeyword(content.at(i), ListOfKeywords[kw]);
        break; // { "usedistances" },
      case 57:
        keywords->addKeyword(content.at(i), ListOfKeywords[kw]);
        break; // { "inversionstrictness" },
      case 58:
        keywords->addKeyword(content.at(i), ListOfKeywords[kw]);
        break; // { "memorypurge" },
      case 59:
        keywords->addKeyword(content.at(i), ListOfKeywords[kw]);
        break; // { "staticbondweighting" },
      case 60:
        keywords->addKeyword(content.at(i), ListOfKeywords[kw]);
        break; // { "staticangleweighting" },
      case 61:
        keywords->addKeyword(content.at(i), ListOfKeywords[kw]);
        break; // { "statictorsionweighting" },
      case 62:
        keywords->addKeyword(content.at(i), ListOfKeywords[kw]);
        break; // { "staticrdcweighting" },
      case 63:
        keywords->addKeyword(content.at(i), ListOfKeywords[kw]);
        break; // { "staticchiralvolumeweighting" },
      case 64:
        keywords->addKeyword(content.at(i), ListOfKeywords[kw]);
        break; // { "staticdistanceweighting" },
      case 65:
        keywords->addKeyword(content.at(i), ListOfKeywords[kw]);
        break; // { "floatingrdcangles" },
      case 66:
        keywords->addKeyword(content.at(i), ListOfKeywords[kw]);
        break; // { "useinitialholonomics" },
      case 67:
        keywords->addKeyword(content.at(i), ListOfKeywords[kw]);
        break; // { "usereducedcovariance" },
      default:
        std::cout << "TITANIA did not save: " << content.at(i) << std::endl;
    }
  }

  return 0;
}

void
InputFile::removeSpecialChars()
{
  char x;
  unsigned int i, j;
  for (i = 0; i < content.size(); ++i)
  {
    //       std::cout << "start in remove special char: " << i << content.at(i)
    //       << std::endl;
    for (j = 0; j < content.at(i).size(); ++j)
    {
      x = content.at(i).at(j);
      if (x == 35)
        break;                       // #
      else if (x < 44 ||             // Special, \" !"#$%&'()*+ \"
               x == 58 ||            // :
               (x > 59 && x < 65) || // \" <=>?@ \"
               (x > 93 && x < 97) || // \" ^_` \"
               x > 122               // \" {|}~DEL \"
      )
      {
        content.at(i).at(j) = 59;
      }                 // ;
      else if (x == 91) // [
      {
        //            content.at(i).at(j) = 59;             // ;
        while (content.at(i).at(++j) != 93)
          ; // ]
            //            content.at(i).at(j) = 59;             // ;
      }
      else if (x == 92 || x == 47) // \ /
      {
        while (j < content.at(i).size() && content.at(i).at(j) != 59)
          ++j;
      }
    }
    //      std::cout << "start in remove special char: " << content.at(i) <<
    //      std::endl;
    ta2csvString(content.at(i));
  }
}

void
ta2csvString(std::string &s)
{
  //   std::cout << "in ta2csvString: " << s << std::endl;
  int i, j;
  std::string tmp;
  for (i = ((int) s.size() - 1); i >= 0; --i)
  {
    if (!isspace(s.at(i)))
      break;
  }

  for (j = 0; j < ((int) s.size()); ++j)
  {
    if (!isspace(s.at(j)))
      break;
  }
  tmp = s.substr(j, i - j + 1);
  s   = tmp;
  for (j = 0; j < ((int) s.size()); ++j)
  {
    if (s.at(j) == 35)
      break;
    if (s.at(j) == 44)
    {
      s.at(j) = 46;
      continue;
    } // replace comma by point
    else if (isspace(s.at(j)))
    {
      int k = j;
      while ((k + 1) < ((int) s.size()) && isspace(s.at(k + 1)))
        ++k;
      if (j == k)
      {
        s.at(j) = 59;
      }
      else
        s.replace(j, (k - j + 1), ";");
    }
  }
  // second parse to remove possible double semicolons due to bad .csv inputs

  for (j = 0; j < ((int) s.size()); ++j)
  {
    if (s.at(j) == 59 && (j + 1) < ((int) s.size()) && s.at(j + 1) == 59)
    {
      unsigned int k = j + 1;
      while ((k + 1) < s.size() && s.at(k + 1) == 59)
        ++k;
      s.replace(j, (k - j + 1), ";");
    }
  }
}

std::string
InputFile::findLabel(std::string s, std::string p)
{
  std::vector<std::string> keys;
  unsigned int i, j;
  i = 0;

  // Collect label fragments
  while (i < s.size())
  {
    j = i;
    while (i < s.size() && s.at(i++) != 59)
      ;
    if (i == s.size())
      ++i;
    keys.push_back(s.substr(j, (i - j - 1)));
  }
  if (keys.size() == 1)
    return keys.at(0);

  // Jump to current keyword
  std::string uf = "";
  for (i = 0; i < unformated.size(); ++i)
  {
    uf = unformated.at(i);
    for (j = 0; j < uf.size(); ++j)
      uf.at(j) = std::tolower(uf.at(j));
    if (uf.find(p) != std::string::npos)
      break; // std::cout << unformated.at(i) << std::endl;
  }
  // Find the label after Keyword
  for (; i < unformated.size(); ++i)
  {
    for (j = 0; j < keys.size(); ++j)
    {
      if (unformated.at(i).find(keys.at(j)) != std::string::npos)
        continue;
      break;
    }
    if (j == keys.size())
      break;
  }

  return unformated.at(i);
}

int
InputFile::includeFile(unsigned int &i)
{
  // Get the filename

  std::string cont = content.at(i);

  unsigned int j, k, l;
  k = cont.size() - 1;
  l = 0;
  while (cont.at(--k) != 59)
    ++l;

  // Link and check the file

  std::string basename = cont.substr(++k, ++l);

  std::string filename;
  if (basename.at(0) == 47)
    filename = basename;
  else
    filename = workingDir + basename;

  std::fstream inc(filename, std::ios::in | std::ios::binary);

  if (!inc.good())
  {
    std::cerr << "ERROR:\tFile " << filename
              << " does not exist. Terminating TITANIA.\n";
    unsigned int uk    = 0;
    std::string errMsg = "";
    std::string tmp    = "";
    this->getKeyLine(filename, errMsg, uk);
    std::cerr << errMsg;
    return ERROR_FILE_DOES_NOT_EXIST;
  }

  // Set the iterator for content and set "FileEnd" marker

  std::vector<std::string>::iterator it;
  k = i;
  if ((i + 1) == content.size())
  {
    content.push_back(FILEEND_HEX);
  }
  else
  {
    ++i;
    content.insert(content.begin() + i, FILEEND_HEX);
  }
  it = (content.begin() + i);

  // Parse the file and save content

  std::string il;

  // content, unformated and commentary have to be updated.
  // start with loading the two other iterators

  std::vector<std::string>::iterator u = unformated.begin();
  std::vector<std::string>::iterator c = commentary.begin();
  bool ac                              = false;
  bool au                              = false;

  for (l = 0, j = 0; l < unformated.size(); ++l)
  {
    if (unformated.at(l).size() == 0)
      continue;
    if (unformated.at(l).at(0) == 35)
    {
      ++j;      // jump to next commentary line;
      continue; // this line cant be the include line
    }
    if (unformated.at(l).find(basename) != (std::string::npos))
      break;
  }

  if ((j) == commentary.size())
    ac = true;
  else if (commentary.size() == 0)
    ac = true;
  else
    c = commentary.begin() + (j);
  if ((l + 1) == unformated.size())
    au = true;
  else
    u = unformated.begin() + (++l);
  while (std::getline(inc, il))
  {
    if (au)
      unformated.push_back(il);
    else
    {
      unformated.insert(u, il);
      u = (unformated.begin() + (++l));
    }

    ta2csvString(il);

    if (il.size() == 0)
    {
      continue;
    }
    else if (il.at(0) == 35)
    {
      if (ac)
        commentary.push_back(il);
      else
      {
        commentary.insert(c, il);

        c = (commentary.begin() + (++j));
      }
      continue;
    }

    content.insert(it, il);
    it = (content.begin() + (++i));
  }

  it = content.begin() + k;
  content.insert(it, FILESTART_HEX);

  this->removeSpecialChars();

  i = k + 1;

  return 0;
}

int
checkLine(std::string s)
{
  unsigned int i;
  i = 0;
  while (i < s.size() && !(s.at(i) == 59 || s.at(i) == 91 || s.at(i) == 93))
    ++i;
  ++i;
  std::string key;
  std::string ikey;
  if (i == s.size())
    key = s;
  else
    key = s.substr(0, (i - 1));
  ikey = s.substr(
      0, (i)); // Some Keywords ( eg. FILESTART_HEX ) have to be case sensitive!
  for (i = 0; i < key.size(); ++i)
    key.at(i) = std::tolower(key.at(i));

  // First parse over key, since more keywords are caseinsensitive

  for (i = 0; i < NumberOfKeywords; ++i)
  {
    if (key.compare(ListOfKeywords[i].keyword))
      continue;
    return i;
  }

  // Second parse with casesensitive keys

  for (i = 0; i < NumberOfKeywords; ++i)
  {
    if (ikey.compare(ListOfKeywords[i].keyword))
      continue;
    return i;
  }

  return 0;
}

void
InputFile::getKeyLine(std::string &key, std::string &error, unsigned int &line)
{
  unsigned int i;
  for (i = line; i < unformated.size(); ++i)
  {
    if (unformated.at(i).find(key) != (std::string::npos))
    {
      line = i;
      error +=
          ("Line " + std::to_string(line) + ": " + unformated.at(line) + "\n");
    }
  }
}

void
InputFile::echo(std::fstream &o)
{
  unsigned int i;
  for (i = 0; i < unformated.size(); ++i)
    o << unformated.at(i) << std::endl;
}

int
InputFile::comms(std::fstream &o)
{
  unsigned int i, j;
  for (i = 0, j = 0; i < commentary.size(); ++i)
  {
    if (commentary.at(i).size() > 2 && commentary.at(i).at(1) == 43)
    {
      o << commentary.at(i).substr(0, commentary.at(i).size() - 1) << std::endl;
      ++j;
    }
  }
  return j;
}
