
#include "FastQ.hpp"

FastQRecord::FastQRecord()
{
  read_name.reserve(256);
  bases.reserve(256);
  qualities.reserve(256);
  extra.reserve(256);
}

void FastQRecord::appendMeta(const std::string &barcode, const std::string &umi) {
  read_name += "_";
  read_name += barcode;
  read_name += "_";
  read_name += umi;
}


int FastQRecord::parseRead(std::istream& fq_stream) {

  std::string readline;

  std::getline(fq_stream, readline);
  std::getline(fq_stream, bases);
  std::getline(fq_stream, qualities);
  std::getline(fq_stream, extra);
  if (fq_stream.eof())
  {  
    return EOF_ERROR;
  }
  
  // TODO check if valid reads can start with another character.

  if ((not read_name.empty()) and read_name[0] != '@') {
    std::cerr << "Error reading in FastQ file" << std::endl;
    return 1;
  }
  
  // Get all the before the space
  const auto split = readline.find(" ");
  read_name = readline.substr(0, split);
  // Get the SPACE and everything after
  read_rest = readline.substr(split);
  
  return 0;
}

int FastQRecord::dumpRead(std::ostream& fq_stream) const
{
  
  fq_stream << read_name << read_rest << "\n"
            << bases << "\n"
            << qualities << "\n"
            << extra << "\n";
  return 0;
}
