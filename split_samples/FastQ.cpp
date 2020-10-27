
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
  // TODO check if valid reads can start with another character.
  if ((not read_name.empty()) and read_name[0] != '@') {
    std::cerr << "Error reading in FastQ file" << std::endl;
    return 1;
  }
  // if (bases.size() < fixed_size) {
  //   std::cerr << "Error reading in FastQ read" << std::endl;
  //   std::cerr << "Read name: " << read_name << std::endl;
  //   return 1;
  // }
  // Get all the before the space
  read_name = read_name.substr(0, readline.find(" "));
  // Get the SPACE and everything after
  read_rest = read_name.substr(readline.find(" "));
  
  return 0;
}

int FastQRecord::dumpRead(std::ostream& fq_stream) const
{ 
  fq_stream << read_name << read_rest << std::endl
            << bases << std::endl
            << qualities << std::endl
            << extra << std::endl;
  return 0;
}
