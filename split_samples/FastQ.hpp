#pragma once
#include <string>
#include <iostream>


class FastQ {
 public:
  std::string read_name;
  std::string bases;
  std::string qualities;
  std::string extra;
  int reserve()
  {
    read_name.reserve(256);
    bases.reserve(256);
    qualities.reserve(256);
    extra.reserve(256);
    return 0;
  }

  int parseRead(std::istream& fq_stream)
  {
    
    
    std::getline(fq_stream, read_name);
    std::getline(fq_stream, bases);
    std::getline(fq_stream, qualities);
    std::getline(fq_stream, extra);
    // TODO check if valid reads can start with another character.
    if ( (not read_name.empty()) and read_name[0] != '@')
    {
      std::cerr << "Error reading in FastQ file" << std::endl;
      return 1;
    }
    return 0;
  }
  int dumpRead(std::ostream& fq_stream){
    
    fq_stream << read_name << std::endl;
    fq_stream << bases << std::endl;
    fq_stream << qualities << std::endl;
    fq_stream << extra << std::endl;
    return 0;
  }
};
