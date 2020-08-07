#pragma once
#include <string>
#include <iostream>


class FastQ {
 public:
  std::string stripped_readname;
  std::string stripped_rest;
  std::string read_name;
  std::string bases;
  std::string qualities;
  std::string extra;
  FastQ()
  {
    stripped_readname.reserve(256);
    read_name.reserve(256);
    bases.reserve(256);
    qualities.reserve(256);
    extra.reserve(256);
  }

  void append_meta(const std::string &barcode, const std::string &umi) {
    read_name = stripped_readname;
    read_name += "_";
    read_name += barcode;
    read_name += "_";
    read_name += umi;
    read_name += " ";
    read_name += stripped_rest;
  }

  void parseMeta() {
  stripped_readname = read_name.substr(read_name.find(" "));
  stripped_rest = read_name.substr(read_name.find(" "));
}

  int parseRead(std::istream& fq_stream, unsigned min_size = 14, unsigned fixed_size = 256)
  {

  std::getline(fq_stream, read_name);
  std::getline(fq_stream, bases);
  std::getline(fq_stream, qualities);
  std::getline(fq_stream, extra);
  // TODO check if valid reads can start with another character.
  if ((not read_name.empty()) and read_name[0] != '@') {
    std::cerr << "Error reading in FastQ file" << std::endl;
    return 1;
  }
  if (bases.size() < fixed_size) {
    std::cerr << "Error reading in FastQ read" << std::endl;
    std::cerr << "Read name: " << read_name << std::endl;
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
