#pragma once
#include <string>
#include <iostream>


class FastQRecord {
 private:
  std::string read_rest;
  std::string read_name;
  std::string bases;
  std::string qualities;
  std::string extra;
 public:
  FastQRecord();
  void appendMeta(const std::string &barcode, const std::string &umi);
  int parseRead(std::istream& fq_stream);
  int dumpRead(std::ostream& fq_stream) const;
  inline const std::string getReadName() const
  {
    const auto& rn = this->read_name;
    return rn.substr(0,rn.find(" "));
  }
  inline const std::string &getReadLine() const
  {
    return read_name;
  }

  inline const std::string& getReadBases() const 
  {
    return bases;
  }
};




