
#include <iostream>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <vector>
#include <unordered_map>
#include <utility>

typedef std::ofstream OutStream;

typedef std::unordered_map<std::string, OutStream*> Streams;

long missed = 0;

Streams barcodes;

int openBarcodeStreams(const std::string& filename)
{
  std::cerr<<"Reading Barcodes from " << filename << std::endl;
  std::ifstream bcFile(filename.c_str());

  for(std::string str; std::getline(bcFile,str);)
  {
    if(str.find(",") == std::string::npos)
    {
      continue;
    }
    std::vector<std::string> fields;
    boost::split(fields, str, boost::is_any_of(","));
    std::string barcode = fields[1];
    std::string sample_name = fields[0];
    sample_name += "_out.sam";
    std::cerr << "File " << sample_name << " for " << barcode << std::endl;
    std::ofstream *fileOut= new std::ofstream(sample_name.c_str(), 
                                               std::ios_base::out | std::ios_base::app);
    barcodes[barcode] = fileOut;
  }
  bcFile.close();
  return 0;
  
}


int main(int argc, char* argv[])
{

  std::ifstream inStream(argv[1], std::ios_base::binary);
  if(openBarcodeStreams(argv[2]))
  {

    std::cerr << "Error in opening barcode file " << std::endl;
  }else
  {
    std::cerr << "Read sample sheet "<< std::endl;
  }
  
  //std::ofstream outStream("test.out.gz", std::ios_base::binary);
  


  unsigned int records = 0;
  std::string decompressedStringBuf;

  std::ifstream ifile(argv[1]);
  long lines = 0;
  std::cerr << std::endl;
  for (std::string str; std::getline(ifile,str);++lines)
  {

    std::vector<std::string> samrecord;
    boost::split(samrecord, str, boost::is_any_of("\t"));
    
    std::vector<std::string> readname;
    boost::split(readname,samrecord[0],boost::is_any_of("_"));
    
    auto it = barcodes.find(readname[1]);
    if(it != barcodes.end())
    {
      *(it->second) << str << std::endl;
    }
    else
    {
      missed++;
    }
    
    if (lines % 100000 == 0)
    {
      std::cerr << "Records:" << lines <<" misses:" <<missed << "\r";
    }
      
  }
  std::cerr << std::endl;
  
  for(auto& b : barcodes)
  {
    auto stream = b.second;
    stream->close();
  }
  
  if (argc == 3)
  {
    

  }

  return 0;

}


