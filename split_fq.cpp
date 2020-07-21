#include <iostream>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <vector>
#include <unordered_map>
#include <utility>
#include <string>



#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>


//#include <boost/iostreams/device/file.hpp>
//#include <boost/iostreams/stream.hpp>


typedef boost::iostreams::filtering_ostream OutStream;


typedef std::unordered_map<std::string, OutStream*> Streams;

typedef std::unordered_map<std::string, std::pair<std::string, unsigned int>> Counts;




Streams barcodes;
Counts reads;
const std::string undetermined = "UNDETERMINED";
std::string unmapped_reads_fq;




OutStream* addStream(const std::string& filename)
{
  std::ofstream *fileOut = new std::ofstream(
      filename.c_str(), 
      std::ios_base::binary | std::ios_base::out);
  
  if (!fileOut->good())
  {
    return NULL;
  }

  OutStream *out = new OutStream();
  out->push(boost::iostreams::gzip_compressor());
  out->push(*fileOut);
  return out;
}


int openBarcodeStreams(const std::string& filename, const std::string& prefix)
{
  std::cerr<<"Reading Barcodes from " << filename << std::endl;
  std::ifstream bcFile(filename.c_str());
  if (!bcFile.good())
  {
    std::cerr << "File " << filename << " could not be opened" << std::endl;
  }
  for(std::string str; std::getline(bcFile,str);)
  {
    
    std::vector<std::string> fields;
    boost::split(fields, str, boost::is_any_of(" ,"));
    std::string sample_name = fields[0];
    std::string barcode = fields[1];
    
    barcodes[barcode] = addStream( prefix + "-" + sample_name + ".fq.gz" );
    reads[barcode] = std::make_pair(sample_name, 0);
    
    std::cerr << "File " << sample_name << " for " << barcode << std::endl;
  }
  bcFile.close();

  barcodes[undetermined] = addStream(undetermined);
  reads[undetermined] = std::make_pair(undetermined, 0);
  
  std::cerr << "Unmapped reads are exported to " << unmapped_reads_fq << std::endl;
  return 0;
}

int parseFastQFile(const std::string& fq_file)
{
  std::ifstream ifile(fq_file, std::ios_base::in | std::ios_base::binary);

  if (!ifile.good())
  {
    std::cerr << "File " << fq_file << " could not be opened" << std::endl;
  }

  unsigned int records = 0;
  
  try
  {
    
    std::string decompressedStringBuf;

    
    boost::iostreams::filtering_istream in;
    in.push(boost::iostreams::gzip_decompressor());
    in.push(ifile);
    
    long lines = 0;
    for (std::string str; std::getline(in,str);++lines)
    {
      //and str.find("_") != std::string::npos
      if ( not (lines % 4 == 0 and !str.empty() and str[0] == '@'))
      {
        continue;
      }
      std::vector<std::string> fields;
      boost::split(fields, str, boost::is_any_of("_"));
      if(fields.size() < 2)
      {
        std::cerr << str << std::endl;
        std::cerr <<  "Line not properly formatted, did UMI tools run? " << std::endl;
        return 1;
      }



      std::string sample_id = fields[1];
      
      auto it = barcodes.find(sample_id);
      auto stream = it == barcodes.end() ? barcodes.at(undetermined) : it->second;

      
      *stream << str << std::endl;

      
      std::string tmp;
      for (int i = 0; i < 3 && std::getline(in,tmp); ++i, lines++)
      {
        *stream << tmp << std::endl;
      }
        
      if (++records % 30000 == 0)
      {
        std::cerr << "written " << records << " so far" << "\r";
      }
      
    }
    std::cerr << std::endl;
  }
  catch (const boost::iostreams::gzip_error& e )
  {
    std::cerr << e.what() << std::endl;
  }

  
  
  std::cerr << "written " << records << " so far" << std::endl;
  
  
  std::cerr << "Flushing buffers" << std::endl;
  for(auto& b : barcodes)
  {
    auto stream = b.second;
    stream->flush();
    boost::iostreams::close(*stream);
  }
  return 0;
  
}

void exportReadsPerSample(const std::string & exp_name)
{
  std::string filename = exp_name + "-report.csv";
  std::ofstream csv(filename.c_str());
  csv << "SampleName,Reads" <<std::endl;
  for (auto& sample : reads) {
    csv << sample.second.first << sample.second.first << std::endl;
  }
  csv.close();
}

int main(int argc, char* argv[])
{
  if (argc != 3)
  {
    std::cerr << "Usage:" << argv[0] << " <UMI_TOOLS input fastq file> <manifest.csv (sample_name, barcode)>" << std::endl;
    return 1;
  }
  std::string fq_file = argv[1];
  
  std::string experiment_name = fq_file;
  experiment_name.erase(
      experiment_name.begin() + experiment_name.find_last_of(".fq.gz"),
      experiment_name.end());

  if (openBarcodeStreams(argv[2], experiment_name)) {
    std::cerr << "Error in loading barcode file " << std::endl;
    return 1;
  }

  if (parseFastQFile(fq_file)) {
    return 1;
  }

  exportReadsPerSample(experiment_name);
  
  //std::ofstream outStream("test.out.gz", std::ios_base::binary);

  return 0;

}


