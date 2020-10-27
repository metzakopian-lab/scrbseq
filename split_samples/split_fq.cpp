#include <iostream>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <vector>
#include <unordered_map>
#include <utility>
#include <string>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>

#include "FastQ.hpp"
#include "io.hpp"

FQStreams barcodes;
SampleCounts reads;

int parseFastQFile(const std::string& fq_file)
{
  std::ifstream ifile(fq_file, std::ios_base::in | std::ios_base::binary);

  if (!ifile.good())
  {
    std::cerr << "File " << fq_file << " could not be opened" << std::endl;
  }
  ifile.close();

  unsigned int records = 0;
  
  try
  {
    
    std::string decompressedStringBuf;

    
    boost::iostreams::filtering_istream in;
    in.push(boost::iostreams::gzip_decompressor(boost::iostreams::zlib::default_window_bits, 4*BUFFER_SIZE));
    in.push(boost::iostreams::file_source(fq_file));
    
    FastQRecord read;
    
    while (in.good())
    {
      //and str.find("_") != std::string::npos
      read.parseRead(in);
      std::vector<std::string> fields;
      boost::split(fields, read.getReadName(), boost::is_any_of("_"));
      if(fields.size() < 2)
      {
        std::cerr << read.getReadName() << std::endl;
        std::cerr <<  "Line not properly formatted, did UMI tools run? " << std::endl;
        return 1;
      }
      
      std::string sample_id = fields[1];
      
      // Increment sample id
      {
        auto read_it = reads.find(sample_id);
        auto& sample_reads = read_it == reads.end() ? reads.at(undetermined).second : read_it->second.second;
        ++sample_reads;
      }

      // Forward record to the relevant file
      auto it = barcodes.find(sample_id);
      auto stream = it == barcodes.end() ? barcodes.at(undetermined) : it->second;
      read.dumpRead(*stream);
      
      // Report progress
      if (++records % 500000 == 0)
      {
        std::cerr << "Read " << records << " records so far" << std::endl;
      }
      
    }
    std::cerr << std::endl;
  }
  catch (const boost::iostreams::gzip_error& e )
  {
    std::cerr << e.what() << std::endl;
  }
  
  std::cerr << "Written " << records << " reads in total!" << std::endl;
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
  for (auto& sample : reads)
  {
    csv << sample.second.first << "," << sample.second.second << std::endl;
  }
  csv.close();
}

int main(int argc, char* argv[])
{
  if (argc != 4)
  {
    std::cerr << "Usage:" << argv[0] << " <UMI_TOOLS input fastq file> <manifest.csv (sample_name, barcode)> <output-prefix>" << std::endl;
    return 1;
  }
  std::cout << "Version 2" << std::endl;
  std::string fq_file = argv[1];
  std::string manifest = argv[2];
  std::string output_prefix = argv[3];
  if (initializeSamples(manifest, output_prefix, barcodes, reads)) {
    std::cerr << "Error in loading barcode file " << std::endl;
    return 1;
  }

  if (parseFastQFile(fq_file)) {
    return 1;
  }

  exportReadsPerSample(output_prefix);
  
  //std::ofstream outStream("test.out.gz", std::ios_base::binary);
  return 0;

}


