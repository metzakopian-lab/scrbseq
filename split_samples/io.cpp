#include <iostream>
#include <vector>

#include <boost/algorithm/string.hpp>

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include "io.hpp"



FQOut* addStream(const std::string& filename)
{
  // http://boost.2283326.n4.nabble.com/Writing-large-binary-files-with-boost-gzip-td3434404.html
  FQOut *out = new FQOut();
  out->push(boost::iostreams::gzip_compressor(boost::iostreams::zlib::default_compression, BUFFER_SIZE));
  out->push(boost::iostreams::file_sink(filename));
  return out;
}


int initializeSamples(const std::string& filename, const std::string& prefix, FQStreams& samples, SampleCounts& counts)
{
  std::string unmapped_reads_fq = "UNDETERMINED.fq.gz";
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
    std::string output_file = prefix + sample_name + ".fq.gz";
    samples[barcode] = addStream(output_file);
    counts[barcode] = std::make_pair(sample_name, 0);
    
    std::cerr << "File " << sample_name << " for " << barcode << std::endl;
  }
  bcFile.close();

  samples[undetermined] = addStream(undetermined);
  counts[undetermined] = std::make_pair(undetermined, 0);
  std::cerr << "Unmapped reads are exported to " << unmapped_reads_fq << std::endl;
  return 0;
}



FQIn* openFastQ(const std::string& fq_file, bool gzip) {
  std::ifstream* ifile = new std::ifstream(fq_file, std::ios_base::in | std::ios_base::binary);
  if (!ifile->good())
  {
    std::cerr << "File " << fq_file << " could not be opened" << std::endl;
  }
  //ifile.close();
     
  auto in = new boost::iostreams::filtering_streambuf<boost::iostreams::input>();
  if (gzip) {
    in->push(boost::iostreams::gzip_decompressor());
  }
  in->push(*ifile);
  
  return new std::istream(in);
}



void exportReadsPerSample(const std::string & exp_name, const SampleCounts& counts)
{
  std::string filename = exp_name + "-report.csv";
  std::ofstream csv(filename.c_str());
  csv << "SampleName,Reads" <<std::endl;
  for (auto& sample : counts)
  {
    csv << sample.second.first << "," << sample.second.second << std::endl;
  }
  csv.close();
}
