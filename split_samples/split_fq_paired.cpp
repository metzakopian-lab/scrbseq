#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include <string>
#include "FastQ.hpp"
#include "io.hpp"

#include <boost/iostreams/filter/gzip.hpp>
#include <boost/program_options.hpp>

FQStreams barcodes;
SampleCounts reads;


std::hash<std::string> hasher;

namespace po = boost::program_options;


// const std::string pattern ="CCCCCCXXXXXXXX";
unsigned barcode_offset = 0;
unsigned barcode_size = 6;
unsigned umi_offset = 6;
unsigned umi_size = 8;


inline bool areMates(FastQRecord& read1, FastQRecord& read2)
{
  return hasher(read1.getReadName()) == hasher(read2.getReadName());
}


// Flush demultiplexed reads
int flush_buffers()
{
  for (auto& b : barcodes)
  {
    auto stream = b.second;
    stream->flush();
    boost::iostreams::close(*stream);
  }

  return 0;
}


int parseFastQFile(const std::string& mRNA_fq_file, const std::string& barcode_fq_file)
{  
  unsigned int records = 0;
  try
  {
    std::string decompressedStringBuf;

    
    auto& mRNA_stream = *openFastQ(mRNA_fq_file);
    auto& barcode_stream = *openFastQ(barcode_fq_file);
   
    FastQRecord mRNA, barcode;
    
    while (mRNA_stream.good() && barcode_stream.good())
    {
      //and str.find("_") != std::string::npos
      auto R1 = mRNA.parseRead(mRNA_stream);
      auto R2 = barcode.parseRead(barcode_stream);
      if( R1 || R2) {
        continue;
      }
      
      if (!areMates(mRNA, barcode)) {
        std::cerr << "Read Name do not match" << std::endl <<
            "mRNA_read" << mRNA.getReadName() << std::endl <<
            "barcode_read" << barcode.getReadName() <<std::endl;
        return 1;
      }
      // TODO check these guys work
      std::string sample_barcode = barcode.getReadBases().substr(barcode_offset, barcode_size);
      std::string umi = barcode.getReadBases().substr(umi_offset, umi_size);
      
      mRNA.appendMeta(sample_barcode, umi);
      // Increment sample id
      {
        auto read_it = reads.find(sample_barcode);
        auto& sample_reads = read_it == reads.end() ? reads.at(undetermined).second : read_it->second.second;
        ++sample_reads;
      }
      
      // Forward record to the relevant file
      auto it = barcodes.find(sample_barcode);
      auto stream = it == barcodes.end() ? barcodes.at(undetermined) : it->second;
      mRNA.dumpRead(*stream);
      
      // Report progress and flush the buffers
      if (++records % 500000 == 0)
      { 
        std::cerr << "Scanned " << records << " reads so far" << std::endl;
        
        // Flush buffers to avoid memory problems
        for (auto& b : barcodes)
        {
          auto stream = b.second;
          stream->flush();
        }
      }
    }
    std::cerr << std::endl;
  }
  catch (const boost::iostreams::gzip_error& e )
  {
    std::cerr << e.what() << std::endl;
    return 1;
  }
  
  std::cerr << "Finished scanning! " << records << " reads in total!" << std::endl;
  flush_buffers();
  std::cerr << "Finished Final Flush" << std::endl;
  
  return 0;
}


int main(int argc, char* argv[])
{
  std::string mRNA_file, barcode_file, manifest, output_prefix;
  

  po::variables_map vm;
    po::options_description desc{"Options"};
    desc.add_options()
        ("help,h", "prints this help message")
        ("manifest,m", po::value<std::string>(&manifest)->required(), "csv containing the sample and barcodes")
        ("rna-file,r", po::value<std::string>(&mRNA_file)->required(), "Gzipped fastq file containing RNA content reads")
        ("barcode-file,b", po::value<std::string>(&barcode_file)->required(), "Gzipped fastq file containing barcode data immediately followed by UMI data")
        ("barcode-length", po::value(&barcode_size)->default_value(6), "Barcode size in nucleotides file reads")
        ("umi-length", po::value(&umi_size)->default_value(8), "UMI length in nucleotides file reads");
    
  try
  {
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
  }
  catch(const po::error &ex){ 
    std::cerr << ex.what() << std::endl;
  }

  umi_offset = barcode_size;
  
  if (initializeSamples(manifest, output_prefix, barcodes, reads)) {
    std::cerr << "Error in loading barcode file " << std::endl;
    return 1;
  }

  
  if (parseFastQFile(mRNA_file, barcode_file)) {
    return 1;
  }

  exportReadsPerSample(output_prefix, reads);
  //std::ofstream outStream("test.out.gz", std::ios_base::binary);
  return 0;
}


