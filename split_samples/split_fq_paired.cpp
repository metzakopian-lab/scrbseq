#include <iostream>
#include <fstream>
#include <vector>

#include <utility>
#include <string>


#include "FastQ.hpp"
#include "io.hpp"


#include <boost/iostreams/filter/gzip.hpp>

FQStreams barcodes;
SampleCounts reads;


std::hash<std::string> hasher;
// const std::string pattern ="CCCCCCXXXXXXXX";
unsigned barcode_offset = 8;
unsigned umi_offset = 14;



inline bool areMates(FastQRecord& read1, FastQRecord& read2)
{
  return hasher(read1.getReadName()) == hasher(read2.getReadName());
}

int parseFastQFile(const std::string& mRNA_fq_file, const std::string& barcode_fq_file)
{  
  unsigned int records = 0;
  try
  {
    std::string decompressedStringBuf;

    
    auto& mRNA_stream = *openFastQ(mRNA_fq_file);
    auto& barcode_stream = *openFastQ(barcode_fq_file);
    
    FastQRecord mRNA_read, barcode_read;

    unsigned umi_size = umi_offset - barcode_offset;
    
    while (mRNA_stream.good() && barcode_stream.good())
    {
      //and str.find("_") != std::string::npos
      if (mRNA_read.parseRead(mRNA_stream) ||
          barcode_read.parseRead(barcode_stream)) {
        continue;
      }
      if (!areMates(mRNA_read, barcode_read)) {
        std::cerr << "Read Name do not match" << std::endl <<
            "mRNA_read" << mRNA_read.getReadName() << std::endl <<
            "barcode_read" << barcode_read.getReadName() <<std::endl;
        return 1;
      }
      // TODO check these guys work
      std::string barcode = barcode_read.getReadBases().substr(barcode_offset);
      std::string umi = barcode_read.getReadBases().substr(barcode_offset, umi_size);
      std::cerr << "Read barcode" << barcode << " with UMI " << umi << std::endl;
      mRNA_read.appendMeta(barcode, umi);
      // Increment sample id
      {
        auto read_it = reads.find(barcode);
        auto& sample_reads = read_it == reads.end() ? reads.at(undetermined).second : read_it->second.second;
        ++sample_reads;
      }
      
      // Forward record to the relevant file
      auto it = barcodes.find(barcode);
      auto stream = it == barcodes.end() ? barcodes.at(undetermined) : it->second;
      mRNA_read.dumpRead(*stream);
      
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
  for (auto& b : barcodes)
  {
    auto stream = b.second;
    stream->flush();
    boost::iostreams::close(*stream);
  }
  return 0;
}



int main(int argc, char* argv[])
{
  if (argc != 4)
  {
    std::cerr << "Usage:" << argv[0] << " <mRNA fastq file> <barcode fastq file> <manifest.csv (sample_name, barcode)> <output-prefix>" << std::endl;
    return 1;
  }
  std::cout << "Version 2" << std::endl;
  std::string mRNA_file = argv[1];
  std::string barcode_file = argv[2];
  
  std::string manifest = argv[2];
  std::string output_prefix = argv[3];
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


