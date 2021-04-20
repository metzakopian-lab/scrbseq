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
      if ( R1 || R2) {
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
        //std::cerr << "Barcode " << barcode.getReadBases() << " "<<sample_barcode << " " << umi << std::endl;
        
        std::cerr << "Scanned " << records << " reads so far" << std::endl;
        flush_buffers();
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

  
  if (argc != 7)
  {
    std::cerr << "Usage:" << argv[0] << "<mRNA fastq file> <barcode fastq file> <manifest.csv (sample_name, barcode)> <output-prefix> <barcode-size> <umi-size>" << std::endl;
    return 1;
  }
  std::cout << "Version 2" << std::endl;
  std::string mRNA_file = argv[1];
  std::string barcode_file = argv[2];
  
  std::string manifest = argv[3];
  std::string output_prefix = argv[4];
  
  if (initializeSamples(manifest, output_prefix, barcodes, reads)) {
    std::cerr << "Error in loading barcode file " << std::endl;
    return 1;
  }

  barcode_size = std::stoi(argv[5]);
  umi_size = std::stoi(argv[6]);
  // Start the UMI straight after the Barcode
  umi_offset = barcode_size;
  
  
  ///std::ios::sync_with_stdio(false);
  if (parseFastQFile(mRNA_file, barcode_file)) {
    return 1;
  }

  exportReadsPerSample(output_prefix, reads);
  //std::ofstream outStream("test.out.gz", std::ios_base::binary);
  return 0;
}


