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
/// 64MB buffer size
#define BUFFER_SIZE 2 << 26 


typedef boost::iostreams::filtering_ostream OutStream;
typedef std::unordered_map<std::string, OutStream*> Streams;
typedef std::unordered_map<std::string, std::pair<std::string, unsigned int>> Counts;


Streams barcodes;
Counts reads;
const std::string undetermined = "UNDETERMINED";
std::string unmapped_reads_fq;

// const std::string pattern ="CCCCCCXXXXXXXX";
unsigned barcode_offset = 8;
unsigned umi_offset = 14;


OutStream* addStream(const std::string& filename)
{
  // http://boost.2283326.n4.nabble.com/Writing-large-binary-files-with-boost-gzip-td3434404.html
  OutStream *out = new OutStream();
  out->push(boost::iostreams::gzip_compressor(boost::iostreams::zlib::default_compression, BUFFER_SIZE));
  out->push(boost::iostreams::file_sink(filename));
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
    std::string output_file = prefix + "-" + sample_name + ".fq.gz";
    barcodes[barcode] = addStream(output_file);
    reads[barcode] = std::make_pair(sample_name, 0);
    
    std::cerr << "File " << sample_name << " for " << barcode << std::endl;
  }
  bcFile.close();

  barcodes[undetermined] = addStream(undetermined);
  reads[undetermined] = std::make_pair(undetermined, 0);
  
  std::cerr << "Unmapped reads are exported to " << unmapped_reads_fq << std::endl;
  return 0;
}

boost::iostreams::filtering_istream* openFastQ(const std::string& fq_file, bool gzip = true) {
  std::ifstream ifile(fq_file, std::ios_base::in | std::ios_base::binary);

  if (!ifile.good())
  {
    std::cerr << "File " << fq_file << " could not be opened" << std::endl;
  }
  ifile.close();

  boost::iostreams::filtering_istream* in = new boost::iostreams::filtering_istream();
  if (gzip) {
    in->push(boost::iostreams::gzip_decompressor(
        boost::iostreams::zlib::default_window_bits, 4 * BUFFER_SIZE));
  }
  in->push(boost::iostreams::file_source(fq_file));
  return in;
}
// TODO
bool areMates(FastQ& read1, FastQ& read2)
{
  // hash strings to cmpr
  auto rname1 = read1.read_name.find(" ");
  //calculate.hash for read1 to rname position
  size_t hash1 = 0;
  auto rname2 = read2.read_name.find(" ");
  size_t hash2 = 0;
  return hash1 == hash2;
}

int parseFastQFile(const std::string& mRNA_fq_file, const std::string& barcode_fq_file)
{
  
  unsigned int records = 0;
  try
  {
    
    std::string decompressedStringBuf;

    
    auto& mRNA_stream = *openFastQ(mRNA_fq_file);
    auto& barcode_stream = *openFastQ(barcode_fq_file);
    
    long lines = 0;
    FastQ mRNA_read, barcode_read;

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
            "mRNA_read" << mRNA_read.read_name << std::endl <<
            "barcode_read" << barcode_read.read_name <<std::endl;
        return 1;
      }
      
      std::string barcode = barcode_read.read_name.substr(barcode_offset);
      std::string umi = barcode_read.read_name.substr(barcode_offset, umi_size);
      mRNA_read.append_meta(barcode, umi);
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
    std::cerr << "Usage:" << argv[0] << " <mRNA fastq file> <barcode fastq file> <manifest.csv (sample_name, barcode)> <output-prefix>" << std::endl;
    return 1;
  }
  std::cout << "Version 2" << std::endl;
  std::string mRNA_file = argv[1];
  std::string barcode_file = argv[2];
  
  std::string manifest = argv[2];
  std::string output_prefix = argv[3];
  if (openBarcodeStreams(manifest, output_prefix)) {
    std::cerr << "Error in loading barcode file " << std::endl;
    return 1;
  }

  if (parseFastQFile(mRNA_file, barcode_file)) {
    return 1;
  }

  exportReadsPerSample(output_prefix);
  
  //std::ofstream outStream("test.out.gz", std::ios_base::binary);
  return 0;

}


