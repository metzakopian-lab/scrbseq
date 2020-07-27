#include <iostream>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <vector>
#include <unordered_map>
#include <utility>
#include <string>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

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
    std::string read_name(256), bases(256), qualities(256), extra(256);
    while (in.good())
    {
      //and str.find("_") != std::string::npos
      std::getline(in, read_name);
      std::getline(in, bases);
      std::getline(in, qualities);
      std::getline(in, extra);
      
      if ( not (lines % 4 == 0 and !str.empty() and read_name[0] == '@'))
      {
        std::cerr << "Error reading in FastQ file" << std::endl;
        return 1;
      }

      
      std::vector<std::string> fields;
      boost::split(fields, read_name, boost::is_any_of("_"));
      if(fields.size() < 2)
      {
        std::cerr << str << std::endl;
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
      
      if (bases.size() > 8)
      {
        *stream << read_name << std::endl;
        *stream << bases << std::endl;
        *stream << qualities << std::endl;
        *stream << extra << std::endl;
      }
      
      // Report progress
      if (++records % 50000 == 0)
      {
        std::cerr << "Read " << records << " records so far" << "\r";
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
  
  std::string fq_file = argv[1];
  std::string manifest = argv[2];
  std::string output_prefix = argv[3];
  if (openBarcodeStreams(manifest, output_prefix)) {
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


