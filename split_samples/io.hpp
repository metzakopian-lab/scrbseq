#pragma once

/// 64MB buffer size
#define BUFFER_SIZE 2 << 26


#include <unordered_map>
#include <boost/iostreams/filtering_stream.hpp>


typedef boost::iostreams::filtering_ostream FQOut;
typedef std::istream FQIn;


typedef std::unordered_map<std::string, FQOut*> FQStreams;
typedef std::unordered_map<std::string, std::pair<std::string, unsigned int>> SampleCounts;

const std::string undetermined = "UNDETERMINED";



FQIn* openFastQ(
    const std::string& fq_file,bool gzip = true);


FQOut* addStream(const std::string& filename);


int initializeSamples(
    const std::string& filename,
    const std::string& prefix,
    FQStreams&,
    SampleCounts&);


void exportReadsPerSample(const std::string & exp_name, const SampleCounts& counts);
