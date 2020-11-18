#include <rocksdb/db.h>

#include <cxxopts.hpp>
#include <indicators.hpp>
#include <iostream>

#include "export.hpp"
#include "rocksdb.hpp"

int main(int argc, char** argv) {
	std::string inputPath;
	std::string outputFile;

	enum class OutFormat { sparsematrix };
	OutFormat outformat = OutFormat::sparsematrix;

	cxxopts::Options options("kconv",
	                         "Convert a raw rocksdb output to various output format");
	options.add_options()("i,input", "input db path",
	                      cxxopts::value<std::string>(inputPath))(
	    "o,output", "output file path",
	    cxxopts::value<std::string>(outputFile)->default_value("kmap_out.csv"))(
	    "f,output-format", "output format (multimap, sparsematrix). Default = sparsematrix",
	    cxxopts::value<std::string>()->default_value("sparsematrix"))("h,help",
	                                                                  "Print usage");

	auto opt = options.parse(argc, argv);

	if (opt.count("help")) {
		std::cout << options.help() << std::endl;
		exit(0);
	}
	if (opt["output-format"].as<std::string>() == "sparsematrix") {
		outformat = OutFormat::sparsematrix;
	} else {
		std::cerr << "Invalid output format \"" << opt["output-format"].as<std::string>()
		          << "\"; available: sparsematrix" << std::endl;
		exit(1);
	}

	RocksDB rdb{};
	rdb.opts.error_if_exists = false;
	rdb.opts.create_if_missing = false;

	rdb.open(inputPath);
	dbToSparseMatrix(rdb, outputFile);
}
