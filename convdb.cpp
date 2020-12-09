#include <cxxopts.hpp>
#include <indicators.hpp>
#include <iostream>

#include "export.hpp"
#include "kmapdb.hpp"

int main(int argc, char** argv) {
	std::string inputPath;
	std::string outputFile;

	enum class OutFormat { sparsematrix, condensed };
	OutFormat outformat = OutFormat::sparsematrix;

	cxxopts::Options options("kconv",
	                         "Convert a raw rocksdb output to various output format");
	options.add_options()("i,input", "input db path",
	                      cxxopts::value<std::string>(inputPath))(
	    "o,output", "output file path",
	    cxxopts::value<std::string>(outputFile)->default_value("kmap_out.csv"))(
	    "f,output-format",
	    "output format (sparsematrix, condensed). Default = sparsematrix",
	    cxxopts::value<std::string>()->default_value("sparsematrix"))("h,help",
	                                                                  "Print usage");

	auto opt = options.parse(argc, argv);

	if (opt.count("help")) {
		std::cout << options.help() << std::endl;
		exit(0);
	}

	if (opt["output-format"].as<std::string>() == "sparsematrix") {
		outformat = OutFormat::sparsematrix;
	} else if (opt["output-format"].as<std::string>() == "condensed") {
		outformat = OutFormat::condensed;
	} else {
		std::cerr << "Invalid output format \"" << opt["output-format"].as<std::string>()
		          << "\"; available: sparsematrix" << std::endl;
		exit(1);
	}

	KmapDB db{};
	db.openRO(inputPath);

	switch (outformat) {
		case OutFormat::condensed:
			// dbToCondensed(db, outputFile);
			break;
		default:
		case OutFormat::sparsematrix:
			dbToSparseMatrix(db, outputFile);
			break;
	}
}
