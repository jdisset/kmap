#include <rocksdb/db.h>

#include <cxxopts.hpp>
#include <indicators.hpp>
#include <iostream>

#include "rocksdb.hpp"
#include "utils.hpp"

int main(int argc, char** argv) {
	RocksDB rdb{};
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
	    cxxopts::value<std::string>()->default_value("sparsematrix"));

	auto opt = options.parse(argc, argv);

	if (opt["output-format"].as<std::string>() == "sparsematrix") {
		outformat = OutFormat::sparsematrix;
	} else {
		std::cerr << "Invalid output format \"" << opt["output-format"].as<std::string>()
		          << "\"; available: sparsematrix" << std::endl;
		exit(1);
	}

	rocksdb::Options opts;
	opts.IncreaseParallelism();
	opts.OptimizeLevelStyleCompaction();

	/*RocksDB rdb(opts);*/
	// rdb.open(inputPath);

	// const int CHUNKSIZE = 500;

	// std::FILE* file = std::fopen(outputFile.c_str(), "w");
	// if (!file) throw std::runtime_error("Error opening output file");

	// size_t estimatedN = 1000;
	////rdb.db->GetProperty("rocksdb.estimate-num-keys", &estimatedN);
	// ProgressBar mainbar{option::BarWidth{50},
	// option::Start{"["},
	// option::Fill{"■"},
	// option::Lead{"■"},
	// option::Remainder{"-"},
	// option::End{"]"},
	// option::PostfixText{"Writing to " + outputFile},
	// option::ShowElapsedTime{true},
	// option::ForegroundColor{Color::yellow},
	// option::MaxProgress{estimatedN + CHUNKSIZE},
	// option::FontStyles{std::vector<FontStyle>{FontStyle::bold}}};

	// rdb.forEach([](const auto* it) {
	// buff += it->key().ToString() + "; " << it->value().ToString() << std::endl;
	/*});*/
}
