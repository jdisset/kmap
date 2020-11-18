#pragma once
#include <rocksdb/db.h>

#include <cxxopts.hpp>
#include <indicators.hpp>
#include <iostream>

#include "rocksdb.hpp"
#include "utils.hpp"

// outputs in sparse matrix format:
// k-mer, sparse matrix indices, sparse matrix values
// the sparse matrix contains the count of every next letter for every dataset
inline void dbToSparseMatrix(RocksDB& rdb, std::string outputFile) {
	const int CHUNKSIZE = 500;

	std::FILE* file = std::fopen(outputFile.c_str(), "w");
	if (!file) throw std::runtime_error("Error opening output file");
	std::string estimatedN = "1000";
	rdb.db->GetProperty("rocksdb.estimate-num-keys", &estimatedN);
	size_t N = std::stoi(estimatedN);
	ProgressBar mainbar{option::BarWidth{50},
	                    option::Start{"["},
	                    option::Fill{"■"},
	                    option::Lead{"■"},
	                    option::Remainder{"-"},
	                    option::End{"]"},
	                    option::PostfixText{"Writing to " + outputFile},
	                    option::ShowElapsedTime{true},
	                    option::ForegroundColor{Color::yellow},
	                    option::MaxProgress{N + CHUNKSIZE},
	                    option::FontStyles{std::vector<FontStyle>{FontStyle::bold}}};

	// write header line
	std::fprintf(file, "kmer, count_mat_indices, count_mat_values\n");

	int c = 0;
	std::string buff;

	rdb.forEach([&](const auto* it) {
		datacount_t d = deserialize(it->value().ToString());
		std::string kmer = it->key().ToString();
		// auto leftpad = raw_seq_t(k - decoded.size(), '[');
		// kmer.insert(kmer.begin(), leftpad.begin(), leftpad.end());
		buff += kmer + "; [";

		std::vector<size_t> allcounts;
		for (const auto& [dataset, counts] : d) {
			assert(counts.size() > 0);
			for (size_t j = 0; j < counts.size() - 1; ++j) {
				if (counts[j] > 0) {  // if this character appears
					buff += "[" + std::to_string(dataset) + "," + std::to_string(j) + "],";
					allcounts.push_back(counts[j]);
				}
			}
		}
		buff.pop_back();  // removes trailing comma
		buff += "]; [";
		for (const auto& c : allcounts) buff += std::to_string(c) + ",";
		buff.pop_back();  // removes trailing comma
		buff += "]\n";
		if (++c > CHUNKSIZE) {
			std::fwrite(buff.c_str(), 1, buff.size(), file);
			buff.clear();
			mainbar.tick(c);
			c = 0;
		}
	});

	std::fwrite(buff.c_str(), 1, buff.size(), file);
	mainbar.tick(c);
	std::fclose(file);
	{  // final progress update
		mainbar.set_option(option::ForegroundColor{Color::green});
		mainbar.set_option(option::PrefixText{"✔ KMap written to " + outputFile});
		mainbar.set_option(option::BarWidth{0});
		mainbar.set_option(option::Fill{""});
		mainbar.set_option(option::Lead{""});
		mainbar.set_option(option::Start{""});
		mainbar.set_option(option::End{""});
		mainbar.set_option(option::ShowPercentage{false});
		mainbar.set_option(option::PostfixText{"                              "});
		mainbar.mark_as_completed();
	}
}
