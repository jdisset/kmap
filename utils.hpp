#pragma once
#include <charconv>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <indicators.hpp>
#include <json.hpp>
#include <memory>
#include <random>
#include <robinhood.hpp>
#include <string_view>
#include <unordered_map>
#include <xhash.hpp>

#include "alphabet.hpp"

using namespace indicators;  // for progress bars
using json = nlohmann::json;
namespace fs = std::filesystem;

// type definitions

using raw_seq_t = std::vector<char>;           // just a normal ascii sequence
using encoded_seq_t = std::vector<int8_t>;     // a sequence encoded into an alphabet
using dataset_t = std::vector<encoded_seq_t>;  // list of encoded sequences

using seqview_t = std::basic_string_view<int8_t>;  // view of encoded sequence

template <typename K, typename V> using umap_t = robin_hood::unordered_map<K, V>;

using kmap_t = umap_t<seqview_t, std::vector<uint64_t>>;      // seqview -> counts
using datacount_t = umap_t<uint64_t, std::vector<uint64_t>>;  // counts per dataset id
using multikmap_t = umap_t<seqview_t, datacount_t>;  // seqview -> counts per dataset

inline datacount_t deserialize(const std::string& s) {
	datacount_t d{};
	auto j = json::parse(s);
	for (const auto& e : j.items())
		d[atoll(e.key().c_str())] = e.value().get<std::vector<uint64_t>>();
	return d;
}

inline std::string serialize(const datacount_t& d) {
	json j{};
	for (const auto& [k, v] : d) j[k] = v;
	return j.dump();
}

inline void addCounts(std::vector<uint64_t>& a, const std::vector<uint64_t>& b) {
	assert(a.size() == b.size());
	for (size_t i = 0; i < a.size(); ++i) a[i] += b[i];
}

// encode / decode. Important for 2 reasons:
// - normalizing uppercase & lowercases
// - encoding of a letter = index in the count array
static const constexpr Alpha alphaMap{};
inline encoded_seq_t encodeSequence(const raw_seq_t& s, const Alphabet& alpha) {
	encoded_seq_t r(s.size() + 2);  // don't forget prepend and append symbols
	r[0] = alphaMap.encode('[', alpha);
	for (size_t i = 1; i < r.size() - 1; ++i) r[i] = alphaMap.encode(s[i - 1], alpha);
	r[r.size() - 1] = alphaMap.encode(']', alpha);
	return r;
}

inline raw_seq_t decodeSequence(const encoded_seq_t& s, const Alphabet& alpha) {
	const auto decode = (alpha == Alphabet::dna) ?
	    [](const int8_t& c) -> char { return alphaMap.decode<Alphabet::dna>(c); } :
	    ((alpha == Alphabet::rna) ?
	     [](const int8_t& c) -> char { return alphaMap.decode<Alphabet::rna>(c); } :
	     [](const int8_t& c) -> char { return alphaMap.decode<Alphabet::protein>(c); });

	raw_seq_t r(s.size() - 2);
	for (size_t i = 0; i < s.size() - 2; ++i) r[i] = decode(s[i + 1]);
	return r;
}

inline raw_seq_t decodeSequenceView(const seqview_t& s, const Alphabet& alpha) {
	const auto decode = (alpha == Alphabet::dna) ?
	    [](const int8_t& c) -> char { return alphaMap.decode<Alphabet::dna>(c); } :
	    ((alpha == Alphabet::rna) ?
	     [](const int8_t& c) -> char { return alphaMap.decode<Alphabet::rna>(c); } :
	     [](const int8_t& c) -> char { return alphaMap.decode<Alphabet::protein>(c); });
	raw_seq_t r(s.size());
	for (size_t i = 0; i < s.size(); ++i) r[i] = decode(s[i]);
	return r;
}

template <typename B>
std::vector<encoded_seq_t> readFasta(const std::string& filepath, const Alphabet& alpha,
                                     B& pbars, size_t barid) {
	fs::path p{filepath};
	std::vector<encoded_seq_t> res;
	std::ifstream file(filepath.c_str());

	size_t charactersRead = 0;
	raw_seq_t currentSeq{};

	auto recordSeq = [&]() {
		if (currentSeq.size() > 0)
			res.push_back(encodeSequence(std::exchange(currentSeq, {}), alpha));
	};

	for (std::string line{}; std::getline(file, line);) {
		if (line[0] == '>')
			recordSeq();
		else {
			charactersRead += line.size();
			currentSeq.insert(currentSeq.end(), line.begin(), line.end());
		}
		if (charactersRead > 10000) pbars[barid].tick(std::exchange(charactersRead, 0));
	}
	recordSeq();
	pbars[barid].tick(charactersRead);

	return res;
}

// returns the number of sequence letters in a fasta file
inline size_t getFastaSize(const fs::path& filepath) {
	std::ifstream file(filepath.c_str());
	size_t nLetters = 0;
	for (std::string line{}; std::getline(file, line);)
		if (line[0] != '>') nLetters += line.size();
	return nLetters;
}

inline auto readPaths(const std::string& inputFile) {
	fs::path filePath{fs::absolute(fs::path(inputFile))};
	auto basePath = filePath.parent_path();
	std::ifstream file(filePath.c_str());

	std::vector<fs::path> paths;
	std::vector<size_t> sizes;

	for (std::string line{}; std::getline(file, line);) {
		fs::path datapath = fs::path(basePath) / line;
		paths.push_back(datapath);
		sizes.push_back(getFastaSize(datapath));
	}

	return std::make_pair(paths, sizes);
}

inline std::pair<std::vector<dataset_t>, size_t> readDatasets(
    const std::vector<fs::path>& paths, const std::vector<size_t>& sizes, size_t offset,
    size_t maxSize, Alphabet alpha, DynamicProgress<ProgressBar>& bars) {
	assert(paths.size() == sizes.size());

	size_t readSize = 0;
	size_t next = offset;
	while (next < sizes.size() && readSize < maxSize) {
		readSize += sizes[next];
		++next;
	}

	auto barid = bars.make_bar(
	    option::BarWidth{50}, option::Start{"["}, option::Fill{"■"}, option::Lead{"■"},
	    option::Remainder{"-"}, option::End{"]"}, option::ShowElapsedTime{true},
	    option::ForegroundColor{Color::yellow}, option::MaxProgress{readSize},
	    option::FontStyles{std::vector<FontStyle>{FontStyle::bold}});

	std::vector<dataset_t> datasets;
	for (size_t i = offset; i < next; ++i) {
		bars[barid].set_option(option::PostfixText{"Loading " + paths[i].string()});
		datasets.push_back(readFasta(paths[i], alpha, bars, barid));
	}
	bars[barid].mark_as_completed();

	return {datasets, next};
}

// generate N random dna sequences of size L
inline std::vector<encoded_seq_t> randomDNASequences(int L, int N) {
	const std::vector<char> letters{'A', 'C', 'G', 'T'};
	std::random_device rd;
	std::mt19937 gen(0);
	std::uniform_int_distribution<> d(0, letters.size() - 1);
	std::vector<encoded_seq_t> allseqs;
	for (int i = 0; i < N; ++i) {
		auto s = std::vector<char>(L);
		for (auto& n : s) n = letters[d(gen)];
		allseqs.push_back(encodeSequence(s, Alphabet::dna));
	}
	return allseqs;
}

inline void dumpMap(const kmap_t& m, int k, Alphabet alpha, std::string outputpath) {
	const int CHUNKSIZE = 10000;

	// std::ofstream file(outputpath);
	// if (file.fail()) throw std::runtime_error("Error opening output file");

	std::FILE* file = std::fopen(outputpath.c_str(), "w");
	if (!file) throw std::runtime_error("Error opening output file");

	ProgressSpinner spinner{option::PostfixText{"Writing to " + outputpath},
	                        option::ForegroundColor{Color::yellow},
	                        option::FontStyles{std::vector<FontStyle>{FontStyle::bold}}};

	auto alphabet = Alpha::getAlphabet(alpha);

	// Title line:
	// KMER, char0, char1, char2, ..., ']'
	// we ignore the last column since it's the begin character. It can't be in the counts
	std::fprintf(file, "%s", "kmer");
	for (size_t i = 0; i < alphabet.size() - 1; ++i) fprintf(file, ", %c", alphabet[i]);
	std::fprintf(file, "%s", "\n");

	int c = 0;
	std::string buff;
	for (const auto& kv : m) {
		auto decoded = decodeSequenceView(kv.first, alpha);
		auto leftpad = raw_seq_t(k - decoded.size(), '[');
		buff.insert(buff.end(), leftpad.begin(), leftpad.end());
		buff.insert(buff.end(), decoded.begin(), decoded.end());

		for (size_t i = 0; i < kv.second.size() - 1; ++i)
			buff += ", " + std::to_string(kv.second[i]);
		buff += "\n";

		if (++c % CHUNKSIZE == 0) {
			std::fwrite(buff.c_str(), 1, buff.size(), file);
			buff.clear();
			spinner.set_progress(static_cast<int>(((float)c / (float)m.size()) * 100));
		}
	}
	std::fwrite(buff.c_str(), 1, buff.size(), file);
	std::fclose(file);
	spinner.set_option(option::ForegroundColor{Color::green});
	spinner.set_option(option::PrefixText{"✔"});
	spinner.set_option(option::ShowSpinner{false});
	spinner.set_option(option::ShowPercentage{false});
	spinner.set_option(option::PostfixText{"All k-mers written!"});
	spinner.mark_as_completed();
}

// merge multiple multikmaps into one multikmap
inline void mergeMultiMaps(multikmap_t& res, std::vector<multikmap_t>& maps,
                           size_t startIndex = 0,
                           DynamicProgress<ProgressBar>* bars = nullptr) {
	size_t N = maps.size();

	size_t barid = 0;
	if (bars) {
		size_t entries = 0;
		for (size_t m = startIndex; m < N; ++m) entries += maps[m].size();
		barid = bars->make_bar(
		    option::BarWidth{50}, option::Start{"["}, option::Fill{"■"}, option::Lead{"■"},
		    option::Remainder{"-"}, option::End{"]"},
		    option::PostfixText{"Merging " + std::to_string(maps.size()) + " K-Maps"},
		    option::ShowElapsedTime{true}, option::ForegroundColor{Color::yellow},
		    option::MaxProgress{entries + 1},
		    option::FontStyles{std::vector<FontStyle>{FontStyle::bold}});
	}

	int ticks = 0;
	for (size_t m = startIndex; m < N; ++m) {
		for (auto& [kmer, map] : maps[m]) {
			if (!res.count(kmer))
				res[kmer] = decltype(map)(std::move(map));
			else {
				for (auto& [dataset, counts] : map) {
					if (res[kmer].count(dataset)) {
						assert(res[kmer].size() == counts.size());
						auto& ref = res[kmer][dataset];
						for (size_t i = 0; i < counts.size(); ++i) ref[i] += counts[i];
					} else {
						res[kmer][dataset] = counts;
					}
				}
			}
			if (++ticks > 1000) {
				if (bars) (*bars)[barid].tick(ticks);
				ticks = 0;
			}
		}
	}
	if (bars) {
		auto& b = (*bars)[barid];
		b.set_option(option::ForegroundColor{Color::green});
		b.set_option(option::PrefixText{"✔ Merged " + std::to_string(maps.size()) +
		                                " kmaps. Cleaning memory..."});
		b.set_option(option::BarWidth{0});
		b.set_option(option::Fill{""});
		b.set_option(option::Lead{""});
		b.set_option(option::Start{""});
		b.set_option(option::End{""});
		b.set_option(option::ShowPercentage{false});
		b.set_option(option::PostfixText{"                              "});
		(*bars)[barid];  // update
	}
	for (size_t m = startIndex; m < N; ++m) maps[m].clear();  // avoid wasting memory
	if (bars) {
		(*bars)[barid].mark_as_completed();
		(*bars)[barid];  // update
	}
}

inline void dumpMultiMap(const multikmap_t& m, const std::vector<std::string>& names,
                         int k, Alphabet alpha, std::string outputpath) {
	const int CHUNKSIZE = 2000;

	std::FILE* file = std::fopen(outputpath.c_str(), "w");
	if (!file) throw std::runtime_error("Error opening output file");

	ProgressBar mainbar{option::BarWidth{50},
	                    option::Start{"["},
	                    option::Fill{"■"},
	                    option::Lead{"■"},
	                    option::Remainder{"-"},
	                    option::End{"]"},
	                    option::PostfixText{"Writing to " + outputpath},
	                    option::ShowElapsedTime{true},
	                    option::ForegroundColor{Color::yellow},
	                    option::MaxProgress{m.size() + 1},
	                    option::FontStyles{std::vector<FontStyle>{FontStyle::bold}}};

	auto alphabet = Alpha::getAlphabet(alpha);

	{  // write header line
		// KMER, dataset, char0, char1, char2, ..., ']'
		// we ignore the last column since it's the begin character. It can't be in the
		// counts
		std::fprintf(file, "%s", "kmer, dataset");
		for (size_t i = 0; i < alphabet.size() - 1; ++i) fprintf(file, ", %c", alphabet[i]);
		std::fprintf(file, "%s", "\n");
	}

	int c = 0;
	std::string buff;

	for (const auto& kv : m) {
		auto decoded = decodeSequenceView(kv.first, alpha);
		auto leftpad = raw_seq_t(k - decoded.size(), '[');
		std::string kmer;
		kmer.insert(kmer.end(), leftpad.begin(), leftpad.end());
		kmer.insert(kmer.end(), decoded.begin(), decoded.end());

		for (const auto& [dataset, counts] : kv.second) {
			buff += kmer;
			buff += ", " + names[dataset];
			for (size_t j = 0; j < counts.size() - 1; ++j) {
				buff += ", " + std::to_string(counts[j]);
			}
			buff += "\n";
		}

		if (++c % CHUNKSIZE == 0) {
			std::fwrite(buff.c_str(), 1, buff.size(), file);
			buff.clear();
			mainbar.tick(CHUNKSIZE);
		}
	}
	std::fwrite(buff.c_str(), 1, buff.size(), file);
	mainbar.tick(c);
	std::fclose(file);

	{  // final spinner update

		mainbar.set_option(option::ForegroundColor{Color::green});
		mainbar.set_option(option::PrefixText{"✔ KMap written to " + outputpath});
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

// outputs in sparse matrix format:
// k-mer, sparse matrix indices, sparse matrix values
// the sparse matrix contains the count of every next letter for every dataset
inline void dumpMultiMap_sparseMatrix(const multikmap_t& m, int k, Alphabet alpha,
                                      std::string outputpath) {
	const int CHUNKSIZE = 500;

	std::FILE* file = std::fopen(outputpath.c_str(), "w");
	if (!file) throw std::runtime_error("Error opening output file");

	ProgressBar mainbar{option::BarWidth{50},
	                    option::Start{"["},
	                    option::Fill{"■"},
	                    option::Lead{"■"},
	                    option::Remainder{"-"},
	                    option::End{"]"},
	                    option::PostfixText{"Writing to " + outputpath},
	                    option::ShowElapsedTime{true},
	                    option::ForegroundColor{Color::yellow},
	                    option::MaxProgress{m.size() + CHUNKSIZE},
	                    option::FontStyles{std::vector<FontStyle>{FontStyle::bold}}};

	auto alphabet = Alpha::getAlphabet(alpha);

	// write header line
	std::fprintf(file, "kmer, count_mat_indices, count_mat_values\n");

	int c = 0;
	std::string buff;

	for (const auto& kv : m) {  // kv: kmer -> vec of char counts per dataset
		auto decoded = decodeSequenceView(kv.first, alpha);
		auto leftpad = raw_seq_t(k - decoded.size(), '[');
		std::string kmer;
		kmer.insert(kmer.end(), leftpad.begin(), leftpad.end());
		kmer.insert(kmer.end(), decoded.begin(), decoded.end());

		buff += kmer + "; [";
		std::vector<size_t> allcounts;
		for (const auto& [dataset, counts] : kv.second) {
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
	}
	std::fwrite(buff.c_str(), 1, buff.size(), file);
	mainbar.tick(c);
	std::fclose(file);

	{  // final progress update
		mainbar.set_option(option::ForegroundColor{Color::green});
		mainbar.set_option(option::PrefixText{"✔ KMap written to " + outputpath});
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
