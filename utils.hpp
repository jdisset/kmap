#pragma once

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <indicators.hpp>
#include <memory>
#include <random>
#include <robinhood.hpp>
#include <unordered_map>
#include <xhash.hpp>

#include "alphabet.hpp"

using raw_seq_t = std::vector<char>;
using encoded_seq_t = std::vector<int8_t>;
using dataset_t = std::vector<encoded_seq_t>;
using namespace indicators;
namespace fs = std::filesystem;

// /!\ CAREFUL
// we are using pointers to the list of sequences to avoid copies
// but it means the references to the sequences MUST BE KEPT VALID
// i.e no insertion/deletions in the sequence datastructure.
struct SeqView {
	const int8_t* begin = nullptr;
	size_t length;
	const int8_t& operator[](size_t i) const { return *(begin + i); }
};

struct seqview_hash {
	size_t operator()(const SeqView& s) const { return xxh::xxhash<64>(s.begin, s.length); }
};

struct seqview_equal {
	bool operator()(const SeqView& s1, const SeqView& s2) const {
		if (s1.length != s2.length) return false;
		for (size_t i = 0; i < s1.length; ++i) {
			if (s1[i] != s2[i]) return false;
		}
		return true;
	}
};

template <typename V>
using seqmap_t = robin_hood::unordered_map<SeqView, V, seqview_hash, seqview_equal>;

using kmap_t = seqmap_t<std::vector<size_t>>;  // SeqView -> counts

using multikmap_t =
    seqmap_t<std::vector<std::vector<size_t>>>;  // SeqView -> counts per dataset

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

inline raw_seq_t decodeSequenceView(const SeqView& s, const Alphabet& alpha) {
	const auto decode = (alpha == Alphabet::dna) ?
	    [](const int8_t& c) -> char { return alphaMap.decode<Alphabet::dna>(c); } :
	    ((alpha == Alphabet::rna) ?
	     [](const int8_t& c) -> char { return alphaMap.decode<Alphabet::rna>(c); } :
	     [](const int8_t& c) -> char { return alphaMap.decode<Alphabet::protein>(c); });
	raw_seq_t r(s.length);
	for (size_t i = 0; i < s.length; ++i) r[i] = decode(s[i]);
	return r;
}

template <typename B>
std::vector<encoded_seq_t> readFasta(const std::string& filepath, const Alphabet& alpha,
                                     B* progress) {
	fs::path p{filepath};

	std::vector<encoded_seq_t> res;

	std::ifstream file(filepath.c_str());
	std::string line;

	size_t totalBytesRead = 0;
	raw_seq_t currentSeq{};

	auto recordSeq = [&]() {
		if (currentSeq.size() > 0) {
			res.push_back(encodeSequence(currentSeq, alpha));
			currentSeq = raw_seq_t();
			progress->tick(totalBytesRead);
			totalBytesRead = 0;
		}
	};

	while (std::getline(file, line)) {
		totalBytesRead += line.size();
		if (line[0] == '>') {
			recordSeq();
		} else {
			currentSeq.insert(currentSeq.end(), line.begin(), line.end());
		}
	}
	recordSeq();

	return res;
}

std::pair<std::vector<dataset_t>, std::vector<std::string>> readDatasets(
    std::string inputFile, Alphabet alpha) {
	// inputFile contains the list of all datasets files, relative to its location

	fs::path filePath{fs::absolute(fs::path(inputFile))};
	auto p = filePath.parent_path();

	std::vector<fs::path> datasets;
	std::ifstream file(filePath.c_str());
	std::string line;

	size_t totalFileSize = 0;
	while (std::getline(file, line)) {
		fs::path datapath{p};
		datapath /= line;
		datasets.push_back(datapath);
		totalFileSize += fs::file_size(datapath);
	}

	std::vector<dataset_t> alldatasets;
	std::vector<std::string> allnames;

	ProgressBar mainbar{option::BarWidth{50},
	                    option::Start{"["},
	                    option::Fill{"■"},
	                    option::Lead{"■"},
	                    option::Remainder{"-"},
	                    option::End{"]"},
	                    option::PostfixText{"Loading datasets"},
	                    option::ShowElapsedTime{true},
	                    option::ForegroundColor{Color::cyan},
	                    option::MaxProgress{totalFileSize},
	                    option::FontStyles{std::vector<FontStyle>{FontStyle::bold}}};

	for (size_t i = 0; i < datasets.size(); ++i) {
		allnames.push_back(datasets[i].filename());
		mainbar.set_option(option::PostfixText{"Loading " + allnames[i]});
		alldatasets.push_back(readFasta(datasets[i].c_str(), alpha, &mainbar));
	}

	mainbar.set_option(option::ForegroundColor{Color::green});
	mainbar.set_option(option::PrefixText{"✔ " + std::to_string(alldatasets.size()) +
	                                      " datasets loaded in memory"});
	mainbar.set_option(option::BarWidth{0});
	mainbar.set_option(option::Fill{""});
	mainbar.set_option(option::Lead{""});
	mainbar.set_option(option::Start{""});
	mainbar.set_option(option::End{""});
	mainbar.set_option(option::ShowPercentage{false});
	mainbar.set_option(option::PostfixText{"                              "});
	mainbar.mark_as_completed();

	return {alldatasets, allnames};
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

inline void collapseMaps(std::vector<kmap_t>& maps) {
	kmap_t& res = maps[0];
	for (size_t m = 1; m < maps.size(); ++m) {
		for (auto& kv : maps[m]) {
			if (res.count(kv.first)) {
				auto& r = res[kv.first];
				assert(res[kv.first].size() == kv.second.size());
				for (size_t i = 0; i < kv.second.size(); ++i) r[i] += kv.second[i];
			} else
				res.insert(std::move(kv));
		}
	}
	maps.resize(1);
}

// merge multiple kmaps into one multikmap
inline multikmap_t mergeMaps(const std::vector<kmap_t>& maps) {
	ProgressSpinner spinner{
	    option::PostfixText{"Merging " + std::to_string(maps.size()) + " K-Maps"},
	    option::ForegroundColor{Color::yellow}, option::ShowElapsedTime{true},
	    option::FontStyles{std::vector<FontStyle>{FontStyle::bold}}};

	multikmap_t res;  // seq -> counts for each dataset
	for (size_t m = 0; m < maps.size(); ++m) {
		for (auto& kv : maps[m]) {
			if (!res.count(kv.first))
				res[kv.first] = std::vector<std::vector<size_t>>(maps.size());
			res[kv.first][m] = kv.second;
		}
		spinner.set_progress(
		    std::min(100, static_cast<int>(((float)m / (float)maps.size()) * 100)));
	}

	{  // spinner update
		spinner.set_option(option::ForegroundColor{Color::green});
		spinner.set_option(option::PrefixText{"✔ " + std::to_string(maps.size()) +
		                                      " kmaps merged (" + std::to_string(res.size()) +
		                                      " entries)"});
		spinner.set_option(option::ShowSpinner{false});
		spinner.set_option(option::ShowPercentage{false});
		spinner.set_option(option::PostfixText{""});
		spinner.mark_as_completed();
	}
	return res;
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
	                    option::ForegroundColor{Color::cyan},
	                    option::MaxProgress{m.size() + 1},
	                    option::FontStyles{std::vector<FontStyle>{FontStyle::bold}}};

	auto alphabet = Alpha::getAlphabet(alpha);

	{  // write header line
		// KMER, dataset, char0, char1, char2, ..., ']'
		// we ignore the last column since it's the begin character. It can't be in the counts
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

		assert(names.size() == kv.second.size());

		for (size_t i = 0; i < kv.second.size(); ++i) {  // for each dataset
			auto& counts = kv.second[i];
			if (counts.size() > 0) {  // not all datasets will have this kmer
				buff += kmer;
				buff += ", " + names[i];
				for (size_t j = 0; j < counts.size() - 1; ++j) {
					buff += ", " + std::to_string(counts[j]);
				}
				buff += "\n";
			}
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
	                    option::ForegroundColor{Color::cyan},
	                    option::MaxProgress{m.size() + 1},
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
		for (size_t i = 0; i < kv.second.size(); ++i) {  // for each dataset
			const auto& counts = kv.second[i];             // counts: = vector of char counts
			if (counts.size() > 0) {                       // if datasets has this kmer
				for (size_t j = 0; j < counts.size() - 1; ++j) {
					if (counts[j] > 0) {  // if this character appears
						buff += "[" + std::to_string(i) + "," + std::to_string(j) + "],";
						allcounts.push_back(counts[j]);
					}
				}
			}
		}
		buff.pop_back();  // removes trailing comma
		buff += "]; [";
		for (const auto& c : allcounts) buff += std::to_string(c) + ",";
		buff.pop_back();  // removes trailing comma
		buff += "]\n";
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
