#pragma once

#include <filesystem>
#include <fstream>
#include <indicators.hpp>
#include <random>
#include <robinhood.hpp>
#include <unordered_map>
#include <xhash.hpp>

#include "alphabet.hpp"

using raw_seq_t = std::vector<char>;
using encoded_seq_t = std::vector<int8_t>;
using namespace indicators;

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

using kmap_t =
    robin_hood::unordered_map<SeqView, std::vector<size_t>, seqview_hash, seqview_equal>;

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

inline std::vector<encoded_seq_t> readFasta(const std::string& filepath,
                                            const Alphabet& alpha) {
	std::filesystem::path p{filepath};
	ProgressSpinner spinner{option::PostfixText{std::string("Reading ") + p.u8string()},
	                        option::ForegroundColor{Color::yellow},
	                        option::FontStyles{std::vector<FontStyle>{FontStyle::bold}}};

	size_t filesize = std::filesystem::file_size(p);

	std::vector<encoded_seq_t> res;

	std::ifstream file(filepath.c_str());
	std::string line;

	size_t totalBytesRead = 0;
	raw_seq_t currentSeq{};

	while (std::getline(file, line)) {
		if (line[0] == '>') {
			totalBytesRead += line.size() + currentSeq.size();
			if (currentSeq.size() > 0) {
				res.push_back(encodeSequence(currentSeq, alpha));
				currentSeq = raw_seq_t();
				spinner.set_progress(
				    static_cast<int>(((float)totalBytesRead / (float)filesize) * 100));
			}
		} else {
			currentSeq.insert(currentSeq.end(), line.begin(), line.end());
		}
	}
	spinner.set_option(option::ForegroundColor{Color::green});
	spinner.set_option(option::PrefixText{"✔"});
	spinner.set_option(option::ShowSpinner{false});
	spinner.set_option(option::ShowPercentage{false});
	spinner.set_option(option::PostfixText{std::string("Loaded ") +
	                                       std::to_string(res.size()) +
	                                       std::string(" sequences from ") + p.u8string()});
	spinner.mark_as_completed();

	return res;
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

inline void dump(const kmap_t& m, int k, Alphabet alpha, std::string outputpath) {
	const int CHUNKSIZE = 10000;
	std::ofstream file(outputpath);
	if (file.fail()) throw std::runtime_error("Error opening output file");
	ProgressSpinner spinner{option::PostfixText{"Writing to " + outputpath},
	                        option::ForegroundColor{Color::yellow},
	                        option::FontStyles{std::vector<FontStyle>{FontStyle::bold}}};
	int c = 0;

	auto alphabet = Alpha::getAlphabet(alpha);
	file << "kmer";
	// we ignore the last column since it's the begin character. It can't be in the counts
	for (size_t i = 0; i < alphabet.size() - 1; ++i) file << ", " << alphabet[i];
	file << std::endl;
	std::stringstream line;
	for (const auto& kv : m) {
		raw_seq_t leftPad{};
		raw_seq_t rightPad{};
		auto decoded = decodeSequenceView(kv.first, alpha);
		if (decoded[0] == '[')
			leftPad = raw_seq_t(k - decoded.size(), '[');
		else if (decoded[decoded.size() - 1] == ']')
			rightPad = raw_seq_t(k - decoded.size(), ']');

		line << leftPad << decoded << rightPad;
		for (size_t i = 0; i < kv.second.size() - 1; ++i) line << ", " << kv.second[i];
		line << "\n";

		if (++c % CHUNKSIZE == 0) {
			file << line.str() << std::flush;
			line = std::stringstream();
			spinner.set_progress(static_cast<int>(((float)c / (float)m.size()) * 100));
		}
	}
	file << line.str() << std::flush;
	spinner.set_option(option::ForegroundColor{Color::green});
	spinner.set_option(option::PrefixText{"✔"});
	spinner.set_option(option::ShowSpinner{false});
	spinner.set_option(option::ShowPercentage{false});
	spinner.set_option(option::PostfixText{"All k-mers written!"});
	spinner.mark_as_completed();
}

inline void collapseMaps(std::vector<kmap_t>& maps) {
	ProgressSpinner spinner{option::PostfixText{"Merging all K-Maps"},
	                        option::ForegroundColor{Color::yellow},
	                        option::FontStyles{std::vector<FontStyle>{FontStyle::bold}}};

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
		spinner.set_progress(static_cast<int>(((float)m / (float)maps.size()) * 100));
	}
	maps.resize(1);

	spinner.set_option(option::ForegroundColor{Color::green});
	spinner.set_option(option::PrefixText{"✔"});
	spinner.set_option(option::ShowSpinner{false});
	spinner.set_option(option::ShowPercentage{false});
	spinner.set_option(option::PostfixText{"All maps merged!"});
	spinner.mark_as_completed();
}
