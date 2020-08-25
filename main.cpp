#include <array>
#include <iostream>
#include <random>
#include <unordered_map>
#include <utility>
#include <vector>

#include "alphabet.hpp"
#include "xhash.hpp"

using raw_seq_t = std::vector<char>;
using encoded_seq_t = std::vector<int8_t>;

// BE VERY CAREFUL
// we are using references to the list of sequences to avoid copies
// but it means the references to the sequences MUST BE KEPT VALID
// i.e no insertion/deletions in the sequence datastructure.
struct SeqView {
	const int8_t* begin = nullptr;
	size_t length;
	const int8_t& operator[](size_t i) const { return *(begin + i); }
};

struct seq_hash_xh {
	size_t operator()(const SeqView& s) const { return xxh::xxhash<64>(s.begin, s.length); }
};

struct seq_equal {
	bool operator()(const SeqView& s1, const SeqView& s2) const {
		if (s1.length != s2.length) return false;
		for (size_t i = 0; i < s1.length; ++i)
			if (s1[i] != s2[i]) return false;
		return true;
	}
};

using kmap_t = std::unordered_map<SeqView, std::vector<size_t>, seq_hash_xh, seq_equal>;

// encode / decode. Important for 2 reasons:
// - normalizing uppercase & lowercases
// - encoding of a letter = index in the count array
static const constexpr Alpha alphaMap{};
encoded_seq_t encodeSequence(const raw_seq_t& s, const Alphabet& alpha) {
	encoded_seq_t r(s.size() + 2);  // don't forget prepend and append symbols
	r[0] = alphaMap.encode('[', alpha);
	for (size_t i = 1; i < r.size() - 1; ++i) {
		auto v = alphaMap.encode(s[i - 1], alpha);
		if (v < 0)
			throw std::runtime_error("Invalid character in sequence");
		else
			r[i] = v;
	}
	r[r.size() - 1] = alphaMap.encode(']', alpha);
	return r;
}

raw_seq_t decodeSequence(const encoded_seq_t& s, const Alphabet& alpha) {
	const auto decode = (alpha == Alphabet::dna) ?
	    [](const int8_t& c) -> char { return alphaMap.decode<Alphabet::dna>(c); } :
	    ((alpha == Alphabet::rna) ?
	     [](const int8_t& c) -> char { return alphaMap.decode<Alphabet::rna>(c); } :
	     [](const int8_t& c) -> char { return alphaMap.decode<Alphabet::protein>(c); });

	raw_seq_t r(s.size() - 2);
	for (size_t i = 0; i < s.size() - 2; ++i) r[i] = decode(s[i + 1]);
	return r;
}

raw_seq_t decodeSequenceView(const SeqView& s, const Alphabet& alpha) {
	const auto decode = (alpha == Alphabet::dna) ?
	    [](const int8_t& c) -> char { return alphaMap.decode<Alphabet::dna>(c); } :
	    ((alpha == Alphabet::rna) ?
	     [](const int8_t& c) -> char { return alphaMap.decode<Alphabet::rna>(c); } :
	     [](const int8_t& c) -> char { return alphaMap.decode<Alphabet::protein>(c); });

	std::cerr << "s ptr = " << std::hex << (void*)s.begin << std::endl;
	std::cerr << "s length = " << s.length << std::endl;
	raw_seq_t r(s.length);
	for (size_t i = 0; i < s.length; ++i) r[i] = decode(s[i]);
	return r;
}

// computeKMap takes a vector of encoded sequences
// first letter of an ecoded sequence is the prepend symbol, last is the append symbol
// there's no need for actually prepending/appending all k symbols
kmap_t computeKMap(std::vector<encoded_seq_t>& seqs, int k, Alphabet alpha) {
	if (k <= 0) throw std::invalid_argument("k must be > 0 ");
	const auto ALPHABET_SIZE = alphaMap.alphabetSizes[alpha];
	kmap_t res;
	for (size_t s = 0; s < seqs.size(); ++s) {
		std::cerr << "computing k maps for sequence " << s << std::endl;
		const int N = seqs[s].size();
		assert(N > 0);
		for (int i = 0; i < (int)N + k - 2; ++i) {
			const size_t left = std::max(0, i - k + 1);
			const size_t right = std::min(N - 1, i);
			const size_t nxt = std::min(N - 1, i + 1);
			const size_t l = (int)right - (int)left + 1;
			std::cerr << "left = " << left;
			std::cerr << " ; right = " << right;
			std::cerr << " ; length = " << l;
			std::cerr << " ; nxt = " << nxt << std::endl;
			std::cerr << ", next = " << alphaMap.decode<Alphabet::dna>(seqs[s][nxt])
			          << std::endl;

			SeqView sv{&seqs[s][left], l};

			if (!res.count(sv)) {
				std::cerr << "doesnt exist" << std::endl;
				std::cerr << "subseq = " << decodeSequenceView(sv, alpha);
				res[sv] = std::vector<size_t>(ALPHABET_SIZE, 0);
				std::cerr << "now c ok" << std::endl;
			}
			res[sv][nxt]++;
		}
	}
	std::cerr << "DONE" << std::endl;
	return res;
}

int main(int, char**) {
	auto s = std::vector<char>{'A', 't', 'a', 't', 'c', 'C', 'g', 't', 'a', 't', 'c', 'C'};
	std::vector<encoded_seq_t> allseqs;
	allseqs.push_back(encodeSequence(s, Alphabet::dna));
	auto res = computeKMap(allseqs, 3, Alphabet::dna);
	for (const auto& kv : res) {
		std::cerr << "Decoding keys..." << std::endl;
		std::cout << decodeSequenceView(kv.first, Alphabet::dna) << " : " << kv.second
		          << std::endl;
	}
	return 0;
}
