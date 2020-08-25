#include <array>
#include <boost/container_hash/hash.hpp>
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

template <typename Container> struct container_hash {
	std::size_t operator()(Container const& c) const {
		return boost::hash_range(c.begin(), c.end());
	}
};

struct seq_hash {
	template <typename Container> size_t operator()(const Container& c) const {
		return xxh::xxhash<64>(c);
	}
};

struct seqview_hash {
	size_t operator()(const SeqView& s) const { return xxh::xxhash<32>(s.begin, s.length); }
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
    std::unordered_map<SeqView, std::vector<size_t>, seqview_hash, seqview_equal>;
// using kmap_t =
// std::unordered_map<encoded_seq_t, std::vector<size_t>, container_hash<encoded_seq_t>>;

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
	raw_seq_t r(s.length);
	for (size_t i = 0; i < s.length; ++i) r[i] = decode(s[i]);
	return r;
}

// computeKMap takes a vector of encoded sequences
// first letter of an ecoded sequence is the prepend symbol, last is the append symbol
// there's no need for actually prepending/appending all k symbols
kmap_t computeKMap(const std::vector<encoded_seq_t>& seqs, int k, Alphabet alpha) {
	if (k <= 0) throw std::invalid_argument("k must be > 0 ");
	const auto ALPHABET_SIZE = alphaMap.alphabetSizes[alpha];
	kmap_t res;
	for (size_t s = 0; s < seqs.size(); ++s) {
		const int N = seqs[s].size();
		assert(N > 0);
		for (int i = 0; i < (int)N + k - 2; ++i) {
			const size_t left = std::max(0, i - k + 1);
			const size_t right = std::min(N - 1, i);
			const size_t nxt = std::min(N - 1, i + 1);
			const size_t l = (int)right - (int)left + 1;
			SeqView sv{&seqs[s][left], l};
			// encoded_seq_t sv(l);
			// for (size_t j = 0; j < l; ++j) sv[j] = seqs[s][j + left];
			if (!res.count(sv)) res[sv] = std::vector<size_t>(ALPHABET_SIZE, 0);
			res[sv][seqs[s][nxt]]++;
		}
	}
	return res;
}

kmap_t mergeMaps(const std::vector<kmap_t>& maps) {
	kmap_t res = maps[0];
	for (size_t m = 1; m < maps.size(); ++m) {
		for (const auto& kv : maps[m]) {
			if (res.count(kv.first)) {
				auto& r = res[kv.first];
				assert(res[kv.first].size() == kv.second.size());
				for (size_t i = 0; i < kv.second.size(); ++i) {
					r[i] += kv.second[i];
				}
			} else
				res.insert(kv);
		}
	}
	return res;
}

int main(int, char**) {
	const std::vector<char> letters{'A', 'C', 'G', 'T'};
	std::random_device rd;   // Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd());  // Standard mersenne_twister_engine seeded with rd()
	std::uniform_int_distribution<> d(0, 3);
	std::vector<encoded_seq_t> allseqs;
	for (int i = 0; i < 1000; ++i) {
		auto s = std::vector<char>(1000);
		for (auto& n : s) n = letters[d(gen)];
		allseqs.push_back(encodeSequence(s, Alphabet::dna));
	}

	std::vector<kmap_t> acc;
	acc.reserve(10);

	auto start = std::chrono::steady_clock::now();
	for (int i = 0; i < 10; ++i) acc.push_back(computeKMap(allseqs, 5, Alphabet::dna));
	auto end = std::chrono::steady_clock::now();
	std::chrono::duration<double> elapsed = end - start;

	auto res = mergeMaps(acc);
	auto end2 = std::chrono::steady_clock::now();
	std::chrono::duration<double> elapsed2 = end2 - end;

	std::cerr << acc.size() << " maps populated in " << elapsed.count() << "s\n";
	std::cerr << res.size() << " entries total. Merged in " << elapsed2.count() << "s\n";
	/*for (const auto& kv : res) {*/
	// std::cerr << "Decoding keys..." << std::endl;
	//// std::cout << decodeSequenceView(kv.first, Alphabet::dna) << " : " << kv.second
	////<< std::endl;
	////
	// std::cout << kv.first << " : " << kv.second
	//<< std::endl;
	/*}*/
	return 0;
}
