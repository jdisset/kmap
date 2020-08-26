#pragma once
#include <unordered_map>

#include "alphabet.hpp"
#include "robinhood.hpp"
#include "xhash.hpp"

using raw_seq_t = std::vector<char>;
using encoded_seq_t = std::vector<int8_t>;

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

inline void collapseMaps(std::vector<kmap_t>& maps) {
	if (maps.size() < 1) return;
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
