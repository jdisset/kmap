#include <iostream>
#include <random>
#include <unordered_map>
#include <utility>
#include <vector>

#include "tinypool.hpp"
#include "utils.hpp"

// computeKMap takes a vector of encoded sequences
// first letter of an ecoded sequence is the prepend symbol, last is the append symbol
// there's no need for actually prepending/appending all k symbols
void computeKMap(const std::vector<encoded_seq_t>& seqs, int k, Alphabet alpha,
                 size_t lowerBound, size_t upperBound, kmap_t& res) {
	if (k <= 0) throw std::invalid_argument("k must be > 0 ");
	const auto ALPHABET_SIZE = alphaMap.alphabetSizes[alpha];
	for (size_t s = lowerBound; s < upperBound; ++s) {
		const int N = seqs[s].size();
		assert(N > 0);
		int prevUnknownCharacter = k + 1;
		for (int i = 0; i < (int)N + k - 2; ++i) {
			const size_t left = std::max(0, i - k + 1);
			const size_t right = std::min(N - 1, i);
			const int nxtchar = seqs[s][std::min(N - 1, i + 1)];
			const size_t l = (int)right - (int)left + 1;

			if (nxtchar < 0)  // unknown character
				prevUnknownCharacter = 0;
			else
				++prevUnknownCharacter;

			if (prevUnknownCharacter > k) {
				SeqView sv{&seqs[s][left], l};
				if (!res.count(sv)) res[sv] = std::vector<size_t>(ALPHABET_SIZE, 0);
				res[sv][nxtchar]++;
			}
		}
	}
}

kmap_t computeKMap(const std::vector<encoded_seq_t>& seqs, int k, Alphabet alpha) {
	kmap_t res;
	computeKMap(seqs, k, alpha, 0, seqs.size(), res);
	return res;
}

int main(int, char**) {
	const size_t NB_SEQUENCES_PER_TASK = 100;
	int K = 10;

	// generate random sequences
	const int L = 100;
	const int N = 10000;
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

	// setup parallel processing
	unsigned int nbThreads = 8;
	TinyPool::ThreadPool tp{nbThreads};
	std::vector<kmap_t> allmaps(nbThreads);
	Alphabet alpha = Alphabet::dna;
	auto start = std::chrono::steady_clock::now();
	for (size_t i = 0; i < allseqs.size(); i += NB_SEQUENCES_PER_TASK) {
		tp.push_work([i, K, alpha, &allseqs, &allmaps](size_t procId) {
			const size_t lower = i;
			const size_t upper = std::min(i + NB_SEQUENCES_PER_TASK, allseqs.size());
			computeKMap(allseqs, K, alpha, lower, upper, allmaps[procId]);
		});
	}
	tp.waitAll();

	auto end = std::chrono::steady_clock::now();
	std::chrono::duration<double> elapsed = end - start;

	collapseMaps(allmaps);
	auto end2 = std::chrono::steady_clock::now();
	std::chrono::duration<double> elapsed2 = end2 - end;

	std::cerr << nbThreads << " maps populated in " << elapsed.count() << "s\n";
	std::cerr << allmaps[0].size() << " entries total. Merged in " << elapsed2.count()
	          << "s\n";
	return 0;
}
