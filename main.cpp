#include <cxxopts.hpp>
#include <indicators.hpp>
#include <iostream>
#include <random>
#include <tinypool.hpp>
#include <unordered_map>
#include <utility>
#include <vector>

#include "utils.hpp"

// computeKMap takes a vector of encoded sequences
// first letter of an ecoded sequence is the prepend symbol, last is the append symbol
// there's no need for actually prepending/appending all k symbols
// we also skip sequences shorter than k
void computeKMap(const std::vector<encoded_seq_t>& seqs, int k, Alphabet alpha,
                 size_t lowerBound, size_t upperBound, kmap_t& res) {
	if (k <= 0) throw std::invalid_argument("k must be > 0 ");
	const auto ALPHABET_SIZE = alphaMap.alphabetSizes[alpha];
	for (size_t s = lowerBound; s < upperBound; ++s) {
		const int N = seqs[s].size();
		assert(N > 0);
		int prevUnknownCharacter = k + 1;
		for (int i = 0; i < (int)N - 1; ++i) {
			const size_t left = std::max(0, i - k + 1);
			const int nxtchar = seqs[s][i + 1];
			const size_t l = i - (int)left + 1;

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

// overload of computeKMap that builds and returns a kmap
kmap_t computeKMap(const std::vector<encoded_seq_t>& seqs, int k, Alphabet alpha) {
	kmap_t res;
	computeKMap(seqs, k, alpha, 0, seqs.size(), res);
	return res;
}

int main(int argc, char** argv) {
	size_t NB_SEQUENCES_PER_TASK = 10;
	int K = 40;
	unsigned int nbThreads = 1;
	bool progress = true;
	Alphabet alpha = Alphabet::dna;
	std::string outputFile;

	cxxopts::Options options(
	    "kmap", "computes all kmers and their associated next character count");

	options.add_options()("k", "k as in K-mers, i.e their size",
	                      cxxopts::value<int>(K)->default_value("30"))(
	    "i,input", "input fasta file path", cxxopts::value<std::string>())(
	    "o,output", "output file path",
	    cxxopts::value<std::string>(outputFile)->default_value("kmap_out.csv"))(
	    "p,progress", "show progress output",
	    cxxopts::value<bool>(progress)->default_value("true"))(
	    "a,alphabet", "alphabet: dna (default), rna or protein",
	    cxxopts::value<std::string>()->default_value("dna"))(
	    "t", "number of threads",
	    cxxopts::value<unsigned int>(nbThreads)->default_value("1"))(
	    "c,chunksize", "number of sequences per subtask. Only useful with multithreaded.",
	    cxxopts::value<size_t>(NB_SEQUENCES_PER_TASK)->default_value("20"))("h,help",
	                                                                        "Print usage");

	auto opt = options.parse(argc, argv);

	if (opt.count("help")) {
		std::cout << options.help() << std::endl;
		exit(0);
	}
	if (!opt.count("input")) {
		std::cerr << "ERROR: An input fasta file is required!!" << std::endl;
		std::cout << options.help() << std::endl;
		exit(1);
	}

	if (opt["alphabet"].as<std::string>() == "dna") {
		alpha = Alphabet::dna;
	} else if (opt["alphabet"].as<std::string>() == "rna") {
		alpha = Alphabet::rna;
	} else if (opt["alphabet"].as<std::string>() == "protein") {
		alpha = Alphabet::protein;
	} else {
		std::cerr << "Invalid alphabet; should be dna, rna or protein" << std::endl;
		exit(1);
	}

	show_console_cursor(false);

	auto allseqs = readFasta(opt["input"].as<std::string>(), alpha);

	// progress bar
	using namespace indicators;

	size_t totalTicks = 0;
	for (const auto& s : allseqs) totalTicks += s.size();

	ProgressBar bar{option::BarWidth{50},
	                option::Start{"["},
	                option::Fill{"■"},
	                option::Lead{"■"},
	                option::Remainder{"-"},
	                option::End{"]"},
	                option::PostfixText{"Computing hashes + counts"},
	                option::ShowElapsedTime{true},
	                option::ForegroundColor{Color::cyan},
	                option::MaxProgress{totalTicks},
	                option::FontStyles{std::vector<FontStyle>{FontStyle::bold}}};

	// setup parallel processing
	TinyPool::ThreadPool tp{nbThreads};
	std::vector<kmap_t> allmaps(nbThreads);

	for (size_t i = 0; i < allseqs.size(); i += NB_SEQUENCES_PER_TASK) {
		tp.push_work(
		    [i, K, NB_SEQUENCES_PER_TASK, alpha, &allseqs, &allmaps, &bar](size_t procId) {
			    const size_t lower = i;
			    const size_t upper = std::min(i + NB_SEQUENCES_PER_TASK, allseqs.size());
			    computeKMap(allseqs, K, alpha, lower, upper, allmaps[procId]);
			    int ticksize = 0;
			    for (size_t j = lower; j < upper; ++j) ticksize += allseqs[j].size();
			    bar.tick(ticksize);
		    });
	}
	tp.waitAll();

	collapseMaps(allmaps);

	dump(allmaps[0], K, alpha, outputFile);

	show_console_cursor(true);
	return 0;
}
