#include <cxxopts.hpp>
#include <indicators.hpp>
#include <iostream>
#include <random>
#include <tinypool.hpp>
#include <unordered_map>
#include <utility>
#include <vector>

#include "utils.hpp"

using namespace indicators;

// computeMultiKMap takes a vector of encoded sequences
// first letter of an ecoded sequence is the prepend symbol, last is the append symbol
// there's no need for actually prepending/appending all k symbols
// we also skip sequences shorter than k
void computeMultiKMap(const std::vector<encoded_seq_t>& seqs, size_t datasetId, int k,
                      Alphabet alpha, size_t lowerBound, size_t upperBound,
                      multikmap_t& res) {
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
				seqview_t sv{&seqs[s][left], l};
				auto& datasetCounts = res[sv][datasetId];
				if (datasetCounts.size() == 0)
					datasetCounts = std::vector<size_t>(ALPHABET_SIZE, 0);
				datasetCounts[nxtchar]++;
			}
		}
	}
}

void computeDatasetMultiKmap(const std::vector<dataset_t>& datasets, size_t datasetId,
                             multikmap_t& outmap, int K, Alphabet alpha, size_t nbThreads,
                             size_t nbSeqPerTask, ProgressBar& mainbar) {
	if (nbThreads > 1) {
		TinyPool::ThreadPool tp{nbThreads};
		std::vector<multikmap_t> allmaps(tp.nThreads);
		for (size_t i = 0; i < datasets[datasetId].size(); i += nbSeqPerTask) {
			tp.push_work([i, K, nbSeqPerTask, alpha, datasetId, nbDatasets = datasets.size(),
			              &datasets, &allmaps, &mainbar](size_t procId) {
				auto& dataset = datasets[datasetId];
				const size_t lower = i;
				const size_t upper = std::min(i + nbSeqPerTask, dataset.size());
				computeMultiKMap(dataset, datasetId, K, alpha, lower, upper, allmaps[procId]);
				int ticksize = 0;
				for (size_t j = lower; j < upper; ++j) ticksize += dataset[j].size();
				mainbar.tick(ticksize);
			});
		}
		tp.waitAll();
		mergeMultiMaps(outmap, allmaps);
	} else {
		computeMultiKMap(datasets[datasetId], datasetId, K, alpha, 0,
		                 datasets[datasetId].size(), outmap);
		for (const auto& d : datasets[datasetId]) mainbar.tick(d.size());
	}
}

int main(int argc, char** argv) {
	size_t nbSeqPerTask = 10;
	int K = 40;
	unsigned int nbSeqThreads = 1;
	unsigned int nbGlobalThreads = 1;
	bool progress = true;
	Alphabet alpha = Alphabet::dna;
	std::string inputFile;
	std::string outputFile;

	enum class OutFormat { multimap, sparsematrix };
	OutFormat outformat = OutFormat::sparsematrix;

	cxxopts::Options options(
	    "kmap", "computes all kmers and their associated next character count");

	options.add_options()("k", "k as in K-mers, i.e their size",
	                      cxxopts::value<int>(K)->default_value("30"))(
	    "i,input",
	    "input file listing each dataset's path, one per line. Paths should be either "
	    "absolute or relative to the location of the input file itself",
	    cxxopts::value<std::string>(inputFile))(
	    "o,output", "output file path",
	    cxxopts::value<std::string>(outputFile)->default_value("kmap_out.csv"))(
	    "f,output-format", "output format (multimap, sparsematrix). Default = sparsematrix",
	    cxxopts::value<std::string>()->default_value("sparsematrix"))(
	    "p,progress", "show progress output",
	    cxxopts::value<bool>(progress)->default_value("true"))(
	    "a,alphabet", "alphabet: dna (default), rna or protein",
	    cxxopts::value<std::string>()->default_value("dna"))(
	    "t,seqthreads",
	    "Number of threads to deploy when processing a dataset. Total number of threads = "
	    "--nparallel x --seqthreads",
	    cxxopts::value<unsigned int>(nbSeqThreads)->default_value("1"))(
	    "n,nparallel",
	    "How many datasets to process in parallel. Total number of threads = "
	    "--nparallel x --seqthreads",
	    cxxopts::value<unsigned int>(nbGlobalThreads)->default_value("1"))(
	    "c,chunksize", "number of sequences per subtask. Only useful with multithreaded.",
	    cxxopts::value<size_t>(nbSeqPerTask)->default_value("20"))("h,help", "Print usage");

	auto opt = options.parse(argc, argv);

	if (opt.count("help")) {
		std::cout << options.help() << std::endl;
		exit(0);
	}

	if (!opt.count("input")) {
		std::cerr << "ERROR: An input file containing the list of datasets is required! Call "
		             "with --help for list of options."
		          << std::endl;
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

	if (opt["output-format"].as<std::string>() == "multimap") {
		outformat = OutFormat::multimap;
	} else if (opt["output-format"].as<std::string>() == "sparsematrix") {
		outformat = OutFormat::sparsematrix;
	} else {
		std::cerr << "Invalid output format \"" << opt["output-format"].as<std::string>()
		          << "\"; should be multimap or sparsematrix" << std::endl;
		exit(1);
	}

	show_console_cursor(false);

	auto [alldatasets, allnames] = readDatasets(inputFile, alpha);

	size_t totalTicks = 0;
	for (const auto& d : alldatasets)
		for (const auto& s : d) totalTicks += s.size();

	ProgressBar mainbar{option::BarWidth{50},
	                    option::Start{"["},
	                    option::Fill{"■"},
	                    option::Lead{"■"},
	                    option::Remainder{"-"},
	                    option::End{"]"},
	                    option::PostfixText{"Hashing datasets"},
	                    option::ShowElapsedTime{true},
	                    option::ForegroundColor{Color::cyan},
	                    option::MaxProgress{totalTicks + 1},
	                    option::FontStyles{std::vector<FontStyle>{FontStyle::bold}}};

	// set up thread pool
	TinyPool::ThreadPool tp{nbGlobalThreads};

	// std::vector<kmap_t> allkmaps(alldatasets.size());
	std::vector<multikmap_t> allkmaps(nbGlobalThreads);

	// compute a kmap (in allkmaps) for each dataset, using the threadpool.
	for (size_t i = 0; i < alldatasets.size(); ++i) {
		tp.push_work([i, K, nbSeqPerTask, nbSeqThreads, alpha, &alldatasets = alldatasets,
		              &allkmaps, &mainbar](size_t threadId) {
			computeDatasetMultiKmap(alldatasets, i, allkmaps[threadId], K, alpha, nbSeqThreads,
			                        nbSeqPerTask, mainbar);
		});
	}
	tp.waitAll();

	mainbar.set_option(option::ForegroundColor{Color::green});
	mainbar.set_option(option::PrefixText{"✔ All datasets hashed"});
	mainbar.set_option(option::BarWidth{0});
	mainbar.set_option(option::Fill{""});
	mainbar.set_option(option::Lead{""});
	mainbar.set_option(option::Start{""});
	mainbar.set_option(option::End{""});
	mainbar.set_option(option::ShowPercentage{false});
	mainbar.set_option(option::PostfixText{"                              "});
	mainbar.mark_as_completed();

	mergeMultiMaps(allkmaps[0], allkmaps, true, 1);
	auto& result = allkmaps[0];

	// output result
	switch (outformat) {
		case OutFormat::multimap:
			dumpMultiMap(result, allnames, K, alpha, outputFile);
			break;

		default:
			std::cerr << "Output format not recognized, using sparsematrix" << std::endl;
		case OutFormat::sparsematrix:
			dumpMultiMap_sparseMatrix(result, K, alpha, outputFile);
			break;
	}

	std::cerr << "Done. Cleaning up memory..." << std::endl;
	show_console_cursor(true);
	return 0;
}
