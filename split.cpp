#include <cxxopts.hpp>
#include <indicators.hpp>
#include <random>
#include <string>
#include <vector>

#include "utils.hpp"

template <typename R>
std::size_t chooseOutput(const std::vector<std::pair<std::string, double>> probs,
                         R& rnd) {
	assert(probs.size() > 0);
	auto normalizedProbs = probs;
	double sum = 0;
	for (const auto& p : probs) sum += p.second;
	if (sum == 0) return 0;
	for (auto& p : normalizedProbs) p.second /= sum;
	std::uniform_real_distribution<double> dChoice(0, 1);
	double choice = dChoice(rnd);
	for (size_t i = 0; i < normalizedProbs.size(); ++i) {
		choice -= normalizedProbs[i].second;
		if (choice <= 0) return i;
	}
	assert(type <= 0);  // should never be reached
	return 0;
}

void splitFasta(const std::string& filepath,
                const std::vector<std::pair<std::string, double>>& outputs) {
	std::random_device rd;
	std::mt19937 gen(rd());

	fs::path p{filepath};
	std::ifstream file(filepath.c_str());
	std::vector<std::ofstream> outFiles;
	for (const auto& o : outputs)
		outFiles.push_back(std::ofstream(fs::path(o.first).c_str()));

	auto saveBlock = [&](const std::string& block) {
		outFiles[chooseOutput(outputs, gen)] << block;
	};

	std::string currentBlock;
	for (std::string line{}; std::getline(file, line);) {
		if (line[0] == '>' || line[0] == '@') {
			saveBlock(currentBlock);
			currentBlock.clear();
		}
		currentBlock += line + "\n";
	}
	saveBlock(currentBlock);
}

int main(int argc, char** argv) {
	std::string inputFile;
	std::vector<std::string> outputStr;
	std::vector<std::pair<std::string, double>> outputs;
	cxxopts::Options options("fasplit", "splits fasta or fastq file in 2");

	options.add_options()("i,input", "input fasta or fastq file",
	                      cxxopts::value<std::string>(inputFile))(
	    "o,output",
	    "one output file path and its ratio, separated by a comma, syntax: "
	    "'filename,ratio' "
	    "where filename is the path to the output fasta file and raio is the percentage "
	    "(a "
	    "real number btwn 0 and 100) of sequences to sample from the input file",
	    cxxopts::value<std::vector<std::string>>(outputStr));

	auto opt = options.parse(argc, argv);
	for (auto& n : outputStr) {
		size_t p = n.find_first_of(',');
		std::string name = n.substr(0, p);
		double ratio = std::stof(n.substr(p + 1, n.size() - p));
		outputs.push_back({name, ratio});
	}

	splitFasta(inputFile, outputs);
	std::cerr << "OK" << std::endl;
	return 0;
}
