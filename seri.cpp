#include <charconv>
#include <iostream>
#include <random>
#include <string_view>

#include "utils.hpp"

#define REQUIRE(expr)                           \
	if (!(expr)) {                                \
		std::cerr << "REQUIRE FAILED" << std::endl; \
		exit(1);                                    \
	}

int main(int, char**) {
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> maxdatasets(0, 10);
	std::uniform_int_distribution<size_t> val(0, std::numeric_limits<uint64_t>::max());
	const int NTEST = 10;
	for (int n = 0; n < NTEST; ++n) {
		datacount_t d;
		int S = maxdatasets(gen);
		const int NLETTERS = maxdatasets(gen);
		for (int i = 0; i < S; ++i) {
			std::vector<uint64_t> counts(NLETTERS);
			for (auto& v : counts) v = val(gen);
			d[i] = counts;
		}
		auto str = serialize(d);
		auto d2 = deserialize(str);
		auto str2 = serialize(d2);
		std::cerr << str << std::endl << std::endl;
		std::cerr << str2 << std::endl;
		REQUIRE(str == str2);
	}

	/* std::string test = "12344|non sense lol wliejf iwew|w elkfh we";*/
	// splt(std::string_view(test), ',',
	//[](const auto& sv, size_t i) { std::cerr << i << ":" << sv << std::endl; });
	// size_t v;
	// std::cerr << firstFromChars(std::string_view(test), '|', v) << std::endl;
	// std::cerr << "v = " << v << std::endl;

	return 0;
}
