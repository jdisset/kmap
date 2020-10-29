#include <chrono>
#include <random>

#include "utils.hpp"

#define CATCH_CONFIG_MAIN
#include <catch.hpp>

bool isEq(const datacount_t& a, const datacount_t& b) {
	if (a.size() != b.size()) return false;
	for (const auto& [k, v] : a)
		if (!b.count(k) || a.at(k) != b.at(k)) return false;
	return true;
}

TEST_CASE("datacount_t back & forth serialization", "[serialization]") {
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> maxdatasets(0, 100);
	std::uniform_int_distribution<size_t> val(0, std::numeric_limits<uint64_t>::max());
	const int NTEST = 100;
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
		REQUIRE(isEq(d, deserialize(str)));
	}
}
