#pragma once
#include <array>
#include <cassert>
#include <sstream>
#include <vector>

// static_for helper
template <std::size_t N> struct num { static const constexpr auto value = N; };
template <class F, std::size_t... Is>
void constexpr static_for(F func, std::index_sequence<Is...>) {
	(func(num<Is>{}), ...);
}
template <std::size_t N, typename F> constexpr void static_for(F func) {
	static_for(func, std::make_index_sequence<N>());
}

template <typename T>
std::ostream& operator<<(std::ostream& output, std::vector<T> const& values) {
	for (auto const& value : values) output << value << " ";
	return output;
}

template <typename T, size_t N>
std::ostream& operator<<(std::ostream& output, std::array<T, N> const& values) {
	for (auto const& value : values) output << value << " ";
	return output;
}

inline std::ostream& operator<<(std::ostream& output, std::vector<int8_t> const& values) {
	for (auto const& value : values) output << (int)value << " ";
	return output;
}

// Alphabet stuff
enum Alphabet { dna, rna, protein, size };
struct Alpha {
	// list of alphabets
	static const constexpr std::tuple<std::array<char, 6>, std::array<char, 6>,
	                                  std::array<char, 22>>
	    alphabets = {{'A', 'C', 'G', 'T', ']', '['},
	                 {'A', 'C', 'G', 'U', ']', '['},
	                 {'A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L',
	                  'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', ']', '['}};

	static const constexpr size_t alphabetSizes[] = {6, 6, 22};

	// array containing encoding char -> int, for each alphabet
	std::array<std::array<int8_t, 128>, Alphabet::size> a;

	constexpr Alpha() : a() {
		static_for<Alphabet::size>([&](auto alpha) {
			for (int i = 0; i < 128; ++i) a[alpha.value][i] = -1;
			for (size_t i = 0; i < std::get<alpha.value>(alphabets).size(); ++i) {
				a[alpha.value][std::get<alpha.value>(alphabets)[i]] = i;
				a[alpha.value][std::get<alpha.value>(alphabets)[i] + 32] = i;  // lowercase
			}
		});
	}

	constexpr inline int8_t encode(const char& c, const Alphabet alpha) const {
		return a[alpha][c];
	}

	template <Alphabet alpha> constexpr inline char decode(const int8_t& c) const {
		if (c < 0) return '*';
		return std::get<alpha>(alphabets)[c];
	}
};
