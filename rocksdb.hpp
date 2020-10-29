#pragma once
#include <rocksdb/db.h>
#include <rocksdb/merge_operator.h>

#include "utils.hpp"

struct RocksDB {
	using Slice = rocksdb::Slice;  // basically a string_view
	std::shared_ptr<rocksdb::DB> db;
	rocksdb::Options opts{};
	rocksdb::WriteOptions writeopts{};

	RocksDB();

	void open(const std::string& path);

	void add(const multikmap_t& kmap, const Alphabet& alpha);
};
