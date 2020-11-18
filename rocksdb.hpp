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
	RocksDB(rocksdb::Options o) : opts(o) {}

	void open(const std::string& path);

	template <typename F> void forEach(F&& f) {
		rocksdb::Iterator* it = db->NewIterator(rocksdb::ReadOptions());
		for (it->SeekToFirst(); it->Valid(); it->Next()) {
			std::forward<F>(f)(it);
		}
		assert(it->status().ok());  // Check for any errors found during the scan
		delete it;
	}

	void add(const multikmap_t& kmap, const Alphabet& alpha,
	         DynamicProgress<ProgressBar>* bars = nullptr);
};
