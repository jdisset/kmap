#include "rocksdb.hpp"
#include "rocksdbmerge.h"

RocksDB::RocksDB() {
	opts.error_if_exists = true;
	opts.create_if_missing = true;
	opts.merge_operator.reset(new kmerMergeOperator);
}

void RocksDB::open(const std::string& path) {
	rocksdb::DB* p;
	auto status = rocksdb::DB::Open(opts, path, &p);
	db.reset(p);
	assert(status.ok());
}

void RocksDB::add(const multikmap_t& kmap, const Alphabet& alpha) {
	for (const auto& [k, v] : kmap) {
		// TODO: fewer copies...
		auto dec = decodeSequenceView(k, alpha);
		std::string ks(dec.begin(), dec.end());
		db->Merge(writeopts, ks, serialize(v));
	}
}
