#include <iostream>

#include "rocksdb.hpp"
#include "rocksdbmerge.h"

RocksDB::RocksDB() {
	opts.error_if_exists = true;
	opts.create_if_missing = true;
	//opts.IncreaseParallelism();
	//opts.OptimizeLevelStyleCompaction();
	opts.merge_operator.reset(new kmerMergeOperator());
}



void RocksDB::open(const std::string& path) {
	rocksdb::DB* p;
	auto status = rocksdb::DB::Open(opts, path.c_str(), &p);
	db.reset(p);
	if (!status.ok()) throw std::runtime_error(status.ToString());
}

void RocksDB::add(const multikmap_t& kmap, const Alphabet& alpha,
                  DynamicProgress<ProgressBar>* bars) {
	PBar p(bars, kmap.size(), "Writing to rocksDB");
	for (const auto& [k, v] : kmap) {
		auto dec = decodeSequenceView(k, alpha);
		Slice ks(&dec[0], dec.size());
		auto s = serialize(v);
		db->Merge(writeopts, ks, s);
		p.step();
	}
	p.completeMsg("Merged batch");
	p.complete();
}
