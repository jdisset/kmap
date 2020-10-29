#pragma once
#include <rocksdb/db.h>
#include <rocksdb/merge_operator.h>

#include "utils.hpp"

// Custom merge operator. When storing counts for a kmer, we will check if
// there are already some counts for this kmer. If yes, we add the current counts for
// new datasets or, if there is an overlap in the new vs existing datasets, we add the
// counts.
class kmerMergeOperator : public rocksdb::AssociativeMergeOperator {
	using Slice = rocksdb::Slice;  // basically a string_view
 public:
	virtual bool Merge([[maybe_unused]] const Slice& key, const Slice* existingValue,
	                   const Slice& value, std::string* newValue,
	                   [[maybe_unused]] rocksdb::Logger* logger) const override;

	virtual const char* Name() const override;
};
