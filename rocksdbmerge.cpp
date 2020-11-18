#include "rocksdbmerge.h"

bool kmerMergeOperator::Merge([[maybe_unused]] const Slice& key,
                              const Slice* existingValue, const Slice& value,
                              std::string* newValue,
                              [[maybe_unused]] rocksdb::Logger* logger) const {

	datacount_t newMap = deserialize(value.ToString());

	if (existingValue) {
		datacount_t existingMap = deserialize((*existingValue).ToString());
		for (const auto& [k, v] : existingMap) addCounts(newMap[k], v);
	}

	*newValue = serialize(newMap);
	return true;  // no error
}

const char* kmerMergeOperator::Name() const { return "kmerMergeOperator"; }
