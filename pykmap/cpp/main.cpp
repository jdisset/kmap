#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cassert>
#include <rocksdb.hpp>
#include <utils.hpp>
#include <vector>

namespace py = pybind11;

using mat2d_t = std::vector<std::vector<uint64_t>>;

struct KmapDB {
	RocksDB rdb;

	KmapDB(const std::string& dbpath) {
		rdb.opts.error_if_exists = false;
		rdb.opts.create_if_missing = false;
		rdb.openRO(dbpath);
	}

	mat2d_t get(const std::string& kmer, size_t NDATASETS) {
		std::string value;
		rdb.db->Get(rocksdb::ReadOptions(), kmer, &value);
		if (value.size() == 0) return mat2d_t();
		datacount_t dcounts = deserialize(value);
		size_t NALPHA = dcounts.begin()->second.size();
		mat2d_t res(NDATASETS, std::vector<uint64_t>(NALPHA, 0));
		for (auto& [d, counts] : dcounts) res[d] = counts;
		return res;
	}
};

PYBIND11_MODULE(pykmapcore, m) {
	m.doc() = "PyKMap core extension";
	py::class_<KmapDB>(m, "KmapDB")
	    .def(py::init<const std::string&>())
	    .def("get", &KmapDB::get);
}
