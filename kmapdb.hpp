#pragma once
#include <sqlite3.h>

#include <filesystem>
#include <sstream>
#include <string>

#include "utils.hpp"
namespace fs = std::filesystem;

struct KmapDB {
	sqlite3 *db = nullptr;
	std::vector<std::string> alphabet;
	size_t alphabetSize = 22;
	KmapDB() {}
	KmapDB(const raw_seq_t &alpha) {
		for (auto a : alpha) alphabet.push_back(std::string(1, a));
		alphabet[alphabet.size() - 1] = "BEGIN_CHAR";
		alphabet[alphabet.size() - 2] = "END_CHAR";
	}

	inline void openW(const std::string &dbfile, bool mergeIfExists = false) {
		if (!mergeIfExists && fs::exists(dbfile))
			throw std::runtime_error("DB file already exists");
		if (sqlite3_open(dbfile.c_str(), &db)) throw std::runtime_error(sqlite3_errmsg(db));
		createTables();
	}

	inline void openRO(const std::string &dbfile) {
		if (!fs::exists(dbfile)) throw std::runtime_error("couldn't find DB file");
		if (sqlite3_open(dbfile.c_str(), &db)) throw std::runtime_error(sqlite3_errmsg(db));
	}

	inline void createTables() {
		if (!db) throw std::invalid_argument("db pointer is null");

		std::string sql =
		    "CREATE TABLE IF NOT EXISTS kmer("
		    "id INTEGER PRIMARY KEY ,"
		    "str TEXT);"
		    "CREATE TABLE IF NOT EXISTS dataset_counts("
		    "id INTEGER PRIMARY KEY,"
		    "kmer_id INTEGER NOT NULL,";

		for (auto a : alphabet) sql += a + " INTEGER NOT NULL DEFAULT 0,";

		sql += "dataset_id INTEGER NOT NULL);";
		exec(sql);
	}

	// inline void createTables() {
	// if (!db) throw std::invalid_argument("db pointer is null");
	// std::string sql =
	//"CREATE TABLE IF NOT EXISTS kmer("
	//"id INTEGER PRIMARY KEY ,"
	//"str TEXT);"
	//"CREATE TABLE IF NOT EXISTS dataset("
	//"id INTEGER PRIMARY KEY,"
	//"path TEXT);"
	//"CREATE TABLE IF NOT EXISTS counts(k"
	//"id INTEGER PRIMARY KEY"
	//"kmer_id INTEGER NOT NULL,"
	//"dataset_id INTEGER NOT NULL,"
	//"letter INTEGER NOT NULL,"
	//"cnt INTEGER NOT NULL DEFAULT 0);";
	// exec(sql);
	//}
	//

	template <typename F> void forEach(F &&f) {
		std::string sql =
		    "SELECT * FROM dataset_counts, kmer WHERE kmer_id = kmer.id ORDER BY str, "
		    "dataset_id;";

		sqlite3_stmt *stmt;
		prepare(sql, &stmt);

		std::string currentKmer;
		int64_t currentDatasetId = -1;
		std::vector<uint64_t> cnts;
		datacount_t dtc;
		std::string kmer;
		while (sqlite3_step(stmt) == SQLITE_ROW) {
			size_t nCols = sqlite3_data_count(stmt);
			size_t alphaSize = nCols - 5;
			kmer = std::string(
			    reinterpret_cast<const char *>(sqlite3_column_text(stmt, nCols - 1)));
			if (kmer != currentKmer) {
				if (dtc.size() > 0) std::forward<F>(f)(currentKmer, dtc);
				currentKmer = kmer;
				dtc = datacount_t();
			}
			int64_t dataset_id = sqlite3_column_int64(stmt, nCols - 3);
			if (dataset_id != currentDatasetId) cnts = std::vector<uint64_t>(alphaSize);
			currentDatasetId = dataset_id;
			for (size_t i = 2; i < alphaSize + 2; ++i) {
				// merges (adds up) if currentDatasetId didnt change
				cnts[i - 2] += sqlite3_column_int64(stmt, i);
			}
			dtc[dataset_id] = cnts;
		}
		if (dtc.size() > 0) std::forward<F>(f)(kmer, dtc);
		sqlite3_finalize(stmt);
	}

	void saveKmer(const std::string &k) {
		std::string sql = "INSERT OR IGNORE INTO kmer(str) VALUES ('" + k + "');";
		exec(sql);
	}

	uint64_t getRowCount() {
		sqlite3_stmt *stmt;
		int rc = sqlite3_prepare_v2(db, "SELECT count(*) FROM dataset_counts;", -1, &stmt, 0);
		rc = sqlite3_step(stmt);
		uint64_t res = sqlite3_column_int64(stmt, 0);
		sqlite3_finalize(stmt);
		return res;
	}

	uint64_t getKmerId(const std::string &k) {
		sqlite3_stmt *stmt;
		int rc = sqlite3_prepare_v2(db, "SELECT id FROM kmer where str=?;", -1, &stmt, 0);
		sqbind(stmt, 1, k);
		rc = sqlite3_step(stmt);
		uint64_t res = sqlite3_column_int64(stmt, 0);
		sqlite3_finalize(stmt);
		return res;
	}

	size_t getLastInsertId() { return sqlite3_last_insert_rowid(db); }

	// inline void add(const multikmap_t &kmap, const Alphabet &alpha, int K,
	// DynamicProgress<ProgressBar> *bars = nullptr) {
	// exec("PRAGMA synchronous = OFF");
	// exec("PRAGMA journal_mode = MEMORY");

	// PBar p(bars, kmap.size() * 2, "Writing to KmapDB");
	// std::vector<uint64_t> kmerIds;
	// kmerIds.reserve(kmap.size());
	//// 1 - save all kmers
	// exec("BEGIN TRANSACTION");
	//{
	// std::string sql =
	//"INSERT INTO kmer(str)"
	//"VALUES (?1);";
	// sqlite3_stmt *stmt;
	// prepare(sql, &stmt);
	// for (const auto &k : kmap) {
	// auto decoded = leftpad(decodeSequenceView(k.first, alpha), K);
	// std::string ks(&decoded[0], decoded.size());
	// bind(stmt, ks);
	// step(stmt);
	// p.step();
	// kmerIds.emplace_back(getLastInsertId());
	//}
	// sqlite3_finalize(stmt);
	//}
	// exec("END TRANSACTION");

	// size_t ki = 0;
	// std::string sql =
	//"INSERT INTO count(kmer_id,dataset_id,letter,cnt)"
	//"VALUES (?1, ?2, ?3, ?4);";
	// sqlite3_stmt *stmt;
	// prepare(sql, &stmt);
	// exec("BEGIN TRANSACTION");
	// for (const auto &[k, v] : kmap) {
	// auto kmerId = kmerIds[ki++];
	// for (const auto &[d, c] : v) {
	// for (unsigned int i = 0; i < c.size(); ++i) {
	// bind(stmt, kmerId, d, i, c[i]);
	// step(stmt);
	//}
	//}
	// p.step();
	//}
	// exec("END TRANSACTION");
	// p.completeMsg("Merged batch");
	// p.complete();
	//}

	inline void add(const multikmap_t &kmap, const Alphabet &alpha, int K,
	                DynamicProgress<ProgressBar> *bars = nullptr) {
		exec("PRAGMA synchronous = OFF");
		exec("PRAGMA journal_mode = MEMORY");

		PBar p(bars, kmap.size() * 3, "Writing to KmapDB");
		std::vector<uint64_t> kmerIds;
		kmerIds.reserve(kmap.size());

		// save all kmers
		exec("BEGIN TRANSACTION");
		{
			std::string sql =
			    "INSERT INTO kmer(str)"
			    "VALUES (?1);";
			sqlite3_stmt *stmt;
			prepare(sql, &stmt);
			for (const auto &k : kmap) {
				auto decoded = leftpad(decodeSequenceView(k.first, alpha), K);
				std::string ks(&decoded[0], decoded.size());
				bind(stmt, ks);
				step(stmt);
				p.step();
				kmerIds.emplace_back(getLastInsertId());
			}
			sqlite3_finalize(stmt);
		}
		exec("END TRANSACTION");

		size_t ki = 0;
		std::string sql = "INSERT INTO dataset_counts(kmer_id,";
		for (auto a : alphabet) sql += a + ",";
		sql += "dataset_id) VALUES (?1";
		for (size_t ia = 0; ia < alphabet.size(); ++ia) sql += ", ?" + std::to_string(ia + 2);
		sql += ",?" + std::to_string(alphabetSize + 2) + ");";

		exec("BEGIN TRANSACTION");
		sqlite3_stmt *stmt;
		prepare(sql, &stmt);
		for (const auto &[k, v] : kmap) {
			auto kmerId = kmerIds[ki++];
			for (const auto &[d, c] : v) {
				sqbind(stmt, 1, kmerId);
				for (unsigned int i = 0; i < c.size(); ++i) sqbind(stmt, i + 2, c[i]);
				sqbind(stmt, alphabetSize + 2, d);
				step(stmt);
			}
			p.step(2);
		}
		exec("END TRANSACTION");
		p.completeMsg("Merged batch");
		p.complete();
	}

	// ---------------------------------------
	//             GENERAL SQL UTILS
	// ---------------------------------------

	static int callbackHandler(void *actualCallback, int argc, char **argv,
	                           char **azColName) {
		auto ptr =
		    (static_cast<std::function<void(int, char **, char **)> *>(actualCallback));
		(*ptr)(argc, argv, azColName);
		return 0;
	}

	inline void exec(std::string sql) {
		if (!db) throw std::invalid_argument("db pointer is null");
		char *err_msg = nullptr;
		int rc = sqlite3_exec(db, sql.c_str(), 0, 0, &err_msg);
		if (rc != SQLITE_OK) {
			std::ostringstream errorMsg;
			errorMsg << "SQL error: " << err_msg << ". \n\nREQ = " << sql << std::endl;
			sqlite3_free(err_msg);
			throw std::invalid_argument(errorMsg.str());
		} else
			sqlite3_free(err_msg);
	}

	template <typename C> void exec(std::string sql, const C &callback) {
		if (!db) throw std::invalid_argument("db pointer is null");
		char *err_msg = nullptr;
		std::function<void(int, char **, char **)> cbackFunc = callback;
		int rc = sqlite3_exec(db, sql.c_str(), callbackHandler, &cbackFunc, &err_msg);
		if (rc != SQLITE_OK) {
			std::ostringstream errorMsg;
			errorMsg << "SQL error: " << err_msg << ". \n\nREQ = " << sql << std::endl;
			sqlite3_free(err_msg);
			throw std::invalid_argument(errorMsg.str());
		} else
			sqlite3_free(err_msg);
	}
	template <typename T, typename = std::enable_if_t<std::is_arithmetic<T>::value>>
	int sqbind(sqlite3_stmt *stmt, int index, T value) {
		if constexpr (std::is_integral<T>::value) {
			if constexpr (sizeof(T) <= sizeof(int))
				return sqlite3_bind_int(stmt, index, value);
			else
				return sqlite3_bind_int64(stmt, index, value);
		}
		return sqlite3_bind_double(stmt, index, static_cast<double>(value));
	}

	int sqbind(sqlite3_stmt *stmt, int index, const std::string &value) {
		return sqlite3_bind_text(stmt, index, value.c_str(), -1, SQLITE_TRANSIENT);
	}

	template <class... Args, size_t... Is>
	int bind_impl(sqlite3_stmt *stmt, std::index_sequence<Is...>, Args &&... args) {
		return (sqbind(stmt, Is + 1, std::forward<Args>(args)) + ...);
	}

	template <class... Args> int bind(sqlite3_stmt *stmt, Args &&... args) {
		int res =
		    bind_impl(stmt, std::index_sequence_for<Args...>{}, std::forward<Args>(args)...);

		return res;
	}

	void prepare(const std::string &sql, sqlite3_stmt **stmt) {
		sqlite3_prepare_v2(db, sql.c_str(), -1, stmt, 0);
	}
	void step(sqlite3_stmt *stmt) {
		sqlite3_step(stmt);
		sqlite3_clear_bindings(stmt);
		sqlite3_reset(stmt);
	}
	~KmapDB() { sqlite3_close(db); }
};
