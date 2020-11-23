set(ROCKSDB_ROOT_DIR ${CMAKE_BINARY_DIR}/thirdparty/rocksdb)

include(FindGit)
find_package(Git REQUIRED)
include(ExternalProject)
ExternalProject_Add(rocksdb
	GIT_REPOSITORY	https://github.com/facebook/rocksdb
	GIT_TAG		    v6.10.2
	GIT_SHALLOW		1
	GIT_PROGRESS	1
	INSTALL_DIR		${ROCKSDB_ROOT_DIR}/bin
    CMAKE_ARGS -DUSE_RTTI=1
		-DCMAKE_INSTALL_PREFIX=${ROCKSDB_ROOT_DIR}
		-DCMAKE_INSTALL_LIBDIR=lib
        -DCMAKE_CXX_STANDARD=17
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
		-DWITH_TESTS=OFF
		-DWITH_GFLAGS=OFF
		-DWITH_BENCHMARK_TOOLS=OFF
        -DWITH_TOOLS=OFF
		-DWITH_BZ2=ON
		-DWITH_LZ4=OFF
		-DWITH_SNAPPY=OFF
		-DWITH_ZLIB=OFF
        -DWITH_ZSTD=OFF
        -DCMAKE_POSITION_INDEPENDENT_CODE=True
		INSTALL_COMMAND $(MAKE) install
		)

ExternalProject_Get_Property(rocksdb BINARY_DIR)

set(ROCKSDB_LIBRARIES
	${ROCKSDB_ROOT_DIR}/lib/librocksdb.a)

set(ROCKSDB_FOUND TRUE)

set(ROCKSDB_INCLUDE_DIR
    ${ROCKSDB_ROOT_DIR}/include)
message(STATUS "Found RocksDB library: ${ROCKSDB_LIBRARIES}")
message(STATUS "Found RocksDB includes: ${ROCKSDB_INCLUDE_DIR}")

mark_as_advanced(
    ROCKSDB_ROOT_DIR
    ROCKSDB_LIBRARIES
    ROCKSDB_INCLUDE_DIR
)
