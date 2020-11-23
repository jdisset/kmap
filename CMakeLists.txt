cmake_minimum_required(VERSION 3.11)

project(kmap)

if( NOT CMAKE_BUILD_TYPE )
	set( CMAKE_BUILD_TYPE Release CACHE STRING
		"Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
		FORCE )
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

# Sanitize modes
set(CMAKE_CXX_FLAGS_SAN " -O1 -g -fsanitize=address,undefined" CACHE STRING "Enable sanitization" FORCE)
mark_as_advanced( CMAKE_CXX_FLAGS_SAN)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -g -Wextra -pedantic")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -g -Wall -Wextra -pedantic")

include(build-rocksdb)
include_directories(${ROCKSDB_INCLUDE_DIR})

find_package(BZip2 REQUIRED)
include_directories(${BZIP2_INCLUDE_DIR})

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/external")

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

add_executable(kmap main.cpp rocksdb.cpp rocksdbmerge.cpp)
add_dependencies(kmap rocksdb)
target_link_libraries(kmap PRIVATE 
	Threads::Threads 
	${ROCKSDB_LIBRARIES} ${LZ4_LIBRARIES} ${ZLIB_LIBRARIES} ${BZIP2_LIBRARIES} ${ZSTD_LIBRARIES} )

add_executable(kconv convdb.cpp rocksdb.cpp rocksdbmerge.cpp)
add_dependencies(kconv rocksdb)
target_link_libraries(kconv PRIVATE 
	Threads::Threads 
	${ROCKSDB_LIBRARIES} ${LZ4_LIBRARIES} ${ZLIB_LIBRARIES} ${BZIP2_LIBRARIES} ${ZSTD_LIBRARIES} )

add_executable(tests test.cpp)
add_dependencies(tests rocksdb)
target_link_libraries(tests PRIVATE 
	Threads::Threads
	${ROCKSDB_LIBRARIES} ${LZ4_LIBRARIES} ${ZLIB_LIBRARIES} ${BZIP2_LIBRARIES} ${ZSTD_LIBRARIES} )

