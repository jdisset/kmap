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
#set(CMAKE_CXX_FLAGS_SAN " -O1 -g -fsanitize=address,undefined" CACHE STRING "Enable sanitization" FORCE)
#mark_as_advanced( CMAKE_CXX_FLAGS_SAN)
#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -g -Wextra -pedantic")
#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -g -Wall -Wextra -pedantic")

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/external")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/external/sqlite")

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
#find_package(SQLite3 REQUIRED)
include_directories(${SQLITE3_INCLUDE_DIRS})

add_executable(kmap main.cpp external/sqlite/sqlite3.c)
#target_link_libraries(kmap ${SQLITE3_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(kmap PRIVATE ${CMAKE_THREAD_LIBS_INIT})

add_executable(kconv convdb.cpp external/sqlite/sqlite3.c)
#target_link_libraries(kconv PRIVATE ${SQLITE3_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(kconv PRIVATE ${CMAKE_THREAD_LIBS_INIT})

#add_executable(tests test.cpp)
#target_link_libraries(tests PRIVATE
