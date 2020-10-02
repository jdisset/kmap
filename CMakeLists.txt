cmake_minimum_required(VERSION 3.1)

project(kmap)

if( NOT CMAKE_BUILD_TYPE )
	set( CMAKE_BUILD_TYPE Release CACHE STRING
		"Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
		FORCE )
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Sanitize modes
set(CMAKE_CXX_FLAGS_SAN " -O1 -g -fsanitize=address,undefined" CACHE STRING "Enable sanitization" FORCE)
mark_as_advanced( CMAKE_CXX_FLAGS_SAN)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -g -Wextra -pedantic")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -g -Wall -Wextra -pedantic")

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/external")

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
add_executable(kmap main.cpp)
target_link_libraries(kmap PRIVATE Threads::Threads)
