cmake_minimum_required (VERSION 3.0)
project(DistanceField)

set(CMAKE_CXX_STANDARD 14)

include_directories(
	t-map-2d
	include
)

add_library( distance_field	
	t-map-2d/tmap2d.hpp
	include/distance_field.hpp
	src/distance_field.cpp
)

add_executable( test_distance_field
	src/bitmap.hpp
	src/main.cpp
)

target_link_libraries( test_distance_field
	distance_field
)
