cmake_minimum_required (VERSION 3.0)
project(DistanceField)

set(CMAKE_CXX_STANDARD 14)

include_directories(
	t-map-2d
)

add_library(distance_field
    distance_field.hpp
    distance_field.cpp
)

add_executable(test_distance_field
	t-map-2d/tmap2d.hpp
    distance_field.hpp
    distance_field.cpp
	bitmap.hpp
    main.cpp
)

