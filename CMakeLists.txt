
project(DistanceField)

set(CMAKE_CXX_STANDARD 11)

add_library(distance_field
    distance_field.hpp
    distance_field.cpp
)

add_executable(test_distance_field
    distance_field.hpp
    distance_field.cpp
	bitmap.hpp
    main.cpp
)

