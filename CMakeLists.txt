
project(DistanceField)

set(CMAKE_CXX_STANDARD 11)

find_package(OpenCV)

add_library(distance_field
    distance_field.hpp
    distance_field.cpp
)

add_executable(test_distance_field
    main.cpp
)

target_link_libraries(test_distance_field
    distance_field
)