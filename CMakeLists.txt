cmake_minimum_required (VERSION 3.0)
project(DistanceField)

set(CMAKE_CXX_STANDARD 14)

include_directories(
	t-map-2d
	include
)

add_library( distance_field	
	src/distance_field.cpp
)


enable_testing()

add_executable( test_distance_field
	EXCLUDE_FROM_ALL
	t-map-2d/tmap2d.hpp
	include/distance_field.hpp
	src/distance_field.cpp	
	test/src/bitmap.hpp
	test/src/main.cpp
)
#target_compile_definitions( test_distance_field PRIVATE DISTANCE_FIELD_DEBUG)

#add_library( test_helper 
#	EXCLUDE_FROM_ALL
#	test/src/bitmap.hpp
#)

add_executable( test_field_error 
	EXCLUDE_FROM_ALL
	test/src/test_field_error.cpp
)
target_link_libraries( test_field_error
#	test_helper 
    distance_field
)
#add_test( field_error_test 
#	test_field_error delta_sweep 0.5 
#)

add_executable( test_bitmap 
	EXCLUDE_FROM_ALL
	
	test/src/bitmap.hpp
	test/src/test_bitmap.cpp
)
#add_test( test_16_ARGB test_bitmap "images/im_16_ARGB.bmp" )
#add_test( test_16_RGB test_bitmap "images/im_16_RGB.bmp" )
#add_test( test_16_XRGB test_bitmap "images/im_16_XRGB.bmp" )
add_test( test_24_RGB test_bitmap "images/im_24_RGB.bmp" )
#add_test( test_32_ARGB test_bitmap "images/im_32_ARGB.bmp" )
#add_test( test_32_XRGB test_bitmap "images/im_32_XRGB.bmp" )
add_test( test_i test_bitmap "images/im_i.bmp" )
add_test( test_i_RLE test_bitmap "images/im_i_RLE.bmp" )

#add_test( argb_bitmap_test test_bitmap "images/32.bmp" )
#add_test( argb_bitmap_test2 test_bitmap "images/32c.bmp" )
#add_test( rgb_bitmap_test test_bitmap "images/24.bmp" )
#add_test( greyscale_bitmap_test test_bitmap "images/8_256.bmp" )
#add_test( palette_bitmap_test test_bitmap "images/8_4.bmp" )
#add_test( mono_bitmap_test test_bitmap "images/1.bmp" )
#add_test( rgb_bitmap_test test_bitmap "images/8_256_rle.bmp" )
#add_test( greyscale_bitmap_test2 test_bitmap "images/felix512.bmp" )
#add_test( mono_bitmap_test2 test_bitmap "images/felix512_mono.bmp" )
#add_test( rle_bitmap_test2 test_bitmap "images/felix512_rle.bmp" )


add_test( win10_1_test test_bitmap "images/win10_1.bmp" )
add_test( win10_24_test test_bitmap "images/win10_24.bmp" )
add_test( win10_4_16_test test_bitmap "images/win10_4_16.bmp" )
add_test( win10_8_256_test test_bitmap "images/win10_8_256.bmp" )
