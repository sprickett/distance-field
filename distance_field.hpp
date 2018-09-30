#pragma once
#include <vector>


class DistanceFieldGenerator
{
public:
	bool operator()(
		const unsigned char *binary_input,
		unsigned char *greyscale_output,
		double max_distance,
		int width,
		int height,
		int input_stride,
		int output_stride)
	{
		return v2(binary_input, greyscale_output, max_distance, width, height, input_stride, output_stride);
	}

	bool v1(
		const unsigned char *binary_input,
		unsigned char *greyscale_output,
		double max_distance,
		int width,
		int height,
		int input_stride,
		int output_stride);

	bool v2(
		const unsigned char *binary_input,
		unsigned char *greyscale_output,
		double max_distance,
		int width,
		int height,
		int input_stride,
		int output_stride);
		
private:
	struct Delta
	{
		typedef short type; 
		type x, y;
	};
	std::vector<Delta> deltas_;

};