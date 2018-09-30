#pragma once
#include <vector>


class DistanceFieldGenerator
{
public:
	bool operator()(
		const unsigned char *binary_input,
		unsigned char *greyscale_output,
		double max_distance,
		unsigned width,
		unsigned height,
		unsigned input_stride,
		unsigned output_stride);

		
private:
	struct Delta
	{
		typedef short type;
		type x, y;
	};
	std::vector<Delta> deltas_;

};