#include "distance_field.hpp"

#include <algorithm>

	template<typename T>
	T clamp(T value, T minimum, T maximum)
	{
			return std::min(std::max(value, minimum), maximum);
	}

	bool DistanceFieldGenerator::v1(
			const unsigned char *binary_input, 
			unsigned char *greyscale_output,
			double max_distance, 
			int width, 
			int height, 
			int input_stride,
			int output_stride)
	{
		constexpr auto max_dim = std::numeric_limits<Delta::type>::min() + 2;
		constexpr Delta infinity = { 0, max_dim };
		const unsigned char* I; // input row pointer
	
		Delta *D; // delta row pointer
		int stride; // delta stride
		int x, y; // position
		unsigned d0, d1; // square distance 
		int dx, dy; // nearest border offset

		// reallocate and initialise deltas		
		stride = width + 2;
		deltas_.clear(); // clear first to avoid any copy
		deltas_.resize(stride*(height+2), infinity); // delta x must be zero
		Delta * D0 = &deltas_[stride + 1];



		// mark immediate borders
		for (y = 0; y < height - 1; ++y)
		{
			I = &binary_input[y*input_stride];
			D = &D0[y*stride];
			for (x = 0; x < width-1; ++x)
			{
				// diagonals 
				if (I[x] != I[x + input_stride + 1])
				{
					D[x] = { 1,1 };
					D[x + stride + 1] = { -1,-1 };
				}
				if (I[x+1] != I[x + input_stride])
				{
					D[x+1] = { -1,1 };
					D[x + stride] = { 1,-1 };
				}

				// horizontals 
				if (I[x] != I[x + 1])
				{
					D[x] = { 1,0 };
					D[x + 1] = { -1,0 };
				}
				if (I[x] != I[x + input_stride])
				{
					D[x] = { 0,1 };
					D[x + stride] = { 0,-1 };
				}

			}
		}

		// handle last pixel
		{
			const unsigned char *I = &binary_input[height*input_stride-1];
			Delta *D = &D0[height*stride - 1];
			if (I[0] != I[-1])
			{
				D[-1] = { 1,0 };
				D[0] = { -1,0 };
			}
			if (I[0] != I[int(0-input_stride)])
			{
				D[0-input_stride] = { 0,1 };
				D[0] = { 0,-1 };
			}
		}

		// forward pass
		for (y = 0; y < height; ++y)
		{
			D = &D0[y*stride];
			for (x = 0; x < width; ++x)
			{
				dx = D[x].x;
				dy = D[x].y;
				d0 = dx * dx + dy * dy;

				dx = D[x - stride - 1].x - 2;
				dy = D[x - stride - 1].y - 2;
				d1 = dx * dx + dy * dy;
				if (d1 < d0)
				{
					d0 = d1;
					D[x].x = dx; 
					D[x].y = dy;
				}
		
				dx = D[x - stride].x;
				dy = D[x - stride].y - 2;
				d1 = dx * dx + dy * dy;
				if (d1 < d0)
				{
					d0 = d1;
					D[x].x = dx;
					D[x].y = dy;
				}

				dx = D[x - stride + 1].x + 2;
				dy = D[x - stride + 1].y - 2;
				d1 = dx * dx + dy * dy;
				if (d1 < d0)
				{
					d0 = d1;
					D[x].x = dx;
					D[x].y = dy;
				}

				dx = D[x - 1].x - 2;
				dy = D[x - 1].y;
				d1 = dx * dx + dy * dy; 
				if (d1 < d0)
				{
					d0 = d1;
					D[x].x = dx;
					D[x].y = dy;
				}
			}
		}

		// backward pass
		for ( y = height - 1; y >= 0; --y) 
		{
			D = &D0[y*stride];
			for (x = width - 1; x >= 0; --x)
			{
				dx = D[x].x;
				dy = D[x].y;
				d0 = dx * dx + dy * dy;

				dx = D[x + stride + 1].x + 2;
				dy = D[x + stride + 1].y + 2;
				d1 = dx * dx + dy * dy;
				if (d1 < d0)
				{
					d0 = d1;
					D[x].x = dx;
					D[x].y = dy;
				}

				dx = D[x + stride].x;
				dy = D[x + stride].y + 2;
				d1 = dx * dx + dy * dy;
				if (d1 < d0)
				{
					d0 = d1;
					D[x].x = dx;
					D[x].y = dy;
				}
				
				dx = D[x + stride - 1].x - 2;
				dy = D[x + stride - 1].y + 2;
				d1 = dx * dx + dy * dy;
				if (d1 < d0)
				{
					d0 = d1;
					D[x].x = dx;
					D[x].y = dy;
				}

				dx = D[x + 1].x + 2;
				dy = D[x + 1].y;
				d1 = dx * dx + dy * dy;
				if (d1 < d0)
				{
					d0 = d1;
					D[x].x = dx;
					D[x].y = dy;
				}
			}
		}
		//unsigned i = 0;
		//for (y = 0; y < height; ++y)
		//{
		//	Delta *D = &D0[y*stride];
		//	for (x = 0; x < width; ++x)
		//	{
		//		printf("%+3d %+3d|", clamp<int>(D[x].x,-99,+99), clamp<int>(D[x].y,-99,99));
		//	}
		//	printf("\n");
		//}
		//printf("\n");



		// mark immediate borders
		double m = 128 / (max_distance * 2.0);
		for (y = 0; y < height; ++y)
		{
			unsigned char* O = &greyscale_output[y*output_stride];
			I = &binary_input[y*input_stride];
			D = &D0[y*stride];
			for (x = 0; x < width; ++x)
			{
				dx = D[x].x;
				dy = D[x].y;
				double dist = sqrt(dx * dx + dy * dy) * 0.5;
				dist = std::min(dist, max_distance);
				dist *= 127.5 / max_distance;
				dist = I[x] != 0 ? dist : -dist;
				dist += 127.5;				
				O[x] = dist;
			}
		}

		return true;
	}

	bool DistanceFieldGenerator::v2(
		const unsigned char *binary_input,
		unsigned char *greyscale_output,
		double max_distance,
		int width,
		int height,
		int input_stride,
		int output_stride)
	{
		constexpr auto max_dim = std::numeric_limits<Delta::type>::min() + 2;
		constexpr Delta infinity = { 0, max_dim }; // delta x must be zero
		const unsigned char* I; // input row pointer

		Delta *D; // delta row pointer
		int stride; // delta stride
		int x, y; // position
		int dx, dy; // nearest border offset
		unsigned d0, d1; // square distance 

		// reallocate and initialise deltas		
		stride = width + 2;
		deltas_.clear(); // clear first to avoid any copy
		deltas_.resize(stride*(height + 2), infinity); 
		Delta * D0 = &deltas_[stride + 1];



		// mark immediate borders
		for (y = 0; y < height - 1; ++y)
		{
			I = &binary_input[y*input_stride];
			D = &D0[y*stride];
			for (x = 0; x < width - 1; ++x)
			{
				// diagonals 
				if (I[x] != I[x + input_stride + 1])
				{
					D[x] = { 1,1 };
					D[x + stride + 1] = { -1,-1 };
				}
				if (I[x + 1] != I[x + input_stride])
				{
					D[x + 1] = { -1,1 };
					D[x + stride] = { 1,-1 };
				}

				// horizontals 
				if (I[x] != I[x + 1])
				{
					D[x] = { 1,0 };
					D[x + 1] = { -1,0 };
				}
				if (I[x] != I[x + input_stride])
				{
					D[x] = { 0,1 };
					D[x + stride] = { 0,-1 };
				}

			}
		}

		// handle last pixel
		{
			const unsigned char *I = &binary_input[height*input_stride - 1];
			Delta *D = &D0[height*stride - 1];
			if (I[0] != I[-1])
			{
				D[-1] = { 1,0 };
				D[0] = { -1,0 };
			}
			if (I[0] != I[int(0 - input_stride)])
			{
				D[0 - input_stride] = { 0,1 };
				D[0] = { 0,-1 };
			}
		}

		// look up table
		struct { int x, y, offset; } const K[] = {
			{-2, -2, -stride - 1, },
			{ 0, -2, -stride, },
			{ 2, -2, -stride + 1, },
			{-2,  0, -1, },

			{ 2,  2, stride + 1, },
			{ 0,  2, stride, },
			{-2,  2, stride - 1, },
			{ 2,  0, 1, },
		};

		// forward pass
		for (y = 0; y < height; ++y)
		{
			D = &D0[y*stride];
			for (x = 0; x < width; ++x)
			{
				dx = D[x].x;
				dy = D[x].y;
				d0 = dx * dx + dy * dy;

				for (int i = 0; i < 4; ++i)
				{
					dx = D[x + K[i].offset].x + K[i].x;
					dy = D[x + K[i].offset].y + K[i].y;
					d1 = dx * dx + dy * dy;
					if (d1 < d0)
					{
						d0 = d1;
						D[x].x = dx;
						D[x].y = dy;
					}
				}
			}
		}

		// backward pass
		for (y = height - 1; y >= 0; --y)
		{
			D = &D0[y*stride];
			for (x = width - 1; x >= 0; --x)
			{
				dx = D[x].x;
				dy = D[x].y;
				d0 = dx * dx + dy * dy;
				for (int i = 4; i < 8; ++i)
				{
					dx = D[x + K[i].offset].x + K[i].x;
					dy = D[x + K[i].offset].y + K[i].y;
					d1 = dx * dx + dy * dy;
					if (d1 < d0)
					{
						d0 = d1;
						D[x].x = dx;
						D[x].y = dy;
					}
				}
			}
		}


		// mark immediate borders
		const double m = 127.5 / (max_distance * 2.0);
		for (y = 0; y < height; ++y)
		{
			unsigned char* O = &greyscale_output[y*output_stride];
			I = &binary_input[y*input_stride];
			D = &D0[y*stride];
			for (x = 0; x < width; ++x)
			{
				dx = D[x].x;
				dy = D[x].y;
				double dist = sqrt(dx * dx + dy * dy) * m;
				dist = std::min(dist, 127.5);
				dist = I[x] != 0 ? dist : -dist;
				dist += 127.5;
				O[x] = dist;
			}
		}

		return true;
	}

