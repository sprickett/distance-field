#include <iostream>
#include <iomanip>
#include <algorithm>

#include "distance_field.hpp"
#include "bitmap.hpp"
#include "tmap2d.hpp"
void drawCircle(TMap<uint8_t>& im, int x, int y, int radius)
{
	int rows = im.height();
	int cols = im.width();

	int r2 = radius * radius;
	int xs = std::max(0, x - radius);
	int ys = std::max(0, y - radius);
	int xe = std::min(cols, x + radius);
	int ye = std::min(rows, y + radius);
	int dy, dx;
	for (int iy = ys; iy < ye; ++iy)
	{
		dy = y - iy;
		dy *= dy;
		uint8_t* B = im.ptr(iy);
		for (int ix = xs; ix < xe; ++ix)
		{
			dx = x - ix;
			if (dx*dx + dy < r2)
				B[ix] ^= 255;
		}
	}

}

const char * toString(const Bitmap::Error& error)
{
	const char * chars = "";
	switch (error)
	{
	case Bitmap::Error::NONE:
		chars = "bitmap::Error::NONE";
		break;
	case Bitmap::Error::FILE_IO_ERROR:
		chars = "bitmap::Error::FILE_IO_ERROR";
		break;
	case Bitmap::Error::UNSUPPORTED_HEADER_TYPE:
		chars = "bitmap::Error::UNSUPPORTED_HEADER_SIZE";
		break;
	case Bitmap::Error::BAD_MAGIC_BYTES:
		chars = "bitmap::Error::BAD_MAGIC_BYTES";
		break;
	case Bitmap::Error::BAD_HEADER:
		chars = "bitmap::Error::BAD_HEADER";
		break;
	case Bitmap::Error::BAD_PALETTE:
		chars = "bitmap::Error::BAD_PALETTE";
		break;
	case Bitmap::Error::BAD_BITS_PER_PIXEL:
		chars = "bitmap::Error::BAD_BITS_PER_PIXEL";
		break;
	case Bitmap::Error::BAD_IMAGE:
		chars = "bitmap::Error::BAD_IMAGE";
		break;
	}
	return chars;
}


int main(int argc, char* argv[])
{
	using namespace distance_field;
	constexpr auto SZ =256;
	TMap < float > distance, check;
	TMap<uint8_t> binary(SZ, SZ);
	binary.setTo(0);
	drawCircle(binary, 0, 0, (SZ * 4) / 5);

	//int bpp;
	//Bitmap bmp;
	//Bitmap::Error result = bmp.load("images/felix256.bmp", binary, bpp);
	//if (result != Bitmap::Error::NONE)
	//{
	//	std::cout << "failed to load bitmap: "<< toString(result) << "\n";
	//	return -1;
	//}
	//if (bpp != 8)
	//{
	//	std::cout << "only 8 bit bitmaps supported: " << toString(result) << "\n";
	//	return -1;
	//}
	
	{ 
		std::cout << "deltaSweep\n";
		TMap < std::pair<int16_t, int16_t> > deltas;
		deltaSweep(binary, deltas);
		signedDistance(binary, deltas, distance);
	}

	//{
	//	TMap < int > signed_squared_distance;
	//	dijkstra(binary, signed_squared_distance);
	//	signedDistance(signed_squared_distance, distance);
	//}

	{ 

		std::cout << "simpleList\n";
		TMap < int > signed_squared_distance;
		simpleList(binary, signed_squared_distance);
		signedDistance(signed_squared_distance, check);
	}

	
	int width = binary.width();
	int height = binary.height();
	if (distance.width() != width || distance.height() != height)
	{
		std::cout << "bad image\n";
		return -1;
	}

	std::cout << "checking error...\n";
	double max_error = -std::numeric_limits<double>::infinity();
	double sum_squares = 0.f;
	int count = 0;
	for (int y = 0; y < height; ++y)
	{
		const float *C = check.ptr(y);
		const float *D = distance.ptr(y);
		for (int x = 0; x < width; ++x)
		{			
			if (D[x] == C[x])
				continue;
			double e = fabs(double(D[x]) - C[x]);
			max_error = std::max(max_error, e);
			sum_squares += e * e;
			++count;
		}
	}
	double size = double(width) * height;
	double variance = sum_squares / size;
	double standard_deviation = std::sqrt(variance);
	
	for (int i = 1; i < argc; ++i)
		std::cout << argv[i] << "\n";
	std::cout << count << " out of " << size << " distances incorrect "
			"(" << std::setprecision(2) << (100.0*count) / size << std::setprecision(5) << "%) \n" <<
			"maximum error: " << max_error << "\n" <<
			"mean square error: " << variance << "\n" <<
			"standard deviation: " << standard_deviation << "\n";

	return standard_deviation > 0.0015;
}