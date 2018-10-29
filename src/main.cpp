#include "distance_field.hpp"
#include "bitmap.hpp"
#include <algorithm>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <functional>

typedef std::chrono::high_resolution_clock Clock;

class BouncingBallImage
{
public:
	struct Ball
	{
		int x, y, vx, vy, r;
	};
	BouncingBallImage(void)
		:width_(1024)
		,height_(1024)
		,data_(width_*height_)
		,balls_({
			{123,701,-7,5,1},
			{456,78,6,-5,31},
			{319,604,-9,7,127},	
			{319,604,8,9,270}
			})
	{
		update();
	}
	int cols(void) const { return width_; }
	int rows(void) const { return height_; }
	const uint8* data(void) const { return data_.data(); }
	void update(void)
	{
		std::fill(data_.begin(), data_.end(), 0);

		int x, y, y2, r2;
		for (auto& b : balls_)
		{
			b.x += b.vx;
			b.y += b.vy;
			x = width_ - b.r;
			y = height_ - b.r;
			r2 = b.r*b.r;


			if (b.x < b.r)
			{
				b.x = 2 * b.r - b.x;
				b.vx = std::abs(b.vx);
			}
			if( b.x > x)
			{
				b.x = 2 * x- b.x;
				b.vx = -std::abs(b.vx);
			}

			
			if (b.y < b.r)
			{
				b.y = 2 * b.r - b.y;
				b.vy = std::abs(b.vy);
			}
			if (b.y > y)
			{
				b.y = 2 * y - b.y;
				b.vy = -std::abs(b.vy);
			}

			for (y = -b.r; y <= b.r; ++y)
			{
				y2 = r2 - y * y;
				uint8* p = data_.data() + (b.y + y) * width_ + b.x;
				for (x = -b.r; x <= b.r; ++x)
					if (x * x < y2)
						p[x] ^= 0xff;
			}
		}
	}
private:
	int width_;
	int height_;
	std::vector<uint8> data_;
	std::vector<Ball> balls_;
};

inline int sign(const int& x)
{
	return (0 < x) - (x < 0);
}

void deltaFieldBruteForce(const TMap<uint8>& binary, TMap < std::pair<int16_t, int16_t> >& field)
{
	int cols = binary.cols();
	int rows = binary.rows();
	field.create(cols,rows);
	int mdy=0, mdx=0, md2 = std::numeric_limits<int>::max();
	for (int y = 0; y < rows; ++y)
	{
		const uint8* B = binary.ptr(y);
		auto* D = field.ptr(y);
		if (y > 0)
		{
			mdx = D[-cols].first;
			mdy = D[-cols].second-1;
			mdy *= sign(mdy);
			mdx *= sign(mdx);
			++mdy;
			++mdx;
			md2 = mdx * mdx + mdy * mdy;
		}
		for (int x = 0; x < cols; ++x)
		{
			int b = B[x];
			if (x > 0)
			{
				mdx = D[x - 1].first-1;
				mdy = D[x - 1].second;
				mdy *= sign(mdy);
				mdx *= sign(mdx);
				++mdy;
				++mdx;
				md2 = mdx * mdx + mdy * mdy ;
			}
		
			for (int yi = 0; yi<rows;++yi)
			{
				int dy = yi - y;
				int dy2 = dy * dy;
				if (dy2 > md2) // row is bust
					continue;
				const uint8* Bi = binary.ptr(yi);
				for (int xi = 0; xi < cols; ++xi)
				{
					if (b == Bi[xi])
						continue;

					int dx = xi - x;
					int d2 = dy2 + dx * dx;
					if (d2 < md2)
					{
						mdx = dx;
						mdy = dy;
						md2 = d2;
					}
				}
			}
			D[x].first = mdx;
			D[x].second = mdy;
		}

	}
	for (int y = 0; y < rows; ++y)
	{
		auto* D = field.ptr(y);
		for (int x = 0; x < cols; ++x)
		{
			mdx = D[x].first;
			mdy = D[x].second;
			D[x].first = mdx * 2 - sign(mdx);
			D[x].second = mdy * 2 - sign(mdy);
		}

	}
}


void draw_dot(TMap<uint8>& im, int x, int y, int)
{
	*im.ptr(x,y) ^= 255;
}
void draw_circle(TMap<uint8>& im, int x, int y, int radius)
{
	int rows = im.rows();
	int cols = im.cols();

	int r2 = radius * radius;
	int xs = std::max(0, x-radius);
	int ys = std::max(0, y-radius);
	int xe = std::min(cols, x + radius);
	int ye = std::min(rows, y + radius);
	int dy, dx;
	for (int iy = ys; iy < ye; ++iy)
	{
		dy = y - iy;
		dy *= dy;
		uint8* B = im.ptr(iy);
		for (int ix = xs; ix < xe; ++ix)
		{
			dx = x - ix;
			if(dx*dx+dy<r2)
				B[ix] ^= 255;
		}
	}

}

void draw_square(TMap<uint8>& im, int x, int y, int side)
{
	int xs = std::max(0, x);
	int ys = std::max(0, y);
	int xe = std::min(im.cols(), x + side);
	int ye = std::min(im.rows(), y + side);
	for (int iy = ys; iy < ye; ++iy)
	{
		uint8* B = im.ptr(iy);
		for (int ix = xs; ix < xe; ++ix)
			B[ix] ^=  255;
	}

}

void draw_diamond(TMap<uint8 >& im, int x, int y, int side)
{
	
	int ys = y - side;	
	int ym = y;
	int ye = y + side;
	int xs = x;
	int xe = x+1;
	int i, ie;
	int s = -1;

	for (y = ys; y<ye; ++y, xs+=s, xe-=s)
	{
		if ((unsigned)y >= im.rows())
			continue;
		uint8* B = im.ptr(y);
		for (i = std::max(0,xs), ie = std::min(im.cols(),xe); i < ie; ++i)
			B[i] ^= 255;
		if (y == ym)
			s = -s;
	}
}

void timeAllTest(const TMap<uint8 >& binary_image, int num_iterations)
{
	TMap<int> sqr_distance;
	TMap < std::pair<int16_t, int16_t> > offsets;
	using namespace std::chrono;
	for (int i = 0; i < num_iterations; ++i)
	{
		auto tp = Clock::now();
		distance_field::dijkstra(binary_image, sqr_distance);
		auto du = Clock::now() - tp;
		std::cout << "dijkstra " << duration_cast<microseconds>(du).count() << "\n";
	}
	for (int i = 0; i < num_iterations; ++i)
	{
		auto tp = Clock::now();
		distance_field::simpleList(binary_image, sqr_distance);
		auto du = Clock::now() - tp;
		std::cout << "simple list " << duration_cast<microseconds>(du).count() << "\n";
	}
	for (int i = 0; i < num_iterations; ++i)
	{
		auto tp = Clock::now();
		distance_field::deltaSweep(binary_image, offsets);
		auto du = Clock::now() - tp;
		std::cout << "delta sweep " << duration_cast<microseconds>(du).count() << "\n";
	}
}

int test(int algorithm_index = 0)
{
	enum {
		W = 256,
		H = W,
		R = (W + H) / 8,
		X = W / 2,
		Y = H /2,
	};	

	TMap < std::pair<int16_t, int16_t> > deltas;
	TMap < std::pair<int16_t, int16_t> > check_deltas;// (W, H);
	TMap<uint8> binary(W,H);
	TMap<float> error(W,H);
	TMap<float> distance;
	TMap<int> sqr_distance;
	binary.setTo(0);

	draw_circle(binary, X, Y, R);
	
	switch (algorithm_index)
	{
	case 0:
		distance_field::dijkstra(binary, sqr_distance);
		distance_field::signedDistance(sqr_distance,distance);
		break;
	case 1:
		distance_field::deltaSweep(binary, deltas);
		distance_field::signedDistance(binary, deltas, distance);
		break;
	case 2:
		distance_field::simpleList(binary, sqr_distance);
		distance_field::signedDistance(sqr_distance, distance);
		break;
	default:
		return -1;
	}

	distance_field::simpleList(binary, sqr_distance);
	distance_field::signedDistance(sqr_distance, error);

	save_bitmap("image.bmp", binary.ptr(), W*H, W, H, W, 8);
	for (int y = 0; y < H; ++y)
	{
		float* D = error.ptr(y);
		uint8* B = binary.ptr(y);
		for (int x = 0; x < W; ++x)
			B[x] = std::min(std::max(0.f, D[x] + 127.5f), 255.f);
	}
	save_bitmap("check.bmp", binary.ptr(), W*H, W, H, W, 8);
	for (int y = 0; y < H; ++y)
	{
		float* D = distance.ptr(y);
		uint8* B = binary.ptr(y);
		for (int x = 0; x < W; ++x)
			B[x] = std::min(std::max(0.f, D[x] + 127.5f), 255.f);
	}
	save_bitmap("distance.bmp", binary.ptr(), W*H, W, H, W, 8);

	for (int y = 0; y < H; ++y)
	{
		float* E = error.ptr(y);
		float* D = distance.ptr(y);
		for (int x = 0; x < W; ++x)
			E[x] -= D[x];
	}
	float min_error = std::numeric_limits<float>::infinity();
	float max_error = -std::numeric_limits<float>::infinity();
	for (int y = 0; y < H; ++y)
	{
		float* E = error.ptr(y);
		for (int x = 0; x < W; ++x)
		{
			max_error = std::max(max_error, E[x]);
			min_error = std::min(min_error, E[x]);
		}
	}
	std::cout << "error range: " << min_error << " - " << max_error <<"\n";
	
	float d = max_error - min_error;
	float m = m > 0 ? 255 / d: 0;
	for (int y = 0; y < H; ++y)
	{
		float* E = error.ptr(y);
		uint8* B = binary.ptr(y);
		for (int x = 0; x < W; ++x)
			B[x] = (E[x] - min_error) * m;
	}
	save_bitmap("error.bmp", binary.ptr(), W*H, W,H,W, 8);

	int ret_val = max_error ? int(std::min(max_error+1.f, float(std::numeric_limits<int>::max()))) : 0;
	return ret_val;
}



int main(int argc, char* argv[])
{
    std::vector<uint8> data;	
	int cols, rows, stride, bits_per_pixel;
	
	const char* save_name = "df.bmp";
	if (argc < 2)
	{
		return test(0);
		TMap<uint8> binary(256,256);
		binary.setTo(0);
		draw_circle(binary, 512, 512, 256);
		timeAllTest(binary, 20);
		return 0;
	}		
	else if (!load_bitmap(argv[1], data, cols, rows, stride, bits_per_pixel))
		std::cout << "Failed to load \"" << argv[1] << "\" bitmap\n";
	else if (bits_per_pixel != 8)
		std::cout << "Only greyscale bitmaps supported\n";
	else
	{
		TMap<uint8> img(cols, rows, data.data(), stride);
		TMap<int> signed_distance_squared;
		distance_field::dijkstra(img, signed_distance_squared);
		for (int y = 0; y < rows; ++y)
		{
			auto* B = img.ptr(y);
			auto* D = signed_distance_squared.ptr(y);
			for (int x = 0; x < cols; ++x)
				B[x] = std::min( std::max( 0.f, 125.5f + (D[x] < 0 ? -std::sqrtf(-D[x]) : std::sqrtf(D[x]))), 255.f);
		}
		if (!save_bitmap(save_name, data.data(), data.size(), cols, rows, stride, bits_per_pixel))
			std::cout << "Failed to save \"" << save_name << "\" bitmap\n";
		else
			std::cout << "Saved \"" << save_name << "\" bitmap\n";
	}

	return 0;
}