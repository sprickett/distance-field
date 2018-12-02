#include "distance_field.hpp"
#include "bitmap.hpp"

#include <algorithm>
#include <numeric>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <functional>
#include <random>
#include <map>


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

inline int lengthSquared(int x, int y)
{
	return x * x + y * y;
}

void deltaFieldBruteForce(const TMap<uint8>& binary, TMap < std::pair<int16_t, int16_t> >& field)
{
	int cols = binary.width();
	int rows = binary.height();
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


void drawDot(TMap<uint8>& im, int x, int y, int)
{
	*im.ptr(x,y) ^= 255;
}
void drawCircle(TMap<uint8>& im, int x, int y, int radius)
{
	int rows = im.height();
	int cols = im.width();

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

void drawSquare(TMap<uint8>& im, int x, int y, int side)
{
	int xs = std::max(0, x);
	int ys = std::max(0, y);
	int xe = std::min(im.width(), x + side);
	int ye = std::min(im.height(), y + side);
	for (int iy = ys; iy < ye; ++iy)
	{
		uint8* B = im.ptr(iy);
		for (int ix = xs; ix < xe; ++ix)
			B[ix] ^=  255;
	}

}

void drawDiamond(TMap<uint8 >& im, int x, int y, int side)
{
	
	int ys = y - side;	
	int ym = y;
	int ye = y + side + 1;
	int xs = x;
	int xe = x+1;
	int i, ie;
	int s = -1;

	for (y = ys; y<ye; ++y, xs+=s, xe-=s)
	{
		if (unsigned(y) >= unsigned(im.height()))
			continue;
		uint8* B = im.ptr(y);
		for (i = std::max(0,xs), ie = std::min(im.width(),xe); i < ie; ++i)
			B[i] ^= 255;
		if (y == ym)
			s = -s;
	}
}

void drawBitmap(TMap<uint8 >& im, const char * filename)
{
	int w, h, s, bpp;
	std::vector<unsigned char> buf;
	load_bitmap(filename, buf, w, h, s, bpp);
	if (bpp != 8)
		return;
	TMap<uint8> bmp(w,h,buf.data(),s );
	w = std::min(bmp.width(), im.width());
	h = std::min(bmp.height(), im.height());
	bmp(0, 0,w,h).copyTo(im(0,0,w,h));
}

void timeAllTest(const TMap<uint8 >& binary_image, int num_iterations)
{
	TMap<int> sqr_distance;
	double sz = binary_image.width()*binary_image.height();
	TMap < std::pair<int16_t, int16_t> > deltas;
	TMap < float > distance;


	using namespace std::chrono;
	for (int i = 0; i < num_iterations; ++i)
	{
		auto tp = Clock::now();
		distance_field::deltaSweep(binary_image, deltas);
		distance_field::signedDistance(binary_image,deltas, distance);
		auto du = Clock::now() - tp;
		std::cout << "delta sweep\t " <<
			duration_cast<microseconds>(du).count() << " us\t " <<
			int(sz/duration<double>(du).count()) << " pixels per second\n";
	}
	using namespace std::chrono;
	for (int i = 0; i < num_iterations; ++i)
	{
		auto tp = Clock::now();
		distance_field::deltaSweepCached(binary_image, deltas, distance);
		auto du = Clock::now() - tp;
		std::cout << "delta sweep cached\t " <<
			duration_cast<microseconds>(du).count() << " us\t " <<
			int(sz / duration<double>(du).count()) << " pixels per second\n";
	}
	for (int i = 0; i < num_iterations; ++i)
	{
		auto tp = Clock::now();
		distance_field::dijkstra(binary_image, sqr_distance);
		distance_field::signedDistance(sqr_distance, distance);
		auto du = Clock::now() - tp;
		std::cout << "dijkstra\t " <<
			duration_cast<microseconds>(du).count() << " us\t " <<
			int(sz/duration<double>(du).count()) << " pixels per second\n";
	}

	// bit too slow...
	//for (int i = 0; i < num_iterations; ++i)
	//{
	//	auto tp = Clock::now();
	//	distance_field::simpleList(binary_image, sqr_distance);
	//	auto du = Clock::now() - tp;
	//	std::cout << "simple list\n " <<
	//		duration_cast<microseconds>(du).count() << " us\n " <<
	//		int(sz/duration<double>(du).count()) << " pixels per second\n";
	//}
}

void saveDebugImage(const char * filename)
{
#ifdef DISTANCE_FIELD_DEBUG
	saveBitmap(filename, distance_field::getDebugImage());
#endif
}

std::vector< std::pair<int,int> > findIntegralDeltas(int sqr_distance)
{
	sqr_distance = abs(sqr_distance);
	int y = int(std::sqrt(sqr_distance));
	std::vector< std::pair<int, int> > deltas;
	for (int x = 0; x <= y; ++x)
	{
		int d = sqr_distance - x * x;
		while (y * y < d)
			++y;
		while (y * y > d)
			--y;
		if (y*y == d)
			deltas.push_back({ x,y });
	}
	return deltas;
}


float test(int algorithm_index = 0)
{
	enum {
		W =256,
		H = W,
		R = (W + H) / 7,
		X = W / 2,
		Y = H /2,
	};	

	TMap < std::pair<int16_t, int16_t> > deltas;
	TMap < std::pair<int16_t, int16_t> > deltas_check;
	TMap<uint8> binary(W,H);
	TMap< std::array<uint8,3 > > bgr(W, H);
	TMap<float> check(W,H);
	TMap<float> distance;
	TMap<int> sqr_distance;
	TMap<int> sqr_distance_check;
	binary.setTo(0);
	drawCircle(binary, X, Y, R);
	drawBitmap(binary, "felix256.bmp");
	//drawDiamond(binary, X, Y, R);

	switch (algorithm_index)
	{
	case 0:
		distance_field::dijkstra(binary, sqr_distance);
		distance_field::signedDistance(sqr_distance,distance);
#ifdef DISTANCE_FIELD_DEBUG
		deltas = distance_field::getDeltaField().clone();
		saveDebugImage("debug_dijkstra.bmp");
#endif
		break;
	case 1:
		distance_field::deltaSweep(binary, deltas);
		distance_field::signedDistance(binary, deltas, distance);
		sqr_distance.create(W, H);
		for (int y = 0; y < H; ++y)
		{
			int* O = sqr_distance.ptr(y);
			auto* I = deltas.ptr(y);
			for (int x = 0; x < W; ++x)
				O[x] = lengthSquared(I[x].first, I[x].second);
		}

#ifdef DISTANCE_FIELD_DEBUG
		deltas = distance_field::getDeltaField().clone();
#endif		
		break;
	case 2:
		distance_field::deltaSweepCached(binary, deltas, distance);
		distance_field::signedDistance(binary, deltas, distance);
		sqr_distance.create(W, H);
		for (int y = 0; y < H; ++y)
		{
			int* O = sqr_distance.ptr(y);
			auto* I = deltas.ptr(y);
			for (int x = 0; x < W; ++x)
				O[x] = lengthSquared(I[x].first, I[x].second);
		}
		break;
	default:
		return std::numeric_limits<float>::infinity();
	}


	distance_field::simpleList(binary, sqr_distance_check);
	distance_field::signedDistance(sqr_distance_check, check);
#ifdef DISTANCE_FIELD_DEBUG
	deltas_check = distance_field::getDeltaField().clone();
	saveDebugImage("debug_list.bmp");
#endif
	saveBitmap("image.bmp", binary);

	for (int y = 0; y < H; ++y)
	{
		float* D = distance.ptr(y);
		uint8* B = binary.ptr(y);
		for (int x = 0; x < W; ++x)
			B[x] = uint8_t(std::min(std::max(0.f, D[x] + 127.5f), 255.f));
	}
	saveBitmap("distance.bmp", binary);

	for (int y = 0; y < H; ++y)
	{
		float* D = check.ptr(y);
		uint8* B = binary.ptr(y);
		for (int x = 0; x < W; ++x)
			B[x] = uint8_t(std::min(std::max(0.f, D[x] + 127.5f), 255.f));
	}
	saveBitmap("distance_check.bmp", binary);

	float min_error = 0;// std::numeric_limits<float>::infinity();
	float max_error = 0;// -std::numeric_limits<float>::infinity();
	float error_ct = 0;
	float error = 0;
	for (int y = 0; y < H; ++y)
	{
		float* E = check.ptr(y);
		float* D = distance.ptr(y);
		for (int x = 0; x < W; ++x)
		{
			E[x] -= D[x];
			error_ct += E[x] != 0.f;
			error += E[x] * E[x];
			max_error = std::max(max_error, E[x]);
			min_error = std::min(min_error, E[x]);
		}
	}
	error /= W * H;
	std::cout << "error: " << error << "\n";
	std::cout << "min: " << min_error << "  max: " << max_error << "  count: " << error_ct  <<"\n";
	
	float d = max_error - min_error;
	float m = d > 0 ? 255 / d: 0;
	for (int y = 0; y < H; ++y)
	{
		float* E = check.ptr(y);
		uint8* B = binary.ptr(y);
		for (int x = 0; x < W; ++x)
			B[x] = uint8_t((E[x] - min_error) * m);
	}
	saveBitmap("error.bmp", binary);


#ifdef DISTANCE_FIELD_DEBUG
	if(error_ct < 64)
	{
		bgr.setTo({ 0,0,0 });
		for (int y = 0; y < H; ++y)
		{
			float* E = check.ptr(y);
			for (int x = 0; x < W; ++x)
			{
				if (!E[x])
					continue;
				auto d = *deltas_check.ptr(y,x);
				for (int yy = 0; yy < H; ++yy)
				{
					auto* D = deltas.ptr(yy);
					auto* C = deltas_check.ptr(yy);
					if (D == C)
						std::cout << "wtf!\n";
					auto* B = bgr.ptr(yy);
					for (int xx = 0; xx < W; ++xx)
					{
						
						if (D[xx] == d)
							B[xx][2] |= 255; 
						if (C[xx] == d)
							B[xx][1] |= 255;
					}
				}
			
			}
			saveBitmap("propagate.bmp", bgr);
		}
	}
	for (int y = 0; y < H; ++y)
	{
		auto* C = deltas_check.ptr(y);
		auto* D = deltas.ptr(y);
		uint8* I = binary.ptr(y);
		for (int x = 0; x < W; ++x)
		{
			uint8 val;
			if (C[x] == D[x])
				val = 0;
			else if (lengthSquared(D[x].first-x*2, D[x].second-y*2) != lengthSquared(C[x].first-x*2, C[x].second-y*2))
				val = 255;
			else
				val = 64;
			I[x] = val;
		}			
	}
	saveBitmap("difference.bmp", binary);
#endif

	return error;
}



template<typename T>
void dump(const TMap<T>& m ,int w =4)
{
	for (int y = 0; y < m.height(); ++y)
	{
		const T* p = m.ptr(y);
		for (int x = 0; x < m.width(); ++x)
			std::cout << std::setw(4) << int(p[x]) << " ";
		std::cout << "\n";
	}
	std::cout << "\n";
}


bool saveAsDistanceField(const char * filename)
{
	std::vector<uint8> data;
	int cols, rows, stride, bits_per_pixel;	
	


	if (!load_bitmap(filename, data, cols, rows, stride, bits_per_pixel))
	{
		std::cout << "Failed to load \"" << filename << "\" bitmap\n";
		return false;
	}

	if (bits_per_pixel != 8)
	{
		std::cout << "Only greyscale bitmaps supported\n";
		return false;
	}

	TMap<uint8> img(cols, rows, data.data(), stride);
	TMap<int> sdsq;
	distance_field::dijkstra(img, sdsq);
	std::transform(sdsq.begin(), sdsq.end(), img.begin(), [](const int& x)->uint8_t {
		float f = float(x);
		f = f < 0 ? -std::sqrtf(-f) : std::sqrtf(f);
		return uint8_t(std::min(std::max(0.f, 125.5f + f), 255.f));
	});

		
	std::string save_name;
	std::getline(std::istringstream(filename), save_name, '.');
	save_name += "_df.bmp";

	if (!saveBitmap(save_name.c_str(), img))
	{
		std::cout << "Failed to save \"" << save_name << "\" bitmap\n";
		return false;
	}

	std::cout << "Saved \"" << save_name << "\" bitmap\n";

	return true;
}




int main(int argc, char* argv[])
{
	if (argc < 2)
	{	
		//test(2);
		//return 0;
		constexpr auto SZ = 512;
		TMap<uint8> binary(SZ,SZ);
		binary.setTo(0);
		drawCircle(binary, SZ, SZ, (SZ * 2) / 3);
		timeAllTest(binary, 20);
		return 0;
	}		
	else
	{
		saveAsDistanceField(argv[1]);
	}
	return 0;
}