#include "distance_field.hpp"
#include "bitmap.hpp"
#include <sstream>
#include <iostream>
#include <iomanip>
#include <chrono>

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
	int width(void) const { return width_; }
	int height(void) const { return height_; }
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
	int width = binary.width();
	int height = binary.height();
	field.create(width,height);
	int mdy=0, mdx=0, md2 = std::numeric_limits<int>::max();
	for (int y = 0; y < height; ++y)
	{
		const uint8* B = binary.ptr(y);
		auto* D = field.ptr(y);
		if (y > 0)
		{
			mdx = D[-width].first;
			mdy = D[-width].second-1;
			mdy *= sign(mdy);
			mdx *= sign(mdx);
			++mdy;
			++mdx;
			md2 = mdx * mdx + mdy * mdy;
		}
		for (int x = 0; x < width; ++x)
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
		
			for (int yi = 0; yi<height;++yi)
			{
				int dy = yi - y;
				int dy2 = dy * dy;
				if (dy2 > md2) // row is bust
					continue;
				const uint8* Bi = binary.ptr(yi);
				for (int xi = 0; xi < width; ++xi)
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
			D[x].first =mdx;
			D[x].second = mdy;
		}

	}
	for (int y = 0; y < height; ++y)
	{
		auto* D = field.ptr(y);
		for (int x = 0; x < width; ++x)
		{
			mdx = D[x].first;
			mdy = D[x].second;
			D[x].first = mdx * 2 - sign(mdx);
			D[x].second = mdy * 2 - sign(mdy);

			//std::cout << D[x].x +x<< " " << D[x].y+y << "\n";
		}

	}
}


void draw_dot(TMap<uint8>& im, int x, int y, int)
{
	*im.ptr(x,y) ^= 255;
}
void draw_circle(TMap<uint8>& im, int x, int y, int radius)
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

void draw_square(TMap<uint8>& im, int x, int y, int side)
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
		if ((unsigned)y >= im.height())
			continue;
		uint8* B = im.ptr(y);
		for (i = std::max(0,xs), ie = std::min(im.width(),xe); i < ie; ++i)
			B[i] ^= 255;
		if (y == ym)
			s = -s;
	}
}

void test3()
{
	enum {
		W = 256,
		H = W,
		R = (W + H) / 8,
		X = W / 2,
		Y = H /2,
	};	

	TMap < std::pair<int16_t, int16_t> > deltas;
	TMap < std::pair<int16_t, int16_t> > check_deltas(W,H);
	TMap<uint8> binary(W,H);
	TMap<float> error(W,H);
	TMap<float> distance;
	TMap<int> sqr_distance;
	binary.setTo(0);

	draw_circle(binary, X, Y, R);

	deltaFieldBruteForce(binary, check_deltas);
	//distance_field::deltaSweep(binary, deltas);
	//distance_field::signedDistance(binary, deltas, distance);
	distance_field::dijkstra(binary, sqr_distance);
	distance_field::signedDistance(sqr_distance,distance);

	distance_field::signedDistance(binary, check_deltas, error);

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
			E[x] = D[x]-E[x];
	}
	float mx = 0;
	for (auto *e = error.ptr(), *E = e+W*H; e<E; ++e)
	{
		if (*e > 16000.f)
			*e = 0;
		mx = std::max(mx, *e);
	}
	std::cout << "max error: " << mx << "\n";
	
	float m = mx > 0 ? 255 / mx : 0;
	for (int y = 0; y < H; ++y)
	{
		float* E = error.ptr(y);
		uint8* B = binary.ptr(y);
		for (int x = 0; x < W; ++x)
			B[x] = E[x] * m;

	}
	save_bitmap("error.bmp", binary.ptr(), W*H, W,H,W, 8);




}

int rand4x4(void)
{
	unsigned r = rand();
	return (r & 3) + ((r >> 2) & 3) + ((r >> 4) & 3) + ((r >> 6) & 3) - 6;
}


int main(int argc, char* argv[])
{
	test3();
	return 0;
	
    std::vector<uint8> data;	
	int width, height, stride, bits_per_pixel;
	
	//const char* save_name = "df.bmp";
	//if (argc < 2)
	//	printf("Usage: %s <path to greyscale bitmap>\n", argv[0]);
	//else if (!load_bitmap(argv[1], data, width, height, stride, bits_per_pixel))
	//	printf("Failed to load \"%s\" bitmap\n", argv[1]);
	//else if (bits_per_pixel != 8)
	//	printf("Only greyscale bitmaps supported\n");
	//else if (!gen.generate(data.data(), data.data(), 128, width, height, stride, stride))
	//	printf("Failed to generate\n");
	//else if (!save_bitmap(save_name, data.data(), data.size(), width, height, stride, bits_per_pixel))
	//	printf("Failed to save \"%s\" bitmap\n", save_name);
	//else
	//	printf("Saved \"%s\" bitmap\n", save_name);

	return 0;
}