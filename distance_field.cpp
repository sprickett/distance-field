#include "distance_field.hpp"

#include <algorithm>
#include <iostream>

using namespace distance_field;

namespace {
	template<typename T>
	T clamp(T value, T minimum, T maximum)
	{
		return std::min(std::max(value, minimum), maximum);
	}
	inline int sign(const int& x)
	{
		return (0 < x) - (x < 0);
	}
	inline int l2Norm(int x, int y)
	{
		return x * x + y * y;
	}

	typedef std::tuple<int, int, int*> dijkstra_node;
	typedef std::pair<int, std::vector< dijkstra_node > > dijkstra_front;
	typedef std::forward_list <  dijkstra_front  > dijkstra_queue;

	void dijkstraInsert(int x, int y, int *pl2, dijkstra_queue& queue)
	{
		int l2 = l2Norm(x, y);
		if (*pl2 <= l2)
			return;
		*pl2 = l2;
		auto i = queue.before_begin();
		auto j = i++;
		while (l2 > i->first) // not checking for end, see init
			j = i++;
		if (l2 < i->first) // insert new l2
			i = queue.emplace_after(j, dijkstra_front(l2, std::vector<dijkstra_node>()));
		i->second.push_back(std::make_tuple(x, y, pl2));
	}

	void dijkstraInit(TMap<int>& l2_buffer, int width, int height, double max_distance)
	{
		l2_buffer.create(width + 2, height + 2);

		// set distance at borders to zero 
		int i;
		int stride = l2_buffer.stride();
		size_t sz = stride * l2_buffer.height();
		auto* D = l2_buffer.ptr();
		for (i = 0; i < stride; ++i)
			D[i] = 0;
		for (i = stride; i < sz; i += stride)
			D[i] = 0;
		for (i = stride * (l2_buffer.height() - 1); i < sz; ++i)
			D[i] = 0;

		double d2 = max_distance * max_distance;
		int max_l2 = d2 < double(std::numeric_limits<int>::max()) ? d2 : std::numeric_limits<int>::max();
		l2_buffer = l2_buffer(1, 1, width, height);
		l2_buffer.setTo(max_l2);
	}

	bool dijkstraSeed(const TMap<unsigned char>& binary_input, TMap<int>& l2_buffer, dijkstra_queue& queue)
	{
		if (!binary_input.ptr())
			return false;

		int x, y;
		int width = l2_buffer.width();
		int height = l2_buffer.height();
		int stride = l2_buffer.stride();
		int input_stride = binary_input.stride();
		const unsigned char* I;
		int * D;

		// mark immediate borders
		for (y = 0; y < height - 1; ++y)
		{
			D = l2_buffer.ptr(y);
			I = binary_input.ptr(y);
			for (x = 0; x < width - 1; ++x)
			{
				// diagonals 
				if (I[x] != I[x + input_stride + 1])
				{
					dijkstraInsert(+1, +1, D + x, queue);
					dijkstraInsert(-1, -1, D + x + stride + 1, queue);
				}
				if (I[x + 1] != I[x + input_stride])
				{
					dijkstraInsert(-1, +1, D + x + 1, queue);
					dijkstraInsert(+1, -1, D + x + stride, queue);
				}
				// horizontal 
				if (I[x] != I[x + 1])
				{
					dijkstraInsert(+1, 0, D + x, queue);
					dijkstraInsert(-1, 0, D + x + 1, queue);
				}
				// vertical
				if (I[x] != I[x + input_stride])
				{
					dijkstraInsert(0, +1, D + x, queue);
					dijkstraInsert(0, -1, D + x + stride, queue);
				}

			}
			// last column
			if (I[x] != I[x + input_stride])
			{
				dijkstraInsert(0, +1, D + x, queue);
				dijkstraInsert(0, -1, D + x + stride, queue);
			}

		}
		// last row
		D = l2_buffer.ptr(y);
		I = binary_input.ptr(y);
		for (x = 0; x < width - 1; ++x)
		{
			if (I[x] != I[x + 1])
			{
				dijkstraInsert(1, 0, D + x, queue);
				dijkstraInsert(-1, 0, D + x + 1, queue);
			}
		}

		return true;
	}

	void dijkstraGenerate(TMap<int>& l2_buffer, dijkstra_queue& queue)
	{
		int x, y;
		int * pl2;
		int l2;
		int stride = l2_buffer.stride();
		std::vector< dijkstra_node > vec;
		while (!queue.empty())
		{
			auto& p = queue.front();
			l2 = p.first;
			std::swap(vec, p.second);
			queue.pop_front(); // p is now undefined!
			for (auto& f : vec)
			{
				x = std::get<0>(f);
				y = std::get<1>(f);
				pl2 = std::get<2>(f);
				if (*pl2 < l2)
					continue; // something better came up
				if (x >= 0) // only try insert if increasing l2
					dijkstraInsert(x + 2, y, pl2 - 1, queue);
				if (y >= 0)
					dijkstraInsert(x, y + 2, pl2 - stride, queue);
				if (x <= 0)
					dijkstraInsert(x - 2, y, pl2 + 1, queue);
				if (y <= 0)
					dijkstraInsert(x, y - 2, pl2 + stride, queue);
			}
		}
	}

	void dijkstraSetSign(const TMap<unsigned char>& binary_input, TMap<int>& l2_buffer)
	{
		int x, y;
		int width = l2_buffer.width();
		int height = l2_buffer.height();
		const unsigned char* I;
		int * D;
		for (y = 0; y < height - 1; ++y)
		{
			D = l2_buffer.ptr(y);
			I = binary_input.ptr(y);
			for (x = 0; x < width - 1; ++x)
				if (I[x])
					D[x] = -D[x];
		}
	}

	bool deltaSweepInit(TMap<std::pair<int16_t, int16_t> >& deltas, int width, int height)
	{
		constexpr int max_dim = std::numeric_limits<int16_t>::max();
		if ((width | height) < 0 || width > max_dim || height > max_dim)
			return false;
		deltas.create(width + 2, height + 2);
		deltas.setTo(std::pair<int16_t, int16_t>(0, std::numeric_limits<int16_t>::min()));
		deltas = deltas(1, 1, width, height);
		return true;
	}

	void deltaSweepSeed(const TMap<uint8_t>&binary, TMap<std::pair<int16_t, int16_t> >& deltas)
	{
		const unsigned char* I;
		std::pair<int16_t, int16_t> *D;
		int x, y;
		int width = binary.width();
		int height = binary.height();
		int binary_stride = binary.stride();
		int delta_stride = deltas.stride();

		// mark immediate borders
		for (y = 0; y < height - 1; ++y)
		{
			I = binary.ptr(y);
			D = deltas.ptr(y);
			for (x = 0; x < width - 1; ++x)
			{
				// diagonals 
				if (I[x] != I[x + binary_stride + 1])
					D[x] = { 1,1 }, D[x + delta_stride + 1] = { -1,-1 };
				if (I[x + 1] != I[x + binary_stride])
					D[x + 1] = { -1,1 }, D[x + delta_stride] = { 1,-1 };
				// horizontal 
				if (I[x] != I[x + 1])
					D[x] = { 1,0 }, D[x + 1] = { -1,0 };
				// vertical
				if (I[x] != I[x + binary_stride])
					D[x] = { 0,1 }, D[x + delta_stride] = { 0,-1 };
			}
			// last column
			if (I[x] != I[x + binary_stride])
				D[x] = { 0,1 }, D[x + delta_stride] = { 0,-1 };

		}
		// last row
		I = binary.ptr(y);
		D = deltas.ptr(y);
		for (x = 0; x < width - 1; ++x)
			if (I[x] != I[x + 1])
				D[x] = { 1,0 }, D[x + 1] = { -1,0 };
	}

	void deltaSweepGenerate(TMap<std::pair<int16_t, int16_t> >& deltas)
	{
		std::pair<int16_t, int16_t> *D;
		int x, y; // position
		unsigned d0, d1; // square distance 
		int dx, dy; // nearest border offset
		int width = deltas.width();
		int height = deltas.height();
		int stride = deltas.stride();


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
		//struct { int x, y, offset; } const K[] = {
		//{-1, -1, -stride_ - 1, },
		//{ 0, -1, -stride_, },
		//{ 1, -1, -stride_ + 1, },
		//{-1,  0, -1, },

		//{ 1,  1, stride_ + 1, },
		//{ 0,  1, stride_, },
		//{-1,  1, stride_ - 1, },
		//{ 1,  0, 1, },
		//};

		// forward pass
		for (y = 0; y < height; ++y)
		{
			D = deltas.ptr(y);
			for (x = 0; x < width; ++x)
			{
				dx = D[x].first;
				dy = D[x].second;
				d0 = dx * dx + dy * dy;

				for (int i = 0; i < 4; ++i)
				{
					dx = D[x + K[i].offset].first + K[i].x;
					dy = D[x + K[i].offset].second + K[i].y;
					d1 = dx * dx + dy * dy;
					if (d1 < d0)
					{
						d0 = d1;
						D[x].first = dx;
						D[x].second = dy;
					}
				}
			}
		}

		// exit early if image is empty 
		if (D[x - 1].second == std::numeric_limits<int16_t>::min())
			return;


		// backward pass
		for (y = height - 1; y >= 0; --y)
		{
			D = deltas.ptr(y);
			for (x = width - 1; x >= 0; --x)
			{
				dx = D[x].first;
				dy = D[x].second;
				d0 = dx * dx + dy * dy;
				for (int i = 4; i < 8; ++i)
				{
					dx = D[x + K[i].offset].first + K[i].x;
					dy = D[x + K[i].offset].second + K[i].y;
					d1 = dx * dx + dy * dy;
					if (d1 < d0)
					{
						d0 = d1;
						D[x].first = dx;
						D[x].second = dy;
					}
				}
			}
		}
	}
}


void distance_field::signedDistance(const TMap<int>& signed_square_distance, TMap<float>& signed_distance)
{
	int width = signed_square_distance.width();
	int height = signed_square_distance.height();
	signed_distance.create(width, height);
	for (int y = 0; y < height; ++y)
	{
		float* F = signed_distance.ptr(y);
		const int* D = signed_square_distance.ptr(y);
		for (int x = 0; x < width; ++x)
		{
			F[x] = D[x] < 0.f ? -std::sqrtf(-D[x]) : std::sqrtf(D[x]);
		}
	}
}

void distance_field::signedDistance(const TMap<unsigned char>& binary, const TMap<std::pair<int16_t, int16_t>>& deltas, TMap<float>& signed_distance)
{
	int width = binary.width();
	int height = binary.height();
	if (deltas.width() != width || deltas.height() != height)
		return;
	signed_distance.create(width, height);
	for (int y = 0; y < height; ++y)
	{
		float* F = signed_distance.ptr(y);
		auto* D = deltas.ptr(y);
		auto* B = binary.ptr(y);
		for (int x = 0; x < width; ++x)
		{
			float d = std::sqrtf(l2Norm(D[x].first, D[x].second));
			F[x] = B[x]? -d : d;
		}
	}
}

void distance_field::deltaSweep(const TMap<unsigned char>& binary_input, TMap< std::pair<int16_t, int16_t> >& delta_output, double max_distance)
{
	if (!deltaSweepInit(delta_output, binary_input.width(), binary_input.height()))
		return;
	deltaSweepSeed(binary_input, delta_output);
	deltaSweepGenerate(delta_output);
}

void distance_field::dijkstra(const TMap<unsigned char>& binary_input, TMap<int>& signed_square_distance_output, double max_distance)
{
	dijkstra_queue queue(1, dijkstra_front(std::numeric_limits<int>::max(), std::vector<dijkstra_node>()));
	dijkstraInit(signed_square_distance_output, binary_input.width(), binary_input.height(), max_distance);
	dijkstraSeed(binary_input, signed_square_distance_output, queue);
	dijkstraGenerate(signed_square_distance_output, queue);
	dijkstraSetSign(binary_input, signed_square_distance_output);
}

bool DRA::init(int width, int height)
{
	constexpr int max_dim = std::numeric_limits<Delta::type>::max();
	if ((width | height) < 0 || width > max_dim || height > max_dim)
		return false;
	width_ = width;
	height_ = height;
	stride_ = width_ + 1;
	size_t sz = stride_ * (height_ + 2) + 1;
	delta_.clear();
	delta_.resize(sz, infinity_delta);
	xy0_ = delta_.data() + stride_ + 1;
	distance_.clear();
	distance_.resize(sz, std::numeric_limits<float>::max());
	d0_ = distance_.data() + stride_ + 1;
	return true;
}

bool DRA::mark(const unsigned char * binary_input, int input_stride)
{
	if (!binary_input)
		return false;

	const unsigned char* I;
	Delta *D;
	float *F;
	int x, y;
	float f2 = sqrtf(2);
	// mark immediate borders
	for (y = 0; y < height_ - 1; ++y)
	{
		I = &binary_input[y*input_stride];
		D = &xy0_[y*stride_];
		F = &d0_[y*stride_];
		for (x = 0; x < width_ - 1; ++x)
		{
			// diagonals 
			if (I[x] != I[x + input_stride + 1])
				D[x] = { 1,1 }, D[x + stride_ + 1] = { -1,-1 }, F[x] = f2, F[x + 1] = f2;
			if (I[x + 1] != I[x + input_stride])
				D[x + 1] = { -1,1 }, D[x + stride_] = { 1,-1 }, F[x] = f2, F[x + 1] = f2;
			// horizontal 
			if (I[x] != I[x + 1])
				D[x] = { 1,0 }, D[x + 1] = { -1,0 }, F[x] = 1, F[x + 1] = 1;
			// vertical
			if (I[x] != I[x + input_stride])
				D[x] = { 0,1 }, D[x + stride_] = { 0,-1 }, F[x] = 1, F[x + 1] = 1;
		}
		// last column
		if (I[x] != I[x + input_stride])
			D[x] = { 0,1 }, D[x + stride_] = { 0,-1 }, F[x] = 1, F[x + 1] = 1;

	}
	// last row
	I = &binary_input[y*input_stride];
	D = &xy0_[y*stride_];
	F = &d0_[y*stride_];
	for (x = 0; x < width_ - 1; ++x)
		if (I[x] != I[x + 1])
			D[x] = { 1,0 }, D[x + 1] = { -1,0 }, F[x] = 1, F[x + 1] = 1;

	return true;
}

void DRA::generate(void)
{
	Delta *D;
	float *F;
	int x, y; // position
	unsigned d0, d1; // square distance 
	int dx, dy; // nearest border offset
	float f2 = sqrt(8);
	float f1 = 2;
	// look up table
	struct { float d;  int x, y, offset; } const K[] = {
		{f2,-2, -2, -stride_ - 1, },
		{f1, 0, -2, -stride_, },
		{f2, 2, -2, -stride_ + 1, },
		{f1,-2,  0, -1, },

		{f2, 2,  2, stride_ + 1, },
		{f1, 0,  2, stride_, },
		{f2,-2,  2, stride_ - 1, },
		{f1, 2,  0, 1, },
	};
	//struct { int x, y, offset; } const K[] = {
	//{-1, -1, -stride_ - 1, },
	//{ 0, -1, -stride_, },
	//{ 1, -1, -stride_ + 1, },
	//{-1,  0, -1, },

	//{ 1,  1, stride_ + 1, },
	//{ 0,  1, stride_, },
	//{-1,  1, stride_ - 1, },
	//{ 1,  0, 1, },
	//};

	// forward pass
	for (y = 0; y < height_; ++y)
	{
		D = xy0_ + y * stride_;
		F = d0_ + y * stride_;
		for (x = 0; x < width_; ++x)
		{
			for (int i = 0; i < 4; ++i)
			{
				if (F[x + K[i].offset] + K[i].d < F[x])
				{
					dx = D[x + K[i].offset].x + K[i].x;
					dy = D[x + K[i].offset].y + K[i].y;
					D[x].x = dx;
					D[x].y = dy;
					F[x] = sqrt(dx*dx + dy * dy);
				}
			}
		}
	}

	// exit early if image is empty 
	if (D[x - 1].y == std::numeric_limits<Delta::type>::min())
		return;


	// backward pass
	for (y = height_ - 1; y >= 0; --y)
	{
		D = xy0_ + y * stride_;
		F = d0_ + y * stride_;
		for (x = width_ - 1; x >= 0; --x)
		{
			for (int i = 4; i < 8; ++i)
			{
				float f = F[x];
				if (F[x + K[i].offset] + K[i].d < f)
				{
					dx = D[x + K[i].offset].x + K[i].x;
					dy = D[x + K[i].offset].y + K[i].y;
					f = sqrt(dx*dx + dy * dy);
				}
				if (f < F[x])
					F[x] = f, D[x].x = dx, D[x].y = dy;
			}
		}
	}
}

bool DRA::operator()(
	const unsigned char *binary_input,
	int width,
	int height,
	int input_stride)
{
	if (!init(width, height) || !mark(binary_input, input_stride))
		return false;
	generate();

	//generate();
	//generate();
	return true;
}