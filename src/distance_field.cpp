#include "distance_field.hpp"

#include <algorithm>
#include <tuple>
#include <vector>
#include <array>
#include <map>
//#include <queue> // priority_queue
#include <iostream>


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
inline int length_squared(int x, int y)
{
	return x * x + y * y;
}
void setSign(const TMap<unsigned char>& binary_input, TMap<int>& l2_buffer)
{
	int x, y;
	int width = l2_buffer.width();
	int height = l2_buffer.height();
	const unsigned char* I;
	int * D;
	for (y = 0; y < height; ++y)
	{
		D = l2_buffer.ptr(y);
		I = binary_input.ptr(y);
		for (x = 0; x < width; ++x)
			if (I[x])
				D[x] = -D[x];
	}
}
#ifdef DISTANCE_FIELD_DEBUG
static TMap<uint8_t> debug_image = TMap<uint8_t>();
static TMap< std::pair<int16_t, int16_t> > debug_delta = TMap< std::pair<int16_t, int16_t> >();
#endif
}

class DijkstraDetail
{
public:
	static void dijkstra(const TMap<unsigned char>& binary_input, TMap<int>& signed_square_distance_output, double max_distance)
	{
		Queue queue;
		initQueue(queue);
		
		init(signed_square_distance_output, binary_input.width(), binary_input.height(), max_distance);
		seed(binary_input, signed_square_distance_output, queue);

		generate(signed_square_distance_output, queue);
		setSign(binary_input, signed_square_distance_output);
	}
private:
	typedef std::tuple<int, int, int*> Node; // x, y, *l2
	typedef std::pair<const int, std::vector< Node > > Front;
	typedef std::map<int, std::vector< Node > > Queue;

#ifdef DISTANCE_FIELD_DEBUG
	static int *pl2_origin;
	static int pl2_stride;
#endif

	static void pop_front(std::map<int, std::vector< Node > > &queue)
	{
		queue.erase(queue.begin());
	}
	static Front& front(std::map<int, std::vector< Node > > &queue)
	{
		return *queue.begin();
	}

	static void initQueue(std::map<int, std::vector< Node > > &queue)
	{
		queue[std::numeric_limits<int>::max()];
	}

	static void insert(int x, int y, int *pl2, std::map<int, std::vector< Node > >& queue)
	{
		int l2 = length_squared(x, y);
		if (*pl2 <= l2)
			return;
		*pl2 = l2;
		queue[l2].push_back(std::make_tuple(x, y, pl2));
#ifdef DISTANCE_FIELD_DEBUG
		int od = pl2 - pl2_origin;
		int ox = od % pl2_stride;
		int oy = od / pl2_stride;
		*debug_delta.ptr(oy,ox) = { ox*2 + x, oy*2+y };
		//*debug_delta.ptr(oy, ox) = { x, y };
#endif
	}


	static void init(TMap<int>& l2_buffer, int width, int height, double max_distance)
	{
		l2_buffer.create(width + 2, height + 2);
		l2_buffer.setTo(0); // just for now...
		// set distance at borders to zero to avoid any distance check
		//l2_buffer(0, 0, l2_buffer.cols(), 1).setTo(0);
		//l2_buffer(0, l2_buffer.rows() - 1, l2_buffer.cols(), 1).setTo(0);
		//l2_buffer(l2_buffer.cols() - 1, 0, 1, l2_buffer.rows()).setTo(0);

		double d2 = max_distance * max_distance;
		int max_l2 = d2 < double(std::numeric_limits<int>::max()) ? d2 : std::numeric_limits<int>::max();
		l2_buffer = l2_buffer(1, 1, width, height);
		l2_buffer.setTo(max_l2);
#ifdef DISTANCE_FIELD_DEBUG
		pl2_origin = l2_buffer.ptr();
		pl2_stride = l2_buffer.stride();
		debug_delta.create(width,height);
		debug_delta.setTo({ 0, 0 });
#endif
	}

	static bool seed(const TMap<unsigned char>& binary_input, TMap<int>& l2_buffer, Queue& queue)
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
					insert(+1, +1, D + x, queue);
					insert(-1, -1, D + x + stride + 1, queue);
				}
				if (I[x + 1] != I[x + input_stride])
				{
					insert(-1, +1, D + x + 1, queue);
					insert(+1, -1, D + x + stride, queue);
				}
				// horizontal 
				if (I[x] != I[x + 1])
				{
					insert(+1, 0, D + x, queue);
					insert(-1, 0, D + x + 1, queue);
				}
				// vertical
				if (I[x] != I[x + input_stride])
				{
					insert(0, +1, D + x, queue);
					insert(0, -1, D + x + stride, queue);
				}

			}
			// last column
			if (I[x] != I[x + input_stride])
			{
				insert(0, +1, D + x, queue);
				insert(0, -1, D + x + stride, queue);
			}

		}
		// last row
		D = l2_buffer.ptr(y);
		I = binary_input.ptr(y);
		for (x = 0; x < width - 1; ++x)
		{
			if (I[x] != I[x + 1])
			{
				insert(+1, 0, D + x, queue);
				insert(-1, 0, D + x + 1, queue);
			}
		}
#ifdef DISTANCE_FIELD_DEBUG
		debugSeeds(binary_input, l2_buffer, queue);
#endif
		return true;
	}

#ifdef DISTANCE_FIELD_DEBUG
	static void debugSeeds(const TMap<unsigned char>& binary_input, TMap<int>& l2_buffer, Queue& queue)
	{
		constexpr auto M1 = 2;
		constexpr auto M2 = M1*2;
		const int cols = binary_input.width();
		const int rows = binary_input.height();
		debug_image.create( cols * M2, rows * M2 );
		debug_image.setTo(128);// = debug_image(0, 0, cols * M2 + 1, rows * M2 + 1);
		

		for (int y = 0; y < rows; ++y)
		{
			auto I = binary_input.ptr(y);
			for (int x = 0; x < cols; ++x)
				debug_image(x*M2, y*M2, M2, M2).setTo(I[x] ? 160 :96);
		}
		auto l2_0 = l2_buffer.ptr();
		unsigned stride = l2_buffer.stride();
		for (auto& v : queue)
		{
			for (auto& t : v.second)
			{
				int* l2 = std::get<2>(t);
				if (v.first < *l2)
					continue;
				int x = l2 - l2_0;
				int y = x / stride;
				x %= stride;
				x = x * M2 + M1;
				y = y * M2 + M1;
				//*debug_image.ptr(y, x) = 0;
				*debug_image.ptr(y + std::get<1>(t)*M1, x + std::get<0>(t)*M1) = 255;
			}
		}
	}
	static void debugQueue(const Queue& queue)
	{
		
		int sum = 0, mx = 0;
		for (auto& p : queue)
		{
			int sz = int(p.second.size());
			mx = std::max(mx, sz);
			sum += sz;
		}
		std::cout << "Dijkstra Queue: size: " << queue.size()
			<< "  mean: " << sum / queue.size()
			<< "  max: " << mx
			<< "\n";
	}
#endif

	static void generate(TMap<int>& l2_buffer, Queue& queue)
	{
		int x, y;
		int * pl2;
		int l2;
		int stride = l2_buffer.stride();
		std::vector< Node > vec;
		while (!queue.empty())
		{
#ifdef DISTANCE_FIELD_DEBUG
			debugQueue(queue);
#endif
			{
				auto& p = front(queue);
				l2 = p.first;
				std::swap(vec, p.second);
			}
			pop_front(queue);
			for (auto& f : vec)
			{
				x = std::get<0>(f);
				y = std::get<1>(f);
				pl2 = std::get<2>(f);
				if (*pl2 < l2)
					continue; // something better came up

				int sx = sign(x);
				int sy = sign(y);

				insert(x + 2 * sx, y + 2 * sy, pl2 - sy * stride - sx, queue);
				insert(x, y + 2 * sy, pl2 - sy * stride, queue);
				insert(x + 2 * sx, y, pl2 - sx, queue);
			}
		}
	}


};

#ifdef DISTANCE_FIELD_DEBUG
int *DijkstraDetail::pl2_origin=nullptr;
int DijkstraDetail::pl2_stride=0;
#endif

class DeltaSweepDetail
{
public:
	typedef std::pair<int16_t, int16_t> delta_type;
	static void deltaSweep(const TMap<unsigned char>& binary_input, TMap< delta_type >& delta_output, double max_distance)
	{
		constexpr size_t KSize = 1;
		if (!init(delta_output, binary_input.width(), binary_input.height(), KSize))
			return;
		seed(binary_input, delta_output);
		generate(delta_output, makeKernel<KSize>(delta_output.stride()));
#ifdef DISTANCE_FIELD_DEBUG
		debug_delta = delta_output.clone();
		for (int y = 0; y < debug_delta.height(); ++y)
		{
			auto* D = debug_delta.ptr(y);
			for (int x = 0; x < debug_delta.width(); ++x)
				D[x].first += x*2, D[x].second += y*2;
		}
#endif
	}

private:
	typedef std::pair< std::pair<int, int>, std::pair<int, int> > kernel_type;
	

	static bool init(TMap<delta_type >& deltas, int width, int height, int hsz)
	{
		constexpr int max_dim = std::numeric_limits<int16_t>::max();
		if ((width | height) < 0 || width > max_dim || height > max_dim)
			return false;
		deltas.create(width + hsz*2, height + hsz*2);
		deltas.setTo({ 0, std::numeric_limits<int16_t>::min() });
		deltas = deltas(hsz,hsz, width, height);
		return true;
	}

	static TMap<delta_type> transpose(const TMap<delta_type>& image)
	{
		const int stride = image.stride();
		const int w = image.height(); // handed!
		const int h = image.width();
		TMap<delta_type> transposed(w+2,h+2);
		transposed.setTo({ 0, std::numeric_limits<int16_t>::min() });
		transposed = transposed(1, 1, w, h);
		for (int y = 0; y < h; ++y)
		{
			auto *O = transposed.ptr(y);
			auto *I = image.ptr(0, y);
			for (int x = 0; x < w; ++x)
				O[x] = I[x*stride];
		}
		return transposed;
	}

	static void seed(const TMap<uint8_t>&binary, TMap<std::pair<int16_t, int16_t> >& deltas)
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

	static void deltaCheck(unsigned& distance, std::pair<int16_t, int16_t>& delta,const std::pair<int16_t, int16_t>& neighbour, int x, int y)
	{
		int dx = neighbour.first + x;
		int dy = neighbour.second + y;
		int dist = length_squared(dx, dy);
		if (dist < distance)
		{
			distance = dist;
			delta.first = dx;
			delta.second = dy;
		}
	}

	template <size_t KSz>
	static void generate(TMap<std::pair<int16_t, int16_t> >& deltas, const std::array<kernel_type, KSz>& K)
	{
		std::pair<int16_t, int16_t> *D;
		int x, y; // position
		unsigned distance; // square distance 	
		int width = deltas.width();
		int height = deltas.height();
	
#ifdef DISTANCE_FIELD_DEBUG
		for (auto& k : K)
		{
			std::cout << k.first.first << ", " << k.first.second << " : " << k.second.first << ", " << k.second.second << "\n";
		}
#endif

		// forward pass
		for (y = 0; y < height; ++y)
		{
			D = deltas.ptr(y);
			// ++ pass
			for (x = 0; x < width; ++x)
			{
				distance = length_squared(D[x].first, D[x].second);
				for (auto& k : K) // we like to think looping this fixed size array gets unrolled
					deltaCheck(distance, D[x], D[x + k.second.first + k.second.second], k.first.first, k.first.second);
			}

			// -+ pass
			for (x = width - 1; x >= 0; --x)
			{
				distance = length_squared(D[x].first, D[x].second);
				auto& k = K.back();
				deltaCheck(distance, D[x], D[x - k.second.first + k.second.second], -k.first.first, k.first.second);
			}
		}

		// exit early if image is empty 
		if (D[width - 1].second == std::numeric_limits<int16_t>::min())
			return;

		// backward pass
		for (y = height - 1; y >= 0; --y)
		{
			D = deltas.ptr(y);
			// -- pass
			for (x = width - 1; x >= 0; --x)
			{
				distance = length_squared(D[x].first, D[x].second);
				for (auto& k : K) 
					deltaCheck(distance, D[x], D[x - k.second.first - k.second.second], -k.first.first, -k.first.second);
			}
			// +- pass
			for (x = 0; x < width; ++x)
			{
				distance = length_squared(D[x].first, D[x].second);
				auto& k = K.back();
				deltaCheck(distance, D[x], D[x + k.second.first - k.second.second], k.first.first, -k.first.second);
			}
		}


	}
	 


	template <size_t HSz>
	static auto makeKernel(int stride)
	{
		constexpr size_t KSz = ((HSz * 2 + 1)*(HSz * 2 + 1) - 1) / 2;
		std::array< kernel_type, KSz > K;
		using std::make_pair;
		size_t i = 0;
		for (int y = -int(HSz); y < 0; ++y)
			for (int x = -int(HSz); x <= int(HSz); ++x, ++i)
				K[i] = make_pair(make_pair(x * 2, y * 2), make_pair(x,y*stride));
		for (int x = -int(HSz); x < 0; ++x, ++i)
			K[i] = make_pair(make_pair(x * 2, 0), make_pair(x, 0));
		return K;
	}

	template <size_t KSz>
	static void showK(const std::array< kernel_type, KSz>& K)
	{
		for (auto& k : K)
		{
			std::cout << k.first.first << ", " << k.first.second << " : " << k.second << "\n";
		}
		std::cout << "\n";
	}

};


class DeltaSweepDetailCached
{
public:
	typedef std::pair<int16_t, int16_t> Delta;


	static void deltaSweep(
			const TMap<unsigned char>& binary_input, 
			TMap< Delta >& delta_output, 
			TMap< float >& distance_output, 
			double max_distance)
	{
		constexpr size_t KSize = 3;
		if (!init(distance_output, delta_output, binary_input.width(), binary_input.height(), KSize))
			return;
		seed(binary_input, distance_output, delta_output);
		generate(distance_output, delta_output, makeKernel<KSize>(delta_output.stride()));

#ifdef DISTANCE_FIELD_DEBUG
		debug_delta = delta_output.clone();
		for (int y = 0; y < debug_delta.height(); ++y)
		{
			auto* D = debug_delta.ptr(y);
			for (int x = 0; x < debug_delta.width(); ++x)
				D[x].first += x * 2, D[x].second += y * 2;
		}
#endif
	}

private:
	struct Neighbour
	{
		Neighbour(void) 
		{};
		Neighbour(int x, int y, int offset, float estimate)
			:x(x), y(y), offset(offset), estimate(estimate)
		{}
		int x;
		int y;
		int offset;
		float estimate;
	};

	static constexpr float root2 = 1.4142135624f;

	static bool init(TMap<float >& distances, TMap< Delta >& deltas, int width, int height, int hsz)
	{
		constexpr int max_dim = std::numeric_limits<int16_t>::max();
		if ((width | height) < 0 || width > max_dim || height > max_dim)
			return false;

		deltas.create(width, height);
		distances.create(width + hsz * 2, height + hsz * 2);
		distances.setTo(std::numeric_limits<float>::infinity());
		distances = distances(hsz, hsz, width, height);
		
		return true;
	}

	static void seed(const TMap<uint8_t>&binary, TMap<float >& distances, TMap< Delta >& deltas)
	{
		const unsigned char* I;
		Delta *V;
		float *D;
		int x, y;
		int width = binary.width();
		int height = binary.height();
		int binary_stride = binary.stride();
		int delta_stride = deltas.stride();
		int distances_stride = distances.stride();

		// mark immediate borders
		for (y = 0; y < height - 1; ++y)
		{
			I = binary.ptr(y);
			V = deltas.ptr(y);
			D = distances.ptr(y);
			for (x = 0; x < width - 1; ++x)
			{
				// diagonals 
				if (I[x] != I[x + binary_stride + 1])
				{
					D[x] = root2;
					V[x] = { 1,1 };
					D[x + distances_stride + 1] = root2;
					V[x + delta_stride + 1] = { -1,-1 };
				}

				if (I[x + 1] != I[x + binary_stride])
				{
					D[x + 1] = root2;
					V[x + 1] = { -1,1 };
					D[x + distances_stride] = root2;
					V[x + delta_stride] = { 1,-1 };
				}
				// horizontal 
				if (I[x] != I[x + 1])
				{
					D[x] = 1.f;
					V[x] = { 1,0 };
					D[x + 1] = 1.f;
					V[x + 1] = { -1,0 };
				}
				// vertical
				if (I[x] != I[x + binary_stride])
				{
					D[x] = 1.f;
					V[x] = { 0,1 };
					D[x + distances_stride] = 1.f;
					V[x + delta_stride] = { 0,-1 };
				}
			}
			// last column
			if (I[x] != I[x + binary_stride])
			{
				D[x] = 1.f;
				V[x] = { 0,1 };
				D[x + distances_stride] = 1.f;
				V[x + delta_stride] = { 0,-1 };
			}

		}
		// last row
		I = binary.ptr(y);
		V = deltas.ptr(y);
		D = distances.ptr(y);
		for (x = 0; x < width - 1; ++x)
		{
			if (I[x] != I[x + 1])
			{
				D[x] = 1.f;
				V[x] = { 1,0 }; 
				D[x + 1] = 1.f;
				V[x + 1] = { -1,0 };
			}
		}
	}

	//static void deltaCheck(unsigned& distance, std::pair<int16_t, int16_t>& delta, float neighbour, int x, int y)
	//{
	//	int dist = length_squared(x, y);
	//	if (dist < distance)
	//	{
	//		distance = dist;
	//		delta.first = dx;
	//		delta.second = dy;
	//	}
	//}

	template <size_t KSz>
	static void generate(TMap<float >& distances, TMap< Delta >& deltas, const std::array<Neighbour, KSz>& K)
	{
		Delta *V;
		float *D;
		int x, y, i; // position
		float distance, estimate; // distance 
		int dx, dy; // nearest border offset
		int width = distances.width();
		int height = distances.height();
		

		// forward pass
		for (y = 0; y < height; ++y)
		{
			V = deltas.ptr(y);
			D = distances.ptr(y);
			// ++ pass
			for (x = 0; x < width; ++x)
			{
				distance = D[x];
				for (auto& k : K) // we like to think looping this fixed size array gets unrolled
				{
					i = x + k.offset;
					estimate = D[i] + k.estimate;
					if (estimate < distance)
					{
						V[x].first = V[i].first + k.x;
						V[x].second = V[i].second + k.y;
						D[x] = distance = std::sqrt(length_squared(V[x].first, V[x].second));
					}
				}					
			}

			// -+ pass
			for (x = width - 1; x >= 0; --x)
			{
				distance = D[x];
				const Neighbour& k = K.back();
				i = x - k.offset; // ! negated
				estimate = D[i] + k.estimate;
				if (estimate < distance)
				{
					V[x].first = V[i].first - k.x; // ! negated
					V[x].second = V[i].second;
					D[x] = distance = std::sqrt(length_squared(V[x].first, V[x].second));
				}
			}
		}

		// TODO
		//// exit early if image is empty 
		//if (D[width - 1].second == std::numeric_limits<int16_t>::min())
		//	return;

		// backward pass
		for (y = height - 1; y >= 0; --y)
		{
			V = deltas.ptr(y);
			D = distances.ptr(y);
			// -- pass
			for (x = width - 1; x >= 0; --x)
			{
				distance = D[x];
				for (auto& k : K)
				{
					i = x - k.offset;
					estimate = D[i] + k.estimate;
					if (estimate < distance)
					{
						V[x].first = V[i].first - k.x;
						V[x].second = V[i].second - k.y;
						D[x] = distance = std::sqrt(length_squared(V[x].first, V[x].second));
					}
				}
			}
			// +- pass
			for (x = 0; x < width; ++x)
			{
				distance = D[x];
				const Neighbour& k = K.back();
				i = x + k.offset; // ! negated
				estimate = D[i] + k.estimate;
				if (estimate < distance)
				{
					V[x].first = V[i].first + k.x; // ! negated
					V[x].second = V[i].second;
					D[x] = distance = std::sqrt(length_squared(V[x].first, V[x].second));
				}
			}
		}


	}



	template <size_t HSz>
	static auto makeKernel(int stride)
	{
		constexpr size_t KSz = ((HSz * 2 + 1)*(HSz * 2 + 1) - 1) / 2;
		std::array< Neighbour, KSz > K;
		using std::make_pair;
		size_t i = 0;
		for (int y = -int(HSz); y < 0; ++y)
			for (int x = -int(HSz); x <= int(HSz); ++x, ++i)
				K[i] = Neighbour(x * 2, y * 2, x + y * stride, sqrtf(length_squared(x * 2, y * 2)));
		for (int x = -int(HSz); x < 0; ++x, ++i)
			K[i] = Neighbour(x * 2, 0, x, sqrtf(length_squared(x * 2, 0)));
		return K;
	}

	template <size_t KSz>
	static void showK(const std::array< Neighbour, KSz>& K)
	{
		for (auto& k : K)
		{
			std::cout << k.first.first << ", " << k.first.second << " : " << k.second << "\n";
		}
		std::cout << "\n";
	}
};


class SimpleListDetail
{
private:
	typedef std::pair<int, int> point;

	static void init(TMap<int>& l2_buffer, int width, int height)
	{
		l2_buffer.create(width , height);
		l2_buffer.setTo(std::numeric_limits<int>::max());
	}

	static void seed(const TMap<uint8_t>&binary, std::vector<point>& points)
	{
		const unsigned char* I;
		int x, y;
		int width = binary.width();
		int height = binary.height();
		int stride = binary.stride();

		// mark immediate borders
		for (y = 0; y < height - 1; ++y)
		{
			I = binary.ptr(y);
			for (x = 0; x < width - 1; ++x)
			{
				// diagonals 
				if (I[x] != I[x + stride + 1] || I[x + 1] != I[x + stride])
					points.push_back(point(x*2+1,y*2+1));
				// horizontal 
				if (I[x] != I[x + 1])
					points.push_back(point(x * 2 + 1, y * 2));
				// vertical
				if (I[x] != I[x + stride])
					points.push_back(point(x * 2, y * 2 + 1));
			}
			// last vertical
			if (I[x] != I[x + stride])
				points.push_back(point(x * 2, y * 2 + 1));

		}
		// last row
		I = binary.ptr(y);
		for (x = 0; x < width - 1; ++x)
			if (I[x] != I[x + 1])
				points.push_back(point(x * 2 + 1, y * 2));
#ifdef DISTANCE_FIELD_DEBUG
		debugSeeds(binary, points);
#endif
	}

#ifdef DISTANCE_FIELD_DEBUG
	static void debugSeeds(const TMap<uint8_t>&binary, std::vector<point>& points)
	{
		constexpr auto M1 = 2;
		constexpr auto M2 = M1 * 2;
		const int cols = binary.width();
		const int rows = binary.height();
		debug_image.create(cols * M2, rows * M2);
		debug_image.setTo(128);
		
		for (int y = 0; y < rows; ++y)
		{
			auto I = binary.ptr(y);
			for (int x = 0; x < cols; ++x)
				debug_image(x*M2, y*M2, M2, M2).setTo(I[x] ? 160 : 96);
		}

		for (auto& p : points)
			*debug_image.ptr(p.second*M1 + M1, p.first*M1+M1) = 255;	
	}
#endif

	static void generate(TMap<int>& l2_buffer, const std::vector<point>& points)
	{
		int* D;
		int x, y; // position
		int d; // square distance 
		int dx, dy; // postion * 2
		int width = l2_buffer.width();
		int height = l2_buffer.height();
	
		// forward pass
		for (y = 0; y < height; ++y)
		{
			dy = y * 2;
			D= l2_buffer.ptr(y);
			for (x = 0; x < width; ++x)
			{
				dx = x * 2;
				d = std::numeric_limits<int>::max();
#ifdef DISTANCE_FIELD_DEBUG
				for (auto& p : points)
				{
					
					int l2 = length_squared(p.first - dx, p.second - dy);
					if (l2 < d)
					{
						d = l2;
						*debug_delta.ptr(y, x) = { p.first , p.second };
						//*debug_delta.ptr(y, x) = { p.first - dx, p.second - dy };
					}
				}
#else	
				for (auto& p : points)
					d = std::min(d, length_squared(p.first - dx, p.second - dy));
#endif		
				D[x] = d;
			}
		}

	}

public:
	static void simpleList(const TMap<unsigned char>& binary_in, TMap< int >& l2_out)
	{
		std::vector<point> points;
		init(l2_out, binary_in.width(), binary_in.height());
#ifdef DISTANCE_FIELD_DEBUG
		debug_delta.create(binary_in.width(), binary_in.height());
		debug_delta.setTo({0, 0});
#endif
		seed(binary_in, points);
		generate(l2_out,points);// , makeKernel<1>(delta_output.stride()));
		setSign(binary_in, l2_out);
	}
};

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
			F[x] = D[x] < 0.f ? -std::sqrt(-D[x]) : std::sqrt(D[x]);
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
			float d = std::sqrtf(length_squared(D[x].first, D[x].second));
			F[x] = B[x]? -d : d;
		}
	}
}

void distance_field::deltaSweep(const TMap<unsigned char>& binary_input, TMap< std::pair<int16_t, int16_t> >& delta_output, double max_distance)
{
	DeltaSweepDetail::deltaSweep(binary_input, delta_output, max_distance);
}

void distance_field::deltaSweepCached(
	const TMap<unsigned char>& binary_input,
	TMap< DeltaSweepDetailCached::Delta >& delta_output,
	TMap< float >& distance_output,
	double max_distance)
{
	DeltaSweepDetailCached::deltaSweep(binary_input, delta_output, distance_output, max_distance);
}

void distance_field::dijkstra(const TMap<unsigned char>& binary_input, TMap<int>& signed_square_distance_output, double max_distance)
{
	DijkstraDetail::dijkstra(binary_input, signed_square_distance_output, max_distance);
}

void distance_field::simpleList(const TMap<unsigned char>& binary_input, 
		TMap<int>& signed_square_distance_output)
{
	SimpleListDetail::simpleList(binary_input, signed_square_distance_output);
}


#ifdef DISTANCE_FIELD_DEBUG
const TMap<uint8_t>& distance_field::getDebugImage(void)
{
	return debug_image;
}
const TMap< std::pair<int16_t, int16_t> >&  distance_field::getDeltaField(void)
{
	return debug_delta;
}
#endif