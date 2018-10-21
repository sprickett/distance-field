#pragma once
#include "tmap2d.hpp"
#include <vector>
#include <memory>
#include <algorithm>
#include <forward_list>
#include <tuple>

namespace distance_field {

	void signedDistance(const TMap<int>& signed_square_distance, TMap<float>& signed_distance);
	void signedDistance(const TMap<unsigned char>& binary, const TMap<std::pair<int16_t, int16_t> >& deltas, TMap<float>& signed_distance);

	void deltaSweep(const TMap<unsigned char>& binary_input,
		TMap<std::pair<int16_t,int16_t> >& delta_output,
		double max_distance = std::numeric_limits<double>::infinity());

	void dijkstra(const TMap<unsigned char>& binary_input, 
			TMap<int>& signed_square_distance_output,
			double max_distance = std::numeric_limits<double>::infinity());



	class DRA
	{
	public:


		bool operator()(
			const unsigned char *binary_input,
			int width,
			int height,
			int input_stride);

		struct Delta
		{
			typedef short type;
			type x, y;
		};

		static constexpr Delta infinity_delta = { 0, std::numeric_limits<Delta::type>::min() }; // delta x must be zero

		int width(void) const { return width_; }
		int height(void) const { return height_; }
		int stride(void) const { return stride_; }
		const Delta *data(void)  const { return xy0_; }
	private:

		bool init(int width, int height);
		bool mark(const unsigned char *binary_input, int input_stride);
		void generate(void);

		int width_;
		int height_;
		int stride_;
		std::vector<Delta> delta_;
		std::vector<float> distance_;
		Delta * xy0_;
		float * d0_;
	};
	

}