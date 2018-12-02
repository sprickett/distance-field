/*
	Copyright (c) 2018, Shaun Prickett
	All rights reserved.

	Redistribution and use in source and binary forms, with or without
	modification, are permitted provided that the following conditions are met:

	* Redistributions of source code must retain the above copyright notice, this
	  list of conditions and the following disclaimer.

	* Redistributions in binary form must reproduce the above copyright notice,
	  this list of conditions and the following disclaimer in the documentation
	  and/or other materials provided with the distribution.

	* Neither the name of the copyright holder nor the names of its
	  contributors may be used to endorse or promote products derived from
	  this software without specific prior written permission.

	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
	AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
	IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
	DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
	FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
	DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
	SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
	CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
	OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
	OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/


#pragma once
#include "tmap2d.hpp"


namespace distance_field {

	void signedDistance(
			const TMap<int>& signed_square_distance, 
			TMap<float>& signed_distance);

	void signedDistance(
			const TMap<unsigned char>& binary,
			const TMap<std::pair<int16_t, int16_t> >& deltas, 
			TMap<float>& signed_distance);

	void deltaSweep(
			const TMap<unsigned char>& binary_input,
			TMap<std::pair<int16_t,int16_t> >& delta_output,
			double max_distance = std::numeric_limits<double>::infinity());

	void deltaSweepCached(
		const TMap<unsigned char>& binary_input,
		TMap< std::pair<int16_t, int16_t> >& delta_output,
		TMap< float >& distance_output,
		double max_distance = std::numeric_limits<double>::infinity());


	void dijkstra(
			const TMap<unsigned char>& binary_input, 
			TMap<int>& signed_square_distance_output,
			double max_distance = std::numeric_limits<double>::infinity());

	void simpleList(
		const TMap<unsigned char>& binary_input,
		TMap<int>& signed_square_distance_output);

#ifdef DISTANCE_FIELD_DEBUG
	const TMap<uint8_t>& getDebugImage(void);
	const TMap< std::pair<int16_t, int16_t> >& getDeltaField(void);
#endif
}