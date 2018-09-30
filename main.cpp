#include "distance_field.hpp"
#include "bitmap.hpp"

int main(int argc, char* argv[])
{
	DistanceFieldGenerator gen;
    std::vector<uint8> data;	
	int width, height, stride, bits_per_pixel;
	
	const char* save_name = "df.bmp";
	if (argc < 2)
		printf("Usage: %s <path to greyscale bitmap>\n", argv[0]);
	else if (!load_bitmap(argv[1], data, width, height, stride, bits_per_pixel))
		printf("Failed to load \"%s\" bitmap\n", argv[1]);
	else if (bits_per_pixel != 8)
		printf("Only greyscale bitmaps supported\n");
	else if (!gen(data.data(), data.data(), 128, width, height, stride, stride))
		printf("Failed to generate\n");
	else if (!save_bitmap(save_name, data, width, height, stride, bits_per_pixel))
		printf("Failed to save \"%s\" bitmap\n", save_name);
	else
		printf("Saved \"%s\" bitmap\n", save_name);

	return 0;
}