
#include "bitmap.hpp"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <string>

std::string toString(const Bitmap::Compression& error)
{
	const char * chars = "";
	switch (error)
	{
	case Bitmap::Compression::RGB:
		chars = "bitmap::Compression::RGB";
		break;
	case Bitmap::Compression::RLE8:
		chars = "bitmap::Compression::RLE8";
		break;
	case Bitmap::Compression::RLE4:
		chars = "bitmap::Compression::RLE4";
		break;
	case Bitmap::Compression::BITFIELDS:
		chars = "bitmap::Compression::BITFIELDS";
		break;
	case Bitmap::Compression::JPEG:
		chars = "bitmap::Compression::JPEG";
		break;
	case Bitmap::Compression::PNG:
		chars = "bitmap::Compression::PNG";
		break;
	case Bitmap::Compression::ALPHABITFIELDS:
		chars = "bitmap::Compression::ALPHABITFIELDS";
		break;
	case Bitmap::Compression::CMYK:
		chars = "bitmap::Compression::CMYK";
		break;
	case Bitmap::Compression::CMYKRLE8:
		chars = "bitmap::Compression::CMYKRLE8";
		break;
	case Bitmap::Compression::CMYKRLE4:
		chars = "bitmap::Compression::CMYKRLE4";
		break;
	}
	return chars;
}
std::string toString(const Bitmap::Error& error)
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

std::string toString(const Bitmap::LogicalColourSpace& lcs)
{
	const char * chars = "";
	switch (lcs)
	{
	case Bitmap::LogicalColourSpace::CALIBRATED_RGB:
		chars = "Bitmap::LogicalColourSpace::CALIBRATED_RGB";
		break;
	case Bitmap::LogicalColourSpace::sRGB:
		chars = "Bitmap::LogicalColourSpace::sRGB";
		break;
	case Bitmap::LogicalColourSpace::WINDOWS_COLOR_SPACE:
		chars = "Bitmap::LogicalColourSpace::WINDOWS_COLOR_SPACE";
		break;
	}
	return chars;
}


int main(int argc, char* argv[])
{
	if (argc < 2)
	{
		std::cout << "Missing path\n";
		return -1;
	}
	
	Bitmap bmp;
	TMap < uint8_t > binary;
	Bitmap::Error result = bmp.load(argv[1], binary);
	const Bitmap::Header &h = bmp.header();

	std::cout << "sizeof(Bitmap::Header): " << sizeof(Bitmap::Header) <<"\n";
	std::cout << "sizeof(Bitmap::InfoHeader): " << sizeof(Bitmap::InfoHeader) << "\n";
	std::cout << "sizeof(Bitmap::V3Header): " << sizeof(Bitmap::InfoHeader)+sizeof(Bitmap::BitFields) << "\n";
	std::cout << "sizeof(Bitmap::V4Header): " << sizeof(Bitmap::InfoHeader) + sizeof(Bitmap::BitFields) +sizeof(Bitmap::ColourSpace) << "\n";
	std::cout << "sizeof(Bitmap::V5Header): " << sizeof(Bitmap::InfoHeader) + sizeof(Bitmap::BitFields) + sizeof(Bitmap::ColourSpace) + sizeof(Bitmap::ColourProfile) << "\n\n";
	std::cout << "File Header\n" <<
			"  file size: " << h.file.file_size << "\n"
			"  file reserved: " << h.file.reserved <<"\n"
			"  data offset: " << h.file.data_offset <<"\n";
	if (h.info.header_size >= 12)
		std::cout << "Core Header\n"
				"  size of header: " << h.info.header_size <<"\n"
				"  width: " << h.info.width <<"\n"
				"  height: " << h.info.height <<"\n"
				"  colour planes: " << h.info.colour_planes <<"\n"
				"  bits per pixel: " << h.info.bits_per_pixel <<"\n";
	if (h.info.header_size >= 40)
		std::cout << "Info Header\n" <<
				"  compression: " << toString(h.info.compression) <<"\n"
				"  size: " << h.info.size <<"\n"
				"  pixel per metre x: " << h.info.pixel_per_metre_x <<"\n"
				"  pixel per metre y: " << h.info.pixel_per_metre_y <<"\n"
				"  num colours: " << h.info.num_colours <<"\n"
				"  num important: " << h.info.num_important << "\n";
	if (h.info.header_size >= 52)
		std::cout << "Bit Fields Header\n" <<
				"  red bits: " << std::hex << h.bf.red_bits <<"\n"
				"  green bits: " << h.bf.green_bits <<"\n"
				"  blue bits: " << h.bf.blue_bits <<"\n";
	if (h.info.header_size >= 56)
		std::cout << 
				"  alpha bits: " << h.bf.alpha_bits << std::dec <<	"\n";
	if (h.info.header_size >= 108)
		std::cout << "Colour Space Header\n" <<
				"  colour space: " << toString(h.cs.colour_space) << " (" << uint32_t(h.cs.colour_space) << ")\n"
				"  red x y z: " << h.cs.red_xyz.x << ", " << h.cs.red_xyz.y << ", " << h.cs.red_xyz.z << "\n"
				"  green x y z: " << h.cs.green_xyz.x << ", " << h.cs.green_xyz.y << ", " << h.cs.green_xyz.z << "\n"
				"  blue x y z: " << h.cs.blue_xyz.x << ", " << h.cs.blue_xyz.y << ", " << h.cs.blue_xyz.z << "\n"
				"  red gamma: " << h.cs.red_gamma << "\n"
				"  green gamma: " << h.cs.green_gamma << "\n"
				"  blue gamma: " << h.cs.blue_gamma << "\n";
	if (h.info.header_size >= 124)
		std::cout << "Colour Profile Header\n" <<
				"  colour space: " << toString(h.cp.colour_space) << " (" <<uint32_t(h.cp.colour_space) << ")\n"
				"  offset: " << h.cp.profile_offset <<"\n"
				"  size: " << h.cp.profile_size <<"\n"
				"  reserved: " << h.cp.reserved <<"\n";
	std::cout << "Palette\n" << std::hex;
	for (auto& c: bmp.palette())
		std::cout <<"  "<< c<<"\n";
	std::cout << std::dec;

	if (result != Bitmap::Error::NONE)
	{
		std::cout << "failed to load bitmap: "<< toString(result) << "\n";
		return -1;
	}
	return 0;
}