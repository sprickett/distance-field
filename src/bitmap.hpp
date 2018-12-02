#pragma once
#include <tmap2d.hpp>
#include <fstream>
#include <vector>
#include <algorithm>
#include <array>
#include <iostream> //!!!!!!!!!!!


typedef int int32;
typedef unsigned char uint8;
typedef unsigned int uint32;
typedef unsigned short uint16;

enum class BitmapCompression : uint32
{
	BI_RGB,				// 	none 	Most common
	BI_RLE8,			//  RLE 8 - bit / pixel 			Can be used only with 8 - bit / pixel bitmaps
	BI_RLE4,			// 	RLE 4 - bit / pixel 			Can be used only with 4 - bit / pixel bitmaps
	BI_BITFIELDS,		//	OS22XBITMAPHEADER : Huffman 1D 	BITMAPV2INFOHEADER : RGB bit field masks,	BITMAPV3INFOHEADER + : RGBA
	BI_JPEG,			// 	OS22XBITMAPHEADER : RLE - 24 	BITMAPV4INFOHEADER + : JPEG image for printing[12]
	BI_PNG,				//									BITMAPV4INFOHEADER + : PNG image for printing[12 
	BI_ALPHABITFIELDS,	// 	RGBA bit field masks 			only Windows CE 5.0 with.NET 4.0 or later
	BI_CMYK,			//	none 							only Windows Metafile CMYK[3]
	BI_CMYKRLE8,		//	RLE - 8 						only Windows Metafile CMYK
	BI_CMYKRLE4,		// 	RLE - 4 						only Windows Metafile CMYK
};

struct BitmapFileHeader
{
	uint32 file_size = 0;	// 02  2 	4 bytes 	The size of the BMP file in bytes
	uint32 reserved = 0;	// 06  6 	4 bytes 	Reserved; actual value depends on the application that creates the image
	uint32 data_offset = 0; // 0A 10 	4 bytes 	The offset, i.e.starting address, of the byte where the bitmap image data(pixel array) can be found.
};

struct BitmapInfoHeader
{
	uint32 header_size = 40;		// 0E 14 4 	the size of this header(40 bytes)
	uint32 width = 0;				// 12 18 4 	the bitmap width in pixels(signed integer)
	uint32 height = 0;				// 16 22 4 	the bitmap height in pixels(signed integer)
	uint16 colour_planes = 1;		// 1A 26 2 	the number of color planes(must be 1)
	uint16 bits_per_pixel = 8;		// 1C 28 2 	the number of bits per pixel, which is the color depth of the image.Typical values are 1, 4, 8, 16, 24 and 32.
	BitmapCompression compression = BitmapCompression::BI_RGB;
	// 1E 30 4 	the compression method being used.See the next table for a list of possible values
	uint32 size = 0;				// 22 34 4 	the image size. This is the size of the raw bitmap data; a dummy 0 can be given for BI_RGB bitmaps.
	uint32 pixel_per_metre_x = 0;		// 26 38 4 	the horizontal resolution of the image. (pixel per metre, signed integer)
	uint32 pixel_per_metre_y = 0;		// 2A 42 4 	the vertical resolution of the image. (pixel per metre, signed integer)
	uint32 num_colours = 0;			// 2E 46 4 	the number of colors in the color palette, or 0 to default to 2n
	uint32 num_important = 0;		// 32 50 4 	the number of important colors used, or 0 when every color is important; generally ignored
};

const uint16 magic_bytes = *(uint16*)"BM";

bool load_bitmap(std::ifstream& f, std::vector<uint8>& data, int& width, int& height, int& stride, int& bits_per_pixel)
{
	uint16 magic = 0;
	BitmapFileHeader bf;
	BitmapInfoHeader bi;
	size_t pos = f.tellg();
	size_t data_pos = 0;
	uint32 i = 0;

	f.read((char *)&magic, sizeof(magic));
	if (magic != magic_bytes)
		goto fail_return;

	f.read((char *)&bf, sizeof(bf));
	f.read((char *)&bi, sizeof(bi));
	if (f.bad() || f.eof()
		|| bi.header_size != 40
		|| bi.colour_planes != 1
		|| bi.compression != BitmapCompression::BI_RGB
		|| bi.num_colours != 256
		)
		goto fail_return;

	width = bi.width;
	height = bi.height;
	bits_per_pixel = bi.bits_per_pixel;
	stride = (((width*bits_per_pixel + 31) / 8) & -4);

	if (bf.file_size - bf.data_offset < bi.size)
		goto fail_return;

	for (i = 0; i < bi.num_colours; ++i)
	{
		//unsigned int c;
		//f.read((char*)&c, sizeof(c));
		//printf("%3u %x\n", i, c);
	}

	data.clear();
	data.resize(bi.size);
	data_pos = pos + bf.data_offset;
	f.seekg(data_pos);
	f.read((char *)&data[0], bi.size);
	if (f.bad())
		goto fail_return;
	return true;
fail_return:
	f.seekg(pos);
	return false;
}
bool load_bitmap(const char * filename, std::vector<uint8>& data, int& width, int& height, int& stride, int& bits_per_pixel)
{
	std::ifstream f(filename, std::ios_base::binary);
	return load_bitmap(f, data, width, height, stride, bits_per_pixel);
}


bool saveBitmap(std::ofstream& f, const TMap<unsigned char>& image, const std::vector<unsigned>& palette = std::vector<unsigned>())
{
	if (image.width() * image.height() <= 0)
		return false;

	std::streampos pos = f.tellp(); // want to rewind under fail
	BitmapFileHeader bf;
	BitmapInfoHeader bi;
	uint32 n;
	   
	bi.width = image.width();
	bi.height = image.height();
	bi.bits_per_pixel = 8;
	bi.compression = BitmapCompression::BI_RGB;
	uint32_t stride = ((bi.width * bi.bits_per_pixel + 31) & -32) / 8u;
	bi.size = stride * bi.height;	
	bi.pixel_per_metre_x = 2835;
	bi.pixel_per_metre_y = 2835;
	size_t palette_size = palette.size();
	bool has_palette = palette_size > 1 && palette_size <= 256;
	bi.num_colours = has_palette? uint32_t(palette_size) : 256;
	bi.num_important = 0;
	bf.data_offset = sizeof(magic_bytes) + sizeof(bf) + sizeof(bi) + bi.num_colours * 4;
	bf.file_size = bf.data_offset + bi.size;
	   
	f.write((const char*)&magic_bytes, sizeof(magic_bytes));
	f.write((const char*)&bf, sizeof(bf));
	f.write((const char*)&bi, sizeof(bi));

	if (has_palette) 
	{
		for (auto n : palette)
			f.write((const char*)&n, sizeof(n));
	}
	else 
	{
		for (n = 0; n < 0x1000000; n += 0x10101)
			f.write((const char*)&n, sizeof(n));
	}

	if (image.stride() == stride)
	{
		f.write((const char*)image.ptr(), bi.size);
	}
	else 
	{
		uint32_t width_in_bytes = (image.width() * bi.bits_per_pixel + 7)/8u;
		uint32_t padding = stride - width_in_bytes;
		n = 0;
		for (unsigned y = 0; y < bi.height; ++y)
		{
			f.write((const char*)image.ptr(y), width_in_bytes);
			f.write((const char*)&n, padding);
		}
	}
	

	if (f.bad())
		goto fail_return;

	return true;
fail_return:
	f.seekp(pos);
	return false;
}

bool saveBitmap(std::ofstream& f, const TMap<std::array<uint8_t,3> >& image)
{
	if (image.width() * image.height() <= 0)
		return false;

	std::streampos pos = f.tellp(); // want to rewind under fail
	BitmapFileHeader bf;
	BitmapInfoHeader bi;
	uint32 n;

	bi.width = image.width();
	bi.height = image.height();
	bi.bits_per_pixel = 24;
	bi.compression = BitmapCompression::BI_RGB;
	uint32_t stride = ((bi.width * bi.bits_per_pixel + 31) & -32) / 8u;
	bi.size = stride * bi.height;
	bi.pixel_per_metre_x = 2835;
	bi.pixel_per_metre_y = 2835;
	bi.num_colours = 0;
	bi.num_important = 0;
	bf.data_offset = sizeof(magic_bytes) + sizeof(bf) + sizeof(bi) + bi.num_colours * 4;
	bf.file_size = bf.data_offset + bi.size;

	f.write((const char*)&magic_bytes, sizeof(magic_bytes));
	f.write((const char*)&bf, sizeof(bf));
	f.write((const char*)&bi, sizeof(bi));


	if (image.stride() == stride)
	{
		f.write((const char*)image.ptr(), bi.size);
	}
	else
	{
		uint32_t width_in_bytes = (image.width() * bi.bits_per_pixel + 7) / 8u;
		uint32_t padding = stride - width_in_bytes;
		n = 0;
		for (unsigned y = 0; y < bi.height; ++y)
		{
			f.write((const char*)image.ptr(y), width_in_bytes);
			f.write((const char*)&n, padding);
		}
	}


	if (f.bad())
		goto fail_return;

	return true;
fail_return:
	f.seekp(pos);
	return false;
}

bool saveBitmap(const char * filename, const TMap<std::array<uint8_t, 3> >& image)
{
	if (image.width() * image.height() <= 0)
		return false;
	std::ofstream fs(filename, std::ios_base::binary);
	return saveBitmap(fs, image);
}

bool saveBitmap(const char * filename, const TMap<unsigned char>& image, const std::vector<unsigned>& palette = std::vector<unsigned>())
{
	if (image.width() * image.height() <= 0)
		return false;
	std::ofstream fs(filename, std::ios_base::binary);
	return saveBitmap(fs, image, palette);
}
