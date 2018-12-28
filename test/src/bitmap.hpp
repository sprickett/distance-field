#pragma once
#include <tmap2d.hpp>
#include <fstream>
#include <vector>
#include <algorithm>
#include <array>
#include <iostream> //!!!!!!!!!!!


class Bitmap
{
public:
	enum class Error : uint32_t
	{	
		NONE,
		FILE_IO_ERROR,
		UNSUPPORTED_HEADER_TYPE,
		BAD_MAGIC_BYTES,
		BAD_HEADER,
		BAD_PALETTE,
		BAD_BITS_PER_PIXEL,
		BAD_IMAGE,
	};


	enum class Compression : uint32_t
	{
		RGB,				// 	none 	Most common
		RLE8,			//  RLE 8 - bit / pixel 			Can be used only with 8 - bit / pixel bitmaps
		RLE4,			// 	RLE 4 - bit / pixel 			Can be used only with 4 - bit / pixel bitmaps
		BITFIELDS,		//	OS22XBITMAPHEADER : Huffman 1D 	BITMAPV2INFOHEADER : RGB bit field masks,	BITMAPV3INFOHEADER + : RGBA
		JPEG,			// 	OS22XBITMAPHEADER : RLE - 24 	BITMAPV4INFOHEADER + : JPEG image for printing[12]
		PNG,				//									BITMAPV4INFOHEADER + : PNG image for printing[12 
		ALPHABITFIELDS,	// 	RGBA bit field masks 			only Windows CE 5.0 with.NET 4.0 or later
		CMYK,			//	none 							only Windows Metafile CMYK[3]
		CMYKRLE8,		//	RLE - 8 						only Windows Metafile CMYK
		CMYKRLE4,		// 	RLE - 4 						only Windows Metafile CMYK
	};

	enum class HeaderSize : uint32_t
	{
		CORE_HEADER = 12,
		INFO_HEADER = 40,
		V2_HEADER = 52,
		V3_HEADER = 56,
		V4_HEADER = 108,
		V5_HEADER = 124, 
	};

	enum class LogicalColourSpace : uint32_t
	{
		CALIBRATED_RGB = 0x00000000,
		sRGB = 0x73524742, // "sRGB"
		WINDOWS_COLOR_SPACE = 0x57696E20 // "Win "
	};



	struct FileHeader
	{
		uint32_t file_size = 0;	// 02  2 	4 bytes 	The size of the BMP file in bytes
		uint32_t reserved = 0;	// 06  6 	4 bytes 	Reserved; actual value depends on the application that creates the image
		uint32_t data_offset = 0; // 0A 10 	4 bytes 	The offset, i.e.starting address, of the byte where the bitmap image data(pixel array) can be found.	
	};
	struct InfoHeader
	{
		uint32_t header_size = 0;// 108;		// 0E 14 4 	the size of this header(40 bytes)
		uint32_t width = 0;				// 12 18 4 	the bitmap width in pixels(signed integer)
		uint32_t height = 0;				// 16 22 4 	the bitmap height in pixels(signed integer)
		uint16_t colour_planes = 1;		// 1A 26 2 	the number of color planes(must be 1)
		uint16_t bits_per_pixel = 0;		// 1C 28 2 	the number of bits per pixel, which is the color depth of the image.Typical values are 1, 4, 8, 16, 24 and 32
		Compression compression = Compression::RGB;	// 1E 30 4 	the compression method being used.See the next table for a list of possible values
		uint32_t size = 0;				// 22 34 4 	the image size. This is the size of the raw bitmap data; a dummy 0 can be given for BI_RGB bitmaps.
		uint32_t pixel_per_metre_x = 0;		// 26 38 4 	the horizontal resolution of the image. (pixel per metre, signed integer)
		uint32_t pixel_per_metre_y = 0;		// 2A 42 4 	the vertical resolution of the image. (pixel per metre, signed integer)
		uint32_t num_colours = 0;			// 2E 46 4 	the number of colors in the color palette, or 0 to default to 2n
		uint32_t num_important = 0;		// 32 50 4 	the number of important colors used, or 0 when every color is important; generally ignored
	};
	struct BitFields
	{
		uint32_t red_bits = 0;// 0xff0000;//red channel bit mask(valid because BI_BITFIELDS is specified)
		uint32_t green_bits = 0;// 0xff00;//green channel bit mask(valid because BI_BITFIELDS is specified)
		uint32_t blue_bits = 0;// 0xff;//blue channel bit mask(valid because BI_BITFIELDS is specified)
		uint32_t alpha_bits = 0;// 0xff000000;//alpha channel bit mask
	};
	struct XYZ
	{
		uint32_t x=0, y=0, z=0;
	};
	struct ColourSpace
	{
		LogicalColourSpace colour_space = LogicalColourSpace::CALIBRATED_RGB; //20 6E 69 57 
		XYZ red_xyz;
		XYZ green_xyz;
		XYZ blue_xyz;
		uint32_t red_gamma = 0; //Unused for LCS "Win " or "sRGB"
		uint32_t green_gamma = 0;
		uint32_t blue_gamma = 0;
	};
	struct ColourProfile
	{
		LogicalColourSpace colour_space;
		uint32_t profile_offset = 0; //Unused for LCS "Win " or "sRGB"
		uint32_t profile_size = 0;
		uint32_t reserved = 0;
	};

	struct Header
	{
		FileHeader file;
		InfoHeader info;
		BitFields bf;
		ColourSpace cs;
		ColourProfile cp;
	};

	static constexpr uint16_t magic_bytes = 'B' + 0x100*'M';





	//Error load(const char * filename, std::vector<uint8_t>& data, int& width, int& height, int& stride, int& bits_per_pixel)
	//{
	//	std::ifstream f(filename, std::ios_base::binary);
	//	return read(f, data, width, height, stride, bits_per_pixel);
	//}

	const Header& header(void) { return h_; }
	const std::vector<uint32_t>& palette(void){ return palette_; }
	const TMap<uint8_t>& image (void) { return data_; }
	void palette(const std::vector<uint32_t>& palette) 
	{
		palette_ = palette;
		is_valid_ &= validatePalette(palette_);
	}

	Error load(std::ifstream& fs , TMap<uint8_t>& image)
	{
		std::cout << "shh\n";
		return read(fs, image);
	}
	Error load(const char * filename, TMap<uint8_t>& image)
	{
		std::cout << "shhp\n";
		std::ifstream fs(filename, std::ios_base::binary);
		if ((fs.bad() | fs.eof() | fs.fail()))
			return  Error::FILE_IO_ERROR;

		return read(fs,image);
	}

	void clear(void)
	{

	}


	//bool saveBitmap(std::ofstream& f, const TMap<std::array<uint8_t, 3> >& image)
	//{
	//	if (image.width() * image.height() <= 0)
	//		return false;
	//
	//	std::streampos pos = f.tellp(); // want to rewind under fail
	//	BitmapHeader h;
	//	uint32_t n;
	//
	//	bi.width = image.width();
	//	bi.height = image.height();
	//	bi.bits_per_pixel = 24;
	//	bi.compression = BitmapCompression::BI_RGB;
	//	uint32_t stride = ((bi.width * bi.bits_per_pixel + 31) & -32) / 8u;
	//	bi.size = stride * bi.height;
	//	bi.pixel_per_metre_x = 2835;
	//	bi.pixel_per_metre_y = 2835;
	//	bi.num_colours = 0;
	//	bi.num_important = 0;
	//	bf.data_offset = sizeof(magic_bytes) + sizeof(bf) + sizeof(bi) + bi.num_colours * 4;
	//	bf.file_size = bf.data_offset + bi.size;
	//
	//	f.write((const char*)&magic_bytes, sizeof(magic_bytes));
	//	f.write((const char*)&bf, sizeof(bf));
	//	f.write((const char*)&bi, sizeof(bi));
	//
	//
	//	if (image.stride() == stride)
	//	{
	//		f.write((const char*)image.ptr(), bi.size);
	//	}
	//	else
	//	{
	//		uint32_t width_in_bytes = (image.width() * bi.bits_per_pixel + 7) / 8u;
	//		uint32_t padding = stride - width_in_bytes;
	//		n = 0;
	//		for (unsigned y = 0; y < bi.height; ++y)
	//		{
	//			f.write((const char*)image.ptr(y), width_in_bytes);
	//			f.write((const char*)&n, padding);
	//		}
	//	}
	//
	//
	//	if (f.bad())
	//		goto fail_return;
	//
	//	return true;
	//fail_return:
	//	f.seekp(pos);
	//	return false;
	//}

	//bool saveBitmap(const char * filename, const TMap<std::array<uint8_t, 3> >& image)
	//{
	//	if (image.width() * image.height() <= 0)
	//		return false;
	//	std::ofstream fs(filename, std::ios_base::binary);
	//	return save(fs, image);
	//}

	Error save(const char * filename, const TMap<unsigned char>& image, int bits_per_pixel = 8, const std::vector<unsigned>& palette = std::vector<unsigned>())
	{
		if (image.width() * image.height() <= 0)
			return Error::BAD_IMAGE;
		std::ofstream fs(filename, std::ios_base::binary);
		return write(fs, image, bits_per_pixel, palette);
	}

private:
	static bool validateBitsPerPixel(int bpp)
	{
		switch (bpp)
		{
		case 1:
		case 4:
		case 8:
		case 16:
		case 24:
		case 32:
			return true;
		default: // unsuported bits per pixel;
			return false;
		}
	}
	static bool validatePalette(const std::vector<uint32_t>& palette)
	{
		return palette.size() <= 256;
	}

	Error read(std::ifstream& f, TMap<uint8_t>& data)
	{
		uint16_t magic = 0;
		size_t pos = f.tellg();
		size_t chk;
		h_ = {};
		palette_.clear();
		data_ = data_(0, 0, 0, 0);
		
		f.read((char *)&magic, sizeof(magic));
		if (magic != magic_bytes)
			return cleanup(f, pos, Error::BAD_MAGIC_BYTES);

		f.read((char *)&h_.file, sizeof(h_.file) + 4);
		if ((f.bad() | f.eof() | f.fail()))
			return cleanup(f, pos, Error::FILE_IO_ERROR);				

		switch (HeaderSize(h_.info.header_size))
		{
		case HeaderSize::CORE_HEADER: 
			// 16 bit width and height
			f.read((char *)&h_.info.width, 2);
			f.seekg(std::ios::cur, 2);
			f.read((char *)&h_.info.height, 2);
			f.seekg(std::ios::cur, 2);
			f.read((char *)&h_.info.colour_planes, 4);
			break;
		case HeaderSize::INFO_HEADER:
			f.read((char *)&h_.info.width, 36); 
			break;
		case HeaderSize::V2_HEADER:
			f.read((char *)&h_.info.width, 48);
			break;
		case HeaderSize::V3_HEADER:
			f.read((char *)&h_.info.width, 52);
			break;
		case HeaderSize::V4_HEADER:
			f.read((char *)&h_.info.width, 36);
			if (h_.info.compression == Bitmap::Compression::BITFIELDS ||
					h_.info.compression == Bitmap::Compression::ALPHABITFIELDS)
			{
				f.read(((char *)&h_.bf), sizeof(h_.bf));
				f.read(((char *)&h_.cs), sizeof(h_.cs));
			}
			else
			{
				f.read(((char *)&h_.cs), sizeof(h_.cs));
				f.read(((char *)&h_.bf), sizeof(h_.bf));
			}
			break;
		case HeaderSize::V5_HEADER:
			f.read((char *)&h_.info.width, 120);
			break;
		default: // unrecognised header size
			return cleanup(f, pos, Error::UNSUPPORTED_HEADER_TYPE);
		}
		if ((f.bad() | f.fail())) // allow eof
			return cleanup(f, pos, Error::FILE_IO_ERROR);

		unsigned stride = (h_.info.width * h_.info.bits_per_pixel + 31) / 32;

		if (!validateBitsPerPixel(h_.info.bits_per_pixel) ||
				h_.info.colour_planes != 1 ||
				h_.info.num_colours > 256 ||
				//h_.v4.size != h_.v4.height * stride ||
				h_.file.file_size < h_.file.data_offset + h_.info.size)
			return cleanup(f, pos, Error::BAD_HEADER);

		palette_.resize(h_.info.num_colours);
		if (!palette_.empty())
			f.read((char*)palette_.data(), palette_.size() * sizeof(palette_.front()));

		if (h_.info.size != h_.info.height * stride)
			data.create(h_.info.size, 1);
		else
			data.create(stride, h_.info.height);
		f.seekg(pos + h_.file.data_offset);
		f.read((char *)data.ptr(), h_.info.size);
	
		if ((f.bad() | f.eof() | f.fail()))
			return cleanup(f, pos, Error::FILE_IO_ERROR);

		return Error::NONE;
	}


	Error write(std::ofstream& f, const TMap<unsigned char>& image, int bits_per_pixel = 8, const std::vector<unsigned>& palette = std::vector<unsigned>())
	{
		if (image.width() * image.height() <= 0)
			return Error::BAD_IMAGE;


		std::streampos pos = f.tellp(); // want to rewind under fail
		uint32_t n;

		h_.info.header_size = sizeof(h_.info);
		h_.info.width = image.width();
		h_.info.height = image.height();
		h_.info.bits_per_pixel = 8;
		h_.info.compression = Compression::RGB;
		uint32_t stride = ((h_.info.width *h_.info.bits_per_pixel + 31) & -32) / 8u;
		h_.info.size = stride * h_.info.height;
		h_.info.pixel_per_metre_x = 2835;
		h_.info.pixel_per_metre_y = 2835;
		size_t palette_size = palette.size();
		bool has_palette = palette_size > 1 && palette_size <= 256;
		h_.info.num_colours = has_palette ? uint32_t(palette_size) : 256;
		h_.info.num_important = 0;
		h_.file.data_offset = int(sizeof(magic_bytes) + sizeof(FileHeader)) + h_.info.header_size + h_.info.num_colours * 4;
		h_.file.file_size = h_.file.data_offset + h_.info.size;

		f.write((const char*)&magic_bytes, sizeof(magic_bytes));
		f.write((const char*)&h_, sizeof(h_));
		if ((f.bad() | f.eof() | f.fail()))
			return Error::FILE_IO_ERROR;
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
			f.write((const char*)image.ptr(), h_.info.size);
		}
		else
		{
			uint32_t width_in_bytes = (image.width() * h_.info.bits_per_pixel + 7) / 8u;
			uint32_t padding = stride - width_in_bytes;
			n = 0;
			for (unsigned y = 0; y < h_.info.height; ++y)
			{
				f.write((const char*)image.ptr(y), width_in_bytes);
				f.write((const char*)&n, padding);
			}
		}

		if ((f.bad() | f.eof() | f.fail()))
			return  Error::FILE_IO_ERROR;
		return  Error::NONE;
	}


	Error cleanup(std::ifstream& f, size_t pos, Error error)
	{
		if (error != Error::NONE)
			f.seekg(pos);
		return error;
	}
	Error cleanup(std::ofstream& f, size_t pos, Error error)
	{
		if (error != Error::NONE)
			f.seekp(pos);
		return error;
	}



	bool is_valid_;
	Header h_;
	TMap<uint8_t> data_;
	std::vector<uint32_t> palette_;

};




