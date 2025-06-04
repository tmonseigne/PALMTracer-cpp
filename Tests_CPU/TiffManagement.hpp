#pragma once

#include <tiffio.h>
#include <vector>
#include <iostream>
#include <stdexcept>

inline void loadTiff3D(const char* filename, std::vector<unsigned short>& data, size_t& width, size_t& height, size_t& planes)
{
	TIFF* tif = TIFFOpen(filename, "r");
	if (!tif) { throw std::runtime_error("Erreur d'ouverture du fichier TIFF"); }

	uint32 w, h;
	uint16 depth;
	size_t frame = 0;

	do {
		TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w);
		TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);
		TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &depth);

		if (depth != 16) { throw std::runtime_error("Le TIFF doit être en 16 bits"); }

		if (frame == 0) {
			width  = w;
			height = h;
		}
		else {
			if (w != width || h != height) {
				TIFFClose(tif);
				throw std::runtime_error("Toutes les images doivent avoir les mêmes dimensions");
			}
		}

		std::vector<unsigned short> buf(w * h);

		for (uint32 row = 0; row < h; ++row) {
			if (TIFFReadScanline(tif, buf.data() + row * w, row) < 0) {
				TIFFClose(tif);
				throw std::runtime_error("Erreur de lecture scanline");
			}
		}

		data.insert(data.end(), buf.begin(), buf.end());
		frame++;
	}
	while (TIFFReadDirectory(tif));

	planes = frame;

	TIFFClose(tif);
}

inline void saveTiff3D(const char* filename, const unsigned short* data, size_t width, size_t height, size_t planes)
{
	TIFF* tif = TIFFOpen(filename, "w");
	if (!tif) { throw std::runtime_error("Erreur à l'ouverture du TIFF pour l'écriture"); }

	for (size_t k = 0; k < planes; ++k) {
		TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, width);
		TIFFSetField(tif, TIFFTAG_IMAGELENGTH, height);
		TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 16);
		TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 1);
		TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, height);
		TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
		TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
		TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);

		for (uint32 row = 0; row < height; ++row) {
			const unsigned short* row_data = data + k * width * height + row * width;
			TIFFWriteScanline(tif, (void*)row_data, row);
		}

		TIFFWriteDirectory(tif);
	}

	TIFFClose(tif);
}
