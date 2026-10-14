#pragma once

#include <stdexcept>
#include <tiffio.h>
#include <vector>
#include "Structs.hpp"

inline void LoadTiff3D(const char* filename, std::vector<uint16_t>& data, palm_tracer::Vec3I& dim)
{
	TIFF* tif = TIFFOpen(filename, "r");
	if (tif == nullptr) { throw std::runtime_error("Erreur d'ouverture du fichier TIFF"); }

	uint32_t w, h;
	uint16_t depth;
	size_t frame = 0;

	do {
		TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w);
		TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);
		TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &depth);

		if (depth != 16) {
			TIFFClose(tif);
			throw std::runtime_error("Le TIFF doit être en 16 bits");
		}

		if (frame == 0) {
			dim.x = w;
			dim.y = h;
		}
		else {
			if (w != dim.x || h != dim.y) {
				TIFFClose(tif);
				throw std::runtime_error("Toutes les images doivent avoir les mêmes dimensions");
			}
		}

		std::vector<uint16_t> buf(static_cast<size_t>(w) * h);

		for (uint32_t row = 0; row < h; ++row) {
			if (TIFFReadScanline(tif, buf.data() + static_cast<size_t>(row) * w, row) < 0) {
				TIFFClose(tif);
				throw std::runtime_error("Erreur de lecture scanline");
			}
		}

		data.insert(data.end(), buf.begin(), buf.end());
		frame++;
	}
	while (TIFFReadDirectory(tif) != 0);

	dim.z = frame;

	TIFFClose(tif);
}

inline void SaveTiff3D(const char* filename, const uint16_t* data, const palm_tracer::Vec3I& dim)
{
	TIFF* tif = TIFFOpen(filename, "w");
	if (tif == nullptr) { throw std::runtime_error("Erreur à l'ouverture du TIFF pour l'écriture"); }

	for (size_t z = 0; z < dim.z; ++z) {
		TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, dim.x);
		TIFFSetField(tif, TIFFTAG_IMAGELENGTH, dim.y);
		TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 16);
		TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 1);
		TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, dim.y);
		TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
		TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
		TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);

		for (uint32_t y = 0; y < dim.y; ++y) {
			const uint16_t* row_data = data + z * dim.x * dim.y + y * dim.x;
			TIFFWriteScanline(tif, const_cast<uint16_t*>(row_data), y);
		}

		TIFFWriteDirectory(tif);
	}

	TIFFClose(tif);
}
