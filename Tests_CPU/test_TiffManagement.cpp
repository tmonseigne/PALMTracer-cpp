#include <gtest/gtest.h>

#include <filesystem>
#include <memory>
#include <vector>

#include "TiffManagement.hpp"

namespace palm_tracer {
/// <summary>Le rejet d'un TIFF 8 bits libère le fichier, notamment pour permettre sa suppression sous Windows.</summary>
TEST(TiffManagementTests, RejectsWrongDepthAndClosesFile)
{
	const auto Directory = std::filesystem::path(__FILE__).parent_path() / "Output";
	std::filesystem::create_directories(Directory);
	const auto Path = Directory / "unsupported-depth.tiff";
	std::unique_ptr<TIFF, decltype(&TIFFClose)> writer(TIFFOpen(Path.string().c_str(), "w"), TIFFClose);
	ASSERT_NE(writer, nullptr);
	TIFFSetField(writer.get(), TIFFTAG_IMAGEWIDTH, uint32_t{1});
	TIFFSetField(writer.get(), TIFFTAG_IMAGELENGTH, uint32_t{1});
	TIFFSetField(writer.get(), TIFFTAG_BITSPERSAMPLE, 8);
	TIFFSetField(writer.get(), TIFFTAG_SAMPLESPERPIXEL, 1);
	TIFFSetField(writer.get(), TIFFTAG_ROWSPERSTRIP, uint32_t{1});
	TIFFSetField(writer.get(), TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
	TIFFSetField(writer.get(), TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
	uint8_t pixel = 42;
	ASSERT_GE(TIFFWriteScanline(writer.get(), &pixel, 0), 0);
	writer.reset();
	std::vector<uint16_t> image;
	Vec3I dim{};
	EXPECT_THROW(LoadTiff3D(Path.string().c_str(), image, dim), std::runtime_error);
	EXPECT_TRUE(std::filesystem::remove(Path));
}
} // namespace palm_tracer
