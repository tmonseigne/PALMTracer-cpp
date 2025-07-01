#include "gtest/gtest.h"

#include <algorithm>

#include "Spline.hpp"
#include "TiffManagement.hpp"

namespace palm_tracer {
//-------------------------------------------------------------------------------------------------
TEST(Spline, BaseCreation) { Spline spline; }

//-------------------------------------------------------------------------------------------------
TEST(Spline, Calibrate)
{
	Spline spline;

	std::vector<uint16_t> input;
	size_t width, height, planes;
	// 1. Charger
	LoadTiff3D("C:/Git/PALMTracer-cpp/Tests_CPU/Input/Beads.tiff", input, height, width, planes);
}
} // namespace palm_tracer
