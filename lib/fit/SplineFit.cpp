#include "SplineFit.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <execution>
#include <numbers>
#include <numeric>
#include <ranges>
#include <vector>

namespace palm_tracer {
void SplineFit::fit(const SegmentationInfo& infos, std::array<double, N_LOCALIZATIONS>& point) { Fitting::fit(infos, point); }
}
