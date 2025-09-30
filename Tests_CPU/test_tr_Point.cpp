#include "gtest/gtest.h"

#include <algorithm>
#include <cmath>
#include <numeric>

#include "Structs.hpp"
#include "tracking/Point.hpp"

namespace palm_tracer {
//-------------------------------------------------------------------------------------------------
TEST(PointTests, Base)
{
	const Point P1;
	std::cout << "P1 : " << P1 << "\n";
	EXPECT_EQ(P1.id_, 0);
	EXPECT_EQ(P1.plane_, 0);
	EXPECT_DOUBLE_EQ(P1.pos_.x, 0.0);
	EXPECT_DOUBLE_EQ(P1.pos_.y, 0.0);
	EXPECT_DOUBLE_EQ(P1.pos_.z, 0.0);
	EXPECT_DOUBLE_EQ(P1.intensity_, 0.0);
	EXPECT_DOUBLE_EQ(P1.surface_, 0.0);

	const Point P2 = Point(1, 2, {.x = 3, .y = 4, .z = 5}, 6, 7);
	std::cout << "P2 : " << P2 << "\n";
	EXPECT_EQ(P2.id_, 1);
	EXPECT_EQ(P2.plane_, 2);
	EXPECT_DOUBLE_EQ(P2.pos_.x, 3);
	EXPECT_DOUBLE_EQ(P2.pos_.y, 4);
	EXPECT_DOUBLE_EQ(P2.pos_.z, 5);
	EXPECT_DOUBLE_EQ(P2.intensity_, 6);
	EXPECT_DOUBLE_EQ(P2.surface_, 7);
	EXPECT_DOUBLE_EQ(P2.GetDistance(P1), std::sqrt(50));

}
} // namespace palm_tracer
