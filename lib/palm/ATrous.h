#pragma once
#include <vector>


class ATrous
{
public:
	size_t height_, width_, size_;
	std::vector<double> tmp_;
	std::vector<std::vector<double>> wDecomp_;
	ATrous() = default;
	~ATrous() = default;

	int ATrousDecomposition2D(const unsigned short* image, size_t height, size_t width, bool* planesToUse, size_t nbPlane);

	//2D transform:
private:
	[[nodiscard]] static constexpr size_t twoToThe(const size_t val) noexcept { return (val == 0) ? 1 : (static_cast<size_t>(1) << val); }
	[[nodiscard]] size_t offsetTest(int i, int j) const;
	[[nodiscard]] size_t offset(size_t i, size_t j) const;
	int testSizes2D(bool* planesToUse, size_t& nbPlane) const;
	void b3SplineFiltering2D(const std::vector<double>& in, std::vector<double>& out, size_t plane);
};
