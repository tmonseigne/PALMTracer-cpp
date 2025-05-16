/*******************************************************************************************************************
								Link class to treat data from segmentation of images 2D (16bit only)
*******************************************************************************************************************/
#pragma once

struct val_2d
{
	double x, y;
};

struct box_2d
{
	size_t xMin, xMax, yMin, yMax;
};

struct ws_fusion_2d
{
	double maxPotentialWithNeighbor, minPotential, maxPotential;
	int x, y, structInContact, structToBecome, newIndex;
	bool kept, locked;
};


struct max_2d
{
	double val;
	size_t x, y;
};

struct variance_2d
{
	double xx, yy, xy;
};

class Info2D
{
public:
	size_t index_;
	size_t peri_;
	double surf_;
	max_2d min_, max_;
	double stdev_;
	val_2d moment_, momentSurf_;
	double mean_;
	double intensity_;
	variance_2d var_;
	val_2d sigma_;
	double angle_;
	box_2d boundingBox_;


	Info2D() = default;
	~Info2D() = default;

	// Copy / Move
	Info2D(Info2D&&) noexcept = default;
	Info2D& operator=(Info2D&&) noexcept = default;
	Info2D(const Info2D&) = default;
	Info2D& operator=(const Info2D&) = default;

	void Reset(const size_t height = 0, const size_t width = 0)
	{
		peri_ = 0;
		surf_ = 0;
		min_ = {.val = 0, .x = 0, .y = 0};
		max_ = {.val = 0, .x = 0, .y = 0};
		stdev_ = 0;
		moment_ = {.x = 0, .y = 0};
		momentSurf_ = {.x = 0, .y = 0};
		mean_ = 0;
		intensity_ = 0;
		var_ = {.xx = 0, .yy = 0, .xy = 0};
		sigma_ = {.x = 0, .y = 0};
		angle_ = 0;
		boundingBox_ = {.xMin = height, .xMax = 0, .yMin = width, .yMax = 0};
	}
};
