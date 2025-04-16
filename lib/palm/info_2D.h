/*******************************************************************************************************************
								Link class to treat data from segmentation of images 2D (16bit only)
*******************************************************************************************************************/
#pragma once

using val2D = struct
{
	double x, y;
};

using box2D = struct
{
	size_t xMin, xMax, yMin, yMax;
};

using wsFusion2D = struct
{
	double maxPotentialWithNeighbor;
	int structureInContact;
	double minPotentialOnTheStructure, maxPotentialOnTheStructure;
	int x, y;
	bool kept, locked;
	int structureToBecaume, newIndex;
};


using maxStruct2D = struct
{
	double val;
	size_t x, y;
};

using varianceStruct2D = struct
{
	double xx, yy, xy;
};

using InfoData2D = struct
{
	size_t peri;
	size_t surf;
	maxStruct2D max;
	val2D moment, momentSurf;
	double mean;
	double intensity;
	varianceStruct2D var;
	val2D sigma;
	double angle;
	size_t index;
	size_t time;
};


class Info2D
{
public:
	size_t peri_;
	double surf_;
	maxStruct2D min_, max_;
	double stdev_;
	val2D moment_, momentSurf_;
	double mean_;
	double intensity_;
	varianceStruct2D var_;
	val2D sigma_;
	double angle_;
	size_t index_;
	box2D boundingBox_;


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
