#pragma once

class Segmentation
{
protected:
	size_t nbInfoStruct_;

public:
	Segmentation() = default;
	virtual ~Segmentation() = default;
	size_t GetNbStructs() const { return nbInfoStruct_; }
};
