/*******************************************************************************************************************
								Pactage for array use
*******************************************************************************************************************/


#ifdef		_ARRAY_TYPE_

#include "NoError.h"

//*****************		CONSTRUCTORS		DESTRUCTOR		ALLOC		**************
ArrayType::ArrayType()
{
	data = NULL;
	nbDim = 0;
	nbPoints = 0;
	sizeAllocated = 0;
	//type=NB_TYPE;
#ifdef INTEGER
	minVal = MINVAL;
	maxVal = MAXVAL;
#endif
	dataAllocatedByConstructor = false;
}


ArrayType::ArrayType(const long i, const long j, const long k, const long l)
{
	data = NULL;
	sizeAllocated = 0;
	const int ret = alloc(i, j, k, l);
	if (ret < 0) { exit(ret); }
}


ArrayType::ArrayType(DataType* _data, const long i, const long j, const long k, const long l)
{
	data = NULL;
	sizeAllocated = 0;
	const int ret = alloc(_data, i, j, k, l);
	if (ret < 0) { exit(ret); }
}


bool ArrayType::isThereData() { return data != NULL; }

int ArrayType::getDim(const long i, const long j, const long k, const long l)
{
	int _nbDim = 4;
	if (l == 1) {
		_nbDim = 3;
		if (k == 1) {
			_nbDim = 2;
			if (j == 1) { _nbDim = 1; }
		}
	}
	return _nbDim;
}

int ArrayType::testDim(const long i, const long j, const long k, const long l)
{
	if (i < 0 || j < 0 || k < 0 || l < 0) {
		nbDim = 0;
		nbPoints = 0;
		data = NULL;
		sizeAllocated = 0;
		return SIZE_ERR;
	}
	return RET_OK;
}

int ArrayType::getDepth() { return 8 * sizeof(DataType); }


int ArrayType::alloc(const long i, const long j, const long k, const long l)
{
	long ret;
	nbDim = getDim(i, j, k, l);

	if ((ret = testDim(i, j, k, l)) < 0) { return ret; }


	if (isThereData()) {
		if (sizeAllocated >= (unsigned long)(i * j * k * l)) {
			nbPoints = i * j * k * l;
			sizeData[0] = i;
			sizeData[1] = j;
			sizeData[2] = k;
			sizeData[3] = l;
		}
		else {
			nbDim = 0;
			nbPoints = 0;
			data = NULL;
			return SIZE_ERR;
		}
	}
	else {
		data = new DataType[i * j * k * l];
		if (data == NULL) {
			nbDim = 0;
			nbPoints = 0;
			sizeAllocated = 0;
			return ALLOC_ERR;
		}
		else {
			sizeAllocated = (unsigned long)(i * j * k * l);
			nbPoints = i * j * k * l;
			sizeData[0] = i;
			sizeData[1] = j;
			sizeData[2] = k;
			sizeData[3] = l;
			dataAllocatedByConstructor = true;
			//type=NB_TYPE;
#ifdef INTEGER
			minVal = MINVAL;
			maxVal = MAXVAL;
#endif
		}
	}
	return RET_OK;
}


int ArrayType::alloc(DataType* _data, const long i, const long j, const long k, const long l)
{
	long ret;
	nbDim = getDim(i, j, k, l);

	if ((ret = testDim(i, j, k, l)) < 0) { return ret; }

	if (isThereData()) {
		if (dataAllocatedByConstructor) {
			nbDim = 0;
			nbPoints = 0;
			sizeAllocated = 0;
			return ALLOC_ERR;
		}
		else {
			nbPoints = (unsigned long)(i * j * k * l);
			sizeAllocated = 0;
			sizeData[0] = i;
			sizeData[1] = j;
			sizeData[2] = k;
			sizeData[3] = l;
			data = _data;
		}
	}
	else {
		sizeAllocated = 0;
		dataAllocatedByConstructor = false;
		data = _data;
		nbPoints = (unsigned long)(i * j * k * l);
		sizeData[0] = i;
		sizeData[1] = j;
		sizeData[2] = k;
		sizeData[3] = l;

#ifdef INTEGER
		minVal = MINVAL;
		maxVal = MAXVAL;
#endif
	}
	return RET_OK;
}

int ArrayType::realloc(const long i, const long j, const long k, const long l)
{
	long ret;
	nbDim = getDim(i, j, k, l);
	if ((ret = testDim(i, j, k, l)) < 0) { return ret; }
	if (isThereData() && !dataAllocatedByConstructor) { return ALLOC_ERR; }
	if (isThereData()) {
		if (sizeAllocated >= (unsigned long)(i * j * k * l)) {
			nbPoints = (unsigned long)(i * j * k * l);
			sizeData[0] = i;
			sizeData[1] = j;
			sizeData[2] = k;
			sizeData[3] = l;
		}
		else {
			delete [] data;
			data = NULL;
			if ((ret = alloc(i, j, k, l)) < 0) { return ret; }
		}
	}
	else if ((ret = alloc(i, j, k, l)) < 0) { return ret; }
	return RET_OK;
}


ArrayType::~ArrayType()
{
	if (dataAllocatedByConstructor) { if (isThereData()) { delete [] data; } }
	else { data = NULL; }
}


//**************		[]		*********************
DataType&ArrayType::operator[](long i) { return data[i]; }

//**************		()		*********************
DataType&ArrayType::operator()(long i) { return data[i]; }

DataType&ArrayType::operator()(long i, long j) { return data[i * sizeData[1] + j]; }

DataType&ArrayType::operator()(long i, long j, long k) { return data[(i * sizeData[1] + j) * sizeData[2] + k]; }

DataType&ArrayType::operator()(long i, long j, long k, long l) { return data[((i * sizeData[1] + j) * sizeData[2] + k) * sizeData[3] + l]; }


//**************		OPERATORS	********************

#define DataType2 unsigned short
#define ArrayType2 USrtArray
#include "operators.cpp"
#undef ArrayType2
#undef DataType2

#define DataType2 float
#define ArrayType2 FltArray
#include "operators.cpp"
#undef ArrayType2
#undef DataType2

#define DataType2 bool
#define BOOL_TYPE1
#define ArrayType2 BlArray
#include "operators.cpp"
#undef BOOL_TYPE1
#undef ArrayType2
#undef DataType2

//**************		++		*********************
#ifndef BOOL_TYPE
ArrayType& ArrayType::operator++()
{
	for (unsigned long i = 0; i < nbPoints; i++) { ++data[i]; }
	return *this;
}
#endif

//**************		--		*********************
#ifndef BOOL_TYPE
ArrayType& ArrayType::operator--()
{
	for (unsigned long i = 0; i < nbPoints; i++) { ++data[i]; }
	return *this;
}
#endif


//*************		standard methods		*****************

int ArrayType::size(int di) { return sizeData[di]; }

void ArrayType::setSize(int di, int valSize)
{
	sizeData[di] = valSize;
	nbDim = getDim(sizeData[0], sizeData[1], sizeData[2], sizeData[3]);
	nbPoints = (unsigned long)(sizeData[0] * sizeData[1] * sizeData[2] * sizeData[3]);
	if (nbPoints > sizeAllocated) { exit(SIZE_ERR); }
}

int ArrayType::dim() { return nbDim; }


#ifndef BOOL_TYPE
DataType ArrayType::min()
{
	if (!isThereData()) { exit(ALLOC_ERR); }
	DataType _min = data[0];
	for (unsigned long i = 0; i < nbPoints; i++) { if (_min > data[i]) { _min = data[i]; } }
	return _min;
}
#endif

#ifndef BOOL_TYPE
DataType ArrayType::max()
{
	if (!isThereData()) { exit(ALLOC_ERR); }
	DataType _max = data[0];
	for (unsigned long i = 0; i < nbPoints; i++) { if (_max < data[i]) { _max = data[i]; } }
	return _max;
}
#endif

#ifndef BOOL_TYPE
double ArrayType::sigma()
{
	double mean = 0.;
	unsigned long i;
	double var = 0.;
	for (i = 0; i < nbPoints; i++) { mean += data[i]; }
	mean /= (double)nbPoints;
	for (i = 0; i < nbPoints; i++) { var += (data[i] - mean) * (data[i] - mean); }
	var /= nbPoints;
	if (var < 0) { return 0.; }
	else { return sqrt(var); }
}
#endif

#ifndef BOOL_TYPE
double ArrayType::power()
{
	double pow = 0.;
	for (unsigned long i = 0; i < nbPoints; i++) { pow += data[i] * data[i]; }
	pow /= (double)nbPoints;
	return pow;
}
#endif


#ifndef BOOL_TYPE
double ArrayType::mean()
{
	double m = 0.;
	for (unsigned long i = 0; i < nbPoints; i++) { m += data[i]; }
	m /= (double)nbPoints;
	return m;
}
#endif

#ifndef BOOL_TYPE
double ArrayType::sum()
{
	double m = 0.;
	for (unsigned long i = 0; i < nbPoints; i++) { m += data[i]; }
	return m;
}
#endif
#endif
