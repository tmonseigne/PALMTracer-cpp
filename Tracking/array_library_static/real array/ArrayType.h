/*******************************************************************************************************************
								Pactage for general array use
*******************************************************************************************************************/

#ifdef		_ARRAY_TYPE_


//LIBERARY
#include <fcntl.h>
#include <io.h>
#include <iostream>
#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include "../../generalParam/generalParam.h"
//extern "C"{


#ifndef _HEADER_FILE_
#define _HEADER_FILE_

//HEADER OF FILES
typedef struct
{
	unsigned long type;
	unsigned long nbDim;
	unsigned long sizeData[4];
} HeaderFile;
#endif


class ArrayType
{
private:
	bool dataAllocatedByConstructor;
	unsigned long sizeData[4];
	unsigned long sizeAllocated;

public:
	unsigned long nbDim;
	unsigned long nbPoints;
	//static unsigned long type;
	DataType* data;
#ifdef INTEGER
	signed long minVal;
	signed long maxVal;
#endif

private:
	int getDim(const long i, const long j, const long k, const long l);
	int testDim(const long i, const long j, const long k, const long l);

public:
	ArrayType();
	ArrayType(const long i, const long j = 1L, const long k = 1L, const long l = 1L);
	ArrayType(DataType* _data, const long i, const long j = 1L, const long k = 1L, const long l = 1L);

	int alloc(const long i, const long j = 1L, const long k = 1L, const long l = 1L);
	int alloc(DataType* _data, const long i, const long j = 1L, const long k = 1L, const long l = 1L);
	int realloc(const long i, const long j = 1L, const long k = 1L, const long l = 1L);

	~ArrayType();

	static unsigned long type() { return (NB_TYPE); }
	int size(int di);
	void setSize(int di, int valSize);
	int dim();
	int getDepth();

	DataType& operator[](long i);
	DataType& operator()(long i);
	DataType& operator()(long i, long j);
	DataType& operator()(long i, long j, long k);
	DataType& operator()(long i, long j, long k, long l);


	//**************		+=		*********************=
#ifndef BOOL_TYPE
	ArrayType& operator+=(const char val);
	ArrayType& operator+=(const unsigned char val);
	ArrayType& operator+=(const short val);
	ArrayType& operator+=(const unsigned short val);
	ArrayType& operator+=(const long val);
	ArrayType& operator+=(const unsigned long val);
	ArrayType& operator+=(const float val);
	ArrayType& operator+=(const double val);

	ArrayType& operator+=(const USrtArray& arrayType);
	ArrayType& operator+=(const LngArray& arrayType);
	ArrayType& operator+=(const FltArray& arrayType);
	ArrayType& operator+=(const DblArray& arrayType);
	ArrayType& operator+=(const BlArray& arrayType);
#endif


	//**************		-=		*********************=
#ifndef BOOL_TYPE
	ArrayType& operator-=(const char val);
	ArrayType& operator-=(const unsigned char val);
	ArrayType& operator-=(const short val);
	ArrayType& operator-=(const unsigned short val);
	ArrayType& operator-=(const long val);
	ArrayType& operator-=(const unsigned long val);
	ArrayType& operator-=(const float val);
	ArrayType& operator-=(const double val);

	ArrayType& operator-=(const USrtArray& arrayType);
	ArrayType& operator-=(const LngArray& arrayType);
	ArrayType& operator-=(const FltArray& arrayType);
	ArrayType& operator-=(const DblArray& arrayType);

#endif

	//**************		*=		*********************=
#ifndef BOOL_TYPE
	ArrayType& operator*=(const char val);
	ArrayType& operator*=(const unsigned char val);
	ArrayType& operator*=(const short val);
	ArrayType& operator*=(const unsigned short val);
	ArrayType& operator*=(const long val);
	ArrayType& operator*=(const unsigned long val);
	ArrayType& operator*=(const float val);
	ArrayType& operator*=(const double val);

	ArrayType& operator*=(const USrtArray& arrayType);
	ArrayType& operator*=(const LngArray& arrayType);
	ArrayType& operator*=(const FltArray& arrayType);
	ArrayType& operator*=(const DblArray& arrayType);
	ArrayType& operator*=(const BlArray& arrayType);

#endif

	//**************		/=		*********************=
#ifndef BOOL_TYPE
	ArrayType& operator/=(const char val);
	ArrayType& operator/=(const unsigned char val);
	ArrayType& operator/=(const short val);
	ArrayType& operator/=(const unsigned short val);
	ArrayType& operator/=(const long val);
	ArrayType& operator/=(const unsigned long val);
	ArrayType& operator/=(const float val);
	ArrayType& operator/=(const double val);

	ArrayType& operator/=(const USrtArray& arrayType);
	ArrayType& operator/=(const LngArray& arrayType);
	ArrayType& operator/=(const FltArray& arrayType);
	ArrayType& operator/=(const DblArray& arrayType);

#endif


	//**************		=		*********************=
	ArrayType& operator=(const char val);
	ArrayType& operator=(const unsigned char val);
	ArrayType& operator=(const short val);
	ArrayType& operator=(const unsigned short val);
	ArrayType& operator=(const long val);
	ArrayType& operator=(const unsigned long val);
	ArrayType& operator=(const float val);
	ArrayType& operator=(const double val);
	ArrayType& operator=(const bool val);

	ArrayType& operator=(const USrtArray& arrayType);
	ArrayType& operator=(const LngArray& arrayType);
	ArrayType& operator=(const FltArray& arrayType);
	ArrayType& operator=(const DblArray& arrayType);
	ArrayType& operator=(const BlArray& arrayType);

#ifndef BOOL_TYPE
	ArrayType& operator++(void);
	ArrayType& operator--(void);
#endif
	bool isThereData(void);

#ifndef BOOL_TYPE
#ifdef min
#undef min
#endif
#ifdef max
#undef max
#endif
	DataType min(void);
	DataType max(void);
	void renormalize(double min, double max);
	void renormalizeZD(double zero, double dilation);
	void renormalizeZD(double zero, double dilation, long min, long max);
#endif


#ifndef BOOL_TYPE
	double sigma(void);
	double power(void);
#endif

#ifndef BOOL_TYPE
	double mean(void);
	double sum(void);
#endif

#ifndef BOOL_TYPE
	void hardThresholdPositive(DataType val);
	void hardThresholdPositive(ArrayType& arrayType, DataType val);
	void softThresholdPositive(DataType val);
	void hardThresholdNegative(DataType val);
	void hardThresholdNegative(ArrayType& arrayType, DataType val);
	void hardThresholdInterval(DataType valMin, DataType valMax, bool between);
	void hardThresholdInterval(ArrayType& arrayType, DataType valMin, DataType valMax, bool between);

#ifndef UNSIGNED
	void hardThresholdPosNeg(DataType val);
	void hardThresholdPosNegContinue(DataType val);
	void softThresholdPosNeg(DataType val);
	void continuousThreshold(DataType val);
#endif
#endif
	void writeArray(wchar_t* fileName);
	void readArray(wchar_t* fileName);
};


//} //end of extern "C"

#endif
