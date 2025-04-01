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
#include "tools/generalParam.h"


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
	static int getDim(long i, long j, long k, long l);
	int testDim(long i, long j, long k, long l);

public:
	ArrayType();
	ArrayType(long i, long j = 1L, long k = 1L, long l = 1L);
	ArrayType(DataType* _data, long i, long j = 1L, long k = 1L, long l = 1L);

	int alloc(long i, long j = 1L, long k = 1L, long l = 1L);
	int alloc(DataType* _data, long i, long j = 1L, long k = 1L, long l = 1L);
	int realloc(long i, long j = 1L, long k = 1L, long l = 1L);

	~ArrayType();

	static unsigned long type() { return NB_TYPE; }
	int size(int di);
	void setSize(int di, int valSize);
	int dim();
	static int getDepth();

	DataType& operator[](long i);
	DataType& operator()(long i);
	DataType& operator()(long i, long j);
	DataType& operator()(long i, long j, long k);
	DataType& operator()(long i, long j, long k, long l);


	//**************		+=		*********************=
#ifndef BOOL_TYPE
	ArrayType& operator+=(char val);
	ArrayType& operator+=(unsigned char val);
	ArrayType& operator+=(short val);
	ArrayType& operator+=(unsigned short val);
	ArrayType& operator+=(long val);
	ArrayType& operator+=(unsigned long val);
	ArrayType& operator+=(float val);
	ArrayType& operator+=(double val);
	ArrayType& operator+=(const USrtArray& arrayType);
	ArrayType& operator+=(const FltArray& arrayType);
	ArrayType& operator+=(const BlArray& arrayType);
#endif


	//**************		-=		*********************=
#ifndef BOOL_TYPE
	ArrayType& operator-=(char val);
	ArrayType& operator-=(unsigned char val);
	ArrayType& operator-=(short val);
	ArrayType& operator-=(unsigned short val);
	ArrayType& operator-=(long val);
	ArrayType& operator-=(unsigned long val);
	ArrayType& operator-=(float val);
	ArrayType& operator-=(double val);
	ArrayType& operator-=(const USrtArray& arrayType);
	ArrayType& operator-=(const FltArray& arrayType);
#endif

	//**************		*=		*********************=
#ifndef BOOL_TYPE
	ArrayType& operator*=(char val);
	ArrayType& operator*=(unsigned char val);
	ArrayType& operator*=(short val);
	ArrayType& operator*=(unsigned short val);
	ArrayType& operator*=(long val);
	ArrayType& operator*=(unsigned long val);
	ArrayType& operator*=(float val);
	ArrayType& operator*=(double val);
	ArrayType& operator*=(const USrtArray& arrayType);
	ArrayType& operator*=(const FltArray& arrayType);
	ArrayType& operator*=(const BlArray& arrayType);
#endif

	//**************		/=		*********************=
#ifndef BOOL_TYPE
	ArrayType& operator/=(char val);
	ArrayType& operator/=(unsigned char val);
	ArrayType& operator/=(short val);
	ArrayType& operator/=(unsigned short val);
	ArrayType& operator/=(long val);
	ArrayType& operator/=(unsigned long val);
	ArrayType& operator/=(float val);
	ArrayType& operator/=(double val);
	ArrayType& operator/=(const USrtArray& arrayType);
	ArrayType& operator/=(const FltArray& arrayType);
#endif

	//**************		=		*********************=
	ArrayType& operator=(char val);
	ArrayType& operator=(unsigned char val);
	ArrayType& operator=(short val);
	ArrayType& operator=(unsigned short val);
	ArrayType& operator=(long val);
	ArrayType& operator=(unsigned long val);
	ArrayType& operator=(float val);
	ArrayType& operator=(double val);
	ArrayType& operator=(bool val);
	ArrayType& operator=(const USrtArray& arrayType);
	ArrayType& operator=(const FltArray& arrayType);
	ArrayType& operator=(const BlArray& arrayType);

#ifndef BOOL_TYPE
	ArrayType& operator++();
	ArrayType& operator--();
#endif
	bool isThereData();

#ifndef BOOL_TYPE
#ifdef min
#undef min
#endif
#ifdef max
#undef max
#endif
	DataType min();
	DataType max();
#endif


#ifndef BOOL_TYPE
	double sigma();
	double power();
#endif

#ifndef BOOL_TYPE
	double mean();
	double sum();
#endif

};


#endif
