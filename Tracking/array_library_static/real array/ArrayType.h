

/*******************************************************************************************************************
								Pactage for general array use
*******************************************************************************************************************/

#ifdef		_ARRAY_TYPE_


//LIBERARY
	#include <stdio.h>
	#include <iostream>
	#include <stdlib.h>
	#include <math.h>
	#include <io.h>
	#include <fcntl.h>
	#include <sys/stat.h>
	#include <memory.h>
	#include "../../DebugMIA/print_debug.h"
	#include "../../generalParam/generalParam.h"
	#include "../../DebugMIA/controlAlloc.h"
//extern "C"{



#ifndef _HEADER_FILE_
#define _HEADER_FILE_
//HEADER OF FILES
typedef struct{
		unsigned long type;
		unsigned long nbDim;
		unsigned long sizeData[4];
}HeaderFile;
#endif


class ArrayType{
private:
		bool dataAllocatedByConstructor; 
		unsigned long sizeData[4];
		unsigned long sizeAllocated;
public:
		unsigned long nbDim;
		unsigned long nbPoints;
		//static unsigned long type;
		DataType	*data;
#ifdef INTEGER
		signed long minVal;
		signed long maxVal;
#endif

private:
		
		#ifdef	DEBUG_AT
		
		const wchar_t* typeName(){ return((wchar_t*)&TYPE_NAME); } 
		
		void writeInfo(const wchar_t* s,const long i, const long j, const long k, const long l);
		void writeInfo(DataType	*_data, const wchar_t* s, const long i, const long j, const long k, const long l);
		#endif
		int getDim(const long i, const long j, const long k, const long l);
		int testDim(const long i, const long j, const long k, const long l);
public:

		ArrayType();
		ArrayType(const long i, const long j=1L, const long k=1L, const long l=1L);
		ArrayType(DataType	*_data, const long i, const long j=1L, const long k=1L, const long l=1L);

		int		alloc(const long i, const long j=1L, const long k=1L, const long l=1L);
		int		alloc(DataType	*_data, const long i, const long j=1L, const long k=1L, const long l=1L);
		int		realloc(const long i, const long j=1L, const long k=1L, const long l=1L);

		~ArrayType();

		static unsigned long type(){return(NB_TYPE);} 
		int size(int di);
		void setSize(int di, int valSize);
		int dim();
		int getDepth();
#ifdef COMPLEX
		void conjugate();
#endif
		DataType &operator[](long i);
		DataType &operator()(long i);
		DataType &operator()(long i,long j);
		DataType &operator()(long i,long j, long k);
		DataType &operator()(long i,long j, long k, long l);


		//**************		+=		*********************=
#ifndef BOOL_TYPE
		ArrayType &operator+=(const char val);
		ArrayType &operator+=(const unsigned char val);
		ArrayType &operator+=(const short val);
		ArrayType &operator+=(const unsigned short val);
		ArrayType &operator+=(const long val);
		ArrayType &operator+=(const unsigned long val);
		ArrayType &operator+=(const float val);
		ArrayType &operator+=(const double val);
#ifdef COMPLEX
		ArrayType &operator+=(const ComplexChr val);
		ArrayType &operator+=(const ComplexSrt val);
		ArrayType &operator+=(const ComplexLng val);
		ArrayType &operator+=(const ComplexFlt val);
		ArrayType &operator+=(const ComplexDbl val);
#endif

		ArrayType &operator+=(const ChrArray &arrayType);
		ArrayType &operator+=(const UChrArray &arrayType);
		ArrayType &operator+=(const SrtArray &arrayType);
		ArrayType &operator+=(const USrtArray &arrayType);
		ArrayType &operator+=(const LngArray &arrayType);
		ArrayType &operator+=(const ULngArray &arrayType);
		ArrayType &operator+=(const FltArray &arrayType);
		ArrayType &operator+=(const DblArray &arrayType);
		ArrayType &operator+=(const BlArray &arrayType);
#ifdef COMPLEX
		ArrayType &operator+=(const CxChrArray &arrayType);
		ArrayType &operator+=(const CxSrtArray &arrayType);
		ArrayType &operator+=(const CxLngArray &arrayType);
		ArrayType &operator+=(const CxFltArray &arrayType);
		ArrayType &operator+=(const CxDblArray &arrayType);
#endif
#endif


		//**************		-=		*********************=
#ifndef BOOL_TYPE
		ArrayType &operator-=(const char val);
		ArrayType &operator-=(const unsigned char val);
		ArrayType &operator-=(const short val);
		ArrayType &operator-=(const unsigned short val);
		ArrayType &operator-=(const long val);
		ArrayType &operator-=(const unsigned long val);
		ArrayType &operator-=(const float val);
		ArrayType &operator-=(const double val);
#ifdef COMPLEX
		ArrayType &operator-=(const ComplexChr val);
		ArrayType &operator-=(const ComplexSrt val);
		ArrayType &operator-=(const ComplexLng val);
		ArrayType &operator-=(const ComplexFlt val);
		ArrayType &operator-=(const ComplexDbl val);
#endif

		ArrayType &operator-=(const ChrArray &arrayType);
		ArrayType &operator-=(const UChrArray &arrayType);
		ArrayType &operator-=(const SrtArray &arrayType);
		ArrayType &operator-=(const USrtArray &arrayType);
		ArrayType &operator-=(const LngArray &arrayType);
		ArrayType &operator-=(const ULngArray &arrayType);
		ArrayType &operator-=(const FltArray &arrayType);
		ArrayType &operator-=(const DblArray &arrayType);
#ifdef COMPLEX
		ArrayType &operator-=(const CxChrArray &arrayType);
		ArrayType &operator-=(const CxSrtArray &arrayType);
		ArrayType &operator-=(const CxLngArray &arrayType);
		ArrayType &operator-=(const CxFltArray &arrayType);
		ArrayType &operator-=(const CxDblArray &arrayType);
#endif
#endif

		//**************		*=		*********************=
#ifndef BOOL_TYPE
		ArrayType &operator*=(const char val);
		ArrayType &operator*=(const unsigned char val);
		ArrayType &operator*=(const short val);
		ArrayType &operator*=(const unsigned short val);
		ArrayType &operator*=(const long val);
		ArrayType &operator*=(const unsigned long val);
		ArrayType &operator*=(const float val);
		ArrayType &operator*=(const double val);
#ifdef COMPLEX
		ArrayType &operator*=(const ComplexChr val);
		ArrayType &operator*=(const ComplexSrt val);
		ArrayType &operator*=(const ComplexLng val);
		ArrayType &operator*=(const ComplexFlt val);
		ArrayType &operator*=(const ComplexDbl val);
#endif

		ArrayType &operator*=(const ChrArray &arrayType);
		ArrayType &operator*=(const UChrArray &arrayType);
		ArrayType &operator*=(const SrtArray &arrayType);
		ArrayType &operator*=(const USrtArray &arrayType);
		ArrayType &operator*=(const LngArray &arrayType);
		ArrayType &operator*=(const ULngArray &arrayType);
		ArrayType &operator*=(const FltArray &arrayType);
		ArrayType &operator*=(const DblArray &arrayType);
		ArrayType &operator*=(const BlArray &arrayType);
#ifdef COMPLEX
		ArrayType &operator*=(const CxChrArray &arrayType);
		ArrayType &operator*=(const CxSrtArray &arrayType);
		ArrayType &operator*=(const CxLngArray &arrayType);
		ArrayType &operator*=(const CxFltArray &arrayType);
		ArrayType &operator*=(const CxDblArray &arrayType);
#endif
#endif

		//**************		/=		*********************=
#ifndef BOOL_TYPE
		ArrayType &operator/=(const char val);
		ArrayType &operator/=(const unsigned char val);
		ArrayType &operator/=(const short val);
		ArrayType &operator/=(const unsigned short val);
		ArrayType &operator/=(const long val);
		ArrayType &operator/=(const unsigned long val);
		ArrayType &operator/=(const float val);
		ArrayType &operator/=(const double val);
#ifdef COMPLEX
		ArrayType &operator/=(const ComplexChr val);
		ArrayType &operator/=(const ComplexSrt val);
		ArrayType &operator/=(const ComplexLng val);
		ArrayType &operator/=(const ComplexFlt val);
		ArrayType &operator/=(const ComplexDbl val);
#endif
		ArrayType &operator/=(const ChrArray &arrayType);
		ArrayType &operator/=(const UChrArray &arrayType);
		ArrayType &operator/=(const SrtArray &arrayType);
		ArrayType &operator/=(const USrtArray &arrayType);
		ArrayType &operator/=(const LngArray &arrayType);
		ArrayType &operator/=(const ULngArray &arrayType);
		ArrayType &operator/=(const FltArray &arrayType);
		ArrayType &operator/=(const DblArray &arrayType);
#ifdef COMPLEX
		ArrayType &operator/=(const CxChrArray &arrayType);
		ArrayType &operator/=(const CxSrtArray &arrayType);
		ArrayType &operator/=(const CxLngArray &arrayType);
		ArrayType &operator/=(const CxFltArray &arrayType);
		ArrayType &operator/=(const CxDblArray &arrayType);
#endif
#endif


		//**************		=		*********************=
		ArrayType &operator=(const char val);
		ArrayType &operator=(const unsigned char val);
		ArrayType &operator=(const short val);
		ArrayType &operator=(const unsigned short val);
		ArrayType &operator=(const long val);
		ArrayType &operator=(const unsigned long val);
		ArrayType &operator=(const float val);
		ArrayType &operator=(const double val);
		ArrayType &operator=(const bool val);
#ifdef COMPLEX
		ArrayType &operator=(const ComplexChr val);
		ArrayType &operator=(const ComplexSrt val);
		ArrayType &operator=(const ComplexLng val);
		ArrayType &operator=(const ComplexFlt val);
		ArrayType &operator=(const ComplexDbl val);
#endif

		ArrayType &operator=(const ChrArray &arrayType);
		ArrayType &operator=(const UChrArray &arrayType);
		ArrayType &operator=(const SrtArray &arrayType);
		ArrayType &operator=(const USrtArray &arrayType);
		ArrayType &operator=(const LngArray &arrayType);
		ArrayType &operator=(const ULngArray &arrayType);
		ArrayType &operator=(const FltArray &arrayType);
		ArrayType &operator=(const DblArray &arrayType);
		ArrayType &operator=(const BlArray &arrayType);
#ifdef COMPLEX
		ArrayType &operator=(const CxChrArray &arrayType);
		ArrayType &operator=(const CxSrtArray &arrayType);
		ArrayType &operator=(const CxLngArray &arrayType);
		ArrayType &operator=(const CxFltArray &arrayType);
		ArrayType &operator=(const CxDblArray &arrayType);
#endif

#ifndef BOOL_TYPE
#ifndef	COMPLEX
		ArrayType &operator++(void);
		ArrayType &operator--(void);
#endif
#endif
		bool	isThereData(void);
		
#ifndef BOOL_TYPE
#ifndef	COMPLEX
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
#endif


#ifndef BOOL_TYPE
#ifndef	COMPLEX
		double	sigma(void);
		double	power(void);
#endif
#endif

#ifndef BOOL_TYPE
#ifndef	COMPLEX
		double	mean(void);
		double	sum(void);
#endif
#endif

#ifndef BOOL_TYPE
#ifndef	COMPLEX		
		void hardThresholdPositive( DataType val );
		void hardThresholdPositive( ArrayType &arrayType, DataType val );
		void softThresholdPositive( DataType val );
		void hardThresholdNegative( DataType val );
		void hardThresholdNegative( ArrayType &arrayType, DataType val );
		void hardThresholdInterval( DataType valMin, DataType valMax, bool between);
		void hardThresholdInterval( ArrayType &arrayType, DataType valMin, DataType valMax, bool between);

#ifndef UNSIGNED
		void hardThresholdPosNeg( DataType val );
		void hardThresholdPosNegContinue( DataType val );
		void softThresholdPosNeg( DataType val );
		void continuousThreshold(DataType val);
#endif
#endif
#endif
		void writeArray(wchar_t* fileName);
		void readArray(wchar_t* fileName);

};


//} //end of extern "C"

#endif

