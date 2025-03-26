#include "stdafx.h"


/*******************************************************************************************************************
								Pactage for array use
*******************************************************************************************************************/


#ifdef		_ARRAY_TYPE_

extern FILE* debugStreamALS;

#ifdef DEBUG_AT
extern ControlAllocArrayType *controlAllocArrayTypePublic;
#endif


//*****************		CONSTRUCTORS		DESTRUCTOR		ALLOC		**************
ArrayType::ArrayType()
{
#ifdef	DEBUG_AT
	openDebugStreamALS();
	fprintf(debugStreamALS, "\t** ArrayType() **, type:%ls, this:0x%p\n", typeName(), this);
	closeDebugStreamALS();
#endif
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

#ifdef DEBUG_AT
	if(controlAllocArrayTypePublic!=NULL)
		controlAllocArrayTypePublic->allocNewArray(this);
#endif
}


ArrayType::ArrayType(const long i, const long j, const long k, const long l)
{
	int ret;
	data = NULL;
	sizeAllocated = 0;
	ret = alloc(i, j, k, l);
	if (ret < 0) {
		openDebugStreamALS();
		fprintf(debugStreamALS, "�� � ERROR this:0x%p constructor: %d � ��\n", this, ret);
		closeDebugStreamALS();
		exit(ret);
	}
#ifdef DEBUG_AT
	if(controlAllocArrayTypePublic!=NULL)
		controlAllocArrayTypePublic->allocNewArray(this);
#endif
}


ArrayType::ArrayType(DataType* _data, const long i, const long j, const long k, const long l)
{
	int ret;
	data = NULL;
	sizeAllocated = 0;
	ret = alloc(_data, i, j, k, l);
	if (ret < 0) {
		openDebugStreamALS();
		fprintf(debugStreamALS, "�� � ERROR this:0x%p constructor: %d � ��\n", this, ret);
		closeDebugStreamALS();
		exit(ret);
	}
#ifdef DEBUG_AT
	if(controlAllocArrayTypePublic!=NULL)
		controlAllocArrayTypePublic->allocNewArray(this);
#endif
}


bool ArrayType::isThereData(void) { return (data != NULL); }

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
	return (_nbDim);
}

#ifdef	DEBUG_AT
void ArrayType::writeInfo(const wchar_t* s, const long i, const long j, const long k, const long l){
	openDebugStreamALS();
	switch(nbDim){
	case 1: 
		fwprintf(debugStreamALS, L"\t** %ls(%d) **, type:%ls, this:0x%p\n", s, i, typeName(), this); 
		break;
	case 2: 
		fwprintf(debugStreamALS, L"\t** %ls(%d, %d) **, type:%ls, this:0x%p\n", s, i,j, typeName(), this); 
		break;
	case 3: 
		fwprintf(debugStreamALS, L"\t** %ls(%d, %d, %d) **, type:%ls, this:0x%p\n", s, i,j,k, typeName(), this); 
		break;
	case 4: 
		fwprintf(debugStreamALS, L"\t** %ls(%d, %d, %d, %d) **, type:%ls, this:0x%p\n", s, i,j,k,l, typeName(), this); 
		break;
	}
	closeDebugStreamALS();
}

void ArrayType::writeInfo(DataType	*_data, const wchar_t* s, const long i, const long j, const long k, const long l){
	openDebugStreamALS();
	switch(nbDim){
	case 1: 
		fwprintf(debugStreamALS, L"\t** %ls(0x%p, %d) **, type:%ls, this:0x%p\n", s, _data, i, typeName(), this); 
		break;
	case 2: 
		fwprintf(debugStreamALS, L"\t** %ls(0x%p, %d, %d) **, type:%ls, this:0x%p\n", s, _data, i,j, typeName(), this); 
		break;
	case 3: 
		fwprintf(debugStreamALS, L"\t** %ls(0x%p, %d, %d, %d) **, type:%ls, this:0x%p\n", s, _data, i,j,k, typeName(), this); 
		break;
	case 4: 
		fwprintf(debugStreamALS, L"\t** %ls(0x%p, %d, %d, %d, %d) **, type:%ls, this:0x%p\n", s, _data, i,j,k,l, typeName(), this); 
		break;
	}
	closeDebugStreamALS();
}
#endif

int ArrayType::testDim(const long i, const long j, const long k, const long l)
{
	if ((i < 0) || (j < 0) || (k < 0) || (l < 0)) {
		nbDim = 0;
		nbPoints = 0;
		data = NULL;
		sizeAllocated = 0;
		return (SIZE_ERR);
	}
	return (RET_OK);
}

int ArrayType::getDepth()
{
#ifdef COMPLEX
		return(16*sizeof(DataTypeSimple));
#else
	return (8 * sizeof(DataType));
#endif
}


int ArrayType::alloc(const long i, const long j, const long k, const long l)
{
	long ret;
	nbDim = getDim(i, j, k, l);
#ifdef	DEBUG_AT
	writeInfo(L"alloc",i,j,k,l);
#endif
	if ((ret = testDim(i, j, k, l)) < 0) return (ret);


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
			return (SIZE_ERR);
		}
	}
	else {
		data = new DataType[i * j * k * l];
		if (data == NULL) {
			nbDim = 0;
			nbPoints = 0;
			sizeAllocated = 0;
			return (ALLOC_ERR);
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
#ifdef	DEBUG_AT

#endif
	return (RET_OK);
}


int ArrayType::alloc(DataType* _data, const long i, const long j, const long k, const long l)
{
	long ret;
	nbDim = getDim(i, j, k, l);
#ifdef	DEBUG_AT
	writeInfo(_data, L"alloc", i,j,k,l);
#endif
	if ((ret = testDim(i, j, k, l)) < 0) return (ret);

	if (isThereData()) {
		if (dataAllocatedByConstructor) {
			nbDim = 0;
			nbPoints = 0;
			sizeAllocated = 0;

			openDebugStreamALS();
			fprintf(debugStreamALS, "�� � ALLOC_ERR in ArrayType::realloc, this:0x%p � ��\n", this);
			fprintf(debugStreamALS, "�� � ArrayType::alloc(data, i,j,... ) can not be used where data have already been allocated by constructor � ��\n", this);
			closeDebugStreamALS();
			return (ALLOC_ERR);
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
		//type=NB_TYPE;
#ifdef INTEGER
		minVal = MINVAL;
		maxVal = MAXVAL;
#endif
	}
	return (RET_OK);
}

int ArrayType::realloc(const long i, const long j, const long k, const long l)
{
	long ret;
	nbDim = getDim(i, j, k, l);
#ifdef	DEBUG_AT
	writeInfo(L"realloc", i,j,k,l);
#endif
	if ((ret = testDim(i, j, k, l)) < 0) return (ret);
	if (isThereData() && !dataAllocatedByConstructor) {
		openDebugStreamALS();
		fprintf(debugStreamALS, "�� � ALLOC_ERR in ArrayType::realloc, this:0x%p � ��\n", this);
		fprintf(debugStreamALS, "�� � ArrayType::realloc can be used if a ArrayType::alloc(data, i,j,... ) have been executed � ��\n", this);
		closeDebugStreamALS();
		return (ALLOC_ERR);
	}
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
			if ((ret = alloc(i, j, k, l)) < 0) return (ret);
		}
	}
	else if ((ret = alloc(i, j, k, l)) < 0) return (ret);
	return (RET_OK);
}


ArrayType::~ArrayType()
{
#ifdef	DEBUG_AT
	long i;
	openDebugStreamALS();
	fprintf(debugStreamALS, "\t** destructor ~ArrayType() (%d units to deallocate (", (dataAllocatedByConstructor ? sizeAllocated:0)); /* ))**, type:%ls, this:0x%p\n",*/
	//		sizeAllocated ,typeName(), this);
	for(i=0;i<(signed)nbDim-1;i++)
		fprintf(debugStreamALS, "%d,", sizeData[i]);
	if(nbDim!=0) fprintf(debugStreamALS, "%d", sizeData[nbDim-1]);
	fprintf(debugStreamALS, ") **, type:%ls, this:0x%p\n", typeName(), this);

	closeDebugStreamALS();
#endif
	if (dataAllocatedByConstructor) { if (isThereData()) delete [] data; }
	else data = NULL;

#ifdef DEBUG_AT
	if(controlAllocArrayTypePublic!=NULL)
		controlAllocArrayTypePublic->deallocArray(this);
#endif
}


//**************		[]		*********************
DataType&ArrayType::operator[](long i)
{
#ifdef	DEBUG_AT
	if((i>=(long)nbPoints) || (i<0)){
		openDebugStreamALS();
		fprintf(debugStreamALS, "�� � ERROR operator[] � ��\n");
		fprintf(debugStreamALS, "\t\t>>> i=%d ,sizeData[0]=%d, ptr=0x%p\n",i, sizeData[0],this);
		closeDebugStreamALS();
		exit(BOUND_ERR);
	}
#endif
	return (data[i]);
}

//**************		()		*********************
DataType&ArrayType::operator()(long i)
{
#ifdef	DEBUG_AT
	if((i>=(long)sizeData[0]) || (i<0) ){
		openDebugStreamALS();
		fprintf(debugStreamALS, "�� � BOUND_ERR operator(%d) max(%d), this:%p � ��\n", 
			i,sizeData[0],this);
		closeDebugStreamALS();
		exit(BOUND_ERR);
	}
	/*if(nbDim!=1){
		openDebugStreamALS();
		fprintf(debugStreamALS, "�� � DIM_ERR operator(i) � ��\n");
		closeDebugStreamALS();
		exit(DIM_ERR);
	}*/
#endif
	return (data[i]);
}

DataType&ArrayType::operator()(long i, long j)
{
#ifdef	DEBUG_AT
	if((i>=(long)sizeData[0]) || (i<0) || (j>=(long)sizeData[1]) || (j<0)){
		openDebugStreamALS();
		fprintf(debugStreamALS, "�� � BOUND_ERR operator(%d,%d) max(%d,%d), this:%p � ��\n", 
			i,j,sizeData[0],sizeData[1],this);
		closeDebugStreamALS();
		exit(BOUND_ERR);
	}
	/*if(nbDim!=2){
		openDebugStreamALS();
		fprintf(debugStreamALS, "�� � DIM_ERR operator(i,j) � ��\n");
		closeDebugStreamALS();
		exit(DIM_ERR);
	}*/
#endif
	return (data[i * sizeData[1] + j]);
}

DataType&ArrayType::operator()(long i, long j, long k)
{
#ifdef	DEBUG_AT
	if((i>=(long)sizeData[0]) || (i<0) || (j>=(long)sizeData[1]) || (j<0) ||
			(k>=(long)sizeData[2]) || (k<0)){
		openDebugStreamALS();
		fprintf(debugStreamALS, "�� � BOUND_ERR operator(%d,%d,%d) max(%d,%d,%d), this:%p � ��\n", 
			i,j,k,sizeData[0],sizeData[1],sizeData[2],this);
		closeDebugStreamALS();
		exit(BOUND_ERR);
	}
	/*if(nbDim!=3){
		openDebugStreamALS();
		fprintf(debugStreamALS, "�� � DIM_ERR operator(i,j,k) � ��\n");
		closeDebugStreamALS();
		exit(DIM_ERR);
	}*/
#endif
	return (data[(i * sizeData[1] + j) * sizeData[2] + k]);
}

DataType&ArrayType::operator()(long i, long j, long k, long l)
{
#ifdef	DEBUG_AT
	if((i>=(long)sizeData[0]) || (i<0) || (j>=(long)sizeData[1]) || (j<0) ||
			(k>=(long)sizeData[2]) || (k<0) || (l>=(long)sizeData[3]) || (l<0)){
		openDebugStreamALS();
		fprintf(debugStreamALS, "�� � BOUND_ERR operator(%d,%d,%d,%d) max(%d,%d,%d,%d), this:0x%p � ��\n", 
			i,j,k,l,sizeData[0],sizeData[1],sizeData[2],sizeData[3],this);
		closeDebugStreamALS();
		exit(BOUND_ERR);
	}
	/*if(nbDim!=4){
		openDebugStreamALS();
		fprintf(debugStreamALS, "�� � DIM_ERR operator(i,j,k,l) � ��\n");
		closeDebugStreamALS();
		exit(DIM_ERR);
	}*/
#endif
	return (data[((i * sizeData[1] + j) * sizeData[2] + k) * sizeData[3] + l]);
}


//**************		OPERATORS	********************

#define DataType2 unsigned char
#define ArrayType2 UChrArray
#include "operators.cpp"
#undef ArrayType2
#undef DataType2

#define DataType2 short
#define ArrayType2 SrtArray
#include "operators.cpp"
#undef ArrayType2
#undef DataType2

#define DataType2 unsigned short
#define ArrayType2 USrtArray
#include "operators.cpp"
#undef ArrayType2
#undef DataType2

#define DataType2 long
#define ArrayType2 LngArray
#include "operators.cpp"
#undef ArrayType2
#undef DataType2

#define DataType2 unsigned long
#define ArrayType2 ULngArray
#include "operators.cpp"
#undef ArrayType2
#undef DataType2

#define DataType2 float
#define ArrayType2 FltArray
#include "operators.cpp"
#undef ArrayType2
#undef DataType2

#define DataType2 double
#define ArrayType2 DblArray
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

#ifdef COMPLEX
	#define COMPLEX2

	#define DataType2 ComplexChr
	#define ArrayType2 CxChrArray
	#include "operators.cpp"
	#undef ArrayType2
	#undef DataType2

	#define DataType2 ComplexSrt
	#define ArrayType2 CxSrtArray
	#include "operators.cpp"
	#undef ArrayType2
	#undef DataType2

	#define DataType2 ComplexLng
	#define ArrayType2 CxLngArray
	#include "operators.cpp"
	#undef ArrayType2
	#undef DataType2

	#define DataType2 ComplexFlt
	#define ArrayType2 CxFltArray
	#include "operators.cpp"
	#undef ArrayType2
	#undef DataType2

	#define DataType2 ComplexDbl
	#define ArrayType2 CxDblArray
	#include "operators.cpp"
	#undef ArrayType2
	#undef DataType2
	
	#undef COMPLEX2
#endif


//**************		++		*********************
#ifndef	COMPLEX
#ifndef BOOL_TYPE
ArrayType& ArrayType::operator++(void)
{
	for (unsigned long i = 0; i < nbPoints; i++) ++data[i];
	return (*this);
}
#endif
#endif

//**************		--		*********************
#ifndef	COMPLEX
#ifndef BOOL_TYPE
ArrayType& ArrayType::operator--(void)
{
	for (unsigned long i = 0; i < nbPoints; i++) ++data[i];
	return (*this);
}
#endif
#endif


//*************		standard methods		*****************


#ifdef	COMPLEX
void ArrayType::conjugate(){
	long l;
	for(l=0;l<(long)nbPoints;l++)
		data[l][1]=-data[l][1];
}
#endif

int ArrayType::size(int di) { return (sizeData[di]); }

void ArrayType::setSize(int di, int valSize)
{
	sizeData[di] = valSize;
	nbDim = getDim(sizeData[0], sizeData[1], sizeData[2], sizeData[3]);
	nbPoints = (unsigned long)(sizeData[0] * sizeData[1] * sizeData[2] * sizeData[3]);
	if (nbPoints > sizeAllocated) {
		openDebugStreamALS();
		fprintf(debugStreamALS, "�� � SIZE_ERR in setSize(int %d di, int %d valSize), this:0x%p � ��\n", di, valSize, this);
		closeDebugStreamALS();
		exit(SIZE_ERR);
	}
}

int ArrayType::dim() { return (nbDim); }


#ifndef BOOL_TYPE
#ifndef	COMPLEX
DataType ArrayType::min(void)
{
	if (!isThereData()) {
		openDebugStreamALS();
		fprintf(debugStreamALS, "�� � ALLOC_ERR in min(void), this:0x%p � ��\n", this);
		closeDebugStreamALS();
		exit(ALLOC_ERR);
	}
	DataType _min = data[0];
	for (unsigned long i = 0; i < nbPoints; i++)
		if (_min > data[i]) _min = data[i];
	return (_min);
}
#endif
#endif

#ifndef BOOL_TYPE
#ifndef	COMPLEX
DataType ArrayType::max(void)
{
	if (!isThereData()) {
		openDebugStreamALS();
		fprintf(debugStreamALS, "�� � ALLOC_ERR in max(void), this:0x%p � ��\n", this);
		closeDebugStreamALS();
		exit(ALLOC_ERR);
	}
	DataType _max = data[0];
	for (unsigned long i = 0; i < nbPoints; i++)
		if (_max < data[i]) _max = data[i];
	return (_max);
}
#endif
#endif

#ifndef BOOL_TYPE
#ifndef	COMPLEX
void ArrayType::renormalize(double _min, double _max)
{
	DataType oldMin, oldMax;
	double val2;
	oldMin = min();
	oldMax = max();
	double val = (_max - _min) / (oldMax - oldMin);

	for (unsigned long i = 0; i < nbPoints; i++) {
		val2 = (data[i] - oldMin) * val + _min;
		data[i] = (DataType)val2;
	}
}
#endif
#endif

#ifndef BOOL_TYPE
#ifndef	COMPLEX
void ArrayType::renormalizeZD(double zero, double dilation)
{
	double val;
	for (unsigned long i = 0; i < nbPoints; i++) {
		val = zero + data[i] * dilation;
		data[i] = (DataType)val;
	}
}
#endif
#endif

#ifndef BOOL_TYPE
#ifndef	COMPLEX
void ArrayType::renormalizeZD(double zero, double dilation, long _min, long _max)
{
	double val;
	for (unsigned long i = 0; i < nbPoints; i++) {
		val = zero + data[i] * dilation;
		if (val < _min) val = (DataType)_min;
		else if (val > _max) val = (DataType)_max;
		data[i] = (DataType)val;
	}
}
#endif
#endif

#ifndef COMPLEX
#ifndef BOOL_TYPE
double ArrayType::sigma(void)
{
	double mean = 0.;
	unsigned long i;
	double var = 0.;
	for (i = 0; i < nbPoints; i++) mean += data[i];
	mean /= (double)nbPoints;
	for (i = 0; i < nbPoints; i++) var += (data[i] - mean) * (data[i] - mean);
	var /= nbPoints;
	if (var < 0) return (0.);
	else return (sqrt(var));
}
#endif
#endif

#ifndef COMPLEX
#ifndef BOOL_TYPE
double ArrayType::power(void)
{
	unsigned long i;
	double pow = 0.;
	for (i = 0; i < nbPoints; i++) pow += (data[i]) * (data[i]);
	pow /= (double)nbPoints;
	return (pow);
}
#endif
#endif


#ifndef COMPLEX
#ifndef BOOL_TYPE
double ArrayType::mean()
{
	unsigned long i;
	double m = 0.;
	for (i = 0; i < nbPoints; i++) m += data[i];
	m /= (double)nbPoints;
	return (m);
}
#endif
#endif

#ifndef COMPLEX
#ifndef BOOL_TYPE
double ArrayType::sum()
{
	unsigned long i;
	double m = 0.;
	for (i = 0; i < nbPoints; i++) m += data[i];
	return (m);
}
#endif
#endif

/*#ifdef COMPLEX
DataType ArrayType::mean(){
	unsigned long i;
	DataType m=(DataTypeComplex)0.;
	for(i=0;i<nbPoints;i++)
		m+=data[i];
	m/=(DataTypeComplex)nbPoints	return(m);
}
#endif

#ifdef COMPLEX
DataType ArrayType::sum(){
	unsigned long i;
	DataType m=(DataTypeComplex)0;
	for(i=0;i<nbPoints;i++)
		m+=data[i];
	return(m);
}
#endif*/


//****************		Threshold		*********************
#ifndef COMPLEX
#ifndef BOOL_TYPE
void ArrayType::hardThresholdPositive(DataType val)
{
	unsigned long i;
	for (i = 0; i < nbPoints; i++)
		if (data[i] < val) data[i] = (DataType)0;
}
#endif
#endif


#ifndef COMPLEX
#ifndef BOOL_TYPE
void ArrayType::hardThresholdPositive(ArrayType& arrayType, DataType val)
{
	unsigned long i;
	for (i = 0; i < nbPoints; i++)
		if (arrayType.data[i] < val) data[i] = (DataType)0;
		else data[i] = arrayType.data[i];
}
#endif
#endif

#ifndef COMPLEX
#ifndef BOOL_TYPE
void ArrayType::hardThresholdNegative(DataType val)
{
	unsigned long i;
	for (i = 0; i < nbPoints; i++)
		if (data[i] > val) data[i] = (DataType)0;
}
#endif
#endif

#ifndef COMPLEX
#ifndef BOOL_TYPE
void ArrayType::hardThresholdNegative(ArrayType& arrayType, DataType val)
{
	unsigned long i;
	for (i = 0; i < nbPoints; i++)
		if (arrayType.data[i] > val) data[i] = (DataType)0;
		else data[i] = arrayType.data[i];
}
#endif
#endif


#ifndef COMPLEX
#ifndef BOOL_TYPE
void ArrayType::softThresholdPositive(DataType val)
{
	unsigned long i;
	for (i = 0; i < nbPoints; i++)
		if (data[i] < val) data[i] = (DataType)0;
		else data[i] -= val;
}
#endif
#endif


#ifndef COMPLEX
#ifndef BOOL_TYPE
void ArrayType::hardThresholdInterval(DataType valMin, DataType valMax, bool between)
{
	unsigned long i;
	if (between) {
		for (i = 0; i < nbPoints; i++)
			if ((data[i] < valMin) || (data[i] > valMax)) data[i] = (DataType)0;
	}
	else
		for (i = 0; i < nbPoints; i++)
			if ((data[i] > valMin) && (data[i] < valMax)) data[i] = (DataType)0;
}
#endif
#endif

#ifndef COMPLEX
#ifndef BOOL_TYPE
void ArrayType::hardThresholdInterval(ArrayType& arrayType, DataType valMin, DataType valMax, bool between)
{
	unsigned long i;
	if (between) {
		for (i = 0; i < nbPoints; i++)
			if ((arrayType.data[i] < valMin) || (arrayType.data[i] > valMax)) data[i] = (DataType)0;
			else data[i] = arrayType.data[i];
	}
	else
		for (i = 0; i < nbPoints; i++)
			if ((arrayType.data[i] > valMin) && (arrayType.data[i] < valMax)) data[i] = (DataType)0;
			else data[i] = arrayType.data[i];
}
#endif
#endif


#ifndef UNSIGNED
#ifndef COMPLEX
#ifndef BOOL_TYPE
void ArrayType::hardThresholdPosNeg(DataType val)
{
	unsigned long i;
	for (i = 0; i < nbPoints; i++)
		if ((data[i] < val) && (data[i] > -val)) data[i] = (DataType)0;
}
#endif
#endif
#endif


#ifndef UNSIGNED
#ifndef COMPLEX
#ifndef BOOL_TYPE
void ArrayType::hardThresholdPosNegContinue(DataType val)
{
	unsigned long i;
	double val_2;
	val_2 = val * val;
	for (i = 0; i < nbPoints; i++) data[i] = (DataType)(data[i] * (1. - (DataType)exp(-data[i] * data[i] * 0.5 / val_2)));
}
#endif
#endif
#endif

#ifndef UNSIGNED
#ifndef COMPLEX
#ifndef BOOL_TYPE
void ArrayType::softThresholdPosNeg(DataType val)
{
	unsigned long i;
	for (i = 0; i < nbPoints; i++)
		if (data[i] < val) {
			if (data[i] > -val) data[i] = (DataType)0;
			else data[i] += (DataType)val;
		}
		else data[i] -= (DataType)val;
}
#endif
#endif
#endif

#ifndef UNSIGNED
#ifndef COMPLEX
#ifndef BOOL_TYPE
void ArrayType::continuousThreshold(DataType val)
{
	unsigned long i;
	double val_2;
	val_2 = val * val;
	for (i = 0; i < nbPoints; i++) data[i] = (DataType)(data[i] * (1. - exp(-data[i] * data[i] * .5 / val_2)));
}
#endif
#endif
#endif

//*******************		input and output		*********************
void ArrayType::writeArray(wchar_t* fileName)
{
	int fileHandle;
	int nbCharWrite;
	HeaderFile header;


	header.nbDim = nbDim;
	header.sizeData[0] = sizeData[0];
	header.sizeData[1] = sizeData[1];
	header.sizeData[2] = sizeData[2];
	header.sizeData[3] = sizeData[3];
	header.type = type();


	if ((fileHandle = _wopen(fileName, O_RDWR | O_CREAT | O_BINARY, S_IWRITE)) != -1) {
		if ((nbCharWrite = write(fileHandle, (void*)&header, sizeof(HeaderFile))) == -1) {
			openDebugStreamALS();
			fprintf(debugStreamALS, "�� � WRITING_FILE_ERR writeArray() � ��\n");
			closeDebugStreamALS();
			exit(WRITING_FILE_ERR);
		}

		//		cout << nbCharWrite << "   " << sizeof(HeaderFile)<< endl;


		if ((nbCharWrite = write(fileHandle, data, sizeof(DataType) * nbPoints)) == -1) {
			openDebugStreamALS();
			fprintf(debugStreamALS, "�� � WRITING_FILE_ERR writeArray() � ��\n");
			closeDebugStreamALS();
			exit(WRITING_FILE_ERR);
		}
		//		cout << nbCharWrite<< " " << sizeof(DataType)*nbPoints<<endl;
#ifdef COMPLEX
		//			cout << data[0][0] << data[0][1] << endl;
#else
		//			cout << data[0]<< endl;
#endif
		close(fileHandle);
	}
	else {
		openDebugStreamALS();
		fprintf(debugStreamALS, "�� � OPENING_FILE_ERR writeArray() � ��\n");
		closeDebugStreamALS();
		exit(OPENING_FILE_ERR);
	}
}


void ArrayType::readArray(wchar_t* fileName)
{
	int fileHandle;


	HeaderFile header;


	if ((fileHandle = _wopen(fileName, O_RDONLY | O_BINARY, S_IREAD)) != -1) {
		int nb;

		if ((nb = read(fileHandle, &header, sizeof(HeaderFile))) == -1) {
			openDebugStreamALS();
			fprintf(debugStreamALS, "�� � READING_FILE_ERR readArray() � ��\n");
			closeDebugStreamALS();
			exit(READING_FILE_ERR);
		}


		switch (header.nbDim) {
			case 1: alloc(header.sizeData[0]);
			case 2: alloc(header.sizeData[0], header.sizeData[1]);
			case 3: alloc(header.sizeData[0], header.sizeData[1], header.sizeData[2]);
			case 4: alloc(header.sizeData[0], header.sizeData[1], header.sizeData[2], header.sizeData[3]);
		}

		if (read(fileHandle, data, sizeof(DataType) * nbPoints) == -1) {
			openDebugStreamALS();
			fprintf(debugStreamALS, "�� � READING_FILE_ERR readArray() � ��\n");
			closeDebugStreamALS();
			exit(READING_FILE_ERR);
		}
		close(fileHandle);
	}
	else {
		openDebugStreamALS();
		fprintf(debugStreamALS, "�� � OPENING_FILE_ERR readArray() � ��\n");
		closeDebugStreamALS();
		exit(OPENING_FILE_ERR);
	}
}


#endif
