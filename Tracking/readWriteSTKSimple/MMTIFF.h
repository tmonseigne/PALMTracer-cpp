/*******************************************************************************************************************
								Class TIFFImage to read and write TIFF and STK for 8 16 & 32bits unsigned integer
*******************************************************************************************************************/


#ifndef __TIFF_IMAGE__
#define __TIFF_IMAGE__

#include <stdlib.h>
#include <string.h>
#include "../DebugMIA/print_debug.h"
#include "../generalParam/NoError.h"

typedef unsigned char TIFFTYPE_BYTE;
typedef char TIFFTYPE_ASCII;
typedef unsigned short TIFFTYPE_SHORT;
typedef unsigned long TIFFTYPE_LONG;

typedef struct
{
	TIFFTYPE_LONG numerator;
	TIFFTYPE_LONG denominator;
} TIFFTYPE_RATIONAL;

typedef char TIFFTYPE_SBYTE;
typedef char TIFFTYPE_UNDEFINED;
typedef short TIFFTYPE_SSHORT;
typedef long TIFFTYPE_SLONG;

typedef struct
{
	TIFFTYPE_SLONG numerator;
	TIFFTYPE_SLONG denominator;
} TIFFTYPE_SRATIONAL;

typedef float TIFFTYPE_FLOAT;
typedef double TIFFTYPE_DOUBLE;

//#define DEBUG_RW_STK

class MMTAGS;


class TAG
{ //TIFF Tag
public:
	TAG(FILE* in);
	TAG(unsigned short _tag, unsigned short _type, unsigned long _count, unsigned long _val, void* _valBuf);
	~TAG();
	unsigned short tag;  //tag num
	unsigned short type; //tag type
	unsigned int count;
	unsigned int val; // or offset if count>1
	void* valBuf;
	unsigned long offset;
	int readTag(FILE* in, MMTAGS* mm);
	int writeTag(FILE* in);
	int writeTagBuf(FILE* in, MMTAGS* mm);

private:
	int readMMUIC1Tag(FILE* in, MMTAGS* mm);
	int readMMUIC2Tag(FILE* in, MMTAGS* mm);
	int readMMUIC3Tag(FILE* in, MMTAGS* mm);
	int readMMUIC4Tag(FILE* in, MMTAGS* mm);
	int writeMMUIC1TagBuf(FILE* in, MMTAGS* mm);
	int writeMMUIC2TagBuf(FILE* in, MMTAGS* mm);
	int writeMMUIC3TagBuf(FILE* in, MMTAGS* mm);
	int writeMMUIC4TagBuf(FILE* in, MMTAGS* mm);
};


class MMTAGS
{ //special MM UIC Tags
public:
	MMTAGS();
	~MMTAGS();
	void** bufMM14;
	char fromMM1OrMM4[48];
	unsigned long writtingOffsetMM1[48];
	unsigned long nbPlanes;
	unsigned long* bufMM2;
	unsigned long* bufMM3;
};


class IFD
{ //image file directory
private:
	bool isMMTag(TAG* _tag);

public:
	IFD(unsigned long _nbPlanes, unsigned long _width, unsigned long _height, unsigned long _depth, unsigned long _samplesPerPixel);
	IFD(unsigned long offsetIFD, FILE* in);
	~IFD();
	TAG** tags;
	unsigned short tagCount;
	MMTAGS* mmTAGS;
	void print();
	TAG* getTag(unsigned long tag, unsigned long& val);
	int addNewTag(TAG* tag);
	int writeAllBufTags(FILE* in);
	int writeIfd(FILE* in);
};


class TIFFImage
{
public:
	TIFFImage();
	TIFFImage(unsigned long nbPlanes, unsigned long width, unsigned long height, unsigned long depth, unsigned long samplesPerPixel);
	~TIFFImage();
	int readImageInfo(char* fileName);
	void printTags();

	void* allocData();
	int readData();
	int readData(void* buffer); //for data already allocated
	int writeImage(char* fileName, void* data);
	int writeImage(char* fileName);
	unsigned long getNbPlanes();
	unsigned long getWidth();
	unsigned long getHeight();
	unsigned long getDepth();           // 8 or 16
	unsigned long getSamplesPerPixel(); //1 (grey levels) or 3 channels (RGB)
	unsigned long getSizeOfData();
	void* getBuffer();
	int addToString(char* stringToAdd);
	unsigned long getFirstPlaneOffset();

private:
	IFD* ifd;
	FILE* in;
	char fileName[1024];
	unsigned long firstPlaneOffset;
	unsigned long nbPlanes;
	unsigned long width;
	unsigned long height;
	unsigned long depth;           // 8 or 16
	unsigned long samplesPerPixel; //1 (grey levels) or 3 channels (RGB)
	unsigned long sizeOfData;
	void* buffer;
};


#endif //__TIFF_IMAGE__
