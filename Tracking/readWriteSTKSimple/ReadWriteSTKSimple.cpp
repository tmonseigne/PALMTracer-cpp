#include "stdafx.h"

/*******************************************************************************************************************
								Read and write TIFF and MM STK for 8 16 & 32bits unsigned integer
*******************************************************************************************************************/

#include "ReadWriteSTKSimple.h"
#include "MMTIFF.h"

int readInfoSTK(char* fileName, unsigned long& depth, unsigned long& samplesPerPixel, unsigned long& nbPlanes, unsigned long& width, unsigned long& height, unsigned long& sizeOfData)
{
	int ret;
	TIFFImage i;
	if ((ret = i.readImageInfo(fileName)) < 0) return (ret);

	depth = i.getDepth();
	samplesPerPixel = i.getSamplesPerPixel();
	nbPlanes = i.getNbPlanes();

	width = i.getWidth();
	height = i.getHeight();
	sizeOfData = i.getSizeOfData();
	return (RET_OK);
}


int readSTK(char* fileName, unsigned char* image)
{
	int ret;
	TIFFImage i;
	if ((ret = i.readImageInfo(fileName)) < 0) return (ret);
	if ((long)i.getSamplesPerPixel() != 1) return (ERR_BAD_FORMAT);
	if ((ret = i.readData(image)) < 0) return (ret);
	return (RET_OK);
}


int readSTK32(char* fileName, unsigned long* image)
{
	int ret;
	TIFFImage i;
	if ((ret = i.readImageInfo(fileName)) < 0) return (ret);
	if ((long)i.getDepth() != 32) return (ERR_BAD_FORMAT);
	if ((long)i.getSamplesPerPixel() != 1) return (ERR_BAD_FORMAT);
	if ((ret = i.readData(image)) < 0) return (ret);
	return (RET_OK);
}


int readSTK16(char* fileName, unsigned short* image)
{
	int ret;
	TIFFImage i;
	if ((ret = i.readImageInfo(fileName)) < 0) return (ret);
	if ((long)i.getDepth() != 16) return (ERR_BAD_FORMAT);
	if ((long)i.getSamplesPerPixel() != 1) return (ERR_BAD_FORMAT);
	if ((ret = i.readData(image)) < 0) return (ret);
	return (RET_OK);
}


int readSTK8(char* fileName, unsigned char* image)
{
	int ret;
	TIFFImage i;
	if ((ret = i.readImageInfo(fileName)) < 0) return (ret);
	if ((long)i.getDepth() != 8) return (ERR_BAD_FORMAT);
	if ((long)i.getSamplesPerPixel() != 1) return (ERR_BAD_FORMAT);
	if ((ret = i.readData(image)) < 0) return (ret);
	return (RET_OK);
}

#define SIZE_BUF	4096


int writeSTK(char* fileNameDest, char* fileNameSource, unsigned char* image)
{
	int ret;
	TIFFImage imTIFF;
	if ((ret = imTIFF.readImageInfo(fileNameSource)) < 0) return (ret);
	// to do
	FILE* source;
	FILE* dest;
	unsigned long readDataCount, writeDataCount;
	source = fopen(fileNameSource, "rb");
	if (source == NULL) { return (ERR_READ_STK); }

	dest = fopen(fileNameDest, "wb");
	if (dest == NULL) { return (ERR_WRITE_STK); }
	static unsigned char buff[SIZE_BUF];
	long stripOffset, sizeTmp;
	int i;
	stripOffset = imTIFF.getFirstPlaneOffset();
	for (i = 0; i < stripOffset; i += SIZE_BUF) {
		if (i + SIZE_BUF >= stripOffset) { sizeTmp = stripOffset - i; }
		else { sizeTmp = SIZE_BUF; }
		if (fread(buff, sizeof(char), sizeTmp, source) != (unsigned long)sizeTmp) { return (ERR_READ_STK); }
		if (fwrite(buff, sizeof(char), sizeTmp, dest) != (unsigned long)sizeTmp) { return (ERR_WRITE_STK); }
	}

	if (fwrite(image, sizeof(char), imTIFF.getSizeOfData(), dest) != (unsigned long)imTIFF.getSizeOfData()) { return (ERR_WRITE_STK); }
	fseek(source, ftell(dest),SEEK_SET);
	readDataCount = SIZE_BUF;
	for (i = 0; (readDataCount == SIZE_BUF) && (readDataCount > 0); i += SIZE_BUF) {
		readDataCount = fread(buff, sizeof(char),SIZE_BUF, source);
		if (readDataCount > 0) { if ((writeDataCount = fwrite(buff, sizeof(char), readDataCount, dest)) != (unsigned long)readDataCount) { return (ERR_WRITE_STK); } }
	}
	fclose(source);
	fclose(dest);
	return (RET_OK);
}


int writeSTK16(char* fileNameDest, char* fileNameSource, unsigned short* image)
{
	int ret;
	TIFFImage imTIFF;
	if ((ret = imTIFF.readImageInfo(fileNameSource)) < 0) return (ret);
	// to do
	FILE* source;
	FILE* dest;
	unsigned long readDataCount, writeDataCount;
	source = fopen(fileNameSource, "rb");
	if (source == NULL) { return (ERR_READ_STK); }

	dest = fopen(fileNameDest, "wb");
	if (dest == NULL) { return (ERR_WRITE_STK); }
	static unsigned char buff[SIZE_BUF];
	long stripOffset, sizeTmp;
	int i;
	stripOffset = imTIFF.getFirstPlaneOffset();
	for (i = 0; i < stripOffset; i += SIZE_BUF) {
		if (i + SIZE_BUF >= stripOffset) { sizeTmp = stripOffset - i; }
		else { sizeTmp = SIZE_BUF; }
		if (fread(buff, sizeof(char), sizeTmp, source) != (unsigned long)sizeTmp) { return (ERR_READ_STK); }
		if (fwrite(buff, sizeof(char), sizeTmp, dest) != (unsigned long)sizeTmp) { return (ERR_WRITE_STK); }
	}

	if (fwrite((unsigned char*)image, sizeof(char), imTIFF.getSizeOfData(), dest) != (unsigned long)imTIFF.getSizeOfData()) { return (ERR_WRITE_STK); }
	fseek(source, ftell(dest),SEEK_SET);
	readDataCount = SIZE_BUF;
	for (i = 0; (readDataCount == SIZE_BUF) && (readDataCount > 0); i += SIZE_BUF) {
		readDataCount = fread(buff, sizeof(char),SIZE_BUF, source);
		if (readDataCount > 0) { if ((writeDataCount = fwrite(buff, sizeof(char), readDataCount, dest)) != (unsigned long)readDataCount) { return (ERR_WRITE_STK); } }
	}
	fclose(source);
	fclose(dest);
	return (RET_OK);
}

int writeSTK32(char* fileNameDest, char* fileNameSource, unsigned long* image)
{
	int ret;
	TIFFImage imTIFF;
	if ((ret = imTIFF.readImageInfo(fileNameSource)) < 0) return (ret);
	// to do
	FILE* source;
	FILE* dest;
	unsigned long readDataCount, writeDataCount;
	source = fopen(fileNameSource, "rb");
	if (source == NULL) { return (ERR_READ_STK); }

	dest = fopen(fileNameDest, "wb");
	if (dest == NULL) { return (ERR_WRITE_STK); }
	static unsigned char buff[SIZE_BUF];
	long stripOffset, sizeTmp;
	int i;
	stripOffset = imTIFF.getFirstPlaneOffset();
	for (i = 0; i < stripOffset; i += SIZE_BUF) {
		if (i + SIZE_BUF >= stripOffset) { sizeTmp = stripOffset - i; }
		else { sizeTmp = SIZE_BUF; }
		if (fread(buff, sizeof(char), sizeTmp, source) != (unsigned long)sizeTmp) { return (ERR_READ_STK); }
		if (fwrite(buff, sizeof(char), sizeTmp, dest) != (unsigned long)sizeTmp) { return (ERR_WRITE_STK); }
	}

	if (fwrite((unsigned char*)image, sizeof(char), imTIFF.getSizeOfData(), dest) != (unsigned long)imTIFF.getSizeOfData()) { return (ERR_WRITE_STK); }
	fseek(source, ftell(dest),SEEK_SET);
	readDataCount = SIZE_BUF;
	for (i = 0; (readDataCount == SIZE_BUF) && (readDataCount > 0); i += SIZE_BUF) {
		readDataCount = fread(buff, sizeof(char),SIZE_BUF, source);
		if (readDataCount > 0) { if ((writeDataCount = fwrite(buff, sizeof(char), readDataCount, dest)) != (unsigned long)readDataCount) { return (ERR_WRITE_STK); } }
	}
	fclose(source);
	fclose(dest);
	return (RET_OK);
}


int writeSTK8(char* fileNameDest, char* fileNameSource, unsigned char* image)
{
	int ret;
	TIFFImage imTIFF;
	if ((ret = imTIFF.readImageInfo(fileNameSource)) < 0) return (ret);
	// to do
	FILE* source;
	FILE* dest;
	unsigned long readDataCount, writeDataCount;
	source = fopen(fileNameSource, "rb");
	if (source == NULL) { return (ERR_READ_STK); }

	dest = fopen(fileNameDest, "wb");
	if (dest == NULL) { return (ERR_WRITE_STK); }
	static unsigned char buff[SIZE_BUF];
	long stripOffset, sizeTmp;
	int i;
	stripOffset = imTIFF.getFirstPlaneOffset();
	for (i = 0; i < stripOffset; i += SIZE_BUF) {
		if (i + SIZE_BUF >= stripOffset) { sizeTmp = stripOffset - i; }
		else { sizeTmp = SIZE_BUF; }
		if (fread(buff, sizeof(char), sizeTmp, source) != (unsigned long)sizeTmp) { return (ERR_READ_STK); }
		if (fwrite(buff, sizeof(char), sizeTmp, dest) != (unsigned long)sizeTmp) { return (ERR_WRITE_STK); }
	}

	if (fwrite((unsigned char*)image, sizeof(char), imTIFF.getSizeOfData(), dest) != (unsigned long)imTIFF.getSizeOfData()) { return (ERR_WRITE_STK); }
	fseek(source, ftell(dest),SEEK_SET);
	readDataCount = SIZE_BUF;
	for (i = 0; (readDataCount == SIZE_BUF) && (readDataCount > 0); i += SIZE_BUF) {
		readDataCount = fread(buff, sizeof(char),SIZE_BUF, source);
		if (readDataCount > 0) { if ((writeDataCount = fwrite(buff, sizeof(char), readDataCount, dest)) != (unsigned long)readDataCount) { return (ERR_WRITE_STK); } }
	}
	fclose(source);
	fclose(dest);
	return (RET_OK);
}
