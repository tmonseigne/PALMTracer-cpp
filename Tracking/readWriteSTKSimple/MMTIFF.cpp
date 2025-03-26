/*******************************************************************************************************************
								Class TIFFImage to read and write TIFF and STK for 8 16 & 32bits unsigned integer
*******************************************************************************************************************/

#include "MMTIFF.h"


#define	TIFFTAG_SUBFILETYPE		254	/* subfile data descriptor */
#define	TIFFTAG_OSUBFILETYPE		255	/* +kind of data in subfile */
#define	TIFFTAG_IMAGEWIDTH		256	/* image width in pixels */
#define	TIFFTAG_IMAGELENGTH		257	/* image height in pixels */
#define	TIFFTAG_BITSPERSAMPLE		258	/* bits per channel (sample) */
#define	TIFFTAG_COMPRESSION		259	/* data compression technique */
#define	TIFFTAG_PHOTOMETRIC		262	/* photometric interpretation */
#define TIFFTAG_IMAGEDESCRIPTION        270     /* info about image */
#define	TIFFTAG_STRIPOFFSETS		273	/* offsets to data strips */
#define	TIFFTAG_SAMPLESPERPIXEL		277	/* samples per pixel */
#define	TIFFTAG_ROWSPERSTRIP		278	/* rows per strip of data */
#define	TIFFTAG_STRIPBYTECOUNTS		279	/* bytes counts for strips */
#define	TIFFTAG_XRESOLUTION		282	/* pixels/resolution in x */
#define	TIFFTAG_YRESOLUTION		283	/* pixels/resolution in y */
#define	TIFFTAG_RESOLUTIONUNIT		296	/* units of resolutions */
#define	TIFFTAG_SOFTWARE		305	/* name & release */
#define	TIFFTAG_DATETIME		306	/* creation date and time */
#define	TIFFTAG_PREDICTOR		317	/* prediction scheme w/ LZW */

#define TIFFTYPE_BYTE_NO		1
#define TIFFTYPE_ASCII_NO		2
#define TIFFTYPE_SHORT_NO		3
#define TIFFTYPE_LONG_NO		4
#define TIFFTYPE_RATIONAL_NO	5
#define TIFFTYPE_SBYTE_NO		6
#define TIFFTYPE_UNDEFINED_NO	7
#define TIFFTYPE_SSHORT_NO		8
#define TIFFTYPE_SLONG_NO		9
#define TIFFTYPE_SRATIONAL_NO	10
#define TIFFTYPE_FLOAT_NO		11
#define TIFFTYPE_DOUBLE_NO		12

#define MMUIC1 33628
#define MMUIC2 33629
#define MMUIC3 33630
#define MMUIC4 33631

extern FILE* debugStreamALS;

long sizeOfTIFFType(unsigned short type)
{
	switch (type) {
		case TIFFTYPE_BYTE_NO: return (sizeof(TIFFTYPE_BYTE));
		case TIFFTYPE_ASCII_NO: return (sizeof(TIFFTYPE_ASCII));
		case TIFFTYPE_SHORT_NO: return (sizeof(TIFFTYPE_SHORT));
		case TIFFTYPE_LONG_NO: return (sizeof(TIFFTYPE_LONG));
		case TIFFTYPE_RATIONAL_NO: return (sizeof(TIFFTYPE_RATIONAL));
		case TIFFTYPE_SBYTE_NO: return (sizeof(TIFFTYPE_SBYTE));
		case TIFFTYPE_UNDEFINED_NO: return (sizeof(TIFFTYPE_UNDEFINED));
		case TIFFTYPE_SSHORT_NO: return (sizeof(TIFFTYPE_SSHORT));
		case TIFFTYPE_SLONG_NO: return (sizeof(TIFFTYPE_SLONG));
		case TIFFTYPE_SRATIONAL_NO: return (sizeof(TIFFTYPE_SRATIONAL));
		case TIFFTYPE_FLOAT_NO: return (sizeof(TIFFTYPE_FLOAT));
		case TIFFTYPE_DOUBLE_NO: return (sizeof(TIFFTYPE_DOUBLE));
		default: return (0);
	}
}


/*******************************************************************************************************************
								Class TAG as specified in the TIFF specification
*******************************************************************************************************************/
TAG::TAG(FILE* in)
{
	fread(&tag, sizeof(short), 1, in);
	fread(&type, sizeof(short), 1, in);
	fread(&count, sizeof(long), 1, in);
	fread(&val, sizeof(long), 1, in);
	valBuf = NULL;
}

TAG::TAG(unsigned short _tag, unsigned short _type, unsigned long _count, unsigned long _val, void* _valBuf)
{
	tag = _tag;
	type = _type;
	count = _count;
	val = _val;
	if (count * sizeOfTIFFType(type) > 4) {
		if (_valBuf == NULL) { valBuf = NULL; }
		else {
			valBuf = malloc(count * sizeOfTIFFType(type));
			memcpy(valBuf, _valBuf, count * sizeOfTIFFType(type));
		}
	}
	else { valBuf = &(val); }
}

TAG::~TAG() { if (count * sizeOfTIFFType(type) > 4) free(valBuf); }


int TAG::writeMMUIC1TagBuf(FILE* in, MMTAGS* mm)
{
	unsigned long off;
	unsigned long tmp;
	long iD;
	for (long iD = 0; iD < 48; iD++) {
		if (mm->bufMM14[iD] != NULL)         //this iD existes
			if (mm->fromMM1OrMM4[iD] == 1) { //the tag iD has been read by UIC1Tag
				off = ftell(in);
				switch (iD) {
					case 4:
					case 5:
					case 16:
					case 17:
					case 21:
					case 22:
					case 23:
					case 24:
					case 38:
					case 39:
					case 42: // avirer
						//case 46:
						mm->writtingOffsetMM1[iD] = off;
						if (fwrite((unsigned long*)(mm->bufMM14[iD]), sizeof(long), 2, in) != 2) return (ERR_WRITE_STK);
						break;
					case 6:
					case 7:
					case 25: mm->writtingOffsetMM1[iD] = off;
						tmp = (unsigned long)strlen((char*)mm->bufMM14[iD]) + 1; //+1 to add the end of string char					
						if (tmp < 5) {
							unsigned long tmp2 = 5;
							if (fwrite(&tmp2, sizeof(long), 1, in) != 1) return (ERR_WRITE_STK);
							if (fwrite((char*)(mm->bufMM14[iD]), sizeof(char), tmp, in) != tmp) return (ERR_WRITE_STK);
							unsigned char tmp3 = 0;
							for (long i = tmp; i < 5; i++) if (fwrite(&tmp3, sizeof(char), 1, in) != 1) return (ERR_WRITE_STK);
						}
						else {
							if (fwrite(&(tmp), sizeof(long), 1, in) != 1) return (ERR_WRITE_STK);
							if (fwrite((char*)(mm->bufMM14[iD]), sizeof(char), tmp, in) != tmp) return (ERR_WRITE_STK);
						}
						break;

						break;
					case 28:
					case 29: mm->writtingOffsetMM1[iD] = off;
						if (fwrite((unsigned long*)(mm->bufMM14[iD]), sizeof(long), mm->nbPlanes * 4, in) != mm->nbPlanes * 4) return (ERR_WRITE_STK);
						break;
					case 40:
					case 46: mm->writtingOffsetMM1[iD] = off;
						if (fwrite((unsigned long*)(mm->bufMM14[iD]), sizeof(long), mm->nbPlanes * 2, in) != mm->nbPlanes * 2) return (ERR_WRITE_STK);
						break;
					case 41: mm->writtingOffsetMM1[iD] = off;
						if (fwrite((unsigned long*)(mm->bufMM14[iD]), sizeof(long), mm->nbPlanes, in) != mm->nbPlanes) return (ERR_WRITE_STK);
						break;
					//case 42:
					case 43:
					case 44:
					case 45: mm->writtingOffsetMM1[iD] = off;
						if (fwrite((unsigned long*)(mm->bufMM14[iD]), sizeof(long), 1, in) != 1) return (ERR_WRITE_STK);
						break;
				}
			}
	}

	count = 0;
	val = ftell(in);
	for (iD = 0; iD < 48; iD++) {            //now write the tag id
		if (mm->bufMM14[iD] != NULL)         //this iD existes
			if (mm->fromMM1OrMM4[iD] == 1) { //the tag iD has been read by UIC1Tag

				switch (iD) {
					case 0:
					case 1:
					case 2:
					case 3:
					case 8:
					case 9:
					case 10:
					case 11:
					case 12:
					case 13:
					case 14:
					case 15:
					case 18:
					case 19:
					case 20:
					case 26:
					case 27:
					case 30:
					case 31:
					case 32:
					case 33:
					case 34:
					case 35:
					case 36:
					case 47: if (fwrite(&iD, sizeof(long), 1, in) != 1) return (ERR_WRITE_STK);
						if (fwrite((unsigned long*)(mm->bufMM14[iD]), sizeof(long), 1, in) != 1) return (ERR_WRITE_STK);
						count++;
						break;
					case 4:
					case 5:
					case 6:
					case 7:
					case 21:
					case 16:
					case 17:
					case 22:
					case 23:
					case 24:
					case 25:
					case 28:
					case 29:
					case 38:
					case 39:
					case 40:
					case 41:
					case 42:
					case 43:
					case 44:
					case 45:
					case 46: if (fwrite(&iD, sizeof(long), 1, in) != 1) return (ERR_WRITE_STK);
						if (fwrite(mm->writtingOffsetMM1 + iD, sizeof(long), 1, in) != 1) return (ERR_WRITE_STK);
						count++;
						break;
				}
			}
	}

	return (RET_OK);
}

int TAG::readMMUIC1Tag(FILE* in, MMTAGS* mm)
{
	long iD;
	long _val;
	long currentPos;
	long tmp;

	if (type == TIFFTYPE_RATIONAL_NO) { //It indicates a stack file created by an obsolete (pre-2.5) version of MetaMorph, deleting the tag
		tag = 0;
		type = 0;
		count = 0;
		val = 0;
		valBuf = NULL;
		return (RET_OK);
	}
	if (fseek(in, val,SEEK_SET) != 0) return (ERR_READ_STK);

	for (long j = 0; j < (long)count; j++) {
		fread(&iD, sizeof(long), 1, in);
		fread(&_val, sizeof(long), 1, in);
		currentPos = ftell(in);


		switch (iD) {
			case 0:
			case 1:
			case 2:
			case 3:
			case 8:
			case 9:
			case 10:
			case 11:
			case 12:
			case 13:
			case 14:
			case 15:
			case 18:
			case 19:
			case 20:
			case 26:
			case 27:
			case 30:
			case 31:
			case 32:
			case 33:
			case 34:
			case 35:
			case 36:
			case 47: if (mm->bufMM14[iD] == NULL) mm->bufMM14[iD] = malloc(sizeof(long));
				if (mm->bufMM14[iD] != NULL) *((unsigned long*)(mm->bufMM14[iD])) = (unsigned long)_val;
				mm->fromMM1OrMM4[iD] = 1;
				break;

			case 4:
			case 5:
			case 16:
			case 17:
			case 21:
			case 22:
			case 23:
			case 24:
			case 38:
			case 39:
			case 42: // avirer
				//case 46:
				if (fseek(in, _val,SEEK_SET) != 0) return (ERR_READ_STK);
				if (mm->bufMM14[iD] == NULL) {
					mm->bufMM14[iD] = malloc(2 * sizeof(long));
					if (mm->bufMM14[iD] != NULL) fread(mm->bufMM14[iD], sizeof(long), 2, in);
					mm->fromMM1OrMM4[iD] = 1;
				}
				break;

			case 6:
			case 7:
			case 25: if (fseek(in, _val,SEEK_SET) != 0) return (ERR_READ_STK);
				fread(&tmp, sizeof(long), 1, in);
				if (mm->bufMM14[iD] == NULL) {
					mm->bufMM14[iD] = malloc(tmp * sizeof(char));
					if (mm->bufMM14[iD] != NULL) fread(mm->bufMM14[iD], sizeof(char), tmp, in);
					mm->fromMM1OrMM4[iD] = 1;
				}
				break;
			case 28:
			case 29: if (fseek(in, _val,SEEK_SET) != 0) return (ERR_READ_STK);
				if (mm->bufMM14[iD] == NULL) {
					mm->bufMM14[iD] = malloc(mm->nbPlanes * 4 * sizeof(long));
					if (mm->bufMM14[iD] != NULL) fread(mm->bufMM14[iD], sizeof(long), mm->nbPlanes * 4, in);
					mm->fromMM1OrMM4[iD] = 1;
				}
				break;
			case 40:
			case 46: if (fseek(in, _val,SEEK_SET) != 0) return (ERR_READ_STK);
				if (mm->bufMM14[iD] == NULL) {
					mm->bufMM14[iD] = malloc(mm->nbPlanes * 2 * sizeof(long));
					if (mm->bufMM14[iD] != NULL) fread(mm->bufMM14[iD], sizeof(long), mm->nbPlanes * 2, in);
					mm->fromMM1OrMM4[iD] = 1;
				}
				break;
			case 41: if (fseek(in, _val,SEEK_SET) != 0) return (ERR_READ_STK);
				if (mm->bufMM14[iD] == NULL) {
					mm->bufMM14[iD] = malloc(mm->nbPlanes * sizeof(long));
					if (mm->bufMM14[iD] != NULL) fread(mm->bufMM14[iD], sizeof(long), mm->nbPlanes, in);
					mm->fromMM1OrMM4[iD] = 1;
				}
				break;
			//case 42:
			case 43:
			case 44:
			case 45: if (fseek(in, _val,SEEK_SET) != 0) return (ERR_READ_STK);
				if (mm->bufMM14[iD] == NULL) {
					mm->bufMM14[iD] = malloc(sizeof(long));
					if (mm->bufMM14[iD] != NULL) fread(mm->bufMM14[iD], sizeof(long), 1, in);
					mm->fromMM1OrMM4[iD] = 1;
				}
				break;
		}
		if (fseek(in, currentPos,SEEK_SET) != 0) return (ERR_READ_STK);
	}
	return (RET_OK);
}


int TAG::writeMMUIC2TagBuf(FILE* in, MMTAGS* mm)
{
	val = ftell(in);
	if (fwrite(mm->bufMM2, sizeof(long), 6 * mm->nbPlanes, in) != 6 * mm->nbPlanes) return (ERR_WRITE_STK);
	return (RET_OK);
}

int TAG::readMMUIC2Tag(FILE* in, MMTAGS* mm)
{
	mm->nbPlanes = count;
	mm->bufMM2 = (unsigned long*)malloc(mm->nbPlanes * 6 * sizeof(long));
	if (mm->bufMM2 == NULL) return (ERR_READ_STK);
	if (fseek(in, val,SEEK_SET) != 0) return (ERR_READ_STK);
	fread(mm->bufMM2, sizeof(long), 6 * mm->nbPlanes, in);
	return (RET_OK);
}

int TAG::writeMMUIC3TagBuf(FILE* in, MMTAGS* mm)
{
	val = ftell(in);
	if (fwrite(mm->bufMM3, sizeof(long), 2 * mm->nbPlanes, in) != 2 * mm->nbPlanes) return (ERR_WRITE_STK);
	return (RET_OK);
}

int TAG::readMMUIC3Tag(FILE* in, MMTAGS* mm)
{
	mm->nbPlanes = count;
	mm->bufMM3 = (unsigned long*)malloc(mm->nbPlanes * 2 * sizeof(long));
	if (mm->bufMM3 == NULL) return (ERR_READ_STK);
	if (fseek(in, val,SEEK_SET) != 0) return (ERR_READ_STK);
	fread(mm->bufMM3, sizeof(long), 2 * mm->nbPlanes, in);
	return (RET_OK);
}

int TAG::writeMMUIC4TagBuf(FILE* in, MMTAGS* mm)
{
	unsigned long tmp;
	val = ftell(in);
	for (short sID = 1; sID < 48; sID++) {
		if (mm->bufMM14[sID] != NULL)         //this iD existes
			if (mm->fromMM1OrMM4[sID] == 4) { //the tag iD has been read by UIC4Tag

				switch (sID) {
					case 0:
					case 1:
					case 2:
					case 3:
					case 8:
					case 9:
					case 10:
					case 11:
					case 12:
					case 13:
					case 14:
					case 15:
					case 18:
					case 19:
					case 20:
					case 26:
					case 27:
					case 30:
					case 31:
					case 32:
					case 33:
					case 34:
					case 35:
					case 36:
					case 47: if (fwrite(&(sID), sizeof(short), 1, in) != 1) return (ERR_WRITE_STK); //writes sID
						if (fwrite((unsigned long*)(mm->bufMM14[sID]), sizeof(long), 1, in) != 1) return (ERR_WRITE_STK);
						break;
					case 4:
					case 5:
					case 16:
					case 17:
					case 21:
					case 22:
					case 23:
					case 24:
					case 38:
					case 39:
					case 42: // avirer
						//case 46:
						if (fwrite(&(sID), sizeof(short), 1, in) != 1) return (ERR_WRITE_STK); //writes sID
						if (fwrite((unsigned long*)(mm->bufMM14[sID]), sizeof(long), 2, in) != 2) return (ERR_WRITE_STK);
						break;
					case 6:
					case 7:
					case 25: if (fwrite(&(sID), sizeof(short), 1, in) != 1) return (ERR_WRITE_STK); //writes sID
						tmp = (unsigned long)strlen((char*)mm->bufMM14[sID]);
						if (tmp < 5) {
							unsigned long tmp2 = 5;
							if (fwrite(&tmp2, sizeof(long), 1, in) != 1) return (ERR_WRITE_STK);
							if (fwrite((char*)(mm->bufMM14[sID]), sizeof(char), tmp, in) != tmp) return (ERR_WRITE_STK);
							char tmp3 = 0;
							for (long i = tmp; i < 5; i++) if (fwrite(&tmp3, sizeof(char), 1, in) != 1) return (ERR_WRITE_STK);
						}
						else {
							if (fwrite(&(tmp), sizeof(long), 1, in) != 1) return (ERR_WRITE_STK);
							if (fwrite((char*)(mm->bufMM14[sID]), sizeof(char), tmp, in) != tmp) return (ERR_WRITE_STK);
						}
						break;
					case 28:
					case 29: if (fwrite(&(sID), sizeof(short), 1, in) != 1) return (ERR_WRITE_STK); //writes sID
						if (fwrite((unsigned long*)(mm->bufMM14[sID]), sizeof(long), mm->nbPlanes * 4, in) != mm->nbPlanes * 4) return (ERR_WRITE_STK);
						break;
					case 40:
					case 46: if (fwrite(&(sID), sizeof(short), 1, in) != 1) return (ERR_WRITE_STK); //writes sID
						if (fwrite((unsigned long*)(mm->bufMM14[sID]), sizeof(long), mm->nbPlanes * 2, in) != mm->nbPlanes * 2) return (ERR_WRITE_STK);
						break;
					case 41: if (fwrite(&(sID), sizeof(short), 1, in) != 1) return (ERR_WRITE_STK); //writes sID
						if (fwrite((unsigned long*)(mm->bufMM14[sID]), sizeof(long), mm->nbPlanes, in) != mm->nbPlanes) return (ERR_WRITE_STK);
						break;
					//case 42:
					case 43:
					case 44:
					case 45: if (fwrite(&(sID), sizeof(short), 1, in) != 1) return (ERR_WRITE_STK); //writes sID
						if (fwrite((unsigned long*)(mm->bufMM14[sID]), sizeof(long), 1, in) != 1) return (ERR_WRITE_STK);
						break;
				}
			}
	}
	tmp = 0; //to stop the tag buffer see MMUIC4 spec
	if (fwrite(&(tmp), sizeof(long), 1, in) != 1) return (ERR_WRITE_STK);
	return (RET_OK);
}

int TAG::readMMUIC4Tag(FILE* in, MMTAGS* mm)
{
	long i;

	short sID;
	long tmp;
	mm->nbPlanes = count;
	if (fseek(in, val,SEEK_SET) != 0) return (ERR_READ_STK);
	for (fread(&sID, sizeof(short), 1, in); sID != 0; fread(&sID, sizeof(short), 1, in)) {
		mm->fromMM1OrMM4[sID] = 4;
		switch (sID) {
			case 0:
			case 1:
			case 2:
			case 3:
			case 8:
			case 9:
			case 10:
			case 11:
			case 12:
			case 13:
			case 14:
			case 15:
			case 18:
			case 19:
			case 20:
			case 26:
			case 27:
			case 30:
			case 31:
			case 32:
			case 33:
			case 34:
			case 35:
			case 36:
			case 47: if (mm->bufMM14[sID] == NULL) mm->bufMM14[sID] = malloc(sizeof(long));
				if (mm->bufMM14[sID] != NULL) fread(mm->bufMM14[sID], sizeof(long), 1, in);
				break;

			case 4:
			case 5:
			case 16:
			case 17:
			case 21:
			case 22:
			case 23:
			case 24:
			case 38:
			case 39:
			case 42: // avirer
				//case 46:
				if (mm->bufMM14[sID] == NULL) mm->bufMM14[sID] = malloc(2 * sizeof(long));
				if (mm->bufMM14[sID] != NULL) fread(mm->bufMM14[sID], sizeof(long), 2, in);
				break;

			case 6:
			case 7:
			case 25: fread(&tmp, sizeof(long), 1, in);
				if (mm->bufMM14[sID] == NULL) mm->bufMM14[sID] = malloc(tmp * sizeof(char));
				fread(mm->bufMM14[sID], sizeof(char), tmp, in);
				break;
			case 28:
			case 29: if (mm->bufMM14[sID] == NULL) mm->bufMM14[sID] = malloc(mm->nbPlanes * 4 * sizeof(long));
				fread(mm->bufMM14[sID], sizeof(long), mm->nbPlanes * 4, in);
				break;
			case 40:
			case 46: if (mm->bufMM14[sID] == NULL) mm->bufMM14[sID] = malloc(mm->nbPlanes * 2 * sizeof(long));
				fread(mm->bufMM14[sID], sizeof(long), mm->nbPlanes * 2, in);
				break;
			case 41: if (mm->bufMM14[sID] == NULL) mm->bufMM14[sID] = malloc(mm->nbPlanes * sizeof(long));
				fread(mm->bufMM14[sID], sizeof(long), mm->nbPlanes, in);
				break;
			//case 42:
			case 43:
			case 44:
			case 45: if (mm->bufMM14[sID] == NULL) mm->bufMM14[sID] = malloc(sizeof(long));
				if (mm->bufMM14[sID] != NULL) fread(mm->bufMM14[sID], sizeof(long), 1, in);
				break;
			case 37:
				//this tag is not read but to get the next tag to FILE poiter is moved
				for (i = 0; i < (long)mm->nbPlanes; i++) {
					long tmp;
					fread(&tmp, sizeof(long), 1, in); // nb of char
					fseek(in, tmp, SEEK_CUR);
				}
				break;
		}
	}
	return (RET_OK);
}

int TAG::readTag(FILE* in, MMTAGS* mm)
{
	int ret;
	switch (tag) {
		case MMUIC1: //defined by UIC
			if ((ret = readMMUIC1Tag(in, mm)) < 0) return (ret);
			break;
		case MMUIC2: //defined by UIC
			if ((ret = readMMUIC2Tag(in, mm)) < 0) return (ret);
			break;
		case MMUIC3: //defined by UIC
			if ((ret = readMMUIC3Tag(in, mm)) < 0) return (ret);
			break;
		case MMUIC4: //defined by UIC
			if ((ret = readMMUIC4Tag(in, mm)) < 0) return (ret);
			break;

		default: //defined by Adobe TIFF
			if (count * sizeOfTIFFType(type) > 4) {
				valBuf = malloc(count * sizeOfTIFFType(type));
				if (valBuf == NULL) return (ERR_READ_STK);
				if (fseek(in, val,SEEK_SET) != 0) return (ERR_READ_STK);
				if (fread(valBuf, sizeOfTIFFType(type), count, in) != 1) return (ERR_READ_STK);
			}
			else { valBuf = &(val); }
			break;
	}
	return (RET_OK);
}


int TAG::writeTagBuf(FILE* in, MMTAGS* mm)
{
	int ret;
	switch (tag) {
		case MMUIC1: //defined by UIC
			if ((ret = writeMMUIC1TagBuf(in, mm)) < 0) return (ret);
			break;
		case MMUIC2: //defined by UIC
			if ((ret = writeMMUIC2TagBuf(in, mm)) < 0) return (ret);
			break;
		case MMUIC3: //defined by UIC
			if ((ret = writeMMUIC3TagBuf(in, mm)) < 0) return (ret);
			break;
		case MMUIC4: //defined by UIC
			if ((ret = writeMMUIC4TagBuf(in, mm)) < 0) return (ret);
			break;

		default: //defined by Adobe TIFF
			if (valBuf != NULL) {
				if (count * sizeOfTIFFType(type) > 4) {
					val = ftell(in);
					if (fwrite(valBuf, sizeOfTIFFType(type), count, in) != count) return (ERR_WRITE_STK);
				}
			}
			break;
	}
	return (RET_OK);
}

int TAG::writeTag(FILE* in)
{
	if (fwrite(&tag, sizeof(short), 1, in) != 1) return (ERR_WRITE_STK);
	if (fwrite(&type, sizeof(short), 1, in) != 1) return (ERR_WRITE_STK);
	if (fwrite(&count, sizeof(long), 1, in) != 1) return (ERR_WRITE_STK);
	if (fwrite(&val, sizeof(long), 1, in) != 1) return (ERR_WRITE_STK);
	return (RET_OK);
}

/*******************************************************************************************************************
								Class IFD as specified in the TIFF specification
*******************************************************************************************************************/
//TAG::TAG(unsigned short _tag, unsigned short _type, unsigned long _count, unsigned long _val, void* _valBuf)
IFD::IFD(unsigned long _nbPlanes, unsigned long _width, unsigned long _height, unsigned long _depth, unsigned long _samplesPerPixel)
{
	tagCount = 20;
	TAG* tagUICMM1 = NULL;
	tags = (TAG**)malloc(tagCount * sizeof(TAG*));
	for (long i = 0; i < tagCount; i++) tags[i] = NULL;

	unsigned int resolution[] = {1, 1};              //default resolution in x & y
	unsigned int uIC2Default[] = {0, 1, 0, 0, 0, 0}; //znum,zdenom,cdate,ctime,mdate,mtime

	tags[0] = new TAG(TIFFTAG_SUBFILETYPE, TIFFTYPE_SHORT_NO, 1, 0, NULL);
	tags[1] = new TAG(TIFFTAG_IMAGEWIDTH, TIFFTYPE_LONG_NO, 1, _width, NULL);
	tags[2] = new TAG(TIFFTAG_IMAGELENGTH,TIFFTYPE_LONG_NO, 1, _height, NULL);
	tags[3] = new TAG(TIFFTAG_BITSPERSAMPLE,TIFFTYPE_SHORT_NO, 1, _depth, NULL);
	tags[4] = new TAG(TIFFTAG_COMPRESSION,TIFFTYPE_SHORT_NO, 1, 1, NULL);
	tags[5] = new TAG(TIFFTAG_PHOTOMETRIC,TIFFTYPE_SHORT_NO, 1, 1, NULL);
	tags[6] = new TAG(TIFFTAG_STRIPOFFSETS,TIFFTYPE_LONG_NO, 1, 0, NULL);
	tags[7] = new TAG(TIFFTAG_SAMPLESPERPIXEL,TIFFTYPE_SHORT_NO, 1, _samplesPerPixel, NULL);
	tags[8] = new TAG(TIFFTAG_ROWSPERSTRIP,TIFFTYPE_LONG_NO, 1, 0, NULL);
	tags[9] = new TAG(TIFFTAG_STRIPBYTECOUNTS,TIFFTYPE_LONG_NO, 1, 0, NULL);
	tags[10] = new TAG(TIFFTAG_XRESOLUTION,TIFFTYPE_RATIONAL_NO, 1, 0, resolution);
	tags[11] = new TAG(TIFFTAG_YRESOLUTION,TIFFTYPE_RATIONAL_NO, 1, 0, resolution);
	tags[12] = new TAG(TIFFTAG_RESOLUTIONUNIT,TIFFTYPE_SHORT_NO, 1, 1, NULL);
	tagCount = 13;
	if (_nbPlanes > 1) { //MM Tags
		tagCount = 14;
		tags[13] = new TAG(MMUIC2, TIFFTYPE_RATIONAL_NO, _nbPlanes, 0, NULL);
		mmTAGS = new MMTAGS();
		mmTAGS->nbPlanes = _nbPlanes;
		mmTAGS->bufMM2 = (unsigned long*)malloc(_nbPlanes * 6 * sizeof(long));
		for (long i = 0; i < (long)_nbPlanes; i++) memcpy(mmTAGS->bufMM2 + 6 * i, uIC2Default, 6 * sizeof(long));
	}
	else { mmTAGS = NULL; }
}

IFD::IFD(unsigned long offsetIFD, FILE* in)
{
	bool mmTagDetected = false;
	TAG* tagUICMM1 = NULL;
	fseek(in, offsetIFD,SEEK_SET);
	fread(&tagCount, sizeof(short), 1, in);

	tags = (TAG**)malloc(tagCount * sizeof(TAG*));
	for (long i = 0; i < tagCount; i++) {
		tags[i] = new TAG(in);
		if (isMMTag(tags[i])) mmTagDetected = true;
	}
	if (mmTagDetected) mmTAGS = new MMTAGS();
	else mmTAGS = NULL;
	//the MM tag MMUIC1 (if it existes) has to be read after MMUIC4 as reading after all other tags
	unsigned long i;
	for (i = 0; i < tagCount; i++) {
		if (tags[i]->tag == MMUIC1) tagUICMM1 = tags[i];
		else tags[i]->readTag(in, mmTAGS);
	}
	if (tagUICMM1 != NULL) tagUICMM1->readTag(in, mmTAGS);
}

IFD::~IFD()
{
	for (long i = 0; i < tagCount; i++) delete tags[i];
	free(tags);
	delete mmTAGS;
}

TAG* IFD::getTag(unsigned long tag, unsigned long& val)
{
	// try to find 'tag' in the image file directory, return the associated value or offset,or 0 in case of failure
	unsigned short i;

	for (i = 0; i < tagCount; i++) {
		if ((tags[i]->tag) == tag) {
			val = tags[i]->val;
			return (tags[i]);
		}
	}
	return (NULL);
}

int IFD::addNewTag(TAG* tag)
{
	TAG** tagsTmp;
	if ((tagsTmp = (TAG**)malloc(tagCount * sizeof(TAG*))) == NULL) return (ALLOC_ERR);

	for (int i = 0; i < tagCount; i++) tagsTmp[i] = tags[i];
	tagCount++;
	free(tags);
	if ((tags = (TAG**)malloc(tagCount * sizeof(TAG*))) == NULL) return (ALLOC_ERR);
	unsigned long i;
	for (i = 0; i < tagCount; i++) tags[i] = tagsTmp[i];
	tags[tagCount - 1] = tag;
	free(tagsTmp);
	return (RET_OK);
}

int IFD::writeAllBufTags(FILE* in)
{
	int ret;
	fseek(in, 0,SEEK_END);
	for (long i = 0; i < tagCount; i++)
		if ((ret = tags[i]->writeTagBuf(in, mmTAGS)) < 0) {
#ifdef DEBUG_RW_STK
			openDebugStreamALS();
			fprintf(debugStreamALS, "The buffer of the tag 0x%x have a problem in writting\n", tags[i]->tag);
			closeDebugStreamALS();
#endif
			return (ret);
		}
	return (RET_OK);
}

int IFD::writeIfd(FILE* in)
{
	int ret;
	unsigned short ifdEnterie;
	unsigned long ifdOffset;
	unsigned long ifdEnd;
	ifdEnterie = (unsigned short)tagCount;
	ifdOffset = ftell(in);
	ifdEnd = 0;
	if (fwrite(&ifdEnterie, sizeof(short), 1, in) != 1) return (ERR_WRITE_STK); //writting the IFD enterie
	for (long i = 0; i < tagCount; i++) if ((ret = tags[i]->writeTag(in)) < 0) return (ret);
	if (fwrite(&ifdEnd, sizeof(long), 1, in) != 1) return (ERR_WRITE_STK);    //to stop the list of tags. 0 specify the there is not more idf.
	fseek(in, 4,SEEK_SET);                                                    //to go at the beginning of the file to set the offset of the first 
	if (fwrite(&ifdOffset, sizeof(long), 1, in) != 1) return (ERR_WRITE_STK); //
	fseek(in, 0,SEEK_END);
	fclose(in);
	return (RET_OK);
}


void IFD::print()
{
	long i, j, k;

	TIFFTYPE_BYTE* val_TIFFTYPE_BYTE;
	TIFFTYPE_ASCII* val_TIFFTYPE_ASCII;
	TIFFTYPE_SHORT* val_TIFFTYPE_SHORT;
	TIFFTYPE_LONG* val_TIFFTYPE_LONG;
	TIFFTYPE_RATIONAL* val_TIFFTYPE_RATIONAL;
	TIFFTYPE_SBYTE* val_TIFFTYPE_SBYTE;
	TIFFTYPE_UNDEFINED* val_TIFFTYPE_UNDEFINED;
	TIFFTYPE_SSHORT* val_TIFFTYPE_SSHORT;
	TIFFTYPE_SLONG* val_TIFFTYPE_SLONG;
	TIFFTYPE_SRATIONAL* val_TIFFTYPE_SRATIONAL;
	TIFFTYPE_FLOAT* val_TIFFTYPE_FLOAT;
	TIFFTYPE_DOUBLE* val_TIFFTYPE_DOUBLE;

	openDebugStreamALS();

	for (i = 0; i < tagCount; i++) {
		fprintf(debugStreamALS, "TAG 0x%x\n", tags[i]->tag);
		fprintf(debugStreamALS, "\t type: %u ", tags[i]->type);
		switch (tags[i]->type) {
			case TIFFTYPE_BYTE_NO: fprintf(debugStreamALS, "BYTE\t");
				break;
			case TIFFTYPE_ASCII_NO: fprintf(debugStreamALS, "ASCII\t");
				break;
			case TIFFTYPE_SHORT_NO: fprintf(debugStreamALS, "SHORT\t");
				break;
			case TIFFTYPE_LONG_NO: fprintf(debugStreamALS, "LONG\t");
				break;
			case TIFFTYPE_RATIONAL_NO: fprintf(debugStreamALS, "RATIONAL\t");
				break;
			case TIFFTYPE_SBYTE_NO: fprintf(debugStreamALS, "SBYTE\t");
				break;
			case TIFFTYPE_UNDEFINED_NO: fprintf(debugStreamALS, "UNDEFINED\t");
				break;
			case TIFFTYPE_SSHORT_NO: fprintf(debugStreamALS, "SSHORT\t");
				break;
			case TIFFTYPE_SLONG_NO: fprintf(debugStreamALS, "SLONG\t");
				break;
			case TIFFTYPE_SRATIONAL_NO: fprintf(debugStreamALS, "SRATIONAL\t");
				break;
			case TIFFTYPE_FLOAT_NO: fprintf(debugStreamALS, "FLOAT\t");
				break;
			case TIFFTYPE_DOUBLE_NO: fprintf(debugStreamALS, "DOUBLE\t");
				break;
		}
		fprintf(debugStreamALS, "count:%u\n", tags[i]->count);

		switch (tags[i]->tag) {
			case MMUIC1:
			case MMUIC4: for (j = 0; j < 48; j++)
					switch (j) {
						case 0:
						case 1:
						case 2:
						case 3:
						case 8:
						case 9:
						case 10:
						case 11:
						case 12:
						case 13:
						case 14:
						case 15:
						case 18:
						case 19:
						case 20:
						case 26:
						case 27:
						case 30:
						case 31:
						case 32:
						case 33:
						case 34:
						case 35:
						case 36:
						case 47: if (mmTAGS->bufMM14[j] != NULL) fprintf(debugStreamALS, "MMTagID:%d\t%d\n", j, *(long*)mmTAGS->bufMM14[j]);
							break;
						case 4:
						case 5:
						case 16:
						case 17:
						case 21:
						case 22:
						case 23:
						case 24:
						case 38:
						case 39:
						case 42: //a virer
							//case 46:
							if (mmTAGS->bufMM14[j] != NULL) fprintf(debugStreamALS, "MMTagID:%d\t(%d, %d)\n", j, *((long*)mmTAGS->bufMM14[j]), *((long*)mmTAGS->bufMM14[j] + 1));
							break;

						case 6:
						case 7:
						case 25: if (mmTAGS->bufMM14[j] != NULL) fprintf(debugStreamALS, "MMTagID:%d\t%hs\n", j, mmTAGS->bufMM14[j]);
							break;
						case 28:
						case 29: if (mmTAGS->bufMM14[j] != NULL) {
								fprintf(debugStreamALS, "MMTagID:%d\tX(*, *) (*, *)\n", j);
								for (k = 0; k < (long)mmTAGS->nbPlanes; k++) fprintf(debugStreamALS, "\t%d %d %d %d\n", *((long*)mmTAGS->bufMM14[j] + k * 4), *((long*)mmTAGS->bufMM14[j] + k * 4 + 1), *((long*)mmTAGS->bufMM14[j] + k * 4 + 2), *((long*)mmTAGS->bufMM14[j] + k * 4 + 3));
							}
							break;
						case 40:
						case 46: if (mmTAGS->bufMM14[j] != NULL) {
								fprintf(debugStreamALS, "MMTagID:%d\tz(*, *)\n", j);
								for (k = 0; k < (long)mmTAGS->nbPlanes; k++) fprintf(debugStreamALS, "\t%d %d\n", *((long*)mmTAGS->bufMM14[j] + k * 2), *((long*)mmTAGS->bufMM14[j] + k * 2 + 1));
							}
							break;
						case 41: if (mmTAGS->bufMM14[j] != NULL) {
								fprintf(debugStreamALS, "MMTagID:%d\tz valid\n", j);
								for (k = 0; k < (long)mmTAGS->nbPlanes; k++) fprintf(debugStreamALS, "\t%d\n", *((long*)mmTAGS->bufMM14[j] + k));
							}
							break;
						//case 42:
						case 43:
						case 44:
						case 45: if (mmTAGS->bufMM14[j] != NULL) fprintf(debugStreamALS, "MMTagID:%d\t%d\n", j, *((long*)mmTAGS->bufMM14[j]));
							break;
					}
				break;
			case MMUIC2: fprintf(debugStreamALS, "MMUIC2 Z distance (*/*), creation date *, time *, mod date *, time *\n");
				for (j = 0; j < (long)tags[i]->count; j++) fprintf(debugStreamALS, "\t(%d, %d), %d, %d, %d, %d\n", mmTAGS->bufMM2[j * 6], mmTAGS->bufMM2[j * 6 + 1], mmTAGS->bufMM2[j * 6 + 2], mmTAGS->bufMM2[j * 6 + 3], mmTAGS->bufMM2[j * 6 + 4], mmTAGS->bufMM2[j * 6 + 5]);
				break;
			case MMUIC3: fprintf(debugStreamALS, "MMUIC3 wavelength (*/*)\n");
				for (j = 0; j < (long)tags[i]->count; j++) fprintf(debugStreamALS, "\t(%d, %d)\n", mmTAGS->bufMM3[j * 2], mmTAGS->bufMM3[j * 2 + 1]);
				break;
			default: if (tags[i]->count * sizeOfTIFFType(tags[i]->type) > 4) fprintf(debugStreamALS, "offset value: 0x%x\n", tags[i]->val);
				else fprintf(debugStreamALS, "value:\n");

				for (j = 0; j < (long)tags[i]->count; j++) {
					if (tags[i]->valBuf != NULL) {
						switch (tags[i]->type) {
							case TIFFTYPE_BYTE_NO: val_TIFFTYPE_BYTE = (TIFFTYPE_BYTE*)tags[i]->valBuf;
								fprintf(debugStreamALS, "%u\t", val_TIFFTYPE_BYTE[j]);
								break;
							case TIFFTYPE_ASCII_NO: val_TIFFTYPE_ASCII = (TIFFTYPE_ASCII*)tags[i]->valBuf;
								if ((val_TIFFTYPE_ASCII[j] != 0x0D) && (val_TIFFTYPE_ASCII[j] != 0x00)) fprintf(debugStreamALS, "%c", val_TIFFTYPE_ASCII[j]);
								break;
							case TIFFTYPE_SHORT_NO: val_TIFFTYPE_SHORT = (TIFFTYPE_SHORT*)tags[i]->valBuf;
								fprintf(debugStreamALS, "%u\t", val_TIFFTYPE_SHORT[j]);
								break;
							case TIFFTYPE_LONG_NO: val_TIFFTYPE_LONG = (TIFFTYPE_LONG*)(tags[i]->valBuf);
								fprintf(debugStreamALS, "%u\t", val_TIFFTYPE_LONG[j]);
								break;
							case TIFFTYPE_RATIONAL_NO: val_TIFFTYPE_RATIONAL = (TIFFTYPE_RATIONAL*)tags[i]->valBuf;
								fprintf(debugStreamALS, "%u/%u\t", val_TIFFTYPE_RATIONAL[j].numerator, val_TIFFTYPE_RATIONAL[j].denominator);
								break;
							case TIFFTYPE_SBYTE_NO: val_TIFFTYPE_SBYTE = (TIFFTYPE_SBYTE*)tags[i]->valBuf;
								fprintf(debugStreamALS, "%d\t", val_TIFFTYPE_SBYTE[j]);
								break;
							case TIFFTYPE_UNDEFINED_NO: val_TIFFTYPE_UNDEFINED = (TIFFTYPE_UNDEFINED*)tags[i]->valBuf;
								fprintf(debugStreamALS, "0x%x\t", val_TIFFTYPE_UNDEFINED[j]);
								break;
							case TIFFTYPE_SSHORT_NO: val_TIFFTYPE_SSHORT = (TIFFTYPE_SSHORT*)tags[i]->valBuf;
								fprintf(debugStreamALS, "%d\t", val_TIFFTYPE_SSHORT[j]);
								break;
							case TIFFTYPE_SLONG_NO: val_TIFFTYPE_SLONG = (TIFFTYPE_SLONG*)tags[i]->valBuf;
								fprintf(debugStreamALS, "%d\t", val_TIFFTYPE_SLONG[j]);
								break;
							case TIFFTYPE_SRATIONAL_NO: val_TIFFTYPE_SRATIONAL = (TIFFTYPE_SRATIONAL*)tags[i]->valBuf;
								fprintf(debugStreamALS, "%d/%d\t", val_TIFFTYPE_SRATIONAL[j].numerator, val_TIFFTYPE_SRATIONAL[j].denominator);
								break;
							case TIFFTYPE_FLOAT_NO: val_TIFFTYPE_FLOAT = (TIFFTYPE_FLOAT*)tags[i]->valBuf;
								fprintf(debugStreamALS, "%lf\t", val_TIFFTYPE_FLOAT[j]);
								break;
							case TIFFTYPE_DOUBLE_NO: val_TIFFTYPE_DOUBLE = (TIFFTYPE_DOUBLE*)tags[i]->valBuf;
								fprintf(debugStreamALS, "%lf\t", val_TIFFTYPE_DOUBLE[j]);
								break;
						}
					}
				}
		}
		fprintf(debugStreamALS, "\n");
	}
	closeDebugStreamALS();
}


bool IFD::isMMTag(TAG* _tag)
{
	if ((_tag->tag == MMUIC1) || (_tag->tag == MMUIC2) || (_tag->tag == MMUIC3) || (_tag->tag == MMUIC4)) return (true);
	else return (false);
}


/*******************************************************************************************************************
								Class MMTAGS as specified in the UIC file format specification
*******************************************************************************************************************/
MMTAGS::MMTAGS()
{
	bufMM14 = (void**)malloc(48 * sizeof(void*));
	for (long j = 0; j < 48; j++) bufMM14[j] = NULL;
	nbPlanes = 0;
	bufMM2 = NULL;
	bufMM3 = NULL;
}

MMTAGS::~MMTAGS()
{
	for (long j = 0; j < 48; j++) free(bufMM14[j]);
	free(bufMM14);
	free(bufMM2);
	free(bufMM3);
}


/*******************************************************************************************************************
								Class TIFFImage
*******************************************************************************************************************/
TIFFImage::TIFFImage()
{
	ifd = NULL;
	in = NULL;
	firstPlaneOffset = 0;
	buffer = NULL;
}

TIFFImage::TIFFImage(unsigned long _nbPlanes, unsigned long _width, unsigned long _height, unsigned long _depth, unsigned long _samplesPerPixel)
{
	ifd = NULL;
	in = NULL;
	firstPlaneOffset = 8;
	buffer = NULL;

	nbPlanes = _nbPlanes;
	width = _width;
	height = _height;
	depth = _depth;
	samplesPerPixel = _samplesPerPixel;
	ifd = new IFD(_nbPlanes, _width, _height, _depth, _samplesPerPixel);
}

TIFFImage::~TIFFImage()
{
	delete ifd;
	free(buffer);
}

int TIFFImage::readImageInfo(char* _fileName)
{
	unsigned short header;
	unsigned long offsetIFD;
	unsigned long stripOffsets;
	TAG* tmpTag;
	strcpy(fileName, _fileName); //wcscpy
	//Reading the hearder file
	in = fopen(fileName, "rb");
	if (in == NULL) { return (ERR_READ_STK); }

	//going to the begining of the file
	if (fseek(in, 0,SEEK_SET) != 0) {
		fclose(in);
		return (ERR_READ_STK);
	}

	//Big or little endian
	if (fread(&header, sizeof(short), 1, in) != 1) {
		fclose(in);
		openDebugStreamALS();
		fprintf(debugStreamALS, "Error %d %hs\n",__LINE__,__FILE__);
		closeDebugStreamALS();
		return (ERR_READ_STK);
	}
	if (header != 0x4949) {
		openDebugStreamALS();
		fprintf(debugStreamALS, "Error file should be little endian tif! %d %hs\n",__LINE__,__FILE__);
		closeDebugStreamALS();
		fclose(in);
		return (ERR_BAD_FORMAT);
	}

	if (fread(&header, sizeof(short), 1, in) != 1) {
		fclose(in);
		openDebugStreamALS();
		fprintf(debugStreamALS, "Error %d %hs\n",__LINE__,__FILE__);
		closeDebugStreamALS();
		return (ERR_READ_STK);
	} //number 42 identifies the file as a TIFF file.

	if (fread(&offsetIFD, sizeof(long), 1, in) != 1) {
		fclose(in);
		openDebugStreamALS();
		fprintf(debugStreamALS, "Error %d %hs\n",__LINE__,__FILE__);
		closeDebugStreamALS();
		return (ERR_READ_STK);
	} //The offset (in bytes) of the first IFD.

	if ((ifd = new IFD(offsetIFD, in)) == NULL) {
		openDebugStreamALS();
		fprintf(debugStreamALS, "Error unable to read the IFD! %d %hs\n",__LINE__,__FILE__);
		closeDebugStreamALS();
		fclose(in);
		return (ERR_BAD_FORMAT);
	}
	if ((tmpTag = ifd->getTag(TIFFTAG_IMAGEWIDTH, width)) == NULL) {
		openDebugStreamALS();
		fprintf(debugStreamALS, "Error unable to read the width! %d %hs\n",__LINE__,__FILE__);
		closeDebugStreamALS();
		fclose(in);
		return (ERR_BAD_FORMAT);
	}
	if ((tmpTag = ifd->getTag(TIFFTAG_IMAGELENGTH, height)) == NULL) {
		openDebugStreamALS();
		fprintf(debugStreamALS, "Error unable to read the height! %d %hs\n",__LINE__,__FILE__);
		closeDebugStreamALS();
		fclose(in);
		return (ERR_BAD_FORMAT);
	}
	if ((tmpTag = ifd->getTag(TIFFTAG_BITSPERSAMPLE, depth)) == NULL) {
		openDebugStreamALS();
		fprintf(debugStreamALS, "Error unable to read the depth! %d %hs\n",__LINE__,__FILE__);
		closeDebugStreamALS();
		fclose(in);
		return (ERR_BAD_FORMAT);
	}
	if (ifd->getTag(TIFFTAG_SAMPLESPERPIXEL, samplesPerPixel) == NULL) { samplesPerPixel = 1; }
	if ((tmpTag = ifd->getTag(TIFFTAG_STRIPOFFSETS, stripOffsets)) == NULL) {
		openDebugStreamALS();
		fprintf(debugStreamALS, "Error unable to read the stripoffsets! %d %hs\n",__LINE__,__FILE__);
		closeDebugStreamALS();
		fclose(in);
		return (ERR_BAD_FORMAT);
	}

	if (tmpTag->count > 1) {
		fseek(in, stripOffsets,SEEK_SET);
		if (!(fread(&firstPlaneOffset, sizeof(unsigned int), 1, in))) {
			openDebugStreamALS();
			fprintf(debugStreamALS, "Error couldn'read plane offset at %d %hs!\n",__LINE__,__FILE__);
			closeDebugStreamALS();
			fclose(in);
			return (ERR_BAD_FORMAT);
		}
	}
	else {
		//only one strip=> tag value
		firstPlaneOffset = stripOffsets;
	}

	if (ifd->mmTAGS != NULL) {
		nbPlanes = ifd->mmTAGS->nbPlanes;
		if (nbPlanes == 0) nbPlanes = 1;
	}
	else nbPlanes = 1;
	sizeOfData = width * height * samplesPerPixel * nbPlanes * depth / 8;

	fclose(in);
	return (RET_OK);
}

unsigned long TIFFImage::getNbPlanes() { return (nbPlanes); }
unsigned long TIFFImage::getWidth() { return (width); }
unsigned long TIFFImage::getHeight() { return (height); }
unsigned long TIFFImage::getDepth() { return (depth); }
unsigned long TIFFImage::getSamplesPerPixel() { return (samplesPerPixel); }
unsigned long TIFFImage::getFirstPlaneOffset() { return (firstPlaneOffset); }

unsigned long TIFFImage::getSizeOfData() { return (sizeOfData); }
void* TIFFImage::getBuffer() { return (buffer); }


void* TIFFImage::allocData()
{
	sizeOfData = width * height * samplesPerPixel * nbPlanes * depth / 8;
	buffer = malloc(sizeOfData);
	if (buffer == NULL) {
		openDebugStreamALS();
		fprintf(debugStreamALS, "Error could not allocate %d bytes! %d %hs\n", sizeOfData,__LINE__,__FILE__);
		closeDebugStreamALS();
		return (NULL);
	}

	return (buffer);
}


int TIFFImage::readData(void* _buffer)
{ //this buffer have to be allocate (getSizeOfData) before calling readData
	sizeOfData = width * height * samplesPerPixel * nbPlanes * depth / 8;
	FILE* in = fopen(fileName, "rb");
	if (in == NULL) { return (ERR_READ_STK); }

	if (fseek(in, firstPlaneOffset,SEEK_SET) != 0) return (ERR_READ_STK);
	int readDataCount;

	if ((readDataCount = fread(_buffer, sizeof(char), sizeOfData, in)) != (long)sizeOfData) {
		openDebugStreamALS();
		fprintf(debugStreamALS, "Error reading the data (%d!=%d bytes)! %d %ls\n", readDataCount, sizeOfData,__LINE__,__FILE__);
		fprintf(debugStreamALS, "firstPlaneOffset %d\n", firstPlaneOffset);
		closeDebugStreamALS();
		fclose(in);
		return (ERR_READ_STK);
	}
	fseek(in, 0,SEEK_END);
	fclose(in);
	in = NULL;
	return (RET_OK);
}

int TIFFImage::readData()
{
	if (buffer == NULL) if (allocData() == NULL) return (ALLOC_ERR);
	return (readData(buffer));
}

void TIFFImage::printTags() { ifd->print(); }

int TIFFImage::writeImage(char* _fileName, void* data)
{
	int ret;
	unsigned long valTmp;
	char head[] = {0x49, 0x49, 0x2a, 0x00};
	unsigned long unkIfdOffset = 0;
	div_t frac;
	TAG* stripOffsets;
	TAG* rowsPerStrip;
	TAG* stripBytesCount;
	unsigned long stripBytesCountTmp;
	unsigned int i;

#ifdef DEBUG_RW_STK
	openDebugStreamALS();
	fprintf(debugStreamALS, "writing image:%ls\n", _fileName);
	closeDebugStreamALS();
#endif

	sizeOfData = width * height * samplesPerPixel * nbPlanes * depth / 8;

	in = fopen(_fileName, "wb");
	if (in == NULL) { return (ERR_WRITE_STK); }

	fseek(in, 0,SEEK_SET); //go to the begining of the file
	//writing the header file
	if (fwrite(head, 1, 4, in) != 4) {
		fclose(in);
		return (ERR_WRITE_STK);
	}
	if (fwrite(&unkIfdOffset, sizeof(unsigned long), 1, in) != 1) {
		fclose(in);
		return (ERR_WRITE_STK);
	}

#ifdef DEBUG_RW_STK
	openDebugStreamALS();
	fprintf(debugStreamALS, "head info written\n");
	closeDebugStreamALS();
#endif

	if (fwrite(data, sizeof(char), sizeOfData, in) != sizeOfData) {
		openDebugStreamALS();
		fprintf(debugStreamALS, "Error writing the data (%d bytes) in the image %ls! %d %hs\n", sizeOfData, _fileName,__LINE__,__FILE__);
		closeDebugStreamALS();
		fclose(in);
		return (ERR_WRITE_STK);
	}

#ifdef DEBUG_RW_STK
	openDebugStreamALS();
	fprintf(debugStreamALS, "data written\n");
	closeDebugStreamALS();
#endif


	//Updates the tags stripOffsets, rowsPerStrip and stripBytesCount.
	if ((stripOffsets = ifd->getTag(TIFFTAG_STRIPOFFSETS, valTmp)) == NULL) {
		fclose(in);
		return (ERR_BAD_FORMAT);
	}
	if ((rowsPerStrip = ifd->getTag(TIFFTAG_ROWSPERSTRIP, valTmp)) == NULL) {
		fclose(in);
		return (ERR_BAD_FORMAT);
	}
	if ((stripBytesCount = ifd->getTag(TIFFTAG_STRIPBYTECOUNTS, valTmp)) == NULL) {
		fclose(in);
		return (ERR_BAD_FORMAT);
	}

	rowsPerStrip->val = (TIFFTYPE_LONG)(8192.0 / (width * samplesPerPixel * depth / 8)); //arbitrary size of strip
	rowsPerStrip->count = 1;
	frac = div((int)height, rowsPerStrip->val); // computes number of strips and how many rows in last strip(if any)
	stripOffsets->count = (frac.rem == 0) ? frac.quot : (frac.quot + 1);
	stripBytesCount->count = stripOffsets->count;
	stripBytesCountTmp = width * samplesPerPixel * rowsPerStrip->val * depth / 8;

	fseek(in, 0,SEEK_END);
	//realloc if the image is new and not reed.
	if (stripOffsets->count * sizeof(TIFFTYPE_LONG) > 4) {
		if (stripOffsets->val == 0) {
			stripBytesCount->valBuf = malloc(stripOffsets->count * sizeof(TIFFTYPE_LONG));
			stripOffsets->valBuf = malloc(stripOffsets->count * sizeof(TIFFTYPE_LONG));
		}
		else {
			stripBytesCount->valBuf = realloc(stripBytesCount->valBuf, stripOffsets->count * sizeof(TIFFTYPE_LONG));
			stripOffsets->valBuf = realloc(stripOffsets->valBuf, stripOffsets->count * sizeof(TIFFTYPE_LONG));
		}
	}

	stripOffsets->val = ftell(in);
	stripBytesCount->val = ftell(in);

	for (i = 0; i < stripOffsets->count; i++) { //writes stripOffsets
		if (stripOffsets->type == TIFFTYPE_LONG_NO) ((TIFFTYPE_LONG*)(stripOffsets->valBuf))[i] = (TIFFTYPE_LONG)(8 + i * stripBytesCountTmp);
		else //TIFFTYPE_SHORT_NO
			((TIFFTYPE_SHORT*)(stripOffsets->valBuf))[i] = (TIFFTYPE_SHORT)(8 + i * stripBytesCountTmp);
	}

	for (i = 0; (long)i < frac.quot; i++) {
		if (stripBytesCount->type == TIFFTYPE_LONG_NO) ((TIFFTYPE_LONG*)(stripBytesCount->valBuf))[i] = (TIFFTYPE_LONG)stripBytesCountTmp;
		else //TIFFTYPE_SHORT_NO
			((TIFFTYPE_SHORT*)(stripBytesCount->valBuf))[i] = (TIFFTYPE_SHORT)stripBytesCountTmp;
	}
	if (frac.rem) { //last strip?
		if (stripBytesCount->type == TIFFTYPE_LONG_NO) ((TIFFTYPE_LONG*)(stripBytesCount->valBuf))[i] = (TIFFTYPE_LONG)(frac.rem * width * samplesPerPixel * depth / 8);
		else //TIFFTYPE_SHORT_NO
			((TIFFTYPE_SHORT*)(stripBytesCount->valBuf))[i] = (TIFFTYPE_SHORT)(frac.rem * width * samplesPerPixel * depth / 8);
	}
#ifdef DEBUG_RW_STK
	openDebugStreamALS();
	fprintf(debugStreamALS, "The tags stripOffsets, rowsPerStrip and stripBytesCount have been updated\n");
	closeDebugStreamALS();
#endif

	if ((ret = ifd->writeAllBufTags(in)) < 0) {
		fclose(in);
		return (ret);
	}

	if ((ret = ifd->writeIfd(in)) < 0) return (ret);

	return (RET_OK);
}

int TIFFImage::addToString(char* stringToAdd)
{
	TAG* tag;
	long i, j, k;
	unsigned long tmp;
	char* strTmp;
	long countStrings;
	char* tmpBuf;
	long strL;
	int ret;
	long countSlashN;
	if ((tag = ifd->getTag(TIFFTAG_IMAGEDESCRIPTION, tmp)) == NULL) {
		//return(ERR_BAD_FORMAT);
		tmpBuf = (char*)malloc(sizeof(char) * __max(2*nbPlanes, 5));
		if (5 > 2 * nbPlanes) for (i = 0; i < 5; i++) tmpBuf[i] = 0;
		else
			for (i = 0; i < (long)nbPlanes; i++) {
				tmpBuf[2 * i] = ' ';
				tmpBuf[2 * i + 1] = 0;
			}

		tag = new TAG(TIFFTAG_IMAGEDESCRIPTION, TIFFTYPE_ASCII_NO, __max(2*nbPlanes, 5), 0, tmpBuf);

		if ((ret = ifd->addNewTag(tag)) < 0) return (ret);
		//free(tmpBuf);
	}

	countStrings = 0;

	for (i = 0; i < (long)tag->count; i++) if (((char*)(tag->valBuf))[i] == 0) countStrings++;


	countSlashN = 0;

	strL = strlen(stringToAdd);
	for (i = 0; i < strL; i++) { if (stringToAdd[i] == '\n') countSlashN++; }

	strL += countSlashN;
	if ((strTmp = (char*)malloc((tag->count + strL * countStrings) * sizeof(char))) == NULL) return (ALLOC_ERR);

	for (k = 0, i = 0; i < (long)tag->count; i++, k++)
		if (((char*)(tag->valBuf))[i] == 0) {
			for (j = 0; j < strL - countSlashN; j++, k++)
				if (stringToAdd[j] == '\n') {
					strTmp[k++] = 0x0D;
					strTmp[k] = 0x0A;
				}
				else strTmp[k] = stringToAdd[j];
			strTmp[k] = 0;
		}
		else strTmp[k] = ((char*)(tag->valBuf))[i];


	tag->count += strL * countStrings;
	free(tag->valBuf);
	if ((tag->valBuf = malloc(tag->count * sizeof(char))) == NULL) return (ALLOC_ERR);
	memcpy(tag->valBuf, strTmp, tag->count * sizeof(char));
	free(strTmp);

	return (RET_OK);
}

int TIFFImage::writeImage(char* fileName) { return (writeImage(fileName, buffer)); }


/*****************************************************************************************************
EXAMPLES:

***1) Open a stack, write an other stack with the same information...

	TIFFImage i;
	int ret;
	if((ret=i.readImageInfo("c:\\stk.stk"))<0) return(ret);
	i.readData();
	i.addToString("coucou c'est moi");

	if((ret=i.writeImage("c:\\stkcpy.stk"))<0) return(ret);
	return(0);


***2) Create a new stack...
	long n=20;
	TIFFImage i=TIFFImage(100, n, 100, 16, 1);

	i.allocData();
	void* data=i.getBuffer();
	i.addToString("coucou c'est moi");
	for(long c=0;c<10000*n;c++)
		((unsigned short *)data)[c]=(unsigned short)rand();
	int ret;
	if((ret=i.writeImage("c:\\test.stk"))<0) return(ret);
	return(0);

*****************************************************************************************************/

#undef	TIFFTAG_SUBFILETYPE
#undef	TIFFTAG_OSUBFILETYPE
#undef	TIFFTAG_IMAGEWIDTH
#undef	TIFFTAG_IMAGELENGTH
#undef	TIFFTAG_BITSPERSAMPLE
#undef	TIFFTAG_COMPRESSION
#undef	TIFFTAG_PHOTOMETRIC
#undef  TIFFTAG_IMAGEDESCRIPTION
#undef	TIFFTAG_STRIPOFFSETS
#undef	TIFFTAG_SAMPLESPERPIXEL
#undef	TIFFTAG_ROWSPERSTRIP
#undef	TIFFTAG_STRIPBYTECOUNTS
#undef	TIFFTAG_XRESOLUTION
#undef	TIFFTAG_YRESOLUTION
#undef	TIFFTAG_RESOLUTIONUNIT
#undef	TIFFTAG_SOFTWARE
#undef	TIFFTAG_DATETIME
#undef	TIFFTAG_PREDICTOR

#undef TIFFTYPE_BYTE_NO
#undef TIFFTYPE_ASCII_NO
#undef TIFFTYPE_SHORT_NO
#undef TIFFTYPE_LONG_NO
#undef TIFFTYPE_RATIONAL_NO
#undef TIFFTYPE_SBYTE_NO
#undef TIFFTYPE_UNDEFINED_NO
#undef TIFFTYPE_SSHORT_NO
#undef TIFFTYPE_SLONG_NO
#undef TIFFTYPE_SRATIONAL_NO
#undef TIFFTYPE_FLOAT_NO
#undef TIFFTYPE_DOUBLE_NO

#undef MMUIC1
#undef MMUIC2
#undef MMUIC3
#undef MMUIC4
