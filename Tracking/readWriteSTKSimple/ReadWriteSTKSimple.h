


/*******************************************************************************************************************
								Read and write TIFF and MM STK for 8 16 & 32bits unsigned integer
*******************************************************************************************************************/


#ifndef __RW_STK__
#define __RW_STK__

		 




int readInfoSTK(char* fileName, unsigned long & depth, unsigned long & samplesPerPixel, unsigned long & nbPlanes, unsigned long & width, unsigned long & height, unsigned long & sizeOfData);
int readSTK(char* fileName, unsigned char *image);
int writeSTK(char* fileNameDest, char* fileNameSource, unsigned char *image);


int readSTK32(char* fileName, unsigned long *image);
int readSTK16(char* fileName, unsigned short *image);
int readSTK8(char* fileName, unsigned char *image);
int writeSTK32(char* fileNameDest, char* fileNameSource, unsigned long *image);
int writeSTK16(char* fileNameDest, char* fileNameSource, unsigned short *image);
int writeSTK8(char* fileNameDest, char* fileNameSource, unsigned char *image);

#endif  //__RW_STK__ 
