#pragma once


#define WIDTH_PORTION 16 // 32
#define HEIGHT_PORTION 16 // 32

#define coeffH0  3. / 8.
#define coeffH1  1. / 4.
#define coeffH2  1. / 16.

#define FIT_BLOCK_SIZE /*32*/64 /*128*/

#define _coeffH0 2. / 3.
#define _coeffH1 1. / 6.

#define THREAD_PORTION 32
#define BLOCK_SIZE 256
#define MAX_NB_STRUCTURES 150/*400*//*800*/// nombre maximum des structures par image
#define MAX_PIXELS_BY_STRUCTURE 50  // nombre maximum de pixels par structure

#define PI 3.1415926535897932384626433832795

#define EPSILON 0.0000001

#define MAX_BLOCKS 65535

//__device__ unsigned short WIDTH/* = 128*/;
//__device__ unsigned short HEIGHT/* = 128*/;
