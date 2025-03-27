/*******************************************************************************************************************
								Tools for printing in a file for debugging
*******************************************************************************************************************/


#include "print_debug.h"

//char debugLine[256];

FILE* debugStream;


void openDebugStream(void) { debugStream = fopen("c:\\debug.out", "a"); }


void closeDebugStream(void) { fclose(debugStream); }
