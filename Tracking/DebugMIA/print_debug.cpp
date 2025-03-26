#include "stdafx.h"

/*******************************************************************************************************************
								Tools for printing in a file for debugging
*******************************************************************************************************************/


#include "print_debug.h"

//char debugLine[256];

FILE *debugStreamALS;



void openDebugStreamALS(void){
	debugStreamALS = _wfopen( L"c:\\debug_ALS.out", L"a" );
}

void emptyALS(void){
	debugStreamALS = _wfopen( L"c:\\debug_ALS.out", L"w" );
	fclose( debugStreamALS );
}


void closeDebugStreamALS(void){
   fclose( debugStreamALS );
}