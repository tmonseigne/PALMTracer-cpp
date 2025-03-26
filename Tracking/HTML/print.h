/*******************************************************************************************************************
									print.cpp to print directly in HTLM files 
*******************************************************************************************************************/

#ifndef __HTML_PRINT__
#define __HTML_PRINT__
class PublicParameters;
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "..\kernel\general.h"
//#include <mbstring.h>
#include "..\generalParam\noError.h"
int printHTML(FILE* file, const wchar_t* format, ...);
int printTableLinesHTML(FILE* stream, const wchar_t* format, ...);
int printTableLinesColorHTML(FILE* stream, const wchar_t* bGColor, const wchar_t* color, const wchar_t* format, ...);
FILE* openHTML(PublicParameters* publicParameters, const wchar_t* title = L"");
void closeHTML(FILE* stream);
int openRowColorHTML(FILE* stream, const wchar_t* bGColor);
int addCellColorHTML(FILE* stream, const wchar_t* color, const wchar_t* st);
int closeRowColorHTML(FILE* stream);
#endif //__HTML_PRINT__
