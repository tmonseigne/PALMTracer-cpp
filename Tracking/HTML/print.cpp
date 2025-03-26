
#include "stdafx.h"


/*******************************************************************************************************************
									print.cpp to print directly in HTLM files 
*******************************************************************************************************************/

#include "print.h"
#define MAX_NB_STRING 1000
void addHTMLTag(wchar_t* buf, const wchar_t* tag, long &j){
	wcscpy(buf+j, tag);
	j+=wcslen(tag);
	//if(buf[j]=='\0')
	//	j--;
}

int printTableLinesHTML(FILE *stream, const wchar_t *format, ...){
	
	va_list v;
	long count;
	long i,j;
	wchar_t buf[MAX_NB_STRING]; 
	wchar_t bufHTML[MAX_NB_STRING];
	
	va_start(v, format);
	count=vswprintf(buf, format, v);
	if(count>MAX_NB_STRING){
		fprintf(stream, "\n\n<br><br><br>Error in \"int printHTML(FILE *file, const wchar_t *format, ...)\"<br>\nToo much charaters in the buffer<br>\n</body>\n</html>");
		return(T_ERR_HTML);
	}

	j=0;
	addHTMLTag(bufHTML, L"\t<tr>\n\t\t<td>",j);
	for(i=0;(j<MAX_NB_STRING)&&(buf[i]!='\0');i++,j++){
		switch(buf[i]){
		case '\t':
			addHTMLTag(bufHTML, L"</td>\n\t\t<td>",j);
			j--;
			break;
		default:
			bufHTML[j]=buf[i];
			break;
		}
	}
	j--;
	addHTMLTag(bufHTML, L"</td>\n\t</tr>\n",j);

	fwprintf(stream, bufHTML);
	va_end(v);
	return(RET_OK);
}

int printTableLinesColorHTML(FILE *stream, const wchar_t *bGColor, const wchar_t *color, const wchar_t *format, ...){ //back ground color and color font
	va_list v;
	long count;
	long i,j;
	wchar_t buf[MAX_NB_STRING]; 
	wchar_t bufHTML[MAX_NB_STRING];
	
	va_start(v, format);
	count=vswprintf(buf, format, v);
	if(count>MAX_NB_STRING){
		fprintf(stream, "\n\n<br><br><br>Error in \"int printHTML(FILE *file, const wchar_t *format, ...)\"<br>\nToo much charaters in the buffer<br>\n</body>\n</html>");
		return(T_ERR_HTML);
	}

	j=0;
	if(bGColor!=NULL){
		addHTMLTag(bufHTML, L"\t<tr bgcolor=\"",j);
		addHTMLTag(bufHTML, bGColor,j);
		addHTMLTag(bufHTML, L"\">\n",j);
	}else{
		addHTMLTag(bufHTML, L"\t<tr>\n",j);
	}

	addHTMLTag(bufHTML, L"\t\t<td>",j);
	if(color!=NULL){
		addHTMLTag(bufHTML, L"<font color=\"",j);
		addHTMLTag(bufHTML, color,j);
		addHTMLTag(bufHTML, L"\">",j);
	}


	for(i=0;(j<MAX_NB_STRING)&&(buf[i]!='\0');i++,j++){
		switch(buf[i]){
		case '\t':
			if(color!=NULL){
				addHTMLTag(bufHTML, L"</font></td><td>",j);
				addHTMLTag(bufHTML, L"<font color=\"",j);
				addHTMLTag(bufHTML, color,j);
				addHTMLTag(bufHTML, L"\">",j);
			}else{
				addHTMLTag(bufHTML, L"</td><td>",j);
			}
			j--;
			break;
		default:
			bufHTML[j]=buf[i];
			break;
		}
	}
	
	if(color!=NULL)
		addHTMLTag(bufHTML, L"</font>",j);
	addHTMLTag(bufHTML, L"</td>\n\t</tr>\n",j);

	fwprintf(stream, bufHTML);
	va_end(v);
	return(RET_OK);
}

int openRowColorHTML(FILE *stream,  const wchar_t *bGColor){ 
	if(bGColor!=NULL){
		fwprintf(stream, L"<table>\n\t<tr bgcolor=\"%ls\">\n\t\t",bGColor);
	}else{
		fwprintf(stream, L"<table>\n\t<tr>\n\t\t");
	}
	return(RET_OK);
}

int addCellColorHTML(FILE *stream,  const wchar_t *color, const wchar_t *st){ 
	if(color!=NULL){
		fwprintf(stream, L"<td><font color=\"%ls\">%ls</font></td>",color,st);
	}else{
		fwprintf(stream, L"<td>%ls</td>",st);
	}
	return(RET_OK);
}

int closeRowColorHTML(FILE *stream){ 
	fwprintf(stream, L"\n\t</tr>\n</table>\n");
	return(RET_OK);
}

int printHTML(FILE *stream, const wchar_t *format, ...){

	va_list v;
	long count;
	long i,j;
	wchar_t buf[MAX_NB_STRING]; 
	wchar_t bufHTML[MAX_NB_STRING];
	bool errorWord=false;
	
	
	va_start(v, format);
	count=vswprintf(buf, format, v);
	if(count>MAX_NB_STRING){
		fwprintf(stream, L"\n\n<br><br><br>Error in \"int printHTML(FILE *file, const wchar_t *format, ...)\"<br>\nToo much charaters in the buffer<br>\n</body>\n</html>");
		return(T_ERR_HTML);
	}


	j=0;

	if(_wcsnicmp(format, L"error", 5)==0)
		errorWord=true;

	if(errorWord)
		addHTMLTag(bufHTML, L"<font color=red>",j);

			
	
	for(i=0;(j<MAX_NB_STRING)&&(buf[i]!='\0');i++,j++){
		switch(buf[i]){
		case '\n':
			addHTMLTag(bufHTML, L"<br>\n",j);
			break;
		default:
			bufHTML[j]=buf[i];
			break;
		}
	}
	if(errorWord){
		j--;
		addHTMLTag(bufHTML, L"</font>\n",j);
	}
	if(bufHTML[j-1]!=L'\0')
		bufHTML[j]=L'\0';
	
	fwprintf(stream, bufHTML);
	va_end(v);



	return(RET_OK);
}


FILE* openHTML(PublicParameters *publicParameters, const wchar_t* title){
	FILE* stream;
	if( (stream=_wfopen(publicParameters->fileRes,L"w"))==NULL){
		wprintf( L"The file %ls can not be opened\n", publicParameters->fileRes);
		return(NULL);
	}

	//HTML 4.01
	//fprintf(stream,"<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01//EN\" \"http://www.w3.org/TR/html4/strict.dtd\">\n");
	fprintf(stream,"<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\" \"http://www.w3.org/TR/html4/loose.dtd\">\n");

	//header
	fprintf(stream,"<html>\n<head>\n");
	fprintf(stream,"\t<META NAME=\"Author\" LANG=\"en\" CONTENT=\"\">\n");
	fprintf(stream,"\t<META NAME=\"Publisher\" CONTENT=\"\">\n");
	fprintf(stream,"\t<META NAME=\"Reply-to\" CONTENT=\"\">\n");
	fprintf(stream,"\t<META NAME=\"Indentifier-URL\" CONTENT=\"\">\n");
	fprintf(stream,"\t<title>\n\t\t%ls\n\t</title>\n\n", title);
	fprintf(stream,"</head>\n<body>\n\n");
	fprintf(stream,"<h1><a name=\"top\">%ls</a></h1>\n\n", title);

	fprintf(stream,"<p><b>This file has been generated by \"%ls\" from software</b><br>\n", publicParameters->fileArg0);
#ifdef DEBUG_AT
	fprintf(stream,"<b>Debug Mode Result</b><br>\n");
#endif
	char buffer [9];
	_strdate( buffer );
	fprintf(stream,"The current date is %s (mm/dd/yy) and ", buffer );
	_strtime( buffer );
	fprintf(stream,"the current time is %s <br>\n", buffer );

	return(stream);
}

void closeHTML(FILE* stream){
	fprintf(stream,"\n</body>\n</html>\n\n");
	fclose(stream);
}

#undef MAX_NB_STRING 
