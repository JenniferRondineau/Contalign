#include <stdio.h>
#include <string.h>
#include "debug.h"
#include "memory.h"


void* _safeMalloc(const char* filename,int line, size_t size)
	{
	void* ptr= malloc(size);
	if(ptr==NULL)
		{
		FATAL("[%s:%d]out of memory",filename,line);
		}
	return ptr;
	}
	
void* _safeCalloc(const char* filename,int line, size_t nmemb,size_t size)
	{
	void* ptr= calloc(nmemb,size);
	if(ptr==NULL)
		{
		FATAL("[%s:%d]out of memory",filename,line);
		}
	return ptr;
	}

void* _safeRealloc(const char* filename,int line, void* ptr0,size_t size)
	{
	void* ptr= realloc(ptr0,size);
	if(ptr==NULL)
		{
		FATAL("[%s:%d]out of memory",filename,line);
		}
	return ptr;
	}
	
char* _safeStrdup(const char* filename,int line, const char* s)
	{
	char* ptr=strdup(s);
	if(ptr==NULL)
		{
		FATAL("[%s:%d]out of memory",filename,line);
		}
	return ptr;
	}

FILE* _safeFOpen(const char* filename,int line, const char* file0,const char* mode)
	{
	FILE* file=NULL;
	LOG("Opening %s with mode=\"%s\"\n",file0,mode);
	file = fopen(file0,mode);
	if ( file== NULL) 
 		{
		FATAL("Cannot open file [%s:%d]\n",filename, line);
	      	}    
	return file; 
	}
	

	
