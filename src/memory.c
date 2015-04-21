#include <stdio.h>
#include <string.h>
#include "memory.h"


void* _safeMalloc(const char* filename,int line, size_t size)
	{
	void* ptr= malloc(size);
	if(ptr==NULL)
		{
		fprintf(stderr,"[%s:%d]out of memory",filename,line);
		exit(EXIT_FAILURE);
		}
	return ptr;
	}
	
void* _safeCalloc(const char* filename,int line, size_t nmemb,size_t size)
	{
	void* ptr= calloc(nmemb,size);
	if(ptr==NULL)
		{
		fprintf(stderr,"[%s:%d]out of memory",filename,line);
		exit(EXIT_FAILURE);
		}
	return ptr;
	}

void* _safeRealloc(const char* filename,int line, void* ptr0,size_t size)
	{
	void* ptr= realloc(ptr0,size);
	if(ptr==NULL)
		{
		fprintf(stderr,"[%s:%d]out of memory",filename,line);
		exit(EXIT_FAILURE);
		}
	return ptr;
	}
	
char* _safeStrdup(const char* filename,int line, const char* s)
	{
	char* ptr=strdup(s);
	if(ptr==NULL)
		{
		fprintf(stderr,"[%s:%d]out of memory",filename,line);
		exit(EXIT_FAILURE);
		}
	return ptr;
	}

FILE* _safeFOpen(const char* filename,int line, const char* file0,const char* mode)
	{
	FILE* file=fopen(file0,mode);
	if ( file== NULL) 
 		{
		fprintf(stderr,"Cannot open file [%s:%d]\n",filename, line);
		exit(EXIT_FAILURE);
	      	}    
	return file; 
	}
	

	
