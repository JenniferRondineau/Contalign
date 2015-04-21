#ifndef MEMORY_H
#define MEMORY_H
#include <stdlib.h>

void* _safeMalloc(const char* filename,int line, size_t size);
void* _safeCalloc(const char* filename,int line, size_t nmemb,size_t size);
void* _safeRealloc(const char* filename,int line, void* ptr0,size_t size);
char* _safeStrdup(const char* filename,int line, const char* s);
FILE* _safeFOpen(const char* filename,int line, const char* file0,const char* mode);

#define safeMalloc(S) _safeMalloc(__FILE__,__LINE__,S)
#define safeCalloc(N,S) _safeCalloc(__FILE__,__LINE__,N,S)
#define safeRealloc(P,S) _safeRealloc(__FILE__,__LINE__,P,S)
#define safeStrdup(C) _safeStrdup(__FILE__,__LINE__,C)
#define safeFOpen(F,M) _safeFOpen(__FILE__,__LINE__,F,M)


#endif
