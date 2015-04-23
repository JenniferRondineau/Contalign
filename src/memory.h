/*

The MIT License (MIT)

Copyright (c) 2015 Jennifer Rondineau


Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

*/
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
