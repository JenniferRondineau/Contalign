#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* main */

int main(int argc, char** argv)
	{
	FILE* input=stdin;
	FILE* output=stdout;
	int c;
	while((c=fgetc(input))!=-1)
	        {
	        fprintf(output,"%c",c);
	        }
	
        return 0;
	}

