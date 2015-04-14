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
#ifndef COUNT_ALIGN_H
#define COUNT_ALIGN_H
#include "debug.h"

/* structure of the program parameters */  
typedef struct Contalign
	{
	samFile *fp; /* initial BAM file */
	samfile_t *out_file; /* final BAM file */
	const char* filename_in; /* name of initial BAM file */
	char* filename_out; /* name of final BAM file */
	FILE* file_fastq; /* Fastq file with unmapped read */
	char* filename_fastq; /* name of fastq file */
	FILE* file; /* report file */
	char* output_report; /* name of output report file */ 
	FILE* full_report; /* full report file */
	char* namefull_report; /* name of full report file */
	bwaidx_t *idx; /* index file */
	char* ref; /* name of reference of contaminant */
	}Contalign,*ContalignPtr;

/* contaminant structure */
typedef struct Contaminants
	{
	char* c_name; /* name of contaminant ex: E.coli */
	float contaminants_count; 
	}Contaminants;
	
/*Fastq structure */
typedef struct Fastq
	{
	char* qname; /* Name of read */ 
	char* read; /* "/1" forward reads ; "/2" reverse reads */
	char* seq; /* sequence */ 
	char* qual; /* the quality values for the sequence */
	Contaminants* contaminant;
	}Fastq;

/* Sample structure */ 
typedef struct SampleName
 	{
 	char* sample_name; /* Name of sample must be unique */
 	long unMap; /* number of unmapped reads by sample*/
 	Fastq* fastq; /*fastq associated with this sample */
 	}Sample;
 
/* Group structure */ 	
typedef struct Group
 	{
 	long unMap; /* number of unmapped reads by read group */
 	char* rgId; /* Name of read group must be unique */
 	Sample* sample; /* Sample associated with this read group */ 
 	}Group;
 	
 	


/* map unmapped reads against a reference of contaminants */
Contaminants *align( Fastq* fastq, Contaminants* contaminant, int fastq_count, int *count_contaminants, bwaidx_t *idx);


#define COUNTALIGN_VERSION "0.1"


#endif



