/*

The MIT License (MIT)

Copyright (c) 2015 Jennifer Rondineau

http://samtools.sourceforge.net/
Authors: Heng Li, Bob Handsaker, Jue Ruan, Colin Hercus, Petr Danecek

https://github.com/lh3/bwa/blob/master/example.c
Authors: Heng Li's

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
 	
 	
 	
#define VERIFY_NOT_NULL(POINTER) do{if(POINTER==NULL) { fprintf(stderr,"Memory Alloc failed File %s Line%d.\n",__FILE__,__LINE__); exit(EXIT_FAILURE);}}while(0)
	
	
									
#define DUMP_FASTQ if ( file_fastq != NULL) {\
		   for (i=0; i < fastq_count; i++)\
			{\
			fprintf(file_fastq, "@%s%s\n%s\n+\n%s\n", fastq[i].qname, fastq[i].read, fastq[i].seq, fastq[i].qual);\
			}\
		   } else \
			{\
			fprintf(stderr,"Unsaved fasta file.\n");\
		   }

/*version of this code */ 
static void version () {printf(" Version : v.1 \n");}

/* map unmapped reads against a reference of contaminants */
Contaminants *align( Fastq* fastq, Contaminants* contaminant, int fastq_count, int *count_contaminants, bwaidx_t *idx);



static int compareGroup(const void *g1, const void *g2)
{
	Group *group1 = (Group *) g1;
	Group *group2 = (Group *) g2;
	return strcmp (group1->rgId, group2->rgId);  //sort Read Group in alphabetical order
}


static int compareContaminant(const void *c1, const void *c2)
{
	Contaminants *contaminant1 = (Contaminants *) c1;
	Contaminants *contaminant2 = (Contaminants *) c2;
	return contaminant2->contaminants_count - contaminant1->contaminants_count; //sort contaminants in descending order
}

