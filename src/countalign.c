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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sam.h"
#include "sam_header.h"
#include "debug.h"
#include <errno.h>  


#define VERIFY_NOT_NULL(POINTER) do{if(POINTER==NULL) { fprintf(stderr,"Memory Alloc failed File %s Line%d.\n",__FILE__,__LINE__); exit(EXIT_FAILURE);}}while(0)

typedef struct SampleName
 	{
 	char* sample_name; /* must be unique */
 	long unMap;
 	}Sample;
 	
typedef struct Group
 	{
 	char* rgId;/* must be unique */
 	Sample* sample;
 	}Group;


static void usage ()
{ 
/* cf isatty */
fprintf(stderr,
"\nUsage : cat file.bam | ./countalign -o file.txt > file.bam \n\
\n\
Description : \n\
\n\
This program reads BAM file and looks for unmapped reads. It released a report containing the total number of reads and the number of unmapped reads by sample. \n\
\n\
Options :\n\
\n\
         -o \033[4mFILE\033[0m   Name of the output file containing the report of unmapped reads.\n\
\n\
         -O \033[4mFILE\033[0m   Name of the output BAM file.\n\
\n\
         -h        Output help and exit immediately. \n\
\n\
         -v        Output version and exit immediately. \n\n");
}



static void version ()
{
printf(" Version : v.1 \n");
}





static int compareGroup(const void *g1, const void *g2)
{
	Group *group1 = (Group *) g1;
	Group *group2 = (Group *) g2;
	return strcmp (group1->rgId, group2->rgId);
}




int main(int argc, char** argv)
	{
	FILE* file=NULL;	
	char* filename_out =NULL; 
	int c;
	char* output_report; 
	bam_hdr_t *header=NULL;
	bam1_t *b = NULL;
	b = bam_init1();
	long nReads=0;
	int n=0, i=0;
	char *rgId=NULL, *sampleName=NULL;	
	Sample* samples=NULL;
	int sample_count=0;
	Group* group=NULL;
	int group_count=0;
	uint8_t *search = NULL; 
	char* Groupsearch = NULL;

	
	// get commande line options
	 while ((c = getopt(argc, argv,"o:O:hv")) != -1) {
	 switch(c)
		 {
		 case 'h': 
		 	{
		 	usage();
		 	return EXIT_SUCCESS;
		 	break;
		 	}
		 case 'v':
		 	{
		 	version();
		 	return EXIT_SUCCESS;
		 	break;
		 	}
		 
		 case 'o': output_report = strdup(optarg); break;
		 case 'O' :
		 	{
		 	filename_out = strdup(optarg);
		 	if( filename_out == NULL )
		 		{
		 		fprintf(stderr,"Cannot alloc memory.\n");
		 		return EXIT_FAILURE;
		 		}
		 	break;
		 	}
		case '?': 
			{
			 fprintf(stderr, "ERROR: option -%c is undefined\n", optopt);
		 	 return EXIT_FAILURE;
		 	 break;
		 	 }
		case ':': 
			{
			fprintf(stderr, "ERROR: option -%c requires an argument\n",optopt);
			return EXIT_FAILURE;
		 	break;
		 	}
		default :
			{
			usage(); 
			return EXIT_FAILURE;
			break;
			}
		 }
	

	 } // si ./ countalign aljsd 

	
	
	samFile *fp = sam_open( "-" ); 
	
		if (fp == NULL) 
			{
        		fprintf(stderr,"Cannot read file. %s.\n",strerror(errno));
        		return EXIT_FAILURE;
			}

	header = sam_hdr_read(fp);

		if( header == NULL)
			{
			fprintf(stderr, "Cannot read header \n");
			return EXIT_FAILURE;
			}


	samfile_t *out_file = samopen((filename_out!=NULL? filename_out : "-"), "wb", header ); 
	
		if (out_file == NULL) 
			{
			fprintf(stderr,"Failed to open output file . %s.\n",strerror(errno));
			return EXIT_FAILURE;
			} 



	void *iter = sam_header_parse2(header->text);
	while ( (iter = sam_header2key_val(iter, "RG","ID","SM",&rgId,&sampleName) )!=NULL)
		{
			Sample* searchSample = NULL;
			
			/* sample must be unique */
			
			for (i=0;i<sample_count;++i)
				{ 
				if (strcmp (samples[i].sample_name, sampleName) == 0)
					{
					searchSample = &samples[i];
					break;
					}
				}
				
			if (searchSample == NULL ) 
				{ 					
				samples = (Sample*)realloc(samples,sizeof(Sample)*(sample_count+1));
				VERIFY_NOT_NULL(samples);
				samples[sample_count].sample_name = strdup(sampleName);
				VERIFY_NOT_NULL(samples[sample_count].sample_name);
				samples[sample_count].unMap=0;
				
				searchSample=&samples[sample_count];
				
				sample_count++;
				}

			 /* group must be unique */							
			
			for (i=0;i<group_count;++i)
			{ 
				if (strcmp (group[i].rgId,rgId) == 0)
				{
				fprintf(stderr,"DUP GROUP");
				return EXIT_FAILURE;
				}
			}
					
								
			group = realloc(group,sizeof(Group)*(group_count+1));	
			VERIFY_NOT_NULL(group);	
			group[group_count].rgId = strdup(rgId);
			VERIFY_NOT_NULL(group[group_count].rgId);
			
			group[group_count].sample = searchSample;  
			group_count++;							
			

			
		}
		

	if ((group_count == 0)  && (sample_count ==0))
		{
		fprintf(stderr,"Error, missing SAM header or @RG tag\n");
		return EXIT_FAILURE;
		}


	qsort( group, group_count, sizeof(Group), compareGroup);


	samfile_t *temp_bam = samopen("temp.bam", "wb", header );
	
		if (temp_bam == NULL) 
			{
			fprintf(stderr,"Failed to open output file . %s.\n",strerror(errno));
			return EXIT_FAILURE;
			} 

	
	DEBUG;
	while(sam_read1(fp, header, b) >= 0)
	{ 
		nReads++;
		if ( (b->core.flag & BAM_FUNMAP ) )
			{
				
			search = bam_aux_get(b, "RG");
			Groupsearch = bam_aux2Z(search);
			
			struct Group key, *res;
			key.rgId = Groupsearch;
			res = bsearch(&Groupsearch, group, group_count, sizeof(Group), compareGroup);
			if (res == NULL)
				{
				fprintf(stderr, "%s unknown read group\n", Groupsearch);
				return EXIT_FAILURE;
				}
			else 
				{
				res->sample->unMap++;
				}
				
			samwrite(temp_bam, b); 	
			}
		samwrite(out_file, b); 
	}
	
	

 	file=fopen(output_report,"w"); 
 		if ( file != NULL) 
 		{
 		for(i=0;i< sample_count;++i)
 			{
      			fprintf(file,"%s\t%lu \n",samples[i].sample_name, samples[i].unMap);
      			}
      		
      		} else 
      			{
	      	 	fprintf(stderr,"Cannot read file. %s.\n",strerror(errno));
	      	 	return EXIT_FAILURE;
      	 		}     	 		
      	 		
      	 		
      	DEBUG;	 		
        fclose(file);
	bam_destroy1(b);
	sam_close(fp); 
	samclose(out_file);
	samclose(temp_bam);
	free(filename_out);


	
	for (i=0; i<group_count; i++)
		{
		free(group[i].rgId);
		}
	free(group);

	
	for (i=0; i<sample_count; i++)
		{
		free(samples[i].sample_name);
		}
	free(samples);
	

	
        return EXIT_SUCCESS;
        }
        
        
        
        // ensuite possible de transformer du BAM en FASTQ par cette commande : samtools bam2fq test.bam > test.fastq
        // le fichier .fastq va pouvoir servir pour faire les alignements avec les réferences des contaminants.
        
