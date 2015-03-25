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




static void usage ()
{ 
/* cf isatty */
fprintf(stderr,
"\nUsage : cat file.bam | ./countalign -o file.txt > file.bam \n\
\n\
Description : \n\
\n\
This program reads BAM file and looks for unmapped reads. It released a report containing the total number of reads and the number of unmapped reads. \n\
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

static int compareSample(const void *s1, const void *s2)
{
	Sample *sample1 = (Sample *) s1;
	Sample *sample2 = (Sample *) s2;
	return strcmp (sample1->sample_name, sample2->sample_name);
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
	long UnMap=0;
	int n=0, i=0, j=0, len=0;
	char *nameResearch=NULL; 
	char *rgId=NULL, *sampleName=NULL;
	
	Sample* samples=NULL;
	int sample_count=0;
	Group* group=NULL;
	int group_count=0;

	char *verifuniq_sample=NULL; 
	int uniq_sample = 0;

	char *verifuniq_group=NULL; 
	int uniq_group = 0;

	
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
			return EXIT_FAILURE;
		
		 }
	
		if (argv[1] == NULL)
		{
			return EXIT_FAILURE;
			}

	 } // si ./ countalign aljsd 

	
	
	samFile *fp = sam_open( "-" ); 
	
		if (fp == NULL) 
		{
        		fprintf(stderr,"Cannot read file. %s.\n",strerror(errno));
        		return EXIT_FAILURE;
		}
	DEBUG;
	header = sam_hdr_read(fp);

	
		if( header == NULL)
		{
			fprintf(stderr, "Cannot read header \n");
			return EXIT_FAILURE;
		}


	
	DEBUG;

	samfile_t *out_file = samopen((filename_out!=NULL? filename_out : "-"), "wb", header ); 
	
		if (out_file == NULL) 
		{
			fprintf(stderr,"Failed to open output file . %s.\n",strerror(errno));
			return EXIT_FAILURE;
		} 

DEBUG;
	


		void *iter = sam_header_parse2(header->text);
		while ( (iter = sam_header2key_val(iter, "RG","ID","SM",&rgId,&sampleName) )!=NULL)
			{
			
				 /* sample must be unique */
				verifuniq_sample = strdup(sampleName);
				
				for (i=0;i<sample_count;++i)
					{ 
					if (strcmp (samples[i].sample_name, verifuniq_sample) == 0)
						{
						uniq_sample=1;
						}
					}
				
				if (uniq_sample == 0 ) 
				{ 					
					samples = realloc(samples,sizeof(Sample)*(n+1));
					VERIFY_NOT_NULL(samples);
					samples[sample_count].sample_name = strdup(sampleName);
					VERIFY_NOT_NULL(samples[sample_count].sample_name);
					sample_count++;
				}


				 /* group must be unique */			
				verifuniq_group = strdup(rgId);
				
				for (i=0;i<group_count;++i)
					{ 
					if (strcmp (group[i].rgId, verifuniq_group) == 0)
						{
						uniq_group=1;
						}
					}
					
				if (uniq_group == 0 ) 
				{ 					
					group = realloc(group,sizeof(Group)*(n+1));	
					VERIFY_NOT_NULL(group);	
					group[group_count].rgId = strdup(rgId);
					VERIFY_NOT_NULL(group[group_count].rgId);
					group_count++;
				}


				DEBUG;				
				
				
			}
		

//qsort(&firstSample, count_samples, sizeof(Sample), compareSample);




DEBUG;
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
			nameResearch = strdup(bam_get_qname(b));
			i=0;
			for (i=0; i < group_count ; i++) 
			{
				len=strlen(group[i].rgId);
									
				if(strncmp ( group[i].rgId, nameResearch, len) == 0)
				{
					fprintf(stderr,"yehhhhh %s \n ", group[i].rgId);				
				}
			}
			//comparer avec ceux dans la structure
			UnMap++;
			//samwrite(temp_bam, b); 	
			}
		samwrite(out_file, b); 
  
	}
	
fprintf(stderr,"le nom recherché %s\n",nameResearch);


	 	/* a la fin, faire un rapport PAR SAMPLE 
 		
 		faudra utiliser 
 		* malloc (nbre echantillons)
 		* bsearch recherche par dichotnomy)q
 		* une structure de données
 	
 		char *sample_name = (char*) malloc(sizeof(char)*(taille+1));
 		char *sample_name = (char*) calloc(truc+1, sizeof(char));
 		free(sample_name);

 		*/
	
 	file=fopen(output_report,"w"); 
 		if ( file != NULL) 
 		{
 		/*for(i=0;i< sample_count;++i)
 			{
      			fprintf(file,"%s\t%lu \n",samples[i].name,samples[i].UnMap);
      			}*/
      			fprintf(file,"%d\n",UnMap);
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

	//cleanup
	//pour tous samples: free(name)
	free(samples);
	//pour tous gpourps: free(name)
	free(group);

	
        return EXIT_SUCCESS;
        }
        
        
        
        // ensuite possible de transformer du BAM en FASTQ par cette commande : samtools bam2fq test.bam > test.fastq
        // le fichier .fastq va pouvoir servir pour faire les alignements avec les réferences des contaminants.
        
