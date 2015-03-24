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




int main(int argc, char** argv)
	{

	FILE* input_file=NULL;
	FILE* file=NULL;	
	char* filename_out =NULL; 
	int c;
	char* output_report; 
	bam_hdr_t *header=NULL;
	bam1_t *b = NULL;
	b = bam_init1();
	long nReads=0;
	long UnMap=0;




	
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


	 } // si ./ countalign aljsd 

	
	
	samFile *fp = sam_open((input_file? input_file : "-")); 
	
		if (fp == NULL) 
		{
        		fprintf(stderr,"Cannot read file. %s.\n",strerror(errno));
        		return EXIT_FAILURE;
		}
	DEBUG;
	header = sam_hdr_read(fp);
	/* ici aller voir les sources de header, recuperer les noms des échantillons, si pas d'échantillons, mettre echantillon par default "SAMPLE" */
	//chercher sam_header.c ; sam_header_line_parse ; 

	
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



	
	sam_header_parse2(header->text);




	samfile_t *temp_bam = samopen("temp.bam", "wb", header );
	
		if (temp_bam == NULL) 
		{
			fprintf(stderr,"Failed to open output file . %s.\n",strerror(errno));
			return EXIT_FAILURE;
		} 



	while(sam_read1(fp, header, b) >= 0)
	{ 
		nReads++;
		if ( (b->core.flag & BAM_FUNMAP ) )
			{
			UnMap++;
			samwrite(temp_bam, b); 	
			}
		samwrite(out_file, b); 
  
	}
	
	
	
 	file=fopen(output_report,"w"); 
 		if ( file != NULL) 
 		{
 		/* a la fin, faire un rapport PAR SAMPLE 
 		
 		faudra utiliser 
 		* malloc (nbre echantillons)
 		* bsearch recherche par dichotnomy)q
 		* une structure de données
 	
 		char *sample_name = (char*) malloc(sizeof(char)*(taille+1));
 		char *sample_name = (char*) calloc(truc+1, sizeof(char));
 		free(sample_name);
 		struct sample_count
 			{
 			char* sample_name;
 			long UnMap;
 			};
 		
 		*/
      			fprintf(file,"%lu \n",UnMap);
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
	

	
        return EXIT_SUCCESS;
        }
        
        
        
        // ensuite possible de transformer du BAM en FASTQ par cette commande : samtools bam2fq test.bam > test.fastq
        // le fichier .fastq va pouvoir servir pour faire les alignements avec les réferences des contaminants.
        
