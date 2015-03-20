#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sam.h"

#include "debug.h"



int main(int argc, char** argv)
	{

	FILE* input_file=NULL;
	FILE* output_file=NULL;
	FILE * fichier;	

	DEBUG;
	int c;
	char* f_in; 



	
	// get commande line options
	 while ((c = getopt(argc, argv,"o:O:")) != -1) {
	 switch(c)
	 {
	 case 'o': f_in = strdup(optarg); break;
	// case 'O' : f_out = strdup(optarg);   break;
	 }
	 }

	
	
 	bam_hdr_t *header=NULL;


	DEBUG;
	
	samFile *fp = sam_open((input_file? input_file : "-")); 
	
	if (fp == NULL) 
	{
        fprintf(stderr,"Cannot read file");
        return 1;
        }
        
	DEBUG;
	header = sam_hdr_read(fp);
	bam1_t *b = NULL;

	
	if( header == NULL)
	{
	fprintf(stderr, "Cannot read header \n");
	}


	DEBUG;
	samfile_t *out_file = samopen((output_file? output_file : "test.bam"), "wb", header ); // reussir à prendre le "> out.bam" au lieu de "sample.bam"

	
	DEBUG;
	
	
	if (out_file == NULL) 
	{
		fprintf(stderr,"Failed to open output file.\n");
	} 

	long nReads=0;
	long nbUnMap=0;
	b = bam_init1();
	
	DEBUG;

	if(NULL == fp)
	{
		fprintf(stderr,"Failed to open input file.\n");
	} 
		      

	

		DEBUG;

		while(sam_read1(fp, header, b) >= 0)
		{ 
			nReads++;
			if ( (b->core.flag & BAM_FUNMAP ) || (b->core.flag & BAM_FMUNMAP ))
				{
				nbUnMap++;
				samwrite(out_file, b); 
				}
  
		}
		DEBUG;
		
	
 		fichier=fopen(f_in,"w+");
      		fprintf(fichier,"Nombre de reads : %lu et le nombre de reads non mappés : %lu \n", nReads, nbUnMap);
        	fclose(fichier);
		DEBUG;

		bam_destroy1(b);
		
		sam_close(fp); 
		samclose(out_file);
		DEBUG;



        return 0;
        }
        
        
        
        // ensuite possible de transformer du BAM en FASTQ par cette commande : samtools bam2fq test.bam > test.fastq
        // le fichier .fastq va pouvoir servir pour faire les alignements avec les réferences des contaminants.
        
