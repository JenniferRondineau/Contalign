#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sam.h"

#include "debug.h"

/* sam_hdr_write(out, header) */

int main(int argc, char** argv)
	{

	FILE* input_file=NULL;
	FILE* output_file=NULL;
	FILE * fichier;	

	DEBUG;
	int c;
	char* f_in; 
	char* f_out; 

	
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
	header = sam_header_read(fp);
	bam1_t *b = NULL;

	if( header == NULL)
	{
	fprintf(stderr, "Try again \n");
	}


	DEBUG;
	samfile_t *out = samopen((output_file? output_file : "sample.bam"), "w", 0 ); // reussir à prendre le "> out.bam" au lieu de "sample.bam"
    	//sam_hdr_write(out, header); 
	DEBUG;

	long nReads=0;
	long nbNoMap=0;
	b = bam_init1();
	
	DEBUG;

	if(NULL == fp)
	{
		fprintf(stderr,"fichier non ouvert\n");
	} 
		      

	

		DEBUG;

		while(sam_read1(fp, header, b) >= 0)
		{ 
			nReads++;
			if ( (b->core.flag & BAM_FUNMAP ) || (b->core.flag & BAM_FMUNMAP ))
				{
				nbNoMap++;

				//samwrite(out, b); Une subtilité pour écrire dans le fichier sam ? allocation dynamique
				}
  
		}
		DEBUG;
		
	

 
 		fichier=fopen(f_in,"w+");
      		fprintf(fichier,"Nombre de reads : %lu et le nombre de reads non mappés : %lu \n", nReads, nbNoMap);
        	fclose(fichier);


		bam_destroy1(b);
		sam_close(fp); 
		samclose(out);
		
	


        return 0;
        }
