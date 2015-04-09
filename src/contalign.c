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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sam.h"
#include "sam_header.h"
#include "debug.h"
#include <errno.h>
#include "bwa.c"
#include "kseq.h" 
#include "bwamem.h"
#include "contalign.h"
#include <getopt.h>

 	
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
static void version () {fputs(" Version : " COUNTALIGN_VERSION "\n",stdout);}

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

#define PRINT_OPTION(SHORT,LONG,TYPE,DEF) \
	fputs("   -" SHORT ",--" LONG "   " TYPE "  " DEF ".\n",stderr)   

static void usage ()
	{ 
	fputs(
	"\nUsage : cat file.bam | ./contalign -r file.fa -s file.fastq -o file.txt > file.bam \n\
	\n\
	Description : \n\
	\n\
	Contalign is a software allowing to read BAM files and looks for unmapped reads. Contalign map unmapped reads against a large reference of contaminants. It releases a report containing the number of unmapped reads by sample and the potential contaminants.  \n\
	\n\
	Options :\n",stderr);
	PRINT_OPTION("o","report","FILE","Name of the output file containing the report of unmapped reads");
	PRINT_OPTION("O","output","FILE","Name of the output BAM file (*.bam)");
	PRINT_OPTION("s","save","FILE","Save FASTQ file (*.fastq)");
	PRINT_OPTION("r","databwa","FILE","Reference file (*.fa)");
	PRINT_OPTION("h","help","FILE","Output help and exit immediately.");
	PRINT_OPTION("v","version","FILE","Output version and exit immediately.");
	}



Contaminants *align( Fastq* fastq, Contaminants* contaminant, int fastq_count, int *count_contaminants, bwaidx_t *idx)
{
	int i=0, n=0, len=0, j=0;
	mem_opt_t *opt;
	opt = mem_opt_init(); // initialize the BWA-MEM parameters to the default values
	for(n=0; n<fastq_count;n++) 
	{
		len= strlen(fastq[n].seq);
		mem_alnreg_v ar;
		ar = mem_align1(opt, idx->bwt, idx->bns, idx->pac, len, fastq[n].seq);  // get all the hits
		for (i = 0; i < ar.n; ++i) // traverse each hit
		{
			mem_aln_t a;
			if (ar.a[i].secondary >= 0) continue; // skip secondary alignments
			a = mem_reg2aln(opt, idx->bns, idx->pac, len, fastq[n].seq, &ar.a[i]); 
			printf("%s%s\t%c\t%s\t%ld\t%d\t\n", fastq[n].qname, fastq[n].read, "+-"[a.is_rev], idx->bns->anns[a.rid].name, (long)a.pos, a.mapq); // print alignment
			Contaminants* searchContaminant = NULL; 
			for (j=0; j< *count_contaminants; j++) 
			{      // search if the contaminant is already registered
				if (strcmp (contaminant[j].c_name, idx->bns->anns[a.rid].name) == 0)
					{
					searchContaminant = &contaminant[j];
					contaminant[j].contaminants_count++;
					break;
					}
			}
			if (searchContaminant == NULL ) 
			{
				contaminant = (Contaminants*)realloc(contaminant,sizeof(Contaminants)*(*count_contaminants+1)); 
				VERIFY_NOT_NULL(contaminant);
				contaminant[*count_contaminants].c_name = strdup(idx->bns->anns[a.rid].name); //Save the new contaminant 
				VERIFY_NOT_NULL(contaminant[*count_contaminants].c_name);
				contaminant[*count_contaminants].contaminants_count=1;
				searchContaminant=&contaminant[*count_contaminants];
				*count_contaminants=*count_contaminants+1;
			}
		}
	// free memory allocated
	free(fastq[n].qname);
	free(fastq[n].seq);
	free(fastq[n].qual);
	free(ar.a);
	}
	free(opt);
return contaminant;
}




int main(int argc, char** argv)
	{
	FILE* file=NULL;
	FILE* file_fastq=NULL;
	int c, i=0, sample_count=0, group_count=0, fastq_count=0;
	int count_contaminants=0, report_contaminant =0 ;
	bam_hdr_t *header=NULL;
	bam1_t *b = NULL;
	b = bam_init1();
	float nReads=0;
	const char *rgId=NULL, *sampleName=NULL;	
	Sample* samples=NULL;
	Contaminants* contaminant=NULL;
	Group* group=NULL;
	Fastq* fastq=NULL;
	uint8_t *search = NULL, *seq= NULL, *qual = NULL;
	char *qname = NULL, *read = NULL, *output_report=NULL, *filename_in=NULL, *filename_out =NULL, *ref=NULL;
	const char *filename_fastq = NULL; 
	int8_t *buf = NULL;
	int32_t qlen =NULL;
	bwaidx_t *idx;
	float pourcent =0;

	
	
	// get commande line options
        int option_index = 0;
        static struct option long_options[] = {
       	    {"help", no_argument, NULL, 'h'},
       	    {"version", no_argument, NULL, 'v' },
            {"save",  required_argument, 0, 's'},
            {"report",  required_argument, 0, 'o'},
            {"databwa",  required_argument, 0, 'r'},
            {"output",  optional_argument, 0, 'O'},
            {"input",  required_argument, 0, 'i'},
            {"report_contaminant", no_argument, 0, 'c'},
            {0,         0,                 0,  0 }
        };
        
	 while ((c = getopt_long(argc, argv,"r:o:O:hvs:i:c", long_options, &option_index)) != -1) {

		 switch(c)
			 {

		 switch(c) {

			case 'h': usage(); return EXIT_SUCCESS; break;
			
			case 'v': version (); return EXIT_SUCCESS; break;
			 
			case 'o': output_report = optarg; break;
			
			case 'i': filename_in = optarg; break;

			case 's': filename_fastq = optarg; break;
				
			case 'r': ref = optarg; break;

			case 'O': filename_out = optarg; break;
			 
			case 'c': report_contaminant = 1; break;

			case '?': fprintf(stderr, "ERROR: option -%c is undefined\n", optopt); return EXIT_FAILURE; break;
			
			case ':': fprintf(stderr, "ERROR: option -%c requires an argument\n",optopt); return EXIT_FAILURE; break;
			 	
			default : { usage(); return EXIT_FAILURE; break; }
			
		}
			 
	 }

	if (filename_fastq != NULL) //open the fastq file, if option "save" is indicated 
			{
			file_fastq=fopen(filename_fastq,"w+");
			if ( file_fastq == NULL) 
 				{
		      	 	fprintf(stderr,"Cannot read file. %s.\n",strerror(errno));
		      	 	return EXIT_FAILURE;
	      	 		}     	 	
			}		
			
	
	samFile *fp = sam_open( filename_in!=NULL? filename_in : "-" ); 
	
		if (fp == NULL) 
			{
        		fprintf(stderr,"Cannot read file. %s.\n",strerror(errno));
        		return EXIT_FAILURE;
			}

	header = sam_hdr_read(fp); //read and copy header of initial BAM file

		if( header == NULL)
			{
			fprintf(stderr, "Cannot read header \n");
			return EXIT_FAILURE;
			}


	samfile_t *out_file = samopen((filename_out!=NULL? filename_out : "-"), "wb", header );  //create the output BAM file with the initial header
	
		if (out_file == NULL) 
			{
			fprintf(stderr,"Failed to open output file . %s.\n",strerror(errno));
			return EXIT_FAILURE;
			} 
		
			

	void *iter = sam_header_parse2(header->text);
	while ( (iter = sam_header2key_val(iter, "RG","ID","SM",&rgId,&sampleName) )!=NULL) //sam_header2key_val() allows to find 'RG' tag in the header of BAM file. 
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
				samples[sample_count].sample_name = strdup(sampleName); //Save the new sample
				VERIFY_NOT_NULL(samples[sample_count].sample_name);
				samples[sample_count].unMap=0;
				
				searchSample=&samples[sample_count]; 
				
				sample_count++;
				}

			 /* group must be unique */							
			
			for (i=0;i<group_count;++i) { 
				if (strcmp (group[i].rgId,rgId) == 0) {
					fprintf(stderr,"DUP GROUP");
					return EXIT_FAILURE;
				}
			}
										
			group = (Group*)realloc(group,sizeof(Group)*(group_count+1));	
			VERIFY_NOT_NULL(group);	
			group[group_count].rgId = strdup(rgId); //Save the new read group
			VERIFY_NOT_NULL(group[group_count].rgId);
			
			group[group_count].sample = searchSample; 
			group_count++;							
		}
		

	if ((group_count == 0)  && (sample_count ==0))
		{
		fprintf(stderr,"Error, missing SAM header or @RG tag\n");
		return EXIT_FAILURE;
		}


	qsort( group, group_count, sizeof(Group), compareGroup);  //sort Read Group in alphabetical order 


	idx = bwa_idx_load(ref, BWA_IDX_ALL); // load the BWA index
		if (NULL == idx) {
			fprintf(stderr, "Index load failed %s.\n",strerror(errno));
			exit(EXIT_FAILURE);
		}
		



	while(sam_read1(fp, header, b) >= 0) // read each reads of BAM file
	{ 
		nReads++;
		if ( (b->core.flag & BAM_FUNMAP ) ) //detect unmapped reads
			{
			Group key;
			if ((b->core.flag & BAM_FMUNMAP)) 
				{
				fprintf(stderr, "Not orphelin \n");
			}
			search = bam_aux_get(b, "RG"); 
			if (search == NULL) {
				fprintf(stderr, "Error, missing @RG tag\n");
				return EXIT_FAILURE;
			}else key.rgId = bam_aux2Z(search); 
			
			struct Group *res;
			res = bsearch(
				&key,
				group,
				group_count,
				sizeof(Group),
				compareGroup
				); 
			if (res == NULL) {
				return EXIT_FAILURE;
			} else 
				{
				res->sample->unMap++; // Count unmapped reads by sample
				}
						
			fastq = (Fastq*)realloc(fastq,sizeof(Fastq)*(fastq_count+1));	
			VERIFY_NOT_NULL(fastq);	
			

			// Write sequence
			qlen = (b->core.l_qseq);
			if (qlen == 0) {
				fprintf(stderr, "Error, missing sequence\n");
				return EXIT_FAILURE;
			}
       			buf=calloc((qlen+1),sizeof(int8_t));       			
       			buf[qlen] = '\0';
       			seq = bam_get_seq(b);
        		for (i = 0; i < qlen; ++i){
           			buf[i] = bam_seqi(seq, i);
           			buf[i] = seq_nt16_str[buf[i]];
           			}

			fastq[fastq_count].seq = strdup((char*)buf);
			VERIFY_NOT_NULL(fastq[fastq_count].seq);
			

			// Write read name
			qname = strdup(bam_get_qname(b));
			if (qname == NULL) {
				fprintf(stderr, "Error, missing read name\n");
				return EXIT_FAILURE;
			}
			if ((b->core.flag & BAM_FREAD1) && !(b->core.flag & BAM_FREAD2)) read="/1";
         		else if ((b->core.flag & BAM_FREAD2) && !(b->core.flag & BAM_FREAD1)) read="/2";
           		
           		fastq[fastq_count].qname = qname;
           		VERIFY_NOT_NULL(fastq[fastq_count].qname);
           		           		
           		fastq[fastq_count].read = read;
           		VERIFY_NOT_NULL(fastq[fastq_count].read);
           		           		
           		
           		// Write quality	
           		qual = bam_get_qual(b);	
           		if (qual == NULL) {
				fprintf(stderr, "Error, missing quality\n");
				return EXIT_FAILURE;
			}
           		buf=calloc((qlen+1),sizeof(int8_t));       			
       			buf[qlen] = '\0';
       			for (i = 0; i < qlen; ++i) buf[i] = 33 + qual[i];


           		fastq[fastq_count].qual = strdup((char*)buf);
           		VERIFY_NOT_NULL(fastq[fastq_count].qual);
           		
			res->sample->fastq = fastq;
			
			free(buf);
			fastq_count++;	

			if(fastq_count == 1000)  
				{
				DUMP_FASTQ; // Write fastq file
				contaminant = align(fastq, contaminant,fastq_count, &count_contaminants, idx); //mapping unmapped reads against reference of contaminants
				fastq_count=0;
				}		
			}
		samwrite(out_file, b); // Write read in the output file
	}
	
        //write and mapping the latest reads
	DUMP_FASTQ; 
	contaminant = align(fastq, contaminant,fastq_count, &count_contaminants, idx); 

	
	qsort( contaminant, count_contaminants, sizeof(Contaminants), compareContaminant); //sort the contaminants in descending order
	
	
	ET SI output_report est NULL ?
	
	// Write the report containing number of unmapped reads by sample and potential contaminants
 	file=fopen(output_report,"w"); 
 		if ( file != NULL) 
 		{
 		for(i=0;i< sample_count;++i)
 			{
      			fprintf(file,"%s\t%lu \n",samples[i].sample_name, samples[i].unMap);
      			}
      		for(i=0;i< count_contaminants;++i)
 			{
 			pourcent = (contaminant[i].contaminants_count / nReads)*100;
      			fprintf(file,"%f%% \t (%.0f/%.0f) \t%s\n", pourcent,contaminant[i].contaminants_count,nReads, contaminant[i].c_name);
      			}
      		} else 
      			{
	      	 	fprintf(stderr,"Cannot read file. %s.\n",strerror(errno));
	      	 	return EXIT_FAILURE;
      	 		}     	 		
 	

        
        // Close files, free and return
       	for (i=0; i<count_contaminants; i++)
		{
		free(contaminant[i].c_name);
		}
	free(contaminant);

	if (file_fastq != NULL) fclose(file_fastq);

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
	
	bwa_idx_destroy(idx);
        fclose(file);
	bam_destroy1(b);
	sam_close(fp);
	samclose(out_file);


	
return EXIT_SUCCESS;
}
        
        
