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
#include <unistd.h>
#include "sam.h"
#include "sam_header.h"
#include "debug.h"
#include <errno.h>
#include "kseq.h" 
#include "bwamem.h"
#include "contalign.h"
#include <getopt.h>


#define VERIFY_NOT_NULL(POINTER) do{if(POINTER==NULL) { fprintf(stderr,"Memory Alloc failed File %s Line%d.\n",__FILE__,__LINE__); exit(EXIT_FAILURE);}}while(0)



/* write fastq in file_fastq */ 						
#define DUMP_FASTQ if ( app->file_fastq != NULL) {\
		   for (i=0; i < fastq_count; i++)\
			{\
			fprintf(app->file_fastq, "@%s%s\n%s\n+\n%s\n", fastq[i].qname, fastq[i].read, fastq[i].seq, fastq[i].qual);\
			}\
		   } else \
			{\
			fprintf(stderr,"Unsaved fasta file.\n");\
		   }


void OpenFile(ContalignPtr app) 
	{
	if (app->filename_fastq != NULL) //open the fastq file, if option "save" is indicated 
		{
		app->file_fastq=fopen(app->filename_fastq,"w+");
		if ( app->file_fastq == NULL) 
 			{
		      	fprintf(stderr,"Cannot open file. %s.\n",strerror(errno));
		      	exit(EXIT_FAILURE);
	      	 	}     	 	
		}		
			
	if (app->namefull_report != NULL) //open the full report, if option "full_report" is indicated
			{
			app->full_report =fopen(app->namefull_report ,"w+");
			if (app->full_report  == NULL) 
 				{
		      	 	fprintf(stderr,"Cannot open file. %s.\n",strerror(errno));
		      	 	exit(EXIT_FAILURE);
	      	 		}     	 	
			}

	if(app->ref == NULL) {
		fprintf(stderr, "Index load failed : %s.\n",strerror(errno));
		exit(EXIT_FAILURE);
	} else { 
		app->idx = bwa_idx_load(app->ref, BWA_IDX_ALL); // load the BWA index
		if (NULL == app->idx) {
			fprintf(stderr, "Index load failed : %s.\n",strerror(errno));
			exit(EXIT_FAILURE);
		}
	}	

	 if (app->output_report != NULL) {
	 	app->file=fopen(app->output_report,"w"); // open the report of unmapped reads and contaminants
		if ( app->file == NULL)
			{
			fprintf(stderr,"Cannot write file. %s.\n",strerror(errno));
			exit(EXIT_FAILURE);
			} 
	}
	
}



int compareGroup(const void *g1, const void *g2)
{
	Group *group1 = (Group *) g1;
	Group *group2 = (Group *) g2;
	return strcmp (group1->rgId, group2->rgId);  //sort Read Group in alphabetical order
}


int compareContaminant(const void *c1, const void *c2)
{
	Contaminants *contaminant1 = (Contaminants *) c1;
	Contaminants *contaminant2 = (Contaminants *) c2;
	return contaminant2->contaminants_count - contaminant1->contaminants_count; //sort contaminants in descending order
}

int compareNameContaminant(const void *c1, const void *c2)
{
	Contaminants *contaminant1 = (Contaminants *) c1;
	Contaminants *contaminant2 = (Contaminants *) c2;
	return strcmp (contaminant1->c_name, contaminant2->c_name); //sort contaminants in alphabetical order
}

Contaminants *align( Fastq* fastq, Contaminants* contaminant, Sample* samples, int fastq_count, int *count_contaminants, int sample_count, ContalignPtr app)
{
	int i=0, n=0, len=0, sample=0;
	mem_opt_t *opt;
	bwaidx_t *idx = app->idx;
	opt = mem_opt_init(); // initialize the BWA-MEM parameters to the default values
	for(sample=0;sample< sample_count;++sample) { 
		for(n=0; n<fastq_count;n++) {
			fastq[n].score = 0 ;
			
			len= strlen(fastq[n].seq);
			mem_alnreg_v ar;
			ar = mem_align1(opt, idx->bwt, idx->bns, idx->pac, len, fastq[n].seq);  // get all the hits
			for (i = 0; i < ar.n; ++i) // traverse each hit
			{
				
				mem_aln_t a;
				if (ar.a[i].secondary >= 0) continue; // skip secondary alignments
				a = mem_reg2aln(opt, idx->bns, idx->pac, len, fastq[n].seq, &ar.a[i]); 
				
				/* Example of report can be obtained after alignment :
				
				int k=0;
				fprintf(stderr,"%s%s\t%c\t%s\t%ld\t%d\t%d\t", fastq[n].qname, fastq[n].read, "+-"[a.is_rev], idx->bns->anns[a.rid].name, (long)a.pos, a.mapq, a.score); // print alignment
				for (k = 0; k < a.n_cigar; ++k) // print CIGAR
					fprintf(stderr,"%d%c", a.cigar[k]>>4, "MIDSH"[a.cigar[k]&0xf]);
					fprintf(stderr,"\n");
				free(a.cigar); // deallocate CIGAR
				
				*/
				

				if (fastq[n].score < a.score ) //compare if this read has been previously mapped to another contaminant
					{	
						
					Contaminants key;
					key.c_name=idx->bns->anns[a.rid].name;

					struct Contaminants *res = NULL; 
					res = bsearch(    // search if the contaminant is already registered
						&key,
						contaminant,
						*count_contaminants, 
						sizeof(Contaminants),
						compareNameContaminant
						);

						if (res != NULL) { // if the contaminant is already registered 
							res->contaminants_count++; 
						}
						else {   // if the contaminant isn't registered 
							contaminant = (Contaminants*)realloc(contaminant,sizeof(Contaminants)*(*count_contaminants+1)); 
							VERIFY_NOT_NULL(contaminant);
							contaminant[*count_contaminants].c_name = idx->bns->anns[a.rid].name; //Save the new contaminant 
							VERIFY_NOT_NULL(contaminant[*count_contaminants].c_name);
							contaminant[*count_contaminants].contaminants_count=1;
							res=&contaminant[*count_contaminants];
							*count_contaminants=*count_contaminants+1;
							qsort( contaminant, *count_contaminants, sizeof(Contaminants), compareNameContaminant); // sort Contaminant in alphabetical order			
						}
			
					fastq[n].contaminant=res;
					fastq[n].score = a.score;
					
				}
				
				//write the full report, if option "full_report" is indicated
				if (app->full_report != NULL && fastq[n].score != 0 ) { 
							fprintf(app->full_report,"%s\t%d\t%s\t%s\t%s\t%s\n", samples[sample].sample_name,samples[sample].fastq[n].original_BAMFLAG,samples[sample].fastq[n].contaminant->c_name,samples[sample].fastq[n].qname, samples[sample].fastq[n].seq, samples[sample].fastq[n].qual);	

				}
		
			}
		// free memory allocated
		free(fastq[n].qname);
		free(fastq[n].seq);
		free(fastq[n].qual);
		free(ar.a);
		}
	}
	free(opt);

return contaminant;
}



void runAppl(ContalignPtr app)
	{
	int i=0, sample_count=0, group_count=0, fastq_count=0, count_contaminants=0 ;
	bam_hdr_t *header=NULL;
	bam1_t *b = NULL;
	b = bam_init1();
	float nReads=0,  pourcent =0;
	const char *rgId=NULL, *sampleName=NULL; 
	Sample* samples=NULL;
	Contaminants* contaminant=NULL;
	Group* group=NULL;
	Fastq* fastq=NULL;
	uint8_t *search = NULL, *seq= NULL, *qual = NULL;
	char *qname = NULL, *read = NULL;
	int8_t *buf = NULL;
	int32_t qlen =NULL;
	
	fprintf(stderr,"Reading from %s...\n",(app->filename_in!=NULL? app->filename_in : "<<stdin>>")); // message indicating the name of the file read
	app->fp = sam_open( app->filename_in!=NULL? app->filename_in : "-" );  // opening bam file
	if (app->fp == NULL) 
		{
		fprintf(stderr,"Cannot read file \"%s\". %s.\n",
			(app->filename_in!=NULL? app->filename_in : "<<stdin>>"),
			strerror(errno)
			);
		exit(EXIT_FAILURE);
		}

	header = sam_hdr_read(app->fp); //read and copy header of initial BAM file

		if( header == NULL)
			{
			fprintf(stderr, "Cannot read header \n");
			exit(EXIT_FAILURE);
			}

	app->out_file = samopen((app->filename_out!=NULL? app->filename_out : "-"), "wb", header );  //create the output BAM file with the initial header
	
		if (app->out_file == NULL) 
			{
			fprintf(stderr,"Failed to open output file . %s.\n",strerror(errno));
			exit(EXIT_FAILURE);
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
					exit(EXIT_FAILURE);
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
		exit(EXIT_FAILURE);
		}

	qsort( group, group_count, sizeof(Group), compareGroup);  //sort Read Group in alphabetical order 

	while(sam_read1(app->fp, header, b) >= 0) // read each reads of BAM file
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
				exit(EXIT_FAILURE);
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
				exit(EXIT_FAILURE);
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
				exit(EXIT_FAILURE);
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
				exit(EXIT_FAILURE);
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
				exit(EXIT_FAILURE);
			}
           		buf=calloc((qlen+1),sizeof(int8_t));       			
       			buf[qlen] = '\0';
       			for (i = 0; i < qlen; ++i) buf[i] = 33 + qual[i];


           		fastq[fastq_count].qual = strdup((char*)buf);
           		VERIFY_NOT_NULL(fastq[fastq_count].qual);
           		
           		fastq[fastq_count].original_BAMFLAG = b->core.flag;
           		
			res->sample->fastq = fastq;
			

			free(buf);
			fastq_count++;	

			if(fastq_count == 1000)  
				{
				DUMP_FASTQ; // Write fastq file
				contaminant = align(fastq, contaminant, samples, fastq_count, &count_contaminants, sample_count, app); //mapping unmapped reads against reference of contaminants
				fastq_count=0;
				}
			

				
			}
		samwrite(app->out_file, b); // Write read in the output file
	}
	
        //write and mapping the latest reads
	DUMP_FASTQ; 
	contaminant = align(fastq, contaminant, samples, fastq_count, &count_contaminants, sample_count, app); 

	
	qsort( contaminant, count_contaminants, sizeof(Contaminants), compareContaminant); //sort the contaminants in descending order



	// Write the report containing number of unmapped reads by sample and potential contaminants
	for(i=0;i< sample_count;++i)
		fprintf(app->file,"%s\t%lu \n",samples[i].sample_name, samples[i].unMap);
	for(i=0;i< count_contaminants;++i)
		{
		pourcent = (contaminant[i].contaminants_count / nReads)*100;
		fprintf(app->file,"%f%% \t (%.0f/%.0f) \t%s\n", pourcent,contaminant[i].contaminants_count,nReads, contaminant[i].c_name);
		}


	// Close files, free and return
	sam_close(app->fp);
	samclose(app->out_file);
	for (i=0; i<count_contaminants; i++)
		{
		free(contaminant[i].c_name);
		}
	free(contaminant);
	
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
	bam_destroy1(b);
	
}











