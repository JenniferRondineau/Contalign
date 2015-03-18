/**

cat file.sam | countalign > count.txt

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include "sam.h"

#include "debug.h"



int main(int argc, char** argv)
	{

	FILE * fichier;
	FILE * rapport; 

	
	rapport=fopen("rapport.bam","w+");
	fclose(rapport);
	DEBUG;
	samfile_t *fp = samopen("/home/bird/packages/samtools/examples/toy.sam", "rb", 0); 
	//samfile_t *out = samopen("rapport.bam", "w", 0); 
	DEBUG;
	bam1_t *b = NULL;

	long nReads=0;
	long nbNoMap=0;
	b = bam_init1();
	
	DEBUG;

	if(NULL == fp)
	{
		fprintf(stderr,"fichier non ouvert\n");
	} 
		 

	

		DEBUG;

		while(samread(fp, b) >= 0)
		{ 
			nReads++;
			if ( b->core.flag == 163 )
				{
				nbNoMap++;
				rapport=fopen("rapport.bam","w+");
				fprintf(rapport, "%d %d %d %d %d %d %d %d ", b->core.flag, b->core.pos, b->core.n_cigar, b->core.tid, b->core.bin, b->core.qual, b->core.l_qname, b->core.l_qseq);
				fclose(rapport);
				}
  
		}
		DEBUG;

 
 		fichier=fopen("count.txt","w+");
      		fprintf(fichier,"Nombre de reads : %lu et le nombre de reads non mappés : %lu \n", nReads, nbNoMap);
        	fclose(fichier);


		bam_destroy1(b);
		samclose(fp); 
		
	


        return 0;
        }
