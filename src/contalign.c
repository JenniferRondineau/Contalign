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


/*version of this code */ 
static void version () {fputs(" Version : " COUNTALIGN_VERSION "\n",stdout);}



#define PRINT_OPTION(SHORT,LONG,TYPE,DEF) \
	fputs("   -" SHORT ",--" LONG "   " TYPE "  " DEF ".\n",stderr)   


static void usage ()
	{ 
	fputs(
	"Version: " COUNTALIGN_VERSION "\n"
	"Usage : cat file.bam | ./contalign -r file.fa -s file.fastq -o file.txt > file.bam \n"
	"\n"
	"Description : \n"
	"\n"
	"Contalign is a software allowing to read BAM files and looks for unmapped reads. Contalign map unmapped reads against a large reference of contaminants. It releases a report containing the number of unmapped reads by sample and the potential contaminants.  \n\
	\n\
	Options :\n",stderr);
	PRINT_OPTION("o","report","FILE","Name of the output file containing the report of unmapped reads");
	PRINT_OPTION("O","output","FILE","Name of the output BAM file (*.bam)");
	PRINT_OPTION("s","save","FILE","Save FASTQ file (*.fastq)");
	PRINT_OPTION("r","databwa","FILE","Reference file (*.fa)");
	PRINT_OPTION("h","help","FILE","Output help and exit immediately.");
	PRINT_OPTION("v","version","FILE","Output version and exit immediately.");
	}



static ContalignPtr ContalignNew()
	{
	return malloc(sizeof(Contalign)) ;
	}
	
 

static void ContalignRelease(ContalignPtr app)  // Close files, free 
	{
	bwa_idx_destroy(app->idx);
	fclose(app->file);	
	if (app->file_fastq != NULL) fclose(app->file_fastq);
	if (app->full_report != NULL) fclose (app->full_report);
	free(app);
	}


int main(int argc, char** argv)
	{
	
	ContalignPtr app = ContalignNew();
	char *type;

	// get commande line options
	int c=0, option_index = 0;
        static struct option long_options[] = {
       	    {"help", no_argument, NULL, 'h'},
       	    {"version", no_argument, NULL, 'v' },
            {"save",  required_argument, 0, 's'},
            {"report",  required_argument, 0, 'o'},
            {"databwa",  required_argument, 0, 'r'},
            {"output",  optional_argument, 0, 'O'},
            {"input",  required_argument, 0, 'I'},
            {"full_report",  required_argument, 0, 'c'},
            {0,         0,                 0,  0 }
        };
        
	 while ((c = getopt_long(argc, argv,"r:o:O:hvs:I:c:", long_options, &option_index)) != -1) {

		 switch(c) {

			case 'h': usage(); return EXIT_SUCCESS;  break;
			
			case 'v': version (); return EXIT_SUCCESS; break;
			 
			case 'o': app->output_report = optarg; break;
			
			case 'I': app->filename_in = optarg; break;

			case 's': app->filename_fastq = optarg; break;
				
			case 'r': app->ref = optarg; break;

			case 'O': app->filename_out = optarg; break;
			
			case 'c': app->namefull_report = optarg; break;
			 
			case '?': fprintf(stderr, "ERROR: option -%c is undefined\n", optopt); return EXIT_FAILURE; break;
			
			case ':': fprintf(stderr, "ERROR: option -%c requires an argument\n",optopt); return EXIT_FAILURE; break;
			 	
			default : usage(); return EXIT_FAILURE; break;
	
		}
	}
	
	
	OpenFile(app) ; 

	if(optind==argc) // if option "--inputfile"
		{
		runAppl(app);	
		}
	else if(optind+1==argc) // only one BAM file or ".list" file 
		{
		app->filename_in=argv[optind]; 
		type = strpbrk(argv[optind], ".");
		if (strcmp(type,".list")==0) { // search if the file format is ".list"
			FILE* list;
			char path[FILENAME_MAX];
			char* fin="\0";
			list=fopen(app->filename_in,"r"); //open file for reading
			if ( list == NULL) {
		      		fprintf(stderr,"Cannot open file. %s.\n",strerror(errno));
		      		exit(EXIT_FAILURE);
	      	 	}
	      	 	while(fgets(path,FILENAME_MAX, list) != 0) // read line by line
            			{ 	 
            			if(path[strlen(path) - 1] == '\n')
            			path[strlen(path) - 1] = '\0';  
            			if(strcmp(path, fin)==0) break;
            			app->filename_in=path;
            			runAppl(app);    // probleme au niveau du nom des contaminants dès qu'il y a plus d'un bam a aligner
            			}	
            		fclose(list);
		} else runAppl(app);

		}
	else
		{
		// list of BAM files
		while(optind<argc)
			{
			app->filename_in=argv[optind];
			runAppl(app);
			++optind;
			}
		}
	
	
	ContalignRelease(app);  // Close files, free 
	
	
return EXIT_SUCCESS;
}
        
        
