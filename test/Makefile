hts.version=1.2.1
samtools.version=1.2

all: test1 test2

test1: 	../src/contalign
	wget -O - "ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/data/HG00096/exome_alignment/HG00096.chrom20.ILLUMINA.bwa.GBR.exome.20120522.bam" | ../src/samtools-${samtools.version}/samtools view -h - | head -n 1000 | LD_LIBRARY_PATH=../src/htslib-${hts.version} $< -r UniVec.fa -o test.txt > /dev/null
	stat test.txt
	echo -n "sample name + number of unmapped read :" 
	tail test.txt | grep "HG" 

test2: 	../src/contalign
	wget -O - "ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/data/HG00096/exome_alignment/HG00096.chrom20.ILLUMINA.bwa.GBR.exome.20120522.bam" | ../src/samtools-${samtools.version}/samtools view -h - | head -n 10000 | LD_LIBRARY_PATH=../src/htslib-${hts.version} $< -r UniVec.fa -o test.txt > /dev/null
	stat test.txt
	echo -n "sample name + number of unmapped read :" 
	tail test.txt | grep "HG" 

../src/contalign: 
	(cd ../src ; make)

clean:
	(cd ../src ; make)
