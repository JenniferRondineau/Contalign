
.PHONY:all test clean

sample.bam= sample.bam
CC=gcc
samtools.version=1.2
samtools.dir=samtools-${samtools.version}
htslib.version=1.2.1
htslib.dir=htslib-${htslib.version}
bwa.version=0.7.12
bwa.dir=bwa-${bwa.version}
CFLAGS= -Wall -g -I ${samtools.dir} -I ${bwa.dir}
LDFLAGS=-L ${samtools.dir} -L ${htslib.dir} -L ${bwa.dir}
sam.libs= -lbam -lhts -lbwa -lm -lz -lpthread
PACKAGES= bwa-${bwa.version}/libbwa.a samtools-${samtools.version}/libam.a htslib-${htslib.version}/libhts.a

all: contalign test 

test: contalign.txt
	@echo -n "Nombre de Reads dans ${sample.bam} : "
	@${samtools.dir}/samtools view -c ${sample.bam}
	@echo -n "Nombre de Reads dans test.bam : "
	@${samtools.dir}/samtools view test.bam | wc -l
	@echo -n "Nombre de Reads NON MAPPES dans ${sample.bam} : "
	@${samtools.dir}/samtools view -f 4 -c ${sample.bam}
	@echo -n "Nombre de Reads NON MAPPES dans $< : "
	@cat $<

contalign.txt : contalign sample reference
	cat ${sample.bam} | LD_LIBRARY_PATH=${htslib.dir} ./contalign -r UniVec.fa -o $@ -s save.fastq > test.bam

contalign : contalign.o 
	${CC} -o $@  ${LDFLAGS} $< ${sam.libs} 

contalign.o : contalign.c debug.h $(PACKAGES)
	${CC} -c -o $@ ${CFLAGS} $< 

bwa-${bwa.version}/libbwa.a:
	rm -rf bwa-${bwa.version}
	wget -O ${bwa.version}.tar.gz "https://github.com/lh3/bwa/archive/${bwa.version}.tar.gz"
	tar xfvz ${bwa.version}.tar.gz
	rm ${bwa.version}.tar.gz
	(cd bwa-${bwa.version} ; make)
	
samtools-${samtools.version}/libam.a:
	rm -rf samtools-${samtools.version}
	wget -O ${samtools.version}.tar.bz2 "https://github.com/samtools/samtools/releases/download/${samtools.version}/samtools-${samtools.version}.tar.bz2"
	tar xvjf ${samtools.version}.tar.bz2
	rm ${samtools.version}.tar.bz2
	(cd samtools-${samtools.version} ; make)
	
htslib-${htslib.version}/libhts.a:
	rm -rf htslib-${htslib.version}
	wget -O ${htslib.version}.tar.gz "https://github.com/samtools/htslib/archive/${htslib.version}.tar.gz"
	tar xfvz ${htslib.version}.tar.gz
	rm ${htslib.version}.tar.gz
	(cd htslib-${htslib.version} ; make)

sample : 
	wget -O - "ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/data/HG00096/exome_alignment/HG00096.chrom20.ILLUMINA.bwa.GBR.exome.20120522.bam" | ../src/samtools-${samtools.version}/samtools view -h - | head -n 50000 > sample.bam

reference : 
	wget -O - "ftp://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec" > UniVec.fa
	bwa-${bwa.version}/bwa index UniVec.fa
	

clean:
	rm -f contalign
	rm -rf $(dir $(TARGETS))
	rm -f htslib-${htslib.version}
	rm -f samtools-${samtools.version}
	rm -f bwa-${bwa.version}


