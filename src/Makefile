
.PHONY:all test clean

sample.bam=/home/bird/masterm2/sample3.bam
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

all: contalign test20150317 test

test: contalign.txt
	@echo -n "Nombre de Reads dans ${sample.bam} : "
	@${samtools.dir}/samtools view -c ${sample.bam}
	@echo -n "Nombre de Reads dans test.bam : "
	@${samtools.dir}/samtools view test.bam | wc -l
	@echo -n "Nombre de Reads NON MAPPES dans ${sample.bam} : "
	@${samtools.dir}/samtools view -f 4 -c ${sample.bam}
	@echo -n "Nombre de Reads NON MAPPES dans $< : "
	@cat $<

contalign.txt : contalign
	cat  ${sample.bam} | LD_LIBRARY_PATH=${htslib.dir} ./contalign -r UniVec.fa -o $@ -s save.fastq > test.bam
	
contalign : contalign.o
	${CC} -o $@  ${LDFLAGS} $< ${sam.libs}

contalign.o : contalign.c debug.h 
	${CC} -c -o $@ ${CFLAGS} $<
	
test20150317: test20150317.c
	${CC} -o $@ $<


bwa-${bwa.version}/libbwa.a:
	wget -O ${bwa.version}.tar.gz "https://github.com/lh3/bwa/archive/${bwa.version}.tar.gz"
	tar xfvz ${bwa.version}.tar.gz
	rm ${bwa.version}.tar.gz
	(cd bwa-${bwa.version}.version} ; make)


samtools-${samtools.version}/libam.a:
	wget -O ${samtools.version}.tar.bz2 "https://github.com/samtools/samtools/releases/download/1.2/samtools-${samtools.version}.tar.bz2"
	tar xvjf ${samtools.version}.tar.bz2
	rm ${samtools.version}.tar.bz2
	(cd samtools-${samtools.version} ; make)


htslib-${htslib.version}/libhts.a:
	wget -O ${htslib.version}.tar.gz "https://github.com/samtools/htslib/archive/${htslib.version}.tar.gz"
	tar xfvz ${htslib.version}.tar.gz
	rm ${htslib.version}.tar.gz
	(cd htslib-${htslib.version} ; make)


clean:
	rm -f contalign
	rm -rf $(dir $(TARGETS))
	rm -f htslib-${htslib.version}
	rm -f samtools-${samtools.version}
	rm -f bwa-${bwa.version}


	
# LD_LIBRARY_PATH=../../packages/samtools/htslib ./contalign
