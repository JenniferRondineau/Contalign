# cible : dependances
#(tab)comment faire cible a partir de dependance
# $@ = nom de la cible
# $< est le nom de la PREMIERE dependance
.PHONY:all test clean
sample.bam=/home/bird/masterm2/sample3.bam
CC=gcc
PACKAGES=../../packages
samtools.dir=${PACKAGES}/samtools
htslib.dir=${samtools.dir}/htslib
bwa.dir=${PACKAGES}/bwa
CFLAGS= -Wall -g -I ${samtools.dir} -I ${bwa.dir}
LDFLAGS=-L ${samtools.dir} -L ${htslib.dir} -L ${bwa.dir}
sam.libs= -lbam -lhts -lbwa -lm -lz -lpthread

all: contalign test20150317 test

test: contalign.txt
	@echo -n "Nombre de Reads dans ${sample.bam} : "
	@${PACKAGES}/samtools/samtools view -c ${sample.bam}
	@echo -n "Nombre de Reads dans test.bam : "
	@${PACKAGES}/samtools/samtools view test.bam | wc -l
	@echo -n "Nombre de Reads NON MAPPES dans ${sample.bam} : "
	@${PACKAGES}/samtools/samtools view -f 4 -c ${sample.bam}
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

clean:
	rm -f contalign
	rm -rf $(dir $(TARGETS))


	
# LD_LIBRARY_PATH=../../packages/samtools/htslib ./contalign
