# cible : dependances
#(tab)comment faire cible a partir de dependance
# $@ = nom de la cible
# $< est le nom de la PREMIERE dependance
.PHONY:all test clean
sample.bam=/home/bird/masterm2/sample.bam
CC=gcc
PACKAGES=../../packages
samtools.dir=${PACKAGES}/samtools
htslib.dir=${samtools.dir}/htslib
CFLAGS= -Wall -g -I ${samtools.dir}
LDFLAGS=-L ${samtools.dir} -L ${htslib.dir} 
sam.libs= -lbam -lhts

all: countalign test20150317 test

test: countalign.txt
	@echo -n "Nombre de Read dans ${sample.bam} : "
	@${PACKAGES}/samtools/samtools view -c ${sample.bam}
	@echo -n "Nombre de Reads dans test.bam : "
	@${PACKAGES}/samtools/samtools view -c test.bam
	@echo -n "Nombre de Reads NON MAPPES dans ${sample.bam} : "
	@${PACKAGES}/samtools/samtools view -f 4 -c ${sample.bam}
	@echo -n "Nombre de Reads NON MAPPES dans $< : "
	@cat $<

countalign.txt : countalign
	cat  ${sample.bam} | LD_LIBRARY_PATH=${htslib.dir} ./countalign -o $@ > test.bam
	
countalign : countalign.o
	${CC} -o $@  ${LDFLAGS} $< ${sam.libs}

countalign.o : countalign.c debug.h
	${CC} -c -o $@ ${CFLAGS} $<
	
test20150317: test20150317.c
	${CC} -o $@ $<

clean:
	rm -f countalign
	rm -rf $(dir $(TARGETS))


	

