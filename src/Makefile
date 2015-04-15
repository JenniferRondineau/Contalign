include config.mk

.PHONY:all clean


CC=gcc
samtools.dir=samtools-${samtools.version}
htslib.dir=${samtools.dir}/htslib-${htslib.version}
bwa.dir=bwa-${bwa.version}
CFLAGS= -Wall -O3 -I ${samtools.dir} -I ${bwa.dir}  -I ${htslib.dir}
LDFLAGS=-L ${samtools.dir} -L ${htslib.dir} -L ${bwa.dir}
sam.libs=  -lm -lz -lpthread
libraries=${bwa.dir}/libbwa.a ${samtools.dir}/libbam.a ${htslib.dir}/libhts.a


all: contalign  

contalign : contalign.o runcontalign.o $(libraries)
	${CC} -o $@ ${DFLAGS} $^ ${sam.libs} 

contalign.o : contalign.c contalign.h debug.h $(libraries)
	${CC} -c -o $@ ${CFLAGS} $< 

runcontalign.o : runcontalign.c contalign.h debug.h $(libraries)
	${CC} -c -o $@ ${CFLAGS} $< 

${bwa.dir}/libbwa.a:
	rm -rf bwa-${bwa.version}
	wget ${WGETFLAGS} -O ${bwa.version}.tar.gz "https://github.com/lh3/bwa/archive/${bwa.version}.tar.gz"
	tar xfvz ${bwa.version}.tar.gz
	rm ${bwa.version}.tar.gz
	(cd $(dir $@) ; make)
	
${samtools.dir}/libbam.a: 
	rm -rf samtools-${samtools.version}
	wget  ${WGETFLAGS} -O samtools-${samtools.version}.tar.bz2 "https://github.com/samtools/samtools/releases/download/${samtools.version}/samtools-${samtools.version}.tar.bz2"
	tar xvjf samtools-${samtools.version}.tar.bz2
	rm samtools-${samtools.version}.tar.bz2
	(cd $(dir $@) ; make)


clean:
	rm -f contalign ${htslib.dir} ${samtools.dir} ${bwa.dir}


