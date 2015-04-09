include ../src/config.mk
sample.bam= sample.bam
reference=UniVec.fa

all: test1 test2



test: contalign.txt ${sample.bam}
	@echo -n "Nombre de Reads dans ${sample.bam} : "
	@${samtools.dir}/samtools view -c ${sample.bam}
	@echo -n "Nombre de Reads dans test.bam : "
	@${samtools.dir}/samtools view test.bam | wc -l
	@echo -n "Nombre de Reads NON MAPPES dans ${sample.bam} : "
	@${samtools.dir}/samtools view -f 4 -c ${sample.bam}
	@echo -n "Nombre de Reads NON MAPPES dans $< : "
	@cat $<

contalign.txt : contalign ${sample.bam} $(addsuffix .bwt,${reference})
	cat ${sample.bam} | LD_LIBRARY_PATH=${htslib.dir} ./contalign -r ${reference} -o $@ -s save.fastq > test.bam
	${samtools.dir}/samtools view -c test.bam

${sample.bam} : 
	wget -O - ${WGETFLAGS} "ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/data/HG00096/exome_alignment/HG00096.chrom20.ILLUMINA.bwa.GBR.exome.20120522.bam" |\
		./samtools-${samtools.version}/samtools view -h - | head -n 50000 |\
		./samtools-${samtools.version}/samtools view -Sb -o $@ -

$(addsuffix .bwt,${reference}): ${reference} bwa-${bwa.version}/libbwa.a
	bwa-${bwa.version}/bwa index $<

${reference} : 
	wget -O $@ "ftp://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec" 
	


test1: 	../src/contalign
	wget ${WGETFLAGS} -O - "ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/data/HG00096/exome_alignment/HG00096.chrom20.ILLUMINA.bwa.GBR.exome.20120522.bam" |\
	../src/samtools-${samtools.version}/samtools view -h - | head -n 1000 | LD_LIBRARY_PATH=../src/htslib-${hts.version} $< -r ../src/UniVec.fa -o test.txt > /dev/null
	stat test.txt
	@echo -n "sample name + number of unmapped read : " 
	@tail test.txt | grep "HG" 
	

test2: 	../src/contalign
	wget -O - "ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/data/HG00096/exome_alignment/HG00096.chrom20.ILLUMINA.bwa.GBR.exome.20120522.bam" | ../src/samtools-${samtools.version}/samtools view -h - | head -n 10000 | LD_LIBRARY_PATH=../src/htslib-${hts.version} $< -r ../src/UniVec.fa -o test.txt > /dev/null
	stat test.txt
	@echo -n "sample name + number of unmapped read : " 
	@tail test.txt | grep "HG" 
	

../src/contalign: 
	(cd ../src ; make)

clean:
	(cd ../src ; make)
