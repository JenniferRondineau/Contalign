include ../src/config.mk
sample.bam= sample.bam
reference_UniVec=UniVec.fa
reference=reference.fa
EUTILS=http://eutils.ncbi.nlm.nih.gov
EFETCH=${EUTILS}/entrez/eutils/efetch.fcgi

all: 
	@echo "select a test"

test_report: contalign.txt ${sample.bam}
	@echo -n "Nombre de Reads dans ${sample.bam} : "
	@../src/samtools-${samtools.version}/samtools view -c ${sample.bam}
	@echo -n "Nombre de Reads dans test.bam : "
	@../src/samtools-${samtools.version}/samtools view test.bam | wc -l
	@echo -n "Nombre de Reads NON MAPPES dans ${sample.bam} : "
	@../src/samtools-${samtools.version}/samtools view -f 4 -c ${sample.bam}
	@echo -n "Nombre de Reads NON MAPPES dans $< : "
	@cat $<

contalign.txt : ../src/contalign ${sample.bam} $(addsuffix .bwt,${reference_UniVec})
	cat ${sample.bam} | $< -r ${reference_UniVec} -o $@ -s save.fastq > test.bam
	../src/samtools-${samtools.version}/samtools view -c test.bam

${sample.bam} : 
	wget -O - ${wgetFLAGS} "ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/data/HG00096/exome_alignment/HG00096.chrom20.ILLUMINA.bwa.GBR.exome.20120522.bam" |\
	../src/samtools-${samtools.version}/samtools view -h - | head -n 50000 |\
	../src/samtools-${samtools.version}/samtools view -Sb -o $@ -

test1: 	../src/contalign $(addsuffix .bwt,${reference_UniVec}) ${sample.bam}
	cat ${sample.bam} |\
	../src/samtools-${samtools.version}/samtools view -h - | head -n 10000 |  $< -r ${reference_UniVec} -o test.txt -c fullreport.txt > /dev/null
	@stat test.txt
	@echo -n "sample name + number of unmapped read : " 
	@tail test.txt | grep "HG" 
	@cat test.txt
	@cat fullreport.txt

test2: 	../src/contalign $(addsuffix .bwt,${reference}) ${sample.bam}
	cat ${sample.bam} |\
	../src/samtools-${samtools.version}/samtools view -h - | head -n 10000 | $< -r ${reference} -o test.txt -c fullreport.txt > /dev/null
	@stat test.txt
	@echo -n "sample name + number of unmapped read : " 
	cat test.txt
	cat fullreport.txt
	
test3: ../src/contalign $(addsuffix .bwt,${reference}) ${sample.bam}
	$< -r ${reference} -o test.txt -s save.fastq -c fullreport.txt sample.bam > /dev/null
	@stat test.txt
	cat test.txt
	tail fullreport.txt
		
test_version: ../src/contalign
	$< -v
	
test_help: ../src/contalign
	$< -h
	
save.fastq: ../src/contalign ${sample.bam} $(addsuffix .bwt,${reference_UniVec})
	cat ${sample.bam} | $< -r ${reference_UniVec} -s $@ > test.bam
	head $@

$(addsuffix .bwt,${reference_UniVec}): ${reference_UniVec}
	../src/bwa-${bwa.version}/bwa index $<
	
${reference_UniVec} : 
	wget -O $@ "ftp://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec" 

$(addsuffix .bwt,${reference}): ${reference}
	../src/bwa-${bwa.version}/bwa index $<

${reference}:
	ftp_proxy=cache.ha.univ-nantes.fr:3128 wget -O $(addsuffix .tmp,$@) "${EFETCH}?db=nucleotide&id=108885074,127763381,13507739,482889214,374352002,749202356,44454,148342,148377268,15828471,751266762,749295773,397655102,193082772,21326584,240047135,26553452,269114774,281306659,291319937,294155300,294660180,304372805,308189587,313664890,313678134,319776738,319801362,325972867,325989358,330370665,331028054,331703020,339320528,344204770,358356470,430768498,47118301,47458835,479336529,526118307,54019969,593776974,593777399,6572414,71893359,71894025,730966571,83319253,90592637,9628880,9629357&rettype=fasta" && \
	mv  $(addsuffix .tmp,$@) $@
	

../src/contalign: 
	(cd ../src ; make)
	

clean:
	(cd ../src ; make)
