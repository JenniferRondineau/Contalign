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
	wget -O - ${WGETFLAGS} "ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/data/HG00096/exome_alignment/HG00096.chrom20.ILLUMINA.bwa.GBR.exome.20120522.bam" |\
	../src/samtools-${samtools.version}/samtools view -h - | head -n 50000 |\
	../src/samtools-${samtools.version}/samtools view -Sb -o $@ -

test1: 	../src/contalign $(addsuffix .bwt,${reference_UniVec}) ${sample.bam}
	cat ${sample.bam} |\
	../src/samtools-${samtools.version}/samtools view -h - | head -n 10000 |  $< -r ${reference_UniVec} -o test.txt > /dev/null
	stat test.txt
	@echo -n "sample name + number of unmapped read : " 
	@tail test.txt | grep "HG" 

test2: 	../src/contalign $(addsuffix .bwt,${reference}) ${sample.bam}
	cat ${sample.bam} |\
	../src/samtools-${samtools.version}/samtools view -h - | head -n 10000 |  $< -r ${reference} -o test.txt > /dev/null
	stat test.txt
	@echo -n "sample name + number of unmapped read : " 
	cat test.txt
		
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

$(reference.fa):
	wget -O - "${EFETCH}?db=nucleotide&id=730966571&rettype=fasta" > $@  # gi|730966571|ref|NC_025822.1| Vibrio phage phi-A318, complete genome
	wget -O - "${EFETCH}?db=nucleotide&id=593776974&rettype=fasta" >> $@ # gi|593776974|ref|NC_023715.1| Yersinia phage Yep-phi, complete genome
	wget -O - "${EFETCH}?db=nucleotide&id=593777399&rettype=fasta" >> $@ # gi|593777399|ref|NC_023718.1| Pseudomonas phage phi_Pto-bp6g, complete genome
	wget -O - "${EFETCH}?db=nucleotide&id=479336529&rettype=fasta" >> $@ # gi|479336529|ref|NC_021062.1| Pseudomonas phage Phi-S1, complete genome
	wget -O - "${EFETCH}?db=nucleotide&id=526118307&rettype=fasta" >> $@ # gi|526118307|ref|NC_021784.1| Phage phi OH2 DNA, complete genome
	wget -O - "${EFETCH}?db=nucleotide&id=358356470&rettype=fasta" >> $@ # gi|358356470|ref|NC_016162.1| Vibrio phage VCY-phi, complete genome
	wget -O - "${EFETCH}?db=nucleotide&id=331028054&rettype=fasta" >> $@ # gi|331028054|ref|NC_015466.1| Roseobacter phage RDJL Phi 1, complete genome
	wget -O - "${EFETCH}?db=nucleotide&id=90592637&rettype=fasta" >> $@ # gi|90592637|ref|NC_007917.1| Clostridium phage phi CD119, complete genome
	wget -O - "${EFETCH}?db=nucleotide&id=281306659&rettype=fasta" >> $@ # gi|281306659|ref|NC_013638.1| Pseudomonas phage phi-2, complete genome
	wget -O - "${EFETCH}?db=nucleotide&id=9629357&rettype=fasta" >> $@ # gi|9629357|ref|NC_001802.1| Human immunodeficiency virus 1, complete genome
	wget -O - "${EFETCH}?db=nucleotide&id=9628880&rettype=fasta" >> $@ # gi|9628880|ref|NC_001722.1| Human immunodeficiency virus 2, complete genome 
	wget -O - "${EFETCH}?db=nucleotide&id=148377268&rettype=fasta" >> $@ # gi|148377268|ref|NC_009497.1|_Mycoplasma_agalactiae_PG2_chromosome, complete sequence
	wget -O - "${EFETCH}?db=nucleotide&id=291319937&rettype=fasta" >> $@ # gi|291319937|ref|NC_013948.1|_Mycoplasma_agalactiae_5632_chromosome, complete sequence
	wget -O - "${EFETCH}?db=nucleotide&id=193082772&rettype=fasta" >> $@ # gi|193082772|ref|NC_011025.1|_Mycoplasma_arthritidis_158L3-1 complete genome
	wget -O - "${EFETCH}?db=nucleotide&id=339320528&rettype=fasta" >> $@ # gi|339320528|ref|NC_015725.1|_Mycoplasma_bovis_Hubei-1, complete genome
	wget -O - "${EFETCH}?db=nucleotide&id=313678134&rettype=fasta" >> $@ # gi|313678134|ref|NC_014760.1|_Mycoplasma_bovis_PG45_clone_MU_clone A2, complete genome
	wget -O - "${EFETCH}?db=nucleotide&id=83319253&rettype=fasta" >> $@ # gi|83319253|ref|NC_007633.1|_Mycoplasma_capricolum_subsp._capricolum ATCC 27343, complete genome
	wget -O - "${EFETCH}?db=nucleotide&id=240047135&rettype=fasta" >> $@ # gi|240047135|ref|NC_012806.1|_Mycoplasma_conjunctivae, complete genome
	wget -O - "${EFETCH}?db=nucleotide&id=294155300&rettype=fasta" >> $@ # gi|294155300|ref|NC_014014.1|_Mycoplasma_crocodyli_MP145, complete genome
	wget -O - "${EFETCH}?db=nucleotide&id=308189587&rettype=fasta" >> $@ # gi|308189587|ref|NC_014552.1|_Mycoplasma_fermentans_JER, complete genome
	wget -O - "${EFETCH}?db=nucleotide&id=319776738&rettype=fasta" >> $@ # gi|319776738|ref|NC_014921.1|_Mycoplasma_fermentans_M64, complete genome
	wget -O - "${EFETCH}?db=nucleotide&id=294660180&rettype=fasta" >> $@ # gi|294660180|ref|NC_004829.2|_Mycoplasma_gallisepticum_str._R(low), complete genome
	wget -O - "${EFETCH}?db=nucleotide&id=108885074&rettype=fasta" >> $@ # gi|108885074|ref|NC_000908.2|_Mycoplasma_genitalium_G37, complete genome
	wget -O - "${EFETCH}?db=nucleotide&id=319801362&rettype=fasta" >> $@ # gi|319801362|emb|FR773153.2|_Mycoplasma_haemofelis_str._Langford_1, complete genome sequence
	wget -O - "${EFETCH}?db=nucleotide&id=269114774&rettype=fasta" >> $@ # gi|269114774|ref|NC_013511.1|_Mycoplasma_hominis_ATCC_23114 chromosome complete genome
	wget -O - "${EFETCH}?db=nucleotide&id=54019969&rettype=fasta" >> $@ # gi|54019969|ref|NC_006360.1|_Mycoplasma_hyopneumoniae_232, complete genome
	wget -O - "${EFETCH}?db=nucleotide&id=72080342&rettype=fasta" >> $@ # gi|72080342|ref|NC_007332.1|_Mycoplasma_hyopneumoniae_7448, complete genome
	wget -O - "${EFETCH}?db=nucleotide&id=71893359&rettype=fasta" >> $@ # gi|71893359|ref|NC_007295.1|_Mycoplasma_hyopneumoniae_J, complete genome
	wget -O - "${EFETCH}?db=nucleotide&id=304372805&rettype=fasta" >> $@ # gi|304372805|ref|NC_014448.1|_Mycoplasma_hyorhinis_HUB-1, complete genome
	wget -O - "${EFETCH}?db=nucleotide&id=313664890&rettype=fasta" >> $@ # gi|313664890|ref|NC_014751.1|_Mycoplasma_leachii_PG50 clone MU clone A8, complete genome
	wget -O - "${EFETCH}?db=nucleotide&id=47458835&rettype=fasta" >> $@ # gi|47458835|ref|NC_006908.1|_Mycoplasma_mobile 163K complete genome
	wget -O - "${EFETCH}?db=nucleotide&id=330370665&rettype=fasta" >> $@ # gi|330370665|ref|NC_015407.1|_Mycoplasma_mycoides_subsp._capri LC str. 95010 plasmid, complete genome
	wget -O - "${EFETCH}?db=nucleotide&id=331703020&rettype=fasta" >> $@ # gi|331703020|ref|NC_015431.1|_Mycoplasma_mycoides_subsp._capri_LC_str. 95010 chromosome, complete sequence
	wget -O - "${EFETCH}?db=nucleotide&id=127763381&rettype=fasta" >> $@ # gi|127763381|ref|NC_005364.2|_Mycoplasma_mycoides_subsp._mycoides_SC_str._PG1 chromosome, complete genome
	wget -O - "${EFETCH}?db=nucleotide&id=26553452&rettype=fasta" >> $@ # gi|26553452|ref|NC_004432.1|_Mycoplasma_penetrans_HF-2_DNA, complete genome
	wget -O - "${EFETCH}?db=nucleotide&id=13507739&rettype=fasta" >> $@ # gi|13507739|ref|NC_000912.1|_Mycoplasma_pneumoniae_M129_chromosome, complete genome
	wget -O - "${EFETCH}?db=nucleotide&id=15828471&rettype=fasta" >> $@ # gi|15828471|ref|NC_002771.1|_Mycoplasma_pulmonis_UAB_CTIP, complete genome
	wget -O - "${EFETCH}?db=nucleotide&id=344204770&rettype=fasta" >> $@ # gi|344204770|ref|NC_015946.1|_Mycoplasma_putrefaciens_KS1, complete genome
	wget -O - "${EFETCH}?db=nucleotide&id=325972867&rettype=fasta" >> $@ # gi|325972867|ref|NC_015155.1|_Mycoplasma_suis_str._Illinois, complete genome
	wget -O - "${EFETCH}?db=nucleotide&id=325989358&rettype=fasta" >> $@ # gi|325989358|ref|NC_015153.1|_Mycoplasma_suis_KI3806 complete genome
	wget -O - "${EFETCH}?db=nucleotide&id=71894025&rettype=fasta" >> $@ # gi|71894025|ref|NC_007294.1|_Mycoplasma_synoviae 53, complete genome
	wget -O - "${EFETCH}?db=nucleotide&id=21326584&rettype=fasta" >> $@ # gi|21326584|ref|NC_003977.1| Hepatitis B virus, complete genome
	wget -O - "${EFETCH}?db=nucleotide&id=47118301&rettype=fasta" >> $@ # gi|47118301|dbj|BA000007.2| Escherichia coli O157:H7 str. Sakai DNA, complete genome
	wget -O - "${EFETCH}?db=nucleotide&id=6572414&rettype=fasta" >> $@ # gi|6572414|emb|Z86099.2| Herpes simplex virus type 2 (strain HG52), complete genome
	wget -O - "${EFETCH}?db=nucleotide&id=430768498&rettype=fasta" >> $@ # gi|430768498|dbj|AB618031.1| Herpes simplex virus (type 1 /strain RH2) DNA, nearly complete genome
	wget -O - "${EFETCH}?db=nucleotide&id=430768498&rettype=fasta" >> $@ # gi|330360|gb|M10593.1|HS4ENVGP Epstein Barr virus major outer envelope glycoprotein genes gp350 and gp220, complete cds
	

../src/contalign: 
	(cd ../src ; make)
	

clean:
	(cd ../src ; make)
