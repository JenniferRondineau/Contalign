# masterm2

Contalign is a software allowing to read BAM files and looks for unmapped reads. 
Contalign map unmapped reads against a large reference of contaminants. It 
releases a report containing the number of unmapped reads by sample and the
potential contaminants. 

# Requirements / Dependencies
Contalign depends on the samtools library <http://samtools.sourceforge.net>, 
the BWA library <http://bio-bwa.sourceforge.net>, and HTSlib <http://www.htslib.org>.
Building him requires samtools and BWA development files to be installes on the 
build machine. 

# Download & install

$ git clone "https://github.com/JenniferRondineau/masterm2.git"

# Compilation 
$ cd masterm2/src
$ make 


