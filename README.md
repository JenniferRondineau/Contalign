# masterm2

Contalign is a software allowing to read BAM files and looks for unmapped reads. 
Contalign map unmapped reads against a large reference of contaminants. It 
releases a report containing the number of unmapped reads by sample and the
potential contaminants. 

# Requirements 

Contalign requires wget <https://www.gnu.org/software/wget/>

# Download & install

```shell
$ git clone "https://github.com/JenniferRondineau/masterm2.git"
```

# Compilation 


```shell
$ cd masterm2/src
$ make
```

to compile contalign.


# Authors

Jennifer Rondineau, bioinformatics student (M1) wrote this code.
 
Heng Li from the Sanger Institute wrote the C version of samtools. 

Bob Handsaker from the Broad Institute implemented the BGZF library. 

Code is also inspired from Heng Li's https://github.com/lh3/bwa/blob/master/example.c

# Report errors

# license
The MIT License (MIT)

Copyright (c) 2015 JenniferRondineau

http://samtools.sourceforge.net/
Authors: Heng Li, Bob Handsaker, Jue Ruan, Colin Hercus, Petr Danecek

https://github.com/lh3/bwa/blob/master/example.c
Authors: Heng Li's

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


