****
kmap
****

Overview
========
The repository provides code to represent sequences files in fasta format in terms of their kmer counts.

### Author
Jean Disset (j.disset@gmail.com)
 
Installation
============

### Requirements 
- gcc >= 9.2.0 
- cmake >= 3.14.1
- libBz2 (available from most package managers)

Kmap depends on rocksDB, which will be automatically pulled from the official repository and compiled when running cmake. RocksDB itslef requires libBz2, which should be available from most system package managers.

### Instructions

Kmap can be built using the included CMakeLists.txt. 
After cloning the repo, run the folling code ::

```
    cd
    mkdir build
    cd build
    cmake ..
    make -j4
```
Processing sequences
====================
kmap starts with a set of sequence files and a lag k.
For each k-mer seen in any of these files, kmap counts how often the k-mer is followed by each letter of a specified alphabet (dna, rna or protein) in each sequence file.
Prior to counting, sequences are prepended with k start symbols [ and suffixed with a ]; kmap includes these symbols in its counting.
These counts are represented in a matrix whose rows represent each sequence file and whose columns represent each letter of the alphabet (in alphabetical order, ] last).
kmap will then output a file of semicolon-separated list of values.
The file will include each kmer seen in the data, followed by the indices and values of the counts matrix in sparse format.

For example, kmapping ::
```
    TAT
    TTT
```
with the two sequences in separate files with k=2 results in ::
```
    [[;[[0,3],[1,3]];[1,1]
    [T;[[0,0],[0,3]];[1,1]
    TA;[[0,3]];[1]
    AT;[[0,4]];[1]
    TT;[[1,3],[1,4]];[1,1]
```

`AT;[[0,4]];[1]` means that the 2-mer `AT` was only seen in the 0-th dataset followed by the letter of the DNA alphabet indexed at 4, the stop symbol ], and this transition was only seen once.

Processing sequence files
=========================
Running `kmap --help` will describe how to use kmap.
The `t, n, c` flags control mutithreading.
The `a` flag takes the alphabet of the sequence, either 'dna', 'rna' or 'protein'.
The `k` flag describes the length of the kmers to be represented or "k" in "k-mer".
The input to kmap, specified by the `i` flag is a text file with the locations of separate sequence files relative to the input file on each line.
The output location is described by the `o` flag.

The processing is done in two steps, first with kmap, then kconv.
As an example, to process protein sequences in sequences1.fa and sequences2.fa at k=3 into kmap.csv, one first must create a file files_list.txt ::
    sequences1.fa
    sequences2.fa
then run kmap ::
    kmap -k 3 -i files_list.txt -o temp_folder -a protein
and finally run kconv ::
    kconv -i temp_folder -o kmap.csv
