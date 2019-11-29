
# kmcEx
K-mers along with their frequency have served as an elementary building block for error correction, repeat detection, multiple sequence alignment, genome assembly, etc., attracting intensive studies in k-mer counting. However, the output of k-mer counters itself is large; very often, it is too large to fit into main memory, leading to highly narrowed usability. We introduce a novel idea of encoding k- mers as well as their frequency, achieving good memory saving and retrieval efficiency. Specifically, we propose a Bloom Filter-like data structure to encode counted k-mers by coupled-bit arrays— one for k-mer representation and the other for frequency encoding. Experiments on five real data sets show that the average memory- saving ratio on all 31-mers is as high as 13.81 as compared with raw input, with 7 hash functions. At the same time, the retrieval time complexity is well controlled (effectively constant), and the false-positive rate is decreased by two orders of magnitude.

# installation 
kmcEx is based on **C++11**，and the installation step is very simple, in the kmodel main directory run 'make'  command, then you can get the executable kmodel.
```
make # run in the main directory 
```
# usage
how to use  the executable kmcEx for testing
```
1. USAGE
     kmcEx [options] <input_file_name> <output_file_name> <working_directory>
     kmcEx [options] <@input_file_names> <output_file_name> <working_directory>
2. OPTIONS
     1) REQUIRED
        input_file_name    - single file in FASTQ format (gziped or not)
        @input_file_names  - file name with list of input files in FASTQ format (gziped or not)
        working_directory  - save temporary files
     2) OPTIONAL
        -k<len>            - k-mer length (default: 31) 
        -t<value>          - total number of threads (default: 4)
        -ci<value>         - exclude k-mers occurring less than <value> times (default: 1)
        -cs<value>         - maximal value of a counter (default: 1023)
        -nh<value>         - number of hash (default: 7)
        -nb<value>         - number of bit array (default: 5)
3. EXAMPLES
     kmcEx -k31 -nh7 -nb5  rs.fastq rs.res /tmp
     kmcEx -k31 -nh7 -nb5  @rs.lst rs.res /tmp
```

# usage demo

```
#!/bin/bash
kn=31
out_dir=/tmp 
ci=1
cs=1023
./kmcEx -ci${ci} -cs${cs} -k${kn} -t8 @sa.lst $out_dir/SA_k${kn}_f${ci}-${cs}.res /tmp
./kmcEx -ci${ci} -cs${cs} -k${kn} -t8 @rs.lst $out_dir/RS_k${kn}_f${ci}-${cs}.res /tmp
./kmcEx -ci${ci} -cs${cs} -k${kn} -t8 @hc14.lst $out_dir/HC14_k${kn}_f${ci}-${cs}.res /tmp
./kmcEx -ci${ci} -cs${cs} -k${kn} -t8 @bi.lst $out_dir/BI_k${kn}_f${ci}-${cs}.res /tmp
./kmcEx -ci${ci} -cs${cs} -k${kn} -t8 @na12878_1.lst $out_dir/NA12878_k${kn}_f${ci}-${cs}.res /tmp
```
the file **rs.lst**  is the list of pathes of fastq files, one line per file, e.g.:
```
/share/share/data/NGS/GAGE/RS/rs_frag_1.fastq
/share/share/data/NGS/GAGE/RS/rs_frag_2.fastq
```
after running the command, in the $out_dir，you will get kmc_database files:
```
RS_k31_f1-1023.res.kmc_pre
RS_k31_f1-1023.res.kmc_suf
```
the model saved  in  /tmp/RS_k31_f1-1023.res  (working_directory)，i.e.:

```
header  km.bin  rest.bin
```


# API usage
create a kmodel and save it to the disk，since the number of kmers having count of 1 (counter=1) is very huge，we distinguish this two situations that counter>1 and counter>=1. you can pass the parameter 'ci' to get different model (if ci=1, then the model uses counter>=1; else the counter>1 is used)
1) create a model, and save it to a disk
```c
#include "kmodel.hpp"
//some arguments
int n_hash = 7, n_bit = 4, ci=1, cs=1023;
//kmc_database is the kmc database file,such as 'RS_k31_f1-1000.res'
string kmc_database = "/tmp/RS_k31_f2-1000.res";
//get a model with ci... 
KModel* kmodel = get_model(ci, cs, n_hash, n_bit);
//initialize the model
kmodel->init_KModel(kmc_database);
//save model to an existing directory "model_dir",
kmodel->save_model("/tmp/rs_f2-1000_model");
```

2) load a model from  a disk, and retrieve the occurrence of kmers
```c
#include "kmodel.hpp"
//some arguments
string save_dir = "/tmp/rs_f2-1000_model";
KModel* kmodel = get_model(save_dir);
//query one kmer 
string kmer = ...;
int occ = kmodel->kmer_to_occ(kmer)
//query multiple kmer with multiple threads
vector<string> kmer_v = ...;
vector<int> out = kmodel->kmer_to_occ(kmer_v);
```

# reference
kmc: https://github.com/refresh-bio/KMC
