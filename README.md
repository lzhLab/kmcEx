# kmcEx
The k-mers along with their frequency have been served as the elementary building block for error correction, repeat detection, multiple sequence alignment, genome assembly, etc, attracting intensive studies in k-mer counting. However, the result of k-mer counters itself is pretty large; very often, it is too large to ﬁt into the main memory, which in turn signiﬁcantly narrows down its usability. To overcome this bottleneck, we introduce a very novel idea of encoding k-mers as well as their frequency achieving ultra memory saving and retrieval eﬃcient. Precisely, we propose a Bloom Filter-like data structure to encode kmers as well as their frequency by coupling-bit arrays— one for k-mer representation and the other for frequency encoding. Experimental results conducted on ﬁve real data sets show that the compression ratio is as high as 7.8 via encoding. Besides, we have achieved constant time complexity in retrieving a k-mer as well as its frequency and two magnitudes smaller false positive rate.

# installation 
kmcEx is based on **C++11**，and the installation step is very simple, in the kmodel main directory run 'make'  command, then you can get the executable kmcEx.
```
make # run in the main directory 
```
# usage
how to use  the executable kmcEx to test
```
1. USAGE
    kmcEx [options] <input_file_name> <output_file_name> <working_directory>
	kmcEx [options] <@input_file_names> <output_file_name> <working_directory>

2. OPTIONS;
	1) REQUIRED
		input_file_name - single file in FASTQ format (gziped or not)
		@input_file_names - file name with list of input files in FASTQ format (gziped or not)
		working_directory -save some temporary files
	2) OPTIONAL
		-k<len> - k-mer length (k from 10 to 256; default: 25
		-t<value> - total number of threads (default: no. of CPU cores)
		-ci<value> - exclude k-mers occurring less than <value> times (default: 2)
		-cx<value> - exclude k-mers occurring more of than <value> times (default: 1e9)
		-cs<value> - maximal value of a counter(default: 1001,exclude)
		-nh<value> - number of hash (default: 6)
		-nb<value> - number of bit array (default: 4)
3. EXAMPLES
	kmcEx -k27 -nh7 -nb4  NA19238.fastq NA.res /tmp
	kmcEx -k27 -nh7 -nb4  @files.lst NA.res /tmp
```

# API usage
create a kmodel and save it to the disk
```C++
#include "kmodel.hpp"
int n_hash=6,n_bit=4; //some parameters
OccuBin* occu_bin = new OccuBin(n_hash);
//input_file is the kmc database file,such as 'NA.res'
occu_bin->init_occ_count(input_file); 
KModel* kmodel = new KModel(occu_bin, n_bit);
kmodel->init_KModel(input_file); //initialize the model
kmodel->save_model(model_dir); //save model to the directory "model_dir"
int occ=kmodel->kmer_to_occ(kmer) //get the occurrence of the kmer
```

load a model from  the disk
```
#include "kmodel.hpp"
KModel* kmodel = new KModel();
kmodel->load_model(model_dir); //load the model from the disk
int occ=kmodel->kmer_to_occ(kmer) //get the occurrence of the kmer
```