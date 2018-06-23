					kmcEx: Kmer count Model tool

VERSION: 1.0
DATE   : April.10, 2018

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
3. EXAMPLES" << endl;
	kmcEx -k27 -nh7 -nb4  NA19238.fastq NA.res /tmp
	kmcEx -k27 -nh7 -nb4  @files.lst NA.res /tmp
	

4.API USAGE
4.1 create a kmodel and save it to the disk
	#include "kmodel.hpp"
	int n_hash=6,n_bit=4; //some parameters
	OccuBin* occu_bin = new OccuBin(n_hash);
	occu_bin->init_occ_count(input_file); //input_file is the kmc database file
	KModel* kmodel = new KModel(occu_bin, n_bit);
	kmodel->init_KModel(input_file); //initialize the model
	kmodel->save_model(model_dir); //save model to the directory "model_dir"
	int occ=kmodel->kmer_to_occ(kmer) //get the occurrence of the kmer

4.2 load a model from the disk
	#include "kmodel.hpp"
	KModel* kmodel = new KModel();
	kmodel->load_model(model_dir); //load model from disk
	int occ=kmodel->kmer_to_occ(kmer) //get the occurrence of the kmer


