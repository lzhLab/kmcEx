//-----------------------------------------------
// Copyright 2016 Guangxi University
// Written by Liang Zhao(S080011@e.ntu.edu.sg)
// Released under the GPL
//-----------------------------------------------
//
// kmcEx - Main driver program
// It's developed based on SGA, originally writen by Jared Simpson (js18@sanger.ac.uk)
//

#include "kmodel.hpp"
#include <stdlib.h>
#include "omp.h"


struct KParams {
	int k = 31; //kmer
	int num_hash = 7; // hashfunction number
	int num_bit = 5; //the number of bit array, default is 5
	int ci = 1; //exclude k - mers occurring less than <value> times(default: 1)
	int cs = 1023; //-cs<value> -maximal value of a counter(default: 1023)
	int t = 4; // -t<value> - total number of threads (default: 4)
	string input_file_name;
	string output_file_name;
	string working_directory = "/tmp";

}kParams;


void read_me() {
	cout << "----------------------------------------------------------------------" << endl;
	cout << "           kmcEx: counted k-mer encoding & decoding                   " << endl;
	cout << "----------------------------------------------------------------------" << endl;
	cout << "VERSION: 1.5" << endl;
	cout << "DATE   : Nov 2nd, 2019" << endl;
	cout << "----------------------------------------------------------------------" << endl << endl;
	cout << "1. USAGE" << endl;
	cout << "     kmcEx [options] <input_file_name> <output_file_name> <working_directory>" << endl;
	cout << "     kmcEx [options] <@input_file_names> <output_file_name> <working_directory>" << endl;
	cout << "2. OPTIONS" << endl;
	cout << "     1) REQUIRED" << endl;
	cout << "        input_file_name    - single file in FASTQ format (gziped or not)" << endl;
	cout << "        @input_file_names  - file name with list of input files in FASTQ format (gziped or not)" << endl;
	cout << "        working_directory  - save temporary files" << endl;
	cout << "     2) OPTIONAL" << endl;
	cout << "        -k<len>            - k-mer length (default: 31)" << endl;
	cout << "        -t<value>          - total number of threads (default: 4)" << endl;
	cout << "        -ci<value>         - exclude k-mers occurring less than <value> times (default: 1)" << endl;
	cout << "        -cs<value>         - maximal value of a counter (default: 1023)" << endl;
	cout << "        -nh<value>         - number of hash (default: 7)" << endl;
	cout << "        -nb<value>         - number of bit array (default: 5)" << endl;
	cout << "3. EXAMPLES" << endl;
	cout << "     kmcEx -k31 -nh7 -nb5  rs.fastq rs.res /tmp" << endl;
	cout << "     kmcEx -k31 -nh7 -nb5  @rs.lst rs.res /tmp" << endl << endl;
}

bool checkFileExist(string path) {
	ifstream fin(path);
	if (!fin)
		return false;
	return true;
}

bool parse_parameters(int argc, char *argv[]) {
	int i;
	if (argc < 4)
		return false;

	for (i = 1; i < argc; ++i)
	{
		if (argv[i][0] != '-')
			break;
		// Number of threads
		if (strncmp(argv[i], "-t", 2) == 0)
			kParams.t = atoi(&argv[i][2]);
		// k-mer length
		else if (strncmp(argv[i], "-k", 2) == 0)
			kParams.k = atoi(&argv[i][2]);
		//hsah number
		else if (strncmp(argv[i], "-nh", 3) == 0)
			kParams.num_hash = atoi(&argv[i][3]);
		//bit number
		else if (strncmp(argv[i], "-nb", 3) == 0)
			kParams.num_bit = atoi(&argv[i][3]);
		// Minimum counter threshold
		else if (strncmp(argv[i], "-ci", 3) == 0)
			kParams.ci = atoi(&argv[i][3]);
		//maximal value of a counter
		else if (strncmp(argv[i], "-cs", 3) == 0)
			kParams.cs = atoi(&argv[i][3]);
	}

	if (argc - i < 3)
		return false;
	i = argc - 3;
	kParams.input_file_name = string(argv[i++]);
	kParams.output_file_name = string(argv[i++]);
	kParams.working_directory = string(argv[i++]);
	if (kParams.input_file_name.length() <= 0) {
		cout << "input_file_name is needed..!";
		return false;
	}
	if (kParams.output_file_name.length() <= 0) {
		cout << "output_file_name is needed..!";
		return false;
	}
	if (kParams.input_file_name.length() <= 0) {
		cout << "working_directory is needed..!";
		return false;
	}
	return true;
}

//void test_create_api(string kmc_db,string save_model) {
//	//some arguments
//	int n_hash = 7, n_bit = 4, ci = 2, cs = 1023;
//	//kmc_database is the kmc database file,such as 'RS_k31_f1-1000.res'
//	string kmc_database = kmc_db;
//	//get a model with ci... 
//	KModel* kmodel = get_model(ci, cs, n_hash, n_bit);
//	//initialize the model
//	kmodel->init(kmc_database);
//	//save model to an existing directory "model_dir",
//	kmodel->save(save_model);
//}



int run(int argc, char* argv[]) {
	bool flag = parse_parameters(argc, argv);
	if (!flag) {
		read_me();
		return -1;
	}

	char cmd[1024];
	sprintf(cmd, "./kmc_api/kmc -k%d -t%d -ci%d -cs%d %s %s %s", kParams.k, kParams.t, kParams.ci, kParams.cs,
		kParams.input_file_name.c_str(), kParams.output_file_name.c_str(), kParams.working_directory.c_str());
	cout << cmd << endl;
	system(cmd);
	cout << endl;

	KModel* kmodel = get_model(kParams.ci, kParams.cs, kParams.num_hash, kParams.num_bit);
	kmodel->init(kParams.output_file_name);
	kmodel->show_kmodel_info();
	//save model to the disk
	string save_dir = kParams.working_directory + "/" + Tools::get_file_name(kParams.output_file_name);
	system(("mkdir -p " + save_dir).c_str()); //new a dirctory - -;
	kmodel->save(save_dir);
}




int _tmain(int argc, char* argv[]) {
	run(argc, argv);
	return 0;
}