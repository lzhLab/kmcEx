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

}Params;


void read_me() {
	cout << "----------------------------------------------------------------------" << endl;
	cout << "           kmcEx: counted k-mer encoding & decoding                   " << endl;
	cout << "----------------------------------------------------------------------" << endl;
	cout << "VERSION: 1.3" << endl;
	cout << "DATE   : Apr 2nd, 2019" << endl;
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
			Params.t = atoi(&argv[i][2]);
		// k-mer length
		else if (strncmp(argv[i], "-k", 2) == 0)
			Params.k = atoi(&argv[i][2]);
		//hsah number
		else if (strncmp(argv[i], "-nh", 3) == 0)
			Params.num_hash = atoi(&argv[i][3]);
		//bit number
		else if (strncmp(argv[i], "-nb", 3) == 0)
			Params.num_bit = atoi(&argv[i][3]);
		// Minimum counter threshold
		else if (strncmp(argv[i], "-ci", 3) == 0)
			Params.ci = atoi(&argv[i][3]);
		//maximal value of a counter
		else if (strncmp(argv[i], "-cs", 3) == 0)
			Params.cs = atoi(&argv[i][3]);
	}

	if (argc - i < 3)
		return false;
	i = argc - 3;
	Params.input_file_name = string(argv[i++]);
	Params.output_file_name = string(argv[i++]);
	Params.working_directory = string(argv[i++]);
	if (Params.input_file_name.length() <= 0) {
		cout << "input_file_name is needed..!";
		return false;
	}
	if (Params.output_file_name.length() <= 0) {
		cout << "output_file_name is needed..!";
		return false;
	}
	if (Params.input_file_name.length() <= 0) {
		cout << "working_directory is needed..!";
		return false;
	}
	return true;
}



int run(int argc, char* argv[]) {
	bool flag = parse_parameters(argc, argv);
	if (!flag) {
		read_me();
		return -1;
	}
	char cmd[1024];
	sprintf(cmd, "./kmc_api/kmc -k%d -t%d -ci%d -cs%d %s %s %s", Params.k, Params.t, Params.ci, Params.cs,
		Params.input_file_name.c_str(), Params.output_file_name.c_str(), Params.working_directory.c_str());
	cout << cmd << endl;
	system(cmd);
	cout << endl;

	//cout << Params.output_file_name << endl;
	cout << "kmcEx status:" << endl;
	KModel* kmodel = get_model(Params.ci, Params.cs, Params.num_hash, Params.num_bit);
	kmodel->init_KModel(Params.output_file_name);
	kmodel->show_kmodel_info();
	//save model to the disk
	string save_dir = Params.working_directory + "/" + Tools::get_file_name(Params.output_file_name);
	system(("mkdir -p " + save_dir).c_str()); //new a dirctory - -;
	kmodel->save_model(save_dir);
	//kmodel->test_model(Params.output_file_name);
	cout << "\n-------\n";
	cout << "   kmcEx model is successfully saved in  :  " << save_dir << "\n\n";
}


void test_kmer_kmodel(string kmer, string model_dir = "save_model") {
	KModel* kmodel = get_model(2, Params.num_hash, Params.num_bit);
	kmodel->load_model(model_dir);
	cout << kmodel->kmer_to_occ(kmer) << endl;
}



int _tmain(int argc, char* argv[]) {
	run(argc, argv);
	return 0;
}
