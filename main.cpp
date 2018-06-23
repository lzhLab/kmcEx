#include "kmodel.hpp"
#include <stdlib.h>


struct KParams {
	int k = 25; //kmer
	int num_hash = 7; // hashfunction number
	int num_bit = 4; //the number of bit array, default is 4
	int ci = 2; //exclude k - mers occurring less than <value> times(default: 2)
	int cx = 1000000000; //-cx<value> -exclude k - mers occurring more of than <value> times(default: 1e9)
	int cs = 1001; //-cs<value> -maximal value of a counter(default: 1001,exclude 1001)
	int t = 0; // -t<value> - total number of threads (default: no. of CPU cores)
	string input_file_name;
	string output_file_name;
	string working_directory = "/tmp";

};
KParams Params;

void read_me() {
	cout << "----------------------------------------------------------------------" << endl;
	cout << "           kmcEx: Kmer count Model tool" << endl;
	cout << "----------------------------------------------------------------------" << endl;
	cout << "VERSION: 1.01" << endl;
	cout << "DATE   : April.10, 2018" << endl;
	cout << "----------------------------------------------------------------------" << endl << endl;
	cout << "1. USAGE" << endl;
	cout << "     kmcEx [options] <input_file_name> <output_file_name> <working_directory>" << endl;
	cout << "     kmcEx [options] <@input_file_names> <output_file_name> <working_directory>" << endl;
	cout << "2. OPTIONS" << endl;
	cout << "	1) REQUIRED" << endl;
	cout << "		input_file_name - single file in FASTQ format (gziped or not)" << endl;
	cout << "		@input_file_names - file name with list of input files in FASTQ format (gziped or not)" << endl;
	cout << "		working_directory -save some temporary files" << endl;
	cout << "	2) OPTIONAL" << endl;
	cout << "		-k<len> - k-mer length (k from 10 to 256; default: 25" << endl;
	cout << "		-t<value> - total number of threads (default: no. of CPU cores)" << endl;
	cout << "		-ci<value> - exclude k-mers occurring less than <value> times (default: 2)" << endl;
	cout << "		-cx<value> - exclude k-mers occurring more of than <value> times (default: 1e9)" << endl;
	cout << "		-cs<value> - maximal value of a counter(default: 1001,exclude)" << endl;
	cout << "		-nh<value> - number of hash (default: 7)" << endl;
	cout << "		-nb<value> - number of bit array (default: 4)" << endl;
	cout << "3. EXAMPLES" << endl;
	cout << "	kmcEx -k27 -nh7 -nb4  rs.fastq rs.res /tmp" << endl;
	cout << "	kmcEx -k27 -nh7 -nb4  @rs.lst rs.res /tmp" << endl << endl;
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
		// Maximum counter threshold
		else if (strncmp(argv[i], "-cx", 3) == 0)
			Params.cx = atoi(&argv[i][3]);
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

void test_fp(KModel* kmodel, string input_file) {
	cout << "testing......." << endl;
	string kmer;
	float counter;
	int occ, t_bin, bin, appearance; //appearance:how many time this kmer appearance in these bit array
	double false_positive = 0, retrieval_error = 0, left_map = 0, last_f = 0;
	CKMCFile kmer_data_base;
	kmer_data_base.OpenForListing(input_file);
	uint32 _kmer_length = kmer_data_base.KmerLength();
	CKmerAPI kmer_object(_kmer_length);

	while (kmer_data_base.ReadNextKmer(kmer_object, counter)) {
		kmer = kmer_object.to_string().c_str();
		if (kmodel->get_kmer_from_map(kmer) > -1) { //can be found in the map
			left_map++;
			continue;
		}
		if (counter >= C_MAX) {//skip the cs in the kmc2
			last_f++;
			continue;
		}
		occ = counter;
		bin = kmodel->get_occu_bin()->occ_to_bin(occ); //occ->bin
		t_bin = kmodel->kmer_to_bin(kmer, appearance);
		if (appearance > 1) {
			false_positive++;
		}
		if (t_bin != bin) //don't consider the result from map
			retrieval_error++;
	}
	kmer_data_base.Close();

	uint64_t bit_kmer_count = kmodel->get_kmer_count() - left_map - last_f;
	printf("false_positive:%0.6f, retrieval_error:%0.6f\n", false_positive / bit_kmer_count, retrieval_error / bit_kmer_count);
}


int _tmain(int argc, char* argv[]) {

	bool flag = parse_parameters(argc, argv);
	if (!flag) {
		read_me();
		return -1;
	}
	char cmd[1024];
	sprintf(cmd, "./kmc_api/kmc -k%d -t%d -ci%d -cx%d -cs%d %s %s %s", Params.k, Params.t, Params.ci, Params.cx, Params.cs,
		Params.input_file_name.c_str(), Params.output_file_name.c_str(), Params.working_directory.c_str());
	cout << cmd << endl;
	system(cmd);
	cout << endl;
	C_MAX = Params.cs;
	OccuBin* occu_bin = new OccuBin(Params.num_hash);
	occu_bin->init_occ_count(Params.output_file_name);
	KModel* kmodel = new KModel(occu_bin, Params.num_bit);
	kmodel->show_header();
	kmodel->init_KModel(Params.output_file_name);
	kmodel->show_header2();
	test_fp(kmodel, Params.output_file_name);
	return 0;
}