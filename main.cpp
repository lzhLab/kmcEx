#include "kmodel.hpp"
#include <stdlib.h>
#include "omp.h"


struct KParams {
	int k = 31; //kmer
	int num_hash = 7; // hashfunction number
	int num_bit = 4; //the number of bit array, default is 4
	int ci = 1; //exclude k - mers occurring less than <value> times(default: 1)
	int cs = 1023; //-cs<value> -maximal value of a counter(default: 1023)
	int t = 4; // -t<value> - total number of threads (default: 4)
	string input_file_name;
	string output_file_name;
	string working_directory = "/tmp";

}Params;


void read_me() {
	cout << "----------------------------------------------------------------------" << endl;
	cout << "           kmcEx: Kmer count Model tool" << endl;
	cout << "----------------------------------------------------------------------" << endl;
	cout << "VERSION: 1.02" << endl;
	cout << "DATE   : JULY.20, 2018" << endl;
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
	cout << "		-k<len> - k-mer length (default: 31)" << endl;
	cout << "		-t<value> - total number of threads (default: 4)" << endl;
	cout << "		-ci<value> - exclude k-mers occurring less than <value> times (default: 1)" << endl;
	cout << "		-cs<value> - maximal value of a counter (default: 1023)" << endl;
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
	cout << "==============test_fp==============" << endl;
	string kmer;
	float counter;
	int occ, t_bin, bin, idx = 0; //appearance:how many time this kmer appearance in these bit array
	double kmodel_false_positive = 0;
	CKMCFile kmer_data_base;
	kmer_data_base.OpenForListing(input_file);
	uint32 _kmer_length = kmer_data_base.KmerLength();
	CKmerAPI kmer_object(_kmer_length);
	uint64 kc_kmodel = kmodel->get_kmodel_kmer_count()-kmodel->get_left_kmer_count();
	//cout << kmodel->get_kmodel_kmer_count() << "-----" << kmodel->get_left_kmer_count() << endl;
	uint64_t kc_bf1 = kmodel->get_once_kmer_count();
	auto start = chrono::high_resolution_clock::now();
	while (kmer_data_base.ReadNextKmer(kmer_object, counter)) {
		kmer_buff[idx].kmer = kmer_object.to_string().c_str();
		kmer_buff[idx++].occ = counter;
		if (idx >= BLOACK_SIZE) {
			#pragma omp parallel for num_threads(THREAD_NUM) reduction(+:kmodel_false_positive) 
			for (int i = 0; i < idx; ++i) {
				if (kmodel->get_kmer_from_map(kmer_buff[i].kmer) > -1) { //can be found in the map
					continue;
				}
				//	continue;
				if (kmodel->kmer_to_bin(kmer_buff[i].kmer) != kmodel->get_occu_bin()->occ_to_bin(kmer_buff[i].occ)) //don't consider the result from map
					++kmodel_false_positive;
			}
			idx = 0;
		}
	}
	//handle the left kmers
	#pragma omp parallel for num_threads(THREAD_NUM) reduction(+:kmodel_false_positive) 
	for (int i = 0; i < idx; ++i) {
		if (kmodel->get_kmer_from_map(kmer_buff[i].kmer) > -1) { //can be found in the map
			continue;
		}
		if (kmodel->kmer_to_bin(kmer_buff[i].kmer) != kmodel->get_occu_bin()->occ_to_bin(kmer_buff[i].occ)) //don't consider the result from map
			++kmodel_false_positive;
	}
	kmer_data_base.Close();
	chrono::duration<double> dur = chrono::high_resolution_clock::now() - start;
	printf("test time:%.4f min\n", dur.count() / 60);
	cout << "FP_count_in_kmodel: " << kmodel_false_positive << endl;
	printf("kmodel_false_positive:%.3e\n\n",kmodel_false_positive / kc_kmodel);
}

void test_fp(KModelOne* kmodel, string input_file) {
	cout << "==============test_fp_one==============" << endl;
	string kmer;
	float counter;
	int occ, t_bin, bin, idx = 0; //appearance:how many time this kmer appearance in these bit array
	double kmodel_false_positive = 0, bf1_false_positive = 0, fn = 0;
	CKMCFile kmer_data_base;
	kmer_data_base.OpenForListing(input_file);
	uint32 _kmer_length = kmer_data_base.KmerLength();
	CKmerAPI kmer_object(_kmer_length);
	uint64 kc_kmodel = kmodel->get_kmodel_kmer_count() - kmodel->get_left_kmer_count();
	uint64_t kc_bf1 = kmodel->get_once_kmer_count();
	auto start = chrono::high_resolution_clock::now();
	while (kmer_data_base.ReadNextKmer(kmer_object, counter)) {
		kmer_buff[idx].kmer = kmer_object.to_string().c_str();
		kmer_buff[idx++].occ = counter;
		if (idx >= BLOACK_SIZE) {
			#pragma omp parallel for num_threads(THREAD_NUM) reduction(+:bf1_false_positive)  reduction(+:kmodel_false_positive) 
			for (int i = 0; i < idx; ++i) {
				if (kmodel->get_kmer_from_map(kmer_buff[i].kmer) > -1) { //can be found in the map
					continue;
				}
				if (kmer_buff[i].occ != 1 && kmodel->check_bloomfilter01(kmer_buff[i].kmer))
					++bf1_false_positive;
				if (kmer_buff[i].occ == 1)
					continue;
				if (kmodel->kmer_to_bin(kmer_buff[i].kmer) != kmodel->get_occu_bin()->occ_to_bin(kmer_buff[i].occ)) //don't consider the result from map
					++kmodel_false_positive;
			}
			idx = 0;
		}
	}
	//handle the left kmers
	#pragma omp parallel for num_threads(THREAD_NUM) reduction(+:bf1_false_positive)  reduction(+:kmodel_false_positive) 
	for (int i = 0; i < idx; ++i) {
		if (kmodel->get_kmer_from_map(kmer_buff[i].kmer) > -1) { //can be found in the map
			continue;
		}
		if (kmer_buff[i].occ != 1 && kmodel->check_bloomfilter01(kmer_buff[i].kmer))
			++bf1_false_positive;
		if (kmer_buff[i].occ == 1)
			continue;
		if (kmodel->kmer_to_bin(kmer_buff[i].kmer) != kmodel->get_occu_bin()->occ_to_bin(kmer_buff[i].occ)) //don't consider the result from map
			++kmodel_false_positive;
	}

	kmer_data_base.Close();
	chrono::duration<double> dur = chrono::high_resolution_clock::now() - start;
	printf("test time:%.4f min\n", dur.count() / 60);
	cout << "FP_count_in_kmodel: " << kmodel_false_positive << endl;
	cout << "FP_count_in_bloomfilter01:" << bf1_false_positive << endl;
	printf("kmodel_false_positive:%.3e\nbloomfilter01_false_positive:%.3e\n\n",
		kmodel_false_positive / kc_kmodel, bf1_false_positive / kc_bf1);
}


//input_file is a testing file
void test_present_absent(string input_file, string save_dir) {
	//cout << "=====test nAbsent nPresent=====" << endl;
	cout << input_file << endl;
	//present include both kmer and frequence,while absent does not
	ifstream fin(input_file);
	if (!fin) {
		cout << "can't open the file " << input_file << endl;
		return;
	}
	vector<string> kmer_v;
	string kmer, frquence;
	//load the testing file into vector
	if ((int)input_file.find("present") > -1) 
	{
		while (fin>> kmer >> frquence){
			kmer_v.push_back(kmer);
		}
	}
	else 
	{
		while (fin >> kmer) {
			kmer_v.push_back(kmer);
		}
	}

	auto start = chrono::high_resolution_clock::now();
	KModel* kmodel = Params.ci > 1 ? new KModel() : new KModelOne();
	kmodel->load_model(save_dir);
	chrono::duration<double> dur = chrono::high_resolution_clock::now() - start;
	cout << "openning time: " << dur.count() << endl;
	
	start = chrono::high_resolution_clock::now();
	uint64 nPresent = 0, nAbsent = 0;
	//vector<int> result_v = kmodel->kmer_to_occ(kmer_v);
	for (auto kmer : kmer_v) {
		if (kmodel->kmer_to_occ(kmer))
			++nPresent;
		else
			++nAbsent;
	}
	dur = chrono::high_resolution_clock::now() - start;
	cout << "query time: " << dur.count() << endl;
	cout << "#nAbsent " << nAbsent << endl;
	cout << "#nPresent " << nPresent << endl << endl;
}

void test_open(string input_file,string absent_path) {
	CKMCFile kmer_data_base;
	kmer_data_base.OpenForListing(input_file);
	uint32 _kmer_length = kmer_data_base.KmerLength();
	CKmerAPI kmer_object(_kmer_length);
	float counter;
	while (kmer_data_base.ReadNextKmer(kmer_object, counter)) {
		cout << counter << endl;
	}
}

void test_load(string input_file,string save_dir = "save_model") {
	KModel* kmodel = get_model(Params.ci, Params.num_hash, Params.num_bit);
	cout << "====test_load=====" << endl;
	auto start = chrono::high_resolution_clock::now();
	kmodel->load_model(save_dir);
	chrono::duration<double> dur = chrono::high_resolution_clock::now() - start;
	cout << "opening time: " << dur.count() <<"\n\n";
	if (Params.ci > 1)
		test_fp(kmodel, input_file);
	else
		test_fp((KModelOne*)kmodel, input_file);
}


void test_check_kmer(string kmc_db, string raw) {
	ifstream fin(raw);
	string kmer;
	float occ1,occ2;
	CKMCFile kmer_data_base;
	kmer_data_base.OpenForListing(kmc_db);
	uint32 _kmer_length = kmer_data_base.KmerLength();
	CKmerAPI kmer_object(_kmer_length);
	while (fin >> kmer >> occ1) {
		kmer_object.from_string(kmer);
		kmer_data_base.CheckKmer(kmer_object, occ2);
			printf("%.0f %.0f\n", occ1, occ2);
	}
	fin.close();
	kmer_data_base.Close();
}


void test_kmc3_fasfq_kmer(string kmc_db, string fastq) {
	ifstream fin(fastq);
	string s1,s2,s3,s4,kmer,kmer2;
	int len,i;
	CKMCFile kmer_data_base;
	kmer_data_base.OpenForRA(kmc_db);
	uint32 _kmer_length = kmer_data_base.KmerLength();
	uint64 occ;
	cout << _kmer_length << endl;
	CKmerAPI kmer_object(_kmer_length);
	while (getline(fin, s1)&& getline(fin, s2)&& getline(fin, s3)&& getline(fin, s4)) {
		if (s4.length() <= 0) break;
		len = s2.length();
		for (i = 0; i < len-_kmer_length; i++) {
			kmer = s2.substr(i, _kmer_length);
			//cout <<"OR " << kmer << endl;
			kmer2 = Tools::get_complementation(kmer);
			kmer = kmer.compare(kmer2) < 0 ? kmer : kmer2;
			kmer_object.from_string(kmer);
			//cout <<"LE " << kmer  << endl;
			if (!kmer_data_base.CheckKmer(kmer_object, occ)) {
				cout << kmer << endl << endl;
			}
		}
	}
	fin.close();
	kmer_data_base.Close();
}

void test_kmc3_present_absent(string kmc_db, string input_file) {
	//cout << "===== kmer3 test nAbsent nPresent=====" << endl;
	cout << input_file << endl;
	//present include both kmer and frequence,while absent does not
	ifstream fin(input_file);
	if (!fin) {
		cout << "can't open the file " << input_file << endl;
		return;
	}
	vector<string> kmer_v;
	string kmer, frquence;
	//load the testing file into vector
	if ((int)input_file.find("present") > -1)
	{
		while (fin >> kmer >> frquence) {
			kmer_v.push_back(kmer);
		}
	}
	else
	{
		while (fin >> kmer) {
			kmer_v.push_back(kmer);
		}
	}

	auto start = chrono::high_resolution_clock::now();
	CKMCFile kmer_data_base;
	kmer_data_base.OpenForRA(kmc_db);
	uint32 _kmer_length = kmer_data_base.KmerLength();
	CKmerAPI kmer_object(_kmer_length);
	chrono::duration<double> dur = chrono::high_resolution_clock::now() - start;
	cout << "openning time: " << dur.count() << endl;

	start = chrono::high_resolution_clock::now();
	uint64 nPresent = 0, nAbsent = 0;
	uint64 count;
	int len = kmer_v.size();
	for (int i = 0; i < len; i++) {
		kmer_object.from_string(kmer_v[i]);
		if (kmer_data_base.CheckKmer(kmer_object, count))
			++nPresent;
		else
			++nAbsent;
	}
	kmer_data_base.Close();
	dur = chrono::high_resolution_clock::now() - start;
	cout << "query time: " << dur.count() << endl;
	cout << "#nAbsent " << nAbsent << endl;
	cout << "#nPresent " << nPresent << endl << endl;
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

	max_counter = Params.cs + 1;
	THREAD_NUM = Params.t;
	cout << Params.output_file_name << endl;
	KModel* kmodel = get_model(Params.ci, Params.num_hash, Params.num_bit);
	kmodel->init_KModel(Params.output_file_name);
	kmodel->show_model_info();
	//save model to the disk
	string save_dir = Params.working_directory + "/" + Tools::get_file_name(Params.output_file_name);
	system(("mkdir -p " + save_dir).c_str()); //new a dirctory - -;
	kmodel->save_model(save_dir);
	cout << "save model successfully.!" << endl;
	Params.ci > 1 ? test_fp(kmodel, Params.output_file_name) : test_fp((KModelOne*)kmodel, Params.output_file_name);
}


int _tmain(int argc, char* argv[]) {
	run(argc, argv);
	return 0;
}