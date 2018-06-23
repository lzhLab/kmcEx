#include "kmodel.hpp"

using namespace std;


int _tmain(int argc, char* argv[]) {
	int n_hash = 7;
	int n_bit = 4;
	string input_file(argv[1]);
	n_hash = atoi(argv[2]);
	cout << "\n**************" << input_file << "**************\n";
	//string model_dir(argv[2]);
	//string flag(argv[3]);
	KModel* kmodel;
	OccuBin* occu_bin = new OccuBin(n_hash);
	occu_bin->init_occ_count(input_file);
	
	return 0;

	kmodel = new KModel(occu_bin, n_bit);
	kmodel->show_header();
	kmodel->init_KModel(input_file);
	kmodel->show_header2();

	//cout << "testing......." << endl;
	//testing.......
	string kmer = "";
	float counter;
	int occ, t_bin, bin, appearance = 1;
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
		t_bin = kmodel->kmer_to_bin(kmer, appearance,bin);
		if (appearance > 1) {
			false_positive++;
		}
		if (t_bin != bin) //don't consider the result from map
			retrieval_error++;
	}
	kmer_data_base.Close();

	uint64_t bit_kmer_count = kmodel->get_kmer_count() - left_map-last_f;
	printf("false_positive:%0.6f, retrieval_error:%0.6f\n", false_positive / bit_kmer_count, retrieval_error / bit_kmer_count);

}