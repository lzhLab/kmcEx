#pragma once
#include <stdint.h>
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include "kmc_api/kmc_file.h"
using namespace std;

int C_MAX = 1001; //the max occurrence,exclude self

//OccuBinMeta
//the bean struct needed by OccBin
//the instance of OccuBinMeta show be a array,which it's index represent occurrence
//e.g OccuBinMeta ocm[10],  ocm[i].occ_bin -> i=occurrence; mean=ocm[i].occ_bin
struct OccuBinMeta {
	int occ_mean; //occ->mean, e.g.100,101,102,103 ->102
	int occ_bin; //occ->bin, e.g.102->64
};


//OccuBin
//OccMap used to translate a occurrence to optimal number index,for example:
//we assume that array A that length is 64 to save occurrence, 1-31 save the orignal and unchanged occurrence,
//and 32-63 save the occurrence translated that we used some compressed method 
class OccuBin
{
public:

	OccuBin(int n_hash = 7) {
		max_bin = (1 << n_hash) - 1;
		this->n_hash = n_hash;
		start_bin = max_bin / 2 + 2; //(1 << n_hash)/2
		occ_bin_meta = new OccuBinMeta[C_MAX];
		occ_count = new int[C_MAX];
		memset(occ_count, 0, sizeof(int)*C_MAX);
	}

	//initialize the occ_count from kmer_count file 
	void init_occ_count(string file_kmer_data_base) {
		float counter;
		int occ;

		CKMCFile kmer_data_base;
		if (!kmer_data_base.OpenForListing(file_kmer_data_base)) {
			cout << "init_occ_count()__can't open the kmer_data_base\n";
			return ;
		}
		kmer_count = kmer_data_base.KmerCount();
		uint32 _kmer_length = kmer_data_base.KmerLength();
		CKmerAPI kmer_object(_kmer_length);
		while (kmer_data_base.ReadNextKmer(kmer_object, counter)) {
			occ = counter;
			if (occ < C_MAX) occ_count[occ]++;
		}
		kmer_data_base.Close();
		init_occ_bin(); 
	}

	int occ_to_bin(int occ) {
		if (occ<start_bin) 
			return occ;
		return occ_bin_meta[occ].occ_bin;
	}

	int bin_to_mean(int occ_bin) {
		if (occ_bin<start_bin) return occ_bin;
		for (int i = start_bin; i<C_MAX; i++) {
			if (occ_bin_meta[i].occ_bin == occ_bin) {
				return occ_bin_meta[i].occ_mean;
			}
		}
		return 0;
	}

	double calc_distance_error() {
		double sum_dis = 0, sum_kmer_count = 0;
		for (int i = 2; i<C_MAX; i++) {
			sum_kmer_count += i*occ_count[i];
			//cout << i << " " << occ_count[i] << " " << occ_bin_meta[i].occ_mean << " " << occ_bin_meta[i].occ_bin << endl;
		}
		for (int i = start_bin; i < C_MAX; i++) {
			sum_dis += abs(occ_bin_meta[i].occ_mean - i)*occ_count[i];
		}
		double dis_error = sum_dis / sum_kmer_count;
		return dis_error;
	}


	int get_n_hash() {
		return n_hash;
	}

	uint64_t get_kmer_count() {
		return kmer_count;
	}

	OccuBinMeta* get_occ_bin_meta() {
		return occ_bin_meta;
	}

private:
	int n_hash;
	int start_bin; //2^(n_hash)/2
	int max_bin; //2^(n_hash)-1
	uint64_t kmer_count;
	OccuBinMeta* occ_bin_meta;
	int* occ_count;

	void init_occ_bin() {
		int sum = 0, i;
		//equal depth: start_bin~(start_bin + start_bin/3)
		int width_start_bin = start_bin + start_bin * 1 / 3;
		int width_start = start_bin + C_MAX * 1 / 3;
		for (i = start_bin; i < width_start; i++) { //start_bin~width_start is the equal depth
			sum += occ_count[i];
		}
		int base = sum / (width_start_bin - start_bin + 1);
		int bin = start_bin, start = start_bin, mean;
		for (i = start_bin, sum = 0; i < width_start; i++) {
			sum += occ_count[i];
			if (!is_more_closer(sum, occ_count[i + 1], base) || i >= width_start - 1) {
				mean = (start + i) / 2;
				set_section_mean(mean, bin, start, i);
				if (bin < width_start_bin) bin++;
				sum = 0;
				start = i + 1; //start=end+1
			}
		}
		//equal width
		width_start_bin = bin; //the bin maybe have not arrive at the width_start_bin
		int bin_capacity = ceil(1.0*(C_MAX - width_start) / (max_bin - width_start_bin + 1));
		for (i = width_start; i < C_MAX - 1; i+= bin_capacity) {
			mean = (2 * i + bin_capacity - 1) / 2;
			int end = i + bin_capacity - 1; //exclude end
			if (end >= C_MAX - 1) { //beyond the C_MAX-1
				end = C_MAX - 1;
				mean = (i+ C_MAX - 1) / 2;
			}
			set_section_mean(mean, width_start_bin, i, end);
			width_start_bin++;
		}
	}

	bool is_more_closer(int sum, int increment, int base) {
		return abs(sum - base) > abs(sum + increment - base);
	}

	void set_section_mean(int mean, int bin, int start, int end) {
		for (int i = start; i <= end; i++) {
			occ_bin_meta[i].occ_mean = mean;
			occ_bin_meta[i].occ_bin = bin;
		}
	}
};
