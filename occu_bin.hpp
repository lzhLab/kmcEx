#pragma once
#include <stdint.h>
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <unordered_map>
#include "kmc_api/kmc_file.h"
using namespace std;


//OccuBinMeta
//the bean struct needed by OccBin
//the instance of OccuBinMeta show be a array,which it's index represent occurrence
//e.g OccuBinMeta ocm[10],  ocm[i].occ_bin -> i=occurrence; mean=ocm[i].occ_bin
class OccuBinMeta {
public:
	uint32_t occ_mean = -1; //occ->mean, e.g.100,101,102 ->101
	uint32_t occ_bin = -1; //occ->bin, e.g.101->65
};


///OccuBin
//[1/4, 1/2, 1/4]
class OccuBin {
public:
	OccuBin(int max_counter,int n_hash = 7) {
		this->max_counter = max_counter;
		this->n_hash = n_hash;
		this->bin_end_index3 = 1 << n_hash;
		this->bin_end_index1 = bin_end_index3 / 4;
		this->bin_end_index2 = bin_end_index1 + bin_end_index3 / 2;
		this->occ_bin_meta = new OccuBinMeta[max_counter];
		//the numbers of bin2 is a half of 1 << n_hash ==1/2
		int bin2_num = bin_end_index3 / 2;
		int bin2_capacity = 3;
		int bin2_start = bin_end_index1;
		for (int i = 0; i < bin2_num; i++) {
			for (int j = 0; j < bin2_capacity; j++) {
				occ_bin_meta[bin2_start + j].occ_mean = bin2_start + 1; //mean value
				occ_bin_meta[bin2_start + j].occ_bin = bin_end_index1 + i;
			}
			bin2_start += bin2_capacity;
		}
		//the numbers of bin2 is 1/4 of 1 << n_hash ==1/4
		int bin3_num = bin_end_index3 / 4;
		int bin3_capacity = (max_counter - bin2_start) / bin3_num;
		for (int i = 0; i <bin3_num; i++) {
			for (int j = 0; j < bin3_capacity; j++) {
				occ_bin_meta[bin2_start + j].occ_mean = (2 * bin2_start + bin3_capacity) / 2; //mean value
				occ_bin_meta[bin2_start + j].occ_bin = bin_end_index2 + i;
			}
			bin2_start += bin3_capacity;
		}
		//handle the left,using the previous reslut
		for (int i = bin2_start; i < max_counter; i++) {
			occ_bin_meta[i].occ_mean = (2 * bin2_start - bin3_capacity) / 2;
			occ_bin_meta[i].occ_bin = bin_end_index3 - 1;
		}
		//bin2mean
		for (int i = bin_end_index1; i<max_counter; i++) {
			bin2mean.insert(make_pair(occ_bin_meta[i].occ_bin, occ_bin_meta[i].occ_mean));
		}
	}


	int occ_to_bin(int occ) {
		if (occ<bin_end_index1)
			return occ;
		return occ_bin_meta[occ].occ_bin;
	}

	uint32_t occ_to_bin(uint32_t occ) {
		if (occ<bin_end_index1)
			return occ;
		return occ_bin_meta[occ].occ_bin;
	}

	uint32_t bin_to_mean(uint32_t occ_bin) {
		if (occ_bin<bin_end_index1) 
			return occ_bin;
		return bin2mean[occ_bin];
	}


	int get_bin_end_index1() {
		return this->bin_end_index1;
	}

	int get_max_counter() {
		return this->max_counter;
	}

	int get_hash_number() {
		return this->n_hash;
	}

	OccuBinMeta* get_occ_bin_meta() {
		return occ_bin_meta;
	}


private:
	int max_counter; //	max_counter = Params.cs + 1;
	int n_hash;
	int bin_end_index1; //2^(n_hash-2)
	int bin_end_index2; //2^(n_hash-2)+2^(n_hash-1)
	int bin_end_index3; //2^(n_hash)
	unordered_map<uint32_t, uint32_t> bin2mean;
	OccuBinMeta* occ_bin_meta;
};
