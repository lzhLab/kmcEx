#pragma once
#include <stdint.h>
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include "kmc_api/kmc_file.h"
using namespace std;

int C_MAX = 1001; //the max occurrence

//OccuBinMeta
//the bean struct needed by OccBin
//the instance of OccuBinMeta show be a array,which it's index represent occurrence
//e.g OccuBinMeta ocm[10],  ocm[i].occ_bin -> i=occurrence; mean=ocm[i].occ_bin
class OccuBinMeta {
public:
	int occ_mean = -1; //occ->mean, e.g.100,101,102 ->101
	int occ_bin = -1; //occ->bin, e.g.101->65
};


///OccuBin
///OccMap used to translate a occurrence to optimal number index,for example:
///we assume that array A that length is 64 to save occurrence, 1-31 save the orignal and unchanged occurrence,
///and 32-63 save the occurrence translated that we used some compressed method 
class OccuBin {
public:
	OccuBin(int n_hash = 7) {
		this->n_hash = n_hash;
		this->bin_end_index3 = 1 << n_hash;
		this->bin_end_index1 = bin_end_index3 / 4;
		this->bin_end_index2 = bin_end_index1 + bin_end_index3 / 2;
		this->occ_bin_meta = new OccuBinMeta[C_MAX];
		//init occ_bin_meta
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
		//cout << bin2_start << endl;
		//printf("%d %d %d\n", bin_end_index1, bin_end_index2, bin_end_index3);
		int bin3_num = bin_end_index3 / 4;
		int bin3_capacity = (C_MAX - bin2_start) / bin3_num;
		for (int i = 0; i <bin3_num; i++) {
			for (int j = 0; j < bin3_capacity; j++) {
				occ_bin_meta[bin2_start + j].occ_mean = (2 * bin2_start + bin3_capacity) / 2; //mean value
				occ_bin_meta[bin2_start + j].occ_bin = bin_end_index2 + i;
			}
			bin2_start += bin3_capacity;
		}
		//handle the left,using the previous reslut
		for (int i = bin2_start; i < C_MAX; i++) {
			occ_bin_meta[i].occ_mean = (2 * bin2_start - bin3_capacity) / 2;
			occ_bin_meta[i].occ_bin = bin_end_index3 - 1;
		}
		//for (int i = 0; i < 1001; i++) {
		//	printf("%d %d %d\n", i, occ_bin_meta[i].occ_mean, occ_bin_meta[i].occ_bin);
		//}
	}

	int occ_to_bin(int occ) {
		if (occ<bin_end_index1)
			return occ;
		return occ_bin_meta[occ].occ_bin;
	}

	int bin_to_mean(int occ_bin) {
		if (occ_bin<bin_end_index1) 
			return occ_bin;
		for (int i = bin_end_index1; i<C_MAX; i++) {
			if (occ_bin_meta[i].occ_bin == occ_bin) {
				return occ_bin_meta[i].occ_mean;
			}
		}
		return 0;  /////////////////
	}


	int get_bin_end_index1() {
		return this->bin_end_index1;
	}

	int get_nhash() {
		return this->n_hash;
	}

	OccuBinMeta* get_occ_bin_meta() {
		return occ_bin_meta;
	}


private:
	int n_hash;
	int bin_end_index1; //2^(n_hash-2)
	int bin_end_index2; //2^(n_hash-2)+2^(n_hash-1)
	int bin_end_index3; //2^(n_hash)
	OccuBinMeta* occ_bin_meta;
};
