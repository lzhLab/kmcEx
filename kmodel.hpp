#pragma once

#ifndef KMODEL_H
#define KMODEL_H

#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <cstdlib>
#include <stdint.h>
#include <stdio.h>
#include <cstring>
#include "kmc_api/kmc_file.h"
#include "occu_bin.hpp"
#include "tools.hpp"
using namespace std;

const int Inf = 0x7fffffff;

//save the model data
struct BitSaveData {
	uint64_t bit_array_length; //bit_array_length when used in hash phase
	uint32_t* hash_seed; //the array to save hash seed
	uint8_t *bit_array_1; //a bit array used to save kmer's occurrence
	uint8_t *bit_array_2; //a tag bit array
};


class KModel
{
public:
	KModel() {}
	KModel(OccuBin* occu_bin, int n_bits) {
		this->occu_bin = occu_bin;
		n_hash = occu_bin->get_n_hash();
		kmer_count = occu_bin->get_kmer_count();
		this->n_bits = n_bits;
		set_bit_array_size(kmer_count);
	}

	// n_bit:how many bit_array used
	void init_KModel(string file_kmer_data_base) {
		int occ32, i, j;
		bits_data = new BitSaveData[n_bits];
		string kmer;
		float counter;
		uint8_t occ8;
		uint64_t key, t_byte_array_size = byte_array_size;
		left_kmer_count = 0;
		for (i = 0; i < n_bits; i++) {
			BitSaveData bit_data;
			bit_data.bit_array_length = t_byte_array_size << 3;
			bit_data.bit_array_1 = new uint8_t[t_byte_array_size];
			bit_data.bit_array_2 = new uint8_t[t_byte_array_size];
			memset(bit_data.bit_array_1, 0, sizeof(uint8_t)*t_byte_array_size);
			memset(bit_data.bit_array_2, 0, sizeof(uint8_t)*t_byte_array_size);
			bit_data.hash_seed = new uint32_t[n_hash];
			//get num_bit hash seed
			for (j = 0; j < n_hash; j++) {
				int index = (i*n_hash + j) % 32;
				bit_data.hash_seed[j] = HashSeeds[index];
			}
			bits_data[i] = bit_data;
			t_byte_array_size = (uint64_t)(t_byte_array_size*0.47); //byte_array_size reduce 1 half
		}

		CKMCFile kmer_data_base;
		if (!kmer_data_base.OpenForListing(file_kmer_data_base)) {
			std::cout << "init_occ_count()__can't open the kmer_data_base\n";
			return;
		}
		uint32 _kmer_length = kmer_data_base.KmerLength();
		CKmerAPI kmer_object(_kmer_length);
		//uint64_t tc = 0;
		while (kmer_data_base.ReadNextKmer(kmer_object, counter)) {
			if (counter >= C_MAX)  //skip the cs in the kmc2
				continue;
			occ32 = counter;
			kmer = kmer_object.to_string().c_str();
			//occ_to_bin,save compressed
			occ8 = (uint8_t)(occu_bin->occ_to_bin(occ32));
			for (i = 0; i < n_bits; i++) {
				if (insert_to_bit_array(kmer, occ8, i)) { //insert success
					c4[i]++;
					break;
				}
			}
			if (i >= n_bits) {//if insert bit array fail,then save it to map
				key = Tools::kmers2uint64(kmer);
				last_map[key] = occ32;
				left_kmer_count++;
			}
		}
		kmer_data_base.Close();
	}

	int kmer_to_occ(string kmer) {
		int occ = get_kmer_from_map(kmer);
		if (occ != -1)
			return occ;
		int count = 0;
		int bin = kmer_to_bin(kmer, count);
		occ = occu_bin->bin_to_mean(bin);
		return occ;
	}

	//alter the char in the pos as the neighbor of kmer
	void shift_one(string kmer, vector<int>  &v_candidates) {
		vector<int> v_bin = find_bitarray(kmer);
		if (v_bin.size() >= 1) {
			v_candidates.push_back(v_bin[0]);
		}
		else if (v_bin.size() <= 0) { // try kmer's complementation
			string nw_kmer_com = Tools::get_complementation(kmer);
			v_bin = find_bitarray(nw_kmer_com);
			if (v_bin.size() >= 1) {
				v_candidates.push_back(v_bin[0]);
			}
		}
	}

	vector<int> get_neighbor_kmer_bin(string kmer) {
		vector<int> v_candidates;
		string s = "ACGT", t1_kmer, t2_kmer;
		int kmer_len = kmer.length();
		t1_kmer = kmer.substr(1); //remove the first char
		t2_kmer = kmer.substr(0, kmer_len - 1);//remove the last char ? -2 kmc2 bug? 
		for (int i = 0; s[i]; i++) { //shift forward
			string nw_kmer = t1_kmer + s[i];
			shift_one(nw_kmer, v_candidates);

		}
		for (int i = 0; s[i]; i++) { //shift back
			string nw_kmer = s[i] + t2_kmer;
			shift_one(nw_kmer, v_candidates);
		}
		return v_candidates;
	}

	//appearance: how many times the k appearances in the different bit array
	int kmer_to_bin(string kmer, int &len_v_bin) {
		vector<int> v_bin = find_bitarray(kmer);
		len_v_bin = v_bin.size();
		if (len_v_bin <= 1) {
			//if (gt_bin != v_bin[0]) printf("%d %d\n", gt_bin, v_bin[0]);
			return v_bin[0]; //only one
		}
		//k can be found in more than one coupling-bit arrays.
		vector<int> v_candidates = get_neighbor_kmer_bin(kmer);
		int len_v_can = v_candidates.size();
		if (len_v_can <= 0)//no any candidates,get the minimum bin
			return Tools::vector_min(v_bin);
		//get the best bin
		int min_dist = 9999, best_bin;
		for (int i = 0; i < len_v_bin; i++) {
			int cur_dist = 0;
			for (int j = 0; j < len_v_can; j++) {
				cur_dist += abs(v_bin[i] - v_candidates[j]);
			}
			if (min_dist > cur_dist) {
				min_dist = cur_dist;
				best_bin = v_bin[i];
			}
		}
		return best_bin;
		//printf("appearance:%d, candidates:%d, gt_bin:%d\n", appearance, len_cand, gt_bin);
		//for (int i = 0; i <appearance; i++)
		//	cout << v_bin[i] << " ";
		//cout << endl;
		//for (int i = 0; i < len_cand; i++)
		//	cout << v_candidates[i] << " ";
		//cout << endl << endl;
	}

	//find in the map
	int get_kmer_from_map(string kmer) {
		uint64_t key = Tools::kmers2uint64(kmer);
		unordered_map<uint64_t, int>::iterator it = last_map.find(key);
		if (it != last_map.end()) {
			return it->second;
		}
		return -1;
	}

	//save model
	void save_model(string save_dir) {
		//1. save model parameters:n_hash,n_bits,kmer_count
		ofstream fout(save_dir + "/param.conf");
		fout << "number_hash " << n_hash << endl;
		fout << "number_bit " << n_bits << endl;
		fout << "kmer_count  " << kmer_count << endl;
		fout << "last_map_size  " << last_map.size() << endl;
		fout << "C_MAX   " << C_MAX << endl;
		fout.close();

		//2. save OccuBin
		FILE *fp_occ_meta = fopen((save_dir + "/occ.bin").c_str(), "wb");
		fwrite(occu_bin->get_occ_bin_meta(), sizeof(OccuBinMeta), C_MAX, fp_occ_meta);
		fclose(fp_occ_meta);

		//3. save KModel
		for (int i = 0; i < n_bits; i++) {
			//3.1 hash_function
			FILE *fp_hash = fopen((save_dir + "/hash_" + std::to_string(i) + ".bin").c_str(), "wb");
			fwrite(bits_data[i].hash_seed, sizeof(uint32_t), n_hash, fp_hash);
			fclose(fp_hash);
			//3.2 bit_array_1
			FILE *fp_bit1 = fopen((save_dir + "/bit1_" + std::to_string(i) + ".bin").c_str(), "wb");
			fwrite(bits_data[i].bit_array_1, sizeof(uint8_t), bits_data[i].bit_array_length >> 3, fp_bit1);
			fclose(fp_bit1);
			//3.3 bit_array_2
			FILE *fp_bit2 = fopen((save_dir + "/bit2_" + std::to_string(i) + ".bin").c_str(), "wb");
			fwrite(bits_data[i].bit_array_2, sizeof(uint8_t), bits_data[i].bit_array_length >> 3, fp_bit2);
			fclose(fp_bit2);
		}

		//4. save last_map
		FILE *fp_map_w = fopen((save_dir + "/last_map.bin").c_str(), "wb");
		unordered_map<uint64_t, int>::iterator iter;
		for (iter = last_map.begin(); iter != last_map.end(); iter++) {
			fwrite(&iter->first, sizeof(uint64_t), 1, fp_map_w);
			fwrite(&iter->second, sizeof(int), 1, fp_map_w);
		}
		fclose(fp_map_w);

	}

	//load model
	void load_model(string save_dir) {
		//1. load model parameters:n_hash,n_bits,kmer_count
		ifstream fin(save_dir + "/param.conf");
		if (!fin) {
			cout << "load_model__cant't open the param.conf !\n";
			return;
		}
		string t_str;
		int last_map_size;
		fin >> t_str >> n_hash;
		fin >> t_str >> n_bits;
		fin >> t_str >> kmer_count;
		fin >> t_str >> last_map_size;
		fin >> t_str >> C_MAX;
		set_bit_array_size(kmer_count);
		fin.close();

		//2. load OccuBin
		occu_bin = new OccuBin(n_hash);
		FILE *fp_occ_meta = fopen((save_dir + "/occ.bin").c_str(), "rb");
		fread(occu_bin->get_occ_bin_meta(), sizeof(OccuBinMeta), C_MAX, fp_occ_meta);
		fclose(fp_occ_meta);
		//3. load KModel
		uint64_t t_byte_array_size = byte_array_size;
		bits_data = new BitSaveData[n_bits];
		for (int i = 0; i < n_bits; i++) {
			//3.1 hash_function
			bits_data[i].hash_seed = new uint32_t[n_hash];
			bits_data[i].bit_array_length = t_byte_array_size << 3;
			FILE *fp_hash = fopen((save_dir + "/hash_" + std::to_string(i) + ".bin").c_str(), "rb");
			fread(bits_data[i].hash_seed, sizeof(uint32_t), n_hash, fp_hash);
			fclose(fp_hash);
			//3.2 bit_array_1
			bits_data[i].bit_array_1 = new uint8_t[t_byte_array_size];
			FILE *fp_bit1 = fopen((save_dir + "/bit1_" + std::to_string(i) + ".bin").c_str(), "rb");
			fread(bits_data[i].bit_array_1, sizeof(uint8_t), t_byte_array_size, fp_bit1);
			fclose(fp_bit1);
			//3.3 bit_array_2
			bits_data[i].bit_array_2 = new uint8_t[t_byte_array_size];
			FILE *fp_bit2 = fopen((save_dir + "/bit2_" + std::to_string(i) + ".bin").c_str(), "rb");
			fread(bits_data[i].bit_array_2, sizeof(uint8_t), t_byte_array_size, fp_bit2);
			fclose(fp_bit2);
			t_byte_array_size = t_byte_array_size*0.47;
		}
		//printf("4. load last_map %d \n",last_map_size);
		//4. load last_map
		uint64_t k;
		int v;
		FILE *fp_map_r = fopen((save_dir + "/last_map.bin").c_str(), "rb");
		for (int i = 0; i < last_map_size; i++) {
			fread(&k, sizeof(uint64_t), 1, fp_map_r);
			fread(&v, sizeof(int), 1, fp_map_r);
			last_map.insert(make_pair(k, v));
		}
		fclose(fp_map_r);
	}

	void show_header() {
		cout << "kmer_count:" << kmer_count << endl;
		cout << "number_bit_array:" << n_bits << endl;
		cout << "byte_array_size:" << byte_array_size << endl;
		cout << "number_hash_function:" << n_hash << endl;
		printf("dist_error:%.6f\n", occu_bin->calc_distance_error());
	}

	void show_header2() {
		uint64_t total_byte_array_size = 0;
		for (int i = 0; i < 4; i++) {
			total_byte_array_size += bits_data[i].bit_array_length >> 3;
		}

		double left = kmer_count;
		for (int i = 0; i < 4; i++) {
			printf("%02d--save_%d--left_%.0f---saveratio_%.3f\n", i + 1, c4[i], left - c4[i], c4[i] / left);
			left -= c4[i];
		}
		cout << "memory used in bit array:" << Tools::filesize_format(total_byte_array_size * 2) << endl;
		cout << "the left kmer is saved to map:" << left_kmer_count << endl;

		//get the dict the number of uint_8 -> sum(0,1)
		int b = 8, n = 1 << b;
		int* dic_num = new int[n];
		for (int i = 0; i<n; i++) {
			uint8_t sum = 0, v = i;
			for (int j = 0; j<b; j++) {
				sum += v & 0x1;
				v >>= 1;
			}
			dic_num[i] = sum;
		}

		double total_c = 0;
		for (int i = 0; i < n_bits; i++) {
			double c2 = 0, byte_size = bits_data[i].bit_array_length >> 3;
			for (uint64_t p = 0; p < byte_size; p++) {
				int idx = bits_data[i].bit_array_2[p];
				c2 += dic_num[idx];
			}
			total_c += c2;
			printf("usage rate of the no.%d bit array:%0.3f\n", i + 1, c2 / bits_data[i].bit_array_length);
		}
		//printf("usage rate of the all bit array:%0.3f\n", total_c / n_bits / bit_array_length);
		delete[] dic_num;
	}

	uint64_t get_kmer_count() {
		return kmer_count;
	}

	int get_left_kmer_count() {
		return left_kmer_count;
	}

	OccuBin* get_occu_bin() {
		return occu_bin;
	}


private:
	OccuBin* occu_bin;
	int n_hash;
	int n_bits;
	int left_kmer_count; //the kmer count that can't save to bit_array and save to last_map
	uint64_t kmer_count;
	uint64_t byte_array_size; //byte size, default is the kmer_count*0.4
	uint64_t bit_array_length; //the real bit table length,bit_array_length=byte_array_size*8
	unordered_map<uint64_t, int> last_map; //save the left kmer
	BitSaveData* bits_data;
	int c4[4] = { 0 }; //how many kmer save in the bit_array

	void set_bit_array_size(uint64_t _kmer_count) {
		//to minimize the false positive -> bit_array_length=kmer_count*n_hash/ln2
		byte_array_size = _kmer_count / 8.0*n_hash / 0.69;
		bit_array_length = byte_array_size << 3; //byte_array_size*8
	}


	//set the position in the bit to one
	void set_bit(uint8_t *bit, uint64_t pos) {
		uint64_t row = pos / 8;
		uint64_t col = pos % 8;
		uint8_t x = 0x1 << (8 - col - 1); //set the bit in the column in the table 
		bit[row] |= x;
	}

	//check the position in the bit is one or zero 
	bool check_bit(const uint8_t *bit, uint64_t pos) {
		uint64_t row = pos / 8;
		uint64_t col = pos % 8;
		return (bit[row] >> (8 - col - 1)) & 0x1;
	}

	//----------------------------------------------------------------------------------
	//if the bit_array_2'value is one and the bit in the occurrence is different from the bit_array_1' value,
	//then we can insert this kmer into bit array,or write it to the disk
	//--------------------------------------------------------------------------------
	bool insert_to_bit_array(std::string kmer, uint8_t occ8, int index) {

		uint64_t* bit1_v = new uint64_t[n_hash];// the binary value of frequency in bit_array_1 
		uint64_t* bit2_v = new uint64_t[n_hash];// the index of bit_array_2
		uint32_t*  hash_seed = bits_data[index].hash_seed;
		uint8_t* bit_array_1 = bits_data[index].bit_array_1;
		uint8_t* bit_array_2 = bits_data[index].bit_array_2;
		uint64_t array_length = bits_data[index].bit_array_length;
		for (int i = 0; i<n_hash; i++) {
			bit1_v[i] = occ8 & 0x1;
			bit2_v[i] = Tools::murmur_hash64(kmer.c_str(), kmer.size(), hash_seed[i]) % array_length;
			occ8 >>= 1;
		}
		bool ok = true; //the flag if can insert this kmer
		for (int i = 0; i<n_hash; i++) {
			uint8_t v1 = check_bit(bit_array_1, bit2_v[i]);
			uint8_t v2 = check_bit(bit_array_2, bit2_v[i]);
			if (v2 && v1 != bit1_v[i]) {
				ok = false;
			}
		}
		if (ok) {
			for (int i = 0; i<n_hash; i++) {
				if (bit1_v[i]) {
					set_bit(bit_array_1, bit2_v[i]);
				}
				set_bit(bit_array_2, bit2_v[i]);
			}
		}
		delete[] bit1_v;
		delete[] bit2_v;
		return ok;
	}

	//travel the bit array to get the occurrence of the kmer
	vector<int> find_bitarray(string kmer) {
		vector<int> v_bin;
		int result = -1;
		uint64_t pos;
		uint8_t* occ_bit = new uint8_t[n_hash];
		for (int i = 0; i< n_bits; i++) {
			bool ok = true;
			for (int j = 0; j< n_hash; j++) {
				pos = Tools::murmur_hash64(kmer.c_str(), kmer.size(), bits_data[i].hash_seed[j]) % bits_data[i].bit_array_length;
				occ_bit[j] = check_bit(bits_data[i].bit_array_1, pos);
				if (!check_bit(bits_data[i].bit_array_2, pos)) {
					ok = false;
				}
			}
			if (ok) {
				result = Tools::bin_to_decimal(occ_bit, n_hash);
				v_bin.push_back(result);
			}
		}
		//find in the map
		//uint64_t key = Tools::kmers2uint64(kmer);
		//unordered_map<uint64_t, int>::iterator it = last_map.find(key);
		//if (it != last_map.end()) {
		//	v_bin.push_back(it->second);
		//}
		delete[] occ_bit;
		return  v_bin;
	}

};

#endif