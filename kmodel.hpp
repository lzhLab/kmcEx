#pragma once

#ifndef KMODEL_H
#define KMODEL_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <stdint.h>
#include <stdio.h>
#include <cstring>
#include "kmc_api/kmc_file.h"
#include "occu_bin.hpp"
#include "tools.hpp"
#include "rest.hpp"
#include <chrono>
#include "omp.h"
#include <atomic> 
using namespace std;

const string BASE_CHAR = "ACGT";
const int BLOACK_SIZE = 1 << 19;

struct KmerBuff {
	string kmer;
	uint32_t occ;
};

//model data struct
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

	KModel(OccuBin* occu_bin, int n_bits, int ci) {
		this->occu_bin = occu_bin;
		this->n_bits = n_bits;
		this->ci = ci;
		cs = occu_bin->get_max_counter() - 1;
		bf_num = ci == 1 ? 1 : 3;
		n_hash = occu_bin->get_hash_number();
		km_back_num_hash = n_hash - 2;
		bf_num_hash = n_hash - 1;
		bf_back_num_hash = n_hash - 2;
	}

	void init(string db_file) {
		//cout << "ci " << ci << endl;
		string kmer;
		uint32 counter;
		CKMCFile kmer_data_base;
		CKmerAPI kmer_object;
		kmc_api(db_file, kmer_data_base, kmer_object);
		km_kmercount = get_km_kmer_count(kmer_data_base, kmer_object);
		init_km_bit(km_kmercount);
		show_header_info(); //message show
		auto start_build = chrono::high_resolution_clock::now();
		while (kmer_data_base.ReadNextKmer(kmer_object, counter)) {
			kmer = kmer_object.to_string();
			if (counter < ci + bf_num)
				push_to_bloomfilter(kmer, counter);
			else
				push_to_array(kmer, counter);
		}
		//handle the last
		push_last_to_array(n_bits);
		push_last_to_bloomfilter();

		//KRestData build
		kld->build();
		delete[] km_buff;
		delete[] bf_buff;
		delete km_buff_num;
		chrono::duration<double> dur = chrono::high_resolution_clock::now() - start_build;
		build_time_cost = dur.count();
	}



	vector<int> kmer_to_occ(vector<string> kmer_v, int t_num = 4) {
		int n = kmer_v.size();
		vector<int> occ_v(n);
#pragma omp parallel for num_threads(t_num) 
		for (int i = 0; i < n; i++) {
			occ_v[i] = kmer_to_occ(kmer_v[i]);
		}
		return occ_v;
	}

	int kmer_to_occ(string kmer, uint32_t r_occ = 0) {
		//1.get the mini kmer
		kmer = Tools::get_min_kmer(kmer);
		//2.check map
		int occ = kld->check_kmer(kmer);
		if (occ != 0) return occ;
		//3.check km_back
		bool is_in_back = check_back_bloomfilter(kmer, km_back, bit_km_back, km_back_num_hash);
		occ = check_all_bf(kmer);
		if (occ != 0 && !is_in_back) return occ;
		//if (occ != 0) return occ;
		if (!is_in_back) return 0;
		//4.check bit
		int bin = kmer_to_bin(kmer, occ);
		occ = occu_bin->bin_to_mean(bin);
		return occ;
	}

	void show_header_info() {
		cout << "KMCEX:" << endl;
		cout << "   kmodel number hash                 :     " << n_hash << endl;
		cout << "   kmodel bit array                   :     " << n_bits << endl;
		cout << "   total kmercount                    :     " << total_kmer_count << endl;
		cout << "   kmercount in blommfilter           :     " << bf_kmercount << endl;
		cout << "   kmercount in kmodel                :     " << km_kmercount << endl;
	}

	void show_kmodel_info() {
		uint64_t total_byte_size = 0;
		uint64_t bf_byte_size = 0;
		for (int i = 0; i < bf_num; i++) {
			bf_byte_size += byte_bf[i] + byte_bf_back[i];
		}
		uint64_t km_byte_size = 0;
		for (int i = 0; i < n_bits; i++) {
			km_byte_size += bits_data[i].bit_array_length >> 3;
		}
		km_byte_size *= 2;
		uint64_t map_byte = kld->get_all_byte_size();
		total_byte_size = bf_byte_size + km_byte_size + map_byte + byte_km_back;
		cout << "   kmercount hash map                 :     " << kld->get_rest_count() << endl;
		cout << "   memory bloomfilter                 :     " << Tools::filesize_format(bf_byte_size) << endl;
		cout << "   memory bit array                   :     " << Tools::filesize_format(km_byte_size) << endl;
		cout << "   memory rest map                    :     " << Tools::filesize_format(map_byte) << endl; //last_map=44byte
		cout << "   total memory                       :     " << Tools::filesize_format(total_byte_size) << endl;
		cout << "   build time cost                    :     " << build_time_cost << endl;

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
		delete[] dic_num;
		//double total_c = 0;
		//for (int i = 0; i < n_bits; i++) {
		//	double c2 = 0, byte_size = bits_data[i].bit_array_length >> 3;
		//	for (uint64_t p = 0; p < byte_size; p++) {
		//		int idx = bits_data[i].bit_array_2[p];
		//		c2 += dic_num[idx];
		//	}
		//	total_c += c2;
		//	printf("   no.%d bit array usage rate          :     %0.3f\n", i + 1, c2 / bits_data[i].bit_array_length);
		//}

	}


	//save model
	void save(string save_dir) {
		////save the header of the model
		ofstream fout(save_dir + "/header");
		fout << "number_hash " << n_hash << endl;
		fout << "number_bit " << n_bits << endl;
		fout << "ci " << ci << endl;
		fout << "cs " << cs << endl;
		fout.close();
		////save km_bin
		//write kmer_counts,bit_bf and bit_bf_back
		FILE *fp_km = fopen((save_dir + "/km.bin").c_str(), "wb");
		fwrite(&km_kmercount, sizeof(uint64_t), 1, fp_km);
		for (int i = 0; i < bf_num; i++) {
			fwrite(&kmer_counts[i], sizeof(uint64_t), 1, fp_km);
		}
		for (int i = 0; i < bf_num; i++) {
			fwrite(bit_bf[i], sizeof(uint8_t), byte_bf[i], fp_km);
			fwrite(bit_bf_back[i], sizeof(uint8_t), byte_bf_back[i], fp_km);
		}
		//write km_back_bloomfilter
		fwrite(km_back, sizeof(uint8_t), byte_km_back, fp_km);
		//write km_bit_array
		for (int i = 0; i < n_bits; i++) {
			//2.1hash_function
			//fwrite(bits_data[i].hash_seed, sizeof(uint32_t), n_hash, fp_km);
			//bit_array_1
			fwrite(bits_data[i].bit_array_1, sizeof(uint8_t), km_byte_size, fp_km);
			//bit_array_2
			fwrite(bits_data[i].bit_array_2, sizeof(uint8_t), km_byte_size, fp_km);
		}
		fclose(fp_km);
		////save last_map
		kld->save_file(save_dir + "/rest.bin");
	}

	//load model
	void load(string save_dir) {
		//load kmer_count
		FILE *fp_km = fopen((save_dir + "/km.bin").c_str(), "rb");
		fread(&km_kmercount, sizeof(uint64_t), 1, fp_km);

		for (int i = 0; i < bf_num; i++) {
			fread(&kmer_counts[i], sizeof(uint64_t), 1, fp_km);
		}
		init_bf_parameter();
		for (int i = 0; i < bf_num; i++) {
			fread(bit_bf[i], sizeof(uint8_t), byte_bf[i], fp_km);
			fread(bit_bf_back[i], sizeof(uint8_t), byte_bf_back[i], fp_km);
		}
		init_km_parameter(km_kmercount);
		//load km_back_bloomfilter
		fread(km_back, sizeof(uint8_t), byte_km_back, fp_km);
		//load KModel
		for (int i = 0; i < n_bits; i++) {
			fread(bits_data[i].bit_array_1, sizeof(uint8_t), km_byte_size, fp_km);
			//bit_array_2
			fread(bits_data[i].bit_array_2, sizeof(uint8_t), km_byte_size, fp_km);
		}
		fclose(fp_km);
		//load last_map
		kld = new KRestData();
		kld->from_file(save_dir + "/rest.bin");
	}

private:
	OccuBin* occu_bin;
	int ci;
	int cs;
	int kmer_length;
	int n_hash;
	int n_bits;
	uint64_t total_kmer_count = 0;
	int bf_num = 1;
	int bf_index[3] = { 1,0,2 };

	//*bloomfilter and k-2 back bloomfilter*//
	uint64_t kmer_counts[3] = { 0 }; 
	uint8_t **bit_bf;
	uint8_t **bit_bf_back;
	uint64_t* byte_bf;
	uint64_t* length_bf;
	uint64_t* byte_bf_back;
	uint64_t* length_bf_back;
	uint64_t bf_kmercount = 0;
	int bf_num_hash = 6;
	int bf_back_num_hash = 5; 

	//*kmodel*//
	double build_time_cost;
	uint64_t km_kmercount = 0;
	uint64_t km_byte_size; //array in byte size
	uint64_t km_bit_size; //km_bit_size=km_byte_size*8

		
	uint8_t* km_back;//k-2_kmer bloomfilter in km
	int km_back_num_hash = 5;
	uint64_t byte_km_back, bit_km_back;

	BitSaveData* bits_data;
	KRestData* kld;

	//**km insertion buffer**//
	KmerBuff **km_buff; //buff array used by thread
	uint32 bucket_size = 1 << 18;
	int* km_buff_num; //the real length of buffer
	int km_buff_size;// bucket_size*this->n_bits;
	uint32_t km_buff_idx = 0;

	//**blommfilter insertion buffer**//
	KmerBuff* bf_buff = new KmerBuff[BLOACK_SIZE]; //buff array used by thread
	uint32_t bf_buff_idx = 0;


	int kmer_to_bin(string kmer, int occ) {
		vector<int> v_bin = find_bitarray(kmer);
		int len_v_bin = v_bin.size();
		if (len_v_bin == 0) {
			return occ; //caused by false positive,may be in one of the bloomfilter
		}
		if (len_v_bin == 1) {
			if (occ) {
				vector<int> v_candidates = get_neighbor_kmer_bin(kmer);
				//if (v_candidates.size() <= 0)return 0;
				int cnt_bf = 0;
				for (auto v : v_candidates)
					if (v < ci + bf_num) cnt_bf++;
				if (cnt_bf >= v_candidates.size() / 2) return occ;
			}
			return v_bin[0];
		}
		//k can be found in more than one coupling-bit arrays.
		vector<int> v_candidates = get_neighbor_kmer_bin(kmer);
		int len_v_can = v_candidates.size();
		if (len_v_can <= 0) { //none candidates=FP
							  //cout << kmer << " " << Tools::get_complementation(kmer) << " " << occ << endl;///
			return 0;
		}
		int min_dist = (2 << 20), best_bin = v_bin[0];
		for (int i = 0; i < len_v_bin; i++) {
			int cur_dist = 0, cur_min = (2 << 20);
			for (int j = 0; j < len_v_can; j++) {
				cur_dist = abs(v_bin[i] - v_candidates[j]);
				if (cur_dist < cur_min) cur_min = cur_dist;
			}
			if (min_dist > cur_min) {
				min_dist = cur_min;
				best_bin = v_bin[i];
			}
		}
		return best_bin;
	}

	//check the neighbor of kmer
	void get_candidates(string kmer, vector<int> &v_candidates) {
		string min_kmer = Tools::get_min_kmer(kmer);
		int v_from_map = kld->check_kmer(min_kmer); //get_kmer_from_map(min_kmer);
		if (v_from_map > 0) {
			v_candidates.push_back(occu_bin->occ_to_bin(v_from_map));
			return;
		}
		int occ = check_all_bf(min_kmer);
		if (occ != 0) {
			v_candidates.push_back(occ);
			return;
		}
		if (check_back_bloomfilter(min_kmer, km_back, bit_km_back, km_back_num_hash)) {
			int v = find_bitarray_one(min_kmer);
			if (v > -1) v_candidates.push_back(v);
		}
	}

	vector<int> get_neighbor_kmer_bin(string kmer) {
		vector<int> v_candidates;
		string t1_kmer, t2_kmer;
		int kmer_len = kmer.length();
		t1_kmer = kmer.substr(1); //remove the first char
		t2_kmer = kmer.substr(0, kmer_len - 1);//remove the last char 
		for (int i = 0; BASE_CHAR[i]; i++) { //shift forward
			string nw_kmer = t1_kmer + BASE_CHAR[i];
			get_candidates(nw_kmer, v_candidates);
		}
		for (int i = 0; BASE_CHAR[i]; i++) { //shift back
			string nw_kmer = BASE_CHAR[i] + t2_kmer;
			get_candidates(nw_kmer, v_candidates);
		}
		return v_candidates;
	}

	int check_all_bf(string kmer) {
		for (int j = 0; j < bf_num; j++) {
			int i = ci == 1 ? j : bf_index[j]; //3 in the first£¨bloomfilter 3 2 4£©
			bool b_bf = check_bloomfilter(kmer, bit_bf[i], length_bf[i], bf_num_hash);
			bool b_bf_back = check_back_bloomfilter(kmer, bit_bf_back[i], length_bf_back[i], bf_back_num_hash);
			if (b_bf && b_bf_back) {
				return i + ci;
			}
		}
		return 0;
	}

	bool check_bloomfilter(std::string kmer, uint8_t* bit_bf, uint64_t bf_length, int num_hash) {
		uint64_t pos;
		int len = kmer.size();
		auto str_kmer = kmer.c_str();
		for (int i = 0; i<num_hash; i++) {
			pos = Tools::murmur_hash64(str_kmer, len, HashSeeds[i]) % bf_length;
			if (check_bit(bit_bf, pos) == 0)
				return false;
		}
		return true;
	}

	//To decrease false positive, check k-2_mer bloomfilter
	bool check_back_bloomfilter(string kmer, uint8_t* bit_bf, uint64_t bf_length, int num_hash) {
		int len = kmer.length();
		string nw_kmer = kmer.substr(1, len - 2);
		return check_bloomfilter(nw_kmer, bit_bf, bf_length, num_hash);
	}

	//open the kmc db
	void kmc_api(string file_kmer_data_base, CKMCFile &kmer_data_base, CKmerAPI &kmer_object) {
		if (!kmer_data_base.OpenForListing(file_kmer_data_base)) {
			cout << "can't open the kmer_data_base " << file_kmer_data_base << endl;
			exit(1);
		}
		kmer_length = kmer_data_base.KmerLength();
		kmer_object = CKmerAPI(kmer_length);
	}

	void init_bf_parameter() {
		bit_bf = new uint8_t*[bf_num];
		bit_bf_back = new uint8_t*[bf_num];
		byte_bf = new uint64_t[bf_num];
		length_bf = new uint64_t[bf_num];
		byte_bf_back = new uint64_t[bf_num];
		length_bf_back = new uint64_t[bf_num];
		for (int i = 0; i < bf_num; i++) {
			//k bloomfilter
			byte_bf[i] = kmer_counts[i] / 5.5 * bf_num_hash;
			length_bf[i] = byte_bf[i] << 3;
			bit_bf[i] = new uint8_t[byte_bf[i]]{ 0 };
			//k-2 back
			byte_bf_back[i] = (kmer_counts[i] >> 3)* bf_back_num_hash;
			length_bf_back[i] = byte_bf_back[i] << 3;
			bit_bf_back[i] = new uint8_t[byte_bf_back[i]]{ 0 };
			bf_kmercount += kmer_counts[i];
		}
	}

	//get the kmer count wihic freq=1 2 3 or 2 3 4 for specific ci
	uint64_t get_km_kmer_count(CKMCFile &kmer_data_base, CKmerAPI kmer_object) {
		uint32 counter;
		while (kmer_data_base.ReadNextKmer(kmer_object, counter)) {
			if (counter < ci + bf_num)
				kmer_counts[counter - ci]++;
		}
		total_kmer_count = kmer_data_base.KmerCount();
		kmer_data_base.RestartListing(); //RestartListing
										 //init bloomfilter
		init_bf_parameter();
		return total_kmer_count - bf_kmercount;
	}

	void init_km_parameter(uint64_t km_kmercount) {
		km_byte_size = (km_kmercount >> 4) * n_hash;
		km_bit_size = km_byte_size << 3;
		byte_km_back = (km_kmercount >> 4)* km_back_num_hash;
		bit_km_back = byte_km_back << 3;
		bits_data = new BitSaveData[n_bits];
		km_back = new uint8_t[byte_km_back]{ 0 };
		for (int i = 0; i < n_bits; i++) {
			BitSaveData bit_data;
			bit_data.bit_array_length = km_bit_size;
			bit_data.bit_array_1 = new uint8_t[km_byte_size]{ 0 };
			bit_data.bit_array_2 = new uint8_t[km_byte_size]{ 0 };
			bit_data.hash_seed = new uint32_t[n_hash];
			//set num_bit hash seed
			for (int j = 0; j < n_hash; j++) {
				int index = (i*n_hash + j) % 128;
				bit_data.hash_seed[j] = HashSeeds[index];
			}
			bits_data[i] = bit_data;
		}
	}

	void init_km_bit(uint64_t km_kmercount) {
		//1. init parameter size 
		init_km_parameter(km_kmercount);
		//2. init km_buff
		km_buff = new KmerBuff*[n_bits];
		km_buff_num = new int[n_bits];
		km_buff_size = bucket_size*n_bits;
		for (int i = 0; i < n_bits; i++) {
			km_buff[i] = new KmerBuff[bucket_size];
			km_buff_num[i] = bucket_size;
		}
		//3.init KRestData
		kld = new KRestData(kmer_length);
	}

	void push_bloomfilter(string kmer, int i) {
		insert_bloomfilter(kmer, bit_bf[i], length_bf[i], bf_num_hash);
		string nw_kmer = kmer.substr(1, kmer.size() - 2);
		insert_bloomfilter(nw_kmer, bit_bf_back[i], length_bf_back[i], bf_back_num_hash);
	}

	void push_to_bloomfilter(string kmer, uint32 occ) {
		bf_buff[bf_buff_idx].kmer = kmer;
		bf_buff[bf_buff_idx++].occ = occ;
		if (bf_buff_idx >= BLOACK_SIZE) {
#pragma omp parallel for num_threads(4)
			for (int i = 0; i < BLOACK_SIZE; ++i) {
				push_bloomfilter(bf_buff[i].kmer, bf_buff[i].occ - ci);
			}
			bf_buff_idx = 0;
		}
	}

	void push_last_to_bloomfilter() {
#pragma omp parallel for num_threads(4)
		for (int i = 0; i < bf_buff_idx; ++i) {
			push_bloomfilter(bf_buff[i].kmer, bf_buff[i].occ - ci);
		}
	}

	void insert_bloomfilter(string kmer, uint8_t* bit_bf, uint64_t bf_length, int num_hash) {
		uint64_t pos;
		int len = kmer.size();
		auto str_kmer = kmer.c_str();
		for (int i = 0; i<num_hash; i++) {
			pos = Tools::murmur_hash64(str_kmer, len, HashSeeds[i]) % bf_length;
			set_bit(bit_bf, pos);
		}
	}

	void push_to_array(string kmer, uint32 occ) {
		int row = km_buff_idx / bucket_size;
		int col = km_buff_idx%bucket_size;
		km_buff[row][col].kmer = kmer;
		km_buff[row][col].occ = occ;
		km_buff_idx++;
		if (km_buff_idx >= km_buff_size) {
			insert_with_thread(km_buff, km_buff_num, n_bits, bucket_size);
			km_buff_idx = 0;
		}
	}

	void push_last_to_array(int n_thread) {
		int row = (km_buff_idx - 1) / bucket_size;
		int col = (km_buff_idx - 1) % bucket_size;
		km_buff_num[row] = col + 1;
		for (int i = row + 1; i < n_thread; i++)
			km_buff_num[i] = 0;
		insert_with_thread(km_buff, km_buff_num, n_thread, bucket_size);
	}

	int reorder_buffer(KmerBuff* a, int n) {
		int il = 0, ir = n - 1;
		while (il < ir) {
			while (il < ir && !a[ir].occ) ir--;
			while (il < ir && a[il].occ) il++;
			if (il < ir) {
				a[il] = a[ir];
				a[ir].occ = 0;
			}
		}
		return a[il].occ ? il + 1 : 0; //il==0 -> All inserted successfully
	}

	//push buffer data into bit_array
	void insert_array(KmerBuff* buff, int index, int &buff_n) {
		for (int c = 0; c < buff_n; c++) {
			bool flag = insert_to_array(buff[c].kmer, occu_bin->occ_to_bin(buff[c].occ), index);
			if (flag) {
				//for kmodel_k-2_mer
				string nw_kmer = buff[c].kmer.substr(1, buff[c].kmer.size() - 2);
				//nw_kmer = Tools::get_min_kmer(nw_kmer);
				insert_bloomfilter(nw_kmer, km_back, bit_km_back, km_back_num_hash);
				buff[c].occ = 0;
			}
		}
		buff_n = reorder_buffer(buff, buff_n);
	}

	void insert_with_thread(KmerBuff** buff, int* buff_real_n, int n_thread, int bucket_size) {
		int i, j;
		//for kmodel
		for (int t = 0; t < n_thread; t++) {
#pragma omp parallel for num_threads(n_thread)
			for (i = 0; i < n_thread; ++i) {
				insert_array(buff[i], (i + t) % n_thread, buff_real_n[i]);
			}
		}
		//save the last items into map
		for (i = 0; i < n_thread; ++i) {
			for (j = 0; j < buff_real_n[i]; ++j) {
				kld->push_back(buff[i][j].kmer, (int)buff[i][j].occ);
			}
			buff_real_n[i] = bucket_size;
		}
	}

	//set the position in the bit to one
	inline void set_bit(uint8_t *bit, uint64_t pos) {
		uint64_t row = pos >> 3;
		uint64_t col = pos & 0x7;
		uint8_t x = 0x1 << (8 - col - 1); //set the bit in the column in the table 
		__sync_fetch_and_or(bit + row, x);//bit[row] |= x;
	}

	//check the position in the bit is one or zero 
	inline bool check_bit(const uint8_t *bit, uint64_t pos) {
		uint64_t row = pos >> 3;
		uint64_t col = pos & 0x7;
		return (bit[row] >> (8 - col - 1)) & 0x1;
	}

	bool insert_to_array(string kmer, uint32_t occ, int index) {
		//make bit1_v and bit2_v class memeber
		uint64_t* bit1_v = new uint64_t[n_hash];// the binary value of frequency in bit_array_1 
		uint64_t* bit2_v = new uint64_t[n_hash];// the index of bit_array_2
		uint32_t*  hash_seed = bits_data[index].hash_seed;
		uint8_t* bit_array_1 = bits_data[index].bit_array_1;
		uint8_t* bit_array_2 = bits_data[index].bit_array_2;
		uint64_t array_length = bits_data[index].bit_array_length;
		for (int i = 0; i < n_hash; i++) {
			bit1_v[i] = occ & 0x1;
			bit2_v[i] = Tools::murmur_hash64(kmer.c_str(), kmer.size(), hash_seed[i]) % array_length;
			occ >>= 1;
		}
		bool ok = true; //the flag if this kmer can be inserted
		for (int i = 0; i < n_hash; i++) {
			uint8_t v1 = check_bit(bit_array_1, bit2_v[i]);
			uint8_t v2 = check_bit(bit_array_2, bit2_v[i]);
			if (v2 && v1 != bit1_v[i]) {
				ok = false;  break;
			}
		}
		if (ok) {
			for (int i = 0; i < n_hash; i++) {
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
				if (result > 0)v_bin.push_back(result);
			}
		}
		delete[] occ_bit;
		return  v_bin;
	}


	//travel the bit array and just get one occurrence of the kmer
	int find_bitarray_one(string kmer) {
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
				if (result != 0) break;
			}
		}
		delete[] occ_bit;
		return  result;
	}
};

KModel* get_model(int ci = 1, int cs = 1023, int num_hash = 7, int num_bit = 5) {
	OccuBin* occu_bin = new OccuBin(cs + 1, num_hash);
	return new KModel(occu_bin, num_bit, ci);
}


KModel* get_model(string save_dir) {
	ifstream fin(save_dir + "/header");
	if (!fin) {
		cout << "load_model: cant't open the header of the model !\n";
		exit(1);
	}
	string t_str;
	int n_hash, n_bits, ci, cs;
	fin >> t_str >> n_hash;
	fin >> t_str >> n_bits;
	fin >> t_str >> ci;
	fin >> t_str >> cs;
	fin.close();
	KModel* km = get_model(ci, cs, n_hash, n_bits);
	km->load(save_dir);
	return km;
}


#endif