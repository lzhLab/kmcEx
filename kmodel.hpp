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
#include <chrono>
#include "omp.h"

using namespace std;


struct KmerBuff {
	string kmer;
	uint32_t occ;
};

const int BLOACK_SIZE = 1048576;
KmerBuff kmer_buff[BLOACK_SIZE];

int THREAD_NUM = 4;

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
		this->n_hash = occu_bin->get_nhash();
		this->n_bits = n_bits;
	}

	// n_bit:how many bit_array used
	void init_KModel(string file_kmer_data_base) {
		string kmer;
		float counter;
		uint8_t occ8;
		int occ32;
		CKMCFile kmer_data_base;
		if (!kmer_data_base.OpenForListing(file_kmer_data_base)) {
			std::cout << "init_occ_count()__can't open the kmer_data_base\n";
			return;
		}
		uint32 _kmer_length = kmer_data_base.KmerLength();
		this->total_kmer_count = kmer_data_base.KmerCount();
		CKmerAPI kmer_object(_kmer_length);
		placeholder(kmer_data_base, kmer_object);
		init_bit_aray();
		show_header_info();	
		uint64_t idx = 0;
		auto start = chrono::high_resolution_clock::now();//////////////
		while (kmer_data_base.ReadNextKmer(kmer_object, counter)) {
			occ32 = counter;
			++occ_count[occ32];
			//kmer = kmer_object.to_string().c_str();
			kmer_buff[idx].kmer = kmer_object.to_string().c_str();
			kmer_buff[idx++].occ = occ32;
			if (idx >= BLOACK_SIZE) {
				#pragma omp parallel for num_threads(THREAD_NUM)
				for (int i = 0; i < idx; ++i) {
					insert_pre(kmer_buff[i].kmer, kmer_buff[i].occ);
				}
				idx = 0;
			}
		}
		#pragma omp parallel for num_threads(THREAD_NUM)
		for (int i = 0; i < idx; ++i) {
			insert_pre(kmer_buff[i].kmer, kmer_buff[i].occ);
		}
		chrono::duration<double> dur = chrono::high_resolution_clock::now() - start;
		cout << "read kmer and init kmodel time: " << dur.count() << endl;
		kmer_data_base.Close();
	}

	//if this kmer can't be found,it will return 0
	virtual int kmer_to_occ(string kmer) {
		int occ = get_kmer_from_map(kmer);
		if (occ != -1)
			return occ;
		int bin = kmer_to_bin(kmer);
		occ = occu_bin->bin_to_mean(bin);
		return occ;
	}

	//if this kmer can't be found,it will return 0
	//t_num: the number of threads
	vector<int> kmer_to_occ(vector<string> kmer_v,int t_num=4) {
		int n = kmer_v.size();
		vector<int> occ_v(n);
		#pragma omp parallel for num_threads(t_num) 
		for (int i = 0; i < n; i++) {
			occ_v[i] = kmer_to_occ(kmer_v[i]);
		}
		return occ_v;
	}

	//alter the char in the pos as the neighbor of kmer
	void shift_one(string kmer, vector<int>  &v_candidates) {
		int v = find_bitarray_one(kmer);
		if (v > -1) {
			v_candidates.push_back(v);
		}
		else { // try kmer's complementation
			string nw_kmer_com = Tools::get_complementation(kmer);
			v = find_bitarray_one(nw_kmer_com);
			if (v >-1) {
				v_candidates.push_back(v);
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

	int kmer_to_bin(string kmer) {
		vector<int> v_bin = find_bitarray(kmer);
		int len_v_bin = v_bin.size();
		if (len_v_bin == 0) {
			//cout << len_v_bin << " : " << kmer << endl;
			//cout << "v_bin.size == 0" << endl;
			return 0;
		}
		if (len_v_bin == 1) {
			return v_bin[0]; //only one
		}
		//k can be found in more than one coupling-bit arrays.
		vector<int> v_candidates = get_neighbor_kmer_bin(kmer);
		int len_v_can = v_candidates.size();
		if (len_v_can <= 0 && len_v_bin > 1)//no any candidates,get the minimum bin
			return Tools::vector_min(v_bin);
		//get the best bin
		int min_dist = 999999, best_bin = v_bin[0];
		for (int i = 0; i < len_v_bin; i++) {
			int cur_dist = 0;
			for (int j = 0; j < len_v_can; j++) {
				cur_dist += abs(v_bin[i] - v_candidates[j]);
			}
			if (min_dist > cur_dist && v_bin[i] != 0) {
				min_dist = cur_dist;
				best_bin = v_bin[i];
			}
		}
		return best_bin;
	}

	//find in the map
	int get_kmer_from_map(string kmer) {
		uint64_t key = Tools::kmers2uint64(kmer);
		unordered_map<uint64_t, uint32_t>::iterator it = last_map.find(key);
		if (it != last_map.end()) {
			return it->second;
		}
		return -1;
	}

	//save model
	virtual void save_model(string save_dir) {
		//1. save model parameters:n_hash,n_bits,kmer_count
		ofstream fout(save_dir + "/param.conf");
		fout << "number_hash " << n_hash << endl;
		fout << "number_bit " << n_bits << endl;
		fout << "total_kmer_count  " << total_kmer_count << endl;
		fout << "once_kmer_count  " << once_kmer_count << endl;
		fout << "last_map_size  " << last_map.size() << endl;
		fout << "C_MAX   " << C_MAX << endl;
		fout.close();
		//2. save OccuBin
		FILE *fp_occ_meta = fopen((save_dir + "/occ.bin").c_str(), "wb");
		fwrite(occu_bin->get_occ_bin_meta(), sizeof(OccuBinMeta), C_MAX, fp_occ_meta);
		fclose(fp_occ_meta);
		//3. save KModel
		FILE *fp_hash = fopen((save_dir + "/hash.bin").c_str(), "wb");
		FILE *fp_bit1 = fopen((save_dir + "/bit1.bin").c_str(), "wb");
		FILE *fp_bit2 = fopen((save_dir + "/bit2.bin").c_str(), "wb");
		for (int i = 0; i < n_bits; i++) {
			//3.1 hash_function
			fwrite(bits_data[i].hash_seed, sizeof(uint32_t), n_hash, fp_hash);
			//3.2 bit_array_1
			fwrite(bits_data[i].bit_array_1, sizeof(uint8_t), bits_data[i].bit_array_length >> 3, fp_bit1);
			//3.3 bit_array_2
			fwrite(bits_data[i].bit_array_2, sizeof(uint8_t), bits_data[i].bit_array_length >> 3, fp_bit2);
		}
		fclose(fp_hash);
		fclose(fp_bit1);
		fclose(fp_bit2);

		//4. save last_map
		FILE *fp_map_w = fopen((save_dir + "/last_map.bin").c_str(), "wb");
		unordered_map<uint64_t, uint32_t>::iterator iter;
		for (iter = last_map.begin(); iter != last_map.end(); iter++) {
			fwrite(&iter->first, sizeof(uint64_t), 1, fp_map_w);
			fwrite(&iter->second, sizeof(uint32_t), 1, fp_map_w);
		}
		fclose(fp_map_w);

		//5.save bloomfilter for one-mer
		//FILE *fp_bloom = fopen((save_dir + "/bloom.bin").c_str(), "wb");
		//fwrite(bit_bf1, sizeof(uint8_t),byte_bf1, fp_bloom);
		//fclose(fp_bloom);
	}

	//load model
	virtual void load_model(string save_dir) {
		//1. load model parameters:n_hash,n_bits,kmer_count
		ifstream fin(save_dir + "/param.conf");
		if (!fin) {
			cout << "load_model: cant't open the param.conf !\n";
			return;
		}
		string t_str;
		fin >> t_str >> n_hash;
		fin >> t_str >> n_bits;
		fin >> t_str >> total_kmer_count;
		fin >> t_str >> once_kmer_count;
		fin >> t_str >> left_kmer_count;
		fin >> t_str >> C_MAX;
		set_bit_array_size();
		fin.close();

		//2. load OccuBin
		occu_bin = new OccuBin(n_hash);
		FILE *fp_occ_meta = fopen((save_dir + "/occ.bin").c_str(), "rb");
		fread(occu_bin->get_occ_bin_meta(), sizeof(OccuBinMeta), C_MAX, fp_occ_meta);
		fclose(fp_occ_meta);
		//3. load KModel
		uint64_t t_byte_array_size = byte_array_size;
		bits_data = new BitSaveData[n_bits];
		FILE *fp_hash = fopen((save_dir + "/hash.bin").c_str(), "rb");
		FILE *fp_bit1 = fopen((save_dir + "/bit1.bin").c_str(), "rb");
		FILE *fp_bit2 = fopen((save_dir + "/bit2.bin").c_str(), "rb");
		for (int i = 0; i < n_bits; i++) {
			//3.1 hash_function
			bits_data[i].hash_seed = new uint32_t[n_hash];
			bits_data[i].bit_array_length = t_byte_array_size << 3;
			fread(bits_data[i].hash_seed, sizeof(uint32_t), n_hash, fp_hash);

			//3.2 bit_array_1
			bits_data[i].bit_array_1 = new uint8_t[t_byte_array_size];
			fread(bits_data[i].bit_array_1, sizeof(uint8_t), t_byte_array_size, fp_bit1);

			//3.3 bit_array_2
			bits_data[i].bit_array_2 = new uint8_t[t_byte_array_size];
			fread(bits_data[i].bit_array_2, sizeof(uint8_t), t_byte_array_size, fp_bit2);
			t_byte_array_size = t_byte_array_size >> 1;
		}
		fclose(fp_hash);
		fclose(fp_bit1);
		fclose(fp_bit2);

		//4. load last_map
		uint64_t k;
		int v;
		FILE *fp_map_r = fopen((save_dir + "/last_map.bin").c_str(), "rb");
		for (int i = 0; i < left_kmer_count; i++) {
			fread(&k, sizeof(uint64_t), 1, fp_map_r);
			fread(&v, sizeof(uint32_t), 1, fp_map_r);
			last_map.insert(make_pair(k, v));
		}
		fclose(fp_map_r);

		////5.load bloomfilter for one-mer
		//FILE *fp_bloom = fopen((save_dir + "/bloom.bin").c_str(), "rb");
		//bit_bf1 = new uint8_t[this->byte_bf1];
		//fread(bit_bf1, sizeof(uint8_t), byte_bf1, fp_bloom);
		//fclose(fp_bloom);
	}

	virtual void show_header_info() {
		cout << "[kmodel_number_hash]:" << n_hash << endl;
		cout << "total_kmer_count:" << total_kmer_count << endl;
		//cout << "kmercount sent into bloomfilter01:" << once_kmer_count << endl;
		//cout << "kmercount sent into kmodel:" << kmodel_kmer_count << endl;
		cout << "kmodel_number_bit_array:" << n_bits << endl;
		cout << "kmodel_byte_array_size:" << byte_array_size << endl;
	}

	void show_model_info() {
		cout << "kmodel_distance_error:" << calc_distance_error() << endl;
		uint64_t total_byte_array_size = 0;
		for (int i = 0; i < 4; i++) {
			total_byte_array_size += bits_data[i].bit_array_length >> 3;
		}
		double left = kmodel_kmer_count; //kmodel_kmer_count
		for (int i = 0; i < n_bits; i++) {
			printf("%02d--save_%d--left_%.0f---saveratio_%.3f\n", i + 1, c_left[i], left - c_left[i], c_left[i] / left);
			left -= c_left[i];
		}
		cout << "the left kmer is saved to map:" << left_kmer_count << endl;
		//cout << "memory used in bloomfilter01:" << Tools::filesize_format(byte_bf1) << endl;
		cout << "memory used in kmodel bit array:" << Tools::filesize_format(total_byte_array_size * 2) << endl;
		cout << "total memory used:" << Tools::filesize_format(total_byte_array_size * 2) << endl;

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


	int get_left_kmer_count() {
		return left_kmer_count;
	}

	OccuBin* get_occu_bin() {
		return occu_bin;
	}

	//the kmer_count that the kmer only appears once
	uint64_t get_once_kmer_count() {
		return this->once_kmer_count;
	};

	//the kmer_count stored in kmodel
	uint64_t get_kmodel_kmer_count() {
		return this->kmodel_kmer_count;
	}


protected:
	OccuBin* occu_bin;
	int n_hash;
	int n_bits;
	int left_kmer_count = 0; //the kmer count that can't save to bit_array and save to last_map
	uint64_t* occ_count; //statastic occurrence
	uint64_t total_kmer_count; //total_kmer_count counted kmercount 
	uint64_t once_kmer_count = 0; //the kmer_count that the kmer only appears once
	uint64_t kmodel_kmer_count; //the kmer_count stored in kmodel
	uint64_t byte_array_size; //byte size
	uint64_t bit_array_length; //the real bit table length,bit_array_length=byte_array_size*8
	unordered_map<uint64_t, uint32_t> last_map; //save the left kmer
	BitSaveData* bits_data;
	int* c_left; //how many kmer saved in the bit_array

	virtual void set_bit_array_size() {
		this->kmodel_kmer_count = total_kmer_count - once_kmer_count;
		//to minimize the false positive -> bit_array_length=kmer_count*n_hash/ln2
		this->byte_array_size = kmodel_kmer_count / 5.52 * n_hash;
		this->bit_array_length = this->byte_array_size << 3;
	}

	//set the position in the bit to one
	inline void set_bit(uint8_t *bit, uint64_t pos) {
		uint64_t row = pos >> 3;
		uint64_t col = pos & 0x7;
		uint8_t x = 0x1 << (8 - col - 1); //set the bit in the column in the table 
		bit[row] |= x;
	}

	//check the position in the bit is one or zero 
	inline  bool check_bit(const uint8_t *bit, uint64_t pos) {
		uint64_t row = pos >> 3;
		uint64_t col = pos & 0x7;
		return (bit[row] >> (8 - col - 1)) & 0x1;
	}

	virtual void placeholder(CKMCFile& kmer_data_base, CKmerAPI kmer_object) {
		//do somthing in child class
		//cout << "placeholder" << endl;
	}

	virtual void init_bit_aray() {
		//init kmodel bit array
		bits_data = new BitSaveData[n_bits];
		c_left = new int[n_bits];
		memset(c_left, 0, sizeof(int)*n_bits);
		set_bit_array_size();
		uint64_t t_byte_array_size = this->byte_array_size;
		for (int i = 0; i < n_bits; i++) {
			BitSaveData bit_data;
			bit_data.bit_array_length = t_byte_array_size << 3;
			bit_data.bit_array_1 = new uint8_t[t_byte_array_size];
			bit_data.bit_array_2 = new uint8_t[t_byte_array_size];
			memset(bit_data.bit_array_1, 0, sizeof(uint8_t)*t_byte_array_size);
			memset(bit_data.bit_array_2, 0, sizeof(uint8_t)*t_byte_array_size);
			bit_data.hash_seed = new uint32_t[n_hash];
			//set num_bit hash seed
			for (int j = 0; j < n_hash; j++) {
				int index = (i*n_hash + j) % 32;
				bit_data.hash_seed[j] = HashSeeds[index];
			}
			bits_data[i] = bit_data;
			t_byte_array_size = t_byte_array_size >> 1; //byte_array_size reduce
		}
		occ_count = new uint64_t[C_MAX];
	}


	virtual void insert_pre(string kmer, int occ) {
		#pragma omp critical
		{
			auto occ8 = (uint8_t)(occu_bin->occ_to_bin(occ));
			int index = insert_to_array(kmer, occ8);//get_insert_idx(kmer, occ8);
			if (index > -1) { // can be inserted
				c_left[index]++;//couting how many kmer save in the bit_array
			}
			else { // store kmer_occ using map
				uint64_t key = Tools::kmers2uint64(kmer);
				last_map.insert(make_pair(key, (uint32_t)occ));
				left_kmer_count++;
			}
		}
	}

	//----------------------------------------------------------------------------------
	//if the bit_array_2'value is one and the bit in the occurrence is different from the bit_array_1' value,
	//then we can insert this kmer into bit array,or write it to the disk
	//--------------------------------------------------------------------------------
	int insert_to_array(std::string kmer, uint8_t occ8) {
		uint64_t* bit1_v = new uint64_t[n_hash];// the binary value of frequency in bit_array_1 
		uint64_t* bit2_v = new uint64_t[n_hash];// the index of bit_array_2
		int idx = -1;
		for (int index = 0; index < this->n_bits; index++) {
			uint8_t occ8_t = occ8;
			uint32_t*  hash_seed = bits_data[index].hash_seed;
			uint8_t* bit_array_1 = bits_data[index].bit_array_1;
			uint8_t* bit_array_2 = bits_data[index].bit_array_2;
			uint64_t array_length = bits_data[index].bit_array_length;
			for (int i = 0; i < n_hash; i++) {
				bit1_v[i] = occ8_t & 0x1;
				bit2_v[i] = Tools::murmur_hash64(kmer.c_str(), kmer.size(), hash_seed[i]) % array_length;
				occ8_t >>= 1;
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
				idx = index; break;
			}
		}
		delete[] bit1_v;
		delete[] bit2_v;
		return idx;
	}

	double calc_distance_error() {
		double sum_dis = 0, sum_kmer_count = 0;
		for (int i = 2; i < C_MAX; i++) {
			sum_kmer_count += i*occ_count[i];
		}
		int start_bin = this->occu_bin->get_bin_end_index1();
		auto occ_bin_meta = this->occu_bin->get_occ_bin_meta();
		for (int i = start_bin; i < C_MAX; i++) {
			sum_dis += abs(occ_bin_meta[i].occ_mean - i)*occ_count[i];
		}
		double dis_error = sum_dis / sum_kmer_count;
		return dis_error;
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
				if (result>0)
					v_bin.push_back(result);
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

class KModelOne: public KModel
{
public:
	KModelOne() {}
	KModelOne(OccuBin* occu_bin, int n_bits): KModel(occu_bin,n_bits)
	{	
	}

	int kmer_to_occ(string kmer) {
		int occ = get_kmer_from_map(kmer);
		if (occ != -1)
			return occ;
		if (check_bloomfilter(kmer, bit_bf1, length_bf1))
			return 1;
		int bin = kmer_to_bin(kmer);
		occ = occu_bin->bin_to_mean(bin);
		return occ;
	}

	bool check_bloomfilter(std::string kmer, uint8_t* bit_bf, uint64_t bf_length) {
		uint64_t pos;
		for (int i = 0; i<bf_n_hash; i++) {
			pos = Tools::murmur_hash64(kmer.c_str(), kmer.size(), HashSeeds[i]) % bf_length;
			if (check_bit(bit_bf, pos) == 0)
				return false;
		}
		return true;
	}

	void insert_pre(string kmer, int occ) {
		if (occ == 1)
			insert_to_bloomfilter(kmer, this->bit_bf1, this->length_bf1);
		else
		{
			#pragma omp critical
			{
				auto occ8 = (uint8_t)(occu_bin->occ_to_bin(occ));
				int index = insert_to_array(kmer, occ8);//get_insert_idx(kmer, occ8);
				if (index > -1) { // can be inserted
					c_left[index]++;//couting how many kmer save in the bit_array
				}
				else { // store kmer_occ using map
					uint64_t key = Tools::kmers2uint64(kmer);
					last_map.insert(make_pair(key, (uint32_t)occ));
					left_kmer_count++;
				}
			}
		}
	}

	//save model
	void save_model(string save_dir) {
		KModel::save_model(save_dir);
		//save bloomfilter for one-mer
		FILE *fp_bloom = fopen((save_dir + "/bloom.bin").c_str(), "wb");
		fwrite(bit_bf1, sizeof(uint8_t),byte_bf1, fp_bloom);
		fclose(fp_bloom);
	}

	//load model
	void load_model(string save_dir) {
		KModel::load_model(save_dir);
		//load bloomfilter for one-mer
		FILE *fp_bloom = fopen((save_dir + "/bloom.bin").c_str(), "rb");
		bit_bf1 = new uint8_t[this->byte_bf1];
		fread(bit_bf1, sizeof(uint8_t), byte_bf1, fp_bloom);
		fclose(fp_bloom);
	}

	bool check_bloomfilter01(std::string kmer) {
		return check_bloomfilter(kmer, this->bit_bf1, this->length_bf1);
	}


protected:
	uint8_t *bit_bf1; //bloomfilter for one occurrence
	uint64_t byte_bf1,length_bf1;
	int bf_n_hash = 7;

	void init_bit_aray() {
		KModel::init_bit_aray();
		//init bloomfilter
		bit_bf1 = new uint8_t[byte_bf1];
		memset(bit_bf1, 0, sizeof(uint8_t)*byte_bf1);
	}

	void set_bit_array_size() {
		KModel::set_bit_array_size();
		this->byte_bf1 = once_kmer_count / 5.52 * bf_n_hash;
		this->length_bf1 = this->byte_bf1 << 3;
	}

	void insert_to_bloomfilter(std::string kmer, uint8_t* &bit_bf, uint64_t bf_length) {
		uint64_t pos;
		int len = kmer.size();
		auto str_kmer = kmer.c_str();
		for (int i = 0; i<bf_n_hash; i++) {
			pos = Tools::murmur_hash64(str_kmer, len, HashSeeds[i]) % bf_length;
			set_bit(bit_bf, pos);
		}
	}

	void placeholder(CKMCFile &kmer_data_base, CKmerAPI kmer_object) {
		cout << "placeholder" << endl;
		float counter;
		auto start = chrono::high_resolution_clock::now();
		//Get the one occurence
		while (kmer_data_base.ReadNextKmer(kmer_object, counter)) {
			if (counter == 1)
				++once_kmer_count;
		}
		chrono::duration<double> dur = chrono::high_resolution_clock::now() - start;
		cout << "read kmer time:" << dur.count() << endl;
		kmer_data_base.RestartListing();
	}

	void show_header_info() {
		cout << "[kmodel_number_hash]:" << n_hash << endl;
		cout << "total_kmer_count:" << total_kmer_count << endl;
		cout << "kmercount sent into bloomfilter01:" << once_kmer_count << endl;
		cout << "kmercount sent into kmodel:" << kmodel_kmer_count << endl;
		cout << "kmodel_number_bit_array:" << n_bits << endl;
		cout << "kmodel_byte_array_size:" << byte_array_size << endl;
	}
};


#endif