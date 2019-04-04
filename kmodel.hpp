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
#include <atomic> 

using namespace std;


struct KmerBuff {
	string kmer;
	uint32_t occ;
};


const int BLOACK_SIZE = 1 << 19;
const int BK = 4096; 



struct BitSaveData {
	uint64_t bit_array_length; 
	uint32_t* hash_seed; 
	uint8_t *bit_array_1; 
	uint8_t *bit_array_2; 
};


class KModel
{
public:
	KModel() {}
	KModel(OccuBin* occu_bin, int n_bits) {
		this->occu_bin = occu_bin;
		this->n_hash = occu_bin->get_nhash();
		this->n_bits = n_bits;
		this->km_back_n_hash = this->n_hash - 2;
	}

	virtual void init_KModel(string file_kmer_data_base) {
		uint32 counter;
		CKMCFile kmer_data_base;
		CKmerAPI kmer_object;
		string kmer;
		init_kmc_api(file_kmer_data_base, kmer_data_base, kmer_object);
		init_bit_aray(kmer_data_base.KmerCount());
		show_header_info();
		init_buff(n_bits);
		auto start = chrono::high_resolution_clock::now();
		while (kmer_data_base.ReadNextKmer(kmer_object, counter)) {
			kmer = kmer_object.to_string();
			push_to_array(kmer, counter);
		}
		push_last_to_array(n_bits);
		chrono::duration<double> dur = chrono::high_resolution_clock::now() - start;
		cout << "   kmcEx model construction time:        " << dur.count() << endl;
		kmer_data_base.Close();
	}

	void set_load_param(uint64_t once_kmer_count, uint64_t kmodel_kmer_count, int last_map_size) {
		this->once_kmer_count = once_kmer_count;
		this->kmodel_kmer_count = kmodel_kmer_count;
		this->last_map_size = last_map_size;
	}

	virtual int kmer_to_occ(string kmer, uint32_t r_occ = 0) {
		kmer = Tools::get_min_kmer(kmer);
		int occ = get_kmer_from_map(kmer);
		if (occ != -1) return occ;
		if (!check_back_bloomfilter(kmer, km_back, length_km_back, km_back_n_hash)) {
			return 0;
		}
		int bin = kmer_to_bin(kmer, false, r_occ);
		occ = occu_bin->bin_to_mean(bin);
		return occ;
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

	int get_kmer_from_map(string kmer) {
		uint64_t key = Tools::kmers2uint64(kmer);
		unordered_map<uint64_t, uint32_t>::iterator it = last_map.find(key);
		if (it != last_map.end()) {
			return it->second;
		}
		return -1;
	}

	virtual void save_model(string save_dir) {
		ofstream fout(save_dir + "/param.conf");
		fout << "number_hash " << n_hash << endl;
		fout << "number_bit " << n_bits << endl;
		fout << "once_kmer_count  " << once_kmer_count << endl;
		fout << "kmodel_kmer_count  " << kmodel_kmer_count << endl;
		fout << "last_map_size  " << last_map.size() << endl;
		fout << "max_counter   " << occu_bin->get_max_counter() << endl;
		fout.close();
		FILE *fp_hash = fopen((save_dir + "/hash.bin").c_str(), "wb");
		FILE *fp_bit1 = fopen((save_dir + "/bit1.bin").c_str(), "wb");
		FILE *fp_bit2 = fopen((save_dir + "/bit2.bin").c_str(), "wb");
		for (int i = 0; i < n_bits; i++) {
			fwrite(bits_data[i].hash_seed, sizeof(uint32_t), n_hash, fp_hash);
			fwrite(bits_data[i].bit_array_1, sizeof(uint8_t), bits_data[i].bit_array_length >> 3, fp_bit1);
			fwrite(bits_data[i].bit_array_2, sizeof(uint8_t), bits_data[i].bit_array_length >> 3, fp_bit2);
		}
		fclose(fp_hash);
		fclose(fp_bit1);
		fclose(fp_bit2);

		FILE *fp_km_back = fopen((save_dir + "/km_back.bin").c_str(), "wb");
		fwrite(km_back, sizeof(uint8_t), byte_km_back, fp_km_back);
		fclose(fp_km_back);

		FILE *fp_map_w = fopen((save_dir + "/last_map.bin").c_str(), "wb");
		unordered_map<uint64_t, uint32_t>::iterator iter;
		uint64_t buff_kmer[BK];
		uint32_t buff_occ[BK], ix = 0;
		for (iter = last_map.begin(); iter != last_map.end(); iter++) {
			buff_kmer[ix] = iter->first;
			buff_occ[ix++] = iter->second;
			if (ix >= BK) {
				fwrite(buff_kmer, sizeof(uint64_t), ix, fp_map_w);
				fwrite(buff_occ, sizeof(uint32_t), ix, fp_map_w);
				ix = 0;
			}
		}
		if (ix) {
			fwrite(buff_kmer, sizeof(uint64_t), ix, fp_map_w);
			fwrite(buff_occ, sizeof(uint32_t), ix, fp_map_w);
		}
		fclose(fp_map_w);
	}

	virtual void load_model(string save_dir) {
		set_bit_array_size(kmodel_kmer_count);
		uint64_t t_byte_array_size = byte_array_size;
		bits_data = new BitSaveData[n_bits];
		FILE *fp_hash = fopen((save_dir + "/hash.bin").c_str(), "rb");
		FILE *fp_bit1 = fopen((save_dir + "/bit1.bin").c_str(), "rb");
		FILE *fp_bit2 = fopen((save_dir + "/bit2.bin").c_str(), "rb");
		for (int i = 0; i < n_bits; i++) {
			bits_data[i].hash_seed = new uint32_t[n_hash];
			bits_data[i].bit_array_length = t_byte_array_size << 3;
			fread(bits_data[i].hash_seed, sizeof(uint32_t), n_hash, fp_hash);
			bits_data[i].bit_array_1 = new uint8_t[t_byte_array_size];
			fread(bits_data[i].bit_array_1, sizeof(uint8_t), t_byte_array_size, fp_bit1);
			bits_data[i].bit_array_2 = new uint8_t[t_byte_array_size];
			fread(bits_data[i].bit_array_2, sizeof(uint8_t), t_byte_array_size, fp_bit2);
		}
		fclose(fp_hash);
		fclose(fp_bit1);
		fclose(fp_bit2);
		FILE *fp_km_back = fopen((save_dir + "/km_back.bin").c_str(), "rb");
		km_back = new uint8_t[byte_km_back];
		fread(km_back, sizeof(uint8_t), byte_km_back, fp_km_back);
		fclose(fp_km_back);
		uint64_t k[BK];
		uint32 v[BK], bin = last_map_size / BK, last = last_map_size % BK;
		FILE *fp_map_r = fopen((save_dir + "/last_map.bin").c_str(), "rb");
		for (int i = 0; i < bin; i++) {
			fread(k, sizeof(uint64_t), BK, fp_map_r);
			fread(v, sizeof(uint32_t), BK, fp_map_r);
			for (int j = 0; j < BK; j++) {
				last_map.insert(pair<uint64_t, uint32_t>(k[j], v[j]));
			}
		}
		if (last) {
			fread(k, sizeof(uint64_t), last, fp_map_r);
			fread(v, sizeof(uint32_t), last, fp_map_r);
			for (int j = 0; j < last; j++) {
				last_map.insert(pair<uint64_t, uint32_t>(k[j], v[j]));
			}
		}
		fclose(fp_map_r);
	}

	void show_header_info() {
		cout << "   Number of hash functions              :  " << n_hash << endl;
		cout << "   Number of coupled-bit arrays          :  " << n_bits << endl;
		cout << "   Number of k-mers having count=1       :  " << once_kmer_count << endl;
		cout << "   Number of k-mers having count>1       :  " << kmodel_kmer_count << endl;
	}

	void show_usage_rate() {
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
			printf("   No.%d coupled-bit array occupancy rate :  %0.3f\n", i + 1, c2 / bits_data[i].bit_array_length);
		}
		delete[] dic_num;
	}

	virtual void show_bloomfilter_info(uint64_t total_byte) {
		cout << "   (k-2)-mer BF size for count>1 k-mers   :  " << Tools::filesize_format(byte_km_back) << endl;
		total_byte += byte_km_back;
		cout << "   Total memory usage                     :  " << Tools::filesize_format(total_byte) << endl;
	}

	void show_kmodel_info() {
		uint64_t total_byte_array_size = 0;
		for (int i = 0; i < n_bits; i++) {
			total_byte_array_size += bits_data[i].bit_array_length >> 3;
		}
		total_byte_array_size *= 2;
		double left = kmodel_kmer_count; 
		uint64_t map_byte = last_map.size() * 44;
		cout << "   Memory usage for coupled-bit arrays   :  " << Tools::filesize_format(total_byte_array_size) << endl;
		cout << "   Memory usage for hash_map             :  " << Tools::filesize_format(map_byte) << endl; 
		total_byte_array_size += map_byte;
		show_bloomfilter_info(total_byte_array_size);
		show_usage_rate();
	}


protected:
	OccuBin* occu_bin;
	int n_hash;
	int n_bits;
	int last_map_size = 0; //the kmer count that can't save to bit_array and save to last_map
	uint64_t once_kmer_count = 0; //the kmer_count that the kmer only appears once
	uint64_t kmodel_kmer_count = 0; //the kmer_count stored in kmodel
	uint64_t byte_array_size; //byte size
	uint64_t bit_array_length; //the real bit table length,bit_array_length=byte_array_size*8
	unordered_map<uint64_t, uint32_t> last_map; //save the left kmer
	BitSaveData* bits_data;
	uint8_t* km_back;
	int km_back_n_hash = 5;
	uint64_t byte_km_back, length_km_back;
	KmerBuff **buff; //buff array used by thread
	uint32 bucket_size = 1 << 18;
	int* buff_n; //the real length of buffer
	int buff_size;// bucket_size*this->n_bits;
	uint64_t buff_idx = 0;
	virtual void get_candidates(string kmer, vector<int> &v_candidates) {
		string min_kmer = Tools::get_min_kmer(kmer);
		int v_from_map = get_kmer_from_map(min_kmer);
		if (v_from_map > -1) {
			v_candidates.push_back(occu_bin->occ_to_bin(v_from_map));
			return;
		}
		if (check_back_bloomfilter(min_kmer, km_back, length_km_back, km_back_n_hash)) {
			int v = find_bitarray_one(min_kmer);
			if (v > -1) v_candidates.push_back(v);
		}
	}

	vector<int> get_neighbor_kmer_bin(string kmer) {
		vector<int> v_candidates;
		string s = "ACGT", t1_kmer, t2_kmer;
		int kmer_len = kmer.length();
		t1_kmer = kmer.substr(1); //remove the first char
		t2_kmer = kmer.substr(0, kmer_len - 1);//remove the last char 
		for (int i = 0; s[i]; i++) { //shift forward
			string nw_kmer = t1_kmer + s[i];
			get_candidates(nw_kmer, v_candidates);
		}
		for (int i = 0; s[i]; i++) { //shift back
			string nw_kmer = s[i] + t2_kmer;
			get_candidates(nw_kmer, v_candidates);
		}
		return v_candidates;
	}

	int kmer_to_bin(string kmer, bool in_bloomfilter, int occ) {
		vector<int> v_bin = find_bitarray(kmer);
		int len_v_bin = v_bin.size();
		if (len_v_bin == 0) {
			return in_bloomfilter ? 1 : 0;
		}
		if (len_v_bin == 1) {
			if (in_bloomfilter) {
				vector<int> v_candidates = get_neighbor_kmer_bin(kmer);
				int cnt_one = 0;
				for (auto v : v_candidates)
					if (v == 1) cnt_one++;
				if (cnt_one >= v_candidates.size() / 2) return 1;
			}
			return v_bin[0]; //80%
		}
		vector<int> v_candidates = get_neighbor_kmer_bin(kmer);
		int len_v_can = v_candidates.size();
		if (len_v_can <= 0) { //none candidates=FP
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

	void insert_bloomfilter(string kmer, uint8_t* bit_bf, uint64_t bf_length, int num_hash) {
		uint64_t pos;
		int len = kmer.size();
		auto str_kmer = kmer.c_str();
		for (int i = 0; i<num_hash; i++) {
			pos = Tools::murmur_hash64(str_kmer, len, HashSeeds[i]) % bf_length;
			set_bit(bit_bf, pos);
		}
	}


	bool check_back_bloomfilter(string kmer, uint8_t* bit_bf, uint64_t bf_length, int num_hash) {
		int len = kmer.length();
		string nw_kmer = kmer.substr(1, len - 2);
		return check_bloomfilter(nw_kmer, bit_bf, bf_length, num_hash);
	}


	void init_kmc_api(string file_kmer_data_base, CKMCFile &kmer_data_base, CKmerAPI &kmer_object) {
		if (!kmer_data_base.OpenForListing(file_kmer_data_base)) {
			cout << "init_occ_count()__can't open the kmer_data_base\n";
			exit(1);
		}
		uint32 _kmer_length = kmer_data_base.KmerLength();
		kmer_object = CKmerAPI(_kmer_length);
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
		return a[il].occ ? il + 1 : 0; 
	}

	void insert_array(KmerBuff* buff, int index, int &buff_n) {
		for (int c = 0; c < buff_n; c++) {
			bool flag = insert_to_array(buff[c].kmer, occu_bin->occ_to_bin(buff[c].occ), index);
			if (flag) {
				string nw_kmer = buff[c].kmer.substr(1, buff[c].kmer.size() - 2);
				insert_bloomfilter(nw_kmer, this->km_back, this->length_km_back, this->km_back_n_hash);
				buff[c].occ = 0;
			}
		}
		buff_n = reorder_buffer(buff, buff_n);
	}

	void insert_with_thread(KmerBuff** buff, int* buff_real_n, int n_thread, int bucket_size) {
		int i, j;
		for (int t = 0; t < n_thread; t++) {
			#pragma omp parallel for num_threads(n_thread)
			for (i = 0; i < n_thread; ++i) {
				insert_array(buff[i], (i + t) % n_thread, buff_real_n[i]);
			}
		}
		for (i = 0; i < n_thread; ++i) {
			for (j = 0; j < buff_real_n[i]; ++j) {
				uint64_t key = Tools::kmers2uint64(buff[i][j].kmer);
				last_map.insert(make_pair(key, buff[i][j].occ));
			}
			buff_real_n[i] = bucket_size;
		}
	}

	virtual void init_buff(int n_thread) {
		buff = new KmerBuff*[n_thread];
		buff_n = new int[n_thread];
		buff_size = bucket_size*n_thread;
		for (int i = 0; i < n_thread; i++) {
			buff[i] = new KmerBuff[bucket_size];
			buff_n[i] = bucket_size;
		}
	}

	void push_to_array(string kmer, uint32 occ) {
		int row = buff_idx / bucket_size;
		int col = buff_idx%bucket_size;
		buff[row][col].kmer = kmer;
		buff[row][col].occ = occ;
		buff_idx++;
		if (buff_idx >= buff_size) {
			insert_with_thread(buff, buff_n, n_bits, bucket_size);
			buff_idx = 0;
		}
	}

	void push_last_to_array(int n_thread) {
		int row = (buff_idx - 1) / bucket_size;
		int col = (buff_idx - 1) % bucket_size;
		buff_n[row] = col + 1;
		for (int i = row + 1; i < n_thread; i++)
			buff_n[i] = 0;
		insert_with_thread(buff, buff_n, n_thread, bucket_size);

		delete[] buff;
		delete buff_n;
	}

	inline void set_bit(uint8_t *bit, uint64_t pos) {
		uint64_t row = pos >> 3;
		uint64_t col = pos & 0x7;
		uint8_t x = 0x1 << (8 - col - 1); //set the bit in the column in the table 
		__sync_fetch_and_or(bit + row, x);
	}

	inline bool check_bit(const uint8_t *bit, uint64_t pos) {
		uint64_t row = pos >> 3;
		uint64_t col = pos & 0x7;
		return (bit[row] >> (8 - col - 1)) & 0x1;
	}

	virtual void set_bit_array_size(uint64 kmer_count) {
		this->kmodel_kmer_count = kmer_count;
		this->byte_array_size = (kmer_count >> 4) * n_hash;
		this->bit_array_length = this->byte_array_size << 3;
		this->byte_km_back = (kmer_count >> 4)* this->km_back_n_hash;
		this->length_km_back = this->byte_km_back << 3;
	}

	virtual void init_bit_aray(uint64 kmer_count) {
		set_bit_array_size(kmer_count);
		bits_data = new BitSaveData[n_bits];
		km_back = new uint8_t[byte_km_back];
		memset(km_back, 0, sizeof(uint8_t)*byte_km_back);
		uint64_t t_byte_array_size = this->byte_array_size;
		for (int i = 0; i < n_bits; i++) {
			BitSaveData bit_data;
			bit_data.bit_array_length = t_byte_array_size << 3;
			bit_data.bit_array_1 = new uint8_t[t_byte_array_size];
			bit_data.bit_array_2 = new uint8_t[t_byte_array_size];
			memset(bit_data.bit_array_1, 0, sizeof(uint8_t)*t_byte_array_size);
			memset(bit_data.bit_array_2, 0, sizeof(uint8_t)*t_byte_array_size);
			bit_data.hash_seed = new uint32_t[n_hash];
			for (int j = 0; j < n_hash; j++) {
				int index = (i*n_hash + j) % 128;
				bit_data.hash_seed[j] = HashSeeds[index];
			}
			bits_data[i] = bit_data;
		}
	}


	bool insert_to_array(string kmer, uint32_t occ, int index) {
		uint64_t* bit1_v = new uint64_t[n_hash];
		uint64_t* bit2_v = new uint64_t[n_hash];
		uint32_t*  hash_seed = bits_data[index].hash_seed;
		uint8_t* bit_array_1 = bits_data[index].bit_array_1;
		uint8_t* bit_array_2 = bits_data[index].bit_array_2;
		uint64_t array_length = bits_data[index].bit_array_length;
		for (int i = 0; i < n_hash; i++) {
			bit1_v[i] = occ & 0x1;
			bit2_v[i] = Tools::murmur_hash64(kmer.c_str(), kmer.size(), hash_seed[i]) % array_length;
			occ >>= 1;
		}
		bool ok = true; 
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

class KModelOne : public KModel
{
public:
	KModelOne() {}
	KModelOne(OccuBin* occu_bin, int n_bits) : KModel(occu_bin, n_bits)
	{
		this->bf_n_hash = occu_bin->get_nhash() - 1;
		this->bf_back_n_hash = occu_bin->get_nhash() - 2;
	}

	void init_buff(int n_thread) {
		KModel::init_buff(n_thread);
		kmer_once_buff = new string[kmer_once_buff_size];
	}

	void push_to_bloom(int t_num = 4) {
#pragma omp parallel for num_threads(t_num) 
		for (int i = 0; i < once_idx; ++i)
		{
			insert_to_bloomfilter(kmer_once_buff[i]);
		}
	}

	void init_KModel(string file_kmer_data_base) {
		uint32 counter;
		int i;
		CKMCFile kmer_data_base;
		CKmerAPI kmer_object;
		init_kmc_api(file_kmer_data_base, kmer_data_base, kmer_object);
		calc_one_kmer_count(kmer_data_base, kmer_object);
		init_bit_aray(kmer_data_base.KmerCount());
		init_buff(n_bits);
		show_header_info(); 
		auto start = chrono::high_resolution_clock::now();
		while (kmer_data_base.ReadNextKmer(kmer_object, counter)) {
			if (counter == 1) {
				kmer_once_buff[once_idx++] = kmer_object.to_string();
				if (once_idx >= kmer_once_buff_size) {
					push_to_bloom();
					once_idx = 0;
				}
			}
			else
				push_to_array(kmer_object.to_string(), counter);
		}
		push_last_to_array(n_bits);
		push_to_bloom();
		delete[] kmer_once_buff;
		chrono::duration<double> dur = chrono::high_resolution_clock::now() - start;
		cout << "   kmcEx model construction time         :  " << dur.count() << endl;
		kmer_data_base.Close();
	}

	int kmer_to_occ(string kmer, uint32_t r_occ = 0) {
		kmer = Tools::get_min_kmer(kmer);
		int occ = get_kmer_from_map(kmer);
		if (occ > -1) {
			return occ;
		}
		bool b_km_back = check_back_bloomfilter(kmer, km_back, length_km_back, km_back_n_hash);
		bool b_bf01 = check_bloomfilter(kmer, bit_bf1, length_bf1, bf_n_hash);
		bool b_bf01_back = check_back_bloomfilter(kmer, bit_bf1_back, length_bf1_back, bf_back_n_hash);
		bool in_bloomfilter = b_bf01 && b_bf01_back;
		if (!b_km_back && !in_bloomfilter) {  
			return 0;
		}
		if (in_bloomfilter && !b_km_back) return 1;
		int bin = kmer_to_bin(kmer, in_bloomfilter, r_occ);
		occ = occu_bin->bin_to_mean(bin);
		return occ;
	}

	void save_model(string save_dir) {
		KModel::save_model(save_dir);
		FILE *fp_bloom = fopen((save_dir + "/bloom.bin").c_str(), "wb");
		fwrite(bit_bf1, sizeof(uint8_t), byte_bf1, fp_bloom);
		fclose(fp_bloom);
		FILE *fp_bloom2 = fopen((save_dir + "/bloom2.bin").c_str(), "wb");
		fwrite(bit_bf1_back, sizeof(uint8_t), byte_bf1_back, fp_bloom2);
		fclose(fp_bloom2);
	}

	void load_model(string save_dir) {
		KModel::load_model(save_dir);
		FILE *fp_bloom = fopen((save_dir + "/bloom.bin").c_str(), "rb");
		bit_bf1 = new uint8_t[byte_bf1];
		fread(bit_bf1, sizeof(uint8_t), byte_bf1, fp_bloom);
		fclose(fp_bloom);
		FILE *fp_bloom2 = fopen((save_dir + "/bloom2.bin").c_str(), "rb");
		bit_bf1_back = new uint8_t[byte_bf1_back];
		fread(bit_bf1_back, sizeof(uint8_t), byte_bf1_back, fp_bloom2);
		fclose(fp_bloom2);
	}

protected:
	uint8_t *bit_bf1; 
	uint8_t *bit_bf1_back;
	uint64_t byte_bf1, length_bf1;
	uint64_t byte_bf1_back, length_bf1_back;
	int bf_n_hash = 6;
	int bf_back_n_hash = 5;

	const uint32 kmer_once_buff_size = 1 << 19;
	string* kmer_once_buff;
	uint32 once_idx = 0;

	void init_bit_aray(uint64 kmer_count) {
		KModel::init_bit_aray(kmer_count);
		bit_bf1 = new uint8_t[byte_bf1];
		bit_bf1_back = new uint8_t[byte_bf1_back];
		memset(bit_bf1, 0, sizeof(uint8_t)*byte_bf1);
		memset(bit_bf1_back, 0, sizeof(uint8_t)*byte_bf1_back);
	}

	void set_bit_array_size(uint64 kmer_count) {
		if (!this->kmodel_kmer_count) 
			this->kmodel_kmer_count = kmer_count - this->once_kmer_count;
		this->byte_array_size = (kmodel_kmer_count >> 4) * n_hash;
		this->bit_array_length = this->byte_array_size << 3;

		this->byte_bf1 = once_kmer_count / 5.5 * bf_n_hash;
		this->length_bf1 = this->byte_bf1 << 3;
		this->byte_bf1_back = (once_kmer_count >> 3)* bf_back_n_hash;
		this->length_bf1_back = this->byte_bf1_back << 3;

		this->byte_km_back = (kmodel_kmer_count >> 4)* this->km_back_n_hash;
		this->length_km_back = this->byte_km_back << 3;
	}

	void insert_to_bloomfilter(string kmer) {
		insert_bloomfilter(kmer, this->bit_bf1, this->length_bf1, this->bf_n_hash);
		string nw_kmer = kmer.substr(1, kmer.size() - 2);
		insert_bloomfilter(nw_kmer, this->bit_bf1_back, this->length_bf1_back, this->bf_back_n_hash);
	}

	void calc_one_kmer_count(CKMCFile &kmer_data_base, CKmerAPI kmer_object) {
		uint32 counter;
		while (kmer_data_base.ReadNextKmer(kmer_object, counter)) {
			if (counter == 1)++once_kmer_count;
		}
		kmer_data_base.RestartListing();
	}

	void get_candidates(string kmer, vector<int> &v_candidates) {
		string min_kmer = Tools::get_min_kmer(kmer);
		int v_from_map = get_kmer_from_map(min_kmer);
		if (v_from_map > -1) {
			v_candidates.push_back(occu_bin->occ_to_bin(v_from_map));
			return;
		}
		bool b_bf01 = check_bloomfilter(min_kmer, bit_bf1, length_bf1, bf_n_hash);
		bool b_bf01_back = check_back_bloomfilter(min_kmer, bit_bf1_back, length_bf1_back, bf_back_n_hash);
		if (b_bf01&&b_bf01_back) {
			v_candidates.push_back(1);
			return;
		}
		if (check_back_bloomfilter(min_kmer, km_back, length_km_back, km_back_n_hash)) {
			int v = find_bitarray_one(min_kmer);
			if (v > -1) v_candidates.push_back(v);
		}
	}

	void show_bloomfilter_info(uint64_t total_byte) {
		cout << "   (k-2)-mer BF size for count>1 k-mers  :  " << Tools::filesize_format(byte_km_back) << endl;
		cout << "   k-mer BF size for count=1 k-mers      :  " << Tools::filesize_format(byte_bf1) << endl;
		cout << "   (k-2)-mer BF size for count=1 k-mers  :  " << Tools::filesize_format(byte_bf1_back) << endl;
		total_byte += byte_km_back + byte_bf1 + byte_bf1_back;
		cout << "   Total memory usage                    :  " << Tools::filesize_format(total_byte) << endl;
	}

};


KModel* get_model(int ci = 1,int cs=1023, int num_hash = 7, int num_bit = 5) {
	OccuBin* occu_bin = new OccuBin(cs + 1, num_hash);
	return ci > 1 ? new KModel(occu_bin, num_bit) : new KModelOne(occu_bin, num_bit);
}

KModel* get_model(string save_dir) {
	ifstream fin(save_dir + "/param.conf");
	if (!fin) {
		cout << "Load_model: cant't open the param.conf !\n";
		exit(1);
	}
	string t_str;
	int n_hash, n_bits, max_counter, last_map_size;
	uint64_t once_kmer_count, kmodel_kmer_count;
	fin >> t_str >> n_hash;
	fin >> t_str >> n_bits;
	fin >> t_str >> once_kmer_count;
	fin >> t_str >> kmodel_kmer_count;
	fin >> t_str >> last_map_size;
	fin >> t_str >> max_counter;
	fin.close();
	int ci = once_kmer_count > 0 ? 1 : 2;
	KModel* km = get_model(ci, max_counter - 1, n_hash, n_bits);
	km->set_load_param(once_kmer_count, kmodel_kmer_count, last_map_size);
	km->load_model(save_dir);
	return km;
}

#endif
