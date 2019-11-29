#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <stdint.h>
#include <chrono>
#include <map>
#include "omp.h"
using namespace std;

struct BinData {
	vector<uint8_t> suffix;
	int count;
};

struct PreData {
	vector<BinData> suffix_array;
};

template <typename UINT>
inline UINT bin2uint(const string &kmer, int start, int end) {
	UINT v = 0;
	for (int i = start; i < end; i++) {
		v <<= 2;
		switch (kmer[i]) {
		case 'A':;       break;
		case 'C': v |= 1; break;
		case 'G': v |= 2; break;
		case 'T': v |= 3;
		}
	}
	return v;
}

//this is for sort
inline bool compare_bin(const BinData &a, const BinData &b) {
	int len = a.suffix.size();
	for (int i = 0; i<len; i++) {
		if (a.suffix[i] != b.suffix[i])
			return a.suffix[i]<b.suffix[i];
	}
	return true;
}

class KRestData {
private:
	int k;
	int pre_len;
	int suf_len;
	int map_size;
	int pre_buffer_size = 0;
	uint64_t suffix_bin_count = 0; 
	int suff_group; 
	uint64_t suff_bin_size; //suff_bin_size=suffix_bin_count*suff_group

	int *hash2index;
	int *pre_buffer;

	int* lows;

	uint8_t* suffix_bin;
	int *count_bin;

	PreData* pd;

	inline int compare_suffix_bin(uint8_t* key, int mid) {
		mid *= suff_group; 
		for (int i = 0; i < suff_group; i++) {
			if (key[i] < suffix_bin[mid + i])
				return 1;
			if (key[i] > suffix_bin[mid + i])
				return -1;
		}
		return 0;
	}

	int get_prefix_len(int k) {
		for (int i = 7; i >= 3; i--) {
			if ((k - i) % 4 == 0)
				return i;
		}
	}

	BinData get_bin_data(const string &kmer, int pre_len, int count = 0) {
		int len = kmer.length(), interval = 4;
		BinData dl;
		for (int i = pre_len; i<len; i += interval) {
			dl.suffix.push_back(bin2uint<uint8_t>(kmer, i, i + interval));
		}
		dl.count = count;
		return dl;
	}

	void stat() {
		for (int i = 0; i < map_size; i++) {
			hash2index[i] = -1;
			int v_count = pd[i].suffix_array.size();
			if (v_count > 0) {
				suffix_bin_count += v_count;
				hash2index[i] = pre_buffer_size++;
			}
		}
	}

	void sort_suffix() {
#pragma omp parallel for num_threads(4)
		for (int i = 0; i<map_size; i++) {
			if (pd[i].suffix_array.size()>0) {
				sort(pd[i].suffix_array.begin(), pd[i].suffix_array.end(), compare_bin);
			}
		}
	}

	//transform the suffix array£¨BinData£©into bin (array)
	void transform() {
		suff_bin_size = suffix_bin_count*suff_group;
		suffix_bin = new uint8_t[suff_bin_size];
		count_bin = new int[suffix_bin_count];
		pre_buffer = new int[++pre_buffer_size]{ 0 };
		int index1 = 1, index2 = 0, index3 = 0;
		for (int i = 0; i<map_size; i++) {
			int v_count = pd[i].suffix_array.size();
			if (v_count>0) {
				pre_buffer[index1++] = pre_buffer[index1 - 1] + v_count;
				for (int j = 0; j<v_count; j++) {
					count_bin[index3++] = pd[i].suffix_array[j].count;
					auto suffix = pd[i].suffix_array[j].suffix;
					for (auto it = suffix.begin(); it != suffix.end(); it++) {
						suffix_bin[index2++] = *it;
					}
				}
			}
		}
	}

public:
	KRestData() {}

	KRestData(int k) {
		this->k = k;
		pre_len = get_prefix_len(k);
		map_size = 1 << (2 * pre_len);
		suf_len = k - pre_len;
		suff_group = suf_len / 4;
		pd = new PreData[map_size];
		hash2index = new int[map_size];
		//cout << "pre_len " << pre_len << endl;
	}

	void push_back(string kmer, int count) {
		int pre_index = bin2uint<uint32_t>(kmer, 0, pre_len);
		BinData db = get_bin_data(kmer, pre_len, count);
		pd[pre_index].suffix_array.push_back(db);
	}

	void build() {
		stat();
		sort_suffix();
		transform();
	}

	void from_file(string file) {
		FILE *fp = fopen(file.c_str(), "rb");
		//k
		fread(&k, sizeof(int), 1, fp);
		//pre_len
		fread(&pre_len, sizeof(int), 1, fp);
		//map size
		fread(&map_size, sizeof(int), 1, fp);
		//pre_count
		fread(&pre_buffer_size, sizeof(int), 1, fp);
		//suffix_bin_size
		fread(&suff_bin_size, sizeof(uint64_t), 1, fp);
		//count_bin_size
		fread(&suffix_bin_count, sizeof(uint64_t), 1, fp);

		//init
		suf_len = k - pre_len;
		suff_group = suf_len / 4;
		hash2index = new int[map_size];
		pre_buffer = new int[pre_buffer_size];
		suffix_bin = new uint8_t[suff_bin_size];
		count_bin = new int[suffix_bin_count];

		//map
		fread(hash2index, sizeof(int), map_size, fp);
		//pre_count_array
		fread(pre_buffer, sizeof(int), pre_buffer_size, fp);
		//suffix bin
		fread(suffix_bin, sizeof(uint8_t), suff_bin_size, fp);
		//count_bin
		fread(count_bin, sizeof(int), suffix_bin_count, fp);
		fclose(fp);
	}

	void save_file(string file) {
		////store to a disk
		FILE *fp = fopen(file.c_str(), "wb");
		//k
		fwrite(&k, sizeof(int), 1, fp);
		//pre_len
		fwrite(&pre_len, sizeof(int), 1, fp);
		//map size
		fwrite(&map_size, sizeof(int), 1, fp);
		//pre_count
		fwrite(&pre_buffer_size, sizeof(int), 1, fp);
		//suffix_bin_size
		fwrite(&suff_bin_size, sizeof(uint64_t), 1, fp);
		//count_bin_size
		fwrite(&suffix_bin_count, sizeof(uint64_t), 1, fp);
		//map
		fwrite(hash2index, sizeof(int), map_size, fp);
		//pre_count_array
		fwrite(pre_buffer, sizeof(int), pre_buffer_size, fp);
		//suffix bin
		fwrite(suffix_bin, sizeof(uint8_t), suff_bin_size, fp);
		//count_bin
		fwrite(count_bin, sizeof(int), suffix_bin_count, fp);
		fclose(fp);
	}

	int check_kmer(string kmer) {
		if (k != kmer.length()) {
			return 0;
		}
		uint8_t* key = new uint8_t[suff_group];
		for (int j = 0, i = pre_len; i<k; i += 4, j++) {
			key[j] = bin2uint<uint8_t>(kmer, i, i + 4);
		}
		int hashindex = bin2uint<uint32_t>(kmer, 0, pre_len);
		int pre_index = hash2index[hashindex];
		if (pre_index < 0) {
			return 0;
		}
		int low = pre_buffer[pre_index], mid;
		int high = pre_buffer[pre_index + 1];
		bool found = false;
		while (low <= high) {
			mid = (low + high) / 2;
			if (compare_suffix_bin(key, mid)>0)
				high = mid - 1;
			else if (compare_suffix_bin(key, mid)<0)
				low = mid + 1;
			else {
				found = true; break;
			}
		}
		delete[] key;
		return found ? count_bin[mid] : 0;
	}

	uint64_t get_rest_count() {
		return suffix_bin_count;
	}

	uint64_t get_all_byte_size() {
		return sizeof(uint8_t)*suff_bin_size + sizeof(int)*suffix_bin_count + sizeof(int)* pre_buffer_size + sizeof(int)*map_size;
	}
};



