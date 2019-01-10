#pragma once
#include <stdint.h>
#include <string>
#include <vector>
#include <map>
#include<fstream>



const uint32_t HashSeeds[32] =
{
	86028157,86941489,87866957,88799033,89743943,90695933,91657987,92632703,
	93617369,94610833,95612849,96628801,97653247,98690023,99740651,100801609,
	101870663,102948731,104041079,105144551,106262993,107388563,108525269,109676867,
	110842351,112019693,113206337,114410599,115627249,116852959,118092629,119349731
};




class Tools
{
public:
	//the hash funtion
	inline static uint64_t murmur_hash64(const void * key, int len, uint32_t seed)
	{
		const uint64_t m = 0xc6a4a7935bd1e995;
		const int r = 47;
		uint64_t h = seed ^ (len * m);
		const uint64_t * data = (const uint64_t *)key;
		const uint64_t * end = data + (len / 8);
		while (data != end) {
			uint64_t k = *data++;

			k *= m;
			k ^= k >> r;
			k *= m;

			h ^= k;
			h *= m;
		}
		const unsigned char * data2 = (const unsigned char*)data;
		switch (len & 7)
		{
		case 7: h ^= uint64_t(data2[6]) << 48;
		case 6: h ^= uint64_t(data2[5]) << 40;
		case 5: h ^= uint64_t(data2[4]) << 32;
		case 4: h ^= uint64_t(data2[3]) << 24;
		case 3: h ^= uint64_t(data2[2]) << 16;
		case 2: h ^= uint64_t(data2[1]) << 8;
		case 1: h ^= uint64_t(data2[0]);
			h *= m;
		};

		h ^= h >> r;
		h *= m;
		h ^= h >> r;
		return h;
	}


	//binary convert into decimal
	static int bin_to_decimal(uint8_t *arr, int length) {
		int x = 0;
		for (int i = length - 1; i >= 0; i--) {
			x <<= 1;
			x |= arr[i];
		}
		return x;
	}

	static uint64_t kmers2uint64(string kmer) {
		uint64_t v = 0;
		for (int i = 0; kmer[i]; i++) {
			v <<= 2;
			char x = toupper(kmer[i]);
			switch (x) {
			case 'A':;       break;
			case 'C': v |= 1; break;
			case 'G': v |= 2; break;
			case 'T': v |= 3;
			}
		}
		return v;
	}

	static uint64_t get_file_size(string file) {
		ifstream in(file);
		if (!in) {
			puts("get_file_size()___can't open the file!");
			exit(1);
		}
		in.seekg(0, ios::end);
		uint64_t size = in.tellg();
		in.close();
		return size;
	}

	static string uint64_to_string(uint64_t u_kmer, int len) {
		string ACGT = "ACGT";
		string kmer(len, 's');
		int v = 0;
		for (int i = len - 1; i >= 0; i--) {
			v = u_kmer & 0x3;
			u_kmer >>= 2;
			kmer[i] = ACGT[v];
		}
		return kmer;
	}

	static string get_file_name(string path) {
		int pos = path.find_last_of('/');
		return path.substr(pos + 1);
	}

	static int str_to_int(char* str) {
		int sum = 0;
		for (int i = 0; str[i]; i++) {
			sum = sum * 10 + str[i] - '0';
		}
		return sum;
	}

	static string filesize_format(uint64_t size) {
		int m = 1024 * 1024;
		return to_string(size / m) + "MB";
	}

	static int vector_min(vector<int> v) {
		int _min = v[0], len = v.size();
		for (int i = 0; i < len; i++) {
			if (_min < v[i]) {
				_min = v[i];
			}
		}
		return _min;
	}

	static uint64_t get_complementation(uint64_t v, int len = 31) {
		uint64_t nw_v = 0, last = 0;
		for (int i = 0; i < len; i++) {
			last = (~(v & 0x3))&(0x3);//get the last two bits 
			nw_v <<= 2; //vacate two bits
			nw_v |= last; //put it into the new value
			v >>= 2; //remove the last two bits
		}
		return nw_v;
	}

	static uint64_t get_min(uint64_t v1, uint64_t v2) {
		return v1 <= v2 ? v1 : v2;
	}


	static uint64_t get_min_com_kmer_uint(string kmer) {
		int len = kmer.length();
		uint64_t u_kmer = kmers2uint64(kmer);
		uint64_t nw_u_kmer = get_complementation(u_kmer, len);
		return get_min(u_kmer, nw_u_kmer);
	}

	static string get_complementation(string kmer) {
		int len = kmer.length();
		uint64_t u_kmer = kmers2uint64(kmer);
		uint64_t nw_u_kmer = get_complementation(u_kmer, len);
		return uint64_to_string(nw_u_kmer, len);
	}

	static string get_min_kmer(string kmer) {
		int len = kmer.length();
		uint64_t u_kmer = kmers2uint64(kmer);
		uint64_t nw_u_kmer = get_complementation(u_kmer, len);
		if (u_kmer <= nw_u_kmer)
			return kmer;
		return uint64_to_string(nw_u_kmer, len);
	}
private:
};



