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

//const uint32_t HashSeeds[32] =
//{
//	93617369,112019693,96628801,119349731,106262993,86941489,98690023,95612849,91657987,99740651,90695933,108525269,
//	97653247,101870663,105144551,86028157,92632703,94610833,110842351,118092629,107388563,114410599,109676867,
//	89743943,115627249,100801609,113206337,104041079,88799033,102948731,116852959,87866957,
//};




class Tools
{
public:
	//the hash funtion
	static uint64_t murmur_hash64(const void * key, int len, uint32_t seed)
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

	static uint64_t kmers2uint64(std::string kmer) {
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

	static uint64_t get_file_size(std::string file) {
		std::ifstream in(file);
		if (!in) {
			puts("get_file_size()___can't open the file!");
			exit(1);
		}
		in.seekg(0, std::ios::end);
		uint64_t size = in.tellg();
		in.close();
		return size;
	}

	static std::string get_file_name(std::string path) {
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

	static std::string filesize_format(uint64_t size) {
		int m = 1024 * 1024;
		return std::to_string(size / m) + "MB";
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

	static std::string get_complementation(std::string kmer) {

		string com = "";
		std::map<char, char> c2c = { { 'A','T' },{ 'C','G' },{ 'T','A' },{ 'G','C' } };
		for (int i = kmer.length() - 1; i >= 0; i--) {
			com += c2c[kmer[i]];
		}
		return com;
	}
private:
};



