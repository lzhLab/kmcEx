KMC_API_DIR = kmc_api

CC      = g++
CFLAGS  = -O3 -m64 -fopenmp -std=c++11 

.cpp.o:
	$(CC) $(CFLAGS) -c $< -o $@
	@chmod +x "kmc_api/kmc"

kmcEx:  $(KMC_API_DIR)/mmer.o  $(KMC_API_DIR)/kmc_file.o $(KMC_API_DIR)/kmer_api.o  main.o
	$(CC) $(CFLAGS) -o $@  main.o  $(KMC_API_DIR)/mmer.o  $(KMC_API_DIR)/kmc_file.o $(KMC_API_DIR)/kmer_api.o

clean:
	-rm -f  $(KMC_API_DIR)/*.o
	-rm -f *.o
	-rm -f kmcEx
