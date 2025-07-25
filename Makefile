# for MacOS with Homebrew
CC = cc
CTAG = -Xpreprocessor -fopenmp -Rpass=loop-vectorize -lomp -I"$(shell brew --prefix libomp)/include" -L"$(shell brew --prefix libomp)/lib" -O3
CSIMDTAG = -Xpreprocessor -fopenmp -Rpass=loop-vectorize -march=native -ffast-math -lomp -I"$(shell brew --prefix libomp)/include" -L"$(shell brew --prefix libomp)/lib" -O3

# for Wisteria Aquarius
# CC = fccpx
# CTAG = -Kfast -Kopenmp
# CSIMDTAG = -Kfast -Kopenmp -Ksimd

default: original_noopt original openmp simd simdv2

# --- original without optimization ---
original_noopt: original_noopt.o
	$(CC) -o original_noopt original_noopt.o -lm
original_noopt.o: original.c
	$(CC) -c original.c -o original_noopt.o

# --- original ---
original: original.o
	$(CC) $(CTAG) -o original original.o -lm
original.o: original.c
	$(CC) $(CTAG) -c original.c -o original.o

# --- openmp ---
openmp: openmp.o
	$(CC) $(CTAG) -o openmp openmp.o -lm
openmp.o: openmp.c
	$(CC) $(CTAG) -c openmp.c -o openmp.o

# --- simd ---
simd: simd.o
	$(CC) $(CSIMDTAG) -o simd simd.o -lm
simd.o: simd.c
	$(CC) $(CSIMDTAG) -c simd.c -o simd.o

# # --- advanced ---
# advanced: advanced.o
# 	$(CC) $(CTAG) -o advanced advanced.o -lm
# advanced.o: advanced.c
# 	$(CC) $(CTAG) -c advanced.c -o advanced.o

simdv2: simdv2.o
	$(CC) $(CSIMDTAG) -o simdv2 simdv2.o -lm
simdv2.o: simdv2.c
	$(CC) $(CSIMDTAG) -c simdv2.c -o simdv2.o

clean:
	rm original.o original openmp.o openmp original_noopt.o original_noopt simd.o simd simdv2.o simdv2
