#!/bin/bash
# Warning: Run this script in a clean directory, otherwise existing files may be deleted!

# Preparing steps: ensure that pong and pophelper are installed

n_replicates=25

# compile Crimp
make -B -C ..
cp ../crimp .

# download CLUMPP along with arabid and chicken datasets:
wget https://rosenberglab.stanford.edu/software/CLUMPP_Linux64.1.1.2.tar.gz
tar -xzf CLUMPP_Linux64.1.1.2.tar.gz

cp CLUMPP_Linux64.1.1.2/CLUMPP .
cp CLUMPP_Linux64.1.1.2/example/arabid/arabid.popfile .
cp CLUMPP_Linux64.1.1.2/example/chicken/chicken.indfile .

# create artificial datasets 'largeK' and 'largeR' (already contained)
#./kmeans_largeK.r
#./kmeans_largeR.r

# shuffle input Q-matrices
./permutate_input.sh $n_replicates

# run benchmarks
./benchmark_arabid.sh $n_replicates
./benchmark_chicken.sh $n_replicates
./benchmark_largeK.sh $n_replicates
./benchmark_largeR.sh $n_replicates

# average results
for f in results_{arabid,chicken,largeK,largeR}.txt; do
    ./summarize.r $f
done
