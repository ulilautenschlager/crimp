#!/bin/bash
# usage ./permutate_input.sh <number of replicates>

n_rep=$1

# store original Q-matrices as separate files
rm -rf input_orig
mkdir -p input_orig/{arabid,chicken,largeK,largeR}

sed 's/.*: *//' chicken.indfile | split -d -a 2 -l 601 - input_orig/chicken/Q
cut -d" " -f2-4 arabid.popfile | head -n 864 | split -d -a 2 -l 96 - input_orig/arabid/Q
split -d -a 2 -l 1001 largeK.txt input_orig/largeK/Q
split -d -a 3 -l 251 largeR.txt input_orig/largeR/Q


# shuffle matrices and columns
rm -rf input_shuffled
for i in $(seq 1 $n_rep); do
    mkdir -p input_shuffled/{arabid,chicken,largeK,largeR}/$i
done

for indir in input_orig/{arabid,chicken,largeK,largeR}; do
    # get number of clusters/columns
    n_col=$(awk 'NR==1 {print NF; exit}' "$indir"/*)
    echo "${indir}/*: $n_col columns"

    for i_rep in $(seq 1 $n_rep); do
        i_mat=0
        
        # process Q-matrices in random order
        for qmat in $(ls "$indir"/* | shuf --random-source=/dev/urandom); do
            # shuffle columns
            awk -v "perm=$(seq 1 $n_col | shuf --random-source=/dev/urandom | tr '\n' ',' | sed 's/,$//')" '
                BEGIN {
                    split(perm, col_order, ",")
                    n = length(col_order)
                }

                {
                    for (i in col_order) {
                        printf "%s%c", $col_order[i], i+0<n ? " " : "\n"
                    }
                }' "$qmat" > "${indir/_orig/_shuffled}/${i_rep}/Q$i_mat"
                
                echo "shuffle $qmat  -->  ${indir/_orig/_shuffled}/${i_rep}/Q$i_mat"
                
            (( i_mat ++ ))
        done
    done
done
