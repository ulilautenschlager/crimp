#!/bin/bash
# usage ./benchmark_largeK.sh <number of replicates>

n_rep=$1

source utils.sh

echo -e "method\tsystem\tuser\telapsed\tswaps\tmax_res_kb\tShannon_ent\tGini_imp\tH\tHprime" > results_largeK.txt


# Crimp (heuristic only)
for obj in gini_imp entropy; do
    for n_iter in 1 5 20 100; do
        method="crimp_heuristic_${obj}_${n_iter}"
        
        for i_rep in $(seq 1 $n_rep); do
        
            # prepare input
            create_crimp_input input_shuffled/largeK/$i_rep > input_crimp.txt
            
            # run Crimp
            if [[ $obj == "gini_imp" ]]; then
                timeout 3600s /usr/bin/time -o time.txt -f "%S\t%U\t%e\t%W\t%M" \
                    ./crimp -Q -n $n_iter input_crimp.txt \
                && echo finished >> time.txt
            else
                timeout 3600s /usr/bin/time -o time.txt -f "%S\t%U\t%e\t%W\t%M" \
                    ./crimp -e -Q -n $n_iter input_crimp.txt \
                && echo finished >> time.txt
            fi
            
            # check if completed
            if grep -q finished time.txt; then
                # evaluate output
                rearrange_columns input_shuffled/largeK/$i_rep input_crimp.txt.permutations > q_matrices_ordered.txt
                scores="$(evaluate_scores q_matrices_ordered.txt)"
                echo -e "${method}\t$(head -n1 time.txt)\t$scores" >> results_largeK.txt
            else
                echo -e "${method}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA" >> results_largeK.txt
                break
            fi

            # clean up
            rm -f input_crimp* time.txt q_matrices_ordered.txt
        done
    done
done


# CLUMPP (heuristics only)
i_alg=2
for alg in Greedy LargeKGreedy; do
    i_obj=1
    for obj in H Hprime; do
        for n_iter in 1 5 20 100; do
            method="CLUMPP_${alg}_${obj}_${n_iter}"
            
            for i_rep in $(seq 1 $n_rep); do
            
                # prepare input
                create_clumpp_input input_shuffled/largeK/$i_rep > input_clumpp.indfile
                sed -e 's/<DATATYPE>/0/' -e "s/<REPEATS>/${n_iter}/" template.paramfile > input_clumpp.paramfile
                
                # run CLUMPP
                 timeout 3600s /usr/bin/time -o time.txt -f "%S\t%U\t%e\t%W\t%M" \
                    ./CLUMPP input_clumpp.paramfile \
                    -i input_clumpp.indfile \
                    -o input_clumpp.outfile \
                    -j input_clumpp.miscfile \
                    -k 100 \
                    -c 1000 \
                    -r 100 \
                    -m $i_alg \
                    -w 0 \
                    -s $i_obj \
                && echo finished >> time.txt
                
                # check if completed
                if grep -q finished time.txt; then
                    # evaluate output
                    sed -n '/The list of permutations/,/The pairwise/p' input_clumpp.miscfile | grep "^[0-9]" > permutations
                    rearrange_columns input_shuffled/largeK/$i_rep permutations > q_matrices_ordered.txt
                    scores="$(evaluate_scores q_matrices_ordered.txt)"
                    echo -e "${method}\t$(head -n1 time.txt)\t$scores" >> results_largeK.txt
                else
                    echo -e "${method}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA" >> results_largeK.txt
                    break
                fi

                # clean up
                rm -f input_clumpp* permutations time.txt q_matrices_ordered.txt
            done
        done
        
        (( i_obj ++ ))
    done
    
    (( i_alg ++ ))
done


# pophelper
for i_rep in $(seq 1 $n_rep); do
    method="pophelper"

    # prepare input
    create_pophelper_input input_shuffled/largeK/$i_rep > input_pophelper.txt
    
    # run pophelper
    timeout 3600s /usr/bin/time -o time.txt -f "%S\t%U\t%e\t%W\t%M" \
        ./pophelper.r \
    && echo finished >> time.txt
    
    # check if completed
    if grep -q finished time.txt; then
        # evaluate output
        rearrange_columns input_shuffled/largeK/$i_rep permutations.csv > q_matrices_ordered.txt
        scores="$(evaluate_scores q_matrices_ordered.txt)"
        echo -e "${method}\t$(head -n1 time.txt)\t$scores" >> results_largeK.txt
    else
        echo -e "${method}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA" >> results_largeK.txt
        break
    fi

    # clean up
    rm -f input_pophelper.txt permutations.csv time.txt q_matrices_ordered.txt
done


# pong
for dist in G jaccard; do
    method="pong_$dist"
    
    for i_rep in $(seq 1 $n_rep); do

        # prepare input
        create_pong_input input_shuffled/largeK/$i_rep > filemap_pong.txt
        
        # run pophelper
        timeout 3600s /usr/bin/time -o time.txt -f "%S\t%U\t%e\t%W\t%M" \
            pong -m filemap_pong.txt --dist_metric $dist -s 0 -o pong_output -f --disable_server -g -v -l largeK_colors_pong.txt \
        && echo finished >> time.txt
        
        # check if completed
        if grep -q finished time.txt; then
            # evaluate output
            cat pong_output/best_alignment_per_K.txt | grep "^Q[0-9]" | sed 's/[ \t*]*$//' | sort -k1,1 | sed 's/^Q[0-9]*\t//' > permutations.csv
            rearrange_columns input_shuffled/largeK/$i_rep permutations.csv > q_matrices_ordered.txt
            scores="$(evaluate_scores q_matrices_ordered.txt)"
            echo -e "${method}\t$(head -n1 time.txt)\t$scores" >> results_largeK.txt
        else
            echo -e "${method}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA" >> results_largeK.txt
            break
        fi

        # clean up
        rm -rf pong_output filemap_pong.txt permutations.csv time.txt q_matrices_ordered.txt
    done
done
