#!/bin/bash
# usage ./benchmark_arabid.sh <number of replicates>

n_rep=$1

source utils.sh

create_arabid_popfile () {
    indir=$1
    for qmat in $(ls "$indir"/Q* | sort); do
        grep -v "^$" $qmat | awk '/[0-9]/ {printf "%d: %s %d\n", NR, $0, (NR==1) ? 10 : 1}'
        echo ""
    done
}

evaluate_scores_arabid() {
    infile=$1
    # calculate mean entropy
    ./crimp -e -n0 -w arabid_weights.txt $infile | awk 'NR==3 {printf "%s\t", $NF}'
    # calculate mean Gini impurity, H and H'
    ./crimp -c -n0 -w arabid_weights.txt $infile | awk 'NR==3 {printf "%s\t", $NF} NR==4 {printf "%s\t%s", $5, $8}' | tr -d ","
}

head -n 95 arabid.popfile | cut -d" " -f5 > arabid_weights.txt

echo -e "method\tsystem\tuser\telapsed\tswaps\tmax_res_kb\tShannon_ent\tGini_imp\tH\tHprime" > results_arabid.txt


# Crimp (heuristic only)
for obj in gini_imp entropy; do
    for n_iter in 1 5 20 100; do
        method="crimp_heuristic_${obj}_${n_iter}"
        
        for i_rep in $(seq 1 $n_rep); do
        
            # prepare input
            create_crimp_input input_shuffled/arabid/$i_rep > input_crimp.txt
            
            # run Crimp
            if [[ $obj == "gini_imp" ]]; then
                timeout 1800s /usr/bin/time -o time.txt -f "%S\t%U\t%e\t%W\t%M" \
                    ./crimp -Q -n $n_iter -w arabid_weights.txt input_crimp.txt \
                && echo finished >> time.txt
            else
                timeout 1800s /usr/bin/time -o time.txt -f "%S\t%U\t%e\t%W\t%M" \
                    ./crimp -e -Q -n $n_iter -w arabid_weights.txt input_crimp.txt \
                && echo finished >> time.txt
            fi
            
            # check if completed
            if grep -q finished time.txt; then
                # evaluate output
                rearrange_columns input_shuffled/arabid/$i_rep input_crimp.txt.permutations > q_matrices_ordered.txt
                scores="$(evaluate_scores_arabid q_matrices_ordered.txt)"
                echo -e "${method}\t$(head -n1 time.txt)\t$scores" >> results_arabid.txt
            else
                echo -e "${method}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA" >> results_arabid.txt
                break
            fi

            # clean up
            rm -f input_crimp* time.txt q_matrices_ordered.txt
        done
    done
done

# Crimp (exhaustive)
for obj in gini_imp entropy; do
    method="crimp_exhaustive_${obj}"
    
    for i_rep in $(seq 1 $n_rep); do
    
        # prepare input
        create_crimp_input input_shuffled/arabid/$i_rep > input_crimp.txt
        
        # run Crimp
        if [[ $obj == "gini_imp" ]]; then
            timeout 1800s /usr/bin/time -o time.txt -f "%S\t%U\t%e\t%W\t%M" \
                ./crimp -x -Q -w arabid_weights.txt input_crimp.txt \
            && echo finished >> time.txt
        else
            timeout 1800s /usr/bin/time -o time.txt -f "%S\t%U\t%e\t%W\t%M" \
                ./crimp -x -e -Q -w arabid_weights.txt input_crimp.txt \
            && echo finished >> time.txt
        fi
        
        # check if completed
        if grep -q finished time.txt; then
            # evaluate output
            rearrange_columns input_shuffled/arabid/$i_rep input_crimp.txt.permutations > q_matrices_ordered.txt
            scores="$(evaluate_scores_arabid q_matrices_ordered.txt)"
            echo -e "${method}\t$(head -n1 time.txt)\t$scores" >> results_arabid.txt
        else
            echo -e "${method}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA" >> results_arabid.txt
            break
        fi

        # clean up
        rm -f input_crimp* time.txt q_matrices_ordered.txt
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
                create_arabid_popfile input_shuffled/arabid/$i_rep > input_clumpp.popfile
                sed -e 's/<DATATYPE>/1/' -e "s/<REPEATS>/${n_iter}/" template.paramfile > input_clumpp.paramfile
                
                # run CLUMPP
                 timeout 1800s /usr/bin/time -o time.txt -f "%S\t%U\t%e\t%W\t%M" \
                    ./CLUMPP input_clumpp.paramfile \
                    -p input_clumpp.popfile \
                    -o input_clumpp.outfile \
                    -j input_clumpp.miscfile \
                    -k 3 \
                    -c 95 \
                    -r 9 \
                    -m $i_alg \
                    -w 1 \
                    -s $i_obj \
                && echo finished >> time.txt
                
                # check if completed
                if grep -q finished time.txt; then
                    # evaluate output
                    sed -n '/The list of permutations/,/The pairwise/p' input_clumpp.miscfile | grep "^[0-9]" > permutations
                    rearrange_columns input_shuffled/arabid/$i_rep permutations > q_matrices_ordered.txt
                    scores="$(evaluate_scores_arabid q_matrices_ordered.txt)"
                    echo -e "${method}\t$(head -n1 time.txt)\t$scores" >> results_arabid.txt
                else
                    echo -e "${method}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA" >> results_arabid.txt
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

# CLUMPP (exhaustive)
i_obj=1
for obj in H Hprime; do
    method="CLUMPP_FullSearch_${obj}"
    
    for i_rep in $(seq 1 $n_rep); do
    
        # prepare input
        create_arabid_popfile input_shuffled/arabid/$i_rep > input_clumpp.popfile
        sed -e 's/<DATATYPE>/1/' -e "s/<REPEATS>/1/" template.paramfile > input_clumpp.paramfile
        
        # run CLUMPP
            timeout 1800s /usr/bin/time -o time.txt -f "%S\t%U\t%e\t%W\t%M" \
            ./CLUMPP input_clumpp.paramfile \
            -p input_clumpp.popfile \
            -o input_clumpp.outfile \
            -j input_clumpp.miscfile \
            -k 3 \
            -c 95 \
            -r 9 \
            -m 1 \
            -w 1 \
            -s $i_obj \
        && echo finished >> time.txt
        
        # check if completed
        if grep -q finished time.txt; then
            # evaluate output
            sed -n '/The list of permutations/,/The pairwise/p' input_clumpp.miscfile | grep "^[0-9]" > permutations
            rearrange_columns input_shuffled/arabid/$i_rep permutations > q_matrices_ordered.txt
            scores="$(evaluate_scores_arabid q_matrices_ordered.txt)"
            echo -e "${method}\t$(head -n1 time.txt)\t$scores" >> results_arabid.txt
        else
            echo -e "${method}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA" >> results_arabid.txt
            break
        fi

        # clean up
        rm -f input_clumpp* permutations time.txt q_matrices_ordered.txt
    done
    
    (( i_obj ++ ))
done
