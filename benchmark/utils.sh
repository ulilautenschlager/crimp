#!/bin/bash

rearrange_columns () {
    indir=$1
    permutations=$2
    
    i=1
    for f in $(ls $indir/Q* | sort); do
        perm=$(sed -n "${i}p" $permutations | sed 's/^[ \t]*//;s/[ \t]$//;s/[ \t]\+/,/g')
        
        awk -v "perm=$perm" '
            BEGIN {
                split(perm, col_order, ",")
                n = length(col_order)
            }

            {
                for (i in col_order) {
                    printf "%s%c", $col_order[i], i+0<n ? " " : "\n"
                }
            }
        ' $f
        
        (( i ++ ))
    done
}

create_pophelper_input () {
    indir=$1
    for qmat in $(ls "$indir"/Q* | sort); do
        grep -v "^$" $qmat | awk '/[0-9]/ {print NR": "$0" 1"}'
        echo ""
    done
}

create_clumpp_input () {
    indir=$1
    for qmat in $(ls "$indir"/Q* | sort); do
        grep -v "^$" $qmat | awk '/[0-9]/ {printf NR" "NR" (0) 1 : "; print $0}'
        echo ""
    done
}

create_crimp_input () {
    indir=$1
    for qmat in $(ls "$indir"/Q* | sort); do
        cat $qmat
        echo ""
    done
}

create_pong_input () {
    indir=$1
    n_col=$(awk 'NR==1 {print NF; exit}' "$indir"/Q*)
    ls "$indir"/Q* | sort | sed "s/$/ $n_col/" | awk '{s=$1; sub(/.*\//,"",s); print s"\t"$2"\t"$1}'
}

evaluate_scores() {
    infile=$1
    # calculate mean entropy
    ./crimp -e -n0 $infile | awk '/^initial mean Shannon/ {printf "%s\t", $NF}'
    # calculate mean Gini impurity, H and H'
    ./crimp -c -n0 $infile | awk '/^initial mean Gini/ {printf "%s\t", $NF} /^CLUMPP scores/ {printf "%s\t%s", $5, $8}' | tr -d ","
}
