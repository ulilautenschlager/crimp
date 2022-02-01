# Crimp

Crimp (cluster relabeling based on impurity minimization) allows to align 
clusters across multiple clusterings, based on their matrices of membership 
coefficients (Q-matrices). This is frequently done to post-process 
outputs from programs like Structure (Pritchard et al. 2000). Similar to CLUMPP 
(Jakobsson & Rosenberg 2007), Crimp's input clusterings have to comprise the 
same number of clusters (K) and clustered objects (C).

## Installation
Download Crimp, e.g. using `git clone https://github.com/ulilautenschlager/crimp`

- Linux: From directory "crimp/", run `make` to compile Crimp.
- Windows: Not yet supported. (However, binaries compiled using MinGW seem to work.)

## Usage
Since Crimp is a command-line tool, you usually need to open a Linux terminal. 
In the following, some familiarity with running command-line tools will be assumed.

Crimp expects one mandatory argument, which is the input file containing 
multiple clusterings (Q-matrices), e.g. `./crimp example.txt`. Additional
options (e.g., `-n 20`) may be placed before the input file. Use `./crimp -h` or
 `crimp.exe -h` to display a complete list and descriptions of available program 
options.


### Input format

Crimp expects a single input file comprising multiple 
whitespace-delimited Q-matrices, separated by one or more empty lines, e.g.

```raw
1  0
1  0
0  1

0  1
0  1
1  0
```

Optionally, each row of a Q-matrix may contain additional information at the 
beginning, ending with ":", e.g.

```raw
individual_1: 1  0
individual_2: 1  0
individual_3: 0  1

individual_1: 0  1
individual_2: 0  1
individual_3: 1  0
```

or like CLUMMP-indfiles:

```raw
1 1 (0) 1 :  1  0
2 2 (0) 1 :  1  0
3 3 (0) 1 :  0  1

1 1 (0) 1 :  0  1
2 2 (0) 1 :  0  1
3 3 (0) 1 :  1  0
```

If present in the input, such information will be ignored and not be part of
Crimp's output.

Note: The clustered objects (rows) have to be identically ordered for all 
Q-matrices!

### Output files

By default, two output files are written: an averaged Q-matrix 
(\<input\>.average) based on the aligned input matrices and the optimized column
permutations (\<input\>.permutations). When using `-r` or `-R`, Crimp also 
rearranges the original Q-matrices and writes them to \<input\>.ordered.

### List of options (see help page)

```raw
  -n RUNS           Number of hillclimbing runs. (default: 10)
                    You may want to increase this value, depending on 
                    the variability of individual runs. If 0 is specified,
                    only the input order is evaluated and no output files
                    are written.
                    
  -s SEED           Initialize pseudorandom number generator for 
                    reproducible optimization (default: 0, i.e. disabled)
                    
  -w WEIGHT_FILE    Row-specific weights can be supplied as separate input file,
                    values must be delimited by spaces, tabs, or newlines.
                    The given values do not have to be normalized.
                    
  -e                Minimize mean Shannon entropy (default: Gini impurity).
  
  -c                Calculate CLUMPP scores H and H'. This is only done for
                    the initial and final solution.

  -h                Print help message and exit.

  -q                Quiet mode 1:
                    Only print final score of each optimization run.
                    If the standard output is redirected into a file 
                    (e.g. "./Crimp example.indfile > test.log"), this option (or -Q)
                    is strongly recommended.
                    
  -Q                Quiet mode 2:
                    Do not report individual optimization runs at all.

  -r                Create an output file for the aligned coefficient matrices.
                    (By default, only files for optimized permutations and
                    averaged coefficients are written.)

  -R                Like -r, but write normalized coefficients (as internally used)
                    where each row sums to approx. one.

  -x                Exhaustive search:
                    Evaluate all (K!)^(R-1) possible permutations.
                    This is only possible for very small problem sizes!

```
