/*
 * Copyright (c) 2021 Ulrich Lautenschlager - All Rights Reserved
 * 
 * Permission is hereby granted, free of charge, to any person 
 * obtaining a copy of this software and associated documentation 
 * files (the "Software"), to deal in the Software without 
 * restriction, including without limitation the rights to use, 
 * copy, modify, merge, publish, distribute, sublicense, and/or sell 
 * copies of the Software, and to permit persons to whom the 
 * Software is furnished to do so, subject to the following 
 * conditions:
 * 
 * The above copyright notice and this permission notice shall be 
 * included in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#if defined(unix) || defined(__unix__) || defined(__unix) || (defined (__APPLE__) && defined (__MACH__))
#include <getopt.h>
#else
#include "getopt.h"
#endif


#define LOG(x) ((x)>0 ? log2(x) : 0) 
#define MAX_LINE_WIDTH 40000
#define MAX_ROW_PREFIX_LENGTH 100
#define MAX_CLUSTERS 5000
#define ERR(file, line, message) fprintf(stderr, "Error in %s, line %d\n%s", file, line, message); if (errno!=0) perror(NULL); exit(EXIT_FAILURE)
#define MAX_CLUSTERS_EXHAUSTIVE 10
#define MAX_CLUSTERINGS_EXHAUSTIVE 40
#define MAX_SOLUTIONS_EXHAUSTIVE 2e9
#define UINT_FIXED_SIZE uint32_t // for more than 4294 clusterings, use uint64_t or decrease precision


// Structure for input data, intermediate results etc., enables more compact function calls.
typedef struct data {
	UINT_FIXED_SIZE *coeff_matrices;    // normalized coefficients (6 digits, represented as int)
	UINT_FIXED_SIZE *coeff_sum_matrix;
	UINT_FIXED_SIZE *diff_vector;
	int *current_permutations;   // current permutations -> rename: permutations
	int n_clusterings;
	int n_clusters;
	int n_items;
	int opt_gini;
	int opt_popfile;
	double *weights;
	double *row_sums;  // C*R row-wise sums of the non-normalized input matrices
	
	// additional information (optional)
	char *row_prefixes;
	int *pop_sizes;
	
	// only used for gini impurity
	double *square_sums;
	double corr_gini;
	
	// only used for shannon entropy
	double *plogp_sums;
	double corr_add;
	double corr_div;
	
	uint32_t prng_state[4];

} DATA;


void free_data(DATA *D) {
	free(D->current_permutations);
	free(D->coeff_matrices);
	free(D->coeff_sum_matrix);
	free(D->diff_vector);
	free(D->plogp_sums);
	free(D->square_sums);
	free(D->weights);
	free(D->row_sums);
	free(D->row_prefixes);
	free(D->pop_sizes);
}


// This function is based on Algorithm T from Knuth 7.2.1.2 (4A, p. 323)
void compute_swap_sequence(int n_clusters, int *swap_list) {
	int d, i, j, k, m, n_swaps;
	
	n_swaps = 1;
	for (i=2; i<=n_clusters; i++) {
		n_swaps *= i;
	}
	
	d = n_swaps/2;
	swap_list[d-1] = 0;
	for (m=3; m<=n_clusters; m++) {
		k = 0;
		d = d/m;
		while (k < n_swaps) {
			k += d;
			for (j=m-2; j>=0; j--) {
				swap_list[k-1] = j;
				k += d;
			}
			swap_list[k-1] += 1;
			k += d;
			for (j=0; j<=m-2; j++) {
				swap_list[k-1] = j;
				k+= d;
			}
		}
	}
	swap_list[n_swaps-1] = 0; // restores original order

	return;
}


// foreign objective functions – naive implemetation, don't use for optimization

double pairwise_matrix_similarity_G(DATA *D, int run1, int run2) {
	
	UINT_FIXED_SIZE coeff1, coeff2;
	double diff, sqrd_diff_sum = 0, denom1 = 0, denom2 = 0;
	int j,k,j1,j2;
	
	for (j=0; j<D->n_clusters; j++) {
		j1 = D->current_permutations[run1 * D->n_clusters + j];
		j2 = D->current_permutations[run2 * D->n_clusters + j];
		
		for (k=0; k<D->n_items; k++) {
			coeff1 = D->coeff_matrices[run1 * D->n_clusters * D->n_items + j1 * D->n_items + k];
			coeff2 = D->coeff_matrices[run2 * D->n_clusters * D->n_items + j2 * D->n_items + k];
			diff = 0.000001 * ((double) coeff1 - (double) coeff2);
			sqrd_diff_sum += diff*diff * D->weights[k] * D->n_items;
			denom1 += pow(0.000001 * (double) coeff1 - 1 /(double) D->n_clusters, 2);
			denom2 += pow(0.000001 * (double) coeff2 - 1 /(double) D->n_clusters, 2);
		}
	}

	return 1 - sqrt(sqrd_diff_sum)/sqrt(sqrt(denom1)*sqrt(denom2));
}

double pairwise_matrix_similarity_G_prime(DATA *D, int run1, int run2) {
	
	UINT_FIXED_SIZE coeff1, coeff2;
	double diff, sqrd_diff_sum = 0;
	int j,k,j1,j2;
	
	for (j=0; j<D->n_clusters; j++) {
		j1 = D->current_permutations[run1 * D->n_clusters + j];
		j2 = D->current_permutations[run2 * D->n_clusters + j];
		
		for (k=0; k<D->n_items; k++) {
			coeff1 = D->coeff_matrices[run1 * D->n_clusters * D->n_items + j1 * D->n_items + k];
			coeff2 = D->coeff_matrices[run2 * D->n_clusters * D->n_items + j2 * D->n_items + k];
			diff = 0.000001 * ((double) coeff1 - (double) coeff2);
			sqrd_diff_sum += diff*diff  * D->weights[k] * D->n_items;
		}
	}

	return 1 - sqrt(sqrd_diff_sum)/sqrt(2 * (double) D->n_items);
}


double avg_similarity_H_prime(DATA *D) {
	
	int i1, i2;
	double sum_pairw_sim = 0;

	for (i1=0; i1<D->n_clusterings-1; i1++) {
		for (i2=i1+1; i2<D->n_clusterings; i2++) {
			sum_pairw_sim += pairwise_matrix_similarity_G_prime(D, i1, i2);
		}
	}

	return 2*sum_pairw_sim/(D->n_clusterings * (D->n_clusterings-1));
}


double avg_similarity_H(DATA *D) {
	
	int i1, i2;
	double sum_pairw_sim = 0;

	for (i1=0; i1<D->n_clusterings-1; i1++) {
		for (i2=i1+1; i2<D->n_clusterings; i2++) {
			sum_pairw_sim += pairwise_matrix_similarity_G(D, i1, i2);
		}
	}

	return 2*sum_pairw_sim/(D->n_clusterings * (D->n_clusterings-1));
}


double mean_kullback_leibler(DATA *D) {
	
	UINT_FIXED_SIZE coeff, coeff_sum;
	double kld = 0;
	int i,j,j1,k;
	
	for (i=0; i < D->n_clusterings; i++) {
		for (j=0; j < D->n_clusters; j++) {
			j1 = D->current_permutations[i * D->n_clusters + j];
			
			for (k=0; k<D->n_items; k++) {
				coeff = D->coeff_matrices[i * D->n_clusters * D->n_items + j1 * D->n_items + k];
				coeff_sum = D->coeff_sum_matrix[j * D->n_items + k];
				kld += (double) coeff * LOG((double) coeff * D->n_clusterings / ((double) coeff_sum)) * D->weights[k];
			}
		}
	}

	return 1e-6 * kld / D->n_clusterings;
}


static inline uint32_t rotl(uint32_t x, int k) {
	return (x << k) | (x >> (32 - k));
}


// xoshiro128** 1.1
uint32_t rand32(uint32_t *state) {

	uint32_t result = rotl(state[1] * 5, 7) * 9;
	uint32_t t = state[1] << 9;

	state[2] ^= state[0];
	state[3] ^= state[1];
	state[1] ^= state[2];
	state[0] ^= state[3];

	state[2] ^= t;

	state[3] = rotl(state[3], 11);

	return result;
}

// return integer from 0,1,...,n-1
static inline uint64_t runif(uint64_t n, uint32_t *prng_state) {
	uint64_t r;
	do {
		r = rand32(prng_state);
	} while (r >= UINT32_MAX - UINT32_MAX % n);
	return r % n;
}

void shuffle(int *arr, int size, uint32_t *prng_state) {
	int i,j,k;
	
	for (i=size-1; i>=0; i--) {
		j = runif(i+1, prng_state);
		k = arr[i];
		arr[i] = arr[j];
		arr[j] = k;
	}
	return;
}

void shuffle2(int *arr1, int *arr2, int size, uint32_t *prng_state) {
	int i,j,tmp1,tmp2;
	
	for (i=size-1; i>=0; i--) {
		j = runif(i+1, prng_state);
		tmp1 = arr1[i];
		arr1[i] = arr1[j];
		arr1[j] = tmp1;

		tmp2 = arr2[i];
		arr2[i] = arr2[j];
		arr2[j] = tmp2;
	}
	return;
}

void random_perm(DATA *D) {
	int i,j;
	for (i=0; i < D->n_clusterings; i++) {
		for (j=0; j < D->n_clusters; j++) {
			D->current_permutations[i * D->n_clusters + j] = j;
		}
		shuffle(D->current_permutations + i * D->n_clusters, D->n_clusters, D->prng_state);
	}
	return;
}


double apply_permutations(DATA *D) {
	int i,j,j_perm,k;
	double plogp_total, square_total, plogp_sum, square_sum, score;
	UINT_FIXED_SIZE coeff_sum;
	
	// for the sake of readability:
	int n_items = D->n_items;
	int n_clusters = D->n_clusters;
	
	// calculate coefficient sums 
	for (j=0; j<n_clusters; j++) {
		for (k=0; k<n_items; k++) {
			coeff_sum = 0;
			for (i=0; i<D->n_clusterings; i++) {
				j_perm = D->current_permutations[i * n_clusters + j];
				coeff_sum += D->coeff_matrices[i * n_clusters * n_items + j_perm * n_items + k];
			}
			D->coeff_sum_matrix[j * n_items + k] = coeff_sum;
		}
	}
	
	if (D->opt_gini) {
		// calculate sums of squares
		for (j=0; j<n_clusters; j++) {
			square_sum = 0;
			for (k=0; k<n_items; k++) {
				coeff_sum = D->coeff_sum_matrix[j * n_items + k];
				square_sum += (double) coeff_sum * (double) coeff_sum * D->weights[k];
			}
			D->square_sums[j] = square_sum;
		}
		
		square_total = 0;
		for (j=0; j<n_clusters; j++) {
			square_total += D->square_sums[j];
		}
		score = 1 - D->corr_gini * square_total;
	}
	else {
		// calculate summed p*log(p)) values
		for (j=0; j<n_clusters; j++) {
			plogp_sum = 0;
			for (k=0; k<n_items; k++) {
				coeff_sum = D->coeff_sum_matrix[j*n_items + k];
				plogp_sum += coeff_sum * LOG((double) coeff_sum) * D->weights[k];
			}
			D->plogp_sums[j] = plogp_sum;
		}
		
		plogp_total = 0;
		for (j=0; j<n_clusters; j++) {
			plogp_total += D->plogp_sums[j];
		}
		score = -plogp_total / D->corr_div + D->corr_add;
	}

	return score;
}

double evaluate_move(DATA *D, int i_run, int cl1, int cl2) {
	int j,k;
	double plogp_total, square_total, score;
	UINT_FIXED_SIZE coeff_sum, *coeff_sum_ptr, *coeff_ptr1, *coeff_ptr2;
	// for the sake of readability:
	int n_items = D->n_items;
	int n_clusters = D->n_clusters;
	
	// map permuted cluster indices to original ones
	int j1 = D->current_permutations[i_run*n_clusters + cl1];
	int j2 = D->current_permutations[i_run*n_clusters + cl2];
	
	coeff_ptr1 = D->coeff_matrices + i_run*n_clusters*n_items + j1*n_items;
	coeff_ptr2 = D->coeff_matrices + i_run*n_clusters*n_items + j2*n_items;
	for (k=0; k<n_items; k++) {
		D->diff_vector[k] = coeff_ptr2[k] - coeff_ptr1[k];
	}
	
	if (D->opt_gini) {
		square_total = 0;
		for (j=0; j<n_clusters; j++) {
			if (j==cl1) {
				coeff_sum_ptr = D->coeff_sum_matrix + j*n_items;
				for (k=0; k<n_items; k++) {
					coeff_sum = coeff_sum_ptr[k] + D->diff_vector[k];
					square_total += (double) coeff_sum * (double) coeff_sum * D->weights[k];
				}
			}
			else if (j==cl2) {
				coeff_sum_ptr = D->coeff_sum_matrix + j*n_items;
				for (k=0; k<n_items; k++) {
					coeff_sum = coeff_sum_ptr[k] - D->diff_vector[k];
					square_total += (double) coeff_sum * (double) coeff_sum * D->weights[k];
				}
			}
			else {
				square_total += D->square_sums[j];
			}
		}
		score = 1 - D->corr_gini * square_total;
	}
	else {
		plogp_total = 0;
		for (j=0; j<n_clusters; j++) {
			if (j==cl1) {
				for (k=0; k<n_items; k++) {
					coeff_sum = D->coeff_sum_matrix[j*n_items + k] + D->diff_vector[k];
					plogp_total += coeff_sum * LOG(coeff_sum) * D->weights[k];
				}
			}
			else if (j==cl2) {
				for (k=0; k<n_items; k++) {
					coeff_sum = D->coeff_sum_matrix[j*n_items + k] - D->diff_vector[k];
					plogp_total += coeff_sum * LOG(coeff_sum) * D->weights[k];
				}
			}
			else {
				plogp_total += D->plogp_sums[j];
			}
		}
		score = -plogp_total / D->corr_div + D->corr_add;
	}
	return score;
}
	
void update_arrays(DATA *D, int i_run, int cl1, int cl2) {
	
	int k,l,start1,start2;
	UINT_FIXED_SIZE coeff_sum;
	// for the sake of readability:
	int n_items = D->n_items;
	int n_clusters = D->n_clusters;
	
	// map permuted cluster indices to original ones
	int j1 = D->current_permutations[i_run*n_clusters + cl1];
	int j2 = D->current_permutations[i_run*n_clusters + cl2];
	
	start1 = i_run*n_clusters*n_items + j1*n_items;
	start2 = i_run*n_clusters*n_items + j2*n_items;
	for (k=0; k<n_items; k++) {
		D->diff_vector[k] = D->coeff_matrices[start2 + k] - D->coeff_matrices[start1 + k];
	}
	
	// swap indices
	l = D->current_permutations[i_run*n_clusters + cl1];
	D->current_permutations[i_run*n_clusters + cl1] = D->current_permutations[i_run*n_clusters + cl2];
	D->current_permutations[i_run*n_clusters + cl2] = l;
	
	if (D->opt_gini) {
		D->square_sums[cl1] = 0;
		D->square_sums[cl2] = 0;
		
		for (k=0; k<n_items; k++) {
			D->coeff_sum_matrix[cl1*n_items + k] += D->diff_vector[k];
			coeff_sum = D->coeff_sum_matrix[cl1*n_items + k];
			D->square_sums[cl1] += (double) coeff_sum * (double) coeff_sum * D->weights[k];
		}

		for (k=0; k<n_items; k++) {
			D->coeff_sum_matrix[cl2*n_items + k] -= D->diff_vector[k];
			coeff_sum = D->coeff_sum_matrix[cl2*n_items + k];
			D->square_sums[cl2] += (double) coeff_sum * (double) coeff_sum * D->weights[k];
		}
	}
	else {
		D->plogp_sums[cl1] = 0;
		D->plogp_sums[cl2] = 0;
		
		for (k=0; k<n_items; k++) {
			D->coeff_sum_matrix[cl1*n_items + k] += D->diff_vector[k];
			coeff_sum = D->coeff_sum_matrix[cl1*n_items + k];
			D->plogp_sums[cl1] += coeff_sum * LOG(coeff_sum) * D->weights[k];
		}

		for (k=0; k<n_items; k++) {
			D->coeff_sum_matrix[cl2*n_items + k] -= D->diff_vector[k];
			coeff_sum = D->coeff_sum_matrix[cl2*n_items + k];
			D->plogp_sums[cl2] += coeff_sum * LOG(coeff_sum) * D->weights[k];
		}
	}

	return;
}

// for exhaustive search:
double evaluate_and_update(DATA *D, int i_run, int cl1, int cl2) {
	
	int j,k,l,start1,start2;
	double plogp_total, square_total, score;
	UINT_FIXED_SIZE coeff_sum;
	// for the sake of readability:
	int n_items = D->n_items;
	int n_clusters = D->n_clusters;
	
	// map permuted cluster indices to original ones
	int j1 = D->current_permutations[i_run*n_clusters + cl1];
	int j2 = D->current_permutations[i_run*n_clusters + cl2];
	
	start1 = i_run*n_clusters*n_items + j1*n_items;
	start2 = i_run*n_clusters*n_items + j2*n_items;
	for (k=0; k<n_items; k++) {
		D->diff_vector[k] = D->coeff_matrices[start2 + k] - D->coeff_matrices[start1 + k];
	}
	
	// swap indices
	l = D->current_permutations[i_run*n_clusters + cl1];
	D->current_permutations[i_run*n_clusters + cl1] = D->current_permutations[i_run*n_clusters + cl2];
	D->current_permutations[i_run*n_clusters + cl2] = l;
	
	if (D->opt_gini) {
		D->square_sums[cl1] = 0;
		D->square_sums[cl2] = 0;
		
		for (k=0; k<n_items; k++) {
			D->coeff_sum_matrix[cl1*n_items + k] += D->diff_vector[k];
			coeff_sum = D->coeff_sum_matrix[cl1*n_items + k];
			D->square_sums[cl1] += (double) coeff_sum * (double) coeff_sum * D->weights[k];
		}

		for (k=0; k<n_items; k++) {
			D->coeff_sum_matrix[cl2*n_items + k] -= D->diff_vector[k];
			coeff_sum = D->coeff_sum_matrix[cl2*n_items + k];
			D->square_sums[cl2] += (double) coeff_sum * (double) coeff_sum * D->weights[k];
		}
		
		square_total = 0;
		for (j=0; j<n_clusters; j++) {
			square_total += D->square_sums[j];
		}
		score = 1 - D->corr_gini * square_total;
	}
	else {
		D->plogp_sums[cl1] = 0;
		D->plogp_sums[cl2] = 0;
		
		for (k=0; k<n_items; k++) {
			D->coeff_sum_matrix[cl1*n_items + k] += D->diff_vector[k];
			coeff_sum = D->coeff_sum_matrix[cl1*n_items + k];
			D->plogp_sums[cl1] += coeff_sum * LOG(coeff_sum) * D->weights[k];
		}

		for (k=0; k<n_items; k++) {
			D->coeff_sum_matrix[cl2*n_items + k] -= D->diff_vector[k];
			coeff_sum = D->coeff_sum_matrix[cl2*n_items + k];
			D->plogp_sums[cl2] += coeff_sum * LOG(coeff_sum) * D->weights[k];
		}
		
		plogp_total = 0;
		for (j=0; j<n_clusters; j++) {
			plogp_total += D->plogp_sums[j];
		}
		score = -plogp_total / D->corr_div + D->corr_add;
	}

	return score;
}


void print_avg(DATA *D, FILE *fp) {
	int j,k;

	for (k=0; k<D->n_items; k++) {
		// print additional info
		fprintf(fp, "%s", D->row_prefixes + k * (MAX_ROW_PREFIX_LENGTH + 2));
		
		for (j=0; j<D->n_clusters; j++) {
			// print coefficients
			fprintf(fp, "%.6lf", D->coeff_sum_matrix[j*D->n_items + k]/(1e6*D->n_clusterings));
			if (j < D->n_clusters-1) {
				fprintf(fp, " ");
			}
		}
		if (D->opt_popfile) {
			fprintf(fp, " %d\n", D->pop_sizes[k]);
		}
		else {
			fprintf(fp, "\n");
		}
	}
}

void print_perm(DATA *D, FILE *fp) {
	int i,j;
	for (i=0; i<D->n_clusterings; i++) {
		for (j=0; j<D->n_clusters; j++) {
			fprintf(fp, "%d%c", 1 + D->current_permutations[i*D->n_clusters + j], (j==D->n_clusters-1 ? '\n' : '\t')); // 1-based indexing
		}
	}
}

void print_coeff_ordered(DATA *D, FILE *fp, int normalized) {
	int i,j,k,l;

	for (i=0; i<D->n_clusterings; i++) {
		for (k=0; k<D->n_items; k++) {
			// print additional info
			fprintf(fp, "%s", D->row_prefixes + k * (MAX_ROW_PREFIX_LENGTH + 2));
			
			// print coefficients
			for (j=0; j<D->n_clusters; j++) {
				l = D->current_permutations[i*D->n_clusters + j];
				if (normalized == 1) {
					fprintf(fp, "%.6lf", D->coeff_matrices[i*D->n_clusters*D->n_items + l*D->n_items + k] /1e6);
				}
				else{
					fprintf(fp, "%.6lf", D->coeff_matrices[i*D->n_clusters*D->n_items + l*D->n_items + k] * D->row_sums[i*D->n_items + k] /1e6);
				}
				if (j < D->n_clusters-1) {
					fprintf(fp, " ");
				}
			}
			if (D->opt_popfile) {
				fprintf(fp, " %d\n", D->pop_sizes[k]);
			}
			else {
				fprintf(fp, "\n");
			}
		}
		fprintf(fp, "\n");
	}
}



// Initialize PRNG, read input into *Data
void initialize_from_files(DATA *Data, char *infile_coeffs, char *infile_weights, int opt_popfile, int opt_popfile_weights, int opt_gini, int opt_seed) {
	int i,j,k,p;
	char line[MAX_LINE_WIDTH];
	double raw_values[MAX_CLUSTERS];
	FILE *infile;
	char *ptr, *coeff_string;
	int n_items, n_clusters, n_clusterings, non_empty_lines, max_clusterings;
	size_t n_bytes_left;

 	double row_sum, sum;
	char *ptr_prefixes;
	char tmp_str[MAX_ROW_PREFIX_LENGTH + 2];  // prefix + trailing space + '\0'
	char tmp_str_reformatted[MAX_ROW_PREFIX_LENGTH + 2];
	
	Data->prng_state[0] = 238258174;
	Data->prng_state[1] = 2790760861;
	Data->prng_state[2] = 3110701451;
	Data->prng_state[3] = 1280583839;

// 	printf("sizeof(UINT_FIXED_SIZE): %d\n", (int) sizeof(UINT_FIXED_SIZE));
	
	// setvbuf(stdout, NULL, _IONBF, 0);
	if ((NULL != infile_weights) && opt_popfile_weights) {
		ERR(__FILE__, __LINE__, "incompatible options\n");
	}

	if (opt_seed) {
		Data->prng_state[1] += opt_seed;
	}
	else {
		Data->prng_state[1] = time(NULL) + clock();
	}
	
	// get input info:

	infile = fopen(infile_coeffs, "r");
	if (NULL == infile) {
		ERR(__FILE__, __LINE__, "can't open input file\n");
	}

	// get number of clusters from first line
	if (NULL == fgets(line, MAX_LINE_WIDTH, infile)) {
		ERR(__FILE__, __LINE__, "can't read input\n");
	}
	rewind(infile);
	
	coeff_string = strrchr(line, ':');
	if (NULL == coeff_string) {
		coeff_string = line;
	}
	else {
		// start parsing after last ":"
		coeff_string += 1;
	}

	n_clusters = 0;
	ptr = strtok(coeff_string, " \t\n\r");
	while (NULL != ptr) {
		n_clusters++;
		ptr = strtok(NULL, " \t\n\r");
	}
	
	if (opt_popfile) {
		n_clusters -= 1;
	}
	
 	if (n_clusters > MAX_CLUSTERS) {
		fprintf(stderr, "Error: too many clusters (max. %d)\n", (int) MAX_CLUSTERS);
		exit(EXIT_FAILURE);
	}
	

	non_empty_lines = 0;
	n_items = 0;
	i = 0;

	while (NULL != fgets(line, MAX_LINE_WIDTH, infile)) {
		if (strspn(line, " \t\r\n") != strlen(line)) {
			non_empty_lines++;
			i++;
		}
		else {
			if (i!=0 && n_items!=0 && n_items!=i) {
				fprintf(stderr, "Error: varying number of rows\n");
				exit(EXIT_FAILURE);
			}
			if (i!=0) {
				n_items = i;
			}
			i = 0;
		}
	}
	if (!feof(infile)) {
		ERR(__FILE__, __LINE__, "can't read input\n");
	}
	rewind(infile);

	n_clusterings = non_empty_lines/n_items;
	printf("%d clusters, %d rows, %d clusterings, %d non-empty lines\n", n_clusters, n_items, non_empty_lines/n_items, non_empty_lines);
	if (n_clusterings < 2 || n_clusters < 2 || 0 == n_items) {
		ERR(__FILE__, __LINE__, "check input!\n");
	}

	// prevent int overflow in coeff_sum_matrix
	max_clusterings = (int) (1e-6 * (pow(2, 8 * (double) sizeof(UINT_FIXED_SIZE)) - 1));
 	if (n_clusterings > max_clusterings) {
		fprintf(stderr, "Error: too many clusterings (max. %d)\nYou may want to define UINT_FIXED_SIZE as uint64_t and recompile.\n", (int) max_clusterings);
		exit(EXIT_FAILURE);
	}

	// allocate memory
	Data->current_permutations = malloc(n_clusterings * n_clusters * sizeof(int));
	Data->coeff_matrices = malloc(n_clusterings * n_items * n_clusters * sizeof(UINT_FIXED_SIZE));
	Data->coeff_sum_matrix = malloc(n_items * n_clusters * sizeof(UINT_FIXED_SIZE));
	Data->diff_vector = malloc(n_items * sizeof(UINT_FIXED_SIZE));
	Data->plogp_sums = malloc(n_clusters * sizeof(double));
	Data->square_sums = malloc(n_clusters * sizeof(double));
	Data->weights = malloc(n_items * sizeof(double));
	Data->row_sums = malloc(n_clusterings * n_items * sizeof(double));
	Data->row_prefixes = calloc(n_items * (MAX_ROW_PREFIX_LENGTH + 2), sizeof(char));
	Data->pop_sizes = malloc(n_items * sizeof(int));
	if (
		NULL == Data->coeff_matrices ||
		NULL == Data->current_permutations ||
		NULL == Data->coeff_sum_matrix ||
		NULL == Data->diff_vector ||
		NULL == Data->plogp_sums ||
		NULL == Data->square_sums ||
		NULL == Data->weights ||
		NULL == Data->row_sums ||
		NULL == Data->row_prefixes ||
		NULL == Data->pop_sizes
	) {
		ERR(__FILE__, __LINE__, "can't allocate memory\n");
	}
	
	Data->n_items = n_items;
	Data->n_clusters = n_clusters;
	Data->n_clusterings = n_clusterings;
	Data->corr_gini = 1e-12 / (n_clusterings * n_clusterings);
	Data->corr_div = 1e6 * n_clusterings;
	Data->corr_add = LOG(1e6 * n_clusterings);
	Data->opt_gini = opt_gini;
	Data->opt_popfile = opt_popfile;
	
	
	// read input file:

	i=0; // clustering
	j=0; // cluster
	k=0; // item

	p=1; // inside matrix?
	
	ptr_prefixes = Data->row_prefixes;

	while (NULL != fgets(line, MAX_LINE_WIDTH, infile)) {
		j=0;
		// empty line -> new clustering
		if (strspn(line, " \t\r\n") == strlen(line)) {
			if (p==1) {
				i++;
				k=0;
				p=0;
			}
		}
		else {
			// read optional info at the begin of each row (everything before ":")
			ptr_prefixes = Data->row_prefixes + k * (MAX_ROW_PREFIX_LENGTH + 2);
			
			n_bytes_left = strcspn(line, ":");
			
			if (line[n_bytes_left] == ':') {
				n_bytes_left += 1; // prefix including ":"
				
				if (n_bytes_left > MAX_ROW_PREFIX_LENGTH) {
					fprintf(stderr, "Error: at the begin of each line, max. %d characters of additional information are allowed\n", (int) MAX_ROW_PREFIX_LENGTH);
					exit(EXIT_FAILURE);
				}
				
				// Warning: the strcats below may cause buffer overflows if improperly modified!
				strncpy(tmp_str, line, n_bytes_left);
				tmp_str[n_bytes_left] = '\0';
				tmp_str_reformatted[0] = '\0';
				
				ptr = strtok(tmp_str, " \t");
				while (NULL != ptr) {
					strcat(tmp_str_reformatted, ptr);
					strcat(tmp_str_reformatted, " ");
					ptr = strtok(NULL, " \t");
				}

				if (i == 0) {
					strncpy(ptr_prefixes, tmp_str_reformatted, MAX_ROW_PREFIX_LENGTH + 2);
				}
				else {
					if (strcmp(ptr_prefixes, tmp_str_reformatted) != 0) {
						fprintf(stderr, "Error: conflicting information at the beginning of the %dth row:\n\tclustering 1: \"%s\"\n\tclustering %d: \"%s\"\n", k+1, ptr_prefixes, i+1, tmp_str_reformatted);
						exit(EXIT_FAILURE);
					}
				}
			}
			
			// read coefficients and, optionally, population sizes
			p=1;
			j=0;
			coeff_string = strrchr(line, ':');
			if (NULL == coeff_string) {
				coeff_string = line;
			}
			else {
				// start parsing after last ":"
				coeff_string += 1;
			}
			row_sum = 0;
			ptr = strtok(coeff_string, " \t\n\r");
			while (NULL != ptr) {
				if (j < n_clusters) {
					raw_values[j] = atof(ptr);
					row_sum += raw_values[j];
				}
				else if ((j == n_clusters) && opt_popfile) {
					if (i == 0) {
						Data->pop_sizes[k] = atoi(ptr);
						if (Data->pop_sizes[k] < 0) {
							ERR(__FILE__, __LINE__, "weights must be non-negative");
						}
					}
					else {
						if (Data->pop_sizes[k] != atoi(ptr)) {
							fprintf(stderr, "Error: conflicting population sizes at the end of the %dth row:\n\tclustering 1: %d\n\tclustering %d: %d\n", k+1, Data->pop_sizes[k], i+1, atoi(ptr));
							exit(EXIT_FAILURE);
						}
					}
				}
				else {
					ERR(__FILE__, __LINE__, "can't read input\n");
				}
				
				ptr = strtok(NULL, " \t\n\r");
				j++;
			}
			
			// normalize coefficients of current row
			for (j=0; j<n_clusters; j++) {
				Data->coeff_matrices[i*n_clusters*n_items + j*n_items + k] = (UINT_FIXED_SIZE) lround(raw_values[j] / row_sum * 1e6);
			}
			// store original row sum
			Data->row_sums[i*n_items + k] = row_sum;
			
			k++;
		}
	}
	if (!feof(infile)) {
		ERR(__FILE__, __LINE__, "can't read input\n");
	}

	fclose(infile);

	/*----------------------------------------------------------------*/
	
	// normalize weights from popfile
	if (opt_popfile_weights) {
		sum = 0;
		for (k=0; k<n_items; k++) {
			sum += (double) Data->pop_sizes[k];
		}
		if (sum == 0) {
			ERR(__FILE__, __LINE__, "sum of weights must be greater than zero\n");
		}
		for (k=0; k<n_items; k++) {
			Data->weights[k] = (double) Data->pop_sizes[k] / sum;
		}
	}
	
	// read weights from separate file:
	
	if (NULL != infile_weights) {
		infile = fopen(infile_weights, "r");
		if (NULL == infile) {
			ERR(__FILE__, __LINE__, "can't open input file\n");
		}
		
		k = 0;
		
		while (NULL != fgets(line, MAX_LINE_WIDTH, infile)) {
			ptr = strtok(line, " \t\n\r");
			while (NULL != ptr) {
				Data->weights[k] = atof(ptr);
				if (Data->weights[k] < 0) {
					ERR(__FILE__, __LINE__, "weights must be non-negative");
				}
				ptr = strtok(NULL, " \t\n\r");
				k++;
				if (k>n_items) {
					ERR(__FILE__, __LINE__, "number of weights doesn't equal number of rows\n");
				}
			}
		}
		if (!feof(infile)) {
			ERR(__FILE__, __LINE__, "can't read input\n");
		}
		
		if (k!=n_items) {
			ERR(__FILE__, __LINE__, "number of weights doesn't equal number of rows\n");
		}
		
		fclose(infile);
		
		sum = 0;
		for (k=0; k<n_items; k++) {
			sum += Data->weights[k];
		}
		if (sum == 0) {
			ERR(__FILE__, __LINE__, "sum of weights must be greater than zero\n");
		}
		for (k=0; k<n_items; k++) {
			Data->weights[k] = Data->weights[k] / sum;
		}
	}
	else if (!opt_popfile_weights) {
		for (k=0; k<n_items; k++) {
			Data->weights[k] = 1 / (double) n_items;
		}
	}

}


void optimize_RRHC(DATA *Data, int opt_gini, int opt_clumpp, int opt_n_runs, int opt_quiet) {
	int i,j,j1,j2,l;
	int run, iter;
	int i_best = 0;
	int j1_best = 0;
	int j2_best = 0;
	int *order_i;
	int *order_j1;
	int *order_j2;
	int n_eval = 0;
	int improved_iter = 0;

	int *best_permutations_so_far, *best_permutations_run;
	
	double current_score, best_score_so_far;
	double candidate_score, best_score_run;
	double H, H_prime;
	
	// for the sake of readability:
	int n_clusters = Data->n_clusters;
	int n_clusterings = Data->n_clusterings;
	int *current_permutations = Data->current_permutations;
	
	int logging_period=1000;
	clock_t t0, t1;
	
	best_permutations_so_far = malloc(n_clusterings * n_clusters * sizeof(int));
	best_permutations_run = malloc(n_clusterings * n_clusters * sizeof(int));
	if (
		NULL == best_permutations_run ||
		NULL == best_permutations_so_far
	) {
		ERR(__FILE__, __LINE__, "can't allocate memory\n");
	}
	
	
	for (i=0; i<n_clusterings; i++) {
		for (j=0; j<n_clusters; j++) {
			current_permutations[i*n_clusters + j] = j;
		}
	}
	memcpy(best_permutations_run, current_permutations, n_clusterings * n_clusters * sizeof(int));
	memcpy(best_permutations_so_far, current_permutations, n_clusterings * n_clusters * sizeof(int));
	
	
	// evaluate input order
	current_score = apply_permutations(Data);
	best_score_so_far = current_score;
	printf("initial mean %s = %.6lf\n", (opt_gini ? "Gini impurity" : "Shannon entropy"), current_score);
	
	if (opt_clumpp) {
		H = avg_similarity_H(Data);
		H_prime = avg_similarity_H_prime(Data);
		printf("CLUMPP scores: H = %.6lf, H\' = %.6lf\n", H, H_prime);
	}
	if (!opt_n_runs) {
		exit(EXIT_SUCCESS);
	}
	if (opt_quiet!=2) {
		printf("neighborhood size: %d\n", n_clusterings * (n_clusters*(n_clusters-1)) / 2);
		printf("------------------\n");
	}
	
	order_i = malloc(n_clusterings * sizeof(int));
	order_j1 = malloc((n_clusters*(n_clusters-1)) / 2 * sizeof(int));
	order_j2 = malloc((n_clusters*(n_clusters-1)) / 2 * sizeof(int));
	if (NULL==order_i || NULL==order_j1 || NULL==order_j2) {
		ERR(__FILE__, __LINE__, "can't allocate memory\n");
	}

	for (i=0; i<n_clusterings; i++) {
		order_i[i]=i;
	}

	l=0;
	for (j1=0; j1<n_clusters-1; j1++) {
		for (j2=j1+1; j2<n_clusters; j2++) {
			order_j1[l] = j1;
			order_j2[l] = j2;
			l++;
		}
	}

	for (run=1; run <= opt_n_runs; run++) {
		random_perm(Data);
		
		best_score_run = apply_permutations(Data);
		if (best_score_run < best_score_so_far) {
			best_score_so_far = best_score_run;
			memcpy(best_permutations_so_far, current_permutations, n_clusterings * n_clusters * sizeof(int));
		}
		
		iter=1;
		
		if (opt_quiet==0) {
			printf("\rrun %d, iteration %d:\t%.6lf     ", run, iter, best_score_run);
		}
		
		while (1) {
			improved_iter = 0;
			shuffle(order_i, n_clusterings, Data->prng_state);

			// evaluate all possible moves
			for (i=0; i<n_clusterings; i++) {
				shuffle2(order_j1, order_j2, (n_clusters*(n_clusters-1)) / 2, Data->prng_state);
				
				if (n_eval==0) {
					t0 = clock();
				}

				for (l=0; l<(n_clusters*(n_clusters-1)) / 2; l++) {
					j1 = order_j1[l];
					j2 = order_j2[l];

					n_eval++;
					candidate_score = evaluate_move(Data, order_i[i], j1, j2);
					
					if (candidate_score < best_score_run) {
						improved_iter = 1;
						i_best = order_i[i];
						j1_best = j1;
						j2_best = j2;
						best_score_run = candidate_score;
						update_arrays(Data, i_best, j1_best, j2_best);
						
						if ((opt_quiet==0) && (n_eval<100)) {
							printf("\rrun %d, iteration %d:\t%.6lf", run, iter, best_score_run);
						}
					}
					
					// measure time for first 100 iterations (note: later iterations may be faster)
					if (opt_quiet==0) {
						if (n_eval==100) {
							t1 = clock();
							if (t1-t0) {
								logging_period = (int) (0.1 * CLOCKS_PER_SEC * 100 / (t1-t0)); // -> ~10-20 updates/sec
							}
							else {
								logging_period = 1000;
							}
						}
						else if (n_eval%logging_period == 0) {
							printf("\rrun %d, iteration %d:\t%.6lf", run, iter, best_score_run);
						}
					}
					
				}
			}
			printf("\rrun %d, iteration %d:\t%.6lf", run, iter, best_score_run);
			iter++;
			if (improved_iter==0) {
				memcpy(best_permutations_run, current_permutations, n_clusterings * n_clusters * sizeof(int));

				if (best_score_run < best_score_so_far) {
					best_score_so_far = best_score_run;
					memcpy(best_permutations_so_far, current_permutations, n_clusterings * n_clusters * sizeof(int));
				}

				break;
			}
		}
		if (opt_quiet==0) {
			printf("\n");
		}
		else if (opt_quiet==1) {
			printf("run %d:\t%.6lf\n", run, best_score_run);
		}
	}
	
	n_eval += opt_n_runs; // account for initial solutions (nitpicking)

	printf("------------------\n%d evaluated solutions\n", n_eval);
	printf("best mean %s = %.6lf\n", (opt_gini ? "Gini impurity" : "Shannon entropy"), best_score_so_far);

	// restore optimal permutations
	memcpy(current_permutations, best_permutations_so_far, n_clusterings * n_clusters * sizeof(int));
	best_score_so_far = apply_permutations(Data);

	free(best_permutations_run);
	free(best_permutations_so_far);
	free(order_i);
	free(order_j1);
	free(order_j2);
}


void optimize_exhaustive(DATA *Data, int opt_gini, int opt_clumpp, int opt_quiet) {
	int i,j;
	int n_swaps, n_solutions;
	int n_eval = 0;
	int *next_swap;
	int *cycle_end;
	
	int *swap_list;
	int *best_permutations_so_far;
	
	double current_score, best_score_so_far;
	double H, H_prime;
	
	// for the sake of readability:
	int n_clusters = Data->n_clusters;
	int n_clusterings = Data->n_clusterings;
	int *current_permutations = Data->current_permutations;
	
	best_permutations_so_far = malloc(n_clusterings * n_clusters * sizeof(int));
	if (NULL == best_permutations_so_far) {
		ERR(__FILE__, __LINE__, "can't allocate memory\n");
	}
	
	for (i=0; i<n_clusterings; i++) {
		for (j=0; j<n_clusters; j++) {
			current_permutations[i*n_clusters + j] = j;
		}
	}
	memcpy(best_permutations_so_far, current_permutations, n_clusterings * n_clusters * sizeof(int));
	
	
	// evaluate input order
	current_score = apply_permutations(Data);
	best_score_so_far = current_score;
	printf("initial mean %s = %.6lf\n", (opt_gini ? "Gini impurity" : "Shannon entropy"), current_score);
	
	if (opt_clumpp) {
		H = avg_similarity_H(Data);
		H_prime = avg_similarity_H_prime(Data);
		printf("CLUMPP scores: H = %.6lf, H\' = %.6lf\n", H, H_prime);
	}
	if (opt_quiet!=2) {
		printf("------------------\n");
	}

	if (n_clusterings > MAX_CLUSTERINGS_EXHAUSTIVE) {
		fprintf(stderr, "Error: too many clusterings for exhaustive mode (max. %d)\n", (int) MAX_CLUSTERINGS_EXHAUSTIVE);
		exit(EXIT_FAILURE);
	}
	if (n_clusters > MAX_CLUSTERS_EXHAUSTIVE) {
		fprintf(stderr, "Error: too many clusters for exhaustive mode (max. %d)\n", (int) MAX_CLUSTERS_EXHAUSTIVE);
		exit(EXIT_FAILURE);
	}
	
	n_swaps = 1;
	for (i=2; i<=n_clusters; i++) {
		n_swaps *= i;
	}
	
	swap_list = malloc(n_swaps * sizeof(int));
	next_swap = malloc(MAX_CLUSTERINGS_EXHAUSTIVE  * sizeof(int));
	cycle_end = malloc(MAX_CLUSTERINGS_EXHAUSTIVE  * sizeof(int));
	
	if (NULL == swap_list ||
		NULL == next_swap ||
		NULL == cycle_end
	) {
		ERR(__FILE__, __LINE__, "can't allocate memory\n");
	}
	
	if (n_swaps > pow((double) MAX_SOLUTIONS_EXHAUSTIVE, (double) 1 / (double) (n_clusterings-1))) {
		fprintf(stderr, "Error: too many possible solutions for exhaustive mode (max. approx. %d)\n", (int) MAX_SOLUTIONS_EXHAUSTIVE);
		exit(EXIT_FAILURE);
	}
	
	compute_swap_sequence(n_clusters, swap_list);
	best_score_so_far = current_score;
	
	for (i=0; i<n_clusterings-1; i++) {
		next_swap[i] = -1;
		cycle_end[i] = n_swaps-2;
	}
	
	n_solutions = 1;
	for (i=0; i<n_clusterings-1; i++) { //
		n_solutions *= n_swaps;
	}
	n_eval = n_solutions;
	
	for (i=0; i<n_solutions; i++) {
		for (j=n_clusterings-2; j>=0; j--) {
			if (next_swap[j] == cycle_end[j]) {
				cycle_end[j] = (cycle_end[j] + n_swaps - 1) % n_swaps;
			}
			else {
				next_swap[j] = (next_swap[j] + 1) % n_swaps;
				//printf("swap clusters %d and %d from clustering %d\n", swap_list[next_swap[j]], swap_list[next_swap[j]]+1, j+1);
				if ((opt_quiet == 0) && (i%10000 == 0)) {
					printf("\revaluate solution %d of %d", i, n_solutions);
				}
				current_score = evaluate_and_update(Data, j+1, swap_list[next_swap[j]], swap_list[next_swap[j]]+1);
				
				if (current_score < best_score_so_far) {
					best_score_so_far = current_score;
					memcpy(best_permutations_so_far, current_permutations, n_clusterings * n_clusters * sizeof(int));
				}
				break;
			}
		}
	}
	printf("\n");

	printf("------------------\n%d evaluated solutions\n", n_eval);
	printf("best mean %s = %.6lf\n", (opt_gini ? "Gini impurity" : "Shannon entropy"), best_score_so_far);

	// restore optimal permutations
	memcpy(current_permutations, best_permutations_so_far, n_clusterings * n_clusters * sizeof(int));
	best_score_so_far = apply_permutations(Data);
	
	free(best_permutations_so_far);
	free(swap_list);
	free(next_swap);
	free(cycle_end);
}


int main(int argc, char **argv) {
	char c;
	DATA Data;
	
	char *opt_weights=NULL;
	int opt_exhaustive=0, opt_clumpp=0, opt_gini=1;
	int opt_n_runs=10, opt_seed=0, opt_quiet=0, opt_write_reordered=0;
	int opt_popfile=0, opt_popfile_weights=0;
	char fname[256];
	FILE *outfile;
	double H, H_prime;
	
	char help_string[] = 
		"Crimp v1.1.0\n"
		"\n"
		"Usage: ./Crimp [OPTION]... FILE\n"
		"\n"
		"Options:\n"
		"\n"
		"  -n RUNS           Number of hillclimbing runs. (default: 10)\n"
		"                    You may want to increase this value, depending on \n"
		"                    the variability of individual runs. If 0 is specified,\n"
		"                    only the input order is evaluated and no output files\n"
		"                    are written.\n"
		"                    \n"
		"  -s SEED           Initialize pseudorandom number generator for \n"
		"                    reproducible optimization (default: 0, i.e. disabled)\n"
		"                    \n"
		"  -w WEIGHT_FILE    Row-specific weights can be supplied as separate input file,\n"
		"                    values must be delimited by spaces, tabs, or newlines.\n"
		"                    The given values do not have to be normalized.\n"
		"                    \n"
		"  -e                Minimize mean Shannon entropy (default: Gini impurity).\n"
		"  \n"
		"  -c                Calculate CLUMPP scores H and H'. This is only done for\n"
		"                    the initial and final solution.\n"
		"\n"
		"  -h                Print help message and exit.\n"
		"                    \n"
		"  -p                Treat input as popfile and ignore last column (population sizes).\n"
		"\n"
		"  -P                Treat input as popfile and use population sizes as weights.\n"
		"                    (incompatible with option -w)\n"
		"\n"
		"  -q                Quiet mode 1:\n"
		"                    Only print final score of each optimization run.\n"
		"                    If the standard output is redirected into a file \n"
		"                    (e.g. \"./Crimp example.indfile > test.log\"), this option (or -Q)\n"
		"                    is strongly recommended.\n"
		"                    \n"
		"  -Q                Quiet mode 2:\n"
		"                    Do not report individual optimization runs at all.\n"
		"\n"
		"  -r                Create an output file for the aligned coefficient matrices.\n"
		"                    (By default, only files for optimized permutations and\n"
		"                    averaged coefficients are written.)\n"
		"\n"
		"  -R                Like -r, but write normalized coefficients (as internally used)\n"
		"                    where each row sums to approx. one.\n"
		"\n"
		"  -x                Exhaustive search:\n"
		"                    Evaluate all (K!)^(R-1) possible permutations.\n"
		"                    This is only possible for very small problem sizes!\n"
		"\n"
	;
	
	setvbuf(stdout, NULL, _IONBF, 0); //
	
	// parse options
	while ((c = getopt(argc, argv, "cehn:pPqQrRs:w:x")) != -1) {
		switch (c) {
			case 'c':
				opt_clumpp = 1;
				break;
			case 'e':
				opt_gini = 0;
				break;
			case 'h':
				printf("%s", help_string);
				exit(EXIT_SUCCESS);
			case 'n':
				opt_n_runs = atoi(optarg);
				break;
			case 'p':
				opt_popfile = 1;
				opt_popfile_weights = 0;
				break;
			case 'P':
				opt_popfile = 1;
				opt_popfile_weights = 1;
				break;
			case 'q':
				opt_quiet = 1;
				break;
			case 'Q':
				opt_quiet = 2;
				break;
			case 'r':
				opt_write_reordered = 1;
				break;
			case 'R':
				opt_write_reordered = 2;
				break;
			case 's':
				opt_seed = atoi(optarg);
				break;
			case 'w':
				opt_weights = optarg;
				break;
			case 'x':
				opt_exhaustive = 1;
				break;
			case '?':
				fprintf(stderr, "Unknown option character %c\n", optopt);
				exit(EXIT_FAILURE);
			default:
				abort ();
		}
	}

	// getopt_port/getopt.c is less flexible than the glibc version 
	if (argc > optind + 1) {
		fprintf(stderr, "Error: invalid option \"%s\". Optional arguments must be specified before the input file.\n", argv[optind+1]);
		exit(EXIT_FAILURE);
	}

	initialize_from_files(&Data, argv[optind], opt_weights, opt_popfile, opt_popfile_weights, opt_gini, opt_seed);
	
	if (opt_exhaustive) {
		optimize_exhaustive(&Data, opt_gini, opt_clumpp, opt_quiet);
	}
	else {
		optimize_RRHC(&Data, opt_gini, opt_clumpp, opt_n_runs, opt_quiet);
	}
	
	if (opt_clumpp) {
		H = avg_similarity_H(&Data);
		H_prime = avg_similarity_H_prime(&Data);
		printf("CLUMPP scores: H = %.6lf, H\' = %.6lf\n", H, H_prime);
	}
	
	strncpy(fname, argv[optind], 240);
	strncat(fname, ".average", 15);
	outfile = fopen(fname, "w");
	print_avg(&Data, outfile);
	fclose(outfile);

	strncpy(fname, argv[optind], 240);
	strncat(fname, ".permutations", 15);
	outfile = fopen(fname, "w");
	print_perm(&Data, outfile);
	fclose(outfile);

	if (opt_write_reordered != 0) {
		strncpy(fname, argv[optind], 240);
		strncat(fname, ".ordered", 15);
		outfile = fopen(fname, "w");
		print_coeff_ordered(&Data, outfile, (opt_write_reordered == 2) ? 1 : 0);
		fclose(outfile);
	}

	free_data(&Data);
}
