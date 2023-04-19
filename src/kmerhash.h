/* reconstructed from VSEARCH authored by Torbjorn Rognes, Frederic Mahe and Tomas Flouri
 * the original source code can be found in https://github.com/torognes/vsearch
 * I have changed some behaviors and add some other functionality for better usage with this project
 */
#ifndef KMER_HASH_H
#define KMER_HASH_H

#include "city.h"
#include "util.h"
#include <stdint.h>
#include <memory.h>
#include <stdlib.h>

// a hash bucket
struct kh_bucket_t{
    unsigned int kmer; // kmer value
    unsigned int pos; // 1-based position, 0 = empty
};

// kmer hash data structure
struct kh_handle_t{
    kh_bucket_t* hash; // buckets to store kmers of one sequence
    unsigned int hash_mask; // 2^(bucket number) - 1, used to fast accurate locate empty bucket
    int size; // bucket size(2*length_of_seq)
    int alloc; // number of buckets allocated
    int maxpos; // just equal length_of_seq
};

// initialize kmer hash
inline struct kh_handle_t* kh_init(){
    kh_handle_t* kh = (kh_handle_t*) malloc(sizeof(kh_handle_t));
    kh->maxpos = 0;
    kh->alloc = 256;
    kh->size = 0;
    kh->hash_mask = kh->alloc - 1; // 0xFF (11111111)
    kh->hash = (kh_bucket_t*) malloc(kh->alloc * sizeof(kh_bucket_t));
    return kh;
}

// free kmer hash
inline void kh_exit(kh_handle_t* kh){
    if(kh->hash) free(kh->hash);
    free(kh);
}

// insert a kmer into an empty hash bucket
inline void kh_insert_kmer(struct kh_handle_t* kh, int k, unsigned int kmer, unsigned int pos){
    /* find free bucket in hash */
    unsigned int j = CityHash64((char*)&kmer, (k+3)/4) & kh->hash_mask;
    while(kh->hash[j].pos) j = (j + 1) & kh->hash_mask;
    kh->hash[j].kmer = kmer;
    kh->hash[j].pos = pos;
}

// insert kmers of a sequence into hash
inline void kh_insert_kmers(struct kh_handle_t* kh, int k, char* seq, int len){
    int kmers = 1 << (2 * k);
    unsigned int kmer_mask = kmers - 1; // 2*k 1bits
    /* reallocate hash table if necessary */
    if(kh->alloc < 2 * len){
        while(kh->alloc < 2 * len) kh->alloc *= 2;
        kh->hash = (kh_bucket_t *)realloc(kh->hash, kh->alloc * sizeof(kh_bucket_t));
    }
    /* init bucket size to 1 */
    kh->size = 1;
    /* update mask to be at least larger than 2 * len */
    while(kh->size < 2 * len) kh->size *= 2;
    kh->hash_mask = kh->size - 1;
    /* update max position */
    kh->maxpos = len;
    /* clear previous buckets */
    memset(kh->hash, 0, kh->size * sizeof(struct kh_bucket_t));
    unsigned int bad = kmer_mask;
    unsigned int kmer = 0;
    char* s = seq;
    for(int pos = 0; pos < len; pos++){
        int c = *s++;
        bad <<= 2ULL;
        bad |= ambig_nuc_mask[c];
        bad &= kmer_mask;

        kmer <<= 2ULL;
        kmer |= nuc_to_2bit[c];
        kmer &= kmer_mask;

        if(!bad){
            /* 1-based pos of start of kmer */
            kh_insert_kmer(kh, k, kmer, pos - k + 1 + 1);
        }
    }
}

inline int kh_find_best_diagonal(struct kh_handle_t* kh, int k, char* seq, int len){
    int diag_counts[kh->maxpos];
    memset(diag_counts, 0, kh->maxpos * sizeof(int));
    int kmers = 1 << (2 * k);
    unsigned int kmer_mask = kmers - 1;
    unsigned int bad = kmer_mask;
    unsigned int kmer = 0;
    char * s = seq + len - 1;
    for(int pos = 0; pos < len; pos++){
        int c = *s--;
        bad <<= 2ULL;
        bad |= ambig_nuc_mask[c];
        bad &= kmer_mask;

        kmer <<= 2ULL;
        kmer |= nuc_to_2bit[nuc_to_cmp[c]];
        kmer &= kmer_mask;

        if(!bad){
            /* find matching buckets in hash */
            unsigned int j = CityHash64((char*)&kmer, (k+3)/4) & kh->hash_mask;
            while(kh->hash[j].pos){
                if(kh->hash[j].kmer == kmer){
                    int fpos = kh->hash[j].pos - 1;
                    int diag = fpos - (pos - k + 1);
                    if(diag >= 0) diag_counts[diag]++;
                }
                j = (j + 1) & kh->hash_mask;
            }
        }
    }
    int best_diag_count = -1;
    int best_diag = -1;
    int good_diags = 0;
    for(int d = 0; d < kh->maxpos - k + 1; d++){
        int diag_len = kh->maxpos - d;
        int minmatch = MAX(1, diag_len - k + 1 - k * MAX(diag_len / 20, 0));
        int c = diag_counts[d];
        if(c >= minmatch) good_diags++;
        if(c > best_diag_count){
            best_diag_count = c;
            best_diag = d;
        }
    }
    if(good_diags == 1) return best_diag;
    else return -1;
}

inline void kh_find_diagonals(struct kh_handle_t* kh, int k, char* seq, int len, int* diags){
    int kmers = 1 << (2 * k);
    unsigned int kmer_mask = kmers - 1;
    unsigned int bad = kmer_mask;
    unsigned int kmer = 0;
    char *s = seq + len - 1;
    for(int pos = 0; pos < len; pos++){
        int c = *s--;
        bad <<= 2ULL;
        bad |= ambig_nuc_mask[c];
        bad &= kmer_mask;
        kmer <<= 2ULL;
        kmer |= nuc_to_2bit[nuc_to_cmp[c]];
        kmer &= kmer_mask;
        if(!bad){
            /* find matching buckets in hash */
            unsigned int j = CityHash64((char*)&kmer, (k+3)/4) & kh->hash_mask;
            while(kh->hash[j].pos){
                if(kh->hash[j].kmer == kmer){
                    int fpos = kh->hash[j].pos - 1;
                    int diag = len + fpos - (pos - k + 1);
                    if(diag >= 0) diags[diag]++;
                }
                j = (j + 1) & kh->hash_mask;
            }
        }
    }
}

#endif
