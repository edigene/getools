/* reconstructed from VSEARCH authored by Torbjorn Rognes, Frederic Mahe and Tomas Flouri
 * the original source code can be found in https://github.com/torognes/vsearch
 * I have changed some behaviors and add some other functionality for better usage with this project
 */
#ifndef DEREP_H_
#define DEREP_H_

#include <stdio.h>
#include <limits.h>
#include <stdint.h>
#include <stdlib.h>
#include <memory.h>
#include <unistd.h>
#include <libgen.h>
#include "htslib/kstring.h"
#include "htslib/khash.h"
#include "options.h"
#include "common.h"
#include "util.h"
#include "krec.h"
#include "city.h"

KHASH_MAP_INIT_STR(es, int32_t)

#define MAX_QUAL_DEP 1000000

// hash bucket
struct bucket_t{
    uint64_t hash; // hash of seq(combined with header hash)
    int64_t seqno_first; // number of first sequence with this seq
    int64_t seqno_last; // number of last sequence with this seq
    int64_t size; // marker, >=1:bucket filled, 0:bucket empty
    char* seq = NULL; // sequence
    char* ext = NULL; // comment sequence
    int len = 0; // sequence length
    uint32_t* qual = NULL; // quality
    uint64_t id = 0; // cluster id
    khash_t(es) *emap = NULL;

    bucket_t(){
        emap = kh_init(es);
    }

    ~bucket_t(){
        if(seq){ free(seq); seq = NULL; }
        if(ext){ free(ext); ext = NULL; }
        if(qual){ free(qual); qual = NULL; }
        if(emap){ 
            khiter_t k;
            for(k = kh_begin(emap); k != kh_end(emap); ++k){
                if(kh_exist(emap, k)){
                    free((void*)kh_key(emap, k));
                }
            }
            kh_destroy(es, emap);
            emap = NULL;
        }
    }

    static void out_bucket_head(FILE* f){
        fprintf(f, "count\tseq\tcomment\n");
    }

    void out_bucket_body(FILE* f){
        if(size) fprintf(f, "%lld\t%s\t%s\n", size, seq, ext);
    }

    void mean_qual(){
        if(qual && size){
            for(int i = 0; i < len; ++i) qual[i] /= MIN(size, MAX_QUAL_DEP);
        }
    }
};

// dedup worker
struct DeReper{
    // external opt
    Options* mOpt;
    // hash
    uint64_t alloc_clusters = 0; // clusters count
    uint64_t hashtablesize = 0; // number of buckets in hashtable
    uint64_t hash_mask = 0; // hashtable location mask
    bucket_t* hashtable = NULL; // hash table
    // seq
    int64_t alloc_seqlen = 0; // allocated sequence length
    char* seq_up = NULL; // normalized sequence(uppercased, T->U)
    char* rc_seq_up = NULL;// reverse complement of seq_up
    // stat
    uint64_t sequencecount = 0;// number of sequences added into hash table
    uint64_t nucleotidecount = 0; // number of bases added into hash table
    int64_t shortest = INT64_MAX; // min length of sequence added into hash table
    int64_t longest = 0; // max length of sequence added into hash table
    uint64_t discarded_short = 0;// too shourt reads not added into hash table
    uint64_t discarded_long = 0; // too long reads not added into hash table
    uint64_t clusters = 0; // unique sequence classes
    double duprate = .0; // dup rate
    int64_t* dupdist = NULL;
    int64_t maxidx = 0;
    int64_t minidx = 0;
    int64_t modidx = 0;
    int64_t quat_dup[4] = {0};
    // options
    int64_t opt_minseqlength = 0; // min seq length to derep
    int64_t opt_maxseqlength = INT_MAX; // max seq length to derep
    int64_t opt_strand = 1; // rev complete considered as dup if opt_strand > 1
    bool summarized = false;

    DeReper(){
        alloc_clusters = 1024;
        alloc_seqlen = 1023;
        hashtablesize = 2 * alloc_clusters;
        hash_mask = hashtablesize - 1;
        hashtable = (bucket_t*)malloc(sizeof(bucket_t)*hashtablesize);
        memset(hashtable, 0, sizeof(bucket_t)*hashtablesize);
        seq_up = (char*)malloc((alloc_seqlen+1)*sizeof(char));
        rc_seq_up = (char*)malloc((alloc_seqlen+1)*sizeof(char));
    }

    DeReper(Options* opt){
        mOpt = opt;
        alloc_clusters = 1024;
        alloc_seqlen = 1023;
        hashtablesize = 2 * alloc_clusters;
        hash_mask = hashtablesize - 1;
        hashtable = (bucket_t*)malloc(sizeof(bucket_t)*hashtablesize);
        memset(hashtable, 0, sizeof(bucket_t)*hashtablesize);
        seq_up = (char*)malloc((alloc_seqlen+1)*sizeof(char));
        rc_seq_up = (char*)malloc((alloc_seqlen+1)*sizeof(char));
    }

    ~DeReper(){
        clean();
        if(dupdist){ free(dupdist); dupdist = NULL; }
    }
        
    void clean(){
        if(hashtable){
            for(uint64_t i = 0; i < hashtablesize; ++i){
                if(hashtable[i].seq){ free(hashtable[i].seq); hashtable[i].seq = NULL; }
                if(hashtable[i].ext){ free(hashtable[i].ext); hashtable[i].ext = NULL; }
                if(hashtable[i].qual){ free(hashtable[i].qual); hashtable[i].qual = NULL; }
                if(hashtable[i].emap){ 
                    khiter_t k;
                    for(k = kh_begin(hashtable[i].emap); k != kh_end(hashtable[i].emap); ++k){
                        if(kh_exist(hashtable[i].emap, k)){
                            free((char*)kh_key(hashtable[i].emap, k));
                        }
                    }
                    kh_destroy(es, hashtable[i].emap);
                    hashtable[i].emap = NULL;
                }
            }
            free(hashtable);
            hashtable = NULL;
        }
        if(seq_up){ free(seq_up); seq_up = NULL; };
        if(rc_seq_up){ free(rc_seq_up); rc_seq_up = NULL; };
    }

    // add one seq to markdup, return true if seq length invalid or dup
    bool derep_rec(const krec1_t* r);
    // add one seq to markdup, return true if seq length invalid or dup
    bool derep_rec(char* seq, int seqlen, char* qual,  bool qrev, char* extseq, int extlen, bool appext = false);
    // convert string to uppercase and replact U by T
    void string_normalize(char* normalized, char* s, unsigned int len);
    // reverse complement nucleotide sequence
    void reverse_complement(char* rc, char* seq, int64_t len);
    // compare two seqs
    int seqcmp(char* a, char* b, int n);
    // rehash to allocate more bins
    void rehash();
    // summary
    void summary();
    // generate json format result
    void reportJSON(kstring_t* s, const char* dh, const char* dm);
    void reportHTML(kstring_t* s);
    // generate tsv format head
    static void tsvHead(kstring_t* s);
    // generate tsv format content
    void tsvBody(kstring_t* s);
};

extern int derep_main(int argc, char** argv);

#endif
