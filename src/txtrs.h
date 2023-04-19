#ifndef VARTYPE_H
#define VARTYPE_H

#include <cstring>
#include <cstdlib>
#include "htslib/kstring.h"

#define C4_A 0
#define C4_C 1
#define C4_G 2
#define C4_T 3
#define C4_U 3
#define C4_N 4

#define SEQARR "ACGTN"
#define REVSEQ "TGCAN"

#define BASE_COMP(a, b) (a + b == 3)

#define X_CODO   0

#define C4_Stop  0
#define C4_Phe   1
#define C4_Leu   2
#define C4_Ser   3
#define C4_Tyr   4
#define C4_Cys   5
#define C4_Trp   6
#define C4_Pro   7
#define C4_His   8
#define C4_Gln   9
#define C4_Arg  10
#define C4_Ile  11
#define C4_Met  12
#define C4_Thr  13
#define C4_Asn  14
#define C4_Lys  15
#define C4_Val  16
#define C4_Ala  17
#define C4_Asp  18
#define C4_Glu  19
#define C4_Gly  20
#define C4_Sec  21 // Selenocysteine
#define C4_Ukn  22 // Unknown (condon contains 'N')
#define ALL_AA_CNT 23

static const char *codon_3letter_names[] = {"*",   "Phe", "Leu", "Ser", "Tyr", "Cys", // 0-5
                                     "Trp", "Pro", "His", "Gln", "Arg", "Ile", // 6-11
                                     "Met", "Thr", "Asn", "Lys", "Val", "Ala", // 12-17
                                     "Asp", "Glu", "Gly", "Sec", // 18-21
                                     "Xaa"};// 22 for unknown


static const char *codon_1letter_names[] = {"*", "F", "L", "S", "Y", "C", // 0-5
                                            "W", "P", "H", "Q", "R", "I", // 6-11
                                            "M", "T", "N", "K", "V", "A", // 12-17
                                            "D", "E", "G", "U", // 18-21A
                                            "X"}; // 22


static const int seq2num_table[256] = {
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
};

// none mitochondrial gene condon matrix
// UGA->termination, AUA->isoleucine
static const int codon_matrix_comm[5][5][5] = {
    { { C4_Lys, C4_Asn, C4_Lys, C4_Asn, C4_Ukn, },
      { C4_Thr, C4_Thr, C4_Thr, C4_Thr, C4_Ukn, },
      { C4_Arg, C4_Ser, C4_Arg, C4_Ser, C4_Ukn, },
      { C4_Ile, C4_Ile, C4_Met, C4_Ile, C4_Ukn, }, 
      { C4_Ukn, C4_Ukn, C4_Ukn, C4_Ukn, C4_Ukn, },
    },

    { { C4_Gln, C4_His, C4_Gln, C4_His, C4_Ukn, },
      { C4_Pro, C4_Pro, C4_Pro, C4_Pro, C4_Ukn, },
      { C4_Arg, C4_Arg, C4_Arg, C4_Arg, C4_Ukn, },
      { C4_Leu, C4_Leu, C4_Leu, C4_Leu, C4_Ukn, }, 
      { C4_Ukn, C4_Ukn, C4_Ukn, C4_Ukn, C4_Ukn, },
    },

    { { C4_Glu, C4_Asp, C4_Glu, C4_Asp, C4_Ukn, },
      { C4_Ala, C4_Ala, C4_Ala, C4_Ala, C4_Ukn, },
      { C4_Gly, C4_Gly, C4_Gly, C4_Gly, C4_Ukn, },
      { C4_Val, C4_Val, C4_Val, C4_Val, C4_Ukn, }, 
      { C4_Ukn, C4_Ukn, C4_Ukn, C4_Ukn, C4_Ukn, },
    },

    { { C4_Stop, C4_Tyr, C4_Stop, C4_Tyr, C4_Ukn,},
      { C4_Ser, C4_Ser, C4_Ser, C4_Ser, C4_Ukn, },
      { C4_Stop, C4_Cys, C4_Trp, C4_Cys, C4_Ukn, },
      { C4_Leu, C4_Phe, C4_Leu, C4_Phe, C4_Ukn, }, 
      { C4_Ukn, C4_Ukn, C4_Ukn, C4_Ukn, C4_Ukn, },
    },
    {
      { C4_Ukn, C4_Ukn, C4_Ukn, C4_Ukn, C4_Ukn, },
      { C4_Ukn, C4_Ukn, C4_Ukn, C4_Ukn, C4_Ukn, },
      { C4_Ukn, C4_Ukn, C4_Ukn, C4_Ukn, C4_Ukn, },
      { C4_Ukn, C4_Ukn, C4_Ukn, C4_Ukn, C4_Ukn, },
      { C4_Ukn, C4_Ukn, C4_Ukn, C4_Ukn, C4_Ukn, },
    },
};

// with mitochondrial gene condon matrix
// UGA->tryptophan, AUA->methionine
static const int codon_matrix_mito[5][5][5] = {
    { { C4_Lys, C4_Asn, C4_Lys, C4_Asn, C4_Ukn, },
      { C4_Thr, C4_Thr, C4_Thr, C4_Thr, C4_Ukn, },
      { C4_Arg, C4_Ser, C4_Arg, C4_Ser, C4_Ukn, },
      { C4_Met, C4_Ile, C4_Met, C4_Ile, C4_Ukn, }, 
      { C4_Ukn, C4_Ukn, C4_Ukn, C4_Ukn, C4_Ukn, },
    },

    { { C4_Gln, C4_His, C4_Gln, C4_His, C4_Ukn, },
      { C4_Pro, C4_Pro, C4_Pro, C4_Pro, C4_Ukn, },
      { C4_Arg, C4_Arg, C4_Arg, C4_Arg, C4_Ukn, },
      { C4_Leu, C4_Leu, C4_Leu, C4_Leu, C4_Ukn, }, 
      { C4_Ukn, C4_Ukn, C4_Ukn, C4_Ukn, C4_Ukn, },
    },

    { { C4_Glu, C4_Asp, C4_Glu, C4_Asp, C4_Ukn, },
      { C4_Ala, C4_Ala, C4_Ala, C4_Ala, C4_Ukn, },
      { C4_Gly, C4_Gly, C4_Gly, C4_Gly, C4_Ukn, },
      { C4_Val, C4_Val, C4_Val, C4_Val, C4_Ukn, }, 
      { C4_Ukn, C4_Ukn, C4_Ukn, C4_Ukn, C4_Ukn, },
    },

    { { C4_Stop, C4_Tyr, C4_Stop, C4_Tyr, C4_Ukn, },
      { C4_Ser, C4_Ser, C4_Ser, C4_Ser, C4_Ukn, },
      { C4_Trp, C4_Cys, C4_Trp, C4_Cys, C4_Ukn, },
      { C4_Leu, C4_Phe, C4_Leu, C4_Phe, C4_Ukn, }, 
      { C4_Ukn, C4_Ukn, C4_Ukn, C4_Ukn, C4_Ukn, },
    },

    {
      { C4_Ukn, C4_Ukn, C4_Ukn, C4_Ukn, C4_Ukn, },
      { C4_Ukn, C4_Ukn, C4_Ukn, C4_Ukn, C4_Ukn, },
      { C4_Ukn, C4_Ukn, C4_Ukn, C4_Ukn, C4_Ukn, },
      { C4_Ukn, C4_Ukn, C4_Ukn, C4_Ukn, C4_Ukn, },
      { C4_Ukn, C4_Ukn, C4_Ukn, C4_Ukn, C4_Ukn, },
    },
};

// DNA level sequence variant  type
enum var_seqs_type{
    var_type_unknow = 0,
    var_type_ref,
    var_type_snp,
    var_type_del,
    var_type_ins,
    var_type_delins,
    var_type_copy, // dup for ins
    var_type_complex,
    var_type_frameshift,
};

static const char *var_seqs_table[9] = {
    "unknown",
    "reftype",
    "snp",
    "del",
    "ins",
    "delins",
    "dup for ins",
    "complex",
    "frameshift",
};

// RNA level function type
enum var_func_type{
    var_is_unknown = 0,
    var_is_reference,
    var_is_intron,
    var_is_noncoding,
    var_is_utr5,
    var_is_utr3,
    var_is_synonymous,
    var_is_missense,
    var_is_nonsense, // stop gained
    var_is_inframe_insertion,
    var_is_inframe_deletion,
    var_is_inframe_delins,
    var_is_frameshift,
    var_is_stop_lost,
    var_is_stop_retained,
    var_is_complex,
    var_is_no_call,
    var_is_transcript_ablation, // whole exome deletion
    var_is_start_lost,
};

static const char *var_type_table[21] = {
    "Unknown",
    "Reference",
    "Intron",
    "Noncoding",
    "Utr5",
    "Utr3",
    "Synonymous",
    "Missense",
    "Nonsense",
    "InframeInsertion",
    "InframeDeletion",
    "InframeDelins",
    "Frameshift",
    "StopLost",
    "StopRetained",
    "Complex",
    "NoCall",
    "TranscriptAblation",
    "StartLost",
    NULL,
    NULL,
};

enum var_splice_type{
    var_is_not_splice = 0,
    var_is_splice_site,
    var_is_splice_donor,
    var_is_splice_acceptor,
};

static const char *splice_type_table[5] = {
    "NotSplice",
    "SpliceSite",
    "SpliceDonor",
    "SpliceAcceptor",
    NULL,
};

static inline int seq2code4(int seq){
    return seq2num_table[seq];
}

static inline int is_stop(const char *codon, int is_mito = 0){
    if(is_mito){
        return codon_matrix_mito[seq2code4(codon[0])][seq2code4(codon[1])][seq2code4(codon[2])] == C4_Stop;
    }else{
        return codon_matrix_comm[seq2code4(codon[0])][seq2code4(codon[1])][seq2code4(codon[2])] == C4_Stop;
    }
}

static inline int codon2aminoid(const char *codon, int ismito){
    if(ismito){
        return codon_matrix_mito[seq2code4(codon[0])][seq2code4(codon[1])][seq2code4(codon[2])];
    }else{
        return codon_matrix_comm[seq2code4(codon[0])][seq2code4(codon[1])][seq2code4(codon[2])];
    }
}

static inline void revc_seq1(char *seq, int l){
    for(int i = 0; i < l/2; ++i){
        char c = REVSEQ[seq2code4(seq[i])];
        seq[i] = REVSEQ[seq2code4(seq[l-1-i])];
        seq[l-1-i] = c;
    }
    if(l & 1){// odd len
        seq[l/2] = REVSEQ[seq2code4(seq[l/2])];
    }
}

static inline char *revc_seq2(const char *seq, int n){
    if(n == 0) return NULL;
    char *rev = (char*)calloc(n + 1, sizeof(char));
    for(int i = 0; i < n; ++i){
        rev[i] = REVSEQ[seq2code4(seq[n-1-i])];
    }
    rev[n] = '\0';
    return rev;
}

static inline bool is_same_seq(const char *a, const char *b, int l){
    for(int i = 0; i < l; ++i){
        if(a[i] != b[i]){
            return false;
        }
    }
    return true;
}

static inline char* nc2aa(const char* seq, int n, bool ismt){
    int c = 0;
    kstring_t ts = {0, 0, 0};
    kstring_t as = {0, 0, 0};
    while(c + 2 < n){
        ts.l = 0;
        kputsn(seq + c, 3, &ts);
        int idx = codon2aminoid(ts.s, ismt);
        ksprintf(&as, "%s", codon_1letter_names[idx]);
        c += 3;
    }
    kputc('\0', &as);
    if(ts.l) free(ts.s);
    return as.s;
}

static inline char* ncr2aa(const char* seq, int rn, bool ismt){
    char* rc = revc_seq2(seq, rn);
    char* as = nc2aa(rc, rn, ismt);
    free(rc);
    return as;
}

#endif
