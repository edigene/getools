/* reconstructed from ksw authored by Attractive Chaos
 * the original source code can be found in https://github.com/attractivechaos/klib
 * I have changed some behaviors and add some other functionality for better usage with this project
 */
#ifndef __AC_KSW4_GETOOLS_H
#define __AC_KSW4_GETOOLS_H

#include <stdint.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <string.h>
#ifdef __aarch64__
#include "sse2neon.h"
#else
#include <emmintrin.h>
#endif
#include <unistd.h>
#include <stdio.h>
#include <zlib.h>
#include "common.h"

#define MAPSTR "MIDNSHP=X"
#ifndef BAM_CIGAR_SHIFT
#define BAM_CIGAR_SHIFT 4
#endif

#ifndef SW_BW
#define SW_BW 50
#endif

#define KSW_FTYPE       uint32_t
#define KSW_FLALN       0x1     // 1st bit
#define KSW_FRALN       0x2     // 2nd bit
#define KSW_FIDLOCMP    0x4     // 3rd bit
#define KSW_FMAYVARSEQ  0x8     // 4th bit
#define KSW_FREFTYPE    0x10    // 5th bit
#define KSW_FMERGED     0x20    // 6th bit
#define KSW_FHIERRLOWQ  0x40    // 7th bit
#define KSW_FLOWFREQ    0x80    // 8th bit
#define KSW_FSNVINRNG   0x100   // 9th bit
#define KSW_FINSINRNG   0x200   // 10th bit
#define KSW_FDELINRNG   0x400   // 11th bit
#define KSW_FDIINRNG    0x800   // 12th bit
#define KSW_FREPSEQR    0x1000  // 13th bit
#define KSW_FSPANSGR    0x2000  // 14th bit
#define KSW_FPAIROLP    0x4000  // 15th bit
#define KSW_FPRIMDIMER  0x8000  // 16th bit
#define KSW_FMISPRIMINN 0x10000 // 17th bit
#define KSW_FMISPRIMOUT 0x20000 // 18th bit
#define KSW_FLOWSURPT   0x40000 // 19th bit
#define KSW_FMANYVARS   0x80000 // 20th bit
#define KSW_FVARINSEQ   0x100000// 21st bit
#define KSW_FRECALLHIT  0x200000// 22nd bit
#define KSW_FRECEXACT   0x400000// 23rd bit
#define KSW_FRECANYHIT  0x800000// 24th bit
#define KSW_SHIFTMS     0x1B    // 28-32bit
#define KSW_MASKMS      0x1F    // 5bit1
#define get_merge_status(x) ((x >> KSW_SHIFTMS) & KSW_MASKMS)
#define set_merge_status(x, m) (x &= 0x7FFFFFF, x |= (((KSW_FTYPE)m) << KSW_SHIFTMS))
#define left_realn_worked(x) (x & KSW_FLALN)
#define right_realn_worked(x) (x & KSW_FRALN)
#define overlap_indel_compatible(x) (x & KSW_FIDLOCMP)
#define is_valid_variant_seq(x) (x & KSW_FMAYVARSEQ)
#define is_reftype_seq(x) (x & KSW_FREFTYPE)

#ifdef __GNUC__
#define LIKELY(x) __builtin_expect((x),1)
#define UNLIKELY(x) __builtin_expect((x),0)
#else
#define LIKELY(x) (x)
#define UNLIKELY(x) (x)
#endif

#ifndef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#endif

#ifndef MAX
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#endif

/** Round an integer to the next closest power-2 integer.
 * @param  x  integer to be rounded (in place)
 * @discussion x will be modified.
 */
#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

#define KSW_XBYTE  0x10000
#define KSW_XSTOP  0x20000
#define KSW_XSUBO  0x40000
#define KSW_XSTART 0x80000

const unsigned char kswge_int2nt[5] = {'A', 'C', 'G', 'T', 'N'};
const unsigned char kswge_nt2int[256] = {
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

/** array index is an ASCII character value from a CIGAR,
 * element value is the corresponding integer opcode between 0 and 8
 */
const uint8_t kswge_encoded_ops[] = {
    0,         0,         0,         0,
    0,         0,         0,         0,
    0,         0,         0,         0,
    0,         0,         0,         0,
    0,         0,         0,         0,
    0,         0,         0,         0,
    0,         0,         0,         0,
    0,         0,         0,         0,
    0 /*   */, 0 /* ! */, 0 /* " */, 0 /* # */,
    0 /* $ */, 0 /* % */, 0 /* & */, 0 /* ' */,
    0 /* ( */, 0 /* ) */, 0 /* * */, 0 /* + */,
    0 /* , */, 0 /* - */, 0 /* . */, 0 /* / */,
    0 /* 0 */, 0 /* 1 */, 0 /* 2 */, 0 /* 3 */,
    0 /* 4 */, 0 /* 5 */, 0 /* 6 */, 0 /* 7 */,
    0 /* 8 */, 0 /* 9 */, 0 /* : */, 0 /* ; */,
    0 /* < */, 7 /* = */, 0 /* > */, 0 /* ? */,
    0 /* @ */, 0 /* A */, 0 /* B */, 0 /* C */,
    2 /* D */, 0 /* E */, 0 /* F */, 0 /* G */,
    5 /* H */, 1 /* I */, 0 /* J */, 0 /* K */,
    0 /* L */, 0 /* M */, 3 /* N */, 0 /* O */,
    6 /* P */, 0 /* Q */, 0 /* R */, 4 /* S */,
    0 /* T */, 0 /* U */, 0 /* V */, 0 /* W */,
    8 /* X */, 0 /* Y */, 0 /* Z */, 0 /* [ */,
    0 /* \ */, 0 /* ] */, 0 /* ^ */, 0 /* _ */,
    0 /* ` */, 0 /* a */, 0 /* b */, 0 /* c */,
    0 /* d */, 0 /* e */, 0 /* f */, 0 /* g */,
    0 /* h */, 0 /* i */, 0 /* j */, 0 /* k */,
    0 /* l */, 0 /* m */, 0 /* n */, 0 /* o */,
    0 /* p */, 0 /* q */, 0 /* r */, 0 /* s */,
    0 /* t */, 0 /* u */, 0 /* v */, 0 /* w */,
    0 /* x */, 0 /* y */, 0 /* z */, 0 /* { */,
    0 /* | */, 0 /* } */, 0 /* ~ */, 0 /*  */
};
typedef struct{
    int qlen, slen;
    uint8_t shift, mdiff, max, size;
    __m128i *qp, *H0, *H1, *E, *Hmax;
} kswq_t;

typedef struct {
	int score; // best score
	int te, qe; // target end and query end [0 based inclusive]
	int score2, te2; // second best score and ending position on the target
	int tb, qb; // target start and query start [0 based inclusive]
        int mm, nvar, nmm, nins, ndel, lsc, tsc;
        uint32_t* cigar;
        int ncigar;
        KSW_FTYPE smask;
} kswr_t;

const kswr_t g_defr = { 0, -1, -1, -1, -1, -1, -1 , 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

struct kswge_opt_t{
    int match = 5;
    int mismatch = 4;
    int gapo = 25;
    int gape = 0;
    char* ref = NULL;
    char* qry = NULL;
    char* qry2 = NULL;
    int bpos = -1;
    int epos = -1;
    int posb = -1;
    int pose = -1;
    int wlen = 60;
};

/**
 * Aligning two sequences
 *
 * @param qlen    length of the query sequence (typically <tlen)
 * @param query   query sequence with 0 <= query[i] < m
 * @param tlen    length of the target sequence
 * @param target  target sequence
 * @param m       number of residue types
 * @param mat     m*m scoring matrix in one-dimension array
 * @param gapo    gap open penalty; a gap of length l cost "-(gapo+l*gape)"
 * @param gape    gap extension penalty
 * @param xtra    extra information (see below)
 * @param qry     query profile (see below)
 *
 * @return        alignment information in a struct; unset values to -1
 *
 * When xtra==0, kswge_align() uses a signed two-byte integer to store a
 * score and only finds the best score and the end positions. The 2nd best
 * score or the start positions are not attempted. The default behavior can
 * be tuned by setting KSW_X* flags:
 *
 *   KSW_XBYTE:  use an unsigned byte to store a score. If overflow occurs,
 *               kswr_t::score will be set to 255
 *
 *   KSW_XSUBO:  track the 2nd best score and the ending position on the
 *               target if the 2nd best is higher than (xtra&0xffff)
 *
 *   KSW_XSTOP:  stop if the maximum score is above (xtra&0xffff)
 *
 *   KSW_XSTART: find the start positions
 *
 * When *qry==NULL, kswge_align() will compute and allocate the query profile
 * and when the function returns, *qry will point to the profile, which can
 * be deallocated simply by free(). If one query is aligned against multiple
 * target sequences, *qry should be set to NULL during the first call and
 * freed after the last call. Note that qry can equal 0. In this case, the
 * query profile will be deallocated in kswge_align().
 */
kswr_t kswge_align(int qlen, uint8_t *query, int tlen, uint8_t *target, int m, const int8_t *mat, int gapo, int gape, int xtra, kswq_t **qry);
kswr_t kswge_align2(int qlen, uint8_t *query, int tlen, uint8_t *target, int m, const int8_t *mat, int o_del, int e_del, int o_ins, int e_ins, int xtra, kswq_t **qry);
kswr_t* kswge_semi_global(int qlen, uint8_t *query, int tlen, uint8_t *target, int m, const int8_t *mat, int gapo, int gape);
kswr_t* kswge_semi_global2(int qlen, uint8_t *query, int tlen, uint8_t *target, int m, const int8_t *mat, int o_del, int e_del, int o_ins, int e_ins);

/**
 * Banded global alignment
 *
 * @param qlen    query length
 * @param query   query sequence with 0 <= query[i] < m
 * @param tlen    target length
 * @param target  target sequence with 0 <= target[i] < m
 * @param m       number of residue types
 * @param mat     m*m scoring mattrix in one-dimension array
 * @param gapo    gap open penalty; a gap of length l cost "-(gapo+l*gape)"
 * @param gape    gap extension penalty
 * @param w       band width
 * @param n_cigar (out) number of CIGAR elements
 * @param cigar   (out) BAM-encoded CIGAR; caller need to deallocate with free()
 *
 * @return        score of the alignment
 */
int kswge_global(int qlen, const uint8_t *query, int tlen, const uint8_t *target, int m, const int8_t *mat, int gapo, int gape, int w, int *n_cigar, uint32_t **cigar);
int kswge_global2(int qlen, const uint8_t *query, int tlen, const uint8_t *target, int m, const int8_t *mat, int o_del, int e_del, int o_ins, int e_ins, int w, int *n_cigar, uint32_t **cigar);

/**
 * Extend alignment
 *
 * The routine aligns $query and $target, assuming their upstream sequences,
 * which are not provided, have been aligned with score $h0. In return,
 * region [0,*qle) on the query and [0,*tle) on the target sequences are
 * aligned together. If *gscore>=0, *gscore keeps the best score such that
 * the entire query sequence is aligned; *gtle keeps the position on the
 * target where *gscore is achieved. Returning *gscore and *gtle helps the
 * caller to decide whether an end-to-end hit or a partial hit is preferred.
 *
 * The first 9 parameters are identical to those in kswge_global()
 *
 * @param h0      alignment score of upstream sequences
 * @param _qle    (out) length of the query in the alignment
 * @param _tle    (out) length of the target in the alignment
 * @param _gtle   (out) length of the target if query is fully aligned
 * @param _gscore (out) score of the best end-to-end alignment; negative if not found
 *
 * @return        best semi-local alignment score
 */
int kswge_extend(int qlen, const uint8_t *query, int tlen, const uint8_t *target, int m, const int8_t *mat, int gapo, int gape, int w, int end_bonus, int zdrop, int h0, int *qle, int *tle, int *gtle, int *gscore, int *max_off);
int kswge_extend2(int qlen, const uint8_t *query, int tlen, const uint8_t *target, int m, const int8_t *mat, int o_del, int e_del, int o_ins, int e_ins, int w, int end_bonus, int zdrop, int h0, int *qle, int *tle, int *gtle, int *gscore, int *max_off);

/** construct a scoring matrix for nucleotides
 *  //  A  C  G  T  N (or other ambiguous code)
 *  //  2 -2 -2 -2  0   A
 *  // -2  2 -2 -2  0   C
 *  // -2 -2  2 -2  0   G
 *  // -2 -2 -2  2  0   T
 *  //  0  0  0  0  0   N (or other ambiguous code)
 * @param match_score match score
 * @param mismatch_score mismatch score
 */
int8_t* kswge_gen_smat(const int match_score, const int mismatch_score);

/** calculate the number of mismatches and modify the cigar string, differentiate matches(=), mimatches(X), softclips(S), as well as update the tb and te in case some bad alignment occured
 * @param ret alignment result
 * @param ref pointer to the reference sequence
 * @param refLen length of reference sequence
 * @param read pointer to the read sequence
 * @param readLen length of the read
 */
void kswge_mark_mismatch (kswr_t* ret, uint8_t* ref, int refLen, uint8_t* read, int readLen);

/** Produce CIGAR 32-bit unsigned integer from CIGAR operation and CIGAR length
 * @param length length of CIGAR
 * @param op_letter CIGAR operation character ('M', 'I', etc)
 * @return 32-bit unsigned integer, representing encoded CIGAR operation and length
 */
static inline uint32_t kswge_gen_cigar (uint32_t length, char op_letter) {
        return (length << BAM_CIGAR_SHIFT) | (kswge_encoded_ops[(int)op_letter]);
}

/** Extract CIGAR operation character from CIGAR 32-bit unsigned integer
 * @param cigar_int 32-bit unsigned integer, representing encoded CIGAR operation and length
 * @return CIGAR operation character ('M', 'I', etc)
 */
static inline char kswge_cigar_opchr(uint32_t cigar_int) {
        return (cigar_int & 0xfU) > 8 ? 'M': MAPSTR[cigar_int & 0xfU];
}

/** Extract length of a CIGAR operation from CIGAR 32-bit unsigned integer
 * @param cigar_int 32-bit unsigned integer, representing encoded CIGAR operation and length
 * @return length of CIGAR operation
 */
static inline uint32_t kswge_cigar_oplen (uint32_t cigar_int) {
        return cigar_int >> BAM_CIGAR_SHIFT;
}

/** Output CIGAR string to stream
 * @param cigar cigars
 * @param len cigar length
 */
static inline void kswge_cigar_output(FILE* f, uint32_t* cigar, int32_t len, const char* title){
    fprintf(f, "%s", title);
    for(int32_t i = 0; i < len; ++i){
        fprintf(f, "%d%c", kswge_cigar_oplen(cigar[i]), kswge_cigar_opchr(cigar[i]));
    }
    fprintf(f, "\n");
}

/** convert a sequence of ACGTN into 01234
 * @param seq characters sequence
 * @return int8_t array rep of seq
 */
static inline uint8_t* kswge_seq2ints(const char* seq, const int len){
    uint8_t* num = (uint8_t*)malloc(len*sizeof(uint8_t));
    for(int i = 0; i < len; ++i) num[i] = kswge_nt2int[(int)seq[i]];
    return num;
}

/* output ksw alignment result */
void kswge_ret_output(FILE* f, const kswr_t* a, const char* rseq, const char* qseq, const unsigned char* table, const char* title, int wlen);

/** right align indel
 * @param ret align result
 * @param qints query int array
 * @param rints ref int array
 * @param beg beginning position to do realign
 * @param end ending position to do realign
 */
void kswge_right_align(kswr_t* ret, uint8_t* qints, uint8_t* rints, int beg=-1, int end=-1);

/** left align indel
 * @param ret align result
 * @param qints query int array
 * @param rints ref int array
 * @param beg beginning position to do realign
 * @param end ending position to do realign
 */
void kswge_left_align(kswr_t* ret, uint8_t* qints, uint8_t* rints, int beg=-1, int end=-1);

/** test whether two alignment compatible, they must be left/right align consistent
 * return -1 if not compatible
 * return 0 if alignment not overlap on ref
 * return 1 if alignment overlap and compatible
 */
int kswge_indel_compatible(kswr_t* aln1, kswr_t* aln2, const char* qseq1, const char* qseq2);

uint32_t* kswge_add_cigar(uint32_t* new_cigar, int32_t* p, int32_t* s, uint32_t length, char op);
uint32_t* store_previous_m(int8_t choice, uint32_t* length_m, uint32_t* length_x, int32_t* p, int32_t* s, uint32_t* new_cigar);

int kswge_main(int argc, char** argv);

inline void kswr_destroy(kswr_t* a){
    if(a->cigar) free(a->cigar);
    free(a);
}

#endif
