/* reconstructed from ksw authored by Attractive Chaos
 * the original source code can be found in https://github.com/attractivechaos/klib
 * I have changed some behaviors and add some other functionality for better usage with this project
 */
#ifndef AC_KREC_H
#define AC_KREC_H

#include <zlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>

#define KS1_SEP_SPACE 0 // isspace(): \t, \n, \v, \f, \r
#define KS1_SEP_TAB   1 // isspace() && !' '
#define KS1_SEP_LINE  2 // line separator: "\n" (Unix) or "\r\n" (Windows)
#define KS1_SEP_MAX   2
#define KS1_BUFFER_SIZE 16384

#define ks1_err(ks) ((ks)->end == -1)
#define ks1_eof(ks) ((ks)->is_eof && (ks)->begin >= (ks)->end)
#define ks1_rewind(ks) ((ks)->is_eof = (ks)->begin = (ks)->end = 0)
#define kseq1_rewind(ks) ((ks)->last_char = (ks)->f->is_eof = (ks)->f->begin = (ks)->f->end = 0)
#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

typedef struct {
    unsigned char *buf;
    int begin, end, is_eof;
    gzFile f;
} kstream1_t;

#ifndef KSTRANDS
#define KSTRANDS "-+"
#endif

#ifndef KSTRING_T
#define KSTRING_T kstring_t
typedef struct kstring_t {
        size_t l, m;
        char *s;
} kstring_t;
#endif

typedef struct {
    int last_char;
        kstream1_t *f;
} kseq1_t;             

typedef struct {
    kstring_t name, comment, seq, qual;
    uint32_t sample:18, off:8, match:1, strand:1, dwhy:4;
} krec1_t;

static inline kstream1_t *ks_init(gzFile f){
    kstream1_t *ks = (kstream1_t*)calloc(1, sizeof(kstream1_t));
    ks->f = f;
    ks->buf = (unsigned char*)malloc(KS1_BUFFER_SIZE);
    return ks;
}
static inline void ks1_destroy(kstream1_t *ks){
    if(ks){
        free(ks->buf); free(ks);
    }
}

static inline int ks1_getc(kstream1_t *ks)
{
    if(ks1_err(ks)) return -3;
    if(ks->is_eof && ks->begin >= ks->end) return -1;
    if(ks->begin >= ks->end){
       ks->begin = 0;
       ks->end = gzread(ks->f, ks->buf, KS1_BUFFER_SIZE);
       if(ks->end == 0) { ks->is_eof = 1; return -1;}
       if(ks->end == -1) { ks->is_eof = 1; return -3;}
    }
    return (int)ks->buf[ks->begin++];
}

static inline int ks1_getuntil2(kstream1_t *ks, int delimiter, kstring_t *str, int *dret, int append){
    int gotany = 0;
    if(dret) *dret = 0;
    str->l = append? str->l : 0;
    for(;;){
        int i;
        if(ks1_err(ks)) return -3;
        if(ks->begin >= ks->end){
            if(!ks->is_eof){
                ks->begin = 0;
                ks->end = gzread(ks->f, ks->buf, KS1_BUFFER_SIZE);
                if(ks->end == 0){ ks->is_eof = 1; break; }
                if(ks->end == -1){ ks->is_eof = 1; return -3; }
            }else break;
        }
        if(delimiter == KS1_SEP_LINE){ 
            for (i = ks->begin; i < ks->end; ++i) 
                if(ks->buf[i] == '\n') break; 
        } else if (delimiter > KS1_SEP_MAX) {
            for (i = ks->begin; i < ks->end; ++i)
                if (ks->buf[i] == delimiter) break;
        } else if (delimiter == KS1_SEP_SPACE) {
            for (i = ks->begin; i < ks->end; ++i)
                if (isspace(ks->buf[i])) break;
        } else if (delimiter == KS1_SEP_TAB) {
            for (i = ks->begin; i < ks->end; ++i)
                if (isspace(ks->buf[i]) && ks->buf[i] != ' ') break; 
        } else i = 0; /* never come to here! */
        if (str->m - str->l < (size_t)(i - ks->begin + 1)) {
            str->m = str->l + (i - ks->begin) + 1;
            kroundup32(str->m);
            str->s = (char*)realloc(str->s, str->m);
        }
        gotany = 1;
        memcpy(str->s + str->l, ks->buf + ks->begin, i - ks->begin); 
        str->l = str->l + (i - ks->begin);
        ks->begin = i + 1;
        if (i < ks->end) {
            if (dret) *dret = ks->buf[i];
            break;
        }
    }
    if (!gotany && ks1_eof(ks)) return -1;
    if (str->s == 0) {
        str->m = 1;
        str->s = (char*)calloc(1, 1);
    } else if (delimiter == KS1_SEP_LINE && str->l > 1 && str->s[str->l-1] == '\r') --str->l; 
    str->s[str->l] = '\0';
    return str->l;
}

static inline int ks1_getuntil(kstream1_t *ks, int delimiter, kstring_t *str, int *dret) { 
    return ks1_getuntil2(ks, delimiter, str, dret, 0);
}

static inline kseq1_t *kseq1_init(gzFile fd){
    kseq1_t *s = (kseq1_t*)calloc(1, sizeof(kseq1_t));
    s->f = ks_init(fd);
    return s;
}

static inline krec1_t *krec1_init(){
    krec1_t *r = (krec1_t*)calloc(1, sizeof(krec1_t)); 
    return r; 
}   
        
static inline void kseq1_destroy(kseq1_t *ks){
    if(!ks) return;
    ks1_destroy(ks->f);
    free(ks);
}

static inline void krec1_destroy(krec1_t *kr)                    
{      
    if(!kr) return;
    if(kr->name.s) {
        free(kr->name.s); 
        kr->name.s = NULL;
    }
    if(kr->comment.s){
        free(kr->comment.s); 
        kr->comment.s = NULL;
    }
    if(kr->seq.s){
        free(kr->seq.s); 
        kr->seq.s = NULL;
    }
    if(kr->qual.s){
        free(kr->qual.s); 
        kr->qual.s = NULL;
    }
    free(kr);
}

static inline void krec1_output(krec1_t *kr, FILE *f){
    if(kr->qual.l == 0){ // fa
        fprintf(f, ">%s\n", kr->name.s);
        fprintf(f, "%s\n", kr->seq.s);
    }else{ // fq
        fprintf(f, "@%s\n", kr->name.s);
        fprintf(f, "%s\n", kr->seq.s);
        fprintf(f, "%c\n", KSTRANDS[kr->strand]);
        fprintf(f, "%s\n", kr->qual.s);
    }
}

static inline void krec1_outgz(krec1_t *kr, gzFile fp){
    if(kr->qual.l == 0){ // fa
        gzprintf(fp, ">%s\n", kr->name.s);
        gzprintf(fp, "%s\n", kr->seq.s);
    }else{ // fq
        gzprintf(fp, "@%s\n", kr->name.s);
        gzprintf(fp, "%s\n", kr->seq.s);
        gzprintf(fp, "%c\n", KSTRANDS[kr->strand]);
        gzprintf(fp, "%s\n", kr->qual.s);
    }
}

/* Return value:
   >=0  length of the sequence (normal)
   -1   end-of-file
   -2   truncated quality string
   -3   error reading stream
 */
static inline int kseq1_read(kseq1_t *seq, krec1_t* rec){ 
    int c,r; 
    kstream1_t *ks = seq->f; 
    if (seq->last_char == 0) { /* then jump to the next header line */ 
        while ((c = ks1_getc(ks)) >= 0 && c != '>' && c != '@'); 
    if (c < 0) return c; /* end of file or error*/ 
        seq->last_char = c; 
    } /* else: the first header char has been read in the previous call */ 
    rec->comment.l = rec->seq.l = rec->qual.l = 0; /* reset all members */ 
    if ((r=ks1_getuntil(ks, 0, &rec->name, &c)) < 0) return r;  /* normal exit: EOF or error */ 
    if (c != '\n') ks1_getuntil(ks, KS1_SEP_LINE, &rec->comment, 0); /* read FASTA/Q comment */ 
    if (rec->seq.s == 0) { /* we can do this in the loop below, but that is slower */ 
        rec->seq.m = 256; 
        rec->seq.s = (char*)malloc(rec->seq.m); 
    } 
    while ((c = ks1_getc(ks)) >= 0 && c != '>' && c != '+' && c != '@' && c != '-') { 
        if (c == '\n') continue; /* skip empty lines */ 
        rec->seq.s[rec->seq.l++] = c; /* this is safe: we always have enough space for 1 char */ 
        ks1_getuntil2(ks, KS1_SEP_LINE, &rec->seq, 0, 1); /* read the rest of the line */ 
    } 
    if (c == '>' || c == '@') seq->last_char = c; /* the first header char has been read */
    if (rec->seq.l + 1 >= rec->seq.m) { /* seq->seq.s[seq->seq.l] below may be out of boundary */ 
        rec->seq.m = rec->seq.l + 2; 
        kroundup32(rec->seq.m); /* rounded to the next closest 2^k */ 
        rec->seq.s = (char*)realloc(rec->seq.s, rec->seq.m); 
    } 
    rec->seq.s[rec->seq.l] = 0;/* null terminated string */ 
    if (c == '+' || c == '-') rec->strand = c == '+' ? 1 : 0;
    else return rec->seq.l; /* FASTA */
    if (rec->qual.m < rec->seq.m) {/* allocate memory for qual in case insufficient */ 
        rec->qual.m = rec->seq.m; 
        rec->qual.s = (char*)realloc(rec->qual.s, rec->qual.m); 
    } 
    while ((c = ks1_getc(ks)) >= 0 && c != '\n'); /* skip the rest of strand line */ 
    if (c == -1) return -2; /* error: no quality string */ 
    while ((c = ks1_getuntil2(ks, KS1_SEP_LINE, &rec->qual, 0, 1) >= 0 && rec->qual.l < rec->seq.l)); 
    if (c == -3) return -3; /* stream error */ 
    seq->last_char = 0;/* we have not come to the next header line */ 
    if (rec->seq.l != rec->qual.l) return -2; /* error: qual string is of a different length */ 
    return rec->seq.l; 
}

#endif
