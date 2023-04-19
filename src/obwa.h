#ifndef OBWA_H
#define OBWA_H

#include "bwalib/bwa.h"
#include "bwalib/bwt.h"
#include "bwalib/sais.h"
#include "bwalib/utils.h"
#include "bwalib/bwamem.h"
#include "krec.h"
#include "htslib/kstring.h"
#include <vector>

// nucleotide compl nucleotide complement table
const unsigned char seq_nt16_cmp[256] = {
    0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15,
   16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  31,
   32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  43,  44,  45,  46,  47,
   48,  49,  50,  51,  52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63,
   64, 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
  'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z',  91,  92,  93,  94,  95,
   64, 'T', 'v', 'G', 'h', 'e', 'f', 'C', 'd', 'i', 'j', 'm', 'l', 'k', 'N', 'o',
  'p', 'q', 'y', 's', 'A', 'a', 'b', 'w', 'x', 'r', 'z', 123, 124, 125, 126, 127
};

// online bwa 
class obwa_t{
    public:
        mem_opt_t* memopt; // align options
        bwaidx_t* bwaidx;   // bwa index
        
        // obwa_t buildor
        obwa_t(){
            bwaidx = 0;
            memopt = mem_opt_init();
            memopt->flag |= MEM_F_SOFTCLIP;
        }
        
        // obwa_t destructor
        ~obwa_t(){
            if(bwaidx){
                bwa_idx_destroy(bwaidx);
            }
            if(memopt){
                free(memopt);
            }
        }

    public:
        
        // pac to bwt, from bwtindex.c bwt_pac2bwt(...) in bwa
        // 'is' algorithm only here, not work for long genomes
        // 50000000 is the threshold used in bwa
        bwt_t* pac2bwt(const uint8_t* pac, int bwt_seq_lenr){
            bwt_t *bwt;
            ubyte_t *buf;
            int i;
            // initialization
            bwt = (bwt_t*)calloc(1, sizeof(bwt_t));
            bwt->seq_len = bwt_seq_lenr; //bwa_seq_len(fn_pac); //dummy
            bwt->bwt_size = (bwt->seq_len + 15) >> 4;
            // prepare sequence
            memset(bwt->L2, 0, 5 * 4);
            buf = (ubyte_t*)calloc(bwt->seq_len + 1, 1);
            for(i = 0; i < (int)bwt->seq_len; ++i){
                buf[i] = pac[i>>2] >> ((3 - (i&3)) << 1) & 3;
                ++bwt->L2[1+buf[i]];
            }
            for(i = 2; i <= 4; ++i) bwt->L2[i] += bwt->L2[i-1];
            // Burrows-Wheeler Transform
            bwt->primary = is_bwt(buf, bwt->seq_len);
            bwt->bwt = (uint32_t*)calloc(bwt->bwt_size, 4);
            for(i = 0; i < (int)bwt->seq_len; ++i) bwt->bwt[i>>4] |= buf[i] << ((15 - (i&15)) << 1);
            free(buf);
            return bwt;
        }
        
        // add an annotation sequence
        bntann1_t* addAnns(const krec1_t* rec, bntann1_t* ann, size_t offset){
            ann->offset = offset;
            ann->name = (char*)malloc(rec->name.l+1); // +1 for \0
            strncpy(ann->name, rec->name.s, rec->name.l+1);
            ann->anno = (char*)malloc(7);
            strcpy(ann->anno, "(null)\0");
            ann->len = rec->seq.l;
            ann->n_ambs = 0; // number of "holes"
            ann->gi = 0; // gi?
            ann->is_alt = 0;
            return ann;
        }

        // add an seq to pac, from bntseq.c add1(...) in bwa
        uint8_t* addSeq2pac(const krec1_t* rec, bntseq_t* bns, uint8_t* pac, int64_t *m_pac, int *m_seqs, int *m_holes, bntamb1_t **q){
            bntann1_t *p;
            int lasts;
            if(bns->n_seqs == *m_seqs){
                *m_seqs <<= 1;
                bns->anns = (bntann1_t*)realloc(bns->anns, *m_seqs * sizeof(bntann1_t));
            }
            p = bns->anns + bns->n_seqs;
            p->name = strndup(rec->name.s, rec->name.l);
            p->anno = (rec->comment.s && rec->comment.l > 0)? strndup(rec->comment.s, rec->comment.l) : strdup("(null)");
            p->gi = 0; p->len = rec->seq.l;
            p->offset = (bns->n_seqs == 0)? 0 : (p-1)->offset + (p-1)->len;
            p->n_ambs = 0;
            for(size_t i = lasts = 0; i < rec->seq.l; ++i){
                int c = nst_nt4_table[(int)rec->seq.s[i]];
                if(c >= 4){ // N
                    if(lasts == rec->seq.s[i]){ // contiguous N
                        ++(*q)->len;
                    }else{
                        if(bns->n_holes == *m_holes){
                            (*m_holes) <<= 1;
                            bns->ambs = (bntamb1_t*)realloc(bns->ambs, (*m_holes) * sizeof(bntamb1_t));
                        }
                        *q = bns->ambs + bns->n_holes;
                        (*q)->len = 1;
                        (*q)->offset = p->offset + i;
                        (*q)->amb = rec->seq.s[i];
                        ++p->n_ambs;
                        ++bns->n_holes;
                    }
                }
                lasts = rec->seq.s[i];
                { // fill buffer
                    if(c >= 4) c = lrand48()&3;
                    if(bns->l_pac == *m_pac){ // double the pac size
                        *m_pac <<= 1;
                        pac = (uint8_t*)realloc(pac, *m_pac/4);
                        memset(pac + bns->l_pac/4, 0, (*m_pac - bns->l_pac)/4);
                    }
                    _set_pac(pac, bns->l_pac, c);
                    ++bns->l_pac;
                }
            }
            ++bns->n_seqs;
            return pac;
        }

        // make pac, from bntseq.c bns_fasta2bntseq(...) in bwa
        uint8_t* makePac(const std::vector<krec1_t*>& v, bool for_only){
            bntseq_t * bns = (bntseq_t*)calloc(1, sizeof(bntseq_t));
            uint8_t *pac = 0;
            int32_t m_seqs, m_holes;
            int64_t m_pac, l;
            bntamb1_t *q;

            bns->seed = 11; // fixed seed for random generator
            m_seqs = m_holes = 8; m_pac = 0x10000;
            bns->anns = (bntann1_t*)calloc(m_seqs, sizeof(bntann1_t));
            bns->ambs = (bntamb1_t*)calloc(m_holes, sizeof(bntamb1_t));
            pac = (uint8_t*) calloc(m_pac/4, 1);
            q = bns->ambs;
            // move through the unaligned sequences
            for(auto& s: v) pac = addSeq2pac(s, bns, pac, &m_pac, &m_seqs, &m_holes, &q);
            if(!for_only){
                // add the reverse complemented sequence
                int64_t ll_pac = (bns->l_pac * 2 + 3) / 4 * 4;
                if(ll_pac > m_pac){
                    m_pac = ll_pac;
                    pac = (uint8_t*)realloc(pac, m_pac/4);
                }
                memset(pac + (bns->l_pac+3)/4, 0, (m_pac - (bns->l_pac+3)/4*4) / 4);
                for(l = bns->l_pac - 1; l >= 0; --l, ++bns->l_pac) _set_pac(pac, bns->l_pac, 3-_get_pac(pac, l));
            }
            bns_destroy(bns);
            return pac;
        }

        // pack to file
        void writePacToFile(const char* file){
            FILE* fp;
            kstring_t ps = {0, 0, 0};
            ksprintf(&ps, "%s.pac", file);
            fp = xopen(ps.s, "wb"); free(ps.s);
            ubyte_t ct;
            err_fwrite(bwaidx->pac, 1, (bwaidx->bns->l_pac>>2) + ((bwaidx->bns->l_pac&3) == 0? 0 : 1), fp);
            if(bwaidx->bns->l_pac % 4 == 0){
                ct = 0;
                err_fwrite(&ct, 1, 1, fp);
            }
            ct = bwaidx->bns->l_pac % 4;
            err_fwrite(&ct, 1, 1, fp);
            err_fflush(fp);
            err_fclose(fp);
        }

        void alignSeq2mem(const char* seq, int len, std::vector<mem_aln_t>& result){
            result.clear();
            if(!bwaidx) return;
            mem_alnreg_v ar = mem_align1(memopt, bwaidx->bwt, bwaidx->bns, bwaidx->pac, len, seq);
            for(size_t i = 0; i < ar.n; ++i){
                mem_aln_t a_aln = mem_reg2aln(memopt, bwaidx->bns, bwaidx->pac, len, seq, &ar.a[i]);
                if(a_aln.score >= memopt->T){
                    result.push_back(a_aln);
                }else{
                    free(a_aln.cigar);
                    free(a_aln.XA);
                }
            }
            if(ar.a) free(ar.a);
        }

        void alignSeq2set(const char* seq, int len, std::vector<mem_aln_t>& result, int maxoff=0xffff){
            result.clear();
            mem_alnreg_v ar = mem_align2(memopt, bwaidx->bwt, bwaidx->bns, bwaidx->pac, len, seq);
            for(size_t i = 0; i < ar.n; ++i){
                if(ar.a[i].rid < 0) continue;
                int64_t pos = ar.a[i].rb < bwaidx->bns->l_pac ? ar.a[i].rb : ar.a[i].re - 1;
                if(pos >= bwaidx->bns->l_pac) continue; // skip rev
                pos -= bwaidx->bns->anns[ar.a[i].rid].offset;
                if(abs(pos-ar.a[i].qb) > maxoff) continue;
                if(abs((ar.a[i].qe - ar.a[i].qb)-abs(ar.a[i].re-ar.a[i].rb)) > maxoff) continue;
                if(abs((ar.a[i].re - ar.a[i].rb)-bwaidx->bns->anns[ar.a[i].rid].len) > maxoff) continue;
                mem_aln_t a_aln = mem_reg2aln(memopt, bwaidx->bns, bwaidx->pac, len, seq, &ar.a[i]);
                if(a_aln.score >= memopt->T){
                    result.push_back(a_aln);
                }else{
                    free(a_aln.cigar);
                    free(a_aln.XA);
                }
            }
            if(ar.a) free(ar.a);
        }

        mem_aln_t besthit2mem(const char* seq, int len){
            mem_aln_t ret;
            ret.cigar = NULL;
            ret.XA = NULL;
            ret.flag |= 4; // BAM_FUNMAP
            mem_alnreg_v ar = mem_align1(memopt, bwaidx->bwt, bwaidx->bns, bwaidx->pac, len, seq);
            if(ar.n == 0) return ret;
            ret = mem_reg2aln(memopt, bwaidx->bns, bwaidx->pac, len, seq, &ar.a[0]);
            for(size_t i = 1; i < ar.n; ++i){
                mem_aln_t a_aln = mem_reg2aln(memopt, bwaidx->bns, bwaidx->pac, len, seq, &ar.a[i]);
                if(a_aln.score > ret.score){
                    free(ret.cigar);
                    free(ret.XA);
                    ret = a_aln;
                }else{
                    free(a_aln.cigar);
                    free(a_aln.XA);
                }
            }
            return ret;
        }

        static void mem_aln_destroy(mem_aln_t& aln){
            if(aln.cigar){ free(aln.cigar); aln.cigar = NULL; }
            if(aln.XA){ free(aln.XA); aln.XA = NULL; }
        }

        // build index from a list of seq
        void buildIndex(const std::vector<krec1_t*>& v){
            if(!v.size()) return;
            // check the integrity of the input data
            for(auto& s: v) assert(s->name.l && s->name.s && s->seq.l && s->seq.s);
            if(bwaidx){ bwa_idx_destroy(bwaidx); bwaidx = NULL; }
            // allocate memory for idx
            bwaidx = (bwaidx_t*)calloc(1, sizeof(bwaidx_t));;
            uint8_t* pac_for = makePac(v, true);
            uint8_t* pac = makePac(v, false);
            size_t tlen = 0;
            for(auto& s: v) tlen += s->seq.l;
            bwt_t *bwt = pac2bwt(pac, 2*tlen);
            free(pac);
            bwt_bwtupdate_core(bwt);
            // build sa from bwt and occ. adds it to bwt struct
            bwt_cal_sa(bwt, 32);
            bwt_gen_cnt_table(bwt);
            // make the bns
            bntseq_t * bns = (bntseq_t*) calloc(1, sizeof(bntseq_t));
            bns->l_pac = tlen;
            bns->n_seqs = v.size();
            bns->seed = 11;
            bns->n_holes = 0;
            // make the anns
            bns->anns = (bntann1_t*)calloc(v.size(), sizeof(bntann1_t));
            size_t offset = 0;
            for(size_t k = 0; k < v.size(); ++k){
                addAnns(v[k],  &bns->anns[k], offset);
                offset += v[k]->seq.l;
            }
            //ambs is "holes", like N bases
            bns->ambs = 0;
            // make the in-memory idx struct
            bwaidx->bwt = bwt;
            bwaidx->bns = bns;
            bwaidx->pac = pac_for;
        }

        void buildIndex(krec1_t* rec){
            std::vector<krec1_t*> recs;
            recs.push_back(rec);
            buildIndex(recs);
        }

        // load external bwt index
        void loadIndex(const char* file){
            bwaidx_t* newIndex = bwa_idx_load(file, BWA_IDX_ALL);
            if(!newIndex){
                fprintf(stderr, "error loading index!!!");
                exit(EXIT_FAILURE);
            }
            if(newIndex){
                bwa_idx_destroy(bwaidx);
            }
            bwaidx = newIndex;
        }
        
        // write index to file
        void writeIndex(const char* file){
            if(!bwaidx) return;
            kstring_t p = {0, 0, 0};
            ksprintf(&p, "%s.bwt", file);
            bwt_dump_bwt(p.s, bwaidx->bwt);
            p.l = 0;
            ksprintf(&p, "%s.sa", file);
            bwt_dump_sa(p.s, bwaidx->bwt);
            bns_dump(bwaidx->bns, file);
            writePacToFile(file);
        }

        // set the gap open penalty
        inline void setGapOpenPenalty(int gapOpenPenalty = 6){
            assert(gapOpenPenalty > 0);
            memopt->o_del = memopt->o_ins = gapOpenPenalty;
        }

        // set the gap extension penalty
        inline void setGapExtendPenalty(int gapExtPenalty = 1){
            assert(gapExtPenalty > 0);
            memopt->e_del = memopt->e_ins = gapExtPenalty;
        }

        // set mismatch penalty
        inline void setMismatchPenalty(int mismatchPenaly = 4){
            assert(mismatchPenaly > 0);
            memopt->b = mismatchPenaly;
        }

        // set the reseed trigger, look for internal seeds inside a seed longer than seedlength * reseed
        inline void setReseedTriger(float reseed = 1.5){
            assert(reseed > 0);
            memopt->split_factor = reseed;
        }

        // set SW alignment bandwidth
        inline void setBandWidth(int width = 100){
            assert(width > 0);
            memopt->w = width;
        }

        // set the SW alignment off-diagonal X-dropoff
        inline void setXDropoff(int dropOff = 100){
            assert(dropOff > 0);
            memopt->zdrop = dropOff;
        }

        // set the 3' clipping penalty
        inline void set3PrimeClipPenalty(int p3ClipPenalty = 5){
            assert(p3ClipPenalty >= 0);
            memopt->pen_clip3 = p3ClipPenalty;
        }

        // set the 5' clipping penalty
        inline void set5PrimeClipPenalty(int p5ClipPenalty = 5){
            assert(p5ClipPenalty >= 0);
            memopt->pen_clip5 = p5ClipPenalty;
        }

        // set seed length
        inline void setMinSeedLength(int seedLen = 19){
            assert(seedLen > 0);
            memopt->min_seed_len = seedLen;
        }

        // get seed length
        inline int getMinSeedLength(){
            return memopt->min_seed_len;
        }

        // set min out score
        inline void setMinOutScore(int minoutscore = 30){
            assert(minoutscore > 0);
            memopt->T = minoutscore;
        }

        // get min out score
        inline int getMinOutScore(){
            return memopt->T;
        }
        
        // set the match score, this should be set first as it will scale penalty options 
        inline void setMatchScore(int matchScore = 1){
            assert(matchScore > 0);
            memopt->b *= matchScore;
            memopt->T *= matchScore;
            memopt->o_del *= matchScore;
            memopt->o_ins *= matchScore;
            memopt->e_del *= matchScore;
            memopt->e_ins *= matchScore;
            memopt->zdrop *= matchScore;
            memopt->pen_clip3 *= matchScore;
            memopt->pen_clip5 *= matchScore;
            memopt->pen_unpaired *= matchScore;
            memopt->a = matchScore;
        }
        
        // fill match score, this should be done after all score set
        inline void fillScoreMatrix(){
            bwa_fill_scmat(memopt->a, memopt->b, memopt->mat);
        }
};

void obwa_usage(char* arg0);
int obwa_main(int argc, char** argv);

#endif
