#include "preidx.h"
#include "htslib/sam.h"

pre_idx_t::pre_idx_t(){
    obwaa = new obwa_t();
}

pre_idx_t::~pre_idx_t(){
    if(obwaa){ delete obwaa; obwaa = NULL; }
    for(auto& obwa: obwas){
        if(obwa){ delete obwa; obwa = NULL;}
    }
}

void pre_idx_t::sets(obwa_t* obwa){
    obwa->setMismatchPenalty(1);
    obwa->setGapOpenPenalty(1);
    obwa->setGapExtendPenalty(1);
    obwa->set3PrimeClipPenalty(0);
    obwa->set5PrimeClipPenalty(0);
    obwa->setMatchScore(1); // after all penalty set
    obwa->fillScoreMatrix(); // after all score/penalty set
}

void pre_idx_t::init(const PrefixList& pfl, int mm, int mo, bool sem, bool dpre){
    maxm = mm;
    maxl = 0;
    sematch = sem;
    droppre = dpre;
    mins = pfl[0]->fseq.size();
    obwas.reserve(pfl.size()*2);
    sidx.reserve(pfl.size()*2);
    strands.reserve(pfl.size()*2);
    lens.reserve(pfl.size()*2);
    boff.reserve(pfl.size()*2);
    seidx.reserve(pfl.size()*2);
    std::vector<krec1_t*> recs;
    for(auto& p: pfl){
        // fwd
        krec1_t* kfwd = NULL;
        krec1_t* krev = NULL;
        if(p->fseq.size()){
            kfwd = krec1_init();
            ksprintf(&kfwd->name, "%d_fwd", p->sidx);
            ksprintf(&kfwd->seq, "%s", p->fseq.c_str());
        }
        // rev
        if(p->rseq.size()){
            krev = krec1_init();
            ksprintf(&krev->name, "%d_rev", p->sidx);
            ksprintf(&krev->seq, "%s", p->rseq.c_str());
        }
        // add to recs, sidx, strand, length, offset, seidx
        if(kfwd){
            recs.push_back(kfwd);
            sidx.push_back(p->sidx);
            strands.push_back(1);
            lens.push_back(p->fseq.size());
            boff.push_back(p->fbo);
            seidx.push_back(p->seidx);
        }
        if(krev){
            recs.push_back(krev);
            sidx.push_back(p->sidx);
            strands.push_back(0);
            lens.push_back(p->rseq.size());
            boff.push_back(p->rbo);
            seidx.push_back(p->seidx);
        }
        // construct partner obwa_t
        if(krev){
            obwa_t* obwar = new obwa_t();
            sets(obwar);
            obwar->buildIndex(krev);
            obwar->setMinSeedLength(MAX(MIN_SEED_LEN, (krev->seq.l - maxm)/(maxm+1)));
            obwar->setMinOutScore(obwar->getMinSeedLength());
            obwas.push_back(obwar);
        }
        if(kfwd){
            obwa_t* obwaf = new obwa_t();
            sets(obwaf);
            obwaf->buildIndex(kfwd);
            obwaf->setMinSeedLength(MAX(MIN_SEED_LEN, (kfwd->seq.l - maxm)/(maxm+1)));
            obwaf->setMinOutScore(obwaf->getMinSeedLength());
            obwas.push_back(obwaf);
        }
        // compute seed length
        if(p->fseq.size()){
            int sl = (p->fseq.size() - maxm)/(maxm+1);
            if(sl < mins) mins = sl;
        }
        if(p->rseq.size()){
            int sl = (p->rseq.size() - maxm)/(maxm+1);
            if(sl < mins) mins = sl;
        }
        // update maxpl
        if(p->fseq.size() > maxl) maxl = p->fseq.size();
        if(p->rseq.size() > maxl) maxl = p->rseq.size();
        // update minpl
        if(p->fseq.size() && p->fseq.size() < minl) minl = p->fseq.size();
        if(p->rseq.size() && p->rseq.size() < minl) minl = p->rseq.size();
    }
    mins = MAX(mins, MIN_SEED_LEN); // at least 6 bp seed length
    obwaa->setMinSeedLength(mins);
    obwaa->setMinOutScore(obwaa->getMinSeedLength());
    obwaa->buildIndex(recs);
    sets(obwaa);
    maxs = recs.size();
    for(auto& kr: recs){ krec1_destroy(kr); kr = NULL; }
    maxl += mo;
    maxo = mo + mm;
}

bool pre_idx_t::match(obwa_t* obwa, krec1_t* r, int& m, int& tid, int& scr){
    m = 0; tid = -1, scr = 0; r->dwhy = MATCH_FAIL_UNKNOWN;
    mem_alnreg_v ar = mem_align2(obwa->memopt, obwa->bwaidx->bwt, obwa->bwaidx->bns, obwa->bwaidx->pac, maxl, r->seq.s);
    if(!ar.n){
        if(ar.a) free(ar.a);
#ifdef PRI_MAIN_FUN
        ++nom;
#endif
        r->dwhy = MATCH_FAIL_ALN;
        return false;
    }
    int hs = MIN_SEED_LEN;
    int hc = 0, hi = 0, nm = 0, oo = 0, xm = 0, xo = 0;
    for(size_t i = 0; i < ar.n; ++i){
        if(ar.a[i].rid < 0 || ar.a[i].rid >= maxs) continue;
        if(ar.a[i].score < hs) break;
        int64_t pos = ar.a[i].rb < obwa->bwaidx->bns->l_pac ? ar.a[i].rb : ar.a[i].re - 1;
        if(pos >= obwa->bwaidx->bns->l_pac) continue; // skip rev
        pos -= obwa->bwaidx->bns->anns[ar.a[i].rid].offset;
        if(ar.a[i].qb > maxo || 
           abs(pos-ar.a[i].qb) > maxo ||
           abs((ar.a[i].qe - ar.a[i].qb)-abs(ar.a[i].re-ar.a[i].rb)) > maxo || 
           abs(abs(ar.a[i].re - ar.a[i].rb) - obwa->bwaidx->bns->anns[ar.a[i].rid].len) > maxo){
            if(ar.a[i].score > hs) xo = 1;
            continue;
        }
        if(ar.a[i].score > hs){
            mem_aln_t a_aln = mem_reg2aln(obwa->memopt, obwa->bwaidx->bns, obwa->bwaidx->pac, maxl, r->seq.s, &ar.a[i]);
            if(a_aln.NM < maxm){
                hs = a_aln.score;
                scr = hs;
                hc = 1;
                hi = i;
                nm = a_aln.NM;
                if(droppre) oo = ar.a[i].qe;
                else{
                    if((a_aln.cigar[0] & 0xf) == 3 || (a_aln.cigar[0] & 0xf) == 4){
                        oo = a_aln.cigar[0] >> 4;
                    }else oo = 0;
                }
            }else{
                xm = 1;
            }
            free(a_aln.cigar);
            if(a_aln.XA) free(a_aln.XA);
        }else if(ar.a[i].score == hs && hs > MIN_SEED_LEN){
            ++hc;
        }
    }
    if(hc == 0){
        if(xm) r->dwhy = MATCH_FAIL_MAXMM;
        else if(xo) r->dwhy = MATCH_FAIL_MAXOFF;
        else r->dwhy = MATCH_FAIL_SCORE;
    }else if(hc > 1){
        r->dwhy = MATCH_FAIL_MULTIPLE;
    }
    if(hc != 1){
        free(ar.a);
        return false;
    }else{ // one and only one is okay
        m = nm;
        tid = ar.a[hi].rid;
        if(droppre) r->off = oo;
        else r->off = oo + boff[ar.a[hi].rid];
        r->sample = sidx[ar.a[hi].rid];
        r->match = true;
        r->strand = strands[ar.a[hi].rid];
        free(ar.a);
        return true;
    }
}

bool pre_idx_t::dmatch(obwa_t* obwa, krec1_t* r, MatResutMap& mrs){
    mrs.clear();
    mem_alnreg_v ar = mem_align2(obwa->memopt, obwa->bwaidx->bwt, obwa->bwaidx->bns, obwa->bwaidx->pac, maxl, r->seq.s);
    if(!ar.n){
        if(ar.a) free(ar.a);
#ifdef PRI_MAIN_FUN
        ++nom;
#endif
        r->dwhy = MATCH_FAIL_ALN; 
        return false;
    }
    int hs = MIN_SEED_LEN, hi = 0, xo = 0, mx = 0;
    for(size_t i = 0; i < ar.n; ++i){
        if(ar.a[i].rid < 0 || ar.a[i].rid >= maxs) continue;
        if(ar.a[i].score < hs) break;
        int64_t pos = ar.a[i].rb < obwa->bwaidx->bns->l_pac ? ar.a[i].rb : ar.a[i].re - 1;
        if(pos >= obwa->bwaidx->bns->l_pac) continue; // skip rev
        pos -= obwa->bwaidx->bns->anns[ar.a[i].rid].offset;
        if(ar.a[i].qb > maxo || 
           abs(pos-ar.a[i].qb) > maxo ||
           abs((ar.a[i].qe - ar.a[i].qb)-abs(ar.a[i].re-ar.a[i].rb)) > maxo || 
           abs(abs(ar.a[i].re - ar.a[i].rb) - obwa->bwaidx->bns->anns[ar.a[i].rid].len) > maxo){
            if(ar.a[i].score > hs) xo = 1;
            continue;
        }
        if(ar.a[i].score > hs){
            mem_aln_t a_aln = mem_reg2aln(obwa->memopt, obwa->bwaidx->bns, obwa->bwaidx->pac, maxl, r->seq.s, &ar.a[i]);
            if(a_aln.NM < maxm){
                hs = a_aln.score;
                hi = i;
                mrs.clear();
                mm_mat_t mr;
                mr.mm = a_aln.NM;
                mr.scr = a_aln.score;
                if(droppre){
                    mr.off = ar.a[i].qe;
                }else{
                    if((a_aln.cigar[0] & 0xf) == 3 || (a_aln.cigar[0] & 0xf) == 4){
                        mr.off = a_aln.cigar[0] >> 4;
                    }else mr.off = 0;
                    mr.off += boff[ar.a[i].rid];
                }
                mrs[ar.a[hi].rid] = mr;
            }else{
                mx = 1;
            }
            free(a_aln.cigar);
            if(a_aln.XA) free(a_aln.XA);
        }else if(ar.a[i].score == hs && hs > MIN_SEED_LEN){
            mem_aln_t a_aln = mem_reg2aln(obwa->memopt, obwa->bwaidx->bns, obwa->bwaidx->pac, maxl, r->seq.s, &ar.a[i]);
            if(a_aln.NM < maxm){
                hi = i;
                mm_mat_t mr;
                mr.mm = a_aln.NM;
                mr.scr = a_aln.score;
                if(droppre){
                    mr.off = ar.a[i].qe;
                }else{
                    if((a_aln.cigar[0] & 0xf) == 3 || (a_aln.cigar[0] & 0xf) == 4){
                        mr.off = a_aln.cigar[0] >> 4;
                    }else mr.off = 0;
                    mr.off += boff[ar.a[i].rid];
                }
                mrs[ar.a[hi].rid] = mr;
            }
            free(a_aln.cigar);
            if(a_aln.XA) free(a_aln.XA);
        }
    }
    free(ar.a);
    if(mrs.size()){
        r->match = true;
        return true;
    }else{
        if(mx) r->dwhy = MATCH_FAIL_MAXMM;
        else if(xo) r->dwhy = MATCH_FAIL_MAXOFF;
        else r->dwhy = MATCH_FAIL_SCORE;
        return false;
    }
}

bool pre_idx_t::match(krec1_t* r1, krec1_t* r2, int& m){
    r1->off = r2->off = 0;
    r1->match = r2->match = false;
    r1->dwhy = r2->dwhy = MATCH_FAIL_UNKNOWN;
    m = 0;
    int tid1 = -1, tid2 = -1, scr = 0;
    if(sematch){
        if(match(obwaa, r1, m, tid1, scr)){
            r2->strand = !r1->strand;
            r2->sample = r1->sample;
            r2->match = true;
            return true;
        }else if(match(obwaa, r2, m, tid2, scr)){
            r1->strand = !r2->strand;
            r1->sample = r2->sample;
            r1->match = true;
            return true;
        }else{
            return false;
        }
    }else{
        int m1 = 0, m2 = 0;
        if(match(obwaa, r1, m1, tid1, scr) && (seidx[tid1] || match(obwas[tid1], r2, m2, tid2, scr))){
            r2->strand = !r1->strand;
            r2->sample = r1->sample;
            m = m1 + m2;
            r1->match = r2->match = true;
            return true;
        }
        if(match(obwaa, r2, m2, tid2, scr) && (seidx[tid2] || match(obwas[tid2], r1, m1, tid1, scr))){
            r1->strand = !r2->strand;
            r1->sample = r2->sample;
            m = m1 + m2;
            r1->match = r2->match = true;
            return true;
        }
        if(r1->match) r1->dwhy = MATCH_FAIL_PEFAIL;
        if(r2->match) r2->dwhy = MATCH_FAIL_PEFAIL;
        return false;
    }
}

bool pre_idx_t::dmatch(krec1_t* r1, krec1_t* r2, int& m){
    r1->off = r2->off = 0;
    r1->match = r2->match = false;
    r1->dwhy = r2->dwhy = MATCH_FAIL_UNKNOWN;
    m = 0;
    MatResutMap msr;
    int mm = 0, tid = -1, scr = 0, hc = 0, nx = 0, hs;
    if(dmatch(obwaa, r1, msr)){
        MatResutMap::iterator hi = msr.begin();
        hs = hi->second.scr;
        for(auto iter = msr.begin(); iter != msr.end(); ++iter){
            if(!seidx[iter->first]){
                if(match(obwas[iter->first], r2, mm, tid, scr)){
                    iter->second.mm += mm;
                    iter->second.scr += scr;
                    if(iter->second.scr > hs) hs = iter->second.scr, hi = iter, hc = 1;
                    else if(iter->second.scr == hs) hi = iter, hc = 2;
                }else{
                    iter->second.scr = -1;
                }
            }else{
                if(iter->second.scr > hi->second.scr) hi = iter, hc = 1;
                else if(iter->second.scr == hi->second.scr) hi = iter, hc = 2;
            }
        }
        if(hi->second.scr < 0) nx = 1;
        if(hi->second.scr > 0 && hc == 1){
            if(strands[hi->first]){
                r1->strand = 1;
                r2->strand = 0;
                r1->off = hi->second.off;
                r1->sample = r2->sample = sidx[hi->first];
            }else{
                r1->strand = 0;
                r2->strand = 1;
                r1->off = hi->second.off;
                r1->sample = r2->sample = sidx[hi->first];
            }
            m = hi->second.mm;
            r1->match = r2->match = true;
            return true;
        }
    }else if(dmatch(obwaa, r2, msr)){
        MatResutMap::iterator hi = msr.begin();
        hs = hi->second.scr;
        for(auto iter = msr.begin(); iter != msr.end(); ++iter){
            if(!seidx[iter->first]){
                if(match(obwas[iter->first], r1, mm, tid, scr)){
                    iter->second.mm += mm;
                    iter->second.scr += scr;
                    if(iter->second.scr > hs) hs = iter->second.scr, hi = iter, hc = 1;
                    else if(iter->second.scr == hs) hi = iter, hc = 2;
                }else{
                    iter->second.scr = -1;
                }
            }else{
                if(iter->second.scr > hi->second.scr) hi = iter, hc = 1;
                else if(iter->second.scr == hi->second.scr) hi = iter, hc = 2;
            }
        }
        if(hi->second.scr < 0) nx = 1;
        if(hi->second.scr > 0 && hc == 1){
            if(strands[hi->first]){
                r1->strand = 0;
                r2->strand = 1;
                r2->off = hi->second.off;
                r1->sample = r2->sample = sidx[hi->first];
            }else{
                r1->strand = 1;
                r2->strand = 0;
                r2->off = hi->second.off;
                r1->sample = r2->sample = sidx[hi->first];
            }
            m = hi->second.mm;
            r1->match = r2->match = true;
            return true;
        }
    }
    if(nx){
        if(r1->match) r1->dwhy = MATCH_FAIL_PEFAIL;
        if(r2->match) r2->dwhy = MATCH_FAIL_PEFAIL;
    }else if(hc == 2){
        if(r1->match) r1->dwhy = MATCH_FAIL_MULTIPLE;
        if(r2->match) r2->dwhy = MATCH_FAIL_MULTIPLE;
    }
    return false;
}

bool pre_idx_t::match(krec1_t* r, int& m){
    int tid = -1, src = 0;
    r->dwhy = MATCH_FAIL_UNKNOWN;
    return match(obwaa, r, m, tid, src);
}

void pri_usage(pre_idx_t* opt, char* arg0){
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage: %s [options]\n\n", arg0);
    fprintf(stderr, "Options: -i FILE input fastq r1\n");
    fprintf(stderr, "         -I FILE input fasta r2\n");
    fprintf(stderr, "         -c FILE primer configure file\n");
    fprintf(stderr, "         -m INT  max mismatch allowed for each match [%d]\n", opt->maxm);
    fprintf(stderr, "         -o INT  max leading offset allowed for each match [%d]\n", opt->maxo);
    fprintf(stderr, "         -s      single match is okay\n");
    fprintf(stderr, "\n");
}

int pri_main(int argc, char** argv){
    pre_idx_t opt;
    if(argc == 1){
        pri_usage(&opt, argv[0]);
        return 0;
    }
    int c = -1;
    while((c = getopt(argc, argv, "i:I:c:m:o:sh")) >= 0){
        switch(c){
            case 'i': opt.inr1 = optarg; break;
            case 'I': opt.inr2 = optarg; break;
            case 'c': opt.cfgf = optarg; break;
            case 'm': opt.maxm = atoi(optarg); break;
            case 'o': opt.maxo = atoi(optarg); break;
            case 's': opt.sematch = true; break;
            case 'h': pri_usage(&opt, argv[0]); return 0; break;
            default: break;
        }
    }
    if(opt.inr1 == NULL || opt.cfgf == NULL){
        fprintf(stderr, "input fq and configure file must be provided\n");
        return 1;
    }

    PrefixList prl;
    file2prl(opt.cfgf, prl);
    opt.init(prl, opt.maxm, opt.maxo, opt.sematch);
    fprintf(stderr, "==========spliter alignment stragety==========\n");
    fprintf(stderr, "      max mismatch allowed: %d\n", opt.maxm);
    fprintf(stderr, "max leading offset allowed: %d\n", opt.maxo);
    fprintf(stderr, "      single match allowed: %d\n", opt.sematch);
    fprintf(stderr, "      min alignment scores: %d\n", opt.obwaa->getMinOutScore());
    fprintf(stderr, " min alignment scores each:\n");
    for(auto& ob: opt.obwas){
        fprintf(stderr, "%d,\n", ob->getMinOutScore());
    }
    fprintf(stderr, "==========spliter alignment stragety==========\n");
    krec1_t* kr1 = krec1_init();
    krec1_t* kr2 = krec1_init();
    gzFile fd1 = gzopen(opt.inr1, "r");
    gzFile fd2 = NULL;
    kseq1_t* ks1 = kseq1_init(fd1);
    kseq1_t* ks2 = NULL;
    if(opt.inr2){
        fd2 =  gzopen(opt.inr2, "r");
        ks2 = kseq1_init(fd2);
    }
    int64_t ttr = 0;
    int64_t spr = 0;
    int m;
    while((kseq1_read(ks1, kr1) >= 0)){
        ++ttr;
        if(ks2){
            kseq1_read(ks2, kr2);
            if(opt.match(kr1, kr2, m)) ++spr;
        }else{
            if(opt.match(kr1, m)) ++spr;
        }
    }
    fprintf(stderr, "TotalReads: %lld\n", ttr);
    fprintf(stderr, "SplitReads: %lld\n", spr);
    fprintf(stderr, "NotAnySMat: %lld\n", opt.nom);
    gzclose(fd1);
    krec1_destroy(kr1); krec1_destroy(kr2);
    kseq1_destroy(ks1);
    if(opt.inr2){
        kseq1_destroy(ks2);
        gzclose(fd2);
    }
    return 0;
}

#ifdef PRI_MAIN_FUN
int main(int argc, char** argv){
    return pri_main(argc, argv);
}
#endif

