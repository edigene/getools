#include "ampidx.h"
#include "htslib/sam.h"

amp_idx_t::amp_idx_t(){
    obwa = new obwa_t();
}

amp_idx_t::~amp_idx_t(){
    if(obwa){ delete obwa; obwa = NULL; }
}

void amp_idx_t::sets(){
    obwa->setMismatchPenalty(1);
    obwa->setGapOpenPenalty(1);
    obwa->setGapExtendPenalty(1);
    obwa->set3PrimeClipPenalty(0);
    obwa->set5PrimeClipPenalty(0);
    obwa->setMatchScore(1); // after all penalty set
    obwa->fillScoreMatrix(); // after all score/penalty set
}

void amp_idx_t::init(const AmpliconList& ampl, bool sem){
    sematch = sem;
    sidx.reserve(ampl.size());
    std::vector<krec1_t*> recs;
    for(size_t i = 0; i < ampl.size(); ++i){
        krec1_t* kr = krec1_init();
        kputs(ampl[i]->aname.c_str(), &kr->name);
        kputs(ampl[i]->aseq.c_str(), &kr->seq);
        recs.push_back(kr);
        sidx.push_back(i);
    }
    obwa->setMinSeedLength(minseed);
    obwa->setMinOutScore(minscore);
    obwa->buildIndex(recs);
    sets();
    maxs = recs.size();
    for(auto& kr: recs){ krec1_destroy(kr); kr = NULL; }
}

bool amp_idx_t::match(krec1_t* r){
    r->dwhy = MATCH_FAIL_UNKNOWN;
    mem_alnreg_v ar = mem_align2(obwa->memopt, obwa->bwaidx->bwt, obwa->bwaidx->bns, obwa->bwaidx->pac, r->seq.l, r->seq.s);
    if(!ar.n){
        if(ar.a) free(ar.a);
        r->dwhy = MATCH_FAIL_ALN;
        return false;
    }
    int hs = 0, hi = 0;
    for(size_t i = 0; i < ar.n; ++i){
        if(ar.a[i].rid < 0 || ar.a[i].rid >= maxs) continue;
        if(ar.a[i].score < hs) break;
        if(ar.a[i].score > hs){
            hs = ar.a[i].score;
            hi = 1;
            r->off = 0;
            r->sample = sidx[ar.a[i].rid];
            r->match = true;
            r->strand = ar.a[i].rb >= obwa->bwaidx->bns->l_pac ? 0 : 1;
        }else if(ar.a[i].score == hs && hs){
            ++hi;
        }
    }
    free(ar.a);
    if(hi == 1) return true;
    else{
        if(hi > 1) r->dwhy = MATCH_FAIL_MULTIPLE;
        else r->dwhy = MATCH_FAIL_SCORE;
        return false;
    }
}

bool amp_idx_t::match(krec1_t* r1, krec1_t* r2){
    r1->off = r2->off = 0;
    r1->match = r2->match = false;
    r1->dwhy = r2->dwhy = MATCH_FAIL_UNKNOWN;
    if(sematch){
        if(match(r1)){
            r2->strand = !r1->strand;
            r2->sample = r1->sample;
            r2->match = true;
            return true;
        }else if(match(r2)){
            r1->strand = !r2->strand;
            r1->sample = r2->sample;
            r1->match = true;
            return true;
        }
        return false;
    }else{
        MatResutMap mr1;
        msecond(r1, mr1);
        if(mr1.empty()) return false;
        MatResutMap mr2;
        msecond(r2, mr2);
        if(mr2.empty()) return false;
        MatResutMap::iterator bi1, bi2;
        int32_t bs = 0;
        for(auto ps1 = mr1.begin(); ps1 != mr1.end(); ++ps1){
            auto ps2 = mr2.find(ps1->first);
            if(ps2 != mr2.end()){
                if(ps1->second.strand != ps2->second.strand){
                    if(ps1->second.scr + ps2->second.scr > bs){
                        bs = ps1->second.scr + ps2->second.scr;
                        bi1 = ps1;
                        bi2 = ps2;
                    }else if(ps1->second.scr + ps2->second.scr == bs){
                        bs = -1;
                    }
                }
            }
        }
        if(bs > 0){
            r1->off = r2->off = 0;
            r1->match = r2->match = true;
            r1->sample = bi1->first;
            r2->sample = bi2->first;
            r1->strand = bi1->second.strand;
            r2->strand = bi2->second.strand;
            return true;
        }
        r1->dwhy = r2->dwhy = MATCH_FAIL_PEFAIL;
    }
    return false;
}

void amp_idx_t::msecond(krec1_t* r, MatResutMap& mrs){
    mrs.clear();
    mem_alnreg_v ar = mem_align2(obwa->memopt, obwa->bwaidx->bwt, obwa->bwaidx->bns, obwa->bwaidx->pac, r->seq.l, r->seq.s);
    if(!ar.n){
        if(ar.a) free(ar.a);
        r->dwhy = MATCH_FAIL_ALN;
        return;
    }
    int hs = 0;
    for(size_t i = 0; i < ar.n; ++i){
        if(ar.a[i].rid < 0 || ar.a[i].rid >= maxs) continue;
        if(ar.a[i].score < hs) break;
        if(ar.a[i].score > hs){
            mrs.clear();
            mm_mat_t mm;
            mm.off = mm.mm = 0;
            mm.scr = ar.a[i].score;
            mm.strand = ar.a[i].rb >= obwa->bwaidx->bns->l_pac ? 0 : 1;
            mrs[sidx[ar.a[i].rid]] = mm;
            hs = ar.a[i].score;
        }else if(ar.a[i].score == hs && hs){
            mm_mat_t mm;
            mm.off = mm.mm = 0;
            mm.scr = ar.a[i].score;
            mm.strand = ar.a[i].rb >= obwa->bwaidx->bns->l_pac ? 0 : 1;
            mrs[sidx[ar.a[i].rid]] = mm;
        }
    }
    free(ar.a);
}
