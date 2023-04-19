#include "spliter.h"

Spliter::Spliter(Options* opt){
    mOpt = opt;
    mAlnSpliter = opt->preidx;
    mAmpSpliter = opt->ampidx;
    mResult = new SplitResult(mOpt);
    mQcR1 = new QcStat(opt);
    mQcR2 = NULL;
    if(mOpt->in2.size()) mQcR2 = new QcStat(opt);
}

Spliter::~Spliter(){
    if(mResult){
        delete mResult; mResult = NULL;
    }
    if(mQcR1){
        delete mQcR1; mQcR1 = NULL;
    }
    if(mQcR2){
        delete mQcR2; mQcR2 = NULL;
    }
}

void Spliter::init(){
    mResult->init();
}

bool Spliter::exactMatch(krec1_t* r, const std::string& p){
    bool match = true;
    for(int off = 0; off <= mOpt->maxoff; ++off){
        match = true;
        for(size_t i = 0; i < p.size(); ++i){
            if(r->seq.s[off+i] != p[i]){
                match = false;
                break;
            }
        }
        if(match){
            r->off = off;
            return true;
        }
    }
    return false;
}

prefix_t* Spliter::exactMatch(krec1_t* r){
    if(r->seq.l < mOpt->minlen) return NULL;
    if(mOpt->useHash){
        for(int off = 0; off <= mOpt->maxoff; ++off){
            auto iter = mOpt->vblens.rbegin();
            uint64_t hk = 0;
            for(size_t i = 0; i < *iter; ++i){
                hk <<= 2;
                hk |= nuc_to_2bit[(int)r->seq.s[i+off]];
            }
            while(iter != mOpt->vblens.rend()){
                auto sit = mOpt->prh[*iter].find(hk);
                if(sit != mOpt->prh[*iter].end()){
                    if(sit->second.size() > 1) return NULL; // multiple other end
                    r->off = off + sit->second[0]->fbo;
                    r->sample = sit->second[0]->sidx;
                    if(mOpt->droppre) r->off += sit->second[0]->fseq.size() - sit->second[0]->fbo;
                    r->strand = sit->second[0]->fwd;
                    return sit->second[0];
                }else{
                    int lbl = *iter;
                    ++iter;
                    if(iter != mOpt->vblens.rend()){
                        hk >>= 2 * (lbl - *iter);
                    }
                }
            }
            
        }
    }else{
        for(int off = 0; off <= mOpt->maxoff; ++off){
            for(auto& bl: mOpt->vblens){
                if(r->seq.l <= bl + off) continue; // skip too short reads
                std::string rseq(r->seq.s + off, r->seq.s + bl);
                auto iter = mOpt->psh[bl].find(rseq);
                if(iter != mOpt->psh[bl].end()){
                    if(iter->second.size() > 1) return NULL; // multiple other end
                    r->off = off + iter->second[0]->fbo;
                    r->sample = iter->second[0]->sidx;
                    r->strand = iter->second[0]->fwd;
                    if(mOpt->droppre) r->off += iter->second[0]->fseq.size() - iter->second[0]->fbo;
                    return iter->second[0];
                }
            }
        }
    }
    return NULL;
}

prefix_t* Spliter::exactMatch(krec1_t* r1, krec1_t* r2){
    if(r1->seq.l < mOpt->minlen || r2->seq.l < mOpt->minlen) return NULL;
    if(mOpt->useHash){
        for(int off = 0; off <= mOpt->maxoff; ++off){
            auto iter = mOpt->vblens.rbegin();
            uint64_t hk = 0;
            for(size_t i = 0; i < *iter; ++i){
                hk <<= 2;
                hk |= nuc_to_2bit[(int)r1->seq.s[i+off]];
            }
            while(iter != mOpt->vblens.rend()){
                auto sit = mOpt->prh[*iter].find(hk);
                if(sit != mOpt->prh[*iter].end()){
                    for(auto& ppt : sit->second){
                        if(exactMatch(r2, ppt->rseq)){
                            r1->off = off + ppt->fbo;
                            r2->off += ppt->rbo;
                            r1->sample = r2->sample = ppt->sidx;
                            r1->strand = ppt->fwd;
                            r2->strand = !r1->strand;
                            if(mOpt->droppre){
                                r1->off += ppt->fseq.size() - ppt->fbo;
                                r2->off += ppt->rseq.size() - ppt->rbo;
                            }
                            return ppt;
                        }
                    }
                    return NULL;
                }else{
                    int lbl = *iter;
                    ++iter;
                    if(iter != mOpt->vblens.rend()){
                        hk >>= 2 * (lbl - *iter);
                    }
                }
            }
        }
    }else{
        for(int off = 0; off <= mOpt->maxoff; ++off){
            for(auto& bl: mOpt->vblens){
                if(r1->seq.l <= bl + off) continue; // skip too short reads
                std::string rseq(r1->seq.s + off, r1->seq.s + bl);
                auto iter = mOpt->psh[bl].find(rseq);
                if(iter != mOpt->psh[bl].end()){
                    for(auto& ppt : iter->second){
                        if(exactMatch(r2, ppt->rseq)){
                            r1->off = off + ppt->fbo;
                            r2->off += ppt->rbo;
                            r1->sample = r2->sample = ppt->sidx;
                            r1->strand = ppt->fwd;
                            r2->strand = !r1->strand;
                            if(mOpt->droppre){
                                r1->off += ppt->fseq.size() - ppt->fbo;
                                r2->off += ppt->rseq.size() - ppt->rbo;
                            }
                            return ppt;
                        }
                    }
                }
            }
        }
    }
    return NULL;
}

void Spliter::splitRead(krec1_t* r){
    r->dwhy = MATCH_FAIL_UNKNOWN;
    if(mOpt->skipfr1 >= r->seq.l){
        r->dwhy = MATCH_FAIL_QC;
        mResult->addDropRead(r, NULL, true);
        return;
    }
    if(!mQcR1->statRead(r)){
        r->dwhy = MATCH_FAIL_QC;
        mResult->addDropRead(r, NULL, true);
        return;
    }
    if(mOpt->skipfr1){
        r->seq.s += mOpt->skipfr1;
        r->seq.l -= mOpt->skipfr1;
        r->qual.s += mOpt->skipfr1;
        r->qual.l -= mOpt->skipfr1;
    }
    if(mAmpSpliter){
        if(mAmpSpliter->match(r)){
            r->match = true;
            mResult->addSplitRead(r, NULL, 0);
        }else{
            mResult->addDropRead(r, NULL, false);
        }
        return;
    }
    // mm disallowed
    r->off = 0;
    r->match = false;
    int m = 0;
    prefix_t* rp = exactMatch(r);
    if(rp){
        r->match = true;
        mResult->addSplitRead(r, NULL, 0);
    }else if(mOpt->maxmm && mAlnSpliter->match(r, m)){
        r->match = true;
        mResult->addSplitRead(r, NULL, m);
    }else{
        mResult->addDropRead(r, NULL, false);
    }
}

void Spliter::splitRead(krec1_t* r1, krec1_t* r2){
    r1->dwhy = r2->dwhy = MATCH_FAIL_UNKNOWN;
    if(mOpt->skipfr1 >= r1->seq.l || mOpt->skipfr2 >= r2->seq.l){
        r1->dwhy = r2->dwhy = MATCH_FAIL_QC;
        mResult->addDropRead(r1, r2, true);
        return;
    }
    if(mQcR1->statRead(r1)){
        if(!mQcR2->statRead(r2)){
            r1->dwhy = r2->dwhy = MATCH_FAIL_QC;
            mResult->addDropRead(r1, r2, true);
            return;
        }
    }else{
        mQcR2->statRead(r2);
        r1->dwhy = r2->dwhy = MATCH_FAIL_QC;
        mResult->addDropRead(r1, r2, true);
        return;
    }
    if(mOpt->skipfr1){
        r1->seq.s += mOpt->skipfr1;
        r1->seq.l -= mOpt->skipfr1;
        r1->qual.s += mOpt->skipfr1;
        r1->qual.l -= mOpt->skipfr1;
    }
    if(mOpt->skipfr2){
        r2->seq.s += mOpt->skipfr2;
        r2->seq.l -= mOpt->skipfr2;
        r2->qual.s += mOpt->skipfr2;
        r2->qual.l -= mOpt->skipfr2;
    }
    if(mAmpSpliter){
        if(mAmpSpliter->match(r1, r2)){
            mResult->addSplitRead(r1, r2, 0);
        }else{
            mResult->addDropRead(r1, r2, false);
        }
        return;
    }
    r1->off = r2->off = 0;
    r1->match = r2->match = false;
    int m;
    prefix_t* rp = NULL;
    if(mOpt->sematch){
        rp = exactMatch(r1);
        if(rp){
            m = 0; 
            r1->match = r2->match = true;
            r2->sample = r1->sample;
            r2->strand = !r1->strand;
            mResult->addSplitRead(r1, r2, m);
        }else{
            rp = exactMatch(r2);
            if(rp){
                m = 0;
                r1->match = r2->match = true;
                r1->sample = r2->sample;
                r1->strand = !r2->strand;
                mResult->addSplitRead(r1, r2, m);
            }else{
                if(mOpt->dmatch){
                    if(mOpt->maxmm && mAlnSpliter->dmatch(r1, r2, m)){
                        r1->match = r2->match = true;
                        mResult->addSplitRead(r1, r2, m);
                    }else{
                        mResult->addDropRead(r1, r2, false);
                    }
                }else{
                    if(mOpt->maxmm && mAlnSpliter->match(r1, r2, m)){
                        r1->match = r2->match = true;
                        mResult->addSplitRead(r1, r2, m);
                    }else{
                        mResult->addDropRead(r1, r2, false);
                    }
                }
            }
        }
    }else{
        rp = exactMatch(r1, r2);
        if(mOpt->maxmm == 0 && (!rp)){
            rp = exactMatch(r2, r1);
        }
        if(rp){
            m = 0;
            r1->match = r2->match = true;
            mResult->addSplitRead(r1, r2, m);
        }else{
            if(mOpt->dmatch){
                if(mOpt->maxmm && mAlnSpliter->dmatch(r1, r2, m)){
                    r1->match = r2->match = true;
                    mResult->addSplitRead(r1, r2, m); 
                }else{
                    mResult->addDropRead(r1, r2, false);
                }
            }else{
                if(mOpt->maxmm && mAlnSpliter->match(r1, r2, m)){
                    r1->match = r2->match = true;
                    mResult->addSplitRead(r1, r2, m);
                }else{
                    mResult->addDropRead(r1, r2, false);
                }
            }
        }
    }
}
